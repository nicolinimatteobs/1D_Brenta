# Confronto di due o più simulazioni 1D dei fittizi 1-4 tra loro e con il 2D: portata e tirante nel tempo in alcune
# sezioni, colmo e tirante massimo lungo il tratto, differenza di portata in S6000 tra le simulazioni, indici in un
# file di testo. Ogni simulazione è una cartella che contiene fittizio_1 … fittizio_4 (es. le uscite con due tempi
# di rilascio lento). Tabelle e geometria dalle righe 1-2 di settings.txt.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/confronto_simulazioni_2D.jl <cartella di uscita> <cartella 1> <etichetta 1> <cartella 2> <etichetta 2> ...
#   es. julia --project=. strumenti/confronto_simulazioni_2D.jl output/Volume_in_out/confronto \
#         output/Volume_in_out/svuotamento_24h "1D, rilascio lento 24 h" output/Volume_in_out/svuotamento_32h "1D, rilascio lento 32 h"
using Plots, Printf, Statistics
const R   = dirname(@__DIR__)
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const OUT = joinpath(R, ARGS[1]); mkpath(OUT)
const SIM = [(joinpath(R, ARGS[j]), ARGS[j+1]) for j in 2:2:length(ARGS)-1]
const COL = ["#eb6834", "#2a78d6", "#1baf7a", "#d55181"]
const SEZ = ["S1000", "S3500", "S4200", "S6000"]
const T0  = 3 * 3600.0                     # tiranti confrontati dalle 3 h (prima il 2D è ancora asciutto)
const ONDE = Dict(1 => [(0.0, 99.0)], 2 => [(0.0, 99.0)], 3 => [(0.0, 99.0)], 4 => [(0.0, 14.5), (14.5, 99.0)])
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(f) if !isempty(strip(l))])
nomesez(xi) = @sprintf("S%04d", round(Int, xi/50)*50)
interp(t, y, τ) = (m = clamp(searchsortedlast(t, τ), 1, length(t) - 1); θ = clamp((τ - t[m]) / (t[m+1] - t[m]), 0, 1); (1 - θ) * y[m] + θ * y[m+1])

G = csvnum(joinpath(R, SET[2])); X = getindex.(G, 3); ZB = getindex.(G, 4); xk = X ./ 1000
isez = [findfirst(xi -> nomesez(xi) == s, X) for s in SEZ]
function leggi2D(f, pref, um)
    L = readlines(f); hdr = strip.(split(L[1], ','))
    M = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
    return M[:, 1], M[:, [findfirst(==("$(pref)_$(nomesez(xi)) [$um]"), hdr) for xi in X]]
end

# 1D: portata nelle sezioni (media delle due interfacce), tirante, scambio con la piana
function leggi1D(D)
    t = vec(mat(joinpath(D, "time.csv"))); Q = mat(joinpath(D, "Q.csv")); H = mat(joinpath(D, "h.csv"))
    L = isfile(joinpath(D, "Q_lat.csv")) ? mat(joinpath(D, "Q_lat.csv")) : zeros(size(H))
    (; t, Qs = 0.5 .* (Q[:, 1:end-1] .+ Q[:, 2:end]), Qin = Q[:, 1], Qout = Q[:, end], H, L)
end
nse(o, s) = 1 - sum((s .- o).^2) / sum((o .- mean(o)).^2)
function discesa(t, q, t0, tf)        # primo istante [h] dopo il colmo dell'onda con portata sotto 50 m³/s
    jj = findall(j -> t0 <= t[j]/3600 < tf, eachindex(t)); isempty(jj) && return NaN
    p = jj[argmax(q[jj])]
    j = findfirst(j -> j > p && t[j]/3600 < tf && q[j] < 50, eachindex(t)); j === nothing ? NaN : t[j]/3600
end
function colmo(t, q, t0, tf)
    jj = findall(j -> t0 <= t[j]/3600 < tf, eachindex(t)); p = jj[argmax(q[jj])]; (q[p], t[p]/3600)
end

io = IOBuffer()
println(io, "Confronto tra simulazioni 1D e 2D, fittizi 1-4: ", join(["$(e) ($(relpath(d, R)))" for (d, e) in SIM], "; "))
println(io, "Portata in S6000: colmo (fittizio_4 per onda), NSE e RMSE agli istanti 2D con portata 2D > 0; anticipo della discesa =")
println(io, "istante 2D − istante 1D in cui la portata in S6000 scende sotto 50 m³/s dopo il colmo (positivo: il 1D scende prima).")
println(io, "Tirante: scarto 1D − 2D allo stesso istante, sezioni 2D bagnate, dalle 3 h (media, scarto quadratico medio, medio assoluto).")
println(io, "Volumi dal Q_lat.csv della simulazione (dopo il riempimento).\n")
pcolmo = Any[]; pdiff = Any[]
for k in 1:4
    t2, Q2 = leggi2D(joinpath(R, "reference/ns_50_fittizio_$(k)_discharge.csv"), "Q", "m3/s")
    _,  E2 = leggi2D(joinpath(R, "reference/ns_50_fittizio_$(k)_wse.csv"), "WSE", "m")
    bag = (Q2 .> 0) .& (E2 .> 0)                                   # 2D bagnato
    H2 = ifelse.(bag, E2 .- permutedims(ZB), NaN)
    S = [leggi1D(joinpath(d, "fittizio_$k")) for (d, _) in SIM]
    th2 = t2 ./ 3600

    # 1-2: portata e tirante nel tempo nelle sezioni scelte
    for (nome, cosa) in (("1_portata_nel_tempo", :Q), ("2_tirante_nel_tempo", :H))
        pp = Any[]
        for (j, (s, i)) in enumerate(zip(SEZ, isez))
            p = plot(; title = "$s (x = $(@sprintf("%.2f", X[i]/1000)) km)", xlabel = "tempo [h]",
                     ylabel = cosa == :Q ? "portata [m³/s]" : "tirante [m]", legend = j == 1 ? :outertop : false, legend_columns = -1, st...)
            if cosa == :Q
                plot!(p, S[1].t ./ 3600, S[1].Qin; color = :gray, linestyle = :dash, linewidth = 1.2, label = "ingresso")
                plot!(p, th2, Q2[:, i]; color = :black, linewidth = 2, label = "2D")
            else
                plot!(p, th2, H2[:, i]; color = :black, linewidth = 2, label = "2D")
            end
            for (n, (sim, (_, e))) in enumerate(zip(S, SIM))
                plot!(p, sim.t ./ 3600, cosa == :Q ? sim.Qs[:, i] : sim.H[:, i]; color = COL[n], linewidth = 1.5, label = e)
            end
            push!(pp, p)
        end
        savefig(plot(pp...; layout = (2, 2), size = (1500, 950), dpi = 130, left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
                     plot_title = (cosa == :Q ? "Portata" : "Tirante") * " nel tempo, fittizio_$k", plot_titlefontsize = 12),
                joinpath(OUT, "$(nome)_fittizio_$k.png"))
    end

    # 3: colmo e tirante massimo lungo il tratto
    a = plot(; title = "fittizio_$k: colmo di portata lungo il tratto", xlabel = "x [km]", ylabel = "portata [m³/s]",
             legend = k == 1 ? :outertop : false, legend_columns = -1, st...)
    plot!(a, xk, vec(maximum(Q2; dims = 1)); color = :black, linewidth = 2, label = "2D")
    b = plot(; title = "fittizio_$k: tirante massimo lungo il tratto", xlabel = "x [km]", ylabel = "tirante [m]", legend = false, st...)
    plot!(b, xk, [maximum(filter(!isnan, H2[:, i]); init = NaN) for i in eachindex(X)]; color = :black, linewidth = 2)
    for (n, (sim, (_, e))) in enumerate(zip(S, SIM))
        plot!(a, xk, vec(maximum(sim.Qs; dims = 1)); color = COL[n], linewidth = 1.5, label = e)
        plot!(b, xk, vec(maximum(sim.H; dims = 1)); color = COL[n], linewidth = 1.5)
    end
    push!(pcolmo, a, b)

    # 4: differenza di portata in S6000 tra le simulazioni e la prima, e tra ciascuna e il 2D
    d = plot(; title = "fittizio_$k: portata in S6000, 1D − 2D", xlabel = "tempo [h]", ylabel = "differenza [m³/s]",
             legend = k == 1 ? :outertop : false, legend_columns = -1, st...)
    hline!(d, [0.0]; color = :black, linewidth = 0.8, label = false)
    ok = Q2[:, end] .> 0
    for (n, (sim, (_, e))) in enumerate(zip(S, SIM))
        q1 = [interp(sim.t, sim.Qout, τ) for τ in t2[ok]]
        plot!(d, th2[ok], q1 .- Q2[ok, end]; color = COL[n], linewidth = 1.5, label = e)
    end
    push!(pdiff, d)

    # indici
    println(io, "fittizio_$k")
    for (n, (sim, (_, e))) in enumerate(zip(S, SIM))
        q1 = [interp(sim.t, sim.Qout, τ) for τ in t2[ok]]
        dh = Float64[]
        for (m, τ) in enumerate(t2), i in eachindex(X)
            (τ < T0 || !bag[m, i]) && continue
            push!(dh, interp(sim.t, view(sim.H, :, i), τ) - H2[m, i])
        end
        dt = diff([0.0; sim.t]); Vu = sum(dt .* vec(sum(max.(-sim.L, 0); dims = 2))); Vr = sum(dt .* vec(sum(max.(sim.L, 0); dims = 2)))
        @printf(io, "  %-28s NSE %.4f, RMSE %.2f m³/s; tirante 1D − 2D: media %+.3f m, scarto q. m. %.3f m, medio ass. %.3f m\n",
                e, nse(Q2[ok, end], q1), sqrt(mean((q1 .- Q2[ok, end]).^2)), mean(dh), sqrt(mean(dh.^2)), mean(abs.(dh)))
        for (o, (t0, tf)) in enumerate(ONDE[k])
            c1, tc1 = colmo(sim.t, sim.Qout, t0, tf); c2, tc2 = colmo(t2, Q2[:, end], t0, tf)
            @printf(io, "  %-28s %scolmo S6000 %.2f m³/s a %.2f h (2D %.2f a %.2f h); anticipo della discesa %.0f min\n", "",
                    length(ONDE[k]) > 1 ? "onda $o: " : "", c1, tc1, c2, tc2, 60 * (discesa(t2, Q2[:, end], t0, tf) - discesa(sim.t, sim.Qout, t0, tf)))
        end
        @printf(io, "  %-28s volumi: uscito verso la piana %.0f m³, rientrato %.0f m³, rimasto nella piana %.0f m³\n", "", Vu, Vr, Vu - Vr)
    end
    println(io)
end
savefig(plot(pcolmo...; layout = (4, 2), size = (1500, 1700), dpi = 120, left_margin = 6Plots.mm, bottom_margin = 4Plots.mm,
             plot_title = "Colmo di portata e tirante massimo lungo il tratto", plot_titlefontsize = 12),
        joinpath(OUT, "3_colmo_e_tirante_massimo_lungo_il_tratto.png"))
savefig(plot(pdiff...; layout = (4, 1), size = (1300, 1500), dpi = 120, left_margin = 6Plots.mm, bottom_margin = 4Plots.mm,
             plot_title = "Differenza di portata 1D − 2D in S6000 (istanti 2D con portata > 0)", plot_titlefontsize = 12),
        joinpath(OUT, "4_differenza_portata_S6000.png"))
s = String(take!(io)); write(joinpath(OUT, "confronto.txt"), s); print(s)
