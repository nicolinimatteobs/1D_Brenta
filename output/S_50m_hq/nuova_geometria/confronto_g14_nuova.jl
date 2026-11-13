# Grafico 14 con la geometria nuova (scalini 2D rifatti): confronto con quello vecchio, senza e con il ponte,
# e controllo con un riferimento coerente (scalini 2D nuovi alla stessa portata, negli stessi istanti).
# Fatto una tantum. Uso, dalla cartella del progetto:
#   julia --project=. confronto_g14_nuova.jl <cartella 1D senza ponte> <cartella 1D con il ponte>
using Plots, Printf
const R = pwd()
const OUT = joinpath(R, "output/S_50m_hq/nuova_geometria")
const SENZA = ARGS[1]; const CON = ARGS[2]
const GRIGIO = "#8a8f98"; const BLU = "#2a78d6"; const ARANCIO = "#eb6834"

csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
G = csvnum(joinpath(R, "input/geometria/Test_Brenta/S_50m/localslope_xcoord.csv"))
x = getindex.(G, 3); zb = getindex.(G, 4); n = length(x); xk = x ./ 1000
nomi = [@sprintf("S%04d", round(Int, xi/50)*50) for xi in x]

# --- grafico 14 già calcolato (file di testo): media sugli scenari di scarto quadratico medio e differenza media ---
function leggi_txt(f)
    rms = Float64[]; dm = Float64[]
    for l in eachline(f)
        occursin(r"^S\d{4}\s", l) || continue
        v = [parse.(Float64, split(s)[2:3]) for s in split(l, '|')[2:end]]
        push!(rms, sum(getindex.(v, 1))/length(v)); push!(dm, sum(getindex.(v, 2))/length(v))
    end
    return (; rms, dm)
end
A = leggi_txt(joinpath(R, "output/S_50m/grafici/errore_tirante_sezioni.txt"))                 # geometria vecchia, senza ponte
E = leggi_txt(joinpath(OUT, "errore_tirante_sezioni_senza_ponte.txt"))                        # geometria nuova, senza ponte
F = leggi_txt(joinpath(OUT, "errore_tirante_sezioni.txt"))                                    # geometria nuova, con il ponte

# --- riferimento coerente: livello degli scalini 2D nuovi alla portata 2D, negli stessi istanti del grafico 14 ---
tab(f) = (L = readlines(f); (strip.(split(L[1], ',')), reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])))
lin(X, Y, q) = (j = clamp(searchsortedlast(X, q), 1, length(X)-1); Y[j] + (Y[j+1]-Y[j])*(q - X[j])/(X[j+1]-X[j]))
SC = map(nomi) do nm
    h, M = tab(joinpath(R, "calcolo_tabelle_daBM/input/S_0050_6000/ns_50_scalini_$nm.csv"))
    q = M[:, findfirst(==("Discharge [m3/s]"), h)]; w = M[:, findfirst(==("Mean wse [m]"), h)]
    k = findall(>(0), q); p = sortperm(q[k]); (q[k][p], w[k][p])
end
function coerente(BASE)
    rms = zeros(n); dm = zeros(n)
    for ev in ("fittizio_1", "fittizio_2", "fittizio_3", "fittizio_4")
        D = joinpath(BASE, ev)
        mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
        t1 = [parse(Float64, strip(l)) for l in eachline(joinpath(D, "time.csv")) if !isempty(strip(l))]
        Q1 = mat("Q.csv"); E1 = mat("eta.csv"); Qs = 0.5 .* (Q1[:, 1:end-1] .+ Q1[:, 2:end])
        function leggi2D(f, pref, u)
            L = readlines(f); hdr = strip.(split(L[1], ','))
            M = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
            return M[:, 1], M[:, [findfirst(==("$(pref)_$(nm) [$u]"), hdr) for nm in nomi]]
        end
        tw, E2 = leggi2D(joinpath(R, "reference/ns_50_$(ev)_wse.csv"), "WSE", "m")
        _,  Q2 = leggi2D(joinpath(R, "reference/ns_50_$(ev)_discharge.csv"), "Q", "m3/s")
        k = findall(τ -> τ > 0 && t1[1] <= τ <= t1[end], tw); E2 = E2[k, :]; Q2 = Q2[k, :]
        it(M, τ) = (i = clamp(searchsortedlast(t1, τ), 1, length(t1) - 1); θ = (τ - t1[i])/(t1[i+1] - t1[i]); (1 - θ) .* M[i, :] .+ θ .* M[i+1, :])
        E1c = reduce(vcat, [permutedims(it(E1, τ)) for τ in tw[k]]); Q1c = reduce(vcat, [permutedims(it(Qs, τ)) for τ in tw[k]])
        m = (E2 .> 1.0) .& (Q2 .> 0.0) .& (abs.(Q1c .- Q2) .<= 0.01 .* Q2)      # stessa selezione del grafico 14
        for j in 1:n
            d = [E1c[i, j] - lin(SC[j][1], SC[j][2], Q2[i, j]) for i in findall(m[:, j])]
            rms[j] += 100sqrt(sum(abs2, d)/length(d))/4; dm[j] += 100sum(d)/length(d)/4
        end
    end
    return (; rms, dm)
end
Ec = coerente(SENZA); Fc = coerente(CON)

default(fontfamily = "sans-serif", titlefontsize = 12, guidefontsize = 11, tickfontsize = 8, legendfontsize = 9,
        linewidth = 2, markersize = 3, markerstrokewidth = 0, framestyle = :box, grid = true, gridalpha = 0.15,
        foreground_color_legend = nothing, background_color_legend = RGBA(1, 1, 1, 0.85))
tk = findall(nm -> parse(Int, nm[2:end]) % 250 == 0, nomi)
function pannello(titolo, serie; etichette = Int[])
    p = plot(; title = titolo, ylabel = "scarto quadratico medio [cm]", xticks = (xk[tk], nomi[tk]), xrotation = 45,
             xlims = (0, 6.1), ylims = (0, 13.5), legend = :topleft, left_margin = 8Plots.mm, bottom_margin = 6Plots.mm)
    vspan!(p, [2.575, 3.025]; color = :gray, alpha = 0.12, linewidth = 0, label = "tratto S2600-S3000")
    vline!(p, [3.025]; color = :gray30, linestyle = :dot, linewidth = 1.5, label = "ponte P17 (tra S3000 e S3050)")
    for (lab, v, c, ls, lw) in serie
        plot!(p, xk, v; label = lab, color = c, linestyle = ls, linewidth = lw, marker = :circle)
    end
    for j in etichette
        annotate!(p, xk[j], serie[end][2][j] + 0.35, text(nomi[j], 8, :black, :bottom))
    end
    return p
end
peggiori(v) = sort(partialsortperm(v, 1:8; rev = true))
p1 = pannello("Grafico 14, media sui quattro scenari. Riferimento: piene 2D dei fittizi (fatte con la mesh vecchia)",
              [("geometria vecchia, senza ponte (output/S_50m)", A.rms, GRIGIO, :solid, 1.5),
               ("geometria nuova, senza ponte", E.rms, BLU, :solid, 2), ("geometria nuova, con il ponte P17", F.rms, ARANCIO, :solid, 2)];
              etichette = peggiori(F.rms))
p2 = pannello("Stessi istanti, riferimento coerente con le tabelle: livello degli scalini 2D nuovi alla stessa portata",
              [("geometria nuova, senza ponte", Ec.rms, BLU, :solid, 2), ("geometria nuova, con il ponte P17", Fc.rms, ARANCIO, :solid, 2)];
              etichette = peggiori(Fc.rms))
xlabel!(p2, "sezione")
savefig(plot(p1, p2; layout = (2, 1), size = (1600, 1050), dpi = 130), joinpath(OUT, "14_confronto_geometria_nuova.png"))

open(joinpath(OUT, "confronto_grafico_14.txt"), "w") do io
    println(io, "Grafico 14 con la geometria nuova: scarto quadratico medio del tirante 1D - 2D a parità di portata (entro l'1%) [cm],")
    println(io, "media sui quattro scenari; tra parentesi la differenza media (negativo: 1D più basso)")
    println(io, "A  = geometria vecchia, senza ponte (output/S_50m);  E = geometria nuova, senza ponte;  F = geometria nuova, con il ponte P17")
    println(io, "     riferimento: piene 2D dei fittizi in reference/, fatte con la mesh vecchia")
    println(io, "Ec, Fc = come E e F, negli stessi istanti, ma il riferimento è il livello degli scalini 2D nuovi alla portata 2D")
    println(io)
    @printf(io, "%-6s %14s %14s %14s | %14s %14s\n", "sez", "A", "E", "F", "Ec", "Fc")
    for j in 1:n
        @printf(io, "%-6s", nomi[j])
        for d in (A, E, F); @printf(io, " %6.1f (%+5.1f)", d.rms[j], d.dm[j]); end
        print(io, " |"); for d in (Ec, Fc); @printf(io, " %6.1f (%+5.1f)", d.rms[j], d.dm[j]); end
        println(io)
    end
    println(io)
    for (a, b, nm) in ((0.05, 2.35, "S0050-S2350"), (2.4, 2.55, "S2400-S2550"), (2.6, 3.0, "S2600-S3000"), (3.05, 3.5, "S3050-S3500"), (3.55, 6.0, "S3550-S6000"), (0.05, 6.0, "tutte le sezioni"))
        kk = findall(xi -> a - 1e-6 <= xi <= b + 1e-6, xk)
        @printf(io, "media su %-17s", nm * ":")
        for (l, d) in (("A", A), ("E", E), ("F", F), ("Ec", Ec), ("Fc", Fc)); @printf(io, "  %-2s %4.1f cm", l, sum(d.rms[kk])/length(kk)); end
        println(io)
    end
end
print(read(joinpath(OUT, "confronto_grafico_14.txt"), String))
