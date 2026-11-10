# Confronto tra più riempimenti (condizioni iniziali e spin-up) dello stesso tratto, dai file IC_*.csv.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/confronto_spinup.jl <cartella dei grafici> <cartella 1> <cartella 2> ...
# es. julia --project=. strumenti/confronto_spinup.jl "output/vaia/prova_IC/confronti IC" output/vaia/prova_IC/Y_min_10 ...
# La geometria (per le ascisse) è quella della riga 2 di settings.txt. Non lancia nessuna simulazione.
using Plots, Printf
const R   = dirname(@__DIR__)
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const GEO = joinpath(R, SET[2])
length(ARGS) >= 3 || error("servono la cartella dei grafici e almeno due cartelle da confrontare")
const OUT  = joinpath(R, ARGS[1]); mkpath(OUT)
const CASI = [joinpath(R, a) for a in ARGS[2:end]]
const COL  = ["#2a78d6", "#eb6834", "#1baf7a", "#d55181", "#7a7a7a"]
const STY  = [:solid, :dash, :dot, :dashdot, :solid]
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

mat(d, f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(d, f)) if !isempty(strip(l))])
G = [parse.(Float64, split(l, ',')) for l in readlines(GEO)[2:end] if !isempty(strip(l))]
x = getindex.(G, 3) ./ 1000; imax = length(x)

D = map(CASI) do d
    t = vec(mat(d, "IC_time.csv")) ./ 3600; Q = mat(d, "IC_Q.csv"); H = mat(d, "IC_h.csv"); E = mat(d, "IC_eta.csv")
    n = length(t)
    dQ = vec(maximum(abs.(Q .- Q[n, :]'); dims = 2))
    (; nome = basename(d), t, Q, H, E, n, dQ)
end
# casi con risultati identici: stessa riga nella legenda
uguali(a, b) = size(a.H) == size(b.H) && a.H == b.H && a.Q == b.Q
etich = [begin
    gemelli = [D[j].nome for j in eachindex(D) if j != i && uguali(D[i], D[j])]
    isempty(gemelli) ? D[i].nome : D[i].nome * " (identico a " * join(gemelli, ", ") * ")"
end for i in eachindex(D)]

# 1) condizioni iniziali: tirante alla partenza (primo passo salvato) e a fine riempimento
p1 = plot(; title = "Tirante alla partenza (primo passo, t = 5 s)", ylabel = "tirante [m]", legend = :topleft, st...)
p2 = plot(; title = "Tirante a fine riempimento", ylabel = "tirante [m]", legend = :topleft, st...)
p3 = plot(; title = "Portata a fine riempimento", ylabel = "portata [m³/s]", xlabel = "distanza [km]", legend = :topleft, st...)
for (i, c) in enumerate(D)
    kw = (color = COL[i], linestyle = STY[i], linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0, label = etich[i])
    plot!(p1, x, c.H[1, :]; kw...); plot!(p2, x, c.H[end, :]; kw...)
    plot!(p3, x, 0.5 .* (c.Q[end, 1:end-1] .+ c.Q[end, 2:end]); kw...)
end
savefig(plot(p1, p2, p3; layout = (3, 1), size = (1000, 950), dpi = 130, left_margin = 6Plots.mm),
        joinpath(OUT, "1_condizioni_iniziali.png"))

# 2) spin-up: portata in uscita, tirante minimo nel tratto, distanza dallo stato finale
q1 = plot(; title = "Portata in uscita durante il riempimento", ylabel = "portata [m³/s]", legend = :topleft, st...)
q2 = plot(; title = "Tirante minimo nel tratto", ylabel = "tirante [m]", legend = :bottomright, st...)
q3 = plot(; title = "Distanza dallo stato finale: max |Q − Q_finale|", ylabel = "[m³/s]", xlabel = "tempo di riempimento [h]",
          yscale = :log10, legend = :bottomleft, st...)
for (i, c) in enumerate(D)
    kw = (color = COL[i], linestyle = STY[i], linewidth = 2, label = etich[i])
    plot!(q1, c.t, c.Q[:, end]; kw...)
    plot!(q2, c.t, vec(minimum(c.H; dims = 2)); kw...)
    m = c.dQ .>= 1e-5; plot!(q3, c.t[m], c.dQ[m]; kw...)
end
savefig(plot(q1, q2, q3; layout = (3, 1), size = (1000, 950), dpi = 130, left_margin = 6Plots.mm),
        joinpath(OUT, "2_spinup.png"))

# 3) differenze rispetto al primo caso, alla partenza e a fine riempimento
r1 = plot(; title = "Differenza di tirante rispetto a $(D[1].nome)", ylabel = "Δ tirante [m]", xlabel = "distanza [km]", legend = :bottom, st...)
hline!(r1, [0.0]; color = :black, linewidth = 0.8, label = false)
for (i, c) in enumerate(D[2:end])
    plot!(r1, x, c.H[1, :] .- D[1].H[1, :]; color = COL[i+1], linestyle = STY[i+1], linewidth = 2, markershape = :circle, markersize = 2.5,
          markerstrokewidth = 0, label = etich[i+1] * ", partenza")
    plot!(r1, x, c.H[end, :] .- D[1].H[end, :]; color = COL[i+1], linewidth = 1, alpha = 0.6, label = etich[i+1] * ", fine riempimento")
end
savefig(plot(r1; size = (1000, 450), dpi = 130, left_margin = 6Plots.mm, bottom_margin = 5Plots.mm), joinpath(OUT, "3_differenze.png"))

io = IOBuffer()
for (i, c) in enumerate(D)
    ok = [maximum(abs.(c.Q[k, 1:imax] .- c.Q[k-1, 1:imax])) < 1e-2 && abs(c.Q[k, end] - c.Q[c.n, end]) < 1e-3 * c.Q[c.n, end] for k in 2:c.n]
    ik = findlast(!, ok); kc = ik === nothing ? 2 : min(ik + 2, c.n)
    @printf(io, "%-10s tirante alla partenza da %.3f a %.3f m | minimo nel riempimento %.3f m | uscita minima %.2f m3/s | criterio soddisfatto a %.2f h | fine a %.2f h\n",
            c.nome, minimum(c.H[1, :]), maximum(c.H[1, :]), minimum(c.H), minimum(c.Q[:, end]), c.t[kc], c.t[end])
end
for c in D[2:end]
    @printf(io, "%s - %s: stato finale, scarto massimo di quota %.4f m\n", c.nome, D[1].nome, maximum(abs.(c.E[end, :] .- D[1].E[end, :])))
end
s = String(take!(io)); write(joinpath(OUT, "riepilogo.txt"), s); println(s)
