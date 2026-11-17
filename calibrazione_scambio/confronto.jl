# Colmo lungo il tratto e portata in S6000, 2D e 1D, fittizi 1-4: 1D con lo sfioro sugli argini (out/sfioro) e una o
# più simulazioni con lo scambio dalla tabella. Stampa colmo, NSE in S6000 e volumi (anche in risultati/<nome>.txt).
# Uso, dalla cartella del progetto:
#   julia --project=. calibrazione_scambio/confronto.jl <nome del file> <cartella in out/> <etichetta> [<cartella> <etichetta> ...]
#   es. julia --project=. calibrazione_scambio/confronto.jl 8_colmo_tarato tarato "1D, tabella ΔQ, soglia 150 000 m³, T = 0,5 h, τ = 1 h"
include(joinpath(@__DIR__, "esegui.jl")); include(joinpath(@__DIR__, "indici.jl"))
using Plots
const DIR = @__DIR__
const NOME = ARGS[1]
const RUNS = [("sfioro", "1D, sfioro sugli argini", "#eb6834"); [(ARGS[j], ARGS[j+1], ["#2a78d6", "#7a7a7a", "#1baf7a"][j ÷ 2]) for j in 2:2:length(ARGS)-1]]
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)
io = IOBuffer()
println(io, "Confronto con il 2D (tabelle S_50m_q2, fittizi 1-4). Colmo in S6000 (fittizio_4: 2ª onda), NSE in S6000, errore medio del")
println(io, "colmo lungo il tratto, volumi usciti e rientrati per tratto (netti in ogni tratto tra due sezioni in tutto l'evento).\n")
@printf(io, "%-10s %-55s %8s %8s %7s %9s %11s %11s %11s %11s\n", "scenario", "simulazione", "colmo2D", "colmo1D", "NSE", "|Δcolmo|", "1D uscito", "1D rientr.", "2D uscito", "2D rientr.")
pp = Any[]
for k in 1:4
    r = REF[k]; j2 = findall(>=(ONDA[k]), r.t2)
    p = plot(; title = "fittizio_$k: colmo lungo il tratto" * (k == 4 ? " (2ª onda)" : ""), xlabel = "distanza [km]", ylabel = "portata al colmo [m³/s]", legend = :bottomleft, st...)
    plot!(p, X ./ 1000, vec(maximum(r.Q2[j2, :]; dims = 1)); color = :black, linewidth = 2, label = "2D")
    h = plot(; title = "fittizio_$k: portata in S6000", xlabel = "tempo [h]", ylabel = "portata [m³/s]", legend = :topright, st...)
    plot!(h, r.t2 ./ 3600, r.Q2[:, end]; color = :black, linewidth = 2, label = "2D")
    for (d, nv, col) in RUNS
        D = joinpath(DIR, "out", d, "fittizio_$k")
        t1 = vec(mat(joinpath(D, "time.csv"))); Q1 = mat(joinpath(D, "Q.csv"))
        Qs = 0.5 .* (Q1[:, 1:end-1] .+ Q1[:, 2:end]); Qs[:, end] .= Q1[:, end]
        plot!(p, X ./ 1000, vec(maximum(Qs[findall(>=(ONDA[k]), t1), :]; dims = 1)); color = col, linewidth = 1.6, label = nv)
        plot!(h, t1 ./ 3600, Q1[:, end]; color = col, linewidth = 1.6, label = nv)
        e = indici(D, k)
        @printf(io, "%-10s %-55s %8.2f %8.2f %7.3f %9.2f %11.0f %11.0f %11.0f %11.0f\n", "fittizio_$k", nv, e.colmo2, e.colmo1, e.nse, e.dcolmo, e.v1_out, e.v1_in, e.v2_out, e.v2_in)
    end
    push!(pp, p, h)
end
s = String(take!(io)); write(joinpath(DIR, "risultati", NOME * ".txt"), s); print(s)
savefig(plot(pp...; layout = (4, 2), size = (1300, 1500), dpi = 120, left_margin = 6Plots.mm, bottom_margin = 4Plots.mm,
             plot_title = "Colmo lungo il tratto e portata in uscita: 2D, sfioro sugli argini e scambio dalla tabella ΔQ", plot_titlefontsize = 12),
        joinpath(DIR, "grafici", NOME * ".png"))
