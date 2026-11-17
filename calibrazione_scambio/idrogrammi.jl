# Idrogrammi 2D e 1D in S3500 (prima del rientro), S4200 (dopo il rientro di S4150) e S6000 (uscita), fittizi 1-4:
# 2D, ingresso, 1D con lo sfioro sugli argini (out/sfioro) e una o più simulazioni con lo scambio dalla tabella.
# Uso, dalla cartella del progetto:
#   julia --project=. calibrazione_scambio/idrogrammi.jl <nome del file> <cartella in out/> <etichetta> [<cartella> <etichetta> ...]
#   es. julia --project=. calibrazione_scambio/idrogrammi.jl 6_idrogrammi_tarato.png tarato "1D, tabella ΔQ, soglia 150 000 m³, T = 0,5 h, τ = 1 h"
using Plots, Printf
const DIR = @__DIR__; const PROGETTO = normpath(joinpath(DIR, ".."))
const NOME = ARGS[1]
const RUNS = [(ARGS[j], ARGS[j+1]) for j in 2:2:length(ARGS)-1]
const COLORI = ["#2a78d6", "#7a7a7a", "#1baf7a", "#d55181"]
const SEZ = ["S3500", "S4200", "S6000"]
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 7)
csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(f) if !isempty(strip(l))])
nomesez(xi) = @sprintf("S%04d", round(Int, xi/50)*50)
x = getindex.(csvnum(joinpath(PROGETTO, "input/geometria/Test_Brenta/S_50m_q2/localslope_xcoord.csv")), 3)
isez = [findfirst(xi -> nomesez(xi) == s, x) for s in SEZ]
VAR = [("sfioro", "1D, sfioro sugli argini", "#eb6834"); [(r, e, COLORI[j]) for (j, (r, e)) in enumerate(RUNS)]]
pp = Any[]
for k in 1:4
    L = readlines(joinpath(PROGETTO, "reference/ns_50_fittizio_$(k)_discharge.csv")); hdr = strip.(split(L[1], ','))
    M = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
    for (s, i) in zip(SEZ, isez)
        p = plot(; title = "fittizio_$k, $s", xlabel = "tempo [h]", ylabel = "portata [m³/s]", legend = (k == 1 && s == "S6000") ? :topright : false, st...)
        for (j, (d, nv, col)) in enumerate(VAR)
            D = joinpath(DIR, "out", d, "fittizio_$k")
            t1 = vec(mat(joinpath(D, "time.csv"))); Q1 = mat(joinpath(D, "Q.csv"))
            j == 1 && plot!(p, t1 ./ 3600, Q1[:, 1]; color = :black, linestyle = :dot, linewidth = 1, label = "ingresso (S0050)")
            plot!(p, t1 ./ 3600, s == "S6000" ? Q1[:, end] : 0.5 .* (Q1[:, i] .+ Q1[:, i+1]); color = col, linewidth = 1.5, label = nv)
        end
        plot!(p, M[:, 1] ./ 3600, M[:, findfirst(==("Q_$s [m3/s]"), hdr)]; color = :black, linewidth = 2, label = "2D")
        push!(pp, p)
    end
end
mkpath(joinpath(DIR, "grafici"))
savefig(plot(pp...; layout = (4, 3), size = (1500, 1450), dpi = 120, left_margin = 6Plots.mm, bottom_margin = 4Plots.mm,
             plot_title = "Idrogrammi 2D e 1D prima del rientro (S3500), dopo il rientro di S4150 (S4200) e in uscita (S6000)", plot_titlefontsize = 12),
        joinpath(DIR, "grafici", NOME))
