# Riepilogo e mappe della taratura del rilascio lento (risultati/rilascio_lento.csv):
#   grafici/9_taratura_rilascio_lento.png e risultati/riepilogo_rilascio_lento.txt
# Uso, dalla cartella del progetto: julia --project=. calibrazione_scambio/grafici_rilascio_lento.jl
using Plots, Printf, Statistics
const D = @__DIR__
L = readlines(joinpath(D, "risultati", "rilascio_lento.csv")); H = split(L[1], ',')
R = [Dict(zip(H, parse.(Float64, split(l, ',')))) for l in L[2:end] if !isempty(strip(l))]
ch(r) = (r["T_lento_h"], r["V_soglia_m3"])
combo = unique(ch.(R)); TLs = sort(unique(first.(combo))); Vs = sort(unique(last.(combo)))
righe(c) = [r for r in R if ch(r) == c]
uno(c) = [r for r in righe(c) if r["onda"] == 1]                      # una riga per scenario per gli indici
nse(c) = mean(r["NSE"] for r in uno(c)); dco(c) = mean(r["dcolmo"] for r in uno(c)); liv(c) = mean(r["mae_livello"] for r in uno(c))
ar(c) = mean(abs(r["anticipo_rientro_min"]) for r in righe(c) if !isnan(r["anticipo_rientro_min"]))
ad(c) = mean(abs(r["anticipo_discesa_min"]) for r in righe(c) if !isnan(r["anticipo_discesa_min"]))
val(c, k, o, n) = (i = findfirst(r -> ch(r) == c && r["scenario"] == k && r["onda"] == o, R); i === nothing ? NaN : R[i][n])
io = IOBuffer()
println(io, "Prova del rilascio lento della parte sotto la soglia (T = 0,5 h e τ = 1 h fissi, tabelle S_50m_q2, fittizi 1-4).")
println(io, "Anticipi in minuti (positivo: il 1D è prima del 2D): rientro = portata in S4200 oltre 67 m³/s; discesa = portata in")
println(io, "S6000 sotto 50 m³/s dopo il colmo. Medie dei valori assoluti su tutte le onde (fittizio_4: due onde).\n")
@printf(io, "%8s %8s %8s %9s %8s %10s %10s   %s\n", "T_lento", "soglia", "NSE med", "|Δcolmo|", "livello", "|rientro|", "|discesa|", "rientro f1/f3/f4a/f4b; colmo f4 2ª onda")
for c in sort(combo; by = c -> -nse(c))
    @printf(io, "%8g %8g %8.3f %9.2f %8.3f %10.0f %10.0f   %s; %.1f\n", c..., nse(c), dco(c), liv(c), ar(c), ad(c),
            join([@sprintf("%.0f", val(c, k, o, "anticipo_rientro_min")) for (k, o) in ((1, 1), (3, 1), (4, 1), (4, 2))], "/"), val(c, 4, 2, "colmo1D"))
end
s = String(take!(io)); write(joinpath(D, "risultati", "riepilogo_rilascio_lento.txt"), s); print(s)
st = (framestyle = :box, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8)
function mappa(f, tit, alto, fmt)
    M = [any(==((t, v)), combo) ? f((t, v)) : NaN for v in Vs, t in TLs]
    p = heatmap(1:length(TLs), 1:length(Vs), M; c = alto ? cgrad(:blues) : cgrad(:blues; rev = true), title = tit, colorbar = false,
                xlabel = "tempo di rilascio lento T_lento [h]", ylabel = "soglia [migliaia di m³]",
                xticks = (1:length(TLs), [@sprintf("%g", t) for t in TLs]), yticks = (1:length(Vs), [@sprintf("%g", v/1000) for v in Vs]), st...)
    ok = filter(!isnan, M); lo, hi = minimum(ok), maximum(ok)
    jb = alto ? argmax(replace(M, NaN => -Inf)) : argmin(replace(M, NaN => Inf))
    for i in eachindex(Vs), j in eachindex(TLs)
        isnan(M[i, j]) && continue
        lum = hi > lo ? (M[i, j] - lo) / (hi - lo) : 0.5
        annotate!(p, j, i, text(Printf.format(Printf.Format(fmt), M[i, j]), 8, (alto ? lum > 0.55 : lum < 0.45) ? :white : :black))
    end
    scatter!(p, [jb[2]], [jb[1]]; marker = :rect, markersize = 20, markercolor = RGBA(0, 0, 0, 0), markerstrokecolor = "#eb6834", markerstrokewidth = 2.5, label = false)
end
savefig(plot(mappa(nse, "NSE medio della portata in S6000", true, "%.3f"), mappa(ar, "anticipo medio del rientro [min]", false, "%.0f"),
             mappa(ad, "anticipo medio della discesa [min]", false, "%.0f"); layout = (1, 3), size = (1650, 520), dpi = 120,
             left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, plot_title = "Rilascio lento della parte sotto la soglia (T = 0,5 h, τ = 1 h; riquadro: il migliore)",
             plot_titlefontsize = 12), joinpath(D, "grafici", "9_taratura_rilascio_lento.png"))
