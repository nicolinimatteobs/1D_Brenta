# Grafici e riepilogo della taratura del serbatoio della piana, da risultati/taratura.csv:
#   grafici/2_taratura_T.png        solo il tempo di svuotamento T (soglia e ritardo nulli), con i volumi scambiati
#   grafici/4_taratura_soglia.png   soglia e T senza ritardo (mappe delle medie sui fittizi 1-4)
#   grafici/6_taratura_ritardo.png  NSE medio per ogni ritardo τ (mappe soglia × T)
#   risultati/riepilogo_taratura.txt
# Uso, dalla cartella del progetto: julia --project=. calibrazione_scambio/grafici_taratura.jl
using Plots, Printf, Statistics
const DIR = @__DIR__; const OUT = joinpath(DIR, "grafici"); mkpath(OUT)
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)
const COL = ["#2a78d6", "#1baf7a", "#eb6834", "#d55181"]
L = readlines(joinpath(DIR, "risultati", "taratura.csv")); H = split(L[1], ',')
R = [Dict(zip(H, parse.(Float64, split(l, ',')))) for l in L[2:end] if !isempty(strip(l))]
chiave(r) = (r["tau_h"], r["V_soglia_m3"], r["T_h"])
val(c, k, n) = (i = findfirst(r -> chiave(r) == c && r["scenario"] == k, R); i === nothing ? NaN : R[i][n])
med(c, n) = mean(val(c, k, n) for k in 1:4)
combo = unique(chiave.(R))
io = IOBuffer()
println(io, "Taratura del serbatoio della piana (scambio 'tabella', tabelle S_50m_q2, fittizi 1-4, affluenti spenti, dt massimo 60 s).")
println(io, "L'acqua uscita arriva al serbatoio dopo il ritardo τ; rientra solo il volume oltre la soglia, al massimo (V − soglia)/T.")
println(io, "Indici: NSE della portata in S6000; |Δcolmo| errore medio del colmo lungo il tratto [m³/s]; livello: errore medio")
println(io, "assoluto del livello allo stesso istante [m] (sezioni 2D bagnate, dalle 3 h). Per il fittizio_4 il colmo è della 2ª onda.\n")

# --- 1) solo T ---
Ts = sort(unique([c[3] for c in combo if c[1] == 0 && c[2] == 0]))
cT(T) = (0.0, 0.0, T)
xp = log10.(Ts .+ 0.125); tk = (xp, [@sprintf("%g", T) for T in Ts])
pan(t, y) = plot(; title = t, xlabel = "tempo di svuotamento T [h]", ylabel = y, xticks = tk, st...)
a = pan("NSE della portata in S6000", "NSE [-]"); b = pan("Colmo in S6000 (fittizio_4: 2ª onda)", "portata [m³/s]")
d = pan("Errore medio del colmo lungo il tratto", "|colmo 1D − 2D| medio [m³/s]"); e = pan("Errore medio del livello (stesso istante)", "|livello 1D − 2D| medio [m]")
f = pan("Volume rientrato / volume uscito (per tratto)", "rientrato / uscito [-]")
for k in 1:4
    s(n) = [val(cT(T), k, n) for T in Ts]; mk = (color = COL[k], marker = :circle, markersize = 3, markerstrokewidth = 0)
    plot!(a, xp, s("NSE"); mk..., label = "fittizio_$k")
    plot!(b, xp, s("colmo1D"); mk..., label = "fittizio_$k 1D")
    hline!(b, [val(cT(Ts[1]), k, "colmo2D")]; color = COL[k], linestyle = :dash, label = k == 1 ? "2D (tratteggio)" : false)
    plot!(d, xp, s("dcolmo"); mk..., label = "fittizio_$k")
    plot!(e, xp, s("mae_livello"); mk..., label = "fittizio_$k")
    plot!(f, xp, s("V1_rientrato_tratti") ./ s("V1_uscito_tratti"); mk..., label = "fittizio_$k 1D")
    hline!(f, [val(cT(Ts[1]), k, "V2_rientrato_tratti") / val(cT(Ts[1]), k, "V2_uscito_tratti")]; color = COL[k], linestyle = :dash, label = k == 1 ? "2D (tratteggio)" : false)
end
mn = [med(cT(T), "NSE") for T in Ts]
plot!(a, xp, mn; color = :black, linewidth = 2.5, label = "media")
plot!(d, xp, [med(cT(T), "dcolmo") for T in Ts]; color = :black, linewidth = 2.5, label = "media")
plot!(e, xp, [med(cT(T), "mae_livello") for T in Ts]; color = :black, linewidth = 2.5, label = "media")
Tb = Ts[argmax(mn)]
g = plot(; title = @sprintf("Volumi scambiati per tratto: 2D e 1D con T = 0 e T = %g h", Tb), ylabel = "volume [migliaia di m³]",
         xticks = (1:4, ["fittizio_$k" for k in 1:4]), legend = :topleft, st...)
for (j, (lab, n1, n2, T, col)) in enumerate([("2D", "V2_uscito_tratti", "V2_rientrato_tratti", Ts[1], :black),
                                             ("1D, T = 0", "V1_uscito_tratti", "V1_rientrato_tratti", 0.0, "#7a7a7a"),
                                             (@sprintf("1D, T = %g h", Tb), "V1_uscito_tratti", "V1_rientrato_tratti", Tb, "#2a78d6")])
    o = (j - 2) * 0.27
    bar!(g, (1:4) .+ o, [val(cT(T), k, n1) for k in 1:4] ./ 1000; bar_width = 0.25, color = col, alpha = 0.35, linecolor = :white, label = lab * ": uscito")
    bar!(g, (1:4) .+ o, [val(cT(T), k, n2) for k in 1:4] ./ 1000; bar_width = 0.25, color = col, linecolor = :white, label = lab * ": rientrato")
end
savefig(plot(a, b, d, e, f, g; layout = (3, 2), size = (1300, 1250), dpi = 120, left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
             plot_title = "Taratura del solo tempo di svuotamento T (soglia e ritardo nulli)", plot_titlefontsize = 12),
        joinpath(OUT, "2_taratura_T.png"))
println(io, "== 1) Solo il tempo di svuotamento T (soglia 0, ritardo 0)")
@printf(io, "%6s %8s %8s %8s %8s %8s %10s %10s\n", "T [h]", "NSE f1", "NSE f2", "NSE f3", "NSE f4", "NSE med", "|Δcolmo|", "livello")
for T in Ts
    @printf(io, "%6g %8.3f %8.3f %8.3f %8.3f %8.3f %10.2f %10.3f\n", T, [val(cT(T), k, "NSE") for k in 1:4]..., med(cT(T), "NSE"), med(cT(T), "dcolmo"), med(cT(T), "mae_livello"))
end
println(io)

# --- mappe soglia × T per un ritardo ---
function mappa(τ, n, tit, alto; Vs = nothing, Ts = nothing)
    Vs = Vs === nothing ? sort(unique([c[2] for c in combo if c[1] == τ])) : Vs
    Ts = Ts === nothing ? sort(unique([c[3] for c in combo if c[1] == τ && c[2] > 0])) : Ts
    M = [med((τ, V, T), n) for V in Vs, T in Ts]
    p = heatmap(1:length(Ts), 1:length(Vs), M; c = alto ? cgrad(:blues) : cgrad(:blues; rev = true), title = tit, colorbar = false,
                xlabel = "tempo di svuotamento T [h]", ylabel = "soglia [migliaia di m³]",
                xticks = (1:length(Ts), [@sprintf("%g", T) for T in Ts]), yticks = (1:length(Vs), [@sprintf("%g", V/1000) for V in Vs]), st...)
    ok = filter(!isnan, M); lo, hi = minimum(ok), maximum(ok)
    jb = alto ? argmax(replace(M, NaN => -Inf)) : argmin(replace(M, NaN => Inf))
    for i in eachindex(Vs), j in eachindex(Ts)
        isnan(M[i, j]) && continue
        lum = hi > lo ? (M[i, j] - lo) / (hi - lo) : 0.5
        annotate!(p, j, i, text(n == "dcolmo" ? @sprintf("%.2f", M[i, j]) : @sprintf("%.3f", M[i, j]), 8, (alto ? lum > 0.55 : lum < 0.45) ? :white : :black))
    end
    scatter!(p, [jb[2]], [jb[1]]; marker = :rect, markersize = 20, markercolor = RGBA(0, 0, 0, 0), markerstrokecolor = "#eb6834", markerstrokewidth = 2.5, label = false)
    return p
end

# --- 2) soglia e T, senza ritardo ---
Vs0 = sort(unique([c[2] for c in combo if c[1] == 0])); Ts0 = sort(unique([c[3] for c in combo if c[1] == 0 && c[2] > 0]))
pp = [mappa(0.0, "NSE", "NSE medio della portata in S6000", true; Vs = Vs0, Ts = Ts0),
      mappa(0.0, "dcolmo", "errore medio del colmo lungo il tratto [m³/s]", false; Vs = Vs0, Ts = Ts0),
      mappa(0.0, "mae_livello", "errore medio del livello [m]", false; Vs = Vs0, Ts = Ts0)]
savefig(plot(pp...; layout = (1, 3), size = (1650, 520), dpi = 120, left_margin = 6Plots.mm, bottom_margin = 6Plots.mm,
             plot_title = "Taratura di soglia e tempo di svuotamento, senza ritardo (medie sui fittizi 1-4; riquadro: il migliore)", plot_titlefontsize = 12),
        joinpath(OUT, "4_taratura_soglia.png"))

# --- 3) ritardo ---
# solo i ritardi con la griglia completa di soglia e T (le prove singole restano nel riepilogo)
VsR = sort(unique([c[2] for c in combo if c[1] > 0])); TsR = sort(unique([c[3] for c in combo if c[1] > 0]))
completo(τ) = all((τ, V, T) in combo for V in VsR, T in TsR)
taus = filter(completo, sort(unique([c[1] for c in combo])))
pp = [mappa(τ, "NSE", @sprintf("ritardo τ = %g h", τ), true; Vs = VsR, Ts = TsR) for τ in taus]
nc = 3; nr = cld(length(pp), nc)
append!(pp, [plot(; framestyle = :none) for _ in 1:(nr * nc - length(pp))])
savefig(plot(pp...; layout = (nr, nc), size = (1500, 470 * nr), dpi = 120, left_margin = 8Plots.mm, bottom_margin = 6Plots.mm,
             plot_title = "Taratura del ritardo τ: NSE medio della portata in S6000 (fittizi 1-4; riquadro: il migliore per ogni τ)", plot_titlefontsize = 12),
        joinpath(OUT, "6_taratura_ritardo.png"))

println(io, "== 2) e 3) Tutte le combinazioni, ordinate per NSE medio (le prime 25)")
@printf(io, "%6s %9s %6s %8s %8s %8s %8s %8s %10s %10s\n", "τ [h]", "soglia", "T [h]", "NSE f1", "NSE f2", "NSE f3", "NSE f4", "NSE med", "|Δcolmo|", "livello")
ord = sort(combo; by = c -> -med(c, "NSE"))
for c in ord[1:min(25, end)]
    @printf(io, "%6g %9g %6g %8.3f %8.3f %8.3f %8.3f %8.3f %10.2f %10.3f\n", c..., [val(c, k, "NSE") for k in 1:4]..., med(c, "NSE"), med(c, "dcolmo"), med(c, "mae_livello"))
end
b2 = combo[argmin([med(c, "dcolmo") for c in combo])]; b3 = combo[argmin([med(c, "mae_livello") for c in combo])]
@printf(io, "\nMigliori: NSE con τ %g h, soglia %g m³, T %g h; colmo lungo il tratto con τ %g h, soglia %g m³, T %g h;\n", ord[1]..., b2...)
@printf(io, "livello con τ %g h, soglia %g m³, T %g h.\n\n", b3...)
println(io, "== Volumi per tratto [m³] con la combinazione migliore per NSE")
@printf(io, "%-10s %11s %11s %11s %11s %11s %9s %9s\n", "scenario", "entrato", "1D uscito", "1D rientr.", "2D uscito", "2D rientr.", "1D rie/us", "2D rie/us")
for k in 1:4
    c = ord[1]
    @printf(io, "%-10s %11.0f %11.0f %11.0f %11.0f %11.0f %9.2f %9.2f\n", "fittizio_$k", val(c, k, "V_entrato"), val(c, k, "V1_uscito_tratti"), val(c, k, "V1_rientrato_tratti"),
            val(c, k, "V2_uscito_tratti"), val(c, k, "V2_rientrato_tratti"), val(c, k, "V1_rientrato_tratti") / val(c, k, "V1_uscito_tratti"),
            val(c, k, "V2_rientrato_tratti") / val(c, k, "V2_uscito_tratti"))
end
s = String(take!(io)); write(joinpath(DIR, "risultati", "riepilogo_taratura.txt"), s); print(s)
