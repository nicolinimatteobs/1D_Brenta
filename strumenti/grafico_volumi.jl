# Conservazione del volume nella simulazione (dopo il riempimento): grafico 15_conservazione_volume.png e
# conservazione_volume.txt in <cartella degli output>/grafici.
# I volumi sono ricostruiti come nel modello (BILANCIO_DI_MASSA): entrato da monte, uscito a valle e scambio con
# la piana sommando portata·dt passo per passo (Q.csv, Q_lat.csv, time.csv); volume in alveo come somma delle aree
# delle tabelle (riga 1 di settings.txt) per la lunghezza delle celle, dal livello di eta.csv (all'inizio, l'ultima
# riga di IC_eta.csv). Senza Q_lat.csv (risultati vecchi) lo scambio è ricavato per differenza e l'errore non c'è.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/grafico_volumi.jl [cartella degli output]
using Plots, Printf
const R   = dirname(@__DIR__)                                    # cartella del progetto
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:15]]
const TAB = joinpath(R, SET[1]); const GEO = joinpath(R, SET[2])
const D   = joinpath(R, length(ARGS) >= 1 ? ARGS[1] : SET[5])
const OUT = joinpath(D, "grafici"); mkpath(OUT)
const BLU = "#2a78d6"; const ARA = "#eb6834"; const VER = "#1baf7a"; const ROS = "#d55181"; const GRI = "#7a7a7a"
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]

# tabelle: per ogni sezione quote, area e pendenza dell'area, come READ_VALUES del modello
T = csvnum(TAB); imax = length(csvnum(GEO))
righe = [[r for r in T if Int(r[1]) == k] for k in 1:imax]
function area(k, η)
    r = righe[k]; ys = getindex.(r, 2)
    z = max(η, ys[1]); j = clamp(searchsortedlast(ys, z), 1, length(ys))
    d = z - ys[j]
    if j == length(ys) && d > 0 && j > 1                          # oltre l'ultima quota: pendenza dell'ultimo tratto
        return r[j][3] + (r[j][3] - r[j-1][3]) / (ys[j] - ys[j-1]) * d
    end
    return r[j][3] + r[j][4] * d
end
x = getindex.(csvnum(GEO), 3); dxb = diff(x)
dxi = [k == 1 ? dxb[1] : k == imax ? dxb[end] : 0.5 * (dxb[k] + dxb[k-1]) for k in 1:imax]   # come SPATIAL_DOMAIN!
Valveo(e) = sum(area(k, e[k]) * dxi[k] for k in 1:imax)

t = vec(mat("time.csv")); Q = mat("Q.csv"); E = mat("eta.csv"); E0 = mat("IC_eta.csv")[end, :]
dt = diff([0.0; t]); th = t ./ 3600
V0 = Valveo(E0); Va = [Valveo(E[n, :]) for n in eachindex(t)] .- V0     # aumento del volume in alveo
Vm = cumsum(dt .* Q[:, 1]); Vv = cumsum(dt .* Q[:, end])
con_lat = isfile(joinpath(D, "Q_lat.csv"))
if con_lat
    L = mat("Q_lat.csv")
    Qu = vec(sum(max.(-L, 0.0); dims = 2)); Qr = vec(sum(max.(L, 0.0); dims = 2))
    Vu = cumsum(dt .* Qu); Vr = cumsum(dt .* Qr)
else
    Qu = zeros(length(t)); Qr = zeros(length(t))
    Vu = Vm .- Vv .- Va; Vr = zeros(length(t))                    # uscito per differenza: errore nullo per costruzione
end
Vp  = Vu .- Vr                                                   # rimasto nella piana (uscito − rientrato)
err = Va .- (Vm .- Vv .- Vp)                                     # errore del bilancio
k3(v) = v ./ 1000

a = plot(; title = "Volumi cumulati dall'inizio della simulazione", xlabel = "tempo [h]", ylabel = "volume [migliaia di m³]", legend = :outerright, st...)
plot!(a, th, k3(Vm); color = :black, linewidth = 2, label = "entrato da monte")
plot!(a, th, k3(Vv); color = BLU, linewidth = 2, label = "uscito a valle")
plot!(a, th, k3(Vu); color = ARA, linewidth = 1.5, label = "uscito verso la piana")
con_lat && plot!(a, th, k3(Vr); color = VER, linewidth = 1.5, label = "rientrato dalla piana")
plot!(a, th, k3(Va); color = ROS, linewidth = 1.5, label = "aumento in alveo")

b = plot(; title = "Entrato = uscito a valle + rimasto nella piana + aumento in alveo", xlabel = "tempo [h]",
         ylabel = "volume [migliaia di m³]", legend = :outerright, st...)
s1 = k3(Vv); s2 = s1 .+ k3(Vp); s3 = s2 .+ k3(Va)
plot!(b, th, s1; fillrange = 0, color = BLU, fillalpha = 0.35, linewidth = 0, label = "uscito a valle")
plot!(b, th, s2; fillrange = s1, color = ARA, fillalpha = 0.35, linewidth = 0, label = "rimasto nella piana")
plot!(b, th, s3; fillrange = s2, color = ROS, fillalpha = 0.6, linewidth = 0, label = "aumento in alveo")
plot!(b, th, k3(Vm); color = :black, linewidth = 2, linestyle = :dash, label = "entrato da monte")

c = plot(; title = @sprintf("Errore del bilancio (a fine simulazione %.2f m³, %.5f%% dell'entrato)", err[end], 100err[end] / Vm[end]),
         xlabel = "tempo [h]", ylabel = "errore [m³]", legend = false, st...)
hline!(c, [0.0]; color = :black, linewidth = 0.8)
plot!(c, th, err; color = :black, linewidth = 1.5)
con_lat || annotate!(c, th[end] / 2, 0.0, text("senza Q_lat.csv lo scambio è ricavato per differenza: errore nullo per costruzione", 8, :center, :bottom, GRI))

d = plot(; title = "Scambio con la piana", xlabel = "tempo [h]", ylabel = "portata [m³/s]", legend = :outerright, st...)
hline!(d, [0.0]; color = :black, linewidth = 0.8, label = false)
plot!(d, th, -Qu; color = ARA, linewidth = 1.5, label = "uscita (−)")
con_lat && plot!(d, th, Qr; color = VER, linewidth = 1.5, label = "rientro (+)")

savefig(plot(a, b, c, d; layout = (2, 2), size = (1500, 900), dpi = 130, left_margin = 6Plots.mm, bottom_margin = 5Plots.mm,
             plot_title = "Conservazione del volume, " * basename(D), plot_titlefontsize = 12),
        joinpath(OUT, "15_conservazione_volume.png"))

io = IOBuffer()
println(io, "Conservazione del volume nella simulazione (dopo il riempimento), $(basename(D)), volumi in m3")
@printf(io, "%-36s %14.1f\n", "volume in alveo all'inizio", V0)
@printf(io, "%-36s %14.1f\n", "volume in alveo alla fine", V0 + Va[end])
@printf(io, "%-36s %14.1f\n", "entrato da monte", Vm[end])
@printf(io, "%-36s %14.1f  (%.2f%%)\n", "uscito a valle", Vv[end], 100Vv[end] / Vm[end])
@printf(io, "%-36s %14.1f  (%.2f%%)\n", "uscito verso la piana", Vu[end], 100Vu[end] / Vm[end])
@printf(io, "%-36s %14.1f  (%.2f%% dell'uscito verso la piana)\n", "rientrato dalla piana", Vr[end], Vu[end] > 0 ? 100Vr[end] / Vu[end] : 0.0)
@printf(io, "%-36s %14.1f  (%.2f%%)\n", "rimasto nella piana", Vp[end], 100Vp[end] / Vm[end])
@printf(io, "%-36s %14.1f\n", "aumento in alveo", Va[end])
@printf(io, "%-36s %14.2f  (%.5f%% dell'entrato; massimo nel tempo %.2f)\n", "errore del bilancio", err[end], 100err[end] / Vm[end], maximum(abs.(err)))
con_lat || println(io, "Q_lat.csv assente: lo scambio è ricavato per differenza.")
s = String(take!(io)); write(joinpath(OUT, "conservazione_volume.txt"), s); print(s)
