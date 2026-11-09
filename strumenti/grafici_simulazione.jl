# Grafici di una simulazione. Legge da settings.txt tabelle, geometria, idrogramma di monte e
# cartella degli output, e salva i grafici in <cartella degli output>/grafici.
# Uso, dalla cartella del progetto:  julia --project=. strumenti/grafici_simulazione.jl
using Plots, Printf
const R     = dirname(@__DIR__)                                  # cartella del progetto
const SET   = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const TAB   = joinpath(R, SET[1])                                # tabelle delle sezioni (riga 1)
const GEO   = joinpath(R, SET[2])                                # geometria (riga 2)
const MONTE = joinpath(R, SET[14])                               # idrogramma di monte (riga 14)
const D     = joinpath(R, SET[22])                               # cartella degli output (riga 22)
const OUT   = joinpath(D, "grafici"); mkpath(OUT)
const BLU = "#2a78d6"; const ARA = "#eb6834"; const VER = "#1baf7a"; const ROS = "#d55181"; const GRI = "#7a7a7a"
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
vec1(f) = [parse(Float64, strip(l)) for l in eachline(joinpath(D, f)) if !isempty(strip(l))]
csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
G = csvnum(GEO)
x = getindex.(G, 3) ./ 1000; zb = getindex.(G, 4); imax = length(x)
T = csvnum(TAB)
ztab = [maximum(r[2] for r in T if r[1] == k) for k in 1:imax]      # quota massima di ogni tabella
P = [parse.(Float64, split(l, ',')) for l in eachline(MONTE) if !isempty(strip(l))]

t  = vec1("time.csv") ./ 3600;  Q = mat("Q.csv");  E = mat("eta.csv");  H = mat("h.csv")
ti = vec1("IC_time.csv") ./ 3600; QI = mat("IC_Q.csv"); EI = mat("IC_eta.csv"); HI = mat("IC_h.csv")
Qs = 0.5 .* (Q[:, 1:end-1] .+ Q[:, 2:end])            # portata nelle sezioni (media delle interfacce)
iout = size(Q, 2)

# 1) condizioni iniziali: stato di partenza e di fine riempimento
p1 = plot(; title = "Condizioni iniziali: quota del pelo libero", ylabel = "quota [m s.l.m.]", legend = :topright, st...)
plot!(p1, x, zb; label = "talweg", color = GRI, linewidth = 1.5, fillrange = minimum(zb) - 1, fillcolor = GRI, fillalpha = 0.15)
plot!(p1, x, EI[1, :]; label = "partenza (moto uniforme)", color = ARA, linewidth = 1.5, linestyle = :dash)
plot!(p1, x, EI[end, :]; label = @sprintf("fine riempimento (t = %.1f h)", ti[end]), color = BLU, linewidth = 2)
p2 = plot(; title = "Tirante", ylabel = "tirante [m]", legend = :topright, st...)
plot!(p2, x, HI[1, :]; label = "partenza", color = ARA, linewidth = 1.5, linestyle = :dash, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
plot!(p2, x, HI[end, :]; label = "fine riempimento", color = BLU, linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
p3 = plot(; title = "Portata a fine riempimento", ylabel = "portata [m³/s]", xlabel = "distanza [km]", legend = false, st...)
plot!(p3, x, 0.5 .* (QI[end, 1:end-1] .+ QI[end, 2:end]); color = BLU, linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
savefig(plot(p1, p2, p3; layout = grid(3, 1, heights = [0.45, 0.3, 0.25]), size = (1000, 900), dpi = 130, left_margin = 6Plots.mm),
        joinpath(OUT, "1_condizioni_iniziali.png"))

# 2) idrogrammi: ingresso, metà tratto e uscita
im = cld(imax, 2)
p = plot(; title = "Idrogrammi di portata", xlabel = "tempo [h]", ylabel = "portata [m³/s]", legend = :topright, st...)
plot!(p, getindex.(P, 1) ./ 3600, getindex.(P, 2); label = "idrogramma di monte (" * basename(MONTE) * ")", color = :black, linestyle = :dash, linewidth = 1.2)
plot!(p, t, Q[:, 1]; label = @sprintf("ingresso (x = %.1f km)", x[1]), color = BLU, linewidth = 2)
plot!(p, t, Qs[:, im]; label = @sprintf("metà tratto (x = %.1f km)", x[im]), color = VER, linewidth = 2)
plot!(p, t, Q[:, iout]; label = @sprintf("uscita (x = %.1f km)", x[end]), color = ARA, linewidth = 2)
riep = IOBuffer()
for (v, n) in ((Q[:, 1], "ingresso"), (Qs[:, im], "metà"), (Q[:, iout], "uscita"))
    k = argmax(v); @printf(riep, "picco %-8s %.2f m3/s a %.2f h\n", n, v[k], t[k])
end
savefig(plot(p; size = (1000, 500), dpi = 130, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), joinpath(OUT, "2_idrogrammi.png"))

# 3) tirante nel tempo in alcune sezioni
sez = unique([1, cld(imax, 4), im, 3*cld(imax, 4), imax])
cols = [BLU, VER, ARA, ROS, GRI]
p = plot(; title = "Tirante nel tempo", xlabel = "tempo [h]", ylabel = "tirante [m]", legend = :topright, st...)
for (j, i) in enumerate(sez)
    plot!(p, t, H[:, i]; label = @sprintf("sezione %d (x = %.1f km)", i, x[i]), color = cols[j], linewidth = 1.8)
end
savefig(plot(p; size = (1000, 500), dpi = 130, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), joinpath(OUT, "3_tiranti_sezioni.png"))

# 4) profili longitudinali: inizio, picco a monte, picco a valle, fine, inviluppo dei massimi
kp1 = argmax(Q[:, 1]); kp2 = argmax(Q[:, iout])
emax = vec(maximum(E; dims = 1))
p = plot(; title = "Profili del pelo libero", xlabel = "distanza [km]", ylabel = "quota [m s.l.m.]", legend = :topright, st...)
plot!(p, x, zb; label = "talweg", color = GRI, linewidth = 1.5, fillrange = minimum(zb) - 1, fillcolor = GRI, fillalpha = 0.15)
plot!(p, x, ztab; label = "quota massima delle tabelle", color = :black, linestyle = :dot, linewidth = 1)
plot!(p, x, emax; label = "inviluppo dei massimi", color = ROS, linewidth = 1.2, linestyle = :dash)
plot!(p, x, E[1, :]; label = @sprintf("t = %.1f h (inizio)", t[1]), color = BLU, linewidth = 1.6)
plot!(p, x, E[kp1, :]; label = @sprintf("t = %.1f h (picco in ingresso)", t[kp1]), color = ARA, linewidth = 1.6)
plot!(p, x, E[kp2, :]; label = @sprintf("t = %.1f h (picco in uscita)", t[kp2]), color = VER, linewidth = 1.6)
plot!(p, x, E[end, :]; label = @sprintf("t = %.1f h (fine)", t[end]), color = :black, linewidth = 1.2)
savefig(plot(p; size = (1100, 550), dpi = 130, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), joinpath(OUT, "4_profili.png"))

# 5) livello massimo meno quota massima della tabella (positivo = valori estrapolati)
sup = emax .- ztab
p = bar(x, sup; title = "Livello massimo meno quota massima della tabella (positivo = estrapolazione)",
        xlabel = "distanza [km]", ylabel = "[m]", legend = false, color = [s > 0 ? ARA : BLU for s in sup], linecolor = :white,
        bar_width = 0.75 * minimum(diff(x)), st...)
hline!(p, [0.0]; color = :black, linewidth = 0.8)
savefig(plot(p; size = (1000, 420), dpi = 130, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), joinpath(OUT, "5_estrapolazione_tabelle.png"))
for i in findall(sup .> 0)
    @printf(riep, "sezione %d (x = %.2f km): livello massimo %.2f m sopra l'ultima quota della tabella\n", i, x[i], sup[i])
end

# 6) mappe di portata e tirante nel tempo e nello spazio
h1 = heatmap(x, t, Qs; title = "Portata [m³/s]", xlabel = "distanza [km]", ylabel = "tempo [h]", color = :Blues, st...)
h2 = heatmap(x, t, H; title = "Tirante [m]", xlabel = "distanza [km]", ylabel = "tempo [h]", color = :Oranges, st...)
savefig(plot(h1, h2; layout = (1, 2), size = (1200, 500), dpi = 130, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm), joinpath(OUT, "6_mappe.png"))

# 7) GIF del passaggio della piena
qmax = 1.05 * maximum(Q)
anim = @animate for k in 1:6:length(t)
    a = plot(x, zb; label = "talweg", color = GRI, fillrange = minimum(zb) - 1, fillcolor = GRI, fillalpha = 0.15,
             ylims = (minimum(zb) - 1, max(maximum(ztab), maximum(emax)) + 1), ylabel = "quota [m s.l.m.]", legend = :topright, st...)
    plot!(a, x, emax; label = "inviluppo dei massimi", color = ROS, linestyle = :dash, linewidth = 1)
    plot!(a, x, E[k, :]; label = "pelo libero", color = BLU, linewidth = 2)
    b = plot(x, Qs[k, :]; ylims = (0, qmax), ylabel = "portata [m³/s]", xlabel = "distanza [km]", legend = false, color = ARA, linewidth = 2, st...)
    plot(a, b; layout = (2, 1), size = (1000, 700), dpi = 100, left_margin = 5Plots.mm,
         plot_title = @sprintf("t = %5.2f h, portata in ingresso %.1f m³/s", t[k], Q[k, 1]), plot_titlefontsize = 11)
end
gif(anim, joinpath(OUT, "7_piena.gif"); fps = 12)

s = String(take!(riep)); write(joinpath(OUT, "riepilogo.txt"), s); println(s); println(readdir(OUT))
