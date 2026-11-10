# Evoluzione del riempimento (spin-up) fino a convergenza: figura e animazione dai file IC_*.csv.
# Legge da settings.txt geometria e cartella degli output, e salva in <cartella degli output>/grafici
# 12_spinup_convergenza.png e 13_spinup_animazione.gif. Non lancia nessuna simulazione.
# Uso, dalla cartella del progetto:  julia --project=. strumenti/grafici_spinup.jl
using Plots, Printf
const R   = dirname(@__DIR__)                                    # cartella del progetto
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const GEO = joinpath(R, SET[2])                                  # geometria (riga 2)
const D   = joinpath(R, length(ARGS) >= 1 ? ARGS[1] : SET[22])   # cartella degli output (riga 22, o primo argomento)
const OUT = joinpath(D, "grafici"); mkpath(OUT)
const BLU = "#2a78d6"; const ARA = "#eb6834"; const VER = "#1baf7a"; const ROS = "#d55181"; const GRI = "#7a7a7a"
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

mat(f)  = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
vec1(f) = [parse(Float64, strip(l)) for l in eachline(joinpath(D, f)) if !isempty(strip(l))]
G = [parse.(Float64, split(l, ',')) for l in readlines(GEO)[2:end] if !isempty(strip(l))]
x = getindex.(G, 3) ./ 1000; imax = length(x)

ts = vec1("IC_time.csv"); t = ts ./ 3600
Q = mat("IC_Q.csv"); E = mat("IC_eta.csv"); H = mat("IC_h.csv")
Qs = 0.5 .* (Q[:, 1:end-1] .+ Q[:, 2:end])
n = length(t); Qf = Q[n, 1]

# Criterio di arresto del modello (a parte la durata minima): |q - qn| < 0,01 e uscita = ingresso entro lo 0,1%
ok = [maximum(abs.(Q[k, 1:imax] .- Q[k-1, 1:imax])) < 1e-2 && abs(Q[k, 1] - Q[k, end]) < 1e-3 * Q[k, 1] for k in 2:n]
kc = findlast(!, ok) + 2; tc = t[kc]

dQ = vec(maximum(abs.(Q .- Q[n, :]'); dims = 2)); dE = vec(maximum(abs.(E .- E[n, :]'); dims = 2))
sQ, sE = round(0.01 * Qf; digits = 2), 0.01  # soglie: 1% della portata finale, 1 cm
tQ = t[findlast(>=(sQ), dQ) + 1]; tE = t[findlast(>=(sE), dE) + 1]
@printf("criterio del modello soddisfatto da t = %.0f s (%.2f h)\n", ts[kc], tc)
@printf("|Q-Qfin| < %.2f da %.2f h; |eta-etafin| < %.2f da %.2f h\n", sQ, tQ, sE, tE)
@printf("uscita minima %.2f m3/s a %.2f h; Q(t=5s) in [%.2f, %.2f]\n", minimum(Q[:, end]), t[argmin(Q[:, end])], extrema(Q[1, :])...)

# Sezioni come in 3_tiranti_sezioni.png
im = cld(imax, 2)
sez = unique([1, cld(imax, 4), im, 3 * cld(imax, 4), imax]); cols = [BLU, VER, ARA, ROS, GRI]
serieQ(i) = i == 1 ? Q[:, 1] : i == imax ? Q[:, end] : Qs[:, i]
nomi(i) = i == 1 ? @sprintf("ingresso (x = %.1f km)", x[i]) : i == imax ? @sprintf("uscita (x = %.1f km)", x[i]) : @sprintf("x = %.1f km", x[i])

# --- Figura riassuntiva ---
a = plot(; title = "Portata nel tempo", ylabel = "portata [m³/s]", legend = :bottomright, xlims = (0, t[end]), st...)
b = plot(; title = "Tirante nel tempo", ylabel = "tirante [m]", legend = :bottomright, xlims = (0, t[end]), ylims = (0, 1.08 * maximum(H[:, sez])), xlabel = "tempo di riempimento [h]", st...)
for (j, i) in enumerate(sez)
    plot!(a, t, serieQ(i); label = nomi(i), color = cols[j], linewidth = 1.8)
    plot!(b, t, H[:, i]; label = nomi(i), color = cols[j], linewidth = 1.8)
end
for p in (a, b)
    vline!(p, [tc]; color = :black, linestyle = :dash, linewidth = 1, label = @sprintf("criterio di arresto soddisfatto (%.2f h)", tc))
end
mask(v) = v .>= 1e-5
c = plot(t[mask(dQ)], dQ[mask(dQ)]; yscale = :log10, ylims = (1e-5, 10), xlims = (0, t[end]), color = BLU, linewidth = 2, label = "massimo sulle interfacce",
         title = "Distanza dallo stato finale: portata", ylabel = "max |Q − Q_finale| [m³/s]", xlabel = "tempo di riempimento [h]", legend = :topright, st...)
hline!(c, [sQ]; color = GRI, linestyle = :dot, label = @sprintf("%.2f m³/s (1%%), da %.2f h", sQ, tQ))
d = plot(t[mask(dE)], dE[mask(dE)]; yscale = :log10, ylims = (1e-5, 10), xlims = (0, t[end]), color = ARA, linewidth = 2, label = "massimo sulle sezioni",
         title = "Distanza dallo stato finale: livello", ylabel = "max |η − η_finale| [m]", xlabel = "tempo di riempimento [h]", legend = :topright, st...)
hline!(d, [sE]; color = GRI, linestyle = :dot, label = @sprintf("1 cm, da %.2f h", tE))
savefig(plot(a, b, c, d; layout = @layout([a{0.34h}; b{0.34h}; [c d]]), size = (1050, 980), dpi = 130, left_margin = 6Plots.mm, bottom_margin = 2Plots.mm,
             plot_title = @sprintf("%s: spin-up fino a convergenza (Q di monte = %.2f m³/s)", basename(D), Qf), plot_titlefontsize = 12),
        joinpath(OUT, "12_spinup_convergenza.png"))

# --- Animazione: passo fitto nel transitorio iniziale, largo verso la convergenza ---
fr = filter(<=(n), unique(vcat(1:3:60, 60:12:720, 720:24:2160, 2160:180:n, n)))
hmax = 1.08 * maximum(H); Qmax = 1.05 * maximum(Q)
hmsf(s) = @sprintf("%d h %02d min %02d s", floor(Int, s / 3600), floor(Int, mod(s, 3600) / 60), round(Int, mod(s, 60)))
function scena(k)
    a = plot(x, H[n, :]; label = "stato finale", color = GRI, linewidth = 3, alpha = 0.55, ylims = (0, hmax), xlims = (0, 1.03 * x[end]),
             ylabel = "tirante [m]", legend = :topright, title = "Tirante lungo il tratto", st...)
    plot!(a, x, H[1, :]; label = "partenza (t = 5 s)", color = ARA, linestyle = :dash, linewidth = 1.2)
    plot!(a, x, H[k, :]; label = "stato corrente", color = BLU, linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
    b = plot(x, fill(Qf, imax); label = "portata di monte", color = GRI, linewidth = 3, alpha = 0.55, ylims = (0, Qmax), xlims = (0, 1.03 * x[end]),
             ylabel = "portata [m³/s]", xlabel = "distanza [km]", legend = :topright, title = "Portata lungo il tratto", st...)
    plot!(b, x, Qs[k, :]; label = "stato corrente", color = BLU, linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
    c = plot(t, Q[:, end]; label = false, color = GRI, linewidth = 1, alpha = 0.5, ylims = (0, Qmax), xlims = (0, t[end]),
             ylabel = "portata [m³/s]", xlabel = "tempo di riempimento [h]", legend = :bottomright, title = @sprintf("Portata in uscita (x = %.1f km)", x[end]), st...)
    plot!(c, t[1:k], Q[1:k, end]; label = false, color = ARA, linewidth = 2)
    vline!(c, [tc]; color = :black, linestyle = :dash, linewidth = 1, label = @sprintf("criterio di arresto soddisfatto (%.2f h)", tc))
    scatter!(c, [t[k]], [Q[k, end]]; color = ARA, markersize = 6, markerstrokecolor = :white, label = false)
    stato = ts[k] >= ts[kc] ? "convergenza raggiunta" : "in convergenza"
    plot(a, b, c; layout = @layout([a{0.42h}; b{0.33h}; c]), size = (1000, 820), dpi = 100, left_margin = 6Plots.mm, bottom_margin = 2Plots.mm,
         plot_title = "Spin-up: t = " * hmsf(ts[k]) * " – " * stato * " (animazione a velocità variabile)", plot_titlefontsize = 12)
end
anim = Animation()
for k in fr; frame(anim, scena(k)); end
ultimo = scena(n)
for _ in 1:24; frame(anim, ultimo); end      # pausa sull'ultimo fotogramma (circa 2 s)
gif(anim, joinpath(OUT, "13_spinup_animazione.gif"); fps = 12)
println("fotogrammi: ", length(fr))
