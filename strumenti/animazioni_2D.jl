# Animazioni della piena con il riferimento 2D: 10_piena_2D.gif (solo 2D) e 11_piena_1D_2D.gif
# (1D e 2D insieme, con la differenza di quota). I dati 2D sono livelli e portate per sezione, un'uscita
# ogni 15 min (reference/ns_50_<evento>_wse.csv e reference/ns_50_<evento>_discharge.csv); le colonne
# sono associate alle sezioni del 1D arrotondando l'ascissa ai 50 m, e il 1D è interpolato linearmente
# negli istanti del 2D. Le sezioni 2D ancora asciutte (livello o portata nulli) non sono disegnate.
# Tabelle e geometria sono quelle delle righe 1-2 di settings.txt: devono essere quelle della simulazione.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/animazioni_2D.jl <cartella degli output> <livelli 2D> <portate 2D>
using Plots, Printf
const R   = dirname(@__DIR__)                                    # cartella del progetto
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const TAB = joinpath(R, SET[1])                                  # tabelle delle sezioni (riga 1)
const GEO = joinpath(R, SET[2])                                  # geometria (riga 2)
length(ARGS) == 3 || error("Uso: animazioni_2D.jl <cartella degli output> <livelli 2D> <portate 2D>")
const D    = joinpath(R, ARGS[1])
const WSE2 = joinpath(R, ARGS[2])
const Q2F  = joinpath(R, ARGS[3])
const OUT  = joinpath(D, "grafici"); mkpath(OUT)
const BLU = "#2a78d6"; const ARA = "#eb6834"; const ROS = "#d55181"; const GRI = "#7a7a7a"
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
G = csvnum(GEO); T = csvnum(TAB)
x = getindex.(G, 3); zb = getindex.(G, 4); imax = length(x); xk = x ./ 1000
ztab = [maximum(r[2] for r in T if r[1] == k) for k in 1:imax]
nomesez(xi) = @sprintf("S%04d", round(Int, xi/50)*50)

# 1D: portata nelle sezioni come media delle interfacce
t1 = [parse(Float64, strip(l)) for l in eachline(joinpath(D, "time.csv")) if !isempty(strip(l))]
Q1 = mat("Q.csv"); E1 = mat("eta.csv"); Qs1 = 0.5 .* (Q1[:, 1:end-1] .+ Q1[:, 2:end])

# 2D: colonne delle sezioni del 1D; la riga a t = 0 è tutta nulla ed è esclusa
function leggi2D(f, pref)
    L = readlines(f); hdr = strip.(split(L[1], ','))
    M = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
    c = [findfirst(==("$(pref)_$(nomesez(xi)) [" * (pref == "WSE" ? "m" : "m3/s") * "]"), hdr) for xi in x]
    any(isnothing, c) && error("colonne mancanti in $f: " * join(nomesez.(x[isnothing.(c)]), ", "))
    return M[:, 1], M[:, c]
end
tw, E2 = leggi2D(WSE2, "WSE"); tq, Q2 = leggi2D(Q2F, "Q")
tw == tq || error("istanti diversi tra livelli e portate 2D")
k2 = findall(τ -> τ > 0 && t1[1] <= τ <= t1[end], tw)          # istanti 2D dentro la simulazione 1D
t2 = tw[k2]; E2 = E2[k2, :]; Q2 = Q2[k2, :]
asciutto = (E2 .<= 1.0) .| (Q2 .<= 0.0)
E2 = ifelse.(asciutto, NaN, E2); Q2 = ifelse.(asciutto, NaN, Q2)

# 1D negli istanti del 2D
function interp_t(M, τ)
    i = clamp(searchsortedlast(t1, τ), 1, length(t1) - 1); θ = (τ - t1[i])/(t1[i+1] - t1[i])
    return (1 - θ) .* M[i, :] .+ θ .* M[i+1, :]
end
E1c = reduce(vcat, [permutedims(interp_t(E1, τ)) for τ in t2])
Q1c = reduce(vcat, [permutedims(interp_t(Qs1, τ)) for τ in t2])
Qin = [interp_t(Q1[:, 1:1], τ)[1] for τ in t2]
dE  = E1c .- E2

nanmax(v) = maximum(filter(!isnan, v))
emax2 = [nanmax(E2[:, i]) for i in 1:imax]
ylz = (minimum(zb) - 1, max(maximum(ztab), maximum(E1), nanmax(E2)) + 1)
ylq = (0, 1.05 * max(maximum(Q1), nanmax(Q2)))
dmax = 1.05 * nanmax(abs.(dE))
ttl(j) = @sprintf("t = %5.2f h", t2[j]/3600)

function quadro10(j)
    a = plot(xk, zb; label = "talweg", color = GRI, fillrange = minimum(zb) - 1, fillcolor = GRI, fillalpha = 0.15,
             ylims = ylz, ylabel = "quota [m s.l.m.]", legend = :topright, st...)
    plot!(a, xk, emax2; label = "inviluppo dei massimi", color = ROS, linestyle = :dash, linewidth = 1)
    plot!(a, xk, E2[j, :]; label = "pelo libero (2D)", color = BLU, linewidth = 2)
    b = plot(xk, Q2[j, :]; ylims = ylq, ylabel = "portata [m³/s]", xlabel = "distanza [km]", legend = false, color = ARA, linewidth = 2, st...)
    plot(a, b; layout = (2, 1), size = (1000, 700), dpi = 100, left_margin = 5Plots.mm,
         plot_title = "Simulazione 2D: " * ttl(j) * @sprintf(", portata in %s %.1f m³/s", nomesez(x[1]), Q2[j, 1]), plot_titlefontsize = 11)
end

function quadro11(j)
    m2 = (color = ARA, linewidth = 2, markershape = :circle, markersize = imax > 60 ? 1.5 : 3, markerstrokewidth = 0)
    a = plot(xk, zb; label = "talweg", color = GRI, fillrange = minimum(zb) - 1, fillcolor = GRI, fillalpha = 0.15,
             ylims = ylz, ylabel = "quota [m s.l.m.]", legend = :topright, st...)
    plot!(a, xk, E1c[j, :]; label = "1D con sfioro", color = BLU, linewidth = 2)
    plot!(a, xk, E2[j, :]; label = "2D", m2...)
    c = plot(xk, dE[j, :]; ylims = (-dmax, dmax), ylabel = "quota 1D − 2D [m]", legend = false, color = :black, linewidth = 1.5, st...)
    hline!(c, [0.0]; color = GRI, linewidth = 0.8)
    b = plot(xk, Q1c[j, :]; ylims = ylq, ylabel = "portata [m³/s]", xlabel = "distanza [km]", label = "1D con sfioro",
             color = BLU, linewidth = 2, legend = :bottomleft, st...)
    plot!(b, xk, Q2[j, :]; label = "2D", m2...)
    plot(a, c, b; layout = grid(3, 1, heights = [0.46, 0.2, 0.34]), size = (1000, 850), dpi = 100, left_margin = 6Plots.mm,
         plot_title = ttl(j) * @sprintf(", portata in ingresso: 1D %.1f m³/s, 2D (%s) %.1f m³/s", Qin[j], nomesez(x[1]), Q2[j, 1]),
         plot_titlefontsize = 11)
end

# 5 fotogrammi al secondo: con un'uscita 2D ogni 15 min, circa 1 h simulata ogni 0,8 s
gif(@animate(for j in eachindex(t2); quadro10(j); end), joinpath(OUT, "10_piena_2D.gif"); fps = 5)
gif(@animate(for j in eachindex(t2); quadro11(j); end), joinpath(OUT, "11_piena_1D_2D.gif"); fps = 5)

v = filter(!isnan, dE); I = argmax(replace(abs.(dE), NaN => -1.0))
@printf("%d fotogrammi; quota 1D - 2D: media %+.3f m, scarto quadratico medio %.3f m, massimo %+.2f m in %s a %.2f h\n",
        length(t2), sum(v)/length(v), sqrt(sum(abs2, v)/length(v)), dE[I], nomesez(x[I[2]]), t2[I[1]]/3600)
