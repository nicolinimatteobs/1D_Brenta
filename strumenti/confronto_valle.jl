# Confronto nella sezione di valle tra il modello 1D e il riferimento 2D: tirante e portata nel tempo.
# Grafico 12_confronto_valle_tirante_portata.png e numeri in confronto_valle.txt, in <cartella degli output>/grafici.
# Il tirante 2D è il livello 2D meno il talweg della geometria 1D (riga 2 di settings.txt); il 1D è h.csv.
# La portata 1D nella sezione è la media delle due interfacce, come negli altri strumenti. Gli istanti in
# cui il 2D è ancora asciutto (livello o portata nulli, 2D che parte da alveo asciutto) non sono disegnati
# e sono esclusi dai numeri. Le colonne 2D sono associate alle sezioni arrotondando l'ascissa ai 50 m.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/confronto_valle.jl <cartella degli output> <livelli 2D> <portate 2D>
using Plots, Printf
const R   = dirname(@__DIR__)                                    # cartella del progetto
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const GEO = joinpath(R, SET[2])                                  # geometria (riga 2)
length(ARGS) == 3 || error("Uso: confronto_valle.jl <cartella degli output> <livelli 2D> <portate 2D>")
const D    = joinpath(R, ARGS[1])
const WSE2 = joinpath(R, ARGS[2])
const Q2F  = joinpath(R, ARGS[3])
const OUT  = joinpath(D, "grafici"); mkpath(OUT)
const BLU = "#2a78d6"; const ARA = "#eb6834"; const GRI = "#7a7a7a"
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
G  = csvnum(GEO); i = length(G)                                  # sezione di valle
xv = G[i][3]; zb = G[i][4]
argini = length(G[i]) >= 6 ? (G[i][5], G[i][6]) : (Inf, Inf)
nome = @sprintf("S%04d", round(Int, xv/50)*50)

# 1D
t1 = [parse(Float64, strip(l)) for l in eachline(joinpath(D, "time.csv")) if !isempty(strip(l))]
Q1 = mat("Q.csv"); H1 = mat("h.csv")
q1 = 0.5 .* (Q1[:, i] .+ Q1[:, i+1]); h1 = H1[:, i]; qin = Q1[:, 1]

# 2D
function colonna(f, n)
    L = readlines(f); hdr = strip.(split(L[1], ','))
    c = findfirst(==(n), hdr); c === nothing && error("colonna $n mancante in $f")
    M = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
    return M[:, 1], M[:, c]
end
tw, e2 = colonna(WSE2, "WSE_$nome [m]"); tq, q2 = colonna(Q2F, "Q_$nome [m3/s]")
tw == tq || error("istanti diversi tra livelli e portate 2D")
k  = findall(τ -> t1[1] <= τ <= t1[end], tw)                     # periodo comune
t2 = tw[k]; e2 = e2[k]; q2 = q2[k]
bagnato = (e2 .> 1.0) .& (q2 .> 0.0)
h2 = ifelse.(bagnato, e2 .- zb, NaN); q2 = ifelse.(bagnato, q2, NaN)

interp(v, τ) = (j = clamp(searchsortedlast(t1, τ), 1, length(t1)-1); v[j] + (v[j+1]-v[j])*(τ - t1[j])/(t1[j+1]-t1[j]))
m   = bagnato
h1c = [interp(h1, τ) for τ in t2[m]]; q1c = [interp(q1, τ) for τ in t2[m]]
dh  = h1c .- h2[m]; dq = q1c .- q2[m]

# grafico
th = t1 ./ 3600; t2h = t2 ./ 3600
m2 = (color = ARA, linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
a = plot(; title = "Sezione $nome: tirante", ylabel = "tirante [m]", legend = :topright, xlims = (0, th[end]), st...)
plot!(a, th, h1; label = "1D con sfioro", color = BLU, linewidth = 2)
plot!(a, t2h, h2; label = "2D", m2...)
for (n, z) in zip(("argine destro", "argine sinistro"), argini)
    isfinite(z) && z - zb < 1.3 * max(maximum(h1), maximum(filter(!isnan, h2))) &&
        hline!(a, [z - zb]; label = n, color = GRI, linestyle = n == "argine destro" ? :dash : :dot, linewidth = 1)
end
b = plot(; title = "Sezione $nome: portata", xlabel = "tempo [h]", ylabel = "portata [m³/s]", legend = :topright, xlims = (0, th[end]), st...)
plot!(b, th, qin; label = "ingresso da monte (1D)", color = GRI, linestyle = :dash, linewidth = 1)
plot!(b, th, q1; label = "1D con sfioro", color = BLU, linewidth = 2)
plot!(b, t2h, q2; label = "2D", m2...)
savefig(plot(a, b; layout = (2, 1), size = (1000, 750), dpi = 130, left_margin = 6Plots.mm),
        joinpath(OUT, "12_confronto_valle_tirante_portata.png"))

# numeri
rms(v) = sqrt(sum(abs2, v)/length(v))
k1 = argmax(h1); k2 = argmax(replace(h2, NaN => -Inf))
j1 = argmax(q1); j2 = argmax(replace(q2, NaN => -Inf))
io = IOBuffer()
@printf(io, "Sezione %s (talweg %.2f m s.l.m.), 1D - 2D, %d istanti 2D bagnati (dati ogni 15 min)\n", nome, zb, count(m))
@printf(io, "tirante massimo: 2D %.3f m a %.2f h | 1D %.3f m a %.2f h | differenza %+.3f m\n",
        h2[k2], t2h[k2], h1[k1], th[k1], h1[k1] - h2[k2])
@printf(io, "tirante: differenza media %+.3f m, scarto quadratico medio %.3f m, massimo %+.3f m\n",
        sum(dh)/length(dh), rms(dh), dh[argmax(abs.(dh))])
@printf(io, "portata massima: 2D %.2f m3/s a %.2f h | 1D %.2f m3/s a %.2f h | differenza %+.2f m3/s\n",
        q2[j2], t2h[j2], q1[j1], th[j1], q1[j1] - q2[j2])
@printf(io, "portata: differenza media %+.2f m3/s, scarto quadratico medio %.2f m3/s\n", sum(dq)/length(dq), rms(dq))
s = String(take!(io)); write(joinpath(OUT, "confronto_valle.txt"), s); print(s)
