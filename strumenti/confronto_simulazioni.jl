# Confronto tra più simulazioni dello stesso tratto (dopo il riempimento), dai file Q.csv, eta.csv, h.csv, time.csv.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/confronto_simulazioni.jl <cartella dei grafici> <cartella 1> <cartella 2> ...
# Le differenze sono calcolate rispetto alla prima cartella. La geometria è quella della riga 2 di settings.txt.
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
x = getindex.(G, 3) ./ 1000; zb = getindex.(G, 4); imax = length(x)

D = map(CASI) do d
    t = vec(mat(d, "time.csv")) ./ 3600; Q = mat(d, "Q.csv"); E = mat(d, "eta.csv"); H = mat(d, "h.csv")
    Qs = 0.5 .* (Q[:, 1:end-1] .+ Q[:, 2:end])
    (; nome = basename(d), t, Q, Qs, E, H)
end
uguali(a, b) = size(a.Q) == size(b.Q) && a.Q == b.Q && a.E == b.E
etich = [begin
    gemelli = [D[j].nome for j in eachindex(D) if j != i && uguali(D[i], D[j])]
    isempty(gemelli) ? D[i].nome : D[i].nome * " (identico a " * join(gemelli, ", ") * ")"
end for i in eachindex(D)]
B = D[1]
stessi_t(c) = length(c.t) == length(B.t) && maximum(abs.(c.t .- B.t)) < 1e-9

# 1) idrogrammi in ingresso, a metà tratto e in uscita
im = cld(imax, 2)
serie(c, k) = k == 1 ? c.Q[:, 1] : k == imax ? c.Q[:, end] : c.Qs[:, k]
pp = Any[]
for (k, nome) in ((1, "ingresso"), (im, "metà tratto"), (imax, "uscita"))
    p = plot(; title = @sprintf("Portata in %s (x = %.1f km)", nome, x[k]), ylabel = "portata [m³/s]", legend = :topright, st...)
    for (i, c) in enumerate(D)
        plot!(p, c.t, serie(c, k); color = COL[i], linestyle = STY[i], linewidth = 2, label = etich[i])
    end
    push!(pp, p)
end
xlabel!(pp[end], "tempo [h]")
savefig(plot(pp...; layout = (3, 1), size = (1000, 950), dpi = 130, left_margin = 6Plots.mm), joinpath(OUT, "4_idrogrammi.png"))

# 2) differenze nel tempo rispetto al primo caso: massimo lungo il tratto di |ΔQ| e |Δη|
d1 = plot(; title = "Differenza massima lungo il tratto rispetto a $(B.nome): portata", ylabel = "max |ΔQ| [m³/s]", legend = :topright, st...)
d2 = plot(; title = "Differenza massima lungo il tratto rispetto a $(B.nome): livello", ylabel = "max |Δη| [m]", xlabel = "tempo [h]", legend = :topright, st...)
for (i, c) in enumerate(D[2:end])
    stessi_t(c) || (@warn "$(c.nome): istanti diversi da $(B.nome), differenze nel tempo saltate"; continue)
    plot!(d1, c.t, vec(maximum(abs.(c.Q .- B.Q); dims = 2)); color = COL[i+1], linestyle = STY[i+1], linewidth = 2, label = etich[i+1])
    plot!(d2, c.t, vec(maximum(abs.(c.E .- B.E); dims = 2)); color = COL[i+1], linestyle = STY[i+1], linewidth = 2, label = etich[i+1])
end
savefig(plot(d1, d2; layout = (2, 1), size = (1000, 700), dpi = 130, left_margin = 6Plots.mm), joinpath(OUT, "5_differenze_nel_tempo.png"))

# 3) inviluppi dei massimi lungo il tratto
e1 = plot(; title = "Inviluppo dei massimi: tirante", ylabel = "tirante [m]", legend = :topleft, st...)
e2 = plot(; title = "Inviluppo dei massimi: portata", ylabel = "portata [m³/s]", xlabel = "distanza [km]", legend = :topleft, st...)
for (i, c) in enumerate(D)
    kw = (color = COL[i], linestyle = STY[i], linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0, label = etich[i])
    plot!(e1, x, vec(maximum(c.H; dims = 1)); kw...); plot!(e2, x, vec(maximum(c.Qs; dims = 1)); kw...)
end
savefig(plot(e1, e2; layout = (2, 1), size = (1000, 700), dpi = 130, left_margin = 6Plots.mm), joinpath(OUT, "6_inviluppi.png"))

io = IOBuffer()
for c in D
    k = argmax(c.Q[:, end])
    @printf(io, "%-10s colmo in uscita %.3f m3/s a %.2f h | tirante massimo nel tratto %.3f m\n", c.nome, c.Q[k, end], c.t[k], maximum(c.H))
end
for c in D[2:end]
    if stessi_t(c)
        @printf(io, "%s - %s: su tutta la simulazione max |ΔQ| %.4f m3/s, max |Δη| %.4f m\n", c.nome, B.nome, maximum(abs.(c.Q .- B.Q)), maximum(abs.(c.E .- B.E)))
    else
        @printf(io, "%s - %s: istanti diversi, differenze punto per punto non calcolate\n", c.nome, B.nome)
    end
end
s = String(take!(io)); write(joinpath(OUT, "riepilogo_simulazioni.txt"), s); println(s)
