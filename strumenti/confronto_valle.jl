# Confronto tra il modello 1D e il riferimento 2D in una o più sezioni: tirante e portata nel tempo.
# - Senza sezioni: sezione di valle, grafico 12_confronto_valle_tirante_portata.png (tirante e portata nel
#   tempo, più il tirante in funzione della portata) e numeri in confronto_valle.txt.
# - Con i nomi di una o più sezioni (es. S1000 S2150 S4000 S6000): grafico 13_confronto_sezioni_tirante_portata.png,
#   una riga per sezione con tirante e portata nel tempo, e numeri in confronto_sezioni.txt.
# Tutto in <cartella degli output>/grafici.
# Il tirante 2D è il livello 2D meno il talweg della geometria 1D (riga 2 di settings.txt); il 1D è h.csv.
# La portata 1D nella sezione è la media delle due interfacce, come negli altri strumenti. Gli istanti in
# cui il 2D è ancora asciutto (livello o portata nulli, 2D che parte da alveo asciutto) non sono disegnati
# e sono esclusi dai numeri. Le colonne 2D sono associate alle sezioni arrotondando l'ascissa ai 50 m.
# Il confronto a parità di portata usa la scala di deflusso degli scalini 2D (moto permanente, da cui vengono
# le tabelle), letta da calcolo_tabelle_daBM/input/S_0050_6000 se c'è il file della sezione: lo scarto è il
# tirante meno quello degli scalini alla stessa portata.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/confronto_valle.jl <cartella degli output> <livelli 2D> <portate 2D> [sezioni...]
using Plots, Printf
const R   = dirname(@__DIR__)                                    # cartella del progetto
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const GEO = joinpath(R, SET[2])                                  # geometria (riga 2)
length(ARGS) >= 3 || error("Uso: confronto_valle.jl <cartella degli output> <livelli 2D> <portate 2D> [sezioni...]")
const D    = joinpath(R, ARGS[1])
const WSE2 = joinpath(R, ARGS[2])
const Q2F  = joinpath(R, ARGS[3])
const OUT  = joinpath(D, "grafici"); mkpath(OUT)
const SCALINI = joinpath(R, "calcolo_tabelle_daBM/input/S_0050_6000/ns_50_scalini_")
const BLU = "#2a78d6"; const ARA = "#eb6834"; const GRI = "#7a7a7a"
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
function tabella(f)
    L = readlines(f); hdr = strip.(split(L[1], ','))
    return hdr, reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
end
G = csvnum(GEO)
nomesez(xi) = @sprintf("S%04d", round(Int, xi/50)*50)

# 1D e 2D, letti una volta sola
t1 = [parse(Float64, strip(l)) for l in eachline(joinpath(D, "time.csv")) if !isempty(strip(l))]
Q1 = mat("Q.csv"); H1 = mat("h.csv"); qin = Q1[:, 1]
hw, Mw = tabella(WSE2); hq, Mq = tabella(Q2F)
Mw[:, 1] == Mq[:, 1] || error("istanti diversi tra livelli e portate 2D")
const k2 = findall(τ -> t1[1] <= τ <= t1[end], Mw[:, 1])        # periodo comune
const t2 = Mw[k2, 1]
interp(v, τ) = (j = clamp(searchsortedlast(t1, τ), 1, length(t1)-1); v[j] + (v[j+1]-v[j])*(τ - t1[j])/(t1[j+1]-t1[j]))
rms(v) = sqrt(sum(abs2, v)/length(v))

# Dati e numeri di una sezione (indice i nella geometria)
function sezione(i)
    xv = G[i][3]; zb = G[i][4]; nome = nomesez(xv)
    argini = length(G[i]) >= 6 ? (G[i][5], G[i][6]) : (Inf, Inf)
    q1 = 0.5 .* (Q1[:, i] .+ Q1[:, i+1]); h1 = H1[:, i]
    cw = findfirst(==("WSE_$nome [m]"), hw); cq = findfirst(==("Q_$nome [m3/s]"), hq)
    (cw === nothing || cq === nothing) && error("colonne 2D della sezione $nome mancanti")
    e2 = Mw[k2, cw]; q2 = Mq[k2, cq]
    m = (e2 .> 1.0) .& (q2 .> 0.0)                               # istanti 2D bagnati
    h2 = ifelse.(m, e2 .- zb, NaN); q2 = ifelse.(m, q2, NaN)
    h1c = [interp(h1, τ) for τ in t2[m]]; q1c = [interp(q1, τ) for τ in t2[m]]
    dh = h1c .- h2[m]; dq = q1c .- q2[m]
    # scala di deflusso degli scalini 2D, se c'è il file della sezione
    fsc = SCALINI * nome * ".csv"
    scalini = nothing
    if isfile(fsc)
        hs, Ms = tabella(fsc)
        qs = Ms[:, findfirst(==("Discharge [m3/s]"), hs)]; zs = Ms[:, findfirst(==("Mean wse [m]"), hs)] .- zb
        w = qs .> 0; p = sortperm(qs[w]); scalini = (qs[w][p], zs[w][p])
    else
        @warn "Scalini 2D non trovati ($fsc): niente confronto a parità di portata in $nome"
    end
    function hsc(q)                                              # tirante degli scalini alla portata q
        Q, H = scalini
        j = clamp(searchsortedlast(Q, q), 1, length(Q)-1)
        return H[j] + (H[j+1]-H[j])*(q - Q[j])/(Q[j+1]-Q[j])
    end
    io = IOBuffer()
    k1 = argmax(h1); kk = argmax(replace(h2, NaN => -Inf)); j1 = argmax(q1); j2 = argmax(replace(q2, NaN => -Inf))
    @printf(io, "Sezione %s (x = %.2f km, talweg %.2f m s.l.m.), 1D - 2D, %d istanti 2D bagnati (dati ogni 15 min)\n",
            nome, xv/1000, zb, count(m))
    @printf(io, "tirante massimo: 2D %.3f m a %.2f h | 1D %.3f m a %.2f h | differenza %+.3f m\n",
            h2[kk], t2[kk]/3600, h1[k1], t1[k1]/3600, h1[k1] - h2[kk])
    @printf(io, "tirante: differenza media %+.3f m, scarto quadratico medio %.3f m, massimo %+.3f m\n",
            sum(dh)/length(dh), rms(dh), dh[argmax(abs.(dh))])
    @printf(io, "portata massima: 2D %.2f m3/s a %.2f h | 1D %.2f m3/s a %.2f h | differenza %+.2f m3/s\n",
            q2[j2], t2[j2]/3600, q1[j1], t1[j1]/3600, q1[j1] - q2[j2])
    @printf(io, "portata: differenza media %+.2f m3/s, scarto quadratico medio %.2f m3/s\n", sum(dq)/length(dq), rms(dq))
    if scalini !== nothing
        d1 = h1c .- hsc.(q1c); d2 = h2[m] .- hsc.(q2[m])
        @printf(io, "a parità di portata, tirante - scalini 2D: 1D media %+.3f m (scarto quadratico medio %.3f m), 2D media %+.3f m (%.3f m)\n",
                sum(d1)/length(d1), rms(d1), sum(d2)/length(d2), rms(d2))
    end
    return (; nome, xv, zb, argini, q1, h1, q2, h2, scalini, testo = String(take!(io)))
end

th = t1 ./ 3600; t2h = t2 ./ 3600
m2 = (color = ARA, linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
function pannelli(s, titolo)
    a = plot(; title = "$titolo: tirante", ylabel = "tirante [m]", legend = :topright, xlims = (0, th[end]), st...)
    plot!(a, th, s.h1; label = "1D con sfioro", color = BLU, linewidth = 2)
    plot!(a, t2h, s.h2; label = "2D", m2...)
    hmax = max(maximum(s.h1), maximum(filter(!isnan, s.h2)))
    for (n, z) in zip(("argine destro", "argine sinistro"), s.argini)
        isfinite(z) && z - s.zb < 1.3 * hmax &&
            hline!(a, [z - s.zb]; label = n, color = GRI, linestyle = n == "argine destro" ? :dash : :dot, linewidth = 1)
    end
    b = plot(; title = "$titolo: portata", ylabel = "portata [m³/s]", legend = :topright, xlims = (0, th[end]), st...)
    plot!(b, th, qin; label = "ingresso da monte (1D)", color = GRI, linestyle = :dash, linewidth = 1)
    plot!(b, th, s.q1; label = "1D con sfioro", color = BLU, linewidth = 2)
    plot!(b, t2h, s.q2; label = "2D", m2...)
    return a, b
end

if length(ARGS) == 3
    # sezione di valle: tirante e portata nel tempo, tirante in funzione della portata
    s = sezione(length(G))
    a, b = pannelli(s, "Sezione di valle $(s.nome)")
    xlabel!(b, "tempo [h]")
    c = plot(; title = "Sezione di valle $(s.nome): tirante a parità di portata", xlabel = "portata [m³/s]",
             ylabel = "tirante [m]", legend = :bottomright, st...)
    s.scalini === nothing || plot!(c, s.scalini...; label = "scalini 2D (moto permanente)", color = :black, linewidth = 1,
                                   linestyle = :dash, markershape = :diamond, markersize = 3, markerstrokewidth = 0)
    plot!(c, s.q1, s.h1; label = "1D con sfioro", color = BLU, linewidth = 2)
    plot!(c, s.q2, s.h2; label = "2D", m2...)
    savefig(plot(a, b, c; layout = (3, 1), size = (1000, 1100), dpi = 130, left_margin = 6Plots.mm),
            joinpath(OUT, "12_confronto_valle_tirante_portata.png"))
    write(joinpath(OUT, "confronto_valle.txt"), s.testo); print(s.testo)
else
    # più sezioni: una riga per sezione, tirante a sinistra e portata a destra
    idx = map(ARGS[4:end]) do n
        i = findfirst(r -> nomesez(r[3]) == n, G)
        i === nothing && error("Sezione $n non trovata nella geometria $GEO")
        i
    end
    S = sezione.(idx)
    P = Any[]
    for s in S
        a, b = pannelli(s, @sprintf("%s (x = %.2f km)", s.nome, s.xv/1000))
        push!(P, a, b)
    end
    xlabel!(P[end-1], "tempo [h]"); xlabel!(P[end], "tempo [h]")
    n = length(S)
    savefig(plot(P...; layout = (n, 2), size = (1500, 380 * n), dpi = 110, left_margin = 6Plots.mm, bottom_margin = 3Plots.mm),
            joinpath(OUT, "13_confronto_sezioni_tirante_portata.png"))
    testo = join(getfield.(S, :testo), "\n")
    write(joinpath(OUT, "confronto_sezioni.txt"), testo); print(testo)
end
