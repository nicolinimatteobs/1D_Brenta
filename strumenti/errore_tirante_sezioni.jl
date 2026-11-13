# Errore del tirante 1D rispetto al 2D, sezione per sezione, negli istanti in cui le due portate coincidono.
# Per ogni scenario (sottocartelle fittizio_* della cartella indicata) e per ogni sezione si prendono gli
# istanti 2D (ogni 15 min) in cui il 2D è bagnato e |Q1D − Q2D| ≤ soglia·Q2D, e si calcola lo scarto
# quadratico medio del tirante 1D − 2D. Il 2D è il riferimento: tirante = livello 2D − talweg della geometria
# 1D (riga 2 di settings.txt). Il 1D è interpolato negli istanti 2D; la sua portata nella sezione è la media
# delle due interfacce, come negli altri strumenti. Le colonne 2D sono associate alle sezioni arrotondando
# l'ascissa ai 50 m; i riferimenti sono reference/ns_50_<scenario>_wse.csv e _discharge.csv.
# Risultati in <cartella>/grafici: 14_errore_tirante_sezioni.png ed errore_tirante_sezioni.txt.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/errore_tirante_sezioni.jl [cartella con gli scenari] [soglia relativa]
# (predefiniti: output/S_50m e 0.01, cioè portate uguali entro l'1%).
using Plots, Printf
const R    = dirname(@__DIR__)                                   # cartella del progetto
const SET  = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const GEO  = joinpath(R, SET[2])                                 # geometria (riga 2)
const BASE = joinpath(R, length(ARGS) >= 1 ? ARGS[1] : "output/S_50m")
const SOGLIA = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.01
const OUT  = joinpath(BASE, "grafici"); mkpath(OUT)
const COL  = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]        # un colore per scenario, sempre nello stesso ordine
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 11, guidefontsize = 10, tickfontsize = 8, legendfontsize = 9)

csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
G = csvnum(GEO); x = getindex.(G, 3); zb = getindex.(G, 4); n = length(x)
nomesez(xi) = @sprintf("S%04d", round(Int, xi/50)*50)
nomi = nomesez.(x)

function scenario(ev)
    D = joinpath(BASE, ev)
    mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
    t1 = [parse(Float64, strip(l)) for l in eachline(joinpath(D, "time.csv")) if !isempty(strip(l))]
    Q1 = mat("Q.csv"); H1 = mat("h.csv"); Qs = 0.5 .* (Q1[:, 1:end-1] .+ Q1[:, 2:end])
    size(H1, 2) == n || error("$ev: $(size(H1, 2)) sezioni negli output, $n nella geometria $GEO")
    function leggi2D(f, pref, u)
        isfile(f) || error("riferimento 2D mancante: $f")
        L = readlines(f); hdr = strip.(split(L[1], ','))
        M = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
        c = [findfirst(==("$(pref)_$(nm) [$u]"), hdr) for nm in nomi]
        any(isnothing, c) && error("colonne mancanti in $f: " * join(nomi[isnothing.(c)], ", "))
        return M[:, 1], M[:, c]
    end
    tw, E2 = leggi2D(joinpath(R, "reference/ns_50_$(ev)_wse.csv"), "WSE", "m")
    tq, Q2 = leggi2D(joinpath(R, "reference/ns_50_$(ev)_discharge.csv"), "Q", "m3/s")
    tw == tq || error("$ev: istanti diversi tra livelli e portate 2D")
    k = findall(τ -> τ > 0 && t1[1] <= τ <= t1[end], tw)        # istanti 2D dentro la simulazione 1D
    E2 = E2[k, :]; Q2 = Q2[k, :]
    function interp_t(M, τ)
        i = clamp(searchsortedlast(t1, τ), 1, length(t1) - 1); θ = (τ - t1[i])/(t1[i+1] - t1[i])
        return (1 - θ) .* M[i, :] .+ θ .* M[i+1, :]
    end
    H1c = reduce(vcat, [permutedims(interp_t(H1, τ)) for τ in tw[k]])
    Q1c = reduce(vcat, [permutedims(interp_t(Qs, τ)) for τ in tw[k]])
    m = (E2 .> 1.0) .& (Q2 .> 0.0) .& (abs.(Q1c .- Q2) .<= SOGLIA .* Q2)   # 2D bagnato e stessa portata
    dh = H1c .- (E2 .- zb')
    cnt = vec(sum(m; dims = 1))
    rmse  = [cnt[j] > 0 ? sqrt(sum(abs2, dh[m[:, j], j])/cnt[j]) : NaN for j in 1:n]
    media = [cnt[j] > 0 ? sum(dh[m[:, j], j])/cnt[j] : NaN for j in 1:n]
    qsel  = [cnt[j] > 0 ? extrema(Q2[m[:, j], j]) : (NaN, NaN) for j in 1:n]
    return (; ev, cnt, rmse, media, qsel)
end

eventi = sort(filter(d -> startswith(d, "fittizio_") && isfile(joinpath(BASE, d, "time.csv")), readdir(BASE)))
isempty(eventi) && error("nessuno scenario fittizio_* in $BASE")
length(eventi) <= length(COL) || error("più di $(length(COL)) scenari: aggiungere colori")
S = scenario.(eventi)

# grafico: una linea per scenario, sezioni sulle x
xk = x ./ 1000
tk = findall(nm -> parse(Int, nm[2:end]) % 250 == 0, nomi)       # etichette ogni 250 m
ytop = 100 * maximum(maximum(filter(!isnan, s.rmse)) for s in S)
p = plot(; title = @sprintf("Errore del tirante 1D rispetto al 2D negli istanti con la stessa portata (entro %g%%)", 100SOGLIA),
         xlabel = "sezione", ylabel = "scarto quadratico medio del tirante 1D − 2D [cm]",
         xticks = (xk[tk], nomi[tk]), xrotation = 45, xlims = (0, xk[end] + 0.1), ylims = (-0.07ytop, 1.1ytop),
         legend = :topleft, size = (1500, 650), dpi = 120, left_margin = 6Plots.mm, bottom_margin = 8Plots.mm, st...)
for (s, c) in zip(S, COL)
    plot!(p, xk, 100 .* s.rmse; label = replace(s.ev, "_" => " "), color = c, linewidth = 2,
          markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
end
# nomi delle 6 sezioni peggiori (sopra i punti) e delle 6 migliori (sotto), per la media sugli scenari
hline!(p, [0.0]; color = :black, linewidth = 0.8, label = "")
rm_medio = [sum(s.rmse[j] for s in S)/length(S) for j in 1:n]
peggiori = sort(partialsortperm(rm_medio, 1:6; rev = true))
migliori = sort(partialsortperm(rm_medio, 1:6))
for (sel, sopra) in ((peggiori, true), (migliori, false))
    for j in sel
        y = sopra ? maximum(100 * s.rmse[j] for s in S) + 0.02ytop : -0.01ytop
        # sezioni vicine entrambe etichettate: una a sinistra e una a destra del punto
        al = (j + 1) in sel ? :right : (j - 1) in sel ? :left : :center
        annotate!(p, xk[j] + (al == :right ? -0.01 : al == :left ? 0.01 : 0.0), y,
                  text(nomi[j], 8, sopra ? :black : "#1a7a3a", sopra ? :bottom : :top, al))
        sopra || scatter!(p, [xk[j]], [maximum(100 * s.rmse[j] for s in S)]; label = "", markershape = :utriangle,
                          markersize = 5, color = "#1a7a3a", markerstrokewidth = 0)
    end
end
savefig(p, joinpath(OUT, "14_errore_tirante_sezioni.png"))

# testo: riepilogo per scenario, sezioni peggiori, tabella completa
io = IOBuffer()
@printf(io, "Scarto quadratico medio del tirante 1D - 2D negli istanti 2D bagnati con |Q1D - Q2D| <= %g%% di Q2D\n", 100SOGLIA)
println(io, "Geometria: $GEO\n")
for s in S
    v = filter(!isnan, s.rmse)
    @printf(io, "%s: istanti per sezione %d-%d; scarto quadratico medio mediana %.1f cm, media %.1f cm, massimo %.1f cm in %s\n",
            s.ev, minimum(s.cnt), maximum(s.cnt), 100sort(v)[(end+1)÷2], 100sum(v)/length(v),
            100maximum(v), nomi[argmax(replace(s.rmse, NaN => -1))])
end
println(io, "\nSezioni peggiori (media sugli scenari dello scarto quadratico medio):")
for j in partialsortperm(rm_medio, 1:15; rev = true)
    @printf(io, "  %s (x = %.2f km): %.1f cm; differenza media %s cm\n", nomi[j], xk[j], 100rm_medio[j],
            join([@sprintf("%+.1f", 100s.media[j]) for s in S], " / "))
end
println(io, "\nSezioni migliori (media sugli scenari dello scarto quadratico medio):")
for j in partialsortperm(rm_medio, 1:15)
    @printf(io, "  %s (x = %.2f km): %.1f cm; differenza media %s cm\n", nomi[j], xk[j], 100rm_medio[j],
            join([@sprintf("%+.1f", 100s.media[j]) for s in S], " / "))
end
println(io, "\nPer sezione: per ogni scenario numero di istanti, scarto quadratico medio e differenza media [cm], portate 2D selezionate [m3/s]")
@printf(io, "%-6s %6s", "sez", "x[km]"); foreach(s -> @printf(io, " | %-28s", s.ev), S); println(io)
for j in 1:n
    @printf(io, "%-6s %6.2f", nomi[j], xk[j])
    foreach(s -> @printf(io, " | %3d %5.1f %+5.1f %5.1f-%5.1f", s.cnt[j], 100s.rmse[j], 100s.media[j], s.qsel[j]...), S)
    println(io)
end
testo = String(take!(io))
write(joinpath(OUT, "errore_tirante_sezioni.txt"), testo)
print(split(testo, "\nPer sezione")[1])
println("\nGrafico: ", joinpath(OUT, "14_errore_tirante_sezioni.png"))
