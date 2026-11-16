# Grafico 14 con argine di S2950 corretto e riferimenti 2D nuovi: stato precedente, senza ponte e con il ponte P17.
# Fatto una tantum. Uso, dalla cartella del progetto:
#   julia --project=. output/S_50m_hq/riferimenti_nuovi/confronto_g14.jl <errore_tirante_sezioni.txt precedente>
using Plots, Printf
const R = pwd()
const OUT = joinpath(R, "output/S_50m_hq/riferimenti_nuovi")
const PREC = ARGS[1]
const GRIGIO = "#8a8f98"; const BLU = "#2a78d6"; const ARANCIO = "#eb6834"

csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
G = csvnum(joinpath(R, "input/geometria/Test_Brenta/S_50m/localslope_xcoord.csv"))
x = getindex.(G, 3); n = length(x); xk = x ./ 1000
nomi = [@sprintf("S%04d", round(Int, xi/50)*50) for xi in x]

# per sezione e scenario: scarto quadratico medio e differenza media [cm]
function leggi_txt(f)
    rms = Vector{Float64}[]; dm = Vector{Float64}[]
    for l in eachline(f)
        occursin(r"^S\d{4}\s", l) || continue
        v = [parse.(Float64, split(s)[2:3]) for s in split(l, '|')[2:end]]
        push!(rms, getindex.(v, 1)); push!(dm, getindex.(v, 2))
    end
    media(M) = [sum(r)/length(r) for r in M]
    return (; rms = media(rms), dm = media(dm), rmsS = rms, dmS = dm)
end
A = leggi_txt(PREC)
E = leggi_txt(joinpath(R, "output/S_50m/grafici/errore_tirante_sezioni.txt"))
F = leggi_txt(joinpath(OUT, "grafici/errore_tirante_sezioni.txt"))
@assert length(A.rms) == n && length(E.rms) == n && length(F.rms) == n

default(fontfamily = "sans-serif", titlefontsize = 12, guidefontsize = 11, tickfontsize = 8, legendfontsize = 9,
        linewidth = 2, markersize = 3, markerstrokewidth = 0, framestyle = :box, grid = true, gridalpha = 0.15,
        foreground_color_legend = nothing, background_color_legend = RGBA(1, 1, 1, 0.85))
tk = findall(nm -> parse(Int, nm[2:end]) % 250 == 0, nomi)
function pannello(titolo, ylab, serie, ylims; etichette = Int[], zero = false)
    p = plot(; title = titolo, ylabel = ylab, xticks = (xk[tk], nomi[tk]), xrotation = 45,
             xlims = (0, 6.1), ylims = ylims, legend = :topleft, left_margin = 8Plots.mm, bottom_margin = 6Plots.mm)
    vspan!(p, [2.575, 3.025]; color = :gray, alpha = 0.12, linewidth = 0, label = "tratto S2600-S3000")
    vline!(p, [3.025]; color = :gray30, linestyle = :dot, linewidth = 1.5, label = "ponte P17 (tra S3000 e S3050)")
    zero && hline!(p, [0.0]; color = :black, linewidth = 1, label = "")
    for (lab, v, c, lw) in serie
        plot!(p, xk, v; label = lab, color = c, linewidth = lw, marker = :circle)
    end
    liv = 0                                    # etichette vicine: alzate a gradini per non sovrapporle
    for (i, j) in enumerate(etichette)
        liv = (i > 1 && xk[j] - xk[etichette[i-1]] < 0.3) ? mod(liv + 1, 3) : 0
        annotate!(p, xk[j], F.rms[j] + 0.35 + 0.7liv, text(nomi[j], 8, :black, :bottom))
    end
    return p
end
peggiori(v, k) = sort(partialsortperm(v, 1:k; rev = true))
p1 = pannello("Scarto quadratico medio del tirante 1D - 2D a parità di portata (entro l'1%), media sui quattro scenari",
              "scarto quadratico medio [cm]",
              [("stato precedente (tabelle, argine e riferimenti vecchi)", A.rms, GRIGIO, 1.5),
               ("senza ponte", E.rms, BLU, 2), ("con il ponte P17", F.rms, ARANCIO, 2)], (0, 15);
              etichette = peggiori(F.rms, 8))
p2 = pannello("Differenza media del tirante 1D - 2D a parità di portata (negativo: 1D più basso)",
              "differenza media [cm]",
              [("senza ponte", E.dm, BLU, 2), ("con il ponte P17", F.dm, ARANCIO, 2)], (-15, 12); zero = true)
xlabel!(p2, "sezione")
savefig(plot(p1, p2; layout = (2, 1), size = (1600, 1050), dpi = 130), joinpath(OUT, "14_confronto_ponte.png"))

open(joinpath(OUT, "confronto_grafico_14.txt"), "w") do io
    println(io, "Grafico 14: scarto quadratico medio del tirante 1D - 2D a parità di portata (entro l'1%) [cm], media sui quattro scenari;")
    println(io, "tra parentesi la differenza media (negativo: 1D più basso)")
    println(io, "P = stato precedente (tabelle, argine di S2950 e riferimenti 2D vecchi, senza ponte, output/S_50m del commit HEAD)")
    println(io, "S = senza ponte (output/S_50m, tabelle e argine nuovi, riferimenti 2D nuovi)")
    println(io, "C = con il ponte P17 (output/S_50m_hq/riferimenti_nuovi, stessi riferimenti)")
    println(io)
    @printf(io, "%-6s %14s %14s %14s\n", "sez", "P", "S", "C")
    for j in 1:n
        @printf(io, "%-6s", nomi[j])
        for d in (A, E, F); @printf(io, " %6.1f (%+5.1f)", d.rms[j], d.dm[j]); end
        println(io)
    end
    println(io)
    for (a, b, nm) in ((0.05, 2.35, "S0050-S2350"), (2.4, 2.55, "S2400-S2550"), (2.6, 3.0, "S2600-S3000"), (3.05, 3.5, "S3050-S3500"), (3.55, 6.0, "S3550-S6000"), (0.05, 6.0, "tutte le sezioni"))
        kk = findall(xi -> a - 1e-6 <= xi <= b + 1e-6, xk)
        @printf(io, "media su %-17s", nm * ":")
        for (l, d) in (("P", A), ("S", E), ("C", F)); @printf(io, "  %s %5.2f cm", l, sum(d.rms[kk])/length(kk)); end
        println(io)
    end
    println(io)
    println(io, "Sezioni S2400-S3050, quattro scenari separati: scarto quadratico medio (differenza media) [cm]")
    println(io, "        |          senza ponte (S)                         |          con il ponte (C)")
    @printf(io, "%-6s |", "sez")
    for blocco in 1:2, s in 1:4; @printf(io, " %12s", "fittizio_$s"); blocco == 1 && s == 4 && print(io, " |"); end
    println(io)
    for j in findall(xi -> 2.4 - 1e-6 <= xi <= 3.05 + 1e-6, xk)
        @printf(io, "%-6s |", nomi[j])
        for (d, sep) in ((E, " |"), (F, ""))
            for s in 1:4; @printf(io, " %5.1f (%+5.1f)", d.rmsS[j][s], d.dmS[j][s]); end
            print(io, sep)
        end
        println(io)
    end
end
print(read(joinpath(OUT, "confronto_grafico_14.txt"), String))
