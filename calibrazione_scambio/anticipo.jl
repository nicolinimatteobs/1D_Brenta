# Istante del rientro: quando la portata in S4200 supera 67 m³/s nel 2D e nel 1D (fittizio_4: le due onde,
# divise a 14,5 h), con l'inizio del rientro nel 1D e il volume uscito prima di allora, per km di tratto.
# Uso, dalla cartella del progetto: julia --project=. calibrazione_scambio/anticipo.jl <cartella in out/> [<cartella> ...]
# Scrive anche risultati/anticipo_<cartella>.txt.
using Printf
const DIR = @__DIR__; const PROGETTO = normpath(joinpath(DIR, ".."))
mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(f) if !isempty(strip(l))])
G = [parse.(Float64, split(l, ',')) for l in readlines(joinpath(PROGETTO, "input/geometria/Test_Brenta/S_50m_q2/localslope_xcoord.csv"))[2:end]]
x = getindex.(G, 3); i = findfirst(xx -> round(Int, xx) == 4200, x)
for run in ARGS
    io = IOBuffer()
    println(io, "Istante del rientro, simulazione out/$run: portata in S4200 oltre 67 m³/s nel 2D e nel 1D (anticipo positivo se il 1D è prima).")
    for k in 1:4
        L = readlines(joinpath(PROGETTO, "reference/ns_50_fittizio_$(k)_discharge.csv")); h = strip.(split(L[1], ','))
        M = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
        t2 = M[:, 1]; q2 = M[:, findfirst(==("Q_S4200 [m3/s]"), h)]
        D = joinpath(DIR, "out", run, "fittizio_$k"); t1 = vec(mat(joinpath(D, "time.csv"))); Q1 = mat(joinpath(D, "Q.csv")); L1 = mat(joinpath(D, "Q_lat.csv"))
        q1 = 0.5 .* (Q1[:, i] .+ Q1[:, i+1])
        for (t0, tf) in (k == 4 ? ((0.0, 14.5), (14.5, 99.0)) : ((0.0, 99.0),))
            f(t, q) = (j = findfirst(jj -> t0 <= t[jj]/3600 < tf && q[jj] >= 67, eachindex(t)); j === nothing ? NaN : t[j]/3600)
            a2, a1 = f(t2, q2), f(t1, q1)
            @printf(io, "fittizio_%d%-9s 2D %5.2f h, 1D %5.2f h, anticipo %4.0f min", k, k == 4 ? (t0 == 0 ? " 1ª onda" : " 2ª onda") : "", a2, a1, 60 * (a2 - a1))
            jr = findfirst(j -> t0 <= t1[j]/3600 < tf && sum(max.(L1[j, :], 0.0)) > 0.5, eachindex(t1))
            if jr === nothing
                println(io, "; nessun rientro nel 1D")
                continue
            end
            dt = diff([0.0; t1])
            vout = [sum(-min.(L1[1:jr, c], 0.0) .* dt[1:jr]) for c in 1:size(L1, 2)]
            z = [sum(vout[c] for c in eachindex(x) if a*1000 <= x[c] < (a+1)*1000) for a in 0:5]
            @printf(io, "; rientro nel 1D da %.2f h, uscito prima per km (0-1 … 5-6): %s migliaia di m³\n", t1[jr]/3600, join([@sprintf("%.0f", v/1000) for v in z], " / "))
        end
    end
    s = String(take!(io)); write(joinpath(DIR, "risultati", "anticipo_$run.txt"), s); print(s)
end
