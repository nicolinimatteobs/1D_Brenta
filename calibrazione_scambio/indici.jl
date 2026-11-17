# Funzioni comuni: lettura dei risultati 1D e dei riferimenti 2D dei fittizi, indici di confronto e volumi.
# Si usa con include("indici.jl") negli script di questa cartella (dopo esegui.jl, che definisce PROGETTO).
using Printf, Statistics

const GEOMETRIA = "S_50m_q2"                         # tabelle con la colonna ΔQ
const T0 = 3 * 3600.0                                 # livelli confrontati dalle 3 h (prima il 2D è ancora asciutto)
const ONDA = Dict(1 => 0.0, 2 => 0.0, 3 => 0.0, 4 => 14.5 * 3600)   # inizio dell'onda per il colmo (fittizio_4: la seconda)

csvnum(f) = [parse.(Float64, split(l, ',')) for l in readlines(f)[2:end] if !isempty(strip(l))]
mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(f) if !isempty(strip(l))])
nomesez(xi) = @sprintf("S%04d", round(Int, xi/50)*50)
trapz(t, y) = sum((t[2:end] .- t[1:end-1]) .* (y[2:end] .+ y[1:end-1]) ./ 2)
interp_t(t1, Y, τ) = (m = clamp(searchsortedlast(t1, τ), 1, length(t1) - 1); θ = (τ - t1[m]) / (t1[m+1] - t1[m]); (1 - θ) .* Y[m, :] .+ θ .* Y[m+1, :])

const X = getindex.(csvnum(joinpath(PROGETTO, "input/geometria/Test_Brenta/$GEOMETRIA/localslope_xcoord.csv")), 3)

# colonne del riferimento 2D (portata o livello) nelle sezioni del 1D
function leggi2D(f, pref)
    L = readlines(f); hdr = strip.(split(L[1], ','))
    M = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
    return M[:, 1], M[:, [findfirst(==("$(pref)_$(nomesez(xi)) [" * (pref == "WSE" ? "m" : "m3/s") * "]"), hdr) for xi in X]]
end

const REF = Dict(k => begin
    t2, Q2 = leggi2D(joinpath(PROGETTO, "reference/ns_50_fittizio_$(k)_discharge.csv"), "Q")
    _, E2 = leggi2D(joinpath(PROGETTO, "reference/ns_50_fittizio_$(k)_wse.csv"), "WSE")
    # volume scambiato in ogni tratto tra due sezioni in tutto l'evento, ∫(Q_i − Q_(i−1)) dt
    v2 = [trapz(t2, Q2[:, i] .- Q2[:, i-1]) for i in 2:length(X)]
    (; t2, Q2, E2, v2)
end for k in 1:4)

"""
    indici(D, k)

Confronto della simulazione nella cartella `D` (fittizio `k`) con il 2D: colmo in S6000 (fittizio_4: seconda
onda), NSE della portata in S6000, errore medio del colmo lungo il tratto, errore medio e medio assoluto del
livello allo stesso istante (sezioni 2D bagnate, dalle 3 h), volumi usciti e rientrati del 1D (lordi, passo
per passo, e netti per tratto) e del 2D (netti per tratto), volume entrato.
"""
function indici(D, k)
    r = REF[k]; imax = length(X)
    t1 = vec(mat(joinpath(D, "time.csv"))); Q1 = mat(joinpath(D, "Q.csv")); E1 = mat(joinpath(D, "eta.csv")); L1 = mat(joinpath(D, "Q_lat.csv"))
    Qs = 0.5 .* (Q1[:, 1:end-1] .+ Q1[:, 2:end]); Qs[:, end] .= Q1[:, end]
    j1 = findall(>=(ONDA[k]), t1); j2 = findall(>=(ONDA[k]), r.t2)
    c1 = vec(maximum(Qs[j1, :]; dims = 1)); c2 = vec(maximum(r.Q2[j2, :]; dims = 1))
    jj = findall(j -> r.Q2[j, end] > 0 && r.t2[j] <= t1[end], eachindex(r.t2))
    q1 = [interp_t(t1, Q1[:, end:end], r.t2[j])[1] for j in jj]; q2 = r.Q2[jj, end]
    nse = 1 - sum(abs2, q1 .- q2) / sum(abs2, q2 .- mean(q2))
    de = Float64[]
    for j in findall(τ -> T0 <= τ <= t1[end], r.t2)
        e1 = interp_t(t1, E1, r.t2[j])
        for i in 1:imax
            (r.E2[j, i] > 1 && r.Q2[j, i] > 0) && push!(de, e1[i] - r.E2[j, i])
        end
    end
    dt = diff([0.0; t1])
    v_out = sum(dt .* vec(sum(min.(L1, 0.0); dims = 2))); v_in = sum(dt .* vec(sum(max.(L1, 0.0); dims = 2)))
    v1 = vec(sum(dt .* L1; dims = 1))[2:end]
    return (colmo1 = c1[end], colmo2 = c2[end], nse, dcolmo = mean(abs.(c1 .- c2)), mae = mean(abs.(de)), bias = mean(de),
            v_out = -v_out, v_in, v1_out = -sum(min.(v1, 0.0)), v1_in = sum(max.(v1, 0.0)),
            v2_out = -sum(min.(r.v2, 0.0)), v2_in = sum(max.(r.v2, 0.0)), v_entrato = trapz(t1, Q1[:, 1]))
end
