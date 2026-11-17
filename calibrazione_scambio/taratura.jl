# Taratura del serbatoio della piana (scambio 'tabella', tabelle S_50m_q2): per ogni combinazione di ritardo τ,
# soglia V e tempo di svuotamento T lancia i fittizi 1-4, li confronta con il 2D e aggiunge una riga per scenario a
# risultati/taratura.csv. Le uscite vanno in out/tau_<τ>_V_<V>_T_<T>/fittizio_N e si tengono solo per le
# combinazioni indicate nel quarto argomento, le altre si cancellano dopo il calcolo degli indici.
# Uso, dalla cartella del progetto:
#   julia --project=. -O2 calibrazione_scambio/taratura.jl <τ [h]> <soglie [m3]> <T [h]> ["τ:V:T,τ:V:T"]
#   es. julia --project=. -O2 calibrazione_scambio/taratura.jl 0.5,1 100000,150000 0.5,1 1:150000:0.5
include(joinpath(@__DIR__, "esegui.jl")); include(joinpath(@__DIR__, "indici.jl"))
lista(s) = parse.(Float64, split(s, ','))
const TAUS, VS, TS = lista(ARGS[1]), lista(ARGS[2]), lista(ARGS[3])
const TIENI = length(ARGS) >= 4 ? [lista(replace(c, ':' => ',')) for c in split(ARGS[4], ',')] : Vector{Float64}[]

fout = joinpath(@__DIR__, "risultati", "taratura.csv"); mkpath(dirname(fout))
isfile(fout) || write(fout, "tau_h,V_soglia_m3,T_h,scenario,colmo1D,colmo2D,NSE,dcolmo,mae_livello,bias_livello,V_uscito,V_rientrato," *
                            "V1_uscito_tratti,V1_rientrato_tratti,V2_uscito_tratti,V2_rientrato_tratti,V_entrato\n")
for τ in TAUS, V in VS, T in TS
    nome = @sprintf("tau_%g_V_%g_T_%g", τ, V, T)
    for k in 1:4
        D = simula(geometria = GEOMETRIA, scenario = k, out = joinpath(@__DIR__, "out", nome, "fittizio_$k"), V = V, T = T, τ = τ)
        e = indici(D, k)
        open(fout, "a") do io
            @printf(io, "%g,%g,%g,%d,%.3f,%.3f,%.4f,%.3f,%.4f,%.4f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f\n", τ, V, T, k, e.colmo1, e.colmo2, e.nse,
                    e.dcolmo, e.mae, e.bias, e.v_out, e.v_in, e.v1_out, e.v1_in, e.v2_out, e.v2_in, e.v_entrato)
        end
        @printf("τ %g h, soglia %g m³, T %g h, fittizio_%d: colmo %.2f (2D %.2f), NSE %.3f, |Δcolmo| %.2f, livello %.3f m\n",
                τ, V, T, k, e.colmo1, e.colmo2, e.nse, e.dcolmo, e.mae)
    end
    [τ, V, T] in TIENI || rm(joinpath(@__DIR__, "out", nome); recursive = true)
end
