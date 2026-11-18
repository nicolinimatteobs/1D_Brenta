include(joinpath(@__DIR__, "esegui.jl")); include(joinpath(@__DIR__, "indici.jl"))
lista(s) = parse.(Float64, split(s, ','))
const I42 = findfirst(xx -> round(Int, xx) == 4200, X)
ONDE(k) = k == 4 ? [(0.0, 14.5), (14.5, 99.0)] : [(0.0, 99.0)]
primo(t, q, t0, tf, cond) = (j = findfirst(jj -> t0 <= t[jj]/3600 < tf && cond(q[jj]), eachindex(t)); j === nothing ? NaN : t[j]/3600)
function discesa(t, q, t0, tf)       # primo istante dopo il colmo dell'onda con portata sotto 50
    jj = findall(j -> t0 <= t[j]/3600 < tf, eachindex(t)); p = jj[argmax(q[jj])]
    j = findfirst(j -> j > p && t[j]/3600 < tf && q[j] < 50, eachindex(t)); j === nothing ? NaN : t[j]/3600
end
# Prova del serbatoio attivo (riga 15 di settings.txt): per ogni combinazione di T e soglia lancia i fittizi 1-4 e
# aggiunge a risultati/attivo.csv indici, volumi e anticipo del rientro (S4200 oltre 67 m³/s) e della discesa (S6000
# sotto 50 m³/s dopo il colmo), per onda. Le uscite in out/ sono cancellate, salvo le combinazioni "T:V" indicate.
# Uso, dalla cartella del progetto:
#   julia --project=. -O2 calibrazione_scambio/attivo.jl <T [h]> <soglie [m3]> <τ [h]> <T_lento [h]> <attivo 0|1> ["T:V,..." da tenere]
#   es. julia --project=. -O2 calibrazione_scambio/attivo.jl 0.25,0.5,1,2,4,8 125000,150000,175000,200000,250000 1 32 1
const TS, VS = lista(ARGS[1]), lista(ARGS[2]); const τ, TL = parse(Float64, ARGS[3]), parse(Float64, ARGS[4])
const ATT = ARGS[5]
const TIENI = length(ARGS) >= 6 ? split(ARGS[6], ',') : String[]
fout = joinpath(@__DIR__, "risultati", "attivo.csv")
isfile(fout) || write(fout, "attivo,T_h,V_soglia_m3,tau_h,T_lento_h,scenario,onda,colmo1D,colmo2D,NSE,dcolmo,mae_livello,V1_uscito,V1_rientrato,V2_uscito,V2_rientrato,anticipo_rientro_min,anticipo_discesa_min\n")
for T in TS, V in VS
    nome = @sprintf("att%s_T_%g_V_%g", ATT, T, V)
    for k in 1:4
        D = simula(geometria = GEOMETRIA, scenario = k, out = joinpath(@__DIR__, "out", nome, "fittizio_$k"), V = V, T = T, τ = τ, T_lento = TL, attivo = ATT == "1")
        e = indici(D, k); r = REF[k]
        t1 = vec(mat(joinpath(D, "time.csv"))); Q1 = mat(joinpath(D, "Q.csv"))
        q42 = 0.5 .* (Q1[:, I42] .+ Q1[:, I42+1])
        for (n, (t0, tf)) in enumerate(ONDE(k))
            ar = 60 * (primo(r.t2, r.Q2[:, I42], t0, tf, >=(67)) - primo(t1, q42, t0, tf, >=(67)))
            ad = 60 * (discesa(r.t2, r.Q2[:, end], t0, tf) - discesa(t1, Q1[:, end], t0, tf))
            open(fout, "a") do io
                @printf(io, "%s,%g,%g,%g,%g,%d,%d,%.3f,%.3f,%.4f,%.3f,%.4f,%.0f,%.0f,%.0f,%.0f,%.1f,%.1f\n", ATT, T, V, τ, TL, k, n, e.colmo1, e.colmo2,
                        e.nse, e.dcolmo, e.mae, e.v1_out, e.v1_in, e.v2_out, e.v2_in, ar, ad)
            end
            @printf("attivo %s, T %g h, soglia %g: fittizio_%d onda %d: colmo %.2f (2D %.2f), NSE %.3f, anticipo rientro %.0f min, discesa %.0f min\n",
                    ATT, T, V, k, n, e.colmo1, e.colmo2, e.nse, ar, ad)
        end
    end
    "$(T):$(V)" in TIENI || rm(joinpath(@__DIR__, "out", nome); recursive = true)
end
