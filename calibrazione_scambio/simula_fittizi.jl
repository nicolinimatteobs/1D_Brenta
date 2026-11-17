# Lancia i fittizi 1-4 con le tabelle S_50m_q2 e uno scambio laterale, con le uscite in out/<nome>/fittizio_N:
# servono per gli idrogrammi e l'anticipo del rientro (sfioro come termine di paragone, combinazioni scelte).
# Uso, dalla cartella del progetto:
#   julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl <nome> <sfioro|tabella|nessuno> [soglia [m3] T [h] τ [h] [T_lento [h] [attivo si|no]]]
#   es. julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl sfioro sfioro
#       julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl tarato_lento tabella 150000 0.5 1 32
#       julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl attivo tabella 150000 0.5 1 32 si
include(joinpath(@__DIR__, "esegui.jl"))
const NOME, SC = ARGS[1], ARGS[2]
const V, T, τ = length(ARGS) >= 5 ? parse.(Float64, ARGS[3:5]) : (0.0, 0.0, 0.0)
const TL = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 0.0
const ATT = length(ARGS) >= 7 && lowercase(ARGS[7]) in ("si", "sì", "1", "true")
for k in 1:4
    simula(geometria = "S_50m_q2", scenario = k, out = joinpath(@__DIR__, "out", NOME, "fittizio_$k"), scambio = SC, V = V, T = T, τ = τ, T_lento = TL, attivo = ATT)
end
