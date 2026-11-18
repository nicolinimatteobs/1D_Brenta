# Lancia il modello (Source/ del progetto) senza toccare il settings.txt del progetto: per ogni simulazione
# crea una cartella di lavoro temporanea con il suo settings.txt, che punta agli input del progetto con
# percorsi assoluti. Le righe non indicate vengono dal settings.txt del progetto (file degli affluenti,
# tempo iniziale); passo massimo 60 s, al massimo 20000 passi, affluenti spenti.
# Si usa con include("esegui.jl") negli script di questa cartella.
const PROGETTO = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PROGETTO, "Source", "Run_Simulation.jl")); using .RUN_SIMULATION

"""
    simula(; geometria, scenario, out, scambio = "tabella", V = 0, T = 0, τ = 0, T_lento = 0, attivo = false)

Simulazione del fittizio `scenario` (1-4) con le tabelle di `input/geometria/Test_Brenta/<geometria>/`,
output in `out`. `scambio` è la riga 10 di settings.txt (sfioro, tabella o nessuno); `V` è la soglia del
serbatoio della piana [m³], `T` il tempo di svuotamento [h], `τ` il ritardo [h] e `T_lento` il tempo di
rilascio lento [h] (righe 11-14); `attivo` è il serbatoio attivo (riga 15).
"""
function simula(; geometria, scenario, out, scambio = "tabella", V = 0, T = 0, τ = 0, T_lento = 0, attivo = false)
    base = [strip(l) for l in readlines(joinpath(PROGETTO, "settings.txt"))[1:15]]
    r = copy(base)
    ass(p) = isabspath(p) ? p : joinpath(PROGETTO, p)
    r[1] = ass("input/geometria/Test_Brenta/$geometria/tabella_sezioni.csv")
    r[2] = ass("input/geometria/Test_Brenta/$geometria/localslope_xcoord.csv")
    r[3] = ass("input/idrogrammi/fittizio_$scenario/piena.csv")
    lowercase(r[4]) == "nessuno" || (r[4] = ass(r[4]))
    r[5] = abspath(out); r[7] = "60"; r[8] = "20000"; r[9] = "no"
    r[10] = scambio; r[11] = string(V); r[12] = string(T); r[13] = string(τ); r[14] = string(T_lento); r[15] = attivo ? "si" : "no"
    dir = mktempdir()
    write(joinpath(dir, "settings.txt"), join(r, "\r\n") * "\r\n")
    cd(() -> RUN_SIMULATION.SEMI_IMPLICIT_BRENTA_2D_1D!(), dir)
    rm(dir; recursive = true)
    return abspath(out)
end
