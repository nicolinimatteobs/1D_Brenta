# Lancia il modello (Source/ del progetto) senza toccare il settings.txt del progetto: per ogni simulazione
# crea una cartella di lavoro temporanea con il suo settings.txt, che punta agli input del progetto con
# percorsi assoluti. Le righe non indicate vengono dal settings.txt del progetto (idrogrammi degli affluenti,
# Ks, tempo iniziale); passo massimo 60 s, al massimo 20000 passi, affluenti spenti.
# Si usa con include("esegui.jl") negli script di questa cartella.
const PROGETTO = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PROGETTO, "Source", "Run_Simulation.jl")); using .RUN_SIMULATION

"""
    simula(; geometria, scenario, out, scambio = "tabella", V = 0, T = 0, τ = 0, T_lento = 0)

Simulazione del fittizio `scenario` (1-4) con le tabelle di `input/geometria/Test_Brenta/<geometria>/`,
output in `out`. `scambio` è la riga 28 di settings.txt (sfioro, tabella o nessuno); `V` è la soglia del
serbatoio della piana [m³], `T` il tempo di svuotamento [h], `τ` il ritardo [h] e `T_lento` il tempo di
rilascio lento [h] (righe 29-32).
"""
function simula(; geometria, scenario, out, scambio = "tabella", V = 0, T = 0, τ = 0, T_lento = 0)
    base = [strip(l) for l in readlines(joinpath(PROGETTO, "settings.txt"))[1:32]]
    r = copy(base)
    ass(p) = isabspath(p) ? p : joinpath(PROGETTO, p)
    r[1] = ass("input/geometria/Test_Brenta/$geometria/tabella_sezioni.csv")
    r[2] = ass("input/geometria/Test_Brenta/$geometria/localslope_xcoord.csv")
    for j in [3:13; 15:21]; r[j] = ass(r[j]); end
    r[14] = ass("input/idrogrammi/fittizio_$scenario/piena.csv")
    r[22] = abspath(out); r[25] = "60"; r[26] = "20000"; r[27] = "no"
    r[28] = scambio; r[29] = string(V); r[30] = string(T); r[31] = string(τ); r[32] = string(T_lento)
    dir = mktempdir()
    write(joinpath(dir, "settings.txt"), join(r, "\r\n") * "\r\n")
    cd(() -> RUN_SIMULATION.SEMI_IMPLICIT_BRENTA_2D_1D!(), dir)
    rm(dir; recursive = true)
    return abspath(out)
end
