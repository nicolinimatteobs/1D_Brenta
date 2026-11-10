# Controllo a secco degli input indicati in settings.txt, senza lanciare la simulazione: legge tabelle,
# geometria e idrogrammi come il modello, verifica la coerenza e stima stato iniziale e livelli al picco.
# Uso, dalla cartella del progetto:  julia --project=. strumenti/controllo_input.jl
const R = dirname(@__DIR__)
cd(R)
module C
    using Printf, Parameters, UnPack
    include(joinpath(dirname(@__DIR__), "Source", "struct.jl"))
    include(joinpath(dirname(@__DIR__), "Source", "module.jl"))
    include(joinpath(dirname(@__DIR__), "Source", "subroutines.jl"))
end
using Printf
p = C.INITIAL_SETTINGS!()
println("tabelle:   ", p.parametri_sezioni, "  (righe ", p.n_row, ", colonne ", p.n_col, ")")
println("geometria: ", p.local_slope, "  (sezioni ", p.imax, ")")
println("output:    ", p.output_dir, "  | Ks ", p.Ks, ", dt massimo ", p.dtmax, " s")
g = C.LOAD_GEOMETRY(p.local_slope, p.imax)
data = permutedims(C.LOAD_SECTION_TABLES(p.parametri_sezioni, p.n_row, p.n_col))
fr, lr = C.SECTION_INDEX(data); ys = copy(data[2, :])
length(lr) == p.imax || error("le tabelle hanno $(length(lr)) sezioni, la geometria $(p.imax)")
all(all(diff(ys[fr[k]:fr[k+1]-1]) .> 0) for k in 1:p.imax) || error("quote non crescenti in qualche tabella")
hydro = C.LOAD_ALL_HYDROGRAPHS!(p); tend = C.HYDROGRAPH_END_TIME(hydro)
Qm = hydro[:caldonazzo].Q
@printf("idrogramma di monte: Q iniziale %.2f, massimo %.2f m3/s; durata simulata %.0f s\n", Qm[1], maximum(Qm), tend)
@printf("pendenze locali: da %.5f a %.5f (|S| minimo %.5f)\n", minimum(g.bx), maximum(g.bx), minimum(abs.(g.bx)))
@printf("talweg e prima quota delle tabelle: scarto massimo %.3f m\n", maximum(abs.(g.b_min .- ys[fr[1:end-1]])))
y0 = [C.UNIFORM_LEVEL(data, ys, fr, lr, i, Qm[1], g.bx[i]; y_min = 0.10) - g.b_min[i] for i in 1:p.imax]
yp = [C.UNIFORM_LEVEL(data, ys, fr, lr, i, maximum(Qm), g.bx[i]; y_min = 0.10) for i in 1:p.imax]
fuori = findall(yp[i] >= ys[fr[i+1]-1] - 1e-6 for i in 1:p.imax)
@printf("tirante iniziale (moto uniforme per %.2f m3/s): da %.3f a %.3f m\n", Qm[1], minimum(y0), maximum(y0))
@printf("altezza delle tabelle sopra il talweg: da %.2f a %.2f m\n", minimum(ys[fr[2:end].-1] .- g.b_min), maximum(ys[fr[2:end].-1] .- g.b_min))
println("sezioni in cui il moto uniforme al picco supera la tabella (valori estrapolati): ", isempty(fuori) ? "nessuna" : fuori)
