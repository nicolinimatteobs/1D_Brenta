# Estrae i livelli del pelo libero delle sezioni da un file di risultati del 2D (BASEMENT, .h5) e li scrive
# in un csv, nello stesso formato di reference/ns_50_fittizio_N_wse.csv: intestazione
# "t [s],WSE_S0050 [m],...,WSE_S6000 [m]", una riga per uscita, separatore virgola, numeri scritti con la
# conversione standard di Julia (sezioni asciutte: 2.2250738585072014e-308, il valore che BASEMENT mette nel file).
# Dove stanno i dati: RESULTS/NodeStrg/StateVar/<uscita> (matrice 30x124 in Julia: variabili x nodestring);
# il livello è la variabile 1; il tempo [s] è l'attributo "t" di ogni dataset. Le colonne sono associate alle
# sezioni PER NOME, con NodeStrg/StrNames (BC_monte, P17_monte, P17_valle e BC_valle sono scartate).
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/estrai_wse_2D.jl <file .h5> <csv di uscita>
# Esempio:
#   julia --project=. strumenti/estrai_wse_2D.jl ../Test_Brenta/results/inputs/ns_50_fittizio_1_results.h5 reference/ns_50_fittizio_1_wse.csv
using HDF5

if length(ARGS) != 2
    println("Uso: julia --project=. strumenti/estrai_wse_2D.jl <file .h5> <csv di uscita>")
    exit(1)
end
h5file, csvfile = ARGS

h5open(h5file, "r") do h
    nomi = String.(strip.(string.(read(h["NodeStrg/StrNames"]))))     # nome di ogni nodestring (colonna)
    cols = findall(n -> occursin(r"^S\d{4}$", n), nomi)               # solo le sezioni, nell'ordine del file
    isempty(cols) && error("nessuna sezione SXXXX in StrNames")
    g = h["RESULTS/NodeStrg/StateVar"]
    uscite = sort(collect(keys(g)))                                   # 0000000, 0000001, ...
    righe = String[]
    for u in uscite
        d = g[u]
        wse = read(d)[1, :]                                           # livello: variabile 1, tutte le nodestring
        length(wse) == length(nomi) || error("uscita $u: $(length(wse)) nodestring invece di $(length(nomi))")
        push!(righe, join([attrs(d)["t"]; wse[cols]], ","))
    end
    intest = "t [s]," * join(["WSE_" * nomi[c] * " [m]" for c in cols], ",")
    write(csvfile, intest * "\n" * join(righe, "\n") * "\n")
    println("$(length(uscite)) righe, $(length(cols)) sezioni scritte in $csvfile")
end
