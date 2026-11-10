using Plots, Plots.Measures, Printf, CSV, DataFrames, Infiltrator

# Caso da elaborare (lanciare dalla cartella 1D_Brenta). La cartella di input contiene i file delle
# sezioni (FilePrefix*nome*Format), sections_dx.csv e bankpoints.csv
# Sezioni ogni 50 m
const InputFolder = "calcolo_tabelle_daBM/input/S_0050_6000"
const OutputFolder = "calcolo_tabelle_daBM/output/S_50m"
const FilePrefix = "ns_50_scalini_"
# Sezioni ogni 200 m
# const InputFolder = "calcolo_tabelle_daBM/input/S_0200_6000"
# const OutputFolder = "calcolo_tabelle_daBM/output/S_200m"
# const FilePrefix = "ns_channel_"
const Format = ".csv"
const g = 9.81      # [m/s^2]

# Limiti dei valori ritenuti ragionevoli, usati solo per i controlli
const iWSE_lim = (0.0005, 0.005)    # [-], la pendenza media del fondo è circa 0.002
const Ks_lim = (10.0, 50.0)         # [m^(1/3)/s]

# Valori usati per ricostruire le righe asciutte
const h_dry = 0.01      # [m], tirante
const Q_dry = 0.1       # [m3/s], portata
const Ks_dry = 20.0     # [m^(1/3)/s], Ks imposto

function slope!(slope, y, x)
    @. slope[begin:end-1] = (y[begin+1:end] - y[begin:end-1]) / (x[begin+1:end] - x[begin:end-1])
    slope[end] = slope[end-1]
end

# Legge i dati di una sezione e "ripara" i valori artificialmente messi a 0 quando è asciutta
function load_section(InputPath)
    data = CSV.read(InputPath, DataFrame)
    wet = data[!, "Wetted area [m2]"] .> 0.0
    zref = data[2, "Reference elevation [m]"]
    data[.!wet, "Reference elevation [m]"] .= zref
    data[.!wet, "Mean bottom elevation [m]"] .= zref
    data[.!wet, "Mean wse [m]"] .= zref
    return data, wet
end

#= RICOSTRUZIONE DELLE RIGHE ASCIUTTE
    La simulazione 2D non può usare portate e livelli troppo bassi per questioni di stabilità,
    quindi i primi scalini possono essere asciutti. Il modello 1D serve per le piene e non vede
    mai portate sotto 1 m3/s, quindi basta che questi valori siano ragionevoli.
    Prima riga: tirante h_dry, portata Q_dry, larghezza e pendenza del primo livello bagnato.
    L'area viene da Manning imponendo Ks_dry, con RH = Ω/W:
        Q = Ks * Ω * (Ω/W)^(2/3) * sqrt(i)  =>  Ω = (Q * W^(2/3) / (Ks * sqrt(i)))^(3/5)
    L'area non è W*h_dry, altrimenti RH sarebbe h_dry e il Ks sarebbe irragionevole.
    Eventuali altre righe asciutte sono interpolate fino al primo livello bagnato.
=#
function rebuild_dry!(data, wet, iWSE)
    k1 = findfirst(wet)
    all(wet[k1:end]) || error("Sono presenti righe asciutte dopo il primo livello bagnato")
    WSE, Q, Ω, W = data[!, "Mean wse [m]"], data[!, "Discharge [m3/s]"], data[!, "Wetted area [m2]"], data[!, "Wetted geometric length [m]"]
    zref = data[begin, "Reference elevation [m]"]
    WSE[begin] = zref + h_dry
    Q[begin]   = Q_dry
    Ω[begin]   = (Q_dry * cbrt(W[k1])^2 / (Ks_dry * sqrt(iWSE[k1])))^(3/5)
    for k in 2:k1-1
        θ = (k - 1) / (k1 - 1)
        WSE[k] = WSE[begin] + θ * (WSE[k1] - WSE[begin])
        Q[k]   = Q[begin] + θ * (Q[k1] - Q[begin])
        Ω[k]   = Ω[begin] + θ * (Ω[k1] - Ω[begin])
    end
    W[begin:k1-1] .= W[k1]
    iWSE[begin:k1-1] .= iWSE[k1]
    return data, iWSE
end

#= CALCOLO DELLA PENDENZA DELLA SUPERFICIE LIBERA
    WSE e wet hanno una riga per ogni scalino di portata e una colonna per ogni sezione,
    con le sezioni ordinate da monte verso valle secondo l'ascissa s. A parità di riga i
    livelli appartengono allo stesso istante, quindi formano il profilo di superficie libera.
    La pendenza è positiva quando il livello scende verso valle ed è calcolata:
    - alle differenze centrate, se le sezioni adiacenti sono entrambe bagnate;
    - con la differenza unilaterale, alle estremità o se una delle adiacenti è asciutta;
    - NaN, se la sezione è asciutta o se lo sono entrambe le adiacenti.
=#
function free_surface_slope(WSE, wet, s)
    nlevels, nsections = size(WSE)
    iWSE = fill(NaN, nlevels, nsections)
    for j in 1:nsections, k in 1:nlevels
        wet[k, j] || continue
        up = (j > 1 && wet[k, j-1]) ? j - 1 : j
        dn = (j < nsections && wet[k, j+1]) ? j + 1 : j
        up == dn && continue
        iWSE[k, j] = (WSE[k, up] - WSE[k, dn]) / (s[dn] - s[up])
    end
    return iWSE
end

#= PENDENZE NON POSITIVE
    Con portate molto basse i livelli di sezioni vicine possono essere quasi uguali e la pendenza
    può risultare nulla o negativa (es. S5950 a 0,86 m3/s con le sezioni ogni 50 m). Il Ks e la
    ricostruzione delle righe asciutte non sarebbero definiti, quindi la pendenza viene sostituita
    con quella del livello bagnato successivo della stessa sezione che sia positiva.
    Restituisce l'elenco delle sostituzioni, da segnalare tra i controlli.
=#
function fix_slope!(iWSE, wet, nomi)
    avvisi = String[]
    nlevels, nsections = size(iWSE)
    for j in 1:nsections, k in 1:nlevels
        (wet[k, j] && !(iWSE[k, j] > 0.0)) || continue
        kn = findfirst(m -> wet[m, j] && iWSE[m, j] > 0.0, k+1:nlevels)
        kn === nothing && error("Sezione $(nomi[j]): nessuna pendenza positiva sopra la riga $k")
        push!(avvisi, @sprintf("%-6s riga %d  iWSE=%.6f non positiva, sostituita con %.6f (riga %d)",
            nomi[j], k, iWSE[k, j], iWSE[k+kn, j], k+kn))
        iWSE[k, j] = iWSE[k+kn, j]
    end
    return avvisi
end

#= CALCOLO DI KS
    Non è banale. La pendenza della superficie libera va calcolata per ogni sezione e per
    ogni valore di WSE. Per semplicità assumo che sia unica per tutte le portate, cioè
    che a ogni WSE corrisponda la sola portata del relativo scalino.
    Il Ks comunque deve cambiare a seconda delle diverse regioni bagnate.
    Valori normali di Ks indicano un approccio fisicamente consistente.
    Se la pendenza non è definita o non è positiva il Ks non è definito (NaN).
    La portata è quella locale della sezione, perché quella che resta nel canale
    principale dopo l'esondazione.
=#
function Ks_equivalent(Q, Ω, RH, iWSE)
    (iWSE > 0.0 && Ω > 0.0) || return NaN
    return sqrt( Q^2 / ((Ω * cbrt(RH)^2 )^2 * iWSE) )
end

# Calcola le caratteristiche geometriche/idrauliche di una sezione, una riga per ogni WSE
function section_table(data, iWSE)
    Q   = data[!, "Discharge [m3/s]"]
    WSE = data[!, "Mean wse [m]"]

    # compila l'area
    Ω   = data[!, "Wetted area [m2]"]
    Ω_s = similar(Ω)
    slope!(Ω_s, Ω, WSE)

    # compila il coefficiente per la correzione della velocità, per ora sempre 1
    β   = ones(length(WSE))
    β_s = similar(β)
    slope!(β_s, β, WSE)

    # compila la larghezza
    W   = data[!, "Wetted geometric length [m]"]
    W_s = similar(W)
    slope!(W_s, W, WSE)

    # compila il raggio idraulico
    RH   = Ω ./ W
    RH_s = similar(RH)
    slope!(RH_s, RH, WSE)

    # calcola il Ks equivalente dell'intera sezione
    Ks = Ks_equivalent.(Q, Ω, RH, iWSE)

    # denominatore del termine di attrito del modello 1D, γ = |u|/jden
    jden   = @. Ks^2 * cbrt(RH)^4 / g
    jden_s = similar(jden)
    slope!(jden_s, jden, WSE)

    return DataFrame(; WSE, Ω, Ω_s, β, β_s, W, W_s, RH, RH_s, iWSE, Ks, jden, jden_s)
end

#= Prepara la tabella di una sezione nel formato letto da READ_VALUES! nel modello 1D: 12 colonne
    fisse. In seconda colonna c'è il livello, non il tirante.
=#
function model_table(table, sezione)
    colonne = ["WSE", "Ω", "Ω_s", "β", "β_s", "jden", "jden_s", "W", "W_s", "RH", "RH_s"]
    return insertcols(table[!, colonne], 1, :sezione => sezione)
end

#= CONTROLLI SUI RISULTATI
    La simulazione prevede un'esondazione e le sezioni riguardano solo il canale principale,
    quindi la portata può variare lungo il tratto: se ne terrà conto nel modello 1D con un
    termine di portata laterale uscente. Per questo la portata di una sezione è confrontata
    solo con i propri scalini precedenti, e deve essere crescente.
    I valori fuori dai limiti sono solo segnalati, non modificati.
=#
function check_section(nome, table, Q)
    avvisi = String[]
    WSE, iWSE, Ks = table[!, "WSE"], table[!, "iWSE"], table[!, "Ks"]
    for k in 2:nrow(table)
        Q[k] > Q[k-1] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q non crescente: %.3f -> %.3f", nome, WSE[k], Q[k-1], Q[k]))
        iWSE_lim[1] <= iWSE[k] <= iWSE_lim[2] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  iWSE=%.5f", nome, WSE[k], Q[k], iWSE[k]))
        Ks_lim[1] <= Ks[k] <= Ks_lim[2] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  Ks=%.2f", nome, WSE[k], Q[k], Ks[k]))
    end
    return avvisi
end

# Crea i grafici di una sezione, mostrandoli e/o salvandoli nella cartella della sezione
function plot_section(table, nome, Folder; dosave::Bool=false, doplot::Bool=false)
    grafici = [
        (:Ω,    "area_bagnata",        "Area bagnata",           "Area [m^2]"),
        (:W,    "larghezza_bagnata",   "Larghezza bagnata",      "Larghezza [m]"),
        (:RH,   "raggio.idraulico",    "Raggio idraulico",       "Raggio [m]"),
        (:iWSE, "pendenza_wse",        "Pendenza sup. libera",   "Pendenza [-]"),
        (:Ks,   "ks",                  "Ks equivalente",         "Ks [m^(1/3)/s]"),
        (:jden, "jden",                "jden",                   "jden [m]"),
    ]
    for (colonna, file, titolo, xlabel) in grafici
        plt = plot(table[!, colonna], table[!, "WSE"];
            title="$titolo - $nome",
            xlabel=xlabel,
            ylabel="Livello idrometrico [m]",
            legend=false,
            size=(450, 600)
        )
        doplot && display(plt)
        dosave && savefig(plt, joinpath(Folder, file*".svg"))
    end
    return nothing
end

function main(; dosave::Bool=false, doplot::Bool=false)
    # leggi il file riguardante tutte le sezioni, ordinandole da monte verso valle
    InputFile = "sections_dx.csv"
    InputPath = joinpath(InputFolder, InputFile)
    sections_all = sort(CSV.read(InputPath, DataFrame), "s_m")

    # leggi i dati delle sezioni del canale principale, cioè quelle elencate in sections_dx.csv
    nomi = String.(sections_all[!, "nome"])
    dati = Dict(nome => load_section(joinpath(InputFolder, FilePrefix*nome*Format)) for nome in nomi)
    # le righe di sezioni diverse devono essere gli stessi scalini: stesso numero di righe e, se c'è
    # la colonna del tempo, stessi istanti
    ref = dati[first(nomi)][1]
    for nome in nomi
        d = dati[nome][1]
        nrow(d) == nrow(ref) || error("La sezione $nome ha un numero di righe diverso dalle altre")
        "t [s]" in names(d) && (d[!, "t [s]"] == ref[!, "t [s]"] || error("La sezione $nome ha istanti diversi dalle altre"))
    end

    # calcola la pendenza della superficie libera per ogni sezione e per ogni WSE
    WSE_all = reduce(hcat, [dati[nome][1][!, "Mean wse [m]"] for nome in nomi])
    wet_all = reduce(hcat, [dati[nome][2] for nome in nomi])
    iWSE_all = free_surface_slope(WSE_all, wet_all, sections_all[!, "s_m"])
    avvisi = fix_slope!(iWSE_all, wet_all, nomi)

    # ciclo sulle sezioni, da monte verso valle
    tabelle = Dict{String,DataFrame}()
    @printf("%-10s %6s %17s %17s\n", "Sezione", "n WSE", "iWSE [-]", "Ks [m^(1/3)/s]")
    for (j, nome) in enumerate(nomi)
        data, wet = dati[nome]
        iWSE = iWSE_all[:, j]

        # sostituisci le righe asciutte con valori ricostruiti
        rebuild_dry!(data, wet, iWSE)

        # il modello 1D vuole livelli crescenti: scarta le righe con livello non superiore al massimo
        # dei precedenti (es. S0050 ogni 50 m, dove negli ultimi scalini il livello scende)
        WSE = data[!, "Mean wse [m]"]
        tieni = [k == 1 || WSE[k] > maximum(WSE[1:k-1]) for k in eachindex(WSE)]
        for k in findall(.!tieni)
            push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  livello non crescente: riga scartata", nome, WSE[k], data[k, "Discharge [m3/s]"]))
        end
        data = data[tieni, :]
        iWSE = iWSE[tieni]

        table = section_table(data, iWSE)
        tabelle[nome] = table

        i_ok = table[!, "iWSE"]
        K_ok = table[!, "Ks"]
        @printf("%-10s %6d %8.5f %8.5f %8.2f %8.2f\n", nome, nrow(table),
            minimum(i_ok), maximum(i_ok), minimum(K_ok), maximum(K_ok))
        append!(avvisi, check_section(nome, table, data[!, "Discharge [m3/s]"]))

        # salva tabella e grafici solo se la tabella è diversa da quella già salvata
        Folder = joinpath(OutputFolder, nome)
        TablePath = joinpath(Folder, "tabella"*Format)
        salva = dosave && !(isfile(TablePath) && isequal(CSV.read(TablePath, DataFrame), table))
        if salva
            mkpath(Folder)
            CSV.write(TablePath, table)
        end
        (salva || doplot) && plot_section(table, nome, Folder; dosave=salva, doplot)
    end

    # tabella unica, con il numero progressivo della sezione da monte verso valle
    tabella_unica = reduce(vcat, [model_table(tabelle[nome], j) for (j, nome) in enumerate(nomi)])
    if dosave
        mkpath(OutputFolder)
        CSV.write(joinpath(OutputFolder, "tabella_sezioni"*Format), tabella_unica)
    end

    # geometria delle sezioni: la pendenza del fondo è calcolata come quella della superficie libera,
    # ma con il segno usato dal modello 1D, negativo quando il fondo scende verso valle
    s = sections_all[!, "s_m"]
    talweg = [dati[nome][1][begin, "Reference elevation [m]"] for nome in nomi]
    local_slope = -vec(free_surface_slope(permutedims(talweg), trues(1, length(nomi)), s))
    # l'ultima sezione è in contropendenza, ma il modello 1D ne usa la pendenza per la velocità
    # di uscita: uso la pendenza media del fondo sull'intero tratto
    local_slope[end] = (talweg[end] - talweg[begin]) / (s[end] - s[begin])
    # quote degli argini destro e sinistro di ogni sezione, dai punti in bankpoints.csv
    bankpoints = CSV.read(joinpath(InputFolder, "bankpoints"*Format), DataFrame)
    z_argine = Dict(String(r["nome"]) => r["z_dtm"] for r in eachrow(bankpoints))
    argine_dx = [z_argine[nome*"_dx"] for nome in nomi]
    argine_sx = [z_argine[nome*"_sx"] for nome in nomi]
    geometria = DataFrame("section" => eachindex(nomi), "local_slope" => local_slope,
        "x_coord_[m]" => s, "talweg_[m]" => talweg,
        "argine_dx_[m]" => argine_dx, "argine_sx_[m]" => argine_sx)
    if dosave
        CSV.write(joinpath(OutputFolder, "localslope_xcoord"*Format), geometria)
    end

    println("\nControlli: $(length(avvisi)) valori fuori dai limiti")
    foreach(println, avvisi)

    @infiltrate false

    return tabelle, tabella_unica, geometria
end

main(; dosave=true)
