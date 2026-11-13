using Plots, Plots.Measures, Printf, CSV, DataFrames, Infiltrator, HDF5

# Caso da elaborare (lanciare dalla cartella 1D_Brenta). La cartella di input contiene i file delle
# sezioni (FilePrefix*nome*Format), centerpoints.csv, bankpoints.csv, la mesh del modello 2D (MeshFile)
# e i suoi risultati (ResultsFile), da cui si calcola β.
# β calcolato dal modello 2D (true) oppure β = 1 (false); con β = 1 mesh e risultati 2D non servono.
# Il modello usa le tabelle con β = 1 (input/geometria/Test_Brenta/S_50m): con β dal 2D i risultati
# cambiano poco, quindi quelle tabelle sono tenute da parte (input/geometria/Test_Brenta/S_50m_beta).
const β_dal_2D = false
# Sezioni ogni 50 m
const InputFolder = "calcolo_tabelle_daBM/input/S_0050_6000"
const OutputFolder = β_dal_2D ? "calcolo_tabelle_daBM/output/S_50m_beta" : "calcolo_tabelle_daBM/output/S_50m"
# Cartella degli input del modello 1D: tabella_sezioni.csv e localslope_xcoord.csv vengono scritti anche qui,
# così il modello usa sempre le tabelle appena calcolate (righe 1-2 di settings.txt)
const ModelFolder = β_dal_2D ? "input/geometria/Test_Brenta/S_50m_beta" : "input/geometria/Test_Brenta/S_50m"
const FilePrefix = "ns_50_scalini_"
const MeshFile = "test_brenta_ns_50_computational-mesh.2dm"
const ResultsFile = "ns_50_scalini_results.h5"
# Il caso a 200 m (input/S_0200_6000, prefisso ns_channel_) non ha mesh e risultati 2D: con il calcolo
# di β non si può più elaborare.
const Format = ".csv"
const g = 9.81      # [m/s^2]

# Limiti dei valori ritenuti ragionevoli, usati solo per i controlli
const iWSE_lim = (0.0005, 0.005)    # [-], la pendenza media del fondo è circa 0.002
const Ks_lim = (10.0, 50.0)         # [m^(1/3)/s]
const β_lim = (1.0, 1.5)            # [-]
const Qrel_lim = 0.10               # scarto relativo massimo tra la portata degli elementi e quella della sezione

# Valori usati per ricostruire le righe asciutte
const h_dry = 0.01      # [m], tirante
const Q_dry = 0.1       # [m3/s], portata
const Ks_dry = 20.0     # [m^(1/3)/s], Ks imposto

#= CONTROLLO DEI FILE DI INPUT
    Prima di elaborare controlla che ci siano tutti i file: quelli fissi della cartella di input (mesh e
    risultati 2D solo se β viene dal 2D) e quelli delle sezioni elencate in centerpoints.csv. Se ne manca
    qualcuno li elenca tutti e si ferma.
=#
function check_input_files()
    mancanti = String[]
    fissi = β_dal_2D ? ("centerpoints.csv", "bankpoints.csv", MeshFile, ResultsFile) : ("centerpoints.csv", "bankpoints.csv")
    for f in fissi
        isfile(joinpath(InputFolder, f)) || push!(mancanti, joinpath(InputFolder, f))
    end
    sezioni = joinpath(InputFolder, "centerpoints.csv")
    if isfile(sezioni)
        for nome in String.(CSV.read(sezioni, DataFrame)[!, "nome"])
            f = joinpath(InputFolder, FilePrefix*nome*Format)
            isfile(f) || push!(mancanti, f)
        end
    end
    isempty(mancanti) && return nothing
    msg = "AVVISO: mancano $(length(mancanti)) file di input (la cartella di lavoro deve essere 1D_Brenta):\n  " *
          join(mancanti, "\n  ")
    printstyled(msg, "\n"; color=:red, bold=true)
    error("File di input mancanti: elaborazione interrotta")
end

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
    Larghezza, pendenza e β sono quelli del primo livello bagnato.
=#
function rebuild_dry!(data, wet, iWSE, β)
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
    β[begin:k1-1] .= β[k1]
    return data, iWSE, β
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

#= LETTURA DELLA MESH (.2dm)
    Restituisce le coordinate dei nodi, i triangoli (nodi in ordine, numerati come nel file) e le
    nodestring (elenco dei nodi da un estremo all'altro). Nel file una nodestring è una riga
    "NS n1 n2 ... -nk nome": l'ultimo nodo ha il segno meno ed è seguito dal nome.
=#
function read_mesh(path)
    nodes = Dict{Int,NTuple{2,Float64}}()
    tri = NTuple{3,Int}[]
    strings = Dict{String,Vector{Int}}()
    cur = Int[]
    for line in eachline(path)
        t = split(line)
        isempty(t) && continue
        if t[1] == "ND"
            nodes[parse(Int, t[2])] = (parse(Float64, t[3]), parse(Float64, t[4]))
        elseif t[1] == "E3T"
            parse(Int, t[2]) == length(tri) + 1 || error("$path: triangoli non numerati in ordine")
            push!(tri, (parse(Int, t[3]), parse(Int, t[4]), parse(Int, t[5])))
        elseif t[1] == "NS"
            for v in t[2:end]
                n = tryparse(Int, v)
                if n === nothing
                    strings[v] = cur
                    cur = Int[]
                else
                    push!(cur, abs(n))
                end
            end
        end
    end
    return nodes, tri, strings
end

#= ELEMENTI A MONTE DI UNA NODESTRING
    Ogni tratto della nodestring è un lato della mesh, condiviso da due triangoli: si prende quello
    dalla parte di monte. Nel modello 2D ogni nodestring ha "upstream_direction" (right o left), cioè
    da che parte sta monte percorrendo la nodestring dal primo all'ultimo nodo. Per ogni tratto
    restituisce il triangolo, la lunghezza del lato e la normale unitaria verso valle.
=#
function upstream_elements(ns, nodes, tri, edges, σ)
    elementi = Tuple{Int,Float64,NTuple{2,Float64}}[]
    for (p, q) in zip(ns[begin:end-1], ns[begin+1:end])
        (xp, yp), (xq, yq) = nodes[p], nodes[q]
        dx, dy = xq - xp, yq - yp
        L = hypot(dx, dy)
        lato = get(edges, minmax(p, q), Int[])
        # σ = +1 se monte è a destra: il baricentro del triangolo di monte ha prodotto vettoriale negativo
        monte = filter(lato) do e
            xc = sum(nodes[k][1] for k in tri[e]) / 3
            yc = sum(nodes[k][2] for k in tri[e]) / 3
            σ * (dx * (yc - yp) - dy * (xc - xp)) < 0
        end
        length(monte) == 1 || error("Nodestring: nessun triangolo (o più di uno) a monte del lato $p-$q")
        push!(elementi, (monte[1], L, (-σ * dy / L, σ * dx / L)))
    end
    return elementi
end

#= COEFFICIENTE DI QUANTITÀ DI MOTO β DAL MODELLO 2D
    Come nel metodo di Engelund (eq. 13 delle slide):
        β = Ω/Q² ∫ u² Y dy
    ma con le velocità del modello 2D al posto di quelle di moto uniforme striscia per striscia.
    Le strisce sono i triangoli appena a monte dei lati della nodestring: tirante Y = quota del pelo
    libero - quota del fondo dell'elemento, portata specifica normale alla sezione qn = q·n, larghezza
    ΔL = lunghezza del lato. Quindi u² Y = qn²/Y e
        Ω = Σ Y ΔL,   Q = Σ qn ΔL,   β = Ω Σ(qn²/Y ΔL) / Q²
    Ω coincide con l'area bagnata dei file delle sezioni: BASEMENT la calcola con gli stessi elementi.
    Q invece no, perché BASEMENT la calcola dai flussi sui lati e non dai valori nelle celle
    (in mediana gli elementi danno il 2% in più). Per β si usa la Q degli elementi, così con
    velocità uniforme vale esattamente 1 e non è mai minore di 1.
    Gli elementi con tirante non superiore al minimo del modello 2D sono asciutti (portata nulla).
    Restituisce β, una riga per scalino e una colonna per sezione (NaN se la sezione è asciutta),
    e gli avvisi sulla portata degli elementi.
=#
function beta_2D(nomi, dati, wet_all)
    nodes, tri, strings = read_mesh(joinpath(InputFolder, MeshFile))
    edges = Dict{Tuple{Int,Int},Vector{Int}}()
    for (e, t) in enumerate(tri), (p, q) in ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
        push!(get!(edges, minmax(p, q), Int[]), e)
    end
    avvisi = String[]
    h5open(joinpath(InputFolder, ResultsFile)) do f
        # i triangoli devono essere le celle del file dei risultati (che numera i nodi da 0)
        topo = read(f["CellsAll/Topology"])
        size(topo, 2) == length(tri) && all(topo[:, e] == collect(tri[e]) .- 1 for e in eachindex(tri)) ||
            error("La mesh $MeshFile non corrisponde a quella di $ResultsFile")
        zb = vec(read(f["CellsAll/BottomEl"]))
        hmin = only(read(f["Parameters/MinWaterDepth"]))
        model = read(f[".config/model"])
        ns_h5 = read(f["NodeStrg/StrNames"])
        uscite = sort(keys(f["RESULTS/CellsAll/HydState"]))
        nlevels = nrow(dati[first(nomi)][1])
        length(uscite) == nlevels || error("$ResultsFile ha $(length(uscite)) uscite, i file delle sezioni $nlevels righe")
        H  = [read(f["RESULTS/CellsAll/HydState/"*u]) for u in uscite]   # quota, qx, qy per cella
        SV = [read(f["RESULTS/NodeStrg/StateVar/"*u]) for u in uscite]   # come i file delle sezioni
        β = fill(NaN, nlevels, length(nomi))
        for (j, nome) in enumerate(nomi)
            haskey(strings, nome) || error("La nodestring $nome non è in $MeshFile")
            js = findfirst(==(nome), ns_h5)
            js === nothing && error("La nodestring $nome non è in $ResultsFile")
            m = match(Regex("\"name\":\"$nome\",\"upstream_direction\":\"(right|left)\""), model)
            m === nothing && error("Direzione di monte della nodestring $nome non trovata in $ResultsFile")
            elementi = upstream_elements(strings[nome], nodes, tri, edges, m[1] == "right" ? 1 : -1)
            data = dati[nome][1]
            for k in 1:nlevels
                # le righe del file della sezione devono essere le uscite del file dei risultati (portata,
                # area e, se bagnata, livello: in quelle asciutte load_section lo ha cambiato)
                stessa = isapprox(SV[k][2, js], data[k, "Discharge [m3/s]"]; rtol=1e-6, atol=1e-6) &&
                         isapprox(SV[k][3, js], data[k, "Wetted area [m2]"]; rtol=1e-6, atol=1e-6) &&
                         (!wet_all[k, j] || isapprox(SV[k][1, js], data[k, "Mean wse [m]"]; rtol=1e-9))
                stessa || error("$nome riga $k: il file della sezione non corrisponde all'uscita $(uscite[k]) di $ResultsFile")
                wet_all[k, j] || continue
                Ω = Qe = M = 0.0
                for (e, L, (nx, ny)) in elementi
                    Y = H[k][1, e] - zb[e]
                    Y > hmin || continue
                    qn = H[k][2, e] * nx + H[k][3, e] * ny
                    Ω += Y * L
                    Qe += qn * L
                    M += qn^2 / Y * L
                end
                # tolleranza per gli arrotondamenti di BASEMENT: a 50 m gli scarti arrivano a 1,4e-4, un
                # elemento sbagliato cambierebbe l'area di qualche percento
                isapprox(Ω, data[k, "Wetted area [m2]"]; rtol=1e-3, atol=1e-3) ||
                    error(@sprintf("%s riga %d: area degli elementi %.4f diversa da quella della sezione %.4f",
                                   nome, k, Ω, data[k, "Wetted area [m2]"]))
                Qe > 0.0 || error("$nome riga $k: portata degli elementi non positiva")
                Qs = data[k, "Discharge [m3/s]"]
                abs(Qe / Qs - 1) <= Qrel_lim || push!(avvisi, @sprintf("%-6s Q=%7.3f  portata degli elementi %7.3f (%+.1f%%)",
                                                                       nome, Qs, Qe, 100 * (Qe / Qs - 1)))
                β[k, j] = Ω * M / Qe^2
            end
        end
        return β, avvisi
    end
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
function section_table(data, iWSE, β)
    Q   = data[!, "Discharge [m3/s]"]
    WSE = data[!, "Mean wse [m]"]

    # compila l'area
    Ω   = data[!, "Wetted area [m2]"]
    Ω_s = similar(Ω)
    slope!(Ω_s, Ω, WSE)

    # compila il coefficiente di quantità di moto: dal modello 2D (beta_2D) o 1, secondo β_dal_2D
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
    WSE, iWSE, Ks, β = table[!, "WSE"], table[!, "iWSE"], table[!, "Ks"], table[!, "β"]
    for k in 2:nrow(table)
        Q[k] > Q[k-1] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q non crescente: %.3f -> %.3f", nome, WSE[k], Q[k-1], Q[k]))
        # β vale almeno 1: sotto, di poco, solo per arrotondamento (un solo elemento bagnato)
        β_lim[1] - 1e-9 <= β[k] <= β_lim[2] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  β=%.3f", nome, WSE[k], Q[k], β[k]))
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
        (:β,    "beta",                "Coefficiente β",         "β [-]"),
        (:jden, "jden",                "jden",                   "jden [m]"),
    ]
    β_dal_2D || filter!(gr -> gr[1] != :β, grafici)    # con β = 1 il grafico di β non serve
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
    check_input_files()

    # leggi il file riguardante tutte le sezioni, ordinandole da monte verso valle
    InputFile = "centerpoints.csv"
    InputPath = joinpath(InputFolder, InputFile)
    sections_all = sort(CSV.read(InputPath, DataFrame), "s_m")

    # leggi i dati delle sezioni del canale principale, cioè quelle elencate in centerpoints.csv
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

    # calcola β per ogni sezione e per ogni WSE dai risultati del modello 2D, oppure β = 1
    if β_dal_2D
        β_all, avvisi_β = beta_2D(nomi, dati, wet_all)
        append!(avvisi, avvisi_β)
    else
        β_all = ones(size(WSE_all))
    end

    # geometria delle sezioni: la pendenza del fondo è calcolata come quella della superficie libera,
    # ma con il segno usato dal modello 1D, negativo quando il fondo scende verso valle
    s = sections_all[!, "s_m"]
    talweg = [dati[nome][1][begin, "Reference elevation [m]"] for nome in nomi]
    local_slope = -vec(free_surface_slope(permutedims(talweg), trues(1, length(nomi)), s))
    # l'ultima sezione è in contropendenza, ma il modello 1D ne usa la pendenza per la portata di
    # uscita: uso la pendenza media del fondo sull'intero tratto
    local_slope[end] = (talweg[end] - talweg[begin]) / (s[end] - s[begin])

    # ciclo sulle sezioni, da monte verso valle
    tabelle = Dict{String,DataFrame}()
    @printf("%-10s %6s %17s %17s %17s\n", "Sezione", "n WSE", "iWSE [-]", "Ks [m^(1/3)/s]", "β [-]")
    for (j, nome) in enumerate(nomi)
        data, wet = dati[nome]
        iWSE = iWSE_all[:, j]
        β = β_all[:, j]

        # sostituisci le righe asciutte con valori ricostruiti
        rebuild_dry!(data, wet, iWSE, β)

        # il modello 1D vuole livelli crescenti: scarta le righe con livello non superiore al massimo
        # dei precedenti (es. S0050 ogni 50 m, dove negli ultimi scalini il livello scende)
        WSE = data[!, "Mean wse [m]"]
        tieni = [k == 1 || WSE[k] > maximum(WSE[1:k-1]) for k in eachindex(WSE)]
        for k in findall(.!tieni)
            push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  livello non crescente: riga scartata", nome, WSE[k], data[k, "Discharge [m3/s]"]))
        end
        data = data[tieni, :]
        iWSE = iWSE[tieni]
        β = β[tieni]

        #= ULTIMA SEZIONE: CONDIZIONE DI USCITA
            Il modello 1D calcola la portata in uscita come moto uniforme, Q = Ω·√(g·jden·S), con il jden
            della tabella dell'ultima sezione e la sua pendenza locale S (la media del fondo, vedi sopra).
            Qui quindi Ks e jden sono ricavati con quella pendenza invece che con la pendenza del pelo
            libero: così, a ogni livello della tabella, il modello dà in uscita la portata degli scalini 2D.
            La pendenza del pelo libero nel 2D tra le ultime due sezioni è molto più ripida (fino a 0,0075
            contro 0,00213), probabilmente per il richiamo dell'uscita del dominio 2D: con quella e il Ks 33
            di settings.txt il tirante in uscita era 10-13 cm più alto del 2D.
        =#
        if j == length(nomi)
            iWSE = fill(abs(local_slope[end]), length(iWSE))
        end

        table = section_table(data, iWSE, β)
        tabelle[nome] = table

        i_ok = table[!, "iWSE"]
        K_ok = table[!, "Ks"]
        @printf("%-10s %6d %8.5f %8.5f %8.2f %8.2f %8.3f %8.3f\n", nome, nrow(table),
            minimum(i_ok), maximum(i_ok), minimum(K_ok), maximum(K_ok), minimum(β), maximum(β))
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
        mkpath(ModelFolder)
        CSV.write(joinpath(OutputFolder, "tabella_sezioni"*Format), tabella_unica)
        CSV.write(joinpath(ModelFolder, "tabella_sezioni"*Format), tabella_unica)
    end

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
        CSV.write(joinpath(ModelFolder, "localslope_xcoord"*Format), geometria)
        println("Tabelle e geometria scritte in $OutputFolder e in $ModelFolder")
    end

    println("\nControlli: $(length(avvisi)) valori fuori dai limiti")
    foreach(println, avvisi)

    @infiltrate false

    return tabelle, tabella_unica, geometria
end

main(; dosave=true)
