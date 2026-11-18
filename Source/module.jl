# Legge settings.txt e restituisce i parametri della simulazione. Le prime 15 righe di settings.txt
# vengono lette per posizione: non vanno spostate. Riga 1 tabelle delle sezioni, riga 2 geometria,
# riga 3 idrogramma di monte, riga 4 file degli affluenti (o "nessuno"), riga 5 cartella degli output,
# righe 6-8 tempo iniziale, passo massimo e numero massimo di passi, riga 9 affluenti attivi (si/no),
# riga 10 scambio laterale con la piana (sfioro, tabella o nessuno), righe 11-13 parametri del serbatoio
# della piana per lo scambio 'tabella': soglia [m3], tempo di svuotamento [h] e ritardo [h], riga 14 tempo
# di rilascio lento della parte sotto la soglia [h], riga 15 serbatoio attivo (si/no: la parte oltre la
# soglia rientra anche quando la tabella non chiede rientri). Numero di righe e di colonne delle tabelle,
# numero di sezioni e numero di affluenti sono ricavati direttamente dai file di input
function INITIAL_SETTINGS!()
    # Legge le prime 15 righe e rimuove spazi
    righe = open("settings.txt", "r") do f
        [strip(readline(f)) for _ in 1:15]
    end
    # Parsing dei dati
    parametri_sezioni = righe[1]
    local_slope       = righe[2]
    file_monte        = righe[3]
    file_affluenti    = righe[4]
    output_dir        = righe[5]
    time              = parse(Float64, righe[6])
    dtmax             = parse(Float64, righe[7])
    Nmax              = round(Int, parse(Float64, righe[8]))  # Numero massimo di passi
    con_affluenti     = PARSE_SI_NO(righe[9], "riga 9 di settings.txt (affluenti attivi)")
    scambio           = PARSE_SCAMBIO(righe[10])
    V_soglia          = PARSE_NON_NEGATIVO(righe[11], "riga 11 di settings.txt (soglia del serbatoio della piana, m3)")
    T_piana           = 3600.0 * PARSE_NON_NEGATIVO(righe[12], "riga 12 di settings.txt (tempo di svuotamento, ore)")
    tau_piana         = 3600.0 * PARSE_NON_NEGATIVO(righe[13], "riga 13 di settings.txt (ritardo, ore)")
    T_lento           = 3600.0 * PARSE_NON_NEGATIVO(righe[14], "riga 14 di settings.txt (rilascio lento, ore)")
    piana_attiva      = PARSE_SI_NO(righe[15], "riga 15 di settings.txt (serbatoio attivo)")

    # Dimensioni lette dai file: righe e colonne delle tabelle, sezioni dal file della geometria
    n_row, n_col = FILE_SIZE(parametri_sezioni)
    imax, _      = FILE_SIZE(local_slope)
    scambio == :tabella && n_col < 14 &&
        error("Riga 10 di settings.txt: lo scambio 'tabella' richiede le colonne ΔQ e ΔQ_s in $parametri_sezioni")
    isfile(file_monte) || error("Riga 3 di settings.txt: idrogramma di monte '$file_monte' non trovato")

    # Affluenti dal file della riga 4, letto solo se sono attivi
    aff_nomi, aff_file, aff_sez = con_affluenti ? READ_TRIBUTARIES(file_affluenti, imax) : (String[], String[], Int[])

    return Parametri(parametri_sezioni, local_slope, file_monte, aff_nomi, aff_file, aff_sez, output_dir,
                     time, dtmax, Nmax, con_affluenti, scambio, V_soglia, T_piana, tau_piana, T_lento,
                     piana_attiva, n_row, imax, n_col)
end

# Legge il file degli affluenti (riga 4 di settings.txt): CSV con l'intestazione nome,idrogramma,sezione e una
# riga per affluente, da monte a valle. idrogramma è il percorso del file (come quelli di settings.txt);
# sezione è il numero della sezione di sbocco (1 … imax) e può restare vuota: allora la sezione è scelta in
# automatico, distribuendo gli affluenti uniformemente lungo il tratto. "nessuno" vuol dire nessun affluente.
# Ogni errore indica la riga del file
function READ_TRIBUTARIES(filename, imax)
    nomi, file, sez = String[], String[], Int[]
    lowercase(filename) in ("nessuno", "") && return nomi, file, sez
    isfile(filename) || error("Riga 4 di settings.txt: file degli affluenti '$filename' non trovato")
    for (n, line) in enumerate(eachline(filename))
        (n == 1 || isempty(strip(line))) && continue   # intestazione e righe vuote
        c = String.(strip.(split(line, ',')))
        length(c) == 2 && push!(c, "")
        length(c) == 3 || error("Riga $n di '$filename': servono 3 colonne (nome,idrogramma,sezione), trovate $(length(c))")
        isempty(c[1]) && error("Riga $n di '$filename': manca il nome dell'affluente")
        isfile(c[2]) || error("Riga $n di '$filename': idrogramma '$(c[2])' dell'affluente $(c[1]) non trovato")
        s = isempty(c[3]) ? 0 : tryparse(Int, c[3])
        (s === nothing || (!isempty(c[3]) && !(1 <= s <= imax))) &&
            error("Riga $n di '$filename': sezione '$(c[3])' dell'affluente $(c[1]) non valida, serve un numero da 1 a $imax o niente")
        push!(nomi, c[1]); push!(file, c[2]); push!(sez, s)
    end
    return nomi, file, sez
end

# Legge un numero non negativo da settings.txt
function PARSE_NON_NEGATIVO(s, dove)
    v = tryparse(Float64, strip(s))
    (v === nothing || v < 0) && error("Valore '$s' non valido per $dove: serve un numero non negativo")
    return v
end

# Legge lo scambio laterale con la piana (riga 10 di settings.txt): sfioro (sfioro libero sugli
# argini, solo uscente), tabella (ΔQ delle tabelle, entrante e uscente) o nessuno
function PARSE_SCAMBIO(s)
    v = lowercase(strip(s))
    v in ("sfioro", "tabella", "nessuno") && return Symbol(v)
    error("Valore '$s' non valido per la riga 10 di settings.txt (scambio laterale): usare sfioro, tabella o nessuno")
end

# Legge un sì/no da settings.txt: accetta si, sì, no, true, false, 1, 0 (maiuscole o minuscole)
function PARSE_SI_NO(s, dove)
    v = lowercase(strip(s))
    v in ("si", "sì", "true", "1") && return true
    v in ("no", "false", "0") && return false
    error("Valore '$s' non valido per $dove: usare si o no")
end

# Conta le righe di dati di un file separato da virgole (esclusa l'intestazione e le righe vuote)
# e le colonne dell'intestazione
function FILE_SIZE(filename)
    open(filename, "r") do f
        n_col = length(split(strip(readline(f)), ','))
        n_row = count(l -> !isempty(strip(l)), eachline(f))
        return n_row, n_col
    end
end

# Crea tutti gli array del modello, inizializzati a zero, e li restituisce per nome.
# Quelli di lunghezza imax stanno sulle sezioni, quelli di lunghezza imax+1 sulle interfacce
function ALLOCATE_ARRAYS!(imax)
    # VETTORI CON VALORI FISSI
    b         = zeros(imax)
    dx_b      = zeros(imax-1)
    x_int     = zeros(imax+1)
    dx_int    = zeros(imax)

    # ARRAYS CON VALORI VARIABILI
    eta       = zeros(imax)
    A_bn      = zeros(imax)
    deta      = zeros(imax)
    rh_b      = zeros(imax)
    rh_int    = zeros(imax+1)
    A_int     = zeros(imax+1)
    A_int_g   = zeros(imax+1)
    A_b       = zeros(imax)
    dA        = zeros(imax)
    q         = zeros(imax+1)
    uf        = zeros(imax+1)
    qstar     = zeros(imax+1)
    γ         = zeros(imax+1)
    aa        = zeros(imax)
    bb        = zeros(imax)
    cc        = zeros(imax)
    di        = zeros(imax)
    ff        = zeros(imax)
    smax      = zeros(imax+1)
    β_b       = zeros(imax)
    β_int     = zeros(imax+1)
    rhs       = zeros(imax)
    den_γ_b   = zeros(imax)
    den_γ_int = zeros(imax+1)
    qn        = zeros(imax+1)

    return (; b, dx_b, x_int, dx_int,
            eta, A_bn, deta, rh_b, rh_int, A_int, A_int_g, A_b, dA,
            q, uf, qstar, γ, aa, bb, cc, di, ff,
            smax, β_b, β_int, rhs, den_γ_b, den_γ_int, qn)

end

# Legge un file di dati separati da virgole e riempie la matrice data (row x col), saltando
# la prima riga di intestazione. Si ferma con un errore se una riga non ha col colonne
function LOAD_DATA!(data, filename, row, col)

    open(filename, "r") do f
        # Salta la prima riga (header)
        readline(f)

        # Leggi riga per riga
        for n in 1:row
            line = strip(readline(f))
            fields = split(line, ',')
            val = parse.(Float64, fields)
            if length(val) != col
                error("Errore: la riga $n del file '$filename' non ha $col colonne")
            end
            data[n, :] = val
        end
    end
    return data
end
