# Legge settings.txt e restituisce i parametri della simulazione. Le righe di settings.txt
# vengono lette per posizione: non vanno spostate. I percorsi degli idrogrammi sono raccolti
# in hydro_files; la riga 25 è la cartella in cui salvare gli output. Le righe 1, 2 e 26 non
# vengono più lette. Numero di righe e di colonne delle tabelle e numero di sezioni sono
# ricavati direttamente dai file di input
function INITIAL_SETTINGS!()
    # Legge esattamente 30 righe e rimuove spazi
    righe = open("settings.txt", "r") do f
        [strip(readline(f)) for _ in 1:30]
    end
    # Parsing dei dati
    parametri_sezioni = righe[3]
    local_slope       = righe[4]
    hydro_files = (brenta1    = righe[5],  brenta2    = righe[6],
                   brenta3    = righe[7],  brenta4    = righe[8],
                   brenta5    = righe[9],  brenta6    = righe[10],
                   brenta7    = righe[11], brenta_new = righe[12],
                   brenta_v1  = righe[13], brenta_v2  = righe[14],
                   brenta_v3  = righe[15], brenta_v4  = righe[16],
                   caldonazzo = righe[17], ceggio     = righe[18],
                   levico     = righe[19], centa      = righe[20],
                   chieppena  = righe[21], grigno     = righe[22],
                   maso       = righe[23], moggio     = righe[24])
    output_dir        = righe[25]
    Ks                = parse(Float64, righe[27])
    time              = parse(Float64, righe[28])
    dtmax             = parse(Float64, righe[29])
    Nmax              = round(Int, parse(Float64, righe[30]))  # Numero massimo di passi

    # Dimensioni lette dai file: righe e colonne delle tabelle, sezioni dal file della geometria
    n_row, n_col = FILE_SIZE(parametri_sezioni)
    imax, _      = FILE_SIZE(local_slope)

    return Parametri(parametri_sezioni, local_slope, hydro_files, output_dir,
                     Ks, time, dtmax, Nmax, n_row, imax, n_col)
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
