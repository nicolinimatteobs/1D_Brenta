# Parametri della simulazione letti da settings.txt; n_row, n_col e imax sono ricavati dai file di input.
# Gli affluenti vengono dal file della riga 4 di settings.txt: nomi, idrogrammi e sezioni di ingresso, da monte a valle
# (vettori vuoti se gli affluenti sono spenti o se il file è "nessuno")
struct Parametri
    parametri_sezioni  ::String
    local_slope        ::String
    file_monte         ::String          # idrogramma di monte
    aff_nomi           ::Vector{String}  # affluenti: nomi
    aff_file           ::Vector{String}  # affluenti: idrogrammi
    aff_sez            ::Vector{Int}     # affluenti: sezione di ingresso (0 = scelta in automatico, distribuzione uniforme)
    output_dir         ::String
    time               ::Float64
    dtmax              ::Float64
    Nmax               ::Int
    con_affluenti      ::Bool
    scambio            ::Symbol   # scambio laterale: :sfioro, :tabella o :nessuno
    V_soglia           ::Float64  # serbatoio della piana (scambio :tabella): volume trattenuto prima del rientro [m3]
    T_piana            ::Float64  # tempo di svuotamento [s]
    tau_piana          ::Float64  # ritardo di percorrenza [s]
    T_lento            ::Float64  # tempo di rilascio lento della parte sotto la soglia [s] (0 = nessun rilascio)
    piana_attiva       ::Bool     # serbatoio attivo: la parte oltre la soglia rientra anche se la tabella non lo chiede
    n_row              ::Int
    imax               ::Int
    n_col              ::Int
end
