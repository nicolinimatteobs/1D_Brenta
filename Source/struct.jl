# Parametri della simulazione letti da settings.txt; n_row, n_col e imax sono ricavati dai file di input.
# hydro_files contiene i percorsi degli idrogrammi, es. params.hydro_files.centa
struct Parametri
    parametri_sezioni  ::String
    local_slope        ::String
    hydro_files        ::NamedTuple
    output_dir         ::String
    Ks                 ::Float64
    time               ::Float64
    dtmax              ::Float64
    Nmax               ::Int
    con_affluenti      ::Bool
    scambio            ::Symbol   # scambio laterale: :sfioro, :tabella o :nessuno
    V_soglia           ::Float64  # serbatoio della piana (scambio :tabella): volume trattenuto prima del rientro [m3]
    T_piana            ::Float64  # tempo di svuotamento [s]
    tau_piana          ::Float64  # ritardo di percorrenza [s]
    T_lento            ::Float64  # tempo di rilascio lento della parte sotto la soglia [s] (0 = nessun rilascio)
    n_row              ::Int
    imax               ::Int
    n_col              ::Int
end
