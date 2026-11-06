# Parametri della simulazione letti da settings.txt; n_row, n_col e imax sono ricavati dai file di input.
# hydro_files contiene i percorsi degli idrogrammi, es. params.hydro_files.centa
struct Parametri
    parametri_sezioni  ::String
    local_slope        ::String
    hydro_files        ::NamedTuple
    ICdischarge        ::String
    Ks                 ::Float64
    time               ::Float64
    dtmax              ::Float64
    Nmax               ::Int
    n_row              ::Int
    imax               ::Int
    n_col              ::Int
end
