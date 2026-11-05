# Parametri della simulazione letti da settings.txt e dal file della struttura.
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
