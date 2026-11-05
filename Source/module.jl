function INITIAL_SETTINGS!()
    # Legge esattamente 30 righe e rimuove spazi
    righe = open("settings.txt", "r") do f
        [strip(readline(f)) for _ in 1:30]
    end
    # Parsing dei dati
    struttura         = righe[1]
    nr_strut          = parse(Int, righe[2])
    parametri_sezioni = righe[3]
    local_slope       = righe[4]
    hydro_brenta1     = righe[5]
    hydro_brenta2     = righe[6]
    hydro_brenta3     = righe[7]
    hydro_brenta4     = righe[8]
    hydro_brenta5     = righe[9]
    hydro_brenta6     = righe[10]
    hydro_brenta7     = righe[11]
    hydro_brenta_new  = righe[12]
    hydro_brenta_v1   = righe[13]
    hydro_brenta_v2   = righe[14]
    hydro_brenta_v3   = righe[15]
    hydro_brenta_v4   = righe[16]
    hydro_caldonazzo  = righe[17]
    hydro_ceggio      = righe[18]
    hydro_levico      = righe[19]
    hydro_centa       = righe[20]
    hydro_chieppena   = righe[21]
    hydro_grigno      = righe[22]
    hydro_maso        = righe[23]
    hydro_moggio      = righe[24]
    ICdischarge       = righe[25]
    Q_output          = righe[26]
    Ks                = parse(Float64, righe[27])
    time              = parse(Float64, righe[28])
    dtmax             = parse(Float64, righe[29])
    temp              = parse(Float64, righe[30]) #max loop in time temporaneo
    # Leggi il file struttura
    # Leggi solo le prime 4 righe del file struttura
    val_rif = open(strip(struttura), "r") do f
        [parse(Int, strip(readline(f))) for _ in 1:4]
    end 

    n_line_per_section = val_rif[1]
    n_row              = val_rif[2]
    imax               = val_rif[3]
    n_col              = val_rif[4]

    # Numero massimo di passi
    Nmax = round(Int, temp)

    return Parametri(
    struttura, nr_strut, parametri_sezioni, local_slope,
    hydro_brenta1, hydro_brenta2, hydro_brenta3, hydro_brenta4,
    hydro_brenta5, hydro_brenta6, hydro_brenta7, hydro_brenta_new,
    hydro_brenta_v1, hydro_brenta_v2, hydro_brenta_v3, hydro_brenta_v4,
    hydro_caldonazzo, hydro_ceggio, hydro_levico, hydro_centa,
    hydro_chieppena, hydro_grigno, hydro_maso, hydro_moggio,
    ICdischarge, Q_output, Ks, time, dtmax, temp, Nmax,
    n_line_per_section, n_row, imax, n_col)
end

#struttura, nr_strut, parametri_sezioni, local_slope, hydro_Brenta,
#    nr_PC, nc_PC, hydro_INTBAC1, hydro_INTBAC2, hydro_INTBAC3, hydro_INTBAC4, hydro_INTBAC5,
#    hydro_INTBAC6, hydro_INTBAC7, hydro_INTBAC8, nr_aff, nc_aff, ICdischarge, Q_output, Ks, time, dtmax, temp, Nmax, 
#    n_line_per_section, n_row, imax, n_col

function ALLOCATE_ARRAYS!(imax, n_col)
    # VETTORI CON VALORI FISSI
    val       = Vector{Float64}(undef, n_col)
    b         = Vector{Float64}(undef, imax)
    dx_b      = Vector{Float64}(undef, imax-1)
    x_int     = Vector{Float64}(undef, imax+1)
    dx_int    = Vector{Float64}(undef, imax)

    # ARRAYS CON VALORI VARIABILI
    eta       = Vector{Float64}(undef, imax)
    A_bn      = Vector{Float64}(undef, imax)
    deta      = Vector{Float64}(undef, imax)
    rh_b      = Vector{Float64}(undef, imax)
    rh_int    = Vector{Float64}(undef, imax+1)
    A_int     = Vector{Float64}(undef, imax+1)
    A_int_g   = Vector{Float64}(undef, imax+1)
    A_b       = Vector{Float64}(undef, imax)
    dA        = Vector{Float64}(undef, imax)
    q         = Vector{Float64}(undef, imax+1)
    uf        = Vector{Float64}(undef, imax+1)
    qstar     = Vector{Float64}(undef, imax+1)
    γ         = Vector{Float64}(undef, imax+1)
    aa        = Vector{Float64}(undef, imax)
    bb        = Vector{Float64}(undef, imax)
    cc        = Vector{Float64}(undef, imax)
    di        = Vector{Float64}(undef, imax)
    ff        = Vector{Float64}(undef, imax)
    smax      = Vector{Float64}(undef, imax+1)
    β_b       = Vector{Float64}(undef, imax)
    β_int     = Vector{Float64}(undef, imax+1)
    rhs       = Vector{Float64}(undef, imax)
    den_γ_b   = Vector{Float64}(undef, imax)
    den_γ_int = Vector{Float64}(undef, imax+1)
    qn        = Vector{Float64}(undef, imax+1)

    return (; val, b, dx_b, x_int, dx_int,
            eta, A_bn, deta, rh_b, rh_int, A_int, A_int_g, A_b, dA,
            q, uf, qstar, γ, aa, bb, cc, di, ff,
            smax, β_b, β_int, rhs, den_γ_b, den_γ_int, qn)

end

function LOAD_DATA!(data, filename, row, col)

    try
        open(filename, "r") do f
            # Salta la prima riga (header)
            readline(f)

            # Leggi riga per riga
            for n in 1:row
                line = strip(readline(f))
                fields = split(line, ',')
                val = parse.(Float64, fields)
                if length(val) != col
                    error("Errore: numero di colonne nella riga $n non corrisponde a $col")
                end
                data[n, :] = val
            end
        end
    catch e
        println("Errore durante la lettura del file '$filename': ", e)
        return nothing
    end
    return data
end
