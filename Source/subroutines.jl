function SPATIAL_DOMAIN!(imax, x_b, dx_b, x_int, dx_int)

    @inbounds for i in 1:imax-1
        dx_b[i] = x_b[i+1] - x_b[i]
    end

    dx_int[1] = dx_b[1]
    
    @inbounds  for i in 2:imax-1
        dx_int[i] = 0.5 * (dx_b[i] + dx_b[i-1])
    end
    dx_int[imax] = dx_b[imax-1]

    @inbounds  for i in 2:imax
        x_int[i] = 0.5 * (x_b[i] + x_b[i-1])
    end

    x_int[1]      = -x_int[2]
    x_int[imax+1] = x_b[imax] + (x_b[imax] - x_int[imax])

    return dx_b, dx_int, x_int
end



# Prima riga di ogni sezione e ultima riga usata. data è la tabella trasposta:
# ogni colonna è una riga di parametri_sezioni, il numero di righe può variare tra le sezioni
function SECTION_INDEX(data)
    n = size(data, 2)
    first_row = [1; [r for r in 2:n if data[1,r] != data[1,r-1]]; n+1]
    return first_row, first_row[1:end-1]
end

# Riga r della sezione con ys[r] <= y < ys[r+1], con dh anche non costante.
# Prova prima l'ultima riga usata (g) e le vicine, altrimenti usa la ricerca binaria
@inline function FIND_ROW(ys, lo, hi, y, g)
    @inbounds begin
        if ys[g] <= y
            (g == hi || y < ys[g+1]) && return g
            g+1 == hi && return hi
            y < ys[g+2] && return g+1
        elseif g > lo && ys[g-1] <= y
            return g-1
        end
        y >= ys[hi] && return hi
        while hi - lo > 1
            m = (lo + hi) >>> 1
            if ys[m] <= y; lo = m; else; hi = m; end
        end
    end
    return lo
end

# Proprietà della sezione k per la quota eta, interpolate sulla tabella trasposta.
# Con y < 0 usa y = 0, oltre il tirante massimo estrapola con l'ultima riga
@inline function READ_VALUES(data, ys, first_row, last_row, k, eta, b_min)
    y = max(eta - b_min, 0.0)
    @inbounds begin
        r = FIND_ROW(ys, first_row[k], first_row[k+1]-1, y, last_row[k])
        last_row[k] = r
        d = y - ys[r]
        return (A   = data[3,r]  + data[4,r]*d,
                β   = data[5,r]  + data[6,r]*d,
                den = data[7,r]  + data[8,r]*d,
                w   = data[9,r]  + data[10,r]*d,
                rh  = data[11,r] + data[12,r]*d)
    end
end

function BARYCENTERS!(imax, data_depth, ys, first_row, last_row, b_min, eta,
                      A_b, β_b, den_γ_b, rh_b)
    @inbounds for i in 1:imax
        g = READ_VALUES(data_depth, ys, first_row, last_row, i, eta[i], b_min[i])
        A_b[i]     = g.A
        β_b[i]     = g.β
        den_γ_b[i] = g.den
        rh_b[i]    = g.rh
    end
    return  A_b, β_b, den_γ_b, rh_b
end

# Legge pendenza, ascissa e talweg di ogni sezione
function LOAD_GEOMETRY(filename, imax)
    data = Array{Float64}(undef, imax, 4)
    LOAD_DATA!(data, filename, imax, 4)
    return (bx = data[:, 2], x_b = data[:, 3], b_min = data[:, 4])
end

# Legge le tabelle delle sezioni: per ogni sezione una riga per ogni tirante
function LOAD_SECTION_TABLES(filename, n_row, n_col)
    data = Array{Float64}(undef, n_row, n_col)
    LOAD_DATA!(data, filename, n_row, n_col)
    return data
end

# Carico idrogrammi da file esterni
function LOAD_HYDROGRAPH!(filename)
    t = Float64[]
    Q = Float64[]

    for (n, line) in enumerate(eachline(filename))
        isempty(strip(line)) && continue
        # Dividi la linea in colonne (virgola o spazi) e converti in Float64
        val = parse.(Float64, split(strip(line), [',', ' ', '\t']; keepempty=false))
        if length(val) != 2
            error("Errore: la riga $n del file '$filename' non ha 2 colonne")
        end
        push!(t, val[1])
        push!(Q, val[2])
    end

    if isempty(t)
        error("Errore: il file '$filename' non contiene dati")
    end

    return (t = t, Q = Q)
end

# Carica gli idrogrammi scritti in settings.txt e presenti su disco, es. hydro[:centa].Q
function LOAD_ALL_HYDROGRAPHS!(params)
    files = (brenta1 = params.hydro_brenta1, brenta2 = params.hydro_brenta2,
             brenta3 = params.hydro_brenta3, brenta4 = params.hydro_brenta4,
             brenta5 = params.hydro_brenta5, brenta6 = params.hydro_brenta6,
             brenta7 = params.hydro_brenta7, brenta_new = params.hydro_brenta_new,
             brenta_v1 = params.hydro_brenta_v1, brenta_v2 = params.hydro_brenta_v2,
             brenta_v3 = params.hydro_brenta_v3, brenta_v4 = params.hydro_brenta_v4,
             caldonazzo = params.hydro_caldonazzo, ceggio = params.hydro_ceggio,
             levico = params.hydro_levico, centa = params.hydro_centa,
             chieppena = params.hydro_chieppena, grigno = params.hydro_grigno,
             maso = params.hydro_maso, moggio = params.hydro_moggio)

    hydro = Dict{Symbol, @NamedTuple{t::Vector{Float64}, Q::Vector{Float64}}}()
    for (name, filename) in pairs(files)
        if isfile(filename)
            hydro[name] = LOAD_HYDROGRAPH!(filename)
        else
            println("Attenzione: idrogramma '$filename' non trovato, non viene caricato")
        end
    end

    return hydro
end

# Tempo finale della simulazione: il minore tra i tempi finali degli idrogrammi
function HYDROGRAPH_END_TIME(hydro)
    t_fin = [v.t[end] for v in values(hydro)]
    tend  = minimum(t_fin)
    if maximum(t_fin) != tend
        println("Attenzione: tempi finali diversi, uso tend = $tend")
    end
    return tend
end

function LAGRANGE_INTERPOLATION!(nr, tk, Qk, time, dt)
    dτ = (tk[nr] - tk[1]) / (nr - 1)
    m = Int(max(2, min(nr - 1, ceil((time + dt) / dτ))))

    t = time + dt
    @inbounds begin
        Qm1, Qm, Qp1 = Qk[m-1], Qk[m], Qk[m+1]
        tm1, tm, tp1 = tk[m-1], tk[m], tk[m+1]
    end

    Q = Qm1 * ((t - tm) * (t - tp1)) / ((tm1 - tm) * (tm1 - tp1)) +
        Qm  * ((t - tm1) * (t - tp1)) / ((tm - tm1) * (tm - tp1)) +
        Qp1 * ((t - tm1) * (t - tm)) / ((tp1 - tm1) * (tp1 - tm))

    return Q
end

function THOMAS_ALGORITHM!(imax, aa, cc, di, ff, deta)
    gn = 0.0
    cn = Vector{Float64}(undef, imax)
    fn = Vector{Float64}(undef, imax)
    cn[1] = cc[1]/di[1]
    fn[1] = ff[1]/di[1]
    @inbounds for i in 2:imax
        gn = di[i] - cn[i-1]*aa[i]
        cn[i] = cc[i]/gn
        fn[i] = (ff[i] - aa[i]*fn[i-1])/gn
    end
    deta[imax] = fn[imax]
    @inbounds for i in (imax-1):-1:1
        deta[i] = fn[i] - cn[i]*deta[i+1]
        if abs(deta[i]) > 100
        #@infiltrate
    end
    end

    
    return (deta)
end
