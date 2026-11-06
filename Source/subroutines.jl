# Costruisce la griglia longitudinale a partire dalle ascisse delle sezioni x_b:
# dx_b = distanza tra sezioni consecutive, x_int = posizione delle interfacce tra le celle,
# dx_int = lunghezza della cella attorno a ogni sezione. Riempie i tre vettori e li restituisce
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



# Individua dove inizia ogni sezione nella tabella trasposta data, in cui ogni colonna è una
# riga di parametri_sezioni. Restituisce first_row (prima riga di ogni sezione, più un valore
# finale n+1) e last_row (ultima riga usata, inizialmente la prima di ogni sezione).
# Il numero di righe può variare tra le sezioni
function SECTION_INDEX(data)
    n = size(data, 2)
    first_row = [1; [r for r in 2:n if data[1,r] != data[1,r-1]]; n+1]
    return first_row, first_row[1:end-1]
end



# Trova la riga r della sezione, compresa tra lo e hi, con ys[r] <= y < ys[r+1], dove ys sono
# le quote del pelo libero, anche con passo non costante. Prova prima l'ultima riga usata (g)
# e le vicine, altrimenti usa la ricerca binaria. Oltre la quota massima restituisce hi
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



# Calcola le proprietà della sezione k per la quota del pelo libero eta, interpolando sulla
# tabella trasposta con le quote ys. Restituisce area A, coefficiente β, den (denominatore di γ),
# larghezza w e raggio idraulico rh. Sotto la prima quota (il fondo) usa la prima riga.
# Oltre la quota massima estrapola dall'ultima riga con la pendenza dell'ultimo tratto della
# tabella (tra le ultime due quote), calcolata dai valori. Aggiorna last_row[k] con la riga trovata
@inline function READ_VALUES(data, ys, first_row, last_row, k, eta)
    @inbounds begin
        lo = first_row[k]
        hi = first_row[k+1] - 1
        z  = max(eta, ys[lo])
        r  = FIND_ROW(ys, lo, hi, z, last_row[k])
        last_row[k] = r
        d = z - ys[r]
        if r == hi && d > 0.0 && hi > lo
            p  = hi - 1
            dy = ys[hi] - ys[p]
            return (A   = data[3,hi]  + (data[3,hi]  - data[3,p]) /dy*d,
                    β   = data[5,hi]  + (data[5,hi]  - data[5,p]) /dy*d,
                    den = data[7,hi]  + (data[7,hi]  - data[7,p]) /dy*d,
                    w   = data[9,hi]  + (data[9,hi]  - data[9,p]) /dy*d,
                    rh  = data[11,hi] + (data[11,hi] - data[11,p])/dy*d)
        end
        return (A   = data[3,r]  + data[4,r]*d,
                β   = data[5,r]  + data[6,r]*d,
                den = data[7,r]  + data[8,r]*d,
                w   = data[9,r]  + data[10,r]*d,
                rh  = data[11,r] + data[12,r]*d)
    end
end

# Aggiorna area, β, denominatore di γ e raggio idraulico in tutte le sezioni,
# leggendoli dalle tabelle in base alla quota eta corrente
function BARYCENTERS!(imax, data_depth, ys, first_row, last_row, eta,
                      A_b, β_b, den_γ_b, rh_b)
    @inbounds for i in 1:imax
        g = READ_VALUES(data_depth, ys, first_row, last_row, i, eta[i])
        A_b[i]     = g.A
        β_b[i]     = g.β
        den_γ_b[i] = g.den
        rh_b[i]    = g.rh
    end
    return  A_b, β_b, den_γ_b, rh_b
end

# Quota del pelo libero di moto uniforme (Manning) per la portata Q nella sezione k, con pendenza S.
# Cerca per bisezione tra fondo + y_min e la quota massima della tabella, sfruttando il fatto che
# la portata di Manning cresce con la quota. Se già con y_min la portata supera Q restituisce
# fondo + y_min: il tirante minimo evita veli troppo sottili (magre, pendenze quasi nulle)
function UNIFORM_LEVEL(data, ys, first_row, last_row, k, Q, Ks, S; y_min = 0.10)
    Qm(e) = (s = READ_VALUES(data, ys, first_row, last_row, k, e); s.A*Ks*sqrt(abs(S))*s.rh^(2/3))
    a = ys[first_row[k]] + y_min
    Qm(a) >= Q && return a
    b = ys[first_row[k+1]-1]
    for _ in 1:60
        m = 0.5*(a + b)
        Qm(m) < Q ? (a = m) : (b = m)
    end
    return 0.5*(a + b)
end



# Legge il file della geometria e restituisce pendenza locale bx, ascissa x_b, talweg b_min e
# quote degli argini destro e sinistro di ogni sezione. I file a 4 colonne non hanno gli argini:
# le quote valgono Inf, quindi non c'è mai sfioro
function LOAD_GEOMETRY(filename, imax)
    _, n_col = FILE_SIZE(filename)
    data = Array{Float64}(undef, imax, n_col)
    LOAD_DATA!(data, filename, imax, n_col)
    argine_dx = n_col >= 6 ? data[:, 5] : fill(Inf, imax)
    argine_sx = n_col >= 6 ? data[:, 6] : fill(Inf, imax)
    return (bx = data[:, 2], x_b = data[:, 3], b_min = data[:, 4], argine_dx, argine_sx)
end



# Portata sfiorata da un argine per unità di lunghezza [m2/s], con la quota eta del pelo libero
# e la quota z_argine della cresta (free overflow, formula 16 di Siviglia et al., ASCE 2009).
# Vale zero se il livello non supera l'argine
@inline function FREE_OVERFLOW(eta, z_argine, Cd, g)
    Y0 = max(eta - z_argine, 0.0)
    return Cd*sqrt(2.0*g*Y0^3)
end



# Legge il file delle tabelle delle sezioni e le restituisce come matrice n_row x n_col:
# per ogni sezione una riga per ogni tirante
function LOAD_SECTION_TABLES(filename, n_row, n_col)
    data = Array{Float64}(undef, n_row, n_col)
    LOAD_DATA!(data, filename, n_row, n_col)
    return data
end



# Legge un singolo idrogramma da file e restituisce i vettori dei tempi t e delle portate Q.
# Accetta colonne separate da virgole, spazi o tabulazioni e salta le righe vuote.
# Si ferma con un errore se una riga non ha 2 colonne o se il file non contiene dati
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



# Carica tutti gli idrogrammi scritti in settings.txt chiamando LOAD_HYDROGRAPH! per ciascuno.
# Restituisce un dizionario, es. hydro[:centa].Q; i file non presenti vengono saltati con un avviso
function LOAD_ALL_HYDROGRAPHS!(params)
    hydro = Dict{Symbol, @NamedTuple{t::Vector{Float64}, Q::Vector{Float64}}}()
    for (name, filename) in pairs(params.hydro_files)
        if isfile(filename)
            hydro[name] = LOAD_HYDROGRAPH!(filename)
        else
            println("Attenzione: idrogramma '$filename' non trovato, non viene caricato")
        end
    end

    return hydro
end



# Calcola il tempo finale della simulazione: il minore tra i tempi finali degli idrogrammi.
# Stampa un avviso se i tempi finali non sono tutti uguali
function HYDROGRAPH_END_TIME(hydro)
    t_fin = [v.t[end] for v in values(hydro)]
    tend  = minimum(t_fin)
    if maximum(t_fin) != tend
        println("Attenzione: tempi finali diversi, uso tend = $tend")
    end
    return tend
end



# Calcola la portata dell'idrogramma al tempo time + dt con l'interpolazione di Lagrange su
# tre punti. Sceglie i punti supponendo un passo temporale costante e tk che parte da zero
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



# Risolve con l'algoritmo di Thomas il sistema tridiagonale di Newton: aa, di e cc sono le
# diagonali inferiore, principale e superiore, ff il termine noto. Scrive la soluzione in deta
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

# Ferma la simulazione se in uno dei vettori passati compare un NaN, indicando dove:
# fase del calcolo, passo, tempo e sezioni o interfacce coinvolte
function CHECK_NAN(fase, nCS, time, IC; vettori...)
    for (nome, v) in vettori
        idx = findall(isnan, v)
        if !isempty(idx)
            error("NaN in $nome dopo $fase: passo $nCS (" * (IC == 1 ? "condizioni iniziali" : "simulazione") *
                  "), t = $time s, posizioni $(idx[1:min(end,10)]) ($(length(idx)) in totale)")
        end
    end
end
