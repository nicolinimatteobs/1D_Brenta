# Crea la cartella degli output e svuota i file di una simulazione precedente.
# Restituisce i percorsi dei file: IC_* per le condizioni iniziali, gli altri per la simulazione
function INIT_OUTPUT(out_dir)
    mkpath(out_dir)
    nomi  = (:IC_Q, :IC_eta, :IC_h, :IC_time, :IC_Q_lat, :Q, :eta, :h, :time, :Q_lat)
    files = NamedTuple{nomi}(Tuple(joinpath(out_dir, string(n, ".csv")) for n in nomi))
    for f in files
        write(f, "")
    end
    return files
end

function WRITE_INITIAL_CONDITIONS!(ICdischarge, imax, x_b, eta,b_min, x_int, q)
    # Scrive le condizioni iniziali in un file
    open(ICdischarge, "w") do io
        for i in 1:imax
            @printf(io, "%12.4f,%12.4f,%12.4f,%12.4f\n", x_b[i], eta[i]- b_min[i], x_int[i], q[i])
        end
        # Riga finale con x_int(imax+1) e q(imax+1)
        @printf(io, "%22s,%1s,%12.4f,%12.4f\n", "", "", x_int[imax+1], q[imax+1])
    end
end

function initialize_Q!(filename::String, imax::Int, x_int::Vector{Float64})
    open(filename, "w") do io
        for i in 1:imax+1
            @printf(io, "%12d ", x_int[i])
        end
        write(io, "\n")
    end
end

# --- Funzione per salvare la portata di ogni istante come colonna ---
function save_Q!(filename::String, q::Vector{Float64})
    open(filename, "a") do io
        for val in q
            @printf(io, "%12.6f ", val)
        end
        write(io, "\n")
    end
end

function save_eta!(filename::String, eta::Vector{Float64})
    open(filename, "a") do io
        for val in eta
            @printf(io, "%12.6f ", val)
        end
        write(io, "\n")
    end
end

function save_h!(filename::String, h::Vector{Float64})
    open(filename, "a") do io
        for val in h
            @printf(io, "%12.6f ", val)
        end
        write(io, "\n")
    end
end

function save_U!(filename::String, U::Vector{Float64})
    open(filename, "a") do io
        for val in U
            @printf(io, "%12.6f ", val)
        end
        write(io, "\n")
    end
end
function save_A!(filename::String, A::Vector{Float64})
    open(filename, "a") do io
        for val in A
            @printf(io, "%12.6f ", val)
        end
        write(io, "\n")
    end
end

function save_fondo!(filename, imax, bx, b_min, eta, x_b, q)
    open(filename, "w") do io
        for i in 1:imax
            @printf(io, "%12.6f %12.6f %12.6f %12.6f %12.6f\n", bx[i], b_min[i], eta[i], x_b[i], q[i])
        end
        write(io, "\n")
    end
end

# --- Funzione per salvare il tempo in un file separato ---
function save_time!(filename::String, time::Float64)
    open(filename, "a") do io
        @printf(io, "%12.6f\n", time)
    end
end

# Controllo del bilancio di massa della simulazione (dopo il riempimento), scritto in filename e stampato.
# Volumi [m3]: in alveo all'inizio e alla fine (aree delle tabelle per la lunghezza delle celle), entrato da
# monte, uscito a valle, entrato dagli affluenti, scambio netto con la piana (positivo se entra).
# errore = variazione del volume in alveo − (monte − valle + affluenti + scambio), rispetto al volume da monte;
# oltre lo 0,1% compare un avviso. Con lo scambio 'tabella' controlla anche il serbatoio della piana: il volume
# che contiene deve essere uguale all'acqua uscita e non ancora rientrata dall'inizio, riempimento compreso.
function BILANCIO_DI_MASSA(filename, scambio, V_ini, V_fin, V_monte, V_valle, V_aff, V_lat, V_piana, V_lat_tutto)
    err   = (V_fin - V_ini) - (V_monte - V_valle + V_aff + V_lat)
    somma = -err   # somma di tutti i termini con il loro segno: deve essere zero
    rel   = V_monte > 0.0 ? somma / V_monte : NaN
    io = IOBuffer()
    println(io, "Bilancio di massa della simulazione (dal termine del riempimento alla fine), volumi in m3")
    println(io, "Somma di tutti i termini con il loro segno (deve essere zero):")
    @printf(io, "  %+18.1f   volume in alveo a fine riempimento\n", V_ini)
    @printf(io, "  %+18.1f   entrato da monte\n", V_monte)
    @printf(io, "  %+18.1f   uscito a valle\n", -V_valle)
    @printf(io, "  %+18.1f   entrato dagli affluenti\n", V_aff)
    @printf(io, "  %+18.1f   scambio con la piana (+ entra, − esce)\n", V_lat)
    @printf(io, "  %+18.1f   volume in alveo alla fine\n", -V_fin)
    @printf(io, "  %18s\n", "------------------")
    @printf(io, "  %+20.3f somma = errore del bilancio (%.2e del volume da monte)\n", somma, rel)
    abs(rel) > 1e-3 && println(io, "ATTENZIONE: errore del bilancio oltre lo 0,1%")
    if scambio == :tabella
        err_p = V_piana + V_lat_tutto
        @printf(io, "%-42s %16.1f\n", "serbatoio della piana alla fine", V_piana)
        @printf(io, "%-42s %16.1f  (deve essere uguale al serbatoio)\n", "uscito e non rientrato, riempimento compreso", -V_lat_tutto)
        abs(err_p) > 1e-6 * max(abs(V_lat_tutto), 1.0) && @printf(io, "ATTENZIONE: il serbatoio della piana non torna (%.3f m3)\n", err_p)
    end
    s = String(take!(io))
    write(filename, s); print(s)
    return err
end
