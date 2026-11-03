#include("module.jl")
#include("struct.jl")

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
