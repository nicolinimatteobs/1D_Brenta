include("module.jl")
#include("struct.jl")

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
     @inbounds for i in 1:imax-1
        dx_b[i] = x_b[i+1] - x_b[i]
    end

    dx_int[1] = dx_b[1]
    @inbounds  for i in 2:imax-1
        dx_int[i] = 0.5 * (dx_b[i] + dx_b[i-1])
    end
    dx_int[imax] = dx_b[imax-1]

    @inbounds for i in 2:imax
        x_int[i] = 0.5 * (x_b[i] + x_b[i-1])
    end

    x_int[1]      = -x_int[2]
    x_int[imax+1] = x_b[imax] + (x_b[imax] - x_int[imax])

    return dx_b, dx_int, x_int
end

function READ_VALUES!(data, k, y, points_per_section, val)

    # Limiti sicuri per i loop
    high = points_per_section * k
    low  = high - (points_per_section - 1)
    #@infiltrate
    for j in low:high
        # Caso esatto
        #@infiltrate
        if abs(data[j,2] - y) < 1e-6
            val[1:12] = data[j, 1:12]
            break

        # Interpolazione lineare (solo se j+1 dentro bounds)
        elseif (data[j,2] < y) && (j + 1 <= size(data,1)) && (data[j+1,2] > y)
            val[1 ]   = data[j,1 ]
            val[2 ]   = y
            val[3 ]   = data[j,3 ] + data[j,4 ]*(y - data[j,2])
            val[5 ]   = data[j,5 ] + data[j,6 ]*(y - data[j,2])
            val[7 ]   = data[j,7 ] + data[j,8 ]*(y - data[j,2])
            val[9 ]   = data[j,9 ] + data[j,10]*(y - data[j,2])
            val[11]   = data[j,11] + data[j,12]*(y - data[j,2])
            break

        # Caso in cui y supera l’ultimo dato disponibile
        elseif (data[j,2] < y) && (j + 1 > size(data,1))
            println("Warning: H = ", y, " out of bounds at section ", data[j,1], ". Using last available data.")
            #@infiltrate
            val[1 ] = data[j,1 ]
            val[2 ] = y
            val[3 ] = data[j,3 ] + data[j,4 ]*(y - data[j,2])
            val[5 ] = data[j,5 ] + data[j,6 ]*(y - data[j,2])
            val[7 ] = data[j,7 ] + data[j,8 ]*(y - data[j,2])
            val[9 ] = data[j,9 ] + data[j,10]*(y - data[j,2])
            val[11] = data[j,11] + data[j,12]*(y - data[j,2])
            break
        end
    end
    return val
end


function BARYCENTERS!(imax, n_line_per_section, data_depth, b_min, eta, 
                      A_b, β_b, den_γ_b, rh_b, val)
    @inbounds @simd for i in 1:imax
        
        yb = eta[i] - b_min[i]
        READ_VALUES!(data_depth, i, yb, n_line_per_section, val)
        A_b[i]     = val[3]
        β_b[i]     = val[5]
        den_γ_b[i] = val[7]
        rh_b[i]    = val[11]
    end
    return  A_b, β_b, den_γ_b, rh_b
end

function LAGRANGE_INTERPOLATION!(nr_PC, tk_PC, Qk_PC, time, dt)
    dτ = (tk_PC[nr_PC] - tk_PC[1]) / (nr_PC - 1)
    m = Int(max(2, min(nr_PC - 1, ceil((time + dt) / dτ))))

    t = time + dt
    @inbounds begin
        Qm1, Qm, Qp1 = Qk_PC[m-1], Qk_PC[m], Qk_PC[m+1]
        tm1, tm, tp1 = tk_PC[m-1], tk_PC[m], tk_PC[m+1]
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