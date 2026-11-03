
#include("module.jl")
#include("struct.jl")
#include("subroutines.jl")
#include("writing.jl")

@views @inbounds function SEMI_IMPLICIT!()
    # --- INIZIALIZZAZIONE STRUTTURE ---
    params     = INITIAL_SETTINGS!()

    # --- CREAZIONE ALIAS PER ACCEDERE AI CAMPI DELLE STRUCT ---
    @unpack imax, n_col, n_line_per_section, Ks, dtmax, Nmax = params

    struttura    = params.struttura;        nr_strut  = params.nr_strut;   parametri_sezioni = params.parametri_sezioni
    local_slope  = params.local_slope;   hydro_Brenta = params.hydro_Brenta;    hydro_Noce = params.hydro_Noce
    nr_PC        = params.nr_PC;          ICdischarge = params.ICdischarge;         Q_output = params.Q_output
    nc_PC        = params.nc_PC;                 time = params.time;                    temp = params.temp
    nr_TN        = params.nr_TN;                n_row = params.n_row;                  nc_TN = params.nc_TN;  

    val, b, bx, b_min, x_b, dx_b, x_int, dx_int, eta, A_bn, deta, rh_b, rh_int,
    A_int, A_int_g, A_b, dA, q, uf, qstar, γ, aa, bb, cc, di, ff, smax, β_b,
    β_int, rhs, den_γ_b, den_γ_int, qn = ALLOCATE_ARRAYS!(imax, n_col)

    # --- CARICAMENTO IDROGRAMMI ---
    data_Brenta = Array{Float64}(undef, nr_PC, nc_PC)
    PROCESS_HYDROGRAPH!(data_Brenta, hydro_Brenta, nr_PC, nc_PC)
    tk_PC = data_Brenta[:, 1]
    Qk_PC = data_Brenta[:, 2]
    data_Brenta = nothing  # Libera memoria

    # --- CARICAMENTO DATI E INIZIALIZZAZIONE DOMINIO ---
    n_row_ifx = imax
    n_col_ifx = 4
    data_ifx = Array{Float64}(undef, n_row_ifx, n_col_ifx)
    LOAD_DATA!(data_ifx, local_slope, n_row_ifx, n_col_ifx)
    bx    = data_ifx[:, 2]
    x_b   = data_ifx[:, 3]
    b_min = data_ifx[:, 4]
    data_ifx = nothing  # Libera memoria

    #@infiltrate

    data_depth = Array{Float64}(undef, n_row, n_col)
    LOAD_DATA!(data_depth, parametri_sezioni, n_row, n_col)

    #@infiltrate
    

    SPATIAL_DOMAIN!(imax ,x_b, dx_b, x_int, dx_int)

    Qin = Qk_PC[1]  # Portata in ingresso iniziale
    ym = 0.0
    exp_dt = 2.0/3.0  # Esponente di dt nel calcolo della portata
    nCS = 0
    solution_counter = 0       
    convergence_counter = 0
    g = 9.81

    #@infiltrate

    @inbounds @simd for i in 1:imax
        #println("Inside loop, i = ", i)
        yl = 0.0
        yr = 5.0
        while abs(yr - yl) > 1e-8
            ym = 0.5 * (yl + yr)

            READ_VALUES!(data_depth, i, ym, n_line_per_section, val)

            qm  = val[3]*Ks*sqrt(abs(bx[i]))*val[11]^(exp_dt)

            dqm = qm - Qin
           
            if dqm < 0.0
                yl = ym
            else
                yr = ym
            end
        end
        eta[i] = b_min[i] + ym
    end

    #@infiltrate
   
    q, A_b, β_b, den_γ_b, rh_b = BARYCENTERS!(imax, n_line_per_section, data_depth, b_min, eta, A_b, β_b, den_γ_b, rh_b, q, Qin, val)  
    
    #@infiltrate# Centroidi

    WRITE_INITIAL_CONDITIONS!(ICdischarge, imax, x_b, eta, x_int, q)  # Scrittura condizioni iniziali

    tend = tk_PC[nr_PC]
    dx   = min(minimum(dx_b), minimum(dx_int))

    #@infiltrate

    initialize_Q!("Output/VAIA/dx_70/Q_Brenta_Vaia_dx_70.csv", imax, x_int)

    open("Output/VAIA/dx_70/time_Brenta_Vaia_dx_70.csv", "w") do io
        write(io, "")   # reset file
    end

    #initialize_Q!("Output/VAIA/dx_70/qn_Vaia_dx_70.csv", imax, x_int)
    qn = copy(q)  # Salva la portata precedente

    #@infiltrate

    @time while nCS < Nmax

    qn .= q  # Salva la portata precedente

        # Aggiorna tol e CFL
        if solution_counter == 0
            tol = 1e-10
            CFL = 0.6
        else
            tol = 1e-10
            CFL = 0.95
        end

        # --- INTERPOLAZIONE PARAMETRI INTERNI ---
        @inbounds for i in 1:imax+1
            if i == 1
                A_int[i]     = A_b[i]
                β_int[i]     = β_b[i]
                den_γ_int[i] = den_γ_b[i]
                rh_int[i]    = rh_b[i]
            elseif i == imax+1
                A_int[i]     = A_b[imax]
                β_int[i]     = β_b[imax]
                den_γ_int[i] = den_γ_b[imax]
                rh_int[i]    = rh_b[imax]
            else
                A_int[i]     = 0.5 * (A_b[i] + A_b[i-1])
                β_int[i]     = 0.5 * (β_b[i] + β_b[i-1])
                den_γ_int[i] = 0.5 * (den_γ_b[i] + den_γ_b[i-1])
                rh_int[i]    = 0.5 * (rh_b[i] + rh_b[i-1])
            end
 
            uf[i]   = q[i] * A_int[i]/(A_int[i]^2 + 1e-3)
            smax[i] = abs(β_int[i] * uf[i])

        end


        amax = maximum(smax)
        dt   = min(dtmax, CFL*dx/(amax + 1e-14))

        #if isnan(dt)
        #    @show dt amax 
        #    @infiltrate
        #end
        


        if time + dt >= tend
            dt = tend - time
        end
        if time >= tend
            break
        end

        # --- Aggiornamento portate qstar ---
        @inbounds for i in 1:imax+1
            if i == 1
                if solution_counter == 0
                    qstar[i] = Qin
                else
                    Qin = LAGRANGE_INTERPOLATION!(nr_PC, tk_PC, Qk_PC, time, dt)
                    qstar[i] = Qin
                end
            elseif i == imax+1
                uout = Ks*sqrt(abs(bx[imax]))*rh_int[i]^exp_dt
                qstar[i] = A_int[i]*uout
            else
                fp = 0.5*(β_int[i]*uf[i]*q[i] + β_int[i+1]*uf[i+1]*q[i+1]) - 0.5*2*max(smax[i], smax[i+1])*(q[i+1]-q[i])
                fm = 0.5*(β_int[i-1]*uf[i-1]*q[i-1] + β_int[i]*uf[i]*q[i]) - 0.5*2*max(smax[i-1], smax[i])*(q[i]-q[i-1])
                qstar[i] = q[i] - dt/dx_b[i-1]*(fp - fm)
            end
            γ[i]       = abs(uf[i])/den_γ_int[i]
            A_int_g[i] = A_int[i]/(1.0 + dt*γ[i])
        end

        #if isnan(qstar)
        #    @show 
        #    @infiltrate
        #end
        # --- Aggiornamento geometrie e preparazione RHS per Thomas ---
        A_bn .= A_b

        @inbounds for i in 1:imax
            if i == 1
                aa[i]  = 0.0 
                cc[i]  = -g*dt^2/(dx_int[i]*dx_b[i])*A_int_g[i+1]
                bb[i]  = -aa[i] - cc[i] 
                rhs[i] = A_bn[i] - dt/dx_int[i]*(qstar[i+1]/(1+dt*γ[i+1]) - qstar[i])
            elseif i == imax
                aa[i]  = -g*dt^2/(dx_int[i]*dx_b[i-1])*A_int_g[i] 
                cc[i]  = 0.0
                bb[i]  = -aa[i] - cc[i] 
                rhs[i] = A_bn[i] - dt/dx_int[i]*(qstar[i+1] - qstar[i]/(1+dt*γ[i]))
            else 
                aa[i]  = -g*dt^2/(dx_int[i]*dx_b[i-1])*A_int_g[i] 
                cc[i]  = -g*dt^2/(dx_int[i]*dx_b[i])*A_int_g[i+1]
                bb[i]  = -aa[i] - cc[i] 
                rhs[i] = A_bn[i] - dt/dx_int[i]*(qstar[i+1]/(1+dt*γ[i+1]) - qstar[i]/(1+dt*γ[i]))
            end
        end

        # --- Ciclo Newton ---
        @inbounds for iNewton in 1:100
            @inbounds for i in 1:imax
                yb = eta[i] - b_min[i]
                READ_VALUES!(data_depth, i, yb, n_line_per_section, val)
                A_bn[i] = val[3]
                dA[i]   = val[9]

                if i == 1
                    ff[i] = A_bn[i] + bb[i]*eta[i] + cc[i]*eta[i+1] - rhs[i]
                elseif i == imax
                    ff[i] = A_bn[i] + aa[i]*eta[i-1] + bb[i]*eta[i] - rhs[i]
                else
                    ff[i] = A_bn[i] + aa[i]*eta[i-1] + bb[i]*eta[i] + cc[i]*eta[i+1] - rhs[i]
                end

                di[i] = dA[i] + bb[i]

                if rh_int[i] + rh_int[i+1] < 1e-3
                    di[i] = 1.0
                    ff[i] = 0.0
                end
            end

            res = sqrt(sum(ff.^2))
            if res < tol
                break
            end

            THOMAS_ALGORITHM!(imax, aa, cc, di, ff, deta)
            eta .= eta .- deta

        end

        # --- Aggiornamento portate finali ---
        @inbounds for i in 1:imax
            if i == 1
                q[i] = qstar[i]
            elseif i == imax
                q[i] = qstar[i]
            else
                q[i] = (qstar[i] - g*A_int[i]*dt/dx_b[i-1]*(eta[i]-eta[i-1]))/(1.0 + dt*γ[i])
            end
        end

        #if isnan(dt)
        #    @show 
        #    @infiltrate
        #end

        # --- Aggiornamento tempo e counter ---
        time += dt
        nCS  += 1
        # --- Salvataggio risultati ---
        if nCS % 500 == 0
            println("Time step: ", nCS, " Time: ", @sprintf("%.4f", time), " dt: ", @sprintf("%.6f", dt))
        end
        if nCS % 100 == 0 #&& solution_counter == 1
            save_Q!("Output/VAIA/dx_70/Q_Brenta_Vaia_dx_70.csv", q)
            #save_Q!("Output/PIENA_31102023/qn_31102023.csv", qn)
            save_time!("Output/VAIA/dx_70/time_Brenta_Vaia_dx_70.csv", time)
        end



        # --- Aggiornamento geometrie finali ---
        @inbounds @simd for i in 1:imax
            yb = eta[i] - b_min[i]
            READ_VALUES!(data_depth, i, yb, n_line_per_section, val)
            A_b[i]     = val[3]
            β_b[i]     = val[5]
            den_γ_b[i] = val[7]
            b[i]       = val[9]
            rh_b[i]    = val[11]
        end



        # --- Controllo stato stazionario ---
        if solution_counter == 0
            convergence_counter = 0
            @inbounds for i in 1:imax
                #@infiltrate
                if abs(q[i] - qn[i]) < 0.01
                    convergence_counter += 1
                end
            end
            if convergence_counter == imax
                println("Converged to steady state at time step ", nCS, ", time = ", time)
                #@infiltrate
                solution_counter = 1
                nCS = 1
                time = 0.0
            end
        end
        
        
        
        
        
        
        
        
    end

end