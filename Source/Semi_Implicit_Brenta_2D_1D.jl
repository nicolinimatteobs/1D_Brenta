@views @inbounds function SEMI_IMPLICIT_BRENTA_2D_1D!()

    # --- INIZIALIZZAZIONE STRUTTURE ---
    params = INITIAL_SETTINGS!()

    # --- CREAZIONE ALIAS PER ACCEDERE AI CAMPI DELLE STRUCT ---
    @unpack imax, n_col, n_row, Ks, dtmax, Nmax,
            parametri_sezioni, local_slope, output_dir, time = params

    # --- ALLOCAZIONE ARRAY ---
    @unpack b, dx_b, x_int, dx_int, eta, A_bn, deta, rh_b, rh_int,
            A_int, A_int_g, A_b, dA, q, uf, qstar, γ, aa, bb, cc, di, ff,
            smax, β_b, β_int, rhs, den_γ_b, den_γ_int, qn = ALLOCATE_ARRAYS!(imax)

    # --- CARICAMENTO DATI SEZIONI E CREAZIONE DOMINIO LONGITUDINALE ---
    @unpack bx, x_b, b_min, argine_dx, argine_sx = LOAD_GEOMETRY(local_slope, imax)
    data_depth = permutedims(LOAD_SECTION_TABLES(parametri_sezioni, n_row, n_col))
    first_row, last_row = SECTION_INDEX(data_depth)
    length(last_row) == imax || error("Le tabelle hanno $(length(last_row)) sezioni, la geometria ($local_slope) ne ha $imax")
    ys         = copy(data_depth[2, :])  # quote del pelo libero (colonna 2); copy: con @views sarebbe una vista
    SPATIAL_DOMAIN!(imax, x_b, dx_b, x_int, dx_int)

    # --- CARICAMENTO IDROGRAMMI ---
    hydro = LOAD_ALL_HYDROGRAPHS!(params)
    tend  = HYDROGRAPH_END_TIME(hydro)

    # --- FILE DI OUTPUT: cartella indicata nella riga 23 di settings.txt ---
    out = INIT_OUTPUT(output_dir)

    # --- ALIAS IDROGRAMMI (tempi e portate) ---
    tk_cald       = hydro[:caldonazzo].t;  Qk_cald       = hydro[:caldonazzo].Q
    tk_brenta1    = hydro[:brenta1].t;     Qk_brenta1    = hydro[:brenta1].Q
    tk_brenta2    = hydro[:brenta2].t;     Qk_brenta2    = hydro[:brenta2].Q
    tk_brenta3    = hydro[:brenta3].t;     Qk_brenta3    = hydro[:brenta3].Q
    tk_brenta4    = hydro[:brenta4].t;     Qk_brenta4    = hydro[:brenta4].Q
    tk_brenta5    = hydro[:brenta5].t;     Qk_brenta5    = hydro[:brenta5].Q
    tk_brenta6    = hydro[:brenta6].t;     Qk_brenta6    = hydro[:brenta6].Q
    tk_brenta_new = hydro[:brenta_new].t;  Qk_brenta_new = hydro[:brenta_new].Q
    tk_brenta_v1  = hydro[:brenta_v1].t;   Qk_brenta_v1  = hydro[:brenta_v1].Q
    tk_brenta_v2  = hydro[:brenta_v2].t;   Qk_brenta_v2  = hydro[:brenta_v2].Q
    tk_brenta_v3  = hydro[:brenta_v3].t;   Qk_brenta_v3  = hydro[:brenta_v3].Q
    tk_brenta_v4  = hydro[:brenta_v4].t;   Qk_brenta_v4  = hydro[:brenta_v4].Q
    tk_ceggio     = hydro[:ceggio].t;      Qk_ceggio     = hydro[:ceggio].Q
    tk_levico     = hydro[:levico].t;      Qk_levico     = hydro[:levico].Q
    tk_centa      = hydro[:centa].t;       Qk_centa      = hydro[:centa].Q
    tk_chieppena  = hydro[:chieppena].t;   Qk_chieppena  = hydro[:chieppena].Q
    tk_grigno     = hydro[:grigno].t;      Qk_grigno     = hydro[:grigno].Q
    tk_maso       = hydro[:maso].t;        Qk_maso       = hydro[:maso].Q
    tk_moggio     = hydro[:moggio].t;      Qk_moggio     = hydro[:moggio].Q
    nr_hydro      = length(tk_cald)  # Numero di righe, uguale per tutti gli idrogrammi


    # ----- PARAMETRI PER RIEMPIMENTO ----
    Qin_idro      = Qk_cald[1]  # Portata in ingresso iniziale
    Qin_IC_target = Qin_idro   # Portata obiettivo per condizioni iniziali
    y_ini          = 0.0    # Tirante minimo iniziale [m] (si parte dal moto uniforme, ma non sotto y_ini)
    T_riempimento  = 21600.0  # Durata della rampa di riempimento [s]
    dt_riempimento = 5.0    # Passo massimo durante la rampa [s]

    Q_brenta1_IC    = Qk_brenta1[1]    # Portata affluente bacino 1 condizioni iniziali
    Q_brenta2_IC    = Qk_brenta2[1]    # Portata affluente bacino 2 condizioni iniziali
    Q_brenta3_IC    = Qk_brenta3[1]    # Portata affluente bacino 3 condizioni iniziali
    Q_brenta4_IC    = Qk_brenta4[1]    # Portata affluente bacino 4 condizioni iniziali
    Q_brenta5_IC    = Qk_brenta5[1]    # Portata affluente bacino 5 condizioni iniziali
    Q_brenta6_IC    = Qk_brenta6[1]    # Portata affluente bacino 6 condizioni iniziali
    Q_brenta_new_IC = Qk_brenta_new[1] # Portata affluente bacino 7 condizioni iniziali
    Q_brenta_v1_IC  = Qk_brenta_v1[1]  # Portata affluente bacino 8 condizioni iniziali
    Q_brenta_v2_IC  = Qk_brenta_v2[1]  # Portata affluente bacino 8 condizioni iniziali
    Q_brenta_v3_IC  = Qk_brenta_v3[1]  # Portata affluente bacino 8 condizioni iniziali
    Q_brenta_v4_IC  = Qk_brenta_v4[1]  # Portata affluente bacino 8 condizioni iniziali
    Q_ceggio_IC     = Qk_ceggio[1]     # Portata affluente Ceggio condizioni iniziali
    Q_levico_IC     = Qk_levico[1]     # Portata affluente Levico condizioni iniziali
    Q_centa_IC      = Qk_centa[1]      # Portata affluente Centa condizioni iniziali
    Q_chieppena_IC  = Qk_chieppena[1]  # Portata affluente Chieppena condizioni iniziali
    Q_grigno_IC     = Qk_grigno[1]     # Portata affluente Grigno condizioni iniziali
    Q_maso_IC       = Qk_maso[1]       # Portata affluente Maso condizioni iniziali
    Q_moggio_IC     = Qk_moggio[1]     # Portata affluente Moggio condizioni iniziali

    # Sezioni di ingresso degli affluenti, distribuite equamente lungo il canale
    n_aff   = 18
    sez_aff = [round(Int, k*imax/(n_aff + 1)) for k in 1:n_aff]

    # Gli affluenti con picco oltre Q_soglia_aff entrano su 3 sezioni (sez_aff e le due vicine),
    # con la portata divisa in parti uguali; gli altri entrano in una sola sezione
    Q_soglia_aff = 10.0  # [m3/s]
    Qmax_aff = maximum.((Qk_brenta1, Qk_brenta2, Qk_brenta3, Qk_brenta4, Qk_brenta5, Qk_brenta6,
                         Qk_brenta_new, Qk_brenta_v1, Qk_brenta_v2, Qk_brenta_v3, Qk_brenta_v4,
                         Qk_ceggio, Qk_levico, Qk_centa, Qk_chieppena, Qk_grigno, Qk_maso, Qk_moggio))
    semi_aff = [Q >= Q_soglia_aff ? 1 : 0 for Q in Qmax_aff]  # sezioni in più per lato
    ncel_aff = 2 .* semi_aff .+ 1                              # sezioni su cui entra ogni affluente
    

    exp_dt        = 2.0/3.0 # Esponente di dt nel calcolo della portata
    nCS           = 1
    IC_counter    = 0       
    IC            = 1       # Indicatore condizioni iniziali attive
    consec_count  = 0 
    g             = 9.81


    # Stato iniziale: in ogni sezione tirante di moto uniforme per la prima portata dell'idrogramma
    # di monte (almeno y_ini), con la portata di moto uniforme locale. Così la portata di monte
    # non cambia durante il riempimento. L'attrito è quello delle tabelle (jden), lo stesso usato
    # dopo in moto vario: con il Ks di settings lo stato iniziale era lontano dall'equilibrio del
    # modello e lo spin-up più lungo
    @inbounds for i in 1:imax
        eta[i] = UNIFORM_LEVEL(data_depth, ys, first_row, last_row, i, Qin_IC_target, bx[i]; y_min = y_ini)
    end
    # Alternativa: velo d'acqua costante di y_ini in tutte le sezioni
    #eta .= b_min .+ y_ini

    @inbounds for i in 1:imax
        sec = READ_VALUES(data_depth, ys, first_row, last_row, i, eta[i])
        A_b[i]  = sec.A
        rh_b[i] = sec.rh
        q[i]    = UNIFORM_DISCHARGE(A_b[i], sec.den, bx[i])
        if i == imax
            q[i+1] += q[i]  # Aggiunta portata in ingresso
        end
    end 
  

    Qin_0 = q[1]  # Portata in ingresso all'inizio della rampa
    Qin   = Qin_0  # Portata in ingresso corrente
    Q_aff = zeros(imax)  # Vettore portate affluenti
    Q_lat = zeros(imax)  # Vettore portate uscenti per esondazione
    Cd    = 0.385        # Coefficiente di efflusso dello sfioro sugli argini

    A_b, β_b, den_γ_b, rh_b = BARYCENTERS!(imax, data_depth, ys, first_row, last_row, eta, A_b, β_b, den_γ_b, rh_b)  
    
    
    dx   = min(minimum(dx_b), minimum(dx_int))

    qn = copy(q)  # Salva la portata precedente

    @time while nCS < Nmax

    qn .= q  # Salva la portata precedente

    if IC == 1
        CFL = 0.6
        tol = 1e-14
    else
        CFL = 0.95
        tol = 1e-14
    end

        # --- INTERPOLAZIONE PARAMETRI INTERNI ---
        A_int[1]     = A_b[1]
        β_int[1]     = β_b[1]
        den_γ_int[1] = den_γ_b[1]
        rh_int[1]    = rh_b[1]

        A_int[imax+1]     = A_b[imax]
        β_int[imax+1]     = β_b[imax]
        den_γ_int[imax+1] = den_γ_b[imax]
        rh_int[imax+1]    = rh_b[imax]

        uf[1]   = q[1] * A_int[1]/(A_int[1]^2 + tol)
        smax[1] = abs(β_int[1] * uf[1])

        uf[imax+1]   = q[imax+1] * A_int[imax+1]/(A_int[imax+1]^2 + tol)
        smax[imax+1] = abs(β_int[imax+1] * uf[imax+1])
        #Threads.@threads :static
        @inbounds @simd for i in 2:imax
            
            # Upwind: area e attrito presi dalla cella a monte invece che dalla media.
            # Evita il dente di sega a monte degli affluenti quando le celle sono molto più lunghe
            # del tratto in cui il livello torna al moto uniforme (circa dx > 0,6*y/S), ma è meno
            # accurato nei tratti profondi
            A_int[i]     = q[i] >= 0.0 ? A_b[i-1] : A_b[i]
            β_int[i]     = 0.5 * (β_b[i] + β_b[i-1])
            den_γ_int[i] = q[i] >= 0.0 ? den_γ_b[i-1] : den_γ_b[i]
            rh_int[i]    = 0.5 * (rh_b[i] + rh_b[i-1])
            # Alternativa con la media (upwind sull'area solo durante la rampa: una cella non cede
            # più acqua di quanta ne ha). Per riattivarla sostituire le righe di A_int e den_γ_int con:
            #A_int[i]     = (IC == 1 && time < T_riempimento) ? (q[i] >= 0.0 ? A_b[i-1] : A_b[i]) : 0.5 * (A_b[i] + A_b[i-1])
            #den_γ_int[i] = 0.5 * (den_γ_b[i] + den_γ_b[i-1])
 
            uf[i]   = q[i] * A_int[i]/(A_int[i]^2 + tol)
            smax[i] = abs(β_int[i] * uf[i])

        end

        amax = maximum(smax)
        dt_lim = (IC == 1 && time < T_riempimento) ? dt_riempimento : dtmax
        dt   = min(dt_lim, CFL*dx/(amax + tol))
        CHECK_NAN("calcolo di dt", nCS, time, IC; dt, uf, smax)

        #@infiltrate

        if time + dt >= tend
            dt = tend - time
        end
        if time >= tend
            break
        end

        # Rampa lineare: le portate salgono al valore iniziale in T_riempimento
        if IC == 1
            r_IC = min(1.0, (time + dt)/T_riempimento)
            Qin  = Qin_0 + r_IC*(Qin_IC_target - Qin_0)
        end
        #@infiltrate
        @inbounds @simd for i in 1:imax+1
            if i == 1
                if IC == 1
                    qstar[i] = Qin
                else
                    Qin_idro = LAGRANGE_INTERPOLATION!(nr_hydro, tk_cald, Qk_cald, time, dt)
                    qstar[i] = Qin_idro
                end
            elseif i == imax+1
                uout = Ks*sqrt(abs(bx[imax]))*(rh_int[i]+tol)^exp_dt
                qstar[i] = A_int[i]*uout
                #@infiltrate
            else
                fp = 0.5*(β_int[i]*uf[i]*q[i] + β_int[i+1]*uf[i+1]*q[i+1]) - 0.5*2*max(smax[i], smax[i+1])*(q[i+1]-q[i])
                fm = 0.5*(β_int[i-1]*uf[i-1]*q[i-1] + β_int[i]*uf[i]*q[i]) - 0.5*2*max(smax[i-1], smax[i])*(q[i]-q[i-1])
                qstar[i] = q[i] - dt/dx_b[i-1]*(fp - fm)
            end
            γ[i]       = abs(uf[i])/(den_γ_int[i] + 1e-14)
            A_int_g[i] = A_int[i]/(1.0 + dt*γ[i])
        end
        CHECK_NAN("calcolo di qstar", nCS, time, IC; qstar, γ, A_int_g)
        # --- Aggiornamento geometrie e preparazione RHS per Thomas ---
        A_bn .= A_b

        aa[1]  = 0.0 
        cc[1]  = -g*dt^2/(dx_int[1]*dx_b[1])*A_int_g[2]
        bb[1]  = -aa[1] - cc[1] 
        rhs[1] = A_bn[1] - dt/dx_int[1]*(qstar[2]/(1+dt*γ[2]) - qstar[1])
        aa[imax]  = -g*dt^2/(dx_int[imax]*dx_b[imax-1])*A_int_g[imax] 
        cc[imax]  = 0.0
        bb[imax]  = -aa[imax] - cc[imax] 
        rhs[imax] = A_bn[imax] - dt/dx_int[imax]*(qstar[imax+1] - qstar[imax]/(1+dt*γ[imax]))
        #Threads.@threads :static
        @inbounds @simd for i in 2:imax-1

            aa[i]  = -g*dt^2/(dx_int[i]*dx_b[i-1])*A_int_g[i] 
            cc[i]  = -g*dt^2/(dx_int[i]*dx_b[i])*A_int_g[i+1]
            bb[i]  = -aa[i] - cc[i] 
            rhs[i] = A_bn[i] - dt/dx_int[i]*(qstar[i+1]/(1+dt*γ[i+1]) - qstar[i]/(1+dt*γ[i]))

        end

        # AFFLUENTI DISATTIVATI (prova su un tratto di Brenta reale): togliere i # dai due blocchi qui sotto per riattivarli
        ####  INSERIMENTO AFFLUENTI CONDIZIONI INIZIALI
        # Ogni affluente entra nella sezione sez_aff e, se ha portata alta, nelle due vicine
        if IC == 1
            @inbounds @simd for i in 1:imax
                if abs(i - sez_aff[1]) <= semi_aff[1]
                    Q_aff[i] += (Q_brenta1_IC) / ncel_aff[1]
                end
                if abs(i - sez_aff[2]) <= semi_aff[2]
                    Q_aff[i] += (Q_brenta2_IC) / ncel_aff[2]
                end
                if abs(i - sez_aff[3]) <= semi_aff[3]
                    Q_aff[i] += (Q_brenta3_IC) / ncel_aff[3]
                end
                if abs(i - sez_aff[4]) <= semi_aff[4]
                    Q_aff[i] += (Q_brenta4_IC) / ncel_aff[4]
                end
                if abs(i - sez_aff[5]) <= semi_aff[5]
                    Q_aff[i] += (Q_brenta5_IC) / ncel_aff[5]
                end
                if abs(i - sez_aff[6]) <= semi_aff[6]
                    Q_aff[i] += (Q_brenta6_IC) / ncel_aff[6]
                end
                # brenta7: da riattivare quando Brenta7.csv viene aggiunto (aumentare n_aff)
                if abs(i - sez_aff[7]) <= semi_aff[7]
                    Q_aff[i] += (Q_brenta_new_IC) / ncel_aff[7]
                end
                if abs(i - sez_aff[8]) <= semi_aff[8]
                    Q_aff[i] += (Q_brenta_v1_IC) / ncel_aff[8]
                end
                if abs(i - sez_aff[9]) <= semi_aff[9]
                    Q_aff[i] += (Q_brenta_v2_IC) / ncel_aff[9]
                end
                if abs(i - sez_aff[10]) <= semi_aff[10]
                    Q_aff[i] += (Q_brenta_v3_IC) / ncel_aff[10]
                end
                if abs(i - sez_aff[11]) <= semi_aff[11]
                    Q_aff[i] += (Q_brenta_v4_IC) / ncel_aff[11]
                end
                if abs(i - sez_aff[12]) <= semi_aff[12]
                    Q_aff[i] += (Q_ceggio_IC) / ncel_aff[12]
                end
                if abs(i - sez_aff[13]) <= semi_aff[13]
                    Q_aff[i] += (Q_levico_IC) / ncel_aff[13]
                end
                if abs(i - sez_aff[14]) <= semi_aff[14]
                    Q_aff[i] += (Q_centa_IC) / ncel_aff[14]
                end
                if abs(i - sez_aff[15]) <= semi_aff[15]
                    Q_aff[i] += (Q_chieppena_IC) / ncel_aff[15]
                end
                if abs(i - sez_aff[16]) <= semi_aff[16]
                    Q_aff[i] += (Q_grigno_IC) / ncel_aff[16]
                end
                if abs(i - sez_aff[17]) <= semi_aff[17]
                    Q_aff[i] += (Q_maso_IC) / ncel_aff[17]
                end
                if abs(i - sez_aff[18]) <= semi_aff[18]
                    Q_aff[i] += (Q_moggio_IC) / ncel_aff[18]
                end
                rhs[i] += dt/dx_int[i] * r_IC * Q_aff[i]
            end
        end
        ###  FINE INSERIMENTO AFFLUENTI CONDIZIONI INIZIALI

    
        ####  INSERIMENTO AFFLUENTI SIMULAZIONE
        # Stesse sezioni delle condizioni iniziali (sez_aff)
        if IC == 0
            @inbounds @simd for i in 1:imax
                if abs(i - sez_aff[1]) <= semi_aff[1]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta1), tk_brenta1, Qk_brenta1, time, dt)) / ncel_aff[1]
                end
                if abs(i - sez_aff[2]) <= semi_aff[2]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta2), tk_brenta2, Qk_brenta2, time, dt)) / ncel_aff[2]
                end
                if abs(i - sez_aff[3]) <= semi_aff[3]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta3), tk_brenta3, Qk_brenta3, time, dt)) / ncel_aff[3]
                end
                if abs(i - sez_aff[4]) <= semi_aff[4]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta4), tk_brenta4, Qk_brenta4, time, dt)) / ncel_aff[4]
                end
                if abs(i - sez_aff[5]) <= semi_aff[5]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta5), tk_brenta5, Qk_brenta5, time, dt)) / ncel_aff[5]
                end
                if abs(i - sez_aff[6]) <= semi_aff[6]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta6), tk_brenta6, Qk_brenta6, time, dt)) / ncel_aff[6]
                end
                # brenta7: da riattivare quando Brenta7.csv viene aggiunto (aumentare n_aff)
                if abs(i - sez_aff[7]) <= semi_aff[7]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta_new), tk_brenta_new, Qk_brenta_new, time, dt)) / ncel_aff[7]
                end
                if abs(i - sez_aff[8]) <= semi_aff[8]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta_v1), tk_brenta_v1, Qk_brenta_v1, time, dt)) / ncel_aff[8]
                end
                if abs(i - sez_aff[9]) <= semi_aff[9]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta_v2), tk_brenta_v2, Qk_brenta_v2, time, dt)) / ncel_aff[9]
                end
                if abs(i - sez_aff[10]) <= semi_aff[10]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta_v3), tk_brenta_v3, Qk_brenta_v3, time, dt)) / ncel_aff[10]
                end
                if abs(i - sez_aff[11]) <= semi_aff[11]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_brenta_v4), tk_brenta_v4, Qk_brenta_v4, time, dt)) / ncel_aff[11]
                end
                if abs(i - sez_aff[12]) <= semi_aff[12]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_ceggio), tk_ceggio, Qk_ceggio, time, dt)) / ncel_aff[12]
                end
                if abs(i - sez_aff[13]) <= semi_aff[13]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_levico), tk_levico, Qk_levico, time, dt)) / ncel_aff[13]
                end
                if abs(i - sez_aff[14]) <= semi_aff[14]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_centa), tk_centa, Qk_centa, time, dt)) / ncel_aff[14]
                end
                if abs(i - sez_aff[15]) <= semi_aff[15]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_chieppena), tk_chieppena, Qk_chieppena, time, dt)) / ncel_aff[15]
                end
                if abs(i - sez_aff[16]) <= semi_aff[16]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_grigno), tk_grigno, Qk_grigno, time, dt)) / ncel_aff[16]
                end
                if abs(i - sez_aff[17]) <= semi_aff[17]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_maso), tk_maso, Qk_maso, time, dt)) / ncel_aff[17]
                end
                if abs(i - sez_aff[18]) <= semi_aff[18]
                    Q_aff[i] += (LAGRANGE_INTERPOLATION!(length(tk_moggio), tk_moggio, Qk_moggio, time, dt)) / ncel_aff[18]
                end
                rhs[i] += dt/dx_int[i] * Q_aff[i]
            end
        end
        ##  FINE INSERIMENTO AFFLUENTI

        ####  PORTATA USCENTE PER ESONDAZIONE
        # Termine sorgente nella sola continuità, come gli affluenti ma con segno opposto: un contributo
        # per ogni argine superato, con lo sfioro lungo tutta la cella (dx_int) e la quota eta del passo precedente
        @inbounds for i in 1:imax
            Q_lat[i] = (FREE_OVERFLOW(eta[i], argine_dx[i], Cd, g) + FREE_OVERFLOW(eta[i], argine_sx[i], Cd, g))*dx_int[i]
            rhs[i]  -= dt/dx_int[i] * Q_lat[i]
        end
        ####  FINE PORTATA USCENTE PER ESONDAZIONE


        CHECK_NAN("termine noto e affluenti", nCS, time, IC; rhs, aa, cc)

        # --- Ciclo Newton ---
        #@infiltrate
        @inbounds for iNewton in 1:100
            @inbounds for i in 1:imax
                sec = READ_VALUES(data_depth, ys, first_row, last_row, i, eta[i])
                A_bn[i] = sec.A
                dA[i]   = sec.w

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
            CHECK_NAN("Thomas, iterazione di Newton $iNewton", nCS, time, IC; ff, di, deta)
            eta .= eta .- deta
        end
        #@infiltrate
        # --- Aggiornamento portate finali ---
        q[1] = qstar[1]
        q[imax+1] = qstar[imax+1]
        #Threads.@threads :static
        @inbounds @simd for i in 2:imax 
    
            q[i] = (qstar[i] - g*A_int[i]*dt/dx_b[i-1]*(eta[i]-eta[i-1]))/(1.0 + dt*γ[i])
          
        end

        CHECK_NAN("aggiornamento delle portate", nCS, time, IC; q)

        # --- Aggiornamento tempo e counter ---
        time += dt
        nCS  += 1
        Q_aff_tot = sum(Q_aff)  # Somma degli affluenti, per il bilancio
        Q_lat_tot = sum(Q_lat)  # Somma delle portate uscenti per esondazione, per il bilancio
        Q_aff .= 0.0  # Reset portate affluenti
        # --- Salvataggio risultati ---
        if nCS % 500 == 0
            println("Time step: ", nCS, " Time: ", @sprintf("%.4f", time), " dt: ", @sprintf("%.6f", dt))
        end
        

        #@infiltrate
        if IC == 1
            save_Q!(out.IC_Q, q)
            save_time!(out.IC_time, time)
            save_eta!(out.IC_eta, eta)
            save_h!(out.IC_h, eta .- b_min)
        else
            save_Q!(out.Q, q)
            save_time!(out.time, time)
            save_eta!(out.eta, eta)
            save_h!(out.h, eta .- b_min)
        end


        # --- Aggiornamento geometrie finali ---
        @inbounds for i in 1:imax
            sec = READ_VALUES(data_depth, ys, first_row, last_row, i, eta[i])
            A_b[i]     = sec.A
            β_b[i]     = sec.β
            den_γ_b[i] = sec.den
            b[i]       = sec.w
            rh_b[i]    = sec.rh
        end

        if IC == 1
            IC_counter = 0
            #Threads.@threads :static
            @inbounds @simd for i in 1:imax
                if abs(q[i] - qn[i]) < 1e-2
                  IC_counter += 1
                end
            end

            # Stazionario: rampa finita, q ferma e ingressi uguali all'uscita
            Q_in_tot = q[1] + Q_aff_tot - Q_lat_tot
            if IC_counter == imax && time >= T_riempimento && abs(Q_in_tot - q[imax+1]) < 1e-3*Q_in_tot
                println("Condizioni iniziali raggiunte al tempo: ", nCS, ", time = ", time)
                IC           = 0 # Disattiva condizioni iniziali
                nCS          = 1
                time         = 0.0
                consec_count = 0
                #@infiltrate
            end
        end 
    end
end
