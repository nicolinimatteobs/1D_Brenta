@views @inbounds function SEMI_IMPLICIT_BRENTA_2D_1D!()

    # --- INIZIALIZZAZIONE STRUTTURE ---
    params = INITIAL_SETTINGS!()

    # --- CREAZIONE ALIAS PER ACCEDERE AI CAMPI DELLE STRUCT ---
    @unpack imax, n_col, n_row, Ks, dtmax, Nmax, con_affluenti, scambio, V_soglia, T_piana, tau_piana, T_lento, piana_attiva,
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
    # con_affluenti (riga 27 di settings.txt): se vero i 18 affluenti entrano nello stato iniziale,
    # nel riempimento e nella simulazione
    Qin_idro       = Qk_cald[1]  # Portata in ingresso iniziale
    Qin_IC_target  = Qin_idro    # Portata di monte durante il riempimento (prima dell'idrogramma)
    y_ini          = 0.0      # Tirante minimo iniziale [m] (si parte dal moto uniforme, ma non sotto y_ini)
    dt_riempimento = 5.0      # Passo massimo durante il riempimento [s]
    T_assestamento = 3600.0   # Durata minima del riempimento [s]
    T_fermo        = 1800.0   # Il livello deve restare fermo per almeno questo tempo [s] ...
    tol_fermo      = 1e-3     # ... cioè cambiare meno di tol_fermo [m] in ogni sezione

    # Idrogrammi degli affluenti, nell'ordine delle sezioni di ingresso (da monte a valle)
    tk_aff = (tk_brenta1, tk_brenta2, tk_brenta3, tk_brenta4, tk_brenta5, tk_brenta6,
              tk_brenta_new, tk_brenta_v1, tk_brenta_v2, tk_brenta_v3, tk_brenta_v4,
              tk_ceggio, tk_levico, tk_centa, tk_chieppena, tk_grigno, tk_maso, tk_moggio)
    Qk_aff = (Qk_brenta1, Qk_brenta2, Qk_brenta3, Qk_brenta4, Qk_brenta5, Qk_brenta6,
              Qk_brenta_new, Qk_brenta_v1, Qk_brenta_v2, Qk_brenta_v3, Qk_brenta_v4,
              Qk_ceggio, Qk_levico, Qk_centa, Qk_chieppena, Qk_grigno, Qk_maso, Qk_moggio)
    Q_aff_IC = [Qk[1] for Qk in Qk_aff]  # Portate degli affluenti nelle condizioni iniziali
    Q_aff_t  = zeros(length(Qk_aff))     # Portate degli affluenti al passo corrente (simulazione)

    # Sezioni di ingresso degli affluenti, distribuite equamente lungo il canale
    n_aff   = 18
    sez_aff = [round(Int, k*imax/(n_aff + 1)) for k in 1:n_aff]

    # Gli affluenti con picco oltre Q_soglia_aff entrano su 3 sezioni (sez_aff e le due vicine),
    # con la portata divisa in parti uguali; gli altri entrano in una sola sezione
    Q_soglia_aff = 10.0  # [m3/s]
    Qmax_aff = maximum.(Qk_aff)
    semi_aff = [Q >= Q_soglia_aff ? 1 : 0 for Q in Qmax_aff]  # sezioni in più per lato
    ncel_aff = 2 .* semi_aff .+ 1                              # sezioni su cui entra ogni affluente


    exp_dt        = 2.0/3.0 # Esponente di dt nel calcolo della portata
    nCS           = 1
    IC_counter    = 0
    IC            = 1       # Indicatore condizioni iniziali attive
    consec_count  = 0
    g             = 9.81


    # Stato iniziale dalle portate cumulate: la portata alle interfacce è quella di monte più gli
    # affluenti entrati a monte (se con_affluenti), e in ogni sezione il tirante è quello di moto
    # uniforme per la portata media della cella (almeno y_ini). L'attrito è quello delle tabelle (jden),
    # lo stesso usato dopo in moto vario. Gli affluenti entrano subito al primo valore, senza rampa:
    # il riempimento serve solo all'assestamento (rigurgiti nei tratti non prismatici, valle degli affluenti)
    Q_aff = zeros(imax)  # Vettore portate affluenti (per cella)
    con_affluenti && TRIBUTARY_INFLOW!(Q_aff, Q_aff_IC, sez_aff, semi_aff, ncel_aff)
    q[1] = Qin_IC_target
    @inbounds for i in 1:imax
        q[i+1] = q[i] + Q_aff[i]
    end
    @inbounds for i in 1:imax
        eta[i] = UNIFORM_LEVEL(data_depth, ys, first_row, last_row, i, 0.5*(q[i] + q[i+1]), bx[i]; y_min = y_ini)
    end
    Q_aff .= 0.0
    # Alternativa: velo d'acqua costante di y_ini in tutte le sezioni
    #eta .= b_min .+ y_ini

    Q_lat = zeros(imax)  # Scambio laterale con la piana per cella [m3/s]: positivo se entra, negativo se esce
    Cd    = 0.385        # Coefficiente di efflusso dello sfioro sugli argini (scambio = :sfioro)
    # Serbatoio virtuale della piana (scambio = :tabella). V_soglia [m3], T_piana [s], tau_piana [s] e
    # T_lento [s] vengono dalle righe 29-32 di settings.txt
    V_piana     = 0.0                 # Volume uscito e non ancora rientrato [m3]
    V_rientrato = 0.0                 # Volume rientrato dall'inizio [m3]
    Vu_tot      = 0.0                 # Volume uscito dall'inizio [m3]
    t_tot       = 0.0                 # Tempo dall'inizio, riempimento compreso (time riparte da 0) [s]
    t_storia    = [0.0]; Vu_storia = [0.0]   # Storia del volume uscito, per il ritardo
    # Rilascio lento: la parte sotto la soglia torna al fiume in T_lento, nelle sezioni di rientro, con pesi
    # proporzionali al ΔQ massimo positivo della tabella di ogni sezione (somma 1)
    w_rientro = [max(0.0, maximum(data_depth[13, first_row[k]:first_row[k+1]-1])) for k in 1:imax]
    w_rientro ./= max(sum(w_rientro), 1e-12)

    # Bilancio di massa della simulazione (dopo il riempimento): volumi entrati da monte, usciti a valle,
    # entrati dagli affluenti e scambiati con la piana, e volume in alveo all'inizio [m3]
    B_monte = 0.0; B_valle = 0.0; B_aff = 0.0; B_lat = 0.0; B_iniziale = 0.0
    B_lat_tutto = 0.0   # scambio netto con la piana dall'inizio, riempimento compreso (per il serbatoio)
    eta_fermo = copy(eta)  # Livello all'inizio della finestra di controllo del riempimento
    t_fermo   = time       # Inizio della finestra
    fermo     = false      # true se nell'ultima finestra il livello è rimasto fermo

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
            # Alternativa con la media (upwind sull'area solo durante il riempimento: una cella non cede
            # più acqua di quanta ne ha). Per riattivarla sostituire le righe di A_int e den_γ_int con:
            #A_int[i]     = IC == 1 ? (q[i] >= 0.0 ? A_b[i-1] : A_b[i]) : 0.5 * (A_b[i] + A_b[i-1])
            #den_γ_int[i] = 0.5 * (den_γ_b[i] + den_γ_b[i-1])

            uf[i]   = q[i] * A_int[i]/(A_int[i]^2 + tol)
            smax[i] = abs(β_int[i] * uf[i])

        end

        amax = maximum(smax)
        dt_lim = IC == 1 ? dt_riempimento : dtmax
        dt   = min(dt_lim, CFL*dx/(amax + tol))
        CHECK_NAN("calcolo di dt", nCS, time, IC; dt, uf, smax)

        #@infiltrate

        if time + dt >= tend
            dt = tend - time
        end
        if time >= tend
            break
        end

        #@infiltrate
        @inbounds @simd for i in 1:imax+1
            if i == 1
                if IC == 1
                    qstar[i] = Qin_IC_target
                else
                    Qin_idro = LAGRANGE_INTERPOLATION!(nr_hydro, tk_cald, Qk_cald, time, dt)
                    qstar[i] = Qin_idro
                end
            elseif i == imax+1
                # Uscita: moto uniforme con l'attrito della tabella dell'ultima sezione (jden) e la sua
                # pendenza locale. In calcolo_tabelle_daBM il jden dell'ultima sezione è ricavato proprio con
                # questa pendenza, quindi la portata in uscita segue la scala di deflusso degli scalini 2D.
                # Prima si usava Manning con il Ks di settings.txt (riga 23, ora non più usato):
                #uout = Ks*sqrt(abs(bx[imax]))*(rh_int[i]+tol)^exp_dt
                #qstar[i] = A_int[i]*uout
                qstar[i] = UNIFORM_DISCHARGE(A_int[i], den_γ_int[i], bx[imax])
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

        ####  INSERIMENTO AFFLUENTI (se con_affluenti)
        # Ogni affluente entra nella sezione sez_aff e, se ha portata alta, nelle due vicine, solo come
        # sorgente nella continuità: nel riempimento con la prima portata dell'idrogramma, in simulazione
        # con quella interpolata al tempo time + dt
        if con_affluenti
            if IC == 1
                TRIBUTARY_INFLOW!(Q_aff, Q_aff_IC, sez_aff, semi_aff, ncel_aff)
            else
                for j in eachindex(Q_aff_t)
                    Q_aff_t[j] = LAGRANGE_INTERPOLATION!(length(tk_aff[j]), tk_aff[j], Qk_aff[j], time, dt)
                end
                TRIBUTARY_INFLOW!(Q_aff, Q_aff_t, sez_aff, semi_aff, ncel_aff)
            end
            @inbounds for i in 1:imax
                rhs[i] += dt/dx_int[i] * Q_aff[i]
            end
        end
        ####  FINE INSERIMENTO AFFLUENTI

        ####  SCAMBIO LATERALE CON LA PIANA (riga 28 di settings.txt)
        # Termine sorgente nella sola continuità, come gli affluenti, esplicito (quota eta del passo
        # precedente). Q_lat è positivo se l'acqua entra, negativo se esce.
        # - sfioro: sfioro libero, solo uscente, un contributo per ogni argine superato, lungo tutta la cella (dx_int);
        # - tabella: ΔQ della tabella della sezione (portata della sezione meno quella della precedente negli
        #   scalini 2D), entrante o uscente; sostituisce lo sfioro, che vi è già compreso. Le uscite sono
        #   quelle della tabella e finiscono in un serbatoio virtuale unico per tutto il tratto (la piana).
        #   I rientri sono quelli della tabella, ma limitati dal serbatoio e ridotti tutti della stessa
        #   frazione: senza limite crescono con il livello, che cresce con i rientri, e il modello crea acqua.
        #   Parametri (righe 29-32 di settings.txt):
        #     V_soglia  la piana trattiene i primi V_soglia m3 (si riempie): rientra solo il volume oltre;
        #     T_piana   tempo di svuotamento: rientro al massimo (V − V_soglia)/T_piana (0: tutto subito);
        #     tau_piana ritardo di percorrenza: l'acqua uscita arriva al serbatoio dopo tau_piana;
        #     T_lento   rilascio lento della parte sotto la soglia, in T_lento, sulle sezioni di rientro
        #               (anche quando la tabella non chiede rientri, es. nel ramo calante; 0: nessun rilascio);
        #   Riga 33, serbatoio attivo (piana_attiva): la parte oltre la soglia rientra comunque al ritmo Q_disp;
        #   quello che la tabella non chiede (ramo calante) va sulle sezioni di rientro con i pesi w_rientro.
        #   Senza, dopo il colmo la tabella smette di chiedere rientri e l'acqua oltre la soglia resta ferma;
        # - nessuno: nessuno scambio.
        if scambio == :sfioro
            @inbounds for i in 1:imax
                Q_lat[i] = -(FREE_OVERFLOW(eta[i], argine_dx[i], Cd, g) + FREE_OVERFLOW(eta[i], argine_sx[i], Cd, g))*dx_int[i]
                rhs[i] += dt/dx_int[i] * Q_lat[i]
            end
        elseif scambio == :tabella
            Q_esce   = 0.0  # portata uscente totale [m3/s]
            Q_rientro = 0.0  # rientro richiesto dalle tabelle [m3/s]
            @inbounds for i in 1:imax
                Q_lat[i] = READ_EXCHANGE(data_depth, ys, first_row, last_row, i, eta[i])
                Q_lat[i] < 0.0 ? (Q_esce -= Q_lat[i]) : (Q_rientro += Q_lat[i])
            end
            # Portata che il serbatoio può restituire in questo passo: volume oltre la soglia diviso per
            # T_piana, e comunque non più dell'acqua disponibile
            if tau_piana > 0.0
                # con il ritardo il volume disponibile è quello uscito fino a tau_piana prima della fine
                # del passo (storia del volume uscito, t_tot conta anche il riempimento)
                Vu_tot += dt*Q_esce
                push!(t_storia, t_tot + dt); push!(Vu_storia, Vu_tot)
                V_arrivato = VOLUME_AT(t_storia, Vu_storia, t_tot + dt - tau_piana)
                V_oltre    = max(V_arrivato - V_rientrato - V_soglia, 0.0)
                Q_disp     = T_piana > 0.0 ? min(V_oltre/T_piana, V_oltre/dt) : V_oltre/dt
                V_sotto    = clamp(V_arrivato - V_rientrato, 0.0, V_soglia)
            else
                V_oltre   = max(V_piana - V_soglia, 0.0)
                Q_disp    = T_piana > 0.0 ? min(V_oltre/T_piana, V_oltre/dt + Q_esce) : V_oltre/dt + Q_esce
                V_sotto   = clamp(V_piana, 0.0, V_soglia)
            end
            f_rientro = Q_rientro > 0.0 ? min(1.0, Q_disp/Q_rientro) : 0.0
            # rilascio lento della parte sotto la soglia: V_sotto/T_lento, sulle sezioni di rientro (w_rientro)
            Q_lento = T_lento > 0.0 ? min(V_sotto/T_lento, V_sotto/dt) : 0.0
            # serbatoio attivo: la parte di Q_disp che la tabella non chiede rientra sulle sezioni di rientro
            Q_extra = piana_attiva ? max(Q_disp - f_rientro*Q_rientro, 0.0) : 0.0
            @inbounds for i in 1:imax
                Q_lat[i] > 0.0 && (Q_lat[i] *= f_rientro)
                T_lento > 0.0 && (Q_lat[i] += Q_lento*w_rientro[i])
                piana_attiva  && (Q_lat[i] += Q_extra*w_rientro[i])
                rhs[i] += dt/dx_int[i] * Q_lat[i]
            end
            # scritto così: con Q_lento = 0 e Q_extra = 0 i conti sono quelli senza rilascio lento e serbatoio attivo
            V_piana     += dt*(Q_esce - f_rientro*Q_rientro - Q_lento - Q_extra)
            V_rientrato += dt*f_rientro*Q_rientro + dt*Q_lento + dt*Q_extra
        end
        ####  FINE SCAMBIO LATERALE


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
        t_tot += dt
        nCS  += 1
        Q_aff_tot = sum(Q_aff)  # Somma degli affluenti, per il bilancio
        Q_lat_tot = sum(Q_lat)  # Scambio laterale netto (positivo se entra), per il bilancio
        Q_aff .= 0.0  # Reset portate affluenti
        # Bilancio di massa: volumi della simulazione (dopo il riempimento) e scambio con la piana di tutto il run
        B_lat_tutto += dt*Q_lat_tot
        if IC == 0
            B_monte += dt*q[1]; B_valle += dt*q[imax+1]; B_aff += dt*Q_aff_tot; B_lat += dt*Q_lat_tot
        end
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
            save_Q!(out.IC_Q_lat, Q_lat)
        else
            save_Q!(out.Q, q)
            save_time!(out.time, time)
            save_eta!(out.eta, eta)
            save_h!(out.h, eta .- b_min)
            save_Q!(out.Q_lat, Q_lat)
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

            # Livello fermo: ogni T_fermo si confronta il livello con quello di inizio finestra
            if time - t_fermo >= T_fermo
                fermo = true
                @inbounds for i in 1:imax
                    abs(eta[i] - eta_fermo[i]) >= tol_fermo && (fermo = false)
                end
                eta_fermo .= eta
                t_fermo    = time
            end

            # Stazionario: durata minima passata, q ferma, livello fermo nell'ultima finestra e ingressi uguali all'uscita
            Q_in_tot = q[1] + Q_aff_tot + Q_lat_tot
            if IC_counter == imax && time >= T_assestamento && fermo && abs(Q_in_tot - q[imax+1]) < 1e-3*Q_in_tot
                println("Condizioni iniziali raggiunte al tempo: ", nCS, ", time = ", time)
                B_iniziale   = sum(A_b .* dx_int)  # volume in alveo all'inizio della simulazione
                IC           = 0 # Disattiva condizioni iniziali
                nCS          = 1
                time         = 0.0
                consec_count = 0
                #@infiltrate
            end
        end
    end

    # --- CONTROLLO DEL BILANCIO DI MASSA ---
    if IC == 1
        println("Bilancio di massa non calcolato: il riempimento non si è concluso")
    else
        BILANCIO_DI_MASSA(joinpath(output_dir, "bilancio.txt"), scambio, B_iniziale, sum(A_b .* dx_int),
                          B_monte, B_valle, B_aff, B_lat, V_piana, B_lat_tutto)
    end
end
