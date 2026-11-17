# Approfondimenti 15-21: sfioro, scambio dalla tabella, serbatoio della piana, fine del riempimento, bilancio di massa,
# file di uscita, costruzione delle tabelle.
import numpy as np
from stile import *

SEMI = "Source/Semi_Implicit_Brenta_2D_1D.jl"
SUB = "Source/subroutines.jl"
WRI = "Source/writing.jl"
CT = "calcolo_tabelle_daBM/src/main.jl"


def tabella_sezione(k):
    T = np.loadtxt(os.path.join(PROGETTO, TABELLE), delimiter=",", skiprows=1)
    return T[T[:, 0] == k]


def uscita(scen, nome, ic=False):
    return np.loadtxt(os.path.join(PROGETTO, USCITE, scen, ("IC_" if ic else "") + nome), ndmin=1)


# =============================================================================================== 15
def p15(D):
    G = np.loadtxt(os.path.join(PROGETTO, GEOMETRIA), delimiter=",", skiprows=1)
    D.h1("A cosa serve")
    D.p("Quando il livello supera la cresta di un argine, l'acqua esce dal canale verso la piana. Con lo scambio "
        "`sfioro` (riga 28 di `settings.txt`) il modello lo rappresenta come uno **stramazzo a parete grossa** lungo "
        "ciascun argine della cella. L'acqua uscita non rientra: è persa per il fiume.")
    D.h1("1. La formula dello sfioro")
    D.p("Su una soglia larga, con carico Y_{0} sopra la cresta, la corrente passa per lo stato critico: tirante "
        "critico 2/3 del carico, velocità √(g·2Y_{0}/3). La portata per unità di lunghezza di argine è quindi")
    D.eq([r"q_l=\frac{2}{3}Y_0\,\sqrt{g\,\frac{2}{3}Y_0}=\left(\frac{2}{3}\right)^{3/2}\sqrt{g}\;Y_0^{3/2}\approx 0.385\,\sqrt{2g\,Y_0^3}",
          r"Y_0=\max\left(\eta-z_{argine},\ 0\right),\qquad C_d=0.385"])
    D.p("È la formula 16 di Siviglia, Stocchino e Colombini (ASCE 2009), in `reference/`. Il coefficiente 0,385 è "
        "proprio (2/3)^{3/2}/√2: sfioro libero, senza rigurgito da valle.")
    def graf(fig):
        ax = fig.add_subplot(111); Y = np.linspace(0, 0.6, 200)
        ax.plot(Y * 100, 0.385 * np.sqrt(2 * 9.81 * Y ** 3), color="#2a78d6", lw=1.6, label="q_l per metro di argine")
        ax.plot(Y * 100, 0.385 * np.sqrt(2 * 9.81 * Y ** 3) * 50, color="#eb6834", lw=1.4, ls="--", label="× 50 m di cella")
        ax.set_xlabel("carico sopra la cresta Y0 [cm]", fontsize=8); ax.set_ylabel("portata [m²/s] e [m³/s]", fontsize=8)
        ax.grid(alpha=0.25); ax.tick_params(labelsize=7); ax.legend(fontsize=7.5)
    D.grafico(graf, "sfioro", "Portata sfiorata in funzione del carico: con 20 cm sopra un argine, una cella di 50 m "
              "perde circa 7,6 m³/s da quel lato.")
    D.h1("2. Dalla formula alla cella")
    D.p("Ogni argine superato contribuisce per tutta la lunghezza della cella, con lo stesso carico:")
    D.eq(r"Q_{lat,i}=-\left[q_l(\eta_i,\,z_{dx,i})+q_l(\eta_i,\,z_{sx,i})\right]\Delta x_i\qquad (<0:\ \mathrm{esce})")
    D.p("Il termine entra solo nella continuità, come gli affluenti ma con segno opposto, ed è **esplicito**: si usa "
        "il livello del passo precedente. È attivo anche nel riempimento e compare nel suo bilancio (l'uscita deve "
        "pareggiare monte + affluenti − sfioro).")
    D.estratto(SUB, "@inline function FREE_OVERFLOW", n=4, titolo=SUB)
    D.estratto(SEMI, "        if scambio == :sfioro", n=5, titolo=SEMI)
    D.p(f"Le quote degli argini vengono da `bankpoints.csv` (un punto del DTM per argine e per sezione). Nel "
        f"Test_Brenta vanno da {min(G[:,4].min(), G[:,5].min()):.2f} a {max(G[:,4].max(), G[:,5].max()):.2f} m s.l.m.; "
        "una geometria a 4 colonne ha argini a `Inf` e non sfiora mai.")
    D.flowchart(dict(dimensioni=(12, 7.8), w=6.2, h=0.95, nodi={
        "a": (6, 7.2, "scambio = sfioro: per ogni cella i", "inizio", 5.6, 0.75),
        "b": (6, 6.0, "Y0_dx = max(η_i − z_dx, 0),  Y0_sx = max(η_i − z_sx, 0)", "proc"),
        "c": (6, 4.8, "q = Cd·√(2g·Y0³) per i due argini (FREE_OVERFLOW)", "funz"),
        "d": (6, 3.6, "Q_lat[i] = −(q_dx + q_sx)·Δx_i", "proc"),
        "e": (6, 2.4, "rhs[i] += Δt/Δx_i · Q_lat[i]", "proc"),
        "f": (6, 1.2, "Q_lat salvato in Q_lat.csv (negativo)", "io", 5.4)},
        archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("e", "f")]),
        didascalia="Flowchart dello sfioro sugli argini.", larghezza=0.8)
    D.h1("3. Limiti")
    D.elenco(["**solo uscita**: l'acqua che nel 2D rientra più a valle nel 1D è persa. Nei fittizi a 50 m il colmo in "
              "uscita si fermava a circa 60 m³/s contro 64-82 del 2D;",
              "**un punto per argine**: la quota della cresta è un solo punto del DTM applicato a 50 m di argine;",
              "**sfioro libero**: non sa che cosa c'è dall'altra parte. Dove la piana è già allagata, un punto basso "
              "dell'argine è un punto di rientro, non di uscita (S3550 nei fittizi);",
              "per questo è stato introdotto lo scambio dalla tabella con il serbatoio della piana (approfondimenti 16-17), "
              "che è la configurazione attuale."])
    D.collegamenti([(16, "16_scambio_tabella.pdf", "Lo scambio dalla tabella ΔQ"), (17, "17_serbatoio_della_piana.pdf", "Il serbatoio della piana"),
                    (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa")])


# =============================================================================================== 16
def p16(D):
    s10, s83 = tabella_sezione(10), tabella_sezione(83)
    D.h1("A cosa serve")
    D.p("Lo scambio `tabella` sostituisce la formula dello sfioro con un'informazione presa direttamente dal modello 2D. "
        "Negli scalini 2D a moto permanente la portata cambia lungo il tratto **solo** per lo scambio con la piana: se "
        "la sezione k ha meno portata della k−1 l'acqua è uscita nel tratto, se ne ha di più è rientrata.")
    D.h1("1. La colonna ΔQ")
    D.eq(r"\Delta Q_k(\eta)=Q^{2D}_k-Q^{2D}_{k-1}\qquad \mathrm{allo\ stesso\ scalino\ di\ portata}")
    D.p("ΔQ < 0 dove l'acqua esce, ΔQ > 0 dove rientra, 0 nella prima sezione e dove una delle due è asciutta. È "
        "calcolata da `calcolo_tabelle_daBM` (approfondimento 21) e scritta nelle colonne 13-14 delle tabelle "
        "`S_50m_q2`, con la sua pendenza rispetto al livello, così si interpola come le altre grandezze. Non servono "
        "né la formula dello sfioro né le quote degli argini.")
    def graf(fig):
        ax = fig.subplots(1, 2)
        for a, s, nome, col in ((ax[0], s10, "S0500 (uscita)", "#eb6834"), (ax[1], s83, "S4150 (rientro)", "#1baf7a")):
            a.plot(s[:, 12], s[:, 1], "-o", ms=2.5, color=col, lw=1.3)
            a.axvline(0, color="#9aa0a6", lw=0.8); a.set_title(nome, fontsize=9)
            a.set_xlabel("ΔQ [m³/s]", fontsize=8); a.grid(alpha=0.25); a.tick_params(labelsize=7)
        ax[0].set_ylabel("livello η [m s.l.m.]", fontsize=8)
    D.grafico(graf, "deltaQ", f"ΔQ(η) nelle tabelle di due sezioni: in S0500 l'acqua esce fino a {s10[:,12].min():.1f} m³/s, "
              f"in S4150 rientra fino a {s83[:,12].max():.1f} m³/s.")
    D.h1("2. READ_EXCHANGE")
    D.p("Legge ΔQ al livello della sezione con la stessa ricerca di riga di `READ_VALUES` (approfondimento 04), ma con "
        "una differenza importante agli estremi: **oltre l'ultima quota usa l'ultimo valore, senza estrapolare**. ΔQ può "
        "cambiare molto tra le ultime righe e una retta prolungata darebbe scambi irrealistici.")
    D.eq(r"\Delta Q(\eta)=\Delta Q_1\ (\eta\leq\eta_1),\qquad \Delta Q_r+\Delta Q'_r\,(\eta-\eta_r)\ (\eta_r\leq\eta<\eta_{r+1}),\qquad \Delta Q_N\ (\eta\geq\eta_N)")
    D.estratto(SUB, "@inline function READ_EXCHANGE", n=11, titolo=SUB)
    D.h1("3. Da solo è instabile")
    D.p("La prima prova usava ΔQ così com'è, con Q_{lat,i} = ΔQ_{i}(η_{i}). Nel `fittizio_1` la portata a 4-5 km è arrivata a "
        "93-95 m³/s con 81 in ingresso e a 6,4 h il calcolo è divergito: nessun NaN, ma il passo temporale è andato a "
        "zero. Il motivo è un **anello di retroazione positivo** nelle sezioni di rientro.")
    D.p("Si consideri una sola cella di rientro, con larghezza W e lunghezza Δx, e una piccola perturbazione del "
        "livello η' attorno all'equilibrio. Dalla continuità:")
    D.eq([r"W\,\Delta x\,\frac{d\eta'}{dt}=\frac{d\Delta Q}{d\eta}\,\eta'+\ldots\quad\Rightarrow\quad \eta'(t)\propto e^{\lambda t},\qquad \lambda=\frac{1}{W\,\Delta x}\,\frac{d\Delta Q}{d\eta}"])
    D.p("Dove l'acqua rientra ΔQ **cresce** con il livello (dΔQ/dη > 0, come in S4150): λ > 0 e la perturbazione cresce "
        "esponenzialmente. Più acqua rientra, più il livello sale, più acqua rientra: il modello **crea acqua**, perché "
        "nulla lega il rientro all'acqua uscita davvero. Nel 2D invece può rientrare solo l'acqua che prima è uscita e "
        "che ha riempito la piana. Il termine di gravità del sistema dei livelli smorza in parte, ma con rientri di "
        "alcuni m³/s non basta.")
    D.idea("le uscite della tabella sono affidabili, i rientri no: devono essere **limitati dal volume disponibile**. "
           "Da qui il serbatoio della piana (approfondimento 17).")
    D.h1("4. Il primo passo del blocco tabella")
    D.p("A ogni passo si legge ΔQ in tutte le celle e si separano le uscite (che vanno al serbatoio) dai rientri "
        "richiesti (che il serbatoio dovrà autorizzare):")
    D.eq(r"Q_{esce}=\sum_i\max\left(-\Delta Q_i,\,0\right),\qquad Q_{rientro}=\sum_i\max\left(\Delta Q_i,\,0\right)")
    D.estratto(SEMI, "        elseif scambio == :tabella", n=7, titolo=SEMI)
    D.flowchart(dict(dimensioni=(13, 8.2), w=5.6, h=0.95, nodi={
        "a": (6.5, 7.6, "READ_EXCHANGE(k, η)", "inizio", 4.4, 0.75),
        "b": (6.5, 6.3, "η ≤ η_1 ?", "dec", 2.8, 1.0),
        "c": (11.2, 6.3, "return ΔQ_1", "fine", 2.6, 0.75),
        "d": (6.5, 4.9, "η ≥ η_N ?", "dec", 2.8, 1.0),
        "e": (11.2, 4.9, "return ΔQ_N\n(senza estrapolare)", "fine", 2.9, 0.9),
        "f": (6.5, 3.5, "FIND_ROW → riga r", "funz", 4.4),
        "g": (6.5, 2.2, "return ΔQ_r + ΔQ'_r·(η − η_r)", "fine", 5.6, 0.8)},
        archi=[("a", "b"), ("b", "c", "sì"), ("b", "d", "no"), ("d", "e", "sì"), ("d", "f", "no"), ("f", "g")]),
        didascalia="Flowchart di `READ_EXCHANGE`.", larghezza=0.8)
    D.collegamenti([(4, "04_tabelle_delle_sezioni.pdf", "Le tabelle"), (17, "17_serbatoio_della_piana.pdf", "Il serbatoio della piana"),
                    (21, "21_costruzione_tabelle.pdf", "Come si calcola ΔQ")])


# =============================================================================================== 17
def p17(D):
    D.h1("A cosa serve")
    D.p("Il serbatoio è il pezzo che rende stabile e fisicamente sensato lo scambio dalla tabella. Rappresenta **tutta "
        "la piana** del tratto come un unico volume virtuale V: le uscite della tabella lo riempiono, i rientri lo "
        "svuotano e non possono restituire più acqua di quella che contiene. Non ha geometria: solo un volume e "
        "quattro parametri più un interruttore (righe 29-33 di `settings.txt`).")
    D.tabella([["Riga", "Parametro", "Simb.", "Significato fisico", "Valore tarato"],
               ["29", "soglia", "V_{s}", "la piana deve riempirsi prima di restituire acqua (avvallamenti, ristagni)", "150 000 m³"],
               ["30", "tempo di svuotamento", "T", "con quale velocità la piana piena restituisce", "0,5 h"],
               ["31", "ritardo", "τ", "tempo per attraversare la piana dall'uscita al rientro", "1 h"],
               ["32", "rilascio lento", "T_{l}", "la parte sotto la soglia torna lentamente (drenaggio)", "32 h"],
               ["33", "serbatoio attivo", "–", "l'acqua oltre la soglia rientra anche quando la tabella non lo chiede", "si"]],
              [0.07, 0.17, 0.09, 0.5, 0.17])
    D.h1("1. Il volume arrivato e il ritardo")
    D.p("L'acqua uscita impiega un tempo τ per raggiungere il serbatoio. Il codice tiene la **storia** del volume uscito "
        "cumulato V_{u}(t) e guarda il valore di τ secondi prima. t conta anche il riempimento (`t_tot`):")
    D.eq([r"V_u(t)=\int_0^t Q_{esce}\,dt,\qquad V_{arr}(t)=V_u(t-\tau)",
          r"V_{disp}=V_{arr}-V_r\qquad (V_r:\ \mathrm{volume\ gi\grave{a}\ rientrato})"])
    D.p("`VOLUME_AT` interpola linearmente la storia (un punto per passo):")
    D.estratto(SUB, "function VOLUME_AT", n=6, titolo=SUB)
    D.h1("2. Quanto può rientrare")
    D.p("Il volume disponibile si divide in due parti, sopra e sotto la soglia:")
    D.eq([r"V_{oltre}=\max\left(V_{disp}-V_s,\ 0\right),\qquad V_{sotto}=\min\left(\max\left(V_{disp},\,0\right),\ V_s\right)",
          r"Q_{disp}=\min\left(\frac{V_{oltre}}{T},\ \frac{V_{oltre}}{\Delta t}\right),\qquad Q_{lento}=\min\left(\frac{V_{sotto}}{T_l},\ \frac{V_{sotto}}{\Delta t}\right)"])
    D.p("Il secondo termine di ogni min impedisce di restituire in un passo più di quanto c'è. Senza ritardo (τ = 0) "
        "si usa V_{piana} al posto di V_{disp} e la portata uscente del passo può rientrare subito "
        "(Q_{disp} ≤ V_{oltre}/Δt + Q_{esce}).")
    D.h1("3. Come si distribuisce")
    D.p("I **rientri richiesti** dalla tabella (le celle con ΔQ > 0) sono ridotti tutti della stessa frazione:")
    D.eq(r"f=\min\left(1,\ \frac{Q_{disp}}{Q_{rientro}}\right)\qquad (f=0\ \mathrm{se}\ Q_{rientro}=0)")
    D.p("Il **rilascio lento** e, con il serbatoio attivo, la parte di Q_{disp} che la tabella non chiede (tipicamente "
        "nel ramo calante, quando i livelli scendono e ΔQ torna negativo) vanno nelle sezioni di rientro con pesi "
        "fissi, proporzionali al massimo ΔQ positivo della tabella di ogni sezione:")
    D.eq([r"Q_{extra}=\max\left(Q_{disp}-f\,Q_{rientro},\ 0\right)\quad (\mathrm{solo\ con\ serbatoio\ attivo})",
          r"w_i=\frac{\max\left(0,\ \max_\eta\Delta Q_i\right)}{\sum_k\max\left(0,\ \max_\eta\Delta Q_k\right)},\qquad \sum_iw_i=1"])
    D.p("Lo scambio finale della cella i, che entra nel termine noto come gli affluenti, è")
    D.eq(r"Q_{lat,i}=\Delta Q_i\ (\Delta Q_i<0)\ \mathrm{oppure}\ f\,\Delta Q_i\ (\Delta Q_i>0),\quad +\ \left(Q_{lento}+Q_{extra}\right)w_i")
    D.h1("4. Il serbatoio conserva il volume")
    D.eq([r"V_{piana}^{n+1}=V_{piana}^{n}+\Delta t\left(Q_{esce}-f\,Q_{rientro}-Q_{lento}-Q_{extra}\right)",
          r"V_r^{n+1}=V_r^{n}+\Delta t\left(f\,Q_{rientro}+Q_{lento}+Q_{extra}\right)"])
    D.p("Sommando Q_{lat,i} su tutte le celle si ottiene −Q_{esce} + fQ_{rientro} + Q_{lento} + Q_{extra} (perché Σw_{i} = 1): "
        "esattamente l'opposto della variazione del serbatoio. **Tutto ciò che esce dal canale entra nel serbatoio e "
        "viceversa**: il volume totale canale + piana si conserva, e il rientro non può mai superare l'acqua uscita. "
        "`BILANCIO_DI_MASSA` lo verifica a fine simulazione (approfondimento 19).")
    D.p("Non può rientrare più acqua di quella uscita: V_{r} ≤ V_{arr} ≤ V_{u}, perché Q_{disp} e Q_{lento} sono "
        "limitati da V_{oltre} e V_{sotto}, che insieme valgono al massimo V_{disp}. È questo vincolo che toglie "
        "l'instabilità dell'approfondimento 16.")
    D.h1("5. Il codice")
    D.estratto(SEMI, "            # Portata che il serbatoio può restituire", "            V_rientrato += dt*f_rientro*Q_rientro", titolo=SEMI)
    D.p("I pesi si calcolano una volta all'inizio:")
    D.estratto(SEMI, "    w_rientro = [max(0.0,", n=2, titolo=SEMI)
    D.flowchart(dict(dimensioni=(15, 17.0), w=6.4, h=0.95, nodi={
        "a": (6.5, 16.4, "scambio = tabella (ogni passo)", "inizio", 5.4, 0.75),
        "b": (6.5, 15.2, "READ_EXCHANGE in ogni cella → ΔQ_i", "funz"),
        "c": (6.5, 13.9, "Q_esce = Σ uscite,  Q_rientro = Σ rientri richiesti", "proc"),
        "d": (6.5, 12.6, "τ > 0 ?", "dec", 2.4, 1.0),
        "e": (2.6, 11.1, "V_u += Δt·Q_esce, storia\nV_arr = VOLUME_AT(t + Δt − τ)\nV_disp = V_arr − V_r", "funz", 5.0, 1.4),
        "f": (10.6, 11.1, "V_disp = V_piana\n(il Q_esce del passo può\nrientrare subito)", "proc", 4.6, 1.4),
        "g": (6.5, 9.4, "V_oltre = max(V_disp − V_s, 0) → Q_disp = min(V_oltre/T, V_oltre/Δt)\nV_sotto = min(max(V_disp, 0), V_s) → Q_lento = min(V_sotto/T_l, V_sotto/Δt)", "proc", 11.6, 1.1),
        "h": (6.5, 7.9, "f = min(1, Q_disp/Q_rientro)", "proc"),
        "i": (6.5, 6.6, "serbatoio attivo ?", "dec", 3.6, 1.0),
        "l": (12.0, 6.6, "Q_extra = max(Q_disp − f·Q_rientro, 0)", "proc", 5.0),
        "m": (6.5, 5.1, "Q_lat,i = ΔQ_i (uscite),  f·ΔQ_i (rientri)\n+ (Q_lento + Q_extra)·w_i", "proc", 6.8, 1.1),
        "n": (6.5, 3.7, "rhs[i] += Δt/Δx_i · Q_lat,i", "proc"),
        "o": (6.5, 2.4, "V_piana += Δt·(Q_esce − f·Q_rientro − Q_lento − Q_extra)\nV_r += Δt·(f·Q_rientro + Q_lento + Q_extra)", "proc", 8.4, 1.1),
        "p": (6.5, 1.1, "fine dello scambio", "fine", 3.6, 0.7)},
        archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d.w", "e.n", "sì", [(2.6, 12.6)]), ("d.e", "f.n", "no", [(10.6, 12.6)]),
               ("e.s", "g.n", None, [(2.6, 10.2), (6.5, 10.2)]), ("f.s", "g.n", None, [(10.6, 10.2), (6.5, 10.2)]),
               ("g", "h"), ("h", "i"), ("i", "l", "sì"), ("i", "m", "no (Q_extra = 0)"), ("l.s", "m.e", None, [(12.0, 5.1)]),
               ("m", "n"), ("n", "o"), ("o", "p")]),
        didascalia="Flowchart dello scambio `tabella` con il serbatoio della piana.")
    t = uscita("fittizio_4", "time.csv"); L = uscita("fittizio_4", "Q_lat.csv"); Qin = uscita("fittizio_4", "Q.csv")[:, 0]
    Li = uscita("fittizio_4", "Q_lat.csv", True); ti = uscita("fittizio_4", "time.csv", True)
    dt = np.diff(np.r_[0, t]); dti = np.diff(np.r_[0, ti])
    V = -(np.sum(dti * Li.sum(1))) - np.cumsum(dt * L.sum(1))
    def graf(fig):
        ax = fig.subplots(2, 1, sharex=True)
        ax[0].plot(t / 3600, -np.minimum(L, 0).sum(1) * -1, color="#eb6834", lw=1.2, label="uscita verso la piana (−)")
        ax[0].plot(t / 3600, np.maximum(L, 0).sum(1), color="#1baf7a", lw=1.2, label="rientro (+)")
        ax[0].axhline(0, color="#9aa0a6", lw=0.7); ax[0].set_ylabel("portata [m³/s]", fontsize=8)
        ax[0].legend(fontsize=7, loc="upper right", ncol=2); ax[0].grid(alpha=0.25); ax[0].tick_params(labelsize=7)
        ax[1].plot(t / 3600, V / 1000, color="#2a78d6", lw=1.4, label="volume nel serbatoio")
        ax[1].axhline(150, color="#b4371f", lw=1, ls="--", label="soglia 150 000 m³")
        ax[1].set_ylabel("volume [migliaia di m³]", fontsize=8); ax[1].set_xlabel("tempo [h]", fontsize=8)
        ax[1].legend(fontsize=7, loc="upper right"); ax[1].grid(alpha=0.25); ax[1].tick_params(labelsize=7)
        fig.set_size_inches(7.6, 5.0)
    D.grafico(graf, "serbatoio", f"`fittizio_4` (due onde) con i valori tarati, ricostruito da `Q_lat.csv`: il "
              f"serbatoio arriva a {V.max()/1000:.0f} mila m³ e a fine evento ne contiene {V[-1]/1000:.0f} mila. "
              "Nella seconda onda la piana è già in parte piena e il rientro arriva prima.")
    D.h1("6. Taratura e risultati")
    D.p("I parametri sono stati tarati sui quattro scenari fittizi confrontando con il 2D (115 combinazioni, circa 480 "
        "simulazioni, tutte stabili e con il bilancio chiuso). Codici, risultati e grafici sono in "
        "`calibrazione_scambio/` con il suo `leggimi.md`. Medie sui quattro scenari:")
    D.tabella([["Scambio", "NSE in S6000", "Errore del colmo lungo il tratto", "Errore del livello"],
               ["sfioro (tabelle q2)", "0,937", "–", "–"],
               ["tabella, serbatoio con T = 0", "0,926", "4,18 m³/s", "5,9 cm"],
               ["tabella, solo T = 4 h", "0,986", "1,86 m³/s", "3,6 cm"],
               ["soglia 150 000 m³, T = 1 h", "0,989", "1,11 m³/s", "3,3 cm"],
               ["soglia, T = 0,5 h, τ = 1 h", "0,991", "1,02 m³/s", "3,1 cm"],
               ["+ rilascio lento 32 h", "0,994", "1,05 m³/s", "3,2 cm"],
               ["**+ serbatoio attivo**", "**0,998**", "**1,08 m³/s**", "**2,7 cm**"]], [0.36, 0.18, 0.28, 0.18])
    D.elenco(["**senza soglia** tutta l'acqua uscita rientra subito: il colmo in uscita resta quello in ingresso. Nel 2D "
              "il volume che resta fuori è simile nei quattro scenari (154-215 mila m³): da qui la soglia;",
              "**il ritardo** porta il rientro al momento giusto: in S4200 la portata supera 67 m³/s 44-50 min prima del 2D "
              "senza ritardo, 1-7 min dopo con τ = 1 h;",
              "**il rilascio lento** svuota la parte sotto la soglia tra due onde: senza, nella seconda onda del "
              "`fittizio_4` il rientro arrivava 112 min prima del 2D;",
              "**il serbatoio attivo** fa rientrare l'acqua oltre la soglia anche nel ramo calante: l'anticipo della "
              "discesa in S6000 scende da 32 a 14 min in media."])
    D.attenzione("un solo serbatoio per tutto il tratto (l'acqua uscita a 0,5 km può rientrare a 4 km), cinque "
                 "parametri tarati e verificati sugli stessi quattro scenari, nessuna piena di verifica indipendente.")
    D.collegamenti([(16, "16_scambio_tabella.pdf", "Lo scambio dalla tabella ΔQ"), (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa"),
                    (2, "02_lettura_settings.pdf", "Le righe 29-33 di settings.txt")])


# =============================================================================================== 18
def p18(D):
    Q = uscita("fittizio_1", "Q.csv", True); E = uscita("fittizio_1", "eta.csv", True); t = uscita("fittizio_1", "time.csv", True)
    D.h1("A cosa serve")
    D.p("Lo stato iniziale (approfondimento 07) è solo una stima: in un tratto non prismatico il moto uniforme sezione "
        "per sezione non è l'equilibrio, e a valle degli affluenti nemmeno le portate cumulate lo sono. Il modello "
        "quindi gira con la portata di monte **ferma al primo valore** (riempimento, `IC = 1`) finché il tratto si "
        "assesta. Il riempimento finisce quando **quattro condizioni** valgono insieme.")
    D.h1("Le quattro condizioni")
    D.eq([r"(1)\quad t\geq T_{ass}=3600\ \mathrm{s}",
          r"(2)\quad |Q^{n+1}_{i}-Q^n_{i}|<0.01\ \mathrm{m^3/s}\quad \mathrm{in\ tutte\ le\ celle}",
          r"(3)\quad |\eta_i(t)-\eta_i(t-T_{f})|<1\ \mathrm{mm}\quad \mathrm{in\ tutte\ le\ sezioni},\qquad T_f=1800\ \mathrm{s}",
          r"(4)\quad \left|Q_{monte}+\sum Q_{aff}+\sum Q_{lat}-Q_{valle}\right|<10^{-3}\left(Q_{monte}+\sum Q_{aff}+\sum Q_{lat}\right)"], numerata=False)
    D.elenco(["(1) una **durata minima**: il solo criterio sulla portata scattava troppo presto (nel canale rettangolare "
              "a 2,76 h con il livello ancora 6 mm lontano dall'equilibrio);",
              "(2) la **portata è ferma** da un passo all'altro;",
              "(3) il **livello è fermo**: ogni 1800 s si confronta con quello di inizio finestra; serve perché il "
              "livello si assesta più lentamente della portata;",
              "(4) il **bilancio**: tutto ciò che entra (monte, affluenti, scambio con la piana) esce a valle entro lo 0,1%."])
    D.estratto(SEMI, "            IC_counter = 0", "                consec_count = 0", titolo=SEMI)
    D.p("Quando le condizioni sono soddisfatte si salva il volume in alveo per il bilancio (`B_iniziale`), si passa a "
        "`IC = 0` e tempo e contatori ripartono da zero: inizia la simulazione con l'idrogramma. `t_tot` invece non "
        "riparte, perché il ritardo del serbatoio della piana deve vedere anche l'acqua uscita nel riempimento.")
    dq = np.abs(np.diff(Q, axis=0)).max(1)
    def graf(fig):
        ax = fig.subplots(1, 2)
        ax[0].semilogy(t[1:] / 3600, dq, color="#2a78d6", lw=1.2); ax[0].axhline(0.01, color="#b4371f", ls="--", lw=1)
        ax[0].set_title("max |Q(n+1) − Q(n)| [m³/s]", fontsize=9); ax[0].set_xlabel("tempo [h]", fontsize=8)
        k = int(round(1800 / (t[1] - t[0]))); w = np.arange(k, len(t), k)
        de = [np.abs(E[j] - E[j - k]).max() * 1000 for j in w]
        ax[1].semilogy(t[w] / 3600, de, "o-", color="#1baf7a", lw=1.2, ms=4); ax[1].axhline(1, color="#b4371f", ls="--", lw=1)
        ax[1].set_title("max variazione del livello in 30 min [mm]", fontsize=9); ax[1].set_xlabel("tempo [h]", fontsize=8)
        for a in ax: a.grid(alpha=0.25, which="both"); a.tick_params(labelsize=7); a.axvline(1, color="#9aa0a6", lw=0.8, ls=":")
    D.grafico(graf, "riempimento", f"Riempimento del `fittizio_1`: {len(t)} passi da 5 s, fine a {t[-1]/3600:.1f} h. "
              "Linee rosse: soglie dei criteri; linea punteggiata: durata minima di 1 h.")
    D.flowchart(dict(dimensioni=(14, 11.2), w=5.4, h=0.95, nodi={
        "a": (6, 10.6, "fine del passo, IC = 1", "inizio", 4.0, 0.75),
        "b": (6, 9.4, "conta le celle con |q − qn| < 0,01", "proc"),
        "c": (6, 8.0, "è passata una finestra\ndi T_fermo = 1800 s ?", "dec", 4.8, 1.2),
        "d": (11.6, 8.0, "fermo = tutte le |Δη| < 1 mm\nnuova finestra", "proc", 4.2, 1.0),
        "e": (6, 6.2, "t ≥ 3600 s  e  tutte le celle ferme\ne  fermo  e  |ingressi − uscita| < 0,1% ?", "dec", 8.2, 1.5),
        "f": (6, 4.3, "B_iniziale = Σ A·Δx\nIC = 0,  time = 0,  nCS = 1", "proc", 5.4, 1.1),
        "g": (6, 3.0, "inizia la simulazione", "fine", 3.8, 0.75),
        "h": (12.2, 6.2, "passo\nsuccessivo", "fine", 2.2, 0.9)},
        archi=[("a", "b"), ("b", "c"), ("c", "d", "sì"), ("c", "e", "no"), ("d.s", "e.n", None, [(11.6, 7.15), (6, 7.15)]),
               ("e", "f", "sì"), ("e", "h", "no"), ("f", "g")]),
        didascalia="Flowchart dei criteri di fine del riempimento.", larghezza=0.9)
    D.h1("Perché non la rampa")
    D.p("Prima il riempimento faceva crescere gli affluenti con una rampa di 6 h, inserendoli uno alla volta. Con le "
        "portate cumulate lo stato finale è identico ma il livello è entro 1 mm dall'equilibrio in 1,7 h invece di "
        "7,4 h nel Test_Brenta e in 4,7 h invece di 10,4 h nel canale rettangolare. Nei fittizi a 50 m il "
        "riempimento finisce a 2,5 h.")
    D.attenzione("la tolleranza assoluta di 0,01 m³/s è poco adatta a portate molto basse; T_{ass}, T_{f} e la tolleranza "
                 "di 1 mm sono provati solo sul Test_Brenta (6 km) e sul canale rettangolare (30 km).")
    D.collegamenti([(7, "07_stato_iniziale.pdf", "Lo stato iniziale"), (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa")])


# =============================================================================================== 19
def p19(D):
    testo = open(os.path.join(PROGETTO, USCITE, "fittizio_1", "bilancio.txt"), encoding="utf-8").read()
    D.h1("A cosa serve")
    D.p("Un modello idraulico deve conservare il volume d'acqua: quello che entra, meno quello che esce, è quello che "
        "resta. `BILANCIO_DI_MASSA` lo verifica a fine simulazione e scrive `bilancio.txt` nella cartella degli output. "
        "È il controllo più semplice per accorgersi di un errore nello schema, nelle tabelle o nello scambio con la piana.")
    D.h1("1. L'equazione del bilancio")
    D.p("Il volume in alveo è la somma delle aree delle tabelle per la lunghezza delle celle. Tra l'inizio della "
        "simulazione (fine del riempimento) e la fine deve valere:")
    D.eq([r"V_{alveo}=\sum_iA_i(\eta_i)\,\Delta x_i",
          r"\varepsilon=\left(V_{alveo}^{fine}-V_{alveo}^{inizio}\right)-\left(V_{monte}-V_{valle}+V_{aff}+V_{lat}\right)"])
    D.p("dove V_{monte}, V_{valle}, V_{aff} e V_{lat} sono le somme passo per passo di Δt·Q (approfondimento 13). Lo schema "
        "scrive la continuità in forma di bilancio di cella (approfondimento 11), quindi ε è zero a meno degli "
        "arrotondamenti e della tolleranza di Newton. Oltre lo 0,1% del volume da monte compare un avviso.")
    D.h1("2. Il controllo del serbatoio")
    D.p("Con lo scambio `tabella` il serbatoio della piana deve contenere esattamente l'acqua uscita e non ancora "
        "rientrata dall'inizio, **riempimento compreso**:")
    D.eq(r"V_{piana}^{fine}=-\sum_{\mathrm{tutti\ i\ passi}}\Delta t\sum_iQ_{lat,i}\qquad (\mathrm{tolleranza}\ 10^{-6}\ \mathrm{relativa})")
    D.p("Il primo membro è aggiornato dal serbatoio (approfondimento 17), il secondo dagli scambi entrati davvero nella "
        "continuità: se sono uguali, nessuna acqua è stata creata o persa nello scambio.")
    D.h1("3. Il codice")
    D.estratto(WRI, "function BILANCIO_DI_MASSA", "    return err", titolo=WRI)
    D.estratto(SEMI, "    # --- CONTROLLO DEL BILANCIO DI MASSA ---", n=8, titolo=SEMI)
    D.h1("4. Un esempio reale")
    D.codice(testo, titolo=USCITE + "/fittizio_1/bilancio.txt")
    import re as _re
    v = {k: float(x) for k, x in _re.findall(r"^(.+?)\s{2,}(-?\d+\.\d)", testo, _re.M)}
    ini = v["volume in alveo all'inizio"]; fin = v["volume in alveo alla fine"]
    mon = v["entrato da monte"]; val = v["uscito a valle"]; lat = v["scambio con la piana (+ entra, − esce)"]
    D.p(f"Variazione in alveo: {fin:.1f} − {ini:.1f} = {fin-ini:.1f} m³. Entrato − uscito: {mon:.1f} − {val:.1f} "
        f"{'+' if lat >= 0 else '−'} {abs(lat):.1f} = {mon-val+lat:.1f} m³. La differenza è sotto il decimo di m³ su 2,7 milioni. Il serbatoio "
        "contiene a fine evento circa 1 000 m³ in più dello scambio della sola simulazione: è l'acqua uscita durante "
        "il riempimento, compresa nel controllo del serbatoio.")
    D.flowchart(dict(dimensioni=(12, 8.8), w=6.4, h=0.95, nodi={
        "a": (6, 8.2, "fine del ciclo temporale", "inizio", 4.2, 0.75),
        "b": (6, 6.9, "il riempimento è finito ?", "dec", 4.2, 1.0),
        "c": (11.0, 6.9, "messaggio:\nbilancio non calcolato", "fine", 2.4, 1.0),
        "d": (6, 5.5, "ε = ΔV_alveo − (monte − valle + aff + lat)", "proc"),
        "e": (6, 4.2, "|ε| > 0,1% del volume da monte ?  → avviso", "proc"),
        "f": (6, 2.9, "scambio = tabella: V_piana = −Σ Δt·ΣQ_lat ?  → avviso", "proc", 7.2),
        "g": (6, 1.6, "scrive e stampa bilancio.txt", "io", 4.8)},
        archi=[("a", "b"), ("b", "c", "no"), ("b", "d", "sì"), ("d", "e"), ("e", "f"), ("f", "g")]),
        didascalia="Flowchart del controllo del bilancio.", larghezza=0.85)
    D.nota("`strumenti/grafico_volumi.jl` (grafico 15) ricostruisce gli stessi volumi passo per passo dagli output e li "
           "disegna nel tempo: entrato = uscito a valle + rimasto nella piana + aumento in alveo.")
    D.collegamenti([(11, "11_sistema_dei_livelli.pdf", "Perché il volume si conserva"), (17, "17_serbatoio_della_piana.pdf", "Il serbatoio della piana"),
                    (20, "20_file_di_uscita.pdf", "I file di uscita")])


# =============================================================================================== 20
def p20(D):
    d = os.path.join(PROGETTO, USCITE, "fittizio_1")
    dim = {f: os.path.getsize(os.path.join(d, f)) for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))}
    D.h1("A cosa serve")
    D.p("I risultati sono scritti in file di testo semplici, uno per grandezza, nella cartella della riga 22 di "
        "`settings.txt`. `INIT_OUTPUT` crea la cartella e svuota i file di una simulazione precedente; le funzioni "
        "`save_*` aggiungono una riga a ogni passo.")
    D.estratto(WRI, "function INIT_OUTPUT", n=9, titolo=WRI)
    D.h1("I file")
    righe = [["File", "Contenuto", "Valori per riga", "Dimensione (`fittizio_1`)"]]
    for f, c, n in [("Q.csv", "portate alle interfacce [m³/s]; l'ultima è l'uscita", "imax + 1"),
                    ("eta.csv", "quota del pelo libero nelle sezioni [m s.l.m.]", "imax"),
                    ("h.csv", "tirante = η − talweg [m]", "imax"),
                    ("time.csv", "tempo dall'inizio della fase [s]", "1"),
                    ("Q_lat.csv", "scambio con la piana [m³/s], + entra, − esce", "imax"),
                    ("IC_*.csv", "gli stessi file per il riempimento", "come sopra"),
                    ("bilancio.txt", "controllo del bilancio di massa (approfondimento 19)", "–")]:
        dd = dim.get(f); righe.append([f"`{f}`", c, n, f"{dd/1e6:.1f} MB" if dd else ("~1/3" if f.startswith("IC") else "1 kB")])
    D.tabella(righe, [0.17, 0.5, 0.15, 0.18])
    D.p("I valori sono separati da spazi, con 6 decimali e larghezza fissa 12; una riga per passo temporale. Il tempo è "
        "in un file separato, così ogni riga di `Q.csv` corrisponde alla stessa riga di `time.csv`. Con 120 sezioni e "
        "3 600 passi ogni file di sezione occupa una decina di MB.")
    D.estratto(WRI, "function save_Q!", n=8, titolo=WRI)
    D.h2("Come leggerli")
    D.codice('# Julia\nusing DelimitedFiles\nQ = readdlm("output/.../Q.csv")        # matrice passi × (imax+1)\nt = vec(readdlm("output/.../time.csv"))\n\n'
             '# Python\nimport numpy as np\nQ = np.loadtxt("output/.../Q.csv")\nQ_sezioni = 0.5*(Q[:, :-1] + Q[:, 1:])   # portata nelle sezioni', titolo="lettura")
    D.p("La portata nelle sezioni, quando serve, si ottiene come media delle due interfacce vicine (così fanno gli "
        "strumenti di confronto).")
    D.h1("Chi li usa")
    D.tabella([["Strumento", "File letti", "Grafici"],
               ["`grafici_simulazione.jl`", "IC_*, Q, eta, h, time", "1-6"],
               ["`confronto_reference.jl`", "Q, eta, Q_lat, time", "7-8"],
               ["`animazioni_2D.jl`", "eta, Q, time", "9"],
               ["`grafici_spinup.jl`", "IC_*", "10-11"],
               ["`confronto_valle.jl`", "Q, h, time", "12-13"],
               ["`errore_tirante_sezioni.jl`", "Q, h, time", "14"],
               ["`grafico_volumi.jl`", "Q, eta, IC_eta, Q_lat, time", "15"],
               ["`confronto_simulazioni_2D.jl`", "Q, h, Q_lat, time", "confronti tra simulazioni"]], [0.38, 0.38, 0.24])
    D.flowchart(dict(dimensioni=(13, 9.0), w=5.4, h=0.95, nodi={
        "a": (6.5, 8.4, "INIT_OUTPUT(output_dir)", "funz", 4.8),
        "b": (6.5, 7.1, "mkpath + 10 file vuoti", "io", 4.6),
        "c": (6.5, 5.8, "a ogni passo", "proc", 3.2),
        "d": (6.5, 4.4, "IC = 1 ?", "dec", 2.6, 1.0),
        "e": (3.2, 2.6, "save_* in IC_Q, IC_eta, IC_h,\nIC_time, IC_Q_lat", "io", 4.6, 1.1),
        "f": (9.8, 2.6, "save_* in Q, eta, h,\ntime, Q_lat", "io", 4.0, 1.1)},
        archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d.w", "e.n", "sì", [(3.2, 4.4)]), ("d.e", "f.n", "no", [(9.8, 4.4)])]),
        didascalia="Creazione dei file e scrittura a ogni passo.", larghezza=0.8)
    D.collegamenti([(13, "13_portate_salvataggio_tabelle.pdf", "Il salvataggio nel ciclo"), (19, "19_bilancio_di_massa.pdf", "bilancio.txt")])


# =============================================================================================== 21
def p21(D):
    T = np.loadtxt(os.path.join(PROGETTO, TABELLE), delimiter=",", skiprows=1)
    D.h1("A cosa serve")
    D.p("Le tabelle del modello (approfondimento 04) non si scrivono a mano: le costruisce `calcolo_tabelle_daBM` dai "
        "risultati del modello 2D BASEMENT a **scalini di portata**. Nel 2D si impone una sequenza di portate costanti "
        "(1, 2, 4 … 100 m³/s), ciascuna mantenuta fino al moto permanente; per ogni sezione (nodestring) BASEMENT "
        "scrive livello medio, area, larghezza e portata. Ogni scalino diventa una riga della tabella.")
    D.p("`julia --project=. calcolo_tabelle_daBM/src/main.jl` dalla cartella del progetto. Le costanti in testa scelgono "
        "gli scalini (`InputFolder`) e dove scrivere: le tabelle vanno sia nell'output dello script sia negli input "
        "del modello (`ModelFolder`).")
    D.h1("1. La pendenza del pelo libero")
    D.p("Per ogni scalino k i livelli delle sezioni formano un profilo. La pendenza in una sezione è calcolata alle "
        "differenze centrate con le vicine bagnate, unilaterali agli estremi:")
    D.eq(r"i_{WSE,j}=\frac{\eta_{j-1}-\eta_{j+1}}{s_{j+1}-s_{j-1}}\qquad(\mathrm{positiva\ se\ il\ livello\ scende\ verso\ valle})")
    D.p("Una pendenza nulla o negativa (livelli quasi uguali a portata molto bassa) è sostituita con quella del livello "
        "bagnato successivo della stessa sezione (`fix_slope!`), altrimenti il K_{s} non sarebbe definito.")
    D.estratto(CT, "function free_surface_slope", n=13, titolo=CT)
    D.h1("2. Il Ks equivalente e jden")
    D.p("Si suppone che in ogni sezione e per ogni scalino valga localmente la formula di moto uniforme con la pendenza "
        "del pelo libero. Il K_{s} che la rende vera è il **K_{s} equivalente** della sezione a quel livello:")
    D.eq([r"Q=K_s\,\Omega\,R_h^{2/3}\sqrt{i_{WSE}}\quad\Rightarrow\quad K_s=\sqrt{\frac{Q^2}{\left(\Omega\,R_h^{2/3}\right)^2 i_{WSE}}},\qquad R_h=\frac{\Omega}{W}",
          r"j_{den}=\frac{K_s^2\,R_h^{4/3}}{g}"])
    D.p("La portata è quella **locale** della sezione, cioè quella che resta nel canale principale dopo eventuali "
        "esondazioni. In questo modo la scabrezza del 1D contiene, livello per livello, tutto ciò che nel 2D rallenta "
        "la corrente: vegetazione, forma irregolare, perdite localizzate.")
    D.estratto(CT, "function Ks_equivalent", n=4, titolo=CT)
    def graf(fig):
        ax = fig.add_subplot(111)
        for k, col in ((20, "#2a78d6"), (60, "#eb6834"), (100, "#1baf7a"), (120, "#1f2328")):
            s = T[T[:, 0] == k]; Rh = s[:, 10]; Ks = np.sqrt(9.81 * s[:, 6] / np.maximum(Rh, 1e-9) ** (4 / 3))
            ax.plot(Ks, s[:, 1] - s[0, 1], "-o", ms=2, lw=1.2, color=col, label=f"sezione {k} (S{50*k:04d})")
        ax.set_xlim(0, 70); ax.set_xlabel("Ks equivalente [m^(1/3)/s]", fontsize=8); ax.set_ylabel("livello − prima quota [m]", fontsize=8)
        ax.grid(alpha=0.25); ax.tick_params(labelsize=7); ax.legend(fontsize=7)
    D.grafico(graf, "ks", "K_{s} equivalente ricavato dalle tabelle di quattro sezioni (K_{s} = √(g·jden/R_{h}^{4/3})). La sezione "
              "120 è l'uscita, con la pendenza media del fondo: il suo K_{s} è circa 39,5.")
    D.h1("3. L'ultima sezione")
    D.p("Nell'ultima sezione K_{s} e j_{den} sono ricavati con la pendenza che il modello usa per la condizione di uscita "
        "(la `local_slope` dell'ultima sezione, pari alla pendenza media del fondo, 0,00213) invece che con la pendenza "
        "del pelo libero, che nel 2D vicino all'uscita è molto ripida (fino a 0,0075). Così la formula di uscita "
        "Q = A√(g j_{den} S) restituisce, a ogni livello, proprio la portata dello scalino 2D: la condizione di valle "
        "del 1D è la scala di deflusso del 2D (approfondimento 09).")
    D.h1("4. Lo scambio con la piana: ΔQ")
    D.eq(r"\Delta Q_{k,j}=Q_{k,j}-Q_{k,j-1}\qquad (0\ \mathrm{nella\ prima\ sezione\ e\ se\ una\ delle\ due\ \grave{e}\ asciutta})")
    D.estratto(CT, "function discharge_difference", n=7, titolo=CT)
    D.h1("5. Righe asciutte e livelli non crescenti")
    D.elenco(["i primi scalini possono essere asciutti in alcune sezioni: la prima riga è ricostruita con tirante 1 cm, "
              "portata 0,1 m³/s e K_{s} = 20 (l'area viene da Manning), le altre fino al primo livello bagnato interpolate "
              "(`rebuild_dry!`). Il modello non scende mai a queste portate, ma i valori devono essere ragionevoli;",
              "il modello vuole quote crescenti: una riga con livello non superiore ai precedenti viene scartata "
              "(in S0050 gli ultimi scalini, dove il livello scende);",
              "le colonne delle pendenze (`_s`) sono calcolate con `slope!`: differenza in avanti, l'ultima uguale alla "
              "penultima. È ciò che rende coerente l'estrapolazione di `READ_VALUES`."])
    D.estratto(CT, "function slope!", n=4, titolo=CT)
    D.h1("6. Il flusso dello script")
    D.flowchart(dict(dimensioni=(14, 17.6), w=7.4, h=0.95, nodi={
        "a": (6.5, 17.0, "main(; dosave = true)", "inizio", 4.0, 0.75),
        "b": (6.5, 15.8, "check_input_files: mancano file ? → errore", "funz"),
        "c": (6.5, 14.5, "centerpoints.csv → sezioni da monte a valle\nload_section per ognuna (Ω, W, Q, WSE)", "io", 7.6, 1.1),
        "d": (6.5, 13.1, "free_surface_slope + fix_slope! (tutte le sezioni)", "funz"),
        "e": (6.5, 11.9, "discharge_difference → ΔQ (dati 2D originali)", "funz"),
        "f": (6.5, 10.7, "β = 1  oppure  beta_2D (mesh e risultati 2D)", "funz"),
        "g": (6.5, 9.5, "talweg e local_slope; ultima = pendenza media del fondo", "proc"),
        "h": (6.5, 8.3, "sezione j (da monte a valle)", "proc", 5.0),
        "i": (6.5, 7.1, "rebuild_dry! · scarta livelli non crescenti", "funz"),
        "l": (6.5, 5.9, "ultima sezione: i_WSE = |local_slope|", "proc"),
        "m": (6.5, 4.7, "section_table: Ω, β, W, Rh, Ks, jden, ΔQ e pendenze", "funz"),
        "n": (6.5, 3.5, "check_section: avvisi fuori dai limiti", "funz"),
        "q": (6.5, 2.2, "altre sezioni ?", "dec", 3.4, 1.0),
        "o": (6.5, 0.8, "tabella_sezioni.csv + localslope_xcoord.csv (argini da\nbankpoints.csv) in output e in ModelFolder", "io", 8.0, 1.1)},
        archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("e", "f"), ("f", "g"), ("g", "h"), ("h", "i"), ("i", "l"),
               ("l", "m"), ("m", "n"), ("n", "q"), ("q.e", "h.e", "sì", [(12.6, 2.2), (12.6, 8.3)]), ("q", "o", "no")]),
        didascalia="Flowchart di `calcolo_tabelle_daBM/src/main.jl`.")
    D.h1("7. β dal modello 2D (facoltativo)")
    D.p("Con `β_dal_2D = true` il coefficiente di quantità di moto è calcolato dalle velocità del 2D sugli elementi a "
        "monte di ogni nodestring, come nel metodo di Engelund:")
    D.eq(r"\beta=\frac{\Omega}{Q^2}\int u^2\,Y\,dy\ \approx\ \frac{\Omega\,\sum_e\left(q_{n,e}^2/Y_e\right)\Delta L_e}{\left(\sum_e q_{n,e}\,\Delta L_e\right)^2},\qquad \Omega=\sum_eY_e\,\Delta L_e")
    D.p("con Y_{e} il tirante dell'elemento, q_{n,e} la portata unitaria normale alla sezione e ΔL_{e} la lunghezza del lato. "
        "Vale 1 con velocità uniforme e in mediana 1,14 nel Test_Brenta, ma l'effetto sul modello è trascurabile: le "
        "tabelle in uso hanno β = 1.")
    D.attenzione("il K_{s} equivalente trascura la variazione dell'altezza cinetica u²/2g tra sezioni vicine (fino a 8,5 cm "
                 "a 50 m): dove il profilo 2D non è monotono il 1D sbaglia di 5-10 cm, altrove 1-2 cm in media.")
    D.collegamenti([(4, "04_tabelle_delle_sezioni.pdf", "Le tabelle nel modello"), (9, "09_parte_esplicita_rusanov.pdf", "La condizione di valle"),
                    (16, "16_scambio_tabella.pdf", "Lo scambio dalla tabella ΔQ")])


PDF = [
    (15, "15_sfioro_argini.pdf", "Lo sfioro sugli argini", "Scambio `sfioro`: stramazzo a parete grossa, `FREE_OVERFLOW`", "§11.1", p15),
    (16, "16_scambio_tabella.pdf", "Lo scambio dalla tabella ΔQ", "Scambio `tabella`: `READ_EXCHANGE` e perché da solo è instabile", "§11.2", p16),
    (17, "17_serbatoio_della_piana.pdf", "Il serbatoio della piana", "Soglia, tempo di svuotamento, ritardo, rilascio lento e serbatoio attivo: `VOLUME_AT`", "§11.2", p17),
    (18, "18_fine_riempimento.pdf", "La fine del riempimento", "I quattro criteri di arresto della fase iniziale", "§12.1", p18),
    (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa", "`BILANCIO_DI_MASSA` e il controllo del serbatoio della piana", "§12.2", p19),
    (20, "20_file_di_uscita.pdf", "I file di uscita", "`INIT_OUTPUT`, `save_*` e come leggere i risultati", "§13", p20),
    (21, "21_costruzione_tabelle.pdf", "La costruzione delle tabelle", "`calcolo_tabelle_daBM`: Ks equivalente, jden, ΔQ, β", "§14", p21),
]
