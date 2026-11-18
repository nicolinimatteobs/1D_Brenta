# Approfondimenti 15-21: sfioro, scambio dalla tabella, serbatoio della piana, fine del riempimento, bilancio di massa,
# file di uscita, costruzione delle tabelle.
import numpy as np
from stile import *

SEMI = "Source/Semi_Implicit_Brenta_2D_1D.jl"
SUB = "Source/subroutines.jl"
WRI = "Source/writing.jl"
CT = "calcolo_tabelle_daBM/src/main.jl"


# =============================================================================================== 15
def p15(D):
    D.h1("A cosa serve")
    D.p("Quando il livello supera la cresta di un argine, l'acqua esce dal canale verso la piana. Con lo scambio "
        "`sfioro` (riga 10 di `settings.txt`) il modello lo rappresenta come uno **stramazzo a parete grossa** lungo "
        "ciascun argine della cella. L'acqua uscita non rientra: è persa per il fiume.")
    D.h1("1. La formula dello sfioro")
    D.p("Su una soglia larga, con carico Y_{0} sopra la cresta, la corrente passa per lo stato critico: tirante "
        "critico 2/3 del carico, velocità √(g·2Y_{0}/3). La portata per unità di lunghezza di argine è quindi")
    D.eq([r"q_l=\frac{2}{3}Y_0\,\sqrt{g\,\frac{2}{3}Y_0}=\left(\frac{2}{3}\right)^{3/2}\sqrt{g}\;Y_0^{3/2}\approx 0.385\,\sqrt{2g\,Y_0^3}",
          r"Y_0=\max\left(\eta-z_{argine},\ 0\right),\qquad C_d=0.385"])
    D.p("È la formula 16 di Siviglia, Stocchino e Colombini (ASCE 2009). Il coefficiente 0,385 è "
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
    D.p("Le quote degli argini sono le colonne 5-6 della geometria (approfondimento 03); `calcolo_tabelle_daBM` le "
        "prende da `bankpoints.csv`, un punto del DTM per argine e per sezione. Una geometria a 4 colonne ha argini a "
        "`Inf` e non sfiora mai.")
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
    D.elenco(["**solo uscita**: l'acqua che nella realtà rientra più a valle nel 1D è persa, quindi a valle delle "
              "esondazioni la portata è sottostimata;",
              "**un punto per argine**: la quota della cresta è un solo punto applicato a tutta la lunghezza della cella;",
              "**sfioro libero**: non sa che cosa c'è dall'altra parte. Dove la piana è già allagata, un punto basso "
              "dell'argine può essere un punto di rientro, non di uscita;",
              "per questo c'è anche lo scambio dalla tabella con il serbatoio della piana (approfondimenti 16-17)."])
    D.collegamenti([(16, "16_scambio_tabella.pdf", "Lo scambio dalla tabella ΔQ"), (17, "17_serbatoio_della_piana.pdf", "Il serbatoio della piana"),
                    (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa")])


# =============================================================================================== 16
def p16(D):
    D.h1("A cosa serve")
    D.p("Lo scambio `tabella` sostituisce la formula dello sfioro con un'informazione presa direttamente dal modello 2D. "
        "Negli scalini 2D a moto permanente la portata cambia lungo il tratto **solo** per lo scambio con la piana: se "
        "la sezione k ha meno portata della k−1 l'acqua è uscita nel tratto, se ne ha di più è rientrata.")
    D.h1("1. La colonna ΔQ")
    D.eq(r"\Delta Q_k(\eta)=Q^{2D}_k-Q^{2D}_{k-1}\qquad \mathrm{allo\ stesso\ scalino\ di\ portata}")
    D.p("ΔQ < 0 dove l'acqua esce, ΔQ > 0 dove rientra, 0 nella prima sezione e dove una delle due è asciutta. È "
        "calcolata da `calcolo_tabelle_daBM` (approfondimento 21) e scritta nelle colonne 13-14 delle tabelle, con la "
        "sua pendenza rispetto al livello, così si interpola come le altre grandezze. Non servono né la formula dello "
        "sfioro né le quote degli argini. Serve invece che le tabelle abbiano 14 colonne.")
    def graf(fig):
        ax = fig.subplots(1, 2); e = np.linspace(100, 103.5, 36)
        uscita = -4.0 * np.maximum(e - 102.5, 0.0) ** 1.5
        rientro = np.where(e > 102.0, 5.0 * np.minimum(e - 102.0, 1.2), 0.0)
        for a, dq, nome, col in ((ax[0], uscita, "sezione di uscita", "#eb6834"), (ax[1], rientro, "sezione di rientro", "#1baf7a")):
            a.plot(dq, e, "-o", ms=2.5, color=col, lw=1.3)
            a.axvline(0, color="#9aa0a6", lw=0.8); a.set_title(nome, fontsize=9)
            a.set_xlabel("ΔQ [m³/s]", fontsize=8); a.grid(alpha=0.25); a.tick_params(labelsize=7)
        ax[0].set_ylabel("livello η [m s.l.m.]", fontsize=8)
    D.grafico(graf, "deltaQ", "Forma tipica di ΔQ(η) (esempio schematico): sotto la quota delle sponde è zero; sopra, "
              "nelle sezioni di uscita diventa negativo, in quelle di rientro positivo e crescente con il livello.")
    D.h1("2. READ_EXCHANGE")
    D.p("Legge ΔQ al livello della sezione con la stessa ricerca di riga di `READ_VALUES` (approfondimento 04), ma con "
        "una differenza importante agli estremi: **oltre l'ultima quota usa l'ultimo valore, senza estrapolare**. ΔQ può "
        "cambiare molto tra le ultime righe e una retta prolungata darebbe scambi irrealistici.")
    D.eq(r"\Delta Q(\eta)=\Delta Q_1\ (\eta\leq\eta_1),\qquad \Delta Q_r+\Delta Q'_r\,(\eta-\eta_r)\ (\eta_r\leq\eta<\eta_{r+1}),\qquad \Delta Q_N\ (\eta\geq\eta_N)")
    D.estratto(SUB, "@inline function READ_EXCHANGE", n=11, titolo=SUB)
    D.h1("3. Da solo è instabile")
    D.p("Usando ΔQ così com'è, con Q_{lat,i} = ΔQ_{i}(η_{i}), il calcolo diverge: a valle dei rientri la portata "
        "supera quella in ingresso e il passo temporale va a zero, anche senza NaN. Il motivo è un **anello di "
        "retroazione positivo** nelle sezioni di rientro.")
    D.p("Si consideri una sola cella di rientro, con larghezza W e lunghezza Δx, e una piccola perturbazione del "
        "livello η' attorno all'equilibrio. Dalla continuità:")
    D.eq([r"W\,\Delta x\,\frac{d\eta'}{dt}=\frac{d\Delta Q}{d\eta}\,\eta'+\ldots\quad\Rightarrow\quad \eta'(t)\propto e^{\lambda t},\qquad \lambda=\frac{1}{W\,\Delta x}\,\frac{d\Delta Q}{d\eta}"])
    D.p("Dove l'acqua rientra ΔQ **cresce** con il livello (dΔQ/dη > 0, figura 1 a destra): λ > 0 e la perturbazione cresce "
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
        "quattro parametri più un interruttore (righe 11-15 di `settings.txt`).")
    D.tabella([["Riga", "Parametro", "Simb.", "Significato fisico", "Unità"],
               ["11", "soglia", "V_{s}", "la piana deve riempirsi prima di restituire acqua (avvallamenti, ristagni)", "m³"],
               ["12", "tempo di svuotamento", "T", "con quale velocità la piana piena restituisce (0 = subito)", "h"],
               ["13", "ritardo", "τ", "tempo per attraversare la piana dall'uscita al rientro", "h"],
               ["14", "rilascio lento", "T_{l}", "la parte sotto la soglia torna lentamente (0 = resta nella piana)", "h"],
               ["15", "serbatoio attivo", "–", "l'acqua oltre la soglia rientra anche quando la tabella non lo chiede", "si/no"]],
              [0.07, 0.19, 0.09, 0.53, 0.12])
    D.p("I valori dipendono dal tratto e dalla sua piana: vanno tarati per ogni caso, confrontando il 1D con un "
        "modello 2D o con misure (portate a valle, istanti dei rientri, volume rimasto nella piana).")
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
    D.h1("6. Che cosa fa ciascun parametro")
    D.elenco(["**senza soglia** (V_{s} = 0) e con T = 0 tutta l'acqua uscita rientra subito: le uscite e i rientri della "
              "tabella si compensano e la piana non lamina. La soglia trattiene il volume che nella realtà resta fuori;",
              "**T** distribuisce nel tempo la restituzione dell'acqua oltre la soglia: con T grande la piana restituisce "
              "lentamente e i colmi a valle calano;",
              "**il ritardo τ** sposta il rientro più tardi, dando il tempo all'acqua di attraversare la piana;",
              "**il rilascio lento** svuota la parte sotto la soglia, per esempio tra due onde di piena: senza, quella "
              "parte resta nella piana per sempre e la seconda onda trova la piana già piena;",
              "**il serbatoio attivo** fa rientrare l'acqua oltre la soglia anche quando la tabella non chiede rientri "
              "(nel ramo calante, quando i livelli scendono). Con `no` quell'acqua resta ferma nel serbatoio."])
    D.attenzione("un solo serbatoio per tutto il tratto: l'acqua uscita in un punto può rientrare in qualunque sezione "
                 "di rientro, anche lontana. Per tratti lunghi, con piane separate, lo schema andrebbe esteso a più serbatoi.")
    D.collegamenti([(16, "16_scambio_tabella.pdf", "Lo scambio dalla tabella ΔQ"), (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa"),
                    (2, "02_lettura_settings.pdf", "Le righe 11-15 di settings.txt")])


# =============================================================================================== 18
def p18(D):
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
    D.elenco(["(1) una **durata minima**: il solo criterio sulla portata può scattare troppo presto, con il livello "
              "ancora qualche millimetro lontano dall'equilibrio;",
              "(2) la **portata è ferma** da un passo all'altro;",
              "(3) il **livello è fermo**: ogni 1800 s si confronta con quello di inizio finestra; serve perché il "
              "livello si assesta più lentamente della portata;",
              "(4) il **bilancio**: tutto ciò che entra (monte, affluenti, scambio con la piana) esce a valle entro lo 0,1%."])
    D.estratto(SEMI, "            IC_counter = 0", "                consec_count = 0", titolo=SEMI)
    D.p("Quando le condizioni sono soddisfatte si salva il volume in alveo per il bilancio (`B_iniziale`), si passa a "
        "`IC = 0` e tempo e contatori ripartono da zero: inizia la simulazione con l'idrogramma. `t_tot` invece non "
        "riparte, perché il ritardo del serbatoio della piana deve vedere anche l'acqua uscita nel riempimento.")
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
    D.p("Nel riempimento il passo è al massimo di 5 s e il CFL è 0,6 (approfondimento 08): conta arrivare a un "
        "equilibrio pulito, non la velocità. La durata dipende dal tratto e dalla portata iniziale: più il tratto è "
        "lungo e lontano dal moto uniforme, più dura.")
    D.attenzione("la tolleranza assoluta di 0,01 m³/s è poco adatta a portate molto basse; T_{ass}, T_{f} e la tolleranza "
                 "di 1 mm sono costanti nel codice: per tratti molto lunghi o portate molto basse vanno verificati.")
    D.collegamenti([(7, "07_stato_iniziale.pdf", "Lo stato iniziale"), (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa")])


# =============================================================================================== 19
def p19(D):
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
    D.h1("4. Che cosa contiene bilancio.txt")
    D.elenco(["volume in alveo all'inizio e alla fine della simulazione e la sua variazione;",
              "volumi entrato da monte, uscito a valle, entrato dagli affluenti e scambiato con la piana (+ entra, − esce);",
              "errore ε in m³ e in percentuale del volume da monte, con l'avviso se supera lo 0,1%;",
              "con lo scambio `tabella`, il volume nel serbatoio della piana e il suo controllo. Il serbatoio contiene "
              "anche l'acqua uscita durante il riempimento, che non compare nei volumi della sola simulazione."])
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
    D.h1("A cosa serve")
    D.p("I risultati sono scritti in file di testo semplici, uno per grandezza, nella cartella della riga 5 di "
        "`settings.txt`. `INIT_OUTPUT` crea la cartella e svuota i file di una simulazione precedente; le funzioni "
        "`save_*` aggiungono una riga a ogni passo.")
    D.estratto(WRI, "function INIT_OUTPUT", n=9, titolo=WRI)
    D.h1("I file")
    righe = [["File", "Contenuto", "Valori per riga"]]
    for f, c, n in [("Q.csv", "portate alle interfacce [m³/s]; l'ultima è l'uscita", "imax + 1"),
                    ("eta.csv", "quota del pelo libero nelle sezioni [m s.l.m.]", "imax"),
                    ("h.csv", "tirante = η − talweg [m]", "imax"),
                    ("time.csv", "tempo dall'inizio della fase [s]", "1"),
                    ("Q_lat.csv", "scambio con la piana [m³/s], + entra, − esce", "imax"),
                    ("IC_*.csv", "gli stessi file per il riempimento", "come sopra"),
                    ("bilancio.txt", "controllo del bilancio di massa (approfondimento 19)", "–")]:
        righe.append([f"`{f}`", c, n])
    D.tabella(righe, [0.2, 0.6, 0.2])
    D.p("I valori sono separati da spazi, con 6 decimali e larghezza fissa 12; una riga per passo temporale. Il tempo è "
        "in un file separato, così ogni riga di `Q.csv` corrisponde alla stessa riga di `time.csv`. Ogni valore occupa "
        "12 caratteri: un file di sezione pesa circa 12·imax byte per passo (con 100 sezioni e 4 000 passi, circa 5 MB).")
    D.estratto(WRI, "function save_Q!", n=8, titolo=WRI)
    D.h2("Come leggerli")
    D.codice('# Julia\nusing DelimitedFiles\nQ = readdlm("output/.../Q.csv")        # matrice passi × (imax+1)\nt = vec(readdlm("output/.../time.csv"))\n\n'
             '# Python\nimport numpy as np\nQ = np.loadtxt("output/.../Q.csv")\nQ_sezioni = 0.5*(Q[:, :-1] + Q[:, 1:])   # portata nelle sezioni', titolo="lettura")
    D.p("La portata nelle sezioni, quando serve, si ottiene come media delle due interfacce vicine (così fanno gli "
        "strumenti di confronto).")
    D.h1("Chi li usa")
    D.p("Gli script di `strumenti/` leggono questi file per i grafici: per esempio `grafici_simulazione.jl` (condizioni "
        "iniziali, idrogrammi, tiranti, profili), `grafici_spinup.jl` (il riempimento, dai file `IC_*`) e "
        "`grafico_volumi.jl` (il bilancio dei volumi nel tempo). Gli strumenti di confronto con un modello 2D "
        "richiedono in più i file di riferimento del 2D. Tabelle e geometria vengono dalle righe 1-2 di "
        "`settings.txt`: devono essere quelle della simulazione che si analizza.")
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
    D.h1("A cosa serve")
    D.p("Le tabelle del modello (approfondimento 04) non si scrivono a mano: le costruisce `calcolo_tabelle_daBM` dai "
        "risultati del modello 2D BASEMENT a **scalini di portata**. Nel 2D si impone una sequenza di portate costanti "
        "che copre l'intervallo di interesse, dalla magra alla piena, ciascuna mantenuta fino al moto permanente; per ogni sezione (nodestring) BASEMENT "
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
    D.p("Il K_{s} equivalente non è un parametro fisico costante: cambia con il livello e da una sezione all'altra, "
        "perché raccoglie tutto quello che il 1D non descrive. Valori molto lontani da quelli tipici dei corsi d'acqua "
        "(per esempio sotto 10 o sopra 60 m^{1/3}/s) segnalano di solito un profilo 2D irregolare in quella sezione.")
    D.h1("3. L'ultima sezione")
    D.p("Nell'ultima sezione K_{s} e j_{den} sono ricavati con la pendenza che il modello usa per la condizione di uscita "
        "(la `local_slope` dell'ultima sezione, pari alla pendenza media del fondo) invece che con la pendenza del pelo "
        "libero, che nel 2D vicino all'uscita del dominio è spesso molto ripida (richiamo della condizione al contorno "
        "del 2D). Così la formula di uscita "
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
              "(succede negli scalini più alti di qualche sezione, dove il livello medio sulla nodestring può scendere);",
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
        "Vale 1 con velocità uniforme e cresce quanto più la velocità è disuniforme nella sezione. Con `β_dal_2D = false` "
        "le tabelle hanno β = 1.")
    D.attenzione("il K_{s} equivalente trascura la variazione dell'altezza cinetica u²/2g tra sezioni vicine. Dove il profilo "
                 "2D non è monotono (il pelo libero risale tra due sezioni) il 1D, che con il solo attrito non può far "
                 "risalire il livello, non riproduce esattamente i livelli del 2D.")
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
