# Approfondimenti 08-14: interfacce e passo, trasporto e contorno, attrito, sistema dei livelli, Newton e Thomas,
# portate e salvataggio, controllo dei NaN.
import numpy as np
from stile import *

SEMI = "Source/Semi_Implicit_Brenta_2D_1D.jl"
SUB = "Source/subroutines.jl"
F1 = USCITE + "/fittizio_1"


def carica_tabelle():
    T = np.loadtxt(os.path.join(PROGETTO, TABELLE), delimiter=",", skiprows=1)
    return [T[T[:, 0] == k] for k in range(1, int(T[:, 0].max()) + 1)]


def leggi(tab, eta, col):
    """Come READ_VALUES: interpolazione con la colonna delle pendenze, estrapolazione oltre l'ultima quota."""
    ys = tab[:, 1]; z = max(eta, ys[0])
    r = min(max(np.searchsorted(ys, z, side="right") - 1, 0), len(ys) - 1)
    d = z - ys[r]
    if r == len(ys) - 1 and d > 0:
        return tab[r, col] + (tab[r, col] - tab[r - 1, col]) / (ys[r] - ys[r - 1]) * d
    return tab[r, col] + tab[r, col + 1] * d


def stato_al_colmo():
    t = np.loadtxt(os.path.join(PROGETTO, F1, "time.csv"))
    Q = np.loadtxt(os.path.join(PROGETTO, F1, "Q.csv")); E = np.loadtxt(os.path.join(PROGETTO, F1, "eta.csv"))
    n = int(Q[:, 0].argmax())
    return t, Q, E, n


# =============================================================================================== 08
def p08(D):
    D.h1("A cosa serve")
    D.p("Ogni passo temporale comincia preparando le grandezze che servono **sulle interfacce** (dove stanno le "
        "portate) a partire da quelle **nelle sezioni** (lette dalle tabelle al livello corrente), e scegliendo "
        "il passo temporale Δt in modo che il trasporto esplicito resti stabile.")
    D.h1("1. Dalle sezioni alle interfacce")
    D.p("Area, β, j_{den} e R_{h} sono note nelle sezioni: `A_b`, `β_b`, `den_γ_b`, `rh_b` (aggiornate alla fine del "
        "passo precedente, approfondimento 13). Sulle interfacce interne si usano due regole diverse:")
    D.eq([r"A_{i-1/2}=A_{i-1}\ \mathrm{se}\ Q_{i-1/2}\geq 0,\qquad A_{i-1/2}=A_{i}\ \mathrm{se}\ Q_{i-1/2}<0\qquad(\mathrm{upwind})",
          r"j_{den,i-1/2}\ \mathrm{come}\ A\ (\mathrm{upwind}),\qquad \beta_{i-1/2}=\frac{\beta_{i-1}+\beta_i}{2},\qquad R_{h,i-1/2}=\frac{R_{h,i-1}+R_{h,i}}{2}"])
    D.p("Agli estremi l'interfaccia prende i valori della sezione vicina. Poi la **velocità** e la **massima "
        "velocità di trasporto** su ogni interfaccia:")
    D.eq(r"u_{i-1/2}=\frac{Q_{i-1/2}\,A_{i-1/2}}{A_{i-1/2}^2+\epsilon}\approx\frac{Q_{i-1/2}}{A_{i-1/2}},\qquad s_{i-1/2}=\left|\beta_{i-1/2}\,u_{i-1/2}\right|,\qquad \epsilon=10^{-14}")
    D.p("La forma Q·A/(A² + ε) evita la divisione per zero in una cella asciutta (A = 0 dà u = 0).")
    D.h2("Perché upwind per area e attrito")
    D.p("Con la media aritmetica (schema centrato) compare un **dente di sega** nel profilo dei livelli nei tratti "
        "poco profondi: livelli alternativamente alti e bassi da una sezione all'altra. Non è un errore dei dati ma la "
        "soluzione stazionaria dello schema centrato quando la cella è lunga rispetto alla distanza in cui il livello "
        "torna al moto uniforme, circa Δx > 0,6·y/S (con y il tirante e S la pendenza). Prendere area e attrito dalla "
        "cella a monte toglie l'oscillazione, al prezzo di un errore del primo ordine nei tratti profondi.")
    D.attenzione("su un canale rettangolare con un velo di 10 cm su 30 km, l'upwind su area **e** attrito insieme è "
                 "instabile e va in NaN; solo su uno dei due, o con la media, è stabile. Sul Test_Brenta (tiranti "
                 "0,4-3 m) funziona. La versione con la media è nel codice, commentata.")
    D.estratto(SEMI, "        @inbounds @simd for i in 2:imax", "            smax[i] = abs(β_int[i] * uf[i])", titolo=SEMI)
    D.h1("2. Il passo temporale")
    D.p("Nello schema semi-implicito il gradiente del pelo libero e l'attrito sono impliciti: la celerità delle onde "
        "di gravità √(g·h) non limita il passo. Resta il limite del **trasporto esplicito** (approfondimento 09): "
        "in un passo l'informazione non deve attraversare più di una cella. È la condizione di Courant-Friedrichs-"
        "Lewy (CFL) sulla sola velocità della corrente:")
    D.eq(r"\Delta t=\min\left(\Delta t_{lim},\ \mathrm{CFL}\,\frac{\Delta x_{min}}{\max_{i}\,s_{i-1/2}+\epsilon}\right),\qquad \Delta x_{min}=\min\left(\min_i\Delta x_{i+1/2},\ \min_i\Delta x_i\right)")
    D.tabella([["", "Riempimento (IC = 1)", "Simulazione (IC = 0)"],
               ["CFL", "0,6", "0,95"],
               ["Δt_{lim}", "`dt_riempimento` = 5 s", "`dtmax` (riga 25), 60 s"],
               ["Scopo", "assestamento senza oscillazioni", "efficienza"]], [0.2, 0.4, 0.4])
    D.p("Infine il passo è accorciato perché la simulazione finisca esattamente a t_{end}, la fine dell'idrogramma più breve.")
    D.estratto(SEMI, "        amax = maximum(smax)", n=12, titolo=SEMI)
    t = np.loadtxt(os.path.join(PROGETTO, F1, "time.csv")); dt = np.diff(np.r_[0, t])
    ti = np.loadtxt(os.path.join(PROGETTO, F1, "IC_time.csv"))
    def graf(fig):
        ax = fig.add_subplot(111)
        ax.plot(t[1:] / 3600, dt[1:], color="#2a78d6", lw=1)
        ax2 = ax.twinx()
        Q = np.loadtxt(os.path.join(PROGETTO, F1, "Q.csv"))
        ax2.plot(t / 3600, Q[:, 0], color="#9aa0a6", lw=1, ls="--")
        ax2.set_ylabel("portata in ingresso [m³/s] (tratteggio)", fontsize=8, color="#5f6368")
        ax.set_xlabel("tempo [h]", fontsize=8); ax.set_ylabel("passo temporale Δt [s]", fontsize=8, color="#2a78d6")
        ax.grid(alpha=0.25); ax.tick_params(labelsize=7); ax2.tick_params(labelsize=7)
    D.grafico(graf, "passo", f"Passo temporale nella simulazione del `fittizio_1` (`svuotamento_32h`): {len(t)} passi, "
              f"Δt da {dt[1:].min():.1f} a {dt.max():.1f} s, in media {dt.mean():.1f} s. Con la piena le velocità "
              "crescono e il passo si accorcia. Nel riempimento il passo è sempre 5 s "
              f"({len(ti)} passi, {ti[-1]/3600:.1f} h).")
    s_med = 0.95 * 50 / dt.mean()
    D.p(f"Con Δx = 50 m e CFL = 0,95, un passo medio di {dt.mean():.1f} s corrisponde a una velocità massima di circa "
        f"{s_med:.1f} m/s. Uno schema esplicito sulle onde di gravità, con h ≈ 3 m, avrebbe avuto bisogno di "
        f"Δt ≈ 0,95·50/(2 + 5,4) ≈ 6 s: il semi-implicito permette passi circa quattro volte più lunghi.")
    D.flowchart(dict(dimensioni=(13, 10.6), w=6.6, h=0.95, nodi={
        "a": (6, 10.0, "Inizio del passo", "inizio", 3.2, 0.75),
        "b": (6, 8.8, "qn = q (portata del passo precedente)", "proc"),
        "c": (6, 7.5, "Interfacce di estremità = valori delle sezioni vicine", "proc"),
        "d": (6, 6.1, "Interfacce interne: A e jden upwind (segno di q)\nβ e Rh come media", "proc", 6.6, 1.1),
        "e": (6, 4.7, "u = q·A/(A² + ε),   s = |β·u|", "proc"),
        "f": (6, 3.4, "Δt = min(Δt_lim, CFL·Δx_min/(max s + ε))", "proc"),
        "g": (6, 2.1, "t + Δt ≥ t_end ?  →  Δt = t_end − t", "dec", 6.4, 1.1),
        "h": (6, 0.8, "verso Q* (approfondimento 09)", "fine", 4.6, 0.7)},
        archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("e", "f"), ("f", "g"), ("g", "h")]),
        didascalia="Flowchart della prima parte del passo temporale.", larghezza=0.85)
    D.collegamenti([(9, "09_parte_esplicita_rusanov.pdf", "La parte esplicita"), (10, "10_attrito_implicito.pdf", "L'attrito implicito"),
                    (13, "13_portate_salvataggio_tabelle.pdf", "Aggiornamento delle tabelle")])


# =============================================================================================== 09
def p09(D):
    D.h1("A cosa serve")
    D.p("Questa parte calcola la portata «provvisoria» Q* su ogni interfaccia: la portata del passo precedente "
        "modificata dal solo **trasporto di quantità di moto**, ∂(βQ²/A)/∂x, calcolato in modo esplicito. Sulle due "
        "interfacce di estremità Q* è invece la **condizione al contorno**: l'idrogramma a monte, il moto uniforme a "
        "valle. Gravità e attrito si aggiungono dopo, in modo implicito.")
    D.h1("1. Il trasporto con il flusso di Rusanov")
    D.p("L'equazione della quantità di moto si integra su una «cella duale» centrata sull'interfaccia i+1/2, "
        "lunga Δx_{i+1/2} e delimitata dalle sezioni i e i+1. Il flusso di quantità di moto F = βuQ si valuta nelle "
        "sezioni, a partire dalle due interfacce vicine:")
    D.eq([r"Q^{*}_{i+1/2}=Q^{n}_{i+1/2}-\frac{\Delta t}{\Delta x_{i+1/2}}\left(F_{i+1}-F_{i}\right)",
          r"F_{i}=\frac{(\beta u Q)_{i-1/2}+(\beta u Q)_{i+1/2}}{2}-\frac{\lambda_i}{2}\left(Q_{i+1/2}-Q_{i-1/2}\right)",
          r"\lambda_i=2\,\max\left(|\beta u|_{i-1/2},\ |\beta u|_{i+1/2}\right)"])
    D.p("Il primo termine di F è la media dei flussi (centrata, del secondo ordine ma instabile da sola); il secondo è "
        "una **diffusione numerica** proporzionale alla massima velocità di propagazione λ, che stabilizza lo schema. "
        "È il flusso di **Rusanov** (o Lax-Friedrichs locale). Per il flusso βQ²/A la velocità di propagazione è "
        "la derivata rispetto a Q:")
    D.eq(r"\lambda=\frac{\partial}{\partial Q}\left(\beta\,\frac{Q^2}{A}\right)=2\,\beta\,u")
    D.p("da cui il fattore 2: nel codice compare come `0.5*2*max(smax[i], smax[i+1])`. La diffusione rende lo schema "
        "del primo ordine e smorza leggermente i fronti, ma garantisce che Q* non sviluppi oscillazioni, purché sia "
        "rispettata la condizione di Courant dell'approfondimento 08.")
    D.h1("2. Condizione al contorno di monte")
    D.p("Sull'interfaccia di ingresso la portata è imposta: nel riempimento il primo valore dell'idrogramma, nella "
        "simulazione il valore interpolato a t + Δt (approfondimento 05):")
    D.eq(r"Q^{*}_{1/2}=Q_{monte}(t^n+\Delta t)")
    D.h1("3. Condizione al contorno di valle")
    D.p("All'uscita non si conosce né il livello né la portata. Si suppone **moto uniforme** nell'ultima sezione e si "
        "calcola la portata con l'area e l'attrito della tabella dell'ultima sezione al passo precedente:")
    D.eq(r"Q^{*}_{imax+1/2}=A_{imax}\,\sqrt{g\,j_{den,imax}\,|S_{imax}|}")
    D.p("S è la `local_slope` dell'ultima sezione, che `calcolo_tabelle_daBM` pone uguale alla **pendenza media del "
        "fondo** (0,00213). Lo stesso script ricava il j_{den} dell'ultima sezione proprio con quella pendenza "
        "(approfondimento 21): per costruzione la formula restituisce, livello per livello, la portata degli scalini "
        "2D in S6000. In altre parole la condizione di valle è la **scala di deflusso** del modello 2D, senza "
        "tarature a mano.")
    D.nota("prima si usava Manning con il K_{s} = 33 della riga 23 di `settings.txt`: a 50 m dava in S6000 tiranti "
           "10-13 cm più alti del 2D, perché nel 2D il pelo libero vicino all'uscita è molto più ripido (richiamo "
           "dell'uscita del dominio). La riga 23 è rimasta solo per non spostare le altre.")
    D.h1("4. Il codice")
    D.estratto(SEMI, "        @inbounds @simd for i in 1:imax+1", "        CHECK_NAN(\"calcolo di qstar\"", titolo=SEMI)
    D.p("Nello stesso ciclo si calcolano γ e l'area «ridotta» A/(1 + Δtγ) di ogni interfaccia, che servono per "
        "l'attrito implicito (approfondimento 10).")
    D.flowchart(dict(dimensioni=(14, 9.0), w=4.2, h=1.0, nodi={
        "a": (7, 8.4, "Per ogni interfaccia i = 1 … imax+1", "inizio", 5.6, 0.75),
        "b": (7, 7.0, "i = 1 ?", "dec", 2.4, 1.0),
        "c": (2.4, 5.4, "IC = 1: Q* = Q_monte(0)\nIC = 0: Q* = LAGRANGE(t + Δt)", "funz", 4.4, 1.1),
        "d": (7, 5.4, "i = imax+1 ?", "dec", 3.0, 1.0),
        "e": (11.6, 3.9, "Q* = UNIFORM_DISCHARGE\n(A, jden, S) ultima sezione", "funz", 4.2, 1.1),
        "f": (7, 3.9, "F a monte e a valle (Rusanov)\nQ* = Q − Δt/Δx·(F⁺ − F⁻)", "proc", 4.6, 1.1),
        "g": (7, 2.2, "γ = |u|/jden,  A_g = A/(1 + Δtγ)", "proc", 5.0),
        "h": (7, 0.8, "CHECK_NAN(qstar, γ, A_int_g)", "funz", 5.0, 0.8)},
        archi=[("a", "b"), ("b.w", "c.n", "sì", [(2.4, 7.0)]), ("b", "d", "no"), ("d.e", "e.n", "sì", [(11.6, 5.4)]),
               ("d", "f", "no"), ("f", "g"), ("c.s", "g.w", None, [(2.4, 2.2)]), ("e.s", "g.e", None, [(11.6, 2.2)]), ("g", "h")]),
        didascalia="Flowchart del calcolo di Q* con le condizioni al contorno.", larghezza=0.9)
    D.attenzione("la portata in uscita usa area e attrito del passo precedente: è esplicita. Il livello nell'ultima "
                 "sezione si adegua nel passo successivo. Con passi di 20-30 s e variazioni lente non dà oscillazioni.")
    D.collegamenti([(5, "05_idrogrammi.pdf", "Gli idrogrammi"), (8, "08_interfacce_e_passo_temporale.pdf", "Interfacce e passo"),
                    (10, "10_attrito_implicito.pdf", "L'attrito implicito"), (21, "21_costruzione_tabelle.pdf", "Il jden dell'ultima sezione")])


# =============================================================================================== 10
def p10(D):
    tabs = carica_tabelle(); t, Q, E, n = stato_al_colmo()
    D.h1("A cosa serve")
    D.p("L'attrito sul fondo e sulle sponde frena la corrente. Nei tratti poco profondi è molto forte: trattarlo in "
        "modo esplicito imporrebbe passi temporali piccolissimi. Il modello lo tratta in modo **implicito**, "
        "sfruttando il fatto che con j_{den} è lineare nella portata.")
    D.h1("1. Il problema dell'attrito esplicito")
    D.p("Considerando solo l'attrito, l'equazione della quantità di moto su un'interfaccia è")
    D.eq(r"\frac{dQ}{dt}=-\gamma\,Q,\qquad \gamma=\frac{|u|}{j_{den}}")
    D.p("che ha soluzione esatta Q(t) = Q_{0}e^{−γt}: la portata decade senza cambiare segno. Con Eulero esplicito, "
        "Q^{n+1} = (1 − Δtγ)Q^{n}: per Δtγ > 1 la portata cambia segno a ogni passo, per Δtγ > 2 le oscillazioni "
        "crescono e il calcolo esplode. Con Eulero implicito invece")
    D.eq(r"Q^{n+1}=Q^n-\Delta t\,\gamma\,Q^{n+1}\quad\Rightarrow\quad Q^{n+1}=\frac{Q^n}{1+\Delta t\,\gamma}")
    D.p("il fattore 1/(1 + Δtγ) è sempre tra 0 e 1: la portata si riduce senza mai cambiare segno, **qualunque** sia "
        "il passo. γ è calcolato con la velocità del passo precedente (linearizzazione), così non servono iterazioni.")
    def graf(fig):
        ax = fig.add_subplot(111)
        x = np.linspace(0, 3, 300)
        ax.plot(x, np.exp(-x), color="#1f2328", lw=2, label="esatto  e^(−Δtγ)")
        ax.plot(x, 1 / (1 + x), color="#2a78d6", lw=1.6, label="implicito  1/(1+Δtγ)")
        ax.plot(x, 1 - x, color="#eb6834", lw=1.6, ls="--", label="esplicito  1−Δtγ")
        ax.axhline(0, color="#9aa0a6", lw=0.8); ax.axhline(-1, color="#9aa0a6", lw=0.6, ls=":")
        ax.set_ylim(-2.1, 1.1); ax.set_xlabel("Δt·γ", fontsize=8); ax.set_ylabel("Q(n+1)/Q(n) per il solo attrito", fontsize=8)
        ax.legend(fontsize=7.5, loc="lower left"); ax.grid(alpha=0.25); ax.tick_params(labelsize=7)
        ax.text(2.05, -1.25, "esplicito instabile\n(|fattore| > 1)", fontsize=7.5, color="#b4371f")
    D.grafico(graf, "attrito", "Fattore di riduzione della portata in un passo dovuto al solo attrito.")
    D.h1("2. Valori di γ nel Test_Brenta")
    g = []
    for i in range(len(tabs)):
        q = 0.5 * (Q[n, i] + Q[n, i + 1]); A = leggi(tabs[i], E[n, i], 2); jd = leggi(tabs[i], E[n, i], 6)
        g.append(abs(q / A) / jd)
    g = np.array(g)
    q0 = [0.5 * (Q[0, i] + Q[0, i + 1]) for i in range(len(tabs))]
    g0 = np.array([abs(q0[i] / leggi(tabs[i], E[0, i], 2)) / leggi(tabs[i], E[0, i], 6) for i in range(len(tabs))])
    dtm = np.diff(np.r_[0, t]).mean()
    D.p(f"Con gli output del `fittizio_1`, all'istante del colmo in ingresso γ va da {g.min():.4f} a {g.max():.4f} 1/s "
        f"(mediana {np.median(g):.4f}); con la portata di base iniziale da {g0.min():.4f} a {g0.max():.4f} 1/s. Con un "
        f"passo di circa {dtm:.0f} s, Δtγ arriva a {dtm*max(g.max(), g0.max()):.2f}: lo schema esplicito sarebbe al limite o "
        "oltre nelle sezioni più scabre o meno profonde, quello implicito non ha problemi.")
    def graf2(fig):
        ax = fig.add_subplot(111); x = np.arange(1, len(g) + 1) * 0.05
        ax.semilogy(x, g0, color="#9aa0a6", lw=1.2, label="portata di base (inizio)")
        ax.semilogy(x, g, color="#2a78d6", lw=1.4, label="colmo in ingresso")
        ax.set_xlabel("x [km]", fontsize=8); ax.set_ylabel("γ = |u|/jden [1/s]", fontsize=8)
        ax.legend(fontsize=7.5); ax.grid(alpha=0.25, which="both"); ax.tick_params(labelsize=7)
    D.grafico(graf2, "gamma", "γ lungo il tratto nel `fittizio_1`, ricostruito dagli output e dalle tabelle.")
    D.h1("3. Nello schema completo")
    D.p("Nella quantità di moto discreta l'attrito si somma al gradiente del pelo libero, anch'esso implicito:")
    D.eq([r"Q^{n+1}_{i+1/2}=Q^{*}_{i+1/2}-g\,A_{i+1/2}\frac{\Delta t}{\Delta x_{i+1/2}}\left(\eta^{n+1}_{i+1}-\eta^{n+1}_i\right)-\Delta t\,\gamma_{i+1/2}\,Q^{n+1}_{i+1/2}",
          r"\Rightarrow\quad Q^{n+1}_{i+1/2}=\frac{Q^{*}_{i+1/2}}{1+\Delta t\,\gamma_{i+1/2}}-\frac{g\,\tilde{A}_{i+1/2}\,\Delta t}{\Delta x_{i+1/2}}\left(\eta^{n+1}_{i+1}-\eta^{n+1}_i\right),\qquad \tilde{A}=\frac{A}{1+\Delta t\,\gamma}"])
    D.p("Compaiono quindi due grandezze, calcolate nel ciclo di Q* (approfondimento 09): il fattore 1/(1 + Δtγ) che "
        "riduce Q*, e l'area ridotta Ã (`A_int_g`), che pesa quanto la differenza di livello riesce a spingere "
        "l'acqua. In una cella molto scabra Ã è piccola: il livello deve salire di più per far passare la stessa portata.")
    D.estratto(SEMI, "            γ[i]       = abs(uf[i])/(den_γ_int[i] + 1e-14)", n=2, titolo=SEMI)
    D.h2("Da dove viene j_{den}")
    D.p("j_{den} = K_{s}²R_{h}^{4/3}/g è nella tabella; il K_{s} è quello **equivalente** che, nel modello 2D a moto "
        "permanente, dà la portata dello scalino con la pendenza del pelo libero misurata (approfondimento 21). Per "
        "questo il modello 1D in moto vario usa la scabrezza del 2D livello per livello, senza un K_{s} unico.")
    D.idea("il denominatore 1 + Δtγ è il cuore della stabilità: l'attrito toglie portata ma non può mai invertirla.")
    D.collegamenti([(9, "09_parte_esplicita_rusanov.pdf", "La parte esplicita"), (11, "11_sistema_dei_livelli.pdf", "Il sistema per i livelli"),
                    (21, "21_costruzione_tabelle.pdf", "Il Ks equivalente")])


# =============================================================================================== 11
def p11(D):
    D.h1("A cosa serve")
    D.p("Qui le due equazioni discrete si uniscono. La portata nuova dipende dai livelli nuovi (quantità di moto "
        "implicita); sostituendola nella continuità si ottiene **un'equazione per cella nei soli livelli** del passo "
        "n+1. Il codice prepara i coefficienti `aa`, `bb`, `cc` e il termine noto `rhs`.")
    D.h1("1. La sostituzione")
    D.p("Dall'approfondimento 10 la portata sull'interfaccia è")
    D.eq(r"Q^{n+1}_{i+1/2}=\frac{Q^{*}_{i+1/2}}{1+\Delta t\,\gamma_{i+1/2}}-\frac{g\,\tilde{A}_{i+1/2}\,\Delta t}{\Delta x_{i+1/2}}\left(\eta_{i+1}-\eta_i\right)")
    D.p("e la continuità della cella i (approfondimento 01) è")
    D.eq(r"A_i(\eta_i)=A^n_i-\frac{\Delta t}{\Delta x_i}\left(Q^{n+1}_{i+1/2}-Q^{n+1}_{i-1/2}\right)+\frac{\Delta t}{\Delta x_i}\left(Q_{aff,i}+Q_{lat,i}\right)")
    D.p("Sostituendo e portando a sinistra i termini con i livelli nuovi:")
    D.eq([r"A_i(\eta_i)+a_i\,\eta_{i-1}+b_i\,\eta_i+c_i\,\eta_{i+1}=r_i",
          r"a_i=-\frac{g\,\Delta t^2}{\Delta x_i\,\Delta x_{i-1/2}}\,\tilde{A}_{i-1/2},\qquad c_i=-\frac{g\,\Delta t^2}{\Delta x_i\,\Delta x_{i+1/2}}\,\tilde{A}_{i+1/2},\qquad b_i=-a_i-c_i",
          r"r_i=A^n_i-\frac{\Delta t}{\Delta x_i}\left(\frac{Q^{*}_{i+1/2}}{1+\Delta t\,\gamma_{i+1/2}}-\frac{Q^{*}_{i-1/2}}{1+\Delta t\,\gamma_{i-1/2}}\right)+\frac{\Delta t}{\Delta x_i}\left(Q_{aff,i}+Q_{lat,i}\right)"])
    D.p("La parte a_{i}η_{i−1} + b_{i}η_{i} + c_{i}η_{i+1} è il **laplaciano discreto** del livello, pesato dalle aree "
        "ridotte: è il termine di gravità, che tende a livellare il pelo libero tra celle vicine. r_{i} contiene "
        "tutto ciò che è noto: il volume del passo precedente, il trasporto esplicito, gli apporti laterali.")
    D.h2("Agli estremi")
    D.p("Nella prima cella non c'è una cella a monte (a_{1} = 0) e la portata in ingresso entra nel termine noto così "
        "com'è, senza il fattore d'attrito: è imposta. Nell'ultima c_{imax} = 0 e la portata in uscita (moto uniforme, "
        "approfondimento 09) entra anch'essa nel termine noto. I livelli agli estremi non sono imposti: si adattano.")
    D.estratto(SEMI, "        A_bn .= A_b", "        ####  INSERIMENTO AFFLUENTI", titolo=SEMI)
    D.p("Il termine noto si completa subito dopo con gli affluenti (approfondimento 06) e con lo scambio con la "
        "piana (approfondimenti 15-17), entrambi nella forma `rhs[i] += dt/dx_int[i]*Q`.")
    D.h1("2. La struttura del sistema")
    D.p("Ogni equazione coinvolge solo il livello della cella e delle due vicine: il sistema è **tridiagonale**. "
        "Scritto per tutte le celle, con A(η) il vettore delle aree:")
    D.eq(r"\mathbf{A}(\boldsymbol{\eta})+\mathbf{T}\,\boldsymbol{\eta}=\mathbf{r},\qquad \mathbf{T}=\mathrm{tridiag}(a_i,\ b_i,\ c_i)")
    def graf(fig):
        ax = fig.add_subplot(111); N = 8
        ax.set_xlim(-0.5, N - 0.5); ax.set_ylim(N - 0.5, -0.5); ax.set_aspect("equal")
        for i in range(N):
            for j in range(N):
                if abs(i - j) <= 1:
                    col = "#2a78d6" if i == j else "#eb6834"
                    ax.add_patch(plt.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9, color=col, alpha=0.85 if i == j else 0.55))
                    ax.text(j, i, "W+b" if i == j else ("a" if j < i else "c"), ha="center", va="center", fontsize=7, color="white")
                else:
                    ax.text(j, i, "0", ha="center", va="center", fontsize=7, color="#9aa0a6")
        ax.set_xticks(range(N)); ax.set_yticks(range(N))
        ax.set_xticklabels([f"η{j+1}" for j in range(N)], fontsize=7); ax.set_yticklabels([f"cella {i+1}" for i in range(N)], fontsize=7)
        ax.tick_params(length=0); ax.set_title("Matrice jacobiana per 8 celle", fontsize=9)
    D.grafico(graf, "matrice", "Struttura della matrice jacobiana del sistema (approfondimento 12): diagonale principale "
              "W_{i} + b_{i}, sotto a_{i}, sopra c_{i}, il resto è zero.", larghezza=0.6)
    D.h2("Proprietà che rendono il calcolo sicuro")
    D.elenco(["a_{i} ≤ 0 e c_{i} ≤ 0 perché aree, g e Δt sono positivi;",
              "b_{i} = −a_{i} − c_{i} = |a_{i}| + |c_{i}|: la somma di ogni riga di T è zero (un livello uniforme non genera flusso);",
              "nella jacobiana la diagonale è W_{i} + b_{i} > |a_{i}| + |c_{i}|, perché la larghezza W_{i} è positiva: la matrice è "
              "**a diagonale strettamente dominante**. Thomas è quindi stabile senza pivoting e il sistema ha una "
              "sola soluzione;",
              "più Δt è grande, più pesa il laplaciano rispetto a W: il passo lungo «accoppia» di più le celle, ma la "
              "dominanza resta."])
    D.h1("3. Conservazione del volume")
    D.p("Moltiplicando la continuità per Δx_{i} e sommando su tutte le celle, i flussi interni si eliminano a coppie "
        "(quello che esce dalla cella i entra nella i+1):")
    D.eq(r"\sum_i\Delta x_i\left(A^{n+1}_i-A^{n}_i\right)=\Delta t\left(Q_{1/2}-Q_{imax+1/2}\right)+\Delta t\sum_i\left(Q_{aff,i}+Q_{lat,i}\right)")
    D.p("Se Newton converge, la variazione del volume in alveo è esattamente uguale a quanto entra meno quanto esce: "
        "è ciò che controlla `BILANCIO_DI_MASSA` (approfondimento 19), con errori di 0,0000%.")
    D.collegamenti([(10, "10_attrito_implicito.pdf", "L'attrito implicito"), (12, "12_newton_thomas.pdf", "Newton + Thomas"),
                    (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa")])


# =============================================================================================== 12
def p12(D):
    tabs = carica_tabelle(); r = tabs[19]
    D.h1("A cosa serve")
    D.p("Il sistema dei livelli (approfondimento 11) è **non lineare**, perché l'area A_{i}(η_{i}) viene dalle tabelle. "
        "Si risolve con il metodo di Newton: a ogni iterazione si linearizza attorno ai livelli correnti e si risolve "
        "un sistema lineare tridiagonale con l'algoritmo di Thomas.")
    D.h1("1. Il metodo di Newton")
    D.p("Si cercano gli zeri del residuo f, una componente per cella:")
    D.eq(r"f_i(\boldsymbol{\eta})=A_i(\eta_i)+a_i\,\eta_{i-1}+b_i\,\eta_i+c_i\,\eta_{i+1}-r_i=0")
    D.p("La jacobiana ha la stessa struttura di T con in più, sulla diagonale, la derivata dell'area, che è la "
        "larghezza in superficie W = dA/dη (letta dalla tabella):")
    D.eq([r"\frac{\partial f_i}{\partial\eta_{i-1}}=a_i,\qquad \frac{\partial f_i}{\partial\eta_i}=W_i(\eta_i)+b_i,\qquad \frac{\partial f_i}{\partial\eta_{i+1}}=c_i",
          r"\mathbf{J}^{(k)}\,\boldsymbol{\delta}^{(k)}=\mathbf{f}^{(k)},\qquad \boldsymbol{\eta}^{(k+1)}=\boldsymbol{\eta}^{(k)}-\boldsymbol{\delta}^{(k)}"])
    D.p("Il punto di partenza è il livello del passo precedente. A(η) è lineare a tratti (approfondimento 04): "
        "Newton trova la soluzione esatta del tratto lineare in cui si trova, e se la soluzione cade in un altro "
        "tratto ci arriva in poche iterazioni. A(η) è anche **convessa** (la larghezza cresce con il livello), "
        "il che rende la convergenza monotona.")
    D.h2("Un esempio su una sola cella")
    D.p("Per vedere il comportamento: nella sezione 20 (S1000) si cerca il livello che dà un'area di 30 m², partendo "
        "dal fondo, con l'interpolazione della tabella e W come derivata.")
    At = 30.0; e = r[0, 1] + 0.05; righe = [["Iterazione", "η [m]", "A(η) − 30 [m²]", "W [m]"]]
    for k in range(7):
        A = leggi(r, e, 2); W = leggi(r, e, 8)
        righe.append([str(k), f"{e:.6f}", f"{A-At:+.3e}", f"{W:.3f}"])
        if abs(A - At) < 1e-12: break
        e = e - (A - At) / W
    D.tabella(righe, [0.2, 0.27, 0.3, 0.23], didascalia="Newton su una sola cella con la tabella vera: pochi passi "
              "per entrare nel tratto giusto, poi il residuo va a zero.")
    D.h2("Celle asciutte")
    D.p("Se nelle due interfacce di una cella il raggio idraulico è quasi nullo (R_{h,i−1/2} + R_{h,i+1/2} < 10^{−3} m), la "
        "cella è asciutta: la sua riga diventa δ_{i} = 0 (diagonale 1, residuo 0) e il livello non cambia. Evita "
        "divisioni per larghezze nulle.")
    D.h2("Il criterio di arresto")
    D.p("Newton si ferma quando la norma euclidea del residuo è sotto `tol` = 10^{−14} m², oppure dopo 100 "
        "iterazioni. La tolleranza è molto stretta, vicina all'errore di arrotondamento sulle aree (decine di m²): "
        "in pratica spesso si fanno tutte le 100 iterazioni anche se il residuo è già trascurabile dopo poche. "
        "È uno dei punti aperti in `CLAUDE.md`: una tolleranza relativa farebbe risparmiare tempo.")
    D.estratto(SEMI, "        # --- Ciclo Newton ---", "            eta .= eta .- deta", titolo=SEMI)
    D.h1("2. L'algoritmo di Thomas")
    D.p("Un sistema tridiagonale con N incognite si risolve con l'eliminazione di Gauss senza riempire la matrice: "
        "costo proporzionale a N invece di N³. Con a_{i} sotto la diagonale, d_{i} sulla diagonale, c_{i} sopra e f_{i} "
        "termine noto:")
    D.eq([r"c'_1=\frac{c_1}{d_1},\qquad f'_1=\frac{f_1}{d_1}",
          r"g_i=d_i-c'_{i-1}\,a_i,\qquad c'_i=\frac{c_i}{g_i},\qquad f'_i=\frac{f_i-a_i\,f'_{i-1}}{g_i}\qquad (i=2\ldots N)",
          r"\delta_N=f'_N,\qquad \delta_i=f'_i-c'_i\,\delta_{i+1}\qquad (i=N-1\ldots 1)"])
    D.p("La prima passata elimina la diagonale inferiore (da monte a valle), la seconda risale da valle a monte. Con la "
        "dominanza diagonale dell'approfondimento 11 ogni g_{i} è positivo e non si divide mai per zero.")
    D.estratto(SUB, "function THOMAS_ALGORITHM!", "    return (deta)", titolo=SUB)
    D.flowchart(dict(dimensioni=(14, 12.2), w=6.2, h=0.95, nodi={
        "a": (6, 11.6, "Newton: η^(0) = livelli del passo precedente", "inizio", 6.8, 0.75),
        "b": (6, 10.3, "Per ogni cella: READ_VALUES(η_i) → A_i, W_i", "funz"),
        "c": (6, 9.0, "f_i = A_i + a_i η_(i−1) + b_i η_i + c_i η_(i+1) − r_i\nd_i = W_i + b_i", "proc", 6.6, 1.1),
        "d": (6, 7.6, "cella asciutta ?  (Rh_int sinistro + destro < 1e-3)", "dec", 6.8, 1.1),
        "e": (12, 7.6, "d_i = 1\nf_i = 0", "proc", 2.4, 1.0),
        "f": (6, 6.1, "‖f‖ < 1e-14 ?", "dec", 3.4, 1.0),
        "g": (12, 6.1, "esci da Newton", "fine", 2.8, 0.75),
        "h": (6, 4.7, "THOMAS_ALGORITHM!(a, c, d, f) → δ", "funz"),
        "i": (6, 3.4, "CHECK_NAN(f, d, δ)", "funz", 4.0, 0.8),
        "l": (6, 2.2, "η = η − δ", "proc", 3.0, 0.8),
        "m": (6, 1.0, "iterazione < 100 ?", "dec", 3.8, 1.0)},
        archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d", "e", "sì"), ("d", "f", "no"), ("e.s", "f.n", None, [(12, 6.85), (6, 6.85)]),
               ("f", "g", "sì"), ("f", "h", "no"), ("h", "i"), ("i", "l"), ("l", "m"),
               ("m.w", "b.w", "sì", [(1.2, 1.0), (1.2, 10.3)]), ("m.e", "g.s", "no", [(12, 1.0)])]),
        didascalia="Flowchart del ciclo di Newton con Thomas.", larghezza=0.9)
    D.collegamenti([(4, "04_tabelle_delle_sezioni.pdf", "Le tabelle (A e W)"), (11, "11_sistema_dei_livelli.pdf", "Il sistema per i livelli"),
                    (13, "13_portate_salvataggio_tabelle.pdf", "Le nuove portate")])


# =============================================================================================== 13
def p13(D):
    D.h1("A cosa serve")
    D.p("Con i nuovi livelli si chiude il passo: si calcolano le nuove portate, si avanza il tempo, si aggiornano i "
        "contatori del bilancio, si scrivono i risultati e si rileggono le tabelle ai nuovi livelli per il passo "
        "successivo.")
    D.h1("1. Le nuove portate")
    D.p("Noti i livelli η^{n+1}, la quantità di moto implicita (approfondimento 10) dà direttamente la portata su ogni "
        "interfaccia interna:")
    D.eq(r"Q^{n+1}_{i-1/2}=\frac{Q^{*}_{i-1/2}-g\,A_{i-1/2}\,\frac{\Delta t}{\Delta x_{i-1/2}}\left(\eta^{n+1}_i-\eta^{n+1}_{i-1}\right)}{1+\Delta t\,\gamma_{i-1/2}},\qquad i=2\ldots imax")
    D.p("Sulle due interfacce di estremità resta il valore imposto: Q_{1/2} = Q*_{1/2} (idrogramma) e Q_{imax+1/2} = Q*_{imax+1/2} "
        "(moto uniforme). È la stessa formula usata per costruire il sistema, quindi le portate sono coerenti con i "
        "livelli: il volume si conserva.")
    D.estratto(SEMI, "        # --- Aggiornamento portate finali ---", "        CHECK_NAN(\"aggiornamento delle portate\"", titolo=SEMI)
    D.h1("2. Tempo, contatori e bilancio")
    D.p("Si avanzano il tempo della fase (`time`, che riparte da zero alla fine del riempimento), il tempo totale "
        "(`t_tot`, che serve al ritardo del serbatoio della piana) e il numero del passo. Si sommano i volumi del "
        "passo per il bilancio di massa: entrato da monte Δt·Q_{1/2}, uscito a valle Δt·Q_{imax+1/2}, affluenti e "
        "scambio laterale. Lo scambio è sommato anche nel riempimento (`B_lat_tutto`), per il controllo del serbatoio.")
    D.eq([r"V_{monte}\leftarrow V_{monte}+\Delta t\,Q_{1/2},\qquad V_{valle}\leftarrow V_{valle}+\Delta t\,Q_{imax+1/2}",
          r"V_{aff}\leftarrow V_{aff}+\Delta t\sum_iQ_{aff,i},\qquad V_{lat}\leftarrow V_{lat}+\Delta t\sum_iQ_{lat,i}"])
    D.estratto(SEMI, "        # --- Aggiornamento tempo e counter ---", "            B_monte += dt*q[1]", titolo=SEMI)
    D.h1("3. Salvataggio")
    D.p("A ogni passo si aggiunge una riga ai file di uscita: con il prefisso `IC_` durante il riempimento, senza "
        "durante la simulazione. Il tirante salvato è livello meno talweg. Il formato è nell'approfondimento 20.")
    D.estratto(SEMI, "            save_Q!(out.IC_Q, q)", "            save_Q!(out.Q_lat, Q_lat)", titolo=SEMI)
    D.h1("4. Le tabelle ai nuovi livelli")
    D.p("Infine si rileggono dalle tabelle area, β, j_{den}, larghezza e raggio idraulico di ogni sezione ai livelli "
        "η^{n+1}. Sono i valori che il passo successivo userà sulle interfacce (approfondimento 08) e come A^{n} "
        "nella continuità. Il volume in alveo per il bilancio è calcolato con le stesse aree.")
    D.estratto(SEMI, "        # --- Aggiornamento geometrie finali ---", n=9, titolo=SEMI)
    D.flowchart(dict(dimensioni=(13, 10.0), w=7.2, h=0.95, nodi={
        "a": (6.5, 9.4, "Livelli η(n+1) da Newton", "inizio", 4.6, 0.75),
        "b": (6.5, 8.2, "Q(n+1) sulle interfacce interne; estremi = Q*", "proc"),
        "c": (6.5, 7.0, "CHECK_NAN(q)", "funz", 3.4, 0.8),
        "d": (6.5, 5.8, "time += Δt,  t_tot += Δt,  nCS += 1", "proc"),
        "e": (6.5, 4.6, "volumi del bilancio (monte, valle, affluenti, scambio)", "proc"),
        "f": (6.5, 3.4, "save_*: IC_*.csv (IC = 1) oppure *.csv (IC = 0)", "io"),
        "g": (6.5, 2.2, "READ_VALUES(η(n+1)) → A, β, jden, W, Rh", "funz"),
        "h": (6.5, 1.0, "IC = 1: criteri di fine riempimento (approf. 18)", "fine", 6.6, 0.75)},
        archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("e", "f"), ("f", "g"), ("g", "h")]),
        didascalia="Flowchart della chiusura del passo temporale.", larghezza=0.85)
    D.attenzione("i file sono aperti e chiusi a ogni passo e a ogni riga (modo `\"a\"`, append). È semplice e robusto "
                 "(se il calcolo si ferma i dati scritti restano), ma con molte sezioni e molti passi pesa sui tempi "
                 "di calcolo.")
    D.collegamenti([(12, "12_newton_thomas.pdf", "Newton + Thomas"), (18, "18_fine_riempimento.pdf", "La fine del riempimento"),
                    (19, "19_bilancio_di_massa.pdf", "Il bilancio di massa"), (20, "20_file_di_uscita.pdf", "I file di uscita")])


# =============================================================================================== 14
def p14(D):
    D.h1("A cosa serve")
    D.p("Un NaN (Not a Number) nasce da operazioni come 0/0 o √(numero negativo) e poi si propaga: in un passo "
        "contamina le celle vicine e in pochi passi tutto il dominio. Senza controlli il programma continuerebbe a "
        "scrivere file pieni di NaN e sarebbe difficile capire **dove** e **quando** è nato. `CHECK_NAN` controlla "
        "i vettori dopo ogni fase del calcolo e, al primo NaN, ferma il programma con un messaggio preciso.")
    D.h1("Il codice")
    D.estratto(SUB, "function CHECK_NAN", n=9, titolo=SUB)
    D.p("Gli argomenti con nome (`vettori...`) permettono di passare quanti vettori si vuole: il nome dell'argomento "
        "compare nel messaggio. `findall(isnan, v)` trova le posizioni; il messaggio ne mostra al massimo 10.")
    D.h1("Dove si controlla")
    D.tabella([["Fase", "Vettori controllati", "Che cosa può generare un NaN"],
               ["calcolo di dt", "`dt`, `uf`, `smax`", "area nulla con portata non nulla, tabelle con valori non validi"],
               ["calcolo di qstar", "`qstar`, `γ`, `A_int_g`", "j_{den} nullo o negativo, idrogramma letto male"],
               ["termine noto e affluenti", "`rhs`, `aa`, `cc`", "affluenti o scambio laterale non validi"],
               ["Thomas, iterazione k di Newton", "`ff`, `di`, `deta`", "diagonale nulla, livello fuori dalle tabelle"],
               ["aggiornamento delle portate", "`q`", "propagazione di uno dei precedenti"]], [0.27, 0.27, 0.46])
    D.h2("Un esempio di messaggio")
    D.codice("ERROR: NaN in deta dopo Thomas, iterazione di Newton 3: passo 1532 (simulazione),\n"
             "t = 36842.1 s, posizioni [57, 58] (2 in totale)", titolo="esempio (non da una simulazione reale)")
    D.p("Il messaggio dà le quattro informazioni che `CLAUDE.md` chiede di riportare quando compare un NaN: la "
        "**fase**, il **passo**, il **tempo** e le **posizioni** (sezioni o interfacce). Le posizioni sono indici: "
        "per le sezioni i corrisponde alla sezione S(50·i), per le interfacce i alla faccia tra le sezioni i−1 e i.")
    D.flowchart(dict(dimensioni=(12, 7.8), w=5.6, h=0.95, nodi={
        "a": (5.5, 7.2, "CHECK_NAN(fase, nCS, time, IC; vettori...)", "inizio", 6.6, 0.75),
        "b": (5.5, 5.9, "per ogni vettore (nome, v)", "proc"),
        "c": (5.5, 4.5, "c'è un NaN in v ?", "dec", 3.8, 1.1),
        "d": (5.5, 2.9, "error(\"NaN in nome dopo fase:\npasso, fase, t, posizioni\")", "fine", 5.6, 1.1),
        "e": (10.2, 4.5, "vettore\nsuccessivo", "proc", 2.2, 1.0),
        "f": (10.2, 2.9, "nessun NaN:\nsi continua", "fine", 2.2, 0.9)},
        archi=[("a", "b"), ("b", "c"), ("c", "d", "sì"), ("c", "e", "no"), ("e.n", "b.e", None, [(10.2, 5.9)]), ("e", "f", "fine")]),
        didascalia="Flowchart di `CHECK_NAN`.", larghezza=0.8)
    D.h1("Che cosa fare quando compare")
    D.numerato(["leggere fase, passo, tempo e posizioni dal messaggio;",
                "guardare nelle uscite (`eta.csv`, `Q.csv`, o gli `IC_*` nel riempimento) le ultime righe scritte "
                "attorno alle posizioni indicate: spesso il livello è sceso sotto il fondo o è salito molto oltre "
                "le tabelle;",
                "controllare la tabella delle sezioni coinvolte (righe asciutte, pendenze, j_{den} nullo);",
                "se nasce nel riempimento, controllare la portata iniziale e lo stato iniziale (approfondimento 07)."])
    D.nota("costa poco: è una scansione dei vettori, molto più leggera di un passo del calcolo. Il caso noto di "
           "divergenza senza NaN è lo scambio dalla tabella senza serbatoio: il passo temporale va a zero ma i "
           "valori restano finiti (approfondimento 16).")
    D.collegamenti([(8, "08_interfacce_e_passo_temporale.pdf", "Il passo temporale"), (12, "12_newton_thomas.pdf", "Newton + Thomas"),
                    (16, "16_scambio_tabella.pdf", "Lo scambio dalla tabella")])


PDF = [
    (8, "08_interfacce_e_passo_temporale.pdf", "Interfacce e passo temporale", "Upwind, velocità di trasporto e condizione di Courant", "§10.1-10.2", p08),
    (9, "09_parte_esplicita_rusanov.pdf", "La parte esplicita", "Il flusso di Rusanov e le condizioni al contorno di monte e di valle", "§10.3", p09),
    (10, "10_attrito_implicito.pdf", "L'attrito implicito", "Perché il fattore 1/(1 + Δtγ) rende lo schema stabile", "§10.4", p10),
    (11, "11_sistema_dei_livelli.pdf", "Il sistema per i livelli", "Dalla quantità di moto e dalla continuità al sistema tridiagonale", "§10.5-10.6", p11),
    (12, "12_newton_thomas.pdf", "Newton + Thomas", "Il sistema non lineare dei livelli e l'algoritmo tridiagonale: `THOMAS_ALGORITHM!`", "§10.7", p12),
    (13, "13_portate_salvataggio_tabelle.pdf", "Portate, salvataggio, tabelle", "La chiusura del passo temporale", "§10.8-10.9", p13),
    (14, "14_controllo_nan.pdf", "Il controllo dei NaN", "`CHECK_NAN`: dove nasce un NaN e come leggerlo", "§10.10", p14),
]
