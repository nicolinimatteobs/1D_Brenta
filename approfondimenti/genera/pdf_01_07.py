# Approfondimenti 01-07: equazioni e schema, settings.txt, griglia, tabelle, idrogrammi, affluenti, stato iniziale.
import numpy as np
from stile import *

SEMI = "Source/Semi_Implicit_Brenta_2D_1D.jl"
SUB = "Source/subroutines.jl"
MOD = "Source/module.jl"


def tabella_sezione(k):
    T = np.loadtxt(os.path.join(PROGETTO, TABELLE), delimiter=",", skiprows=1)
    return T[T[:, 0] == k]


def geometria():
    return np.loadtxt(os.path.join(PROGETTO, GEOMETRIA), delimiter=",", skiprows=1)


# =============================================================================================== 01
def p01(D):
    D.h1("1. Che cosa risolve il modello")
    D.p("Il modello descrive la propagazione di una piena lungo un corso d'acqua schematizzato come una linea: "
        "in ogni sezione trasversale si considerano un solo livello del pelo libero η e una sola portata Q. "
        "Le incognite dipendono quindi da una sola coordinata spaziale x (l'ascissa lungo l'asse del fiume) "
        "e dal tempo t. Le equazioni sono quelle di **de Saint-Venant**, che esprimono la conservazione "
        "del volume e della quantità di moto in un tronco di canale di lunghezza infinitesima dx.")
    D.h2("Ipotesi")
    D.elenco(["**pressione idrostatica**: le accelerazioni verticali sono trascurabili, il livello è orizzontale lungo la sezione;",
              "**pendenza piccola** del fondo, così che la quota del pelo libero η = z_{b} + h (talweg più tirante);",
              "**velocità media di sezione** u = Q/A; la non uniformità della velocità nella sezione è riassunta dal "
              "coefficiente di quantità di moto β (Boussinesq);",
              "**attrito** come in moto uniforme (Gauckler-Strickler), con il coefficiente di scabrezza equivalente "
              "ricavato dal modello 2D e scritto nelle tabelle;",
              "gli **apporti laterali** (affluenti e scambio con la piana) entrano nella sola continuità, senza quantità di moto propria."])
    D.h2("Equazione di continuità")
    D.p("Il volume d'acqua nel tronco dx cambia nel tempo per la differenza tra la portata che entra e quella che "
        "esce e per gli apporti laterali q_{lat} [m²/s], positivi se l'acqua entra nel canale:")
    D.eq(r"\frac{\partial A}{\partial t}+\frac{\partial Q}{\partial x}=q_{lat}")
    D.p("A(η) è l'area bagnata, funzione del livello attraverso la forma della sezione (qui: attraverso le tabelle). "
        "Nel modello q_{lat} raccoglie gli affluenti (approfondimento 06) e lo scambio con la piana "
        "(approfondimenti 15-17).")
    D.h2("Equazione della quantità di moto")
    D.p("La variazione della portata nel tempo è dovuta al trasporto di quantità di moto, alla spinta idrostatica "
        "(gradiente del pelo libero) e all'attrito sul contorno:")
    D.eq(r"\frac{\partial Q}{\partial t}+\frac{\partial}{\partial x}\left(\beta\,\frac{Q^2}{A}\right)+g\,A\,\frac{\partial \eta}{\partial x}+g\,A\,J=0")
    D.p("J è la pendenza d'attrito. Con la formula di Gauckler-Strickler, scabrezza K_{s} [m^{1/3}/s] e raggio idraulico R_{h}:")
    D.eq(r"J=\frac{u\,|u|}{K_s^2\,R_h^{4/3}},\qquad u=\frac{Q}{A}")
    D.p("Il codice non usa K_{s} e R_{h} separatamente ma la loro combinazione, letta direttamente dalle tabelle:")
    D.eq(r"j_{den}=\frac{K_s^2\,R_h^{4/3}}{g}\quad\Rightarrow\quad g\,A\,J=\frac{|u|}{j_{den}}\,Q=\gamma\,Q,\qquad \gamma=\frac{|u|}{j_{den}}\;[\mathrm{s}^{-1}]")
    D.p("Il termine d'attrito diventa così **lineare in Q** con un coefficiente γ che dipende dallo stato: è la forma "
        "che permette di trattarlo in modo implicito senza iterazioni (approfondimento 10).")
    D.nota("Il moto uniforme si ottiene annullando le derivate: g A J = g A S, cioè u² = g·j_{den}·S e "
           "Q = A·√(g·j_{den}·S). È la formula di `UNIFORM_DISCHARGE`, usata per lo stato iniziale e per la "
           "condizione al contorno di valle (approfondimenti 07 e 09).")

    D.h1("2. L'idea dello schema semi-implicito")
    D.p("Le onde di piena in un fiume sono lente (si muovono con la velocità della corrente, pochi m/s), ma le "
        "equazioni contengono anche **onde di gravità** con celerità c = √(g·h), più veloci. Uno schema "
        "completamente esplicito dovrebbe rispettare la condizione di Courant sulla celerità più alta, "
        "Δt ≤ Δx/(|u| + c): con Δx = 50 m e h = 3 m sarebbe Δt ≤ 7 s circa. Uno schema completamente implicito "
        "non ha limiti ma richiede di risolvere ogni passo un grande sistema non lineare.")
    D.p("Lo schema semi-implicito (nella famiglia di Casulli) prende il meglio dei due:")
    D.elenco(["**esplicito** il termine convettivo ∂(βQ²/A)/∂x, che è non lineare ma lento: il passo è limitato solo "
              "dalla velocità della corrente (approfondimenti 08-09);",
              "**implicito** il gradiente del pelo libero g A ∂η/∂x, responsabile delle onde di gravità: così il passo "
              "non è limitato dalla celerità √(g h) (approfondimento 11);",
              "**implicito** l'attrito γ Q, che può essere molto forte nei tratti poco profondi (approfondimento 10)."])
    D.p("Sostituendo la quantità di moto discretizzata nella continuità si ottiene un sistema nei soli livelli del "
        "nuovo passo. È non lineare perché A dipende da η (le sezioni non sono rettangolari), ed è "
        "**tridiagonale** perché ogni cella parla solo con le due vicine: si risolve con il metodo di Newton e, "
        "a ogni iterazione, con l'algoritmo di Thomas (approfondimento 12).")
    D.h2("La griglia sfalsata")
    D.p("I livelli η_{i} sono definiti nelle **sezioni** (centri delle celle, i = 1 … imax), le portate Q_{i+1/2} "
        "nelle **interfacce** tra due celle (imax + 1 valori, la prima è l'ingresso da monte e l'ultima l'uscita). "
        "Nel codice l'interfaccia i+1/2 ha indice `i+1`: `q[1]` è l'ingresso, `q[imax+1]` l'uscita. Il dettaglio "
        "è nell'approfondimento 03.")
    D.h2("Le equazioni discrete in sintesi")
    D.p("Con n il passo corrente, n+1 il nuovo, Δx_{i} la lunghezza della cella e Δx_{i+1/2} la distanza tra le sezioni:")
    D.eq([r"Q^{*}_{i+1/2}=Q^n_{i+1/2}-\frac{\Delta t}{\Delta x_{i+1/2}}\left(F_{i+1}-F_i\right)",
          r"Q^{n+1}_{i+1/2}=\frac{Q^{*}_{i+1/2}-g\,A_{i+1/2}\,\frac{\Delta t}{\Delta x_{i+1/2}}\left(\eta^{n+1}_{i+1}-\eta^{n+1}_{i}\right)}{1+\Delta t\,\gamma_{i+1/2}}",
          r"A_i(\eta_i^{n+1})=A_i^n-\frac{\Delta t}{\Delta x_i}\left(Q^{n+1}_{i+1/2}-Q^{n+1}_{i-1/2}\right)+\frac{\Delta t}{\Delta x_i}\left(Q_{aff,i}+Q_{lat,i}\right)"])
    D.p("F è il flusso convettivo di Rusanov (approfondimento 09), Q* la portata «provvisoria» che contiene solo il "
        "trasporto, Q_{aff} e Q_{lat} le portate degli affluenti e dello scambio con la piana nella cella [m³/s]. "
        "La terza equazione, con la seconda sostituita dentro, è il sistema per i livelli (approfondimento 11).")
    D.idea("trasporto esplicito, gravità e attrito impliciti; il volume si conserva esattamente perché la continuità è "
           "scritta in forma di bilancio di cella: ciò che esce da una cella entra nella vicina.")

    D.h1("3. Il flowchart del programma")
    D.p("Il programma parte da `main.jl`, che include il modulo `RUN_SIMULATION` e chiama "
        "`SEMI_IMPLICIT_BRENTA_2D_1D!()`. La funzione legge tutto, prepara lo stato iniziale e poi esegue un unico "
        "ciclo temporale che serve prima per il **riempimento** (`IC = 1`, la portata di monte resta al primo valore "
        "fino all'equilibrio) e poi per la **simulazione** vera (`IC = 0`, si segue l'idrogramma).")
    D.flowchart(dict(dimensioni=(15, 21.5), w=6.2, h=0.95, nodi={
        "a": (5.5, 20.8, "main.jl → SEMI_IMPLICIT_BRENTA_2D_1D!()", "inizio", 7.0, 0.8),
        "b": (5.5, 19.6, "INITIAL_SETTINGS!\nlegge le 33 righe di settings.txt", "funz", 6.2),
        "c": (5.5, 18.3, "ALLOCATE_ARRAYS! · LOAD_GEOMETRY\nLOAD_SECTION_TABLES · SECTION_INDEX", "funz", 6.2),
        "d": (5.5, 17.0, "SPATIAL_DOMAIN!: griglia sfalsata\nLOAD_ALL_HYDROGRAPHS! · INIT_OUTPUT", "funz", 6.2),
        "e": (5.5, 15.7, "Stato iniziale: portate cumulate +\nUNIFORM_LEVEL in ogni sezione", "proc", 6.2),
        "f": (5.5, 14.4, "nCS < Nmax  e  t < t_end ?", "dec", 5.4, 1.1),
        "g": (5.5, 13.0, "Grandezze alle interfacce (upwind),\npasso Δt = min(Δt_max, CFL·Δx/max|βu|)", "proc", 6.2),
        "h": (5.5, 11.7, "Q* esplicito (Rusanov) + condizioni al\ncontorno; attrito implicito → A/(1+Δtγ)", "proc", 6.2),
        "i": (5.5, 10.4, "Coefficienti del sistema tridiagonale\n(aa, bb, cc, rhs)", "proc", 6.2),
        "l": (5.5, 9.1, "Affluenti e scambio con la piana\nnel termine noto: rhs += Δt/Δx·Q", "proc", 6.2),
        "m": (5.5, 7.8, "Newton + THOMAS_ALGORITHM!\nnuovi livelli η al passo n+1", "funz", 6.2),
        "n": (5.5, 6.5, "Nuove portate Q (passo n+1); salvataggio;\ntabelle ai nuovi livelli (READ_VALUES)", "proc", 6.2),
        "o": (5.5, 5.1, "IC = 1 ed equilibrio raggiunto ?", "dec", 5.4, 1.1),
        "p": (5.5, 3.7, "IC = 0, t = 0: inizia la simulazione", "proc", 6.2),
        "q": (12.2, 14.4, "BILANCIO_DI_MASSA\nbilancio.txt", "funz", 3.6),
        "r": (12.2, 12.9, "Fine", "fine", 2.4, 0.8),
        "n1": (12.2, 9.1, "CHECK_NAN dopo ogni fase:\nun NaN ferma il calcolo\ncon fase, passo, tempo, posizioni", "nota", 3.9, 1.3),
        "n2": (12.2, 5.1, "Durante il riempimento la\nportata di monte resta al primo\nvalore dell'idrogramma", "nota", 3.9, 1.3),
    }, archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("e", "f"), ("f", "g", "sì"), ("g", "h"), ("h", "i"),
              ("i", "l"), ("l", "m"), ("m", "n"), ("n", "o"), ("o", "p", "sì"),
              ("o.w", "f.w", "no", [(1.6, 5.1), (1.6, 14.4)]), ("p.w", "f.w", None, [(1.6, 3.7), (1.6, 14.4)]),
              ("f.e", "q.w", "no"), ("q", "r")]),
        didascalia="Flowchart di `SEMI_IMPLICIT_BRENTA_2D_1D!`: lettura, stato iniziale e ciclo temporale, comune a "
                   "riempimento e simulazione. I riquadri con le barre laterali sono chiamate a funzioni.")
    D.h2("Dove si trova ogni passo")
    D.tabella([["Passo", "Codice", "Approf."],
               ["Lettura della configurazione", "`INITIAL_SETTINGS!` (module.jl)", "02"],
               ["Geometria e griglia", "`LOAD_GEOMETRY`, `SPATIAL_DOMAIN!`", "03"],
               ["Tabelle delle sezioni", "`LOAD_SECTION_TABLES`, `READ_VALUES`, `FIND_ROW`", "04"],
               ["Idrogrammi", "`LOAD_HYDROGRAPH!`, `LAGRANGE_INTERPOLATION!`", "05"],
               ["Affluenti", "`TRIBUTARY_INFLOW!`", "06"],
               ["Stato iniziale", "`UNIFORM_LEVEL`, `UNIFORM_DISCHARGE`", "07"],
               ["Interfacce e passo temporale", "ciclo `while`, prima parte", "08"],
               ["Trasporto e condizioni al contorno", "calcolo di `qstar`", "09"],
               ["Attrito", "`γ`, `A_int_g`", "10"],
               ["Sistema per i livelli", "`aa`, `bb`, `cc`, `rhs`", "11"],
               ["Newton + Thomas", "`THOMAS_ALGORITHM!`", "12"],
               ["Portate, salvataggio, tabelle", "fine del ciclo", "13"],
               ["Controllo dei NaN", "`CHECK_NAN`", "14"],
               ["Scambio con la piana", "`FREE_OVERFLOW`, `READ_EXCHANGE`, `VOLUME_AT`", "15-17"],
               ["Fine del riempimento", "criteri di arresto", "18"],
               ["Bilancio di massa", "`BILANCIO_DI_MASSA`", "19"],
               ["File di uscita", "`INIT_OUTPUT`, `save_*`", "20"],
               ["Costruzione delle tabelle", "`calcolo_tabelle_daBM`", "21"]], [0.36, 0.48, 0.16])
    D.h2("Il punto di ingresso")
    D.estratto("main.jl", 1, titolo="main.jl")
    D.estratto("Source/Run_Simulation.jl", 1, titolo="Source/Run_Simulation.jl")
    D.collegamenti([(2, "02_lettura_settings.pdf", "Lettura di settings.txt"),
                    (8, "08_interfacce_e_passo_temporale.pdf", "Grandezze alle interfacce e passo temporale"),
                    (11, "11_sistema_dei_livelli.pdf", "Il sistema per i livelli")])


# =============================================================================================== 02
def p02(D):
    D.h1("A cosa serve")
    D.p("`INITIAL_SETTINGS!` legge `settings.txt` e restituisce una struttura `Parametri` con tutto quello che serve "
        "alla simulazione: i percorsi dei file di input, la cartella degli output, i parametri numerici e le "
        "scelte di modello (affluenti, scambio con la piana). È la prima funzione chiamata dal solutore.")
    D.p("Il file è letto **per posizione**: conta il numero di riga, non un'etichetta. Per questo le prime 33 righe "
        "non vanno mai spostate; dopo ci sono una legenda e vecchie configurazioni, che il modello non legge.")
    D.h1("Il contenuto delle 33 righe")
    D.tabella([["Riga", "Contenuto", "Uso nel codice"],
               ["1", "Tabelle delle sezioni (`tabella_sezioni.csv`)", "`parametri_sezioni` → approf. 04"],
               ["2", "Geometria (`localslope_xcoord.csv`)", "`local_slope` → approf. 03"],
               ["3-13", "Idrogrammi degli affluenti Brenta1-6, BrentaNuovo, BrentaVecchio1-4", "`hydro_files` → approf. 05-06"],
               ["14", "Idrogramma da monte (chiave `caldonazzo`)", "condizione di monte"],
               ["15-21", "Affluenti Ceggio, Levico, Centa, Chieppena, Grigno, Maso, Moggio", "`hydro_files`"],
               ["22", "Cartella degli output", "`output_dir` → approf. 20"],
               ["23", "Ks (non più usato, la riga deve restare)", "`Ks`"],
               ["24", "Tempo iniziale [s]", "`time`"],
               ["25", "Passo temporale massimo [s]", "`dtmax` → approf. 08"],
               ["26", "Numero massimo di passi", "`Nmax`"],
               ["27", "Affluenti attivi (si/no)", "`con_affluenti`"],
               ["28", "Scambio con la piana: sfioro, tabella, nessuno", "`scambio` → approf. 15-17"],
               ["29", "Serbatoio: soglia [m³]", "`V_soglia`"],
               ["30", "Serbatoio: tempo di svuotamento T [h]", "`T_piana` (in s)"],
               ["31", "Serbatoio: ritardo τ [h]", "`tau_piana` (in s)"],
               ["32", "Serbatoio: rilascio lento [h]", "`T_lento` (in s)"],
               ["33", "Serbatoio attivo (si/no)", "`piana_attiva`"]], [0.08, 0.55, 0.37])
    D.h1("Flowchart")
    D.flowchart(dict(dimensioni=(14, 13.5), w=5.6, h=0.9, nodi={
        "a": (5, 12.9, "INITIAL_SETTINGS!()", "inizio", 4.2, 0.75),
        "b": (5, 11.7, "Legge le prime 33 righe e toglie gli spazi", "io"),
        "c": (5, 10.5, "Righe 1-22: percorsi (stringhe)\nrighe 3-21 raccolte in hydro_files", "proc"),
        "d": (5, 9.2, "Righe 23-26: parse(Float64, …)\nNmax arrotondato a intero", "proc"),
        "e": (5, 7.9, "Riga 27 e 33: PARSE_SI_NO\nriga 28: PARSE_SCAMBIO", "funz"),
        "f": (5, 6.6, "Righe 29-32: PARSE_NON_NEGATIVO\n(ore × 3600 → secondi)", "funz"),
        "g": (5, 5.3, "FILE_SIZE: righe e colonne delle tabelle,\nnumero di sezioni dalla geometria", "funz"),
        "h": (5, 3.9, "scambio = tabella e\nmeno di 14 colonne ?", "dec", 4.8, 1.2),
        "i": (11, 3.9, "error(…)\nservono le colonne ΔQ", "fine", 3.2, 0.9),
        "l": (5, 2.4, "return Parametri(…)", "fine", 3.6, 0.75),
        "n": (11, 7.3, "Un valore non valido\n(testo al posto di un numero,\nnumero negativo, si/no errato)\nferma subito il programma\ncon il numero di riga", "nota", 4.0, 2.0),
    }, archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("e", "f"), ("f", "g"), ("g", "h"), ("h", "i", "sì"),
              ("h", "l", "no")]), didascalia="Flowchart di `INITIAL_SETTINGS!`.")
    D.h1("Il codice")
    D.estratto(MOD, "function INITIAL_SETTINGS!", "        error(\"Riga 28", titolo=MOD + " (prima parte)")
    D.estratto(MOD, "    return Parametri(", n=4, titolo=MOD + " (fine)")
    D.p("Tre funzioni piccole controllano i valori e danno un messaggio che dice **quale riga** è sbagliata:")
    D.estratto(MOD, "# Legge un numero non negativo", "# Conta le righe di dati", titolo=MOD)
    D.h2("FILE_SIZE: le dimensioni dai file")
    D.p("Il numero di sezioni, il numero di righe delle tabelle e il numero di colonne non sono scritti in "
        "`settings.txt` ma ricavati dai file: le colonne dall'intestazione, le righe contando quelle non vuote. "
        "Così non è possibile che configurazione e dati siano incoerenti.")
    D.estratto(MOD, "function FILE_SIZE", n=8, titolo=MOD)
    D.h2("La struttura Parametri")
    D.p("I valori sono raccolti in una `struct` immutabile; nel solutore si estraggono con la macro `@unpack` "
        "(pacchetto UnPack), che crea variabili locali con lo stesso nome dei campi.")
    D.estratto("Source/struct.jl", 1, titolo="Source/struct.jl")
    D.h1("Perché la conversione delle unità")
    D.p("Nelle righe 30-32 i tempi del serbatoio della piana si scrivono in **ore**, più comodi da leggere e da "
        "tarare; il modello lavora in secondi, quindi sono moltiplicati per 3600 alla lettura. La soglia (riga 29) "
        "è già in m³. Tutte e quattro devono essere numeri non negativi.")
    D.attenzione("gli strumenti di `strumenti/` leggono solo le prime 26 righe di `settings.txt`; tabelle e "
                 "geometria vengono dalle righe 1-2, quindi devono essere quelle della simulazione che si analizza.")
    D.collegamenti([(1, "01_equazioni_e_schema.pdf", "Equazioni e schema"),
                    (3, "03_geometria_e_griglia.pdf", "Geometria e griglia"),
                    (17, "17_serbatoio_della_piana.pdf", "Il serbatoio della piana")])


# =============================================================================================== 03
def p03(D):
    G = geometria()
    D.h1("A cosa serve")
    D.p("Il fiume è descritto da imax sezioni trasversali (nel Test_Brenta 120, ogni 50 m). Il file della geometria "
        "(riga 2 di `settings.txt`) dà per ogni sezione l'ascissa, la quota del fondo (talweg), la pendenza locale "
        "e le quote degli argini. `LOAD_GEOMETRY` lo legge e `SPATIAL_DOMAIN!` costruisce la **griglia sfalsata**: "
        "le celle attorno alle sezioni e le interfacce tra le celle.")
    D.h1("Il file della geometria")
    D.tabella([["Colonna", "Nome", "Significato", "Variabile"],
               ["1", "`section`", "numero progressivo da monte", "–"],
               ["2", "`local_slope`", "pendenza del fondo, negativa verso valle", "`bx`"],
               ["3", "`x_coord_[m]`", "ascissa lungo l'asse", "`x_b`"],
               ["4", "`talweg_[m]`", "quota del punto più basso", "`b_min`"],
               ["5-6", "`argine_dx_[m]`, `argine_sx_[m]`", "quote delle creste arginali (facoltative)", "`argine_dx`, `argine_sx`"]],
              [0.12, 0.3, 0.38, 0.2])
    D.p(f"Nel Test_Brenta a 50 m: x da {G[0,2]:.2f} a {G[-1,2]:.0f} m, talweg da {G[0,3]:.2f} a {G[-1,3]:.2f} m s.l.m., "
        f"pendenza media del fondo {(G[-1,3]-G[0,3])/(G[-1,2]-G[0,2]):.5f}. Con 4 colonne gli argini valgono "
        "`Inf` e non c'è mai sfioro.")
    D.estratto(SUB, "function LOAD_GEOMETRY", n=8, titolo=SUB)
    D.h1("La griglia sfalsata")
    D.p("Livelli e portate non sono definiti negli stessi punti. Il **livello** η_{i} sta nella sezione i, che è il "
        "centro della cella i; la **portata** sta sulle **interfacce** tra le celle. In questo modo la differenza "
        "di livello tra due sezioni vicine, che spinge l'acqua, è centrata proprio dove si calcola la portata, e la "
        "differenza tra le portate di due interfacce è centrata nella cella di cui si calcola il volume.")
    def disegno(fig):
        ax = fig.add_subplot(111); ax.set_xlim(-0.3, 6.3); ax.set_ylim(-1.5, 1.7); ax.axis("off")
        xs = [0.5, 1.7, 3.0, 4.4, 5.6]
        xi = [0.0] + [0.5 * (xs[k] + xs[k + 1]) for k in range(4)] + [6.2]
        for k, x in enumerate(xi):
            ax.plot([x, x], [-0.5, 0.5], color="#eb6834", lw=1.6)
            lab = "q[1]\ningresso" if k == 0 else ("q[imax+1]\nuscita" if k == 5 else f"q[{k+1}]")
            ax.text(x, 0.62, lab, ha="center", va="bottom", fontsize=8, color="#b4371f")
        ax.plot([xi[0], xi[-1]], [0, 0], color="#9aa0a6", lw=0.8)
        for k, x in enumerate(xs):
            ax.plot(x, 0, "o", color="#2a78d6", ms=7)
            ax.text(x, -0.28, f"η[{k+1}]", ha="center", va="top", fontsize=8.5, color="#1f2328")
        ax.annotate("", xy=(xs[1], -0.95), xytext=(xs[2], -0.95), arrowprops=dict(arrowstyle="<->", color="#2a78d6"))
        ax.text(0.5 * (xs[1] + xs[2]), -1.05, "dx_b[2]: distanza tra le sezioni", ha="center", va="top", fontsize=8, color="#2a78d6")
        ax.annotate("", xy=(xi[2], 1.32), xytext=(xi[3], 1.32), arrowprops=dict(arrowstyle="<->", color="#b4371f"))
        ax.text(0.5 * (xi[2] + xi[3]), 1.4, "dx_int[3]: lunghezza della cella 3", ha="center", va="bottom", fontsize=8, color="#b4371f")
        ax.text(3.1, -1.45, "monte  →  valle", ha="center", fontsize=8, color="#5f6368")
    D.grafico(disegno, "griglia", "Griglia sfalsata con 5 sezioni: livelli nelle sezioni (punti blu), portate nelle "
              "interfacce (trattini arancioni). Le interfacce sono a metà tra due sezioni.")
    D.p("In formule, con x_{i} l'ascissa della sezione i:")
    D.eq([r"\Delta x_{i+1/2}=x_{i+1}-x_i \qquad (\mathtt{dx\_b[i]},\ i=1\ldots imax-1)",
          r"x_{i-1/2}=\frac{x_{i-1}+x_i}{2} \qquad (\mathtt{x\_int[i]},\ i=2\ldots imax)",
          r"\Delta x_i=\frac{\Delta x_{i-1/2}+\Delta x_{i+1/2}}{2},\qquad \Delta x_1=\Delta x_{3/2},\qquad \Delta x_{imax}=\Delta x_{imax-1/2}"])
    D.p("Le due celle di estremità hanno la lunghezza della distanza con la sezione vicina: la prima e l'ultima "
        "interfaccia sono quindi «fuori» dal tratto, alla stessa distanza delle altre. Per le sezioni equidistanti "
        "del Test_Brenta tutte le celle sono lunghe 50 m.")
    D.h2("Due lunghezze, due usi")
    D.elenco(["`dx_int[i]` (Δx_{i}) divide i flussi di volume nella **continuità** della cella i e moltiplica le aree "
              "per il volume in alveo (bilancio di massa, approfondimento 19);",
              "`dx_b[i]` (Δx_{i+1/2}) divide le differenze di livello e di flusso nella **quantità di moto** "
              "dell'interfaccia i+1/2, cioè nel calcolo di Q* e di Q^{n+1} (approfondimenti 09-11);",
              "la più piccola delle due, `dx`, entra nella condizione di Courant (approfondimento 08)."])
    D.h1("Il codice")
    D.estratto(SUB, "function SPATIAL_DOMAIN!", "    return dx_b, dx_int, x_int", titolo=SUB)
    D.h1("La pendenza locale")
    D.p("`bx` (colonna 2) non serve al moto vario: le pendenze del fondo sono già dentro i livelli delle tabelle. "
        "Serve solo per il **moto uniforme**: lo stato iniziale (`UNIFORM_LEVEL`, approfondimento 07) e la portata "
        "in uscita nell'ultima sezione (approfondimento 09). Per questo l'ultima sezione usa la pendenza media del "
        "fondo su tutto il tratto, scelta in `calcolo_tabelle_daBM` (approfondimento 21).")
    n_contro = int((G[:-1, 1] > 0).sum())
    D.attenzione(f"a 50 m le pendenze locali sono molto irregolari: {n_contro} sezioni hanno il fondo in "
                 "contropendenza. Il codice usa sempre il valore assoluto, quindi in quelle sezioni lo stato iniziale "
                 "è solo una prima stima, che il riempimento corregge.")
    D.collegamenti([(4, "04_tabelle_delle_sezioni.pdf", "Le tabelle delle sezioni"),
                    (11, "11_sistema_dei_livelli.pdf", "Il sistema per i livelli"),
                    (15, "15_sfioro_argini.pdf", "Lo sfioro sugli argini")])


# =============================================================================================== 04
def p04(D):
    k = 20; r = tabella_sezione(k); G = geometria()
    D.h1("A cosa serve")
    D.p("Il modello non conosce la forma geometrica delle sezioni. Per ogni sezione ha una **tabella** che, al variare "
        "della quota del pelo libero η, dà le grandezze idrauliche: area bagnata A, coefficiente β, il termine "
        "d'attrito j_{den}, la larghezza in superficie W e il raggio idraulico R_{h}. Le tabelle vengono dal modello 2D "
        "(approfondimento 21) e sono lette a ogni passo e a ogni iterazione di Newton da `READ_VALUES`, la funzione "
        "più chiamata del programma.")
    D.h1("Il formato")
    D.tabella([["Col.", "Nome", "Significato", "Unità"],
               ["1", "`sezione`", "numero della sezione (righe raggruppate per sezione)", "–"],
               ["2", "`WSE`", "quota del pelo libero η (non il tirante), crescente", "m s.l.m."],
               ["3-4", "`Ω`, `Ω_s`", "area bagnata e sua pendenza dA/dη", "m², m"],
               ["5-6", "`β`, `β_s`", "coefficiente di quantità di moto e pendenza", "–, 1/m"],
               ["7-8", "`jden`, `jden_s`", "K_{s}²·R_{h}^{4/3}/g e pendenza", "m, –"],
               ["9-10", "`W`, `W_s`", "larghezza bagnata e pendenza", "m, –"],
               ["11-12", "`RH`, `RH_s`", "raggio idraulico e pendenza", "m, –"],
               ["13-14", "`ΔQ`, `ΔQ_s`", "scambio con la piana dagli scalini 2D (facoltative)", "m³/s, m²/s"]],
              [0.1, 0.2, 0.52, 0.18])
    D.p(f"Le tabelle di `S_50m_q2` hanno da 44 a 52 righe per sezione, una per ogni scalino di portata del 2D. "
        f"Per esempio la sezione {k} (S1000) va da η = {r[0,1]:.2f} a {r[-1,1]:.2f} m: {len(r)} righe, area da "
        f"{r[0,2]:.2f} a {r[-1,2]:.1f} m², j_{{den}} da {r[0,6]:.2f} a {r[-1,6]:.1f} m.")
    def graf(fig):
        ax = fig.subplots(1, 3, sharey=True)
        for a, (c, nome, u) in zip(ax, [(2, "area A", "m²"), (8, "larghezza W", "m"), (6, "j_den", "m")]):
            a.plot(r[:, c], r[:, 1], "-o", ms=2.5, color="#2a78d6", lw=1.2)
            a.set_xlabel(f"{nome} [{u}]", fontsize=8); a.grid(alpha=0.25); a.tick_params(labelsize=7)
        ax[0].set_ylabel("livello η [m s.l.m.]", fontsize=8)
        fig.suptitle("Tabella della sezione 20 (S1000)", fontsize=9)
    D.grafico(graf, "tabella_S1000", "Tre colonne della tabella di S1000 (`S_50m_q2`): ogni punto è uno scalino di "
              "portata del modello 2D.")
    D.h1("L'interpolazione")
    D.p("Tra due righe r e r+1 con η_{r} ≤ η < η_{r+1} ogni grandezza φ è interpolata linearmente. La pendenza è "
        "già scritta nella tabella (colonne pari), così l'interpolazione costa una moltiplicazione e una somma:")
    D.eq(r"\varphi(\eta)=\varphi_r+\varphi'_r\,(\eta-\eta_r),\qquad \varphi'_r=\frac{\varphi_{r+1}-\varphi_r}{\eta_{r+1}-\eta_r}")
    D.p("Agli estremi:")
    D.elenco(["**sotto la prima quota** (il fondo) si usa la prima riga: η è sostituito da max(η, η_{1});",
              "**oltre l'ultima quota** si estrapola con la pendenza dell'ultimo tratto, calcolata dai valori delle ultime "
              "due righe:"])
    D.eq(r"\varphi(\eta)=\varphi_N+\frac{\varphi_N-\varphi_{N-1}}{\eta_N-\eta_{N-1}}\,(\eta-\eta_N),\qquad \eta>\eta_N")
    D.p("L'estrapolazione serve perché durante le piene più forti il livello può superare l'ultimo scalino 2D (al colmo, "
        "in S2100-S2150, di qualche centimetro): fermarsi all'ultima riga renderebbe l'area costante e farebbe "
        "crescere il livello senza limite. Prolungare la retta tiene area e larghezza continue.")
    D.p("La larghezza W = dA/dη è anche la **derivata** che serve a Newton (approfondimento 12): la colonna 9 è la "
        "larghezza bagnata del 2D, che coincide con la pendenza dell'area a meno degli arrotondamenti.")
    D.h1("La ricerca della riga: FIND_ROW")
    D.p("Per sapere in quale intervallo cade η bisogna cercare la riga. Una ricerca lineare costerebbe O(N) a ogni "
        "chiamata; la ricerca binaria O(log N). `FIND_ROW` sfrutta il fatto che tra un passo e l'altro il livello "
        "cambia poco: prova prima l'**ultima riga usata** g (salvata in `last_row`) e le due vicine, e solo se η "
        "è più lontano usa la ricerca binaria. Quasi sempre la risposta è immediata.")
    D.flowchart(dict(dimensioni=(15, 12.2), w=4.6, h=0.9, nodi={
        "a": (5.5, 11.6, "FIND_ROW(ys, lo, hi, y, g)", "inizio", 4.6, 0.75),
        "b": (5.5, 10.2, "ys[g] ≤ y ?", "dec", 3.6, 1.1),
        "c": (5.5, 8.6, "g = hi  oppure  y < ys[g+1] ?", "dec", 5.0, 1.1),
        "d": (11.5, 8.6, "return g", "fine", 2.4, 0.75),
        "e": (5.5, 7.0, "g+1 = hi ?", "dec", 3.4, 1.0),
        "e2": (11.5, 7.0, "return hi", "fine", 2.4, 0.75),
        "f": (5.5, 5.5, "y < ys[g+2] ?", "dec", 3.6, 1.0),
        "f2": (11.5, 5.5, "return g+1", "fine", 2.4, 0.75),
        "h": (1.6, 8.6, "ys[g−1] ≤ y ?", "dec", 2.9, 1.0),
        "h2": (1.6, 7.0, "return g−1", "fine", 2.4, 0.75),
        "i": (5.5, 3.9, "y ≥ ys[hi] ?", "dec", 3.4, 1.0),
        "i2": (11.5, 3.9, "return hi", "fine", 2.4, 0.75),
        "l": (5.5, 2.4, "ricerca binaria tra lo e hi", "proc", 4.4),
        "m": (5.5, 1.1, "return lo", "fine", 2.4, 0.75),
    }, archi=[("a", "b"), ("b", "c", "sì"), ("c", "d", "sì"), ("c", "e", "no"), ("e", "e2", "sì"), ("e", "f", "no"),
              ("f", "f2", "sì"), ("f", "i", "no"), ("b.w", "h.n", "no", [(1.6, 10.2)]), ("h", "h2", "sì"),
              ("h.w", "i.w", "no", [(0.4, 8.6), (0.4, 3.9)]), ("i", "i2", "sì"), ("i", "l", "no"), ("l", "m")]),
        didascalia="Flowchart di `FIND_ROW`: prima i tentativi vicini all'ultima riga usata, poi la ricerca binaria.")
    D.estratto(SUB, "@inline function FIND_ROW", n=17, titolo=SUB)
    D.h1("READ_VALUES")
    D.estratto(SUB, "@inline function READ_VALUES", n=24, titolo=SUB)
    D.h2("La tabella trasposta")
    D.p("Il file ha una riga per livello e una colonna per grandezza; il solutore lo **traspone** "
        "(`data_depth = permutedims(...)`, colonna = riga del file). Julia memorizza le matrici per colonne, quindi i "
        "valori di uno stesso livello sono contigui in memoria e vengono letti insieme. `SECTION_INDEX` trova dove "
        "inizia ogni sezione (`first_row`), così la sezione k occupa le colonne da `first_row[k]` a `first_row[k+1]-1`.")
    D.estratto(SUB, "function SECTION_INDEX", n=5, titolo=SUB)
    D.h2("Un esempio numerico")
    j = 20; eta = 0.5 * (r[j, 1] + r[j + 1, 1])
    A = r[j, 2] + r[j, 3] * (eta - r[j, 1])
    D.p(f"In S1000, per η = {eta:.3f} m la riga trovata è la {j+1}ª (η_{{r}} = {r[j,1]:.3f} m, A_{{r}} = {r[j,2]:.3f} m², "
        f"pendenza {r[j,3]:.3f} m): A = {r[j,2]:.3f} + {r[j,3]:.3f}·{eta-r[j,1]:.3f} = **{A:.3f} m²**; "
        f"la riga successiva ha A = {r[j+1,2]:.3f} m².")
    D.h1("Le colonne ΔQ")
    D.p("Le tabelle a 14 colonne hanno in più lo scambio con la piana ΔQ(η), letto da una funzione separata, "
        "`READ_EXCHANGE`, che **non** estrapola oltre l'ultima quota (approfondimento 16). Le tabelle a 12 colonne "
        "restano valide: `LOAD_SECTION_TABLES` aggiunge due colonne di zeri.")
    D.estratto(SUB, "function LOAD_SECTION_TABLES", n=7, titolo=SUB)
    D.collegamenti([(12, "12_newton_thomas.pdf", "Newton + Thomas (dove si usano A e W)"),
                    (16, "16_scambio_tabella.pdf", "Lo scambio dalla tabella ΔQ"),
                    (21, "21_costruzione_tabelle.pdf", "Come si costruiscono le tabelle")])


# =============================================================================================== 05
def p05(D):
    with open(os.path.join(PROGETTO, "input/idrogrammi/fittizio_1/piena.csv")) as f:
        I = np.array([[float(v) for v in l.split(",")] for l in f if l.strip()])
    D.h1("A cosa serve")
    D.p("Gli idrogrammi danno la portata nel tempo: quello di monte (riga 14 di `settings.txt`) è la condizione al "
        "contorno di ingresso, gli altri 18 sono gli affluenti. Sono tabelle di valori ogni Δτ (900 s negli scenari "
        "fittizi), mentre il modello avanza con passi Δt irregolari di 20-30 s: serve un modo per conoscere la "
        "portata a qualunque istante. `LAGRANGE_INTERPOLATION!` la calcola con un polinomio di secondo grado.")
    D.h1("Lettura")
    D.p("`LOAD_HYDROGRAPH!` legge due colonne, tempo [s] e portata [m³/s], separate da virgole, spazi o tabulazioni, "
        "senza intestazione. `LOAD_ALL_HYDROGRAPHS!` le carica tutte in un dizionario con le chiavi dei nomi "
        "(`hydro[:caldonazzo]`, `hydro[:centa]`, …). La durata simulata è il tempo finale **più breve**:")
    D.eq(r"t_{end}=\min_j\,t^{(j)}_{N_j}")
    D.estratto(SUB, "function HYDROGRAPH_END_TIME", n=8, titolo=SUB)
    D.p(f"L'idrogramma del `fittizio_1` ha {len(I)} valori ogni {I[1,0]-I[0,0]:.0f} s, da {I[0,1]:.2f} m³/s a "
        f"un colmo di {I[:,1].max():.2f} m³/s a {I[I[:,1].argmax(),0]/3600:.0f} h, per {I[-1,0]/3600:.0f} h.")
    D.h1("L'interpolazione di Lagrange su tre punti")
    D.p("La portata serve al tempo t = t^{n} + Δt (fine del passo: la condizione di monte è imposta nel nuovo stato). "
        "Si scelgono tre punti consecutivi dell'idrogramma attorno a t e si usa il polinomio di grado 2 che passa "
        "per essi:")
    D.eq([r"Q(t)=Q_{m-1}\,L_{m-1}(t)+Q_m\,L_m(t)+Q_{m+1}\,L_{m+1}(t)",
          r"L_{m-1}=\frac{(t-t_m)(t-t_{m+1})}{(t_{m-1}-t_m)(t_{m-1}-t_{m+1})},\quad L_m=\frac{(t-t_{m-1})(t-t_{m+1})}{(t_m-t_{m-1})(t_m-t_{m+1})},\quad L_{m+1}=\frac{(t-t_{m-1})(t-t_m)}{(t_{m+1}-t_{m-1})(t_{m+1}-t_m)}"])
    D.p("Ogni L vale 1 nel suo punto e 0 negli altri due, quindi il polinomio riproduce esattamente i valori "
        "dell'idrogramma nei punti dati. L'indice centrale m si ricava supponendo il passo dell'idrogramma "
        "**costante** e il primo tempo **nullo**:")
    D.eq(r"\Delta\tau=\frac{t_N-t_1}{N-1},\qquad m=\min\left(N-1,\ \max\left(2,\ \left\lceil\frac{t^n+\Delta t}{\Delta\tau}\right\rceil\right)\right)")
    D.p("Rispetto all'interpolazione lineare, la parabola segue la curvatura dell'idrogramma vicino al colmo e non "
        "lo «taglia»: con valori ogni 15 minuti la differenza al colmo è piccola ma sistematica.")
    def graf(fig):
        ax = fig.add_subplot(111)
        i0 = int(I[:, 1].argmax())
        sel = I[i0 - 4:i0 + 5]
        ax.plot(sel[:, 0] / 3600, sel[:, 1], "o", color="#1f2328", ms=5, label="valori dell'idrogramma")
        tt = np.linspace(sel[0, 0], sel[-1, 0], 400)
        dτ = (I[-1, 0] - I[0, 0]) / (len(I) - 1); N = len(I)
        def lag(t):
            m = int(min(N - 1, max(2, np.ceil(t / dτ)))) - 1
            t0, t1, t2 = I[m - 1:m + 2, 0]; q0, q1, q2 = I[m - 1:m + 2, 1]
            return (q0 * (t - t1) * (t - t2) / ((t0 - t1) * (t0 - t2)) + q1 * (t - t0) * (t - t2) / ((t1 - t0) * (t1 - t2))
                    + q2 * (t - t0) * (t - t1) / ((t2 - t0) * (t2 - t1)))
        ax.plot(tt / 3600, [lag(t) for t in tt], color="#2a78d6", lw=1.5, label="Lagrange su tre punti (modello)")
        ax.plot(tt / 3600, np.interp(tt, I[:, 0], I[:, 1]), "--", color="#eb6834", lw=1.2, label="interpolazione lineare")
        ax.set_xlabel("tempo [h]", fontsize=8); ax.set_ylabel("portata [m³/s]", fontsize=8)
        ax.grid(alpha=0.25); ax.tick_params(labelsize=7); ax.legend(fontsize=7, loc="lower center")
    D.grafico(graf, "lagrange", "Interpolazione dell'idrogramma del `fittizio_1` attorno al colmo, con la stessa "
              "scelta dei punti del modello.")
    D.estratto(SUB, "function LAGRANGE_INTERPOLATION!", n=16, titolo=SUB)
    D.flowchart(dict(dimensioni=(12, 7.6), w=6.4, h=0.95, nodi={
        "a": (6, 7.0, "LAGRANGE_INTERPOLATION!(nr, tk, Qk, time, dt)", "inizio", 7.4, 0.75),
        "b": (6, 5.75, "Δτ = (t_N − t_1)/(N − 1)\nm = ⌈(time + dt)/Δτ⌉ limitato a [2, N − 1]", "proc"),
        "c": (6, 4.35, "Punti (t_{m−1}, Q_{m−1}), (t_m, Q_m), (t_{m+1}, Q_{m+1})", "proc"),
        "d": (6, 3.0, "Q = Σ Q_k · L_k(time + dt)", "proc"),
        "e": (6, 1.8, "return Q", "fine", 2.6, 0.75)},
        archi=[("a", "b"), ("b", "c"), ("c", "d"), ("d", "e")]),
        didascalia="Flowchart di `LAGRANGE_INTERPOLATION!`.", larghezza=0.8)
    D.attenzione("l'indice m è calcolato supponendo che l'idrogramma parta da t = 0 e abbia passo costante. Con un "
                 "file a passo irregolare i tre punti scelti non sarebbero quelli attorno a t. Per idrogrammi diversi "
                 "dai fittizi conviene controllarlo. Oltre l'ultimo punto la parabola estrapola: per questo la "
                 "simulazione si ferma a t_{end}.")
    D.h1("Dove si usa")
    D.elenco(["nella simulazione, per la portata in ingresso `qstar[1]` (approfondimento 09);",
              "per i 18 affluenti, se la riga 27 è `si` (approfondimento 06);",
              "**non** durante il riempimento: lì la portata di monte resta al primo valore dell'idrogramma."])
    D.collegamenti([(6, "06_affluenti.pdf", "Gli affluenti"), (9, "09_parte_esplicita_rusanov.pdf", "Condizioni al contorno")])


# =============================================================================================== 06
def p06(D):
    D.h1("A cosa serve")
    D.p("Lungo il Brenta entrano 18 affluenti. Il modello li tratta come **sorgenti puntuali di volume**: la loro "
        "portata entra nella continuità delle celle in cui sfociano, senza quantità di moto propria (la corrente "
        "laterale arriva in genere quasi perpendicolare all'asse e la sua spinta lungo l'asse è trascurabile). Si "
        "accendono e spengono con la riga 27 di `settings.txt`; negli scenari fittizi sono spenti.")
    D.h1("Dove entrano")
    D.p("Non avendo ancora le posizioni reali, gli affluenti sono distribuiti uniformemente lungo il tratto:")
    D.eq(r"i_{aff}^{(k)}=\mathrm{round}\left(\frac{k\cdot imax}{19}\right),\qquad k=1\ldots 18")
    D.p("Gli affluenti più grandi, con picco di almeno 10 m³/s, sono distribuiti su **tre** sezioni (quella di "
        "sbocco e le due vicine), gli altri su una sola, così una portata grande non si concentra in una cella "
        "lunga 50 m:")
    D.eq(r"Q_{aff,i}=\sum_{k:\ |i-i_{aff}^{(k)}|\leq s_k}\frac{Q^{(k)}(t)}{2\,s_k+1},\qquad s_k=1\ \mathrm{se}\ \max_t Q^{(k)}\geq 10\ \mathrm{m^3/s},\ \mathrm{altrimenti}\ 0")
    D.estratto(SEMI, "    # Sezioni di ingresso degli affluenti", n=10, titolo=SEMI)
    D.estratto(SUB, "function TRIBUTARY_INFLOW!", n=9, titolo=SUB)
    D.h1("Come entrano nelle equazioni")
    D.p("Il termine sorgente compare solo nella continuità della cella i. Integrato sulla cella e sul passo:")
    D.eq(r"A_i^{n+1}=A_i^n-\frac{\Delta t}{\Delta x_i}\left(Q_{i+1/2}^{n+1}-Q_{i-1/2}^{n+1}\right)+\frac{\Delta t}{\Delta x_i}\,Q_{aff,i}")
    D.p("Nel codice si somma al termine noto `rhs` del sistema per i livelli (approfondimento 11), dopo averne "
        "calcolato la parte idrodinamica. Durante il riempimento la portata è il primo valore di ogni idrogramma, "
        "in simulazione quella interpolata a t + Δt con Lagrange (approfondimento 05):")
    D.estratto(SEMI, "        ####  INSERIMENTO AFFLUENTI", "        ####  FINE INSERIMENTO AFFLUENTI", titolo=SEMI)
    D.flowchart(dict(dimensioni=(13, 9.6), w=5.4, h=0.95, nodi={
        "a": (5.5, 9.0, "Inserimento affluenti (ogni passo)", "inizio", 5.6, 0.75),
        "b": (5.5, 7.7, "con_affluenti ?", "dec", 3.4, 1.0),
        "c": (5.5, 6.2, "IC = 1 ?", "dec", 2.6, 1.0),
        "d": (2.4, 4.6, "Q^{(k)} = primo valore\ndell'idrogramma", "proc", 3.9),
        "e": (8.6, 4.6, "Q^{(k)} = LAGRANGE_\nINTERPOLATION!(t + Δt)", "funz", 3.9),
        "f": (5.5, 3.1, "TRIBUTARY_INFLOW!: Q_aff in ogni cella", "funz"),
        "g": (5.5, 1.8, "rhs[i] += Δt/Δx_i · Q_aff[i]", "proc"),
        "h": (11.3, 7.7, "nessun affluente", "fine", 2.6, 0.75)},
        archi=[("a", "b"), ("b", "c", "sì"), ("b", "h", "no"), ("c.w", "d.n", "sì", [(2.4, 6.2)]), ("c.e", "e.n", "no", [(8.6, 6.2)]),
               ("d.s", "f.w", None, [(2.4, 3.1)]), ("e.s", "f.e", None, [(8.6, 3.1)]), ("f", "g")]),
        didascalia="Flowchart dell'inserimento degli affluenti a ogni passo.", larghezza=0.9)
    D.h1("Nello stato iniziale e nel bilancio")
    D.p("Lo stato iniziale parte dalle **portate cumulate**: la portata all'interfaccia i+1/2 è quella di monte più "
        "tutti gli affluenti entrati più a monte, Q_{i+1/2} = Q_{i−1/2} + Q_{aff,i} (approfondimento 07). Il "
        "volume entrato dagli affluenti è sommato passo per passo e compare in `bilancio.txt` "
        "(approfondimento 19); il riempimento si ferma solo quando l'uscita è uguale a monte + affluenti "
        "(approfondimento 18).")
    D.attenzione("le posizioni uniformi sono provvisorie. Affluenti e scambio con la piana insieme non sono ancora "
                 "stati provati.")
    D.collegamenti([(5, "05_idrogrammi.pdf", "Gli idrogrammi"), (7, "07_stato_iniziale.pdf", "Lo stato iniziale"),
                    (11, "11_sistema_dei_livelli.pdf", "Il sistema per i livelli")])


# =============================================================================================== 07
def p07(D):
    G = geometria(); k = 20; r = tabella_sezione(k)
    D.h1("A cosa serve")
    D.p("Prima della piena il fiume deve trovarsi in uno stato di equilibrio con la portata iniziale. Il modello "
        "lo costruisce in due tempi: prima una **stima** sezione per sezione con il moto uniforme (questo "
        "approfondimento), poi un **riempimento** in cui il modello gira con portata costante finché tutto si "
        "assesta (approfondimento 18).")
    D.h1("Le portate cumulate")
    D.p("La portata di partenza all'interfaccia di monte è il primo valore dell'idrogramma; scendendo verso valle "
        "si aggiungono gli affluenti (se attivi). In ogni sezione si usa la media delle due interfacce:")
    D.eq([r"Q_{1/2}=Q_{monte}(t=0),\qquad Q_{i+1/2}=Q_{i-1/2}+Q_{aff,i}", r"\overline{Q}_i=\frac{Q_{i-1/2}+Q_{i+1/2}}{2}"])
    D.h1("Il moto uniforme con l'attrito delle tabelle")
    D.p("In moto uniforme il livello è parallelo al fondo e la pendenza d'attrito è uguale a quella del fondo, J = S. "
        "Con l'attrito scritto come nel modello (approfondimento 01):")
    D.eq(r"g\,A\,S=\frac{|u|}{j_{den}}\,Q\quad\Rightarrow\quad u^2=g\,j_{den}\,S\quad\Rightarrow\quad Q_u(\eta)=A(\eta)\,\sqrt{g\,j_{den}(\eta)\,|S|}")
    D.p("Poiché j_{den} = K_{s}²R_{h}^{4/3}/g, questa è la formula di Gauckler-Strickler, Q = A K_{s} R_{h}^{2/3} √S, con "
        "il K_{s} equivalente della tabella. È esattamente l'attrito usato in moto vario: lo stato iniziale è quindi "
        "coerente con il resto del calcolo.")
    D.estratto(SUB, "UNIFORM_DISCHARGE(A, den, S)", n=1, titolo=SUB)
    D.h1("La bisezione: UNIFORM_LEVEL")
    D.p("Serve il livello η tale che Q_{u}(η) = Q̄_{i}. Non c'è una formula inversa perché A e j_{den} sono "
        "tabellati, quindi si cerca per **bisezione** nell'intervallo tra il fondo (più un tirante minimo y_{min}) "
        "e l'ultima quota della tabella:")
    D.eq([r"a_0=\eta_{fondo}+y_{min},\qquad b_0=\eta_{N}",
          r"m=\frac{a+b}{2}:\quad Q_u(m)<\overline{Q}\ \Rightarrow\ a=m,\qquad \mathrm{altrimenti}\ b=m",
          r"|\eta-\eta^{*}|\leq\frac{b_0-a_0}{2^{60}}\approx 10^{-18}\ \mathrm{m}\quad (60\ \mathrm{iterazioni})"])
    D.p("Se già al tirante minimo la portata di moto uniforme supera Q̄, la funzione restituisce quel livello. Nel "
        "solutore y_{min} = `y_ini` = 0, quindi nessun tirante minimo.")
    D.flowchart(dict(dimensioni=(13, 11.2), w=5.4, h=0.95, nodi={
        "a": (5.5, 10.6, "UNIFORM_LEVEL(…, k, Q, S; y_min)", "inizio", 5.6, 0.75),
        "b": (5.5, 9.4, "a = fondo + y_min", "proc"),
        "c": (5.5, 8.0, "Q_u(a) ≥ Q ?", "dec", 3.2, 1.0),
        "d": (11.0, 8.0, "return a", "fine", 2.2, 0.75),
        "e": (5.5, 6.6, "b = ultima quota della tabella", "proc"),
        "f": (5.5, 5.3, "m = (a + b)/2;  READ_VALUES(m)\nQ_u(m) = A√(g·jden·|S|)", "proc", 5.6, 1.0),
        "g": (5.5, 3.8, "Q_u(m) < Q ?", "dec", 3.2, 1.0),
        "h": (2.4, 2.4, "a = m", "proc", 2.2, 0.8),
        "i": (8.6, 2.4, "b = m", "proc", 2.2, 0.8),
        "l": (11.0, 5.3, "dopo 60 iterazioni:\nreturn (a + b)/2", "fine", 3.4, 1.0)},
        archi=[("a", "b"), ("b", "c"), ("c", "d", "sì"), ("c", "e", "no"), ("e", "f"), ("f", "g"),
               ("g.w", "h.n", "sì", [(2.4, 3.8)]), ("g.e", "i.n", "no", [(8.6, 3.8)]),
               ("h.w", "f.w", None, [(0.6, 2.4), (0.6, 5.3)]), ("i.s", "f.w", None, [(8.6, 1.4), (0.3, 1.4), (0.3, 5.3)]),
               ("f.e", "l.w", "fine")]),
        didascalia="Flowchart di `UNIFORM_LEVEL`.", larghezza=0.95)
    D.estratto(SUB, "function UNIFORM_LEVEL", n=11, titolo=SUB)
    def graf(fig):
        ax = fig.add_subplot(111)
        S = abs(G[k - 1, 1])
        Q = r[:, 2] * np.sqrt(9.81 * r[:, 6] * S)
        ax.plot(Q, r[:, 1], "-o", ms=2.5, color="#2a78d6", lw=1.3, label="Q_u(η) dalla tabella")
        Qt = 5.07; e = np.interp(Qt, Q, r[:, 1])
        ax.axvline(Qt, color="#eb6834", lw=1, ls="--"); ax.axhline(e, color="#eb6834", lw=1, ls="--")
        ax.plot(Qt, e, "o", color="#eb6834", label=f"Q = {Qt} m³/s → η ≈ {e:.2f} m")
        ax.set_xlabel("portata di moto uniforme [m³/s]", fontsize=8); ax.set_ylabel("livello η [m s.l.m.]", fontsize=8)
        ax.set_title(f"S1000, pendenza locale |S| = {S:.5f}", fontsize=9); ax.grid(alpha=0.25); ax.tick_params(labelsize=7)
        ax.legend(fontsize=7, loc="lower right")
    D.grafico(graf, "moto_uniforme", "Scala di deflusso di moto uniforme in S1000 costruita con la tabella: lo stato "
              "iniziale con la portata di base del `fittizio_1` (5,07 m³/s) è il punto arancione.")
    D.attenzione("la bisezione funziona se Q_{u}(η) cresce con il livello. Con il K_{s} equivalente variabile delle "
                 "tabelle non è garantito in ogni sezione, ma sul Test_Brenta non ha dato problemi. Nelle sezioni in "
                 "contropendenza si usa |S|, quindi la stima può essere lontana dall'equilibrio: il riempimento la corregge.")
    D.h1("Il codice nel solutore")
    D.estratto(SEMI, "    Q_aff = zeros(imax)  # Vettore portate affluenti", n=10, titolo=SEMI)
    D.p("Poi `BARYCENTERS!` legge dalle tabelle area, β, j_{den} e R_{h} ai livelli iniziali, e il ciclo temporale parte "
        "in modalità riempimento.")
    D.h2("Perché non il velo costante e non la rampa")
    D.p("In passato si partiva da un velo d'acqua costante e si faceva crescere la portata con una rampa di 6 ore. "
        "Partire dalle portate cumulate con il moto uniforme porta all'equilibrio in circa 1,7 h invece di 7,4 h "
        "nel Test_Brenta, con lo stesso stato finale. La storia delle prove è in `CLAUDE.md`.")
    D.collegamenti([(4, "04_tabelle_delle_sezioni.pdf", "Le tabelle"), (9, "09_parte_esplicita_rusanov.pdf", "La condizione di valle"),
                    (18, "18_fine_riempimento.pdf", "La fine del riempimento")])


PDF = [
    (1, "01_equazioni_e_schema.pdf", "Equazioni e schema numerico", "Le equazioni di de Saint-Venant, l'idea dello schema semi-implicito e il flowchart del programma", "§2 e §5", p01),
    (2, "02_lettura_settings.pdf", "Lettura di settings.txt", "`INITIAL_SETTINGS!`, `PARSE_*`, `FILE_SIZE` e la struttura `Parametri`", "§6", p02),
    (3, "03_geometria_e_griglia.pdf", "Geometria e griglia sfalsata", "`LOAD_GEOMETRY` e `SPATIAL_DOMAIN!`", "§7.1-7.2", p03),
    (4, "04_tabelle_delle_sezioni.pdf", "Le tabelle delle sezioni", "Formato, interpolazione ed estrapolazione: `READ_VALUES`, `FIND_ROW`, `SECTION_INDEX`", "§7.3-7.4", p04),
    (5, "05_idrogrammi.pdf", "Gli idrogrammi", "Lettura e interpolazione di Lagrange: `LOAD_HYDROGRAPH!`, `LAGRANGE_INTERPOLATION!`", "§8.1", p05),
    (6, "06_affluenti.pdf", "Gli affluenti", "Sorgenti di volume nella continuità: `TRIBUTARY_INFLOW!`", "§8.2 e §10.6", p06),
    (7, "07_stato_iniziale.pdf", "Lo stato iniziale", "Portate cumulate e moto uniforme: `UNIFORM_DISCHARGE`, `UNIFORM_LEVEL`", "§9", p07),
]
