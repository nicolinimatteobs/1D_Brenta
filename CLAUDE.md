# 1D_Brenta: contesto del progetto

Modello idraulico 1D semi-implicito in Julia per la propagazione delle piene nel fiume Brenta,
sviluppato per una tesi di laurea magistrale (UNITN). Il codice, i commenti, i nomi e le risposte
sono in **italiano**. Repository: `github.com/nicolinimatteobs/1D_Brenta`, si lavora su `main`
o su branch personali che poi vengono uniti.

## Regole di lavoro

- **Non modificare il codice del modello senza il consenso dell'utente.** Prima di un cambiamento
  non banale spiega cosa vuoi fare e perché, poi aspetta l'ok. Per prove ed esperimenti usa
  **copie temporanee** del progetto (cartella di lavoro fuori dal repository), non i file in `Source/`.
- **Non far partire simulazioni del progetto senza un comando esplicito.** Le prove in copie
  temporanee vanno bene.
- **Commit e merge li fa l'utente**, a meno che non lo chieda espressamente.
- Se compare un NaN il modello si ferma da solo (`CHECK_NAN`): riportare fase, passo, tempo e posizioni.
- I file `.jl` e `settings.txt` usano fine riga **CRLF**: mantienila quando modifichi i file
  (attenzione a `sed -i`, che può convertirli in LF).
- Per i grafici si usa `Plots.jl`, con testi in italiano e unità di misura.

## Come si usa

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"     # la prima volta
julia --project=. -O2 main.jl                          # condizioni iniziali + simulazione
julia --project=. strumenti/controllo_input.jl         # controllo a secco degli input, nessuna simulazione
julia --project=. strumenti/grafici_simulazione.jl     # grafici in <output>/grafici
julia --project=. strumenti/confronto_reference.jl     # confronto con il 2D in reference/
```

Tutti i comandi vanno lanciati dalla cartella del progetto: i percorsi in `settings.txt` sono relativi.
Una simulazione del Test_Brenta (30 sezioni, 24 h) richiede circa 15 s.

## Struttura

| Percorso | Contenuto |
|---|---|
| `main.jl` | Avvia `RUN_SIMULATION.SEMI_IMPLICIT_BRENTA_2D_1D!()` |
| `Source/Run_Simulation.jl` | Modulo che include gli altri file |
| `Source/Semi_Implicit_Brenta_2D_1D.jl` | Il solutore: riempimento (condizioni iniziali) e simulazione |
| `Source/subroutines.jl` | Lettura tabelle e geometria, `READ_VALUES`, `UNIFORM_LEVEL`, idrogrammi, Thomas, `CHECK_NAN` |
| `Source/module.jl` | `INITIAL_SETTINGS!` (lettura di `settings.txt`), `FILE_SIZE`, `ALLOCATE_ARRAYS!`, `LOAD_DATA!` |
| `Source/struct.jl` | `struct Parametri` |
| `Source/writing.jl` | `INIT_OUTPUT` e funzioni di salvataggio |
| `settings.txt` | Configurazione della simulazione (vedi sotto) |
| `input/geometria/<caso>/` | Tabelle delle sezioni e geometria |
| `input/idrogrammi/<evento>/` | Idrogrammi `t [s], Q [m3/s]` |
| `output/<caso>/` | Risultati (`IC_*.csv`, `*.csv`) e grafici |
| `reference/` | Risultati del modello 2D usati come riferimento (`piena_ns_discharge.csv`) |
| `calcolo_tabelle_daBM/` | Script `src/main.jl` che costruisce le tabelle 1D dai risultati 2D (BASEMENT) a scalini di portata |
| `strumenti/` | Script di controllo, grafici e confronto |

## settings.txt

Letto **per posizione** (27 righe, non spostarle; dalla riga 28 in poi sono note e vecchie configurazioni):

| Riga | Contenuto |
|---|---|
| 1 | Tabelle delle sezioni (es. `input/geometria/Test_Brenta/tabella_sezioni.csv`) |
| 2 | Geometria (es. `input/geometria/Test_Brenta/localslope_xcoord.csv`) |
| 3-14 | Idrogrammi affluenti Brenta1-6, Brenta7, BrentaNuovo, BrentaVecchio1-4 |
| 15 | Idrogramma in ingresso da monte (chiave interna `caldonazzo`) |
| 16-22 | Idrogrammi affluenti Ceggio, Levico, Centa, Chieppena, Grigno, Maso, Moggio |
| 23 | Cartella degli output (creata se non esiste) |
| 24 | Ks |
| 25 | Tempo iniziale |
| 26 | Passo temporale massimo [s] |
| 27 | Numero massimo di passi |

- Numero di righe e colonne delle tabelle e numero di sezioni **sono ricavati dai file**
  (`FILE_SIZE`): non serve più il vecchio `struttura_output.txt`.
- Tutti gli idrogrammi devono esistere tranne `Brenta7.csv`, non ancora usato (viene saltato con un avviso).
- La durata simulata è il tempo finale **più breve** tra gli idrogrammi letti.

## Formato degli input

**Tabelle delle sezioni** (CSV con intestazione, 12 colonne, righe raggruppate per sezione e con
quote crescenti): `sezione, WSE, Ω, Ω_s, β, β_s, jden, jden_s, W, W_s, RH, RH_s`.
- La colonna 2 è la **quota del pelo libero** (non il tirante). Il passo tra le quote può essere non costante.
- Le colonne `_s` sono le pendenze rispetto alla quota, usate per l'interpolazione lineare.
- `jden = Ks²·RH^(4/3)/g` è il denominatore del termine d'attrito (γ = |u|/jden). La scabrezza usata
  in moto vario è quindi **quella contenuta nelle tabelle**: in quelle da `calcolo_tabelle_daBM` è
  il Ks equivalente ricavato dal 2D, livello per livello. Il Ks di `settings.txt` serve solo per lo
  stato iniziale (`UNIFORM_LEVEL`) e per la condizione di moto uniforme all'uscita.

**Geometria**: `section, local_slope, x_coord_[m], talweg_[m]`, una riga per sezione. La pendenza
locale è negativa verso valle (nel codice si usa `abs`).

**Output** (separati da spazi, una riga per passo temporale):
- `Q.csv`: portate alle **interfacce** (imax+1 valori: l'ultimo è l'uscita);
- `eta.csv`, `h.csv`: quota del pelo libero e tirante nelle **sezioni** (imax valori);
- `time.csv`: tempo [s];
- `IC_*`: gli stessi file per la fase di riempimento.

## Il metodo numerico

- Schema **semi-implicito**: avvezione esplicita (flusso tipo Rusanov), attrito implicito,
  equazione della quota risolta con Newton e algoritmo di Thomas (tridiagonale).
- **Lettura delle tabelle** (`READ_VALUES`): le tabelle sono trasposte (`data[colonna, riga]`, più
  veloce in Julia) e la riga si cerca con `FIND_ROW`, che parte dall'ultima riga usata e altrimenti
  usa la ricerca binaria. Interpolazione lineare con le colonne `_s`. Sotto la prima quota si usa la
  prima riga; **oltre l'ultima quota si estrapola con la pendenza dell'ultimo tratto** (calcolata
  dai valori delle ultime due righe).
- **Grandezze alle interfacce**: attualmente **upwind** per area e attrito (`A_int`, `den_γ_int`
  presi dalla cella a monte), media per β e Rh. La versione con la media è nel codice, commentata.
  - Con la media e celle lunghe rispetto a y/S (circa dx > 0,6·y/S) compare un **dente di sega** nei
    tratti poco profondi: è la soluzione stazionaria dello schema centrato, non un errore di input.
  - L'upwind lo elimina ma è meno accurato nei tratti profondi.
  - **Attenzione:** sul canale rettangolare di prova senza affluenti (velo di circa 10 cm su 30 km) l'upwind
    su area e attrito **insieme** è instabile e va in NaN. Upwind solo sull'area, solo sull'attrito
    o con la media è stabile. Sul Test_Brenta (tiranti di 0,4-3 m) funziona.

## Condizioni iniziali (riempimento)

1. In ogni sezione si parte dal **tirante di moto uniforme** (Manning, pendenza locale) per la prima
   portata dell'idrogramma di monte, calcolato da `UNIFORM_LEVEL` per bisezione sulla tabella,
   con un tirante minimo `y_ini = 0.10` m. La vecchia partenza (velo costante) è commentata.
2. Rampa di `T_riempimento = 21600` s con passo massimo `dt_riempimento = 5` s: la portata di monte
   va al primo valore dell'idrogramma (con la partenza da moto uniforme resta già costante) e gli
   affluenti salgono da 0 al loro primo valore (`r_IC`).
3. Il riempimento finisce quando la rampa è conclusa, la portata non cambia più (`|q − qn| < 0,01`
   in tutte le celle) e l'uscita è uguale agli ingressi entro lo 0,1%. Poi tempo e contatori
   ripartono da zero e inizia la simulazione vera.

In un tratto **non prismatico** il moto uniforme sezione per sezione non è l'equilibrio: durante il
riempimento il tratto si assesta (rigurgiti) e trattiene acqua. Sul Test_Brenta il tirante medio
passa da 0,60 a 0,85 m. Per partire più vicini all'equilibrio servirebbe un profilo di moto permanente.

## Affluenti

- 18 affluenti distribuiti uniformemente: `sez_aff = round(k*imax/19)`.
- Quelli con picco ≥ 10 m³/s entrano su 3 sezioni (`semi_aff`, `ncel_aff`), gli altri su una.
- Entrano solo come sorgente nella continuità (`rhs += dt/dx * Q_aff`), senza quantità di moto.
- **Ora sono commentati** (blocchi "INSERIMENTO AFFLUENTI" nel solutore) per il test sul tratto reale
  senza affluenti. Per riattivarli basta togliere i `#` dai due blocchi.

## Stato attuale e risultati

**Test_Brenta** (configurazione corrente di `settings.txt`):
- tratto di 6 km, 30 sezioni ogni 200 m (S0200-S6000);
- tabelle da `calcolo_tabelle_daBM`;
- idrogramma `input/idrogrammi/prova_tratto_brenta/piena.csv` (5,07 → 83,34 m³/s a 7 h, durata 24 h);
- Ks 33, dt massimo 60 s, upwind, nessun affluente;
- output in `output/tratto_brenta/`.

Risultati:
- **Riempimento:** finisce a 6 h, senza NaN.
- **Il modello 1D quasi non lamina:** colmo da 83,34 a 83,15 m³/s, con un ritardo di circa 33 min sui 6 km.
- **Confronto con il 2D** (`reference/piena_ns_discharge.csv`, sezione S6000):
  - colmo 2D: 73,09 m³/s a 10,75 h;
  - colmo 1D: 83,15 m³/s a 7,53 h (+14%, 3,2 h prima);
  - NSE 0,73.
- **Nel 2D** il colmo cala lungo il tratto e ha un ritardo di 225 min. L'idrogramma ha un plateau
  intorno a 50 m³/s tra 5,5 e 8 h, e circa il 5% del volume non arriva in S6000 entro 24 h.
  Il quadro è compatibile con un'esondazione nelle golene.
- **Nel 1D** al colmo il livello supera la quota massima delle tabelle nelle sezioni 2-17, di 5-24 cm:
  lì i valori sono estrapolati e non c'è nessuna golena.
- Le colonne `Q_P17_*` del riferimento (colmo 56,67 m³/s) non sono ancora state interpretate.

**Canale rettangolare di prova** (`input/geometria/rettangolare_L30km_dx200`, 30 km, B = 10 m,
pendenza 0,003, Ks 33, idrogrammi `vaia` con 18 affluenti), risultati in `output/rettangolare/`:
- **Confronto dx** 20 / 100 / 200 m: picchi entro lo 0,05%; dx 100-200 m bastano. Tempi di calcolo:
  dx 20 circa 390 s, dx 100 circa 30 s, dx 200 circa 13 s.
- **Magre e tabelle:** con portate molto basse le tabelle con larghezza e β nulli sul fondo (sotto
  10 cm) fanno fallire Newton e non conservano la massa. Servono tabelle fisiche e fitte vicino
  al fondo. Resoconto completo in `output/rettangolare/resoconto_magra/resoconto.md`.

## Problemi aperti e prossimi passi possibili

- **Tabelle più alte, con le golene**, per evitare l'estrapolazione e rappresentare la laminazione
  vista nel 2D: è il punto principale per avvicinare il 1D al riferimento.
- **Upwind:** l'instabilità con area e attrito insieme su veli sottili e tratti lunghi è da
  risolvere, oppure va scelto un upwind solo parziale.
- **Stato iniziale vicino all'equilibrio** nei tratti non prismatici: profilo di moto permanente,
  o almeno la pendenza tra i talweg al posto di `local_slope`.
- **Criterio di arresto del riempimento:** la tolleranza assoluta 0,01 m³/s su `|q − qn|` è poco
  adatta alle portate molto basse.
- **Varie:**
  - Newton fa sempre 100 iterazioni (tolleranza 1e-14);
  - `THOMAS_ALGORITHM!` alloca due vettori a ogni chiamata;
  - `Brenta7.csv` non esiste ancora.
