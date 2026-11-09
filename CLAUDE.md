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
- Nel repository i file sono salvati con fine riga **LF** (su Windows possono comparire CRLF nella
  copia locale): mantieni quello che trovi nel file.
- **`CLAUDE.md` va tenuto aggiornato** a mano a mano che si lavora.
- Per i grafici si usa `Plots.jl`, con testi in italiano e unità di misura.

## Come si usa

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"     # la prima volta
julia --project=. -O2 main.jl                          # condizioni iniziali + simulazione
julia --project=. strumenti/controllo_input.jl         # controllo a secco degli input, nessuna simulazione
julia --project=. strumenti/grafici_simulazione.jl     # grafici in <output>/grafici
julia --project=. strumenti/grafici_spinup.jl          # figura e animazione del riempimento in <output>/grafici
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
| `reference/` | Risultati del modello 2D usati come riferimento (`ns_channel_piena_discharge.csv`) |
| `calcolo_tabelle_daBM/` | Script `src/main.jl` che costruisce le tabelle 1D dai risultati 2D (BASEMENT) a scalini di portata |
| `strumenti/` | Script di controllo, grafici e confronto |
| `idromaker/main.jl` | Crea gli idrogrammi degli scenari fittizi in `input/idrogrammi/fittizio_*/piena.csv` (`Centa.csv` è solo un termine di paragone) |

## settings.txt

Letto **per posizione** (26 righe, non spostarle; dalla riga 27 in poi sono note e vecchie configurazioni):

| Riga | Contenuto |
|---|---|
| 1 | Tabelle delle sezioni (es. `input/geometria/Test_Brenta/tabella_sezioni.csv`) |
| 2 | Geometria (es. `input/geometria/Test_Brenta/localslope_xcoord.csv`) |
| 3-13 | Idrogrammi affluenti Brenta1-6, BrentaNuovo, BrentaVecchio1-4 |
| 14 | Idrogramma in ingresso da monte (chiave interna `caldonazzo`) |
| 15-21 | Idrogrammi affluenti Ceggio, Levico, Centa, Chieppena, Grigno, Maso, Moggio |
| 22 | Cartella degli output (creata se non esiste) |
| 23 | Ks |
| 24 | Tempo iniziale |
| 25 | Passo temporale massimo [s] |
| 26 | Numero massimo di passi |

- Numero di righe e colonne delle tabelle e numero di sezioni **sono ricavati dai file**
  (`FILE_SIZE`): non serve più il vecchio `struttura_output.txt`.
- Tutti gli idrogrammi indicati devono esistere. `Brenta7.csv`, non ancora usato, non è più in `settings.txt`.
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
locale è negativa verso valle (nel codice si usa `abs`). Dopo queste possono esserci
`argine_dx_[m], argine_sx_[m]`, quote assolute degli argini (in `calcolo_tabelle_daBM` vengono da
`input/bankpoints.csv`). I file a 4 colonne restano validi: gli argini valgono `Inf` e non c'è sfioro.

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

## Esondazione (portata uscente)

- Sfioro libero sugli argini, formula 16 di `reference/siviglia_stocchino_colombini_ASCE2009.pdf`:
  `q_l = Cd·√(2g·Y₀³)` [m²/s], con `Cd = 0.385` e `Y₀` tirante sopra la cresta (`FREE_OVERFLOW`).
- Un contributo per ogni argine superato: `Q_lat = (q_dx + q_sx)·dx_int`, cioè lo sfioro interessa
  tutta la lunghezza della cella con lo stesso `Y₀`.
- Entra solo nella continuità, come gli affluenti ma con segno opposto (`rhs -= dt/dx * Q_lat`), ed è
  esplicito (`eta` del passo precedente). È attivo anche nel riempimento e compare nel suo bilancio.
- L'acqua uscita non rientra: non ci sono golene né sfioro rigurgitato.
- `Q_lat` non viene salvato: si ricava da `eta.csv` con la stessa formula.

## Stato attuale e risultati

**Test_Brenta** (configurazione corrente di `settings.txt`):
- tratto di 6 km, 30 sezioni ogni 200 m (S0200-S6000);
- tabelle da `calcolo_tabelle_daBM`;
- idrogramma `input/idrogrammi/prova_tratto_brenta/piena.csv` (5,07 → 83,34 m³/s a 7 h, durata 24 h);
- Ks 33, dt massimo 60 s, upwind, nessun affluente;
- output in `output/tratto_brenta_sfioro/`, con lo sfioro sugli argini (configurazione corrente);
- `output/tratto_brenta/` conserva i risultati **senza sfioro**: non va sovrascritta. Per rifarli
  serve una geometria a 4 colonne e un'altra cartella alla riga 22 di `settings.txt`.

**Con lo sfioro sugli argini** (`output/tratto_brenta_sfioro/`; i risultati nell'elenco più sotto
sono invece **senza sfioro**):
- colmo in S6000 da 83,15 a 70,49 m³/s (−15%), a 7,47 h invece che a 7,55 h: il ritardo quasi non cambia;
- esce il 4,4% del volume entrato (circa 119 000 m³), quasi tutto da S0800 e S1000 sull'argine destro;
- portata uscente massima 12,8 m³/s; sfiorano solo S0600, S0800, S1000, S1200 e S2000;
- bilancio di massa chiuso e nessuna oscillazione con dt = 60 s;
- il livello supera ancora la quota massima delle tabelle in 5 sezioni (3, 6, 7, 10, 11), al massimo di 2,3 cm;
- grafico del confronto con e senza sfioro: `output/tratto_brenta/grafici/10_effetto_sfioro.png`
  (fatto una tantum sulla prova temporanea: `strumenti/grafici_simulazione.jl` non lo rigenera);
- animazioni della piena 2D in `output/tratto_brenta_sfioro/grafici/`: `10_piena_2D.gif` (solo 2D) e
  `11_piena_1D_2D.gif` (1D con sfioro e 2D insieme, con la differenza di quota). Fatte una tantum:
  livelli e portate 2D vengono da `../Test_Brenta/results/inputs/ns_channel_piena_results.h5`
  (30 sezioni, un'uscita ogni 900 s), perché `calcolo_tabelle_daBM/input/sections_ns_channel_S0200_S6000/`
  contiene gli **scalini di portata** per le tabelle, non la piena. Differenza di quota 1D − 2D:
  scarto quadratico medio 6 cm, massimo 27 cm (x = 5,8 km, t = 11,25 h, 1D più basso);
- evoluzione del riempimento in `output/tratto_brenta_sfioro/grafici/`: `12_spinup_convergenza.png` e
  `13_spinup_animazione.gif` (dagli `IC_*.csv`, con `strumenti/grafici_spinup.jl`). I criteri di arresto su portata e
  bilancio sono soddisfatti da 2,95 h; il riempimento dura 6 h solo per `T_riempimento`.

Risultati senza sfioro:
- **Riempimento:** finisce a 6 h, senza NaN.
- **Il modello 1D quasi non lamina:** colmo da 83,34 a 83,15 m³/s, con un ritardo di circa 33 min sui 6 km.
- **Confronto con il 2D** (`reference/ns_channel_piena_discharge.csv`, sezione S6000), con il
  colmo 2D di 66,28 m³/s a 9,00 h:

  | Modello 1D | Colmo [m³/s] | Istante [h] | Differenza colmo | Volume | RMSE [m³/s] | NSE |
  |---|---|---|---|---|---|---|
  | Senza sfioro | 83,15 | 7,53 | +25,5% | +5,7% | 6,05 | 0,933 |
  | Con sfioro | 70,49 | 7,45 | +6,4% | +1,1% | 2,98 | 0,984 |

  Grafici: `output/tratto_brenta/grafici/` senza sfioro, `output/tratto_brenta_sfioro/grafici/` con sfioro.
- **Nel 2D** il colmo cala soprattutto nei primi 2,4 km (da 80,6 a 66,3 m³/s), poi resta costante, con
  un ritardo di 120 min dalla prima all'ultima sezione (28 min nel 1D con sfioro). L'idrogramma in
  S6000 ha un plateau intorno a 66 m³/s fino a circa 10,5 h. In S6000 transita il 94,5% del volume
  entrato, più del 93,2% di S3000: nei tratti a valle l'acqua sembra rientrare nel canale. Il 1D non
  ha questo rientro, e dopo le 10 h la sua portata è più bassa del 2D fino a 8 m³/s.
- **Nel 1D** al colmo il livello supera la quota massima delle tabelle nelle sezioni 2-17, di 3-24 cm:
  lì i valori sono estrapolati e non c'è nessuna golena.
- Nel riferimento attuale le colonne `Q_P17_*` hanno lo stesso colmo di S3000 (66,29 m³/s a 8,75 h).

**Canale rettangolare di prova** (`input/geometria/rettangolare_L30km_dx200`, 30 km, B = 10 m,
pendenza 0,003, Ks 33, idrogrammi `vaia` con 18 affluenti), risultati in `output/rettangolare/`:
- **Confronto dx** 20 / 100 / 200 m: picchi entro lo 0,05%; dx 100-200 m bastano. Tempi di calcolo:
  dx 20 circa 390 s, dx 100 circa 30 s, dx 200 circa 13 s.
- **Magre e tabelle:** con portate molto basse le tabelle con larghezza e β nulli sul fondo (sotto
  10 cm) fanno fallire Newton e non conservano la massa. Servono tabelle fisiche e fitte vicino
  al fondo. Resoconto completo in `output/rettangolare/resoconto_magra/resoconto.md`.

## Scenari fittizi

Idrogrammi di monte creati da `idromaker/main.jl` (lanciare `julia --project=. idromaker/main.jl`
dalla cartella del progetto): Gumbel con portata di base 5 m³/s e passo 900 s. Valori provati su una
copia temporanea, 1D con sfioro:

| Scenario | Caratteristiche | Durata | Colmo [m³/s] | Colmo a S6000 | Ore con sfioro | Volume sfiorato |
|---|---|---|---|---|---|---|
| `fittizio_1` | Riferimento (uguale a `prova_tratto_brenta`) | 24 h | 83,34 a 7 h | 70,49 | 4,3 | 118 940 m³ (4,4%) |
| `fittizio_2` | Come l'1, portata più bassa | 24 h | 74,83 a 7 h | 69,50 | 3,2 | 37 091 m³ (1,5%) |
| `fittizio_3` | Come l'1, colmo mantenuto per 6 h | 30 h | 83,34 da 7 a 13 h | 70,50 | 10,3 | 396 194 m³ (8,8%) |
| `fittizio_4` | Due colmi come l'1, a 7 e 22 h | 39 h | 83,34 | 70,49 | 8,7 | 239 275 m³ (4,6%) |

Nel 1D l'esondazione inizia intorno a 66 m³/s (S1000). Il colmo a valle cambia poco tra gli scenari
perché lo sfioro limita la portata a circa 70 m³/s.

## Problemi aperti e prossimi passi possibili

- **Esondazione, il punto principale per avvicinare il 1D al riferimento.** Includere le golene nelle
  tabelle sarebbe la strategia migliore, ma non è banale, soprattutto dove l'acqua in golena arriva da
  monte invece di uscire in quella sezione. Si procede quindi, un passo alla volta, con una **portata
  laterale uscente** (vedi "Esondazione"). I dati 2D di riferimento non sono ancora definitivi: per
  ora interessa l'effetto del termine, non far combaciare gli idrogrammi.
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

## Considerazione finale per lunedì
I risultati sono incoraggianti, il modello 1D si è avvicinato molto a quello che riusciamo a misurare con il 2D. Non sono convinto che le le informazioni estratte dal 2D siano effettivamente utilizzate in maniera appropriata nell'1D, per lo meno però possiamo fare un confronto.  
Ci sono ulteriori aspetti che non abbiamo ancora considerato:
- Unire gli affluenti allo sfioro;
- Nel 2D, considerando tutto il fiume, si vede chiaramente che l'acqua tende a rientrare dopo un certo punto. Attualmente non abbiamo questa dinamica, perchè l'acqua che facciamo uscire poi non rientra più. La S6000 mostra un copmportamento che fa pensare che l'acqua sia rientrata, perchè non solo il picco si è ridotto, ma è anche in ritardo. Dobbiamo cercare di capire come includere questo aspetto in maniera intelligiente.