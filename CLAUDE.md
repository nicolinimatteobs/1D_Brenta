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
- **`CLAUDE.md` descrive lo stato committato del progetto** e va aggiornato insieme a ogni commit che
  cambia codice, input, configurazione o risultati. Il lavoro in corso non committato non va descritto.
  In fondo, "Ultimo commit descritto" indica fino a dove arriva la descrizione: per aggiornarla si
  parte da `git log <quel commit>..HEAD`.
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
Gli script dei grafici accettano argomenti facoltativi, per lavorare su una cartella diversa da quella
della riga 22 di `settings.txt` senza modificarlo:

```bash
julia --project=. strumenti/grafici_simulazione.jl output/fittizio_1 input/idrogrammi/fittizio_1/piena.csv
julia --project=. strumenti/grafici_spinup.jl      output/fittizio_1
julia --project=. strumenti/confronto_reference.jl output/fittizio_1 reference/ns_50_fittizio_1_discharge.csv
```

`confronto_reference.jl` confronta solo il periodo comune tra simulazione e riferimento e, sezione per
sezione, scarta gli istanti in cui il riferimento ha portata nulla (2D che parte da alveo asciutto).
Se il riferimento non ha ancora raggiunto il colmo quando finisce, il colmo non viene confrontato.
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
| `input/geometria/<caso>/` | Tabelle delle sezioni e geometria; per il Test_Brenta due spaziature, `Test_Brenta/S_200m/` e `Test_Brenta/S_50m/` |
| `input/idrogrammi/<evento>/` | Idrogrammi `t [s], Q [m3/s]` |
| `output/<caso>/` | Risultati (`IC_*.csv`, `*.csv`) e grafici |
| `reference/` | Risultati del modello 2D usati come riferimento (`ns_channel_piena_discharge.csv`) |
| `calcolo_tabelle_daBM/` | Script `src/main.jl` che costruisce le tabelle 1D dai risultati 2D (BASEMENT) a scalini di portata (vedi "Costruzione delle tabelle") |
| `strumenti/` | Script di controllo, grafici e confronto |
| `idromaker/main.jl` | Crea gli idrogrammi degli scenari fittizi in `input/idrogrammi/fittizio_*/piena.csv` (`Centa.csv` è solo un termine di paragone) |

## settings.txt

Letto **per posizione** (27 righe, non spostarle; dopo ci sono la legenda e vecchie configurazioni):

| Riga | Contenuto |
|---|---|
| 1 | Tabelle delle sezioni (es. `input/geometria/Test_Brenta/S_200m/tabella_sezioni.csv`) |
| 2 | Geometria (es. `input/geometria/Test_Brenta/S_200m/localslope_xcoord.csv`) |
| 3-13 | Idrogrammi affluenti Brenta1-6, BrentaNuovo, BrentaVecchio1-4 |
| 14 | Idrogramma in ingresso da monte (chiave interna `caldonazzo`) |
| 15-21 | Idrogrammi affluenti Ceggio, Levico, Centa, Chieppena, Grigno, Maso, Moggio |
| 22 | Cartella degli output (creata se non esiste) |
| 23 | Ks (solo condizione al contorno di valle) |
| 24 | Tempo iniziale |
| 25 | Passo temporale massimo [s] |
| 26 | Numero massimo di passi |
| 27 | Affluenti attivi: `si` o `no` (accetta anche `sì`, `true`/`false`, `1`/`0`) |

- Numero di righe e colonne delle tabelle e numero di sezioni **sono ricavati dai file**
  (`FILE_SIZE`): non serve più il vecchio `struttura_output.txt`.
- Tutti gli idrogrammi indicati devono esistere, anche con gli affluenti spenti. `Brenta7.csv`, non
  ancora usato, non è più in `settings.txt`.
- La durata simulata è il tempo finale **più breve** tra gli idrogrammi letti.
- Il `settings.txt` committato punta al canale rettangolare con l'idrogramma `fittizio_2`, output in
  `output/fittizio_2` e affluenti spenti: controllare righe 1-2, 14, 22 e 27 prima di lanciare.

## Formato degli input

**Tabelle delle sezioni** (CSV con intestazione, 12 colonne, righe raggruppate per sezione e con
quote crescenti): `sezione, WSE, Ω, Ω_s, β, β_s, jden, jden_s, W, W_s, RH, RH_s`.
- La colonna 2 è la **quota del pelo libero** (non il tirante). Il passo tra le quote può essere non costante.
- Le colonne `_s` sono le pendenze rispetto alla quota, usate per l'interpolazione lineare.
- `jden = Ks²·RH^(4/3)/g` è il denominatore del termine d'attrito (γ = |u|/jden). La scabrezza usata
  in moto vario è quindi **quella contenuta nelle tabelle**: in quelle da `calcolo_tabelle_daBM` è
  il Ks equivalente ricavato dal 2D, livello per livello. Anche lo stato iniziale usa l'attrito delle
  tabelle (`UNIFORM_DISCHARGE`, `Q = A·√(g·jden·S)`). Il Ks di `settings.txt` serve ormai solo per
  la condizione di moto uniforme all'uscita, che quindi non è coerente con le tabelle (nell'ultima
  sezione del Test_Brenta la tabella dà circa 18 contro 33).

**Geometria**: `section, local_slope, x_coord_[m], talweg_[m]`, una riga per sezione. La pendenza
locale è negativa verso valle (nel codice si usa `abs`). Dopo queste possono esserci
`argine_dx_[m], argine_sx_[m]`, quote assolute degli argini (in `calcolo_tabelle_daBM` vengono da
`bankpoints.csv` nella cartella di input del caso). I file a 4 colonne restano validi: gli argini valgono `Inf` e non c'è sfioro.

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

1. **Partenza dalle portate cumulate.** La portata alle interfacce è la prima portata dell'idrogramma
   di monte più gli affluenti entrati più a monte (solo se la riga 27 è `si`). In ogni sezione il
   tirante è quello di **moto uniforme** per la portata media della cella, con la pendenza locale e
   l'attrito delle tabelle (`UNIFORM_LEVEL`, bisezione con `UNIFORM_DISCHARGE`), e non sotto
   `y_ini = 0.0` m (quindi nessun tirante minimo). La vecchia partenza (velo costante) è commentata.
   - Prima si usava Manning con il Ks di `settings.txt` (33), mentre nelle tabelle del Test_Brenta
     il Ks ai tiranti iniziali va da 1,9 a 50 (mediana 16): nel fittizio_1 al primo passo la portata
     crollava fino a 0,19 m³/s. Con l'attrito delle tabelle resta tra 3,98 e 6,43 m³/s.
   - La bisezione suppone che la portata di moto uniforme cresca con la quota: con il Ks variabile
     delle tabelle non è garantito, ma sul Test_Brenta ha funzionato.
2. **Nessuna rampa:** portata di monte e affluenti sono al loro primo valore da subito. Passo massimo
   `dt_riempimento = 5` s.
3. **Arresto** quando valgono tutte insieme:
   - sono passati almeno `T_assestamento = 3600` s;
   - la portata non cambia più (`|q − qn| < 0,01` in tutte le celle);
   - il livello è fermo: ogni `T_fermo = 1800` s si confronta con quello di inizio finestra, e deve
     essere cambiato meno di `tol_fermo = 1 mm` in tutte le sezioni;
   - l'uscita è uguale agli ingressi (monte + affluenti − sfioro) entro lo 0,1%.

   Poi tempo e contatori ripartono da zero e inizia la simulazione vera.

In un tratto **non prismatico** il moto uniforme sezione per sezione non è l'equilibrio: durante il
riempimento il tratto si assesta (rigurgiti) e trattiene acqua. Anche a valle degli affluenti
la partenza cumulata non è l'equilibrio. Per partire ancora più vicini servirebbe un profilo di moto permanente.

**Perché così** (prove in copie temporanee, grafici in `output/IC/Test_Brenta_vaia/` e
`output/IC/rettangolare_vaia/`, con gli affluenti `vaia` e portata di monte `vaia/Caldonazzo.csv`):
- Prima c'era una rampa di 6 h in cui gli affluenti salivano da 0 al primo valore. Inserirli uno alla
  volta allunga il riempimento di 3-4 volte e con dt 60 s va in NaN o non converge.
- Con le portate cumulate lo stato finale è identico a quello con la rampa. Il livello è entro 1 mm
  dall'equilibrio a 1,73 h invece di 7,35 h (Test_Brenta) e a 4,67 h invece di 10,39 h (rettangolare):
  circa 5,6-5,7 h in meno (`7_assestamento.png`, `8_rampa_e_assestamento.gif`).
- Il solo criterio su portata e bilancio scattava troppo presto: nel rettangolare a 2,76 h, con il
  livello ancora a 6 mm dall'equilibrio. Da qui la durata minima e il controllo sul livello.
- Con il codice attuale il riempimento finisce a 2,5 h nel Test_Brenta con affluenti e nel fittizio_1,
  a 5,5 h nel rettangolare con affluenti, con lo stato finale entro 0,03 mm dall'equilibrio. Le
  simulazioni successive coincidono con quelle fatte con la rampa (scarti sotto 1 mm e 0,03 m³/s).
- Nel Test_Brenta con affluenti al primo passo l'uscita arriva a 50 m³/s (contro 36): la condizione
  al contorno di valle usa il Ks 33 di `settings.txt`, la tabella nell'ultima sezione circa 18.

## Affluenti

- 18 affluenti distribuiti uniformemente: `sez_aff = round(k*imax/19)`.
- Quelli con picco ≥ 10 m³/s entrano su 3 sezioni (`semi_aff`, `ncel_aff`), gli altri su una.
- Entrano solo come sorgente nella continuità (`rhs += dt/dx * Q_aff`), senza quantità di moto.
  `TRIBUTARY_INFLOW!` (`subroutines.jl`) distribuisce le portate dei 18 affluenti sulle celle; in
  simulazione si interpolano con `LAGRANGE_INTERPOLATION!` al tempo `time + dt`.
- Si accendono e spengono con la **riga 27 di `settings.txt`** (`con_affluenti` in `Parametri`), che
  vale per stato iniziale, riempimento e simulazione. Non serve toccare il codice.

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
- tratto di 6 km, 30 sezioni ogni 200 m (S0200-S6000), in `input/geometria/Test_Brenta/S_200m/`;
  tutti i risultati qui sotto sono con questa geometria (quella a 50 m non è ancora stata simulata);
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
  (30 sezioni, un'uscita ogni 900 s), perché `calcolo_tabelle_daBM/input/S_0200_6000/`
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

**Risultati salvati e confronto con il 2D.** Le simulazioni 1D con sfioro sono in `output/fittizio_1`
… `fittizio_4`, ciascuna con la sua cartella `grafici/`: grafici 1-7 e GIF di `grafici_simulazione.jl`,
8-9 e `confronto_S6000.txt` di `confronto_reference.jl`, 12-13 di `grafici_spinup.jl`. Mancano i
grafici 10-11 (piena 2D), che richiedono i livelli 2D: per i fittizi ci sono solo le portate.

I riferimenti 2D `reference/ns_50_fittizio_*_discharge.csv` hanno colonne ogni 50 m (S0050-S6000)
e sono diversi da `ns_channel_piena_discharge.csv`:
- il 2D **parte da alveo asciutto**: in S6000 la portata è nulla fino a 2,5 h, quando arriva il fronte,
  mentre il 1D parte con la portata di base;
- tutti e quattro coprono ormai tutta la piena (`fittizio_1` è rimasto quello del commit `7ba5b92`;
  rilanciato con il codice di quel commit, il 1D dà gli stessi file byte per byte).

| Scenario | Riferimento fino a | Colmo in S6000, 2D → 1D | RMSE [m³/s] | NSE | Volume 1D − 2D |
|---|---|---|---|---|---|
| `fittizio_1` | 24 h (completo) | 68,14 a 9,75 h → 70,49 a 7,45 h (+3,5%, 138 min prima) | 3,22 | 0,981 | +0,7% |
| `fittizio_2` | 24 h (completo) | 64,16 a 8,00 h → 69,50 a 7,48 h (+8,3%, 31 min prima) | 2,30 | 0,988 | +2,9% |
| `fittizio_3` | 30 h (completo) | 81,71 a 13,50 h → 70,50 a 8,68 h (−13,7%, 289 min prima) | 6,14 | 0,954 | −5,3% |
| `fittizio_4` | 39 h (completo) | 1ª onda 68,14 a 9,75 h → 70,49 a 7,45 h (+3,4%, 138 min prima); 2ª onda 78,54 a 24,75 h → 70,49 a 22,40 h (−10,2%, 141 min prima) | 4,32 | 0,967 | −1,9% |

Nel `fittizio_1` il ritardo del colmo dalla prima all'ultima sezione è di 165 min nel 2D e 28 nel 1D.
Nel `fittizio_4` la prima onda del 2D è uguale al `fittizio_1` (68,14 m³/s a 9,75 h in S6000). Nella
seconda, a tratto già bagnato, il colmo 2D cala fino a 66-67 m³/s nei primi 3,4 km e poi risale fino a
78-80 m³/s tra 3,6 e 4,2 km: è lo stesso rientro del `fittizio_3`. Il 1D dà 70,49 m³/s in entrambe le
onde. Il ritardo del colmo dalla prima all'ultima sezione è di 165 min nel 2D e 27-28 nel 1D, per
entrambe le onde. **Attenzione:** `confronto_reference.jl` cerca un solo colmo su tutta la simulazione,
quindi per il `fittizio_4` il grafico 9 e il "ritardo di 1065 min" di `confronto_S6000.txt` mescolano le
due onde (a monte il massimo è nella prima, a valle nella seconda). I valori per onda qui sopra sono
calcolati a parte, dividendo a 14,5 h.

Nel `fittizio_2` il colmo 2D cala da 74,3 a 64,2 m³/s nei primi 2,4 km e poi resta costante, senza
risalite: con questa portata non si vede il rientro del `fittizio_3`. Il 1D lamina meno (69,5 m³/s da
1,2 km in poi) e il ritardo del colmo dalla prima all'ultima sezione è di 60 min nel 2D e 29 nel 1D.
Anche qui il confronto è stato rifatto senza rilanciare la simulazione (vedi sotto).

Nel `fittizio_3` il 2D in S6000 resta a 62-66 m³/s fino a circa 9,5 h (sotto il 1D), poi sale fino a
81,7 m³/s, quasi tutta la portata in ingresso (83,3), e la mantiene fino a circa 15,5 h. Il colmo 2D
lungo il tratto risale da 67 a 82 m³/s tra 3,4 e 4,2 km. Con il picco lungo le aree esondate sembrano
riempirsi e l'acqua rientra nel canale a valle: il 1D non lo riproduce, perché l'acqua sfiorata non
rientra. Il ritardo del colmo dalla prima all'ultima sezione è di 360 min nel 2D e 56 nel 1D.
Il confronto è stato rifatto senza rilanciare la simulazione: con il codice attuale (affluenti attivi,
`y_ini = 0`) e `settings.txt` sul canale rettangolare, un nuovo run di `fittizio_3` non sarebbe
confrontabile con il 2D. Anche `confronto_reference.jl` legge la geometria dalla riga 2 di
`settings.txt`: con la configurazione attuale va lanciato da una copia con la geometria del Test_Brenta.

## Costruzione delle tabelle (`calcolo_tabelle_daBM`)

`julia --project=. calcolo_tabelle_daBM/src/main.jl` dalla cartella del progetto. Il caso si sceglie con
le costanti in testa al file (`InputFolder`, `OutputFolder`, `FilePrefix`):

| Caso | Input | Prefisso dei file | Output | Sezioni |
|---|---|---|---|---|
| 50 m (attivo) | `input/S_0050_6000/` | `ns_50_scalini_` | `output/S_50m/` | 120, S0050-S6000 |
| 200 m | `input/S_0200_6000/` | `ns_channel_` | `output/S_200m/` | 30, S0200-S6000 |

- Ogni cartella di input contiene i file delle sezioni, `sections_dx.csv` (nome e ascissa) e
  `bankpoints.csv` (quote degli argini). Le sezioni BC e P17 non sono usate.
- I file a 50 m non hanno la colonna `t [s]`: il controllo che le sezioni abbiano gli stessi scalini
  usa il numero di righe, più il tempo quando c'è.
- In output: una cartella per sezione (`tabella.csv` e grafici, riscritti solo se la tabella cambia),
  `tabella_sezioni.csv` e `localslope_xcoord.csv`, da copiare in `input/geometria/Test_Brenta/S_<passo>m/`.
- **Correzioni automatiche**, segnalate tra gli avvisi a fine esecuzione (a 200 m non intervengono e
  i risultati sono identici a prima):
  - una pendenza del pelo libero nulla o negativa viene sostituita con quella del livello bagnato
    successivo della stessa sezione (a 50 m solo S5950 a 0,86 m³/s: −3,2·10⁻⁵ → 1,2·10⁻⁴); altrimenti
    Ks e `jden` sarebbero NaN e la ricostruzione delle righe asciutte fallirebbe;
  - le righe con livello non superiore al massimo dei precedenti vengono scartate, perché il modello
    vuole quote crescenti (a 50 m le ultime due di S0050, a 94,8 e 99,9 m³/s, dove il livello scende).

**Geometria a 50 m** (`input/geometria/Test_Brenta/S_50m/`, 22 livelli per sezione, 20 in S0050):
- nessun NaN e valori positivi; Ks mediano 20,4 (20,9 a 200 m), con la stessa proporzione di valori
  fuori da 10-50; restano Ks di 3-9 in alcune sezioni tra 0,6 e 1,25 km a livelli intermedi e 139
  nell'ultima riga di S0100;
- **pendenze del fondo molto più irregolari**: 23 sezioni in contropendenza e 6 con |pendenza| sotto
  2·10⁻⁴, fino a 7·10⁻⁶ in S1800 (l'ultima sezione usa sempre la pendenza media, −0,00213). Lo stato
  iniziale usa il moto uniforme con `|local_slope|`: a 5,07 m³/s il tirante supera 1,5 m in 7 sezioni
  (massimo 1,32 m a 200 m) e in S1800 il livello andrebbe oltre la tabella. Prima di simulare conviene
  valutare una pendenza più regolare per lo stato iniziale.

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
  adatta alle portate molto basse; `T_assestamento`, `T_fermo` e `tol_fermo` sono provati solo su
  Test_Brenta (6 km) e rettangolare (30 km).
- **Condizione al contorno di valle** con il Ks di `settings.txt`, non coerente con le tabelle.
- **Varie:**
  - Newton fa sempre 100 iterazioni (tolleranza 1e-14);
  - `THOMAS_ALGORITHM!` alloca due vettori a ogni chiamata;
  - `strumenti/controllo_input.jl` calcola il tirante iniziale con la sola portata di monte e un
    minimo di 10 cm, non con la partenza cumulata.

## Considerazione finale per lunedì
I risultati sono incoraggianti, il modello 1D si è avvicinato molto a quello che riusciamo a misurare con il 2D. Non sono convinto che le le informazioni estratte dal 2D siano effettivamente utilizzate in maniera appropriata nell'1D, per lo meno però possiamo fare un confronto.  
Ci sono ulteriori aspetti che non abbiamo ancora considerato:
- Unire gli affluenti allo sfioro;
- Nel 2D, considerando tutto il fiume, si vede chiaramente che l'acqua tende a rientrare dopo un certo punto. Attualmente non abbiamo questa dinamica, perchè l'acqua che facciamo uscire poi non rientra più. La S6000 mostra un copmportamento che fa pensare che l'acqua sia rientrata, perchè non solo il picco si è ridotto, ma è anche in ritardo. Dobbiamo cercare di capire come includere questo aspetto in maniera intelligiente.

## Ultimo commit descritto

`21b8f7a` ("gestione affluenti e IC cumulata su prima portata di idrogramma"): partenza dalle portate
cumulate senza rampa, nuovo criterio di arresto del riempimento, affluenti dalla riga 27 di
`settings.txt`, indentazione sistemata. Comprende `eac5cbf` (grafici delle prove in `output/IC/`) e
`3b3020b`. Strumenti utili per i riempimenti: `strumenti/confronto_spinup.jl`
(`<cartella dei grafici> <cartella 1> <cartella 2> ...`) e `strumenti/confronto_simulazioni.jl`
(stessi argomenti, per le simulazioni).
