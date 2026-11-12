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
julia --project=. strumenti/grafici_simulazione.jl output/S_50m/fittizio_1 input/idrogrammi/fittizio_1/piena.csv
julia --project=. strumenti/grafici_spinup.jl      output/S_50m/fittizio_1
julia --project=. strumenti/confronto_reference.jl output/S_50m/fittizio_1 reference/ns_50_fittizio_1_discharge.csv
julia --project=. strumenti/animazioni_2D.jl       output/S_50m/fittizio_1 reference/ns_50_fittizio_1_wse.csv reference/ns_50_fittizio_1_discharge.csv
julia --project=. strumenti/confronto_valle.jl     output/S_50m/fittizio_1 reference/ns_50_fittizio_1_wse.csv reference/ns_50_fittizio_1_discharge.csv
julia --project=. strumenti/confronto_valle.jl     output/S_50m/fittizio_1 reference/ns_50_fittizio_1_wse.csv reference/ns_50_fittizio_1_discharge.csv S1000 S2150 S4000 S6000
```

Numerazione dei grafici in `<output>/grafici`:

| N. | File | Script |
|---|---|---|
| 1-6 | condizioni iniziali, idrogrammi, tiranti, profili, estrapolazione, mappe | `grafici_simulazione.jl` |
| 7-8 | `7_confronto_S6000.png`, `8_colmo_lungo_il_tratto.png` (e `confronto_S6000.txt`) | `confronto_reference.jl` |
| 9 | `9_piena_1D_2D.gif` | `animazioni_2D.jl` |
| 10-11 | `10_spinup_convergenza.png`, `11_spinup_animazione.gif` | `grafici_spinup.jl` |
| 12 | `12_confronto_valle_tirante_portata.png` (e `confronto_valle.txt`) | `confronto_valle.jl` |
| 13 | `13_confronto_sezioni_tirante_portata.png` (e `confronto_sezioni.txt`) | `confronto_valle.jl` con i nomi delle sezioni |

La GIF della piena del solo 1D (vecchio 7) e quella del solo 2D (vecchio 10) sono state tolte. Le cartelle
`output/tratto_brenta*`, `output/rettangolare` e `output/vaia` hanno ancora la numerazione vecchia. Il
grafico 13 e il terzo pannello del 12 ci sono solo in `output/S_50m/` (le altre cartelle hanno il 12 a
due pannelli).

`animazioni_2D.jl` (argomenti obbligatori) crea `9_piena_1D_2D.gif` (1D e 2D sovrapposti, con la
differenza di quota): il 1D è interpolato negli istanti del 2D (ogni 15 min) e le sezioni 2D ancora
asciutte non sono disegnate. Stampa media, scarto quadratico medio e massimo della differenza di quota,
ma il massimo cade di solito sul fronte 2D che arriva (prime 2,5 h) e non è significativo.

`confronto_valle.jl` (tre argomenti obbligatori, poi eventuali nomi di sezioni) confronta 1D e 2D:
- senza sezioni fa il grafico 12 nell'ultima sezione: tirante e portata nel tempo, più il tirante in
  funzione della portata con la scala di deflusso degli scalini 2D (da `calcolo_tabelle_daBM/input/S_0050_6000/`);
- con i nomi di una o più sezioni fa il grafico 13: una riga per sezione, tirante e portata nel tempo
  affiancati (in `output/S_50m` S1000, S2150, S4000 e S6000, cioè 8 pannelli).
- Il tirante 2D è il livello 2D meno il talweg della geometria 1D; gli istanti 2D asciutti sono esclusi.
- Nei file di testo c'è anche lo scarto **a parità di portata**: tirante meno quello degli scalini 2D
  alla stessa portata, per il 1D e per il 2D. È il confronto giusto per i livelli (vedi "Livelli 1D e 2D
  a parità di portata"): allo stesso istante pesano anche lo sfasamento e la portata diversa.

`confronto_reference.jl` confronta solo il periodo comune tra simulazione e riferimento e, sezione per
sezione, scarta gli istanti in cui il riferimento ha portata nulla (2D che parte da alveo asciutto).
Se il riferimento non ha ancora raggiunto il colmo quando finisce, il colmo non viene confrontato.
Le colonne del riferimento sono associate alle sezioni arrotondando l'ascissa ai 50 m (con le sezioni
a 200 m dà le stesse colonne di prima).
Tutti gli strumenti leggono tabelle e geometria dalle righe 1-2 di `settings.txt`: devono essere quelle
della simulazione che si sta analizzando.
Una simulazione del Test_Brenta richiede pochi secondi: circa 15 s con 30 sezioni, compresa la
compilazione; con 120 sezioni (50 m) 2-3 s di calcolo per 24-39 h simulate.

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
| `input/geometria/<caso>/` | Tabelle delle sezioni e geometria; per il Test_Brenta `Test_Brenta/S_200m/`, `Test_Brenta/S_50m/` (quella in uso, β = 1) e `Test_Brenta/S_50m_beta/` (β dal 2D, tenuta da parte). Quelle a 50 m le scrive direttamente `calcolo_tabelle_daBM` |
| `input/idrogrammi/<evento>/` | Idrogrammi `t [s], Q [m3/s]` |
| `output/<caso>/` | Risultati (`IC_*.csv`, `*.csv`) e grafici; i fittizi sono raggruppati per geometria in `output/S_200m/`, `output/S_50m/` e `output/S_50m_beta/` (prova con β dal 2D) |
| `reference/` | Risultati del modello 2D usati come riferimento: portate (`ns_channel_piena_discharge.csv`, `ns_50_fittizio_N_discharge.csv`) e, per i fittizi, livelli (`ns_50_fittizio_N_wse.csv`, colonne `WSE_S0050 [m]` … `WSE_S6000 [m]`, estratti da `../Test_Brenta/results/inputs/tmp/fittizi/ns_50_fittizio_N_results.h5`, gruppo `RESULTS/NodeStrg/StateVar`, colonna 1) |
| `calcolo_tabelle_daBM/` | Script `src/main.jl` che costruisce le tabelle 1D dai risultati 2D (BASEMENT) a scalini di portata e le scrive anche negli input del modello (vedi "Costruzione delle tabelle") |
| `strumenti/` | Script di controllo, grafici e confronto |
| `idromaker/main.jl` | Crea gli idrogrammi degli scenari fittizi in `input/idrogrammi/fittizio_*/piena.csv` (`Centa.csv` è solo un termine di paragone) |

## settings.txt

Letto **per posizione** (27 righe, non spostarle; dopo ci sono la legenda e vecchie configurazioni):

| Riga | Contenuto |
|---|---|
| 1 | Tabelle delle sezioni (es. `input/geometria/Test_Brenta/S_50m/tabella_sezioni.csv`) |
| 2 | Geometria (es. `input/geometria/Test_Brenta/S_50m/localslope_xcoord.csv`) |
| 3-13 | Idrogrammi affluenti Brenta1-6, BrentaNuovo, BrentaVecchio1-4 |
| 14 | Idrogramma in ingresso da monte (chiave interna `caldonazzo`) |
| 15-21 | Idrogrammi affluenti Ceggio, Levico, Centa, Chieppena, Grigno, Maso, Moggio |
| 22 | Cartella degli output (creata se non esiste) |
| 23 | Ks: **non più usato** (l'uscita usa la tabella dell'ultima sezione), ma la riga deve restare e contenere un numero |
| 24 | Tempo iniziale |
| 25 | Passo temporale massimo [s] |
| 26 | Numero massimo di passi |
| 27 | Affluenti attivi: `si` o `no` (accetta anche `sì`, `true`/`false`, `1`/`0`) |

- Numero di righe e colonne delle tabelle e numero di sezioni **sono ricavati dai file**
  (`FILE_SIZE`): non serve più il vecchio `struttura_output.txt`.
- Tutti gli idrogrammi indicati devono esistere, anche con gli affluenti spenti. `Brenta7.csv`, non
  ancora usato, non è più in `settings.txt`.
- La durata simulata è il tempo finale **più breve** tra gli idrogrammi letti.
- Il `settings.txt` committato è la configurazione di lavoro: geometria `input/geometria/Test_Brenta/S_50m/`,
  idrogramma `fittizio_1`, output in `output/S_50m/fittizio_1`, dt massimo 60 s, affluenti spenti.
  Lanciare `main.jl` così riproduce byte per byte `output/S_50m/fittizio_1` (verificato). Per gli
  altri fittizi cambiare righe 14 e 22. Prima di lanciare controllare sempre righe 1-2, 14, 22 e 27:
  gli strumenti leggono tabelle e geometria dalle righe 1-2.

## Formato degli input

**Tabelle delle sezioni** (CSV con intestazione, 12 colonne, righe raggruppate per sezione e con
quote crescenti): `sezione, WSE, Ω, Ω_s, β, β_s, jden, jden_s, W, W_s, RH, RH_s`.
- La colonna 2 è la **quota del pelo libero** (non il tirante). Il passo tra le quote può essere non costante.
- Le colonne `_s` sono le pendenze rispetto alla quota, usate per l'interpolazione lineare.
- `jden = Ks²·RH^(4/3)/g` è il denominatore del termine d'attrito (γ = |u|/jden). La scabrezza usata
  in moto vario è quindi **quella contenuta nelle tabelle**: in quelle da `calcolo_tabelle_daBM` è
  il Ks equivalente ricavato dal 2D, livello per livello. Anche lo stato iniziale e la condizione di
  uscita usano l'attrito delle tabelle (`UNIFORM_DISCHARGE`, `Q = A·√(g·jden·S)`).
- **Ultima sezione:** il suo `jden` non descrive l'attrito ma la condizione di uscita. In
  `calcolo_tabelle_daBM` il Ks è ricavato con la pendenza che il modello usa all'uscita (la `local_slope`
  dell'ultima sezione, cioè la pendenza media del fondo), così `Q = A·√(g·jden·S)` riproduce la scala di
  deflusso degli scalini 2D (vedi "Condizione al contorno di valle").
- `β` è il coefficiente di quantità di moto: 1 nelle tabelle in uso, calcolato dal 2D in quelle di
  `S_50m_beta` (vedi "Costruzione delle tabelle").

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
- **β** entra nel flusso convettivo (`β·u·q`) e nella velocità d'onda di Rusanov (`β·|u|`): con β > 1
  il passo temporale si accorcia (a 50 m da circa 25 a 21-23 s).

### Condizione al contorno di valle

- La portata in uscita è esplicita, di moto uniforme, con l'attrito della tabella dell'ultima sezione:
  `qstar[imax+1] = UNIFORM_DISCHARGE(A_int, den_γ_int, bx[imax])`, cioè `Q = A·√(g·jden·S)`, con area e
  `jden` dell'ultima sezione al passo precedente e S = `local_slope` dell'ultima sezione.
- Siccome `calcolo_tabelle_daBM` ricava il `jden` dell'ultima sezione proprio con quella pendenza, la
  portata in uscita segue la **scala di deflusso degli scalini 2D** in S6000, senza tarature a mano.
- **Prima** si usava Manning con il Ks della riga 23 di `settings.txt` (33):
  `Q = A·Ks·√S·Rh^(2/3)` (la riga è rimasta, commentata). A 50 m dava in S6000 tiranti 10-13 cm più alti
  del 2D.
- **Attenzione:** le tabelle delle altre geometrie (`S_200m`, rettangolare, `vaia`) hanno ancora il
  `jden` dell'ultima sezione ricavato con la pendenza del pelo libero. Con il codice attuale la loro
  uscita non è più quella dei risultati salvati (a 200 m la tabella dà un Ks di circa 18 invece di 33).
  Prima di riusarle vanno ricostruite come quelle a 50 m.

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
  al contorno di valle usava il Ks 33 di `settings.txt`, la tabella nell'ultima sezione circa 18
  (prova fatta prima della condizione di valle attuale).

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

La configurazione di lavoro attuale è quella a 50 m (vedi "Scenari fittizi con la geometria a 50 m").
I risultati di questa sezione e degli scenari a 200 m sono **precedenti** alla condizione di valle
attuale: sono fatti con il Ks 33 di `settings.txt` all'uscita.

**Test_Brenta a 200 m:**
- tratto di 6 km, 30 sezioni ogni 200 m (S0200-S6000), in `input/geometria/Test_Brenta/S_200m/`;
  tutti i risultati qui sotto sono con questa geometria (con quella a 50 m ci sono solo gli scenari
  fittizi, vedi "Scenari fittizi con la geometria a 50 m");
- tabelle da `calcolo_tabelle_daBM`;
- idrogramma `input/idrogrammi/prova_tratto_brenta/piena.csv` (5,07 → 83,34 m³/s a 7 h, durata 24 h);
- Ks 33 all'uscita, dt massimo 60 s, upwind, nessun affluente;
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

**Risultati salvati e confronto con il 2D.** Le simulazioni 1D con sfioro con la **geometria a 200 m**
(fatte con il codice di `7ba5b92`, partenza da moto uniforme con rampa di 6 h) sono in
`output/S_200m/fittizio_1` … `fittizio_4`, ciascuna con la sua cartella `grafici/` (grafici 1-12, vedi
"Come si usa"). Lo stesso vale per `output/S_50m/`.

**Tempi di calcolo** (`log_tempi.txt` in ogni cartella dei fittizi, a 200 e a 50 m). Rimisurati su copie
temporanee con il codice che ha prodotto ciascun risultato, aggiungendo solo le stampe dei tempi: i
risultati rilanciati sono identici byte per byte a quelli salvati. i9-9900K, Julia 1.12.6, 1 thread,
`-O2`, seconda esecuzione nello stesso processo (cioè senza compilazione):

| Geometria | Condizioni iniziali | Simulazione |
|---|---|---|
| 200 m (rampa, 7ba5b92) | 6,0 h simulate, 4321 passi, 0,51-0,63 s | 0,27-0,32 s (24-39 h) |
| 50 m (21b8f7a, Ks 33 all'uscita) | 2,5 h simulate, 1801 passi, 0,57-0,58 s | 1,20 s (24 h) - 2,01 s (39 h) |

A 50 m la simulazione costa 4-7 volte di più (4 volte le sezioni e passo da 60 a circa 28 s); le
condizioni iniziali costano circa lo stesso perché il nuovo arresto dimezza il tempo simulato. La
prima esecuzione, con la compilazione, aggiunge 0,1-0,2 s. I `log_tempi.txt` di `output/S_50m` si
riferiscono ai run con il Ks 33 all'uscita: con la condizione di valle attuale i passi sono circa il 2%
in più (3519 invece di 3453 nel `fittizio_1`), i tempi non sono stati rimisurati.

**Differenza di quota 1D − 2D** (dalle 3 h in poi, dopo l'arrivo del fronte 2D, su tutte le sezioni;
a 50 m con la condizione di valle attuale):

| Scenario | 200 m: media / scarto quadratico medio | 50 m: media / scarto quadratico medio |
|---|---|---|
| `fittizio_1` | −0,002 / 0,069 m | −0,051 / 0,109 m |
| `fittizio_2` | +0,007 / 0,055 m | −0,032 / 0,078 m |
| `fittizio_3` | −0,023 / 0,103 m | −0,105 / 0,207 m |
| `fittizio_4` | −0,013 / 0,083 m | −0,069 / 0,143 m |

A 50 m, a valle di 3 km il 1D sta circa 25 cm sotto il 2D: è l'acqua uscita da S2950 e S3550 che nel
1D non rientra (vedi sotto). Con il Ks 33 all'uscita le medie erano 3-4 mm meno negative, perché
l'uscita alzava i livelli degli ultimi chilometri.

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
quindi per il `fittizio_4` il grafico 8 e il "ritardo di 1065 min" di `confronto_S6000.txt` mescolano le
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
Il confronto è stato rifatto senza rilanciare la simulazione: con il codice attuale (partenza
cumulata, `y_ini = 0`, condizione di valle dalla tabella) un nuovo run a 200 m non riprodurrebbe i
risultati salvati. Gli strumenti leggono la geometria dalle righe 1-2 di `settings.txt`, che ora sono a
50 m: per i risultati a 200 m vanno lanciati da una copia con la geometria `S_200m`.

### Scenari fittizi con la geometria a 50 m

Risultati in `output/S_50m/fittizio_1` … `fittizio_4`, con grafici 1-13, animazioni e confronti.
Codice attuale: partenza dalle portate cumulate, nuovo arresto del riempimento, **condizione di valle
dalla tabella dell'ultima sezione**, β = 1, affluenti spenti (riga 27 `no`), sfioro attivo, dt
massimo 60 s. È la configurazione del `settings.txt` committato (con le righe 14 e 22 sullo scenario).

- Nessun NaN; riempimento finito a 2,5 h in tutti e quattro; bilancio di massa chiuso entro lo 0,15%;
  il livello non supera mai la quota massima delle tabelle (a 200 m la superava di 1-2 cm). Le
  pendenze del fondo irregolari non hanno dato problemi nello stato iniziale.
- Rispetto ai run precedenti con il Ks 33 all'uscita cambiano solo i livelli degli ultimi chilometri
  (in S6000 il tirante massimo scende da 1,89 a 1,73 m): portate, colmi, NSE e sfioro restano uguali
  entro 0,01 m³/s e qualche m³, perché l'uscita non influenza la parte a monte dove si esonda.
- **In uscita la portata si ferma a circa 55 m³/s in tutti gli scenari**, contro circa 70 a 200 m:

| Scenario | Colmo in S6000: 2D | 1D 200 m | 1D 50 m | NSE 50 m (200 m) | Volume sfiorato 50 m (200 m) |
|---|---|---|---|---|---|
| `fittizio_1` | 68,14 | 70,49 | 55,01 (−19%) | 0,935 (0,981) | 406 232 m³, 14,9% (4,4%) |
| `fittizio_2` | 64,16 | 69,50 | 54,49 (−15%) | 0,965 (0,988) | 266 779 m³, 10,8% (1,5%) |
| `fittizio_3` | 81,71 | 70,50 | 55,01 (−33%) | 0,785 (0,954) | 1 018 174 m³, 22,5% (8,8%) |
| `fittizio_4` (2ª onda) | 78,54 | 70,49 | 55,01 (−30%) | 0,871 (0,967) | 817 610 m³, 15,6% (4,6%) |

- **Fino a 2,9 km il colmo 1D segue quasi esattamente il 2D** (circa 68 m³/s nel `fittizio_1`), molto
  meglio che a 200 m. Poi due sezioni che a 200 m non c'erano fanno uscire molta acqua:
  - **S2950**: argine destro a 426,21 m, solo 2,33 m sopra il talweg (il più basso del tratto). Nel
    `fittizio_1` da sola sfiora 191 000 m³, quasi metà del totale, e il colmo scende a 58 m³/s;
  - **S3550**: argine sinistro a 425,11 m, altri 52 000 m³; il colmo scende a 55 m³/s.
- Nel 2D proprio lì il colmo **risale** da 65,5 a 68,1 m³/s: l'acqua rientra. Nel 1D esce e non torna,
  quindi da 3,6 km in poi resta circa 13 m³/s sotto. È lo stesso rientro dei `fittizio_3` e
  `fittizio_4`, qui molto più evidente: con le sezioni a 50 m pesa ancora di più la mancanza del rientro.
- Sfiorano 12 sezioni tra 0,25 e 3,55 km (9 nel `fittizio_2`), contro 5 tra 0,6 e 2 km a 200 m.
- Anche qui, per il `fittizio_4`, grafico 8 e ritardo in `confronto_S6000.txt` mescolano le due onde.

### Tirante e portata nelle sezioni (grafici 12 e 13)

Tirante 1D − scalini 2D **a parità di portata**, media nel tempo sugli istanti 2D bagnati (50 m,
condizione di valle attuale; per confronto il 2D sta entro ±0,4 cm dagli scalini in tutte le sezioni):

| Sezione | `fittizio_1` | `fittizio_2` | `fittizio_3` | `fittizio_4` |
|---|---|---|---|---|
| S1000 | −0,5 cm | −0,2 cm | −1,0 cm | −0,5 cm |
| S2150 | −0,6 cm | −0,5 cm | −1,1 cm | −0,9 cm |
| S4000 | −2,6 cm | −2,6 cm | −2,6 cm | −2,8 cm |
| S6000 | +0,1 cm | +0,2 cm | +0,1 cm | +0,1 cm |

- In S6000 con il Ks 33 all'uscita lo stesso scarto era +11 / +12 cm.
- **Allo stesso istante** in S1000 e S2150 1D e 2D coincidono quasi (tirante massimo 1-4 cm più basso,
  portata massima entro 0,8 m³/s). In S4000 e S6000 il 1D sta 14-64 cm sotto al colmo, ma solo perché
  ha meno portata (55 contro 64-82 m³/s): è l'acqua sfiorata a S2950 e S3550 che nel 1D non rientra.
- Il tirante massimo in S6000 (1D − 2D) va da −14 a −38 cm, con scarto quadratico medio 6-19 cm; a
  200 m (Ks 33) da −8 a +15 cm, 7-8 cm.

### Livelli 1D e 2D a parità di portata (questione chiusa)

Le tabelle vengono dal 2D, quindi ci si aspetterebbe che a parità di portata il 1D dia gli stessi livelli
del 2D. Con il Ks 33 all'uscita, nel grafico 12 non succedeva neanche sotto i 40 m³/s, dove non c'è
esondazione e le portate coincidono. Analisi fatta su copie temporanee. **Conclusione:** lo scarto
grande veniva dalla condizione di valle, ora corretta; quello che resta nel tratto è di 1-2 cm in
media e dipende dal Ks equivalente delle tabelle, non da β.

**Prima di tutto: confrontare a parità di portata, non allo stesso istante.** Sul ramo crescente il
1D arriva prima del 2D (laminazione e ritardo minori), quindi allo stesso istante è più alto anche se la
sua relazione livello-portata fosse esatta. Parte di ciò che si vede nel grafico 12 è solo questo
sfasamento.

**Verifica 1: piena dei fittizi**, istanti con portata 2D tra 5 e 40 m³/s (`fittizio_1`), in ogni sezione:
- livello della piena 2D − livello degli scalini 2D alla stessa portata: a 50 m in media 0,5 cm.
  **Il 2D è coerente con le tabelle** e non c'è effetto del moto vario (isteresi);
- livello 1D − livello degli scalini alla portata 1D: in media 1,7 cm nelle sezioni interne, ma
  **+9 cm in S6000** (+8 cm sul 2D). Lo scarto nasce nell'ultima sezione;
- a 200 m lo scarto medio tra 1D e 2D è di 3 cm, ma in parte viene dal riferimento: le tabelle a
  200 m sono costruite dagli scalini di `ns_channel`, mentre la piena di riferimento viene da `ns_50`, e le
  due simulazioni 2D non danno la stessa relazione livello-portata (fino a 11 cm in S4200 e S0800).

**Verifica 2: moto permanente.** Il 1D è stato portato all'equilibrio con portata costante di 20 e
40 m³/s (fase di riempimento) e il profilo è stato confrontato con gli scalini 2D alla stessa portata:

| Prova | Tirante in S6000 | Scarto medio sul tratto | Scarto massimo |
|---|---|---|---|
| 50 m, Ks 33 | +10 / +13 cm | 1,8-2,0 cm | +15 / +18 cm (S5900) |
| 50 m, Ks d'uscita tarato (40,2 / 39,4) | 0 | 1,8 cm | +9 / +10 cm (S5900) |
| 50 m, uscita dalla tabella (attuale) | 0 / 0 cm | 1,75 / 1,84 cm | +9 / +10 cm (S5900) |
| 50 m, uscita dalla tabella e β dal 2D | 0 / 0 cm | 1,69 / 1,71 cm | +9 cm (S5900) |
| 200 m, Ks 33 | 0 / +1 cm | 1,7 cm | −5 / −7 cm |

La condizione di valle attuale dà lo stesso profilo del Ks d'uscita tarato a mano, ma la taratura
viene dalla tabella.

**Causa 1: la condizione di uscita (lo scarto in S6000).**
- Nell'ultima sezione il livello non lo decide la tabella ma la portata di moto uniforme
  `Q = A·Ks·√i·Rh^(2/3)`, con il Ks di `settings.txt` (33) e la pendenza media del fondo (0,00213).
- Nel 2D, verso la fine del tratto, il pelo libero scende molto più ripido: tra S5950 e S6000 cala di
  26 cm in 50 m, cioè circa 0,005. Probabilmente è il richiamo dell'uscita del dominio 2D.
- Con Ks 33 e 0,00213 il 1D dava quindi tiranti troppo alti: a 50 m 8-13 cm tra 10 e 40 m³/s.
- Il Ks d'uscita che riproduce il livello 2D in S6000 è circa 39,5 a 50 m, quasi costante da 15 a
  80 m³/s. A 200 m è 32,7-33,9: lì Ks 33 è giusto per caso, e per questo il problema non si vedeva.
- L'effetto risale per qualche centinaio di metri: con il Ks tarato lo scarto medio sugli ultimi 2 km
  scende da 2,0-2,3 a 1,8-1,9 cm.
- **Corretta:** ora la portata in uscita usa il `jden` dell'ultima tabella, ricavato con la pendenza
  d'uscita (Ks circa 39,5): vedi "Condizione al contorno di valle". In S6000 il 1D sta sulla scala di
  deflusso degli scalini 2D (+0,1 / +0,2 cm nei fittizi) e il Ks della riga 23 non serve più.

**Causa 2: il Ks equivalente delle tabelle (lo scarto che resta nel tratto, 1-2 cm in media).**
- `calcolo_tabelle_daBM` ricava il Ks dalla pendenza del pelo libero calcolata tra le due sezioni
  vicine, a monte e a valle, e trascura la variazione dell'altezza cinetica `u²/2g`.
- Il 1D invece mette l'attrito sull'interfaccia, con `jden` della cella a monte (upwind), e risolve
  anche il termine convettivo. A moto permanente impone quindi tra due sezioni un dislivello che non
  coincide con quello del 2D da cui è stato ricavato il Ks.
- A 50 m `u²/2g` cambia in media di 2 cm tra sezioni vicine, fino a 8,5 cm.
- Il profilo 2D del livello medio sulla linea di sezione (`Mean wse`) è irregolare: tra S5850 e S5900
  cala di 39 cm, tra S5900 e S5950 risale di 1 cm.
- A moto permanente il 1D non può far risalire il pelo libero con il solo attrito, e in questi punti
  sbaglia di 5-10 cm (S5900, S2400, S2950). Gli errori non si sommano lungo il tratto: in media si
  compensano.
- Per ridurlo servirebbe ricavare il Ks dalla pendenza dell'energia (pelo libero più `u²/2g`), oppure
  tarare `jden` direttamente sullo schema discreto del 1D. Prima andrebbero rese più regolari le
  sezioni in cui il profilo 2D non è monotono. Per ora non si fa: 1-2 cm sono accettabili.

**Ipotesi 3 (scartata): β = 1.** Le tabelle usavano β = 1 per ipotesi. β è stato calcolato dal 2D come
nel metodo di Engelund (vedi "Costruzione delle tabelle"): vale in mediana 1,14, tra 1,0 e 1,7 sopra
5 m³/s, e cambia molto tra sezioni vicine (1,05-1,4).
- A moto permanente lo scarto medio cala solo di 0,5-1,4 mm (tabella sopra); localmente il livello
  si sposta di ±2-3 cm (S3350-S3400, S3750), in un verso o nell'altro, e in media si compensa.
- Nei fittizi (`output/S_50m_beta/`, fatti con il Ks 33 all'uscita e confrontati con i run di allora):
  colmo in S6000 +0,2-0,3 m³/s, NSE +0,001-0,005, sfioro −1%, livelli spostati in media di 0,5 cm
  (al massimo −4 cm in S0050 e +3,6 cm in S3350). Grafici di confronto in
  `output/S_50m_beta/confronto_beta/fittizio_N.png`.
- **Si resta con β = 1.** Le tabelle con β sono in `input/geometria/Test_Brenta/S_50m_beta/` (già con
  la condizione di valle attuale), le simulazioni in `output/S_50m_beta/` no.

### Rientro dell'acqua nel 2D

Nel 2D in alcune sezioni la portata cresce rispetto alla sezione subito a monte: l'acqua uscita più a
monte rientra. Analisi delle portate `reference/ns_50_fittizio_*_discharge.csv`, ogni 50 m. Una sezione
è segnata quando cresce il colmo (più di 0,3 m³/s) o il volume transitato in tutto l'evento (più dello
0,3% del volume entrato).

Il segno più sicuro è una crescita **mentre la piena sta ancora salendo** (portata da monte in
aumento). In quella fase il tratto si riempie e trattiene acqua, quindi senza ingressi laterali la
portata cala andando verso valle: se cresce, l'acqua arriva dal lato. Quando la piena cala, invece, il
tratto restituisce l'acqua immagazzinata e la portata cresce verso valle anche senza rientro.

| Sezioni | Crescita del colmo [m³/s] | Scenari | Note |
|---|---|---|---|
| S0400 | +1,1 / +1,8 | tutti | anche con la piena in salita; dopo le uscite di S0150-S0350 |
| S0950 (e S1100) | +0,4 / +0,8 | tutti | |
| S1600-S1700 | fino a +3,1 (S1700) | tutti (S1600 non nel `fittizio_2`) | anche con la piena in salita; dopo le uscite di S1000-S1400 |
| S3000 | +1,4 / +1,6 | tutti | anche con la piena in salita; il 2D perde 25-66 mila m³ in S2950 e li recupera 50 m dopo |
| S3550 (e S3600) | +2,3 / +4,6 | 1, 3, 4 | nel `fittizio_2` S3550 perde acqua |
| S2700, S3700-S3750, S4000 | fino a +1,3 | 3, 4 | |
| **S4150** | +6,5 / +7,2 | 3, 4 (nel 1 solo +11 000 m³ di volume) | il rientro maggiore, fino a 167 000 m³; il colmo sale da circa 74 a 83 m³/s |
| S4250-S4700, S5650-S5700 | fino a +2,5 | 3, 4 | |

Nel `fittizio_4` crescono anche con la piena in salita S0850 e S0950. Rientra circa metà dell'acqua
persa nel `fittizio_1` (175 000 m³ su 356 000) e circa l'80% nel `fittizio_3` (879 000 su 1 093 000).
Nel 1D a 50 m, invece, S2950 da sola sfiora 191 000 m³ che non tornano.

Due avvertenze:
- le portate 2D sono misurate su linee di sezione: un aumento può anche voler dire che l'acqua ha
  aggirato la sezione precedente passando per la piana. Per il 1D è comunque acqua che esce e torna;
- l'abbinamento tra uscite e rientri è dedotto dalla vicinanza delle sezioni, non tracciato.

## Costruzione delle tabelle (`calcolo_tabelle_daBM`)

`julia --project=. calcolo_tabelle_daBM/src/main.jl` dalla cartella del progetto. Lavora solo sul caso
a 50 m (`input/S_0050_6000/`, prefisso `ns_50_scalini_`, 120 sezioni S0050-S6000). L'interruttore
`β_dal_2D` in testa al file sceglie β:

| `β_dal_2D` | β | Output dello script | Input del modello (scritti direttamente) |
|---|---|---|---|
| `false` (attivo) | 1 | `calcolo_tabelle_daBM/output/S_50m/` | `input/geometria/Test_Brenta/S_50m/` |
| `true` | dal 2D | `calcolo_tabelle_daBM/output/S_50m_beta/` | `input/geometria/Test_Brenta/S_50m_beta/` |

- `tabella_sezioni.csv` e `localslope_xcoord.csv` sono scritti sia nell'output dello script sia negli
  input del modello (`ModelFolder`), quindi sono sempre uguali: non vanno più copiati a mano.
  Rilanciato con i dati attuali lo script non cambia nessun file.
- **Controllo dei file all'avvio:** se manca un file di input (`sections_dx.csv`, `bankpoints.csv`,
  i file delle sezioni e, con β dal 2D, mesh e risultati 2D) lo script li elenca tutti in rosso e si ferma.
- Il caso a 200 m (`input/S_0200_6000/`, prefisso `ns_channel_`) non è più configurato nello script:
  le sue tabelle in `input/geometria/Test_Brenta/S_200m/` sono quelle vecchie (vedi "Condizione al
  contorno di valle").
- La cartella di input contiene i file delle sezioni, `sections_dx.csv` (nome e ascissa),
  `bankpoints.csv` (quote degli argini) e, per β, la mesh `test_brenta_ns_50_computational-mesh.2dm`
  e i risultati `ns_50_scalini_results.h5` del 2D a scalini. Le sezioni BC e P17 non sono usate.
- I file a 50 m non hanno la colonna `t [s]`: il controllo che le sezioni abbiano gli stessi scalini
  usa il numero di righe, più il tempo quando c'è.
- In output: una cartella per sezione (`tabella.csv` e grafici, riscritti solo se la tabella cambia;
  `beta.svg` solo con β dal 2D), `tabella_sezioni.csv` e `localslope_xcoord.csv`.
- **Ultima sezione:** Ks e `jden` sono ricavati con la pendenza d'uscita del modello (`local_slope`
  dell'ultima sezione, pendenza media del fondo 0,00213) invece che con la pendenza del pelo libero
  (in S6000 tra 0,0014 e 0,0075, molto ripida per il richiamo dell'uscita del dominio 2D). Così il
  modello in uscita riproduce le portate degli scalini 2D. Il Ks che ne esce è circa 39,5 tra 15 e 80
  m³/s (7 nella riga asciutta ricostruita, 19 a 0,85 m³/s; nessun avviso). Nella tabella della sezione
  la colonna `iWSE` riporta quella pendenza.

**β dal modello 2D** (`beta_2D`, solo con `β_dal_2D = true`), come nel metodo di Engelund
(`β = Ω/Q² ∫ u² Y dy`) ma con le velocità del 2D:
- dalla mesh si leggono nodi, triangoli e nodestring; ogni nodestring segue i lati della mesh
  (3-12 lati per sezione, lunghi 2-5 m, solo il canale principale);
- per ogni lato si prende il triangolo **a monte** (dalla parte indicata da `upstream_direction` nella
  configurazione del 2D, `right` per tutte le sezioni); tirante Y = quota del pelo libero − fondo
  dell'elemento, portata normale `qn = q·n`, larghezza ΔL = lunghezza del lato. Elementi con tirante
  non oltre il minimo del 2D (0,05 m) esclusi;
- `Ω = Σ Y ΔL`, `Q = Σ qn ΔL`, `β = Ω·Σ(qn²/Y·ΔL)/Q²`. Ω coincide con l'area dei file delle sezioni (scarti
  sotto 1,4·10⁻⁴, arrotondamenti di BASEMENT): BASEMENT usa gli stessi elementi. Q no (in mediana il 2%
  in più, in S5900 il 12-30%), perché BASEMENT la calcola dai flussi sui lati: si usa la Q degli
  elementi, così β ≥ 1 e vale 1 con velocità uniforme;
- controlli bloccanti: la mesh deve corrispondere ai risultati, le righe dei file delle sezioni alle
  uscite dell'`.h5`, l'area degli elementi a quella delle sezioni. Avvisi: portata degli elementi oltre
  il 10% da quella della sezione (47, quasi tutti sotto 5 m³/s, e S5900), β fuori da 1-1,5 (12, quasi
  tutti a circa 1 m³/s; sopra 1,5 a 5-10 m³/s solo S3550, S3800, S4500, S5900);
- righe asciutte ricostruite con il β del primo livello bagnato;
- risultato: β mediano 1,14, fino a 4,2 nella prima riga bagnata (circa 1 m³/s). Il pacchetto `HDF5`
  è nell'ambiente del progetto. L'effetto sul modello è trascurabile (vedi "Ipotesi 3" in "Livelli 1D
  e 2D a parità di portata").
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
  2·10⁻⁴, fino a 7·10⁻⁶ in S1800 (l'ultima sezione usa sempre la pendenza media, −0,00213, che è anche
  quella della condizione di uscita). Lo stato
  iniziale usa il moto uniforme con `|local_slope|`: a 5,07 m³/s il tirante supera 1,5 m in 7 sezioni
  (massimo 1,32 m a 200 m) e in S1800 il livello andrebbe oltre la tabella. Negli scenari fittizi il
  riempimento lo ha comunque assestato senza problemi (2,5 h); per portate iniziali diverse resta da
  valutare una pendenza più regolare per lo stato iniziale.

## Problemi aperti e prossimi passi possibili

- **Esondazione, il punto principale per avvicinare il 1D al riferimento.** Includere le golene nelle
  tabelle sarebbe la strategia migliore, ma non è banale, soprattutto dove l'acqua in golena arriva da
  monte invece di uscire in quella sezione. Si procede quindi, un passo alla volta, con una **portata
  laterale uscente** (vedi "Esondazione"). I dati 2D di riferimento non sono ancora definitivi: per
  ora interessa l'effetto del termine, non far combaciare gli idrogrammi. La corrispondenza dei
  livelli a parità di portata è invece considerata **chiusa** (vedi "Livelli 1D e 2D a parità di portata").
- **Upwind:** l'instabilità con area e attrito insieme su veli sottili e tratti lunghi è da
  risolvere, oppure va scelto un upwind solo parziale.
- **Stato iniziale vicino all'equilibrio** nei tratti non prismatici: profilo di moto permanente,
  o almeno la pendenza tra i talweg al posto di `local_slope`.
- **Criterio di arresto del riempimento:** la tolleranza assoluta 0,01 m³/s su `|q − qn|` è poco
  adatta alle portate molto basse; `T_assestamento`, `T_fermo` e `tol_fermo` sono provati solo su
  Test_Brenta (6 km) e rettangolare (30 km).
- **Rientro dell'acqua sfiorata: il prossimo passo.** Nel 2D l'acqua uscita rientra più a valle (vedi
  "Rientro dell'acqua nel 2D"); nel 1D esce e non torna, ed è ormai la causa principale dello scarto
  con il 2D (da 3 km in poi il 1D ha circa 13 m³/s e 25 cm in meno). L'idea è includerlo **senza
  modellare la piana alluvionale**. Ragionamenti fatti finora, niente è deciso:
  - il meccanismo per far rientrare l'acqua c'è già: come gli affluenti, una portata laterale entrante
    nella sola continuità (`rhs += dt/dx·Q`), senza quantità di moto;
  - servono tre cose per ogni uscita: **dove** rientra (sezione o gruppo di sezioni a valle), **quanto**
    (frazione del volume uscito: nel 2D circa metà nel `fittizio_1`, 80% nel `fittizio_3`; il resto va
    perso altrove) e **quando** (ritardo);
  - possibili schemi, dal più semplice: (1) rientro con ritardo puro, `Q_r(t) = α·Q_lat(t − τ)`;
    (2) un **serbatoio virtuale** per ogni zona di uscita (solo un volume, senza geometria della piana)
    che si svuota verso la sezione di rientro come serbatoio lineare, `Q_r = V/T`; (3) come il (2) ma
    con lo svuotamento che dipende dal livello del fiume nella sezione di rientro (sfioro inverso, con
    la stessa formula dello sfioro), così l'acqua torna solo quando il fiume si abbassa;
  - nel 2D ci sono due comportamenti diversi: rientri a breve distanza già con la piena in salita
    (S0400, S1600-S1700, S3000, cioè ritardo piccolo e poco accumulo) e il grande rientro di S4150 nei
    `fittizio_3` e 4, che arriva solo dopo che le aree esondate si sono riempite (serve un accumulo).
    Lo schema (1) basterebbe per i primi, per il secondo serve almeno il (2) o il (3);
  - parametri (α, τ, T, abbinamento uscita-rientro) da ricavare dalle portate 2D sezione per sezione
    (volumi persi e recuperati, istanti), oggi abbinati solo per vicinanza; il bilancio di massa del
    riempimento e della simulazione dovrà includere i volumi nei serbatoi;
  - il caso S2950 → S3000 (il 2D perde 25-66 mila m³ e li recupera 50 m dopo) sembra più acqua che
    aggira la linea di sezione che una vera esondazione: dipende da come sono definite le nodestring
    (punto seguente).
- **Definizione delle nodestring del 2D.** Oggi ogni nodestring va da sponda a sponda del solo canale
  principale (15-25 m, 3-12 lati di mesh): tabelle, Ks e portate 2D di riferimento descrivono solo il
  canale, e per il 1D tutto ciò che supera gli argini è acqua uscita. L'alternativa è **prolungare la
  nodestring nella piana fino a dove la pendenza del terreno cambia e l'acqua sceglie una direzione
  anziché l'altra** (uno spartiacque locale): la parte di piana che scorre verso valle insieme al fiume
  farebbe parte della sezione, quella oltre lo spartiacque sarebbe davvero persa. Ragionamenti:
  - pro: l'acqua che aggira una sezione e rientra subito (S2950 → S3000 e simili) diventerebbe portata
    della sezione, non uscita più rientro; le tabelle includerebbero la golena attiva (sezione
    composta) e una parte del rientro sarebbe compresa senza modellare la piana; lo sfioro resterebbe
    solo dove l'acqua se ne va davvero;
  - contro e dubbi: serve rifare la simulazione 2D a scalini (e i riferimenti dei fittizi) con le nuove
    nodestring; il punto in cui fermarle va definito con un criterio (DTM, linee di displuvio,
    direzione delle velocità 2D agli scalini) e può cambiare con la portata; le nodestring non devono
    incrociarsi e devono restare circa perpendicolari al moto; nelle sezioni composte il Ks
    equivalente è più variabile, β diventa più importante (canale e golena hanno velocità molto
    diverse: il calcolo di β dal 2D è già pronto) e l'ipotesi di pelo libero orizzontale lungo la
    sezione regge meno; quote di sfioro (argini) e tabelle andrebbero ridefinite agli estremi nuovi;
  - una via di mezzo: prolungare solo le nodestring dove il 2D mostra aggiramenti e rientri a breve
    distanza.
- **Affluenti e sfioro** insieme: non ancora provati.
- **Ks equivalente delle tabelle** ricavato dalla pendenza del pelo libero senza l'altezza cinetica:
  lascia scarti di 1-2 cm in media e fino a 10 cm dove il profilo 2D è irregolare. Per ora accettato.
- **Geometrie non aggiornate alla condizione di valle attuale:** `S_200m`, rettangolare e `vaia` hanno
  ancora l'ultima tabella con il Ks dalla pendenza del pelo libero; vanno ricostruite prima di
  riusarle. Anche le simulazioni in `output/S_50m_beta/` sono state fatte con il Ks 33 all'uscita.
- **Varie:**
  - Newton fa sempre 100 iterazioni (tolleranza 1e-14);
  - `THOMAS_ALGORITHM!` alloca due vettori a ogni chiamata;
  - `strumenti/controllo_input.jl` calcola il tirante iniziale con la sola portata di monte e un
    minimo di 10 cm, non con la partenza cumulata, e stampa ancora il Ks della riga 23 come se fosse usato;
  - `exp_dt` e il campo `Ks` di `Parametri` restano nel codice ma non sono più usati.

## Ultimo commit descritto

`e034518` ("aggiunto calcolo di beta"), più le modifiche successive non ancora committate quando
questo file è stato scritto (12/11/2026), da committare insieme a esso:
- condizione di valle dalla tabella dell'ultima sezione (`Source/Semi_Implicit_Brenta_2D_1D.jl`,
  commenti in `Source/module.jl`), Ks della riga 23 non più usato;
- `calcolo_tabelle_daBM/src/main.jl`: interruttore `β_dal_2D` (β = 1 di nuovo in uso), Ks dell'ultima
  sezione con la pendenza d'uscita, scrittura diretta negli input del modello; tabelle rigenerate in
  `calcolo_tabelle_daBM/output/S_50m*` e `input/geometria/Test_Brenta/S_50m*`;
- `settings.txt` sulla configurazione a 50 m (`fittizio_1`, affluenti spenti);
- fittizi a 50 m rifatti in `output/S_50m/` con i grafici 1-13; `strumenti/confronto_valle.jl` con il
  grafico 13 e il confronto a parità di portata.

Commit descritti in precedenza: fino a `21b8f7a`, poi `c753289`-`373ba48` (geometria a 50 m, risultati
dei fittizi, grafici), `d7ee4e6` e `282eaf2` (analisi dei livelli e della condizione di valle),
`e034518` (calcolo di β e prova `S_50m_beta`). Strumenti utili per i riempimenti:
`strumenti/confronto_spinup.jl` (`<cartella dei grafici> <cartella 1> <cartella 2> ...`) e
`strumenti/confronto_simulazioni.jl` (stessi argomenti, per le simulazioni; non distingue cartelle con
lo stesso nome e salta le differenze se gli istanti sono diversi).
