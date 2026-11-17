# Calibrazione dello scambio laterale con la piana (soglia, T, ritardo, rilascio lento, serbatoio attivo)

Codici, risultati e grafici della taratura dei parametri del **serbatoio della piana**, usato dallo scambio
laterale `tabella` (riga 28 di `settings.txt`). Tutte le prove sono sui fittizi 1-4 con le tabelle
`input/geometria/Test_Brenta/S_50m_q2/` (scalini 2D a 2 m³/s, con le colonne ΔQ e ΔQ_s), affluenti spenti,
dt massimo 60 s, confronto con i riferimenti 2D `reference/ns_50_fittizio_N_*.csv`.

## Il metodo

- `calcolo_tabelle_daBM` aggiunge alle tabelle **ΔQ = Q della sezione − Q della sezione precedente** negli
  scalini 2D (moto permanente), riga per riga: negativo dove l'acqua esce verso la piana, positivo dove rientra.
- Nel modello (`scambio = tabella`) ΔQ è letto dalla tabella al livello della sezione ed entra nella sola
  continuità, come gli affluenti; sostituisce lo sfioro sugli argini.
- **Senza limiti i rientri crescono con il livello, che cresce con i rientri**: nel `fittizio_1` la portata
  a 4-5 km arrivava a 95 m³/s con 81 in ingresso e il modello divergeva a 6,4 h. Per questo:
  - le **uscite** sono quelle della tabella e finiscono in un **serbatoio virtuale unico** per tutto il tratto;
  - i **rientri** sono quelli della tabella (solo dove ΔQ > 0), ridotti tutti della stessa frazione in modo da
    non superare quanto il serbatoio può restituire.
- Parametri del serbatoio (righe 29-33 di `settings.txt`):

  | Riga | Parametro | Significato |
  | 29 | soglia [m³] | la piana trattiene i primi m³ usciti (si riempie): rientra solo il volume oltre la soglia |
  | 30 | T [h] | tempo di svuotamento: rientro al massimo (V − soglia)/T; 0 = tutto subito |
  | 31 | τ [h] | ritardo di percorrenza: l'acqua uscita arriva al serbatoio dopo τ |
  | 32 | T_lento [h] | rilascio lento: la parte sotto la soglia torna al fiume in T_lento, nelle sezioni di rientro (pesi: ΔQ massimo positivo della tabella), anche quando la tabella non chiede rientri; 0 = resta nella piana |
  | 33 | serbatoio attivo (si/no) | la parte oltre la soglia rientra comunque al ritmo (V − soglia)/T: quello che la tabella non chiede (ramo calante) va nelle sezioni di rientro con gli stessi pesi del rilascio lento |
- A fine simulazione il modello scrive `bilancio.txt` nella cartella degli output: volumi in alveo, entrati,
  usciti e scambiati, errore del bilancio e, con lo scambio `tabella`, il controllo del serbatoio della piana.
  In tutte le prove l'errore è 0,0000% e il serbatoio torna al m³.

## Risultati

Valori tarati: **soglia 150 000 m³, T = 0,5 h, τ = 1 h, rilascio lento 32 h, serbatoio attivo** (NSE medio in S6000
massimo; il rilascio lento è stato tarato dopo, con gli altri tre fissi, e il serbatoio attivo per ultimo, ritarando
T e soglia). Medie sui fittizi 1-4:

| Simulazione | NSE in S6000 | errore del colmo lungo il tratto | errore del livello |
|---|---|---|---|
| sfioro sugli argini | 0,937 | – | – |
| tabella, serbatoio senza parametri (T = 0) | 0,926 | 4,18 m³/s | 5,9 cm |
| tabella, solo T = 4 h | 0,986 | 1,86 m³/s | 3,6 cm |
| tabella, soglia 150 000 m³, T = 1 h | 0,989 | 1,11 m³/s | 3,3 cm |
| tabella, soglia 150 000 m³, T = 0,5 h, τ = 1 h | 0,991 | 1,02 m³/s | 3,1 cm |
| tabella, come sopra con rilascio lento in 32 h | 0,994 | 1,05 m³/s | 3,2 cm |
| **tabella, come sopra con il serbatoio attivo** | **0,998** | **1,08 m³/s** | **2,7 cm** |

- Colmo in S6000 (2D → 1D tarato con rilascio lento): 69,3 → 69,9; 63,6 → 63,4; 81,7 → 81,8; 78,6 → 77,9 m³/s
  (`fittizio_4`: seconda onda). Con lo sfioro era sempre circa 60.
- **Senza soglia** (grafici 1-3) tutta l'acqua uscita rientra subito: il colmo in uscita resta quello in ingresso.
  Con il solo T il rientro è graduale e troppo presto, mentre nel 2D arriva a scatto, quando la piana si è riempita.
- Nel 2D il volume che resta fuori (uscito − rientrato, per tratto) è simile nei quattro scenari, 154-215 mila m³:
  da qui l'idea della soglia.
- **Il ritardo** cambia poco gli indici medi (le combinazioni migliori sono tutte tra 0,989 e 0,991) ma porta il
  rientro al momento giusto: in S4200 la portata supera 67 m³/s 44-50 min prima del 2D senza ritardo, 1-7 min
  dopo con τ = 1 h (`risultati/anticipo_*.txt`).
- Il minimo è largo: T da 0,5 a 2 h, soglia 100-150 mila m³ e τ 0,5-1,5 h danno NSE medi tra 0,989 e 0,991.
- **Senza rilascio lento** il serbatoio non si svuota mai sotto la soglia: alla seconda onda del `fittizio_4` la
  piana risulta già piena e il rientro arriva 112 min prima del 2D; nel ramo calante il 1D scende 39 min prima
  (in media, portata in S6000 sotto 50 m³/s).
- **Rilascio lento** (griglia T_lento 4-64 h × soglia 75-175 mila m³, con T e τ fissi, grafici 9-11): con 32 h e
  soglia 150 000 m³ l'anticipo della seconda onda scende da 112 a 23 min e quello della discesa da 39 a 32 min;
  nelle prime onde il rientro arriva 4-11 min dopo il 2D. Nel `fittizio_2` ora rientrano 47 000 m³ (2D 51 000).
  Con un rilascio più veloce (4 h) la discesa è quasi esatta (6-8 min) ma il rientro arriva 48 min tardi: la
  piana si svuota troppo presto e va riempita di nuovo. Tempi oltre 32 h non migliorano.
- **Serbatoio attivo** (riga 33; griglia T 0,25-8 h × soglia 125-250 mila m³, τ 1 h e rilascio lento 32 h fissi,
  `risultati/attivo.csv`, `riepilogo_attivo.txt`, grafici 12-13). Senza, dopo il colmo la tabella smette di chiedere
  rientri e l'acqua oltre la soglia resta ferma nel serbatoio: il 1D scendeva 32 min prima del 2D. Il migliore resta
  T = 0,5 h e soglia 150 000 m³: NSE medio da 0,994 a 0,998, errore del livello da 3,2 a 2,7 cm, anticipo della
  discesa in S6000 da 32 a 14 min (per scenario 36 → 20, 10 → 10, 36 → 4, 36 → 20 e 39 → 15 min), errore
  sull'istante del rientro in S4200 da 12 a 7 min. Colmo in S6000 69,9 / 63,4 / 81,8 / 77,3 m³/s (2D 69,3 / 63,6 /
  81,7 / 78,6). Con T o soglia più grandi il serbatoio restituisce acqua prima del colmo e i colmi calano.

## Limiti

- Nel ramo calante il 1D scende ancora prima del 2D: circa 30 min senza serbatoio attivo, 14 min in media con
  (20 min nel `fittizio_1`, in parte per l'invaso del canale che nel 1D manca anche a monte del rientro).
- Il 1D fa uscire il 10-35% di acqua in meno del 2D (volumi netti per tratto).
- Un solo serbatoio per tutto il tratto: l'acqua uscita a 0,5 km può rientrare a 4 km. Per il Brenta intero
  serviranno serbatoi per zona o parametri tratto per tratto.
- **Taratura e verifica sugli stessi quattro scenari**, senza una piena di controllo indipendente.

## File

| File | Contenuto |
|---|---|
| `esegui.jl` | `simula(...)`: lancia il modello (`Source/` del progetto) in una cartella temporanea con il suo `settings.txt`, senza toccare quello del progetto |
| `indici.jl` | Lettura del 2D, indici di confronto (colmo, NSE, errori di colmo e livello) e volumi scambiati |
| `taratura.jl` | Lancia le combinazioni di τ, soglia e T e aggiunge le righe a `risultati/taratura.csv` |
| `simula_fittizi.jl` | Lancia i fittizi 1-4 con uno scambio e dei parametri, uscite in `out/<nome>/` |
| `grafici_taratura.jl` | Grafici 2, 4, 6 e `risultati/riepilogo_taratura.txt` da `risultati/taratura.csv` |
| `confronto.jl` | Colmo lungo il tratto e portata in S6000 (grafici 1 e 8), con indici e volumi |
| `idrogrammi.jl` | Idrogrammi in S3500, S4200 e S6000 (grafici 3, 5, 7) |
| `anticipo.jl` | Istante del rientro in S4200, 2D e 1D (`risultati/anticipo_<nome>.txt`) |
| `rilascio_lento.jl` | Lancia le combinazioni di T_lento e soglia (T e τ fissi) e aggiunge le righe a `risultati/rilascio_lento.csv`, con l'anticipo di rientro e discesa per onda |
| `grafici_rilascio_lento.jl` | Grafico 9 e `risultati/riepilogo_rilascio_lento.txt` |
| `attivo.jl` | Lancia le combinazioni di T e soglia con il serbatoio attivo (τ e T_lento fissi) e aggiunge le righe a `risultati/attivo.csv`, con l'anticipo di rientro e discesa per onda |
| `riepilogo_attivo.jl` | `risultati/riepilogo_attivo.txt` da `risultati/attivo.csv` |
| `risultati/taratura.csv` | Una riga per combinazione e scenario: 73 combinazioni (11 di solo T, 25 di soglia e T, 36 con il ritardo, più 1 prova) |
| `risultati/rilascio_lento.csv` | Una riga per combinazione, scenario e onda: 21 combinazioni di T_lento e soglia (T 0,5 h, τ 1 h) |
| `risultati/attivo.csv` | Una riga per combinazione, scenario e onda: 20 combinazioni di T e soglia con il serbatoio attivo più la versione senza (τ 1 h, T_lento 32 h) |
| `out/` | Uscite delle simulazioni (non nel repository: si rigenerano con gli script) |

Grafici in `grafici/`:

| N. | File | Contenuto |
|---|---|---|
| 1 | `1_serbatoio_T0_confronto.png` | Serbatoio senza parametri (T = 0) contro sfioro e 2D |
| 2 | `2_taratura_T.png` | Solo il tempo di svuotamento T, con i volumi scambiati |
| 3 | `3_idrogrammi_T.png` | Idrogrammi con T = 4 h senza soglia |
| 4 | `4_taratura_soglia.png` | Mappe di soglia e T, senza ritardo |
| 5 | `5_idrogrammi_soglia.png` | Idrogrammi con soglia 150 000 m³ e T = 1 h |
| 6 | `6_taratura_ritardo.png` | Mappe di soglia e T per ogni ritardo |
| 7 | `7_idrogrammi_tarato.png` | Idrogrammi con i valori tarati |
| 8 | `8_colmo_tarato.png` | Colmo lungo il tratto e portata in S6000 con i valori tarati |
| 9 | `9_taratura_rilascio_lento.png` | Mappe di T_lento e soglia: NSE, anticipo del rientro e della discesa |
| 10 | `10_idrogrammi_rilascio_lento.png` | Idrogrammi con i valori tarati, senza e con rilascio lento in 32 h |
| 11 | `11_colmo_rilascio_lento.png` | Colmo lungo il tratto e portata in S6000 con il rilascio lento |
| 12 | `12_idrogrammi_attivo.png` | Idrogrammi con i valori tarati, senza e con il serbatoio attivo |
| 13 | `13_colmo_attivo.png` | Colmo lungo il tratto e portata in S6000 con il serbatoio attivo |

## Come rifare tutto

Dalla cartella del progetto (ogni simulazione dura circa 15 s):

```bash
julia --project=. -O2 calibrazione_scambio/taratura.jl 0 0 0,0.25,0.5,1,2,3,4,6,8,12,24
julia --project=. -O2 calibrazione_scambio/taratura.jl 0 50000,100000,150000,200000,250000 0.25,0.5,1,2,4
julia --project=. -O2 calibrazione_scambio/taratura.jl 0.5,1,1.5,2 50000,100000,150000 0.5,1,2
julia --project=. calibrazione_scambio/grafici_taratura.jl
julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl sfioro sfioro
julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl solo_T0 tabella 0 0 0
julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl solo_T4 tabella 0 4 0
julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl soglia_T1 tabella 150000 1 0
julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl tarato tabella 150000 0.5 1
julia --project=. calibrazione_scambio/confronto.jl 1_serbatoio_T0_confronto solo_T0 "1D, tabella ΔQ, serbatoio con T = 0"
julia --project=. calibrazione_scambio/idrogrammi.jl 3_idrogrammi_T.png solo_T4 "1D, tabella ΔQ, T = 4 h, senza soglia"
julia --project=. calibrazione_scambio/idrogrammi.jl 5_idrogrammi_soglia.png soglia_T1 "1D, tabella ΔQ, soglia 150 000 m³, T = 1 h"
julia --project=. calibrazione_scambio/idrogrammi.jl 7_idrogrammi_tarato.png tarato "1D, tabella ΔQ, soglia 150 000 m³, T = 0,5 h, τ = 1 h"
julia --project=. calibrazione_scambio/confronto.jl 8_colmo_tarato tarato "1D, tabella ΔQ, soglia 150 000 m³, T = 0,5 h, τ = 1 h"
julia --project=. calibrazione_scambio/anticipo.jl soglia_T1 tarato
julia --project=. -O2 calibrazione_scambio/rilascio_lento.jl 0,8 150000 0.5 1
julia --project=. -O2 calibrazione_scambio/rilascio_lento.jl 4,8,16,32 75000,100000,125000,150000 0.5 1
julia --project=. -O2 calibrazione_scambio/rilascio_lento.jl 48,64 150000,175000 0.5 1
julia --project=. calibrazione_scambio/grafici_rilascio_lento.jl
julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl tarato_lento tabella 150000 0.5 1 32
julia --project=. calibrazione_scambio/idrogrammi.jl 10_idrogrammi_rilascio_lento.png tarato "1D, tarato senza rilascio lento" tarato_lento "1D, tarato con rilascio lento in 32 h"
julia --project=. calibrazione_scambio/confronto.jl 11_colmo_rilascio_lento tarato_lento "1D, tabella ΔQ, valori tarati con rilascio lento in 32 h"
julia --project=. calibrazione_scambio/anticipo.jl tarato_lento
julia --project=. -O2 calibrazione_scambio/attivo.jl 0.5 150000 1 32 0
julia --project=. -O2 calibrazione_scambio/attivo.jl 0.5,1,2,4,8 150000,200000,250000 1 32 1
julia --project=. -O2 calibrazione_scambio/attivo.jl 0.25,0.5 125000,175000 1 32 1
julia --project=. -O2 calibrazione_scambio/attivo.jl 0.25 150000 1 32 1
julia --project=. calibrazione_scambio/riepilogo_attivo.jl
julia --project=. -O2 calibrazione_scambio/simula_fittizi.jl attivo tabella 150000 0.5 1 32 si
julia --project=. calibrazione_scambio/idrogrammi.jl 12_idrogrammi_attivo.png tarato_lento "1D, tarato con rilascio lento (senza serbatoio attivo)" attivo "1D, tarato con serbatoio attivo"
julia --project=. calibrazione_scambio/confronto.jl 13_colmo_attivo attivo "1D, tabella ΔQ, valori tarati con serbatoio attivo"
julia --project=. calibrazione_scambio/anticipo.jl attivo
```

`taratura.jl` e `rilascio_lento.jl` aggiungono righe ai loro file senza controllare i doppioni: per rifare tutto
da capo, cancellare prima i file. Le tarature sono state fatte con prototipi dello stesso codice, con i parametri
passati da variabili d'ambiente: con `Source/` attuale i risultati sono identici byte per byte senza rilascio
lento (controllato su quattro casi) e uguali entro l'ultima cifra stampata (10⁻⁶) con il rilascio lento. Il serbatoio attivo è stato provato con una copia
di `Source/`: con quello attuale (riga 33) i risultati dei quattro fittizi sono identici byte per byte, e con la riga 33
a `no` sono identici a quelli di prima.
