# Calibrazione dello scambio laterale con la piana (soglia, T, ritardo)

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
- Parametri del serbatoio (righe 29-31 di `settings.txt`):

  | Riga | Parametro | Significato |
  | 29 | soglia [m³] | la piana trattiene i primi m³ usciti (si riempie): rientra solo il volume oltre la soglia |
  | 30 | T [h] | tempo di svuotamento: rientro al massimo (V − soglia)/T; 0 = tutto subito |
  | 31 | τ [h] | ritardo di percorrenza: l'acqua uscita arriva al serbatoio dopo τ |

## Risultati

Valori tarati: **soglia 150 000 m³, T = 0,5 h, τ = 1 h** (NSE medio in S6000 massimo). Medie sui fittizi 1-4:

| Simulazione | NSE in S6000 | errore del colmo lungo il tratto | errore del livello |
|---|---|---|---|
| sfioro sugli argini | 0,937 | – | – |
| tabella, serbatoio senza parametri (T = 0) | 0,926 | 4,18 m³/s | 5,9 cm |
| tabella, solo T = 4 h | 0,986 | 1,86 m³/s | 3,6 cm |
| tabella, soglia 150 000 m³, T = 1 h | 0,989 | 1,11 m³/s | 3,3 cm |
| **tabella, soglia 150 000 m³, T = 0,5 h, τ = 1 h** | **0,991** | **1,02 m³/s** | **3,1 cm** |

- Colmo in S6000 (2D → 1D tarato): 69,3 → 70,1; 63,6 → 63,2; 81,7 → 81,9; 78,6 → 78,8 m³/s (`fittizio_4`:
  seconda onda). Con lo sfioro era sempre circa 60.
- **Senza soglia** (grafici 1-3) tutta l'acqua uscita rientra subito: il colmo in uscita resta quello in ingresso.
  Con il solo T il rientro è graduale e troppo presto, mentre nel 2D arriva a scatto, quando la piana si è riempita.
- Nel 2D il volume che resta fuori (uscito − rientrato, per tratto) è simile nei quattro scenari, 154-215 mila m³:
  da qui l'idea della soglia.
- **Il ritardo** cambia poco gli indici medi (le combinazioni migliori sono tutte tra 0,989 e 0,991) ma porta il
  rientro al momento giusto: in S4200 la portata supera 67 m³/s 44-50 min prima del 2D senza ritardo, 1-7 min
  dopo con τ = 1 h (`risultati/anticipo_*.txt`).
- Il minimo è largo: T da 0,5 a 2 h, soglia 100-150 mila m³ e τ 0,5-1,5 h danno NSE medi tra 0,989 e 0,991.

## Limiti

- **Seconda onda del `fittizio_4`**: il rientro arriva ancora circa 110 min prima. Il serbatoio non si svuota
  mai sotto la soglia, quindi alla seconda onda la piana risulta già piena; nel 2D fra le due onde si svuota in
  parte. Servirebbe uno svuotamento lento della parte sotto la soglia.
- Nel ramo calante il 1D scende ancora 0,5-1 h prima del 2D.
- Nel `fittizio_2` escono 136 000 m³, meno della soglia, e non rientra niente (nel 2D il 25%); il colmo in uscita
  è comunque giusto. In generale il 1D fa uscire il 10-30% di acqua in meno del 2D.
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
| `risultati/taratura.csv` | Una riga per combinazione e scenario: 73 combinazioni (11 di solo T, 25 di soglia e T, 36 con il ritardo, più 1 prova) |
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
```

`taratura.jl` aggiunge righe a `risultati/taratura.csv` senza controllare i doppioni: per rifare tutto da capo,
cancellare prima il file. Le tarature sono state fatte con un prototipo dello stesso codice, con i parametri
passati da variabili d'ambiente: con `Source/` attuale i risultati sono identici byte per byte (controllato su
due combinazioni).
