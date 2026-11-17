# Riepilogo della prova del serbatoio attivo: risultati/riepilogo_attivo.txt da risultati/attivo.csv.
# Uso, dalla cartella del progetto: julia --project=. calibrazione_scambio/riepilogo_attivo.jl
using Printf, Statistics
const D = @__DIR__
L = readlines(joinpath(D, "risultati", "attivo.csv")); H = split(L[1], ',')
R = [Dict(zip(H, parse.(Float64, split(l, ',')))) for l in L[2:end] if !isempty(strip(l))]
ch(r) = (r["attivo"], r["T_h"], r["V_soglia_m3"])
combo = unique(ch.(R))
righe(c) = [r for r in R if ch(r) == c]
uno(c) = [r for r in righe(c) if r["onda"] == 1]
m(c, n) = mean(r[n] for r in uno(c))
am(c, n) = (v = [abs(r[n]) for r in righe(c) if !isnan(r[n])]; isempty(v) ? NaN : mean(v))
val(c, k, o, n) = (i = findfirst(r -> ch(r) == c && r["scenario"] == k && r["onda"] == o, R); i === nothing ? NaN : R[i][n])
io = IOBuffer()
println(io, "Serbatoio attivo (la parte oltre la soglia rientra comunque al ritmo (V − soglia)/T) contro la versione attuale.")
println(io, "Ritardo 1 h e rilascio lento 32 h fissi, tabelle S_50m_q2, fittizi 1-4. Anticipi in minuti (positivo: 1D prima del 2D):")
println(io, "rientro = S4200 oltre 67 m³/s; discesa = S6000 sotto 50 m³/s dopo il colmo; medie dei valori assoluti su tutte le onde.\n")
@printf(io, "%6s %6s %8s %8s %9s %8s %9s %9s   %s\n", "attivo", "T [h]", "soglia", "NSE med", "|Δcolmo|", "livello", "|rientro|", "|discesa|",
        "colmo S6000 f1/f2/f3/f4 (2D 69,3/63,6/81,7/78,6); discesa f1/f2/f3/f4a/f4b")
for c in sort(combo; by = c -> -m(c, "NSE"))
    @printf(io, "%6s %6g %8g %8.3f %9.2f %8.3f %9.0f %9.0f   %s; %s\n", c[1] == 1 ? "sì" : "no", c[2], c[3], m(c, "NSE"), m(c, "dcolmo"),
            m(c, "mae_livello"), am(c, "anticipo_rientro_min"), am(c, "anticipo_discesa_min"),
            join([@sprintf("%.1f", val(c, k, k == 4 ? 2 : 1, "colmo1D")) for k in 1:4], "/"),
            join([@sprintf("%.0f", val(c, k, o, "anticipo_discesa_min")) for (k, o) in ((1, 1), (2, 1), (3, 1), (4, 1), (4, 2))], "/"))
end
s = String(take!(io)); write(joinpath(D, "risultati", "riepilogo_attivo.txt"), s); print(s)
