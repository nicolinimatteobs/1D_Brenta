# Confronto tra le portate del modello 1D e un riferimento 2D in reference/ (colonne Q_S<ascissa>,
# dati ogni 15 min; a t = 0 vale zero ed è escluso). Idrogramma nell'ultima sezione e colmo in tutte
# le sezioni; le colonne sono associate alle sezioni arrotondando l'ascissa ai 50 m.
# Tutti i confronti sono fatti sul periodo comune: se il riferimento finisce prima della simulazione,
# il colmo si cerca solo fino alla sua fine e, nelle sezioni in cui il riferimento non ha ancora
# raggiunto il colmo (massimo all'ultimo istante disponibile), il colmo non viene confrontato.
# Nel grafico del colmo lungo il tratto sono segnate con il nome le sezioni del 1D che sfiorano in un tratto
# (predefinito tra 2,5 e 4,0 km): il volume uscito da ogni argine è ricavato da eta.csv con la formula del
# modello (sfioro libero, Cd = 0,385, lungo tutta la cella) ed è riportato anche nel file di testo. Serve
# una geometria con le quote degli argini (6 colonne), altrimenti non viene segnato niente.
# Grafici e numeri in <cartella degli output>/grafici.
# Uso, dalla cartella del progetto:
#   julia --project=. strumenti/confronto_reference.jl [cartella degli output] [file di riferimento] [km iniziale] [km finale]
# Senza argomenti usa la cartella della riga 22 di settings.txt e reference/ns_channel_piena_discharge.csv.
using Plots, Printf
const R   = dirname(@__DIR__)                                    # cartella del progetto
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:26]]
const GEO = joinpath(R, SET[2])                                  # geometria (riga 2)
const D   = joinpath(R, length(ARGS) >= 1 ? ARGS[1] : SET[22])   # cartella degli output
const REF = joinpath(R, length(ARGS) >= 2 ? ARGS[2] : joinpath("reference", "ns_channel_piena_discharge.csv"))
const TRATTO = length(ARGS) >= 4 ? (parse(Float64, ARGS[3]), parse(Float64, ARGS[4])) : (2.5, 4.0)   # tratto con i nomi [km]
const OUT = joinpath(D, "grafici"); mkpath(OUT)
const BLU = "#2a78d6"; const ARA = "#eb6834"; const GRI = "#7a7a7a"
st = (framestyle = :box, gridalpha = 0.25, titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8)

mat(f) = reduce(vcat, [permutedims(parse.(Float64, split(l))) for l in eachline(joinpath(D, f)) if !isempty(strip(l))])
t  = [parse(Float64, strip(l)) for l in eachline(joinpath(D, "time.csv")) if !isempty(strip(l))]
Q  = mat("Q.csv"); Qs = 0.5 .* (Q[:, 1:end-1] .+ Q[:, 2:end])          # portata nelle sezioni
G  = [parse.(Float64, split(l, ',')) for l in readlines(GEO)[2:end] if !isempty(strip(l))]
x  = getindex.(G, 3)

L   = readlines(REF)
hdr = strip.(split(L[1], ','))
Rf  = reduce(vcat, [permutedims(parse.(Float64, split(l, ','))) for l in L[2:end] if !isempty(strip(l))])
tr  = Rf[:, 1]
col(n) = findfirst(==(n), hdr)
nref = splitext(basename(REF))[1]

interp(tt, v, s) = (i = clamp(searchsortedlast(tt, s), 1, length(tt)-1); v[i] + (v[i+1]-v[i])*(s - tt[i])/(tt[i+1]-tt[i]))
ok   = (tr .>= t[1]) .& (tr .<= t[end])                          # istanti del riferimento nel periodo comune
tfin = min(tr[end], t[end])                                      # fine del periodo comune
parz = tr[end] < t[end] - 1.0                                    # riferimento più corto della simulazione
km_ok = t .<= tfin                                               # istanti del modello nel periodo comune
vol(tt, v) = sum(0.5 .* (v[1:end-1] .+ v[2:end]) .* diff(tt))
function confronto(i, nome)
    # istanti validi: periodo comune e riferimento già bagnato (nei 2D che partono da alveo asciutto
    # la portata vale 0 finché non arriva il fronte, e il 1D parte invece con la portata di base)
    m  = ok .& (Rf[:, col(nome)] .> 0)
    count(m) < 2 && return (; valido = false, tt = tr[m], qr = Rf[m, col(nome)], colmo_ok = false)
    tt = tr[m]; qr = Rf[m, col(nome)]; qm = [interp(t, Qs[:, i], s) for s in tt]
    e  = qm .- qr
    kr = argmax(qr)
    km = argmax(Qs[km_ok, i])
    nse = 1 - sum(e.^2) / sum((qr .- sum(qr)/length(qr)).^2)
    return (; valido = true, tt, qr, qm, e, pr = qr[kr], tpr = tt[kr], pm = Qs[km, i], tpm = t[km],
            colmo_ok = kr < length(qr),                          # colmo del riferimento raggiunto prima della sua fine
            rmse = sqrt(sum(e.^2)/length(e)), emax = maximum(abs.(e)), nse,
            vr = vol(tt, qr), vm = vol(tt, qm))
end
nomecol(xi) = @sprintf("Q_S%04d [m3/s]", round(Int, xi/50)*50)

nome_u = nomecol(x[end]); c = confronto(length(x), nome_u)
p1 = plot(; title = "Sezione $(nome_u[3:7]): portata", ylabel = "portata [m³/s]", legend = :topright, st...)
c.valido && plot!(p1, c.tt ./ 3600, c.qr; label = "riferimento 2D ($nref)", color = :black, linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
plot!(p1, t ./ 3600, Qs[:, end]; label = "modello 1D", color = BLU, linewidth = 2)
plot!(p1, tr[ok] ./ 3600, Rf[ok, col("Q_BC_monte [m3/s]")]; label = "ingresso da monte (2D)", color = GRI, linestyle = :dash, linewidth = 1)
parz && vline!(p1, [tfin/3600]; color = ARA, linestyle = :dot, linewidth = 1.2, label = @sprintf("fine del riferimento (%.2f h)", tfin/3600))
p2 = plot(; title = "Differenza modello - riferimento", xlabel = "tempo [h]", ylabel = "ΔQ [m³/s]", legend = false,
          xlims = (0, t[end]/3600), st...)
hline!(p2, [0.0]; color = :black, linewidth = 0.8)
c.valido && plot!(p2, c.tt ./ 3600, c.e; color = ARA, linewidth = 1.8, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
c.valido || annotate!(p2, t[end]/7200, 0.0, text("il riferimento è asciutto in questa sezione per tutto il periodo disponibile", 9, :center, :bottom, GRI))
savefig(plot(p1, p2; layout = grid(2, 1, heights = [0.65, 0.35]), size = (1000, 750), dpi = 130, left_margin = 6Plots.mm),
        joinpath(OUT, "7_confronto_S6000.png"))

pr = Float64[]; pm = Float64[]; tpr = Float64[]; tpm = Float64[]; cok = Bool[]
for i in eachindex(x)
    n = nomecol(x[i]); col(n) === nothing && error("colonna mancante nel riferimento: $n")
    ci = confronto(i, n)
    if ci.valido
        push!(pr, ci.pr); push!(pm, ci.pm); push!(tpr, ci.tpr); push!(tpm, ci.tpm)
    else
        km = argmax(Qs[km_ok, i]); push!(pr, NaN); push!(pm, Qs[km, i]); push!(tpr, NaN); push!(tpm, t[km])
    end
    push!(cok, ci.colmo_ok)
end
prv = ifelse.(cok, pr, NaN); tprv = ifelse.(cok, tpr, NaN)       # colmi del riferimento non raggiunti: esclusi

# Sfioro del 1D sezione per sezione: volume uscito dall'argine destro e dal sinistro in tutta la simulazione
Vdx = zeros(length(x)); Vsx = zeros(length(x))
if all(r -> length(r) >= 6, G)
    E1  = mat("eta.csv"); dxb = diff(x)
    dxc = [k == 1 ? dxb[1] : k == length(x) ? dxb[end] : 0.5*(dxb[k] + dxb[k-1]) for k in eachindex(x)]   # lunghezza della cella
    qsf(η, z) = (Y = max(η - z, 0.0); 0.385*sqrt(2*9.81*Y^3))   # sfioro libero per unità di lunghezza [m²/s]
    Vdx .= [vol(t, qsf.(E1[:, k], G[k][5]) .* dxc[k]) for k in eachindex(x)]
    Vsx .= [vol(t, qsf.(E1[:, k], G[k][6]) .* dxc[k]) for k in eachindex(x)]
end
sfiora = findall(k -> TRATTO[1] <= x[k]/1000 <= TRATTO[2] && Vdx[k] + Vsx[k] > 1.0, eachindex(x))
nomesez(xi) = @sprintf("S%04d", round(Int, xi/50)*50)
lato(k) = Vdx[k] > 1.0 && Vsx[k] > 1.0 ? "due argini" : Vdx[k] > 1.0 ? "argine destro" : "argine sinistro"
km_it(v) = replace(@sprintf("%g", v), "." => ",")
nota = parz ? @sprintf(" (fino a %.2f h, fine del riferimento)", tfin/3600) : ""
q1 = plot(; title = "Portata al colmo lungo il tratto" * nota, ylabel = "portata [m³/s]", legend = isempty(sfiora) ? :topright : :bottomleft, st...)
plot!(q1, x ./ 1000, prv; label = "riferimento 2D", color = :black, linewidth = 2, markershape = :circle, markersize = 3, markerstrokewidth = 0)
plot!(q1, x ./ 1000, pm; label = "modello 1D", color = BLU, linewidth = 2, markershape = :circle, markersize = 3, markerstrokewidth = 0)
if !isempty(sfiora)
    # nomi delle sezioni del 1D che sfiorano nel tratto indicato, a destra del punto
    scatter!(q1, x[sfiora] ./ 1000, pm[sfiora]; label = "1D: sezioni con sfioro tra $(km_it(TRATTO[1])) e $(km_it(TRATTO[2])) km",
             markershape = :diamond, markersize = 7, color = ARA, markerstrokecolor = :white, markerstrokewidth = 1)
    alto = maximum(filter(!isnan, vcat(prv, pm))) - minimum(filter(!isnan, vcat(prv, pm)))
    for k in sfiora
        annotate!(q1, x[k]/1000 + 0.05, pm[k] + 0.015alto,
                  text(@sprintf("%s (%s, %.0f mila m³ usciti)", nomesez(x[k]), lato(k), (Vdx[k] + Vsx[k])/1000), 8, :left, :bottom, :black))
    end
end
q2 = plot(; title = "Istante del colmo" * nota, xlabel = "distanza [km]", ylabel = "tempo [h]", legend = :topleft, st...)
plot!(q2, x ./ 1000, tprv ./ 3600; label = "riferimento 2D (dati ogni 15 min)", color = :black, linewidth = 2, markershape = :circle, markersize = 3, markerstrokewidth = 0)
plot!(q2, x ./ 1000, tpm ./ 3600; label = "modello 1D", color = BLU, linewidth = 2, markershape = :circle, markersize = 3, markerstrokewidth = 0)
if !all(cok)
    annotate!(q1, x[end]/1000, maximum(pm), text("sezioni senza punto nero: il riferimento non ha ancora raggiunto il colmo", 8, :right, :top, GRI))
end
savefig(plot(q1, q2; layout = (2, 1), size = (1000, 700), dpi = 130, left_margin = 6Plots.mm), joinpath(OUT, "8_colmo_lungo_il_tratto.png"))

io = IOBuffer()
@printf(io, "Confronto %s (modello 1D - riferimento %s), dati ogni 15 min\n", nome_u[3:7], basename(REF))
parz && @printf(io, "ATTENZIONE: il riferimento arriva a %.0f s (%.2f h), la simulazione a %.0f s (%.2f h): tutti i valori sono sul periodo comune\n",
                tr[end], tr[end]/3600, t[end], t[end]/3600)
if !c.valido
    @printf(io, "nessun confronto in %s: il riferimento vi è asciutto (portata nulla) per tutto il periodo disponibile\n", nome_u[3:7])
else
    c.tt[1] > tr[ok][1] + 1 && @printf(io, "il riferimento in %s è asciutto fino a %.2f h: il confronto parte da lì\n", nome_u[3:7], c.tt[1]/3600)
    if c.colmo_ok
        @printf(io, "colmo riferimento %.2f m3/s a %.2f h | modello %.2f m3/s a %.2f h | differenza %+.2f m3/s (%+.2f%%), %+.0f min\n",
                c.pr, c.tpr/3600, c.pm, c.tpm/3600, c.pm - c.pr, 100(c.pm - c.pr)/c.pr, (c.tpm - c.tpr)/60)
    else
        @printf(io, "colmo non confrontabile: il riferimento non lo ha ancora raggiunto (massimo disponibile %.2f m3/s a %.2f h; modello nel periodo comune %.2f m3/s a %.2f h)\n",
                c.pr, c.tpr/3600, c.pm, c.tpm/3600)
    end
    @printf(io, "RMSE %.3f m3/s, scarto massimo %.3f m3/s, NSE %.4f (%d istanti)\n", c.rmse, c.emax, c.nse, length(c.e))
    @printf(io, "volume transitato nel periodo confrontato: riferimento %.0f m3, modello %.0f m3 (%+.2f%%)\n", c.vr, c.vm, 100(c.vm - c.vr)/c.vr)
end
if all(cok)
    @printf(io, "colmo lungo il tratto: differenza massima %.2f m3/s; ritardo del colmo dalla prima all'ultima sezione: riferimento %.0f min, modello %.0f min\n",
            maximum(abs.(pm .- pr)), (tpr[end]-tpr[1])/60, (tpm[end]-tpm[1])/60)
else
    @printf(io, "colmo lungo il tratto: confrontabile in %d sezioni su %d (nelle altre il riferimento non lo ha ancora raggiunto)", count(cok), length(cok))
    any(cok) && @printf(io, "; differenza massima %.2f m3/s", maximum(abs.(pm[cok] .- pr[cok])))
    println(io)
end
if sum(Vdx) + sum(Vsx) > 1.0
    @printf(io, "sfioro del 1D in tutta la simulazione (da eta.csv, formula del modello): %.0f m3 da %d sezioni\n",
            sum(Vdx) + sum(Vsx), count(>(1.0), Vdx .+ Vsx))
    @printf(io, "sezioni del 1D che sfiorano tra %s e %s km: %s\n", km_it(TRATTO[1]), km_it(TRATTO[2]), isempty(sfiora) ? "nessuna" :
            join([@sprintf("%s %.0f m3 (%s)", nomesez(x[k]), Vdx[k] + Vsx[k], lato(k)) for k in sfiora], ", "))
end
s = String(take!(io)); write(joinpath(OUT, "confronto_S6000.txt"), s); println(s)
