# Confronto tra le portate del modello 1D e il riferimento 2D reference/piena_ns_discharge.csv
# (colonne Q_S0200 ... Q_S6000, una ogni 200 m, dati ogni 15 min; a t = 0 vale zero ed è escluso).
# Idrogramma nell'ultima sezione e colmo in tutte le sezioni; le colonne sono associate alle
# sezioni arrotondando l'ascissa ai 200 m. Grafici e numeri in <cartella degli output>/grafici.
# Uso, dalla cartella del progetto:  julia --project=. strumenti/confronto_reference.jl
using Plots, Printf
const R   = dirname(@__DIR__)                                    # cartella del progetto
const SET = [strip(l) for l in readlines(joinpath(R, "settings.txt"))[1:27]]
const GEO = joinpath(R, SET[2])                                  # geometria (riga 2)
const D   = joinpath(R, SET[23])                                 # cartella degli output (riga 23)
const REF = joinpath(R, "reference", "piena_ns_discharge.csv")
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

interp(tt, v, s) = (i = clamp(searchsortedlast(tt, s), 1, length(tt)-1); v[i] + (v[i+1]-v[i])*(s - tt[i])/(tt[i+1]-tt[i]))
ok = tr .>= t[1]
vol(tt, v) = sum(0.5 .* (v[1:end-1] .+ v[2:end]) .* diff(tt))
function confronto(i, nome)
    qr = Rf[ok, col(nome)]; qm = [interp(t, Qs[:, i], s) for s in tr[ok]]
    e  = qm .- qr
    kr = argmax(qr); km = argmax(Qs[:, i])
    nse = 1 - sum(e.^2) / sum((qr .- sum(qr)/length(qr)).^2)
    return (; qr, qm, e, pr = qr[kr], tpr = tr[ok][kr], pm = Qs[km, i], tpm = t[km],
            rmse = sqrt(sum(e.^2)/length(e)), emax = maximum(abs.(e)), nse,
            vr = vol(tr[ok], qr), vm = vol(tr[ok], qm))
end
nomecol(xi) = @sprintf("Q_S%04d [m3/s]", round(Int, xi/200)*200)

nome_u = nomecol(x[end]); c = confronto(length(x), nome_u)
p1 = plot(; title = "Sezione $(nome_u[3:7]): portata", ylabel = "portata [m³/s]", legend = :topright, st...)
plot!(p1, tr[ok] ./ 3600, c.qr; label = "riferimento 2D (piena_ns_discharge)", color = :black, linewidth = 2, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
plot!(p1, t ./ 3600, Qs[:, end]; label = "modello 1D", color = BLU, linewidth = 2)
plot!(p1, tr[ok] ./ 3600, Rf[ok, col("Q_BC_monte [m3/s]")]; label = "ingresso da monte", color = GRI, linestyle = :dash, linewidth = 1)
p2 = plot(; title = "Differenza modello - riferimento", xlabel = "tempo [h]", ylabel = "ΔQ [m³/s]", legend = false, st...)
hline!(p2, [0.0]; color = :black, linewidth = 0.8)
plot!(p2, tr[ok] ./ 3600, c.e; color = ARA, linewidth = 1.8, markershape = :circle, markersize = 2.5, markerstrokewidth = 0)
savefig(plot(p1, p2; layout = grid(2, 1, heights = [0.65, 0.35]), size = (1000, 750), dpi = 130, left_margin = 6Plots.mm),
        joinpath(OUT, "8_confronto_S6000.png"))

pr = Float64[]; pm = Float64[]; tpr = Float64[]; tpm = Float64[]
for i in eachindex(x)
    n = nomecol(x[i]); col(n) === nothing && error("colonna mancante nel riferimento: $n")
    ci = confronto(i, n); push!(pr, ci.pr); push!(pm, ci.pm); push!(tpr, ci.tpr); push!(tpm, ci.tpm)
end
q1 = plot(; title = "Portata al colmo lungo il tratto", ylabel = "portata [m³/s]", legend = :topright, st...)
plot!(q1, x ./ 1000, pr; label = "riferimento 2D", color = :black, linewidth = 2, markershape = :circle, markersize = 3, markerstrokewidth = 0)
plot!(q1, x ./ 1000, pm; label = "modello 1D", color = BLU, linewidth = 2, markershape = :circle, markersize = 3, markerstrokewidth = 0)
q2 = plot(; title = "Istante del colmo", xlabel = "distanza [km]", ylabel = "tempo [h]", legend = :topleft, st...)
plot!(q2, x ./ 1000, tpr ./ 3600; label = "riferimento 2D (dati ogni 15 min)", color = :black, linewidth = 2, markershape = :circle, markersize = 3, markerstrokewidth = 0)
plot!(q2, x ./ 1000, tpm ./ 3600; label = "modello 1D", color = BLU, linewidth = 2, markershape = :circle, markersize = 3, markerstrokewidth = 0)
savefig(plot(q1, q2; layout = (2, 1), size = (1000, 700), dpi = 130, left_margin = 6Plots.mm), joinpath(OUT, "9_colmo_lungo_il_tratto.png"))

io = IOBuffer()
@printf(io, "Confronto %s (modello 1D - riferimento piena_ns_discharge.csv), dati ogni 15 min da t = %.0f s\n", nome_u[3:7], tr[ok][1])
@printf(io, "colmo riferimento %.2f m3/s a %.2f h | modello %.2f m3/s a %.2f h | differenza %+.2f m3/s (%+.2f%%), %+.0f min\n",
        c.pr, c.tpr/3600, c.pm, c.tpm/3600, c.pm - c.pr, 100(c.pm - c.pr)/c.pr, (c.tpm - c.tpr)/60)
@printf(io, "RMSE %.3f m3/s, scarto massimo %.3f m3/s, NSE %.4f\n", c.rmse, c.emax, c.nse)
@printf(io, "volume transitato: riferimento %.0f m3, modello %.0f m3 (%+.2f%%)\n", c.vr, c.vm, 100(c.vm - c.vr)/c.vr)
@printf(io, "colmo lungo il tratto: differenza massima %.2f m3/s; ritardo del colmo dalla prima all'ultima sezione: riferimento %.0f min, modello %.0f min\n",
        maximum(abs.(pm .- pr)), (tpr[end]-tpr[1])/60, (tpm[end]-tpm[1])/60)
s = String(take!(io)); write(joinpath(OUT, "confronto_S6000.txt"), s); println(s)
