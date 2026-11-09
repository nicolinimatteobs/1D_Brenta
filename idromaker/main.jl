using Printf, Plots, Plots.Measures, Infiltrator, DataFrames, CSV

gumbel(t, μ, β) = @. (1/β) * exp( -(t - μ)/β - exp(-(t - μ)/β) )

# Gumbel con il picco mantenuto per una durata T: la salita è la stessa, poi la portata resta
# al colmo da μ a μ+T e la discesa è quella della Gumbel spostata in avanti di T
gumbel_plateau(t, μ, β, T) = gumbel(@.(min(t, max(μ, t - T))), μ, β)

const ParentDir = normpath(joinpath(@__DIR__, ".."))

function main(; dosave=false)
    centa = CSV.read("idromaker/Centa.csv", DataFrame)

    plt1 = plot(
        (centa.Time[11:31].-36000.0)./3600, centa.Discharge[11:31];
        legend=true,
        label="Centa",
        xlabel="tempo [h]",
        ylabel="portata [m³/s]",
        # xlim=(3.5e4, 1.1e5), # (11, 31)
    ); display(plt1)


    # Parametri della piena di riferimento (scenario 1)
    dt       = 900.0            # s, passo degli idrogrammi
    durata   = 86400.0          # s
    μ        = 25200
    β        = 10800
    α        = 2300000
    baseflow = 5.0              # m^3/s

    # Scenario 2: come l'1 con portata più bassa. Il colmo è di circa 75 m³/s: nel modello 1D
    # l'esondazione inizia intorno a 66 m³/s (S1000)
    α_bassa  = 2050000

    # Scenario 3: come l'1 con il picco mantenuto per T_picco. La durata aumenta di T_picco
    T_picco  = 6*3600.0         # s

    # Scenario 4: due picchi distanziati di Δ_picchi. La durata aumenta di Δ_picchi.
    # Il secondo picco si somma alla coda del primo: α2 è ridotto perché il colmo resti uguale
    Δ_picchi = 15*3600.0        # s
    μ2       = μ + Δ_picchi
    α2       = α * (1 - gumbel(μ2, μ, β)/gumbel(μ2, μ2, β))

    t1 = collect(0.0:dt:durata)
    t3 = collect(0.0:dt:durata+T_picco)
    t4 = collect(0.0:dt:durata+Δ_picchi)

    scenari = [
        ("fittizio_1", t1, baseflow .+ α * gumbel(t1, μ, β)),
        ("fittizio_2", t1, baseflow .+ α_bassa * gumbel(t1, μ, β)),
        ("fittizio_3", t3, baseflow .+ α * gumbel_plateau(t3, μ, β, T_picco)),
        ("fittizio_4", t4, baseflow .+ α * gumbel(t4, μ, β) .+ α2 * gumbel(t4, μ2, β)),
    ]

    # percorso output (Da "approccio_2D_1D")
    test    = "input"
    name    = "piena.csv"

    for (evento, time, discharge) in scenari
        plot!(plt1, time./3600, discharge; label=evento)

        OutputPath = joinpath(ParentDir, test, "idrogrammi", evento, name)
        output = DataFrame(time=time, discharge=round.(discharge, digits=2))

        if dosave
            mkpath(dirname(OutputPath))
            CSV.write(OutputPath, output; header=false)
        end
    end
    display(plt1)

    return nothing
end

main(dosave=true)
