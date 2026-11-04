using Plots, Plots.Measures, Printf, CSV, DataFrames, Infiltrator

const InputFolderParent = "input"
const InputFolderSections = "sections_split_S0200_S6000"
const OutputFolder = "output"

function slope!(slope, y, x) 
    @. slope[begin:end-1] = (y[begin+1:end] - y[begin:end-1]) / (x[begin+1:end] - x[begin:end-1]) 
    slope[end] = slope[end-1]
end

@views function main()
    # leggi il file riguardante tutte le sezioni
    InputFile = "sections_dx.csv"
    InputPath = joinpath(InputFolderParent, InputFolderSections, InputFile)
    sections_all = CSV.read(InputPath, DataFrame)
    size_sections_all = size(sections_all)
    nrow    = size_sections_all[1]
    dxs     = Array{Float64}(undef, nrow-1)
    iWSE_vect    = Array{Float64}(undef, nrow-1)
    dxs     = sections_all[!, "dist_next"]


    InputFiles = ["scalini_ns_split_S4000"]
    Format = ".csv"

    # Predispongo alcune cose per il ciclo sulle sezioni
    
    # leggi i dati
    InputPath = joinpath(InputFolderParent, InputFolderSections, InputFiles[begin]*Format)
    data = CSV.read(InputPath, DataFrame)
    
    # "ripara" i primi valori artificialmente messi a 0
    data[begin, "Reference elevation [m]"] = data[2, "Reference elevation [m]"]
    data[begin, "Mean bottom elevation [m]"] = data[2, "Reference elevation [m]"]
    data[begin, "Mean wse [m]"] = data[2, "Reference elevation [m]"]
    
    # Alloca i vettori delle dimensioni corrette:
    size_data = size(data)
    nrow, ncol = size_data[1], size_data[2]
    WSE     = Array{Float64}(undef, nrow)
    # WSE_s   = Array{Float64}(undef, nrow)
    Ω       = Array{Float64}(undef, nrow)
    Ω_s     = Array{Float64}(undef, nrow)
    W       = Array{Float64}(undef, nrow)
    W_s     = Array{Float64}(undef, nrow)
    RH      = Array{Float64}(undef, nrow)
    RH_s    = Array{Float64}(undef, nrow)
    Ks      = Array{Float64}(undef, nrow)
    jden    = Array{Float64}(undef, nrow)
    jden_s  = Array{Float64}(undef, nrow)

    # compila la water surface elevation
    WSE                 = data[!, "Mean wse [m]"]
    # WSE_s[begin:end-1]  = 

    # compila l'area
    Ω   = data[!, "Wetted area [m2]"]
    slope!(Ω_s, Ω, WSE)
    # display(plot(WSE, Ω, label="Ω"))
    plt1 = plot(Ω, WSE; 
        title="Area bagnata",
        xlabel="Area [m^2]",
        ylabel="Livello idrometrico [m]",
        legend=false,
        size=(450, 600)
    )
    display(plt1)
    savefig(plt1, joinpath(OutputFolder, "area_bagnata.svg"))

    # compila la larghezza
    W = data[!, "Wetted geometric length [m]"]
    slope!(W_s, W, WSE)
    # display(plot(WSE, W, label="W"))
    plt2 = plot(W, WSE; 
        title="Larghezza bagnata",
        xlabel="Larghezza [m]",
        ylabel="Livello idrometrico [m]",
        legend=false,
        size=(450, 600)
    )
    display(plt2)
    savefig(plt2, joinpath(OutputFolder, "larghezza_bagnata.svg"))

    # compila il raggio idraulico
    RH = Ω ./ W
    RH[begin] = 0.0
    slope!(RH_s, RH, WSE)
    # display(plot(WSE, RH, label="RH"))
    plt3 = plot(RH, WSE; 
        title="Raggio idraulico",
        xlabel="Raggio [m]",
        ylabel="Livello idrometrico [m]",
        legend=false,
        size=(450, 600)
    )
    display(plt3)
    savefig(plt3, joinpath(OutputFolder, "raggio.idraulico.svg"))

    # Calcola il Ks equivalente dell'intera sezione
    Q = data[!, "Discharge [m3/s]"]
    iWSE = 0.00246
    Ks = @. sqrt( Q / ((Ω * cbrt(RH)^2 )^2 - iWSE) )
    display(plot(WSE, Ks, label="Ks"))

    #= CALCOLO DI KS
        Non è banale. la pendenza della superficie libera va calcolata per ogni valore di 
        WSE. Per semplicità intanto assumo che sia unica per tutte le portate.
        Il Ks comunque deve cambiare a seconda delle diverse regioni bagnate. 
        Valori normali di Ks indicano un approccio fisicamente consistente.
    =#

    jden = @. 1.0/(Ω * Ks * cbrt(RH)^2)^2

    slope!(jden_s, jden, WSE)
    display(plot(WSE, jden, label="jden"))

    @infiltrate false

    return nothing
end

main()