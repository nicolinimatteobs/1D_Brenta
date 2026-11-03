module RUN_SIMULATION_RIEMPIMENTO
    using Printf, Parameters, UnPack, Infiltrator, LinearAlgebra, BenchmarkTools, Base.Threads
    Threads.nthreads()

    include("subroutines.jl")
    #include("module.jl")
    #include("struct.jl")
    include("Semi_Implicit_Riempimento.jl")
    include("writing.jl")

    export SEMI_IMPLICIT_RIEMPIMENTO!
end