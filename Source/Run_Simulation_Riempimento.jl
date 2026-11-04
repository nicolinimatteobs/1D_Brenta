module RUN_SIMULATION_RIEMPIMENTO
    using Printf, Parameters, UnPack, Infiltrator, LinearAlgebra, BenchmarkTools, Base.Threads
    Threads.nthreads()

    include("struct.jl")
    include("module.jl")
    include("subroutines.jl")
    include("writing.jl")
    include("Semi_Implicit_Riempimento.jl")

    export SEMI_IMPLICIT_RIEMPIMENTO!
end
