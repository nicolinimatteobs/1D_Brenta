module RUN_SIMULATION
    using Printf, Parameters, UnPack, Infiltrator, LinearAlgebra, BenchmarkTools

    include("subroutines.jl")
    #include("module.jl")
    #include("struct.jl")
    include("Semi_Implicit_Brenta_2D_1D.jl")
    include("writing.jl")

    export SEMI_IMPLICIT!
end