include("Source/Run_Simulation_Riempimento.jl")
#include("Source/Run_Simulation.jl")

using .RUN_SIMULATION_RIEMPIMENTO
#using .RUN_SIMULATION

RUN_SIMULATION_RIEMPIMENTO.SEMI_IMPLICIT_RIEMPIMENTO!()
#RUN_SIMULATION.SEMI_IMPLICIT!()


