# MyNeurosciencePackage.jl
module MyNeurosciencePackage

include(".\\BasicFunctions.jl")
using .BasicFunctions
export heaviside,rect_linear

include(".\\NeuronalStructures.jl")
using .NeuronalStructures
export simulation_parameters

include(".\\DynamicsFunction.jl")
using .DynamicsFunction
export simulation_step!,time_step


include(".\\NetworkConstruction.jl")
using .NetworkConstruction
export create_network

#include(".\\Plasticity.jl")
#using .Plasticity


end