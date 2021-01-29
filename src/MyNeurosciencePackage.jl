# MyNeurosciencePackage.jl
module MyNeurosciencePackage

include(".\\BasicFunctions.jl")
using .BasicFunctions

include(".\\NeuronalStructures.jl")
using .NeuronalStructures


include(".\\DynamicsFunction.jl")
using .DynamicsFunction

include(".\\NetworkConstruction.jl")
using .NetworkConstruction

include(".\\Plasticity.jl")
using .Plasticity

# test
end