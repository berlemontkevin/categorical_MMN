# MyNeurosciencePackage.jl
module MyNeurosciencePackage

include(".\\BasicFunctions.jl")
using .BasicFunctions
export heaviside,rect_linear, sigmoid

include(".\\NeuronalStructures.jl")
using .NeuronalStructures
export simulation_parameters, dendrites_param_sigmoid, dend_sigmoid, soma_PC
export neural_integrator, microcircuit
export get_dynamics, save_dynamics
export vip_cell, pv_cell, sst_cell
export gaba_syn, ampa_syn, nmda_syn


include(".\\DynamicsFunction.jl")
using .DynamicsFunction
export simulation_step!,time_step
export dendrite_input_output!
export sum_input!, current_to_frequency!, update_firing!, update_dend!
export full_time_dynamics

include(".\\NetworkConstruction.jl")
using .NetworkConstruction
export create_network

#include(".\\Plasticity.jl")
#using .Plasticity


end