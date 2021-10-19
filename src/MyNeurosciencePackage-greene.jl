# MyNeurosciencePackage-greene.jl
module MyNeurosciencePackage-greene

include("./BasicFunctions.jl")
using .BasicFunctions
export heaviside,rect_linear, sigmoid, create_process!, OU_process, update_process!

include("./NeuronalStructures.jl")
using .NeuronalStructures
export simulation_parameters, dendrites_param_sigmoid, dend_sigmoid, soma_PC
export neural_integrator, microcircuit
export get_dynamics, save_dynamics
export vip_cell, pv_cell, sst_cell
export gaba_syn, ampa_syn, nmda_syn
export bump_structure, layer_bump

include("./DynamicsFunction.jl")
using .DynamicsFunction
export simulation_step!,time_step
export dendrite_input_output!
export sum_input!, current_to_frequency, update_firing!, update_dend!
export full_time_dynamics, current_synapses!, synapse_derivative, update_adaptation!
export update_syn!


include("./NetworkConstruction.jl")
using .NetworkConstruction
export create_network, construct_two_local_microcircuit_integrator
export create_deterministic_oddball, create_deterministic_oscillations_oddball
export construct_two_local_microcircuit_integrator_full_param
export orientation_kernel, create_layer_bump

include("./PlotFunctions.jl")
using .PlotFunctions
export plot_local_circuit, plot_local_circuit_synapses


include("./DataAnalysis.jl")
using .DataAnalysis
export get_mean_firing_rate, compute_MMN_oddball, compute_MMN_time

#include(".\\Plasticity.jl")
#using .Plasticity


end