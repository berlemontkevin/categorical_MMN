module NeuronalStructures




##############################################################################################################
# Abstract types
###############################################################################################################
module AbstractNeuronalTypes
using Parameters
export eq_diff_method, neuron, dendrite, pyr_cell, synapses, interneuron


"""
    Type: eq_diff_method

Defines the type of method that will be used to solve the differential equations

...
# Types currently implemented
- Euler method
....
"""
abstract type eq_diff_method end

"""
    Type: neuron

Abstract type for all the type of neurons

...
# Neurons currently implemented
- dendrites
- pyramidal cells
- interneurons
....
"""
abstract type neuron end



"""
    Type: synapses
    
Abstract type for the types of synapses
"""
abstract type synapses end


"""
    Type: dendrite
    
Abstract type for the types of dendrites
"""
abstract type dendrite <: neuron end

"""
    Type: pyr_cell
    
Abstract type for the types of pyramidal cells
"""
abstract type pyr_cell <: neuron end

"""
    Type: interneuron
    
Abstract type for the types of interneuron
"""
abstract type interneuron <: neuron end

end
using .AbstractNeuronalTypes

###############################################################################################################
# Differential equation method
###############################################################################################################
module EqDiffMethod

using ..AbstractNeuronalTypes
using Parameters
export euler_method

"""
    Struct: Euler method

Defines the parameters for euler method. dt is the time step and record indicates if the simulation is going to record all the time steps or not.
"""
@with_kw struct euler_method <: eq_diff_method
    dt::Float64
    record::Bool = true
end

end
using .EqDiffMethod
###############################################################################################################
# Synapses
###############################################################################################################

module Synapses

using ..AbstractNeuronalTypes

using Parameters
export syn_connection

"""
    Struct: syn_connection

Define the structure of a synaptic connexion. It consists of a pre-synaptic neuron and a Float (or vector if history needs to be saved) representing the synaptic weight.
It will be an element of the neuron struct.
"""
struct syn_connection <: synapses
    pre_syn::neuron
    w::Vector{Float64}
end

end
using.Synapses

###############################################################################################################
# Neuronal Models
###############################################################################################################
module NeuronalModels
export rectified_linear_neurons

using ..AbstractNeuronalTypes, ..Synapses

using Parameters


"""
    Struct: rectified_linear_neurons

It consists in the model of a neuron as a rectified_linear neuron

...
# Arguments
- Isyn::Float64 = synaptic input toward the neuron
- τ:;Float64 = time constant of the neuron
- neurons_list::Vector{syn_connection} = list of all presynaptic connexions going to this neuron
- fr::Float64 = firing rate of the neuron (vector if history needs to be saved)
- Ibg::Float64 background input toward this neuron (vector if history needs to be saved)
- objective_rate::Float64 = objective rate of the neural population. If <0, it means that it shouldn't be taken into account
#TODO introduire un test automatique lors de la construction
....
"""
@with_kw struct rectified_linear_neurons <: neuron
    Isyn::Vector{Float64} = [0.0] 
    τ::Float64 = 2.0 
    neurons_list::Vector{syn_connection} = syn_connection[] 
    fr::Vector{Float64} = [0.0] 
    Ibg::Vector{Float64} = [0.0] 
    objective_rate::Float64 = -1.0
end

end

using .NeuronalModels


#######################
# Hertag & Sprekler 2020
###########################
module Hertag2020_Structures

using Parameters
using ..AbstractNeuronalTypes, ..Synapses

export soma_hertag, dendrites_hertag


"""
    Structure: Soma hertag <: pyr cell

Defines the soma following (https://elifesciences.org/articles/57541)

...
# Arguments
- TauE::Float64 = time constant of the soma en seconde
- fr::Vector{Float64} = vector of the history of firing rate of the soma
- Isyn::Vector{Float64} = vector of the history of synaptic inputs to the soma
- Ibg::Vector{Float64} = vector of the history of background input to the soma
- Itot::Vector{Float64} = vector of the total synaptic input 
- LamdaD::Float64 = leaking conductivity of the dendrites
- LambdaE::Float64 = leaking ocnductivity of the soma
- Thetar::Float64 =  rheobase of the neuron
- objective rate::Float64 = objective rate of this neuron (not used if negative)
- dendrites list::Vector{dendrite} =  list of the dendrites connected to this pyr_cell
- neurons list::Vector{syn connection} = vector of presynaptic connections to this soma
...
"""
@with_kw struct soma_hertag  <: pyr_cell
    τE::Float64 = 0.060 
    fr::Vector{Float64} = [0.0] 
    Isyn::Vector{Float64} = [0.0]
    Ibg::Vector{Float64} = [0.0] 
    Itot::Vector{Float64} = [0.0] 
    λD::Float64 = 0.27
    ID::Vector{Float64} = [0.0]
    λE::Float64 = 0.31
    Θr::Float64 = 14.0 
    IE::Vector{Float64} = [0.0]


    objective_rate::Float64 = -1.0 

    
    dendrites_list::Vector{dendrite} = dendrite[] 
    neurons_list::Vector{syn_connection} = syn_connection[]


end


"""
    Struct: dendriteshertag

...
# Arguments
- neurons list::Vector{syn connection} = list of the presynaptic neurons connected to the dendrite
- Thetac::Float64 =  minimal input to produce a calcium spike
- Itot::Vector{Float64} =  total synaptically generated input in the dendrites
- gain calcium::Float64= calcium gain of the calcium spike
- Isyn::Vector{Float64} = total synaptic input to the dendrite
- pyr target::Vector{soma hertag} =  target of the dendrites: pyramidal cells
- Ibg::Vector{Float64} = background input to the dendrite
- calcium spike::Vector{Float64} =  input resulting from calclium spike
...
"""
@with_kw  struct dendrites_hertag <: dendrite
    neurons_list::Vector{syn_connection} = syn_connection[]
    Θc::Float64 = 28.0 
    Itot::Vector{Float64} = [0.0]
    gain_calcium::Float64 = 7 # Hz
    Isyn::Vector{Float64} = [0.0]
    pyr_target::Vector{soma_hertag} = soma_hertag[] 
    Ibg::Vector{Float64} = [0.0]
    calcium_spike::Vector{Float64} = [0.0]

end

end

using .Hertag2020_Structures



module attractor_network
using Parameters
using ..AbstractNeuronalTypes, ..EqDiffMethod

export wong_wang_cell

@with_kw struct wong_wang_cell <: neuron
    r::Array{Float64} = [0.0]
    s::Array{Float64} = [0.0] # synaptic gatign variable
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    Ioutput::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Array{Float64} = [0.3255]
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Array{Float64} = [0.0]
    γ::Float64 = 0.641
    a::Float64 = 270.0
    b::Float64 = 108.0
    d::Float64 = 0.154

end



end
using .attractor_network



###############################################################################################################
# Neural network
###############################################################################################################
module NeuralNetwork
using Parameters
using ..AbstractNeuronalTypes, ..EqDiffMethod
using ..attractor_network
export neural_population, neural_motif,wong_wang_network, simulation_parameters


"""
    Struct: neural_population

This tructure defined what consists in a neural population
#TODO Un example de construction de neural population
...
# Arguments
- N = number of neurons in the population
- type neuron::String = specific type of the neurons of the population (like rectified lineae)
- list neurons::Vector{neurons} = list of the neurons of this population (only one type of neuron for now)
- type global::Stirng = type of the neurons composing this neural population (interneuon, pyr_cells, dendrites)
...
"""
struct neural_population
    N::Int64 
    type_neuron::String 
    list_neurons::Vector{neuron}
    type_global::String 
end



"""
    Struct: neural_motif

Defines the structure of what constitutes a neural network

...
# Arguments
- list pop : Vector of all the neural populaitons of the neural motif
- m prob = matrix of the probability of connectivity between all the populations
- m weights = matrix of the synaptic weights
...
"""
struct neural_motif 
    list_pop::Vector{neural_population} 
    m_prob::Array{Float64} 
    m_weights::Array{Float64} 

end



@with_kw struct simulation_parameters
    dt::Float64 = 0.0005
    Tfin::Float64 = 2.0
    Tstimdebut::Float64 = 1.0
    
end

@with_kw struct wong_wang_network 
  noiseamp::Float64  =0.02# noise amp of the OU process
  i0::Float64 = 0.3255 #resting state of the OU process
  jn11::Float64 =0.2609# Synaptic strength unit 1 to 1
  jn21::Float64 = 0.0497# Synaptic strength unit 2 to 1
  jn12::Float64 = 0.0497# Synaptic strength unit 1 to 2
  jn22::Float64 = 0.2609# Synaptic strength unit 2 to 2
  tnmda::Float64 =100.0* 0.001# time constant of NMDA receptors
  tampa::Float64 = 2.0 * 0.001# time constant of AMPA receptors
  threshold::Float64# threhsold of the decision
  list_units::Array{wong_wang_cell} = []

end




end
using .NeuralNetwork



module RateDendrites
using Parameters

using ..AbstractNeuronalTypes, ..EqDiffMethod, ...BasicFunctions
export dendrites_param, soma_Sean2020, dend_Sean2020
# inputs to the rate model will be the time average total conductance of exc and inh synapses


@with_kw struct dendrites_param
    # param of the dendrites from :
    c1::Float64 = 0.12 #* 0.001
    c2::Float64 = 0.13624 #* 0.001
    c3::Float64 = 7.0
    c4::Float64 = 0.0 #* 0.001
    c5::Float64 = 0.00964 #* 0.001
    c6::Float64 = 0.02 #* 0.001 #nA units

end


@with_kw mutable struct soma_Sean2020 <: pyr_cell

    param_c = dendrites_param()
    Idendexc::Float64 = 0.0
    Idendinh::Float64 = 0.0
    f_I_curve = zeros(3)
    r::Array{Float64} = [0.0]
    Iinput::Array{Float64} = [0.0]
    a::Float64 = 0.135 * 1000
    b::Float64 = 54
    c::Float64 = 0.308 #secondes
    den::dendrite #dendrite connected to the soma
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 310.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()

end


@with_kw mutable struct dend_Sean2020 <:dendrite
    param_c = dendrites_param()
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    Ioutput::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 30.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()

end


end

using .RateDendrites



module local_circuit
# let's start with only one dendrites
using Parameters
using ..AbstractNeuronalTypes, ..EqDiffMethod, ...BasicFunctions, ..NeuralNetwork
export pv_cell, sst_cell, vip_cell, gaba_syn, ampa_syn,nmda_syn, local_circuit_interneuron
abstract type excitatory_synapse <: synapses end
abstract type local_circuit_interneuron <: interneuron end
export microcircuit
    
## for now interneurons
#TODO better choice of name (maybe having an abstract pv and so type): Note that everything is static for now (one layer)
@with_kw mutable struct pv_cell <: local_circuit_interneuron
    r::Array{Float64} = [0.0] #rate of the neuron
    c_I::Float64 = 330.0
    r0::Float64 = -95.0
    Iinput::Array{Float64} = [0.0]
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 300.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
end


@with_kw mutable struct sst_cell <: local_circuit_interneuron
  r::Array{Float64} = [0.0] #rate of the neuron
    c_I::Float64 = 132.0
    r0::Float64 = -33.0
    Iinput::Array{Float64} = [0.0]
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 300.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
end


@with_kw mutable struct vip_cell <: local_circuit_interneuron
    r::Array{Float64} = [0.0] #rate of the neuron
    c_I::Float64 = 132.0
    r0::Float64 = -33.0 #Hz
    Iinput::Array{Float64} = [0.0]
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 300.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
end


@with_kw mutable struct gaba_syn <: synapses
    τ::Float64 = 5 * 0.001 #s (except for dendrites)
    neuron_pre::neuron
    neuron_post::neuron
    γ::Float64 =  2
    g::Float64
    s::Array{Float64} = [0.0]
    dt::Float64=0.0005

end

@with_kw mutable struct ampa_syn <: synapses
    τ::Float64 = 2 * 0.001 #s
    neuron_pre::neuron
    neuron_post::neuron
    γ::Float64 = 5
    g::Float64
    s::Array{Float64} = [0.0]
    dt::Float64=0.0005

end

@with_kw mutable struct nmda_syn <: synapses
    neuron_pre::neuron
    neuron_post::neuron
    γ::Float64 = 0.641 * 2
    τ::Float64 = 60 * 0.001 #s
    g::Float64
    s::Array{Float64} = [0.0]
    dt::Float64=0.0005

end



@with_kw mutable struct microcircuit
    list_soma::Array{pyr_cell}= []
    list_sst::Array{sst_cell} = []
    list_vip::Array{vip_cell} = []
    list_pv::Array{pv_cell} = []
    list_dend::Array{dendrite} = []
    nn::Array{wong_wang_network} = []

end


end

using .local_circuit






# pour les synapses il faut garder le rate qui arrive
# Il faut donc une structure separee? 
# 
# 


end