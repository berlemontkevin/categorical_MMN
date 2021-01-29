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

###############################################################################################################
# Neural network
###############################################################################################################
module NeuralNetwork
using Parameters
using ..AbstractNeuronalTypes, ..EqDiffMethod
export neural_population, neural_motif


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

end
using .NeuralNetwork







end