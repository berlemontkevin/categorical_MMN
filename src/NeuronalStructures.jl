module NeuronalStructures

using StaticArrays
export simulation_parameters

export dendrites_param_sigmoid, dend_sigmoid, soma_PC
export dend_sigmoid, dendrite
export neural_integrator, microcircuit
export get_dynamics,save_dynamics
export vip_cell, pv_cell, sst_cell
export gaba_syn, ampa_syn, nmda_syn
export sst_cell_with_adaptation

export neuron

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
using Parameters, StaticArrays
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

using Parameters,StaticArrays
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

using Parameters,StaticArrays


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




module attractor_network
using Parameters, StaticArrays
using ..AbstractNeuronalTypes, ..EqDiffMethod
using ...BasicFunctions

export wong_wang_cell#, neural_integrator, dynamique_integrator

@with_kw struct wong_wang_cell <: neuron
    r::Vector{Float64} = [0.0]
    s::Vector{Float64} = [0.0] # synaptic gatign variable
    Iexc::Vector{Float64} = [0.0]
    Iinh::Vector{Float64} = [0.0]
    Ioutput::Vector{Float64} = [0.0]
    list_syn::Vector{synapses} = []
    Ibg::Vector{Float64} = [0.3255]
    Itot::Vector{Float64} = [0.0]
    Inoise::Vector{Float64} = [0.0]
    Istim::Vector{Float64} = [0.0]
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
using Parameters, StaticArrays
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
    name::String
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
    m_prob::Vector{Float64} 
    m_weights::Vector{Float64} 

end



@with_kw struct simulation_parameters
    dt::Float64 = 0.0005
    Tfin::Float64 = 20.0
    Tinit::Float64 = 2.0
    TISI::Float64 = 1.0    
    Tstimduration::Float64=0.5
    #Stimulus::Dict{String,Vector{Float64, floor(Int,(Tinit+Tfin)/dt)}} = Dict()
    


end

@with_kw struct wong_wang_network 
  noiseamp::Float64  =0.009# noise amp of the OU process
  i0::Float64 = 0.3255 #resting state of the OU process
  jn11::Float64 =0.2609# Synaptic strength unit 1 to 1
  jn21::Float64 = 0.0497# Synaptic strength unit 2 to 1
  jn12::Float64 = 0.0497# Synaptic strength unit 1 to 2
  jn22::Float64 = 0.2609# Synaptic strength unit 2 to 2
  tnmda::Float64 =100.0* 0.001# time constant of NMDA receptors
  tampa::Float64 = 2.0 * 0.001# time constant of AMPA receptors
  threshold::Float64# threhsold of the decision
  list_units::Vector{wong_wang_cell} = []
  EEself::Float64 = 11.0*0.0001
  EEcross::Float64 = 11.0*0.0001/2.0


end




end
using .NeuralNetwork



module RateDendrites
using Parameters,StaticArrays

using ..AbstractNeuronalTypes, ..EqDiffMethod, ...BasicFunctions



end

using .RateDendrites



module local_circuit
# let's start with only one dendrites
using Parameters, StaticArrays
using ..AbstractNeuronalTypes, ..attractor_network, ..EqDiffMethod, ...BasicFunctions, ..NeuralNetwork, ..RateDendrites
export dendrites_param, soma_Sean2020, dend_Sean2020, dendrites_param_sigmoid, dend_sigmoid, adaptation_variables
export dynamique_soma_PC, dynamique_dend_sigmoid
export soma_PC, neural_integrator
export pv_cell, sst_cell, vip_cell, gaba_syn, ampa_syn,nmda_syn, local_circuit_interneuron
export local_circuit_interneuron_with_adaptation, local_circuit_interneuron_without_adaptation, dynamique_variables_interneurons, dynamique_variables_synapses
abstract type excitatory_synapse <: synapses end
abstract type local_circuit_interneuron <: interneuron end
abstract type local_circuit_interneuron_with_adaptation <: local_circuit_interneuron end
abstract type local_circuit_interneuron_without_adaptation <: local_circuit_interneuron end
abstract type dynamique_variables_interneurons end
abstract type dynamique_variables_synapses end

export microcircuit
export save_dynamics, get_dynamics
export sst_cell_with_adaptation    


using DataFrames, DrWatson#,CSV#,CSVFiles




# inputs to the rate model will be the time average total conductance of exc and inh synapses


@with_kw struct facilitation_parameters
    U::Float64 = 0.2
    τ::Float64 = 1.5 # secondes
    dt::Float64 = 0.0005
end

@with_kw struct depression_parameters
    fD::Float64 = 0.2
    τ::Float64 = 0.5 # secondes
    dt::Float64 = 0.0005
end


@with_kw mutable struct state_of_synapse
    s::Float64 = 0.0
    u::Float64 = 1.0/2.5
    d::Float64 = 1.0/2.5
    fr_pre::Float64 = 0.0
    ds::Float64 = 0.0
end

@with_kw mutable struct dynamique_gaba_syn <: dynamique_variables_synapses
    s::Float64 = 0.0
    u::Float64 = 1.0/2.5
    d::Float64 = 1.0/2.5
end

@with_kw struct gaba_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    τ::Float64 = 5 * 0.001 #s (except for dendrites)
    #neuron_pre::neuron
    #neuron_post::neuron
    γ::Float64 =  2
    g::Float64
    #s::MVector{1,Float64} = @MVector [0.0]
    dt::Float64 = 0.0005
    facilitation::Bool = false
    #u::MVector{1,Float64} = @MVector [1.0/2.5]
    mult::Float64 = 2.5 # account for the fact that u is below 1
    f_param::facilitation_parameters = facilitation_parameters()
    d_param::depression_parameters = depression_parameters()
    #d::MVector{1,Float64} = @MVector [1.0/2.5]
    depression::Bool = false
    name::String
    s_save::Vector{Float64} = [0.0]
    u_save::Vector{Float64}=[1.0/2.5]
    d_save::Vector{Float64}=[1.0/2.5]
end


@with_kw struct global_synapse <: synapses
    dynamique_variables::dynamique_gaba_syn = dynamique_gaba_syn()


end



@with_kw struct param_synapse
    f_param = facilitation_parameters()
    d_param = depression_parameters()
    facilitation::Bool = false
    depression::Bool = false


end



@with_kw mutable struct dynamique_ampa_syn <: dynamique_variables_synapses
    s::Float64 = 0.0
    u::Float64 = 1.0/2.5
    d::Float64 = 1.0/2.5
end

@with_kw struct ampa_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    τ::Float64 = 2 * 0.001 #s
   # neuron_pre::neuron
   # neuron_post::neuron
    γ::Float64 = 5
    g::Float64
    #s::MVector{1,Float64} = @MVector [0.0]
    dt::Float64 = 0.0005
    facilitation::Bool = false
    #u::MVector{1,Float64} = @MVector [1.0/2.5]
    mult::Float64 = 2.5 # account for the fact that u is below 1
    f_param::facilitation_parameters = facilitation_parameters()
    d_param::depression_parameters = depression_parameters()
    #d::MVector{1,Float64} = @MVector [1.0/2.5]
    depression::Bool = false
    name::String
    s_save::Vector{Float64} = [0.0]
    u_save::Vector{Float64}=[1.0/2.5]
    d_save::Vector{Float64}=[1.0/2.5]
end


@with_kw mutable struct dynamique_nmda_syn <: dynamique_variables_synapses
    s::Float64 = 0.0
    u::Float64 = 1.0/2.5
    d::Float64 = 1.0/2.5
end



@with_kw struct nmda_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    #neuron_pre::neuron
    #neuron_post::neuron
    γ::Float64 = 0.641 * 2
    τ::Float64 = 60 * 0.001 #s
    g::Float64
    #s::MVector{1,Float64} = @MVector [0.0]
    dt::Float64 = 0.0005
    facilitation::Bool = false
    #u::MVector{1,Float64} = @MVector [1.0/2.5]
    mult::Float64 = 2.5 # account for the fact that u is below 1
    f_param::facilitation_parameters = facilitation_parameters()
    d_param::depression_parameters = depression_parameters()
    #d::MVector{1,Float64} = @MVector [0.3]
    depression::Bool = false
    name::String
    s_save::Vector{Float64} = [0.0]
    u_save::Vector{Float64}=[1.0/2.5]
    d_save::Vector{Float64}=[1.0/2.5]
end


@with_kw struct dendrites_param
    # param of the dendrites from :
    c1::Float64 = 0.12 #* 0.001
    c2::Float64 = 0.13624 #* 0.001
    c3::Float64 = 7.0
    c4::Float64 = 0.0 #* 0.001
    c5::Float64 = 0.00964 #* 0.001
    c6::Float64 = 0.02 #* 0.001 #nA units

end

@with_kw struct adaptation_variables
    #TODO ut it in another module
    sA::Vector{Float64} = [0.0]
    sA_save::Vector{Float64} = [0.0]

    τA::Float64 = 0.1 #seconds as for facilitation and 0.05
    gA::Float64 = -0.01 # value to test 0.5, was -0.1 before
end


struct dendrites_param_sigmoid
     # param of the dendrites from : (units to precise)
     c1::Float64 
     c2::Float64 
     c3::Float64 
     c4::Float64 
     c5::Float64 
     c6::Float64 
end

#TODO define a constant param per neuron
@with_kw mutable struct dynamique_dend_sigmoid
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ioutput::Float64 = 0.0
    Ibg::Float64 = 100.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0
end

@with_kw struct dend_sigmoid <:dendrite
    param_c::dendrites_param_sigmoid = [0.0,0.0,0.0,0.0,0.0,0.0]
    #Iexc::MVector{1,Float64} = @MVector [0.0]
    #Iinh::MVector{1,Float64} = @MVector [0.0]
    #Ioutput::MVector{1,Float64} = @MVector [0.0]
    list_syn::Vector{synapses} = []

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()

    Ibg::Float64 = 100.0 * 0.001
    #Itot::MVector{1,Float64} = @MVector [0.0]
    Inoise::Vector{Float64} = [0.0]
    #Istim::MVector{1,Float64} = @MVector [0.0]
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
    name::String = "temp"
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ioutput_save::Vector{Float64} = [0.0]
    Itot_save::Vector{Float64} = [0.0]
    Inoise_save::Vector{Float64} = [0.0]
    dynamique_variables::dynamique_dend_sigmoid = dynamique_dend_sigmoid()
  end


@with_kw mutable struct dynamique_soma_PC
    Idendexc::Float64 = 0.0
    Idendinh::Float64 = 0.0
    r::Float64 = 0.0
    Iinput::Float64 = 0.0
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0
    Ibg::Float64 = 0.15#0.35#310.0 * 0.001

end


@with_kw struct soma_PC <: pyr_cell
    dynamique_variables::dynamique_soma_PC = dynamique_soma_PC()
  #  Idendexc::MVector{1,Float64} = @MVector [0.0]
   # Idendinh::MVector{1,Float64} = @MVector [0.0]
    f_I_curve = zeros(3)
   # r::MVector{1,Float64} = @MVector [0.0]
   # Iinput::MVector{1,Float64} = @MVector [0.0]
    a::Float64 = 135.0
    b::Float64 = 54.0
    c::Float64 = 0.308 #secondes
    den::dend_sigmoid #dendrite connected to the soma
    #Iexc::MVector{1,Float64} = @MVector [0.0]
    #Iinh::MVector{1,Float64} = @MVector [0.0]
    list_syn::Vector{synapses} = Vector{synapses}()

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()

    Ibg::Float64 = 0.15#0.35#310.0 * 0.001
    #Itot::MVector{1,Float64} = @MVector [0.0]
    Inoise::Vector{Float64} = [0.0]
    #Istim::MVector{1,Float64} = @MVector [0.0]
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
    adaptation::adaptation_variables = adaptation_variables()
    adaptation_boolean::Bool = false # boolean of adaptation or not
    name::String
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ioutput_save::Vector{Float64} = [0.0]
    Itot_save::Vector{Float64} = [0.0]
    Inoise_save::Vector{Float64} = [0.0]
    r_save::Vector{Float64} = [0.0]

end

@with_kw mutable struct soma_Sean2020 <: pyr_cell

    param_c = dendrites_param()
    Idendexc::Float64 = 0.0
    Idendinh::Float64 = 0.0
    f_I_curve = zeros(3)
    r::Float64 = 0.0
    Iinput::Float64 = 0.0
    a::Float64 = 0.135 * 1000
    b::Float64 = 54
    c::Float64 = 0.308 #secondes
    den::dendrite #dendrite connected to the soma
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    list_syn::Vector{synapses} = Vector{synapses}()
    Ibg::Float64 = 330.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Vector{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ioutput_save::Vector{Float64} = [0.0]
    Itot_save::Vector{Float64} = [0.0]
    Inoise_save::Vector{Float64} = [0.0]
    r_save::Vector{Float64} = [0.0]
end


@with_kw mutable struct dend_Sean2020 <:dendrite
    param_c = dendrites_param()
    Iexc::Vector{Float64} = [0.0]
    Iinh::Vector{Float64} = [0.0]
    Ioutput::Vector{Float64} = [0.0]
    list_syn::Vector{synapses} = Vector{synapses}()
    Ibg::Float64 = 30.0 * 0.001
    Itot::Vector{Float64} = [0.0]
    Inoise::Vector{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()

end


## for now interneurons

@with_kw mutable struct dynamique_pv_cell <: dynamique_variables_interneurons
    r::Float64 = 0.0		
    Iinput::Float64 = 0.0
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 300.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0

end

#TODO better choice of name (maybe having an abstract pv and so type): Note that everything is static for now (one layer)
# @with_kw struct pv_cell <: local_circuit_interneuron_without_adaptation
#     dynamique_variables::dynamique_pv_cell = dynamique_pv_cell()
#     #r::MVector{1,Float64} = @MVector [0.0]		
#     c_I::Float64 = 330.0
#     r0::Float64 = -95.0
#     #Iinput::MVector{1,Float64} = @MVector [0.0]
#     #Iexc::MVector{1,Float64} = @MVector [0.0]
#     #Iinh::MVector{1,Float64} = @MVector [0.0]
#     list_syn::Vector{synapses} = Vector{synapses}()
#     Ibg::Float64 = 300.0 * 0.001
#     #Itot::MVector{1,Float64} = @MVector [0.0]
#     Inoise::Vector{Float64} = [0.0]
#     #Istim::MVector{1,Float64} = @MVector [0.0]
#     dt::Float64 = 0.0005
#     τ::Float64 = 0.002
#     OU_process::OU_process = OU_process()
#     adaptation::adaptation_variables = adaptation_variables()
#     adaptation_boolean = false # boolean of adaptation or not
#     name::String
#     #Iexc_save::Vector{Float64} = [0.0]
#     #Iinh_save::Vector{Float64} = [0.0]
#     #Ioutput_save::Vector{Float64} = [0.0]
#     #Itot_save::Vector{Float64} = [0.0]
#     #Inoise_save::Vector{Float64} = [0.0]
#     r_save::Vector{Float64} = [0.0]
# end


@with_kw struct pv_cell <: local_circuit_interneuron
    dynamique_variables::dynamique_pv_cell = dynamique_pv_cell()
    #r::MVector{1,Float64} = @MVector [0.0]		
    c_I::Float64 = 330.0
    r0::Float64 = -95.0
    #Iinput::MVector{1,Float64} = @MVector [0.0]
    #Iexc::MVector{1,Float64} = @MVector [0.0]
    #Iinh::MVector{1,Float64} = @MVector [0.0]
    list_syn::Vector{synapses} = Vector{synapses}()

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()


    Ibg::Float64 = 300.0 * 0.001
    #Itot::MVector{1,Float64} = @MVector [0.0]
    Inoise::Vector{Float64} = [0.0]
    #Istim::MVector{1,Float64} = @MVector [0.0]
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
    adaptation::adaptation_variables = adaptation_variables()
    adaptation_boolean::Bool = true # boolean of adaptation or not
    name::String
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ioutput_save::Vector{Float64} = [0.0]
    Itot_save::Vector{Float64} = [0.0]
    Inoise_save::Vector{Float64} = [0.0]
    r_save::Vector{Float64} = [0.0]


end


@with_kw mutable struct dynamique_sst_cell <: dynamique_variables_interneurons
    r::Float64 = 0.0		
    Iinput::Float64 = 0.0
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 300.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0

end

# @with_kw struct sst_cell <: local_circuit_interneuron_without_adaptation
#     dynamique_variables::dynamique_sst_cell = dynamique_sst_cell()
#     #r::MVector{1,Float64} = @MVector [0.0]		
#     c_I::Float64 = 132.0
#     r0::Float64 = -33.0
#     #Iinput::MVector{1,Float64} = @MVector [0.0]
#     #Iexc::MVector{1,Float64} = @MVector [0.0]
#     #Iinh::MVector{1,Float64} = @MVector [0.0]
#     list_syn::Vector{synapses} = Vector{synapses}()
#     Ibg::Float64 = 300.0 * 0.001
#     #Itot::MVector{1,Float64} = @MVector [0.0]
#     Inoise::Vector{Float64} = [0.0]
#     #Istim::MVector{1,Float64} = @MVector [0.0]
#     dt::Float64 = 0.0005
#     τ::Float64 = 0.002
#     adaptation::adaptation_variables = adaptation_variables()
#     adaptation_boolean = false # boolean of adaptation or not
#     OU_process::OU_process = OU_process()
#     name::String
#     Iexc_save::Vector{Float64} = [0.0]
#     Iinh_save::Vector{Float64} = [0.0]
#     Ioutput_save::Vector{Float64} = [0.0]
#     Itot_save::Vector{Float64} = [0.0]
#     Inoise_save::Vector{Float64} = [0.0]
#     r_save::Vector{Float64} = [0.0]
# end

@with_kw struct sst_cell <: local_circuit_interneuron
    dynamique_variables::dynamique_sst_cell = dynamique_sst_cell()
    #r::MVector{1,Float64} = @MVector [0.0]		
    c_I::Float64 = 132.0
    r0::Float64 = -33.0
    #Iinput::MVector{1,Float64} = @MVector [0.0]
    #Iexc::MVector{1,Float64} = @MVector [0.0]
    #Iinh::MVector{1,Float64} = @MVector [0.0]
    list_syn::Vector{synapses} = Vector{synapses}()
    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()


    Ibg::Float64 = 300.0 * 0.001
    #Itot::MVector{1,Float64} = @MVector [0.0]
    Inoise::Vector{Float64} = [0.0]
    #Istim::MVector{1,Float64} = @MVector [0.0]
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    adaptation::adaptation_variables = adaptation_variables()
    adaptation_boolean::Bool = true # boolean of adaptation or not
    OU_process::OU_process = OU_process()
    name::String
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ioutput_save::Vector{Float64} = [0.0]
    Itot_save::Vector{Float64} = [0.0]
    Inoise_save::Vector{Float64} = [0.0]
    r_save::Vector{Float64} = [0.0]
    
end

@with_kw mutable struct dynamique_vip_cell <: dynamique_variables_interneurons
    r::Float64 = 0.0		
    Iinput::Float64 = 0.0
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 300.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0

end

# @with_kw struct vip_cell <: local_circuit_interneuron_without_adaptation
#     dynamique_variables::dynamique_vip_cell = dynamique_vip_cell()
#     #r::MVector{1,Float64} = @MVector [0.0]		
#     c_I::Float64 = 132.0
#     r0::Float64 = -33.0 #Hz
#     #Iinput::MVector{1,Float64} = @MVector [0.0]
#     #Iexc::MVector{1,Float64} = @MVector [0.0]
#     #Iinh::MVector{1,Float64} = @MVector [0.0]
#     list_syn::Vector{synapses} = Vector{synapses}()
#     Ibg::Float64 = 300.0 * 0.001
#     #Itot::MVector{1,Float64} = @MVector [0.0]
#     Inoise::Vector{Float64} = [0.0]
#     #Istim::MVector{1,Float64} = @MVector [0.0]
#     dt::Float64 = 0.0005
#     τ::Float64 = 0.002
#     adaptation::adaptation_variables = adaptation_variables()
#     adaptation_boolean = false # boolean of adaptation or not
#     OU_process::OU_process = OU_process()
#     name::String
#     Iexc_save::Vector{Float64} = [0.0]
#     Iinh_save::Vector{Float64} = [0.0]
#     Ioutput_save::Vector{Float64} = [0.0]
#     Itot_save::Vector{Float64} = [0.0]
#     Inoise_save::Vector{Float64} = [0.0]
#     r_save::Vector{Float64} = [0.0]
# end


@with_kw struct vip_cell <: local_circuit_interneuron
    dynamique_variables::dynamique_vip_cell = dynamique_vip_cell()
    #r::MVector{1,Float64} = @MVector [0.0]		
    c_I::Float64 = 132.0
    r0::Float64 = -33.0 #Hz
    #Iinput::MVector{1,Float64} = @MVector [0.0]
    #Iexc::MVector{1,Float64} = @MVector [0.0]
    #Iinh::MVector{1,Float64} = @MVector [0.0]
    list_syn::Vector{synapses} = Vector{synapses}()

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()

    Ibg::Float64 = 300.0 * 0.001
    #Itot::MVector{1,Float64} = @MVector [0.0]
    Inoise::Vector{Float64} = [0.0]
    #Istim::MVector{1,Float64} = @MVector [0.0]
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    adaptation::adaptation_variables = adaptation_variables()
    adaptation_boolean::Bool = true # boolean of adaptation or not
    OU_process::OU_process = OU_process()
    name::String
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ioutput_save::Vector{Float64} = [0.0]
    Itot_save::Vector{Float64} = [0.0]
    Inoise_save::Vector{Float64} = [0.0]
    r_save::Vector{Float64} = [0.0]
end

@with_kw mutable struct dynamique_integrator
    r::Float64 = 0.0		
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 310.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0
end


 # will be in simulation parameters later for the time dt
@with_kw struct neural_integrator <: neuron
    dynamique_variables::dynamique_integrator = dynamique_integrator()
    τ::Float64 = 0.8 # test d'une cste de temps
    dt::Float64 = 0.0005
    α::Float64 = 0.9  #0.8 before
    a::Float64 = 135.0
    b::Float64 = 54.0
    c::Float64 = 0.308 #secondes
    r::MVector{1,Float64} = @MVector [0.0]		
    Iexc::MVector{1,Float64} = @MVector [0.0]
    Iinh::MVector{1,Float64} = @MVector [0.0]
   # list_syn::Vector{synapses} = []

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()


    Ibg::Float64 = 310.0 * 0.001
    Itot::MVector{1,Float64} = @MVector [0.0]
    Inoise::Vector{Float64} = [0.0]
    Istim::MVector{1,Float64} = @MVector [0.0]
    OU_process::OU_process = OU_process()
    name::String
    r_save::Vector{Float64}= [0.0]		
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ibg_save::Float64 = 310.0 * 0.001
    Itot_save::Vector{Float64} = [0.0]
end


@with_kw mutable struct microcircuit
    list_soma::Vector{pyr_cell}= Vector{pyr_cell}()
    list_sst::Vector{local_circuit_interneuron} = []
    list_vip::Vector{vip_cell} = []
    list_pv::Vector{pv_cell} = []
    list_dend::Vector{dendrite} = []
    nn::Vector{wong_wang_network} = []
    list_integrator::Vector{neural_integrator} = []
    name::String = "microciruit"
end




function save_dynamics(c::microcircuit,notebook_file::String, save_parameters)
    # funciton that save the dynamics of the network into a csv

    df = DataFrame()

    
    for pop in [c.list_pv, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop

                df[!,string("r-",n.name)] = round.(n.r,digits=6)
                df[!,string("Itot-",n.name)] = round.(n.Itot,digits=6)

            for s in n.list_syn
                df[!,string("s-",s.name)] = round.(s.s,digits=6)

            end
        end

    end
    for pop in [c.list_dend]
        for n in pop
            df[!,string("Ioutput-",n.name)] = round.(n.Ioutput, digits=6)

           
        end
    end


    wsave(datadir("simulations",notebook_file, savename(save_parameters, "csv")),df)



end




function get_dynamics(notebook_file::String, save_parameters)
    # funciton that retrieves the dynamics from the csv

 #   df = wload(datadir("simulations",notebook_file, savename(save_parameters, "csv")))
    df = CSV.read(datadir("simulations",notebook_file, savename(save_parameters, "csv")), DataFrame)
    dtemp = wload(datadir("simulations",notebook_file, savename(save_parameters, "jld2")))
    c = dtemp["circuit"]

    for pop in [c.list_pv, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop

            temp_name = string("r-",n.name)
            n.r = df[:,temp_name]
            temp_name = string("Itot-",n.name)
            n.Itot = df[:,temp_name]

            for s in n.list_syn
                temp_name = string("s-",s.name)
                s.s = df[:,temp_name]
            end
        end

    end

    for pop in [c.list_dend]
        for n in pop

            temp_name = string("Ioutput-",n.name)
            n.Ioutput = df[:,temp_name]
           
        end
    end
    return c

end


end

using .local_circuit






# pour les synapses il faut garder le rate qui arrive
# Il faut donc une structure separee? 
# 
# 


end