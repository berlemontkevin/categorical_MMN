module NeuronalStructures

using StaticArrays
using ..BasicFunctions

export simulation_parameters

export dendrites_param_sigmoid, dend_sigmoid, soma_PC
export dend_sigmoid, dendrite
export neural_integrator, microcircuit
export get_dynamics,save_dynamics
export vip_cell, pv_cell, sst_cell
export gaba_syn, ampa_syn, nmda_syn
export sst_cell_with_adaptation
export bump_structure, layer_bump
export parameters_bump_attractor, parameters_inter_microcircuit, parameters_interneurons, parameters_microcircuit, parameters_syn_strength_microcircuit
export neuron

##############################################################################################################
# Abstract types
###############################################################################################################
module AbstractNeuronalTypes
export neuron, dendrite, pyr_cell, synapses, interneuron
export excitatory_synapse, local_circuit_interneuron
export dynamique_variables_interneurons, dynamique_variables_synapses
export eq_diff_method

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


"""
    Type: excitatory_synapse
    
Abstract type for the excitatory synapses
"""
abstract type excitatory_synapse <: synapses end

"""
    Type: local_circuit_interneuron
    
Abstract type for the local circuits with many interneurons types
"""
abstract type local_circuit_interneuron <: interneuron end

"""
    Type: dynamique_variables_interneurons
    
Abstract type for the dynamical variables (memory allocation) of the interneurons
"""
abstract type dynamique_variables_interneurons end

"""
    Type: dynamique_variables_synapses    
Abstract type for the dynamical variables (memory allocation) of the synapses
"""
abstract type dynamique_variables_synapses end

abstract type eq_diff_method end



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

...
# Arguments
- dt::Float64 = time step of the Euler resolution
- record::Bool = does the dynamics record the activity of the networks
...
    """
@with_kw struct euler_method <: eq_diff_method
    dt::Float64 = 0.0005
    record::Bool = true
end

end
using .EqDiffMethod

# for now everything is always solved with Euler Method






###############################################################################################################
# Synapses
###############################################################################################################

module Synapses

using ..AbstractNeuronalTypes, ...BasicFunctions

using Parameters,StaticArrays



export facilitation_parameters, depression_parameters
export state_of_synapse
export gaba_syn, ampa_syn, nmda_syn


"""
    Struct: facilitation_parameters
    
Parameters for facilitaiton as synaptic plasticity

    ...
    # Arguments
    - U::Float64 = 
    - τ::Float64 = time constant of the facilitation
    ...

"""
@with_kw struct facilitation_parameters
    U::Float64 = 0.2
    τ::Float64 = 1.5 # secondes
end



"""
    Struct: depression_parameters
    
Parameters for depression as synaptic plasticity

    ...
    # Arguments
    - fD::Float64 = 
    - τ::Float64 = time constant of the depression
    ...

"""
@with_kw struct depression_parameters
    fD::Float64 = 0.2
    τ::Float64 = 0.5 # secondes
end

"""
    Struct: state_of_synapse
    
Mutable structure that contains all the dynamically allocated variables for the synaptic equations
    ...
    # Arguments
    - s::Float64 = synaptic gating variables 
    - u::Float64 = facilitation variable starting point
    - d::Float64 = depression variable starting point
    - fr_pre::Float64 = starting point of the firing rate of the presynaptic neuron
    - ds::Float64 = starting point for the derivative of the gating variable
    ...

"""
@with_kw mutable struct state_of_synapse <: dynamique_variables_synapses
    s::Float64 = 0.0
    u::Float64 = 1.0 / 2.5
    d::Float64 = 1.0 / 2.5
    fr_pre::Float64 = 0.0
    ds::Float64 = 0.0
end


"""
    Struct: gaba_syn

    Structure for the synapses GABA 
    ...
        # Arguments
        - dynamique_variables::state_of_synapse = 
    ...
"""
@with_kw struct gaba_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    τ::Float64 = 5 * 0.001 # s (except for dendrites)
    γ::Float64 =  2
    g::Float64
    dt::Float64 = 0.0005
    facilitation::Bool = false
    mult::Float64 = 2.5 # account for the fact that u is below 1
    f_param::facilitation_parameters = facilitation_parameters()
    d_param::depression_parameters = depression_parameters()
    depression::Bool = false
    name::String
    s_save::Vector{Float64} = [0.0]
    u_save::Vector{Float64} = [1.0 / 2.5]
    d_save::Vector{Float64} = [1.0 / 2.5]
end

@with_kw struct ampa_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    τ::Float64 = 2 * 0.001 # s
    γ::Float64 = 5
    g::Float64
    dt::Float64 = 0.0005
    facilitation::Bool = false
    mult::Float64 = 2.5 # account for the fact that u is below 1
    f_param::facilitation_parameters = facilitation_parameters()
    d_param::depression_parameters = depression_parameters()
    depression::Bool = false
    name::String
    s_save::Vector{Float64} = [0.0]
    u_save::Vector{Float64} = [1.0 / 2.5]
    d_save::Vector{Float64} = [1.0 / 2.5]
end


@with_kw struct nmda_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    γ::Float64 = 0.641 * 2
    τ::Float64 = 60 * 0.001 # s
    g::Float64
    dt::Float64 = 0.0005
    facilitation::Bool = false
    mult::Float64 = 2.5 # account for the fact that u is below 1
    f_param::facilitation_parameters = facilitation_parameters()
    d_param::depression_parameters = depression_parameters()
    depression::Bool = false
    name::String
    s_save::Vector{Float64} = [0.0]
    u_save::Vector{Float64} = [1.0 / 2.5]
    d_save::Vector{Float64} = [1.0 / 2.5]
end


end
using .Synapses

###############################################################################################################
# Neuronal Models
###############################################################################################################
module NeuronalModels
export rectified_linear_neurons

using ..AbstractNeuronalTypes, ..Synapses, ...BasicFunctions

using Parameters,StaticArrays


export dend_sigmoid

"""
    Struct: rectified_linear_neurons

It consists in the model of a neuron as a rectified_linear neuron

...
# Arguments
- Isyn::Float64 = synaptic input toward the neuron
- τ::Float64 = time constant of the neuron
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
    fr::Vector{Float64} = [0.0] 
    Ibg::Vector{Float64} = [0.0] 
    objective_rate::Float64 = -1.0
end


@with_kw struct dendrites_param_sigmoid
    # param of the dendrites from : (units to precise)
    c1::Float64 = 0.12
    c2::Float64 = -7.0
    c3::Float64 = -0.482
    c4::Float64 = 0.00964
    c5::Float64 = 0.19624
    c6::Float64 = 0.0
end

@with_kw mutable struct dynamique_dend_sigmoid
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ioutput::Float64 = 0.0
    Ibg::Float64 = 100.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0
end

@with_kw struct dend_sigmoid <: dendrite
    param_c::dendrites_param_sigmoid = dendrites_param_sigmoid()
    list_syn::Vector{synapses} = []

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()

    Ibg::Float64 = 100.0 * 0.001
    Inoise::Vector{Float64} = [0.0]
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






@with_kw mutable struct adaptation_variables
    sA::Float64 = 0.0
    sA_save::Vector{Float64} = [0.0]
    τA::Float64 = 0.1 # seconds as for facilitation and 0.05
    gA::Float64 = -0.01 # value to test 0.5, was -0.1 before
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
    Ibg::Float64 = 0.15
end



@with_kw struct soma_PC <: pyr_cell
    dynamique_variables::dynamique_soma_PC = dynamique_soma_PC()
    f_I_curve = zeros(3)
    a::Float64 = 135.0
    b::Float64 = 54.0
    c::Float64 = 0.308 # secondes
    den::dend_sigmoid = dend_sigmoid() # dendrite connected to the soma
    list_syn::Vector{synapses} = Vector{synapses}()

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()

    Ibg::Float64 = 0.15# 0.35#310.0 * 0.001
    Inoise::Vector{Float64} = [0.0]
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
    
    preferred_stim::Float64 = 0.0 
    
end
export soma_PC


@with_kw mutable struct dynamique_rectified_linear <: dynamique_variables_interneurons
    r::Float64 = 0.0		
    Iinput::Float64 = 0.0
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 300.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0

end


@with_kw struct pv_cell <: local_circuit_interneuron
    dynamique_variables::dynamique_rectified_linear = dynamique_rectified_linear()
    c_I::Float64 = 330.0
    r0::Float64 = -95.0
    list_syn::Vector{synapses} = Vector{synapses}()

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()

    Ibg::Float64 = 300.0 * 0.001
    Inoise::Vector{Float64} = [0.0]
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
export pv_cell


@with_kw struct sst_cell <: local_circuit_interneuron
    dynamique_variables::dynamique_rectified_linear = dynamique_rectified_linear()
    c_I::Float64 = 132.0
    r0::Float64 = -33.0
    list_syn::Vector{synapses} = Vector{synapses}()
    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()


    Ibg::Float64 = 300.0 * 0.001
    Inoise::Vector{Float64} = [0.0]
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
export sst_cell



@with_kw struct vip_cell <: local_circuit_interneuron
    dynamique_variables::dynamique_rectified_linear = dynamique_rectified_linear()
    c_I::Float64 = 132.0
    r0::Float64 = -33.0 # Hz
    list_syn::Vector{synapses} = Vector{synapses}()

    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()

    Ibg::Float64 = 300.0 * 0.001
    Inoise::Vector{Float64} = [0.0]
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
export vip_cell



 # will be in simulation parameters later for the time dt
 @with_kw struct neural_integrator <: neuron
    dynamique_variables::dynamique_rectified_linear = dynamique_rectified_linear()
    τ::Float64 = 0.8 
    dt::Float64 = 0.0005
    α::Float64 = 0.9  # 0.8 before
    a::Float64 = 135.0
    b::Float64 = 54.0
    c::Float64 = 0.308 # secondes
    
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
    r_save::Vector{Float64} = [0.0]		
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ibg_save::Float64 = 310.0 * 0.001
    Itot_save::Vector{Float64} = [0.0]
end
export neural_integrator

end

using .NeuronalModels



###############################################################################################################
# Simulation types and parameters
###############################################################################################################
module Simulations

using Parameters
export  simulation_parameters




@with_kw struct simulation_parameters
    
    Tfin::Float64 = 20.0
    Tinit::Float64 = 2.0
    TISI::Float64 = 1.0    
    Tstimduration::Float64 = 0.5
    # Stimulus::Dict{String,Vector{Float64, floor(Int,(Tinit+Tfin)/dt)}} = Dict()
    
end


export parameters_bump_attractor, parameters_inter_microcircuit, parameters_interneurons, parameters_microcircuit, parameters_syn_strength_microcircuit

"""
    Type: parameters_microcircuit

List of all the parameters used for the simulations for the microcircuit
...
# Functions used in
- construct_two_local_microcircuit_integrator
....
"""
@with_kw struct parameters_microcircuit
    dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0)
    sst_adaptation::Bool = true
    soma_adaptation::Bool = true
    pv_to_soma_depression::Bool = true
    soma_to_vip_facilitation::Bool = true
    sst_to_vip_facilitation::Bool = true
    soma_to_sst_facilitation::Bool = true
    vip_to_sst_facilitation::Bool = true
    soma_to_pv_depression::Bool = true
    int_to_vip_depression::Bool = true
    int_to_pv_depression::Bool = true
    int_to_sst_facilitation::Bool = true
    int_to_dend_depression::Bool = true
    integrator_tc::Float64 = 0.8
    time_tot::Int = 1000
    noise::Bool = true
    int_to_sst_connection::Bool = true
    top_down_to_interneurons::Vector{Float64} = [0.47, 0.31, 0.22]
    preferred::Float64 = 0.0

end


@with_kw struct parameters_syn_strength_microcircuit
    gaba_sst_to_dend::Float64 = -0.09
    gaba_pv_to_soma::Float64 = -0.01
    nmda_soma_to_soma::Float64 = 0.18
    nmda_soma_to_vip::Float64 = -0.058
    gaba_sst_to_vip::Float64 = -0.1
    nmda_soma_to_sst::Float64 = 0.0435
    gaba_vip_to_sst::Float64 = -0.05
    nmda_soma_to_pv::Float64 = 0.17
    gaba_sst_to_pv::Float64 = -0.17
    gaba_pv_to_pv::Float64 = -0.18
    nmda_soma_to_int::Float64 = 0.15
    nmda_int_to_dend::Float64 = 0.4

end


@with_kw struct parameters_interneurons

    Ibg_vip::Float64 = 0.25
    Ibg_sst::Float64 = 0.25
    Ibg_pv::Float64 = 0.29
    top_down_to_interneurons::Vector{Float64} = [0.47, 0.31, 0.22]


end


@with_kw struct parameters_inter_microcircuit
    cross_int_to_vip_depression::Bool = true
    cross_int_to_pv_depression::Bool = true
    cross_int_to_sst_facilitation::Bool = true
    cross_int_to_dend_depression::Bool = true
    cross_soma_to_sst_facilitation::Bool = true
end

@with_kw struct parameters_bump_attractor
    num_circuits::Int = 128
end



@with_kw struct bump_structure
    num_circuits::Int64 = 128
    σ::Float64 = 43.2

end
export bump_structure


end
using .Simulations






###############################################################################################################
# Neural network
###############################################################################################################
module NeuralNetwork
using Parameters, StaticArrays
using ..AbstractNeuronalTypes, ..Simulations, ..NeuronalModels



@with_kw struct microcircuit{T <: pyr_cell, D <: dendrite}
    soma::T 
    sst::sst_cell
    vip::vip_cell
    pv::pv_cell
    dend::D
    integrator::neural_integrator
    name::String = "microciruit"
end
export microcircuit


@with_kw struct layer_bump{T <: pyr_cell, D <: dendrite}
    bump_param::bump_structure
    list_microcircuit::SVector{128,microcircuit{T,D}}
    eq_diff_method::eq_diff_method = euler_method()
end
export layer_bump




end
using .NeuralNetwork





module local_circuit
# let's start with only one dendrites
using Parameters, StaticArrays
using ..AbstractNeuronalTypes, ...BasicFunctions, ..NeuralNetwork
export dendrites_param, soma_Sean2020, dend_Sean2020, dendrites_param_sigmoid, dend_sigmoid, adaptation_variables
export dynamique_soma_PC, dynamique_dend_sigmoid
export soma_PC, neural_integrator
export pv_cell, sst_cell, vip_cell, gaba_syn, ampa_syn,nmda_syn, local_circuit_interneuron
export local_circuit_interneuron_with_adaptation, local_circuit_interneuron_without_adaptation, dynamique_variables_interneurons, dynamique_variables_synapses




export save_dynamics, get_dynamics
export sst_cell_with_adaptation    
export bump_structure, layer_bump



using DataFrames, DrWatson# ,CSV#,CSVFiles




# inputs to the rate model will be the time average total conductance of exc and inh synapses









@with_kw struct dendrites_param
    # param of the dendrites from :
    c1::Float64 = 0.12 # * 0.001
    c2::Float64 = 0.13624 # * 0.001
    c3::Float64 = 7.0
    c4::Float64 = 0.0 # * 0.001
    c5::Float64 = 0.00964 # * 0.001
    c6::Float64 = 0.02 # * 0.001 #nA units

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
    c::Float64 = 0.308 # secondes
    den::dendrite # dendrite connected to the soma
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


@with_kw mutable struct dend_Sean2020 <: dendrite
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



@with_kw mutable struct dynamique_integrator
    r::Float64 = 0.0		
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 310.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0
end







function save_dynamics(c::microcircuit, notebook_file::String, save_parameters)
    # funciton that save the dynamics of the network into a csv

    df = DataFrame()

    
    for pop in [c.list_pv, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop

                df[!,string("r-", n.name)] = round.(n.r, digits=6)
                df[!,string("Itot-", n.name)] = round.(n.Itot, digits=6)

            for s in n.list_syn
                df[!,string("s-", s.name)] = round.(s.s, digits=6)

            end
        end

    end
    for pop in [c.list_dend]
        for n in pop
            df[!,string("Ioutput-", n.name)] = round.(n.Ioutput, digits=6)

           
        end
    end


    wsave(datadir("simulations", notebook_file, savename(save_parameters, "csv")), df)



end




function get_dynamics(notebook_file::String, save_parameters)
    # funciton that retrieves the dynamics from the csv

 #   df = wload(datadir("simulations",notebook_file, savename(save_parameters, "csv")))
    df = CSV.read(datadir("simulations", notebook_file, savename(save_parameters, "csv")), DataFrame)
    dtemp = wload(datadir("simulations", notebook_file, savename(save_parameters, "jld2")))
    c = dtemp["circuit"]

    for pop in [c.list_pv, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop

            temp_name = string("r-", n.name)
            n.r = df[:,temp_name]
            temp_name = string("Itot-", n.name)
            n.Itot = df[:,temp_name]

            for s in n.list_syn
                temp_name = string("s-", s.name)
                s.s = df[:,temp_name]
            end
        end

    end

    for pop in [c.list_dend]
        for n in pop

            temp_name = string("Ioutput-", n.name)
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