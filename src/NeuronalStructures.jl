module NeuronalStructures

using StaticArrays
using ..BasicFunctions

export simulation_parameters

export dendrites_param_sigmoid, dend_sigmoid, soma_PC
export dend_sigmoid, dendrite
export neural_integrator, microcircuit
export get_dynamics,save_dynamics
export vip_cell, pv_cell, sst_cell, ngfc_cell
export gaba_syn, ampa_syn, nmda_syn
export sst_cell_with_adaptation
export layer_bump
export parameters_bump_attractor, parameters_inter_microcircuit, parameters_interneurons, parameters_microcircuit, parameters_syn_strength_microcircuit
export neuron, euler_method
export adaptation_variables
export get_firing_rate

export ring_microcircuit, ring_integrator, neural_network


##############################################################################################################
# Abstract types
###############################################################################################################
module AbstractNeuronalTypes
export neuron, dendrite, pyr_cell, synapses, interneuron
export excitatory_synapse, local_circuit_interneuron
export dynamique_variables_interneurons, dynamique_variables_synapses
export eq_diff_method

"""
Abstract type for all the type of neurons

# Neurons currently implemented

- dendrites
- pyramidal cells
- interneurons
"""
abstract type neuron end



"""
  
Abstract type for the types of synapses
"""
abstract type synapses end


"""
Abstract type for the types of dendrites
"""
abstract type dendrite <: neuron end

"""
Abstract type for the types of pyramidal cells
"""
abstract type pyr_cell <: neuron end

"""
Abstract type for the types of interneuron
"""
abstract type interneuron <: neuron end


"""
Abstract type for the excitatory synapses
"""
abstract type excitatory_synapse <: synapses end

"""
Abstract type for the local circuits with many interneurons types
"""
abstract type local_circuit_interneuron <: interneuron end

"""
Abstract type for the dynamical variables (memory allocation) of the interneurons
"""
abstract type dynamique_variables_interneurons end

"""
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
Defines the parameters for euler method. dt is the time step and record indicates if the simulation is going to record all the time steps or not.

# Fields

- dt::Float64 = time step of the Euler resolution
- record::Bool = does the dynamics record the activity of the networks
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
Parameters for facilitaiton as synaptic plasticity

# Fields

- `U::Float64` = 
- `τ::Float64` = time constant of the facilitation

"""
@with_kw struct facilitation_parameters
    U::Float64 = 0.2
    τ::Float64 = 1.5 # secondes
end



"""
Parameters for depression as synaptic plasticity

# Fields

- `fD::Float64` = 
- `τ::Float64` = time constant of the depression

"""
@with_kw struct depression_parameters
    fD::Float64 = 0.2
    τ::Float64 = 0.2 # secondes
end

"""
   
Mutable structure that contains all the dynamically allocated variables for the synaptic equations

# Fields
    
- `s::Float64` = synaptic gating variables 
- `u::Float64` = facilitation variable starting point
- `d::Float64` = depression variable starting point
- `fr_pre::Float64` = starting point of the firing rate of the presynaptic neuron
- `ds::Float64` = starting point for the derivative of the gating variable
    
"""
@with_kw mutable struct state_of_synapse <: dynamique_variables_synapses
    s::Float64 = 0.0
    u::Float64 = 1.0 / 2.5
    d::Float64 = 1.0 / 2.5
    fr_pre::Float64 = 0.0
    ds::Float64 = 0.0
end


"""
Structure for the synapses GABA 


...
# Fields
- `dynamique_variables::state_of_synapse` = dynamical allocation of the variables
- `τ::Float64` = time constant of the synapse_derivative
- `γ::Float64` = 
- `g::Float64` = 
- `facilitation::Bool` = variable that indicates if the synapse is subject to facilitation
- `mult::Float64` = 
- `f_param`::facilitation_parameters = parameters for the facilitation of the synapse
- `d_param`::depression_parameters = parameters for the depression of the synapse
- `depression`::Bool = variable that indicates if the synapse is subject to depression
- `name::String` = Name of the synapses (TODO the structure of the names)
- `s_save::Vector{Float64}`` = Vector that saves the gating variable of the synapse
- `u_save::Vector{Float64}`` = Vector that saves the facilitaiton variable of the synapse
- `d_save::Vector{Float64}`` = Vector that saves the depression variable of the synapse
...
"""
@with_kw struct gaba_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    τ::Float64 = 5 * 0.001 # s (except for dendrites)
    γ::Float64 =  2
    g::Float64
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


"""
Structure for the synapses AMPA

...
# Fields
- `dynamique_variables::state_of_synapse` = dynamical allocation of the variables
- `τ::Float64` = time constant of the synapse_derivative
- `γ::Float64` = 
- `g::Float64` = 
- `facilitation::Bool` = variable that indicates if the synapse is subject to facilitation
- `mult::Float64` = 
- `f_param`::facilitation_parameters = parameters for the facilitation of the synapse
- `d_param`::depression_parameters = parameters for the depression of the synapse
- `depression`::Bool = variable that indicates if the synapse is subject to depression
- `name::String` = Name of the synapses (TODO the structure of the names)
- `s_save::Vector{Float64}`` = Vector that saves the gating variable of the synapse
- `u_save::Vector{Float64}`` = Vector that saves the facilitaiton variable of the synapse
- `d_save::Vector{Float64}`` = Vector that saves the depression variable of the synapse
...
"""
@with_kw struct ampa_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    τ::Float64 = 2 * 0.001 # s
    γ::Float64 = 5
    g::Float64
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



"""
Structure for the synapses NMDA

 
...
# Fields
- `dynamique_variables::state_of_synapse` = dynamical allocation of the variables
- `τ::Float64` = time constant of the synapse_derivative
- `γ::Float64` = 
- `g::Float64` = 
- `facilitation::Bool` = variable that indicates if the synapse is subject to facilitation
- `mult::Float64` = 
- `f_param`::facilitation_parameters = parameters for the facilitation of the synapse
- `d_param`::depression_parameters = parameters for the depression of the synapse
- `depression`::Bool = variable that indicates if the synapse is subject to depression
- `name::String` = Name of the synapses (TODO the structure of the names)
- `s_save::Vector{Float64}`` = Vector that saves the gating variable of the synapse
- `u_save::Vector{Float64}`` = Vector that saves the facilitaiton variable of the synapse
- `d_save::Vector{Float64}`` = Vector that saves the depression variable of the synapse
...
"""
@with_kw struct nmda_syn <: synapses
    dynamique_variables::state_of_synapse = state_of_synapse()
    γ::Float64 = 0.641 * 2
    τ::Float64 = 60 * 0.001 # s
    g::Float64
    facilitation::Bool = false
    mult::Float64 = 2.5 # account for the fact that u is below 1
    f_param::facilitation_parameters = facilitation_parameters()
    d_param::depression_parameters = depression_parameters()
    depression::Bool = false
    name::String
    s_save::Vector{Float64} = [0.5]
    u_save::Vector{Float64} = [1.0 / 2.5]
    d_save::Vector{Float64} = [1.0 / 2.5]
end


end
using .Synapses

###############################################################################################################
# Neuronal Models
###############################################################################################################
module NeuronalModels

using ..AbstractNeuronalTypes, ..Synapses, ...BasicFunctions

using Parameters,StaticArrays


export dend_sigmoid


"""
Parameters for the dendrite with the sigmoid transfer function

``I_{output} = f(I_{exc},I_{inh}) ``
with ``f`` a sigmoid function whose argument is

`` y = \\frac{I_{exc} - c_2 * I_{inh} + c_6}{c_3 * I_{inh} + c_4}``

`` I_{output} = c_1 * (-0.5 + σ(y)) + c_5   ``

"""
@with_kw struct dendrites_param_sigmoid
    c1::Float64 = 0.12
    c2::Float64 = -7.0
    c3::Float64 = -0.482
    c4::Float64 = 0.00964
    c5::Float64 = 0.19624
    c6::Float64 = 0.0
end
export dendrites_param_sigmoid

"""
Dynamical variables of the dendrite
"""
@with_kw mutable struct dynamique_dend_sigmoid
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ioutput::Float64 = 0.0
    Ibg::Float64 = 100.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0
end

"""
Complete structure for the sigmoid dendrite

...
# Fields 
- `param_c::dendrites_param_sigmoid` = param of the dendrites (transfer function)
- `list_syn_pre_gaba::Vector{gaba_syn}` =  liste of GABA synapse where the unit is the presynaptic
- `list_syn_post_gaba::Vector{gaba_syn}` = =  liste of GABA synapse where the unit is the postsynaptic
- `list_syn_pre_ampa::Vector{ampa_syn}` = liste of AMPA synapse where the unit is the presynaptic
- `list_syn_post_ampa::Vector{ampa_syn}` = liste of AMPA synapse where the unit is the postsynaptic
- `list_syn_pre_nmda::Vector{nmda_syn}` = liste of NMDA synapse where the unit is the presynaptic
- `list_syn_post_nmda::Vector{nmda_syn}` = liste of NMDA synapse where the unit is the postsynaptic
- `Ibg::Float64` = Background current to the dendrite
- `τ::Float64` = Time constant of the dendrite
- `OU_process::OU_process` = Ornstein Uhlenbeck process for the noise to the dendrite
- `name::String` = name of the dendrite
- `Iexc_save::Vector{Float64}` = Array to save Iexc to the dendrite
- `Iinh_save::Vector{Float64}`
- `Ioutput_save::Vector{Float64}` 
- `Itot_save::Vector{Float64}` 
- `Inoise_save::Vector{Float64}` 
- `dynamique_variables::dynamique_dend_sigmoid` =  dynamical variables of the dendrite
...

"""
@with_kw struct dend_sigmoid <: dendrite
    param_c::dendrites_param_sigmoid = dendrites_param_sigmoid()
    
    list_syn_pre_gaba::Vector{gaba_syn} = Vector{gaba_syn}()
    list_syn_post_gaba::Vector{gaba_syn} = Vector{gaba_syn}()

    list_syn_pre_ampa::Vector{ampa_syn} = Vector{ampa_syn}()
    list_syn_post_ampa::Vector{ampa_syn} = Vector{ampa_syn}()

    list_syn_pre_nmda::Vector{nmda_syn} = Vector{nmda_syn}()
    list_syn_post_nmda::Vector{nmda_syn} = Vector{nmda_syn}()

    Ibg::Float64 = 100.0 * 0.001
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





"""
Mutable structure for the dynamical variables of the synaptic adaptation

...
# Fields
-`sA::Float64 `
-`sA_save::Vector{Float64}`
-`τA::Float64` Time constant of adaptation (100ms by default)
-`gA::Float64` 
...

"""
@with_kw mutable struct adaptation_variables
    sA::Float64 = 0.0
    sA_save::Vector{Float64} = [0.0]
    τA::Float64 = 0.2 # seconds as for facilitation and 0.05
    gA::Float64 = -0.01 # value to test 0.5, was -0.1 before
end
export adaptation_variables


"""
Mutable structure with all the variables for the soma
"""
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


"""
Structure for the soma with one sigmoid dendrite
"""
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
    r_save::Vector{Float64} = [9.6]
    
    preferred_stim::Float64 = 0.0 
    
end
export soma_PC


"""
Mutable structure for the dynamical variables of a rectified linear neuron
"""
@with_kw mutable struct dynamique_rectified_linear <: dynamique_variables_interneurons
    r::Float64 = 0.0		
    Iinput::Float64 = 0.0
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 329.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0

end


"""
Structure of a pv cell
"""
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
    r_save::Vector{Float64} = [11.4]


end
export pv_cell


"""
Structure of a sst cell
"""
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
    r_save::Vector{Float64} = [2.94]
    
end
export sst_cell


"""
Structure of a vip cell
"""
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
    r_save::Vector{Float64} = [10.64]
end
export vip_cell

"""
Structure of a ngfc cell. For now it has the same parameters as sst and pv
"""
@with_kw struct ngfc_cell <: local_circuit_interneuron
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
    r_save::Vector{Float64} = [2.11]
end
export ngfc_cell


"""
Structure of a neural integrator
"""
 @with_kw struct neural_integrator <: neuron
    dynamique_variables::dynamique_rectified_linear = dynamique_rectified_linear()
    τ::Float64 = 0.2#0.8 
    dt::Float64 = 0.0005
    α::Float64 = 0.5 #0.9  # 0.8 before
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
    OU_process::OU_process = OU_process()
    name::String = "temp"
    r_save::Vector{Float64} = [0.43]		
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
using ..NeuronalModels


"""
Parameters of the simulations
"""
@with_kw struct simulation_parameters
    
    Tfin::Float64 = 20.0
    Tinit::Float64 = 2.0
    TISI::Float64 = 1.0    
    Tstimduration::Float64 = 0.5
    current::Dict{String, Vector{Float64}} = Dict()
    # Stimulus::Dict{String,Vector{Float64, floor(Int,(Tinit+Tfin)/dt)}} = Dict()
    
end


export parameters_bump_attractor, parameters_inter_microcircuit, parameters_interneurons, parameters_microcircuit, parameters_syn_strength_microcircuit

"""
List of all the parameters used for the simulations for the microcircuit

...
# Functions used in:
- construct_two_local_microcircuit_integrator

....
"""
@with_kw mutable struct parameters_microcircuit

    τA::Float64 = 0.1 #ms
    gA::Float64 = -0.01
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
    preferred::Float64 = 0.0
    ngfc_to_dend_depression = true #depression comes fron Abs et al 2018

end



"""
Parameters for the synaptic strengths of a microcircuit
"""
@with_kw mutable struct parameters_syn_strength_microcircuit
    gaba_sst_to_dend::Float64 = -0.09
    gaba_pv_to_soma::Float64 = -0.01
    nmda_soma_to_soma::Float64 = 0.18
    nmda_soma_to_vip::Float64 = 0.058
    gaba_sst_to_vip::Float64 = -0.1
    nmda_soma_to_sst::Float64 = 0.0435
    gaba_vip_to_sst::Float64 = -0.17
    nmda_soma_to_pv::Float64 = 0.17
    gaba_sst_to_pv::Float64 = -0.17
    gaba_pv_to_pv::Float64 = -0.18
    nmda_soma_to_int::Float64 = 0.15
    nmda_int_to_dend::Float64 = 0.4
    gaba_sst_to_ngfc::Float64 = -0.4
    gaba_ngfc_to_dend::Float64 = -0.09

end



"""
Parameters specific to interneurons. THis current version, assume full connection to NGFC (probably to be modulated by the synaptic strength)
"""
@with_kw mutable struct parameters_interneurons

    Ibg_vip::Float64 = 0.29
    Ibg_sst::Float64 = 0.25
    Ibg_pv::Float64 = 0.29
    Ibg_ngfc::Float64 = 0.25
    top_down_to_interneurons::Vector{Float64} = 0.1.*[0.66, 0.34, 0.0,0.25] 
    # with 0.01 the normalization factor of top-down

end

"""
Parameters specific to cross-talk between 2 microcircuits
"""
@with_kw struct parameters_inter_microcircuit
    cross_int_to_vip_depression::Bool = true
    cross_int_to_pv_depression::Bool = true
    cross_int_to_sst_facilitation::Bool = true
    cross_int_to_dend_depression::Bool = true
    cross_soma_to_sst_facilitation::Bool = true
end


"""
Parameters of a ring model
"""
@with_kw struct parameters_bump_attractor
    num_circuits::Int = 128
    σ::Float64 = 43.2
    Jminus::Float64 = -0.5
    Jplus::Float64 = 1.43 #nA
end




end
using .Simulations






###############################################################################################################
# Neural network
###############################################################################################################
module NeuralNetwork
using Parameters, StaticArrays
using ..AbstractNeuronalTypes, ..Simulations
using ..NeuronalModels
using ..EqDiffMethod


abstract type ring_layer end

"""
Strucutre of a microcircuit (3 interneurons, soma, dendrite and integrator)
"""
@with_kw struct microcircuit{T <: pyr_cell, D <: dendrite}
    soma::T 
    sst::sst_cell
    vip::vip_cell
    pv::pv_cell
    dend::D
    ngfc::ngfc_cell
    integrator::neural_integrator = neural_integrator()
    name::String = "microciruit"
    eq_diff_method::euler_method = euler_method()

end
export microcircuit

"""
Structure of a ring model of linear integrator
"""
@with_kw struct ring_integrator <: ring_layer
    neurons::SVector{64,neural_integrator}
    name::String = "ring_integrator"
    eq_diff_method::euler_method = euler_method()
    parameters::parameters_bump_attractor = parameters_bump_attractor()
    parameters_microcircuit::parameters_microcircuit = parameters_microcircuit()

end
export ring_integrator

"""
Structure of a ring model of microcircuit
"""
@with_kw struct ring_microcircuit{T <: pyr_cell, D <: dendrite, E <: eq_diff_method} <: ring_layer
    neurons::SVector{128,microcircuit{T,D}}
    name::String = "ring_microcircuit"
    eq_diff_method::E = euler_method()
    parameters::parameters_bump_attractor = parameters_bump_attractor()
    parameters_m::parameters_microcircuit = parameters_microcircuit()
    parameters_i::parameters_interneurons = parameters_interneurons()

end
export ring_microcircuit

"""
Structure of the full network
"""
@with_kw struct neural_network{T <: pyr_cell, D <: dendrite, E <: eq_diff_method}
    ring_integrator::ring_integrator
    list_microcircuit::ring_microcircuit{T,D,E}
    name::String = "neural_network"
    eq_diff_method::E = euler_method()
   
end
export neural_network

"""
Structure of ring layer
"""
@with_kw struct layer_bump{T <: pyr_cell, D <: dendrite, E <: eq_diff_method}
    bump_param::parameters_bump_attractor
    list_microcircuit::SVector{128,microcircuit{T,D}}
    eq_diff_method::E 
end
export layer_bump



end
using .NeuralNetwork



module get_functions

using Parameters, StaticArrays
using ..AbstractNeuronalTypes, ..Simulations
using ..NeuronalModels
using ..EqDiffMethod
using ..NeuralNetwork

"""
    get_firing_rate(neuron,time)

get the firing rate of a neuron at time t
"""
function get_firing_rate(neuron::T where T <: neuron, time::Int64)
    return neuron.r_save[time]
end


"""
    get_current_to_bump(layer_bump, time, name of the neuron)

get the current to the bump at time time
"""
function get_current_to_bump(layer_bump::layer_bump{T , D, E} where {T <: pyr_cell ,D <:dendrite, E <: eq_diff_method}, time::Int64, name::String)
    return layer_bump.sim.current[name][time]
end
export get_current_to_bump



"""
    get_firing_rate(layer_bump, time, name of the neuron)

get the firing rate of all the neurons of this type acroos the bump, at time time
"""
function get_firing_rate(layer_bump::layer_bump{T , D, E} where {T <: pyr_cell ,D <:dendrite, E <: eq_diff_method}, time::Int64, name::String)
   
    size = layer_bump.bump_param.num_circuits
   
    firing_rate = zeros(size)

    time = floor(Int,time)

    for i=1:size
        if name =="soma"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit[i].soma, time)
        elseif name =="vip"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit[i].vip, time)
        elseif name =="pv"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit[i].pv, time)
        elseif name =="sst"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit[i].sst, time)
        elseif name =="integrator"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit[i].integrator, time)
        elseif name =="ngfc"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit[i].ngfc, time)
        end

    end

    return firing_rate
end

function get_firing_rate(layer_bump::neural_network{T , D, E} where {T <: pyr_cell ,D <:dendrite, E <: eq_diff_method}, time::Float64, name::String)
   
    size = layer_bump.list_microcircuit.parameters.num_circuits
   
    firing_rate = zeros(size)

    time = floor(Int,time)

    for i=1:size
        if name =="soma"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].soma, time)
        elseif name =="vip"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].vip, time)
        elseif name =="pv"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].pv, time)
        elseif name =="sst"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].sst, time)
        elseif name =="integrator"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].integrator, time)
        elseif name =="ngfc"
            firing_rate[i] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].ngfc, time)
        end

    end

    return firing_rate
end

function get_firing_rate(layer_bump::ring_integrator, time::Float64)
   
    size = 64#layer_bump.parameters.num_circuits
   
    firing_rate = zeros(size)

    time = floor(Int,time)

    for i=1:size
     
            firing_rate[i] = get_firing_rate(layer_bump.neurons[i], time)
      

    end

    return firing_rate
end
export get_firing_rate


"""
    get_firing_rate(neuron, timerange)

get the firing of a neuron during a window of timerange    
"""
function get_firing_rate(neuron::T where T <: neuron, timerange::StepRange{Int64, Int64})
    return neuron.r_save[timerange]
end

"""
    get_firing_rate(layer_bump, timerange, name of the neuron)

get the firing rate of this neuron type acroos the bump, during a window of timerange
"""
function get_firing_rate(layer_bump::layer_bump{T , D, E} where {T <: pyr_cell ,D <:dendrite, E <: eq_diff_method}, timerange::StepRange{Int64, Int64}, name::String)
   
    size = layer_bump.bump_param.num_circuits
   
    firing_rate = zeros(size,length(timerange))

    for i=1:size
        if name =="soma"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit[i].soma, timerange)
        elseif name =="vip"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit[i].vip, timerange)
        elseif name =="pv"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit[i].pv, timerange)
        elseif name =="sst"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit[i].sst, timerange)
        elseif name =="integrator"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit[i].integrator, timerange)
        elseif name =="ngfc"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit[i].ngfc, timerange)
        end

    end

    return firing_rate
end

function get_firing_rate(layer_bump::neural_network{T , D, E} where {T <: pyr_cell ,D <:dendrite, E <: eq_diff_method}, timerange::StepRange{Int64, Int64}, name::String)
   
    size = layer_bump.list_microcircuit.parameters.num_circuits
   
    firing_rate = zeros(size,length(timerange))

    for i=1:size
        if name =="soma"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].soma, timerange)
        elseif name =="vip"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].vip, timerange)
        elseif name =="pv"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].pv, timerange)
        elseif name =="sst"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].sst, timerange)
        elseif name =="integrator"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].integrator, timerange)
        elseif name =="ngfc"
            firing_rate[i,:] = get_firing_rate(layer_bump.list_microcircuit.neurons[i].ngfc, timerange)
        end

    end

    return firing_rate
end


function get_firing_rate(layer_bump::ring_integrator, timerange::StepRange{Int64, Int64})
   
    size = 64#layer_bump.parameters.num_circuits
   
    firing_rate = zeros(size,length(timerange))

    for i=1:size
     
            firing_rate[i,:] = get_firing_rate(layer_bump.neurons[i], timerange)
       
    end

    return firing_rate
end


"""
    get_firing_rate(layer_bump, timerange, name of the neuron, indice of the microcircuit)

get the firing rate of a neuron type in a specific microcircuit, during a window of timerange
"""
function get_firing_rate(layer_bump::layer_bump{T , D, E} where {T <: pyr_cell ,D <:dendrite, E <: eq_diff_method}, timerange::StepRange{Int64, Int64}, name::String, indice::Int)
   
    if name =="soma"
        return get_firing_rate(layer_bump.list_microcircuit[indice].soma, timerange)
    elseif name =="vip"
        return get_firing_rate(layer_bump.list_microcircuit[indice].vip, timerange)
    elseif name =="pv"
        return get_firing_rate(layer_bump.list_microcircuit[indice].pv, timerange)
    elseif name =="sst"
        return get_firing_rate(layer_bump.list_microcircuit[indice].sst, timerange)
    elseif name =="integrator"
        return get_firing_rate(layer_bump.list_microcircuit[indice].integrator, timerange)
    elseif name =="ngfc"
        return get_firing_rate(layer_bump.list_microcircuit[indice].ngfc, timerange)
    end
end

function get_firing_rate(layer_bump::neural_network{T , D, E} where {T <: pyr_cell ,D <:dendrite, E <: eq_diff_method}, timerange::StepRange{Int64, Int64}, name::String, indice::Int)
   
    if name =="soma"
        return get_firing_rate(layer_bump.list_microcircuit.neurons[indice].soma, timerange)
    elseif name =="vip"
        return get_firing_rate(layer_bump.list_microcircuit.neurons[indice].vip, timerange)
    elseif name =="pv"
        return get_firing_rate(layer_bump.list_microcircuit.neurons[indice].pv, timerange)
    elseif name =="sst"
        return get_firing_rate(layer_bump.list_microcircuit.neurons[indice].sst, timerange)
    elseif name =="integrator"
        return get_firing_rate(layer_bump.list_microcircuit.neurons[indice].integrator, timerange)
    elseif name =="ngfc"
        return get_firing_rate(layer_bump.list_microcircuit.neurons[indice].ngfc, timerange)
    end
end

function get_firing_rate(layer_bump::ring_integrator, timerange::StepRange{Int64, Int64},indice::Int)
   
   
        return get_firing_rate(layer_bump.neurons[indice], timerange)
    
end
end

using .get_functions

export get_current_to_bump, get_firing_rate



end