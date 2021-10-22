


module DynamicsFunction

using Parameters, StaticArrays


using ..BasicFunctions

using ..NeuronalStructures



export simulation_step!,time_step!

export dendrite_input_output!, update_firing!, update_syn!


export sum_input!, current_to_frequency!, update_dend!, full_time_dynamics!
export current_synapses!, synapse_derivative!, update_adaptation!




module DendritesFunctions

using StaticArrays, Parameters
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.Synapses
using ...NeuronalStructures.Simulations
using ...BasicFunctions

"""
    dendrite_input_output!(dend)

Function that compute the output current of the dendrite with a sigmoid function

Equations are found in [`dend_sigmoid`](@ref)

"""
function dendrite_input_output!(d::dend_sigmoid)
    @unpack c1, c2, c3, c4, c5, c6 = d.param_c
    
    @fastmath d.dynamique_variables.Ioutput = c1 * (-0.5 + sigmoid((d.dynamique_variables.Iexc - c2 * d.dynamique_variables.Iinh + c6 ) / (c3 * d.dynamique_variables.Iinh + c4))) + c5
    
   push!(d.Ioutput_save, d.dynamique_variables.Ioutput)

end
export dendrite_input_output!


end
using .DendritesFunctions



module NeuronalFunctions
using StaticArrays, Parameters
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.Synapses
using ...NeuronalStructures.Simulations

using ...BasicFunctions

export sum_input!, current_to_frequency!, update_adaptation!

"""
    update_adaptation!(neuron, euler_method)

Update the adaptation variable of a neuron according to:

''\\frac{ds_A}{dt} = -s_A/\\tau_A + r  ''
"""
function update_adaptation!(n::T where T <: neuron, euler::euler_method)
    @fastmath   n.adaptation.sA += euler.dt * (-n.adaptation.sA / n.adaptation.τA + n.dynamique_variables.r)
    push!(n.adaptation.sA_save, n.adaptation.sA)
    nothing
end

"""
    sum_input!(interneuron, time, euler_method)

Update the total current arriving into an interneuron at time `time`. Apply the adaptation before if needed.

The total current is updated into `neuron.dynamique_variables.Itot`
"""
function sum_input!(n::T where {T <: local_circuit_interneuron}, index::Int64, euler::euler_method)
 
    if n.adaptation_boolean
        update_adaptation!(n, euler)
       
        @inbounds n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.adaptation.gA * n.adaptation.sA[1]
       push!(n.Itot_save , n.dynamique_variables.Itot)
    else
        @inbounds n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim 
    end
  

end

"""
    sum_input!(soma, time, euler_method)
    
Apply to a soma of a PC
"""
function sum_input!(n::soma_PC, index::Int64, euler::euler_method)

    if n.adaptation_boolean
        update_adaptation!(n,euler)
      
        @inbounds  n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.den.dynamique_variables.Ioutput + n.adaptation.gA * n.adaptation.sA[1]
         push!(n.Itot_save , n.dynamique_variables.Itot)

    else
        @inbounds  n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.den.dynamique_variables.Ioutput
        push!(n.Itot_save , n.dynamique_variables.Itot)
    end
    return
end

"""
    sum_input!(neural_integrator, time)

Apply to a neural integrator
"""
function sum_input!(n::neural_integrator, index::Int64)
    n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim 
    push!(n.Itot_save , n.dynamique_variables.Itot )
end


"""
    current_to_frequency!(soma)

Update the firing rate of the neuron according to the input current and the transfer function

``r = \\frac{a*I_{tot} - b}{1.0 - \\exp{-c*(a*I_{tot} - b)}}``
"""
function current_to_frequency!(neuron::soma_PC)

   @fastmath neuron.dynamique_variables.r = (neuron.a * neuron.dynamique_variables.Itot - neuron.b) / (1.0 - exp(-neuron.c * (neuron.a * neuron.dynamique_variables.Itot - neuron.b)))
   nothing
end

"""
    current_to_frequency!(neuronal_integrator)

``r = \\frac{a*I_{tot} - b}{1.0 - \\exp{-c*(a*I_{tot} - b)}}``
"""
function current_to_frequency!(neuron::neural_integrator)
    
    @fastmath neuron.dynamique_variables.r = neuron.α .* neuron.dynamique_variables.r .+ (neuron.a * neuron.dynamique_variables.Itot - neuron.b) / (1.0 - exp(-neuron.c * (neuron.a * neuron.dynamique_variables.Itot - neuron.b)))
    
    nothing

end

"""
    current_to_frequency!(interneuron)

``r = [c_I * I_{tot} + r-0]_+``
"""
function current_to_frequency!(neuron::T where {T <: local_circuit_interneuron})
 
    neuron.dynamique_variables.r = rect_linear(neuron.c_I * neuron.dynamique_variables.Itot + neuron.r0)
    nothing

end




end
using .NeuronalFunctions


module SynapsesFunctions

using Parameters, StaticArrays



using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.Synapses
using ...NeuronalStructures.Simulations


using ..DendritesFunctions
using ..NeuronalFunctions

using ...BasicFunctions

export synapse_derivative!, current_synapses!

"""
    synapse_derivative!(ampa_syn, euler_method)

Update the gating variable of the synapse according ot the dynamical equations following Euler method

**If Facilitation**

``\\frac{du}{dt} = (U -u)/\\tau_F - U*(1-u) r_{pre}  ``

``\\frac{ds}{dt} = -s/τ + \\alpha γ u r_{pre}     ``

with α a multiplicator factor.

**If Depression**

``\\frac{dd}{dt} = (1 -d)/\\tau_D - d*(1-f_D) r_{pre}  ``

``\\frac{ds}{dt} = -s/τ + \\alpha γ d r_{pre}     ``

with α a multiplicator factor.

**Otherwise** 

``\\frac{ds}{dt} = -s/τ + γ r_{pre}     ``
"""
function synapse_derivative!(s::ampa_syn, euler::euler_method)
    if s.facilitation
        @fastmath  s.dynamique_variables.u += euler.dt * ((s.f_param.U - s.dynamique_variables.u) / s.f_param.τ + s.f_param.U * (1.0 - s.dynamique_variables.u) * s.dynamique_variables.fr_pre)
        push!(s.u_save,s.dynamique_variables.u)
        @fastmath s.dynamique_variables.ds =  -s.dynamique_variables.s / s.τ + s.γ * s.mult * s.dynamique_variables.u * s.dynamique_variables.fr_pre

        s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds
    elseif s.depression
        @fastmath  s.dynamique_variables.d += euler.dt * ((1.0 - s.dynamique_variables.d) / s.f_param.τ - s.dynamique_variables.d * (1.0 - s.d_param.fD) * s.dynamique_variables.fr_pre)
        push!(s.d_save,s.dynamique_variables.d)
         @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * s.mult * s.dynamique_variables.d * s.dynamique_variables.fr_pre

         s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds
    else
         @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * s.dynamique_variables.fr_pre

         s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds
    end
    push!(s.s_save,s.dynamique_variables.s)
nothing
end

"""
    synapse_derivative!(gaba_syn, euler_method)

Update the gating variable of the synapse according ot the dynamical equations following Euler method

**If Facilitation**

``\\frac{du}{dt} = (U -u)/\\tau_F - U*(1-u) r_{pre}  ``

``\\frac{ds}{dt} = -s/τ + \\alpha γ u r_{pre}     ``

with α a multiplicator factor.

**If Depression**

``\\frac{dd}{dt} = (1 -d)/\\tau_D - d*(1-f_D) r_{pre}  ``

``\\frac{ds}{dt} = -s/τ + \\alpha γ d r_{pre}     ``

with α a multiplicator factor.

**Otherwise** 

``\\frac{ds}{dt} = -s/τ + γ r_{pre}     ``
"""
function synapse_derivative!(s::gaba_syn, euler::euler_method)
    if s.facilitation
        @fastmath  s.dynamique_variables.u += euler.dt * ((s.f_param.U - s.dynamique_variables.u) / s.f_param.τ + s.f_param.U * (1.0 - s.dynamique_variables.u) * s.dynamique_variables.fr_pre)
        push!(s.u_save,s.dynamique_variables.u)
        @fastmath s.dynamique_variables.ds =  -s.dynamique_variables.s / s.τ + s.γ * s.mult * s.dynamique_variables.u * s.dynamique_variables.fr_pre

        s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds
    elseif s.depression
        @fastmath  s.dynamique_variables.d += euler.dt * ((1.0 - s.dynamique_variables.d) / s.f_param.τ - s.dynamique_variables.d * (1.0 - s.d_param.fD) * s.dynamique_variables.fr_pre)
        push!(s.d_save,s.dynamique_variables.d)
         @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * s.mult * s.dynamique_variables.d * s.dynamique_variables.fr_pre

         s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds
    else
         @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * s.dynamique_variables.fr_pre

         s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds
    end
    push!(s.s_save,s.dynamique_variables.s)
nothing
end

"""
    synapse_derivative!(nmda_syn, euler_method)

Update the gating variable of the synapse according ot the dynamical equations following Euler method

**If Facilitation**

``\\frac{du}{dt} = (U -u)/\\tau_F - U*(1-u) r_{pre}  ``

``\\frac{ds}{dt} = -s/τ + \\alpha γ u (1 - s) r_{pre}     ``

with α a multiplicator factor.

**If Depression**

``\\frac{dd}{dt} = (1 -d)/\\tau_D - d*(1-f_D) r_{pre}  ``

``\\frac{ds}{dt} = -s/τ + \\alpha γ d (1 - s) r_{pre}     ``

with α a multiplicator factor.

**Otherwise** 

``\\frac{ds}{dt} = -s/τ + γ (1 - s) r_{pre}     ``
"""
function synapse_derivative!(s::nmda_syn, euler::euler_method)
    if s.facilitation
        @fastmath s.dynamique_variables.u += euler.dt * ((s.f_param.U - s.dynamique_variables.u) / s.f_param.τ + s.f_param.U * (1.0 - s.dynamique_variables.u) * s.dynamique_variables.fr_pre)
       push!(s.u_save,s.dynamique_variables.u)
       @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * (1.0 - s.dynamique_variables.s) * s.dynamique_variables.u * s.mult * s.dynamique_variables.fr_pre
        
       s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds

    elseif s.depression
        @fastmath s.dynamique_variables.d += euler.dt * ((1.0 - s.dynamique_variables.d) / s.f_param.τ - s.dynamique_variables.d * (1.0 - s.d_param.fD) * s.dynamique_variables.fr_pre)
        push!(s.d_save,s.dynamique_variables.d)

        @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * (1.0 - s.dynamique_variables.s)* s.mult * s.dynamique_variables.d * s.dynamique_variables.fr_pre
    
        s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds
    else
        @fastmath s.dynamique_variables.ds =  -s.dynamique_variables.s / s.τ + s.γ * (1.0 - s.dynamique_variables.s) * s.dynamique_variables.fr_pre
        
        s.dynamique_variables.s += euler.dt * s.dynamique_variables.ds
    end
    push!(s.s_save,s.dynamique_variables.s)

    nothing

end

"""
    synapse_derivation!(neuron, euler_method)

Update all the synapse of a neuron
"""
function synapse_derivative!(n::neuron, euler::euler_method)

    @simd for s in n.list_syn_post_ampa
        synapse_derivative!(s, euler)
   
    end
    @simd for s in n.list_syn_post_gaba
        synapse_derivative!(s, euler)
    
    end

    @simd for s in n.list_syn_post_nmda
        synapse_derivative!(s, euler)
    
    end
    nothing
end

"""
    synapse_derivative!(microcircuit{soma_PC,dend_sigmoid}, euler_method)

Update all the synapse of a microcircuit
"""
function synapse_derivative!(lc::microcircuit{soma_PC,dend_sigmoid}, euler::euler_method)
    synapse_derivative!(lc.pv, euler)
    synapse_derivative!(lc.dend, euler)

    synapse_derivative!(lc.vip, euler)

    synapse_derivative!(lc.sst, euler)
    synapse_derivative!(lc.soma, euler)
    synapse_derivative!(lc.integrator, euler)
    nothing
end


"""
    current_synapses!(interneuron, dict_current, time_index)

Compute the sum of the synaptic current as well as the stimulus input from the matching name in the dicitonnary
"""
function current_synapses!(n::T where {T <: local_circuit_interneuron}, d::Dict{String,Vector{Float64}}, index::Int64)
        n.dynamique_variables.Iexc = 0.0
        n.dynamique_variables.Iinh = 0.0

        @simd for s in n.list_syn_post_nmda
            @fastmath n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
        end



        @simd for s in n.list_syn_post_gaba
            @fastmath n.dynamique_variables.Iinh += s.g * s.dynamique_variables.s
   
        end

 


        @simd for s in n.list_syn_post_ampa
            @fastmath n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
      
        end



    push!(n.Iexc_save , n.dynamique_variables.Iexc)
    push!(n.Iinh_save , n.dynamique_variables.Iinh)
    n.dynamique_variables.Istim = d[n.name][index]

    nothing

end

"""
    current_synapses!(soma_PC, dict_current, time_index)
"""
function current_synapses!(n::soma_PC, d::Dict{String,Vector{Float64}}, index::Int64)
        n.dynamique_variables.Iexc = 0.0
        n.dynamique_variables.Iinh = 0.0

        @simd for s in n.list_syn_post_nmda
            @fastmath n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
        end



        @simd for s in n.list_syn_post_gaba
            @fastmath n.dynamique_variables.Iinh += s.g * s.dynamique_variables.s
   
        end

 


        @simd for s in n.list_syn_post_ampa
            @fastmath n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
      
        end



    push!(n.Iexc_save , n.dynamique_variables.Iexc)
    push!(n.Iinh_save , n.dynamique_variables.Iinh)
    @inbounds  n.dynamique_variables.Istim = d[n.name][index]
   
    
    nothing

end

"""
    current_synapses!(neural_integrator, dict_current, time_index)
"""
function current_synapses!(n::neural_integrator, d::Dict{String,Vector{Float64}}, index::Int64)
        n.dynamique_variables.Iexc = 0.0
        n.dynamique_variables.Iinh = 0.0

        @simd for s in n.list_syn_post_nmda
            @fastmath n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
        end



        @simd for s in n.list_syn_post_gaba
            @fastmath n.dynamique_variables.Iinh += s.g * s.dynamique_variables.s
   
        end

        @simd for s in n.list_syn_post_ampa
            @fastmath n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
      
        end



    push!(n.Iexc_save , n.dynamique_variables.Iexc)
    push!(n.Iinh_save , n.dynamique_variables.Iinh)
    n.dynamique_variables.Istim = d[n.name][index]
   
   nothing
end

"""
    current_synapses!(microcircuit, dict_current, time_index)

Apply 'current_synapses!' on all the neurons of a microcircuit
"""
function current_synapses!(lc::microcircuit{soma_PC, dend_sigmoid}, d::Dict{String,Vector{Float64}}, index::Int64)
   
    current_synapses!(lc.soma,d,index)
    current_synapses!(lc.sst,d,index)
    current_synapses!(lc.pv,d,index)
    current_synapses!(lc.vip,d,index)
    current_synapses!(lc.integrator,d,index)

    nothing

end

"""
    current_synapses!(dend_sigmoid, dict_current, time_index)
"""
function current_synapses!(n::dend_sigmoid, d::Dict{String,Vector{Float64}}, index::Int64)
   
    n.dynamique_variables.Iexc = 0.0
    n.dynamique_variables.Iinh = 0.0

    @simd for s in n.list_syn_post_nmda
        n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
    end
   

    @simd for s in n.list_syn_post_gaba
        n.dynamique_variables.Iinh += s.g * s.dynamique_variables.s
    end

    @simd for s in n.list_syn_post_ampa
        n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
    end


    push!(n.Iexc_save , n.dynamique_variables.Iexc)
    push!(n.Iinh_save , n.dynamique_variables.Iinh)
    n.dynamique_variables.Istim = d[n.name][index]
   
   nothing
end

end
using .SynapsesFunctions


module NumericalIntegration

using Parameters, StaticArrays

using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.Synapses
using ...NeuronalStructures.Simulations


using ..DendritesFunctions
using ..NeuronalFunctions
using ..SynapsesFunctions

using ...BasicFunctions


export  update_firing!, synapse_derivative!, current_synapses!, update_dend!
export update_firing!
export update_syn!, update_adaptation!



"""
    update_firing!(neuron, euler_method)

Update the firing rate of a neuron according to Euler method

Note that the function first update the fequency according to [`current_to_frequency!`](@ref) then save the result.
"""
function update_firing!(n::T where {T <: neuron}, euler::euler_method)
    current_to_frequency!(n)
    @fastmath  temp = euler.dt / n.τ * (-n.r_save[end] + n.dynamique_variables.r)
    push!(n.r_save, n.r_save[end])
    n.r_save[end] += temp
    n.dynamique_variables.r = n.r_save[end]
    nothing
end


"""
    update_dend!(dend_sigmoid, dict_currentm time_index)

Update the synaptic current towards the dendrite, the total current and the output current that the dendrite generates.
"""
function update_dend!(d::dend_sigmoid, dic::Dict{String,Vector{Float64}}, index::Int64)
    current_synapses!(d, dic, index)
    @inbounds d.dynamique_variables.Iexc += d.OU_process.noise[index] + d.dynamique_variables.Istim + d.dynamique_variables.Ibg

    d.dynamique_variables.Itot = d.dynamique_variables.Iexc + d.dynamique_variables.Iinh
    push!(d.Itot_save, d.dynamique_variables.Itot)
    
    dendrite_input_output!(d)
    nothing
end



"""
    update_syn!(neuron)

Update the presynaptic firing rates of all synapses of a neuron
"""
function update_syn!(n::T where T <: neuron)
    for s in n.list_syn_pre_ampa
        s.dynamique_variables.fr_pre = n.dynamique_variables.r
    end

    for s in n.list_syn_pre_gaba
        s.dynamique_variables.fr_pre = n.dynamique_variables.r
    end
    
    for s in n.list_syn_pre_nmda
        s.dynamique_variables.fr_pre = n.dynamique_variables.r
    end
    nothing
end

"""
    time_step!(microcircuit{soma_PC,dend_sigmoid}, simulation_parameters, time_index)

Update all the variables of a microcircuit for a time step
"""
function time_step!(c::microcircuit{soma_PC,dend_sigmoid}, sim::simulation_parameters, index::Int64)

    synapse_derivative!(c, c.eq_diff_method)
    update_dend!(c.dend, sim.current, index)

    current_synapses!(c, sim.current, index)
    
    
    sum_input!(c.soma, index, c.eq_diff_method)
    sum_input!(c.sst, index, c.eq_diff_method)
    sum_input!(c.vip, index, c.eq_diff_method)
    sum_input!(c.pv, index, c.eq_diff_method)
    sum_input!(c.integrator, index)
    
     


    current_to_frequency!(c.soma)
    current_to_frequency!(c.vip)
    current_to_frequency!(c.sst)
    current_to_frequency!(c.pv)
    current_to_frequency!(c.integrator)
        
        

    update_firing!(c.soma, c.eq_diff_method)
    update_firing!(c.vip, c.eq_diff_method)
    update_firing!(c.sst, c.eq_diff_method)
    update_firing!(c.pv, c.eq_diff_method)
    update_firing!(c.integrator, c.eq_diff_method)
        


    update_syn!(c.soma)
    update_syn!(c.vip)
    update_syn!(c.sst)
    update_syn!(c.pv)
    update_syn!(c.integrator)
    nothing
end

"""
    time_step!(SVector{128,microcircuit{soma_PC,dend_sigmoid}}, simulation_parameters, time_index)

    Update all the variables of a ring model for a time step
"""
function time_step!(l_c::SVector{128,microcircuit{soma_PC,dend_sigmoid}}, sim::simulation_parameters, index::Int64)

    for c in l_c
        time_step!(c, sim, index)
    end
    nothing
end
export time_step!


"""
    full_time_dynamics!(layer_bump{soma_PC}, simulation_parameters)

Compute the dynamics on the full time of the simulaiton
"""
function full_time_dynamics!(l_c::layer_bump{soma_PC}, sim::simulation_parameters)

    list_time = 0.0:l_c.eq_diff_method.dt:(l_c.eq_diff_method.dt * length(d[l_c.list_microcircuit[1].soma.name]))
    for index = 2:length(list_time[1:end-1]) 
         time_step!(l_c.list_microcircuit, sim, index)
    end
    nothing

end
export full_time_dynamics!

"""
    full_time_dynamics!(microcircuit, simualtion_parameters)
"""
function full_time_dynamics!(l_c::Vector{microcircuit{soma_PC, dend_sigmoid}}, sim::simulation_parameters)

    list_time = 0.0:l_c[1].eq_diff_method.dt:(l_c[1].eq_diff_method.dt * length(sim.current[l_c[1].soma.name]))
    for index = 2:length(list_time[1:end-1]) 
    for i=1:length(l_c)
         time_step!(l_c[i], sim, index)
    end
end
    nothing

end

end
using .NumericalIntegration


end