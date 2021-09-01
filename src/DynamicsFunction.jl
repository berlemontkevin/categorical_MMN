


module DynamicsFunction

using Parameters, StaticArrays


using ..BasicFunctions

using ..NeuronalStructures

export simulation_step!,time_step

export dendrite_input_output!, update_firing!, update_syn!


export sum_input!, current_to_frequency, update_dend!, full_time_dynamics
export current_synapses!, synapse_derivative, update_adaptation!




module GeneralDynamicalFunctions

using StaticArrays


using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.attractor_network


using ...BasicFunctions

export update_fr!, update_syn_current!

####################
# Function general for neurons
#####################

"""
    Function: update_fr!(n::rectified_linear_neurons,euler_m::euler_method)

Updates  the firing rate of a rectified linear neuron. Takes the differential equation method in the arguments
"""
function update_fr!(n::rectified_linear_neurons, euler_m::euler_method)
    @. n.fr += euler_m.dt / n.τ * (-n.fr + n.Isyn )
    @. n.fr = rect_linear(n.fr)
    push!(n.fr_save, n.fr[1])

end


"""
    Function: update_syn_current!(n::neuron)

Updates  the synaptic current toward a neuron. IMPORTANT: the backgrounbd current is included in Isyn
let's note that the synaptic weights have a SIGN
"""
function update_syn_current!(n::neuron)

    n.Isyn = @MVector [0.0]
    @simd for syn in n.neurons_list
        @. n.Isyn += syn.pre_syn.fr * syn.w
    end
    @. n.Isyn += n.Ibg 
    push!(n.Isyn_save, n.Isyn[1])
end



end

using .GeneralDynamicalFunctions




module local_circuit_functions

using Parameters, StaticArrays

# using ...NeuronalStructures
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.RateDendrites
using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.attractor_network


using ..GeneralDynamicalFunctions


using ...BasicFunctions


export dendrite_input_output!, sum_input!, current_to_frequency, update_firing!, synapse_derivative, current_synapses!, update_dend!
export update_firing!
export update_syn!, update_adaptation!

function dendrite_input_output!(d::dend_sigmoid)
    @unpack c1, c2, c3, c4, c5, c6 = d.param_c
    
    @fastmath d.dynamique_variables.Ioutput = c1 * (-0.5 + sigmoid((d.dynamique_variables.Iexc - c2 * d.dynamique_variables.Iinh + c6 ) / (c3 * d.dynamique_variables.Iinh + c4))) + c5
    
   # push!(d.Ioutput_save, d.dynamique_variables.Ioutput)

end


function sum_input!(n::T where {T <: local_circuit_interneuron}, index::Int64)
 
    if n.adaptation_boolean
        update_adaptation!(n)
       
        @inbounds n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.adaptation.gA * n.adaptation.sA[1]
   #     push!(n.Itot_save , n.dynamique_variables.Itot)
    else
        @inbounds n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim 
    end
  

end


function sum_input!(n::soma_PC, index::Int64)

    if n.adaptation_boolean
        update_adaptation!(n)
      
        @inbounds  n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.den.dynamique_variables.Ioutput + n.adaptation.gA * n.adaptation.sA[1]
       # push!(n.Itot_save , n.dynamique_variables.Itot)

    else
        @inbounds  n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.den.dynamique_variables.Ioutput
    #   #  push!(n.Itot_save , n.dynamique_variables.Itot)
    end
    return
end



function sum_input!(n::neural_integrator, index::Int64)
    @inbounds n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim 
    # push!(n.Itot_save , n.dynamique_variables.Itot )
end


function current_to_frequency(neuron::pyr_cell)
    # r is the rate of the neuron
    # a = neuron.a
    # b = neuron.b
    # d = neuron.c
    
    # input_current = neuron.Itot

   @fastmath neuron.dynamique_variables.r = (neuron.a * neuron.dynamique_variables.Itot - neuron.b) / (1.0 - exp(-neuron.c * (neuron.a * neuron.dynamique_variables.Itot - neuron.b)))
   nothing
end


function current_to_frequency(neuron::neural_integrator)
    # a = neuron.a
    # b = neuron.b
    # d = neuron.c
    # for now, the integrtor is on the f-I curve
    
    @fastmath neuron.dynamique_variables.r = neuron.α .* neuron.dynamique_variables.r .+ (neuron.a * neuron.dynamique_variables.Itot - neuron.b) / (1.0 - exp(-neuron.c * (neuron.a * neuron.dynamique_variables.Itot - neuron.b)))
    
    nothing

end


function current_to_frequency(neuron::T where {T <: local_circuit_interneuron})
    # r is the rate of the neuron
    # c_I = neuron.c_I
    # r0 = neuron.r0
    # input_current = neuron.Itot

    neuron.dynamique_variables.r = rect_linear(neuron.c_I * neuron.dynamique_variables.Itot + neuron.r0)
    nothing

end



function update_firing!(n::T where {T <: neuron})
    # update firing rates following Euler law
    @fastmath  temp = n.dt / n.τ * (-n.r_save[end] + n.dynamique_variables.r)

    push!(n.r_save, n.r_save[end])

    n.r_save[end] += temp

    n.dynamique_variables.r = n.r_save[end]
    #push!(n.r_save, n.dynamique_variables.r)
    nothing
end


function synapse_derivative(s::ampa_syn)
    if s.facilitation
        @fastmath  s.dynamique_variables.u += s.f_param.dt * ((s.f_param.U - s.dynamique_variables.u) / s.f_param.τ + s.f_param.U * (1.0 - s.dynamique_variables.u) * s.dynamique_variables.fr_pre)
        # push!(s.u_save,s.dynamique_variables.u)
        @fastmath s.dynamique_variables.ds =  -s.dynamique_variables.s / s.τ + s.γ * s.mult * s.dynamique_variables.u * s.dynamique_variables.fr_pre

        s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
    elseif s.depression
        @fastmath  s.dynamique_variables.d += s.f_param.dt * ((1.0 - s.dynamique_variables.d) / s.f_param.τ - s.dynamique_variables.d * (1.0 - s.d_param.fD) * s.dynamique_variables.fr_pre)
        # push!(s.d_save,s.dynamique_variables.d)
         @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * s.mult * s.dynamique_variables.d * s.dynamique_variables.fr_pre

         s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
    else
         @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * s.dynamique_variables.fr_pre

         s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
    end
nothing
end


function synapse_derivative(s::gaba_syn)
    # paper Sean: Inhibition onto the dendrites is slower than inhibition elsewhere (ali and Thomson 2008)
    s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ  + s.γ * s.dynamique_variables.fr_pre

    s.dynamique_variables.s += s.dt * s.dynamique_variables.ds

nothing
end


function synapse_derivative(s::nmda_syn)
    if s.facilitation
        @fastmath s.dynamique_variables.u += s.f_param.dt * ((s.f_param.U - s.dynamique_variables.u) / s.f_param.τ + s.f_param.U * (1.0 - s.dynamique_variables.u) * s.dynamique_variables.fr_pre)
       # push!(s.u_save,s.dynamique_variables.u)
       @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * (1.0 - s.dynamique_variables.s) * s.dynamique_variables.u * s.mult * s.dynamique_variables.fr_pre
        
       s.dynamique_variables.s += s.dt * s.dynamique_variables.ds

    elseif s.depression
        @fastmath s.dynamique_variables.d += s.f_param.dt * ((1.0 - s.dynamique_variables.d) / s.f_param.τ - s.dynamique_variables.d * (1.0 - s.d_param.fD) * s.dynamique_variables.fr_pre)
      #  push!(s.d_save,s.dynamique_variables.d)

        @fastmath s.dynamique_variables.ds = -s.dynamique_variables.s / s.τ + s.γ * s.mult * s.dynamique_variables.d * s.dynamique_variables.fr_pre
    
        s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
    else
        @fastmath s.dynamique_variables.ds =  -s.dynamique_variables.s / s.τ + s.γ * (1.0 - s.dynamique_variables.s) * s.dynamique_variables.fr_pre
        
        s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
    end
nothing

end

function current_synapses!(n::neuron, d::Dict{String,Vector{Float64}}, index::Int64)
    # compute the sum of syn currents
    # separe in two the currents (due to the dendrites)
    # for n in ln 

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



   # push!(n.Iexc_save , n.dynamique_variables.Iexc)
   # push!(n.Iinh_save , n.dynamique_variables.Iinh)
        @inbounds  n.dynamique_variables.Istim = d[n.name][index]
   
    # end
    return

end

function add_current!(gain::Float64, s_dyn::Float64, Iexc::Float64, Iinh::Float64)

    if gain < 0.0
        Iinh += gain * s_dyn
    else
        Iexc += gain * s_dyn
    end

    return
end


function current_synapses!(n::dend_sigmoid)
   
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


    # push!(n.Iexc_save , n.dynamique_variables.Iexc)
    # push!(n.Iinh_save , n.dynamique_variables.Iinh)
    return
end


function update_dend!(d::dend_sigmoid, index::Int64)
    current_synapses!(d)
    @inbounds d.dynamique_variables.Iexc += d.OU_process.noise[index] + d.dynamique_variables.Istim + d.dynamique_variables.Ibg

    d.dynamique_variables.Itot = d.dynamique_variables.Iexc + d.dynamique_variables.Iinh
   # push!(d.Itot_save, d.dynamique_variables.Itot)
    
    dendrite_input_output!(d)
end


function update_adaptation!(n::T where T <: neuron)
    @fastmath @inbounds  @. n.adaptation.sA += n.dt * (-n.adaptation.sA / n.adaptation.τA + n.dynamique_variables.r)
  #  push!(adapt.sA_save, adapt.sA)

end


function update_syn!(n::T where T <: neuron)

    temp = n.dynamique_variables.r
    for s in n.list_syn_pre_ampa
        s.dynamique_variables.fr_pre = temp
    end

    for s in n.list_syn_pre_gaba
        s.dynamique_variables.fr_pre = temp
    end
    
    for s in n.list_syn_pre_nmda
        s.dynamique_variables.fr_pre = temp
    end
    return
end

end

using .local_circuit_functions

module time_dynamics


using ...NeuronalStructures.attractor_network

using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using Parameters, StaticArrays
using ...NeuronalStructures.AbstractNeuronalTypes

using ...NeuronalStructures.RateDendrites
using ...NeuronalStructures.local_circuit

using ..local_circuit_functions
# using ..attractor_network_functions

using ...BasicFunctions


export time_step, full_time_dynamics


### TODO: change the funciton to update the variable s of the synapses


function time_step(l_c::Vector{microcircuit}, sim::simulation_parameters, d::Dict{String,Vector{Float64}}, index::Int64)

    for c in l_c
        
            for n in c.list_pv
                @simd for s in n.list_syn_post_ampa
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
                #     push!(s.s_save , s.dynamique_variables.s)

                end
                @simd for s in n.list_syn_post_gaba
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end

            
                @simd for s in n.list_syn_post_nmda
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end
            end
            for n in c.list_dend
                @simd for s in n.list_syn_post_ampa
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
                #     push!(s.s_save , s.dynamique_variables.s)

                end
                @simd for s in n.list_syn_post_gaba
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end

            
                @simd for s in n.list_syn_post_nmda
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end
            end
            for n in c.list_sst
                @simd for s in n.list_syn_post_ampa
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
                #     push!(s.s_save , s.dynamique_variables.s)

                end
                @simd for s in n.list_syn_post_gaba
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end

            
                @simd for s in n.list_syn_post_nmda
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end
            end
            for n in c.list_vip
                @simd for s in n.list_syn_post_ampa
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
                #     push!(s.s_save , s.dynamique_variables.s)

                end
                @simd for s in n.list_syn_post_gaba
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end

            
                @simd for s in n.list_syn_post_nmda
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end
            end
            for n in c.list_soma
                @simd for s in n.list_syn_post_ampa
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
                #     push!(s.s_save , s.dynamique_variables.s)

                end
                @simd for s in n.list_syn_post_gaba
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end

            
                @simd for s in n.list_syn_post_nmda
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end
            end
            for n in c.list_integrator
                @simd for s in n.list_syn_post_ampa
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
                #     push!(s.s_save , s.dynamique_variables.s)

                end
                @simd for s in n.list_syn_post_gaba
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end

            
                @simd for s in n.list_syn_post_nmda
                    synapse_derivative(s)
                    # s.dynamique_variables.s += s.dt * s.dynamique_variables.ds
               #     push!(s.s_save , s.dynamique_variables.s)

                end
            end
        


    # update dend
    for dend in c.list_dend
        update_dend!(dend, index)
    end
    # update neurons
    for temp_n in c.list_soma
        current_synapses!(temp_n, d, index)
    end
    for temp_n in c.list_vip
        current_synapses!(temp_n, d, index)
    end
    for temp_n in c.list_sst
        current_synapses!(temp_n, d, index)
    end
    for temp_n in c.list_pv
        current_synapses!(temp_n, d, index)
    end
    for temp_n in c.list_integrator
        current_synapses!(temp_n, d, index)
    end
    
    for temp_n in c.list_soma
        sum_input!(temp_n, index)
    end
    for temp_n in c.list_sst
        sum_input!(temp_n, index)
    end
    for temp_n in c.list_vip
        sum_input!(temp_n, index)
    end
    for temp_n in c.list_pv
        sum_input!(temp_n, index)
    end
    for temp_n in c.list_integrator
        sum_input!(temp_n, index)
    end
     


        for temp_n in c.list_soma
            current_to_frequency(temp_n)
        end
        for temp_n in c.list_vip
            current_to_frequency(temp_n)
        end
        for temp_n in c.list_sst
            current_to_frequency(temp_n)
        end
        for temp_n in c.list_pv
            current_to_frequency(temp_n)
        end
        for temp_n in c.list_integrator
            current_to_frequency(temp_n)
        end
        

        for temp_n in c.list_soma
            update_firing!(temp_n)
        end
        
        for temp_n in c.list_vip
            update_firing!(temp_n)
        end
        
        for temp_n in c.list_sst
            update_firing!(temp_n)
        end
        
        for temp_n in c.list_pv
            update_firing!(temp_n)
        end
        
        for temp_n in c.list_integrator
            update_firing!(temp_n)
        end


        for temp_n in c.list_soma
            update_syn!(temp_n)
        end
        for temp_n in c.list_vip
            update_syn!(temp_n)
        end
        for temp_n in c.list_sst
            update_syn!(temp_n)
        end
        for temp_n in c.list_pv
            update_syn!(temp_n)
        end
        for temp_n in c.list_integrator
            update_syn!(temp_n)
        end
 

    end
    return
end


function full_time_dynamics(l_c::Vector{microcircuit}, sim::simulation_parameters, d::Dict{String,Vector{Float64}})
    # TODO: to write better. compute the whole time dynamics of the microcuit
    list_time = 0.0:sim.dt:(sim.dt * length(d[l_c[1].list_soma[1].name]))
   

    for index = 2:length(list_time[2:end]) 
         time_step(l_c, sim, d, index)
    end
    # TODO: do a nice return function for this dynamics
    return
end

# TODO
   
end
using .time_dynamics


end