


module DynamicsFunction

using Parameters, StaticArrays



using ..BasicFunctions

export simulation_step!,time_step

export dendrite_input_output!, update_firing!


export sum_input!, current_to_frequency!, update_dend!, full_time_dynamics
export current_synapses!

module GeneralDynamicalFunctions

using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.attractor_network

using StaticArrays

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
    # if euler_m.record
        @. n.fr += euler_m.dt / n.τ * (-n.fr + n.Isyn )
        @. n.fr = rect_linear(n.fr)
    # else
    #     n.fr[end] += euler_m.dt / n.τ * (-n.fr + n.Isyn[end] )
    #     n.fr[end] = rect_linear(n.fr[end])
    # end
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
    push!(n.Isyn_save,n.Isyn[1])
end



end

using .GeneralDynamicalFunctions




module local_circuit_functions

using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using Parameters, StaticArrays
using ...NeuronalStructures.AbstractNeuronalTypes

using ...NeuronalStructures.RateDendrites
using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.attractor_network


using ..GeneralDynamicalFunctions


using ...BasicFunctions


export dendrite_input_output!, sum_input!, current_to_frequency!, update_firing!, synapse_derivative, current_synapses!, update_dend!
export update_firing!


function dendrite_input_output!(d::dend_Sean2020)
    # dendridic function
    @unpack  c1,c2,c3,c4,c5,c6 =  d.param_c 
    @fastmath β = c5*exp(-d.Iinh[end]/c6)

    @fastmath push!(d.Ioutput , c1*(tanh((d.Iexc[end] +c3*d.Iinh[end] + c4)/β)) + c2)
end

#TODO : struct d'une dendrite prenant en argument le type de dendrite

function dendrite_input_output!(d::dend_sigmoid)
    @unpack c1, c2, c3, c4, c5, c6 = d.param_c
 #   @fastmath y = (d.dynamique_variables.Iexc .- c2.*d.dynamique_variables.Iinh .+ c6 )/(c3.*d.dynamique_variables.Iinh .+ c4)
    
    @fastmath d.dynamique_variables.Ioutput = c1.*(-0.5 .+ sigmoid((d.dynamique_variables.Iexc .- c2.*d.dynamique_variables.Iinh .+ c6 )/(c3.*d.dynamique_variables.Iinh .+ c4))).+c5
    
   # push!(d.Ioutput_save, d.dynamique_variables.Ioutput)



end


# need an abstract type for the neurons
function sum_input!(n::T where {T<: local_circuit_interneuron_without_adaptation} )
    @inbounds push!(n.Inoise, n.OU_process.noise[end]) #TODO pre-create noise terms

        n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.dynamique_variables.Inoise + n.dynamique_variables.Ibg + n.dynamique_variables.Istim 
     #   push!(n.Itot_save , n.dynamique_variables.Itot)
end

function sum_input!(n::T where {T<: local_circuit_interneuron_without_adaptation}, index::Int64 )
  #  push!(n.Inoise, n.OU_process.noise[index]) #TODO pre-create noise terms

  @inbounds n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim 
     #   push!(n.Itot_save , n.dynamique_variables.Itot)
end


# function sum_input!(n::T where {T<: local_circuit_interneuron_with_adaptation} )
#    push!(n.Inoise, n.OU_process.noise[end]) #TODO pre-create noise terms

    
#         update_adaptation!(n, n.dt)
#         Iadapt = n.adaptation.gA * n.adaptation.sA[1]
#         n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.dynamique_variables.Inoise + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + Iadapt
#     #    push!(n.Itot_save , n.dynamique_variables.Itot)

  

# end

function sum_input!(n::T where {T<: local_circuit_interneuron_with_adaptation}, index::Int64 )
  #  push!(n.Inoise, n.OU_process.noise[index]) #TODO pre-create noise terms

    
        update_adaptation!(n, n.dt)
       
     #   Iadapt = 
        @inbounds n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.adaptation.gA * n.adaptation.sA[1]
   #     push!(n.Itot_save , n.dynamique_variables.Itot)

  

end

function sum_input!(n::soma_Sean2020)
    push!(n.Inoise, n.OU_process.noise[end])
    n.Itot = n.Iinh + n.Iexc + n.Inoise[end] + n.Ibg + n.Istim + n.den.Ioutput
    push!(n.Itot_save , n.dynamique_variables.Itot)
end


function sum_input!(n::soma_PC)
    push!(n.Inoise, n.OU_process.noise[end])

    if n.adaptation_boolean
        update_adaptation!(n, n.dt)
     
        #Iadapt = n.adaptation.gA * n.adaptation.sA
        n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.dynamique_variables.Inoise + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.den.dynamique_variables.Ioutput + n.adaptation.gA * n.adaptation.sA[1]
      #  push!(n.Itot_save , n.dynamique_variables.Itot)

    else
        n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.dynamique_variables.Inoise + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.den.dynamique_variables.Ioutput
     #   push!(n.Itot_save , n.dynamique_variables.Itot)
    end
end

function sum_input!(n::soma_PC,index::Int64)
   # push!(n.Inoise, n.OU_process.noise[index])

    #if n.adaptation_boolean
        update_adaptation!(n, n.dt)
      
        #Iadapt = n.adaptation.gA * n.adaptation.sA
        @inbounds  n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.den.dynamique_variables.Ioutput + n.adaptation.gA * n.adaptation.sA[1]
       # push!(n.Itot_save , n.dynamique_variables.Itot)

    # else
    #     @inbounds  n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim + n.den.dynamique_variables.Ioutput
    #   #  push!(n.Itot_save , n.dynamique_variables.Itot)
    # end
    return
end



function sum_input!(n::neural_integrator)
    n.dynamique_variables.Inoise = n.OU_process.noise[end]
    n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.dynamique_variables.Inoise + n.dynamique_variables.Ibg + n.dynamique_variables.Istim 
 #   push!(n.Itot_save , n.dynamique_variables.Itot )
end


function sum_input!(n::neural_integrator, index::Int64)
    #push!(n.Inoise, n.OU_process.noise[index])
    @inbounds n.dynamique_variables.Itot = n.dynamique_variables.Iinh + n.dynamique_variables.Iexc + n.OU_process.noise[index] + n.dynamique_variables.Ibg + n.dynamique_variables.Istim 
    #push!(n.Itot_save , n.dynamique_variables.Itot )
end


function current_to_frequency!(neuron::pyr_cell)
    # r is the rate of the neuron
    #a = neuron.a
    #b = neuron.b
    #d = neuron.c
    
    #input_current = neuron.Itot

    return @fastmath (neuron.a*neuron.dynamique_variables.Itot - neuron.b)/(1.0 - exp(-neuron.c*(neuron.a*neuron.dynamique_variables.Itot - neuron.b)))

end


function current_to_frequency!(neuron::neural_integrator)
    #a = neuron.a
    #b = neuron.b
    #d = neuron.c
    # for now, the integrtor is on the f-I curve
    
    return  @fastmath neuron.α.*neuron.dynamique_variables.r .+ (neuron.a*neuron.dynamique_variables.Itot - neuron.b)/(1.0 - exp(-neuron.c*(neuron.a*neuron.dynamique_variables.Itot - neuron.b)))
    
end


function current_to_frequency!(neuron::T where {T<: local_circuit_interneuron})
    # r is the rate of the neuron
    #c_I = neuron.c_I
    #r0 = neuron.r0
    #input_current = neuron.Itot

    return  rect_linear(neuron.c_I*neuron.dynamique_variables.Itot + neuron.r0)
end



function update_firing!(n::neuron)
    # update firing rates following Euler law
    @fastmath  n.dynamique_variables.r += n.dt/n.τ*(-n.dynamique_variables.r + current_to_frequency!(n))
    push!(n.r_save,n.dynamique_variables.r)
end


function synapse_derivative(s::ampa_syn)
    if s.facilitation
        @fastmath  s.dynamique_variables.u += s.f_param.dt*((s.f_param.U - s.dynamique_variables.u)/s.f_param.τ + s.f_param.U*(1.0 - s.dynamique_variables.u)*s.neuron_pre.dynamique_variables.r)
        #push!(s.u_save,s.dynamique_variables.u)
        return  -s.dynamique_variables.s/s.τ + s.γ*s.mult*s.dynamique_variables.u*s.neuron_pre.dynamique_variables.r
    elseif s.depression
        @fastmath  s.dynamique_variables.d += s.f_param.dt*((1.0 - s.dynamique_variables.d)/s.f_param.τ - s.dynamique_variables.d*(1.0 - s.d_param.fD)*s.neuron_pre.dynamique_variables.r)
        #push!(s.d_save,s.dynamique_variables.d)
        return  -s.dynamique_variables.s/s.τ + s.γ*s.mult*s.dynamique_variables.d*s.neuron_pre.dynamique_variables.r
    else
        return  -s.dynamique_variables.s/s.τ + s.γ*s.neuron_pre.dynamique_variables.r
    end

end


function synapse_derivative(s::gaba_syn)
    #paper Sean: Inhibition onto the dendrites is slower than inhibition elsewhere (ali and Thomson 2008)
    return  @fastmath -s.dynamique_variables.s/s.τ + s.γ*s.neuron_pre.dynamique_variables.r

end


function synapse_derivative(s::nmda_syn)
    if s.facilitation
        @fastmath s.dynamique_variables.u += s.f_param.dt*((s.f_param.U - s.dynamique_variables.u)/s.f_param.τ + s.f_param.U*(1.0 - s.dynamique_variables.u)*s.neuron_pre.dynamique_variables.r)
       # push!(s.u_save,s.dynamique_variables.u)
        
        return  @fastmath -s.dynamique_variables.s/s.τ + s.γ*(1.0 - s.dynamique_variables.s)*s.dynamique_variables.u*s.mult*s.neuron_pre.dynamique_variables.r
    elseif s.depression
        @fastmath s.dynamique_variables.d += s.f_param.dt*((1.0 - s.dynamique_variables.d)/s.f_param.τ - s.dynamique_variables.d*(1.0 - s.d_param.fD)*s.neuron_pre.dynamique_variables.r)
      #  push!(s.d_save,s.dynamique_variables.d)

        return  @fastmath -s.dynamique_variables.s/s.τ + s.γ*s.mult*s.dynamique_variables.d*s.neuron_pre.dynamique_variables.r
    else
        return  @fastmath -s.dynamique_variables.s/s.τ + s.γ*(1.0 - s.dynamique_variables.s)*s.neuron_pre.dynamique_variables.r
    end

end

function current_synapses!(ln::Vector{N} where N <: neuron,d::Dict{String, Vector{Float64}},index::Int64)
    #compute the sum of syn currents
    # separe in two the currents (due to the dendrites)
for n in ln

    n.dynamique_variables.Iexc =0.0
    n.dynamique_variables.Iinh =0.0

    @simd for s in n.list_syn
        if s.g <0.0
             n.dynamique_variables.Iinh += s.g * s.dynamique_variables.s
        else
             n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
        end
    end

   # push!(n.Iexc_save , n.dynamique_variables.Iexc)
   # push!(n.Iinh_save , n.dynamique_variables.Iinh)
   @inbounds  n.dynamique_variables.Istim = d[n.name][index]
end
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


function current_synapses!(n::dend_sigmoid, list_syn::Vector{T} where T<: synapses)
    #compute the sum of syn currents
    # separe in two the currents (due to the dendrites)

    n.dynamique_variables.Iexc =0.0
    n.dynamique_variables.Iinh =0.0

     for s in list_syn
        # if s.g < 0.0
        #    n.dynamique_variables.Iinh += s.g * s.dynamique_variables.s
        # else
        #     n.dynamique_variables.Iexc += s.g * s.dynamique_variables.s
        # end
        add_current!(s.g, s.dynamique_variables.s, n.dynamique_variables.Iexc, n.dynamique_variables.Iinh)
    end

    # push!(n.Iexc_save , n.dynamique_variables.Iexc)
    # push!(n.Iinh_save , n.dynamique_variables.Iinh)
    return
end



function update_dend!(d::dend_Sean2020)
    current_synapses!(d)
    update_process!(d.OU_process)
    @. d.Iexc += d.Inoise[end] + d.Istim + d.Ibg
    dendrite_input_output!(d)
end

function update_dend!(d::dend_sigmoid)
    current_synapses!(d, d.list_syn)
   # update_process!(d.OU_process)
    d.dynamique_variables.Iexc += d.dynamique_variables.Inoise + d.dynamique_variables.Istim + d.dynamique_variables.Ibg
    d.dynamique_variables.Itot = d.dynamique_variables.Iexc + d.dynamique_variables.Iinh
 #   push!(d.Itot_save, d.dynamique_variables.Itot)
    dendrite_input_output!(d)
end



function update_dend!(d::dend_sigmoid, index::Int64)
    current_synapses!(d, d.list_syn)
    #update_process!(d.OU_process)
    @inbounds d.dynamique_variables.Iexc += d.OU_process.noise[index] + d.dynamique_variables.Istim + d.dynamique_variables.Ibg
    d.dynamique_variables.Itot = d.dynamique_variables.Iexc + d.dynamique_variables.Iinh
   # push!(d.Itot_save, d.dynamique_variables.Itot)
    dendrite_input_output!(d)
end


function update_adaptation!(n::T where T <:neuron, dt::Float64)
    #adapt = n.adaptation
    @fastmath @inbounds  n.adaptation.sA[1] += dt*(-n.adaptation.sA[1]/n.adaptation.τA + n.dynamique_variables.r)
  #  push!(adapt.sA_save, adapt.sA)

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
#using ..attractor_network_functions

using ...BasicFunctions


export time_step, full_time_dynamics

function OU_process(n::T where T <: neuron)
    update_process!(n.OU_process)
end

### TODO: change the funciton to update the variable s of the synapses

function time_step(c::microcircuit,sim::simulation_parameters)

    for pop in [c.list_pv, c.list_dend, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop
            @simd for s in n.list_syn
                s.dynamique_variables.s += s.dt * synapse_derivative(s)
           #     push!(s.s_save , s.dynamique_variables.s)

            end
        end

    end


    #update dend
    update_dend!.(c.list_dend)

    #update neurons
    current_synapses!.(c.list_soma)
    current_synapses!.(c.list_vip)
    current_synapses!.(c.list_sst)
    current_synapses!.(c.list_pv)
    current_synapses!.(c.list_integrator)
    OU_process.(c.list_soma)
    OU_process.(c.list_sst)
    OU_process.(c.list_vip)
    OU_process.(c.list_pv)
    OU_process.(c.list_integrator)

    current_to_frequency!.(c.list_soma)
    current_to_frequency!.(c.list_sst)
    current_to_frequency!.(c.list_vip)
    current_to_frequency!.(c.list_pv)
    current_to_frequency!.(c.list_integrator)


    sum_input!.(c.list_soma)
    sum_input!.(c.list_sst)
    sum_input!.(c.list_vip)
    sum_input!.(c.list_pv)
    sum_input!.(c.list_integrator)


    update_firing!.(c.list_soma)
    update_firing!.(c.list_vip)
    update_firing!.(c.list_sst)
    update_firing!.(c.list_pv)
    update_firing!.(c.list_integrator)

    @simd for nn in c.nn
      #  for n in nn.list_units
       #     current_synapses!(n)
        
        #end
      update_s!(nn,sim)
    end
   
    # for nn in c.nn
    #     for n in nn.list_units
    #        @simd for s in n.list_syn

    #             s.s = s.s+ s.dt * synapse_derivative(s)
    #     #    push!(s.s_save, s.s)
    #         end
    #     end
   # update_s!(nn,sim)
   # end
end

function time_step(l_c::Vector{microcircuit},sim::simulation_parameters, d::Dict{String, Vector{Float64}},index::Int64)

    for c in l_c
    for pop in [c.list_pv, c.list_dend, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop
           @simd for s in n.list_syn
                
             s.dynamique_variables.s +=  s.dt * synapse_derivative(s)
         #     push!(s.s_save,s.dynamique_variables.s)
            end
        end

    end


    #update dend
    update_dend!.(c.list_dend,index)

    #update neurons
    current_synapses!(c.list_soma,d,index)
    current_synapses!(c.list_vip,d,index)
    current_synapses!(c.list_sst,d,index)
    current_synapses!(c.list_pv,d,index)
    current_synapses!(c.list_integrator,d,index)
    #OU_process.(c.list_soma)
    #OU_process.(c.list_sst)
    #OU_process.(c.list_vip)
    #OU_process.(c.list_pv)
    #OU_process.(c.list_integrator)

    #current_to_frequency!.(c.list_soma)
    #current_to_frequency!.(c.list_sst)
    #current_to_frequency!.(c.list_vip)
    #current_to_frequency!.(c.list_pv)
    #current_to_frequency!.(c.list_integrator)


    sum_input!.(c.list_soma,index)
    sum_input!.(c.list_sst,index)
    sum_input!.(c.list_vip,index)
    sum_input!.(c.list_pv,index)
    sum_input!.(c.list_integrator,index)


    update_firing!.(c.list_soma)
    update_firing!.(c.list_vip)
    update_firing!.(c.list_sst)
    update_firing!.(c.list_pv)
    update_firing!.(c.list_integrator)


end
end

# function full_time_dynamics(c::microcircuit, sim::simulation_parameters)
#     # TODO: to write better. compute the whole time dynamics of the microcuit
#     list_time = 0.0:sim.dt:sim.Tfin
#     list_stim = zeros(length(list_time))
#     temp_t = 0.0
#     for i=1:length(list_stim)-1
#         temp_t +=sim.dt

#         if temp_t>sim.Tstimduration && i>sim.Tstimdebut/sim.dt
#             if list_stim[i] == 0.0
#                 list_stim[i] = 0.2
#             else
#                 list_stim[i] = 0.0
#             end
#             temp_t = 0.0

#         end
#         list_stim[i+1] = list_stim[i]
#     end


#     for (index,t) in enumerate(list_time[2:end])
       
#         c.list_soma[1].Istim = sim.Stimulus[index] #nA
#         time_step(c,sim)
#     end
#     # TODO: do a nice return function for this dynamics

# end


function full_time_dynamics(l_c::Vector{microcircuit}, sim::simulation_parameters, d::Dict{String, Vector{Float64}})
    # TODO: to write better. compute the whole time dynamics of the microcuit
    list_time = 0.0:sim.dt:(sim.dt*length(d[l_c[1].list_soma[1].name]))
   
 # TODO a dictionnary of stimlus that can be access through the names of the neurons

    for (index,t) in enumerate(list_time[2:end])
        time_step(l_c,sim,d,index)
    end
    # TODO: do a nice return function for this dynamics

end

# TODO

# TODO

    
end
using .time_dynamics


end