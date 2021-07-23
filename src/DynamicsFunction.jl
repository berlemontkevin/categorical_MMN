


module DynamicsFunction

using Parameters



using ..BasicFunctions

export simulation_step!,time_step

export dendrite_input_output!, update_firing!


export sum_input!, current_to_frequency!, update_dend!, full_time_dynamics


module GeneralDynamicalFunctions

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
    # if euler_m.record
        n.fr += euler_m.dt / n.τ * (-n.fr + n.Isyn )
        n.fr = rect_linear(n.fr)
    # else
    #     n.fr[end] += euler_m.dt / n.τ * (-n.fr + n.Isyn[end] )
    #     n.fr[end] = rect_linear(n.fr[end])
    # end
    push!(n.fr_save, n.fr)

end


"""
    Function: update_syn_current!(n::neuron)

Updates  the synaptic current toward a neuron. IMPORTANT: the backgrounbd current is included in Isyn
let's note that the synaptic weights have a SIGN
"""
function update_syn_current!(n::neuron)

    n.Isyn = 0.0
    for syn in n.neurons_list
        n.Isyn += syn.pre_syn.fr * syn.w
    end
    n.Isyn += n.Ibg 
    push!(n.Isyn_save,n.Isyn)
end



end

using .GeneralDynamicalFunctions




module local_circuit_functions

using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using Parameters
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
    β = c5*exp(-d.Iinh[end]/c6)

    push!(d.Ioutput , c1*(tanh((d.Iexc[end] +c3*d.Iinh[end] + c4)/β)) + c2)
end

#TODO : struct d'une dendrite prenant en argument le type de dendrite

function dendrite_input_output!(d::dend_sigmoid)
    @unpack c1, c2, c3, c4, c5, c6 = d.param_c
    y = (d.Iexc .- c2.*d.Iinh .+ c6 )/(c3.*d.Iinh .+ c4)
    
    d.Ioutput = c1.*(-0.5 .+ sigmoid(y)).+c5
    
    push!(d.Ioutput_save, d.Ioutput )



end


# need an abstract type for the neurons
function sum_input!(n::T where {T<: local_circuit_interneuron} )
    push!(n.Inoise, n.OU_process.noise[end]) #TODO pre-create noise terms

    if n.adaptation_boolean
        update_adaptation!(n)

        Iadapt = n.adaptation.gA * n.adaptation.sA
        n.Itot = n.Iinh + n.Iexc + n.Inoise[end] + n.Ibg + n.Istim + Iadapt
        push!(n.Itot_save , n.Itot)

    else
        n.Itot = n.Iinh + n.Iexc + n.Inoise[end] + n.Ibg + n.Istim 
        push!(n.Itot_save , n.Itot)

    end

end

function sum_input!(n::soma_Sean2020)
    push!(n.Inoise, n.OU_process.noise[end])
    n.Itot = n.Iinh + n.Iexc + n.Inoise[end] + n.Ibg + n.Istim + n.den.Ioutput
    push!(n.Itot_save , n.Itot)

end
function sum_input!(n::soma_PC)
    push!(n.Inoise, n.OU_process.noise[end])

    if n.adaptation_boolean
        update_adaptation!(n)

        Iadapt = n.adaptation.gA * n.adaptation.sA
        n.Itot = n.Iinh + n.Iexc + n.Inoise[end] + n.Ibg + n.Istim + n.den.Ioutput+ Iadapt
        push!(n.Itot_save , n.Itot)

    else
        n.Itot = n.Iinh + n.Iexc + n.Inoise[end] + n.Ibg + n.Istim + n.den.Ioutput
        push!(n.Itot_save , n.Itot)
    end
end

function sum_input!(n::neural_integrator)
    push!(n.Inoise, n.OU_process.noise[end])
    n.Itot = n.Iinh + n.Iexc + n.Inoise[end] + n.Ibg + n.Istim 
    push!(n.Itot_save , n.Itot )
end

function current_to_frequency!(neuron::pyr_cell)
    # r is the rate of the neuron
    a = neuron.a
    b = neuron.b
    d = neuron.c
    
    input_current = neuron.Itot

    return (a*input_current - b)/(1 - exp(-d*(a*input_current - b)))

end


function current_to_frequency!(neuron::neural_integrator)
    a = neuron.a
    b = neuron.b
    d = neuron.c
    # for now, the integrtor is on the f-I curve
    
    return neuron.α.*neuron.r .+ (a*neuron.Itot - b)/(1 - exp(-d*(a*neuron.Itot - b)))
    
end


function current_to_frequency!(neuron::T where {T<: local_circuit_interneuron})
    # r is the rate of the neuron
    c_I = neuron.c_I
    r0 = neuron.r0
    input_current = neuron.Itot

    return rect_linear(c_I*input_current + r0)
end



function update_firing!(n::neuron)
    # update firing rates following Euler law
    n.r = n.r + n.dt/n.τ*(-n.r + current_to_frequency!(n))
    push!(n.r_save,n.r)
end


function synapse_derivative(s::ampa_syn)
    if s.facilitation
        s.u += s.f_param.dt*((s.f_param.U - s.u)/s.f_param.τ + s.f_param.U*(1.0 - s.u)*s.neuron_pre.r)
        push!(s.u_save,s.u)
        return -s.s/s.τ + s.γ*s.mult*s.u*s.neuron_pre.r
    elseif s.depression
        s.d += s.f_param.dt*((1.0 - s.d)/s.f_param.τ - s.d*(1.0 - s.d_param.fD)*s.neuron_pre.r)
        push!(s.d_save,s.d)
        return -s.s/s.τ + s.γ*s.mult*s.d*s.neuron_pre.r
    else
        return -s.s/s.τ + s.γ*s.neuron_pre.r
    end

end


function synapse_derivative(s::gaba_syn)
    #paper Sean: Inhibition onto the dendrites is slower than inhibition elsewhere (ali and Thomson 2008)
    return -s.s/s.τ + s.γ*s.neuron_pre.r

end


function synapse_derivative(s::nmda_syn)
    if s.facilitation
        s.u += s.f_param.dt*((s.f_param.U - s.u)/s.f_param.τ + s.f_param.U*(1.0 - s.u)*s.neuron_pre.r)
        push!(s.u_save,s.u)

        return -s.s/s.τ + s.γ*(1.0 - s.s)*s.u*s.mult*s.neuron_pre.r
    elseif s.depression
        s.d += s.f_param.dt*((1.0 - s.d)/s.f_param.τ - s.d*(1.0 - s.d_param.fD)*s.neuron_pre.r)
        push!(s.d_save,s.d)

        return -s.s/s.τ + s.γ*s.mult*s.d*s.neuron_pre.r
    else
        return -s.s/s.τ + s.γ*(1.0 - s.s)*s.neuron_pre.r
    end

end

function current_synapses!(ln::Vector{N} where N <: neuron,d::Dict{String, Vector{Float64}},index::Int64)
    #compute the sum of syn currents
    # separe in two the currents (due to the dendrites)
for n in ln

    n.Iexc = 0.0
    n.Iinh = 0.0

    for s in n.list_syn
        if s.g <0.0
            n.Iinh += s.g * s.s
        else
            n.Iexc += s.g * s.s
        end
    end

    push!(n.Iexc_save , n.Iexc)
    push!(n.Iinh_save , n.Iinh)
    n.Istim = d[n.name][index]
end


end
function current_synapses!(n::dend_sigmoid)
    #compute the sum of syn currents
    # separe in two the currents (due to the dendrites)

    n.Iexc = 0.0
    n.Iinh = 0.0

    for s in n.list_syn
        if s.g <0.0
            n.Iinh += s.g * s.s
        else
            n.Iexc += s.g * s.s
        end
    end

    push!(n.Iexc_save , n.Iexc)
    push!(n.Iinh_save , n.Iinh)

end



function update_dend!(d::dend_Sean2020)
    current_synapses!(d)
    update_process!(d.OU_process)
    d.Iexc += d.Inoise[end] + d.Istim + d.Ibg
    dendrite_input_output!(d)
end

function update_dend!(d::dend_sigmoid)
    current_synapses!(d)
    update_process!(d.OU_process)
    d.Iexc += d.Inoise[end] + d.Istim + d.Ibg
    d.Itot=d.Iexc + d.Iinh
    push!(d.Itot_save, d.Itot)
    dendrite_input_output!(d)
end


function update_adaptation!(n::neuron)
    adapt = n.adaptation
    adapt.sA += n.dt*(-adapt.sA/adapt.τA + n.r)
    push!(adapt.sA_save, adapt.sA)

end


end
using .local_circuit_functions

module attractor_network_functions

using ...NeuronalStructures.attractor_network
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork


using ..GeneralDynamicalFunctions


using ...BasicFunctions

export update_s!
# TODO: rewrite the main funcitons to overwrite them
function current_to_frequency!(n::wong_wang_cell)
    # taken from Abott and Chance
    temp = n.a * n.Itot[end] - n.b
    push!(n.r , temp / (1.0 - exp(-n.d * temp)))
end
    
    
function update_s!(network::wong_wang_network,simu::simulation_parameters)
    
    
      
        unit1 = network.list_units[1]
        unit2 = network.list_units[2]
    
        # update syn current
    for n in network.list_units
        push!(n.Iexc , 0.0)
        push!(n.Iinh , 0.0)
    
        for s in n.list_syn
            if s.g <0.0
                n.Iinh[end] += s.g * s.neuron_pre.r[end]
            else
                n.Iexc[end] += s.g * s.neuron_pre.r[end]
            end
        end
    end
    # for now no modif of Istim TODO (ameliorer)
       # push!(unit1.Istim,unit1.Istim[end])
        #push!(unit2.Istim,unit2.Istim[end])


        push!(unit1.Ibg, unit1.Iexc[end] + unit1.Iinh[end])
        push!(unit2.Ibg, unit2.Iexc[end] + unit2.Iinh[end])
    
        push!(unit1.Itot , network.jn11*unit1.s[end] - network.jn21*unit2.s[end] + unit1.Ibg[end] + unit1.Istim[end]  + unit1.Inoise[end])
        push!(unit2.Itot , network.jn22*unit2.s[end] - network.jn12*unit1.s[end] + unit2.Ibg[end] + unit2.Istim[end] + unit2.Inoise[end])
      
        # update noise
        dummy_time = simu.dt/network.tampa
        push!(unit1.Inoise, unit1.Inoise[end] + (dummy_time)*(-unit1.Inoise[end]) + sqrt(dummy_time)*network.noiseamp*randn())
        push!(unit2.Inoise, unit2.Inoise[end] + (dummy_time)*(-unit2.Inoise[end]) + sqrt(dummy_time)*network.noiseamp*randn())
      
        # update S variable and rate
       current_to_frequency!(unit1)
       current_to_frequency!(unit2)
       
        push!(unit1.s, unit1.s[end] + simu.dt*(-unit1.s[end]/network.tnmda + (1.0-unit1.s[end])*unit1.γ*unit1.r[end]))
        push!(unit2.s, unit2.s[end] + simu.dt*(-unit2.s[end]/network.tnmda + (1.0-unit2.s[end])*unit2.γ*unit2.r[end]))
      
      
      
end
    

end

using .attractor_network_functions



module time_dynamics


using ...NeuronalStructures.attractor_network

using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using Parameters
using ...NeuronalStructures.AbstractNeuronalTypes

using ...NeuronalStructures.RateDendrites
using ...NeuronalStructures.local_circuit

using ..local_circuit_functions
using ..attractor_network_functions

using ...BasicFunctions


export time_step, full_time_dynamics

function OU_process(n::neuron)
    update_process!(n.OU_process)
end

### TODO: change the funciton to update the variable s of the synapses

function time_step(c::microcircuit,sim::simulation_parameters)

    for pop in [c.list_pv, c.list_dend, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop
            for s in n.list_syn
                s.s = s.s + s.dt * synapse_derivative(s)
                push!(s.s_save , s.s)

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

    for nn in c.nn
      #  for n in nn.list_units
       #     current_synapses!(n)
        
        #end
      update_s!(nn,sim)
    end
   
    for nn in c.nn
        for n in nn.list_units
            for s in n.list_syn

                s.s = s.s+ s.dt * synapse_derivative(s)
            push!(s.s_save, s.s)
            end
        end
   # update_s!(nn,sim)
    end
end

function time_step(l_c::Vector{microcircuit},sim::simulation_parameters, d::Dict{String, Vector{Float64}},index::Int64)

    for c in l_c
    for pop in [c.list_pv, c.list_dend, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop
            for s in n.list_syn
                
              s.s = s.s + s.dt * synapse_derivative(s)
              push!(s.s_save,s.s)
            end
        end

    end


    #update dend
    update_dend!.(c.list_dend)

    #update neurons
    current_synapses!(c.list_soma,d,index)
    current_synapses!(c.list_vip,d,index)
    current_synapses!(c.list_sst,d,index)
    current_synapses!(c.list_pv,d,index)
    current_synapses!(c.list_integrator,d,index)
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

#     for nn in c.nn
#       #  for n in nn.list_units
#        #     current_synapses!(n)
        
#         #end
#       update_s!(nn,sim)
#     end
   
#     for nn in c.nn
#         for n in nn.list_units
#             for s in n.list_syn

#                 s.s = s.s+ s.dt * synapse_derivative(s)
#             push!(s.s_save,s.s)
#             end
#         end
#    # update_s!(nn,sim)
#     end
end
end

function full_time_dynamics(c::microcircuit, sim::simulation_parameters)
    # TODO: to write better. compute the whole time dynamics of the microcuit
    list_time = 0.0:sim.dt:sim.Tfin
    list_stim = zeros(length(list_time))
    temp_t = 0.0
    for i=1:length(list_stim)-1
        temp_t +=sim.dt

        if temp_t>sim.Tstimduration && i>sim.Tstimdebut/sim.dt
            if list_stim[i] == 0.0
                list_stim[i] = 0.2
            else
                list_stim[i] = 0.0
            end
            temp_t = 0.0

        end
        list_stim[i+1] = list_stim[i]
    end


    for (index,t) in enumerate(list_time[2:end])
       
        c.list_soma[1].Istim = sim.Stimulus[index] #nA
        time_step(c,sim)
    end
    # TODO: do a nice return function for this dynamics

end


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