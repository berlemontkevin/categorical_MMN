


module DynamicsFunction

using Parameters



using ..BasicFunctions

export simulation_step!,time_step





module GeneralDynamicalFunctions

using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork

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
    if euler_m.record
        push!(n.fr, n.fr[end])
        n.fr[end] += euler_m.dt / n.τ * (-n.fr[end] + n.Isyn[end] )
        n.fr[end] = rect_linear(n.fr[end])
    else
        n.fr[end] += euler_m.dt / n.τ * (-n.fr[end] + n.Isyn[end] )
        n.fr[end] = rect_linear(n.fr[end])
    end

end


"""
    Function: update_syn_current!(n::neuron)

Updates  the synaptic current toward a neuron. IMPORTANT: the backgrounbd current is included in Isyn
let's note that the synaptic weights have a SIGN
"""
function update_syn_current!(n::neuron)

    n.Isyn[end] = 0.0
    for syn in n.neurons_list
        n.Isyn[end] += syn.pre_syn.fr[end] * syn.w[end]
    end
    n.Isyn[end] += n.Ibg[end] 

end

end
using .GeneralDynamicalFunctions

module Hertag2020_Functions
###################################
# Functions Hertag
###################################


using ...NeuronalStructures.Hertag2020_Structures
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork


using ..GeneralDynamicalFunctions


using ...BasicFunctions

export simulation_step!

"""
    Function: compute_calcium_spike!(dend::dendrites_hertag,euler_m::euler_method)

Compute the calcium spike of a dendrite if a calcium dendritic event occurs 9following hertag 2020). Computed using euler_method
"""
function compute_calcium_spike!(dend::dendrites_hertag, euler_m::euler_method)
    if euler_m.record 
        push!(dend.calcium_spike, dend.gain_calcium * heaviside(dend.Itot[end] - dend.Θc))
    else
        dend.calcium_spike[end] = dend.gain_calcium * heaviside(dend.Itot[end] - dend.Θc)
    end
end

"""
    Function: update_current!(d::dendrites_hertag,euler_m::euler_method)

Update the total current arriving at a dendrite following Hertag 2020 model for dendrites
"""
function update_current!(d::dendrites_hertag, euler_m::euler_method)
    if euler_m.record 
        p = d.pyr_target[1]

        push!(d.Itot,  p.λE * p.Isyn[end] + (1.0 - p.λD) * d.Isyn[end])
    else
        p = d.pyr_target[1]
        d.Itot[end] = p.λE * p.Isyn[end] + (1.0 - p.λD) * d.Isyn[end]
    end
end


"""
    Function: update_fr!(pyr_cell::soma_hertag,euler_m::euler_method)

Update the firing rate of a pyramidal cell following Hertag 2020 model
"""
function GeneralDynamicalFunctions.update_fr!(pyr_cell::soma_hertag, euler_m::euler_method)
    if euler_m.record 
        push!(pyr_cell.fr, pyr_cell.fr[end] + euler_m.dt / pyr_cell.τE * (-pyr_cell.fr[end] + rect_linear(pyr_cell.Itot[end] - pyr_cell.Θr)))
    else
        pyr_cell.fr[end] += euler_m.dt / pyr_cell.τE * (-pyr_cell.fr[end] + rect_linear(pyr_cell.Itot[end] - pyr_cell.Θr))

    end
end


"""
    Function: update_current!(pyr_cell::soma_hertag)

Update the total current of a pyramidal cell following Hertag 2020 model
"""
function update_current!(pyr_cell::soma_hertag, euler_m::euler_method)
    if euler_m.record 
        push!(pyr_cell.Itot, 0.0)
    else
        pyr_cell.Itot[end] = 0.0
        
    end
    
    for d in pyr_cell.dendrites_list
        pyr_cell.Itot[end] += pyr_cell.λD .* rect_linear.(d.Isyn[end] .+ d.calcium_spike[end]) 
        
    end
    pyr_cell.Itot[end] += (1.0 .- pyr_cell.λE) .* pyr_cell.Isyn[end]
end


"""
    Function: background_objective!(n::soma_hertag)
    
This function sets up the bacground current needed to achieve the objective of a soma modeled following Hertag 2020
"""
function background_objective!(n::soma_hertag, euler_m::euler_method)
    temp2 = 0.0
    for syn in n.neurons_list
        temp2 += syn.pre_syn.objective_rate * syn.w[end]
    end
    d = n.dendrites_list[1]
    temp = 0.0
    for syn in d.neurons_list
        temp += syn.pre_syn.objective_rate * syn.w[end]
    end
    temp = max(temp, 0.0)
    if euler_m.record 
        push!(n.Ibg,    (n.objective_rate + n.Θr - (n.λD * temp)) / (1.0 - n.λE) - temp2)
    else
        n.Ibg[end] = (n.objective_rate + n.Θr - (n.λD * temp)) / (1.0 - n.λE) - temp2
    end
end





"""
    Function: simulation_step!(nn::neural_motif,euler_m::euler_method,type_network=="Hertag 2020")


...
This function computes one dt step for a neural motif. 
...
"""
function simulation_step!(nn::neural_motif, euler_m::euler_method, type_network="Hertag 2020")



    if type_network == "Hertag 2020"
        N_pop = length(nn.list_pop)
        list_types = String[]
        for i = 1:N_pop
            push!(list_types, nn.list_pop[i].type_global)
        end
    
        dendrites_pop = findall(x -> x == "dendrites", list_types)
        no_dendrites_pop = findall(x -> x != "dendrites", list_types)
        pyr_cells_pop = findall(x -> x == "pyr_cells", list_types)

    
        for np in nn.list_pop[no_dendrites_pop]
            for n in np.list_neurons
                update_syn_current!(n) # of everyone

            end
        end

        for np in nn.list_pop[dendrites_pop]
            for n in np.list_neurons
                update_syn_current!(n) # of everyone
                update_current!(n, euler_m) # all dendrites
                compute_calcium_spike!(n, euler_m) # compute all calcium spikes)
            end
        end

        for np in nn.list_pop[pyr_cells_pop]
            for n in np.list_neurons
                update_current!(n,euler_m) # all pyr cells
            end
        end

        for np in nn.list_pop[no_dendrites_pop]
            for n in np.list_neurons
                update_fr!(n, euler_m) # update all firing rates

            end
        end
    end

end

end

using .Hertag2020_Functions



module local_circuit_functions

using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.NeuralNetwork
using Parameters
using ...NeuronalStructures.AbstractNeuronalTypes

using ...NeuronalStructures.RateDendrites
using ...NeuronalStructures.local_circuit

using ..GeneralDynamicalFunctions


using ...BasicFunctions


export dendrite_input_output!, sum_input, current_to_frequency!, update_firing!, synapse_derivative, current_synapses!, update_dend!

function dendrite_input_output!(d::dend_Sean2020)
    # dendridic function
    @unpack  c1,c2,c3,c4,c5,c6 =  d.param_c 
    β = c5*exp(-d.Iinh[end]/c6)

    push!(d.Ioutput , c1*(tanh((d.Iexc[end] +c3*d.Iinh[end] + c4)/β)) + c2)
end

# need an abstract type for the neurons
function sum_input(n::T where {T<: local_circuit_interneuron} )
    push!(n.Itot , n.Iinh[end] + n.Iexc[end] + n.Inoise[end] + n.Ibg + n.Istim )
end

function sum_input(n::soma_Sean2020)
    push!(n.Itot , n.Iinh[end] + n.Iexc[end] + n.Inoise[end] + n.Ibg + n.Istim + n.den.Ioutput[end])

end


function current_to_frequency!(neuron::pyr_cell)
    # r is the rate of the neuron
    a = neuron.a
    b = neuron.b
    d = neuron.c
    
    input_current = neuron.Itot[end]

    return (a*input_current - b)/(1 - exp(-d*(a*input_current - b)))

end


function current_to_frequency!(neuron::T where {T<: local_circuit_interneuron})
    # r is the rate of the neuron
    c_I = neuron.c_I
    r0 = neuron.r0
    input_current = neuron.Itot[end]

    return rect_linear(c_I*input_current + r0)
end

function update_firing!(n::neuron)
    # update firing rates following Euler law
    push!(n.r , n.r[end] + n.dt/n.τ*(-n.r[end] + current_to_frequency!(n)))

end


function synapse_derivative(s::ampa_syn)

    return -s.s[end]/s.τ + s.γ*s.neuron_pre.r[end]

end


function synapse_derivative(s::gaba_syn)
    #paper Sean: Inhibition onto the dendrites is slower than inhibition elsewhere (ali and Thomson 2008)
    return -s.s[end]/s.τ + s.γ*s.neuron_pre.r[end]

end


function synapse_derivative(s::nmda_syn)

    return -s.s[end]/s.τ + s.γ*(1.0 - s.s[end])*s.neuron_pre.r[end]

end

function current_synapses!(n::neuron)
    #compute the sum of syn currents
    # separe in two the currents (due to the dendrites)

    push!(n.Iexc , 0.0)
    push!(n.Iinh , 0.0)

    for s in n.list_syn
        if s.g <0.0
            n.Iinh[end] += s.g * s.s[end]
        else
            n.Iexc[end] += s.g * s.s[end]
        end
    end
end



function update_dend!(d::dend_Sean2020)
    current_synapses!(d)
    update_process!(d.OU_process)
    d.Iexc[end] += d.Inoise[end] + d.Istim + d.Ibg
    dendrite_input_output!(d)
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


export time_step

function OU_process(n::neuron)
    update_process!(n.OU_process)
end

function time_step(c::microcircuit,sim::simulation_parameters)

    #update dend
    update_dend!.(c.list_dend)

    #update neurons
    current_synapses!.(c.list_soma)
    current_synapses!.(c.list_vip)
    current_synapses!.(c.list_sst)
    current_synapses!.(c.list_pv)
    OU_process.(c.list_soma)
    OU_process.(c.list_sst)
    OU_process.(c.list_vip)
    OU_process.(c.list_pv)
    current_to_frequency!.(c.list_soma)
    current_to_frequency!.(c.list_sst)
    current_to_frequency!.(c.list_vip)
    current_to_frequency!.(c.list_pv)


    sum_input.(c.list_soma)
    sum_input.(c.list_sst)
    sum_input.(c.list_vip)
    sum_input.(c.list_pv)


    update_firing!.(c.list_soma)
    update_firing!.(c.list_vip)
    update_firing!.(c.list_sst)
    update_firing!.(c.list_pv)

    for nn in c.nn
      #  for n in nn.list_units
       #     current_synapses!(n)
        
        #end
    update_s!(nn,sim)
    end
for pop in [c.list_pv, c.list_dend, c.list_sst, c.list_vip, c.list_soma]
    for n in pop
        for s in n.list_syn
        push!(s.s , s.s[end]+ s.dt * synapse_derivative(s))
        end
    end

end

for nn in c.nn
    for n in nn.list_units
        for s in n.list_syn

            push!(s.s , s.s[end]+ s.dt * synapse_derivative(s))
            
        end
    end
   # update_s!(nn,sim)
end
end


    
end
using .time_dynamics


end