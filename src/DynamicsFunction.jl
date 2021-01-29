


module DynamicsFunction
# using DrWatson
# @quickactivate "1-project-categorical-MMN"

# include(srcdir("NeuronalStructures.jl"))
# push!(LOAD_PATH,srcdir())
# using Parameters
using Parameters



using ..BasicFunctions

export simulation_step!


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
    # in tis diff eq, the total FR is rectified
    # the record variables indicates if n.fr is a list or just the last firing rate
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
"""
function update_syn_current!(n::neuron)
    # update the synaptic current towards one neuron
    # Rule: The sign of weights is taken into account in w
    # Background current is added at the end
    n.Isyn[end] = 0.0
    for syn in n.neurons_list
        n.Isyn[end] += syn.pre_syn.fr[end] * syn.w[end]
    end
    n.Isyn[end] += n.Ibg[end] # add the background input to the syn current

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
    # function that updates the total synaptic input generated in the dendrites
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
    Function: simulation_step!(nn::neural_network,euler_m::euler_method,type_network=="Hertag 2020")

    # TODO faire la doc
"""
function simulation_step!(nn::neural_motif, euler_m::euler_method, type_network="Hertag 2020")

    # cette function va s'appliquer des qu'il y a des dendrites et des pyr cells 
    # le modele ici est specifique de hertag pour calcium spike
    # TODO inclure la methode de l'equa diff

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

end