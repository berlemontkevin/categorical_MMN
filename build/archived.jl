

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




module Hertag2020_Construction

using ...NeuronalStructures.NeuralNetwork
using ..GeneralConstruction

export construct_network_hertag!

"""
    Function: construct_network_hertag!(nn::neural_network)

Takes a neural network as argument. Modify it with everything constructed.
Following Hertag model, One to one dendrites connectivity
"""
function construct_network_hertag!(nn::neural_motif, normalized=true)

    N_pop = length(nn.list_pop)
    N_neurons = zeros(N_pop)
    for i = 1:N_pop
        N_neurons[i] = nn.list_pop[i].N
    end

    list_types = String[]
    for i = 1:N_pop
        push!(list_types, nn.list_pop[i].type_global)
    end


    for np in nn.list_pop
        construct_neuron!(np)
    end

    connecting_dendrites_one_to_one!(nn.list_pop[findfirst(isequal("pyr_cells"), list_types)].list_neurons,
    nn.list_pop[findfirst(isequal("dendrites"), list_types)].list_neurons)
    

    if normalized
        for i = 1:N_pop
            for j = 1:N_pop
                if nn.m_prob[i,j] != 0.0
                    temp = nn.m_weights[i,j] ./ (nn.m_prob[i,j] .* nn.list_pop[j].N)
                    connecting_two_populations!(nn.list_pop[i].list_neurons, nn.list_pop[j].list_neurons, nn.m_prob[i,j], temp)
                end
            end
        end

    end

end

end
using .Hertag2020_Construction
