
module NetworkConstruction


module GeneralConstruction

using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuralNetwork

using ...NeuronalStructures.Hertag2020_Structures

using ...NeuronalStructures.NeuronalModels

using ...NeuronalStructures.Synapses

export connecting_dendrites_one_to_one!, connecting_two_populations!
export construct_neuron!

"""
    Function: connecting_dendrites_one_to_one!(pyr_list::Vector{pyr_cell}, dend_list::Vector{dendrites})

This function connects each pyramidal cell to an unique corresponding dendrite
"""
function connecting_dendrites_one_to_one!(pyr_list::Vector{T}, dend_list::Vector{T}) where {T <: neuron,Y <: neuron}
    # this function connects the dendrites to only one pyr cell
    for i = 1:length(pyr_list)
        push!(pyr_list[i].dendrites_list, dend_list[i])
        push!(dend_list[i].pyr_target, pyr_list[i])

    end
end


"""
    Function: connecting_two_populations!(post_list::Array{T,1},pre_list::Array{Y,1},prob::Float64,weight::Float64) where {T<:neuron,Y<:neuron}

This function connects two neural populations together

...
# Arguments
- post_list = list ofp ostsynaptic neurons
- pre_list = list of presynaptic neurons
- prob: prob of connexions from pre-neurons to post-neurons
- weight = synaptic weight of this connection
...
"""
function connecting_two_populations!(post_list::Array{T,1}, pre_list::Array{Y,1}, prob::Float64, weight::Float64) where {T <: neuron,Y <: neuron}
    for i = 1:length(post_list)
        for j = 1:length(pre_list)
            r = rand()
            if r < prob
                push!(post_list[i].neurons_list,
                syn_connection(pre_list[j], [weight]))
            end
        end
    end

end



"""
    Function: construct_neuron!(np::neural_population)

Constructs a list of neuron using the parameters of a neural populaiton
"""
function construct_neuron!(np::neural_population)
    if np.type_neuron == "rectified linear neurons"
        for i = 1:np.N
            push!(np.list_neurons, rectified_linear_neurons())
        end

    elseif np.type_neuron == "soma hertag"
        for i = 1:np.N
            push!(np.list_neurons, soma_hertag())
        end

    elseif np.type_neuron == "dendrite hertag"
        for i = 1:np.N
            push!(np.list_neurons, dendrites_hertag())
        end


    end




end
end

using .GeneralConstruction




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


end