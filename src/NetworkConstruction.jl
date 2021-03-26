
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


module attractor_network_construction

using ...NeuronalStructures.NeuralNetwork
using ..GeneralConstruction
using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.attractor_network
 export construct_attractor_network

function construct_attractor_network()
## TODO : add some parameters as dictionnary


    e1 = wong_wang_cell()
    e2 = wong_wang_cell()
    e1.Istim[1] = 5.2*0.0001*30.0*1.1
    e1.Istim[1] = 5.2*0.0001*30.0*1.1
    e2.Istim[1] = 5.2*0.0001*30.0*0.9
    ww_network = wong_wang_network(threshold = 15.0)
    push!(ww_network.list_units,e1)
    push!(ww_network.list_units,e2)
    return ww_network
end

end

using .attractor_network_construction



module local_microcircuit_network

using ...NeuronalStructures.NeuralNetwork
using ..GeneralConstruction
using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.attractor_network
using ...NeuronalStructures.RateDendrites

using ..attractor_network_construction

export create_network
## TODO: un dictionnairedes param
function construct_local_microcircuit(c::microcircuit)
    # construct the microcircuit in sean case

    vip1 = vip_cell()
    sst1 = sst_cell()
    pv1 = pv_cell()
    dend1 = dend_Sean2020()
    E1 = soma_Sean2020(den=dend1)
    
    vip2 = vip_cell()
    sst2 = sst_cell()
    dend2 = dend_Sean2020()
    E2 = soma_Sean2020(den=dend2)
    
    
    
    push!(dend1.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst1,neuron_post=dend1, g = -0.09))
    push!(E1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E1, g = -0.001))
    push!(E1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=E1, g = 0.18))
    
    push!(vip1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=vip1, g = 0.058))
    push!(vip1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=vip1, g = -0.1))
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst1, g = 0.0435))
    push!(sst1.list_syn, gaba_syn(neuron_pre=vip1,neuron_post=sst1, g = -0.05))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=pv1, g = 0.0435))
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=pv1, g = -0.17))
    push!(pv1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=pv1, g = -0.18))
    
    #connect to E2
    push!(dend2.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst2,neuron_post=dend2, g = -0.09))
    push!(E2.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E2, g = -0.001))
    push!(E2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=E2, g = 0.18))
    
    push!(vip2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=vip2, g = 0.058))
    push!(vip2.list_syn, gaba_syn(neuron_pre=sst2,neuron_post=vip2, g = -0.1))
    
    push!(sst2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=sst2, g = 0.0435))
    push!(sst2.list_syn, gaba_syn(neuron_pre=vip2,neuron_post=sst2, g = -0.05))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E2,neuron_post=pv1, g = 0.0435))
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst2,neuron_post=pv1, g = -0.17))
    
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E2,neuron_post=sst1, g = 0.0435))
    push!(sst2.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst2, g = 0.0435))

    
    push!(c.list_dend,dend1)
    push!(c.list_dend,dend2)

    push!(c.list_soma,E1)
    push!(c.list_soma,E2)

    push!(c.list_vip,vip1)
    push!(c.list_vip,vip2)

    push!(c.list_sst,sst1)
    push!(c.list_sst,sst2)

    push!(c.list_pv,pv1)



end

function connect_areas(ww_network::wong_wang_network,c::microcircuit)

    e1 = ww_network.list_units[1]
    e2 = ww_network.list_units[2]

    E1 = c.list_soma[1]
    E2 = c.list_soma[2]
    push!(c.nn, ww_network )

    push!(e1.list_syn,nmda_syn(neuron_pre=E1,neuron_post=e1, g = 0.3255/0.6))
    push!(e2.list_syn,nmda_syn(neuron_pre=E1,neuron_post=e2, g = 0.3255))

    push!(e2.list_syn,nmda_syn(neuron_pre=E2,neuron_post=e2, g = 0.3255/0.6))
    push!(e1.list_syn,nmda_syn(neuron_pre=E2,neuron_post=e2, g = 0.3255))

end



function create_network()

    c = microcircuit()
    construct_local_microcircuit(c)
    ww = construct_attractor_network()
    connect_areas(ww,c)
    return c
end

end
using .local_microcircuit_network




end