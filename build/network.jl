using DrWatson
@quickactivate "1-project-categorical-MMN"

#loading packages
using Parameters
include(srcdir("basic_functions.jl"))
include(srcdir("structures.jl"))
###################

#### This file will define the type and functions used to create the networks


struct neural_population
    N::Int # number of neurons
    type_neuron::String # types of the neurons in this population
    list_neurons::Vector{neuron}
    type_global::String # can be pyr_cells, dendrites or interneurons for now
end

struct neural_network
    # this structure will contain the informations to construct the network
    list_pop::Vector{neural_population} #list of neural population
    m_prob::Array{N,N} # matrix of the probabilities of connectivity
    m_weights::Array{N,N} # matrix of the synaptic weights

end

function construct_neuron!(np::neural_population)
    # this function will construct a neural populaiton depending on the type of neurons
    if np.type_neuron == "rectified linear neurons"
        for i=1:np.N
            push!(np.list_neurons, rectified_linear_neurons())
        end

    end

end


function connecting_dendrites_one_to_one!(pyr_list::Vector{pyr_cell}, dend_list::Vector{dendrites})
    # this function connects the dendrites to only one pyr cell
    for i=1:length(pyr_list)
        push!(pyr_list[i].dendrites_list, dend_list[i])
    end
end

function connecting_two_populations!(post_list::Array{T,1},pre_list::Array{Y,1},prob::Float64,weight::Float64) where {T<:neuron,Y<:neuron}
    for i=1:length(post_list)
        for j=1:length(pre_list)
            r = rand()
            if r<prob
                push!(post_list[i].neurons_list,
                syn_connection(pre_syn = pre_list[j],w = weight))
            end
        end
    end

end



function construct_network_hertag!(nn::neural_network)

    N_pop = length(n.list_pop)
    N_neurons = zeros(N_pop)
    for i=1:N_pop
        N_neurons[i] = n.list_pop[i].N
    end

    list_types = Vector{String}
    for i=1:N_pop
        push!(list_types,nn.list_pop[i].type_global)
    end


    for np in nn.list_pop
        construct_neuron!(np)
    end

    connecting_dendrites_one_to_one!(nn.listpop[findfirst(isequal("pyr_cells"),list_types)].listneurons,
    nn.listpop[findfirst(isequal("dendrites"),list_types)].listneurons)
    
    for i=1:N_pop
        for j=1:N_pop
            connecting_two_populations!(np.list_pop[i].list_neurons,np.list_pop[j].list_neurons,nn.m_prob[i,j],nn.m_weights[i,j])
        end
    end

end



@with_kw struct syn_probabilities
    pEP::Float64
    pDS::Float64
    pPE::Float64
    pPS::Float64
    pPP::Float64
    pPV::Float64
    pSE::Float64
    pSV::Float64
    pVE::Float64
    pVS::Float64

end

@with_kw struct syn_probabilities_hertag
    pEP::Float64
    pDE::Float64
    pDS::Float64
    pPE::Float64
    pPS::Float64
    pPP::Float64
    pPV::Float64
    pSE::Float64
    pSV::Float64
    pVE::Float64
    pVS::Float64

end


@with_kw struct syn_weights_hertag
    wEP::Float64
    wDE::Float64
    wDS::Float64
    wPE::Float64
    wPS::Float64
    wPP::Float64
    wPV::Float64
    wSE::Float64
    wSV::Float64
    wVE::Float64
    wVS::Float64

end


function construct_network(mat::syn_probabilities)
    NE = 70
    NPV = 10
    NSST = 10 
    NVIP = 10

    pv_list = pv_neurons[]
    for i = 1:NPV
        push!(pv_list, pv_neurons())
    end

    dendrites_list = dendrites[]
    for i = 1:NE
        push!(dendrites_list, dendrites())
    end

    # for now, initialisaiton with only one dendrite
    pyr_list = pyr_cells[]
    for i = 1:NE
        push!(pyr_list, pyr_cells(dendrites_list = [dendrites_list[i]]))
    end

    sst_list = som_neurons[]
    for i = 1:NSST
        push!(sst_list, som_neurons())
    end

    vip_list = vip_neurons[]
    for i = 1:NVIP
        push!(vip_list, vip_neurons())
    end


  
    # assign the syn_probabilities
    ## pyr neurons
    assign_lists!(pyr_list,pv_list,mat.pEP)

    ##pv neurons
    assign_lists!(pv_list,pyr_list,mat.pPE)
    assign_lists!(pv_list,pv_list,mat.pPP)
    assign_lists!(pv_list,sst_list,mat.pPS)
    assign_lists!(pv_list,vip_list,mat.pPV)

    ##sst neurons
    assign_lists!(sst_list,pyr_list,mat.pSE)
    assign_lists!(sst_list,vip_list,mat.pSV)

    ## vip neurons
    assign_lists!(vip_list,pyr_list,mat.pVE)
    assign_lists!(vip_list,sst_list,mat.pVS)

    ## dendrites
    assign_lists!(dendrites_list,sst_list,mat.pDS)


    return network(pv_list,pyr_list,sst_list,vip_list,dendrites_list)
end


function construct_network_hertag(mat::syn_probabilities_hertag,s::syn_weights_hertag)
    NE = 70
    NPV = 10
    NSST = 10 
    NVIP = 10

    pv_list = pv_neurons_hertag[]
    for i = 1:NPV
        push!(pv_list, pv_neurons_hertag())
    end

    dendrites_list = dendrites_hertag[]
    for i = 1:NE
        push!(dendrites_list, dendrites_hertag())
    end

    # for now, initialisaiton with only one dendrite
    pyr_list = soma_hertag[]
    for i = 1:NE
        push!(pyr_list, soma_hertag(dendrites_list = [dendrites_list[i]]))
        push!(dendrites_list[i].pyr_target, pyr_list[i])
    end

    sst_list = som_neurons_hertag[]
    for i = 1:NSST
        push!(sst_list, som_neurons_hertag())
    end

    vip_list = vip_neurons_hertag[]
    for i = 1:NVIP
        push!(vip_list, vip_neurons_hertag())
    end


  
    # assign the syn_probabilities
    ## pyr neurons
    assign_lists!(pyr_list,pv_list,mat.pEP,s.wEP)

    ##pv neurons
    assign_lists!(pv_list,pyr_list,mat.pPE,s.wPE)
    assign_lists!(pv_list,pv_list,mat.pPP,s.wPP)
    assign_lists!(pv_list,sst_list,mat.pPS,s.wPS)
    assign_lists!(pv_list,vip_list,mat.pPV,s.wPV)

    ##sst neurons
    assign_lists!(sst_list,pyr_list,mat.pSE,s.wSE)
    assign_lists!(sst_list,vip_list,mat.pSV,s.wSV)

    ## vip neurons
    assign_lists!(vip_list,pyr_list,mat.pVE,s.wVE)
    assign_lists!(vip_list,sst_list,mat.pVS,s.wVS)

    ## dendrites
    assign_lists!(dendrites_list,sst_list,mat.pDS,s.wDS)
    assign_lists!(dendrites_list,pyr_list,mat.pDE,s.wDE)


    return network_hertag(pv_list,pyr_list,sst_list,vip_list,dendrites_list)
end


function assign_lists!(post_list::Array{T,1},pre_list::Array{Y,1},prob::Float64) where {T<:neuron,Y<:neuron}
    for i=1:length(post_list)
        for j=1:length(pre_list)
            r = rand()
            if r<prob
                push!(post_list[i].neurons_list,
                syn_connection(pre_syn = pre_list[j],w = 0.0))
            end
        end
    end

end


function assign_lists!(post_list::Array{T,1},pre_list::Array{Y,1},prob::Float64,weight::Float64) where {T<:neuron,Y<:neuron}
    for i=1:length(post_list)
        for j=1:length(pre_list)
            r = rand()
            if r<prob
                push!(post_list[i].neurons_list,
                syn_connection(pre_syn = pre_list[j],w = weight))
            end
        end
    end

end
