using DrWatson
@quickactivate "1-project-categorical-MMN"

# loading packages
using Parameters
include(srcdir("structures.jl"))


include(srcdir("basic_functions.jl"))


# This file will define the functions that are applied  to the structures

## Dendrites functions

# TODO: what to do with yang dendrites model
# function update_VD!(d::dendrites)
#     d.VD = 30.0 * (1.0 + tanh((d.gE - d.g12)/d.β)) + d.V0 + d.EL
# end

# function update_g12!(d::dendrites)
#     d.g12 = d.bg * (d.gLD + d.gI)
# end

# function update_β!(d::dendrites)
#     d.β = d.k * exp(d.gI/d.γ)
# end

# function update_snmda(d::dendrites,syn_nmda::syn_NMDA_dendrites)
#     d.snmda = 1.0 - 1.0/(1 + 0.001*d.rE * syn_nmda.τNMDAdec*syn_nmda.αNMDA*syn_nmda.τNMDArise)
# end

# function update_gE!(d::dendrites)
#     d.gI = d.snmda * d.gE
# end
# function update_gI!(d::dendrites,sGABA::syn_GABA_dendrites)
#     # compute gI from the GABA properties
#     d.gI = d.rI * sGABA.τGABA * sGABA.gGABAmax 
# end


function compute_calcium_spike!(dend::dendrites_hertag, euler_m::euler_method)
    # function that computes the calcium dendritic event
    if euler_m.record == true
        push!(dend.calcium_spike, dend.gain_calcium * heaviside(dend.Itot[end] - dend.Θc))
    else
        dend.calcium_spike[end] = dend.gain_calcium * heaviside(dend.Itot[end] - dend.Θc)
    end
end


function update_current!(d::dendrites_hertag, euler_m::euler_method)
    # function that updates the total synaptic input generated in the dendrites
    if euler_m.record == true
        push!(d.Itot,  p.λE * p.Isyn[end] + (1.0 - p.λD) * d.Isyn[end])
    else
        p = d.pyr_target[1]
        d.Itot[end] = p.λE * p.Isyn[end] + (1.0 - p.λD) * d.Isyn[end]
    end
end



################

## Functions generique neuron structure

function update_syn_current!(n::neuron)
    # update the synaptic current towards one neuron
    # Rule: The sign of weights is taken into account in w
    # Background current is added at the end
    n.Isyn = 0.0
    for syn in n.neurons_list
        n.Isyn += syn.pre_syn.fr * syn.w
    end
    n.Isyn += n.Ibg # add the background input to the syn current

end





###################


## Firing rates updates

function update_fr!(n::rectified_linear_neurons, euler_m::euler_method)
    # in tis diff eq, the total FR is rectified
    # the record variables indicates if n.fr is a list or just the last firing rate
    if euler_method.record
        push!(n.fr, n.fr[end])
        n.fr[end] += euler_m.dt / n.τ * (-n.fr[end] + n.Isyn[end] )
        n.fr[end] = rect_linear(n.fr[end])
    else
        n.fr[end] += euler_m.dt / n.τ * (-n.fr[end] + n.Isyn[end] )
        n.fr[end] = rect_linear(n.fr[end])
    end

end

 
function update_fr!(pr::pyr_cells)
    # function that updates the firing rate for the pyramidal cells
    # it updates the somatic FR
    # TODO preciser le type de pyramidal cells
    temp = max(0.0, pr.I + pr.param_firingrate.β) / pr.param_firingrate.γ
    pr.r = exp(fr_soma.α * log(temp))
end


function update_fr!(pyr_cell::soma_hertag, euler_m::euler_method)
    # update the fr of the soma for Hertag model
    pyr_cell.fr = pyr_cell.fr + euler_m.dt / pyr_cell.τE * (-pyr_cell.fr + rect_linear(pyr_cell.Itot - pyr_cell.Θr))

end




###################


## Pyramidal cells functions 
function update_I_DS!(p::pyr_cells)
    temp = 0.0
    for i in p.dendrites_list
        temp = temp + i.VD
    end
    temp = temp / length(p.dendrites_list)
    p.I = 8.0 * (temp + 55.0 )
end



function update_current!(pyr_cell::soma_hertag)
    # function that updates the update the total somatic input toward the soma
    d = pyr_cell.dendrites_list[1] # only 1 dendrite connected
    pyr_cell.Itot = pyr_cell.λD * rect_linear(d.Isyn + d.calcium_spike) + 
    (1.0 - pyr_cell.λE) * pyr_cell.Isyn

end


####################


## Setting up the network parameters



function normalisation_weights(S_weights::syn_weights_hertag, P_weights::syn_probabilities_hertag)
    # function that normalizes all the weights to be indpt of N_neurons
    # TODO function a retravailler
    S_weights.wEP = S_weights.wEP / (P_weights.pEP * 10)
    S_weights.wDS = S_weights.wDS / (P_weights.pDS * 10)
    S_weights.wDE = S_weights.wDE / (P_weights.pDE * 70)
    S_weights.wPE = S_weights.wPE / (P_weights.pPE * 70)
    S_weights.wPP = S_weights.wPP / (P_weights.pPP * 10)
    S_weights.wPS = S_weights.wPS / (P_weights.pPS * 10)
    S_weights.wPV = S_weights.wPV / (P_weights.pPV * 10)
    S_weights.wSE = S_weights.wSE / (P_weights.pSE * 70)
    S_weights.wSV = S_weights.wSV / (P_weights.pSV * 10)
    S_weights.wVE = S_weights.wVE / (P_weights.pVE * 70)
    S_weights.wVS = S_weights.wVS / (P_weights.pVS * 10)




end

function background_objective!(n::soma_hertag)
    # Set the objective background input for pyramidal cells following Hertag et al. equations
    temp2 = 0.0
    for syn in n.neurons_list
        temp2 += syn.pre_syn.objective_rate * syn.w
    end
    d = n.dendrites_list[1]
    temp = 0.0
    for syn in d.neurons_list
        temp += syn.pre_syn.objective_rate * syn.w
    end
    temp = max(temp, 0.0)
    n.Ibg = (n.objective_rate + n.Θr - (n.λD * temp)) / (1.0 - n.λE) - temp2
end

function background_objective!(n::rectified_linear_neurons)
    # Set the objective backgroun input for rectified linear neurons
    temp2 = 0.0
    for syn in n.neurons_list
        temp2 += syn.pre_syn.objective_rate * syn.w
    end
    n.Ibg = n.objective_rate - temp2
end


###############################



## Constructing the neural network

"""
    construct_neuron!(np::neural_population)

This function constructs a neural population depending on the arguments of np

...
# Arguments
- `np.N::Integer`: the number of neurons to construct
- `np.type_neuron=rectified linear neurons`: the neurons of this population will be rectified linear neurons
...
"""
function construct_neuron!(np::neural_population)
    if np.type_neuron == "rectified linear neurons"
        for i = 1:np.N
            push!(np.list_neurons, rectified_linear_neurons())
        end

    end

end


function connecting_dendrites_one_to_one!(pyr_list::Vector{pyr_cell}, dend_list::Vector{dendrites})
    # this function connects the dendrites to only one pyr cell
    for i = 1:length(pyr_list)
        push!(pyr_list[i].dendrites_list, dend_list[i])
    end
end

function connecting_two_populations!(post_list::Array{T,1}, pre_list::Array{Y,1}, prob::Float64, weight::Float64) where {T <: neuron,Y <: neuron}
    for i = 1:length(post_list)
        for j = 1:length(pre_list)
            r = rand()
            if r < prob
                push!(post_list[i].neurons_list,
                syn_connection(pre_syn=pre_list[j], w=weight))
            end
        end
    end

end



function construct_network_hertag!(nn::neural_network)

    N_pop = length(n.list_pop)
    N_neurons = zeros(N_pop)
    for i = 1:N_pop
        N_neurons[i] = n.list_pop[i].N
    end

    list_types = Vector{String}
    for i = 1:N_pop
        push!(list_types, nn.list_pop[i].type_global)
    end


    for np in nn.list_pop
        construct_neuron!(np)
    end

    connecting_dendrites_one_to_one!(nn.listpop[findfirst(isequal("pyr_cells"), list_types)].listneurons,
    nn.listpop[findfirst(isequal("dendrites"), list_types)].listneurons)
    
    for i = 1:N_pop
        for j = 1:N_pop
            connecting_two_populations!(np.list_pop[i].list_neurons, np.list_pop[j].list_neurons, nn.m_prob[i,j], nn.m_weights[i,j])
        end
    end

end


###################################


## simulation temporelles

function simulation_step!(nn::neural_network, euler_m::euler_method, type_network == "Hertag 2020")

    # cette function va s'appliquer des qu'il y a des dendrites et des pyr cells 
    # le modele ici est specifique de hertag pour calcium spike

    if type_network == "Hertag 2020"
        list_types = Vector{String}
        for i = 1:N_pop
            push!(list_types, nn.list_pop[i].type_global)
        end
    
        dendrites_pop = findall("dendrites", list_types)
        no_dendrites_pop = findall(not("dendrites"), list_types)
        pyr_cells_pop = findall("pyr_cells", list_types)


        for np in nn.list_pop[no_dendrites_pop]
            for n in np
                update_syn_current!(n) # of everyone

            end
        end

        for np in nn.list_pop[dendrites_pop]
            for n in np
                update_syn_current!(n) # of everyone
                update_current!(n) # all dendrites
                compute_calcium_spike!(n,euler_m) # compute all calcium spikes)
            end
        end

        for np in nn.list_pop[pyr_cells_pop]
            for n in np
                update_current!(n) # all pyr cells
            end
        end

        for np in nn.list_pop[no_dendrites_pop]
            for n in np
                update_fr!(n, euler_m) # update all firing rates

            end
        end
    end

end




#####################