
module local_circuit
# let's start with only one dendrites
using Parameters, StaticArrays
using ..AbstractNeuronalTypes, ...BasicFunctions, ..NeuralNetwork
export dendrites_param, soma_Sean2020, dend_Sean2020, dendrites_param_sigmoid, dend_sigmoid, adaptation_variables
export dynamique_soma_PC, dynamique_dend_sigmoid
export soma_PC, neural_integrator
export pv_cell, sst_cell, vip_cell, gaba_syn, ampa_syn,nmda_syn, local_circuit_interneuron
export local_circuit_interneuron_with_adaptation, local_circuit_interneuron_without_adaptation, dynamique_variables_interneurons, dynamique_variables_synapses




export save_dynamics, get_dynamics
export sst_cell_with_adaptation    
export bump_structure, layer_bump



using DataFrames, DrWatson# ,CSV#,CSVFiles




# inputs to the rate model will be the time average total conductance of exc and inh synapses









@with_kw struct dendrites_param
    # param of the dendrites from :
    c1::Float64 = 0.12 # * 0.001
    c2::Float64 = 0.13624 # * 0.001
    c3::Float64 = 7.0
    c4::Float64 = 0.0 # * 0.001
    c5::Float64 = 0.00964 # * 0.001
    c6::Float64 = 0.02 # * 0.001 #nA units

end







@with_kw mutable struct soma_Sean2020 <: pyr_cell

    param_c = dendrites_param()
    Idendexc::Float64 = 0.0
    Idendinh::Float64 = 0.0
    f_I_curve = zeros(3)
    r::Float64 = 0.0
    Iinput::Float64 = 0.0
    a::Float64 = 0.135 * 1000
    b::Float64 = 54
    c::Float64 = 0.308 # secondes
    den::dendrite # dendrite connected to the soma
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    list_syn::Vector{synapses} = Vector{synapses}()
    Ibg::Float64 = 330.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Vector{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
    Iexc_save::Vector{Float64} = [0.0]
    Iinh_save::Vector{Float64} = [0.0]
    Ioutput_save::Vector{Float64} = [0.0]
    Itot_save::Vector{Float64} = [0.0]
    Inoise_save::Vector{Float64} = [0.0]
    r_save::Vector{Float64} = [0.0]
end


@with_kw mutable struct dend_Sean2020 <: dendrite
    param_c = dendrites_param()
    Iexc::Vector{Float64} = [0.0]
    Iinh::Vector{Float64} = [0.0]
    Ioutput::Vector{Float64} = [0.0]
    list_syn::Vector{synapses} = Vector{synapses}()
    Ibg::Float64 = 30.0 * 0.001
    Itot::Vector{Float64} = [0.0]
    Inoise::Vector{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64 = 0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()

end


## for now interneurons


@with_kw mutable struct dynamique_sst_cell <: dynamique_variables_interneurons
    r::Float64 = 0.0		
    Iinput::Float64 = 0.0
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 300.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0

end



@with_kw mutable struct dynamique_vip_cell <: dynamique_variables_interneurons
    r::Float64 = 0.0		
    Iinput::Float64 = 0.0
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 300.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0

end



@with_kw mutable struct dynamique_integrator
    r::Float64 = 0.0		
    Iexc::Float64 = 0.0
    Iinh::Float64 = 0.0
    Ibg::Float64 = 310.0 * 0.001
    Itot::Float64 = 0.0
    Inoise::Float64 = 0.0
    Istim::Float64 = 0.0
end







function save_dynamics(c::microcircuit, notebook_file::String, save_parameters)
    # funciton that save the dynamics of the network into a csv

    df = DataFrame()

    
    for pop in [c.list_pv, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop

                df[!,string("r-", n.name)] = round.(n.r, digits=6)
                df[!,string("Itot-", n.name)] = round.(n.Itot, digits=6)

            for s in n.list_syn
                df[!,string("s-", s.name)] = round.(s.s, digits=6)

            end
        end

    end
    for pop in [c.list_dend]
        for n in pop
            df[!,string("Ioutput-", n.name)] = round.(n.Ioutput, digits=6)

           
        end
    end


    wsave(datadir("simulations", notebook_file, savename(save_parameters, "csv")), df)



end




function get_dynamics(notebook_file::String, save_parameters)
    # funciton that retrieves the dynamics from the csv

 #   df = wload(datadir("simulations",notebook_file, savename(save_parameters, "csv")))
    df = CSV.read(datadir("simulations", notebook_file, savename(save_parameters, "csv")), DataFrame)
    dtemp = wload(datadir("simulations", notebook_file, savename(save_parameters, "jld2")))
    c = dtemp["circuit"]

    for pop in [c.list_pv, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
        for n in pop

            temp_name = string("r-", n.name)
            n.r = df[:,temp_name]
            temp_name = string("Itot-", n.name)
            n.Itot = df[:,temp_name]

            for s in n.list_syn
                temp_name = string("s-", s.name)
                s.s = df[:,temp_name]
            end
        end

    end

    for pop in [c.list_dend]
        for n in pop

            temp_name = string("Ioutput-", n.name)
            n.Ioutput = df[:,temp_name]
           
        end
    end
    return c

end


end

using .local_circuit




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
