
module NetworkConstruction
using DrWatson
using JLD2
using Parameters

using ..BasicFunctions

export create_network, construct_two_local_microcircuit_integrator, create_deterministic_oscillations_oddball

export construct_two_local_microcircuit_integrator_full_param
export create_deterministic_oddball
export orientation_kernel, create_layer_bump


module GeneralConstruction

using ...NeuronalStructures.AbstractNeuronalTypes
using ...NeuronalStructures.NeuralNetwork

# using ...NeuronalStructures.Hertag2020_Structures

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
    e1.Istim[1] = 0.3255
    e2.Istim[1] = 0.3255
    ww_network = wong_wang_network(threshold=15.0)
    push!(ww_network.list_units, e1)
    push!(ww_network.list_units, e2)
    return ww_network
end

end

using .attractor_network_construction



module local_microcircuit_network

using Base:list_append!!
using ...BasicFunctions

using ...NeuronalStructures.NeuralNetwork
using ..GeneralConstruction
using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.attractor_network
using ...NeuronalStructures.RateDendrites

using ..attractor_network_construction

using ...NeuronalStructures.parameters_simulations


using DrWatson

export create_network, construct_two_local_microcircuit_integrator, construct_two_local_microcircuit_integrator_full_param, construct_one_local_microcircuit_integrator_full_param!, construct_one_local_microcircuit_integrator!, construct_two_local_microcircuit_integrator, construct_one_local_microcircuit_integrator
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
    
    
    
    push!(dend1.list_syn, gaba_syn(τ=10 * 0.001, neuron_pre=sst1, neuron_post=dend1, g=-0.09))
    push!(E1.list_syn, gaba_syn(neuron_pre=pv1, neuron_post=E1, g=-0.001))
    push!(E1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=E1, g=0.18))
    
    push!(vip1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=vip1, g=0.058))
    push!(vip1.list_syn, gaba_syn(neuron_pre=sst1, neuron_post=vip1, g=-0.1))
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=sst1, g=0.0435))
    push!(sst1.list_syn, gaba_syn(neuron_pre=vip1, neuron_post=sst1, g=-0.05))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=pv1, g=0.0435))
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst1, neuron_post=pv1, g=-0.17))
    push!(pv1.list_syn, gaba_syn(neuron_pre=pv1, neuron_post=pv1, g=-0.18))
    
    # connect to E2
    push!(dend2.list_syn, gaba_syn(τ=10 * 0.001, neuron_pre=sst2, neuron_post=dend2, g=-0.09))
    push!(E2.list_syn, gaba_syn(neuron_pre=pv1, neuron_post=E2, g=-0.001))
    push!(E2.list_syn, nmda_syn(neuron_pre=E2, neuron_post=E2, g=0.18))
    
    push!(vip2.list_syn, nmda_syn(neuron_pre=E2, neuron_post=vip2, g=0.058))
    push!(vip2.list_syn, gaba_syn(neuron_pre=sst2, neuron_post=vip2, g=-0.1))
    
    push!(sst2.list_syn, nmda_syn(neuron_pre=E2, neuron_post=sst2, g=0.0435))
    push!(sst2.list_syn, gaba_syn(neuron_pre=vip2, neuron_post=sst2, g=-0.05))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E2, neuron_post=pv1, g=0.0435))
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst2, neuron_post=pv1, g=-0.17))
    
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E2, neuron_post=sst1, g=0.0435))
    push!(sst2.list_syn, nmda_syn(neuron_pre=E1, neuron_post=sst2, g=0.0435))

    
    push!(c.list_dend, dend1)
    push!(c.list_dend, dend2)

    push!(c.list_soma, E1)
    push!(c.list_soma, E2)

    push!(c.list_vip, vip1)
    push!(c.list_vip, vip2)

    push!(c.list_sst, sst1)
    push!(c.list_sst, sst2)

    push!(c.list_pv, pv1)



end

function connect_areas(ww_network::wong_wang_network, c::microcircuit)

    e1 = ww_network.list_units[1]
    e2 = ww_network.list_units[2]

    E1 = c.list_soma[1]
    E2 = c.list_soma[2]
    push!(c.nn, ww_network)

    push!(e1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=e1, g=ww_network.EEself))
    push!(e2.list_syn, nmda_syn(neuron_pre=E1, neuron_post=e2, g=ww_network.EEcross))

    push!(e2.list_syn, nmda_syn(neuron_pre=E2, neuron_post=e2, g=ww_network.EEself))
    push!(e1.list_syn, nmda_syn(neuron_pre=E2, neuron_post=e2, g=ww_network.EEcross))

    sst1 = c.list_sst[1]
    sst2 = c.list_sst[2]    

    push!(sst1.list_syn, nmda_syn(neuron_pre=e1, neuron_post=sst1, g=0.22))
    push!(sst2.list_syn, nmda_syn(neuron_pre=e2, neuron_post=sst2, g=0.22))


end

function connect_areas(ww_network::wong_wang_network, c::microcircuit, temp_topdowm::Float64)
    # TODO
    e1 = ww_network.list_units[1]
    e2 = ww_network.list_units[2]

    E1 = c.list_soma[1]
    E2 = c.list_soma[2]
    push!(c.nn, ww_network)

    push!(e1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=e1, g=ww_network.EEself))
    push!(e2.list_syn, nmda_syn(neuron_pre=E1, neuron_post=e2, g=ww_network.EEcross))

    push!(e2.list_syn, nmda_syn(neuron_pre=E2, neuron_post=e2, g=ww_network.EEself))
    push!(e1.list_syn, nmda_syn(neuron_pre=E2, neuron_post=e2, g=ww_network.EEcross))

    push!(E1.den.list_syn, nmda_syn(neuron_pre=e1, neuron_post=E1.den, g=0.2))
    push!(E2.den.list_syn, nmda_syn(neuron_pre=e2, neuron_post=E2.den, g=0.2))



    vip1 = c.list_vip[1]
    vip2 = c.list_vip[2]    

    push!(vip1.list_syn, nmda_syn(neuron_pre=e1, neuron_post=vip1, g=temp_topdowm))
    push!(vip2.list_syn, nmda_syn(neuron_pre=e2, neuron_post=vip2, g=temp_topdowm))


end

function create_network()

    c = microcircuit()
    construct_local_microcircuit(c)
    ww = construct_attractor_network()
    connect_areas(ww, c)
    return c
end

function create_network(t::Float64)
    # TODO
    c = microcircuit()
    construct_local_microcircuit(c)
    ww = construct_attractor_network()
    connect_areas(ww, c, t)
    return c
end

function construct_local_microcircuit_integrator(c::microcircuit; dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation=false, adaptation=false, depression=false, td_to_vip=true)
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name



    # Important note for now, there can only be adaptation at the soma

    vip1 = vip_cell(name="vipcell1")
    sst1 = sst_cell(name="sstcell1")
    pv1 = pv_cell(name="pvcell1")
    dend1 = dend_sigmoid(param_c=dend_param, name="dend1")
    E1 = soma_PC(den=dend1, adaptation_boolean=adaptation, name="ecell1")
    integrator1 = neural_integrator(name="integrator1")
    
    
    push!(dend1.list_syn, gaba_syn(τ=10 * 0.001, neuron_pre=sst1, neuron_post=dend1, g=-0.09, name="sst1-to-dend1"))# 0.09
    push!(E1.list_syn, gaba_syn(neuron_pre=pv1, neuron_post=E1, g=-0.005, name="pv1-to-ecell1"))# -0.001
    push!(E1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=E1, g=0.18, name="ecell1-to-ecell1"))
    
    push!(vip1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=vip1, g=0.058, name="ecell1-to-vip1"))
    push!(vip1.list_syn, gaba_syn(neuron_pre=sst1, neuron_post=vip1, g=-0.1, name="sst1-to-vip1"))
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=sst1, g=0.0435, name="ecell1-to-sst1"))
    push!(sst1.list_syn, gaba_syn(neuron_pre=vip1, neuron_post=sst1, g=-0.05, name="ss1-to-dend1"))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=pv1, g=0.08435, name="ecell1-to-pv1"))# 0.0435
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst1, neuron_post=pv1, g=-0.17, name="sst1-to-pv1"))# -0.17
    push!(pv1.list_syn, gaba_syn(neuron_pre=pv1, neuron_post=pv1, g=-0.18, name="pv1-to-pv1"))
    
       
    if td_to_vip
        push!(vip1.list_syn, nmda_syn(neuron_pre=integrator1, neuron_post=vip1, g=0.22, depression=depression, name="integrator1-to-vip1"))
    else 
        push!(sst1.list_syn, nmda_syn(neuron_pre=integrator1, neuron_post=sst1, g=0.22, facilitation=facilitation, name="integrator1-to-sst1"))
    end
    push!(integrator1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=integrator1, g=0.15, name="ecell1-to-integrator1"))

    push!(c.list_dend, dend1)
  
    push!(c.list_soma, E1)
   
    push!(c.list_vip, vip1)
  
    push!(c.list_sst, sst1)
  
    push!(c.list_pv, pv1)

    push!(c.list_integrator, integrator1)



end




function create_network(keyword::String; dend_param=dendrites_param_sigmoid(0.20, -7.0, -0.482, 0.00964, 0.11, 0.1), facilitation=false, adaptation=false, save=false, save_parameters=Dict(), notebook_file="", depression=false)
    if keyword == "integrator"
        c = microcircuit()
        construct_local_microcircuit_integrator(c; dend_param, facilitation,adaptation, depression)
        return c
   # elseif keyword == "two-microcircuit"
    #    list_microcircuit = construct_two_local_microcircuit_integrator(; dend_param, facilitation,adaptation, depression)
     #   return list_microcircuit
    end

    if save
        # saving of the circuit without any dynamics
        temp_dict = Dict("circuit" => c)
        wsave(datadir("simulations", notebook_file, DrWatson.savename(save_parameters, "jld2")), temp_dict)
    end

   
end




function construct_one_local_microcircuit_integrator!(c::microcircuit; dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation=false, adaptation=false, depression=false, td_to_vip=true,integrator_tc=0.8, time_tot=1000)
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name
    # has to set up the tests on adding adaptation and depression and facilitation

    vip1 = vip_cell(name=string(c.name, "-", "vipcell1"), Ibg=0.25, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

    sst1 = sst_cell(name=string(c.name, "-", "sstcell1"), Ibg=0.25, adaptation_boolean=true, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

    pv1 = pv_cell(name=string(c.name, "-", "pvcell1"), Ibg=0.29, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))
    dend1 = dend_sigmoid(param_c=dend_param, name=string(c.name, "-", "dend1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))
    E1 = soma_PC(den=dend1, adaptation_boolean=adaptation, name=string(c.name, "-", "ecell1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))
    integrator1 = neural_integrator(τ=integrator_tc, name=string(c.name, "-", "integrator1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))
    
    create_process!(vip1.OU_process)
    create_process!(sst1.OU_process)
    create_process!(pv1.OU_process)
    create_process!(dend1.OU_process)
    create_process!(E1.OU_process)
    create_process!(integrator1.OU_process)

    dend1_gaba = gaba_syn(τ=10 * 0.001, g=-0.09, name=string(c.name, "-", "sst1-to-dend1"))
    push!(dend1.list_syn_post_gaba, dend1_gaba)# 0.09
    push!(sst1.list_syn_pre_gaba, dend1_gaba)# 0.09

    E1_gabe = gaba_syn(g=-0.001, depression=true, name=string(c.name, "-", "pv1-to-ecell1"))
    push!(E1.list_syn_post_gaba, E1_gabe)# -0.001
    push!(pv1.list_syn_pre_gaba, E1_gabe)# -0.001


    E1_nmda = nmda_syn(g=0.18, name=string(c.name, "-", "ecell1-to-ecell1"))
    push!(E1.list_syn_pre_nmda, E1_nmda)
    push!(E1.list_syn_post_nmda, E1_nmda)



    # nmda_syn(neuron_pre=E1, neuron_post=vip1, g=0.058, depression=true, name=string(c.name, "-", "ecell1-to-vip1"))
    vip1_nmda = nmda_syn(g=0.058, depression=true, name=string(c.name, "-", "ecell1-to-vip1"))
    push!(vip1.list_syn_post_nmda, vip1_nmda)
    push!(E1.list_syn_pre_nmda, vip1_nmda)



    vip1_gaba = gaba_syn(g=-0.1, facilitation=true, name=string(c.name, "-", "sst1-to-vip1"))
    push!(vip1.list_syn_post_gaba, vip1_gaba)
    push!(sst1.list_syn_pre_gaba, vip1_gaba)
    



    sst1_nmda = nmda_syn(g=0.0435, facilitation=true, name=string(c.name, "-", "ecell1-to-sst1"))
    push!(sst1.list_syn_post_nmda, sst1_nmda)
    push!(E1.list_syn_pre_nmda, sst1_nmda)


    sstvip_gaba = gaba_syn(g=-0.05, facilitation=true, name=string(c.name, "-", "ss1-to-dend1"))
    push!(sst1.list_syn_post_gaba, sstvip_gaba)
    push!(vip1.list_syn_pre_gaba, sstvip_gaba)



    pvE1_nmda = nmda_syn(g=0.04435, depression=true, name=string(c.name, "-", "ecell1-to-pv1"))
    push!(pv1.list_syn_post_nmda, pvE1_nmda)# 0.0435
    push!(E1.list_syn_pre_nmda, pvE1_nmda)# 0.0435


    pvsst_gaba = gaba_syn(g=-0.17, name=string(c.name, "-", "sst1-to-pv1"))
    push!(pv1.list_syn_post_gaba, pvsst_gaba)# -0.17
    push!(sst1.list_syn_pre_gaba, pvsst_gaba)# -0.17

    pvpv_gaba = gaba_syn(g=-0.18, name=string(c.name, "-", "pv1-to-pv1"))
    push!(pv1.list_syn_post_gaba, pvpv_gaba)
    push!(pv1.list_syn_pre_gaba, pvpv_gaba)

       

    vipint_nmda = nmda_syn(g=0.47, depression=depression, name=string(c.name, "-", "integrator1-to-vip1"))
    push!(vip1.list_syn_post_nmda, vipint_nmda)
    push!(integrator1.list_syn_pre_nmda, vipint_nmda)


    pvint_nmda = nmda_syn(g=0.31, depression=depression, name=string(c.name, "-", "integrator1-to-pv1"))
    push!(pv1.list_syn_post_nmda, pvint_nmda)
    push!(integrator1.list_syn_pre_nmda, pvint_nmda)


    sstint_nmda = nmda_syn(g=0.22, facilitation=facilitation, name=string(c.name, "-", "integrator1-to-sst1"))
    push!(sst1.list_syn_post_nmda, sstint_nmda)
    push!(integrator1.list_syn_pre_nmda, sstint_nmda)


    intE_nmda = nmda_syn(g=0.15, name=string(c.name, "-", "ecell1-to-integrator1"))
    push!(integrator1.list_syn_post_nmda, intE_nmda)
    push!(E1.list_syn_pre_nmda, intE_nmda)


    dendint_nmda = nmda_syn(g=0.4, depression=true, name=string(c.name, "-", "integrator1-to-dend1"))
    push!(dend1.list_syn_post_nmda, dendint_nmda)	
    push!(integrator1.list_syn_pre_nmda, dendint_nmda)	

        
    push!(c.list_dend, dend1)
  
    push!(c.list_soma, E1)
   
    push!(c.list_vip, vip1)
  
    push!(c.list_sst, sst1)
  
    push!(c.list_pv, pv1)

    push!(c.list_integrator, integrator1)



end

function construct_one_local_microcircuit_integrator_full_param!(c::microcircuit; dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), sst_adaptation=true, soma_adaptation=true, pv_to_soma_depression=true, soma_to_vip_facilitation=true, sst_to_vip_facilitation=true, soma_to_sst_facilitation=true, vip_to_sst_facilitation=true, soma_to_pv_depression=true, int_to_vip_depression=true, int_to_pv_depression=true, int_to_sst_facilitation=true, int_to_dend_depression=true,integrator_tc=0.8, time_tot=1000, noise=true, preferred=0.0, int_to_sst_connection=true, top_down_to_interneurons=[0.47, 0.31, 0.22])
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name
    # has to set up the tests on adding adaptation and depression and facilitation
    if noise
        vip1 = vip_cell(name=string(c.name, "-", "vipcell1"), Ibg=0.25, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        sst1 = sst_cell(name=string(c.name, "-", "sstcell1"), Ibg=0.25, adaptation_boolean=sst_adaptation, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        pv1 = pv_cell(name=string(c.name, "-", "pvcell1"), Ibg=0.29, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        dend1 = dend_sigmoid(param_c=dend_param, name=string(c.name, "-", "dend1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        E1 = soma_PC(den=dend1, adaptation_boolean=soma_adaptation, name=string(c.name, "-", "ecell1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)), preferred_stim=preferred)
        integrator1 = neural_integrator(τ=integrator_tc, name=string(c.name, "-", "integrator1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))
    else
        vip1 = vip_cell(name=string(c.name, "-", "vipcell1"), Ibg=0.25, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        sst1 = sst_cell(name=string(c.name, "-", "sstcell1"), Ibg=0.25, adaptation_boolean=sst_adaptation, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        pv1 = pv_cell(name=string(c.name, "-", "pvcell1"), Ibg=0.30, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        dend1 = dend_sigmoid(param_c=dend_param, name=string(c.name, "-", "dend1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        E1 = soma_PC(den=dend1, adaptation_boolean=soma_adaptation, name=string(c.name, "-", "ecell1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0), preferred_stim=preferred)
        integrator1 = neural_integrator(τ=integrator_tc, name=string(c.name, "-", "integrator1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))
    end
    create_process!(vip1.OU_process)
    create_process!(sst1.OU_process)
    create_process!(pv1.OU_process)
    create_process!(dend1.OU_process)
    create_process!(E1.OU_process)
    create_process!(integrator1.OU_process)

    dend1_gaba = gaba_syn(τ=10 * 0.001, g=-0.09, name=string(c.name, "-", "sst1-to-dend1"))
    push!(dend1.list_syn_post_gaba, dend1_gaba)# 0.09
    push!(sst1.list_syn_pre_gaba, dend1_gaba)# 0.09

    E1_gabe = gaba_syn(g=-0.01, depression=pv_to_soma_depression, name=string(c.name, "-", "pv1-to-ecell1"))
    push!(E1.list_syn_post_gaba, E1_gabe)# -0.001
    push!(pv1.list_syn_pre_gaba, E1_gabe)# -0.001


    E1_nmda = nmda_syn(g=0.18, name=string(c.name, "-", "ecell1-to-ecell1"))
    push!(E1.list_syn_pre_nmda, E1_nmda)
    push!(E1.list_syn_post_nmda, E1_nmda)



    # nmda_syn(neuron_pre=E1, neuron_post=vip1, g=0.058, depression=true, name=string(c.name, "-", "ecell1-to-vip1"))
    vip1_nmda = nmda_syn(g=0.058, depression=soma_to_vip_facilitation, name=string(c.name, "-", "ecell1-to-vip1"))
    push!(vip1.list_syn_post_nmda, vip1_nmda)
    push!(E1.list_syn_pre_nmda, vip1_nmda)



    vip1_gaba = gaba_syn(g=-0.1, facilitation=sst_to_vip_facilitation, name=string(c.name, "-", "sst1-to-vip1"))
    push!(vip1.list_syn_post_gaba, vip1_gaba)
    push!(sst1.list_syn_pre_gaba, vip1_gaba)
    



    sst1_nmda = nmda_syn(g=0.0435, facilitation=soma_to_sst_facilitation, name=string(c.name, "-", "ecell1-to-sst1"))
    push!(sst1.list_syn_post_nmda, sst1_nmda)
    push!(E1.list_syn_pre_nmda, sst1_nmda)


    sstvip_gaba = gaba_syn(g=-0.05, facilitation=vip_to_sst_facilitation, name=string(c.name, "-", "ss1-to-dend1"))
    push!(sst1.list_syn_post_gaba, sstvip_gaba)
    push!(vip1.list_syn_pre_gaba, sstvip_gaba)



    pvE1_nmda = nmda_syn(g=0.17, depression=soma_to_pv_depression, name=string(c.name, "-", "ecell1-to-pv1"))
    push!(pv1.list_syn_post_nmda, pvE1_nmda)# 0.0435
    push!(E1.list_syn_pre_nmda, pvE1_nmda)# 0.0435


    pvsst_gaba = gaba_syn(g=-0.17, name=string(c.name, "-", "sst1-to-pv1"))
    push!(pv1.list_syn_post_gaba, pvsst_gaba)# -0.17
    push!(sst1.list_syn_pre_gaba, pvsst_gaba)# -0.17

    pvpv_gaba = gaba_syn(g=-0.18, name=string(c.name, "-", "pv1-to-pv1"))
    push!(pv1.list_syn_post_gaba, pvpv_gaba)
    push!(pv1.list_syn_pre_gaba, pvpv_gaba)

       

   

    if int_to_sst_connection
        sstint_nmda = nmda_syn(g=top_down_to_interneurons[3], facilitation=int_to_sst_facilitation, name=string(c.name, "-", "integrator1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sstint_nmda)
        push!(integrator1.list_syn_pre_nmda, sstint_nmda)

        vipint_nmda = nmda_syn(g=top_down_to_interneurons[1], depression=int_to_vip_depression, name=string(c.name, "-", "integrator1-to-vip1"))
        push!(vip1.list_syn_post_nmda, vipint_nmda)
        push!(integrator1.list_syn_pre_nmda, vipint_nmda)
    
    
        pvint_nmda = nmda_syn(g=top_down_to_interneurons[1], depression=int_to_pv_depression, name=string(c.name, "-", "integrator1-to-pv1"))
        push!(pv1.list_syn_post_nmda, pvint_nmda)
        push!(integrator1.list_syn_pre_nmda, pvint_nmda)

        sst1_nmda = nmda_syn(g=0.0435, facilitation=soma_to_sst_facilitation, name=string(c.name, "-", "ecell1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sst1_nmda)
        push!(E1.list_syn_pre_nmda, sst1_nmda)

    else
        vipint_nmda = nmda_syn(g=0.59, depression=int_to_vip_depression, name=string(c.name, "-", "integrator1-to-vip1"))
        push!(vip1.list_syn_post_nmda, vipint_nmda)
        push!(integrator1.list_syn_pre_nmda, vipint_nmda)
    
        sstint_nmda = nmda_syn(g=0.15, facilitation=int_to_sst_facilitation, name=string(c.name, "-", "integrator1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sstint_nmda)
        push!(integrator1.list_syn_pre_nmda, sstint_nmda)

        pvint_nmda = nmda_syn(g=0.36, depression=int_to_pv_depression, name=string(c.name, "-", "integrator1-to-pv1"))
        push!(pv1.list_syn_post_nmda, pvint_nmda)
        push!(integrator1.list_syn_pre_nmda, pvint_nmda)

        sst1_nmda = nmda_syn(g=0.0535, facilitation=soma_to_sst_facilitation, name=string(c.name, "-", "ecell1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sst1_nmda)
        push!(E1.list_syn_pre_nmda, sst1_nmda)

    end

    intE_nmda = nmda_syn(g=0.15, name=string(c.name, "-", "ecell1-to-integrator1"))
    push!(integrator1.list_syn_post_nmda, intE_nmda)
    push!(E1.list_syn_pre_nmda, intE_nmda)


    dendint_nmda = nmda_syn(g=0.4, depression=int_to_dend_depression, name=string(c.name, "-", "integrator1-to-dend1"))
    push!(dend1.list_syn_post_nmda, dendint_nmda)	
    push!(integrator1.list_syn_pre_nmda, dendint_nmda)	

        
    push!(c.list_dend, dend1)
  
    push!(c.list_soma, E1)
   
    push!(c.list_vip, vip1)
  
    push!(c.list_sst, sst1)
  
    push!(c.list_pv, pv1)

    push!(c.list_integrator, integrator1)



end


function construct_one_local_microcircuit_integrator!(c::microcircuit; param_microcircuit = parameters_microcircuit(),  param_syn_strength_microcircuit = parameters_syn_strength_microcircuit() )
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name
    # has to set up the tests on adding adaptation and depression and facilitation

    @unpack dend_param, sst_adaptation, soma_adaptation, pv_to_soma_depression, soma_to_vip_facilitation, sst_to_vip_facilitation, soma_to_sst_facilitation, vip_to_sst_facilitation, soma_to_pv_depression, int_to_vip_depression, int_to_pv_depression, int_to_sst_facilitation, int_to_dend_depression, integrator_tc, time_tot, noise, int_to_sst_connection, top_down_to_interneurons, preferred = param_microcircuit


    @unpack   gaba_sst_to_dend, gaba_pv_to_soma, nmda_soma_to_soma, nmda_soma_to_vip,
    gaba_sst_to_vip, nmda_soma_to_sst, gaba_vip_to_sst, nmda_soma_to_pv, gaba_sst_to_pv, gaba_pv_to_pv, nmda_soma_to_sst, nmda_soma_to_int, nmda_int_to_dend = param_syn_strength_microcircuit




    if noise
        vip1 = vip_cell(name=string(c.name, "-", "vipcell1"), Ibg=0.25, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        sst1 = sst_cell(name=string(c.name, "-", "sstcell1"), Ibg=0.25, adaptation_boolean=sst_adaptation, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        pv1 = pv_cell(name=string(c.name, "-", "pvcell1"), Ibg=0.29, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        dend1 = dend_sigmoid(param_c=dend_param, name=string(c.name, "-", "dend1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        E1 = soma_PC(den=dend1, adaptation_boolean=soma_adaptation, name=string(c.name, "-", "ecell1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)), preferred_stim=preferred)
        integrator1 = neural_integrator(τ=integrator_tc, name=string(c.name, "-", "integrator1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))
    else
        vip1 = vip_cell(name=string(c.name, "-", "vipcell1"), Ibg=0.25, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        sst1 = sst_cell(name=string(c.name, "-", "sstcell1"), Ibg=0.25, adaptation_boolean=sst_adaptation, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        pv1 = pv_cell(name=string(c.name, "-", "pvcell1"), Ibg=0.29, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        dend1 = dend_sigmoid(param_c=dend_param, name=string(c.name, "-", "dend1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        E1 = soma_PC(den=dend1, adaptation_boolean=soma_adaptation, name=string(c.name, "-", "ecell1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0), preferred_stim=preferred)
        integrator1 = neural_integrator(τ=integrator_tc, name=string(c.name, "-", "integrator1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))
    end
    create_process!(vip1.OU_process)
    create_process!(sst1.OU_process)
    create_process!(pv1.OU_process)
    create_process!(dend1.OU_process)
    create_process!(E1.OU_process)
    create_process!(integrator1.OU_process)

    dend1_gaba = gaba_syn(τ=10 * 0.001, g = gaba_sst_to_dend, name=string(c.name, "-", "sst1-to-dend1"))
    push!(dend1.list_syn_post_gaba, dend1_gaba)# 0.09
    push!(sst1.list_syn_pre_gaba, dend1_gaba)# 0.09

    E1_gabe = gaba_syn(g= gaba_pv_to_soma, depression=pv_to_soma_depression, name=string(c.name, "-", "pv1-to-ecell1"))
    push!(E1.list_syn_post_gaba, E1_gabe)# -0.001
    push!(pv1.list_syn_pre_gaba, E1_gabe)# -0.001


    E1_nmda = nmda_syn(g= nmda_soma_to_soma, name=string(c.name, "-", "ecell1-to-ecell1"))
    push!(E1.list_syn_pre_nmda, E1_nmda)
    push!(E1.list_syn_post_nmda, E1_nmda)



    vip1_nmda = nmda_syn(g= nmda_soma_to_vip, depression=soma_to_vip_facilitation, name=string(c.name, "-", "ecell1-to-vip1"))
    push!(vip1.list_syn_post_nmda, vip1_nmda)
    push!(E1.list_syn_pre_nmda, vip1_nmda)



    vip1_gaba = gaba_syn(g= gaba_sst_to_vip, facilitation=sst_to_vip_facilitation, name=string(c.name, "-", "sst1-to-vip1"))
    push!(vip1.list_syn_post_gaba, vip1_gaba)
    push!(sst1.list_syn_pre_gaba, vip1_gaba)
    



    sst1_nmda = nmda_syn(g= nmda_soma_to_sst, facilitation=soma_to_sst_facilitation, name=string(c.name, "-", "ecell1-to-sst1"))
    push!(sst1.list_syn_post_nmda, sst1_nmda)
    push!(E1.list_syn_pre_nmda, sst1_nmda)


    sstvip_gaba = gaba_syn(g= gaba_vip_to_sst, facilitation=vip_to_sst_facilitation, name=string(c.name, "-", "ss1-to-dend1"))
    push!(sst1.list_syn_post_gaba, sstvip_gaba)
    push!(vip1.list_syn_pre_gaba, sstvip_gaba)



    pvE1_nmda = nmda_syn(g= nmda_soma_to_pv, depression=soma_to_pv_depression, name=string(c.name, "-", "ecell1-to-pv1"))
    push!(pv1.list_syn_post_nmda, pvE1_nmda)
    push!(E1.list_syn_pre_nmda, pvE1_nmda)


    pvsst_gaba = gaba_syn(g= gaba_sst_to_pv, name=string(c.name, "-", "sst1-to-pv1"))
    push!(pv1.list_syn_post_gaba, pvsst_gaba)
    push!(sst1.list_syn_pre_gaba, pvsst_gaba)

    pvpv_gaba = gaba_syn(g= gaba_pv_to_pv, name=string(c.name, "-", "pv1-to-pv1"))
    push!(pv1.list_syn_post_gaba, pvpv_gaba)
    push!(pv1.list_syn_pre_gaba, pvpv_gaba)

       

   

    if int_to_sst_connection
        sstint_nmda = nmda_syn(g=top_down_to_interneurons[3], facilitation=int_to_sst_facilitation, name=string(c.name, "-", "integrator1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sstint_nmda)
        push!(integrator1.list_syn_pre_nmda, sstint_nmda)

        vipint_nmda = nmda_syn(g=top_down_to_interneurons[1], depression=int_to_vip_depression, name=string(c.name, "-", "integrator1-to-vip1"))
        push!(vip1.list_syn_post_nmda, vipint_nmda)
        push!(integrator1.list_syn_pre_nmda, vipint_nmda)
    
    
        pvint_nmda = nmda_syn(g=top_down_to_interneurons[1], depression=int_to_pv_depression, name=string(c.name, "-", "integrator1-to-pv1"))
        push!(pv1.list_syn_post_nmda, pvint_nmda)
        push!(integrator1.list_syn_pre_nmda, pvint_nmda)

        sst1_nmda = nmda_syn(g= nmda_soma_to_sst, facilitation=soma_to_sst_facilitation, name=string(c.name, "-", "ecell1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sst1_nmda)
        push!(E1.list_syn_pre_nmda, sst1_nmda)

    else
        vipint_nmda = nmda_syn(g=0.59, depression=int_to_vip_depression, name=string(c.name, "-", "integrator1-to-vip1"))
        push!(vip1.list_syn_post_nmda, vipint_nmda)
        push!(integrator1.list_syn_pre_nmda, vipint_nmda)
    
        sstint_nmda = nmda_syn(g=0.15, facilitation=int_to_sst_facilitation, name=string(c.name, "-", "integrator1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sstint_nmda)
        push!(integrator1.list_syn_pre_nmda, sstint_nmda)

        pvint_nmda = nmda_syn(g=0.36, depression=int_to_pv_depression, name=string(c.name, "-", "integrator1-to-pv1"))
        push!(pv1.list_syn_post_nmda, pvint_nmda)
        push!(integrator1.list_syn_pre_nmda, pvint_nmda)

        sst1_nmda = nmda_syn(g=0.0535, facilitation=soma_to_sst_facilitation, name=string(c.name, "-", "ecell1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sst1_nmda)
        push!(E1.list_syn_pre_nmda, sst1_nmda)

    end

    intE_nmda = nmda_syn(g= nmda_soma_to_int, name=string(c.name, "-", "ecell1-to-integrator1"))
    push!(integrator1.list_syn_post_nmda, intE_nmda)
    push!(E1.list_syn_pre_nmda, intE_nmda)


    dendint_nmda = nmda_syn(g= nmda_int_to_dend, depression=int_to_dend_depression, name=string(c.name, "-", "integrator1-to-dend1"))
    push!(dend1.list_syn_post_nmda, dendint_nmda)	
    push!(integrator1.list_syn_pre_nmda, dendint_nmda)	

        
    push!(c.list_dend, dend1)
  
    push!(c.list_soma, E1)
   
    push!(c.list_vip, vip1)
  
    push!(c.list_sst, sst1)
  
    push!(c.list_pv, pv1)

    push!(c.list_integrator, integrator1)



end


function construct_one_local_microcircuit_integrator(name::String; param_microcircuit = parameters_microcircuit(),  param_syn_strength_microcircuit = parameters_syn_strength_microcircuit() )
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name
    # has to set up the tests on adding adaptation and depression and facilitation

    @unpack dend_param, sst_adaptation, soma_adaptation, pv_to_soma_depression, soma_to_vip_facilitation, sst_to_vip_facilitation, soma_to_sst_facilitation, vip_to_sst_facilitation, soma_to_pv_depression, int_to_vip_depression, int_to_pv_depression, int_to_sst_facilitation, int_to_dend_depression, integrator_tc, time_tot, noise, int_to_sst_connection, top_down_to_interneurons, preferred = param_microcircuit


    @unpack   gaba_sst_to_dend, gaba_pv_to_soma, nmda_soma_to_soma, nmda_soma_to_vip,
    gaba_sst_to_vip, nmda_soma_to_sst, gaba_vip_to_sst, nmda_soma_to_pv, gaba_sst_to_pv, gaba_pv_to_pv, nmda_soma_to_sst, nmda_soma_to_int, nmda_int_to_dend = param_syn_strength_microcircuit




    if noise
        vip1 = vip_cell(name=string(name, "-", "vipcell1"), Ibg=0.25, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        sst1 = sst_cell(name=string(name, "-", "sstcell1"), Ibg=0.25, adaptation_boolean=sst_adaptation, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        pv1 = pv_cell(name=string(name, "-", "pvcell1"), Ibg=0.29, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        dend1 = dend_sigmoid(param_c=dend_param, name=string(name, "-", "dend1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))

        E1 = soma_PC(den=dend1, adaptation_boolean=soma_adaptation, name=string(name, "-", "ecell1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)), preferred_stim=preferred)
        integrator1 = neural_integrator(τ=integrator_tc, name=string(name, "-", "integrator1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot)))
    else
        vip1 = vip_cell(name=string(name, "-", "vipcell1"), Ibg=0.25, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        sst1 = sst_cell(name=string(name, "-", "sstcell1"), Ibg=0.25, adaptation_boolean=sst_adaptation, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        pv1 = pv_cell(name=string(name, "-", "pvcell1"), Ibg=0.29, OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        dend1 = dend_sigmoid(param_c=dend_param, name=string(name, "-", "dend1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))

        E1 = soma_PC(den=dend1, adaptation_boolean=soma_adaptation, name=string(name, "-", "ecell1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0), preferred_stim=preferred)
        integrator1 = neural_integrator(τ=integrator_tc, name=string(name, "-", "integrator1"), OU_process=BasicFunctions.OU_process(noise=zeros(time_tot), σ=0.0))
    end
    create_process!(vip1.OU_process)
    create_process!(sst1.OU_process)
    create_process!(pv1.OU_process)
    create_process!(dend1.OU_process)
    create_process!(E1.OU_process)
    create_process!(integrator1.OU_process)

    dend1_gaba = gaba_syn(τ=10 * 0.001, g = gaba_sst_to_dend, name=string(name, "-", "sst1-to-dend1"))
    push!(dend1.list_syn_post_gaba, dend1_gaba)# 0.09
    push!(sst1.list_syn_pre_gaba, dend1_gaba)# 0.09

    E1_gabe = gaba_syn(g= gaba_pv_to_soma, depression=pv_to_soma_depression, name=string(name, "-", "pv1-to-ecell1"))
    push!(E1.list_syn_post_gaba, E1_gabe)# -0.001
    push!(pv1.list_syn_pre_gaba, E1_gabe)# -0.001


    E1_nmda = nmda_syn(g= nmda_soma_to_soma, name=string(name, "-", "ecell1-to-ecell1"))
    push!(E1.list_syn_pre_nmda, E1_nmda)
    push!(E1.list_syn_post_nmda, E1_nmda)



    vip1_nmda = nmda_syn(g= nmda_soma_to_vip, depression=soma_to_vip_facilitation, name=string(name, "-", "ecell1-to-vip1"))
    push!(vip1.list_syn_post_nmda, vip1_nmda)
    push!(E1.list_syn_pre_nmda, vip1_nmda)



    vip1_gaba = gaba_syn(g= gaba_sst_to_vip, facilitation=sst_to_vip_facilitation, name=string(name, "-", "sst1-to-vip1"))
    push!(vip1.list_syn_post_gaba, vip1_gaba)
    push!(sst1.list_syn_pre_gaba, vip1_gaba)
    



    sst1_nmda = nmda_syn(g= nmda_soma_to_sst, facilitation=soma_to_sst_facilitation, name=string(name, "-", "ecell1-to-sst1"))
    push!(sst1.list_syn_post_nmda, sst1_nmda)
    push!(E1.list_syn_pre_nmda, sst1_nmda)


    sstvip_gaba = gaba_syn(g= gaba_vip_to_sst, facilitation=vip_to_sst_facilitation, name=string(name, "-", "ss1-to-dend1"))
    push!(sst1.list_syn_post_gaba, sstvip_gaba)
    push!(vip1.list_syn_pre_gaba, sstvip_gaba)



    pvE1_nmda = nmda_syn(g= nmda_soma_to_pv, depression=soma_to_pv_depression, name=string(name, "-", "ecell1-to-pv1"))
    push!(pv1.list_syn_post_nmda, pvE1_nmda)
    push!(E1.list_syn_pre_nmda, pvE1_nmda)


    pvsst_gaba = gaba_syn(g= gaba_sst_to_pv, name=string(name, "-", "sst1-to-pv1"))
    push!(pv1.list_syn_post_gaba, pvsst_gaba)
    push!(sst1.list_syn_pre_gaba, pvsst_gaba)

    pvpv_gaba = gaba_syn(g= gaba_pv_to_pv, name=string(name, "-", "pv1-to-pv1"))
    push!(pv1.list_syn_post_gaba, pvpv_gaba)
    push!(pv1.list_syn_pre_gaba, pvpv_gaba)

       

   

    if int_to_sst_connection
        sstint_nmda = nmda_syn(g=top_down_to_interneurons[3], facilitation=int_to_sst_facilitation, name=string(name, "-", "integrator1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sstint_nmda)
        push!(integrator1.list_syn_pre_nmda, sstint_nmda)

        vipint_nmda = nmda_syn(g=top_down_to_interneurons[1], depression=int_to_vip_depression, name=string(name, "-", "integrator1-to-vip1"))
        push!(vip1.list_syn_post_nmda, vipint_nmda)
        push!(integrator1.list_syn_pre_nmda, vipint_nmda)
    
    
        pvint_nmda = nmda_syn(g=top_down_to_interneurons[1], depression=int_to_pv_depression, name=string(name, "-", "integrator1-to-pv1"))
        push!(pv1.list_syn_post_nmda, pvint_nmda)
        push!(integrator1.list_syn_pre_nmda, pvint_nmda)

        sst1_nmda = nmda_syn(g= nmda_soma_to_sst, facilitation=soma_to_sst_facilitation, name=string(name, "-", "ecell1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sst1_nmda)
        push!(E1.list_syn_pre_nmda, sst1_nmda)

    else
        vipint_nmda = nmda_syn(g=0.59, depression=int_to_vip_depression, name=string(name, "-", "integrator1-to-vip1"))
        push!(vip1.list_syn_post_nmda, vipint_nmda)
        push!(integrator1.list_syn_pre_nmda, vipint_nmda)
    
        sstint_nmda = nmda_syn(g=0.15, facilitation=int_to_sst_facilitation, name=string(name, "-", "integrator1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sstint_nmda)
        push!(integrator1.list_syn_pre_nmda, sstint_nmda)

        pvint_nmda = nmda_syn(g=0.36, depression=int_to_pv_depression, name=string(name, "-", "integrator1-to-pv1"))
        push!(pv1.list_syn_post_nmda, pvint_nmda)
        push!(integrator1.list_syn_pre_nmda, pvint_nmda)

        sst1_nmda = nmda_syn(g=0.0535, facilitation=soma_to_sst_facilitation, name=string(name, "-", "ecell1-to-sst1"))
        push!(sst1.list_syn_post_nmda, sst1_nmda)
        push!(E1.list_syn_pre_nmda, sst1_nmda)

    end

    intE_nmda = nmda_syn(g= nmda_soma_to_int, name=string(name, "-", "ecell1-to-integrator1"))
    push!(integrator1.list_syn_post_nmda, intE_nmda)
    push!(E1.list_syn_pre_nmda, intE_nmda)


    dendint_nmda = nmda_syn(g= nmda_int_to_dend, depression=int_to_dend_depression, name=string(name, "-", "integrator1-to-dend1"))
    push!(dend1.list_syn_post_nmda, dendint_nmda)	
    push!(integrator1.list_syn_pre_nmda, dendint_nmda)	

        
    c = microcircuit{soma_PC, dend_sigmoid}(soma = E1, sst = sst1, vip = vip1, pv = pv1, dend = dend1, integrator = integrator1, name = name )
    return c

end



function construct_two_local_microcircuit_integrator_full_param(;dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), sst_adaptation=true, soma_adaptation=true, pv_to_soma_depression=true, soma_to_vip_facilitation=true, sst_to_vip_facilitation=true, soma_to_sst_facilitation=true, vip_to_sst_facilitation=true, soma_to_pv_depression=true, int_to_vip_depression=true, int_to_pv_depression=true, int_to_sst_facilitation=true, int_to_dend_depression=true,integrator_tc=0.8, time_tot=1000,cross_int_to_vip_depression=true, cross_int_to_pv_depression=true, cross_int_to_sst_facilitation=true, cross_int_to_dend_depression=true ,cross_soma_to_sst_facilitation=true, noise=true, int_to_sst_connection=true, top_down_to_interneurons=[0.47, 0.31, 0.22])


    c1 = microcircuit(name="microcircuit1")
    c2 = microcircuit(name="microcircuit2")
    
    construct_one_local_microcircuit_integrator_full_param!(c1; dend_param, sst_adaptation , soma_adaptation , pv_to_soma_depression, soma_to_vip_facilitation, sst_to_vip_facilitation, soma_to_sst_facilitation , vip_to_sst_facilitation, soma_to_pv_depression , int_to_vip_depression, int_to_pv_depression, int_to_sst_facilitation, int_to_dend_depression,integrator_tc=integrator_tc, time_tot=time_tot, noise=noise, int_to_sst_connection=int_to_sst_connection, top_down_to_interneurons=top_down_to_interneurons)
    construct_one_local_microcircuit_integrator_full_param!(c2; dend_param, sst_adaptation , soma_adaptation , pv_to_soma_depression, soma_to_vip_facilitation, sst_to_vip_facilitation, soma_to_sst_facilitation , vip_to_sst_facilitation, soma_to_pv_depression , int_to_vip_depression, int_to_pv_depression, int_to_sst_facilitation, int_to_dend_depression,integrator_tc=integrator_tc, time_tot=time_tot, noise=noise, int_to_sst_connection=int_to_sst_connection, top_down_to_interneurons=top_down_to_interneurons)


    sst1so2 = nmda_syn(g=0.0435, facilitation=cross_soma_to_sst_facilitation, name=string("ecell2-to-sst1"))
    push!(c1.list_sst[1].list_syn_post_nmda, sst1so2)
    push!(c2.list_soma[1].list_syn_pre_nmda, sst1so2)

    sst2so1 = nmda_syn(g=0.0435, facilitation=cross_soma_to_sst_facilitation, name=string("ecell1-to-sst2"))
    push!(c2.list_sst[1].list_syn_post_nmda, sst2so1)
    push!(c1.list_soma[1].list_syn_pre_nmda, sst2so1)

    dend1int2 = nmda_syn(g=0.1, depression=cross_int_to_dend_depression, name=string("integrator2-to-ecell1"))
    push!(c1.list_dend[1].list_syn_post_nmda, dend1int2)
    push!(c2.list_integrator[1].list_syn_pre_nmda, dend1int2)

    dend2int1 = nmda_syn(g=0.1, depression=cross_int_to_dend_depression, name=string("integrator1-to-ecell2"))
    push!(c2.list_dend[1].list_syn_post_nmda, dend2int1)
    push!(c1.list_integrator[1].list_syn_pre_nmda, dend2int1)

        
        
        
    list_microcircuit = [c1,c2]

    return list_microcircuit
end

function construct_two_local_microcircuit_integrator(;param_microcircuit = parameters_microcircuit(), param_inter_microcircuit = parameters_inter_microcircuit())


    @unpack dend_param, sst_adaptation, soma_adaptation, pv_to_soma_depression, soma_to_vip_facilitation, sst_to_vip_facilitation, soma_to_sst_facilitation, vip_to_sst_facilitation, soma_to_pv_depression, int_to_vip_depression, int_to_pv_depression, int_to_sst_facilitation, int_to_dend_depression, integrator_tc, time_tot, noise, int_to_sst_connection, top_down_to_interneurons = param_microcircuit

    @unpack cross_int_to_vip_depression, cross_int_to_pv_depression, cross_int_to_sst_facilitation, cross_int_to_dend_depression, cross_soma_to_sst_facilitation = param_inter_microcircuit

    # c1 = microcircuit(name="microcircuit1")
    # c2 = microcircuit(name="microcircuit2")
    
    c1 = construct_one_local_microcircuit_integrator("microcircuit1"; param_microcircuit = param_microcircuit, param_inter_microcircuit = param_inter_microcircuit)

    c2 = construct_one_local_microcircuit_integrator("microcircuit2"; param_microcircuit = param_microcircuit, param_inter_microcircuit = param_inter_microcircuit)


    # sst1so2 = nmda_syn(g=0.0435, facilitation=cross_soma_to_sst_facilitation, name=string("ecell2-to-sst1"))
    # push!(c1.list_sst[1].list_syn_post_nmda, sst1so2)
    # push!(c2.list_soma[1].list_syn_pre_nmda, sst1so2)

    # sst2so1 = nmda_syn(g=0.0435, facilitation=cross_soma_to_sst_facilitation, name=string("ecell1-to-sst2"))
    # push!(c2.list_sst[1].list_syn_post_nmda, sst2so1)
    # push!(c1.list_soma[1].list_syn_pre_nmda, sst2so1)

    # dend1int2 = nmda_syn(g=0.1, depression=cross_int_to_dend_depression, name=string("integrator2-to-ecell1"))
    # push!(c1.list_dend[1].list_syn_post_nmda, dend1int2)
    # push!(c2.list_integrator[1].list_syn_pre_nmda, dend1int2)

    # dend2int1 = nmda_syn(g=0.1, depression=cross_int_to_dend_depression, name=string("integrator1-to-ecell2"))
    # push!(c2.list_dend[1].list_syn_post_nmda, dend2int1)
    # push!(c1.list_integrator[1].list_syn_pre_nmda, dend2int1)

    sst1so2 = nmda_syn(g=0.0435, facilitation=cross_soma_to_sst_facilitation, name=string("ecell2-to-sst1"))
    push!(c1.sst.list_syn_post_nmda, sst1so2)
    push!(c2.soma.list_syn_pre_nmda, sst1so2)

    sst2so1 = nmda_syn(g=0.0435, facilitation=cross_soma_to_sst_facilitation, name=string("ecell1-to-sst2"))
    push!(c2.sst.list_syn_post_nmda, sst2so1)
    push!(c1.soma.list_syn_pre_nmda, sst2so1)

    dend1int2 = nmda_syn(g=0.1, depression=cross_int_to_dend_depression, name=string("integrator2-to-ecell1"))
    push!(c1.dend.list_syn_post_nmda, dend1int2)
    push!(c2.integrator.list_syn_pre_nmda, dend1int2)

    dend2int1 = nmda_syn(g=0.1, depression=cross_int_to_dend_depression, name=string("integrator1-to-ecell2"))
    push!(c2.dend.list_syn_post_nmda, dend2int1)
    push!(c1.integrator.list_syn_pre_nmda, dend2int1)
        
        
    list_microcircuit = [c1,c2]

    return list_microcircuit
end

# TODO add a list of microcircuit as argument to the dynamics functions

function construct_two_local_microcircuit_integrator(; dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation=false, adaptation=false, depression=false, td_to_vip=true, τ_integrator=0.8, time_tot=1000)
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name


    c1 = microcircuit(name="microcircuit1")
    c2 = microcircuit(name="microcircuit2")
    
    construct_one_local_microcircuit_integrator!(c1; dend_param, facilitation, adaptation, td_to_vip, depression, integrator_tc=τ_integrator, time_tot=time_tot)
    construct_one_local_microcircuit_integrator!(c2; dend_param, facilitation, adaptation, td_to_vip, depression, integrator_tc=τ_integrator, time_tot=time_tot)

# nmda_syn(neuron_pre=c2.list_soma[1], neuron_post=c1.list_sst[1], g=0.0435, facilitation=true, name=string("ecell2-to-sst1"))
    sst1so2 = nmda_syn(g=0.0435, facilitation=true, name=string("ecell2-to-sst1"))
    push!(c1.list_sst[1].list_syn_post_nmda, sst1so2)
    push!(c2.list_soma[1].list_syn_pre_nmda, sst1so2)

# nmda_syn(neuron_pre=c1.list_soma[1], neuron_post=c2.list_sst[1], g=0.0435, facilitation=true, name=string("ecell1-to-sst2"))
    sst2so1 = nmda_syn(g=0.0435, facilitation=true, name=string("ecell1-to-sst2"))
    push!(c2.list_sst[1].list_syn_post_nmda, sst2so1)
    push!(c1.list_soma[1].list_syn_pre_nmda, sst2so1)

# nmda_syn(neuron_pre=c2.list_integrator[1], neuron_post=c1.list_dend[1], g=0.1, depression=true, name=string("integrator2-to-ecell1"))
    dend1int2 = nmda_syn(g=0.1, depression=true, name=string("integrator2-to-ecell1"))
    push!(c1.list_dend[1].list_syn_post_nmda, dend1int2)
    push!(c2.list_integrator[1].list_syn_pre_nmda, dend1int2)

# nmda_syn(neuron_pre=c1.list_integrator[1], neuron_post=c2.list_dend[1], g=0.1, depression=true, name=string("integrator1-to-ecell2"))
    dend2int1 = nmda_syn(g=0.1, depression=true, name=string("integrator1-to-ecell2"))
    push!(c2.list_dend[1].list_syn_post_nmda, dend2int1)
    push!(c1.list_integrator[1].list_syn_pre_nmda, dend2int1)

        
        
        
    list_microcircuit = [c1,c2]

    return list_microcircuit
end



end
using .local_microcircuit_network



module global_network


using StaticArrays
using Base:list_append!!
using ...BasicFunctions
using ...NeuronalStructures.NeuralNetwork
using ..GeneralConstruction
using ...NeuronalStructures.local_circuit
using ...NeuronalStructures.attractor_network
using ...NeuronalStructures.RateDendrites
using ..attractor_network_construction
using ..local_microcircuit_network
using ...NeuronalStructures.parameters_simulations
using DrWatson

export orientation_kernel, create_layer_bump

function orientation_kernel(x::Float64, y::Float64, param::bump_structure)

    min_dif = abs(x - y)
    if min_dif > 180.0
        min_dif -= 360.0
    end

    return exp(-min_dif.^2 / param.σ.^2)

end

# 43.2 degrees in engel 2015 paper
# 128 units initially, equally spaced
function connect_two_microcircuit(c1::microcircuit, c2::microcircuit, bump_param::bump_structure)
    
#     temp = orientation_kernel(c1.list_soma[1].preferred_stim, c2.list_soma[1].preferred_stim, bump_param)
#     # TODO automatic version with facilitation
#     e1tosst2 = nmda_syn(g=0.0135 * temp, facilitation=true, name=string("ecell2-to-sst1"))
#     push!(c2.list_sst[1].list_syn_post_nmda, e1tosst2)
#     push!(c1.list_soma[1].list_syn_pre_nmda, e1tosst2)

#     e2tosst1 = nmda_syn(g=0.0135 * temp, facilitation=true, name=string("ecell2-to-sst1"))
#     push!(c1.list_sst[1].list_syn_post_nmda, e2tosst1)
#     push!(c2.list_soma[1].list_syn_pre_nmda, e2tosst1)

#     dend1int2 = nmda_syn(g=0.1 * temp, depression=true, name=string("integrator2-to-ecell1"))
#     push!(c1.list_dend[1].list_syn_post_nmda, dend1int2)
#     push!(c2.list_integrator[1].list_syn_pre_nmda, dend1int2)

# # nmda_syn(neuron_pre=c1.list_integrator[1], neuron_post=c2.list_dend[1], g=0.1, depression=true, name=string("integrator1-to-ecell2"))
#     dend2int1 = nmda_syn(g=0.1 * temp, depression=true, name=string("integrator1-to-ecell2"))
#     push!(c2.list_dend[1].list_syn_post_nmda, dend2int1)
#     push!(c1.list_integrator[1].list_syn_pre_nmda, dend2int1)

temp = orientation_kernel(c1.soma.preferred_stim, c2.soma.preferred_stim, bump_param)
# TODO automatic version with facilitation
e1tosst2 = nmda_syn(g=0.0135 * temp, facilitation=true, name=string("ecell2-to-sst1"))
push!(c2.sst.list_syn_post_nmda, e1tosst2)
push!(c1.soma.list_syn_pre_nmda, e1tosst2)

e2tosst1 = nmda_syn(g=0.0135 * temp, facilitation=true, name=string("ecell2-to-sst1"))
push!(c1.sst.list_syn_post_nmda, e2tosst1)
push!(c2.soma.list_syn_pre_nmda, e2tosst1)

dend1int2 = nmda_syn(g=0.1 * temp, depression=true, name=string("integrator2-to-ecell1"))
push!(c1.dend.list_syn_post_nmda, dend1int2)
push!(c2.integrator.list_syn_pre_nmda, dend1int2)

# nmda_syn(neuron_pre=c1.list_integrator[1], neuron_post=c2.list_dend[1], g=0.1, depression=true, name=string("integrator1-to-ecell2"))
dend2int1 = nmda_syn(g=0.1 * temp, depression=true, name=string("integrator1-to-ecell2"))
push!(c2.dend.list_syn_post_nmda, dend2int1)
push!(c1.integrator.list_syn_pre_nmda, dend2int1)

    



end
    


function create_layer_bump(bump_param::bump_structure; dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), sst_adaptation=true, soma_adaptation=true, pv_to_soma_depression=true, soma_to_vip_facilitation=true, sst_to_vip_facilitation=true, soma_to_sst_facilitation=true, vip_to_sst_facilitation=true, soma_to_pv_depression=true, int_to_vip_depression=true, int_to_pv_depression=true, int_to_sst_facilitation=true, int_to_dend_depression=true,integrator_tc=0.8, time_tot=1000, noise=true, preferred=0.0)

    layer1 = layer_bump(bump_param=bump_param)

    for i = 1:bump_param.num_circuits
        c1 = microcircuit{soma_PC}(name="microcircuit$i")
        construct_one_local_microcircuit_integrator_full_param!(c1; dend_param, sst_adaptation , soma_adaptation , pv_to_soma_depression, soma_to_vip_facilitation, sst_to_vip_facilitation, soma_to_sst_facilitation , vip_to_sst_facilitation, soma_to_pv_depression , int_to_vip_depression, int_to_pv_depression, int_to_sst_facilitation, int_to_dend_depression,integrator_tc=integrator_tc, time_tot=time_tot, noise=noise, preferred=i * 360.0 / bump_param.num_circuits)
        push!(layer1.list_microcircuit, c1)
    end


    for i = 1:(bump_param.num_circuits - 1)
        for j = (i + 1):bump_param.num_circuits
            connect_two_microcircuit(layer1.list_microcircuit[i], layer1.list_microcircuit[j], bump_param)
        end
    end

    return layer1
end




function create_layer_bump(bump_param::bump_structure, num_circuits::Int64; param_microcircuit = parameters_microcircuit(), param_inter_microcircuit = parameters_inter_microcircuit(), param_bump_attractor = parameters_bump_attractor())


    # layer1 = layer_bump{soma_PC, dend_sigmoid}(bump_param=bump_param)

    # layer_microcircuit_list = Vector{microcircuit{soma_PC, dend_sigmoid}}()

    # # num_circuits = bump_param.num_circuits

    # for i = 1:num_circuits
    # #     c1 = microcircuit{soma_PC, dend_sigmoid}(name="microcircuit$i")

    #     c1 = construct_one_local_microcircuit_integrator("microcircuit$i";param_microcircuit = parameters_microcircuit(time_tot = param_microcircuit.time_tot, preferred = i * 360.0 / bump_param.num_circuits))
    #     push!(layer_microcircuit_list, c1)
    # end
    temp = @SVector [ construct_one_local_microcircuit_integrator("microcircuit$i";param_microcircuit = parameters_microcircuit(time_tot = param_microcircuit.time_tot, preferred = i * 360.0 / bump_param.num_circuits)) for i=1:128]
    #temp = SVector{128}(layer_microcircuit_list)
    layer1 = layer_bump{soma_PC, dend_sigmoid}(bump_param=bump_param, list_microcircuit = temp)


    for i = 1:(num_circuits - 1)
        for j = (i + 1):num_circuits
            # connect_two_microcircuit(layer_microcircuit_list[i], layer_microcircuit_list[j], bump_param)
            connect_two_microcircuit(layer1.list_microcircuit[i], layer1.list_microcircuit[j], bump_param)
        end
    end

    #  temp = @SVector [layer_microcircuit_list[i] for i=1:128]
    # #temp = SVector{128}(layer_microcircuit_list)
    # layer1 = layer_bump{soma_PC, dend_sigmoid}(bump_param=bump_param, list_microcircuit = temp)

    return layer1
end


end
using .global_network



module Stimuli

using Parameters

export create_deterministic_oddball, create_deterministic_oscillations_oddball

function create_deterministic_oddball(param::Dict{String,Float64}, dt=0.0005)
    # this function will take into account ISI, strength, number of repetition

    @unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = param



    array_stim_submodule1 =  Vector{Float64}()
    array_stim_submodule2 =  Vector{Float64}()
    
    temp_index_stim = 0
    while temp_index_stim < number_repetitions
            
        
        array_stim_submodule1 = [array_stim_submodule1; zeros(round(Int, Tinter / dt)) ; value_stim_f1 * ones(round(Int, Tstim / dt))]
                
        array_stim_submodule2 = [array_stim_submodule2 ; zeros(round(Int, Tinter / dt)) ; value_stim_f2 * ones(round(Int, Tstim / dt))]
                
        temp_index_stim += 1
                    
    end
    array_stim_submodule1 = [array_stim_submodule1  ; zeros(round(Int, Tinter / dt)) ; value_stim_f2 * ones(round(Int, Tstim / dt))]
                
    array_stim_submodule2 = [array_stim_submodule2 ; zeros(round(Int, Tinter / dt)) ; value_stim_f1 * ones(round(Int, Tstim / dt))]

    array_stim_submodule1 = [0.0 * ones(round(Int, initial_time / dt));array_stim_submodule1;0.0 * ones(round(Int, initial_time / dt))]
    array_stim_submodule2 = [0.0 * ones(round(Int, initial_time / dt));array_stim_submodule2; 0.0 * ones(round(Int, initial_time / dt))]


    return array_stim_submodule1, array_stim_submodule2

end

function create_deterministic_oscillations_oddball(param::Dict{String,Float64}, dt=0.0005)
    # this function will take into account ISI, strength, number of repetition
    # this function will create an oddball paradigm with oscillations

    @unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = param



    array_stim_submodule1 =  Vector{Float64}()
    array_stim_submodule2 =  Vector{Float64}()
    
    temp_index_stim = 0
    while temp_index_stim < number_repetitions
            
        
        array_stim_submodule1 = [array_stim_submodule1; zeros(round(Int, Tinter / dt)) ; value_stim_f1 * ones(round(Int, Tstim / dt));zeros(round(Int, Tinter / dt)) ; value_stim_f2 * ones(round(Int, Tstim / dt))]
                
        array_stim_submodule2 = [array_stim_submodule2 ; zeros(round(Int, Tinter / dt)) ; value_stim_f2 * ones(round(Int, Tstim / dt)); zeros(round(Int, Tinter / dt)) ; value_stim_f1 * ones(round(Int, Tstim / dt))]
                
        temp_index_stim += 1
                    
    end
    array_stim_submodule1 = [array_stim_submodule1  ; zeros(round(Int, Tinter / dt)) ; value_stim_f2 * ones(round(Int, Tstim / dt))]
                
    array_stim_submodule2 = [array_stim_submodule2 ; zeros(round(Int, Tinter / dt)) ; value_stim_f1 * ones(round(Int, Tstim / dt))]

    array_stim_submodule1 = [0.0 * ones(round(Int, initial_time / dt));array_stim_submodule1;0.0 * ones(round(Int, initial_time / dt))]
    array_stim_submodule2 = [0.0 * ones(round(Int, initial_time / dt));array_stim_submodule2; 0.0 * ones(round(Int, initial_time / dt))]


    return array_stim_submodule1, array_stim_submodule2

end


end
using .Stimuli



end