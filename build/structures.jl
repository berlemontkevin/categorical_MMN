using Parameters


#####
# Robert 2016
#########

"""
Type: eq_diff_method

Defines the type of method that will be used to solve the differential equations

...
# Types currently implemented
- Euler method
....
"""
abstract type eq_diff_method end


abstract type neuron end

abstract type synapses end

abstract type dendrite <: neuron end
abstract type pyr_cell <: neuron end
abstract type interneuron <: neuron end







@with_kw struct euler_method <: eq_diff_method
    dt::Float64
    record::Bool = true
end



@with_kw struct syn_connection <: synapses
    #structure d'une connexion synaptique (mutable pour l'instant)
    pre_syn::neuron
    w::Vector{Float64}
end


@with_kw struct rectified_linear_neurons <: neuron
    # structure for the VIP neurons
    Isyn::Vector{Float64} = [0.0] # synaptic input toward this neuron
    τ::Float64 = 2.0 #ms #for now same for all interneurons
    neurons_list::Vector{syn_connection} = syn_connection[] # list of all presynaptic connections
    fr::Vector{Float64} = [0.0] # current firing rate of the neuron
    Ibg::Vector{Float64} = [0.0] #background input

end


# @with_kw mutable struct som_neurons <: interneuron
#     # som interneurons as rate units
#     Irh::Float64 = 40.0 # pA # It's the rheobase of the som neuron
#     β::Float64 = 90.0 #Hz/nA # f-I curve slope # taken from Spike timing dependent plasticity 2007
#     Isyn::Float64 = 0.0 # synaptic current toward the neuron
#     neurons_list::Array{syn_connection,1} = syn_connection[] # list of connected neurons
#     τ::Float64 = 2.0 #ms #for now same for all interneurons
#     Ibg::Float64 = 0.0 #background input

# end

# @with_kw mutable struct vip_neurons <: interneuron
#     # structure for the VIP neurons
#     Isyn::Float64 = 0.0
#     τ::Float64 = 2.0 #ms #for now same for all interneurons
#     neurons_list::Array{syn_connection,1} = syn_connection[]
#     fr::Float64 = 0.0
#     Ibg::Float64 = 0.0 #background input

# end

# @with_kw mutable struct pv_neurons <: interneuron
#     # structure of PV neurons
#     Isyn::Float64 = 0.0
#     τ::Float64 = 2.0 #ms #for now same for all interneurons
#     neurons_list::Array{syn_connection,1} = syn_connection[]
#     fr::Float64 = 0.0
#     Ibg::Float64 = 0.0 #background input

# end










@with_kw mutable struct dendrites <: neuron
    # dendrites structures
    β::Float64 = 0.0
    V0::Float64 = 0.78 #mV
    EL::Float64 = 0.0
    bg::Float64 = 5.56
    g12::Float64 = 0.0
    gLD::Float64 = 4.0 #leak conductance of the dendrite
    k::Float64 = 9.64 #nS
    γ::Float64 = 6.54 #nS
    gE::Float64 = 0.0# excitatory input to the dendrites
    rI::Float64 = 0.0# inhibitory input to the dendrites
    gI::Float64 = 0.0
    VD::Float64 = 0.0
    snmda::Float64=0.0
    rE::Float64 = 30.0
    neurons_list::Array{syn_connection,1} = syn_connection[]
end

@with_kw struct syn_GABA_dendrites <: synapses
    # structure of GABA synapses into dendrites
    τGABA::Float64 = 0.020 #s
    gGABAmax::Float64 = 4.0 #nS
end

@with_kw struct syn_NMDA_dendrites <: synapses
    # structure of GABA synapses into dendrites
    τNMDAdec::Float64 = 100 #ms
    τNMDArise::Float64 = 2 #ms
    αNMDA::Float64 = 0.3
end





@with_kw mutable struct firing_rate_soma
    # parameters of the firing rate function of the soma
    # if the function is a power law
    α::Float64 = 2.89 # exposant
    β::Float64 = 174.386 # somme dans le max
    γ::Float64 = 45.16 # normalisation

end

@with_kw mutable struct pyr_cells <: neuron
    # pyr cells as rate units
    gE::Float64 = 0.0
    gI::Float64 = 0.0
    r::Float64 = 0.0# firing rate of the soma
    I::Float64 = 0.0# total input current
    ID::Float64 = 0.0# dendrites input to soma (averaged)
    Isyn::Float64 = 0.0 #synaptic current
    dendrites_list::Array{dendrites} = dendrites[] # list of the dendrites connected to this pyr_cell
    neurons_list::Array{syn_connection,1} = syn_connection[]
    param_firingrate::firing_rate_soma = firing_rate_soma()
    Ibg::Float64 = 0.0 #background input

end








# function update_fr!(pr::pyr_cells)
#     # function that updates the firing rate for the pyramidal cells
#     # it updates the somatic FR
#     temp = max(0.0, pr.I + pr.param_firingrate.β)/pr.param_firingrate.γ
#     pr.r = exp(fr_soma.α * log(temp))
# end

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

# function update_I_DS!(p::pyr_cells)
#     temp = 0.0
#     for i in p.dendrites_list
#         temp = temp + i.VD
#     end
#     temp = temp/length(p.dendrites_list)
#     p.I = 8.0 * (temp + 55.0 )
# end


# function update_fr!(som::som_neurons,euler_m::euler_method)
#     # in tis diff eq, the total FR is rectified
#     som.syn += rect_linear(euler_m.dt/som.τ * (-som.fr + som.Isyn + som.Iext))
# end

# function update_fr!(pv::pv_neurons,euler_m::euler_method)
#     # in tis diff eq, the total FR is rectified
#     pv.syn += rect_linear(euler_m.dt/pv.τ * (-pv.fr + pv.Isyn + pv.Iext))
# end

# function update_fr!(vip::vip_neurons,euler_m::euler_method)
#     vip.syn += rect_linear(euler_m.dt/vip.τ * (-vip.fr + vip.Isyn + vip.Iext))
# end


# function update_syn_current!(n::neuron)
#     # update the synaptic current towards one neuron
#     # Rule: The sign of weights is taken into account in w
#     n.Isyn = 0.0
#     for syn in n.neurons_list
#         n.Isyn += syn.pre_syn.fr * syn.w
#     end
#     n.Isyn += n.Ibg #add the background input to the syn current

# end



# function simulation(net::network)
#     # function that perform a loop of the simulation
#     # take a whole network as argument


# end


#############################################################
################# Hertag Sprekler 2020
#############################################################


@with_kw struct objective_rate
    fr::Float64

end


@with_kw struct soma_hertag  <: pyr_cell
    # somatic structure of Hertag 2020
    τE::Float64 = 0.060 #ms
    fr::Vector{Float64} = [0.0] # firing rate of somatic compartment
    Isyn::Vector{Float64} = [0.0]# synaptic input to the soma
    Ibg::Vector{Float64} = [0.0] # background input
    Itot::Vector{Float64} = [0.0] #total synaptic input
    λD::Float64 = 0.27
    ID::Vector{Float64} = [0.0]# total synaptic input into dendrites
    λE::Float64 = 0.31
    Θr::Float64 = 14.0 # Hz rheobase of the neuron
    IE::Vector{Float64} = [0.0]# total synaptic inputs into soma


    objective_rate::Float64 = 1.0 #Hz

    
    dendrites_list::Vector{dendrite} = dendrites[] # list of the dendrites connected to this pyr_cell
    neurons_list::Vector{syn_connection} = syn_connection[]


end

@with_kw  struct dendrites_hertag <: dendrite
    # structure of hertag dendrites
    neurons_list::Array{syn_connection,1} = syn_connection[]
    Θc::Float64 = 28.0 # Hz minimal input to produce a calcium spike
    Itot::Vector{Float64} = [0.0]# total synaptically generated input in the dendrites
    gain_calcium::Float64= 7 # Hz
    Isyn::Vector{Float64} = [0.0]# synaptic input (correspond to ID)
    pyr_target::Array{soma_hertag} = soma_hertag[] # target of the dendrites
    Ibg::Vector{Float64} = [0.0]#background input to the dendrite
    calcium_spike::Vector{Float64} = [0.0]# input from calclium spike

end


@with_kw struct som_neurons_hertag <: interneuron
    # som interneurons as rate units
    Isyn::Vector{Float64} = [0.0]
    τ::Float64 = 0.002 #ms #for now same for all interneurons
    neurons_list::Array{syn_connection,1} = syn_connection[]
    fr::Vector{Float64} = [0.0]
    Ibg::Vector{Float64} = [0.0] #background input
    objective_rate::Float64 = 2.0 #Hz


end

@with_kw  struct vip_neurons_hertag <: interneuron
    # structure for the VIP neurons
    Isyn::Vector{Float64} = [0.0]
    τ::Float64 = 0.002 #ms #for now same for all interneurons
    neurons_list::Array{syn_connection,1} = syn_connection[]
    fr::Vector{Float64} = [0.0]
    Ibg::Vector{Float64} = [0.0] #background input
    objective_rate::Float64 = 4.0 #Hz


end

@with_kw struct pv_neurons_hertag <: interneuron
    # structure of PV neurons
    Isyn::Vector{Float64} = [0.0]
    τ::Float64 = 0.002 #ms #for now same for all interneurons
    neurons_list::Array{syn_connection,1} = syn_connection[]
    fr::Vector{Float64} = [0.0]
    Ibg::Vector{Float64} = [0.0] #background input
    objective_rate::Float64 = 2.0 #Hz


end



function update_fr!(n::rectified_linear_neurons,euler_m::euler_method)
    # in tis diff eq, the total FR is rectified
    # the record variables indicates if n.fr is a list or just the last firing rate
    if record
        push!(n.fr,n.fr[end])
        n.fr[end] += euler_m.dt/n.τ * (-n.fr[end] + n.Isyn[end] )
        n.fr[end] = rect_linear(n.fr[end])
    else
        n.fr[end] += euler_m.dt/n.τ * (-n.fr[end] + n.Isyn[end] )
        n.fr[end] = rect_linear(n.fr[end])
    end

end





function update_current!(pyr_cell::soma_hertag)
    # function that updates the update the total somatic input toward the soma
    d = pyr_cell.dendrites_list[1] #only 1 dendrite connected
    pyr_cell.Itot = pyr_cell.λD * rect_linear(d.Isyn + d.calcium_spike) + 
    (1.0 - pyr_cell.λE)*pyr_cell.Isyn

end

function compute_calcium_spike!(dend::dendrites_hertag)
    # function that computes the calcium dendritic event
    dend.calcium_spike = dend.gain_calcium * heaviside(dend.Itot - dend.Θc)
    end


function update_current!(d::dendrites_hertag)
    # function that updates the total synaptic input generated in the dendrites
    p = d.pyr_target[1]
    d.Itot = p.λE * p.Isyn + (1.0 - p.λD) * d.Isyn
end

mutable struct network_hertag
    #structure of the whole network
    pv_list::Array{pv_neurons_hertag,1}
    pc_list::Array{pyr_cell,1}
    sst_list::Array{som_neurons_hertag,1}
    vip_list::Array{vip_neurons_hertag,1}
    dendrites_list::Array{dendrite,1}

end

##TODO; need to change the structure to be able to be muttable 
function simulation_step!(netw::network_hertag,euler_m::euler_method)

    for n in netw.pv_list
        update_syn_current!(n) #of everyone
    end
    for n in netw.pc_list
        update_syn_current!(n) #of everyone
    end
    for n in netw.sst_list
        update_syn_current!(n) #of everyone
    end
    for n in netw.vip_list
        update_syn_current!(n) #of everyone
    end
    for n in netw.dendrites_list
        update_syn_current!(n) #of everyone
        update_current!(n) #all dendrites
    compute_calcium_spike!(n) # compute all calcium spikes)
    end

    for n in netw.pc_list
        update_current!(n) # all pyr cells
    end


    for n in netw.pv_list
        update_fr!(n,euler_m) # update all firing rates
    end
    for n in netw.pc_list
        update_fr!(n,euler_m) # update all firing rates
    end
    for n in netw.sst_list
        update_fr!(n,euler_m) # update all firing rates
    end
    for n in netw.vip_list
        update_fr!(n,euler_m) # update all firing rates
    end
    


end




######################################################


######################################################


######################################################

#####
# Hertag 2020
#########

#rate dendrites
# time averaged voltage
# firing rate soma
# 



# une classe PC qui a tout en computes
# pour les connexions il y aura une liste qui donnera les indices des connexions
# les constantes qui ne bougent pas seront surement a mettre en constantes globales (optimisation plus tard)


# a verifier un peu mieux ce aque font leurs dendrites mais c'est bizarre

####
# Dendrites Yang 2016
# ######

# @with_kw mutable struct pyr_cell_Yang2016
#     V0::Float64=0.78 #mV
#     EL::Float64= -70 #mV
#     β::Float64
#     g12::Float64
#     bg::Float64 = 5.56
#     k::Float64 = 9.64 #nS
#     γ::Float64 = 6.54 #nS
#     Gc::Float64 = 8 #nS
#     r::Float64 # rate of the soma
#     VD::Float64 # voltage dendritic compartment
#     gI::Float64 # conductance of inhibitory neurons
#     gE::Float64 # conductance of excitatory neurons
#     Ereset::Float64 = 0.0


# end


# function rate_soma!(pyr_cell_Yang2016)
#     # function that computes the firing rate of the soma 
#     temp = max(0.0, I + 174.86)/45.16
#     pyr_cell_Yang2016.r = exp(2.89 * log(temp))
# end

# function update_g12!(pyr_cell_Yang2016)
#     # update g12 for the dendrites potential
#     pyr_cell_Yang2016.g12 = pyr_cell_Yang2016.bg * (pyr_cell_Yang2016.gLD + pyr_cell_Yang2016.gI)
# end


# function update_VD(pyr_cell_Yang2016)
#     # update the voltage average of dedrites compartment
#     temp = tanh((pyr_cell_Yang2016.gE - pyr_cell_Yang2016.g12)/pyr_cell_Yang2016.β)
#     pyr_cell_Yang2016.VD = 30 * (1.0 + temp) + pyr_cell_Yang2016.V0 + pyr_cell_Yang2016.EL
# end

# function update_β!()
#     # update β depending on the inhibitory conductance
#     pyr_cell_Yang2016.β = pyr_cell_Yang2016.k * exp(pyr_cell_Yang2016.gI / γ)
# end







# function normalisation_weights(S_weights::syn_weights_hertag , P_weights::syn_probabilities_hertag)
#     # ici on suppose un reseau 70 10 10 10
#     S_weights.wEP = S_weights.wEP/(P_weights.pEP*10)
#     S_weights.wDS = S_weights.wDS/(P_weights.pDS*10)
#     S_weights.wDE = S_weights.wDE/(P_weights.pDE*70)
#     S_weights.wPE = S_weights.wPE/(P_weights.pPE*70)
#     S_weights.wPP = S_weights.wPP/(P_weights.pPP*10)
#     S_weights.wPS = S_weights.wPS/(P_weights.pPS*10)
#     S_weights.wPV = S_weights.wPV/(P_weights.pPV*10)
#     S_weights.wSE = S_weights.wSE/(P_weights.pSE*70)
#     S_weights.wSV = S_weights.wSV/(P_weights.pSV*10)
#     S_weights.wVE = S_weights.wVE/(P_weights.pVE*70)
#     S_weights.wVS=S_weights.wVS/(P_weights.pVS*10)




# end



struct neural_population
    N::Int64 # number of neurons
    type_neuron::String # types of the neurons in this population
    list_neurons::Vector{neuron}
    type_global::String # can be pyr_cells, dendrites or interneurons for now
end

struct neural_network
    # this structure will contain the informations to construct the network
    list_pop::Vector{neural_population} #list of neural population
    m_prob::Array{Float64} # matrix of the probabilities of connectivity
    m_weights::Array{Float64} # matrix of the synaptic weights

end
