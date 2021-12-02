module NetworkConstruction
using DrWatson
using JLD2
using Parameters

using ..BasicFunctions

export create_network, construct_two_local_microcircuit_integrator, create_deterministic_oscillations_oddball

export construct_two_local_microcircuit_integrator_full_param
export create_deterministic_oddball
export orientation_kernel, create_layer_bump

export construct_one_local_microcircuit_integrator
export connect_two_microcircuit!





module local_microcircuit_network

using Base:list_append!!
using ...BasicFunctions

using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.Synapses

using ...NeuronalStructures.Simulations


using DrWatson

export construct_two_local_microcircuit_integrator, construct_one_local_microcircuit_integrator

"""
    construct_one_local_microcircuit_integrator(name, parameters, synaptic_strengths)

Return a microcircuit constucted according to the following arguments:

- `name::String` = name of the microcircuit
- `param_micorcircuit` = parameters of the microcircuit
- `param_syn_strength_microcircuit` = parameters of the strengths of the synapses

"""
function construct_one_local_microcircuit_integrator(name::String; param_microcircuit = parameters_microcircuit(),  param_syn_strength_microcircuit = parameters_syn_strength_microcircuit() , param_interneurons = parameters_interneurons())

    @unpack τA, gA, dend_param, sst_adaptation, soma_adaptation, pv_to_soma_depression, soma_to_vip_facilitation, sst_to_vip_facilitation, soma_to_sst_facilitation, vip_to_sst_facilitation, soma_to_pv_depression, int_to_vip_depression, int_to_pv_depression, int_to_sst_facilitation, int_to_dend_depression, integrator_tc, time_tot, noise, int_to_sst_connection,  preferred, ngfc_to_dend_depression = param_microcircuit


    @unpack   gaba_sst_to_dend, gaba_pv_to_soma, nmda_soma_to_soma, nmda_soma_to_vip,
    gaba_sst_to_vip, nmda_soma_to_sst, gaba_vip_to_sst, nmda_soma_to_pv, gaba_sst_to_pv, gaba_pv_to_pv, nmda_soma_to_sst, nmda_soma_to_int, nmda_int_to_dend, gaba_sst_to_ngfc, gaba_ngfc_to_dend = param_syn_strength_microcircuit

    @unpack Ibg_vip, Ibg_sst, Ibg_pv, Ibg_ngfc, top_down_to_interneurons = param_interneurons



    if noise

        # create the OU process of the cells
        vip_OU = create_OU_process(time_tot)
        pv_OU = create_OU_process(time_tot)
        sst_OU = create_OU_process(time_tot)
        pc_OU = create_OU_process(time_tot)
        dend_OU = create_OU_process(time_tot)
        integrator_OU = create_OU_process(time_tot)
        ngfc_OU = create_OU_process(time_tot)
        
    else
        # create the OU process of the cells

        vip_OU = create_OU_process(time_tot, σ = 0.0)
        pv_OU = create_OU_process(time_tot, σ = 0.0)
        sst_OU = create_OU_process(time_tot, σ = 0.0)
        pc_OU = create_OU_process(time_tot, σ = 0.0)
        dend_OU = create_OU_process(time_tot, σ = 0.0)
        integrator_OU = create_OU_process(time_tot, σ = 0.0)
        ngfc_OU = create_OU_process(time_tot, σ = 0.0)

    end

    vip1 = vip_cell(name=string(name, "-", "vipcell1"), Ibg=Ibg_vip, OU_process=vip_OU)

    sst1 = sst_cell(name=string(name, "-", "sstcell1"), Ibg=Ibg_sst, adaptation_boolean=sst_adaptation, OU_process=sst_OU, adaptation = adaptation_variables(gA=gA, τA = τA))

    pv1 = pv_cell(name=string(name, "-", "pvcell1"), Ibg=Ibg_pv, OU_process=pv_OU)

    ngfc1 = ngfc_cell(name=string(name, "-", "ngfccell1"), Ibg=Ibg_ngfc, OU_process=ngfc_OU)

    dend1 = dend_sigmoid(param_c=dend_param, name=string(name, "-", "dend1"), OU_process=dend_OU)

    E1 = soma_PC(den=dend1, adaptation_boolean=soma_adaptation, name=string(name, "-", "ecell1"), OU_process=pc_OU, preferred_stim=preferred, adaptation = adaptation_variables(gA=gA, τA = τA))

    integrator1 = neural_integrator(τ=integrator_tc, name=string(name, "-", "integrator1"), OU_process=integrator_OU)


    dend1_gaba = gaba_syn(τ=10 * 0.001, g = gaba_sst_to_dend, name=string(name, "-", "sst1-to-dend1"))
    push!(dend1.list_syn_post_gaba, dend1_gaba)
    push!(sst1.list_syn_pre_gaba, dend1_gaba)

    E1_gabe = gaba_syn(g= gaba_pv_to_soma, depression=pv_to_soma_depression, name=string(name, "-", "pv1-to-ecell1"))
    push!(E1.list_syn_post_gaba, E1_gabe)
    push!(pv1.list_syn_pre_gaba, E1_gabe)


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

       
    ngfcsst_gaba = gaba_syn(g= gaba_sst_to_ngfc, name=string(name, "-", "sst1-to-ngfc1"))
    push!(ngfc1.list_syn_post_gaba, ngfcsst_gaba)
    push!(sst1.list_syn_pre_gaba, ngfcsst_gaba)

    dendngfc_gaba = gaba_syn(g= gaba_ngfc_to_dend, name=string(name, "-", "ngfc1-to-dend1"), depression=ngfc_to_dend_depression)
    push!(dend1.list_syn_post_gaba, dendngfc_gaba)
    push!(ngfc1.list_syn_pre_gaba, dendngfc_gaba)

    ngfcint_nmda = nmda_syn(g = top_down_to_interneurons[4], name=string(name, "-", "integrator1-to-ngfc1"))
    push!(ngfc1.list_syn_post_nmda, ngfcint_nmda)
    push!(integrator1.list_syn_pre_nmda, ngfcint_nmda)


    if int_to_sst_connection
        # sstint_nmda = nmda_syn(g=top_down_to_interneurons[3], facilitation=int_to_sst_facilitation, name=string(name, "-", "integrator1-to-sst1"))
        # push!(sst1.list_syn_post_nmda, sstint_nmda)
        # push!(integrator1.list_syn_pre_nmda, sstint_nmda)

        vipint_nmda = nmda_syn(g=top_down_to_interneurons[1], depression=int_to_vip_depression, name=string(name, "-", "integrator1-to-vip1"))
        push!(vip1.list_syn_post_nmda, vipint_nmda)
        push!(integrator1.list_syn_pre_nmda, vipint_nmda)
    
    
        pvint_nmda = nmda_syn(g=top_down_to_interneurons[2], depression=int_to_pv_depression, name=string(name, "-", "integrator1-to-pv1"))
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

        
    c = microcircuit{soma_PC, dend_sigmoid}(soma = E1, sst = sst1, vip = vip1, pv = pv1, ngfc = ngfc1, dend = dend1, integrator = integrator1, name = name )
    return c

end


"""
    construct_two_local_microcircuit_integrator(;parameters, synaptic_strengths)

Construct two inter-connected interneuron microcircuits
"""
function construct_two_local_microcircuit_integrator(;param_microcircuit = parameters_microcircuit(), param_syn_microcircuit = parameters_syn_strength_microcircuit(), param_inter_microcircuit = parameters_inter_microcircuit())


    @unpack dend_param, sst_adaptation, soma_adaptation, pv_to_soma_depression, soma_to_vip_facilitation, sst_to_vip_facilitation, soma_to_sst_facilitation, vip_to_sst_facilitation, soma_to_pv_depression, int_to_vip_depression, int_to_pv_depression, int_to_sst_facilitation, int_to_dend_depression, integrator_tc, time_tot, noise, int_to_sst_connection, top_down_to_interneurons = param_microcircuit

    @unpack cross_int_to_vip_depression, cross_int_to_pv_depression, cross_int_to_sst_facilitation, cross_int_to_dend_depression, cross_soma_to_sst_facilitation = param_inter_microcircuit

    
    c1 = construct_one_local_microcircuit_integrator("microcircuit1"; param_microcircuit = param_microcircuit, param_syn_strength_microcircuit = param_syn_microcircuit)

    c2 = construct_one_local_microcircuit_integrator("microcircuit2"; param_microcircuit = param_microcircuit, param_syn_strength_microcircuit = param_syn_microcircuit)


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

end
using .local_microcircuit_network



module global_network

# TODO to make better these modules
using StaticArrays
using Base:list_append!!
using ...BasicFunctions
using ...NeuronalStructures.NeuralNetwork
using ..local_microcircuit_network
using ...NeuronalStructures.Simulations
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.Synapses
using ...NeuronalStructures.EqDiffMethod

using DrWatson

export orientation_kernel, create_layer_bump
export connect_two_microcircuit!

"""
    orientation_kernel(x,y,parameters_bump_attractor)

Compute the orientation of a ring model:

``K = \\exp{-(x-y)^2/\\sigma^2 }``

Note that in `Engel 2015`, ``\\sigma = 43.2`` degrees and there are 128 units uniformly distributed 
"""
function orientation_kernel(x::Float64, y::Float64, param::parameters_bump_attractor)
    min_dif = abs(x - y)
    if min_dif > 180.0
        min_dif -= 360.0
    end

    return exp(-min_dif.^2 / param.σ.^2)

end


"""
    connect_two_microcircuit!(c1,c2,parameters_bump_attractor, parameters_inter_microcircuit)

Connect to microcircuit within a ring model together
"""
function connect_two_microcircuit!(c1::microcircuit, c2::microcircuit, bump_param::parameters_bump_attractor,param_inter_microcircuit::parameters_inter_microcircuit)
        

    @unpack cross_int_to_vip_depression, cross_int_to_pv_depression, cross_int_to_sst_facilitation, cross_int_to_dend_depression, cross_soma_to_sst_facilitation = param_inter_microcircuit

    temp = orientation_kernel(c1.soma.preferred_stim, c2.soma.preferred_stim, bump_param)
    e1tosst2 = nmda_syn(g=0.00435 * temp, facilitation=cross_soma_to_sst_facilitation, name=string("ecell2-to-sst1"))
    push!(c2.sst.list_syn_post_nmda, e1tosst2)
    push!(c1.soma.list_syn_pre_nmda, e1tosst2)

    e2tosst1 = nmda_syn(g=0.00435 * temp, facilitation=cross_soma_to_sst_facilitation, name=string("ecell2-to-sst1"))
    push!(c1.sst.list_syn_post_nmda, e2tosst1)
    push!(c2.soma.list_syn_pre_nmda, e2tosst1)

    dend1int2 = nmda_syn(g=0.1 * temp, depression=cross_int_to_dend_depression, name=string("integrator2-to-ecell1"))
    push!(c1.dend.list_syn_post_nmda, dend1int2)
    push!(c2.integrator.list_syn_pre_nmda, dend1int2)

    dend2int1 = nmda_syn(g=0.1 * temp, depression=cross_int_to_dend_depression, name=string("integrator1-to-ecell2"))
    push!(c2.dend.list_syn_post_nmda, dend2int1)
    push!(c1.integrator.list_syn_pre_nmda, dend2int1)

    

    nothing

end
    
"""
    create_parameters_microcircuit(i::Int64, param_inter_microcircuit::parameters_inter_microcircuit)

create the parameters of a microcircuit for stimuli index i
"""
function create_parameters_microcircuit(i::Int64, param_microcircuit::parameters_microcircuit, bump_param::parameters_bump_attractor)
    param_microcircuit.preferred = i*360.0/bump_param.num_circuits

    return param_microcircuit
end


"""
    create_layer_bump(parameters_bump_attractor, num_circuits; parameters_microcircuit, parameters_inter_microcircuit)

Return a ring model composed of `N` microcircuits
    """
function create_layer_bump(bump_param::parameters_bump_attractor, num_circuits::Int64; param_microcircuit = parameters_microcircuit(), param_inter_microcircuit = parameters_inter_microcircuit(), param_syn_strength_microcircuit = parameters_syn_strength_microcircuit())


    temp = @SVector [ construct_one_local_microcircuit_integrator("microcircuit$i";param_microcircuit = create_parameters_microcircuit(i,param_microcircuit,bump_param), param_syn_strength_microcircuit = param_syn_strength_microcircuit) for i=1:128]
    layer1 = layer_bump{soma_PC, dend_sigmoid,euler_method}(bump_param=bump_param, list_microcircuit = temp, eq_diff_method = euler_method())


    for i = 1:(num_circuits - 1)
        for j = (i + 1):num_circuits
            connect_two_microcircuit!(layer1.list_microcircuit[i], layer1.list_microcircuit[j], bump_param, param_inter_microcircuit)
        end
    end

    return layer1
end


end
using .global_network



module Stimuli

using Parameters
using ...NeuronalStructures.NeuralNetwork
using ...NeuronalStructures.Simulations
using ...NeuronalStructures.NeuronalModels
using ...NeuronalStructures.Synapses
using ...NeuronalStructures.EqDiffMethod
using ...NeuronalStructures.AbstractNeuronalTypes

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

function create_deterministic_oddball(layer_bump::layer_bump{soma_PC, dend_sigmoid, euler_method},param::Dict{String,Float64})
    # this function will take into account ISI, strength, number of repetition

    @unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = param


    dt = layer_bump.eq_diff_method.dt


    total_time = round(Int,((intitial_time + number_repetitions*(Tinter+Tstim) + Tstim)/dt))

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


"""
create the stimulation protocol for the oddball paradigm for one neuron
"""
function create_deterministic_oddball(neuron::T where T <: neuron, param::Dict{String,Float64}, euler_method::eq_diff_method)
    # this function will take into account ISI, strength, number of repetition

    @unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = param

    dt = euler_method.dt



end

"""

create the simulation currents from a deterministic oddball
"""
function create_deterministic_oddball(layer_bump::layer_bump{soma_PC, dend_sigmoid, euler_method},param::Dict{String,Float64})
    # this function will take into account ISI, strength, number of repetition

    @unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = param

    dt = layer_bump.eq_diff_method.dt

    total_time = round(Int,((intitial_time + number_repetitions*(Tinter+Tstim) + Tstim)/dt))



end


end
using .Stimuli



end