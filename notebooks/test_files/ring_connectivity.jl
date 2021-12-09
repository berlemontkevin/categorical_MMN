using DrWatson
	
quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
push!(LOAD_PATH,srcdir())

using Revise
using MyNeurosciencePackage

using CairoMakie
using Parameters, JLD2, Statistics


param_bump = parameters_bump_attractor()
param_micro = parameters_microcircuit()
param_inter = parameters_interneurons()
param_syn = parameters_syn_strength_microcircuit()
param_int_m = parameters_inter_microcircuit()

ring_m = create_ring_microcircuit(param_bump, param_micro, param_inter)

intra_connectivity!(ring_m, param_micro, param_syn)

ring_a = create_ring_attractor(param_bump, param_micro)

ring_m.neurons[1].sst

ring_a.neurons[1]

inter_connectivity!(ring_a, param_syn, param_int_m, param_bump)

ring_a.neurons[1].list_syn_post_nmda[1].g

temp_strength = zeros(63)
for i=1:63
    temp_strength[i] = ring_a.neurons[1].list_syn_post_nmda[i].g
end

scatter(temp_strength)

inter_connectivity!(ring_m, param_syn, param_int_m, param_bump)

temp_str_soma = zeros(127)
ring_m.neurons[1].soma
for i=1:127
    t = floor(Int, i+5)
    # print(t)
    temp_str_soma[i] = ring_m.neurons[1].soma.list_syn_pre_nmda[t].g
end

scatter(temp_str_soma)


inter_connectivity!(ring_m, ring_a, param_syn,param_inter, param_int_m, param_bump,param_micro)


ring_m.neurons[1].ngfc
strength_ngfc = zeros(64)
for i=1:64
    strength_ngfc[i] = ring_m.neurons[1].ngfc.list_syn_post_nmda[i].g
end
scatter(strength_ngfc)

nn = neural_network{soma_PC, dend_sigmoid, euler_method}(ring_integrator = ring_a, list_microcircuit = ring_m)

nn.ring_integrator


nn2 = create_two_layers()




sim = simulation_parameters()

ring_simu = Dict{String,Vector{Float64}}()
time_tot = 4500

nn = create_two_layers(param_syn = parameters_syn_strength_microcircuit(nmda_soma_to_int = 0.01, nmda_soma_to_pv = 0.17, nmda_soma_to_vip = 0.35, gaba_vip_to_sst = -0.05, gaba_ngfc_to_dend = -0.15, gaba_sst_to_ngfc = -0.4, nmda_soma_to_sst =  0.0435, nmda_int_to_dend = 0.05*0.4),
param_micro = parameters_microcircuit(time_tot = time_tot, soma_to_pv_depression = false, sst_to_vip_facilitation = true, int_to_vip_depression = true, vip_to_sst_facilitation = true, soma_to_sst_facilitation = true, ngfc_to_dend_depression = true, int_to_dend_depression = true,noise =true))


stim = 50.0
stimlist = vcat(zeros(500),ones(1000), zeros(1000), ones(1000), zeros(1000))#, ones(1000), zeros(1000), ones(1000), zeros(1000),ones(1000), zeros(1000))

for i=1:128
   
   sim.current[string("microcircuit$i","-","ecell1")] = 0.2.*stimlist.*orientation_kernel(stim,nn.list_microcircuit.neurons[i].soma.preferred_stim, nn.list_microcircuit.parameters)
   sim.current[string("microcircuit$i","-","dend1")] = 0.0.*stimlist
   sim.current[string("microcircuit$i","-","vipcell1")] = 0.0.*stimlist
   sim.current[string("microcircuit$i","-","sstcell1")] = 0.0.*stimlist
   sim.current[string("microcircuit$i","-","pvcell1")] = 0.0.*stimlist
   sim.current[string("microcircuit$i","-","ngfccell1")] = 0.0.*stimlist


   
end
for i = 1:64
    sim.current[string("integrator","-","$i")] = 0.0.*stimlist

end


full_time_dynamics!(nn, sim)
@code_warntype full_time_dynamics!(nn, sim)
@code_warntype time_step!(nn, sim,10)

bump_animation(nn,2,projectdir("notebooks","animation\\"),"test_ring_connectivity",".mp4")

lines(get_firing_rate(nn,1:1:4500, "soma", 20))

lines(get_firing_rate(nn,1:1:4500, "pv", 20))

lines(get_firing_rate(nn,1:1:4500, "sst", 20))

lines(get_firing_rate(nn,1:1:4500, "vip", 20))

lines(get_firing_rate(nn,1:1:4500, "ngfc", 20))

lines(get_firing_rate(nn.ring_integrator,1:1:4500,  18))

a=nn.ring_integrator.neurons[18].Iexc_save
lines(a[1:2100])

b = nn.list_microcircuit.neurons[18].dend.Itot_save
lines(b)


c = nn.list_microcircuit.neurons[18].soma.Ioutput_save
lines(c)