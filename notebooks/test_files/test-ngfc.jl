using DrWatson
	
quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
push!(LOAD_PATH,srcdir())

using Revise
using MyNeurosciencePackage
using CairoMakie
using Parameters, JLD2, Statistics	
using DataFrames, PlutoUI
#using LsqFit
using ColorSchemes, StaticArrays

set_theme!(theme_light())

sim = simulation_parameters()
ring_simu = Dict{String,Vector{Float64}}()
time_tot = 15000
bump = create_layer_bump(parameters_bump_attractor(),128, param_microcircuit = parameters_microcircuit(time_tot = time_tot, soma_to_pv_depression = false, sst_to_vip_facilitation = false, int_to_vip_depression = true, vip_to_sst_facilitation = true, soma_to_sst_facilitation = true, ngfc_to_dend_depression = true, int_to_dend_depression = true), param_syn_strength_microcircuit = parameters_syn_strength_microcircuit(nmda_soma_to_int = 0.15, nmda_soma_to_pv = 0.17, nmda_soma_to_vip = 0.3, gaba_vip_to_sst = -0.05, gaba_ngfc_to_dend = -0.15, gaba_sst_to_ngfc = -0.4 ))
stim = 50.0

stimlist = vcat(zeros(5000),ones(1000), zeros(1000), ones(1000), zeros(1000), ones(1000), zeros(1000), ones(1000), zeros(1000),ones(1000), zeros(1000))



for i=1:128
   
   sim.current[string("microcircuit$i","-","ecell1")] = 0.2.*stimlist.*orientation_kernel(stim,bump.list_microcircuit[i].soma.preferred_stim, bump.bump_param)
   sim.current[string("microcircuit$i","-","dend1")] = 0.0.*stimlist
   sim.current[string("microcircuit$i","-","vipcell1")] = 0.0.*stimlist
   sim.current[string("microcircuit$i","-","sstcell1")] = 0.0.*stimlist
   sim.current[string("microcircuit$i","-","pvcell1")] = 0.0.*stimlist
   sim.current[string("microcircuit$i","-","ngfccell1")] = 0.0.*stimlist

   sim.current[string("microcircuit$i","-","integrator1")] = 0.0.*stimlist

   
end

full_time_dynamics!(bump, sim)

# need to do the maths of why vip is too low asfiring rates
lines(get_firing_rate(bump,1:1:7000, "soma", 20))

lines(get_firing_rate(bump,1:1:15000, "pv", 20))

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax,get_firing_rate(bump,1:1:15000, "pv", 20))

lines!(ax,get_firing_rate(bump,1:1:15000, "sst", 20))
fig

lines!(ax,get_firing_rate(bump,1:1:15000, "vip", 20))
lines!(ax,get_firing_rate(bump,1:1:15000, "soma", 20))
lines!(ax,get_firing_rate(bump,1:1:15000, "ngfc", 20))

fig

lines(get_firing_rate(bump,1:1:15000, "integrator", 20))

####
# Seems that VIP and SST are not necessarily active at the same time during ISO and Cross stimuli

# depression from ngfc to dendrits seems to be the reason why non adaptation of soma firing rates
# the only interesting fact would be that it allows ap erfect consideratoin of a learinng rule as the network is almost perfectly stable upon presentation of a succession of Stimuli
# TODO: think about what kind of learning rule would be nice to have
####

######
# How to find the starting state of the network without any stimulation
#Probbaly to chekc better solutions
# maybe synapses to add???
# understand why this script doesn't work
####

using JuMP
using Ipopt

model = Model(Ipopt.Optimizer)

# define the variables for the firing rates of the network
@variable(model, rsoma >= 0.0)
@variable(model, rint >= 0.0)
@variable(model, rngfc >= 0.0)
@variable(model, rpv >= 0.0)
@variable(model, rvip >= 0.0)
@variable(model, rsst >= 0.0)

# variable for the current going into the network
@variable(model, Idend)
@variable(model, Isoma)
@variable(model, Iint)
@variable(model, Ipv)
@variable(model, Ingfc)
@variable(model, Ivip)
@variable(model, Isst)

# variables for the synapses of the network
@variable(model, syn_soma)
@variable(model, syn_int)
@variable(model, syn_pv)
@variable(model, syn_ngfc)
@variable(model, syn_vip)
@variable(model, syn_sst)

# constraintes on the synaptic rates of the network
γnmda = 0.641 * 2.0
τnmda = 60.0*0.001

γgaba =  2.0
τgaba =  0.005

@NLconstraint(model, syn_soma == τnmda * γnmda*rsoma/(1.0 + τnmda * γnmda*rsoma)) 
@NLconstraint(model, syn_int == τnmda * γnmda*rint/(1.0 + τnmda * γnmda*rint))

@NLconstraint(model, syn_pv == τgaba * γgaba * rpv)
@NLconstraint(model, syn_ngfc == τgaba * γgaba * rngfc)
@NLconstraint(model, syn_vip == τgaba * γgaba * rvip)
@NLconstraint(model, syn_sst == τgaba * γgaba * rsst)

@NLconstraint(model, Ipv == 0.29 + 0.17*syn_soma - 0.18*syn_pv - 0.17*syn_sst + 0.34*syn_int)
@NLconstraint(model, Ivip == 0.29 + 0.058*syn_soma  - 0.1*syn_sst + 0.66*syn_int)
@NLconstraint(model, Isst == 0.25 + 0.0635*syn_soma - 0.17*syn_vip)
@NLconstraint(model, Ingfc == 0.25 + 0.5*syn_int - 0.2*syn_sst)
@NLconstraint(model, Iint == 0.31 + 0.15*syn_soma)
@NLconstraint(model, Isoma == 0.31 + 0.15*syn_soma + Idend - 0.1*syn_pv)
@NLconstraint(model, Idend == 0.4*syn_int -0.09*syn_sst -0.2*syn_ngfc + 0.1)

a = 135.0
b = 54.0
c = 0.308

@NLconstraint(model, rsoma == (a*Isoma - b)/(1.0 - exp(-c*(a*Isoma-b))))
@NLconstraint(model, rint == Iint/0.9)
@NLconstraint(model, rngfc == 132.0*Ingfc - 33.0 )
@NLconstraint(model, rpv == 330.0*Ipv - 95.0)
@NLconstraint(model, rvip == 132.0*Ivip - 33.0)
@NLconstraint(model, rsst == 132.0*Isst - 33.0)


JuMP.optimize!(model)
print(model)


solution_summary(model)


value(rsoma)
value(rint)
value(rngfc)
value(rsst)
value(rvip)
value(rpv)

value(syn_soma)
value(syn_int)

value(syn_pv)
value(syn_ngfc)
value(syn_vip)
value(syn_sst)

value(Idend)