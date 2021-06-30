### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ a2d39c60-d06a-11eb-1364-c1e29d8f626f
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())

end

# ╔═╡ 1334881e-d06b-11eb-3ec8-13f9d3037c39
begin
	using Revise
	using MyNeurosciencePackage
	using Plots
	plotlyjs()
	using Parameters, JLD2, JuMP
end

# ╔═╡ e5810c90-d06b-11eb-26d3-e378a54e0bf7


# ╔═╡ 16aa4df0-d06b-11eb-3d89-6fafd70abb11
begin
	

#TODO add name field to microcircuit
function construct_one_local_microcircuit_integrator!(c::microcircuit; dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation = false, adaptation = false, depression = false, td_to_vip = true)
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name

    vip1 = vip_cell(name = string(c.name, "-","vipcell1"),Ibg = 0.25)
    sst1 = sst_cell(name = string(c.name, "-","sstcell1"), Ibg = 0.25)
    pv1 = pv_cell(name = string(c.name, "-","pvcell1"), Ibg = 0.27)
    dend1 = dend_sigmoid(param_c = dend_param,name = string(c.name, "-","dend1"))
    E1 = soma_PC(den=dend1, adaptation_boolean = adaptation,name = string(c.name, "-","ecell1"))
    integrator1 = neural_integrator(name = string(c.name, "-","integrator1"))
    
    
    push!(dend1.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst1,neuron_post=dend1, g = -0.09, name = string(c.name, "-","sst1-to-dend1")))#0.09
    push!(E1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E1, g = -0.005, name = string(c.name, "-","pv1-to-ecell1")))#-0.001
    push!(E1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=E1, g = 0.18, name = string(c.name, "-","ecell1-to-ecell1")))
    
    push!(vip1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=vip1, g = 0.058, name = string(c.name, "-","ecell1-to-vip1")))
    push!(vip1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=vip1, g = -0.1, name = string(c.name, "-","sst1-to-vip1")))
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst1, g = 0.0435, name = string(c.name, "-","ecell1-to-sst1")))
    push!(sst1.list_syn, gaba_syn(neuron_pre=vip1,neuron_post=sst1, g = -0.05, name = string(c.name, "-","ss1-to-dend1")))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=pv1, g = 0.08435, name = string(c.name, "-","ecell1-to-pv1")))#0.0435
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=pv1, g = -0.17, name = string(c.name, "-","sst1-to-pv1")))#-0.17
    push!(pv1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=pv1, g = -0.18, name = string(c.name, "-","pv1-to-pv1")))
    
       
   
        push!(vip1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=vip1, g = 0.47, depression = depression, name = string(c.name, "-","integrator1-to-vip1")))
    push!(pv1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=pv1, g = 0.31, depression = depression, name = string(c.name, "-","integrator1-to-pv1")))
    
		
		push!(sst1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=sst1, g = 0.22, facilitation = facilitation, name = string(c.name, "-","integrator1-to-sst1")))
    push!(integrator1.list_syn, nmda_syn(neuron_pre = E1, neuron_post = integrator1, g = 0.15, name = string(c.name, "-","ecell1-to-integrator1")))

    push!(c.list_dend,dend1)
  
    push!(c.list_soma,E1)
   
    push!(c.list_vip,vip1)
  
    push!(c.list_sst,sst1)
  
    push!(c.list_pv,pv1)

    push!(c.list_integrator,integrator1)



end

#TODO add a list of microcircuit as argument to the dynamics functions

function construct_two_local_microcircuit_integrator(; dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation = false, adaptation = false, depression = false, td_to_vip = true)
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name


    c1 = microcircuit(name="microcircuit1")
    c2 = microcircuit(name = "microcircuit2")
    
    construct_one_local_microcircuit_integrator!(c1; dend_param, facilitation, adaptation, td_to_vip, depression)
    construct_one_local_microcircuit_integrator!(c2; dend_param, facilitation, adaptation, td_to_vip, depression)


    push!(c1.list_sst[1].list_syn, nmda_syn(neuron_pre=c2.list_soma[1],neuron_post=c1.list_sst[1], g = 0.0435, name = string("ecell2-to-sst1")))
    push!(c2.list_sst[1].list_syn, nmda_syn(neuron_pre=c1.list_soma[1],neuron_post=c2.list_sst[1], g = 0.0435, name = string("ecell1-to-sst2")))


    list_microcircuit = [c1,c2]

    return list_microcircuit
end

	
	
end

# ╔═╡ 9c147c04-fc96-4873-97eb-1bf65839ee93
begin
lc = construct_two_local_microcircuit_integrator()
sim = simulation_parameters()
t_tot = Int64(sim.Tfin/sim.dt)
	
end

# ╔═╡ 570c0cd7-8469-448a-a6eb-5accefed980e
temp_dict=Dict(
"microcircuit1-vipcell1" =>zeros(t_tot), 
"microcircuit1-sstcell1" =>zeros(t_tot),
"microcircuit1-pvcell1" =>zeros(t_tot),
"microcircuit1-ecell1" =>zeros(t_tot),
"microcircuit1-integrator1" =>zeros(t_tot),
"microcircuit2-vipcell1" =>zeros(t_tot),
"microcircuit2-sstcell1" =>zeros(t_tot),
"microcircuit2-pvcell1" =>zeros(t_tot),
"microcircuit2-ecell1" =>zeros(t_tot),
"microcircuit2-integrator1" =>zeros(t_tot)
)

# ╔═╡ 4f7a9670-077c-488c-be66-264db4616673
begin
	
	function create_stimu(stim_strength::Float64,sim::simulation_parameters; Tsim = 0.500, Tinterval=1) #in s
		
		Tstop = sim.Tfin
		dt = sim.dt
		
		
		
		
		dt_pattern = Int((Tsim+Tinterval)/dt)
		temp_stim = vcat(0.0.*ones(Int(Tinterval/dt)), stim_strength.*ones(Int(Tsim/dt)))
		dummy = 0
		stim = copy(temp_stim)
		
		while dummy < Int(Tstop/dt)
		
			stim = vcat(stim,temp_stim)
			dummy += dt_pattern
			
		end
		print(dummy)
		# only take the t_tot number of elements
		return stim[1:Int(Tstop/dt)]
		
	end
	
	
	
	
	
end

# ╔═╡ 368c79d3-707c-4847-a7f0-a87758447d97
sim

# ╔═╡ 970c91ff-ad74-4d1a-9f5a-c62bf0a3dd27
begin
	
	temp_dict["microcircuit1-ecell1"] = create_stimu(0.1,sim)
	temp_dict["microcircuit2-ecell1"] = create_stimu(0.1,sim)


end

# ╔═╡ c939e024-db53-4a4b-a756-0739e54c4a23
full_time_dynamics(lc,sim,temp_dict)

# ╔═╡ f0ae35a9-5e24-4f71-801f-6ec945bcad5a
begin
	plot(lc[1].list_soma[1].r)
	plot!(lc[2].list_soma[1].r)
	
end

# ╔═╡ 1afa4083-8c7f-4512-8bea-c8b95fd5f015
begin
	plot(lc[1].list_sst[1].r)
	plot!(lc[2].list_sst[1].r)
	
end

# ╔═╡ ca8f5903-655e-4f3a-b552-baebd978fba9
lc[1].list_sst[1].list_syn[1]

# ╔═╡ 7fb6cc17-7c63-4b5e-b352-b7111982db4b
begin
	plot(lc[1].list_integrator[1].r)
	plot!(lc[2].list_integrator[1].r)
	
end

# ╔═╡ Cell order:
# ╠═a2d39c60-d06a-11eb-1364-c1e29d8f626f
# ╠═1334881e-d06b-11eb-3ec8-13f9d3037c39
# ╠═e5810c90-d06b-11eb-26d3-e378a54e0bf7
# ╠═16aa4df0-d06b-11eb-3d89-6fafd70abb11
# ╠═9c147c04-fc96-4873-97eb-1bf65839ee93
# ╠═570c0cd7-8469-448a-a6eb-5accefed980e
# ╠═4f7a9670-077c-488c-be66-264db4616673
# ╠═368c79d3-707c-4847-a7f0-a87758447d97
# ╠═970c91ff-ad74-4d1a-9f5a-c62bf0a3dd27
# ╠═c939e024-db53-4a4b-a756-0739e54c4a23
# ╠═f0ae35a9-5e24-4f71-801f-6ec945bcad5a
# ╠═1afa4083-8c7f-4512-8bea-c8b95fd5f015
# ╠═ca8f5903-655e-4f3a-b552-baebd978fba9
# ╠═7fb6cc17-7c63-4b5e-b352-b7111982db4b
