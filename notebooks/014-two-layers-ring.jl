### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 75f6a42e-5832-11ec-3134-6398493813e0
begin
using DrWatson
	
quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
push!(LOAD_PATH,srcdir())

end

# ╔═╡ 9fa0ceaa-845f-47a2-8089-f7e134c813a8
begin

using Revise
using MyNeurosciencePackage
using PlutoUI
using CairoMakie
using Parameters, JLD2, Statistics

end

# ╔═╡ 6d814b54-e441-4ca5-a9ad-ec076a792c2b
md"""
This notebook investigate the construction of a two layers neural network and its behavior when presented to a somatic input.

Two things to note:
- The upper layer is a ring constructed with integrator networks
- Parameters for the ring are from Engel et al 2015



"""

# ╔═╡ 219906fa-edc7-46f3-b6b6-723604fa899e
begin 
sim = simulation_parameters()
time_tot = 10500

stim = 50.0
stimlist = vcat(zeros(500),ones(1000), zeros(1000), ones(1000), zeros(1000), ones(1000), zeros(1000), ones(1000), zeros(1000),ones(1000), zeros(1000))


end

# ╔═╡ 5886bcce-ab44-4380-a747-ad7f7907765e
md"""
Parameters of the network

"""

# ╔═╡ 87c58b84-d9aa-4790-b415-0168f5f0d47b
begin 

param_syn_sim = parameters_syn_strength_microcircuit(
	nmda_soma_to_int = 0.01, 
	nmda_soma_to_pv = 0.17, 
	nmda_soma_to_vip = 0.35, 
	gaba_vip_to_sst = -0.05, 
	gaba_ngfc_to_dend = -0.15, 
	gaba_sst_to_ngfc = -0.4, 
	nmda_soma_to_sst =  0.0435, 
	nmda_int_to_dend = 0.05*0.4)

param_micro_sim = parameters_microcircuit(
	time_tot = time_tot, 
	soma_to_pv_depression = false, 
	sst_to_vip_facilitation = true,
	int_to_vip_depression = true, 
	vip_to_sst_facilitation = true,
	soma_to_sst_facilitation = true, 
	ngfc_to_dend_depression = true, 
	int_to_dend_depression = true,
	noise =true)
end

# ╔═╡ 98bb0372-4773-49ec-bc91-5058c67a58ee
begin

	
ring_simu = Dict{String,Vector{Float64}}()

nn = create_two_layers(param_syn = param_syn_sim,param_micro = param_micro_sim)



end

# ╔═╡ e00e68b4-03c5-4ad9-aa9e-0153446d0571
begin

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

end

# ╔═╡ db5dbd35-e226-4780-af76-9a1d98629b96
full_time_dynamics!(nn, sim)

# ╔═╡ 2cac8391-1d96-464a-9af2-d5a54171d60f
bump_animation(nn,5,projectdir("notebooks","animation\\"),"test_ring_connectivity",".mp4")

# ╔═╡ 7872700d-d360-4e2c-9a0f-29ca94eb4452
LocalResource(projectdir("notebooks","animation","test_ring_connectivity.mp4"))


# ╔═╡ cb02bbf4-3ec7-468a-9b71-9ab2af4e11cf
begin

	fig = Figure()
	ax1 = Axis(fig[1,1])

	lines!(ax1,1:1:64,get_firing_rate(nn.ring_integrator,1400.0))
	lines!(ax1,1:1:64,get_firing_rate(nn.ring_integrator,3400.0))
	lines!(ax1,1:1:64,get_firing_rate(nn.ring_integrator,5400.0))
	lines!(ax1,1:1:64,get_firing_rate(nn.ring_integrator,7400.0))
	lines!(ax1,1:1:64,get_firing_rate(nn.ring_integrator,9400.0))

		fig
end

# ╔═╡ dc42957d-8279-443d-8355-211b73322fd9


# ╔═╡ ba3ae29a-6f06-444c-ba94-e656537df2b6


# ╔═╡ Cell order:
# ╟─75f6a42e-5832-11ec-3134-6398493813e0
# ╟─9fa0ceaa-845f-47a2-8089-f7e134c813a8
# ╟─6d814b54-e441-4ca5-a9ad-ec076a792c2b
# ╠═219906fa-edc7-46f3-b6b6-723604fa899e
# ╠═5886bcce-ab44-4380-a747-ad7f7907765e
# ╠═87c58b84-d9aa-4790-b415-0168f5f0d47b
# ╠═98bb0372-4773-49ec-bc91-5058c67a58ee
# ╠═e00e68b4-03c5-4ad9-aa9e-0153446d0571
# ╠═db5dbd35-e226-4780-af76-9a1d98629b96
# ╠═2cac8391-1d96-464a-9af2-d5a54171d60f
# ╠═7872700d-d360-4e2c-9a0f-29ca94eb4452
# ╠═cb02bbf4-3ec7-468a-9b71-9ab2af4e11cf
# ╠═dc42957d-8279-443d-8355-211b73322fd9
# ╠═ba3ae29a-6f06-444c-ba94-e656537df2b6
