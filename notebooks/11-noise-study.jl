### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ ef519000-1f99-11ec-2c53-7db5988d8d93
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())
end

# ╔═╡ 2f9182be-b36a-420f-b34e-adaeeb23c47f
begin
	using Revise
	using MyNeurosciencePackage
	using CairoMakie
	using Parameters, JLD2, Statistics	
	using DataFrames, PlutoUI
	#using LsqFit
	using ColorSchemes, StaticArrays
	
	set_theme!(theme_light())

	
end

# ╔═╡ 3fd4fa19-28dc-4372-8a89-f8d5527f47d0
params = Dict(
	"number_repetitions" => 6.0,
    "τ_integrator" => 2.0,
	"value_stim_f1" => 0.2,
	"Tinter" => 1.0,#
	"Tstim" => 0.5,
	"Tfin" => 1.0,
	"value_stim_f2" => 0.0,    
	"initial_time" => 4.0
)

# ╔═╡ f51ae86a-a170-4cbc-923e-661a327c6d52
begin
	function run_simu1_script(params; sst_adaptation = true, soma_adaptation = true, pv_to_soma_depression = true, soma_to_vip_facilitation = true, sst_to_vip_facilitation = true, soma_to_sst_facilitation = true, vip_to_sst_facilitation = true, soma_to_pv_depression = true, int_to_vip_depression = true, int_to_pv_depression = true, int_to_sst_facilitation = true, int_to_dend_depression = true,integrator_tc=0.8, time_tot=1000,cross_int_to_vip_depression = true, cross_int_to_pv_depression = true, cross_int_to_sst_facilitation = true, cross_int_to_dend_depression = true ,cross_soma_to_sst_facilitation = true)
		
		stim_1, stim_2 = create_deterministic_oddball(params)
	
	@unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

	t_tot = length(stim_1)
	
	lc = construct_two_local_microcircuit_integrator_full_param(dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0),sst_adaptation=sst_adaptation, soma_adaptation=soma_adaptation, pv_to_soma_depression=pv_to_soma_depression, soma_to_vip_facilitation=soma_to_vip_facilitation, sst_to_vip_facilitation=sst_to_vip_facilitation, soma_to_sst_facilitation=soma_to_sst_facilitation, vip_to_sst_facilitation=vip_to_sst_facilitation, soma_to_pv_depression=soma_to_pv_depression, int_to_vip_depression=int_to_vip_depression, int_to_pv_depression=int_to_pv_depression, int_to_sst_facilitation=int_to_sst_facilitation, int_to_dend_depression=int_to_dend_depression,integrator_tc = τ_integrator, time_tot=t_tot,cross_int_to_vip_depression=cross_int_to_vip_depression, cross_int_to_pv_depression=cross_int_to_pv_depression, cross_int_to_sst_facilitation=cross_int_to_sst_facilitation, cross_int_to_dend_depression=cross_int_to_dend_depression ,cross_soma_to_sst_facilitation=cross_soma_to_sst_facilitation)
            
            
    sim = simulation_parameters(Tstimduration=Tstim, TISI=Tinter, Tfin=Tfin, Tinit=2000)
        
    
            
            # For now dictionnary is define manually
    oddball = Dict(
    "microcircuit1-vipcell1" => zeros(t_tot), 
    "microcircuit1-sstcell1" => zeros(t_tot),
    "microcircuit1-pvcell1" => zeros(t_tot),
    "microcircuit1-ecell1" => zeros(t_tot),
    "microcircuit1-integrator1" => zeros(t_tot),
    "microcircuit2-vipcell1" => zeros(t_tot),
    "microcircuit2-sstcell1" => zeros(t_tot),
    "microcircuit2-pvcell1" => zeros(t_tot),
    "microcircuit2-ecell1" => zeros(t_tot),
    "microcircuit2-integrator1" => zeros(t_tot)
    )
	
	oddball["microcircuit1-ecell1"] = stim_1
    oddball["microcircuit2-ecell1"] = stim_2
        
	
	full_time_dynamics(lc, sim, oddball)
	return lc,oddball,sim
	end
end

# ╔═╡ 7274e5c9-90ec-4a86-a293-640a44fdb1f3
begin
	lc,oddball,sim = run_simu1_script(params)
fig = plot_local_circuit(lc, sim, oddball, params)
end

# ╔═╡ 7199c541-e95e-4767-a3b9-a148c85f47ec
compute_MMN_oddball(lc,params)

# ╔═╡ f5e72495-b68c-4e2d-93c4-80cdba6b1ebb
md"""
Average MMN on a noisy network
"""

# ╔═╡ 423a0eba-6f56-4407-a090-a4530019cd45
begin
	function compute_average_MMN(tot_simu::Int64, params)
		list_MMN = Float64[]
		for i=1:tot_simu
		lc,oddball,sim = run_simu1_script(params)
		push!(list_MMN, compute_MMN_oddball(lc,params))
		end
		return list_MMN
	end
	
	list_MMN = compute_average_MMN(1000,params)
	
	
	
end

# ╔═╡ ee7a807f-a0ab-41dc-b524-7a00bfe90339
density(list_MMN)

# ╔═╡ 45e0232d-2ea8-49dd-8446-5e22c99931b1
mean(list_MMN)

# ╔═╡ Cell order:
# ╠═ef519000-1f99-11ec-2c53-7db5988d8d93
# ╠═2f9182be-b36a-420f-b34e-adaeeb23c47f
# ╠═3fd4fa19-28dc-4372-8a89-f8d5527f47d0
# ╠═f51ae86a-a170-4cbc-923e-661a327c6d52
# ╠═7274e5c9-90ec-4a86-a293-640a44fdb1f3
# ╠═7199c541-e95e-4767-a3b9-a148c85f47ec
# ╠═f5e72495-b68c-4e2d-93c4-80cdba6b1ebb
# ╠═423a0eba-6f56-4407-a090-a4530019cd45
# ╠═ee7a807f-a0ab-41dc-b524-7a00bfe90339
# ╠═45e0232d-2ea8-49dd-8446-5e22c99931b1
