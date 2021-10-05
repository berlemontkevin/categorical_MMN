### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ a7532200-1b11-11ec-139c-c5fe3272aee2
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())
end

# ╔═╡ b0aa5d53-00ff-4b8c-a8e5-905eb25a20e8
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

# ╔═╡ ac0516db-654d-48d9-8881-18990b915be3
begin
	# Plots variables
	noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")

end

# ╔═╡ 69efc5c1-3279-4c53-ba7b-8a5ce25265bc

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

# ╔═╡ a37d8a1c-1ee3-4206-9c9d-d2032c8fe1bf
begin
	function run_simu1_script(params; sst_adaptation = true, soma_adaptation = true, pv_to_soma_depression = true, soma_to_vip_facilitation = true, sst_to_vip_facilitation = true, soma_to_sst_facilitation = true, vip_to_sst_facilitation = true, soma_to_pv_depression = true, int_to_vip_depression = true, int_to_pv_depression = true, int_to_sst_facilitation = true, int_to_dend_depression = true,integrator_tc=0.8, time_tot=1000,cross_int_to_vip_depression = true, cross_int_to_pv_depression = true, cross_int_to_sst_facilitation = true, cross_int_to_dend_depression = true ,cross_soma_to_sst_facilitation = true, noise=true)
		
		stim_1, stim_2 = create_deterministic_oddball(params)
	
	@unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

	t_tot = length(stim_1)
	
	lc = construct_two_local_microcircuit_integrator_full_param(dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0),sst_adaptation=sst_adaptation, soma_adaptation=soma_adaptation, pv_to_soma_depression=pv_to_soma_depression, soma_to_vip_facilitation=soma_to_vip_facilitation, sst_to_vip_facilitation=sst_to_vip_facilitation, soma_to_sst_facilitation=soma_to_sst_facilitation, vip_to_sst_facilitation=vip_to_sst_facilitation, soma_to_pv_depression=soma_to_pv_depression, int_to_vip_depression=int_to_vip_depression, int_to_pv_depression=int_to_pv_depression, int_to_sst_facilitation=int_to_sst_facilitation, int_to_dend_depression=int_to_dend_depression,integrator_tc = τ_integrator, time_tot=t_tot,cross_int_to_vip_depression=cross_int_to_vip_depression, cross_int_to_pv_depression=cross_int_to_pv_depression, cross_int_to_sst_facilitation=cross_int_to_sst_facilitation, cross_int_to_dend_depression=cross_int_to_dend_depression ,cross_soma_to_sst_facilitation=cross_soma_to_sst_facilitation,noise = noise)
            
            
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

# ╔═╡ 37518325-f6e2-497a-9072-972cdc1c4c34
begin
	lc,oddball,sim = run_simu1_script(params; noise=false)
fig = plot_local_circuit(lc, sim, oddball, params)
end

# ╔═╡ 9e5adf20-3ad5-4cd4-8283-91918ed2eb5e
md""" ## No depression on int to dend


"""

# ╔═╡ 5ae76c69-d9a3-46af-a513-84943be17e13
begin
	lccid,oddballcid,simcid = run_simu1_script(params, cross_int_to_dend_depression = false)
	figcid = plot_local_circuit(lccid, simcid, oddballcid, params)
end

# ╔═╡ e8e16bac-720f-4696-a0a9-aaaa1865b181


# ╔═╡ dc1f33e4-fdf6-4010-b6a1-87446d3a0080
begin
	lciv,oddballiv,simiv = run_simu1_script(params, int_to_vip_depression = false)
	figiv = plot_local_circuit(lciv, simiv, oddballiv, params)
end

# ╔═╡ ba55401f-c545-4841-ade7-fe7851468038


# ╔═╡ b1a7eb26-4118-4e9e-b864-4347b3ede0c2
begin
	lcif,oddballif,simif = run_simu1_script(params, int_to_sst_facilitation = false, int_to_vip_depression = false)
	figif = plot_local_circuit(lcif, simif, oddballif, params)
end

# ╔═╡ 652c1ee0-fc3f-4bdf-aa7f-4e1c52c7b550


# ╔═╡ 65247c10-4c25-47a4-bd22-92ced66b4182


# ╔═╡ Cell order:
# ╠═a7532200-1b11-11ec-139c-c5fe3272aee2
# ╠═b0aa5d53-00ff-4b8c-a8e5-905eb25a20e8
# ╠═ac0516db-654d-48d9-8881-18990b915be3
# ╠═69efc5c1-3279-4c53-ba7b-8a5ce25265bc
# ╟─a37d8a1c-1ee3-4206-9c9d-d2032c8fe1bf
# ╠═37518325-f6e2-497a-9072-972cdc1c4c34
# ╠═9e5adf20-3ad5-4cd4-8283-91918ed2eb5e
# ╠═5ae76c69-d9a3-46af-a513-84943be17e13
# ╠═e8e16bac-720f-4696-a0a9-aaaa1865b181
# ╠═dc1f33e4-fdf6-4010-b6a1-87446d3a0080
# ╠═ba55401f-c545-4841-ade7-fe7851468038
# ╠═b1a7eb26-4118-4e9e-b864-4347b3ede0c2
# ╠═652c1ee0-fc3f-4bdf-aa7f-4e1c52c7b550
# ╠═65247c10-4c25-47a4-bd22-92ced66b4182
