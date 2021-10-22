### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ a7532200-1b11-11ec-139c-c5fe3272aee2
begin
	using DrWatson
	
	#quickactivate("/mnt/c/Users/kevin/OneDrive/1-Research_projects/1-project-categorical-MMN","1-project-categorical-MMN")]
	@quickactivate "1-project-categorical-MMN"
	push!(LOAD_PATH,srcdir())
end

# ╔═╡ b0aa5d53-00ff-4b8c-a8e5-905eb25a20e8
begin
	using Revise
	using MyNeurosciencePackage
	using CairoMakie
	using Parameters, Statistics	
	using DataFrames, PlutoUI
	#using LsqFit
	using ColorSchemes, StaticArrays
	
	set_theme!(theme_light())

	
end

# ╔═╡ 65247c10-4c25-47a4-bd22-92ced66b4182
begin
	using JLD, JLD2
	# tot_MMN_wt_int_to_sst = generate_totMMN(int_to_sst_facilitation = false, noise = false)
	# save(datadir("temp\\script10-temp1.jld"),"tot_MMN_wt_int_to_sst ",tot_MMN_wt_int_to_sst )

# 	 tot_MMN_wt_int_to_vip = generate_totMMN(noise = false, int_to_vip_depression = false)
# 	 save(datadir("temp\\script10-temp2.jld"),"tot_MMN_wt_int_to_vip",tot_MMN_wt_int_to_vip )
	
# 		 tot_MMN_wt_td = generate_totMMN(sst_adaptation = false, noise=false, int_to_sst_facilitation = false, int_to_dend_depression = false, int_to_vip_depression = false, cross_int_to_dend_depression = false)
# 	 save(datadir("temp\\script10-temp3.jld"),"tot_MMN_wt_td",tot_MMN_wt_td )
	
	
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
    "τ_integrator" => 4.0,
	"value_stim_f1" => 0.2,
	"Tinter" => 1.0,#
	"Tstim" => 0.5,
	"Tfin" => 1.0,
	"value_stim_f2" => 0.0,    
	"initial_time" => 4.0
)

# ╔═╡ a37d8a1c-1ee3-4206-9c9d-d2032c8fe1bf

begin
	function run_simu1_script(params;param_microcircuit = parameters_microcircuit(), param_syn_microcircuit = parameters_syn_strength_microcircuit(), param_inter_microcircuit = parameters_inter_microcircuit())
		
		stim_1, stim_2 = create_deterministic_oddball(params)
	
	@unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

	t_tot = length(stim_1)
	param_microcircuit.time_tot = t_tot
	lc = construct_two_local_microcircuit_integrator(param_microcircuit = param_microcircuit, param_syn_microcircuit = param_syn_microcircuit, param_inter_microcircuit = param_inter_microcircuit)
            
    
            
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
    "microcircuit2-integrator1" => zeros(t_tot),
	"microcircuit1-dend1" => zeros(t_tot),
	"microcircuit2-dend1" => zeros(t_tot)
    )
	
	oddball["microcircuit1-ecell1"] = stim_1
    oddball["microcircuit2-ecell1"] = stim_2
        

    sim = simulation_parameters(Tstimduration=Tstim, TISI=Tinter, Tfin=Tfin, Tinit=2000,current = oddball)
	full_time_dynamics!(lc, sim)
	return lc,oddball,sim
	end
end

# ╔═╡ 37518325-f6e2-497a-9072-972cdc1c4c34
begin
	param_microcircuit = parameters_microcircuit(noise = false,top_down_to_interneurons = [0.54,0.36,0.12])
	lc,oddball,sim = run_simu1_script(params;param_microcircuit = param_microcircuit)
fig_defualt = plot_local_circuit(lc, sim, oddball, params)
end

# ╔═╡ b40a7097-53e9-4e87-a6c1-af647810989a
begin
	
	param_microcircuit_test = parameters_microcircuit(noise = false,top_down_to_interneurons = [0.54,0.36,0.12], sst_to_vip_facilitation = false,vip_to_sst_facilitation = false)
	
	syn_strength = parameters_syn_strength_microcircuit(gaba_vip_to_sst = -0.17)
	
	
	lc_test,oddball_test,sim_test = run_simu1_script(params;param_microcircuit = param_microcircuit_test, param_syn_microcircuit = syn_strength)
	
fig_test = plot_local_circuit(lc_test, sim_test, oddball_test, params)
end

# ╔═╡ 9e5adf20-3ad5-4cd4-8283-91918ed2eb5e
md""" ## No depression on int to dend


"""

# ╔═╡ 5ae76c69-d9a3-46af-a513-84943be17e13
begin
	lccid,oddballcid,simcid = run_simu1_script(params, cross_int_to_dend_depression = false,noise = false)
	figcid = plot_local_circuit(lccid, simcid, oddballcid, params)
end

# ╔═╡ e8e16bac-720f-4696-a0a9-aaaa1865b181


# ╔═╡ dc1f33e4-fdf6-4010-b6a1-87446d3a0080
begin
	lciv,oddballiv,simiv = run_simu1_script(params, int_to_vip_depression = false,noise=false)
	figiv = plot_local_circuit(lciv, simiv, oddballiv, params)
end

# ╔═╡ ba55401f-c545-4841-ade7-fe7851468038


# ╔═╡ b1a7eb26-4118-4e9e-b864-4347b3ede0c2
begin
	lcif,oddballif,simif = run_simu1_script(params, sst_adaptation = false, noise=false, int_to_sst_facilitation = false, int_to_dend_depression = false, int_to_vip_depression = false, cross_int_to_dend_depression = false)
	figif = plot_local_circuit(lcif, simif, oddballif, params)
end

# ╔═╡ 260e2503-ba5a-42a2-927f-d5969b107268
begin
	lcif2,oddballif2,simif2 = run_simu1_script(params, noise=false, soma_to_sst_facilitation = false, soma_to_vip_facilitation = false, soma_to_pv_depression = false)
	figif2 = plot_local_circuit(lcif2, simif2, oddballif2, params)
	# TODO full figure
end

# ╔═╡ 0ef65380-d7ab-471c-a56b-2c1db1d3fa40


# ╔═╡ 652c1ee0-fc3f-4bdf-aa7f-4e1c52c7b550
md""" # Figures for Lab Meeting


"""

# ╔═╡ be78b3fd-c711-4ee3-b71e-5c4d53e74ce7
begin
	# TODO: do the MMN effect as a function of ISI for example (and \tau int)
	# Do i want the number of repetitions as a parameters? Maybe just one, 3 and 6?
	# Want a function that compute MMN effect?
	
	
	
	function generate_totMMN(;int_to_sst_facilitation = true, noise = false, int_to_vip_depression = true,sst_adaptation = false,int_to_dend_depression = true,  cross_int_to_dend_depression = true)
		
			num_list = collect(1.0:2.0:11.0)
τ_list = collect(0.1:0.4:5.0)
Tinter_list = collect(0.5:0.5:10.0)
	
tot_MMN = zeros(length(num_list),length(τ_list), length(Tinter_list))
	for (i,num) in enumerate(num_list)
		for (k,τ) in enumerate(τ_list)
			for (j,T) in enumerate(Tinter_list)
			
				params_list = Dict(
	"number_repetitions" => num,
    "τ_integrator" => τ,
	"value_stim_f1" => 0.2,
	"Tinter" => T,#
	"Tstim" => 0.5,
	"Tfin" => 1.0,
	"value_stim_f2" => 0.0,    
	"initial_time" => 4.0
)
					
		lc,oddball,sim = run_simu1_script(params_list;int_to_sst_facilitation = int_to_sst_facilitation, noise = noise, int_to_vip_depression = int_to_vip_depression,sst_adaptation = sst_adaptation,int_to_dend_depression = int_to_dend_depression,  cross_int_to_dend_depression = cross_int_to_dend_depression)
		tot_MMN[i,k,j] =compute_MMN_oddball(lc,params_list)
	end
				
		end
		
	end
	
	return tot_MMN
	end
	
end

# ╔═╡ 96547ab8-642b-4c83-a854-976ba9863197
begin
	
	dic_tot = load(datadir("temp/script10-temp3.jld"))
	load_tot_MMN = dic_tot["tot_MMN_wt_td"]
	# dic_tot = load(datadir("temp\\script10-temp2.jld"))
	# load_tot_MMN = dic_tot["tot_MMN_wt_int_to_vip"]
	num_list = collect(1.0:2.0:11.0)
τ_list = collect(0.1:0.4:5.0)
Tinter_list = collect(0.5:0.5:10.0)
end

# ╔═╡ 15496953-ac1f-4ea2-b19a-6a90b369355b
begin
	# in this oddball context, the exact time of computing the effect is known
	
	fig_temp = Figure()
	index_τ = 3
	τ_value = τ_list[index_τ]
	ax_temp = Axis(fig_temp[1,1],title="MMN vs ISI - τ = $τ_value s - Network without TD facilitation/depression",xlabel="ISI (s)",ylabel="Δ firing rate (Hz)")
	for (i,n) in enumerate(num_list)

		lines!(ax_temp,Tinter_list,load_tot_MMN[i,index_τ,:],linewidth=3,label="Nbr rep = $n")
	end
	fig_temp
	fig_temp[1, 2] = Legend(fig_temp, ax_temp, "Legend", framevisible = false)
	fig_temp
	save(plotsdir("notebook10-MMN_wt_td.png"),fig_temp)
end

# ╔═╡ 139c02d6-df65-4e88-bf80-55c304d82351


# ╔═╡ af4ea3f3-a9df-4e0d-8e0a-c369ced9c502
begin
	lcif3,oddballif3,simif3 = run_simu1_script(params, noise=false, int_to_sst_connection = false, int_to_pv_depression = true, int_to_vip_depression = true, int_to_dend_depression = true, soma_to_sst_facilitation = true)
	figif3 = plot_local_circuit(lcif3, simif3, oddballif3, params)
	# TODO full figure

	# influence of the strength of sst

	fig = Figure()
	colors = to_colormap(:viridis, length(collect(0.1:0.4:5.0)))

	ax = Axis(fig[1,1])
	jnum=0
	influence_sst = zeros(length(collect(0.1:0.4:5.0)), length(collect(0.5:0.5:7.0)), length(collect(0.0:0.1:1.0)))
	knum = 0
	for k=0.5:0.5:7.0
		knum+=1
		jnum=0
	for j=0.1:0.4:5.0
		jnum +=1
		inum =0
		for i=0.0:0.1:1.0
			inum +=1
			params_list = Dict(
    	"number_repetitions" => 7.0,#collect(1.0:1.0:10.0),
        "τ_integrator" => j,#collect(0.1:0.4:5.0),
    	"value_stim_f1" => 0.2,
    	"Tinter" => k,#1.0,#collect(0.5:0.5:7.0),#
    	"Tstim" => 0.5,
    	"Tfin" => 1.0,
    	"value_stim_f2" => 0.0,    
    	"initial_time" => 5.0,
        "sst_strength" => i#collect(0.0:0.1:1.0)
    )
	tempname = savename(params_list)
	temp = load(datadir("sims","script4","oddball_task_$tempname.jld2"))

	influence_sst[jnum,knum,inum] = temp["MMN effect"]
	end
	#scatter!(ax,collect(0.0:0.1:1.0),influence_sst, label="tau = $j", color=colors[jnum])
end
end
#fig[1, 2] = Legend(fig, ax, "Legend", framevisible = false)
perc = collect(0.0:0.1:1.0).*0.22
ax_hm, hm_rep = heatmap(ax,collect(0.5:0.5:7.0),perc,influence_sst[3,:,:])
Colorbar(fig[1, 2],hm_rep,label="MMN effect (Hz)")
#SupTitle_rep = Label(heatmap_fixed_rep[0,1:2], "Heatmap of MMN - Number of repetitions = $value_heatmap_rep ")
ax_hm.xlabel= "τ integrator"
ax_hm.ylabel = "Tinter (s)"
fig


fig_2d = Figure()
ax2d = Axis(fig_2d[1,1],xlabel = "Percentage of int to sst", ylabel = "MMn effect ( HZ)", title="ISI = 1s")
colors = to_colormap(:viridis, length(collect(0.1:0.8:5.0)))
for jnum=1:length(collect(0.1:0.8:5.0))
	j = collect(0.1:0.8:5.0)[jnum]
scatter!(ax2d,perc,influence_sst[jnum,2,:], label="tau = $j (s)", color=colors[jnum], markersize=15)
end
fig_2d[1, 2] = Legend(fig_2d, ax2d, "Legend", framevisible = false)

fig_2d





end

# ╔═╡ 12bac34b-0af1-49e2-b43d-6b5921099898
begin
	
	#ax_fixed_τ = Axis(heatmap_fixed_τ[1,1])
	index_heatmap = 3
	value_heatmap = τ_list[index_heatmap]
	heatmap_fixed_τ = Figure()
	
	ax_hm,hm = heatmap(heatmap_fixed_τ[1,1],num_list,Tinter_list, load_tot_MMN[:,index_heatmap,:])
	Colorbar(heatmap_fixed_τ[1, 2],hm,label="MMN effect (Hz)")
	SupTitle = Label(heatmap_fixed_τ[0,1:2], "Heatmap of MMN - τ integrator = $value_heatmap s")
	ax_hm.xlabel = "Number of repetitions"
	ax_hm.ylabel = "Tinter (s)"
	heatmap_fixed_τ
	
	save(plotsdir("notebook9-heatmap_τ=$value_heatmap.png"),heatmap_fixed_τ)
	heatmap_fixed_τ
end

# ╔═╡ Cell order:
# ╟─a7532200-1b11-11ec-139c-c5fe3272aee2
# ╟─b0aa5d53-00ff-4b8c-a8e5-905eb25a20e8
# ╠═ac0516db-654d-48d9-8881-18990b915be3
# ╠═69efc5c1-3279-4c53-ba7b-8a5ce25265bc
# ╠═a37d8a1c-1ee3-4206-9c9d-d2032c8fe1bf
# ╠═37518325-f6e2-497a-9072-972cdc1c4c34
# ╠═b40a7097-53e9-4e87-a6c1-af647810989a
# ╠═9e5adf20-3ad5-4cd4-8283-91918ed2eb5e
# ╠═5ae76c69-d9a3-46af-a513-84943be17e13
# ╠═e8e16bac-720f-4696-a0a9-aaaa1865b181
# ╠═dc1f33e4-fdf6-4010-b6a1-87446d3a0080
# ╠═ba55401f-c545-4841-ade7-fe7851468038
# ╠═b1a7eb26-4118-4e9e-b864-4347b3ede0c2
# ╠═260e2503-ba5a-42a2-927f-d5969b107268
# ╠═af4ea3f3-a9df-4e0d-8e0a-c369ced9c502
# ╠═0ef65380-d7ab-471c-a56b-2c1db1d3fa40
# ╠═652c1ee0-fc3f-4bdf-aa7f-4e1c52c7b550
# ╠═be78b3fd-c711-4ee3-b71e-5c4d53e74ce7
# ╠═65247c10-4c25-47a4-bd22-92ced66b4182
# ╠═96547ab8-642b-4c83-a854-976ba9863197
# ╠═15496953-ac1f-4ea2-b19a-6a90b369355b
# ╠═12bac34b-0af1-49e2-b43d-6b5921099898
# ╠═139c02d6-df65-4e88-bf80-55c304d82351
