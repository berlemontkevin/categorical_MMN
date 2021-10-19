### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 517c9c3e-0cce-11ec-0a91-07fda58ce84d
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())
end

# ╔═╡ 9fd2c4b4-975b-45b4-9901-5c22180ce984
begin
	using Revise
	using MyNeurosciencePackage
	using CairoMakie
	using Parameters, JLD2, Statistics	
	using DataFrames, PlutoUI
	using LsqFit
	using ColorSchemes, StaticArrays
	
	set_theme!(theme_light())

	
end

# ╔═╡ e9d59c50-5173-4674-bcfb-97dcde7b7035
begin
	using JLD
	#save(datadir("temp\\script9-temp1.jld"),"tot_MMN",tot_MMN)

	
end

# ╔═╡ 859449c2-238e-430c-aebc-0d226314df0b
begin
	# Plots variables
	noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")

end

# ╔═╡ eab23a48-c873-408b-91da-5c208a82f117

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

# ╔═╡ 70649557-de48-424d-ada1-12e6c5540417
begin
	function run_simu1_script(params)
		stim_1, stim_2 = create_deterministic_oddball(params)
	
	@unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

	t_tot = length(stim_1)
	
	lc = construct_two_local_microcircuit_integrator(facilitation=true, depression=true, adaptation=true,  τ_integrator = τ_integrator ,time_tot = t_tot)
            
            
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

# ╔═╡ b3bbd648-8b82-414b-ae39-5cf70f70f06a
begin
	lc,oddball,sim = run_simu1_script(params)
	fig = plot_local_circuit(lc, sim, oddball, params)
	save(plotsdir("notebook9-oddballparadigm.png"), fig) # output size = 800 x 600 pixels
fig
end

# ╔═╡ fd27c550-f629-48ef-b804-71319d9532ec
lc[1].list_dend[1].list_syn_post_nmda

# ╔═╡ 6a0d37bf-f458-4273-909d-4f06e440b0fe
fig_syn = plot_local_circuit_synapses(lc, sim, oddball, params)

# ╔═╡ 707cae76-fdf3-4d6b-856b-b57a65e97ec8
begin
	num_list = collect(1.0:2.0:11.0)
τ_list = collect(0.1:0.4:5.0)
Tinter_list = collect(0.5:0.5:10.0)
	
tot_MMN = zeros(length(num_list),length(τ_list), length(Tinter_list))
end

# ╔═╡ 2bc2f601-a36d-498c-b208-a5099c9e2d59
begin
	# TODO: do the MMN effect as a function of ISI for example (and \tau int)
	# Do i want the number of repetitions as a parameters? Maybe just one, 3 and 6?
	# Want a function that compute MMN effect?
	
	
	
	
# 	for (i,num) in enumerate(num_list)
# 		for (k,τ) in enumerate(τ_list)
# 			for (j,T) in enumerate(Tinter_list)
			
# 				params_list = Dict(
# 	"number_repetitions" => num,
#     "τ_integrator" => τ,
# 	"value_stim_f1" => 0.2,
# 	"Tinter" => T,#
# 	"Tstim" => 0.5,
# 	"Tfin" => 1.0,
# 	"value_stim_f2" => 0.0,    
# 	"initial_time" => 4.0
# )
					
# 		lc,oddball,sim = run_simu1_script(params_list)

# 	tot_MMN[i,k,j] =compute_MMN_oddball(lc,oddball,params_list)
# 	end
				
# 		end
		
# 	end
	
	
	
end

# ╔═╡ 4907dfcb-8679-4b0f-b4d8-d340160096ef
begin 
	dict_tot_MMN = load(datadir("temp\\script9-temp1.jld"))
	load_tot_MMN = dict_tot_MMN["tot_MMN"]
end

# ╔═╡ 9435048e-f62c-40cc-8fe1-622f2ee3fe7f
begin
	# in this oddball context, the exact time of computing the effect is known
	
	fig_temp = Figure()
	index_τ = 3
	τ_value = τ_list[index_τ]
	ax_temp = Axis(fig_temp[1,1],title="MMN vs ISI - τ = $τ_value s",xlabel="ISI",ylabel="Δ firing rate (Hz)")
	for (i,n) in enumerate(num_list)

		lines!(ax_temp,Tinter_list,load_tot_MMN[i,index_τ,:],linewidth=3,label="rep = $n")
	end
	fig_temp
	fig_temp[1, 2] = Legend(fig_temp, ax_temp, "Legend", framevisible = false)
	fig_temp

end

# ╔═╡ f488eb2a-c827-4ba2-b74d-7d7d5d42b5ad
begin
	# let's fix the ISI to 1 sec and the number of rep to 3
	#τ integrator could be a parameter to vary
	# 
	
		fig_temp_τ = Figure()
	index_Tinter = 2
	Tinter_value = Tinter_list[index_Tinter]
	ax_temp_τ = Axis(fig_temp_τ[1,1],title="MMN vs τ - Tinter = $Tinter_value s",xlabel="τ",ylabel="Δ firing rate (Hz)")
	for (i,n) in enumerate(num_list)

		lines!(ax_temp_τ,τ_list,load_tot_MMN[i,:,index_Tinter],linewidth=3,label="rep = $n")
	end
	fig_temp_τ[1, 2] = Legend(fig_temp_τ, ax_temp_τ, "Legend", framevisible = false)
	fig_temp_τ

end

# ╔═╡ fe8683fe-b658-4088-b27c-7b72b4626997
md"""
Index of τ for the heatmap: 
$(@bind index_heatmap PlutoUI.Slider(1:length(τ_list),default=3, show_value=true))

"""

# ╔═╡ 715e4e5c-28b9-4cb9-94fa-1644970903af
begin
	
	#ax_fixed_τ = Axis(heatmap_fixed_τ[1,1])
#	index_heatmap = 3
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

# ╔═╡ b616c49e-7051-45bf-8515-88e19f5c7b4b
md"""
Index of Tinter for the heatmap: 
$(@bind index_heatmap_Tinter PlutoUI.Slider(1:length(Tinter_list),default=2, show_value=true))

"""

# ╔═╡ 3aed38b3-fa55-4f55-8e42-023bf59d35f8
begin
	

	value_heatmap_Tinter = Tinter_list[index_heatmap_Tinter]
	heatmap_fixed_Tinter = Figure()
	
	ax_hm_Tinter,hm_Tinter = heatmap(heatmap_fixed_Tinter[1,1],num_list,τ_list, load_tot_MMN[:,:,index_heatmap_Tinter])
	Colorbar(heatmap_fixed_Tinter[1, 2],hm_Tinter,label="MMN effect (Hz)")
	SupTitle_Tinter = Label(heatmap_fixed_Tinter[0,1:2], "Heatmap of MMN - Tinter integrator = $value_heatmap_Tinter s")
	ax_hm_Tinter.xlabel = "Number of repetitions"
	ax_hm_Tinter.ylabel = "τ integrator"
	save(plotsdir("notebook9-heatmap_Tinter=$value_heatmap_Tinter.png"),heatmap_fixed_Tinter)
	heatmap_fixed_Tinter
	
	
	
end

# ╔═╡ 8d0df090-562e-433e-a9b8-e95934a6876d
md"""
Index of number of repetitions for the heatmap: 
$(@bind index_heatmap_rep PlutoUI.Slider(1:length(num_list),default=2, show_value=true))

"""

# ╔═╡ 68878b41-c4ee-4260-bc38-34449d3b8dbb
begin
	

	value_heatmap_rep = num_list[index_heatmap_rep]
	heatmap_fixed_rep = Figure()
	
	ax_hm_rep,hm_rep = heatmap(heatmap_fixed_rep[1,1],τ_list,Tinter_list, load_tot_MMN[index_heatmap_rep,:,:])
	Colorbar(heatmap_fixed_rep[1, 2],hm_rep,label="MMN effect (Hz)")
	SupTitle_rep = Label(heatmap_fixed_rep[0,1:2], "Heatmap of MMN - Number of repetitions = $value_heatmap_rep ")
	ax_hm_rep.xlabel= "τ integrator"
	ax_hm_rep.ylabel = "Tinter (s)"
	
	save(plotsdir("notebook9-heatmap_rep=$value_heatmap_rep.png"),heatmap_fixed_rep)
	heatmap_fixed_rep
	
	
	
end

# ╔═╡ f2a23925-983d-4469-bf3f-84ab1204efa8
begin
	# this figure represents the time delayu of MMN in the model
	temp = compute_MMN_time(lc,oddball, params)
	lines(temp)
end

# ╔═╡ 43b7c57c-8d07-4590-bae0-a93d92cb2478
temp[50]

# ╔═╡ 0c9e9a3b-c3cf-4bcc-b3ad-f8c6f6103d17
md""" # Cutting off some connectivities

"""

# ╔═╡ cd30528d-997d-4fcb-b8df-115e01ab5ee3


# ╔═╡ b1f796c2-a937-4f5b-bf2b-6225361f85ba


# ╔═╡ 9655a270-7533-4fdd-a368-37cfe2fe4fd9


# ╔═╡ 5d8a17e1-db3e-4997-873c-035f09295fb1


# ╔═╡ 7f957bd0-b28b-4327-9ccd-d578f729a013
begin
	function run_simu2_script(params)
		stim_1, stim_2 = create_deterministic_oscillations_oddball(params)
	
	@unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

	t_tot = length(stim_1)
	
	lc = construct_two_local_microcircuit_integrator(facilitation=true, depression=true, adaptation=true,  τ_integrator = τ_integrator ,time_tot = t_tot)
            
            
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

# ╔═╡ 3939c616-87e0-4dc4-897a-c4ada55208f7
begin
	lc2,oddball2,sim2 = run_simu2_script(params)
fig2 = plot_local_circuit(lc2, sim2, oddball2, params)
end

# ╔═╡ 81b87579-e082-4162-923c-86877df27222
fig2_syn = plot_local_circuit_synapses(lc2, sim2, oddball2, params)


# ╔═╡ f4d06694-8af9-4b69-9b9b-8153230d6612


# ╔═╡ de6f0572-6331-4338-8edf-29b09fe8ea74


# ╔═╡ d7ca103f-b799-4fb6-a8ec-fce2eda6fb97


# ╔═╡ 1ec8e014-7840-4a63-b2db-fa91258903ef
begin
	
	fig_stim = Figure()
	ax_stim = Axis(fig_stim[1,1])
	
	lines!(ax_stim,stim_1)
	lines!(ax_stim,stim_2)

	
	fig_stim
	
end

# ╔═╡ b0a25481-144f-4be7-848f-0adac1bed1fe
length(lc[1].list_soma[1].r_save)

# ╔═╡ 330e7ebc-5a91-4759-b759-0377bb9e4cfd
length(stim_1)

# ╔═╡ eb5d26c0-51b3-4c29-8657-dda970a5a587


# ╔═╡ Cell order:
# ╠═517c9c3e-0cce-11ec-0a91-07fda58ce84d
# ╠═9fd2c4b4-975b-45b4-9901-5c22180ce984
# ╠═859449c2-238e-430c-aebc-0d226314df0b
# ╠═eab23a48-c873-408b-91da-5c208a82f117
# ╠═70649557-de48-424d-ada1-12e6c5540417
# ╠═b3bbd648-8b82-414b-ae39-5cf70f70f06a
# ╠═fd27c550-f629-48ef-b804-71319d9532ec
# ╠═6a0d37bf-f458-4273-909d-4f06e440b0fe
# ╠═707cae76-fdf3-4d6b-856b-b57a65e97ec8
# ╠═2bc2f601-a36d-498c-b208-a5099c9e2d59
# ╠═e9d59c50-5173-4674-bcfb-97dcde7b7035
# ╠═4907dfcb-8679-4b0f-b4d8-d340160096ef
# ╠═9435048e-f62c-40cc-8fe1-622f2ee3fe7f
# ╠═f488eb2a-c827-4ba2-b74d-7d7d5d42b5ad
# ╟─fe8683fe-b658-4088-b27c-7b72b4626997
# ╠═715e4e5c-28b9-4cb9-94fa-1644970903af
# ╟─b616c49e-7051-45bf-8515-88e19f5c7b4b
# ╠═3aed38b3-fa55-4f55-8e42-023bf59d35f8
# ╟─8d0df090-562e-433e-a9b8-e95934a6876d
# ╟─68878b41-c4ee-4260-bc38-34449d3b8dbb
# ╠═f2a23925-983d-4469-bf3f-84ab1204efa8
# ╠═43b7c57c-8d07-4590-bae0-a93d92cb2478
# ╠═0c9e9a3b-c3cf-4bcc-b3ad-f8c6f6103d17
# ╠═cd30528d-997d-4fcb-b8df-115e01ab5ee3
# ╠═b1f796c2-a937-4f5b-bf2b-6225361f85ba
# ╠═9655a270-7533-4fdd-a368-37cfe2fe4fd9
# ╠═5d8a17e1-db3e-4997-873c-035f09295fb1
# ╠═7f957bd0-b28b-4327-9ccd-d578f729a013
# ╠═3939c616-87e0-4dc4-897a-c4ada55208f7
# ╠═81b87579-e082-4162-923c-86877df27222
# ╠═f4d06694-8af9-4b69-9b9b-8153230d6612
# ╠═de6f0572-6331-4338-8edf-29b09fe8ea74
# ╠═d7ca103f-b799-4fb6-a8ec-fce2eda6fb97
# ╠═1ec8e014-7840-4a63-b2db-fa91258903ef
# ╠═b0a25481-144f-4be7-848f-0adac1bed1fe
# ╠═330e7ebc-5a91-4759-b759-0377bb9e4cfd
# ╠═eb5d26c0-51b3-4c29-8657-dda970a5a587
