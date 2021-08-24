### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 63c5c372-10ac-401e-81bf-70c929e3b313
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())

end

# ╔═╡ 9b053399-c4f2-46b5-9890-f1b538c2bfac
begin
	using Revise
	using MyNeurosciencePackage
	using CairoMakie
	using Parameters, JLD2, Statistics	
	using DataFrames, PlutoUI
	
	
	set_theme!(theme_light())

	
end

# ╔═╡ 33aa235e-fec6-11eb-2e65-9b8f1169abae
md""" # Analysis of the integrator time constant and MMN time constant


From data of script 3

"""

# ╔═╡ 4607fe2b-ebe7-4c1a-9f50-302aa0254987
begin
	# Plots variables
	noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")

end

# ╔═╡ 5098df57-da12-4d3b-9452-0080994b3c18
 begin
		using LsqFit
	using ColorSchemes
    fig1 = Figure( font = noto_sans)
	
	
	value_stim = 0.3
	end_tc = 2.1
	freq = 0.8
	
	list_MMN_tc = Float64[]
	
    ax1 = fig1[1, 1:2] = Axis(fig1, title = "Mean firing rate vs ISI for stim $value_stim", xlabel="Inter stim", ylabel = "Δ mean firing rate")
	
	
	tc_list = collect(0.1:0.1:end_tc)
	
	colors = to_colormap(:viridis, length(tc_list))

	fig1_list = zeros(length(tc_list),length(collect(0.2:0.3:15.0)))
	
	
	for j=1:length(tc_list)
		
		tc_int = tc_list[j]
		
		fig1_dict = Dict(
	"f1" => freq,
	"value_stim_f1" => value_stim,
	"Tinter" => collect(0.2:0.3:15.0),
	"Tstim" => 0.5,
	"Tfin" => 300.0,
	"value_stim_f2" => 0.0,
    "τ_integrator" => tc_int
)
	
	for (i,dict_temp) in enumerate(dict_list(fig1_dict))
		
		data = load(datadir("sims","script3",savename("oddball_task",dict_temp,"jld2")))

	mean_value = get_mean_firing_rate(data,"microcircuit1-ecell1","current-to-e1",dict_temp["value_stim_f1"])
		
		fig1_list[j,i]= mean_value
	end
	
		
		fig1_list[j,:] = (-fig1_list[j,:] .+ fig1_list[j,end])./fig1_list[j,end]
		
	

	
		
	
	@. model(x, p) = p[3] - p[1]*tanh(x*p[2])#p[1]*exp(-x*p[2])
	
	p0 = [1.0,0.1,25.0]
	
	fit = curve_fit(model, collect(0.2:0.3:15.0), fig1_list[j,:], p0)
		
		
	push!(list_MMN_tc, 1.0/coef(fit)[2])
		
	#lines!(ax1,collect(0.2:0.3:15.0),model(collect(0.2:0.3:15.0),coef(fit)),linewidth=2, colormap=ColorSchemes.lajolla, linestyle=:dot)
			lines!(ax1,collect(0.2:0.3:15.0),fig1_list[j,:],linewidth=5,color=colors[j], label = "tc int - $tc_int")

		
	end
	
	fig1[1, 3] = Legend(fig1, ax1, "Legend", framevisible = false)

	
	fig1
end

# ╔═╡ 3ca24a85-e209-403e-af2d-d671c2394e63
begin
	# Fitted time constants as function of integrator time constant
	lines(tc_list,list_MMN_tc,linewidth=4)
	
	
end

# ╔═╡ 4e8f9c4b-c7df-4b43-a3d2-4022d98e6c85
begin
	# cell to plot all the time constants of MMN
	
	  fig_tc = Figure( font = noto_sans)
	
	
	value_stim_list = collect(0.2:0.1:0.4)
	value_freq = collect(0.7:0.1:0.9)
	end_tc2 = 3.0
	tc_list2 = collect(0.1:0.5:end_tc2)

	
	

	
    ax_tc = fig_tc[1, 1:2] = Axis(fig_tc, title = "Mean firing rate vs ISI", xlabel="TC integrator", ylabel = "Δ mean firing rate")
	

	list_MMN = zeros(length(value_stim_list),length(value_freq),length(tc_list2))

list_map= [:reds, :greens, :blues]
	
	
	for q=1:length(value_stim_list)
		
		value_tc = value_stim_list[q]
		
		for w=1:length(value_freq)
	
			freq_tc = value_freq[w]
			
	figtc_list = zeros(length(tc_list2),length(collect(0.2:0.6:15.0)))
	
				colors_list = to_colormap(list_map[q], length(value_freq))

	
	for j=1:length(tc_list2)
		
		tc_int = tc_list[j]
		
		fig1_dict = Dict(
	"f1" => freq_tc,
	"value_stim_f1" => value_tc,
	"Tinter" => collect(0.2:0.6:15.0),
	"Tstim" => 0.5,
	"Tfin" => 300.0,
	"value_stim_f2" => 0.0,
    "τ_integrator" => tc_int
)
	
	for (i,dict_temp) in enumerate(dict_list(fig1_dict))
		
		data = load(datadir("sims","script3",savename("oddball_task",dict_temp,"jld2")))

	mean_value = get_mean_firing_rate(data,"microcircuit1-ecell1","current-to-e1",dict_temp["value_stim_f1"])
		
		figtc_list[j,i]= mean_value
	end
	
		
				figtc_list[j,:] = (-figtc_list[j,:] .+ figtc_list[j,end])./figtc_list[j,end]
		
	

	
		
	
	@. model(x, p) = p[3] - p[1]*tanh(x*p[2])#p[1]*exp(-x*p[2])
	
	p0 = [1.0,0.1,25.0]
	
	fit = curve_fit(model, collect(0.2:0.6:15.0), figtc_list[j,:], p0)
		
		
	list_MMN[q,w,j] = 1.0/coef(fit)[2]
		

		
	end
			lines!(ax_tc, tc_list2, list_MMN[q,w,:], label="freq =$freq_tc and value_stim = $value_tc", linewidth = 3, color = colors_list[w])
			
		end
	end
		fig_tc[1, 3] = Legend(fig_tc, ax_tc, "Legend", framevisible = false)

	fig_tc
end

# ╔═╡ f1421b5d-e88e-49fb-82d7-bf7289fd4207


# ╔═╡ be96e4fa-65cb-4fad-8c9d-a71fa99adaef


# ╔═╡ 005e59fe-4f25-4f3e-8712-e8a5c761fb64


# ╔═╡ 4a53013c-91ca-453f-9b0d-26096ccf822a


# ╔═╡ Cell order:
# ╠═33aa235e-fec6-11eb-2e65-9b8f1169abae
# ╠═63c5c372-10ac-401e-81bf-70c929e3b313
# ╠═9b053399-c4f2-46b5-9890-f1b538c2bfac
# ╠═4607fe2b-ebe7-4c1a-9f50-302aa0254987
# ╠═5098df57-da12-4d3b-9452-0080994b3c18
# ╠═3ca24a85-e209-403e-af2d-d671c2394e63
# ╠═4e8f9c4b-c7df-4b43-a3d2-4022d98e6c85
# ╠═f1421b5d-e88e-49fb-82d7-bf7289fd4207
# ╠═be96e4fa-65cb-4fad-8c9d-a71fa99adaef
# ╠═005e59fe-4f25-4f3e-8712-e8a5c761fb64
# ╠═4a53013c-91ca-453f-9b0d-26096ccf822a
