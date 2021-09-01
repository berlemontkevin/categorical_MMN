### A Pluto.jl notebook ###
# v0.15.1

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
	using LsqFit
	using ColorSchemes
	
	set_theme!(theme_light())

	
end

# ╔═╡ 33aa235e-fec6-11eb-2e65-9b8f1169abae
md""" # Analysis of the integrator time constant and MMN time constant


From data of script 3


The idea is to try what kind of impact does frequency and stimulus strength have on the time constant of the decreasing MMN.


"""

# ╔═╡ 4607fe2b-ebe7-4c1a-9f50-302aa0254987
begin
	# Plots variables
	noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")

end

# ╔═╡ 6bac87c9-84cc-441b-80e0-cfa5ddddd874
md"""
The analysis of MMN time constant is the following. 

The different parameters are IS, τ integrator, value of the stimlus, frequency. For every set of parameters, the simulation is run for 6mins. 



TODO: graphic that explains how to measure the mean firing rate.




"""

# ╔═╡ af717f0e-2c0b-4691-87c9-9e13c8c097ee
md"""
**What does the mean firing rate as function of ISI looks like?**

"""

# ╔═╡ e975c7bd-ed9c-4fbc-bc88-b4fb66a846cf
md"""

Value of stimuli:
$(@bind value_stim PlutoUI.Slider(0.1 : 0.1 : 0.6, show_value=true) )

Frequencies:
$(@bind freq PlutoUI.Slider(0.5 : 0.1 : 0.9, show_value=true) )

"""

# ╔═╡ 92677b31-eb38-437f-bbb8-8d12d7840adb
begin
fig_ISI = Figure( font = noto_sans)
	
	
	# value_stim = 0.3
	end_tc = 3.0
	# freq = 0.8
	
	τstep = 0.4

	
    ax_ISI_1 = fig_ISI[1, 1] = Axis(fig_ISI, title = "Mean firing rate vs ISI for stim $value_stim and f=$freq", xlabel="ISI", ylabel = "Δ mean firing rate")
ax_ISI_2 = fig_ISI[2, 1] = Axis(fig_ISI, title = "Mean firing rate vs ISI for subsystem2", xlabel="ISI", ylabel = "Δ mean firing rate")	


	
	tc_list = collect(0.1:τstep:end_tc)
	
	colors = to_colormap(:viridis, length(tc_list))

	subsystem1 = zeros(length(tc_list),length(collect(0.2:0.3:15.0)))
	subsystem2 = zeros(length(tc_list),length(collect(0.2:0.3:15.0)))
	
	
	for j=1:length(tc_list)
		
		tc_int = tc_list[j]
		
		fig1_dict = Dict(
	"f1" => freq,
	"value_stim_f1" => value_stim,
	"Tinter" => collect(0.2:0.3:15.0),
	"Tstim" => 0.5,
	"Tfin" => 360.0,
	"value_stim_f2" => 0.0,
    "τ_integrator" => tc_int
)
	
	for (i,dict_temp) in enumerate(dict_list(fig1_dict))
		
		data = load(datadir("sims","script3",savename("oddball_task",dict_temp,"jld2")))

	mean_value = get_mean_firing_rate(data,"microcircuit1-ecell1","current-to-e1",dict_temp["value_stim_f1"])
	mean_value_sys2 = get_mean_firing_rate(data,"microcircuit2-ecell1","current-to-e2",dict_temp["value_stim_f1"])	
			
			
		subsystem1[j,i]= mean_value# - mean_value_sys2
		subsystem2[j,i]= mean_value_sys2
	end
	
		
	

	
		
	lb = [-Inf, 0.01, 0.0]
	ub = [0.0, Inf, Inf]
	@. model(x, p) = p[1]*exp(-x*p[2]) + p[3]
	p0 = [-4.0,0.02,25.0]
	
	
	ind = map(!,isnan.(subsystem1[j,:]))

	fit = curve_fit(model, collect(0.2:0.3:15.0)[ind], subsystem1[j,ind], p0, lower = lb, upper = ub)

	
			lines!(ax_ISI_1,collect(0.2:0.3:15.0),subsystem1[j,:],linewidth=5,color=colors[j], label = "τ int - $tc_int")
	
	
			lines!(ax_ISI_2,collect(0.2:0.3:15.0),subsystem2[j,:],linewidth=5,color=colors[j])

		
	end
	
	fig_ISI[1:2, 2] = Legend(fig_ISI, ax_ISI_1, "Legend", framevisible = false)
	
	fig_ISI
end

# ╔═╡ 864bc534-db2d-4bb7-962e-d348656e0adf
md"""
First thing that can be seen is that the mean firing rate between the two subnetworks is different. Moreover,  the increase of firing rate seens to be exponential in all the cases.

Important point is that for the MMN effect, what matters is the difference between firing rate and not just the firing rate of one pyramidal cell.



"""

# ╔═╡ b9e8a257-6b81-4c5b-9c43-bae728885e55
md"""
Now due to noise (low frequency of apparition of some stimuli), in the next part I fit both of tese curves before doing the differences between them.





Nota bene: longer simulations will be run on the cluster to obtain the true dynamical behavior.
"""

# ╔═╡ 5098df57-da12-4d3b-9452-0080994b3c18
 begin
		
    fig1 = Figure( font = noto_sans)
	
	x = nothing
	
	
	list_MMN_tc = Float64[]
	
    ax1 = fig1[1, 1:2] = Axis(fig1, title = "Mean firing rate vs ISI for stim $value_stim", xlabel="Inter stim", ylabel = "Δ mean firing rate")
	
	
	
	fig1_list = zeros(length(tc_list),length(collect(0.2:0.3:15.0)))
	value_sys2 = zeros(length(tc_list),length(collect(0.2:0.3:15.0)))
	
	
	for j=1:length(tc_list)
		
		tc_int = tc_list[j]
		
		fig1_dict = Dict(
	"f1" => freq,
	"value_stim_f1" => value_stim,
	"Tinter" => collect(0.2:0.3:15.0),
	"Tstim" => 0.5,
	"Tfin" => 360.0,
	"value_stim_f2" => 0.0,
    "τ_integrator" => tc_int
)
	
	for (i,dict_temp) in enumerate(dict_list(fig1_dict))
		
		data = load(datadir("sims","script3",savename("oddball_task",dict_temp,"jld2")))

	mean_value = get_mean_firing_rate(data,"microcircuit1-ecell1","current-to-e1",dict_temp["value_stim_f1"])
	mean_value_sys2 = get_mean_firing_rate(data,"microcircuit2-ecell1","current-to-e2",dict_temp["value_stim_f1"])	
			
			
		fig1_list[j,i]= mean_value# - mean_value_sys2
		value_sys2[j,i] = mean_value_sys2
	end
	
		
		fig1_list[j,:] = (fig1_list[j,:])# .+ fig1_list[j,end])./fig1_list[j,end]
		
	

	
		
	lb = [-Inf, 0.01, 0.0]
	ub = [0.0, Inf, Inf]

	@. model(x, p) = p[1]*exp(-x*p[2]) + p[3]#p[1]*exp(-x*p[2])p[3] - p[1]*tanh(x*p[2])
	
	p0 = [-4.0,0.02,25.0]
				ind = map(!,isnan.(fig1_list[j,:]))

	fit = curve_fit(model, collect(0.2:0.3:15.0)[ind], fig1_list[j,ind], p0, lower = lb, upper = ub)
		ind = map(!,isnan.(value_sys2[j,:]))
		# x = md"temp $ind"
	fit2 = curve_fit(model, collect(0.2:0.3:15.0)[ind], value_sys2[j,ind], p0, lower = lb, upper = ub)
# a = coef(fit)[2]
# 		lines!(ax1,collect(0.2:0.3:15.0),fig1_list[j,:],linewidth=5,color=:red, label = "tc int - $a")

	#push!(list_MMN_tc, 1.0/coef(fit)[2])
		fig1_list[j,:] = model(collect(0.2:0.3:15.0), coef(fit)) .-model(collect(0.2:0.3:15.0), coef(fit2))
	#lines!(ax1,collect(0.2:0.3:15.0),model(collect(0.2:0.3:15.0),coef(fit)),linewidth=2, colormap=ColorSchemes.lajolla, linestyle=:dot)
			lines!(ax1,collect(0.2:0.3:15.0),fig1_list[j,:],linewidth=5,color=colors[j], label = "tc int - $tc_int")

		
	end
	
	fig1[1, 3] = Legend(fig1, ax1, "Legend", framevisible = false)
 #leg1 = Colorbar(fig1, ax1, label = "stim value", ticklabelsize = 14,
    #            labelpadding = 5, width = 10)
	
	fig1
	# x
end

# ╔═╡ 4e8f9c4b-c7df-4b43-a3d2-4022d98e6c85
begin
	# cell to plot all the time constants of MMN
	
	  fig_tc = Figure( font = noto_sans)
	
	
	value_stim_list = collect(0.1:0.1:0.6)
	value_freq = collect(0.5:0.1:0.8)
	end_tc2 = 2.2
	tc_list2 = collect(0.1:0.2:end_tc2)

	
	average_list = zeros(length(value_stim_list),length(value_freq))
	average_percent_list = zeros(length(value_stim_list),length(value_freq))
	
	
	
	
    ax_tc = fig_tc[1, 1:2] = Axis(fig_tc, title = "Mean firing rate vs ISI", xlabel="TC integrator", ylabel = "Δ mean firing rate")
	

	Δ_list = zeros(length(value_stim_list),length(value_freq),length(tc_list2))

list_map= [:reds, :greens, :blues, :grays, :thermal, :viridis, :grays]
	
			for w=1:length(value_freq)

	for q=1:length(value_stim_list)
		
		value_tc = value_stim_list[q]
		
	
			freq_tc = value_freq[w]
			
	figtc_list = zeros(length(tc_list2),length(collect(0.2:0.3:15.0)))
	
				colors_list = to_colormap(list_map[w], length(value_stim_list))

	
	for j=1:length(tc_list2)
		
		tc_int = tc_list2[j]
		
		fig1_dict = Dict(
	"f1" => freq_tc,
	"value_stim_f1" => value_tc,
	"Tinter" => collect(0.2:0.3:15.0),
	"Tstim" => 0.5,
	"Tfin" => 360.0,
	"value_stim_f2" => 0.0,
    "τ_integrator" => tc_int
)
	
	for (i,dict_temp) in enumerate(dict_list(fig1_dict))
		
		data = load(datadir("sims","script3",savename("oddball_task",dict_temp,"jld2")))

	mean_value = get_mean_firing_rate(data,"microcircuit1-ecell1","current-to-e1",dict_temp["value_stim_f1"])
		
		figtc_list[j,i]= mean_value
	end
		
						average_list[q,w] += figtc_list[j,1]/length(tc_list2) 
				
				figtc_list[j,:] = (-figtc_list[j,:] .+ figtc_list[j,end])./figtc_list[j,end]
		
	

	
		
	
	@. model(x, p) = p[3] - p[1]*tanh(x*p[2])#p[1]*exp(-x*p[2])
	
	p0 = [1.0,0.1,25.0]
	
	fit = curve_fit(model, collect(0.2:0.3:15.0), figtc_list[j,:], p0)
		
		
	Δ_list[q,w,j] = 1.0/coef(fit)[2]
		

		
	end
							average_percent_list[q,w] = mean(figtc_list[:,1]) 

			lines!(ax_tc, tc_list2, Δ_list[q,w,:], label="freq =$freq_tc and value_stim = $value_tc", linewidth = 3, color = colors_list[q])
			
		end
	end
		fig_tc[1, 3] = Legend(fig_tc, ax_tc, "Legend", framevisible = false)

	fig_tc
end

# ╔═╡ 09d8d74f-a8e2-494b-b6e6-06694426666c
begin
	 fig_color = Figure( font = noto_sans)
     ax_color1 = Axis(fig_color, xlabel = "τ integrator", ylabel = "τ MMN effect", title = "Freq of 0.5 ")
     blues_list = to_colormap(:blues, length(value_stim_list))
	for i=1:length(value_stim_list)
     lines!(ax_color1, tc_list2, Δ_list[i,1,:], color= blues_list[i], linewidth = 3)
	end
	fig_color[1, 1] = ax_color1
	
	 ax_color2 = Axis(fig_color,xlabel = "τ integrator", ylabel = "τ MMN effect", title = "Freq of 0.6 ")
     reds_list = to_colormap(:reds, length(value_stim_list))
	for i=1:length(value_stim_list)
     lines!(ax_color2, tc_list2, Δ_list[i,2,:], color= reds_list[i], linewidth = 3)
	end
	fig_color[1, 2] = ax_color2
	
	
	 ax_color3 = Axis(fig_color,xlabel = "τ integrator", ylabel = "τ MMN effect", title = "Freq of 0.7 ")
     greens_list = to_colormap(:greens, length(value_stim_list))
	for i=1:length(value_stim_list)
     lines!(ax_color3, tc_list2, Δ_list[i,3,:], color= greens_list[i], linewidth = 3)
	end
	fig_color[2, 1] = ax_color3
	
	
	 ax_color4 = Axis(fig_color,xlabel = "τ integrator", ylabel = "τ MMN effect", title = "Freq of 0.8 ")
     grays_list = to_colormap(:grays, length(value_stim_list))
	for i=1:length(value_stim_list)
		lines!(ax_color4, tc_list2, Δ_list[i,4,:], color= grays_list[end-i+1], linewidth = 3)
	end
	fig_color[2, 2] = ax_color4
	
	  # cbar_fig_color = Colorbar(fig_color, value_stim_list, label = "color value", ticklabelsize = 14, labelpadding = 5, width = 10)
	
	
	fig_color
	
	
end

# ╔═╡ f1421b5d-e88e-49fb-82d7-bf7289fd4207
md"""
### Heatmap of the effect

"""

# ╔═╡ be96e4fa-65cb-4fad-8c9d-a71fa99adaef
begin
	 fig = Figure( font = noto_sans)
     ax = Axis(fig, title = "Heatmap for stim = 0.1 ")
     hmap = heatmap!(ax,value_freq, tc_list2, Δ_list[1,:,:], colormap = :blues)
	 cbar = Colorbar(fig, hmap, label = "Δ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig[1, 1] = ax
    fig[1, 2] = cbar
	
	ax2 = Axis(fig, title = "Heatmap for stim = 0.2 ")
     hmap2 = heatmap!(ax2,value_freq, tc_list2, Δ_list[2,:,:], colormap = :blues)
	 cbar2 = Colorbar(fig, hmap2, label = "Δ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig[1, 3] = ax2
    fig[1, 4] = cbar2
	
	ax3 = Axis(fig, title = "Heatmap for stim = 0.3 ")
     hmap3 = heatmap!(ax3,value_freq, tc_list2, Δ_list[3,:,:], colormap = :blues)
	 cbar3 = Colorbar(fig, hmap3, label = "Δ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig[2, 1] = ax3
    fig[2, 2] = cbar3
	
	ax4 = Axis(fig, title = "Heatmap for stim = 0.4 ")
     hmap4 = heatmap!(ax4,value_freq, tc_list2, Δ_list[4,:,:], colormap = :blues)
	 cbar4 = Colorbar(fig, hmap4, label = "Δ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig[2, 3] = ax4
    fig[2, 4] = cbar4
	
	fig
end

# ╔═╡ 005e59fe-4f25-4f3e-8712-e8a5c761fb64
md"""
#### Strenght of the effect
"""

# ╔═╡ 63a0d410-1c90-413a-9d09-e800c71eff55
begin
	 fig_strength = Figure( font = noto_sans)
     ax_strength = Axis(fig_strength, title = "Heatmap for stim = 0.1 ")
     hmap_strength = heatmap!(ax_strength,value_stim_list, value_freq, average_list[:,:], colormap = :blues)
	 cbar_strength = Colorbar(fig_strength, hmap_strength, label = "Δ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig_strength[1, 1] = ax_strength
    fig_strength[1, 2] = cbar_strength
	
# 	ax_strength2 = Axis(fig_strength, title = "Heatmap for stim = 0.2 ")
#      hmap_strength2 = heatmap!(ax_strength2,value_freq, value_stim_list, average_list[:,:], colormap = :blues)
# 	 cbar_strength2 = Colorbar(fig_strength, hmap_strength2, label = "Δ values", width = 15,           ticksize=15, tickalign = 1, height = Relative(3.55/4))
# 	fig_strength[1, 3] = ax_strength2
#     fig_strength[1, 4] = cbar_strength2
	
# 	ax_strength3 = Axis(fig_strength, title = "Heatmap for stim = 0.3 ")
#      hmap_strength3 = heatmap!(ax3,value_freq, value_stim_list, average_list[:,:], colormap = :blues)
# 	 cbar_strength3 = Colorbar(fig_strength, hmap_strength3, label = "Δ values", width = 15,
#                 ticksize=15, tickalign = 1, height = Relative(3.55/4))
# 	fig_strength[2, 1] = ax_strength3
#     fig_strength[2, 2] = cbar_strength3
	
# 	ax4 = Axis(fig_strength, title = "Heatmap for stim = 0.4 ")
#      hmap_strength4 = heatmap!(ax_strength4,value_freq, value_stim_list, average_list[:,:], colormap = :blues)
# 	 cbar_strength4 = Colorbar(fig_strength, hmap_strength4, label = "Δ values", width = 15,
#                 ticksize=15, tickalign = 1, height = Relative(3.55/4))
# 	fig[2, 3] = ax_strength4
#     fig[2, 4] = cbar_strength4
	
	fig_strength
end

# ╔═╡ 4a53013c-91ca-453f-9b0d-26096ccf822a
lines(average_list[:,3])

# ╔═╡ Cell order:
# ╟─33aa235e-fec6-11eb-2e65-9b8f1169abae
# ╟─63c5c372-10ac-401e-81bf-70c929e3b313
# ╟─9b053399-c4f2-46b5-9890-f1b538c2bfac
# ╠═4607fe2b-ebe7-4c1a-9f50-302aa0254987
# ╠═6bac87c9-84cc-441b-80e0-cfa5ddddd874
# ╟─af717f0e-2c0b-4691-87c9-9e13c8c097ee
# ╟─e975c7bd-ed9c-4fbc-bc88-b4fb66a846cf
# ╟─92677b31-eb38-437f-bbb8-8d12d7840adb
# ╠═864bc534-db2d-4bb7-962e-d348656e0adf
# ╠═b9e8a257-6b81-4c5b-9c43-bae728885e55
# ╠═5098df57-da12-4d3b-9452-0080994b3c18
# ╠═4e8f9c4b-c7df-4b43-a3d2-4022d98e6c85
# ╟─09d8d74f-a8e2-494b-b6e6-06694426666c
# ╠═f1421b5d-e88e-49fb-82d7-bf7289fd4207
# ╠═be96e4fa-65cb-4fad-8c9d-a71fa99adaef
# ╠═005e59fe-4f25-4f3e-8712-e8a5c761fb64
# ╠═63a0d410-1c90-413a-9d09-e800c71eff55
# ╠═4a53013c-91ca-453f-9b0d-26096ccf822a
