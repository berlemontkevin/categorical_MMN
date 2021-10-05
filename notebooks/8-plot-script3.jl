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
	
	quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
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
$(@bind freq PlutoUI.Slider(0.5 : 0.1 : 0.9, show_value=true, default = 0.6) )

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

	ISI_list = vcat(collect(0.2:0.3:15.0),collect(16.0:1.0:20.0))

	
	tc_list = vcat(collect(0.1:τstep:end_tc),collect(3.1:3.0:16.0))#collect(0.1:τstep:end_tc)
	
	colors = to_colormap(:viridis, length(tc_list))

	subsystem1 = zeros(length(tc_list),length(ISI_list))
	subsystem2 = zeros(length(tc_list),length(ISI_list))
	
	
	for j=1:length(tc_list)
		
		tc_int = tc_list[j]
		
		fig1_dict = Dict(
	"f1" => freq,
	"value_stim_f1" => value_stim,
	"Tinter" => ISI_list,
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
	
		
	

	
	
	
	
	
	
			lines!(ax_ISI_1,ISI_list,subsystem1[j,:],linewidth=5,color=colors[j], label = "τ int - $tc_int")
	
	
			lines!(ax_ISI_2,ISI_list,subsystem2[j,:],linewidth=5,color=colors[j])

		
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
	
		
	lb = [-Inf, 0.01, 0.0]
	ub = [0.0, Inf, Inf]#[0.0, Inf, Inf]
	@. model(x, p) = p[1]*exp(-x*p[2]) + p[3]#p[1]*tanh(x*p[2]) + p[3]#p[1]*exp(-x*p[2]) + p[3]
	p0 = [-1.0,0.02,25.0]
		
    fig_Δ = Figure( font = noto_sans)
	
	
	
    ax_Δ = fig_Δ[1, 1:2] = Axis(fig_Δ, title = "Δ Mean firing rate vs ISI for stim $value_stim and f=$freq", xlabel="Inter stim", ylabel = "Δ mean firing rate")
	
	for j=1:length(tc_list)
		
	subsystem1_ind = map(!,isnan.(subsystem1[j,:]))
	subsystem2_ind = map(!,isnan.(subsystem2[j,:]))

	fit1 = curve_fit(model, ISI_list[subsystem1_ind], subsystem1[j,subsystem1_ind], p0, lower = lb, upper = ub)
		
	fit2 = curve_fit(model, ISI_list[subsystem2_ind], subsystem2[j,subsystem2_ind], p0, lower = lb, upper = ub)
		
		lines!(ax_Δ,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
		
	end
	
	
	fig_Δ
end

# ╔═╡ c14c58d2-ff7f-434c-ab4e-64df6dd1b020
md"""
Let's now try to compare different frequencies together to see some behavior trend.


$(@bind value_stim_2 PlutoUI.Slider(0.1:0.1:0.6, show_value = true))
"""

# ╔═╡ a21e9cd7-99b5-4cbc-8a7a-889cb0a19f85
begin
	
fig_comparison_freq = Figure( font = noto_sans)
	
	
	# value_stim = 0.3
	end_tc_fig = 3.0
	

	
ax_comparison_freq_1 = fig_comparison_freq[1, 1] = Axis(fig_comparison_freq, title = "Δ Mean firing rate vs ISI for f=0.6", xlabel="ISI", ylabel = "Δ mean firing rate")
	
	
ax_comparison_freq_2 = fig_comparison_freq[2, 1] = Axis(fig_comparison_freq, title = "Δ Mean firing rate vs ISI for f=0.7", xlabel="ISI", ylabel = "Δ mean firing rate")
	
	ax_comparison_freq_3 = fig_comparison_freq[1, 2] = Axis(fig_comparison_freq, title = "Δ Mean firing rate vs ISI for f=0.8", xlabel="ISI", ylabel = "Δ mean firing rate")
	
	
ax_comparison_freq_4 = fig_comparison_freq[2, 2] = Axis(fig_comparison_freq, title = "Δ Mean firing rate vs ISI for f=0.9", xlabel="ISI", ylabel = "Δ mean firing rate")	



	
	tc_list_fig = collect(0.1:τstep:end_tc_fig)
	
	colors_fig = to_colormap(:viridis, length(tc_list))

	subsystem1_fig = zeros(length(tc_list),length(ISI_list))
	subsystem2_fig = zeros(length(tc_list),length(ISI_list))
	
	
	for (index,freq_fig) in enumerate(collect(0.6:0.1:0.9))
	
	
	for j=1:length(tc_list)
		
		tc_int = tc_list[j]
		
		fig1_dict = Dict(
	"f1" => freq_fig,
	"value_stim_f1" => value_stim_2,
	"Tinter" => ISI_list,
	"Tstim" => 0.5,
	"Tfin" => 360.0,
	"value_stim_f2" => 0.0,
    "τ_integrator" => tc_int
)
	
	for (i,dict_temp) in enumerate(dict_list(fig1_dict))
		
		data = load(datadir("sims","script3",savename("oddball_task",dict_temp,"jld2")))

	mean_value = get_mean_firing_rate(data,"microcircuit1-ecell1","current-to-e1",dict_temp["value_stim_f1"])
	mean_value_sys2 = get_mean_firing_rate(data,"microcircuit2-ecell1","current-to-e2",dict_temp["value_stim_f1"])	
			
			
		subsystem1_fig[j,i]= mean_value# - mean_value_sys2
		subsystem2_fig[j,i]= mean_value_sys2
	end
	
		subsystem1_ind = map(!,isnan.(subsystem1_fig[j,:]))
	subsystem2_ind = map(!,isnan.(subsystem2_fig[j,:]))

	fit1 = curve_fit(model, ISI_list[subsystem1_ind], subsystem1_fig[j,subsystem1_ind], p0, lower = lb, upper = ub)
		
	fit2 = curve_fit(model, ISI_list[subsystem2_ind], subsystem2_fig[j,subsystem2_ind], p0, lower = lb, upper = ub)
		
			if index==1
		lines!(ax_comparison_freq_1,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
	
			elseif index ==2
						lines!(ax_comparison_freq_2,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
	
			elseif index ==3
						lines!(ax_comparison_freq_3,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
	
			else
						lines!(ax_comparison_freq_4,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
	
			end
	
	
		end
	end
	
	fig_comparison_freq

end

# ╔═╡ 0ffb5306-a695-433b-805b-b665cb5019a9
md"""
Different points have to be noted:
- The general behavior of the network is consistent across frequencies and strength of stimuli
- There seem to be two time constants in the network. The one of the increase of firing rate of the E-cells with ISI. And the one of the MMN effect.
- Interesting feature of the MMN effect, the effect seems stronger for a specific ISI. It is necessary to check if this is observed in the litterature.
-


Let's now see if the value of the stimulus has a kind of impact on the behavior of these curves.

"""

# ╔═╡ 56faff41-67f3-4100-a28c-ab0b04d8486f
md"""
Frequency value:
$(@bind freq_fig_stim PlutoUI.Slider(0.5:0.1:0.9, show_value = true, default = 0.6))
"""

# ╔═╡ c0c213dd-6a51-41ce-8d92-30f0109ec047
begin
	
fig_comparison_stim = Figure( font = noto_sans)
	
	
	

	
ax_comparison_stim_1 = fig_comparison_stim[1, 1] = Axis(fig_comparison_stim, title = "Δ Mean firing rate vs ISI for stim=0.1", xlabel="ISI", ylabel = "Δ mean firing rate")
	
	
ax_comparison_stim_2 = fig_comparison_stim[2, 1] = Axis(fig_comparison_stim, title = "Δ Mean firing rate vs ISI for stim=0.2", xlabel="ISI", ylabel = "Δ mean firing rate")
	
	ax_comparison_stim_3 = fig_comparison_stim[1, 2] = Axis(fig_comparison_stim, title = "Δ Mean firing rate vs ISI for stim=0.3", xlabel="ISI", ylabel = "Δ mean firing rate")
	
	
ax_comparison_stim_4 = fig_comparison_stim[2, 2] = Axis(fig_comparison_stim, title = "Δ Mean firing rate vs ISI for stim=0.4", xlabel="ISI", ylabel = "Δ mean firing rate")	
ax_comparison_stim_5 = fig_comparison_stim[3, 1] = Axis(fig_comparison_stim, title = "Δ Mean firing rate vs ISI for stim=0.5", xlabel="ISI", ylabel = "Δ mean firing rate")	
ax_comparison_stim_6 = fig_comparison_stim[3, 2] = Axis(fig_comparison_stim, title = "Δ Mean firing rate vs ISI for stim=0.6", xlabel="ISI", ylabel = "Δ mean firing rate")	



	

	subsystem1_fig_stim = zeros(length(tc_list),length(ISI_list))
	subsystem2_fig_stim = zeros(length(tc_list),length(ISI_list))
	
	
	fitΔ = zeros(length(collect(0.1:0.1:0.6)),length(tc_list),length(ISI_list))
	
	
	for (index,stim_fig) in enumerate(collect(0.1:0.1:0.6))
	
	
	for j=1:length(tc_list)
		
		tc_int = tc_list[j]
		
		fig1_dict = Dict(
	"f1" => freq_fig_stim,
	"value_stim_f1" => stim_fig,
	"Tinter" => ISI_list,
	"Tstim" => 0.5,
	"Tfin" => 360.0,
	"value_stim_f2" => 0.0,
    "τ_integrator" => tc_int
)
	
	for (i,dict_temp) in enumerate(dict_list(fig1_dict))
		
		data = load(datadir("sims","script3",savename("oddball_task",dict_temp,"jld2")))

	mean_value = get_mean_firing_rate(data,"microcircuit1-ecell1","current-to-e1",dict_temp["value_stim_f1"])
	mean_value_sys2 = get_mean_firing_rate(data,"microcircuit2-ecell1","current-to-e2",dict_temp["value_stim_f1"])	
			
			
		subsystem1_fig[j,i]= mean_value# - mean_value_sys2
		subsystem2_fig[j,i]= mean_value_sys2
	end
	
		subsystem1_ind = map(!,isnan.(subsystem1_fig[j,:]))
	subsystem2_ind = map(!,isnan.(subsystem2_fig[j,:]))

	fit1 = curve_fit(model, ISI_list[subsystem1_ind], subsystem1_fig[j,subsystem1_ind], p0, lower = lb, upper = ub)
		
	fit2 = curve_fit(model, ISI_list[subsystem2_ind], subsystem2_fig[j,subsystem2_ind], p0, lower = lb, upper = ub)
				fitΔ[index,j,:] = model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2))

			if index==1
		lines!(ax_comparison_stim_1,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
				
				
			elseif index ==2
						lines!(ax_comparison_stim_2,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])

			elseif index ==3
						lines!(ax_comparison_stim_3,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
	
			elseif index == 4
						lines!(ax_comparison_stim_4,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
				
				elseif index == 5
						lines!(ax_comparison_stim_5,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
				
				elseif index == 6
						lines!(ax_comparison_stim_6,ISI_list,model(ISI_list, coef(fit1)) .-model(ISI_list, coef(fit2)),linewidth=5,color=colors[j])
	
			end
	
	
		end
	end
	
	fig_comparison_stim

end

# ╔═╡ f671bc10-befc-427c-8690-f0e2e0677710
md"""
It seems that an increased stimulus flattens the MMN effect. This is however the direct firing rate difference and not the percentage.

**Question: What kind of function is this?**

"""

# ╔═╡ f0a54ab8-baa0-482b-bf63-6728a3e49b5e
md"""
#### Analysis of the firing rates

"""

# ╔═╡ 17e690c9-c56a-463e-a01b-9386678aaf41
md"""
The next plot looks at the time constant of the exponential fit to the variation of mean firing rate with respect to the ISI. 

It's worth noting that this analysis just cnsider the time constant and not the exact value at the end.


The gradient of colors represent the strengths of the stimuli.

"""

# ╔═╡ 4e8f9c4b-c7df-4b43-a3d2-4022d98e6c85
begin
	# cell to plot all the time constants of MMN
	
	  fig_tc = Figure( font = noto_sans)
	
	
	value_stim_list = collect(0.1:0.1:0.6)
	value_freq = collect(0.6:0.1:0.9)
	end_tc2 = 3.0
	tc_list2 = collect(0.1:0.2:end_tc2)

	
	average_list = zeros(length(value_stim_list),length(value_freq))
	average_percent_list = zeros(length(value_stim_list),length(value_freq))
	
	
	
	
    ax_tc = fig_tc[1, 1:2] = Axis(fig_tc, title = "Mean firing rate vs ISI", xlabel="TC integrator", ylabel = "Δ mean firing rate")
	

	Δ_list = zeros(length(value_stim_list),length(value_freq),length(tc_list))

list_map= [:reds, :greens, :blues, :grays, :thermal, :viridis, :grays]
	
			for w=1:length(value_freq)

	for q=1:length(value_stim_list)
		
		value_tc = value_stim_list[q]
		
	
			freq_tc = value_freq[w]
			
	figtc_list = zeros(length(tc_list),length(ISI_list))
	
				colors_list = to_colormap(list_map[w], length(value_stim_list))

	
	for j=1:length(tc_list)
		
		tc_int = tc_list[j]
		
		fig1_dict = Dict(
	"f1" => freq_tc,
	"value_stim_f1" => value_tc,
	"Tinter" => ISI_list,
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
		
						average_list[q,w] += figtc_list[j,1]/length(tc_list) 
				
# 				 figtc_list[j,:] = (-figtc_list[j,:] .+ figtc_list[j,end])./figtc_list[j,end]
		
	

	
		
	
#  	@. model(x, p) = p[3] + p[1]*exp(-x*p[2])
	
#  	p0 = [1.0,0.1,5.0]
	
	fit = curve_fit(model, ISI_list, figtc_list[j,:], p0, lower = lb, upper = ub)
		
		
	Δ_list[q,w,j] = 1.0/coef(fit)[2]
		

		
	end
							average_percent_list[q,w] = mean(figtc_list[:,1]) 

			lines!(ax_tc, tc_list, Δ_list[q,w,:], label="freq =$freq_tc and value_stim = $value_tc", linewidth = 3, color = colors_list[q])
			
		end
	end
		fig_tc[1, 3] = Legend(fig_tc, ax_tc, "Legend", framevisible = false)

	# fig_tc
end

# ╔═╡ 09d8d74f-a8e2-494b-b6e6-06694426666c
begin
	 fig_color = Figure( font = noto_sans)
     ax_color1 = Axis(fig_color, xlabel = "τ integrator", ylabel = "τ firing rates S1", title = "Freq of 0.6 ")
     blues_list = to_colormap(:blues, length(value_stim_list))
	for i=1:length(value_stim_list)
     lines!(ax_color1, tc_list, Δ_list[i,1,:], color= blues_list[i], linewidth = 3)
	end
	fig_color[1, 1] = ax_color1
	
	 ax_color2 = Axis(fig_color,xlabel = "τ integrator", ylabel = "τ firing rates S1", title = "Freq of 0.7 ")
     reds_list = to_colormap(:reds, length(value_stim_list))
	for i=1:length(value_stim_list)
     lines!(ax_color2, tc_list, Δ_list[i,2,:], color= reds_list[i], linewidth = 3)
	end
	fig_color[1, 2] = ax_color2
	
	
	 ax_color3 = Axis(fig_color,xlabel = "τ integrator", ylabel = "τ firing rates S1", title = "Freq of 0.8 ")
     greens_list = to_colormap(:greens, length(value_stim_list))
	for i=1:length(value_stim_list)
     lines!(ax_color3, tc_list, Δ_list[i,3,:], color= greens_list[i], linewidth = 3)
	end
	fig_color[2, 1] = ax_color3
	
	
	 ax_color4 = Axis(fig_color,xlabel = "τ integrator", ylabel = "τ firing rates S1", title = "Freq of 0.9 ")
     grays_list = to_colormap(:grays, length(value_stim_list))
	for i=1:length(value_stim_list)
		lines!(ax_color4, tc_list, Δ_list[i,4,:], color= grays_list[end-i+1], linewidth = 3)
	end
	fig_color[2, 2] = ax_color4
	
	  # cbar_fig_color = Colorbar(fig_color, value_stim_list, label = "color value", ticklabelsize = 14, labelpadding = 5, width = 10)
	
	
	fig_color
	
	
end

# ╔═╡ f1421b5d-e88e-49fb-82d7-bf7289fd4207
md"""
**Heatmap of the effect**

"""

# ╔═╡ be96e4fa-65cb-4fad-8c9d-a71fa99adaef
begin
	 fig = Figure( font = noto_sans)
     ax = Axis(fig, title = "Heatmap for stim = 0.1 ", xlabel="Frequencies of stim", ylabel="τ integrator")
     hmap = heatmap!(ax,value_freq, tc_list, Δ_list[1,:,:], colormap = :blues)
	 cbar = Colorbar(fig, hmap, label = "τ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig[1, 1] = ax
    fig[1, 2] = cbar
	
	ax2 = Axis(fig, title = "Heatmap for stim = 0.2 ")
     hmap2 = heatmap!(ax2,value_freq, tc_list, Δ_list[2,:,:], colormap = :blues)
	 cbar2 = Colorbar(fig, hmap2, label = "τ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig[1, 3] = ax2
    fig[1, 4] = cbar2
	
	ax3 = Axis(fig, title = "Heatmap for stim = 0.3 ")
     hmap3 = heatmap!(ax3,value_freq, tc_list, Δ_list[3,:,:], colormap = :blues)
	 cbar3 = Colorbar(fig, hmap3, label = "τ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig[2, 1] = ax3
    fig[2, 2] = cbar3
	
	ax4 = Axis(fig, title = "Heatmap for stim = 0.4 ")
     hmap4 = heatmap!(ax4,value_freq, tc_list, Δ_list[4,:,:], colormap = :blues)
	 cbar4 = Colorbar(fig, hmap4, label = "τ values", width = 15,
                ticksize=15, tickalign = 1, height = Relative(3.55/4))
	fig[2, 3] = ax4
    fig[2, 4] = cbar4
	
	fig
end

# ╔═╡ d6f5b0d0-03c6-4aca-bb45-aaa95b06a035


# ╔═╡ f45547ef-6189-4013-89b6-0ce1c9d09cd7


# ╔═╡ 5780ab0b-08ca-45ac-a613-fbfd2504f1aa


# ╔═╡ d6b1c61d-6e95-4ad2-9492-c6e0d28a76f8


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

# ╔═╡ 1aa2b36a-3504-4b50-a850-001b7b23c7f0
md""" ### Analysis of the Δ of firing rate

The question is *how to fit the difference of mean firing rates between subnetworks?*

Two ideas come to mind:
- An inversed mexican hat
- A sum of exp andtanh functions

At first look both these functions would give access to a time constant and allow for a detailed comparison between parameters.


"""

# ╔═╡ 239ba79c-8f87-4d95-8e1d-2d05124168a0
begin
	
	figtemp = Figure()
	axtemp = Axis(figtemp)
	
	@. hat_fit(x,p) = p[1]*exp(-(x - p[3])^2/p[4]) + p[2]#*x + p[5]#*exp(-x*p[5]) + p[6]
	
	@. exp_tanh_fit(x,p) = p[1] + p[2]*exp(-x*p[3]) + p[4]*tanh((x-p[6])*p[5])
	
	
	pat = [-1.0,-1.0,1.0,2.0,1.0,1.0]
	pexp = [-10.0,1.0,1.0,10.0,10.0,4.0]
	
	
	lb_hat = [-10.0,-10.0,0.0]
	ub_hat = []
	
	
	lb_exptanh = [-19.5,0.0,0.01,0.0,0.01,-000.0]
	ub_exptanh = [20.0,10.0,50.0,40.0,10.0,10.0]
	
	# TODO better code to have all the data before fitting
	
	temp_data = fitΔ[3,10,:]
	
		fit_hat = curve_fit(hat_fit, ISI_list, temp_data, pat)
	
		fit_exptanh = curve_fit(exp_tanh_fit, ISI_list, temp_data, pexp)

	lines!(axtemp,ISI_list, hat_fit(ISI_list, coef(fit_hat)), color=:red)
	lines!(axtemp,ISI_list, exp_tanh_fit(ISI_list, coef(fit_exptanh)), linewidth = 4)
	lines!(axtemp,ISI_list, temp_data,linewidth = 3)

	
	figtemp[1,1] = axtemp
	figtemp
	
	
	
	
end

# ╔═╡ f121b6ea-f4b3-4bff-9a9c-7b62e214802c
coef(fit_exptanh)

# ╔═╡ 47320f6f-111d-4e3b-acd1-f1613c3ad03e
coef(fit_hat)

# ╔═╡ 5b1d3f55-e152-4226-b50f-109b650260b9
begin
	
end

# ╔═╡ Cell order:
# ╠═33aa235e-fec6-11eb-2e65-9b8f1169abae
# ╠═63c5c372-10ac-401e-81bf-70c929e3b313
# ╠═9b053399-c4f2-46b5-9890-f1b538c2bfac
# ╠═4607fe2b-ebe7-4c1a-9f50-302aa0254987
# ╠═6bac87c9-84cc-441b-80e0-cfa5ddddd874
# ╟─af717f0e-2c0b-4691-87c9-9e13c8c097ee
# ╟─e975c7bd-ed9c-4fbc-bc88-b4fb66a846cf
# ╟─92677b31-eb38-437f-bbb8-8d12d7840adb
# ╟─864bc534-db2d-4bb7-962e-d348656e0adf
# ╟─b9e8a257-6b81-4c5b-9c43-bae728885e55
# ╟─5098df57-da12-4d3b-9452-0080994b3c18
# ╟─c14c58d2-ff7f-434c-ab4e-64df6dd1b020
# ╟─a21e9cd7-99b5-4cbc-8a7a-889cb0a19f85
# ╟─0ffb5306-a695-433b-805b-b665cb5019a9
# ╟─56faff41-67f3-4100-a28c-ab0b04d8486f
# ╟─c0c213dd-6a51-41ce-8d92-30f0109ec047
# ╟─f671bc10-befc-427c-8690-f0e2e0677710
# ╟─f0a54ab8-baa0-482b-bf63-6728a3e49b5e
# ╟─17e690c9-c56a-463e-a01b-9386678aaf41
# ╟─4e8f9c4b-c7df-4b43-a3d2-4022d98e6c85
# ╟─09d8d74f-a8e2-494b-b6e6-06694426666c
# ╟─f1421b5d-e88e-49fb-82d7-bf7289fd4207
# ╟─be96e4fa-65cb-4fad-8c9d-a71fa99adaef
# ╠═d6f5b0d0-03c6-4aca-bb45-aaa95b06a035
# ╠═f45547ef-6189-4013-89b6-0ce1c9d09cd7
# ╠═5780ab0b-08ca-45ac-a613-fbfd2504f1aa
# ╠═d6b1c61d-6e95-4ad2-9492-c6e0d28a76f8
# ╠═005e59fe-4f25-4f3e-8712-e8a5c761fb64
# ╠═63a0d410-1c90-413a-9d09-e800c71eff55
# ╠═4a53013c-91ca-453f-9b0d-26096ccf822a
# ╠═1aa2b36a-3504-4b50-a850-001b7b23c7f0
# ╠═239ba79c-8f87-4d95-8e1d-2d05124168a0
# ╠═f121b6ea-f4b3-4bff-9a9c-7b62e214802c
# ╠═47320f6f-111d-4e3b-acd1-f1613c3ad03e
# ╠═5b1d3f55-e152-4226-b50f-109b650260b9
