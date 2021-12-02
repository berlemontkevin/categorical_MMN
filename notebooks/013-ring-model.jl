### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ a12d8f60-3113-11ec-3810-95a5aa96c35e
begin
		using DrWatson
	
	quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())
end

# ╔═╡ 98dd4b78-6e43-4ce7-a599-d7feda16b334
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

# ╔═╡ 517913f8-0766-4001-86b0-9a6cbf8ab7d9
using BenchmarkTools

# ╔═╡ c040262f-ad0e-450c-8b3f-56a7fa3ad0f6
md""" #Plotting some dynamics

"""

# ╔═╡ 9a988fd9-e261-4439-9c56-417a5b390589
	function plot_bump_indice(test_layer,indice)
	    fig = Figure()
	temp_fr = zeros(6,128)
	    ax = fig[1, 1] = Axis(fig, title = "Firing rates")

	for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].soma.r_save[indice]
		temp_fr[2,i] = test_layer.list_microcircuit[i].pv.r_save[indice]
		temp_fr[3,i] = test_layer.list_microcircuit[i].sst.r_save[indice]
		temp_fr[4,i] = test_layer.list_microcircuit[i].vip.r_save[indice]
		temp_fr[6,i] = test_layer.list_microcircuit[i].integrator.r_save[indice]
		temp_fr[5,i] = test_layer.list_microcircuit[i].ngfc.r_save[indice]
		
	end
	
	lines!(ax,temp_fr[1,:],linewidth=4,color=(:blue,0.5),label="E")
	lines!(ax,temp_fr[2,:],linewidth=4,color=(:orange,0.5),label="pv")

	lines!(ax,temp_fr[3,:],linewidth=4,color=(:red,0.5),label="sst")

	lines!(ax,temp_fr[4,:],linewidth=4,color=(:green,0.5),label="vip")
	lines!(ax,temp_fr[5,:],linewidth=4,color=(:black,0.5),label="ngfc")

	lines!(ax,temp_fr[6,:],linewidth=4,color=(:purple,0.5),label="integrator")
	ylims!(ax,0,55)
	fig[1, 2] = Legend(fig, ax)
		
		
	

return fig
		
	end

# ╔═╡ c163f2e4-5ef8-4849-a0ae-6ce0749e6c3b
string(projectdir("notebooks","animation"),"test_functions",".mp4")

# ╔═╡ f51db6b0-61e2-42b5-8682-9448ac5877ee
LocalResource(projectdir("notebooks","animation","test_functions.mp4"))

# ╔═╡ d2adecdc-9ecf-4570-9ed8-5c5e8313aeea
md""" ## Benchmarking the ring model


"""

# ╔═╡ 4c680a9d-1211-451f-94dc-333c7bc652b1
 #microcircuit{soma_PC,dend_sigmoid}(name="microcircuit1")

# ╔═╡ b0f8ec9b-fe26-413e-926c-583324500446
with_terminal() do
	
	@time create_layer_bump(parameters_bump_attractor(),128; param_microcircuit = parameters_microcircuit())
	
	
end

# ╔═╡ 7fc2e01f-64fd-4307-b3fe-06be97af5ef2
begin
	# TODO function that generates currents for all the cells
	 sim = simulation_parameters()
	 ring_simu = Dict{String,Vector{Float64}}()
	 bump = create_layer_bump(parameters_bump_attractor(),128, param_microcircuit = parameters_microcircuit(time_tot = 1000))
	stim = 50.0
	
	for i=1:128
		
		sim.current[string("microcircuit$i","-","ecell1")] = 0.2.*ones(1000).*orientation_kernel(stim,bump.list_microcircuit[i].soma.preferred_stim, bump.bump_param)
		sim.current[string("microcircuit$i","-","dend1")] = 0.0.*ones(1000)
		sim.current[string("microcircuit$i","-","vipcell1")] = 0.0.*ones(1000)
		sim.current[string("microcircuit$i","-","sstcell1")] = 0.0.*ones(1000)
		sim.current[string("microcircuit$i","-","pvcell1")] = 0.0.*ones(1000)
		sim.current[string("microcircuit$i","-","ngfccell1")] = 0.0.*ones(1000)

		sim.current[string("microcircuit$i","-","integrator1")] = 0.0.*ones(1000)

		
	end
end

# ╔═╡ fb06e280-c95a-4a59-99ff-6806c8410f59
	fig_1000 = plot_bump_indice(bump,4000)


# ╔═╡ 312cba43-d235-4cdb-8225-d81a437c199f
get_firing_rate(bump,1000, "ngfc")

# ╔═╡ 2c07be88-f460-4cf9-b54c-84eb9d5877b4
lines(bump.list_microcircuit[20].ngfc.r_save)

# ╔═╡ 9575c9d3-4435-4f38-b351-4ccf076f9927
begin
	function load_bump_E(time;test_layer = bump, sim = sim)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].soma.r_save[floor(Int,time[])]
		

	end
		return temp_fr[1,:]
	end
	
	
	function load_bump_pv(time;test_layer = bump, sim = sim)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[2,i] = test_layer.list_microcircuit[i].pv.r_save[floor(Int,time[])]

	end
		return temp_fr[2,:]
	end
	
		function load_bump_sst(time;test_layer = bump, sim = sim)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[3,i] = test_layer.list_microcircuit[i].sst.r_save[floor(Int,time[])]

	end
		return temp_fr[3,:]
	end
	
		function load_bump_vip(time;test_layer = bump, sim = sim)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[4,i] = test_layer.list_microcircuit[i].vip.r_save[floor(Int,time[])]

	end
		return temp_fr[4,:]
	end
	
		function load_bump_integrator(time;test_layer = bump, sim = sim)
		temp_fr = zeros(5,128)

	for i=1:128
		
		temp_fr[5,i] = test_layer.list_microcircuit[i].integrator.r_save[floor(Int,time[])]

	end
		return temp_fr[5,:]
	end
	
	
		function load_bump_current(time;test_layer = bump, sim = sim)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[1,i] = sim.current[string("microcircuit$i","-","ecell1")][floor(Int,time[])]
	

	end
		return temp_fr[1,:]
	end
end

# ╔═╡ 487ed0d5-07cb-47f5-b051-05f8d662903d
begin
	time = Node(1.0)

	
	
xs = range(1, 128, length=128)

ys_E = @lift(load_bump_E.($time))
ys_pv = @lift(load_bump_pv.($time))
ys_sst = @lift(load_bump_sst.($time))
ys_vip = @lift(load_bump_vip.($time))
ys_integrator = @lift(load_bump_integrator.($time))
#ys_current = @lift(load_bump_current.($time))

fig_animation = Figure()
axis_animation = Axis(fig_animation[1,1],title = @lift("t = $(round(0.5.*$time, digits = 1)) ms"))
#xis_animation2 = Axis(fig_animation[2,1], yticklabelcolor = :red, yaxisposition = :right)

ylims!(axis_animation,0,60)
lines!(axis_animation,xs, ys_E, color = :blue, linewidth = 4)
lines!(axis_animation,xs, ys_pv, color = :orange, linewidth = 4)
lines!(axis_animation,xs, ys_sst, color = :red, linewidth = 4)
lines!(axis_animation,xs, ys_vip, color = :green, linewidth = 4)
lines!(axis_animation,xs, ys_integrator, color = :purple, linewidth = 4)

	
	#lines!(axis_animation2,xs, ys_current, color = :black, linewidth = 4)
#ylims!(axis_animation2,0,0.5)

# scatter!(xs, ys_2, color = :red, markersize = 15)

framerate = 50
t_tot = length(bump.list_microcircuit[1].soma.r_save)-1
timestamps = range(1, t_tot, step=2)
# todo the real framerate
	
	
record(fig_animation, "animation\\bump_animation.mp4", timestamps;
        framerate = framerate) do t
    time[] = floor(Int,t)
end
	

LocalResource(projectdir("notebooks","animation","bump_animation.mp4"))
	
	
end

# ╔═╡ faec8df9-5f1f-4254-8c45-17dd3ce778c8
bump_animation(bump,2,projectdir("notebooks","animation\\"),"test_functions",".mp4")

# ╔═╡ c26e8deb-b971-4610-b9fb-1bb873ccb9ee
with_terminal() do
	
	@time full_time_dynamics!(bump, sim)
	
end

# ╔═╡ b1446125-06ce-4cf6-9661-53000d02dfc4
with_terminal() do
	
	@code_warntype full_time_dynamics!(bump, sim)
	
end

# ╔═╡ d4d3eb7e-6bbf-46f7-af66-87ada476c2f7
with_terminal() do
	
	@code_warntype time_step!(bump.list_microcircuit, sim,10)
	
end

# ╔═╡ 4f272847-bf74-436b-a63b-99766bea496c
@benchmark time_step!(bump.list_microcircuit, sim, 100)

# ╔═╡ b60cd5ba-3d80-410e-96a7-a4d0874bcf69
@benchmark full_time_dynamics!(bump, sim)

# ╔═╡ c34f1c5f-6175-4b7c-a8b1-118100160eee
@benchmark synapse_derivative!(bump.list_microcircuit[1].soma, bump.eq_diff_method)

# ╔═╡ aea0f11c-9fae-4e62-8c9f-c5c8c9e27ce8
with_terminal() do
	
	@code_warntype synapse_derivative!(bump.list_microcircuit[1].soma, bump.eq_diff_method)
	
end

# ╔═╡ Cell order:
# ╠═a12d8f60-3113-11ec-3810-95a5aa96c35e
# ╠═98dd4b78-6e43-4ce7-a599-d7feda16b334
# ╠═c040262f-ad0e-450c-8b3f-56a7fa3ad0f6
# ╠═9a988fd9-e261-4439-9c56-417a5b390589
# ╠═fb06e280-c95a-4a59-99ff-6806c8410f59
# ╠═312cba43-d235-4cdb-8225-d81a437c199f
# ╠═2c07be88-f460-4cf9-b54c-84eb9d5877b4
# ╠═9575c9d3-4435-4f38-b351-4ccf076f9927
# ╠═487ed0d5-07cb-47f5-b051-05f8d662903d
# ╠═c163f2e4-5ef8-4849-a0ae-6ce0749e6c3b
# ╠═faec8df9-5f1f-4254-8c45-17dd3ce778c8
# ╠═f51db6b0-61e2-42b5-8682-9448ac5877ee
# ╠═d2adecdc-9ecf-4570-9ed8-5c5e8313aeea
# ╠═517913f8-0766-4001-86b0-9a6cbf8ab7d9
# ╠═4c680a9d-1211-451f-94dc-333c7bc652b1
# ╠═b0f8ec9b-fe26-413e-926c-583324500446
# ╠═7fc2e01f-64fd-4307-b3fe-06be97af5ef2
# ╠═c26e8deb-b971-4610-b9fb-1bb873ccb9ee
# ╠═b1446125-06ce-4cf6-9661-53000d02dfc4
# ╠═d4d3eb7e-6bbf-46f7-af66-87ada476c2f7
# ╠═4f272847-bf74-436b-a63b-99766bea496c
# ╠═b60cd5ba-3d80-410e-96a7-a4d0874bcf69
# ╠═c34f1c5f-6175-4b7c-a8b1-118100160eee
# ╠═aea0f11c-9fae-4e62-8c9f-c5c8c9e27ce8
