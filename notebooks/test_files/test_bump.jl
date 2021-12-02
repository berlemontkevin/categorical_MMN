### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ a9602c60-260f-11ec-26fa-6b4383161e4e
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())
end

# ╔═╡ 52c79c6b-69b0-4fa8-a299-eabe688a889c
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

# ╔═╡ 2027b524-64e9-4015-8231-cf4a5689628a
begin
	 ring_simu = Dict{String,Vector{Float64}}()
	test_layer = create_layer_bump(bump_structure(),time_tot=18000, noise = false)
	stim = 50.0
	
	for i=1:128
		
		ring_simu[string("microcircuit$i","-","ecell1")] = [zeros(5000);0.2*ones(1000)*orientation_kernel(stim,test_layer.list_microcircuit[i].list_soma[1].preferred_stim, test_layer.bump_param);zeros(2000);0.2*ones(1000)*orientation_kernel(stim,test_layer.list_microcircuit[i].list_soma[1].preferred_stim, test_layer.bump_param);zeros(2000);0.2*ones(1000)*orientation_kernel(stim,test_layer.list_microcircuit[i].list_soma[1].preferred_stim, test_layer.bump_param);zeros(2000);0.2*ones(1000)*orientation_kernel(50.0,test_layer.list_microcircuit[i].list_soma[1].preferred_stim, test_layer.bump_param);zeros(2000);0.2*ones(1000)*orientation_kernel(200.0,test_layer.list_microcircuit[i].list_soma[1].preferred_stim, test_layer.bump_param)]
	end
	
	
end

# ╔═╡ c84c3027-ba6e-4b5f-813a-facc867371d0
test_layer.list_microcircuit[10].list_sst[1]

# ╔═╡ f059dba5-10c7-4635-bc2a-26af524816a4
test_layer.list_microcircuit[10].list_sst[1].list_syn_post_nmda

# ╔═╡ e106b40b-9328-42f7-93e8-88accbb48d9a
begin
temp_list = Float64[]

for i=1:127
		push!(temp_list,test_layer.list_microcircuit[50].list_sst[1].list_syn_post_nmda[i+2].g)
		
end
	
end

# ╔═╡ 543156de-5ac3-4401-bdb9-fce040aedb69
lines(temp_list)

# ╔═╡ 55da7210-97b5-4f46-9e9c-5aeea5d9be1e
temp_list

# ╔═╡ d78fc3c4-79ca-46e7-ad1f-8aadeacd5cbd
ring_simu

# ╔═╡ a5926cde-47d4-4250-b4da-093e8a01d261
begin
	temp_current = zeros(128)

	for i=1:128
		temp_current[i] =ring_simu[string("microcircuit$i","-","ecell1")][end]
		
	end

	lines(temp_current)
	
end

# ╔═╡ faecca48-f8dc-4495-95ea-1041cc319e49
 sim = simulation_parameters()

# ╔═╡ 053795ec-8df7-45af-b6c0-a2c896f54bca
full_time_dynamics(test_layer, sim, ring_simu)

# ╔═╡ 5c4e252a-49c4-4296-8bad-d768b4dedf9d
lines(test_layer.list_microcircuit[20].list_pv[1].r_save)

# ╔═╡ 5adcf2e2-2679-4985-91fa-35aaa0b0d55a
begin
	
	function plot_bump_indice(indice)
	    fig = Figure()
	temp_fr = zeros(5,128)
	    ax = fig[1, 1] = Axis(fig, title = "Firing rates")

	for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].list_soma[1].r_save[indice]
		temp_fr[2,i] = test_layer.list_microcircuit[i].list_pv[1].r_save[indice]
		temp_fr[3,i] = test_layer.list_microcircuit[i].list_sst[1].r_save[indice]
		temp_fr[4,i] = test_layer.list_microcircuit[i].list_vip[1].r_save[indice]
		temp_fr[5,i] = test_layer.list_microcircuit[i].list_integrator[1].r_save[indice]

	end
	
	lines!(ax,temp_fr[1,:],linewidth=4,color=(:blue,0.5),linestyle=:dash,label="E")
	lines!(ax,temp_fr[2,:],linewidth=4,color=(:orange,0.5),linestyle=:dash,label="pv")

	lines!(ax,temp_fr[3,:],linewidth=4,color=(:red,0.5),linestyle=:dash,label="sst")

	lines!(ax,temp_fr[4,:],linewidth=4,color=(:green,0.5),linestyle=:dash,label="vip")
	lines!(ax,temp_fr[5,:],linewidth=4,color=(:purple,0.5),linestyle=:dash,label="integrator")
	ylims!(ax,0,55)
	fig[1, 2] = Legend(fig, ax)
		
		
		for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].list_soma[1].r_save[17050]
		temp_fr[2,i] = test_layer.list_microcircuit[i].list_pv[1].r_save[17050]
		temp_fr[3,i] = test_layer.list_microcircuit[i].list_sst[1].r_save[17050]
		temp_fr[4,i] = test_layer.list_microcircuit[i].list_vip[1].r_save[17050]
		temp_fr[5,i] = test_layer.list_microcircuit[i].list_integrator[1].r_save[17050]

	end
			lines!(ax,temp_fr[1,:],linewidth=4,color=:blue,label="E")
	lines!(ax,temp_fr[2,:],linewidth=4,color=:orange,label="pv")

	lines!(ax,temp_fr[3,:],linewidth=4,color=:red,label="sst")

	lines!(ax,temp_fr[4,:],linewidth=4,color=:green,label="vip")
	lines!(ax,temp_fr[5,:],linewidth=4,color=:purple,label="integrator")

return fig
		
	end
	
	
	fig_5050 = plot_bump_indice(5050)
	
end

# ╔═╡ b64c68e9-4472-4d46-a579-a3e722d3abc1
fig_8050 = plot_bump_indice(8050)

# ╔═╡ fd49b0e7-d91b-4f54-908e-f8a643873bb8
fig_11050 = plot_bump_indice(11050)

# ╔═╡ 299027cd-9655-4005-a2c1-ec229221b8a2
fig_14050 = plot_bump_indice(14050)

# ╔═╡ 8e1bff84-19fb-416c-b8fd-dd02fcd36af4
begin
	function load_bump_E(time)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].list_soma[1].r_save[floor(Int,time[])]
		temp_fr[2,i] = test_layer.list_microcircuit[i].list_pv[1].r_save[floor(Int,time[])]
		temp_fr[3,i] = test_layer.list_microcircuit[i].list_sst[1].r_save[floor(Int,time[])]
		temp_fr[4,i] = test_layer.list_microcircuit[i].list_vip[1].r_save[floor(Int,time[])]
		temp_fr[5,i] = test_layer.list_microcircuit[i].list_integrator[1].r_save[floor(Int,time[])]

	end
		return temp_fr[1,:]
	end
	
	
	function load_bump_pv(time)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].list_soma[1].r_save[floor(Int,time[])]
		temp_fr[2,i] = test_layer.list_microcircuit[i].list_pv[1].r_save[floor(Int,time[])]
		temp_fr[3,i] = test_layer.list_microcircuit[i].list_sst[1].r_save[floor(Int,time[])]
		temp_fr[4,i] = test_layer.list_microcircuit[i].list_vip[1].r_save[floor(Int,time[])]
		temp_fr[5,i] = test_layer.list_microcircuit[i].list_integrator[1].r_save[floor(Int,time[])]

	end
		return temp_fr[2,:]
	end
	
		function load_bump_sst(time)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].list_soma[1].r_save[floor(Int,time[])]
		temp_fr[2,i] = test_layer.list_microcircuit[i].list_pv[1].r_save[floor(Int,time[])]
		temp_fr[3,i] = test_layer.list_microcircuit[i].list_sst[1].r_save[floor(Int,time[])]
		temp_fr[4,i] = test_layer.list_microcircuit[i].list_vip[1].r_save[floor(Int,time[])]
		temp_fr[5,i] = test_layer.list_microcircuit[i].list_integrator[1].r_save[floor(Int,time[])]

	end
		return temp_fr[3,:]
	end
	
		function load_bump_vip(time)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].list_soma[1].r_save[floor(Int,time[])]
		temp_fr[2,i] = test_layer.list_microcircuit[i].list_pv[1].r_save[floor(Int,time[])]
		temp_fr[3,i] = test_layer.list_microcircuit[i].list_sst[1].r_save[floor(Int,time[])]
		temp_fr[4,i] = test_layer.list_microcircuit[i].list_vip[1].r_save[floor(Int,time[])]
		temp_fr[5,i] = test_layer.list_microcircuit[i].list_integrator[1].r_save[floor(Int,time[])]

	end
		return temp_fr[4,:]
	end
	
		function load_bump_integrator(time)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[1,i] = test_layer.list_microcircuit[i].list_soma[1].r_save[floor(Int,time[])]
		temp_fr[2,i] = test_layer.list_microcircuit[i].list_pv[1].r_save[floor(Int,time[])]
		temp_fr[3,i] = test_layer.list_microcircuit[i].list_sst[1].r_save[floor(Int,time[])]
		temp_fr[4,i] = test_layer.list_microcircuit[i].list_vip[1].r_save[floor(Int,time[])]
		temp_fr[5,i] = test_layer.list_microcircuit[i].list_integrator[1].r_save[floor(Int,time[])]

	end
		return temp_fr[5,:]
	end
	
	
		function load_bump_current(time)
		temp_fr = zeros(5,128)

	for i=1:128
		temp_fr[1,i] = ring_simu[string("microcircuit$i","-","ecell1")][floor(Int,time[])]
		temp_fr[2,i] = test_layer.list_microcircuit[i].list_pv[1].r_save[floor(Int,time[])]
		temp_fr[3,i] = test_layer.list_microcircuit[i].list_sst[1].r_save[floor(Int,time[])]
		temp_fr[4,i] = test_layer.list_microcircuit[i].list_vip[1].r_save[floor(Int,time[])]
		temp_fr[5,i] = test_layer.list_microcircuit[i].list_integrator[1].r_save[floor(Int,time[])]

	end
		return temp_fr[1,:]
	end
	
end

# ╔═╡ 92f8283f-9038-48e2-a953-a239caf49b71
begin
	time = Node(4500.0)

	
	
m1 = maximum(test_layer.list_microcircuit[18].list_soma[1].r_save[5050:6000])
m2 = maximum(test_layer.list_microcircuit[18].list_soma[1].r_save[8050:9000])
m3 = maximum(test_layer.list_microcircuit[18].list_soma[1].r_save[11050:12000])

xs = range(1, 128, length=128)

ys_E = @lift(load_bump_E.($time))
ys_pv = @lift(load_bump_pv.($time))
ys_sst = @lift(load_bump_sst.($time))
ys_vip = @lift(load_bump_vip.($time))
ys_integrator = @lift(load_bump_integrator.($time))
ys_current = @lift(load_bump_current.($time))

	# ys_2 = @lift(cos.(xs .- $time) .+ 3)
fig_animation = Figure()
axis_animation = Axis(fig_animation[1,1],title = @lift("t = $(round(0.5.*$time, digits = 1)) ms"))
#xis_animation2 = Axis(fig_animation[2,1], yticklabelcolor = :red, yaxisposition = :right)

ylims!(axis_animation,0,60)
lines!(axis_animation,xs, ys_E, color = :blue, linewidth = 4)
lines!(axis_animation,xs, ys_pv, color = :orange, linewidth = 4)
lines!(axis_animation,xs, ys_sst, color = :red, linewidth = 4)
lines!(axis_animation,xs, ys_vip, color = :green, linewidth = 4)
lines!(axis_animation,xs, ys_integrator, color = :purple, linewidth = 4)
hlines!(axis_animation,m1)
hlines!(axis_animation,m2)
hlines!(axis_animation,m3)

	
	#lines!(axis_animation2,xs, ys_current, color = :black, linewidth = 4)
#ylims!(axis_animation2,0,0.5)

# scatter!(xs, ys_2, color = :red, markersize = 15)

framerate = 50
t_tot = length(test_layer.list_microcircuit[1].list_soma[1].r_save)
timestamps = range(4500, t_tot, step=4)
# todo the real framerate
	
	
record(fig_animation, "animation\\bump_animation.mp4", timestamps;
        framerate = framerate) do t
    time[] = floor(Int,t)
end
	

LocalResource(projectdir("notebooks","test_files","animation","bump_animation.mp4"))
	
	
end

# ╔═╡ 6850c6dc-f65c-4d00-a4e1-c3381b338112

lines(test_layer.list_microcircuit[18].list_soma[1].r_save)


# ╔═╡ 1dd228c1-86cf-41f7-b1b2-281a2e3aa97e
lines(test_layer.list_microcircuit[23].list_sst[1].r_save)


# ╔═╡ 4bcff9f8-141a-4a68-81a6-3f9aa7444074


# ╔═╡ ff3acfae-6051-439a-b69e-b657ece1c79c


# ╔═╡ Cell order:
# ╠═a9602c60-260f-11ec-26fa-6b4383161e4e
# ╠═52c79c6b-69b0-4fa8-a299-eabe688a889c
# ╠═c84c3027-ba6e-4b5f-813a-facc867371d0
# ╠═f059dba5-10c7-4635-bc2a-26af524816a4
# ╠═e106b40b-9328-42f7-93e8-88accbb48d9a
# ╠═543156de-5ac3-4401-bdb9-fce040aedb69
# ╠═55da7210-97b5-4f46-9e9c-5aeea5d9be1e
# ╠═2027b524-64e9-4015-8231-cf4a5689628a
# ╠═d78fc3c4-79ca-46e7-ad1f-8aadeacd5cbd
# ╠═a5926cde-47d4-4250-b4da-093e8a01d261
# ╠═faecca48-f8dc-4495-95ea-1041cc319e49
# ╠═053795ec-8df7-45af-b6c0-a2c896f54bca
# ╠═5c4e252a-49c4-4296-8bad-d768b4dedf9d
# ╠═5adcf2e2-2679-4985-91fa-35aaa0b0d55a
# ╠═b64c68e9-4472-4d46-a579-a3e722d3abc1
# ╠═fd49b0e7-d91b-4f54-908e-f8a643873bb8
# ╠═299027cd-9655-4005-a2c1-ec229221b8a2
# ╠═8e1bff84-19fb-416c-b8fd-dd02fcd36af4
# ╠═92f8283f-9038-48e2-a953-a239caf49b71
# ╠═6850c6dc-f65c-4d00-a4e1-c3381b338112
# ╠═1dd228c1-86cf-41f7-b1b2-281a2e3aa97e
# ╠═4bcff9f8-141a-4a68-81a6-3f9aa7444074
# ╠═ff3acfae-6051-439a-b69e-b657ece1c79c
