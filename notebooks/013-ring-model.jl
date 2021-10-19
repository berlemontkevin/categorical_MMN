### A Pluto.jl notebook ###
# v0.16.1

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

# ╔═╡ d2adecdc-9ecf-4570-9ed8-5c5e8313aeea
md""" ## Trying to benchmark the ring model


"""

# ╔═╡ 4c680a9d-1211-451f-94dc-333c7bc652b1
 microcircuit{soma_PC,dend_sigmoid}(name="microcircuit1")

# ╔═╡ b0f8ec9b-fe26-413e-926c-583324500446
with_terminal() do
	
	@time create_layer_bump(bump_structure(); param_microcircuit = parameters_microcircuit())
	
	
end

# ╔═╡ 7fc2e01f-64fd-4307-b3fe-06be97af5ef2
begin
	
	 sim = simulation_parameters()
	 ring_simu = Dict{String,Vector{Float64}}()
	 bump = create_layer_bump(bump_structure(), param_microcircuit = parameters_microcircuit(time_tot = 1000))
	stim = 50.0
	
	for i=1:128
		
		ring_simu[string("microcircuit$i","-","ecell1")] = 0.2.*ones(1000).*orientation_kernel(stim,bump.list_microcircuit[i].list_soma[1].preferred_stim, bump.bump_param)
	end
end

# ╔═╡ c26e8deb-b971-4610-b9fb-1bb873ccb9ee
with_terminal() do
	
	@time full_time_dynamics(bump, sim, ring_simu)
	
end

# ╔═╡ b1446125-06ce-4cf6-9661-53000d02dfc4
with_terminal() do
	
	@code_warntype full_time_dynamics(bump, sim, ring_simu)
	
end

# ╔═╡ d4d3eb7e-6bbf-46f7-af66-87ada476c2f7
with_terminal() do
	
	@code_warntype time_step(bump.list_microcircuit, sim, ring_simu,10)
	
end

# ╔═╡ 4f272847-bf74-436b-a63b-99766bea496c
@benchmark time_step(bump.list_microcircuit, sim, ring_simu,100)

# ╔═╡ b60cd5ba-3d80-410e-96a7-a4d0874bcf69
@benchmark full_time_dynamics(bump, sim, ring_simu)

# ╔═╡ a2d81277-3e2b-485c-bfc5-c2749a073f09
with_terminal() do
	
	@time time_step(bump.list_microcircuit, sim, ring_simu,100)
	
end

# ╔═╡ c34f1c5f-6175-4b7c-a8b1-118100160eee
@benchmark synapse_derivative(bump.list_microcircuit[1].list_soma[1])

# ╔═╡ 239afd2b-20a3-4b62-a203-fd27fc232a48
@benchmark synapse_derivative(bump.list_microcircuit[1])

# ╔═╡ aea0f11c-9fae-4e62-8c9f-c5c8c9e27ce8
with_terminal() do
	
	@code_warntype synapse_derivative(bump.list_microcircuit[1].list_soma[1])
	
end

# ╔═╡ Cell order:
# ╠═a12d8f60-3113-11ec-3810-95a5aa96c35e
# ╠═98dd4b78-6e43-4ce7-a599-d7feda16b334
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
# ╠═a2d81277-3e2b-485c-bfc5-c2749a073f09
# ╠═c34f1c5f-6175-4b7c-a8b1-118100160eee
# ╠═239afd2b-20a3-4b62-a203-fd27fc232a48
# ╠═aea0f11c-9fae-4e62-8c9f-c5c8c9e27ce8
