### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 83e0f230-d394-11eb-10d1-533506e356d5
begin
		using DrWatson
	
	quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())
	using CSV, DataFrames, Query
	using Makie, GLMakie
	
end

# ╔═╡ 1d0effdf-59a9-42cd-b32c-fc2929271de8
md"""
# Brief analyze of Kennedy's data

From: Markov paper with the 29-area data (Cerebral Cortex) and Sean dopamine paper

The following figure represents the FLN matrix for the 40 areas.
"""

# ╔═╡ bd78b913-7944-4c66-9d68-c6acbb196177
FLNdata = DataFrame(CSV.File("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\_research\\walsh-40-region-kennedy-hierarchy-master\\FLNdata.csv", header=false))

# ╔═╡ ce6dda35-ce40-4b22-bd7d-913f871d60db
begin
	
	areas_by_hier = ["V1",
"V2",
"V4",
"3",
"1",
"V6",
"DP",
"MT",
"8m",
"8l",
"2",
"5",
"TEO",
"F1",
"F4",
"STPc",
"7a",
"F3",
"F5",
"46d",
"9/46v",
"9/46d",
"10",
"7m",
"TEpd",
"PBr",
"LIP",
"32",
"F2",
"9",
"25",
"7b",
"STPi",
"45A",
"F7",
"ProM",
"'8B",
"24c",
"STPr",
"OPRO"]
	
end

# ╔═╡ ab8693a5-7449-47aa-8a33-9f1ce74ae6fe
begin
	#pyplot()
	set_theme!(theme_light())

	df_FLN=Matrix(FLNdata)
	
	
	f = Figure()
	ax = Axis(f[1, 1],xticksize=1.0)
	hm = heatmap!(log10.(df_FLN))
	ax.xticks = (1:1:40,areas_by_hier)
	ax.xticklabelrotation = 2pi/8
	ax.yticks = (1:1:40,areas_by_hier)
	ax.xlabel="Projecting from"
	ax.ylabel = "Projecting to"
	Colorbar(f[1, 2], hm)

	f
	
	
end

# ╔═╡ 261a9569-4803-4783-ace4-a8097ebc3a95
md""" ## Extraction of connections to MT

"""

# ╔═╡ 0c9fe8c1-4dae-4bf0-b863-537d4b923f17
begin
	
	to_MT = DataFrame(From = areas_by_hier, Connection_strength = df_FLN[findfirst(isequal("MT"),areas_by_hier),:])
end

# ╔═╡ 5013ed1c-0d80-47da-b26b-e74d0994aa79
md"""
The previous table shows the connection strengths toward MT. The following highlights the strongest ones.
"""

# ╔═╡ 51a94217-f74a-4a3b-bb32-1a5b2807a1b9
begin 
	
to_MT_significant = to_MT |>
    @filter(_.Connection_strength > 0.001) |>
    @orderby_descending(_.Connection_strength) |>
    DataFrame
end

# ╔═╡ Cell order:
# ╠═83e0f230-d394-11eb-10d1-533506e356d5
# ╠═1d0effdf-59a9-42cd-b32c-fc2929271de8
# ╟─bd78b913-7944-4c66-9d68-c6acbb196177
# ╟─ce6dda35-ce40-4b22-bd7d-913f871d60db
# ╠═ab8693a5-7449-47aa-8a33-9f1ce74ae6fe
# ╠═261a9569-4803-4783-ace4-a8097ebc3a95
# ╠═0c9fe8c1-4dae-4bf0-b863-537d4b923f17
# ╟─5013ed1c-0d80-47da-b26b-e74d0994aa79
# ╠═51a94217-f74a-4a3b-bb32-1a5b2807a1b9
