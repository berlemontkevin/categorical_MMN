### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 9b98eac0-94df-11eb-3158-190d8196de7c
begin
	using DrWatson
	quickactivate(@__DIR__)
end

# ╔═╡ 929c9bf8-f561-4deb-94e6-174234ec7b86
begin
	using DifferentialEquations
	using Makie, Plots
	
end

# ╔═╡ afadaa9a-f12b-48fe-b303-224be531ed4c


# ╔═╡ 24f653b6-d511-4808-8cf3-1c6d0736a6b7


# ╔═╡ 955cdcfa-95aa-45f8-8f62-68c7768a5237


# ╔═╡ 70b0cddf-0f5c-41f7-9fb6-1e9eedffe4d7
pyplot()

# ╔═╡ 1f11ecc4-1264-469f-ac4a-2d178e1311bc
begin
	
	f(x,y)= 0.2609.*x .- 0.0497.*y .+ 0.3255
	h(x) = (270.0.*x .- 108.0)./(1.0 .- exp.(-0.1540 .*(270.0.*x .- 108.0)))
	g(x,y) = x./0.1 .- (1.0 .-x).*0.641.*h(f(x,y))./1000
end

# ╔═╡ 8222eb70-1d13-4d17-bd65-28ddead751f0
begin
s2 = 0.2
s1 = 0.001:0.01:1.0
	



end

# ╔═╡ 4cc57ec3-ee67-4fbe-bfcb-786b7a4d159d


# ╔═╡ 8c7d6785-1aff-4cf4-8b71-81691663e597
	plot(h.(f.(s1,s2)))

# ╔═╡ 945c7003-a8b7-415d-833b-b88d54228365


# ╔═╡ Cell order:
# ╠═9b98eac0-94df-11eb-3158-190d8196de7c
# ╠═929c9bf8-f561-4deb-94e6-174234ec7b86
# ╠═afadaa9a-f12b-48fe-b303-224be531ed4c
# ╠═24f653b6-d511-4808-8cf3-1c6d0736a6b7
# ╠═955cdcfa-95aa-45f8-8f62-68c7768a5237
# ╠═70b0cddf-0f5c-41f7-9fb6-1e9eedffe4d7
# ╠═1f11ecc4-1264-469f-ac4a-2d178e1311bc
# ╠═8222eb70-1d13-4d17-bd65-28ddead751f0
# ╠═4cc57ec3-ee67-4fbe-bfcb-786b7a4d159d
# ╠═8c7d6785-1aff-4cf4-8b71-81691663e597
# ╠═945c7003-a8b7-415d-833b-b88d54228365
