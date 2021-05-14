### A Pluto.jl notebook ###
# v0.14.5

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

# ╔═╡ 9f89ae64-8d89-4e7d-9f21-b19bb731116d
begin
	# Loading packages
	using Revise
	using Parameters
	push!(LOAD_PATH,srcdir())
	
	using Plots
	plotlyjs()
	using MyNeurosciencePackage
	
end

# ╔═╡ a83fc540-ac1e-11eb-026a-698bbfe91259
begin
	using DrWatson
	quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN","1-project-categorical-MMN")
	
end

# ╔═╡ 92a20122-737d-4050-ac75-47219d611ee6
begin
	using PlutoUI
	@bind stim_current Slider(0.0001:0.1:1.0)
end

# ╔═╡ 45a55b69-e856-4663-9631-6e916ff656d9


# ╔═╡ 463f0ee2-8a46-4d7d-9da1-b267654598f6
dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.15624, 0.0)

# ╔═╡ d62ac590-c7c4-40fa-83f4-8d9dd8e79b45
dend = dend_sigmoid(param_c = dend_param)

# ╔═╡ 13bbb829-14e6-4750-8b2d-d67c19719d84
md""" # Example of a dendritic synapse


"""

# ╔═╡ a9a4daba-1d46-4414-b7d5-845657e3e7b4
begin
	
	soma = soma_PC(den = dend)
	dend.Iinh[end] = -0.8
	update_dend!(dend)
	
	
end

# ╔═╡ 2988f99e-6daf-4585-b746-3d045237eb55
begin

	sum_input!(soma)
	update_firing!(soma)
	
end

# ╔═╡ 470d2d4f-d780-4f96-8c3f-aa7cb4cc1748
dend

# ╔═╡ 6fa7c0d1-d07f-4833-bfe8-31eb267cfe28


# ╔═╡ 760a7129-21e7-48cd-a08d-d165f3675a49


# ╔═╡ ac3d412b-10f7-4aa9-be99-90c0f81454da
md""" # Simple integrator neuron


"""

# ╔═╡ 21d23481-35b4-4c74-96c6-0a938a885c74
begin
	
	integrator = neural_integrator(α=0.1)
	
	for i=1:100
		integrator.Iexc[end] = 0.3
		sum_input!(integrator)
		update_firing!(integrator)
	end
	
	
	
	
end

# ╔═╡ 6dd247a0-2002-4238-9c1f-c1161a11d4c4
plot(0.0:integrator.τ:100*integrator.τ,integrator.r)

# ╔═╡ 4e1c06e1-3109-4b85-85ed-ccd25d87c6aa
plot(integrator.Itot)

# ╔═╡ 9f391c16-04c7-4a3c-b9f7-0d12f6bc48e0
md""" # All the synapses equations

 let's assume that only facilitation exists in this local microcircuit

`` \frac{du}{dt} = \frac{U - u}{τ^ u} + U(1-u) r_E ``

After that, two short term plasticity mechanisms are possible.

For NMDA synapses:
`` \frac{ds}{dt} = - \frac{s}{τ} + u(1-s)γr_E ``

For AMPA synapse:
``\frac{ds}{dt} = - \frac{s}{τ} +u γr_E ``

parameters are ``U = 0.2, τ = 1.5``

All the currents are defined by ``g \times s``


"""

# ╔═╡ 4542719f-3737-4339-bfe7-6e3c0300fc79
md""" # Construction of the network


"""

# ╔═╡ 105ddeef-202b-45e9-bc4c-dd8890398012
begin
	
	function constant_stimulus(sim::simulation_parameters,value::Float64)
		sim.Stimulus = value.*ones(Int(sim.Tfin/sim.dt))
		
	end
	
	
	function alternating_stimulus(sim::simulation_parameters,value::Float64)
		 list_stim = zeros(Int(sim.Tfin/sim.dt))
    	temp_t = 0.0
		sim.Tstimduration = 0.5
		
   	 for i=1:length(list_stim)-1
        temp_t +=sim.dt

        if temp_t>sim.Tstimduration && i>sim.Tstimdebut/sim.dt
            if list_stim[i] == 0.0
                list_stim[i] = value
            else
                list_stim[i] = 0.0
            end
            temp_t = 0.0

        end
        list_stim[i+1] = list_stim[i]
    end
		sim.Stimulus = list_stim
		
	end
	
end

# ╔═╡ 6e2684ca-6061-471c-9a74-5b55a8224911
begin
	# plotting function
	
	function plot_circuit(c::microcircuit, sim::simulation_parameters;title="")
		PC = c.list_soma[1]
		PV = c.list_pv[1]
		SST = c.list_sst[1]
		VIP = c.list_vip[1]
		Integrator = c.list_integrator[1]
		Dend = c.list_dend[1]
		liste_time = 0.0:sim.dt:sim.Tfin
		temp =Integrator.list_syn[1].s[end]
		fig = plot()
		plot!(fig,liste_time,PC.r,label="PC")
		plot!(fig,liste_time,PV.r,label="PV")
		plot!(fig,liste_time,SST.r,label="SST")
		plot!(fig,liste_time,VIP.r,label="VIP")
		plot!(fig,liste_time, Integrator.r, label="integrator $temp")
		title!(fig,title)
		return fig
	end
	
	
	function make_dynamics(keywords::Array{String},sim::simulation_parameters;facilitation = false)
		#todo replace keywords by something nicer
		if keywords[2] == "constant"
			constant_stimulus(sim,stim_current)

		elseif keywords[2] == "alternating"
				alternating_stimulus(sim,stim_current)

		else
			warn("Don't forget to add a type of input to your neuron. Keywords[2]")
		end
		
		circuit = create_network(keywords[1]; facilitation)
		full_time_dynamics(circuit, sim)
		
		return circuit

	end
	
	
end

# ╔═╡ d6d45f7c-5da9-4c09-841e-b6185e5f40de
0.12/0.005/2.0

# ╔═╡ 45b87a22-dacf-432b-b631-6267066a76c2
begin
	
	simulation_param = simulation_parameters()
	
	circuit_integrator = make_dynamics(["integrator","constant"],simulation_param)#create_network("integrator")
	
	
	#constant_stimulus(simulation_param,stim_current)
	
	#circuit = create_network()

	#full_time_dynamics(circuit_integrator, simulation_param)
	
	#full_time_dynamics(circuit, simulation_param)
	
	
end

# ╔═╡ e0b1d724-9f81-49ab-a559-1df99ad0bd18
plot(circuit_integrator.list_dend[1].list_syn[1].s)

# ╔═╡ f2fcfd9f-232d-4e9a-8057-a6d751dd7bf1
plot(circuit_integrator.list_vip[1].list_syn[2].s)

# ╔═╡ 8774e37d-1177-41f7-a171-82e1ef775f5c
plot(circuit_integrator.list_dend[1].Iinh)

# ╔═╡ d2673b6f-1d68-4216-bead-0222a655fbed
fig_constant_wt_adaptation = plot_circuit(circuit_integrator,simulation_param, title ="Constant stim and without facilitation")

# ╔═╡ 43aafa43-43e7-4ca7-ab6a-14104ed10a30
begin
	
	
	circuit_integrator_facilitation = make_dynamics(["integrator","constant"],simulation_param; facilitation = true)
	

	
	
end

# ╔═╡ 80e3502e-e5ee-402d-9f6b-56836c3653a2
fig_constant_with_facilitation = plot_circuit(circuit_integrator_facilitation,simulation_param, title ="Constant stim and with facilitation")

# ╔═╡ a24b91f3-9d8b-4c13-a7ef-c6ab1aef698e
md"""
Study of a non constant stimulus on the network
"""

# ╔═╡ 3229a7c2-0d3b-4d5f-a991-5aff5a9b77c8
begin
	circuit_integrator_alternating = make_dynamics(["integrator","alternating"],simulation_param)
		
	fig_alternating_wt_adaptation = plot_circuit(circuit_integrator_alternating,simulation_param, title ="Alternating stim and without facilitation")
end

# ╔═╡ 2d407f57-c901-4f07-812e-2d6e5c7d035c
begin
	circuit_integrator_alternating_facilitation = make_dynamics(["integrator","alternating"],simulation_param; facilitation = true)
		
	fig_alternating_with_facilitation = plot_circuit(circuit_integrator_alternating_facilitation,simulation_param, title ="Alternating stim and with facilitation")
end

# ╔═╡ 5bc03f96-b5b8-40bd-9537-e7c3cfec4d11
md""" # Adaptation of pyramidal cells

"""

# ╔═╡ 90904108-a34d-450a-8ae6-6f1eeaf75b32
begin
	circuit_integrator_alternating_adaptation = create_network("integrator",facilitation=false,adaptation = true)
		
	circuit_integrator_alternating_adaptation.list_soma[1].adaptation.gA = -0.01
	circuit_integrator_alternating_adaptation.list_soma[1].adaptation.τA = 0.1

	
	alternating_stimulus(simulation_param,0.4)

	full_time_dynamics(circuit_integrator_alternating_adaptation, simulation_param)
	
	fig_alternating_with_adaptation = plot_circuit(circuit_integrator_alternating_adaptation,simulation_param, title ="Alternating stim and with adaptation")
end

# ╔═╡ 7397c618-5651-4799-a06a-028c144ba068
plot(0.0:simulation_param.dt:simulation_param.Tfin,circuit_integrator_alternating_adaptation.list_soma[1].adaptation.sA)

# ╔═╡ d8d2c7c2-5771-48f3-8127-28a708705bb9
begin
	circuit_integrator_constant_adaptation = create_network("integrator",facilitation=false,adaptation = true)
		
	circuit_integrator_constant_adaptation.list_soma[1].adaptation.gA = -0.01
		circuit_integrator_constant_adaptation.list_soma[1].adaptation.τA = 1.0

	
	constant_stimulus(simulation_param,0.2)

	full_time_dynamics(circuit_integrator_constant_adaptation, simulation_param)
	
	fig_constant_with_adaptation = plot_circuit(circuit_integrator_constant_adaptation,simulation_param, title ="constant stim and with adaptation")
end

# ╔═╡ f4555ff6-6b0b-40c4-a313-8f4da70eafab
plot(0.0:simulation_param.dt:simulation_param.Tfin,circuit_integrator_constant_adaptation.list_soma[1].adaptation.sA)

# ╔═╡ e4bce060-6041-4e79-b75b-26826d8f11f8
plot(0.0:simulation_param.dt:simulation_param.Tfin,circuit_integrator_alternating_adaptation.list_soma[1].adaptation.sA)

# ╔═╡ 2bce631a-c71f-457d-ac33-180b7cf5aefa
plot(circuit_integrator_alternating_adaptation.list_soma[1].OU_process.noise)

# ╔═╡ 41aca281-92f6-4025-bf40-93e076196b2a
plot(circuit_integrator_alternating_adaptation.list_soma[1].Itot)

# ╔═╡ Cell order:
# ╠═a83fc540-ac1e-11eb-026a-698bbfe91259
# ╠═9f89ae64-8d89-4e7d-9f21-b19bb731116d
# ╠═45a55b69-e856-4663-9631-6e916ff656d9
# ╠═463f0ee2-8a46-4d7d-9da1-b267654598f6
# ╠═d62ac590-c7c4-40fa-83f4-8d9dd8e79b45
# ╠═13bbb829-14e6-4750-8b2d-d67c19719d84
# ╠═a9a4daba-1d46-4414-b7d5-845657e3e7b4
# ╠═2988f99e-6daf-4585-b746-3d045237eb55
# ╠═470d2d4f-d780-4f96-8c3f-aa7cb4cc1748
# ╠═6fa7c0d1-d07f-4833-bfe8-31eb267cfe28
# ╠═760a7129-21e7-48cd-a08d-d165f3675a49
# ╠═ac3d412b-10f7-4aa9-be99-90c0f81454da
# ╠═21d23481-35b4-4c74-96c6-0a938a885c74
# ╠═6dd247a0-2002-4238-9c1f-c1161a11d4c4
# ╠═4e1c06e1-3109-4b85-85ed-ccd25d87c6aa
# ╠═9f391c16-04c7-4a3c-b9f7-0d12f6bc48e0
# ╠═4542719f-3737-4339-bfe7-6e3c0300fc79
# ╠═105ddeef-202b-45e9-bc4c-dd8890398012
# ╠═6e2684ca-6061-471c-9a74-5b55a8224911
# ╠═92a20122-737d-4050-ac75-47219d611ee6
# ╠═e0b1d724-9f81-49ab-a559-1df99ad0bd18
# ╠═f2fcfd9f-232d-4e9a-8057-a6d751dd7bf1
# ╠═8774e37d-1177-41f7-a171-82e1ef775f5c
# ╠═d6d45f7c-5da9-4c09-841e-b6185e5f40de
# ╠═45b87a22-dacf-432b-b631-6267066a76c2
# ╠═d2673b6f-1d68-4216-bead-0222a655fbed
# ╠═43aafa43-43e7-4ca7-ab6a-14104ed10a30
# ╠═80e3502e-e5ee-402d-9f6b-56836c3653a2
# ╠═a24b91f3-9d8b-4c13-a7ef-c6ab1aef698e
# ╠═3229a7c2-0d3b-4d5f-a991-5aff5a9b77c8
# ╠═2d407f57-c901-4f07-812e-2d6e5c7d035c
# ╠═5bc03f96-b5b8-40bd-9537-e7c3cfec4d11
# ╠═90904108-a34d-450a-8ae6-6f1eeaf75b32
# ╠═7397c618-5651-4799-a06a-028c144ba068
# ╠═d8d2c7c2-5771-48f3-8127-28a708705bb9
# ╠═f4555ff6-6b0b-40c4-a313-8f4da70eafab
# ╠═e4bce060-6041-4e79-b75b-26826d8f11f8
# ╠═2bce631a-c71f-457d-ac33-180b7cf5aefa
# ╠═41aca281-92f6-4025-bf40-93e076196b2a
