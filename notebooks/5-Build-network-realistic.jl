### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 713aba3e-c85d-11eb-3fea-ed5a6366e7aa
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())

end

# ╔═╡ 601088b5-c4ee-44ce-9e51-13a01b13d7f4
begin
	using Revise
	using MyNeurosciencePackage
	using Plots
	plotlyjs()
	using Parameters, JLD2, JuMP
end

# ╔═╡ 7ca9beda-fc04-401b-b73c-1033ec029d08
begin
	
	
function construct_local_microcircuit_integrator_temp(c::microcircuit; dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation = false, adaptation = false, depression = false, td_to_vip = true)
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name

    vip1 = vip_cell(name = "vipcell1", Ibg = 0.25)
    sst1 = sst_cell(name = "sstcell1", Ibg = 0.2)
    pv1 = pv_cell(name = "pvcell1", Ibg = 0.27)
    dend1 = dend_sigmoid(param_c = dend_param,name = "dend1")
    E1 = soma_PC(den=dend1, adaptation_boolean = false,name = "ecell1")
    integrator1 = neural_integrator(name = "integrator1")
    
    
    push!(dend1.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst1,neuron_post=dend1, g = -0.09, name = "sst1-to-dend1"))#0.09
    push!(E1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E1, g = -0.005, name = "pv1-to-ecell1"))#-0.001
    push!(E1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=E1, g = 0.18, name = "ecell1-to-ecell1"))
    
    push!(vip1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=vip1, g = 0.058, name = "ecell1-to-vip1"))
    push!(vip1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=vip1, g = -0.1, name = "sst1-to-vip1"))
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst1, g = 0.1435, name = "ecell1-to-sst1"))
    push!(sst1.list_syn, gaba_syn(neuron_pre=vip1,neuron_post=sst1, g = -0.05, name = "ss1-to-dend1"))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=pv1, g = 0.09435, name = "ecell1-to-pv1"))#0.08435tpicola 
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=pv1, g = -0.17, name = "sst1-to-pv1"))#-0.17
    push!(pv1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=pv1, g = -0.18, name = "pv1-to-pv1"))
    
       
       push!(vip1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=vip1, g = 0.47, depression = depression, name = "integrator1-to-vip1"))
    
        push!(sst1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=sst1, g = 0.22, facilitation = facilitation, name = "integrator1-to-sst1"))
		
		        push!(pv1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=pv1, g = 0.31, depression = true, name = "integrator1-to-pv1"))
		
   
		push!(integrator1.list_syn, nmda_syn(neuron_pre = E1, neuron_post = integrator1, g = 0.25, name = "ecell1-to-integrator1"))

    push!(c.list_dend,dend1)
  
    push!(c.list_soma,E1)
   
    push!(c.list_vip,vip1)
  
    push!(c.list_sst,sst1)
  
    push!(c.list_pv,pv1)

    push!(c.list_integrator,integrator1)



end

	
	
end

# ╔═╡ bf2ec59a-6ac1-4d89-aaa9-e9fd8e10718f
begin 	
	function plot_circuit(c::microcircuit, sim::simulation_parameters;title="")
		PC = c.list_soma[1]
		PV = c.list_pv[1]
		SST = c.list_sst[1]
		VIP = c.list_vip[1]
		Integrator = c.list_integrator[1]
		Dend = c.list_dend[1]
		liste_time = 0.0:sim.dt:sim.Tfin
		
		fig = plot(linewidth=2)
		plot!(fig,liste_time,PC.r,label="PC", lw =2)
		plot!(fig,liste_time,PV.r,label="PV", lw =2)
		plot!(fig,liste_time,SST.r,label="SST", lw =2)
		plot!(fig,liste_time,VIP.r,label="VIP", lw =2)
		plot!(fig,liste_time, Integrator.r, label="integrator", lw =2)
		title!(fig,title)
		xlabel!("Time (s)")
		ylabel!("Firing rate (Hz)")
		return fig
	end

end

# ╔═╡ d3d7bfae-1f70-4ff8-afa1-0c83bc421224
begin
	function constant_stimulus(sim::simulation_parameters,value::Float64)
    sim.Stimulus = value.*ones(Int(sim.Tfin/sim.dt))
    
	end
	
	
function alternating_stimulus(sim::simulation_parameters,value::Float64)
     list_stim = zeros(Int(sim.Tfin/sim.dt))
    temp_t = 0.0
    
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

# ╔═╡ 48aab8eb-19b8-4688-8bea-548d7e480c7f
begin
	c = microcircuit()
    construct_local_microcircuit_integrator_temp(c; depression = true)
	
	sim = simulation_parameters()
    alternating_stimulus(sim,0.1)
    full_time_dynamics(c, sim)

	
end

# ╔═╡ b383e632-3341-495c-9b72-161bd79229c9
begin
	
	plot_circuit(c,sim)
	
end

# ╔═╡ de3fbc1e-6285-40b3-8121-9cf3229741b7
begin
	c_cst = microcircuit()
    construct_local_microcircuit_integrator_temp(c_cst; depression = true)
	
	sim_cst = simulation_parameters()
    constant_stimulus(sim_cst,0.2)
    full_time_dynamics(c_cst, sim_cst)

		plot_circuit(c_cst,sim_cst)

end

# ╔═╡ 07ee78b4-9208-4272-a6d6-33b94261d918
begin
	plot(c.list_vip[1].list_syn[3].d)
	plot!(c_cst.list_vip[1].list_syn[3].d)
	
end

# ╔═╡ 923a9a2a-3b0a-45ba-8da4-7b40a66ba53c
begin
	plot(0.0:sim.dt:sim.Tfin,c.list_vip[1].list_syn[3].s)
	plot!(0.0:sim.dt:sim.Tfin,c_cst.list_vip[1].list_syn[3].s)
	
end

# ╔═╡ 1cf5e336-7e2a-4ec5-8392-7e9223abb7a2
begin
	c_a = microcircuit()
    construct_local_microcircuit_integrator_temp(c_a; facilitation = true, td_to_vip =  false)
	
	sim_a = simulation_parameters()
    alternating_stimulus(sim_a,0.2)
    full_time_dynamics(c_a, sim_a)

	plot(c_a.list_sst[1].list_syn[3].s)

end

# ╔═╡ c10b5f77-bb84-4ff7-805c-13f17a9ecd17
plot(0.0:sim.dt:sim.Tfin,c_a.list_sst[1].list_syn[3].u)


# ╔═╡ 32481654-b508-4af2-b886-17a33c69623e
plot(0.0:sim.dt:sim.Tfin,c.list_sst[1].list_syn[2].s)


# ╔═╡ f316ff15-d8f2-4bf9-a97e-82b4092f2af8
plot(0.0:sim.dt:sim.Tfin,c.list_dend[1].list_syn[1].s)



# ╔═╡ 2708263a-941c-4575-9bc5-c46d057a1bf8
plot(0.0:sim.dt:sim.Tfin,c.list_dend[1].Itot)

# ╔═╡ fe917d95-06e7-4ae7-9f8a-19c0035c502e
plot(0.0:sim.dt:sim.Tfin,c_a.list_dend[1].Itot)

# ╔═╡ Cell order:
# ╠═713aba3e-c85d-11eb-3fea-ed5a6366e7aa
# ╠═601088b5-c4ee-44ce-9e51-13a01b13d7f4
# ╠═7ca9beda-fc04-401b-b73c-1033ec029d08
# ╠═bf2ec59a-6ac1-4d89-aaa9-e9fd8e10718f
# ╠═d3d7bfae-1f70-4ff8-afa1-0c83bc421224
# ╠═48aab8eb-19b8-4688-8bea-548d7e480c7f
# ╠═b383e632-3341-495c-9b72-161bd79229c9
# ╠═de3fbc1e-6285-40b3-8121-9cf3229741b7
# ╠═07ee78b4-9208-4272-a6d6-33b94261d918
# ╠═923a9a2a-3b0a-45ba-8da4-7b40a66ba53c
# ╠═1cf5e336-7e2a-4ec5-8392-7e9223abb7a2
# ╠═c10b5f77-bb84-4ff7-805c-13f17a9ecd17
# ╠═32481654-b508-4af2-b886-17a33c69623e
# ╠═f316ff15-d8f2-4bf9-a97e-82b4092f2af8
# ╠═2708263a-941c-4575-9bc5-c46d057a1bf8
# ╠═fe917d95-06e7-4ae7-9f8a-19c0035c502e
