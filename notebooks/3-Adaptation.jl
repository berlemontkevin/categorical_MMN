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
	using PlutoUI

end

# ╔═╡ a83fc540-ac1e-11eb-026a-698bbfe91259
begin
	using DrWatson
	quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN","1-project-categorical-MMN")
	
end

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

# ╔═╡ 92a20122-737d-4050-ac75-47219d611ee6
md""" Value of the current sent to the soma:

0.1 to 1.0 $@bind stim_current Slider(0.1:0.1:1.0, default = 0.5) 

"""


# ╔═╡ 62256341-cfe2-455c-97df-d5e174e58347
md""" Value of the time constant of adaptation and facilitation:

0.1 to 2.0 $@bind τ Slider(0.1:0.1:2.0, default = 0.5)


"""

# ╔═╡ de41db49-4c77-479c-85e4-5a830c3c0a70
md""" Value of the inter stimulus interval:

0.1 to 1.0 $@bind T_inter Slider(0.1:0.2:1.0, default = 0.5)


"""

# ╔═╡ df760430-5d76-43fa-90ba-cafb32638321
md"""
The chosen parameters are:

stim-current = $stim_current

``\tau`` = $τ

Tinter = $T_inter



"""

# ╔═╡ b213aced-a8a5-4f0a-b8fb-944b58042180
md""" ### Constant current to the network

"""

# ╔═╡ 9254d341-8763-47ee-90b3-5f1bbb015395
md""" ## Alternating current to the soma


"""

# ╔═╡ 3707a172-d950-46f0-8bae-5b6d0f262166
md""" ## Analyze of alternating currents

$@bind value_duration NumberField(0:0.08, default=0.02)


let's look at the variation of the $value_duration 
s duration peak through time.

"""

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
	
	
	function make_dynamics(keywords::Array{String},sim::simulation_parameters;facilitation = false,stim=stim_current)
		#todo replace keywords by something nicer
		if keywords[2] == "constant"
			constant_stimulus(sim,stim)

		elseif keywords[2] == "alternating"
				alternating_stimulus(sim,stim)

		else
			warn("Don't forget to add a type of input to your neuron. Keywords[2]")
		end
		
		circuit = create_network(keywords[1]; facilitation)
		full_time_dynamics(circuit, sim)
		
		return circuit

	end
	
	
	
	function analyze_alternating(network::microcircuit,sim::simulation_parameters)
		list_values = []
		PC=network.list_soma[1]
		dt = sim.dt
		stimd = sim.Tstimduration
		temp_dt = sim.Tstimduration + value_duration
		while temp_dt <= sim.Tfin
			push!(list_values, PC.r[Int(round(temp_dt/dt))])
			temp_dt+=2.0*sim.Tstimduration
		end
		
		return list_values
	end
	
	
	using LsqFit
	
	model_sig(x, p) = p[1].*sigmoid.(-(x.-p[4])./p[2]) .+ p[3]#p[1].*exp.(-x./p[2]) .+ p[3]
model_exp(x, p) = p[1].*exp.(-(x)./p[2]) .+ p[3]#p[1].*exp.(-x./p[2]) .+ p[3]

function get_param_fit(c::microcircuit, sim::simulation_parameters,model)
    lb = [1e-9,0.1,1.0,1e-9]
    #lb = [0.0,0.0,0.0,0.0]
       values_fit = analyze_alternating(c,sim)
    ub = [2.0*maximum(values_fit),6.0,maximum(values_fit),30.0]
    p0 = [0.5*maximum(values_fit),0.5,0.5*maximum(values_fit),1]

       xdata_stim = 1.5*sim.Tstimduration:2.0*sim.Tstimduration:sim.Tfin
    fit = curve_fit(model, xdata_stim, values_fit, p0, lower = lb, upper = ub)

		return xdata_stim, fit.param,values_fit,maximum(values_fit)
end


	
	
end

# ╔═╡ 45b87a22-dacf-432b-b631-6267066a76c2
begin
	
 	simulation_param_constant = simulation_parameters()
	
 	circuit_integrator = make_dynamics(["integrator","constant"],simulation_param_constant)
	
 	circuit_integrator_facilitation = make_dynamics(["integrator","constant"],simulation_param_constant; facilitation = true)
	
	

d_adapt=Dict("facilitation" => false,
"adaptation" => true,
            "stimlist" =>stim_current,
                "τ" => τ,
                "Tinterstim" => T_inter
                )
temp_load=    wload(datadir("simulations","notebook-3", savename(d_adapt, "jld2")))
circuit_adapt = temp_load["circuit"]
simulation_param = temp_load["simulation_param"]


d_control=Dict("facilitation" => false,
"adaptation" => false,
            "stimlist" =>stim_current,
                "τ" => τ,
                "Tinterstim" => T_inter
                )
temp_load_control=    wload(datadir("simulations","notebook-3", savename(d_control, "jld2")))
circuit_control = temp_load_control["circuit"]
	
d_facilitation=Dict("facilitation" => true,
"adaptation" => false,
            "stimlist" =>stim_current,
                "τ" => τ,
                "Tinterstim" => T_inter
                )
temp_load_facilitation=    wload(datadir("simulations","notebook-3", savename(d_facilitation, "jld2")))
circuit_facilitation = temp_load_facilitation["circuit"]

	
	
	
	
	
end

# ╔═╡ ba159d20-f0f7-4068-9a5b-24d1102d05a5
fig_constant_wt_adaptation = plot_circuit(circuit_integrator,simulation_param, title ="Constant stim and without facilitation")

# ╔═╡ 80e3502e-e5ee-402d-9f6b-56836c3653a2
fig_constant_with_facilitation = plot_circuit(circuit_integrator_facilitation,simulation_param, title ="Constant stim and with facilitation")

# ╔═╡ d8d2c7c2-5771-48f3-8127-28a708705bb9
begin
	circuit_integrator_constant_adaptation = create_network("integrator",facilitation=false,adaptation = true)
		
	circuit_integrator_constant_adaptation.list_soma[1].adaptation.gA = -0.01
	circuit_integrator_constant_adaptation.list_soma[1].adaptation.τA = τ

	
	constant_stimulus(simulation_param_constant,stim_current)

	full_time_dynamics(circuit_integrator_constant_adaptation, simulation_param_constant)
	
	fig_constant_with_adaptation = plot_circuit(circuit_integrator_constant_adaptation,simulation_param_constant, title ="constant stim and with adaptation")
end

# ╔═╡ 3229a7c2-0d3b-4d5f-a991-5aff5a9b77c8
fig_alternating_wt_adaptation = 	plot_circuit(circuit_control,simulation_param, title ="Alternating stim and without facilitation")


# ╔═╡ 2d407f57-c901-4f07-812e-2d6e5c7d035c
	fig_alternating_with_facilitation = plot_circuit(circuit_facilitation,simulation_param, title ="Alternating stim and with facilitation")


# ╔═╡ 6b82843b-1ad6-4e52-b44c-c3bfc3074d11
	fig_alternating_adaptation = plot_circuit(circuit_adapt,simulation_param, title ="Alternating stim and with adaptation")

# ╔═╡ 90904108-a34d-450a-8ae6-6f1eeaf75b32
begin
	function plot_adaptation_alternating(τA::Float64)
		
		circuit_integrator_alternating_adaptation = create_network("integrator",facilitation=false,adaptation = true)
		
	circuit_integrator_alternating_adaptation.list_soma[1].adaptation.gA = -0.01
	circuit_integrator_alternating_adaptation.list_soma[1].adaptation.τA = τA

	
	alternating_stimulus(simulation_param,0.4)

	full_time_dynamics(circuit_integrator_alternating_adaptation, simulation_param)
	
	fig_alternating_with_adaptation = plot_circuit(circuit_integrator_alternating_adaptation,simulation_param, title ="Alternating stim and with adaptation")
	return fig_alternating_with_adaptation
	end
end

# ╔═╡ 71ac3678-66d6-4f8c-ae0c-6f5792000352
begin
	values_facilitation = analyze_alternating(circuit_facilitation,simulation_param)
	
	values_wtfacilitation = analyze_alternating(circuit_control,simulation_param)
	
		values_adaptation = analyze_alternating(circuit_adapt,simulation_param)
	
	
	scatter(values_facilitation,color=:green,label="Facilitation")
	scatter!(values_wtfacilitation,color=:blue,label="Control")
	scatter!(values_adaptation,color = :red,label="adaptation")
	xlabel!("Number of peaks")
	ylabel!("Firing rate")
end

# ╔═╡ 1963ddc5-8a3d-4fe4-89fa-51497afdd24c
md"""
Let's fit control and adaptation curve by a sigmoid. For the adaptation case it can be done but not for all parameters due to the look of the dynamics.


`` f(x) = p1 * sigmoid(-(x-p_2)/p_4) + p_3 ``

"""

# ╔═╡ ef413bdf-bc21-4b02-b9f2-f4eda91288a6
begin
	
	
d_adapt_fast=Dict("facilitation" => false,
"adaptation" => true,
            "stimlist" =>0.5,
                "τ" => 0.1,
                "Tinterstim" => 0.5
                )
temp_load_fast=    wload(datadir("simulations","notebook-3", savename(d_adapt_fast, "jld2")))
circuit_adapt_fast = temp_load_fast["circuit"]
simulation_param_fast = temp_load_fast["simulation_param"]

	
	
d_control_fast=Dict("facilitation" => false,
"adaptation" => false,
            "stimlist" =>0.5,
                "τ" => 0.1,
                "Tinterstim" => 0.5
                )
temp_control_fast=    wload(datadir("simulations","notebook-3", savename(d_control_fast, "jld2")))
circuit_control_fast = temp_control_fast["circuit"]
	
	
d_f_fast=Dict("facilitation" => true,
"adaptation" => false,
            "stimlist" =>0.5,
                "τ" => 1.5,
                "Tinterstim" => 0.5
                )
temp_f_fast=    wload(datadir("simulations","notebook-3", savename(d_f_fast, "jld2")))
circuit_f_fast = temp_f_fast["circuit"]
	
	values_facilitation_fast = analyze_alternating(circuit_f_fast,simulation_param_fast)
	
	values_wtfacilitation_fast = analyze_alternating(circuit_control_fast,simulation_param_fast)
	
	values_adaptation_fast = analyze_alternating(circuit_adapt_fast,simulation_param_fast)
	
	
	scatter(values_facilitation_fast,color=:green,label="Facilitation")
	scatter!(values_wtfacilitation_fast,color=:blue,label="Control")
	scatter!(values_adaptation_fast,color = :red,label="adaptation")
	xlabel!("Number of peaks")
	ylabel!("Firing rate")
	
	
    xdata_stim, fit_param,m,n = get_param_fit(circuit_f_fast,simulation_param_fast,model_sig) 
	plot!(xdata_stim,model_sig(xdata_stim,fit_param),color=:green,linewidth=2)
	
	
    xdata_stim, fit_param_adapt,ma,na = get_param_fit(circuit_adapt_fast,simulation_param_fast,model_sig) 
	plot!(xdata_stim,model_sig(xdata_stim,fit_param_adapt),color=:red,linewidth=2)
	
	
    xdata_stim, fit_param_control,mc,nc = get_param_fit(circuit_control_fast,simulation_param_fast,model_sig) 
	plot!(xdata_stim,model_sig(xdata_stim,fit_param_control),color=:blue,linewidth=2)
end

# ╔═╡ ba0d669c-1acc-4866-823f-d075fb7f52bb
md"""
$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p1-p3-fixedstim-facilitation.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p1-p3-fixedstim-control.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p2-fixedstim-facilitation.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p2-fixedstim-control.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p4-fixedstim-facilitation.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p4-fixedstim-control.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p1-p3-fixedTinter-facilitation.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p1-p3-fixedTinter-control.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\3-adaptation\\fig-p1-p3-fixedtau-fixedTinter.png"))


"""

# ╔═╡ Cell order:
# ╠═9f89ae64-8d89-4e7d-9f21-b19bb731116d
# ╠═a83fc540-ac1e-11eb-026a-698bbfe91259
# ╟─9f391c16-04c7-4a3c-b9f7-0d12f6bc48e0
# ╟─4542719f-3737-4339-bfe7-6e3c0300fc79
# ╟─105ddeef-202b-45e9-bc4c-dd8890398012
# ╟─6e2684ca-6061-471c-9a74-5b55a8224911
# ╟─92a20122-737d-4050-ac75-47219d611ee6
# ╟─62256341-cfe2-455c-97df-d5e174e58347
# ╟─de41db49-4c77-479c-85e4-5a830c3c0a70
# ╟─df760430-5d76-43fa-90ba-cafb32638321
# ╟─b213aced-a8a5-4f0a-b8fb-944b58042180
# ╟─ba159d20-f0f7-4068-9a5b-24d1102d05a5
# ╟─80e3502e-e5ee-402d-9f6b-56836c3653a2
# ╟─d8d2c7c2-5771-48f3-8127-28a708705bb9
# ╟─45b87a22-dacf-432b-b631-6267066a76c2
# ╟─9254d341-8763-47ee-90b3-5f1bbb015395
# ╟─3229a7c2-0d3b-4d5f-a991-5aff5a9b77c8
# ╟─2d407f57-c901-4f07-812e-2d6e5c7d035c
# ╟─6b82843b-1ad6-4e52-b44c-c3bfc3074d11
# ╟─90904108-a34d-450a-8ae6-6f1eeaf75b32
# ╟─3707a172-d950-46f0-8bae-5b6d0f262166
# ╟─71ac3678-66d6-4f8c-ae0c-6f5792000352
# ╟─1963ddc5-8a3d-4fe4-89fa-51497afdd24c
# ╟─ef413bdf-bc21-4b02-b9f2-f4eda91288a6
# ╠═ba0d669c-1acc-4866-823f-d075fb7f52bb
