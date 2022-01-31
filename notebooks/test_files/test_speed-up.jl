### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ b2a49532-eb00-11eb-21eb-2f173ae5d57b
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())

end

# ╔═╡ a0fa119b-08b1-44f8-b253-71fccf78f2cc
begin
	using Revise
	using MyNeurosciencePackage
	using CairoMakie
	using Parameters, JLD2, JuMP, Statistics
end

# ╔═╡ d436dc2b-107b-4a3b-a1c6-2c5969fc7714
begin
	

#TODO add name field to microcircuit
function construct_one_local_microcircuit_integrator!(c::microcircuit; dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation = false, adaptation = false, depression = false, td_to_vip = true)
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name

    vip1 = vip_cell(name = string(c.name, "-","vipcell1"),Ibg = 0.25)
    sst1 = sst_cell(name = string(c.name, "-","sstcell1"), Ibg = 0.25,adaptation_boolean = true)
    pv1 = pv_cell(name = string(c.name, "-","pvcell1"), Ibg = 0.29)
    dend1 = dend_sigmoid(param_c = dend_param,name = string(c.name, "-","dend1"))
    E1 = soma_PC(den=dend1, adaptation_boolean = adaptation,name = string(c.name, "-","ecell1"))
    integrator1 = neural_integrator(name = string(c.name, "-","integrator1"))
    
    
    push!(dend1.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst1,neuron_post=dend1, g = -0.09, name = string(c.name, "-","sst1-to-dend1")))#0.09
    push!(E1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E1, g = -0.001,depression = true, name = string(c.name, "-","pv1-to-ecell1")))#-0.001
    push!(E1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=E1, g = 0.18, name = string(c.name, "-","ecell1-to-ecell1")))
    
    push!(vip1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=vip1, g = 0.058, depression = true, name = string(c.name, "-","ecell1-to-vip1")))
    push!(vip1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=vip1, g = -0.1, facilitation = true, name = string(c.name, "-","sst1-to-vip1")))
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst1, g = 0.0435, facilitation = true, name = string(c.name, "-","ecell1-to-sst1")))
    push!(sst1.list_syn, gaba_syn(neuron_pre=vip1,neuron_post=sst1, g = -0.05, facilitation = true, name = string(c.name, "-","ss1-to-dend1")))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=pv1, g = 0.04435, depression = true, name = string(c.name, "-","ecell1-to-pv1")))#0.0435
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=pv1, g = -0.17, name = string(c.name, "-","sst1-to-pv1")))#-0.17
    push!(pv1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=pv1, g = -0.18, name = string(c.name, "-","pv1-to-pv1")))
    
       
   
        push!(vip1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=vip1, g = 0.47, depression = depression, name = string(c.name, "-","integrator1-to-vip1")))
    push!(pv1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=pv1, g = 0.31, depression = depression, name = string(c.name, "-","integrator1-to-pv1")))
    
		
		push!(sst1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=sst1, g = 0.22, facilitation = facilitation, name = string(c.name, "-","integrator1-to-sst1")))
		
    push!(integrator1.list_syn, nmda_syn(neuron_pre = E1, neuron_post = integrator1, g = 0.15, name = string(c.name, "-","ecell1-to-integrator1")))

	    push!(dend1.list_syn, nmda_syn(neuron_pre = integrator1, neuron_post = dend1, g = 0.4,depression = true, name = string(c.name, "-","integrator1-to-dend1")))	
		
		
    push!(c.list_dend,dend1)
  
    push!(c.list_soma,E1)
   
    push!(c.list_vip,vip1)
  
    push!(c.list_sst,sst1)
  
    push!(c.list_pv,pv1)

    push!(c.list_integrator,integrator1)



end

#TODO add a list of microcircuit as argument to the dynamics functions

function construct_two_local_microcircuit_integrator(; dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation = false, adaptation = false, depression = false, td_to_vip = true)
    # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name


    c1 = microcircuit(name="microcircuit1")
    c2 = microcircuit(name = "microcircuit2")
    
    construct_one_local_microcircuit_integrator!(c1; dend_param, facilitation, adaptation, td_to_vip, depression)
    construct_one_local_microcircuit_integrator!(c2; dend_param, facilitation, adaptation, td_to_vip, depression)


    push!(c1.list_sst[1].list_syn, nmda_syn(neuron_pre=c2.list_soma[1],neuron_post=c1.list_sst[1], g = 0.0435, facilitation = true, name = string("ecell2-to-sst1")))
    push!(c2.list_sst[1].list_syn, nmda_syn(neuron_pre=c1.list_soma[1],neuron_post=c2.list_sst[1], g = 0.0435, facilitation = true, name = string("ecell1-to-sst2")))

   push!(c1.list_dend[1].list_syn, nmda_syn(neuron_pre=c2.list_integrator[1],neuron_post=c1.list_dend[1], g = 0.1, depression = true, name = string("integrator2-to-ecell1")))
     push!(c2.list_dend[1].list_syn, nmda_syn(neuron_pre=c1.list_integrator[1],neuron_post=c2.list_dend[1], g = 0.1, depression = true, name = string("integrator1-to-ecell2")))

		
		
		
    list_microcircuit = [c1,c2]

    return list_microcircuit
end

	
	
end

# ╔═╡ 371793e6-9bfe-4daf-a26b-0aa89c82b167
temp_dict=Dict(
"microcircuit1-vipcell1" =>zeros(t_tot), 
"microcircuit1-sstcell1" =>zeros(t_tot),
"microcircuit1-pvcell1" =>zeros(t_tot),
"microcircuit1-ecell1" =>zeros(t_tot),
"microcircuit1-integrator1" =>zeros(t_tot),
"microcircuit2-vipcell1" =>zeros(t_tot),
"microcircuit2-sstcell1" =>zeros(t_tot),
"microcircuit2-pvcell1" =>zeros(t_tot),
"microcircuit2-ecell1" =>zeros(t_tot),
"microcircuit2-integrator1" =>zeros(t_tot)
)

# ╔═╡ f02a080a-e259-4ec8-83af-ec896297ba8b
begin
	
	function create_stimu(stim_strength::Float64,sim::simulation_parameters; Tsim = 0.500, Tinterval=1) #in s
		
		Tstop = sim.Tfin
		dt = sim.dt
		
		
		
		
		dt_pattern = Int((Tsim+Tinterval)/dt)
		temp_stim = vcat(0.0.*ones(Int(Tinterval/dt)), stim_strength.*ones(Int(Tsim/dt)))
		dummy = 0
		stim = copy(temp_stim)
		
		while dummy < Int(Tstop/dt)
		
			stim = vcat(stim,temp_stim)
			dummy += dt_pattern
			
		end
		print(dummy)
		# only take the t_tot number of elements
		return stim[1:Int(Tstop/dt)]
		
	end
	
	
	
	
	
end

# ╔═╡ 08f3d0a8-7897-4c3f-9815-12f5d83eb688
sim

# ╔═╡ 2bec0a29-96e0-4f5b-9324-d5429ffbe713
begin
	
	#needs to construct a general array with all the stim. The question is, how to make it automatic??? Symbolic Julia?
	function generate_oddball_task(param::Dict{String,Float64})#f1::Float64, value_stim_f1::Float64; Tinter = 1.0, Tstim = 0.5, Tfin = 60.0, value_stim_f2 = 0.0)
	
	@unpack f1, value_stim_f1,Tinter, Tstim, Tfin, value_stim_f2 = param
		
	f2 = 1-f1
	dt = 0.0005
	t_tot = Int64(Tfin/dt)
		
		
		# create the network with all adapt, depression and so on (TODO the main fig)
	lc = construct_two_local_microcircuit_integrator(facilitation=true, depression = true, adaptation = true)
		
		
	sim = simulation_parameters(Tstimduration = Tstim, TISI = Tinter, Tfin=Tfin, Tinit = 2000)
	
		t_tot +=2000 # Adding a resting state time

		
		#For now dictionnary is define manually
		oddball=Dict(
"microcircuit1-vipcell1" =>zeros(t_tot), 
"microcircuit1-sstcell1" =>zeros(t_tot),
"microcircuit1-pvcell1" =>zeros(t_tot),
"microcircuit1-ecell1" =>zeros(t_tot),
"microcircuit1-integrator1" =>zeros(t_tot),
"microcircuit2-vipcell1" =>zeros(t_tot),
"microcircuit2-sstcell1" =>zeros(t_tot),
"microcircuit2-pvcell1" =>zeros(t_tot),
"microcircuit2-ecell1" =>zeros(t_tot),
"microcircuit2-integrator1" =>zeros(t_tot)
)
	
	array_stim_f1 = zeros(floor(Int,Tfin/dt))
	array_stim_f2 = zeros(floor(Int,Tfin/dt))

	nbr_stim = floor(Int,Tfin/(Tinter + Tstim))
	temp_index_stim = 0
	while temp_index_stim < nbr_stim-1
		
		if rand()<f1
		array_stim_f1[floor(Int,1+temp_index_stim*(Tinter+Tstim)/dt):floor(Int,(temp_index_stim+1)*(Tinter+Tstim)/dt)] = [zeros(floor(Int,Tinter/dt)) ; value_stim_f1*ones(floor(Int,Tstim/dt))]
			
		array_stim_f2[floor(Int,1+temp_index_stim*(Tinter+Tstim)/dt):floor(Int,(temp_index_stim+1)*(Tinter+Tstim)/dt)] = [zeros(floor(Int,Tinter/dt)) ; value_stim_f2*ones(floor(Int,Tstim/dt))]
		else
			
			array_stim_f1[floor(Int,1+temp_index_stim*(Tinter+Tstim)/dt):floor(Int,(temp_index_stim+1)*(Tinter+Tstim)/dt)] = [zeros(floor(Int,Tinter/dt)) ; value_stim_f2*ones(floor(Int,Tstim/dt))]
			
			array_stim_f2[floor(Int,1+temp_index_stim*(Tinter+Tstim)/dt):floor(Int,(temp_index_stim+1)*(Tinter+Tstim)/dt)] = [zeros(floor(Int,Tinter/dt)) ; value_stim_f1*ones(floor(Int,Tstim/dt))]
		end
			
			temp_index_stim +=1
				
	end
	# for now, the end of the array is 0
	
	
	
	# of course it will be necessary to transform it in a function
	
	oddball["microcircuit1-ecell1"] = [0.1*ones(2000);array_stim_f1]
	oddball["microcircuit2-ecell1"] = [0.1*ones(2000);array_stim_f2]
	
		full_time_dynamics(lc,sim,oddball)
	
	times_stim = collect((2000+Int((Tinter+0.01)/dt)):Int((Tstim+Tinter)/dt):Int(Tfin/dt+2000))
			
			d = Dict()
			 for c in lc
  			  for pop in [c.list_pv, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
			# for now dend is not saved
					for n in pop
					d[n.name] = n.r_save[times_stim]
					
				end
				end
				
			end
			
			d["current-to-e1"] = oddball["microcircuit1-ecell1"][times_stim]
			d["current-to-e2"] = oddball["microcircuit2-ecell1"][times_stim]
			d["times_step"] = times_stim
	#Need to save all firing rates, all times and currents
			
			
# 	true_times_stim = times_stim[oddball["microcircuit1-ecell1"][times_stim] .== value_stim_f1]
			
# 		true_times_stim_f2 = times_stim[oddball["microcircuit1-ecell2"][times_stim] .== value_stim_f1]
	#@time full_time_dynamics(lc,sim,oddball)
	
			return d
			
		end
end

# ╔═╡ 17177775-89fe-43cf-a454-8919dbf5065d
param = Dict(
	"f1" => 0.8,
	"value_stim_f1" => 0.2,
	"Tinter" => 1.0,
	"Tstim" => 0.5,
	"Tfin" => 10.0,
	"value_stim_f2" => 0.0
)

# ╔═╡ 3c2b0efc-df33-4a5f-9128-bcfff869dbba
d=generate_oddball_task(param)

# ╔═╡ 00216ed7-8797-4ba7-adfe-4ac59be8c26e
scatter(d["microcircuit1-ecell1"])

# ╔═╡ b69bad9d-300a-42f9-997b-fd022a2d43d6
scatter(d["current-to-e1"])

# ╔═╡ c94a3def-357f-4692-9344-43429c0672d0
data,file = produce_or_load(datadir("sims","script1"), param, generate_oddball_task,prefix="oddball_task")

# ╔═╡ bbcdb7e3-e0f0-4b5c-acab-4364b857ab05
@unpack data

# ╔═╡ 4ba2d264-90da-4fe0-aa5b-44bc74647a6d


# ╔═╡ 31eef824-0da0-402e-8f54-d2b8aedb81bd
begin
	fig = plot_local_circuit(lc,sim,oddball)
	fig
end

# ╔═╡ Cell order:
# ╠═b2a49532-eb00-11eb-21eb-2f173ae5d57b
# ╠═a0fa119b-08b1-44f8-b253-71fccf78f2cc
# ╠═d436dc2b-107b-4a3b-a1c6-2c5969fc7714
# ╠═371793e6-9bfe-4daf-a26b-0aa89c82b167
# ╠═f02a080a-e259-4ec8-83af-ec896297ba8b
# ╠═08f3d0a8-7897-4c3f-9815-12f5d83eb688
# ╠═2bec0a29-96e0-4f5b-9324-d5429ffbe713
# ╠═17177775-89fe-43cf-a454-8919dbf5065d
# ╠═3c2b0efc-df33-4a5f-9128-bcfff869dbba
# ╠═00216ed7-8797-4ba7-adfe-4ac59be8c26e
# ╠═b69bad9d-300a-42f9-997b-fd022a2d43d6
# ╠═c94a3def-357f-4692-9344-43429c0672d0
# ╠═bbcdb7e3-e0f0-4b5c-acab-4364b857ab05
# ╠═4ba2d264-90da-4fe0-aa5b-44bc74647a6d
# ╠═31eef824-0da0-402e-8f54-d2b8aedb81bd
