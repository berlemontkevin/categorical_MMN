### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ a2d39c60-d06a-11eb-1364-c1e29d8f626f
begin
	using DrWatson
	
	quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
	push!(LOAD_PATH,srcdir())

end

# ╔═╡ 1334881e-d06b-11eb-3ec8-13f9d3037c39
begin
	using Revise
	using MyNeurosciencePackage
	using CairoMakie
	using Parameters, JLD2, JuMP, Statistics
end

# ╔═╡ e5810c90-d06b-11eb-26d3-e378a54e0bf7
md""" # Setup of the network

- Adaptation to pyramidal cells
- E to SST should be facilitated (Lee et al 2013) (Long range)
- E to VIP long range are depression (Karnani et al 2016)
- Int to dendrite = depression
- E to SST local should be facilitated
- VIP and SST facilitating
- PV to E depression
- E to PV depression
- E to VIP depression

"""

# ╔═╡ 16aa4df0-d06b-11eb-3d89-6fafd70abb11
begin
	

# #TODO add name field to microcircuit
# function construct_one_local_microcircuit_integrator!(c::microcircuit; dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation = false, adaptation = false, depression = false, td_to_vip = true)
#     # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name

#     vip1 = vip_cell(name = string(c.name, "-","vipcell1"),Ibg = 0.25)
#     sst1 = sst_cell(name = string(c.name, "-","sstcell1"), Ibg = 0.25,adaptation_boolean = true)
#     pv1 = pv_cell(name = string(c.name, "-","pvcell1"), Ibg = 0.29)
#     dend1 = dend_sigmoid(param_c = dend_param,name = string(c.name, "-","dend1"))
#     E1 = soma_PC(den=dend1, adaptation_boolean = adaptation,name = string(c.name, "-","ecell1"))
#     integrator1 = neural_integrator(name = string(c.name, "-","integrator1"))
    
    
#     push!(dend1.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst1,neuron_post=dend1, g = -0.09, name = string(c.name, "-","sst1-to-dend1")))#0.09
#     push!(E1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E1, g = -0.001,depression = true, name = string(c.name, "-","pv1-to-ecell1")))#-0.001
#     push!(E1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=E1, g = 0.18, name = string(c.name, "-","ecell1-to-ecell1")))
    
#     push!(vip1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=vip1, g = 0.058, depression = true, name = string(c.name, "-","ecell1-to-vip1")))
#     push!(vip1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=vip1, g = -0.1, facilitation = true, name = string(c.name, "-","sst1-to-vip1")))
    
#     push!(sst1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst1, g = 0.0435, facilitation = true, name = string(c.name, "-","ecell1-to-sst1")))
#     push!(sst1.list_syn, gaba_syn(neuron_pre=vip1,neuron_post=sst1, g = -0.05, facilitation = true, name = string(c.name, "-","ss1-to-dend1")))
    
#     push!(pv1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=pv1, g = 0.04435, depression = true, name = string(c.name, "-","ecell1-to-pv1")))#0.0435
#     push!(pv1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=pv1, g = -0.17, name = string(c.name, "-","sst1-to-pv1")))#-0.17
#     push!(pv1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=pv1, g = -0.18, name = string(c.name, "-","pv1-to-pv1")))
    
       
   
#         push!(vip1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=vip1, g = 0.47, depression = depression, name = string(c.name, "-","integrator1-to-vip1")))
#     push!(pv1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=pv1, g = 0.31, depression = depression, name = string(c.name, "-","integrator1-to-pv1")))
    
		
# 		push!(sst1.list_syn,nmda_syn(neuron_pre=integrator1,neuron_post=sst1, g = 0.22, facilitation = facilitation, name = string(c.name, "-","integrator1-to-sst1")))
		
#     push!(integrator1.list_syn, nmda_syn(neuron_pre = E1, neuron_post = integrator1, g = 0.15, name = string(c.name, "-","ecell1-to-integrator1")))

# 	    push!(dend1.list_syn, nmda_syn(neuron_pre = integrator1, neuron_post = dend1, g = 0.4,depression = true, name = string(c.name, "-","integrator1-to-dend1")))	
		
		
#     push!(c.list_dend,dend1)
  
#     push!(c.list_soma,E1)
   
#     push!(c.list_vip,vip1)
  
#     push!(c.list_sst,sst1)
  
#     push!(c.list_pv,pv1)

#     push!(c.list_integrator,integrator1)



# end

# #TODO add a list of microcircuit as argument to the dynamics functions

# function construct_two_local_microcircuit_integrator(; dend_param = dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation = false, adaptation = false, depression = false, td_to_vip = true)
#     # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name


#     c1 = microcircuit(name="microcircuit1")
#     c2 = microcircuit(name = "microcircuit2")
    
#     construct_one_local_microcircuit_integrator!(c1; dend_param, facilitation, adaptation, td_to_vip, depression)
#     construct_one_local_microcircuit_integrator!(c2; dend_param, facilitation, adaptation, td_to_vip, depression)


#     push!(c1.list_sst[1].list_syn, nmda_syn(neuron_pre=c2.list_soma[1],neuron_post=c1.list_sst[1], g = 0.0435, facilitation = true, name = string("ecell2-to-sst1")))
#     push!(c2.list_sst[1].list_syn, nmda_syn(neuron_pre=c1.list_soma[1],neuron_post=c2.list_sst[1], g = 0.0435, facilitation = true, name = string("ecell1-to-sst2")))

#    push!(c1.list_dend[1].list_syn, nmda_syn(neuron_pre=c2.list_integrator[1],neuron_post=c1.list_dend[1], g = 0.1, depression = true, name = string("integrator2-to-ecell1")))
#      push!(c2.list_dend[1].list_syn, nmda_syn(neuron_pre=c1.list_integrator[1],neuron_post=c2.list_dend[1], g = 0.1, depression = true, name = string("integrator1-to-ecell2")))

		
		
		
#     list_microcircuit = [c1,c2]

#     return list_microcircuit
# end

	
	
end

# ╔═╡ 4f7a9670-077c-488c-be66-264db4616673
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

# ╔═╡ 970c91ff-ad74-4d1a-9f5a-c62bf0a3dd27


# ╔═╡ 316e1a4b-3490-4c82-ba75-e2741ffccaf7
begin
	
	#needs to construct a general array with all the stim. The question is, how to make it automatic??? Symbolic Julia?
	
	f1 = 0.2
	f2 = 1-f1
	value_stim_f1 = 0.2
	value_stim_f2 = 0.0
	Tfin = 80.0
	dt = 0.0005
	Tinter = 1.0
	Tstim = 0.5
	t_tot = Int64(Tfin/dt)
	lc = construct_two_local_microcircuit_integrator(facilitation=true, depression = true, adaptation = true)
	sim = simulation_parameters(Tfin=10.0)

	sim.Tfin = Tfin
	t_tot +=2000
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
	
	
end

# ╔═╡ 570c0cd7-8469-448a-a6eb-5accefed980e
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

# ╔═╡ 368c79d3-707c-4847-a7f0-a87758447d97
sim

# ╔═╡ 7c4c482e-2397-476e-925a-1bbc295d2ec7
begin
	fig = plot_local_circuit(lc,sim,oddball)
	fig
end

# ╔═╡ a50b968f-df2a-4388-96dd-3b2c26296038
md""" ## Plan to analyze the data


The idea would be to extract the firing rate 10ms after the stim (for example)
Plotting the mean firing rate w/r the frequency


"""

# ╔═╡ 338dfbbf-3701-467a-8eb5-0f77ca7bfd93
begin
	
	times_stim = 4100:3000:160000
	
	true_times_stim = times_stim[oddball["microcircuit1-ecell1"][times_stim] .== 0.2]
	
	
	mean(lc[1].list_soma[1].r[true_times_stim])
	
	
	scatter(true_times_stim,lc[1].list_soma[1].r[true_times_stim])
	
end

# ╔═╡ 173b5c1a-c05d-4c3c-81da-1680fda339ed
begin
	
	#needs to construct a general array with all the stim. The question is, how to make it automatic??? Symbolic Julia?
	function test_freq()
		fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98)
            )
    ax1 = fig[1, 1] = Axis(fig, title = "Local circuit 1")
	list_mean = Float64[]	
	ax2 = fig[2, 1] = Axis(fig, title = "Mean of local circuit 1")
	
		for value_stim_f1 in 0.1:0.1:0.9
	for f1 in 0.1:0.1:0.9
	f2 = 1-f1
	#value_stim_f1 = 0.2
	value_stim_f2 = 0.0
	Tfin = 120.0
	dt = 0.0005
	Tinter = 1.0
	Tstim = 0.5
	t_tot = Int64(Tfin/dt)
	lc = construct_two_local_microcircuit_integrator(facilitation=true, depression = true, adaptation = true)
	sim = simulation_parameters(Tfin=10.0)

	sim.Tfin = Tfin
	t_tot +=2000
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
	
	times_stim = 4100:3000:200000
			
	
	true_times_stim = times_stim[oddball["microcircuit1-ecell1"][times_stim] .== value_stim_f1]
	
	
	#mean(lc[1].list_soma[1].r[true_times_stim])
	
	
	scatter!(ax1,true_times_stim.*sim.dt,lc[1].list_soma[1].r[true_times_stim],label = "Frequency = $f1",color=(:blue, f1))
			
			push!(list_mean,mean(lc[1].list_soma[1].r[true_times_stim]))
		end
		axislegend(ax1)
		
		#lines!(ax2, 0.1:0.1:0.9,list_mean, linewidth=2)
		end
			return fig,list_mean
	end
end

# ╔═╡ 509551a9-b49b-4fbe-9995-3b7c1b500be4
fig_test,test_mean = test_freq()

# ╔═╡ 3a546d30-5343-4517-8f9f-db408f436eff
a = Dict("List mean r"=> test_mean)

# ╔═╡ af4e8f2d-c209-4f13-b85a-adb0d2cf880a
save(datadir("sims","temp","freq_stim_list.jld2"),a)

# ╔═╡ f26ea19f-ee30-4c52-b97d-50d5df9cd58b
scatter(test_mean)

# ╔═╡ Cell order:
# ╠═a2d39c60-d06a-11eb-1364-c1e29d8f626f
# ╠═1334881e-d06b-11eb-3ec8-13f9d3037c39
# ╠═e5810c90-d06b-11eb-26d3-e378a54e0bf7
# ╠═16aa4df0-d06b-11eb-3d89-6fafd70abb11
# ╠═570c0cd7-8469-448a-a6eb-5accefed980e
# ╠═4f7a9670-077c-488c-be66-264db4616673
# ╠═368c79d3-707c-4847-a7f0-a87758447d97
# ╠═970c91ff-ad74-4d1a-9f5a-c62bf0a3dd27
# ╠═316e1a4b-3490-4c82-ba75-e2741ffccaf7
# ╠═7c4c482e-2397-476e-925a-1bbc295d2ec7
# ╠═a50b968f-df2a-4388-96dd-3b2c26296038
# ╠═338dfbbf-3701-467a-8eb5-0f77ca7bfd93
# ╠═173b5c1a-c05d-4c3c-81da-1680fda339ed
# ╠═509551a9-b49b-4fbe-9995-3b7c1b500be4
# ╠═3a546d30-5343-4517-8f9f-db408f436eff
# ╠═af4e8f2d-c209-4f13-b85a-adb0d2cf880a
# ╠═f26ea19f-ee30-4c52-b97d-50d5df9cd58b
