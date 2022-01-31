using DrWatson
@quickactivate "1-project-categorical-MMN"
push!(LOAD_PATH, srcdir())

#script that computes the data for the phase diagram of ISI and integrator time constant

	using MyNeurosciencePackage
	using Parameters, JLD2, JuMP, Statistics


# function construct_one_local_microcircuit_integrator!(c::microcircuit; dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation=false, adaptation=false, depression=false, td_to_vip=true,integrator_tc = 0.8, time_tot = 1000)
#         # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name
    
#     vip1 = vip_cell(name=string(c.name, "-", "vipcell1"), Ibg=0.25, OU_process = OU_process(noise=zeros(time_tot)))
#     sst1 = sst_cell(name=string(c.name, "-", "sstcell1"), Ibg=0.25, adaptation_boolean=true), OU_process = OU_process(noise=zeros(time_tot))
#     pv1 = pv_cell(name=string(c.name, "-", "pvcell1"), Ibg=0.29, OU_process = OU_process(noise=zeros(time_tot)))
#     dend1 = dend_sigmoid(param_c=dend_param, name=string(c.name, "-", "dend1"), OU_process = OU_process(noise=zeros(time_tot)))
#     E1 = soma_PC(den=dend1, adaptation_boolean=adaptation, name=string(c.name, "-", "ecell1"), OU_process = OU_process(noise=zeros(time_tot)))
#     integrator1 = neural_integrator(τ = integrator_tc, name=string(c.name, "-", "integrator1"), OU_process = OU_process(noise=zeros(time_tot)))
        
#     create_process!(vip1.OU_process)
#     create_process!(sst1.OU_process)
#     create_process!(pv1.OU_process)
#     create_process!(dend1.OU_process)
#     create_process!(E1.OU_process)
#     create_process!(integrator1.OU_process)

#     push!(dend1.list_syn, gaba_syn(τ=10 * 0.001, neuron_pre=sst1, neuron_post=dend1, g=-0.09, name=string(c.name, "-", "sst1-to-dend1")))# 0.09
#     push!(E1.list_syn, gaba_syn(neuron_pre=pv1, neuron_post=E1, g=-0.001, depression=true, name=string(c.name, "-", "pv1-to-ecell1")))# -0.001
#     push!(E1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=E1, g=0.18, name=string(c.name, "-", "ecell1-to-ecell1")))
        
#     push!(vip1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=vip1, g=0.058, depression=true, name=string(c.name, "-", "ecell1-to-vip1")))
#     push!(vip1.list_syn, gaba_syn(neuron_pre=sst1, neuron_post=vip1, g=-0.1, facilitation=true, name=string(c.name, "-", "sst1-to-vip1")))
        
#     push!(sst1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=sst1, g=0.0435, facilitation=true, name=string(c.name, "-", "ecell1-to-sst1")))
#     push!(sst1.list_syn, gaba_syn(neuron_pre=vip1, neuron_post=sst1, g=-0.05, facilitation=true, name=string(c.name, "-", "ss1-to-dend1")))
        
#     push!(pv1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=pv1, g=0.04435, depression=true, name=string(c.name, "-", "ecell1-to-pv1")))# 0.0435
#     push!(pv1.list_syn, gaba_syn(neuron_pre=sst1, neuron_post=pv1, g=-0.17, name=string(c.name, "-", "sst1-to-pv1")))# -0.17
#     push!(pv1.list_syn, gaba_syn(neuron_pre=pv1, neuron_post=pv1, g=-0.18, name=string(c.name, "-", "pv1-to-pv1")))
        
           
       
#     push!(vip1.list_syn, nmda_syn(neuron_pre=integrator1, neuron_post=vip1, g=0.47, depression=depression, name=string(c.name, "-", "integrator1-to-vip1")))
#     push!(pv1.list_syn, nmda_syn(neuron_pre=integrator1, neuron_post=pv1, g=0.31, depression=depression, name=string(c.name, "-", "integrator1-to-pv1")))
        
            
#     push!(sst1.list_syn, nmda_syn(neuron_pre=integrator1, neuron_post=sst1, g=0.22, facilitation=facilitation, name=string(c.name, "-", "integrator1-to-sst1")))
            
#     push!(integrator1.list_syn, nmda_syn(neuron_pre=E1, neuron_post=integrator1, g=0.15, name=string(c.name, "-", "ecell1-to-integrator1")))
    
#     push!(dend1.list_syn, nmda_syn(neuron_pre=integrator1, neuron_post=dend1, g=0.4, depression=true, name=string(c.name, "-", "integrator1-to-dend1")))	
            
            
#     push!(c.list_dend, dend1)
      
#     push!(c.list_soma, E1)
       
#     push!(c.list_vip, vip1)
      
#     push!(c.list_sst, sst1)
      
#     push!(c.list_pv, pv1)
    
#     push!(c.list_integrator, integrator1)
    
    
    
# end
    
#     # TODO add a list of microcircuit as argument to the dynamics functions
    
# function construct_two_local_microcircuit_integrator(; dend_param=dendrites_param_sigmoid(0.12, -7.0, -0.482, 0.00964, 0.19624, 0.0), facilitation=false, adaptation=false, depression=false, td_to_vip=true, τ_integrator = 0.8, time_tot = 1000)
#         # once the network will be setup with more microcuicuit. It will be necessary to add the microcircuit in the name
    
    
#     c1 = microcircuit(name="microcircuit1")
#     c2 = microcircuit(name="microcircuit2")
        
#     construct_one_local_microcircuit_integrator!(c1; dend_param, facilitation, adaptation, td_to_vip, depression, integrator_tc = τ_integrator, time_tot = time_tot)
#     construct_one_local_microcircuit_integrator!(c2; dend_param, facilitation, adaptation, td_to_vip, depression, integrator_tc = τ_integrator, time_tot = time_tot)
    
    
#     push!(c1.list_sst[1].list_syn, nmda_syn(neuron_pre=c2.list_soma[1], neuron_post=c1.list_sst[1], g=0.0435, facilitation=true, name=string("ecell2-to-sst1")))
#     push!(c2.list_sst[1].list_syn, nmda_syn(neuron_pre=c1.list_soma[1], neuron_post=c2.list_sst[1], g=0.0435, facilitation=true, name=string("ecell1-to-sst2")))
    
#     push!(c1.list_dend[1].list_syn, nmda_syn(neuron_pre=c2.list_integrator[1], neuron_post=c1.list_dend[1], g=0.1, depression=true, name=string("integrator2-to-ecell1")))
#     push!(c2.list_dend[1].list_syn, nmda_syn(neuron_pre=c1.list_integrator[1], neuron_post=c2.list_dend[1], g=0.1, depression=true, name=string("integrator1-to-ecell2")))
    
            
            
            
#     list_microcircuit = [c1,c2]
    
#     return list_microcircuit
# end
    

    	# needs to construct a general array with all the stim. The question is, how to make it automatic??? Symbolic Julia?
	function generate_oddball_task(param::Dict{String,Float64})# f1::Float64, value_stim_f1::Float64; Tinter = 1.0, Tstim = 0.5, Tfin = 60.0, value_stim_f2 = 0.0)
	
    @unpack f1, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator = param
            
    f2 = 1 - f1
    dt = 0.0005
    t_tot = Int64(Tfin / dt)
    t_tot += 2000 # Adding a resting state time

            
            # create the network with all adapt, depression and so on (TODO the main fig)
    lc = construct_two_local_microcircuit_integrator(facilitation=true, depression=true, adaptation=true,  τ_integrator = τ_integrator ,time_tot = t_tot)
            
            
    sim = simulation_parameters(Tstimduration=Tstim, TISI=Tinter, Tfin=Tfin, Tinit=2000)
        
    
            
            # For now dictionnary is define manually
    oddball = Dict(
    "microcircuit1-vipcell1" => zeros(t_tot), 
    "microcircuit1-sstcell1" => zeros(t_tot),
    "microcircuit1-pvcell1" => zeros(t_tot),
    "microcircuit1-ecell1" => zeros(t_tot),
    "microcircuit1-integrator1" => zeros(t_tot),
    "microcircuit2-vipcell1" => zeros(t_tot),
    "microcircuit2-sstcell1" => zeros(t_tot),
    "microcircuit2-pvcell1" => zeros(t_tot),
    "microcircuit2-ecell1" => zeros(t_tot),
    "microcircuit2-integrator1" => zeros(t_tot)
    )
        
    array_stim_f1 = zeros(round(Int, Tfin / dt))
    array_stim_f2 = zeros(round(Int, Tfin / dt))
    
    nbr_stim = round(Int, Tfin / (Tinter + Tstim))
    temp_index_stim = 0
    while temp_index_stim < nbr_stim - 1
            
        if rand() < f1
            #print(round(Int, (temp_index_stim + 1) * (Tinter + Tstim) / dt))
            array_stim_f1[round(Int, 1 + temp_index_stim * (Tinter + Tstim) / dt):round(Int, (temp_index_stim + 1) * (Tinter + Tstim) / dt)] = [zeros(round(Int, Tinter / dt)) ; value_stim_f1 * ones(round(Int, Tstim / dt))]
                
            array_stim_f2[round(Int, 1 + temp_index_stim * (Tinter + Tstim) / dt):round(Int, (temp_index_stim + 1) * (Tinter + Tstim) / dt)] = [zeros(round(Int, Tinter / dt)) ; value_stim_f2 * ones(round(Int, Tstim / dt))]
        else
                
            array_stim_f1[round(Int, 1 + temp_index_stim * (Tinter + Tstim) / dt):round(Int, (temp_index_stim + 1) * (Tinter + Tstim) / dt)] = [zeros(round(Int, Tinter / dt)) ; value_stim_f2 * ones(round(Int, Tstim / dt))]
                
            array_stim_f2[round(Int, 1 + temp_index_stim * (Tinter + Tstim) / dt):round(Int, (temp_index_stim + 1) * (Tinter + Tstim) / dt)] = [zeros(round(Int, Tinter / dt)) ; value_stim_f1 * ones(round(Int, Tstim / dt))]
        end
                
        temp_index_stim += 1
                    
    end
        # for now, the end of the array is 0
        
        
        
        # of course it will be necessary to transform it in a function
        
    oddball["microcircuit1-ecell1"] = [0.1 * ones(2000);array_stim_f1]
    oddball["microcircuit2-ecell1"] = [0.1 * ones(2000);array_stim_f2]
        
    full_time_dynamics(lc, sim, oddball)
        
    times_stim = collect((2000 + round(Int,(Tinter + 0.01) / dt)):round(Int,(Tstim + Tinter) / dt):round(Int,Tfin / dt + 2000))
                
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
        # Need to save all firing rates, all times and currents
                
                
    # 	true_times_stim = times_stim[oddball["microcircuit1-ecell1"][times_stim] .== value_stim_f1]
                
    # 		true_times_stim_f2 = times_stim[oddball["microcircuit1-ecell2"][times_stim] .== value_stim_f1]
        # @time full_time_dynamics(lc,sim,oddball)
        
    return d
                
end



function map_from_dict(dict::Dict{String,Float64})
    produce_or_load(datadir("sims", "script3"), dict, generate_oddball_task, prefix="oddball_task")
end
    

general_args = Dict(
    "τ_integrator" => collect(0.1:3.0:16.0),#collect(0.1:0.2:3.0),
	"f1" => collect(0.5:0.1:0.9),#min was 0.5
	"value_stim_f1" => collect(0.1:0.1:0.6),
	"Tinter" => vcat(collect(0.2:0.3:15.0),collect(16.0:1.0:20.0)),#
	"Tstim" => 0.5,
	"Tfin" => 360.0,
	"value_stim_f2" => 0.0    
)
dicts = dict_list(general_args)
map(map_from_dict, dicts)


# general_args = Dict(
#     "τ_integrator" => collect(0.1:0.3:3.0) or 0.2,
# 	"f1" => collect(0.5:0.1:0.9),
# 	"value_stim_f1" => collect(0.1:0.1:0.6),
# 	"Tinter" => collect(0.2:0.3:15.0),
# 	"Tstim" => 0.5,
# 	"Tfin" => 300.0,
# 	"value_stim_f2" => 0.0    
# )