using DrWatson
@quickactivate "1-project-categorical-MMN"
push!(LOAD_PATH, srcdir())
using Revise
#script that computes the data for the phase diagram of ISI and integrator time constant

	using MyNeurosciencePackage
	using Parameters, JLD2, JuMP, Statistics



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
                
    # d = Dict()
    # for c in lc
    #     for pop in [c.list_pv, c.list_sst, c.list_vip, c.list_soma, c.list_integrator]
    #             # for now dend is not saved
    #         for n in pop
    #             d[n.name] = n.r_save[times_stim]
                        
    #         end
    #     end
                    
    # end
                
    # d["current-to-e1"] = oddball["microcircuit1-ecell1"][times_stim]
    # d["current-to-e2"] = oddball["microcircuit2-ecell1"][times_stim]
    # d["times_step"] = times_stim
        # Need to save all firing rates, all times and currents
                
                
    # 	true_times_stim = times_stim[oddball["microcircuit1-ecell1"][times_stim] .== value_stim_f1]
                
    # 		true_times_stim_f2 = times_stim[oddball["microcircuit1-ecell2"][times_stim] .== value_stim_f1]
        # @time full_time_dynamics(lc,sim,oddball)
        
    return lc, sim, oddball
                
end



function map_from_dict(dict::Dict{String,Float64})
    produce_or_load(datadir("sims", "script3"), dict, generate_oddball_task, prefix="oddball_task")
end
    

general_args = Dict(
    "τ_integrator" => 0.1,
	"f1" => 0.4,
	"value_stim_f1" => 0.1,
	"Tinter" => 0.2,
	"Tstim" => 0.5,
	"Tfin" => 300.0,
	"value_stim_f2" => 0.0    
)
dicts = dict_list(general_args)
generate_oddball_task(general_args)

@time generate_oddball_task(general_args)

@code_warntype generate_oddball_task(general_args)


using BenchmarkTools

@benchmark generate_oddball_task(general_args)


#  1.132586 seconds (1.59 M allocations: 91.013 MiB, 2.09% gc time, 99.71% compilation time)
# 1.151863 seconds (11.20 M allocations: 448.816 MiB, 3.22% gc time)
# 1.116177 seconds (10.44 M allocations: 378.299 MiB, 4.55% gc time)
#  1.089667 seconds (10.10 M allocations: 324.466 MiB, 3.28% gc time)

using Traceur

lc, sim, oddball =  generate_oddball_task(general_args)

lc[1]

@trace sum_input!(lc[1].list_soma[1],1)

@trace sum_input!(lc[1].list_vip[1],1)


@benchmark     update_adaptation!(lc[1].list_soma[1])
@benchmark     sum_input!(lc[1].list_soma[1],1)
@benchmark current_to_frequency(lc[1].list_soma[1])
@benchmark update_dend!(lc[1].list_dend,10)
@benchmark update_firing!(lc[1].list_sst[1])
@benchmark update_syn!(lc[1].list_sst[1])
@benchmark synapse_derivative(lc[1].list_vip[1].list_syn_post_nmda[1])

function temp(l_c::Vector{microcircuit}, sim::simulation_parameters, d::Dict{String,Vector{Float64}}, index::Int64)
         time_step(l_c, sim, oddball, 10)
        #  full_time_dynamics(l_c, sim, d)
   
end
temp(lc)
@benchmark temp(lc, sim, oddball, 10)


synapse_derivative(lc[1].list_vip[1].list_syn_post_gaba[1])
lc[1].list_vip[1].list_syn_post_gaba[1].dynamique_variables
@code_warntype synapse_derivative(lc[1].list_vip[1].list_syn_post_gaba[1])
@trace synapse_derivative(lc[1].list_vip[1].list_syn_post_gaba[1])


@trace update_adaptation!(lc[1].list_soma[1], lc[1].list_soma[1].dt)

@trace current_to_frequency(lc[1].list_sst[1])
@code_warntype current_to_frequency(lc[1].list_sst[1])


@trace update_dend!(lc[1].list_dend[1],1)

@trace update_dend!(lc[1].list_dend[1])

@trace update_firing!(lc[1].list_sst[1])
@code_warntype update_firing!(lc[1].list_sst[1])


@code_warntype synapse_derivative(lc[1].list_dend[1].list_syn_post_nmda[1])
@trace synapse_derivative(lc[1].list_dend[1].list_syn_post_nmda[1])
@trace synapse_derivative(lc[1].list_soma[1].list_syn_post_gaba[1])


@trace current_synapses!(lc[1].list_dend[1])
@code_warntype current_synapses!(lc[1].list_dend[1])

@code_warntype update_syn!(lc[1].list_soma[1])
@trace update_syn!(lc[1].list_soma[1])


@trace current_synapses!(lc[1].list_dend[1],lc[1].list_dend[1].list_syn)

@trace sigmoid(1.0)

@trace exp(1.0)

@trace generate_oddball_task(general_args)


@code_warntype current_synapses!(lc[1].list_dend[1],lc[1].list_dend[1].list_syn)





mutable struct dyn_variables
    A::Float64
    B::Float64
    C::Float64
 end

 function add_term(s::dyn_variables)
   s.A= 0.0
 end

 let s = dyn_variables(1.0, 2.0, 3.0)
     @benchmark add_term($s)
 end