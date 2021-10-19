# scp command to copy the results
# 

using DrWatson
@quickactivate "1-project-categorical-MMN"
push!(LOAD_PATH, srcdir())

using Distributed

addprocs(3)

@everywhere using Distributed, DrWatson
@everywhere @quickactivate "1-project-categorical-MMN"
@everywhere push!(LOAD_PATH, srcdir())
#script that computes the data for the phase diagram of ISI and integrator time constant

@everywhere	using MyNeurosciencePackage
@everywhere	using Parameters, JLD2, Statistics


    

    	# needs to construct a general array with all the stim. The question is, how to make it automatic??? Symbolic Julia?
        @everywhere	function generate_oddball_task(param::Dict{String,Float64})# f1::Float64, value_stim_f1::Float64; Tinter = 1.0, Tstim = 0.5, Tfin = 60.0, value_stim_f2 = 0.0)
	
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



# function map_from_dict(dict::Dict{String,Float64})
#     produce_or_load(datadir("sims", "script3"), dict, generate_oddball_task, prefix="oddball_task")
# end
    

@everywhere function my_own_map_from_dict(dict::Dict{String,Float64})
    temp = generate_oddball_task(dict)
    tempname = savename(dict)
    save("/scratch/kb3856/1-project-categorical-MMN/data/sims/script3/oddball_task_$tempname.jld2",temp)

end



general_args = Dict(
    "τ_integrator" => 20.0,#collect(0.1:0.2:3.0),
	"f1" => collect(0.5:0.1:0.9),#min was 0.5
	"value_stim_f1" => collect(0.1:0.1:0.6),
	"Tinter" => vcat(collect(0.2:0.3:15.0),collect(16.0:1.0:20.0)),#
	"Tstim" => 0.5,
	"Tfin" => 360.0,
	"value_stim_f2" => 0.0    
)
dicts = dict_list(general_args)
pmap(my_own_map_from_dict, dicts)


# general_args = Dict(
#     "τ_integrator" => collect(0.1:0.3:3.0) or 0.2,
# 	"f1" => collect(0.5:0.1:0.9),
# 	"value_stim_f1" => collect(0.1:0.1:0.6),
# 	"Tinter" => collect(0.2:0.3:15.0),
# 	"Tstim" => 0.5,
# 	"Tfin" => 300.0,
# 	"value_stim_f2" => 0.0    
# )