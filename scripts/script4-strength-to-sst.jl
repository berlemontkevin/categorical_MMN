# scp command to get the data
# scp kb3856@greene.hpc.nyu.edu: /scratch/kb3856/1-project-categorical-MMN/data/sims/script4 /mnt/c/Users/kevin/OneDrive/1-Research_projects/1-project-categorical-MMN/data/sims

using DrWatson
@quickactivate "1-project-categorical-MMN"
push!(LOAD_PATH, srcdir())

using Distributed
addprocs(5)

@everywhere using DrWatson, Distributed

@everywhere @quickactivate "1-project-categorical-MMN"
@everywhere push!(LOAD_PATH, srcdir())


@everywhere using MyNeurosciencePackage
@everywhere using Parameters, JLD2, Statistics



@everywhere function run_simu_oddball(params::Dict{String,Float64})
    stim_1, stim_2 = create_deterministic_oddball(params)

@unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time, sst_strength = params

t_tot = length(stim_1)

td_to_interneurons = [0.47,0.31,0.22.*sst_strength]./sum([0.47,0.31,0.22.*sst_strength])


lc = construct_two_local_microcircuit_integrator_full_param(noise=false,  integrator_tc = τ_integrator ,time_tot = t_tot, top_down_to_interneurons = td_to_interneurons)
        
        
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

oddball["microcircuit1-ecell1"] = stim_1
oddball["microcircuit2-ecell1"] = stim_2
    

full_time_dynamics(lc, sim, oddball)
temp = compute_MMN_oddball(lc,params)

return Dict("MMN effect" => temp, "params" => params)
end




params_list = Dict(
    	"number_repetitions" => collect(1.0:1.0:10.0),
        "τ_integrator" => collect(0.1:0.4:5.0),
    	"value_stim_f1" => 0.2,
    	"Tinter" => collect(0.5:0.5:7.0),#
    	"Tstim" => 0.5,
    	"Tfin" => 1.0,
    	"value_stim_f2" => 0.0,    
    	"initial_time" => 5.0,
        "sst_strength" => collect(0.0:0.1:1.0)
    )



@everywhere function my_own_map_from_dict(dict::Dict{String,Float64})
        temp = run_simu_oddball(dict)
        tempname = savename(dict)
        #save("/scratch/kb3856/1-project-categorical-MMN/data/sims/script4/oddball_task_$tempname.jld2",temp)
        save(datadir("sims","script4","oddball_task_$tempname.jld2"),temp)
    end


dicts = dict_list(params_list)

pmap(my_own_map_from_dict, dicts)
