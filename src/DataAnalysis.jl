module DataAnalysis

using Statistics
using ..NeuronalStructures

export get_mean_firing_rate, compute_MMN_oddball, compute_MMN_time

function get_mean_firing_rate(data::Dict{String,Any},unit::String,current::String,value::Float64)
    # get 

    indices = data[current].==value
    value_mean = mean(data[unit][indices])

    return value_mean
end



function compute_MMN_oddball(lc::Vector{microcircuit{soma_PC,dend_sigmoid}}, sim::Dict{String,Float64})
    
    time_oddball = sim["initial_time"] +sim["Tinter"]+ (sim["number_repetitions"])*(sim["Tinter"]+sim["Tstim"]) + 0.025
    time_frequent = sim["initial_time"] +sim["Tinter"]+ (sim["number_repetitions"]-1)*(sim["Tinter"]+sim["Tstim"]) +0.025


    fr_frequent = lc[1].soma.r_save[round(Int,time_frequent/0.0005)]
    fr_oddball = lc[2].soma.r_save[round(Int,time_oddball/0.0005)]

    return fr_oddball - fr_frequent

end


function compute_MMN_time(lc::Vector{microcircuit{soma_PC,dend_sigmoid}}, sim::Dict{String,Float64})
    time_oddball_d = sim["initial_time"] +sim["Tinter"]+ (sim["number_repetitions"])*(sim["Tinter"]+sim["Tstim"]) 
    time_frequent_d = sim["initial_time"] +sim["Tinter"]+ (sim["number_repetitions"]-1)*(sim["Tinter"]+sim["Tstim"]) 

    time_oddball_f = sim["initial_time"] +sim["Tinter"]+ (sim["number_repetitions"])*(sim["Tinter"]+sim["Tstim"]) +sim["Tstim"]
    time_frequent_f = sim["initial_time"] +sim["Tinter"]+ (sim["number_repetitions"]-1)*(sim["Tinter"]+sim["Tstim"]) +sim["Tstim"]
    fr_frequent = lc[1].soma.r_save[round(Int,time_frequent_d/0.0005):round(Int,time_frequent_f/0.0005)]
    fr_oddball = lc[2].soma.r_save[round(Int,time_oddball_d/0.0005):round(Int,time_oddball_f/0.0005)]

    return fr_oddball .- fr_frequent

end


module bump_analysis





end
using .bump_analysis



end