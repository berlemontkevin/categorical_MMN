module DataAnalysis

using Statistics
export get_mean_firing_rate

function get_mean_firing_rate(data::Dict{String,Any},unit::String,current::String,value::Float64)
    # get 

    indices = data[current].==value
    value_mean = mean(data[unit][indices])

    return value_mean
end

# TODO: mean firing rate when stim is presented
# TODO: mena firing rate difference between two populaitons



end