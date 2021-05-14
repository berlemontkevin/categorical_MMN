module BasicFunctions

    using Parameters 
    export rect_linear, sigmoid
    export heaviside
    export f_I_Abott
    export OU_process, update_process!

"""
    Function: rect_linear

    This function return 0 or x depending on the sign of x
"""
function rect_linear(x::Float64)
    temp = max(0, x)
    return temp
end


"""
    FUnction: sigmoid

    This function returns the sigmoid of x
"""
function sigmoid(x::Float64)
    return 1.0./(1.0.+ exp.(-x))
end

"""
    Function: heaviside(x::Float64)

Compute the heaviside function for x
"""
function heaviside(x::Float64)
    if x > 0.0
        temp = 1.0
    else
        temp = 0.0
    end
    return temp
end


"""
    Function: f_I_Abott(V::Float64, τ::Float64)

Return the rate fo the neuron depeinding on the voltage
Constructed following Abott Chance 2005
"""
function f_I_Abott(V::Float64, τ::Float64)
    Vth = -50
    Vτ = -60
    v = 1.0
    temp = 1.0 - exp(-(V - Vth) / v)
    temp2 = (V - Vth) / (τ * (Vth - Vτ))

    return temp2 / temp
end


@with_kw mutable struct OU_process
    τ::Float64 = 0.002
    dt::Float64 = 0.0005
    σ::Float64 = 5 * 0.001
    noise::Array{Float64} = [0.0]
end


function update_process!(s_process::OU_process)
    push!(s_process.noise, s_process.noise[end])
    s_process.noise[end] +=  s_process.dt / s_process.τ * (-s_process.noise[end]+ sqrt(s_process.τ*s_process.σ*s_process.σ)*randn())
   
end











end