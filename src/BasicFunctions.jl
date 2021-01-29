module BasicFunctions

    export rect_linear
    export heaviside
    export f_I_Abott

"""
    Function: rect_linear

    This function return 0 or x depending on the sign of x
"""
function rect_linear(x::Float64)
    temp = max(0, x)
    return temp
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


end