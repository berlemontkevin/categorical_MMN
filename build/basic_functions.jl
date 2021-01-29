module Basic_functions

    export rect_linear
    export heaviside
    export f_I_Abott


function rect_linear(x::Float64)
    # this function takes max(0,x)
    temp = max(0,x)
    return temp
end

function heaviside(x::Float64)
    # compute the heaviside function
    if x>0.0
        temp = 1.0
    else
        temp = 0.0
    end
    return temp
end

function f_I_Abott(V::Float64,τ::Float64)
    # return the rate depending on thevoltage
    # Shape following Abott Chance 2005: parameters from Garcia 2017
    # Need to have the time constant of the neuron as argument
    Vth = -50
    Vτ = -60
    v = 1.0
    temp = 1.0 - exp(-(V-Vth)/v)
    temp2 = (V - Vth)/(τ*(Vth - Vτ))

    return temp2/temp
end


end