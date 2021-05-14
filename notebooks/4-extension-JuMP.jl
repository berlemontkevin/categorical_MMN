# this file is an extension of the notebook 4-adaptation theoretical

using JuMP


# test of JuMP module
using GLPK
model = Model(GLPK.Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@objective(model, Min, 12x + 20y)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
print(model)
optimize!(model)
@show termination_status(model)
@show primal_status(model)
@show dual_status(model)
@show objective_value(model)
@show value(x)
@show value(y)
@show shadow_price(c1)
@show shadow_price(c2)

# using NLP 

import NLopt

model = Model(NLopt.Optimizer)
set_optimizer_attribute(model, "algorithm", :LD_MMA)

a1 = 2
b1 = 0
a2 = -1
b2 = 1

@variable(m, x1)
@variable(m, x2 >= 0)

@NLobjective(m, Min, sqrt(x2))
@NLconstraint(m, x2 >= (a1*x1+b1)^3)
@NLconstraint(m, x2 >= (a2*x1+b2)^3)

set_start_value(x1, 1.234)
set_start_value(x2, 5.678)

optimize!(model)


@show value(x1)

println("got ", objective_value(model), " at ", [value(x1), value(x2)])


# Ipopt test
using Ipopt
rocket = Model(Ipopt.Optimizer)
set_silent(rocket)

h_0 = 1    # Initial height
v_0 = 0    # Initial velocity
m_0 = 1    # Initial mass
g_0 = 1    # Gravity at the surface

T_c = 3.5  # Used for thrust
h_c = 500  # Used for drag
v_c = 620  # Used for drag
m_c = 0.6  # Fraction of initial mass left at end

c     = 0.5 * sqrt(g_0 * h_0)  # Thrust-to-fuel mass
m_f   = m_c * m_0              # Final mass
D_c   = 0.5 * v_c * m_0 / g_0  # Drag scaling
T_max = T_c * g_0 * m_0        # Maximum thrust

n = 800    # Time steps
@variables(rocket, begin
    Δt ≥ 0, (start = 1/n) # Time step
    # State variables
    v[1:n] ≥ 0            # Velocity
    h[1:n] ≥ h_0          # Height
    m_f ≤ m[1:n] ≤ m_0    # Mass
    # Control variables
    0 ≤ T[1:n] ≤ T_max    # Thrust
end)
@objective(rocket, Max, h[n])

fix(v[1], v_0; force = true)
fix(h[1], h_0; force = true)
fix(m[1], m_0; force = true)
fix(m[n], m_f; force = true)

@NLexpressions(rocket, begin
    # Drag(h,v) = Dc v^2 exp( -hc * (h - h0) / h0 )
    drag[j = 1:n], D_c * (v[j]^2) * exp(-h_c * (h[j] - h_0) / h_0)
    # Grav(h)   = go * (h0 / h)^2
    grav[j = 1:n], g_0 * (h_0 / h[j])^2
    # Time of flight
    t_f, Δt * n
end)


for j in 2:n
    # h' = v
    # Rectangular integration
    # @NLconstraint(rocket, h[j] == h[j - 1] + Δt * v[j - 1])
    # Trapezoidal integration
    @NLconstraint(rocket, h[j] == h[j - 1] + 0.5 * Δt * (v[j] + v[j - 1]))
    # v' = (T-D(h,v))/m - g(h)
    # Rectangular integration
    # @NLconstraint(
    #     rocket,
    #     v[j] == v[j - 1] + Δt *((T[j - 1] - drag[j - 1]) / m[j - 1] - grav[j - 1])
    # )
    # Trapezoidal integration
    @NLconstraint(
        rocket,
        v[j] == v[j-1] +
            0.5 * Δt * (
                (T[j] - drag[j] - m[j] * grav[j]) / m[j] +
                (T[j - 1] - drag[j - 1] - m[j - 1] * grav[j - 1]) / m[j - 1]
            )
    )
    # m' = -T/c
    # Rectangular integration
    # @NLconstraint(rocket, m[j] == m[j - 1] - Δt * T[j - 1] / c)
    # Trapezoidal integration
    @NLconstraint(rocket, m[j] == m[j - 1] - 0.5 * Δt * (T[j] + T[j-1]) / c)
end

println("Solving...")
status = optimize!(rocket)

println("Max height: ", objective_value(rocket))


############################
# Set up of the model of the local microcircuit
############################
using JuMP
using Ipopt
function f(a)
    if a > 0.0
        return a
    else
        return 0.0
    end

end

function sigmoid(x)
    return 1.0./(1.0.+ exp.(-x))
end

function g(Iexc,Iinh)
    c1 = 0.2
    c2 =  -7.0
    c3 = -0.482
    c4 = 0.00964
    c5 = 0.11
    c6 = 0.1
    y = (Iexc - c2*Iinh + c6 )/(c3*Iinh + c4)
    
    return  c1*(-0.5 + sigmoid(y))+c5 
end


function test_opti(Istim::Float64)
    ### constant values
    τ_AMPA = 0.002
    τ_GABA = 0.005
    τ_int = 0.8
    τ_GABAdend = 0.01
    τ_NMDA = 0.06

    r0PV = -95.0
    r0SST = -33.0
    r0VIP = -33.0

    cPV = 330.0
    cSST = 132.0
    cVIP = 132.0

    c1 = 0.2
    c2 =  -7.0
    c3 = -0.482
    c4 = 0.00964
    c5 = 0.11
    c6 = 0.1

    τadapt = 1.0

    Ibg_soma =  0.15
    Ibg_dend = 0.1
    Ibg_pv = 0.3
    Ibg_sst = 0.2
    I_bg_vip = 0.3
    Ibg_int=0.310

    a = 135.0
    b = 54.0
    d = 0.308

    γnmda = 0.641*2
    γgaba = 2.0


    g_ssttodend = -0.09
    g_pvtopc = -0.005
    g_pctopc = 0.18
    g_pctovip = 0.058
    g_ssttovip = -0.1
    g_pctosst = 0.0435
    g_viptosst = -0.05
    g_pctopv = 0.08435
    g_ssttopv = -0.17
    g_pvtopv = -0.18
    g_inttosst = 0.22
    g_pctoint = 0.15



    α = -1.0 +0.8
    ######################


    ####
    # variables declarations

    microcircuit = Model(Ipopt.Optimizer)

    @variable(microcircuit, Rsst >= 0.0)
    @variable(microcircuit, Rpv >= 0.0)
    @variable(microcircuit, Rvip >= 0.0)
    @variable(microcircuit, Rpc >= 0.0)
    @variable(microcircuit, Rint >= 0.0)

    @variable(microcircuit, sssttodend >= 0.0)
    @variable(microcircuit, spc >= 0.0)
    @variable(microcircuit, ssst >= 0.0)
    @variable(microcircuit, svip >= 0.0)
    @variable(microcircuit, spv >= 0.0)
    @variable(microcircuit, sint >= 0.0)

    @variable(microcircuit, Isoma)
    @variable(microcircuit, Idend )
    @variable(microcircuit, Iexc >= 0.0)
    @variable(microcircuit, Iinh <= 0.0)
    @variable(microcircuit, Ivip)
    @variable(microcircuit, Isst)
    @variable(microcircuit, Ipv)
    @variable(microcircuit,Iint)
    ############


    register(microcircuit, :f, 1, f; autodiff = true)
    register(microcircuit, :sigmoid, 1, sigmoid; autodiff = true)
    register(microcircuit, :g, 2, g; autodiff = true)
    ##############
    # Constraints for solving





    #s Constraints
    @NLconstraint(microcircuit, cssst, ssst == τ_GABA*γgaba*Rsst)
    @NLconstraint(microcircuit, csvip, svip == τ_GABA*γgaba*Rvip)
    @NLconstraint(microcircuit, cspv, spv == τ_GABA*γgaba*Rpv)
    @NLconstraint(microcircuit, cssstdend, sssttodend == τ_GABAdend*γgaba*Rsst)

    @NLconstraint(microcircuit, cspc, 0.0 == -spc/τ_NMDA + γnmda*(1.0 - spc)*Rpc)
    @NLconstraint(microcircuit, csint, 0.0 == -sint/τ_NMDA + γnmda*(1.0 - sint)*Rint)

    # firing rate constraints
    @NLconstraint(microcircuit, crsst, -Rsst + f(cSST*Isst + r0SST)==0.0)
    @NLconstraint(microcircuit, crpv, -Rpv + f(cPV*Ipv + r0PV)==0.0)
    @NLconstraint(microcircuit, crvip, -Rvip + f(cVIP*Ivip + r0VIP) ==0.0)
    @NLconstraint(microcircuit, csoma, -Rpc + (a*Isoma-b)/(1.0 -exp(-d*(a*Isoma -b))) == 0.0)
    @NLconstraint(microcircuit,cint, α*Rint + (a*Iint-b)/(1.0 -exp(-d*(a*Iint -b))) == 0.0)

    #current Constraints


    @NLconstraint(microcircuit, cIdend, Idend == g(Iexc,Iinh))
    @NLconstraint(microcircuit, cIexcdend, Iexc == Ibg_dend)
    @NLconstraint(microcircuit, cIinhdend,  Iinh == g_ssttodend*sssttodend)
    @NLconstraint(microcircuit, cIsst, Isst == Ibg_sst + g_viptosst*svip + g_inttosst*sint + g_pctosst*spc)
    @NLconstraint(microcircuit, cIpv, Ipv == Ibg_pv + g_pctopv*spc + g_ssttopv*ssst +g_pvtopv*spv)
    @NLconstraint(microcircuit, cIvip, Ivip == I_bg_vip + g_ssttovip*ssst +g_pctovip*spc)
    @NLconstraint(microcircuit, cIsoma, Isoma == Ibg_soma + Istim + Idend + g_pvtopc*spv +g_pctopc*spc)
    @NLconstraint(microcircuit, cIint, Iint == g_pctoint*spc + Ibg_int)



    ###########


    optimize!(microcircuit)

return value(Rpc),value(Rint),value(Rvip),value(Rsst),value(Rpv)

end 


Rpc,Rint,Rsst,Rvip,Rpv = test_opti(0.1)

test_opti(0.2)