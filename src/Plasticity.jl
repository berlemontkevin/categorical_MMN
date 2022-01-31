module Plasticity

using SpecialFunctions
using QuadGK
using Cubature

const γ=0.57721

function kappa(f)
    # compute the kappa function of Graupner and brunel 2012
    return exp(-f*γ)/gamma(f)

end


function shot_noise(c,f)
    # recursive function to compute the integral
    if c<1
        return kappa(f)*(c^(f-1))

    else
        a,b=hquadrature(t -> shot_noise(t-1,f)*t^(-f),1.0,c,abstol=1e-8)
        return c^(f-1)*(kappa(f) - f * a)
    end

end


function cum_int_shot_noise(Cpre,Cpost,Ctot,fpre,fpost)

a,b = quadgk(t -> shot_noise(t/Cpre,fpre)*shot_noise((Ctot - t)/Cpost,fpost),0.0,Ctot)

return a
end


function time_alpha(θx,fpre,fpost,Cpre,Cpost)
# compute the fraciton of time above a thrshold
# voir si la borne max est a changee
return quadgk(t -> cum_int_shot_noise(Cpre,Cpost,t,fpre,fpost),θx,10.0)

end


function up_down_weights(θd, θp, fpre, fpost, Cpre, Cpost, γp, γd, σ_carre,τ)
    # compute up and down weights

αd = time_alpha(θd,τ*fpre,τ*fpost,Cpre,Cpost)
αp = time_alpha(θp,τ*fpre,τ*fpost,Cpre,Cpost)
Γp = γp * αp
Γd = γd * αd

ρbar = Γp/(Γp + Γd)

σp_carre = σ_carre * (αp + αd) / (Γp + Γd)

U = 0.5 * (1.0 + erf(- (0.5 - ρbar + ρbar * exp(-(Γp + Γd)/τ))/(sqrt(σp_carre* (1.0 - exp(-2.0*(Γp + Γd)/τ))))))

D = 0.5 * (1.0 - erf(- (0.5 - ρbar + (ρbar-1.0) * exp(-(Γp + Γd)/τ))/(sqrt(σp_carre* (1.0 - exp(-2.0*(Γp + Γd)/τ))))))

    return U,D

end



function int_shot_noise(c,f)
    a,b = quadgk(t -> shot_noise(t,f),c,10.0)

    return a
end




end


Cpre = 0.56
Cpost  = 1.24
θd = 1.0
θp = 1.3
γd = 331.9
γp = 725.08
σ_carre = 3.35*3.35
τ = 0.34636

shot_noise(7.5,1.0)


a,b = hquadrature(t -> shot_noise(t,1.0),0.0,5.0,abstol=1e-8)


int_shot_noise(1.0,1.0)




cum_int_shot_noise(Cpre,Cpost,1.0,10.0,10.0)




time_alpha(θd,10.0,10.0,Cpre,Cpost)

U,D = up_down_weights(θd, θp, 10.0, 10.0, Cpre, Cpost, γp, γd, σ_carre,τ)


module STDP


function shot_noise(f)
    list_shot = Float64[]
    list_c = Float64[]

    κf = exp(-f*γ)/gamma(f)

    for c in 0.0:0.001:0.999
        push!(list_shot,κf*(c^(f-1)))
        push!(list_c, c)
    end

    c=1.0
    while list_shot[end]>1e-5 || c <3.0

        temp = sum(0.001.*list_shot[1:length(list_c[1000:end])].*(list_c[1000:end] ).^(-f))

        push!(list_shot,c^(f-1)*(κf - f * temp))
        push!(list_c, c)
        c+=0.001
    end

    return list_shot,list_c
end



function cum_int_shot_noise(Cpre,Cpost,Ctot,fpre,fpost)

    a,b = quadgk(t -> shot_noise(t/Cpre,fpre)*shot_noise((Ctot - t)/Cpost,fpost),0.0,Ctot)
    
    return a
    end
    
    
    function time_alpha(θx,fpre,fpost,Cpre,Cpost)
    # compute the fraciton of time above a thrshold
    # voir si la borne max est a changee
    return quadgk(t -> cum_int_shot_noise(Cpre,Cpost,t,fpre,fpost),θx,10.0)
    
    end
    


end


a,b = shot_noise(10.0)

using Makie
using CairoMakie

figure, axis, scatterobject = scatter(b,a)
figure

d = [sum(a[Int(t):end].*0.001) for t in 1:length(b)]
figure, axis, scatterobject = scatter(b[1:10:end],d[1:10:end])

sum(a[1:end].*0.001)