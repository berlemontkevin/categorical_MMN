### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 30a76def-2a04-412d-89c9-383a106ed7ae
begin
	# test of JUMP module
	using JuMP
	
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
	
	
end

# ╔═╡ ff667110-b363-11eb-0f66-37c2f5cf3240
md""" #Theoretical study of adaptation and facilitation in a local circuit


The goal of this notebook is to study the impact of adaptation and facilitation on the local microcircuit. At first the model will be subject to constant input.




At the steady state the equations of the network are:

`` \tau_{AMPA} \frac{d R_{SST}}{dt} = -R + \left[ c I_{SST} + r^0\right]_+``

`` \tau_{AMPA} \frac{d R_{PV}}{dt} = -R + \left[ c I_{PV} + r^0\right]_+``

`` \tau_{AMPA} \frac{d R_{VIP}}{dt} = -R + \left[ c I_{VIP} + r^0\right]_+``

`` \tau_{AMPA} \frac{d R_{PC}}{dt} = - R_{PC} + \frac{a I_{PC} + b}{1 - \exp(-d (a I_{PC} - b))}``

`` I _{PC} = I_{soma} + I_{dend}``

`` \tau_{int} \frac{d r_{int}}{dt} = - \alpha r_{int} + I_{int} (t)``

`` I_{soma} = I_{stim} + g_{pvtopc} s_{pvtopc} + g_{pctopc} s_{pctopc}``

`` I_{pv} = I_{bgpv} + g_{pctopv} s_{pctopv} + g_{ssttopv} s_{ssttopv}``

`` I_{sst} = I_{bgsst} + g_{viptosst} s_{viptosst} + g_{inttosst} s_{inttosst}``

`` I_{vip} = I_{bgvip} + g_{ssttovip} s_{ssttovip} ``

`` I_{int} = g_{pctoint} s_{pctoint}``

`` I_{dend} = g(I_{exc},I_{inh})``


`` I_{exc,dend} = I_{bg,dend} ``

`` I_{inh,dend} = g_{ssttodend} s_{ssttodend}``

## GABA synapses

`` \frac{d s_{pvtopc}}{dt} = -s_{pvtopc}/\tau_{} + \gamma R_{PV}` ``

`` \frac{d s_{viptosst}}{dt} = -s_{viptosst}/\tau_{} + \gamma R_{vip}` ``

`` \frac{d s_{ssttovip}}{dt} = -s_{ssttovip}/\tau_{} + \gamma R_{sst}` ``

`` \frac{d s_{ssttodend}}{dt} = -s_{ssttodend}/\tau_{} + \gamma R_{sst}` ``

`` \frac{d s_{ssttopv}}{dt} = -s_{ssttopv}/\tau_{} + \gamma R_{sst}` ``



## NMDA synapses

`` \frac{d s_{pctoint}}{dt} = -s_{pctoint}/\tau_{} + \gamma (1 - s_{pctoint}) R_{PC}``

`` \frac{d s_{pctopc}}{dt} = -s_{pctopc}/\tau_{} + \gamma (1 - s_{pctopc}) R_{pc}` ``

`` \frac{d s_{pctosst}}{dt} = -s_{pctosst}/\tau_{} + \gamma (1 - s_{pctosst}) R_{pc}` ``

`` \frac{d s_{pctopv}}{dt} = -s_{pctopv}/\tau_{} + \gamma (1 - s_{pctopv}) R_{pc}` ``

`` \frac{d s_{pctovip}}{dt} = -s_{pctovip}/\tau_{} + \gamma (1 - s_{pctovip}) R_{pc}` ``

`` \frac{d s_{inttosst}}{dt} = -s_{inttosst}/\tau_{} + \gamma (1 - s_{inttosst}) R_{int}``

"""

# ╔═╡ 30bfa65f-3d93-408a-8a2f-078852e589d1


# ╔═╡ Cell order:
# ╠═ff667110-b363-11eb-0f66-37c2f5cf3240
# ╠═30a76def-2a04-412d-89c9-383a106ed7ae
# ╠═30bfa65f-3d93-408a-8a2f-078852e589d1
