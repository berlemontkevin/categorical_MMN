### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 30bfa65f-3d93-408a-8a2f-078852e589d1
begin 
	using PlutoUI
param=Dict(
"τ_AMPA" => 0.002,
"τ_GABA" => 0.005,
"τ_int" => 0.8,
"τ_GABAdend" => 0.01,
"τ_NMDA" => 0.06,

"r0PV" => -95.0,
"r0SST" => -33.0,
"r0VIP" => -33.0,

"cPV" =>330.0,
"cSST" => 132.0,
"cVIP" => 132.0,

"c1" => 0.2,
"c2" =>  -7.0,
"c3" => -0.482,
"c4" => 0.00964,
"c5" => 0.11,
"c6" => 0.1,

"τadapt" => 1.0,

"Ibg_soma" =>  0.15,
"Ibg_dend" => 0.1,
"Ibg_pv" => 0.3,
"Ibg_sst" => 0.2,
"I_bg_vip" => 0.3,
"Ibg_int"=>0.310,

"a" => 135.0,
"b" => 54.0,
"d" => 0.308,

"γnmda" =>0.641*2,
"γgaba" => 2.0,


"g_ssttodend" => -0.09,
"g_pvtopc" => -0.005,
"g_pctopc" => 0.18,
"g_pctovip" => 0.058,
"g_ssttovip" => -0.1,
"g_pctosst" => 0.0435,
"g_viptosst" => -0.05,
"g_pctopv" => 0.08435,
"g_ssttopv" => -0.17,
"g_pvtopv" => -0.18,
"g_inttosst" => 0.22,
"g_pctoint" => 0.15,

"α" => -0.2,


"τA" => 0.1,
"gA" => -0.5,

"τ_facilitation" => 1.5,
"U_inttosst" => 0.2

)
end

# ╔═╡ 7b9f7bf0-dec7-4d03-8c70-4583153c27c8
md""" # Schematic of the model


$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\A1-Homemade\\model-with-integrator.png"))


"""

# ╔═╡ ff667110-b363-11eb-0f66-37c2f5cf3240
md""" # Theoretical study of adaptation and facilitation in a local circuit


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

# ╔═╡ 30a76def-2a04-412d-89c9-383a106ed7ae
md""" ## Adaptation and Facilitation

`` \frac{du}{dt} = \frac{U - u}{τ^ u} + U(1-u) r_E ``

`` \frac{ds}{dt} = - \frac{s}{τ} + u(1-s)γr_E ``


Only applies from integrator to SST for now.



`` \frac{dsA}{dt} = -sA/\tau + r(t) ``
with ``sA`` the adatation variable.

"""

# ╔═╡ 56ea37d7-8ba9-4af7-95d5-34f4f730f32c
md"""
All these equations can be solved theoretically in order to study more easily the influence of the different parameters on the solutions, when subject to a constant input.

Analysis of the influence of time constant for adaptation and facilitation let me decide to keep the one of Sean's paper for now, meaning ``\tau_{facilitation} = 1.5`` and ``\tau_A = 0.1``.


"""

# ╔═╡ 97701ee4-3254-47db-8d74-d3625e110bc4
md"""
$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\plots\\4-Theoretical-adaptation\\figΔ-firing-rate-stimlist.png"))

The x-axis represents the input then to the soma from 0 to 1 nA. The dots represent the difference of firing rate between between facilitaiton (green), or adaptation (red), and the microcircuit without it.

One can note that if for VIP and PV adaptation and facilitation have a similar impact, for the three other subtypes it's not the case.


Moreover, a strong stimulus leads to PC and the integrator to be similar to the standard case.

For SST, however, it always lead to an increase in firing rate, which is logical due to the facilitation mechanism.





"""

# ╔═╡ 2cfc1ce6-db13-4894-8d9c-3bc41ecf9019
md"""
The case of non constant current can't really be solved theoretically as too many variables lead to the steady state. However, it is possible to study this problem numerically
"""

# ╔═╡ c3a021a0-aecd-4a90-a69b-3db35b293139


# ╔═╡ Cell order:
# ╟─7b9f7bf0-dec7-4d03-8c70-4583153c27c8
# ╟─ff667110-b363-11eb-0f66-37c2f5cf3240
# ╟─30a76def-2a04-412d-89c9-383a106ed7ae
# ╟─56ea37d7-8ba9-4af7-95d5-34f4f730f32c
# ╟─30bfa65f-3d93-408a-8a2f-078852e589d1
# ╟─97701ee4-3254-47db-8d74-d3625e110bc4
# ╟─2cfc1ce6-db13-4894-8d9c-3bc41ecf9019
# ╠═c3a021a0-aecd-4a90-a69b-3db35b293139
