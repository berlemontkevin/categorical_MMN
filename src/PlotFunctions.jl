module PlotFunctions

using CairoMakie
using Parameters
using ..NeuronalStructures

export plot_local_circuit, plot_local_circuit_synapses



function plot_local_circuit(lc::Array{microcircuit},sim::simulation_parameters,list_current, params::Dict{String,Float64})
		# list_current being the list of currents that are sent to the networks

    @unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
            resolution = (2000, 1000), font = noto_sans)
    ax1 = fig[1, 1] = Axis(fig, title = "Submodule 1")

    time = 0.0:sim.dt:(sim.dt*(length(list_current[lc[1].list_soma[1].name])-1))
    
    max_yaxis = maximum(lc[1].list_soma[1].r_save)

    pc1 = lines!(ax1,time,lc[1].list_soma[1].r_save, color=:blue, linewidth=2,label="PC")
    vip1 = lines!(ax1,time,lc[1].list_vip[1].r_save, color=:green, linewidth=2,label="VIP")
    pv1 = lines!(ax1,time,lc[1].list_pv[1].r_save, color=:orange, linewidth=2,label="PV")
    sst1 = lines!(ax1,time,lc[1].list_sst[1].r_save, color=:red, linewidth=2,label="SST")
    int1 = lines!(ax1,time,lc[1].list_integrator[1].r_save, color=:purple,linewidth=2)
    
    
    
    ax2 = fig[1, 2] = Axis(fig, title = "Submodule 2")

    pc2 = lines!(ax2,time,lc[2].list_soma[1].r_save, color=:blue, linewidth=2,label="PC")
    vip2 = lines!(ax2,time,lc[2].list_vip[1].r_save, color=:green, linewidth=2,label="VIP")
    pv2 = lines!(ax2,time,lc[2].list_pv[1].r_save, color=:orange, linewidth=2,label="PV")
    sst2 = lines!(ax2,time,lc[2].list_sst[1].r_save, color=:red, linewidth=2,label="SST")
    int2 = lines!(ax2,time,lc[2].list_integrator[1].r_save, color=:purple,linewidth=2)

    
#     ax3 = fig[3,1] = Axis(fig, title = "Currents to Ecells")
    
#     cpc1 = lines!(ax3, list_current[lc[1].list_soma[1].name])
#     cpc2 = lines!(ax3, list_current[lc[2].list_soma[1].name])
    
    leg = fig[1, end+1] = Legend(fig,
                            [pc1,vip1,pv1,sst1,int1],
                            ["PC", "VIP", "PV", "SST","Integrator"])
    
    ylims!(ax1,0,max_yaxis + 10.0)
    ylims!(ax2,0,max_yaxis + 10.0)
    


    stim_time = [(initial_time + Tinter*i+Tstim*(i-1), 5.0 + max_yaxis) for i=1:(number_repetitions+1)]
  
    for i=1:number_repetitions
        text!(ax1, "A", position = (stim_time[round(Int,i)]), color=:red,  textsize = 30.0)
        text!(ax2, "A", position = (stim_time[round(Int,i)]), color=:red,  textsize = 30.0)
    end
    text!(ax1, "B", position = stim_time[end], color=:blue,  textsize = 30.0)
    text!(ax2, "B", position = stim_time[end], color=:blue,  textsize = 30.0)

    return fig
end



function plot_local_circuit_synapses(lc::Array{microcircuit},sim::simulation_parameters,list_current, params::Dict{String,Float64})
    # list_current being the list of currents that are sent to the networks

@unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
        resolution = (2000, 1000), font = noto_sans)
ax1 = fig[1, 1] = Axis(fig, title = "Submodule 1")

time = 0.0:sim.dt:(sim.dt*(length(list_current[lc[1].list_soma[1].name])-1))

 max_yaxis = maximum(lc[1].list_integrator[1].list_syn_pre_nmda[1].s_save)

syn_int_to_vip1 = lines!(ax1,time,lc[1].list_integrator[1].list_syn_pre_nmda[1].s_save, color=:blue, linewidth=2,label="Syn Int to VIP")
syn_int_to_sst1 = lines!(ax1,time,lc[1].list_integrator[1].list_syn_pre_nmda[3].s_save, color=:red, linewidth=2,label="Syn Int to VIP")
#syn_int_to_pv1 = lines!(ax1,time,lc[1].list_soma[1].Iexc_save, color=:blue, linewidth=2,label="Syn Int to VIP")

# vip1 = lines!(ax1,time,lc[1].list_vip[1].r_save, color=:green, linewidth=2,label="VIP")
# pv1 = lines!(ax1,time,lc[1].list_pv[1].r_save, color=:orange, linewidth=2,label="PV")
# sst1 = lines!(ax1,time,lc[1].list_sst[1].r_save, color=:red, linewidth=2,label="SST")
# int1 = lines!(ax1,time,lc[1].list_integrator[1].r_save, color=:purple,linewidth=2)



 ax2 = fig[1, 2] = Axis(fig, title = "Submodule 2")

 syn_int_to_vip2 = lines!(ax2,time,lc[2].list_integrator[1].list_syn_pre_nmda[1].s_save, color=:blue, linewidth=2,label="Syn Int to VIP")
 syn_int_to_sst2 = lines!(ax2,time,lc[2].list_integrator[1].list_syn_pre_nmda[3].s_save, color=:red, linewidth=2,label="Syn Int to VIP")
 # pv2 = lines!(ax2,time,lc[2].list_pv[1].r_save, color=:orange, linewidth=2,label="PV")
# sst2 = lines!(ax2,time,lc[2].list_sst[1].r_save, color=:red, linewidth=2,label="SST")
# int2 = lines!(ax2,time,lc[2].list_integrator[1].r_save, color=:purple,linewidth=2)


#     ax3 = fig[3,1] = Axis(fig, title = "Currents to Ecells")

#     cpc1 = lines!(ax3, list_current[lc[1].list_soma[1].name])
#     cpc2 = lines!(ax3, list_current[lc[2].list_soma[1].name])

# leg = fig[1, end+1] = Legend(fig,
#                         [pc1,vip1,pv1,sst1,int1],
#                         ["PC", "VIP", "PV", "SST","Integrator"])

# ylims!(ax1,0,max_yaxis + 10.0)
# ylims!(ax2,0,max_yaxis + 10.0)



# stim_time = [(initial_time + Tinter*i+Tstim*(i-1), 5.0 + max_yaxis) for i=1:(number_repetitions+1)]

# for i=1:number_repetitions
#     text!(ax1, "A", position = (stim_time[round(Int,i)]), color=:red,  textsize = 30.0)
#     text!(ax2, "A", position = (stim_time[round(Int,i)]), color=:red,  textsize = 30.0)
# end
# text!(ax1, "B", position = stim_time[end], color=:blue,  textsize = 30.0)
# text!(ax2, "B", position = stim_time[end], color=:blue,  textsize = 30.0)

return fig
end









end