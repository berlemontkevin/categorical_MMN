module PlotFunctions

using CairoMakie
using ..NeuronalStructures

export plot_local_circuit



function plot_local_circuit(lc::Array{microcircuit},sim::simulation_parameters,list_current)
		# list_current being the list of currents that are sent to the networks
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
            resolution = (1000, 1400), font = noto_sans)
    ax1 = fig[1, 1] = Axis(fig, title = "Local circuit 1")

    time = 0.0:sim.dt:(sim.dt*length(list_current[lc[1].list_soma[1].name]))
    
    
    pc1 = lines!(ax1,time,lc[1].list_soma[1].r_save, color=:blue, linewidth=2,label="PC")
    vip1 = lines!(ax1,time,lc[1].list_vip[1].r_save, color=:green, linewidth=2,label="VIP")
    pv1 = lines!(ax1,time,lc[1].list_pv[1].r_save, color=:orange, linewidth=2,label="PV")
    sst1 = lines!(ax1,time,lc[1].list_sst[1].r_save, color=:red, linewidth=2,label="SST")
    int1 = lines!(ax1,time,lc[1].list_integrator[1].r_save, color=:purple,linewidth=2)
    
    
    
    ax2 = fig[2, 1] = Axis(fig, title = "Local circuit 2")

    pc2 = lines!(ax2,time,lc[2].list_soma[1].r_save, color=:blue, linewidth=2,label="PC")
    vip2 = lines!(ax2,time,lc[2].list_vip[1].r_save, color=:green, linewidth=2,label="VIP")
    pv2 = lines!(ax2,time,lc[2].list_pv[1].r_save, color=:orange, linewidth=2,label="PV")
    sst2 = lines!(ax2,time,lc[2].list_sst[1].r_save, color=:red, linewidth=2,label="SST")
    int2 = lines!(ax2,time,lc[2].list_integrator[1].r_save, color=:purple,linewidth=2)

    
    ax3 = fig[3,1] = Axis(fig, title = "Currents to Ecells")
    
    cpc1 = lines!(ax3, list_current[lc[1].list_soma[1].name])
    cpc2 = lines!(ax3, list_current[lc[2].list_soma[1].name])
    
    leg = fig[1, end+1] = Legend(fig,
                            [pc1,vip1,pv1,sst1,int1],
                            ["PC", "VIP", "PV", "SST","Integrator"])
    
    ylims!(ax1,0,40.0)
    ylims!(ax2,0,40.0)
    
    return fig
end









end