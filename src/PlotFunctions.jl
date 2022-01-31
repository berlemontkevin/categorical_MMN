module PlotFunctions

using CairoMakie
using Parameters
using ..NeuronalStructures

export plot_local_circuit, plot_local_circuit_synapses



function plot_local_circuit(lc::Vector{microcircuit{soma_PC,dend_sigmoid}},sim::simulation_parameters,current, params::Dict{String,Float64})
		# current being the list of currents that are sent to the networks

    @unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
            resolution = (2000, 1000), font = noto_sans)
    ax1 = fig[1, 1] = Axis(fig, title = "Submodule 1",xlabel="Time (s)", ylabel="Firing rate (Hz)")

    time = 0.0:lc[1].eq_diff_method.dt:(lc[1].eq_diff_method.dt*(length(current[lc[1].soma.name])-1))
    
    max_yaxis = maximum(lc[1].soma.r_save)

    pc1 = lines!(ax1,time,lc[1].soma.r_save, color=:blue, linewidth=2,label="PC")
    vip1 = lines!(ax1,time,lc[1].vip.r_save, color=:green, linewidth=2,label="VIP")
    pv1 = lines!(ax1,time,lc[1].pv.r_save, color=:orange, linewidth=2,label="PV")
    sst1 = lines!(ax1,time,lc[1].sst.r_save, color=:red, linewidth=2,label="SST")
    int1 = lines!(ax1,time,lc[1].integrator.r_save, color=:purple,linewidth=4)
    
    
    
    ax2 = fig[1, 2] = Axis(fig, title = "Submodule 2",xlabel="Time (s)", ylabel="Firing rate (Hz)")

    pc2 = lines!(ax2,time,lc[2].soma.r_save, color=:blue, linewidth=2,label="PC")
    vip2 = lines!(ax2,time,lc[2].vip.r_save, color=:green, linewidth=2,label="VIP")
    pv2 = lines!(ax2,time,lc[2].pv.r_save, color=:orange, linewidth=2,label="PV")
    sst2 = lines!(ax2,time,lc[2].sst.r_save, color=:red, linewidth=2,label="SST")
    int2 = lines!(ax2,time,lc[2].integrator.r_save, color=:purple,linewidth=4)

    
    leg = fig[1, end+1] = Legend(fig,
                            [pc1,vip1,pv1,sst1,int1],
                            ["PC", "VIP", "PV", "SST","Integrator"],labelsize = 30)
    
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


function plot_local_circuit(lc::Array{microcircuit},sim::simulation_parameters,current)
    # current being the list of currents that are sent to the networks

    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
            resolution = (2000, 1000), font = noto_sans)
    ax1 = fig[1, 1] = Axis(fig, title = "Submodule 1",xlabel="Time (s)", ylabel="Firing rate (Hz)")

    time = 0.0:sim.dt:(sim.dt*(length(current[lc[1].soma[1].name])-1))

    max_yaxis = maximum(lc[1].soma[1].r_save)

    pc1 = lines!(ax1,time,lc[1].soma[1].r_save, color=:blue, linewidth=2,label="PC")
    vip1 = lines!(ax1,time,lc[1].vip[1].r_save, color=:green, linewidth=2,label="VIP")
    pv1 = lines!(ax1,time,lc[1].pv[1].r_save, color=:orange, linewidth=2,label="PV")
    sst1 = lines!(ax1,time,lc[1].sst[1].r_save, color=:red, linewidth=2,label="SST")
    int1 = lines!(ax1,time,lc[1].integrator[1].r_save, color=:purple,linewidth=4)



    ax2 = fig[1, 2] = Axis(fig, title = "Submodule 2",xlabel="Time (s)", ylabel="Firing rate (Hz)")

    pc2 = lines!(ax2,time,lc[2].soma[1].r_save, color=:blue, linewidth=2,label="PC")
    vip2 = lines!(ax2,time,lc[2].vip[1].r_save, color=:green, linewidth=2,label="VIP")
    pv2 = lines!(ax2,time,lc[2].pv[1].r_save, color=:orange, linewidth=2,label="PV")
    sst2 = lines!(ax2,time,lc[2].sst[1].r_save, color=:red, linewidth=2,label="SST")
    int2 = lines!(ax2,time,lc[2].integrator[1].r_save, color=:purple,linewidth=4)


    #     ax3 = fig[3,1] = Axis(fig, title = "Currents to Ecells")

    #     cpc1 = lines!(ax3, current[lc[1].soma[1].name])
    #     cpc2 = lines!(ax3, current[lc[2].soma[1].name])

    leg = fig[1, end+1] = Legend(fig,
                            [pc1,vip1,pv1,sst1,int1],
                            ["PC", "VIP", "PV", "SST","Integrator"],labelsize = 30)

    ylims!(ax1,0,max_yaxis + 10.0)
    ylims!(ax2,0,max_yaxis + 10.0)





return fig
end


function plot_local_circuit_synapses(lc::Array{microcircuit},sim::simulation_parameters,current, params::Dict{String,Float64})
        # current being the list of currents that are sent to the networks

    @unpack number_repetitions, value_stim_f1, Tinter, Tstim, Tfin, value_stim_f2, τ_integrator, initial_time = params

    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
            resolution = (2000, 1000), font = noto_sans)
    ax1 = fig[1, 1] = Axis(fig, title = "Submodule 1")

    time = 0.0:sim.dt:(sim.dt*(length(current[lc[1].soma[1].name])-1))

    max_yaxis = maximum(lc[1].integrator[1].syn_pre_nmda[1].s_save)

    syn_int_to_vip1 = lines!(ax1,time,lc[1].integrator[1].syn_pre_nmda[1].s_save, color=:blue, linewidth=2,label="Syn Int to VIP")
    syn_int_to_sst1 = lines!(ax1,time,lc[1].integrator[1].syn_pre_nmda[3].s_save, color=:red, linewidth=2,label="Syn Int to VIP")
    #syn_int_to_pv1 = lines!(ax1,time,lc[1].soma[1].Iexc_save, color=:blue, linewidth=2,label="Syn Int to VIP")

    # vip1 = lines!(ax1,time,lc[1].vip[1].r_save, color=:green, linewidth=2,label="VIP")
    # pv1 = lines!(ax1,time,lc[1].pv[1].r_save, color=:orange, linewidth=2,label="PV")
    # sst1 = lines!(ax1,time,lc[1].sst[1].r_save, color=:red, linewidth=2,label="SST")
    # int1 = lines!(ax1,time,lc[1].integrator[1].r_save, color=:purple,linewidth=2)



    ax2 = fig[1, 2] = Axis(fig, title = "Submodule 2")

    syn_int_to_vip2 = lines!(ax2,time,lc[2].integrator[1].syn_pre_nmda[1].s_save, color=:blue, linewidth=2,label="Syn Int to VIP")
    syn_int_to_sst2 = lines!(ax2,time,lc[2].integrator[1].syn_pre_nmda[3].s_save, color=:red, linewidth=2,label="Syn Int to VIP")
    # pv2 = lines!(ax2,time,lc[2].pv[1].r_save, color=:orange, linewidth=2,label="PV")
    # sst2 = lines!(ax2,time,lc[2].sst[1].r_save, color=:red, linewidth=2,label="SST")
    # int2 = lines!(ax2,time,lc[2].integrator[1].r_save, color=:purple,linewidth=2)


    #     ax3 = fig[3,1] = Axis(fig, title = "Currents to Ecells")

    #     cpc1 = lines!(ax3, current[lc[1].soma[1].name])
    #     cpc2 = lines!(ax3, current[lc[2].soma[1].name])

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




module bump_plots
using CairoMakie
using ..NeuronalStructures


"""
    bump_animation(layer_bump,time_tot, framerate, save_path, save_name, save_format)
"""
function bump_animation(layer_bump, frame_step, save_path, save_name, save_format; start_time = 1.0)

    time = Node(start_time)

    x_range = range(1,length(layer_bump.list_microcircuit), length = length(layer_bump.list_microcircuit))

    y_pc = @lift(get_firing_rate(layer_bump, $time, "soma"))
    y_vip = @lift(get_firing_rate(layer_bump, $time, "vip"))
    y_pv = @lift(get_firing_rate(layer_bump, $time, "pv"))
    y_sst = @lift(get_firing_rate(layer_bump, $time, "sst"))
    y_int = @lift(get_firing_rate(layer_bump, $time, "integrator"))
    

    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
            resolution = (2000, 1000))
    ax1 = fig[1, 1] = Axis(fig, title = @lift("time = $(round(0.5.*$time, digits = 1)) ms"))

    ylims!(ax1,0,60.0)
    lines!(ax1, x_range, y_pc, color=:blue, linewidth=2,label="PC")
    lines!(ax1, x_range, y_vip, color=:green, linewidth=2,label="VIP")
    lines!(ax1, x_range, y_pv, color=:orange, linewidth=2,label="PV")
    lines!(ax1, x_range, y_sst, color=:red, linewidth=2,label="SST")
    lines!(ax1, x_range, y_int, color=:purple, linewidth=2,label="Integrator")

    time_tot = length(layer_bump.list_microcircuit[1].soma.r_save)

    timestamps = range(1, time_tot, step=floor(Int, frame_step))

    record(fig, string(save_path,save_name, save_format), timestamps; framerate = 50) do t
        time[] = floor(Int,t)
    end
    return nothing
end


"""
    bump_animation(layer_bump,time_tot, framerate, save_path, save_name, save_format)
"""
function bump_animation(nn::neural_network{soma_PC, dend_sigmoid, euler_method}, frame_step, save_path, save_name, save_format; start_time = 1.0)

    time = Node(start_time)

    x_range = range(1,length(nn.list_microcircuit.neurons), length = length(nn.list_microcircuit.neurons))

    y_pc = @lift(get_firing_rate(nn, $time, "soma"))
    y_vip = @lift(get_firing_rate(nn, $time, "vip"))
    y_pv = @lift(get_firing_rate(nn, $time, "pv"))
    y_sst = @lift(get_firing_rate(nn, $time, "sst"))
    y_ngfc = @lift(get_firing_rate(nn, $time, "ngfc"))
    y_int = @lift(get_firing_rate(nn.ring_integrator, $time))
    

    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
            resolution = (2000, 1000))
    ax1 = fig[1, 1] = Axis(fig, title = @lift("time = $(round(0.5.*$time, digits = 1)) ms"))

    ylims!(ax1,0,50.0)
    lines!(ax1, x_range, y_pc, color=:blue, linewidth=4,label="PC")
    lines!(ax1, x_range, y_vip, color=:green, linewidth=4,label="VIP")
    lines!(ax1, x_range, y_pv, color=:orange, linewidth=4,label="PV")
    lines!(ax1, x_range, y_sst, color=:red, linewidth=4,label="SST")
    lines!(ax1, x_range, y_ngfc, color=:black, linewidth=4,label="NGFC")
    lines!(ax1, 1:2:128, y_int, color=:purple, linewidth=4,label="Integrator")

    time_tot = length(nn.list_microcircuit.neurons[1].soma.r_save)

    timestamps = range(1, time_tot, step=floor(Int, frame_step))

    fig[1, 2] = Legend(fig, ax1, "Legend", framevisible = false)

    record(fig, string(save_path,save_name, save_format), timestamps; framerate = 50) do t
        time[] = floor(Int,t)
    end
    return nothing
end
export bump_animation

end

using .bump_plots

export bump_animation




end