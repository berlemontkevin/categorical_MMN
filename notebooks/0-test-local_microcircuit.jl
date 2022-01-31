using DrWatson
@quickactivate "1-project-categorical-MMN"

using Revise
#loading packages
using Parameters

push!(LOAD_PATH,srcdir())

using MyNeurosciencePackage



#######################
# Adding some functions



######################


## Create the network
microcircuit = create_network()

simu = simulation_parameters()

microcircuit.nn


input_d1 = [0.0,0.1,0.0,0.1,0.1,0.1,0.0]
input_d2 = [0.0,0.0,0.1,0.0,0.0,0.0,0.1]
duration_inputs = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
duration_break = 1.0


for i=1:length(input_d1)

    microcircuit.list_soma[1].Istim = input_d1[i]
    microcircuit.list_soma[2].Istim = input_d2[i]

   
    for t=1:Int(duration_inputs[i]/simu.dt)
        for nn in microcircuit.nn
            for n in nn.list_units
        push!(n.Istim, 0.3255)
       
        end
    end
        time_step(microcircuit,simu)

    end

    microcircuit.list_soma[1].Istim = 0.0
    microcircuit.list_soma[2].Istim = 0.0
    for t=1:Int(duration_break/simu.dt)
        for nn in microcircuit.nn
            for n in nn.list_units
        push!(n.Istim, 0.0)
       
        end
    end
        time_step(microcircuit,simu)

    end
end

#TODO integrer le courant de backgroun dans la transmission en g

#########################
# Plots (TODO en packge)

using Plots
plotlyjs()

Ttot = sum(duration_inputs) + duration_break*length(duration_inputs)
plot(0.0:simu.dt:Ttot,microcircuit.list_soma[1].r)
plot!(0.0:simu.dt:Ttot,microcircuit.list_soma[2].r)

plot(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[1].s)
plot!(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[2].s)

plot(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[1].r)
plot!(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[2].r)


plot(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[1].Ibg)
plot!(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[2].Ibg)

plot(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[1].Inoise)
plot!(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[2].Inoise)


plot(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[1].Istim)
plot!(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[2].Istim)

plot(0.0:simu.dt:Ttot,microcircuit.list_sst[1].r)
plot!(0.0:simu.dt:Ttot,microcircuit.list_sst[2].r)

plot(0.0:simu.dt:Ttot,microcircuit.list_sst[1].Iexc)
plot!(0.0:simu.dt:Ttot,microcircuit.list_sst[2].Iexc)


##################################

### Investigation influence du feedback top Down

microcircuit = create_network(0.42)
microcircuit_sans = create_network(0.0)

function temporal_simu(microcircuit)
    simu = simulation_parameters()



    input_d1 = [0.0,0.1,0.0,0.1,0.1,0.1,0.0]
    input_d2 = [0.0,0.0,0.1,0.0,0.0,0.0,0.1]
    duration_inputs = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
    duration_break = 1.0
    for i=1:length(input_d1)

    microcircuit.list_soma[1].Istim = input_d1[i]
    microcircuit.list_soma[2].Istim = input_d2[i]

   
    for t=1:Int(duration_inputs[i]/simu.dt)
        for nn in microcircuit.nn
            for n in nn.list_units
        push!(n.Istim, 0.3255)
       
        end
    end
        time_step(microcircuit,simu)

    end

    microcircuit.list_soma[1].Istim = 0.0
    microcircuit.list_soma[2].Istim = 0.0
    for t=1:Int(duration_break/simu.dt)
        for nn in microcircuit.nn
            for n in nn.list_units
        push!(n.Istim, 0.0)
       
        end
    end
        time_step(microcircuit,simu)

    end
    end
end

temporal_simu(microcircuit)

temporal_simu(microcircuit_sans)

plot(0.0:simu.dt:Ttot,microcircuit.nn[1].list_units[1].r,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.nn[1].list_units[1].r,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_soma[1].r,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_soma[1].r,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_vip[1].r,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_vip[1].r,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_soma[1].Iexc,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_soma[1].Iexc,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_soma[1].Iinh,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_soma[1].Iinh,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_dend[1].Iinh,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_dend[1].Iinh,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_dend[1].Iexc,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_dend[1].Iexc,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_dend[1].Ioutput,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_dend[1].Ioutput,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_soma[1].Itot,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_soma[1].Itot,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_pv[1].r,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_pv[1].r,linewidth=2)

plot(0.0:simu.dt:Ttot,microcircuit.list_sst[1].r,linewidth=2)
plot!(0.0:simu.dt:Ttot,microcircuit_sans.list_sst[1].r,linewidth=2)

###########################################