using DrWatson
@quickactivate "1-project-categorical-MMN"

using Revise
#loading packages
using Parameters

push!(LOAD_PATH,srcdir())

using MyNeurosciencePackage




#############
# Test construction modele 2 areas


c = MyNeurosciencePackage.NetworkConstruction.local_microcircuit_network.create_network()


simu = MyNeurosciencePackage.NeuronalStructures.NeuralNetwork.simulation_parameters()

c.nn

for i=1:1000
    #MyNeurosciencePackage.DynamicsFunction.time_dynamics.time_step(c,simu)
    #    push!(temp,E1.Itot)
    time_step(c,simu)
end


c.list_dend[1].Istim = 0.1
for i=1:1000
    MyNeurosciencePackage.DynamicsFunction.time_dynamics.time_step(c,simu)
#    push!(temp,E1.Itot)
end

using Plots
gr()
plotlyjs()

plot(c.list_soma[1].r)


plot(c.nn[1].list_units[1].r)

#



### Definition of Hertag parameters
Hertag_struct = syn_probabilities_hertag(pEP = 0.6,
pDS = 0.55,
pDE = 0.1,
pPE = 0.45,
pPP = 0.5,
pPS = 0.6,
pPV = 0.5,
pSE = 0.35,
pSV = 0.5,
pVE = 0.1,
pVS=0.45)

function eqn8_9_hertag!(VP,MP,VE,s::syn_weights_hertag)
    # compute the synaptic weights following Hertag paper
    s.wPS = -(VP + abs(s.wVS)*MP - (1-s.wPP)/(abs(s.wEP)) * VE) # gain = 0.07
    s.wPV = -(abs(s.wSV)*abs(s.wPS) + (1-s.wSV*s.wVS)*MP)
    return s.wPS, s.wPV
end


Hertag_weights = syn_weights_hertag(wEP = -2.8,
wDS = -3.5,
wDE = 0.42,
wPE = 1.5,
wPP = -0.1,
wPS = -0.3,
wPV = -0.6,
wSE = 1.0,
wSV = -0.6,
wVE = 1.0,
wVS=-0.5)

eqn8_9_hertag!(1.0,0.0,1.0,Hertag_weights)

Hertag_weights.wPS

Hertag_weights.wPV

normalisation_weights(Hertag_weights,Hertag_struct)
test = construct_network_hertag(Hertag_struct,Hertag_weights)

for n in test.pc_list
    n.Ibg = 28 #Hz
end
for n in test.pv_list
    n.Ibg = 2 #Hz
end
for n in test.vip_list
    n.Ibg = 2 #Hz
end
for n in test.sst_list
    n.Ibg = 2 #Hz
end


function backgroun_input_objective!(netw::network_hertag)

    for n in netw.pc_list
        background_objective!(n)
    end

    for n in netw.pv_list
        background_objective!(n)
    end

    for n in netw.vip_list
        background_objective!(n)
    end

    for n in netw.sst_list
        background_objective!(n)
    end



end


euler = euler_method(dt=0.0002)
backgroun_input_objective!(test)
temp_pyr = []
temp_sst = []
temp_vip = []
temp_pv = []
temp_Itot = []


for i=1:1000
simulation_step!(test,euler)
push!(temp_pyr,test.pc_list[1].fr)
push!(temp_sst,test.sst_list[1].fr)
push!(temp_vip,test.vip_list[5].fr)
push!(temp_pv,test.pv_list[1].fr)
push!(temp_Itot, test.pc_list[1].Itot)


end


nlist = []
for d in test.dendrites_list
    push!(nlist,length(d.neurons_list))
end
histogram(nlist)

using Plots
plotlyjs()


plot(temp_Itot)

plot(temp_pyr)
plot!(temp_sst)
plot(temp_pv)
plot(temp_vip)

temp_vip
temp_pv
temp_sst





test.pv_list[4]

length(test.pv_list[1].neurons_list)

