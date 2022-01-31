# test creating hertag network


using DrWatson
@quickactivate "1-project-categorical-MMN"

# loading packages
using Parameters
using Revise

push!(LOAD_PATH,srcdir())


import MyNeurosciencePackage.NeuronalStructures: neuron, dendrite, pyr_cell, interneuron, synapses, eq_diff_method

import MyNeurosciencePackage.NeuronalStructures: neural_population, neural_motif

import MyNeurosciencePackage.NetworkConstruction: construct_network_hertag!

import MyNeurosciencePackage.NeuronalStructures: euler_method

import MyNeurosciencePackage.DynamicsFunction: simulation_step!
# Order of neural populaiton
# \         Pyr Cells, Dendrites PV, SST, VIP
# Pyr Cells
# Dendrites
# PV
# SST
# VIP

connection_matrix = [[0.0  0.0  0.6  0.0   0.0] 
[0.1  0.0  0.0  0.55  0.0] 
[0.45  0.0  0.5  0.6  0.5] 
[0.35  0.0  0.0  0.0  0.5] 
[0.1  0.0  0.0  0.45  0.0]]

connection_weights = [[0.0 0.0 -2.8 0.0  0.0] 
[0.42 0.0 0.0 -3.5 0.0] 
[1.5 0.0 -0.1 -0.3 -0.6] 
[1.0 0.0 0.0 0.0 -0.6] 
[1.0 0.0 0.0 -0.5 0.0]]


# Construction of the circuit motif
pyr_cell_pop = neural_population(70,
                                "soma hertag",
                                pyr_cell[],
                                "pyr_cells")

dendrites_pop = neural_population(70,
                                "dendrite hertag",
                                dendrite[],
                                "dendrites")

pv_pop = neural_population(10,
                                "rectified linear neurons",
                                neuron[],
                                "interneurons")

vip_pop = neural_population(10,
                                "rectified linear neurons",
                                neuron[],
                                "interneurons")

sst_pop = neural_population(10,
                                "rectified linear neurons",
                                neuron[],
                                "interneurons")



hertag_motif = neural_motif([pyr_cell_pop, dendrites_pop, pv_pop, sst_pop, vip_pop],
connection_matrix,
 connection_weights
)



construct_network_hertag!(hertag_motif)




for n in pyr_cell_pop.list_neurons
    push!(n.Ibg ,28) #Hz
end
for n in pv_pop.list_neurons
    push!(n.Ibg,2) #Hz
end
for n in vip_pop.list_neurons
    push!(n.Ibg,2) #Hz
end
for n in sst_pop.list_neurons
    push!(n.Ibg,2) #Hz
end

euler = euler_method(dt=0.0002)

for i=1:10000
simulation_step!(hertag_motif,euler)
end


@time simulation_step!(hertag_motif,euler)


### test de vitesse au propre
using Plots
plotlyjs()
plot(hertag_motif.list_pop[1].list_neurons[10].fr)


hertag_motif.list_pop[1].list_neurons[10].Ibg