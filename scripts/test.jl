using DrWatson
@quickactivate "1-project-categorical-MMN"

push!(LOAD_PATH,srcdir())

using MyNeurosciencePackage



using Plots
plotlyjs()



include(srcdir("structures.jl"))


pyr_cell = pyr_cells_hertag2020(rE = 0.0,
I = 0.0,
ID = 0.0,
IE= 28.0,
I0D =0.0,
c=0.0)




function simulation(dt,t_final)
    # perform the temporal simulation of the dendrite
    for i=1:t_final
    update_fr_soma!(pyr_cell,dt)
    update_current!(pyr_cell)
    compute_calcium_spike!(pyr_cell)
    update_I0D!(pyr_cell)
end

end

dt = 0.5
t_final = 100

simulation(dt,t_final)
pyr_cell

t_final = 500

save_r = []
for i=1:200
I = 1.0*i
pyr_cell = pyr_cells_hertag2020(rE = 0.0,
I = 0.0,
ID = 0.0,
IE= I,
I0D =0.0,
c=0.0)
simulation(dt,t_final)
push!(save_r,pyr_cell.rE)
end

plot(save_r)