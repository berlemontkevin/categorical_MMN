using DrWatson
@quickactivate "1-project-categorical-MMN"


using Plots
plotlyjs()



include(srcdir("structures.jl"))



Gabeprop = syn_GABA_dendrites()

d = dendrites(EL=-70.0)

pcell = pyr_cells(dendrites_list=[d])


ilist = LinRange(-200,50,50)

soma_prop = firing_rate_soma()

FRlist = []
for i in ilist
    pcell.I =i
    update_fr_soma!(pcell,soma_prop)
    push!(FRlist, pcell.r)
end


plot(ilist,FRlist)


d.rI = 0.0 #Hz

gelist = LinRange(0.0,50.0,40)
nmda = syn_NMDA_dendrites()

VDlist = []
for i in gelist
    d.gE =i
    update_gI!(d,Gabeprop)
    update_β!(d)
    update_snmda(d,nmda)
    #update_gE!(d)
    update_g12!(d)
    update_VD!(d)
    push!(VDlist,d.VD)
end

plot(gelist, VDlist)
