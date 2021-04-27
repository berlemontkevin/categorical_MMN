using DrWatson
@quickactivate "1-project-categorical-MMN"

using Revise
#loading packages
using Parameters

push!(LOAD_PATH,srcdir())

using Plots
plotlyjs()
using MyNeurosciencePackage

######################
# Larkum 2019 model
###########################

const g = 1.2*sqrt(1000)
const Ithr = 0.5 # nA
const D = 0.5

FI(Isoma,Idend) = g * sqrt(rect_linear(Isoma + D*Idend - Ithr))


Isoma = 0.0:0.01:1.0
Idend = 0.0:0.01:1.0

FIrate = zeros(length(Isoma),length(Idend))
for i=1:length(Isoma)
    for j=1:length(Idend)

        FIrate[i,j] = FI(Isoma[i],Idend[j])
    end
end

plot(Isoma,Idend,FIrate,st=:surface)
plot!(xlabel="Isoma")
plot!(ylabel="Idend")


fig = plot()
for i=1:20:100
    plot!(fig,Idend,FIrate[i,:],label=Isoma[i],linewidth=2)
end
fig

figIsoma = plot()
for i=1:20:100
    plot!(figIsoma,Isoma,FIrate[:,i],label=Idend[i],linewidth=2)
end
figIsoma


###############
# Seans model
#################


FIsean(Iexc,Iinh) = 0.12*tanh((Iexc + 7.0*Iinh + 0.0)/(0.00964*exp(-Iinh/0.02))) + 0.13624

AbottChance(x) = (135*x -54)/(1.0-exp(-0.308*(135*x-54)))

Iexc=0.03
plot(Iinh,FIsean.(Iexc,Iinh),linewidth=2,label = Iexc)



Iexc = 0.03
Iinh = -0.5:0.001:0.0
plot(Iinh,AbottChance.(0.31.+FIsean.(Iexc,Iinh)),linewidth=2,label = Iexc)

Iexc = 0.05
plot!(Iinh,AbottChance.(0.31.+FIsean.(Iexc,Iinh)),linewidth=2,label=Iexc)

Iexc = 0.1
plot!(Iinh,AbottChance.(0.31.+FIsean.(Iexc,Iinh)),linewidth=2,label=Iexc)

Iexc = 0.0
plot!(Iinh,AbottChance.(0.31.+FIsean.(Iexc,Iinh)),linewidth=2,label=Iexc)


################### 
# Spiking calcium event
#####################