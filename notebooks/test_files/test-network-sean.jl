using DrWatson
@quickactivate "1-project-categorical-MMN"

#loading packages
using Parameters

push!(LOAD_PATH,srcdir())


"""
    Type: eq_diff_method

Defines the type of method that will be used to solve the differential equations

...
# Types currently implemented
- Euler method
....
"""
abstract type eq_diff_method end

"""
    Type: neuron

Abstract type for all the type of neurons

...
# Neurons currently implemented
- dendrites
- pyramidal cells
- interneurons
....
"""
abstract type neuron end



"""
    Type: synapses
    
Abstract type for the types of synapses
"""
abstract type synapses end


"""
    Type: dendrite
    
Abstract type for the types of dendrites
"""
abstract type dendrite <: neuron end

"""
    Type: pyr_cell
    
Abstract type for the types of pyramidal cells
"""
abstract type pyr_cell <: neuron end

"""
    Type: interneuron
    
Abstract type for the types of interneuron
"""
abstract type interneuron <: neuron end

@with_kw mutable struct OU_process
    τ::Float64 = 0.002
dt::Float64 = 0.0005
σ::Float64 = 5 * 0.001
end


@with_kw struct dendrites_param
    c1::Float64 = 0.12 #* 0.001
    c2::Float64 = 0.13624 #* 0.001
    c3::Float64 = 7.0
    c4::Float64 = 0.0 #* 0.001
    c5::Float64 = 0.00964 #* 0.001
    c6::Float64 = 0.02 #* 0.001 #nA units

end


@with_kw mutable struct soma <: pyr_cell

    param_c = dendrites_param()
    Idendexc::Float64 = 0.0
    Idendinh::Float64 = 0.0
    f_I_curve = zeros(3)
    r::Array{Float64} = [0.0]
    Iinput::Array{Float64} = [0.0]
    a::Float64 = 0.135 * 1000
    b::Float64 = 54
    c::Float64 = 0.308 #secondes
    den::dendrite #dendrite connected to the soma
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 310.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()

end


@with_kw mutable struct dend <:dendrite
    param_c = dendrites_param()
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    Ioutput::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 30.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()

end


@with_kw mutable struct pv_cell <: interneuron
    r::Array{Float64} = [0.0] #rate of the neuron
    c_I::Float64 = 330.0
    r0::Float64 = -95.0
    Iinput::Array{Float64} = [0.0]
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 300.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
end


@with_kw mutable struct sst_cell <: interneuron
  r::Array{Float64} = [0.0] #rate of the neuron
    c_I::Float64 = 132.0
    r0::Float64 = -33.0
    Iinput::Array{Float64} = [0.0]
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 300.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
end


@with_kw mutable struct vip_cell <: interneuron
    r::Array{Float64} = [0.0] #rate of the neuron
    c_I::Float64 = 132.0
    r0::Float64 = -33.0 #Hz
    Iinput::Array{Float64} = [0.0]
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Float64 = 300.0 * 0.001
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Float64 = 0.0
    dt::Float64=0.0005
    τ::Float64 = 0.002
    OU_process::OU_process = OU_process()
end

#let's do for now without a;; the dopamine gradient


@with_kw struct simulation_parameters
    dt::Float64 = 0.0005
    Tfin::Float64 = 2.0
    Tstimdebut::Float64 = 1.0
    
end


@with_kw struct wong_wang_cell <: neuron
    r::Array{Float64} = [0.0]
    s::Array{Float64} = [0.0] # synaptic gatign variable
    Iexc::Array{Float64} = [0.0]
    Iinh::Array{Float64} = [0.0]
    Ioutput::Array{Float64} = [0.0]
    list_syn::Array{synapses} = []
    Ibg::Array{Float64} = [0.3255]
    Itot::Array{Float64} = [0.0]
    Inoise::Array{Float64} = [0.0]
    Istim::Array{Float64} = [0.0]
    γ::Float64 = 0.641
    a::Float64 = 270.0
    b::Float64 = 108.0
    d::Float64 = 0.154

end

@with_kw struct wong_wang_network 
    noiseamp::Float64  =0.02# noise amp of the OU process
  i0::Float64 = 0.3255 #resting state of the OU process
  jn11::Float64 =0.2609# Synaptic strength unit 1 to 1
  jn21::Float64 = 0.0497# Synaptic strength unit 2 to 1
  jn12::Float64 = 0.0497# Synaptic strength unit 1 to 2
  jn22::Float64 = 0.2609# Synaptic strength unit 2 to 2
  tnmda::Float64 =100.0* 0.001# time constant of NMDA receptors
  tampa::Float64 = 2.0 * 0.001# time constant of AMPA receptors
  threshold::Float64# threhsold of the decision
  list_units::Array{wong_wang_cell} = []

end


function dendrite_input_output!(d::dend)
    @unpack  c1,c2,c3,c4,c5,c6 =  d.param_c 
    β = c5*exp(-d.Iinh[end]/c6)

    push!(d.Ioutput , c1*(tanh((d.Iexc[end] +c3*d.Iinh[end] + c4)/β)) + c2)
end

function dendrite_input_output!(a,b)
    @unpack  c1,c2,c3,c4,c5,c6 =  dend1.param_c 
    β = c5*exp.(-b./c6)

    return  c1.*(tanh.((a +c3.*b .+ c4)./β)) .+ c2
end


function sum_input(n::interneuron)
    push!(n.Itot , n.Iinh[end] + n.Iexc[end] + n.Inoise[end] + n.Ibg + n.Istim )
end

function sum_input(n::soma)
    push!(n.Itot , n.Iinh[end] + n.Iexc[end] + n.Inoise[end] + n.Ibg + n.Istim + n.den.Ioutput[end])

end


function current_to_frequency!(neuron::pyr_cell)
    # r is the rate of the neuron
    a = neuron.a
    b = neuron.b
    d = neuron.c
    
    input_current = neuron.Itot[end]

    return (a*input_current - b)/(1 - exp(-d*(a*input_current - b)))

end


function current_to_frequency!(neuron::interneuron)
    # r is the rate of the neuron
    c_I = neuron.c_I
    r0 = neuron.r0
    input_current = neuron.Itot[end]

    return rect_linear(c_I*input_current + r0)
end

function update_firing!(n::neuron)
    # update firing rates following Euler law
    push!(n.r , n.r[end] + n.dt/n.τ*(-n.r[end] + current_to_frequency!(n)))

end


function synapse_derivative(s::ampa_syn)

    return -s.s[end]/s.τ + s.γ*s.neuron_pre.r[end]

end


function synapse_derivative(s::gaba_syn)
    #paper Sean: Inhibition onto the dendrites is slower than inhibition elsewhere (ali and Thomson 2008)
    return -s.s[end]/s.τ + s.γ*s.neuron_pre.r[end]

end


function synapse_derivative(s::nmda_syn)

    return -s.s[end]/s.τ + s.γ*(1.0 - s.s[end])*s.neuron_pre.r[end]

end

function current_synapses!(n::neuron)
    #compute the sum of syn currents
    # separe in two the currents (due to the dendrites)

    push!(n.Iexc , 0.0)
    push!(n.Iinh , 0.0)

    for s in n.list_syn
        if s.g <0.0
            n.Iinh[end] += s.g * s.s[end]
        else
            n.Iexc[end] += s.g * s.s[end]
        end
    end
end
function OU_process(n::neuron)
    s = n.OU_process
    push!(n.Inoise , n.Inoise[end] + s.dt / s.τ * (-n.Inoise[end] + sqrt(s.τ*s.σ*s.σ)*randn()))
   
end

d_param = dendrites_param()


vip1 = vip_cell()
sst1 = sst_cell()
pv1 = pv_cell()
dend1 = dend()
E1 = soma(den=dend1)

vip2 = vip_cell()
sst2 = sst_cell()
dend2 = dend()
E2 = soma(den=dend2)



push!(dend1.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst1,neuron_post=dend1, g = -0.09))
push!(E1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E1, g = -0.001))
push!(E1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=E1, g = 0.18))

push!(vip1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=vip1, g = 0.058))
push!(vip1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=vip1, g = -0.1))

push!(sst1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst1, g = 0.0435))
push!(sst1.list_syn, gaba_syn(neuron_pre=vip1,neuron_post=sst1, g = -0.05))

push!(pv1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=pv1, g = 0.0435))
push!(pv1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=pv1, g = -0.17))
push!(pv1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=pv1, g = -0.18))

#connect to E2
push!(dend2.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst2,neuron_post=dend2, g = -0.09))
push!(E2.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E2, g = -0.001))
push!(E2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=E2, g = 0.18))

push!(vip2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=vip2, g = 0.058))
push!(vip2.list_syn, gaba_syn(neuron_pre=sst2,neuron_post=vip2, g = -0.1))

push!(sst2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=sst2, g = 0.0435))
push!(sst2.list_syn, gaba_syn(neuron_pre=vip2,neuron_post=sst2, g = -0.05))

push!(pv1.list_syn, nmda_syn(neuron_pre=E2,neuron_post=pv1, g = 0.0435))
push!(pv1.list_syn, gaba_syn(neuron_pre=sst2,neuron_post=pv1, g = -0.17))


push!(sst1.list_syn, nmda_syn(neuron_pre=E2,neuron_post=sst1, g = 0.0435))
push!(sst2.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst2, g = 0.0435))



function update_dend!(d::dend)
    current_synapses!(d)
    OU_process(d)
    d.Iexc[end] += d.Inoise[end] + d.Istim + d.Ibg
    dendrite_input_output!(d)
end

p = OU_process()

function time_step(c::microcircuit,sim::simulation_parameters)

    #update dend
    update_dend!.(c.list_dend)

    #update neurons
    current_synapses!.(c.list_soma)
    current_synapses!.(c.list_vip)
    current_synapses!.(c.list_sst)
    current_synapses!.(c.list_pv)
    OU_process.(c.list_soma)
    OU_process.(c.list_sst)
    OU_process.(c.list_vip)
    OU_process.(c.list_pv)
    current_to_frequency!.(c.list_soma)
    current_to_frequency!.(c.list_sst)
    current_to_frequency!.(c.list_vip)
    current_to_frequency!.(c.list_pv)


    sum_input.(c.list_soma)
    sum_input.(c.list_sst)
    sum_input.(c.list_vip)
    sum_input.(c.list_pv)


    update_firing!.(c.list_soma)
    update_firing!.(c.list_vip)
    update_firing!.(c.list_sst)
    update_firing!.(c.list_pv)

    for nn in c.nn
        for n in nn.list_units
            current_synapses!(n)
        
        end
    update_s!(nn,simu)
    end
for pop in [c.list_pv, c.list_dend, c.list_sst, c.list_vip, c.list_soma]
    for n in pop
        for s in n.list_syn
        push!(s.s , s.s[end]+ s.dt * synapse_derivative(s))
        end
    end

end

for nn in c.nn
    for n in nn.list_units
        for s in n.list_syn

            push!(s.s , s.s[end]+ s.dt * synapse_derivative(s))
            
        end
    end

end
end


function rect_linear(x::Float64)
    temp = max(0, x)
    return temp
end


@with_kw mutable struct microcircuit
    list_soma::Array{soma}= []
    list_sst::Array{sst_cell} = []
    list_vip::Array{vip_cell} = []
    list_pv::Array{pv_cell} = []
    list_dend::Array{dend} = []
    nn::Array{wong_wang_network} = []

end

circuit1 = microcircuit()

# TODO add parameters structure
function construct_modele_sean(c::microcircuit)
    # construct the microcircuit in sean case

    vip1 = vip_cell()
    sst1 = sst_cell()
    pv1 = pv_cell()
    dend1 = dend()
    E1 = soma(den=dend1)
    
    vip2 = vip_cell()
    sst2 = sst_cell()
    dend2 = dend()
    E2 = soma(den=dend2)
    
    
    
    push!(dend1.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst1,neuron_post=dend1, g = -0.09))
    push!(E1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E1, g = -0.001))
    push!(E1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=E1, g = 0.18))
    
    push!(vip1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=vip1, g = 0.058))
    push!(vip1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=vip1, g = -0.1))
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst1, g = 0.0435))
    push!(sst1.list_syn, gaba_syn(neuron_pre=vip1,neuron_post=sst1, g = -0.05))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E1,neuron_post=pv1, g = 0.0435))
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst1,neuron_post=pv1, g = -0.17))
    push!(pv1.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=pv1, g = -0.18))
    
    #connect to E2
    push!(dend2.list_syn, gaba_syn(τ = 10*0.001, neuron_pre=sst2,neuron_post=dend2, g = -0.09))
    push!(E2.list_syn, gaba_syn(neuron_pre=pv1,neuron_post=E2, g = -0.001))
    push!(E2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=E2, g = 0.18))
    
    push!(vip2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=vip2, g = 0.058))
    push!(vip2.list_syn, gaba_syn(neuron_pre=sst2,neuron_post=vip2, g = -0.1))
    
    push!(sst2.list_syn, nmda_syn(neuron_pre=E2,neuron_post=sst2, g = 0.0435))
    push!(sst2.list_syn, gaba_syn(neuron_pre=vip2,neuron_post=sst2, g = -0.05))
    
    push!(pv1.list_syn, nmda_syn(neuron_pre=E2,neuron_post=pv1, g = 0.0435))
    push!(pv1.list_syn, gaba_syn(neuron_pre=sst2,neuron_post=pv1, g = -0.17))
    
    
    push!(sst1.list_syn, nmda_syn(neuron_pre=E2,neuron_post=sst1, g = 0.0435))
    push!(sst2.list_syn, nmda_syn(neuron_pre=E1,neuron_post=sst2, g = 0.0435))

    
    push!(c.list_dend,dend1)
    push!(c.list_dend,dend2)

    push!(c.list_soma,E1)
    push!(c.list_soma,E2)

    push!(c.list_vip,vip1)
    push!(c.list_vip,vip2)

    push!(c.list_sst,sst1)
    push!(c.list_sst,sst2)

    push!(c.list_pv,pv1)


    e1 = wong_wang_cell()
    e2 = wong_wang_cell()
    e1.Istim[1] = 5.2*0.0001*30.0*1.1
    e1.Istim[1] = 5.2*0.0001*30.0*1.1
    e2.Istim[1] = 5.2*0.0001*30.0*0.9
    ww_network = wong_wang_network(threshold = 15.0)
    push!(c.nn, ww_network )
    push!(ww_network.list_units,e1)
    push!(ww_network.list_units,e2)

    push!(e1.list_syn,nmda_syn(neuron_pre=E1,neuron_post=e1, g = 0.3255/0.6))
    push!(e2.list_syn,nmda_syn(neuron_pre=E1,neuron_post=e2, g = 0.3255))

    push!(e2.list_syn,nmda_syn(neuron_pre=E2,neuron_post=e2, g = 0.3255/0.6))
    push!(e1.list_syn,nmda_syn(neuron_pre=E2,neuron_post=e2, g = 0.3255))


end


construct_modele_sean(circuit1)


circuit1




function simulation(param_simu:simulation_parameters)


end




temp = []
for i=1:1000
time_step()
#push!(temp,E1.r)
end



temp


using Plots

plotlyjs()


plot(E1.r)
plot!(E2.r)

plot(E1.Inoise)

plot(dend1.Ioutput)

plot(dendrite_input_output!.(dend1.Iexc,dend1.Iinh))

plot(dend1.Iexc)

plot(dend1.Iinh)

plot(pv1.r)

plot(temp)

circuit1.list_dend[1].Istim = 0.1
for i=1:1000
    time_step(circuit1,simu)
#    push!(temp,E1.Itot)
end
plot(circuit1.list_soma[1].r)
plot!(circuit1.list_soma[2].r)


plot(circuit1.nn[1].list_units[1].Iexc)

plot!(circuit1.nn[1].list_units[2].Iexc)


plot(circuit1.nn[1].list_units[1].r)
plot!(circuit1.nn[1].list_units[2].r)

plot(circuit1.nn[1].list_units[1].list_syn[1].s)

plot!(circuit1.nn[1].list_units[2].list_syn[1].s)


circuit1.list_dend[2].Istim = 0.1
circuit1.list_dend[1].Istim = 0.0

for i=1:1000
    time_step(circuit1,simu)
#    push!(temp,E1.Itot)
end
plot(circuit1.list_soma[1].r)
plot!(circuit1.list_soma[2].r)









dend1.Istim = 0.0
dend2.Istim = 0.1
for i=1:1000
    time_step()
#    push!(temp,E1.Itot)
end
plot(E1.r)
plot!(E2.r)


plot(dend1.Ioutput)
plot(dend2.Ioutput)



simu = simulation_parameters()

e1 = wong_wang_cell()
e2 = wong_wang_cell()
e1.Istim[1] = 5.2*0.0001*30.0*1.1
e1.Istim[1] = 5.2*0.0001*30.0*1.1
e2.Istim[1] = 5.2*0.0001*30.0*0.9

ww_network = wong_wang_network(threshold = 15.0)
push!(ww_network.list_units,e1)
push!(ww_network.list_units,e2)


function current_to_frequency!(n::wong_wang_cell)
# taken from Abott and Chance
temp = n.a * n.Itot[end] - n.b
push!(n.r , temp / (1.0 - exp(-n.d * temp)))
end


function update_s!(network::wong_wang_network,simu::simulation_parameters)


  
    unit1 = network.list_units[1]
    unit2 = network.list_units[2]

    # update syn current

    push!(unit1.Ibg, unit1.Iexc[end] + unit1.Iinh[end])
    push!(unit2.Ibg, unit2.Iexc[end] + unit2.Iinh[end])

    push!(unit1.Itot , network.jn11*unit1.s[end] - network.jn21*unit2.s[end] + unit1.Ibg[end] + unit1.Istim[end]  + unit1.Inoise[end])
    push!(unit2.Itot , network.jn22*unit2.s[end] - network.jn12*unit1.s[end] + unit2.Ibg[end] + unit2.Istim[end] + unit2.Inoise[end])
  
    # update noise
    dummy_time = simu.dt/network.tampa
    push!(unit1.Inoise, unit1.Inoise[end] + (dummy_time)*(-unit1.Inoise[end]) + sqrt(dummy_time)*network.noiseamp*randn())
    push!(unit2.Inoise, unit2.Inoise[end] + (dummy_time)*(-unit2.Inoise[end]) + sqrt(dummy_time)*network.noiseamp*randn())
  
    # update S variable and rate
   current_to_frequency!(unit1)
   current_to_frequency!(unit2)
   
    push!(unit1.s, unit1.s[end] + simu.dt*(-unit1.s[end]/network.tnmda + (1.0-unit1.s[end])*unit1.γ*unit1.r[end]))
    push!(unit2.s, unit2.s[end] + simu.dt*(-unit2.s[end]/network.tnmda + (1.0-unit2.s[end])*unit2.γ*unit2.r[end]))
  
  
  
end


function update_cells!(n::wong_wang_cell,sim::simulation_parameters)
# update eveything in the wong wang cell

push!(n.s , n.s[end] + sim.dt * (-n.s[end]/n.τ + (1.0 - n.s[end])*n.γ*n.r[end]))  
current_to_frequency!(n)
end
 


  

  # TODO: do the wang wong network alone

function time_simu(n::wong_wang_network,simu::simulation_parameters)

for i in 1:Int(simu.Tfin/simu.dt)

    update_s!(n,simu)

end
end

time_simu(ww_network,simu)



using Plots
plotlyjs()


plot(ww_network.list_units[1].r)
plot!(ww_network.list_units[2].r)