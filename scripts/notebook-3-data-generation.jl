using DrWatson

quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN","1-project-categorical-MMN")
using Revise
using Parameters
push!(LOAD_PATH,srcdir())
using JLD2
using MyNeurosciencePackage


using Plots
plotlyjs()


using JuMP
using Ipopt


function constant_stimulus(sim::simulation_parameters,value::Float64)
    sim.Stimulus = value.*ones(Int(sim.Tfin/sim.dt))
    
end


function alternating_stimulus(sim::simulation_parameters,value::Float64)
     list_stim = zeros(Int(sim.Tfin/sim.dt))
    temp_t = 0.0
    
    for i=1:length(list_stim)-1
    temp_t +=sim.dt

    if temp_t>sim.Tstimduration && i>sim.Tstimdebut/sim.dt
        if list_stim[i] == 0.0
            list_stim[i] = value
        else
            list_stim[i] = 0.0
        end
        temp_t = 0.0

    end
    list_stim[i+1] = list_stim[i]
end
    sim.Stimulus = list_stim
    
end


function make_dynamics(keywords::Array{String},sim::simulation_parameters;facilitation = false,stim=0.0, adaptation = false, τ = 0.1,save = false, save_parameters = Dict(), notebook_file = "" )
    #todo replace keywords by something nicer
    if keywords[2] == "constant"
        constant_stimulus(sim,stim)

    elseif keywords[2] == "alternating"
            alternating_stimulus(sim,stim)

    else
      #  warn("Don't forget to add a type of input to your neuron. Keywords[2]")
    end
    
    circuit = create_network(keywords[1]; facilitation,adaptation,save = save, save_parameters = save_parameters, notebook_file = notebook_file)
    circuit.list_soma[1].adaptation.τA = τ
    circuit.list_sst[1].list_syn[3].f_param.τ = τ

    full_time_dynamics(circuit, sim)
    
    return circuit

end




function analyze_alternating(network::microcircuit,sim::simulation_parameters)
    list_values = []
    PC=network.list_soma[1]
    dt = sim.dt
    temp_dt = sim.Tstimduration + sim.Tstimduration/2
    while temp_dt <= sim.Tfin
        push!(list_values, PC.r[Int(round(temp_dt/dt))])
        temp_dt+=2.0*sim.Tstimduration
    end
    
    return list_values
end



allparams= Dict("facilitation" => [true, false],
"adaptation" => [true,false],
            "stimlist" => collect(0.0:0.1:1.0),
                "τ" => collect(0.1:0.1:2.0),
                "Tinterstim" => collect(0.1:0.2:1.0)
                )

dicts = dict_list(allparams)



function makesim(d::Dict;save = false, save_parameters = Dict(), notebook_file = "")
    simulation_param = simulation_parameters()

    simulation_param.Tstimduration = d["Tinterstim"]
    fa = d["facilitation"]
    a=d["adaptation"]


circuit = make_dynamics(["integrator","alternating"],simulation_param; facilitation = fa,adaptation = a, stim = d["stimlist"], τ=d["τ"],save = save, save_parameters = save_parameters, notebook_file = notebook_file)
circuitdict = Dict("circuit"=>circuit, "simulation_param" => simulation_param)
    return circuitdict
end


# function makesim_temp(d::Dict)
#     simulation_param = simulation_parameters()
#     simulation_param.Tfin=simulation_param.dt
#     simulation_param.Tstimduration = d["Tinterstim"]
#     fa = d["facilitation"]
#     a=d["adaptation"]


#     circuit = make_dynamics(["integrator","alternating"],simulation_param; facilitation = fa,adaptation = a, stim = d["stimlist"], τ=d["τ"])
#     circuitdict = Dict("circuit"=>circuit, "simulation_param" => simulation_param)
#     return circuitdict
# end



for (i, d) in enumerate(dicts)
    f = makesim(d;save=true, notebook_file="notebook-3",save_parameters = d)
    save_dynamics(f["circuit"],"notebook-3", d)

end




###############
# Generation de certaines figures
####################	
	

############
# Fit functions
#############
using LsqFit
	
model_sig(x, p) = p[1].*sigmoid.(-(x.-p[4])./p[2]) .+ p[3]#p[1].*exp.(-x./p[2]) .+ p[3]
model_exp(x, p) = p[1].*exp.(-(x)./p[2]) .+ p[3]#p[1].*exp.(-x./p[2]) .+ p[3]

function get_param_fit(c::microcircuit, sim::simulation_parameters,model)
    lb = [1e-9,0.1,1.0,1e-9]
    #lb = [0.0,0.0,0.0,0.0]
       values_fit = analyze_alternating(c,sim)
    ub = [2.0*maximum(values_fit),6.0,maximum(values_fit),30.0]
    p0 = [0.5*maximum(values_fit),0.5,0.5*maximum(values_fit),1]

       xdata_stim = 1.5*sim.Tstimduration:2.0*sim.Tstimduration:sim.Tfin
    fit = curve_fit(model, xdata_stim, values_fit, p0, lower = lb, upper = ub)

		return xdata_stim, fit.param,values_fit,maximum(values_fit)
end





function fit_sig(x,y)
    m = Model(Ipopt.Optimizer)
    JuMP.register(m, :sigmoid, 1, sigmoid, autodiff=true)
    JuMP.register(m, :model_sig, 2, model_sig, autodiff=true)

    #m = Model(solver=ClpSolver())
    @variable(m, p[1:4] >= 1e-9)
    @variable(m, t >= 0.)
    #@NLconstraint(m, p[2] >= 1e-9);
    @NLconstraint(m, [i=1:length(y)],  p[1]*1/(1+exp(-(x[i]-p[4])/(-p[2]))) + p[3]- y[i] <= t);
    @NLconstraint(m, [i=1:length(y)], -t <= p[1]*1/(1+exp(-(x[i]-p[4])/(-p[2]))) + p[3] - y[i]);
    @NLobjective(m, Min, t)
    optimize!(m)
    return JuMP.value.(p)
end



# adaptation matrice
adaptation_params= Dict("facilitation" => [false],
"adaptation" => [true],
            "stimlist" => 0.5,#collect(0.1:0.1:1.0),
                "τ" => 0.1,
                "Tinterstim" => collect(0.1:0.2:1.0)
                )

adaptation_dict = dict_list(adaptation_params)

τ_list = collect(0.2:0.1:2.0)
stim_list = collect(0.1:0.1:1.0)


array_adaptation = zeros(4,length(stim_list),length(adaptation_params["Tinterstim"]),length(τ_list))
array_facilitation = zeros(4,length(stim_list),length(adaptation_params["Tinterstim"]),length(τ_list))
array_control = zeros(4,length(stim_list),length(adaptation_params["Tinterstim"]),length(τ_list))

for (k,stim) in enumerate(stim_list)
for (j,τ) in enumerate(τ_list)
for (i, d_temp) in enumerate(adaptation_dict)
    d_temp["τ"] = τ
    d_temp["stimlist"] = stim

    temp=savename(d_temp,"jld2")
    temp_load=    wload(datadir("simulations","notebook-3", savename(d_temp, "jld2")))
    circuit = temp_load["circuit"]
    simulation_param = temp_load["simulation_param"]
    fit_param=zeros(4)
    try
    xdata_stim, fit_param,m,n = get_param_fit(circuit,simulation_param,model_sig) 
    catch
        fit_param[1] = 0
        fit_param[2] = 0
        fit_param[3] = 0
        fit_param[4] = 0
    end

    p1 = fit_param[1]
    p2 = fit_param[2]
    p3 = fit_param[3]
    p4 = fit_param[4]

    array_adaptation[:,k,i,j] = fit_param

    d_control = copy(d_temp)
    d_control["adaptation"] = false

    temp=savename(d_control,"jld2")
    temp_load=    wload(datadir("simulations","notebook-3", savename(d_control, "jld2")))
    circuit2 = temp_load["circuit"]
    
    simulation_param2 = temp_load["simulation_param"]
    fit_param_control = zeros(4)
    try
    xdata_stim, fit_param_control,m_control,n = get_param_fit(circuit2,simulation_param2,model_sig) 
    catch
        
    end
    array_control[:,k,i,j] = fit_param_control
 

    d_facilitation = copy(d_temp)
    d_facilitation["adaptation"] = false
    d_facilitation["facilitation"] = true

    temp=savename(d_facilitation,"jld2")
    temp_load=    wload(datadir("simulations","notebook-3", savename(d_facilitation, "jld2")))
    circuit = temp_load["circuit"]
    simulation_param = temp_load["simulation_param"]
    fit_param_f = zeros(4)

    try
    xdata_stim, fit_param_f,m_f,n = get_param_fit(circuit,simulation_param,model_sig) 
    catch
    end
    
    array_facilitation[:,k,i,j] = fit_param_f
  
end
end
end

T_list = collect(0.1:0.2:1.0)



fig = plot()
heatmap!(fig,τ_list,T_list,array_facilitation[1,6,:,:]./array_facilitation[3,6,:,:])
xlabel!(fig,"τ (s)")
ylabel!(fig,"T interstim")
title!(fig,"p1/p3 for the facilitation network")
savefig(fig,plotsdir("3-Adaptation","fig-p1-p3-fixedstim-facilitation.png"))


fig_control = plot()
heatmap!(fig_control,τ_list,T_list,array_control[1,6,:,:]./array_control[3,6,:,:])
xlabel!(fig_control,"τ (s)")
ylabel!(fig_control,"T interstim")
title!(fig_control,"p1/p3 for the control network")
savefig(fig_control,plotsdir("3-Adaptation","fig-p1-p3-fixedstim-control.png"))


fig = plot()
heatmap!(fig,τ_list,T_list,array_facilitation[2,6,:,:])
xlabel!(fig,"τ (s)")
ylabel!(fig,"T interstim")
title!(fig,"p2 for the facilitation network")
savefig(fig,plotsdir("3-Adaptation","fig-p2-fixedstim-facilitation.png"))



fig_control = plot()
heatmap!(fig_control,τ_list,T_list,array_control[2,6,:,:])
xlabel!(fig_control,"τ (s)")
ylabel!(fig_control,"T interstim")
title!(fig_control,"p2 for the control network")
savefig(fig_control,plotsdir("3-Adaptation","fig-p2-fixedstim-control.png"))


fig = plot()
heatmap!(fig,τ_list,T_list,array_facilitation[4,6,:,:])
xlabel!(fig,"τ (s)")
ylabel!(fig,"T interstim")
title!(fig,"p4 for the facilitation network")
savefig(fig,plotsdir("3-Adaptation","fig-p4-fixedstim-facilitation.png"))



fig_control = plot()
heatmap!(fig_control,τ_list,T_list,array_control[4,6,:,:])
xlabel!(fig_control,"τ (s)")
ylabel!(fig_control,"T interstim")
title!(fig_control,"p4 for the control network")
savefig(fig_control,plotsdir("3-Adaptation","fig-p4-fixedstim-control.png"))


##########
# at tInterstim fixed
###########


fig = plot()
heatmap!(fig,τ_list,stim_list,array_facilitation[1,:,3,:]./array_facilitation[3,:,3,:])
xlabel!(fig,"τ (s)")
ylabel!(fig,"Stimlist")
title!(fig,"p1/p3 for the facilitation network")
savefig(fig,plotsdir("3-Adaptation","fig-p1-p3-fixedTinter-facilitation.png"))



fig = plot()
heatmap!(fig,τ_list,stim_list,array_control[1,:,3,:]./array_control[3,:,3,:])
xlabel!(fig,"τ (s)")
ylabel!(fig,"Stimlist")
title!(fig,"p1/p3 for the control network")
savefig(fig,plotsdir("3-Adaptation","fig-p1-p3-fixedTinter-control.png"))




#########
# at tau fixed
###########
fig = plot()
plot!(fig,stim_list,array_adaptation[1,:,3,1]./array_adaptation[3,:,3,1],w=2,color=:red,label="adapt")
plot!(fig,stim_list,array_facilitation[1,:,3,1]./array_facilitation[3,:,3,1],w=2,color=:green,label="facilitaiton")
plot!(fig,stim_list,array_control[1,:,3,1]./array_control[3,:,3,1],w=2,color=:blue,label="control")
xlabel!(fig,"Stimlist")
title!(fig,"p1/p3 at fixed τ and Tinter")
savefig(fig,plotsdir("3-Adaptation","fig-p1-p3-fixedtau-fixedTinter.png"))
















heatmap(τ_list,T_list,array_adaptation[1,6,:,:])
heatmap(τ_list,T_list,array_adaptation[3,5,:,:])
heatmap(τ_list,T_list,array_adaptation[2,5,:,:])
heatmap(τ_list,collect(0.3:0.2:1.0),array_adaptation[4,5,:,:])


heatmap(τ_list,T_list,array_adaptation[1,6,:,:]./array_adaptation[3,6,:,:])
heatmap(τ_list,T_list,array_control[1,6,:,:]./array_control[3,6,:,:])

heatmap(τ_list,collect(0.3:0.2:1.0),array_facilitation[4,5,:,:])

heatmap(τ_list,collect(0.3:0.2:1.0),array_control[4,5,:,:])


heatmap(τ_list,collect(0.1:0.2:1.0),array_p1_facilitation)

heatmap(τ_list,collect(0.1:0.1:1.0),array_p2_facilitation)

heatmap(τ_list,collect(0.1:0.1:1.0),array_p2_adaptation)


scatter(array_p1_adaptation)
scatter!(array_p1_facilitation)

scatter(array_p2_adaptation)
scatter!(array_p2_facilitation)

scatter(array_p3_adaptation)
scatter!(array_p3_facilitation)

scatter(array_p4_adaptation)
scatter!(array_p4_facilitation)



scatter(array_p3_adaptation)
scatter(array_p2_adaptation)


#The idea is to compare the relative p1 and p3 with respect to the absolute maximum that would be achieved by the neurons
# Want to have a 2D map of the results (stim, time constant) first

# issue is that adaptation really doesn't work! meaning I should juste have specific examples but the 2D map doesn't work for it.






d=Dict("facilitation" => false,
"adaptation" => true,
            "stimlist" =>0.6,
                "τ" => 0.1,
                "Tinterstim" => 0.3#collect(0.1:0.2:1.0)
                )
temp=savename(d,"jld2")
    
temp_load=    wload(datadir("simulations","notebook-3", savename(d, "jld2")))
circuit = temp_load["circuit"]
simulation_param = temp_load["simulation_param"]

xdata_stim, fit_param = get_param_fit(circuit,simulation_param,model_sig)


using Plots
plotlyjs()
plot(xdata_stim, model_sig(xdata_stim,fit_param))

function get_param_fit_temp(c::microcircuit, sim::simulation_parameters,model)

    lb = [0.0,0.0,0.0,0.0]
       values_fit = analyze_alternating(c,sim)
    ub = [2.0*maximum(values_fit),100.0,maximum(values_fit),Inf]
    p0 = [0.5*maximum(values_fit),1,0.5*maximum(values_fit),1]

       xdata_stim = 1.5*sim.Tstimduration:2.0*sim.Tstimduration:sim.Tfin
    fit = curve_fit(model, xdata_stim, values_fit, p0, lower = lb, upper = ub)

       return xdata_stim, fit.param,maximum(values_fit),values_fit
end

xdata_stim, fit_param,mtemp, values_temp = get_param_fit_temp(circuit,simulation_param,model_sig)

values_fit = analyze_alternating(circuit,simulation_param)
scatter(values_fit)
xdata_stim,y,values_temp,t = get_param_fit(circuit,simulation_param,model_sig)


plot(circuit.list_soma[1].r)

scatter(xdata_stim,values_temp)
plot!(xdata_stim, model_sig(xdata_stim,y))
plot!(xdata_stim, model_exp(xdata_stim,fit_param))
