using DrWatson
quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN","1-project-categorical-MMN")
push!(LOAD_PATH,srcdir())

using Revise
using MyNeurosciencePackage
using Parameters
using JLD2

###########
# test des noms de colonnes comme variables et symbols
##########


s = Symbol("foo")

temp = "foo2"
s2 = :($temp)

s2 = Symbol(temp)

s2 == :

using DataFrames



df = DataFrame()

df[!,temp] = [0]


df[!,s2]


#############
# test des sauvegardes des donnees
###########

allparams= Dict("facilitation" => false,
"adaptation" => false,
            "stimlist" => 0.0,
                "τ" => 0.1,
                "Tinterstim" => 0.1
                )

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
                

function makesim(d::Dict;save = false, save_parameters = Dict(), notebook_file = "")
                    simulation_param = simulation_parameters()
                
                    simulation_param.Tstimduration = d["Tinterstim"]
                    fa = d["facilitation"]
                    a=d["adaptation"]
                
                
  circuit = make_dynamics(["integrator","alternating"],simulation_param; facilitation = fa,adaptation = a, stim = d["stimlist"], τ=d["τ"],save = save, save_parameters = save_parameters, notebook_file = notebook_file)
  circuitdict = Dict("circuit"=>circuit, "simulation_param" => simulation_param)
                    return circuitdict
                end


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

f = makesim(allparams; save=true, notebook_file="notebook-3",save_parameters = allparams)

using DataFrames, CSV
save_dynamics(f["circuit"],"notebook-3", allparams)


temp = wload(datadir("simulations","notebook-3",savename(allparams,"jld2")))



temp2 = get_dynamics("notebook-3", allparams)


df = wload(datadir("simulations","notebook-3", savename(allparams, "csv")))

using Plots
plotlyjs()


plot(temp2.list_soma[1].r)