using DrWatson
@quickactivate "1-project-categorical-MMN"
using Documenter

push!(LOAD_PATH,srcdir())

include(srcdir("MyNeurosciencePackage.jl"))
# include(srcdir("NeuronalStructures.jl"))
# include(srcdir("DynamicsFunction.jl"))
# include(srcdir("NetworkConstruction.jl"))
# include(srcdir("BasicFunctions.jl"))
# include(srcdir("Plasticity.jl"))


# using MyNeurosciencePackage


using MyNeurosciencePackage.NeuronalStructures, MyNeurosciencePackage.BasicFunctions, 
MyNeurosciencePackage.DynamicsFunction, MyNeurosciencePackage.NetworkConstruction



makedocs(format = Documenter.HTML(prettyurls = false),
modules = [NeuronalStructures, DynamicsFunction, NetworkConstruction],
sitename="My Documentation",
authors ="kevin Berlemont",
pages =[
    "Home" => "index.md",
    "Structures" => "NeuronalStructures.md",
    "Dynamical Functions" => "DynamicsFunction.md",
    "Construction of the networks" => "NetworkConstruction.md",
    "API" => "api.md"
])