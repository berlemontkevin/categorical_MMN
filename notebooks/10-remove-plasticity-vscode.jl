using DrWatson
	
quickactivate("C:\\Users\\kevin\\OneDrive\\1-Research_projects\\1-project-categorical-MMN","1-project-categorical-MMN")
push!(LOAD_PATH,srcdir())

using Revise

using MyNeurosciencePackage


	using WGLMakie
	using Parameters, JLD2, Statistics	
	using DataFrames, PlutoUI
	#using LsqFit
	using ColorSchemes, StaticArrays
	
	set_theme!(theme_light())


lines(1:10)