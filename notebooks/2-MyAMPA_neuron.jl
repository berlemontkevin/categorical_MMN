### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 2d2f220f-a6e4-45e2-823b-d2be87566d8e
begin
	using Plots, PlutoUI
	plotlyjs()
end

# ╔═╡ dfb4ae70-a37a-11eb-0310-dd0f3c12ba94
md"""
# Fit of dendritic neuron

In order to fit the mean dendritic voltage to the inhibitory and excitatory inputs, I use the following function:

`` y= \frac{(I_{exc} - (\beta * I_{inh} + \delta))}{(\zeta * I_{inh} + \mu)}``

`` <V_D> = \alpha(-0.5 +\frac{1}{1 + \exp(-y)}) + \gamma``

A few key points:
- Instead of the exponential, I just use a normalization by the inhibitory current. Fitting are similar (even better) and the formula seems more natural
- Using the sigmoind instead of the hyperbolic tangent leads to an accurate fit and remove the low excitation issue
- In average, the increase of firing rate dues to increase of inhibition was only happening at very low excitatory rates


$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\data\\sims\\code John\\fit-sigmoid_function.png"))

$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\data\\sims\\code John\\fit_inverse_sigmoid.png"))
"""

# ╔═╡ 86acf60a-1338-4487-a2f6-d70306655404
begin
	
	function sigmoid(x)
		
		return 1.0./(1.0.+exp.(-x))
	end
	
	function dendritic_voltage(p, exc, inh)
		# in mV
		mid_point = (inh .+ p[6]).*p[1]
		spread = p[2] .* inh .+ p[3]
		y = 2.0.*p[5].*(-0.5 .+ sigmoid.((exc .- mid_point)./spread)) .+p[4] .-70.0
		return y
	end


	param_fit = [ 6.48453644,  1.29942558,  8.32439742, 30.44449538, 26.47887547,
        6.50400446]
	
	
	inh_constant = 0.08 #multiplciative ocnstant to look at currents
	exc_constant = 7.82 # same for exc constant
	
	inh_list = LinRange(0.0,100,11)*inh_constant
	exc_list = LinRange(0.25,10.0,20) * exc_constant
	
	
	vd = zeros(length(exc_list),length(inh_list))
	for i=1:length(exc_list)
		vd[i,:] = dendritic_voltage(param_fit, exc_list[i], inh_list)
	end
	plot(exc_list,vd,color=:blues, line_z=(1:length(inh_list))', linewidth=3, colorbar=false)
	plot!(legend = false)
	xlabel!("gE tot")
	ylabel!("Mean dendritic voltage")
	title!("VD for many inhibitory rates")
end

# ╔═╡ 715344bd-caa4-447d-b0ba-aa0536407a05
begin 
	plot(inh_list,vd',color=:reds, line_z=(1:length(exc_list))', linewidth=3, 		colorbar=false)
	plot!(legend = false)
	xlabel!("gI tot")
	ylabel!("Mean dendritic voltage")
	title!("VD for many excitatory inputs")

end

# ╔═╡ b99525e1-651b-406c-aa6d-3f737d0606e8
begin 
	plot(inh_list,vd[1:3,:]',color=:reds, line_z=(1:length(3))', linewidth=3, 		colorbar=false)
	plot!(legend = false)
	xlabel!("gI tot")
	ylabel!("Mean dendritic voltage")
	title!("VD for many excitatory inputs at low exc")

end

# ╔═╡ 7a78bdea-339b-4fca-9cc7-f04ced919834
md"""
Now it's time to look at the 2D f-I curve of this kind of dendritic neuron to be sure it's behaving accordingly to the few experimental data


The first thing is to see if the behavior of the somatic firing rate is non-linear with the assymetry in dendritic and somatic inputs.

Due to the assymetry in excitation and inhibition, it's necessary to work at fixed inhibition toward the distal dendrites.

"""

# ╔═╡ ca3fad61-e712-46c8-a84a-34ed6fbd40d7
begin
	
	function firing_rate(vd,Isoma)
		Idend = 8.0.*(vd .+ 55.0)
		Itot = Isoma .+ Idend
		fr = ((Itot.+74.86).*(Itot.>-74.86)./45.16).^2.89 
	end
	
	
	Inh = 4.0 #let's define the constant inhibitory input sent to the dendrite
	
	vd_list = zeros(length(exc_list))
	vd_list = dendritic_voltage(param_fit, exc_list, Inh)
	Isoma_list = LinRange(-150.0,50.0,100)
	
	fIcurve_2D = zeros(length(vd_list),length(Isoma_list))
	for i=1:length(Isoma_list)
		fIcurve_2D[:,i] = firing_rate(vd_list,Isoma_list[i])
	end
	Idend_list = 8.0.*(vd_list .+ 55.0)
	plot(Isoma_list,Idend_list,fIcurve_2D,st=:surface)
	xlabel!("Current to soma")
	ylabel!("Current from dendrite")
end

# ╔═╡ 613026c8-9ec1-408a-8e31-77eaed5aa8cd
md"""
Now let's study a simplified model of the dendritic model. The mean dendritic voltage is directly inserted into current from dendrite to soma just to ismplified the writing.

`` y= \frac{(I_{exc} - (\beta * I_{inh} + \delta))}{(\zeta * I_{inh} + \mu)}``

`` I_{dend \rightarrow soma} = \alpha(-0.5 +\frac{1}{1 + \exp(-y)}) + \gamma``

"""

# ╔═╡ 9620592f-1c3a-4911-93d5-d7b11f33e210
begin
	
	FIsean(Iexc,Iinh) = 0.12.*tanh((Iexc .+ 7.0.*Iinh .+ 0.0)./(0.00964.*exp.(-Iinh./0.02))) .+ 0.13624
 
	AbottChance(x) = (135*x -54)/(1.0-exp(-0.308*(135*x-54)))



	FIseansigmoid(Iexc,Iinh) = 0.12.*2.0.*(-0.5.+sigmoid.((Iexc .+ 7.0.*Iinh .+ 0.0)./(0.00964.*exp.(-Iinh./0.02)))) .+ 0.13624
	
	FIseansigmoidv2(Iexc,Iinh) = 0.20.*1.0.*(-0.5.+sigmoid.((Iexc .+ 7.0.*Iinh .+ 0.1)./(-0.00964.*50.0.*Iinh.+0.00964))) .+ 0.11#19624
end

# ╔═╡ 485b4d95-6795-4cc5-9f68-64fe73a3400d
	@bind Iexc_cst Slider(0.0:0.1:1.0)

# ╔═╡ 0cb4322b-1e60-4a49-a4bc-f7875dd937ce
begin

	Iinh = -0.5:0.001:0.0
	plot(Iinh,AbottChance.(0.4.+FIsean.(Iexc_cst,Iinh)),linewidth=2,label = Iexc_cst)
	plot!(Iinh,AbottChance.(0.4.+FIseansigmoid.(Iexc_cst,Iinh)),linewidth=2,label = "sigmoid")
	plot!(Iinh,AbottChance.(0.4.+FIseansigmoidv2.(Iexc_cst,Iinh)),linewidth=2,label = Iexc_cst)
	#plot!(Iinh,AbottChance.(0.01.+0.0.*Iinh),linewidth=2,label = Iexc_cst)
	
end

# ╔═╡ 783d1513-4710-44dd-b05f-7661169003c0
begin
	using LsqFit
	m(t,p) = 1.0./(-p[1]*t.+p[2])
	p0 = [0.5, 0.001]
	tdata = Iinh
	ydata = 1.0./(0.00964.*exp.(-Iinh./0.02))
	
	fit = curve_fit(m, tdata, ydata, p0)
	
	pfit = fit.param
	
	plot(Iinh,1.0./(-pfit[1].*Iinh.+pfit[2]),linewidth=2,label = "fit")
	plot!(Iinh,1.0./(0.00964.*exp.(-Iinh./0.02)),linewidth=2,label = "original")
	
	
end

# ╔═╡ 8a2fe1c2-ebdf-4e92-9550-34099617434c
begin

	plot(Iinh,1.0./(-0.00964.*50.0.*Iinh.+0.00964),linewidth=2,label = "lin")
	plot!(Iinh,1.0./(0.00964.*exp.(-Iinh./0.02)),linewidth=2,label = Iexc_cst)
	
	
end

# ╔═╡ 024f8b3e-175f-4551-81eb-809cd0fa7908
begin

	plot(Iinh,(FIsean.(Iexc_cst,Iinh)),linewidth=2,label = Iexc_cst)
	plot!(Iinh,(FIseansigmoid.(Iexc_cst,Iinh)),linewidth=2,label = Iexc_cst)
	plot!(Iinh,(FIseansigmoidv2.(Iexc_cst,Iinh)),linewidth=2,label = Iexc_cst)
	#plot!(Iinh,AbottChance.(0.01.+0.0.*Iinh),linewidth=2,label = Iexc_cst)
	
end

# ╔═╡ 094b97dd-fecb-42ff-be1b-fdb1ff12b83e
md"""
2D fI curve for the sigmoid rate model

	
"""

# ╔═╡ 281e0487-1770-44d8-80d2-2939e5953698
@bind cst_Iinh Slider(-0.5:0.05:0.0)

# ╔═╡ cf8fd2d7-457f-49a3-bce3-3d58269ca125
begin
	#cst_Iinh = -0.1
	soma_list = 0.0:0.01:0.6
	exc_dend_list = 0.0:0.01:1.7
	
	sigmoid_2D = zeros(length(soma_list),length(exc_dend_list))
	
	for i=1:length(soma_list)
		for j=1:length(exc_dend_list)
			sigmoid_2D[i,j] = AbottChance(soma_list[i] + FIseansigmoidv2.(exc_dend_list[j],cst_Iinh))
		end
	end
	
	plot(exc_dend_list,soma_list,sigmoid_2D,st=:surface, c=:blues)
	xlabel!("Current to dend")
	ylabel!("Current to soma")
	title!("Firing rate (Hz)")
	
end

# ╔═╡ 0a043d7a-c3b3-4f4d-bb79-dbe73eb12640
md"""
Fitting the exponential funciton in Sean's paper and the linear sigmoid one
"""

# ╔═╡ 5d7545ab-29a9-4efc-bd62-f2a0ba6d27e5
begin
	# Fitting the 2D dendritic current
	Fitsigmoid(p,x) = p[1].*(-0.5.+sigmoid.((x[1] .+ p[2].*x[2] .+ 0.0)./(-p[3].*Iinh.+p[4]))) .+ p[5]
	
	tdata2D = [exc_dend_list, Iinh]
	#ydata = 
	
end

# ╔═╡ b1b410ba-1984-4cc2-b86c-4bdbd7442081
begin 
	
	Sigmoid2D_withparam(p,Iexc,Iinh) =  0.12.*1.0.*(-0.5.+sigmoid.((Iexc .+ 7.0.*Iinh .+ 0.0)./(-p[1].*Iinh.+p[2]))) .+ 0.15624
	
	plot(Iinh,AbottChance.(0.4.+Sigmoid2D_withparam(pfit,Iexc_cst,Iinh)),linewidth=2,label = Iexc_cst)
	
	
end

# ╔═╡ c45193ca-dac7-4954-8f96-6bffedb491cd
@bind cst_Iinh_fit Slider(-0.5:0.05:0.0)

# ╔═╡ dceee7cd-0793-4a0a-b843-82b91f41f9e6
md"""
$(LocalResource("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN\\_research\\0-fig-divers\\larkum2019.png"))
"""

# ╔═╡ 89109d9c-4b6d-4ee5-8dbe-4239e1d14565
begin
	#cst_Iinh = -0.1
	#soma_list = 0.0:0.01:0.6
	#e#xc_dend_list = 0.0:0.01:1.7
	
	sigmoid_2DFit = zeros(length(soma_list),length(exc_dend_list))
	
	for i=1:length(soma_list)
		for j=1:length(exc_dend_list)
			sigmoid_2DFit[i,j] = AbottChance(soma_list[i] + Sigmoid2D_withparam(pfit,exc_dend_list[j],cst_Iinh_fit))
		end
	end
	
	plot(exc_dend_list,soma_list,sigmoid_2DFit,st=:surface, c=:blues)
	xlabel!("Current to dend")
	ylabel!("Current to soma")
	title!("Firing rate for fit(Hz)")
	
end

# ╔═╡ c2ac1c3e-e4f9-4fe7-af3d-e5f00ea9ba09
begin
	current_2DFit = zeros(length(Iinh),length(exc_dend_list))
	
	for i=1:length(Iinh)
		for j=1:length(exc_dend_list)
			current_2DFit[i,j] = Sigmoid2D_withparam(pfit,exc_dend_list[j],Iinh[i])
		end
	end
	
	plot(exc_dend_list,Iinh,current_2DFit,st=:surface, c=:blues)
	xlabel!("Current Exc to dend")
	ylabel!("Current Inh to dend")
	title!("Current sent to soma (sigmoid)")
	
end

# ╔═╡ 267bd8a8-b4af-43e7-8034-6b47cb46de2d
begin
	current_2DNoFit = zeros(length(Iinh),length(exc_dend_list))
	
	for i=1:length(Iinh)
		for j=1:length(exc_dend_list)
			current_2DNoFit[i,j] = FIseansigmoidv2(exc_dend_list[j],Iinh[i])
		end
	end
	
	plot(exc_dend_list,Iinh,current_2DNoFit,st=:surface, c=:blues)
	xlabel!("Current Exc to dend")
	ylabel!("Current Inh to dend")
	title!("Current sent to soma (sigmoid no exact fit)")
	
end

# ╔═╡ 8211c947-c072-4c98-bbe7-353aef99f10b
begin
	current_2DSean = zeros(length(Iinh),length(exc_dend_list))
	
	for i=1:length(Iinh)
		for j=1:length(exc_dend_list)
			current_2DSean[i,j] = FIsean(exc_dend_list[j],Iinh[i])
		end
	end
	
	plot(exc_dend_list,Iinh,current_2DSean,st=:surface, c=:blues)
	xlabel!("Current Exc to dend")
	ylabel!("Current Inh to dend")
	title!("Current sent to soma (Sean)")
	
end

# ╔═╡ b8a93ac9-e576-41ea-927a-7791531eefc7
md"""
To conclude:
- The sigmoid function can be fitted succesfully
- Mathematicla analysis explains why we don't have the issues of the previous function
- In addition the main points in Sean's parameters are kept

If the goal is to have a simpler function (with implicit dendritic voltage use), then I suggest the following (parameters taken from fit and adaptation of Robert or Sean papers):

`` I_{dend \rightarrow soma} = p[1] * \left(-0.5 + \sigma \left( \frac{I_{exc} + p[2] I_{inh}}{p[3] * I_{inh} + p[4]}\right) \right) + p[5]``

with ``\sigma`` the sigmoid function.

Parameters are:
- `` p[1] = 0.12 ``
- `` p[1] = 7.0``
- `` p[1] = -0.482``
- `` p[4] = 0.00964``
- `` p[5] = 0.19624``



The firing rate function is then the standard Abott and Chance one.
"""

# ╔═╡ 50e93fcc-b359-4393-8038-666671536e3b


# ╔═╡ Cell order:
# ╟─dfb4ae70-a37a-11eb-0310-dd0f3c12ba94
# ╟─2d2f220f-a6e4-45e2-823b-d2be87566d8e
# ╟─86acf60a-1338-4487-a2f6-d70306655404
# ╟─715344bd-caa4-447d-b0ba-aa0536407a05
# ╟─b99525e1-651b-406c-aa6d-3f737d0606e8
# ╟─7a78bdea-339b-4fca-9cc7-f04ced919834
# ╟─ca3fad61-e712-46c8-a84a-34ed6fbd40d7
# ╠═613026c8-9ec1-408a-8e31-77eaed5aa8cd
# ╠═9620592f-1c3a-4911-93d5-d7b11f33e210
# ╠═485b4d95-6795-4cc5-9f68-64fe73a3400d
# ╠═0cb4322b-1e60-4a49-a4bc-f7875dd937ce
# ╟─8a2fe1c2-ebdf-4e92-9550-34099617434c
# ╟─024f8b3e-175f-4551-81eb-809cd0fa7908
# ╠═094b97dd-fecb-42ff-be1b-fdb1ff12b83e
# ╠═281e0487-1770-44d8-80d2-2939e5953698
# ╟─cf8fd2d7-457f-49a3-bce3-3d58269ca125
# ╟─0a043d7a-c3b3-4f4d-bb79-dbe73eb12640
# ╠═783d1513-4710-44dd-b05f-7661169003c0
# ╠═5d7545ab-29a9-4efc-bd62-f2a0ba6d27e5
# ╠═b1b410ba-1984-4cc2-b86c-4bdbd7442081
# ╠═c45193ca-dac7-4954-8f96-6bffedb491cd
# ╠═dceee7cd-0793-4a0a-b843-82b91f41f9e6
# ╠═89109d9c-4b6d-4ee5-8dbe-4239e1d14565
# ╠═c2ac1c3e-e4f9-4fe7-af3d-e5f00ea9ba09
# ╟─267bd8a8-b4af-43e7-8034-6b47cb46de2d
# ╠═8211c947-c072-4c98-bbe7-353aef99f10b
# ╠═b8a93ac9-e576-41ea-927a-7791531eefc7
# ╠═50e93fcc-b359-4393-8038-666671536e3b
