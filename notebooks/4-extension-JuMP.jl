# this file is an extension of the notebook 4-adaptation theoretical

############################
# Set up of the model of the local microcircuit
############################
using DrWatson
quickactivate("C:\\Users\\kevin\\Documents\\MyDocuments\\3-NYU\\1-project-categorical-MMN","1-project-categorical-MMN")

using JuMP
using Ipopt

function f(a)
    if a > 0.0
        return a
    else
        return 0.0
    end

end

function sigmoid(x)
    return 1.0./(1.0.+ exp.(-x))
end

function g(Iexc,Iinh)
    c1 = 0.2
    c2 =  -7.0
    c3 = -0.482
    c4 = 0.00964
    c5 = 0.11
    c6 = 0.1
    y = (Iexc - c2*Iinh + c6 )/(c3*Iinh + c4)
    
    return  c1*(-0.5 + sigmoid(y))+c5 
end

param=Dict(
"τ_AMPA" => 0.002,
"τ_GABA" => 0.005,
"τ_int" => 0.8,
"τ_GABAdend" => 0.01,
"τ_NMDA" => 0.06,

"r0PV" => -95.0,
"r0SST" => -33.0,
"r0VIP" => -33.0,

"cPV" =>330.0,
"cSST" => 132.0,
"cVIP" => 132.0,

"c1" => 0.2,
"c2" =>  -7.0,
"c3" => -0.482,
"c4" => 0.00964,
"c5" => 0.11,
"c6" => 0.1,

"τadapt" => 1.0,

"Ibg_soma" =>  0.15,
"Ibg_dend" => 0.1,
"Ibg_pv" => 0.3,
"Ibg_sst" => 0.2,
"I_bg_vip" => 0.3,
"Ibg_int"=>0.310,

"a" => 135.0,
"b" => 54.0,
"d" => 0.308,

"γnmda" =>0.641*2,
"γgaba" => 2.0,


"g_ssttodend" => -0.09,
"g_pvtopc" => -0.005,
"g_pctopc" => 0.18,
"g_pctovip" => 0.058,
"g_ssttovip" => -0.1,
"g_pctosst" => 0.0435,
"g_viptosst" => -0.05,
"g_pctopv" => 0.08435,
"g_ssttopv" => -0.17,
"g_pvtopv" => -0.18,
"g_inttosst" => 0.22,
"g_pctoint" => 0.15,

"α" => -0.2,


"τA" => 0.1,
"gA" => -0.1,

"τ_facilitation" => 1.5,
"U_inttosst" => 0.2

)




function standard_circuit_solution(Istim::Float64, param)
    ### constant values
    τ_AMPA = param["τ_AMPA"]
    τ_GABA = param["τ_GABA"]
    τ_int = param["τ_int"]
    τ_GABAdend = param["τ_GABAdend"]
    τ_NMDA = param["τ_NMDA"]

    r0PV = param["r0PV"]
    r0SST = param["r0SST"]
    r0VIP = param["r0VIP"]

    cPV = param["cPV"]
    cSST = param["cSST"]
    cVIP = param["cVIP"]

    c1 = param["c1"]
    c2 =  param["c2"]
    c3 = param["c3"]
    c4 = param["c4"]
    c5 = param["c5"]
    c6 = param["c6"]

    τadapt = param["τadapt"]

    Ibg_soma =  param["Ibg_soma"]
    Ibg_dend = param["Ibg_dend"]
    Ibg_pv = param["Ibg_pv"]
    Ibg_sst = param["Ibg_sst"]
    I_bg_vip = param["I_bg_vip"]
    Ibg_int= param["Ibg_int"]

    a = param["a"]
    b = param["b"]
    d = param["d"]

    γnmda = param["γnmda"]
    γgaba = param["γgaba"]


    g_ssttodend = param["g_ssttodend"]
    g_pvtopc = param["g_pvtopc"]
    g_pctopc = param["g_pctopc"]
    g_pctovip = param["g_pctovip"]
    g_ssttovip = param["g_ssttovip"]
    g_pctosst = param["g_pctosst"]
    g_viptosst = param["g_viptosst"]
    g_pctopv = param["g_pctopv"]
    g_ssttopv = param["g_ssttopv"]
    g_pvtopv = param["g_pvtopv"]
    g_inttosst = param["g_inttosst"]
    g_pctoint = param["g_pctoint"]



    α = param["α"]
    ######################


    ####
    # variables declarations

    microcircuit = Model(Ipopt.Optimizer)

    @variable(microcircuit, Rsst >= 0.0)
    @variable(microcircuit, Rpv >= 0.0)
    @variable(microcircuit, Rvip >= 0.0)
    @variable(microcircuit, Rpc >= 0.0)
    @variable(microcircuit, Rint >= 0.0)

    @variable(microcircuit, sssttodend >= 0.0)
    @variable(microcircuit, spc >= 0.0)
    @variable(microcircuit, ssst >= 0.0)
    @variable(microcircuit, svip >= 0.0)
    @variable(microcircuit, spv >= 0.0)
    @variable(microcircuit, sint >= 0.0)

    @variable(microcircuit, Isoma)
    @variable(microcircuit, Idend )
    @variable(microcircuit, Iexc >= 0.0)
    @variable(microcircuit, Iinh <= 0.0)
    @variable(microcircuit, Ivip)
    @variable(microcircuit, Isst)
    @variable(microcircuit, Ipv)
    @variable(microcircuit,Iint)
    ############


    register(microcircuit, :f, 1, f; autodiff = true)
    register(microcircuit, :sigmoid, 1, sigmoid; autodiff = true)
    register(microcircuit, :g, 2, g; autodiff = true)
    ##############
    # Constraints for solving





    #s Constraints
    @NLconstraint(microcircuit, cssst, ssst == τ_GABA*γgaba*Rsst)
    @NLconstraint(microcircuit, csvip, svip == τ_GABA*γgaba*Rvip)
    @NLconstraint(microcircuit, cspv, spv == τ_GABA*γgaba*Rpv)
    @NLconstraint(microcircuit, cssstdend, sssttodend == τ_GABAdend*γgaba*Rsst)

    @NLconstraint(microcircuit, cspc, 0.0 == -spc/τ_NMDA + γnmda*(1.0 - spc)*Rpc)
    @NLconstraint(microcircuit, csint, 0.0 == -sint/τ_NMDA + γnmda*(1.0 - sint)*Rint)

    # firing rate constraints
    @NLconstraint(microcircuit, crsst, -Rsst + f(cSST*Isst + r0SST)==0.0)
    @NLconstraint(microcircuit, crpv, -Rpv + f(cPV*Ipv + r0PV)==0.0)
    @NLconstraint(microcircuit, crvip, -Rvip + f(cVIP*Ivip + r0VIP) ==0.0)
    @NLconstraint(microcircuit, csoma, -Rpc + (a*Isoma-b)/(1.0 -exp(-d*(a*Isoma -b))) == 0.0)
    @NLconstraint(microcircuit,cint, α*Rint + (a*Iint-b)/(1.0 -exp(-d*(a*Iint -b))) == 0.0)

    #current Constraints


    @NLconstraint(microcircuit, cIdend, Idend == g(Iexc,Iinh))
    @NLconstraint(microcircuit, cIexcdend, Iexc == Ibg_dend)
    @NLconstraint(microcircuit, cIinhdend,  Iinh == g_ssttodend*sssttodend)
    @NLconstraint(microcircuit, cIsst, Isst == Ibg_sst + g_viptosst*svip + g_inttosst*sint + g_pctosst*spc)
    @NLconstraint(microcircuit, cIpv, Ipv == Ibg_pv + g_pctopv*spc + g_ssttopv*ssst +g_pvtopv*spv)
    @NLconstraint(microcircuit, cIvip, Ivip == I_bg_vip + g_ssttovip*ssst +g_pctovip*spc)
    @NLconstraint(microcircuit, cIsoma, Isoma == Ibg_soma + Istim + Idend + g_pvtopc*spv +g_pctopc*spc)
    @NLconstraint(microcircuit, cIint, Iint == g_pctoint*spc + Ibg_int)



    ###########


    optimize!(microcircuit)

return value(Rpc),value(Rint),value(Rvip),value(Rsst),value(Rpv)
end 



function adaptation_circuit_solution(Istim::Float64, param)
    ### constant values
    τ_AMPA = param["τ_AMPA"]
    τ_GABA = param["τ_GABA"]
    τ_int = param["τ_int"]
    τ_GABAdend = param["τ_GABAdend"]
    τ_NMDA = param["τ_NMDA"]

    r0PV = param["r0PV"]
    r0SST = param["r0SST"]
    r0VIP = param["r0VIP"]

    cPV = param["cPV"]
    cSST = param["cSST"]
    cVIP = param["cVIP"]

    c1 = param["c1"]
    c2 =  param["c2"]
    c3 = param["c3"]
    c4 = param["c4"]
    c5 = param["c5"]
    c6 = param["c6"]

    τadapt = param["τadapt"]

    Ibg_soma =  param["Ibg_soma"]
    Ibg_dend = param["Ibg_dend"]
    Ibg_pv = param["Ibg_pv"]
    Ibg_sst = param["Ibg_sst"]
    I_bg_vip = param["I_bg_vip"]
    Ibg_int= param["Ibg_int"]

    a = param["a"]
    b = param["b"]
    d = param["d"]

    γnmda = param["γnmda"]
    γgaba = param["γgaba"]


    g_ssttodend = param["g_ssttodend"]
    g_pvtopc = param["g_pvtopc"]
    g_pctopc = param["g_pctopc"]
    g_pctovip = param["g_pctovip"]
    g_ssttovip = param["g_ssttovip"]
    g_pctosst = param["g_pctosst"]
    g_viptosst = param["g_viptosst"]
    g_pctopv = param["g_pctopv"]
    g_ssttopv = param["g_ssttopv"]
    g_pvtopv = param["g_pvtopv"]
    g_inttosst = param["g_inttosst"]
    g_pctoint = param["g_pctoint"]



    α = param["α"]

    τA = param["τA"]
    gA = param["gA"]
    ######################


    ####
    # variables declarations

    microcircuit = Model(Ipopt.Optimizer)

    @variable(microcircuit, Rsst >= 0.0)
    @variable(microcircuit, Rpv >= 0.0)
    @variable(microcircuit, Rvip >= 0.0)
    @variable(microcircuit, Rpc >= 0.0)
    @variable(microcircuit, Rint >= 0.0)

    @variable(microcircuit, sssttodend >= 0.0)
    @variable(microcircuit, spc >= 0.0)
    @variable(microcircuit, ssst >= 0.0)
    @variable(microcircuit, svip >= 0.0)
    @variable(microcircuit, spv >= 0.0)
    @variable(microcircuit, sint >= 0.0)

    @variable(microcircuit, Isoma)
    @variable(microcircuit, Idend )
    @variable(microcircuit, Iexc >= 0.0)
    @variable(microcircuit, Iinh <= 0.0)
    @variable(microcircuit, Ivip)
    @variable(microcircuit, Isst)
    @variable(microcircuit, Ipv)
    @variable(microcircuit,Iint)

    @variable(microcircuit,sA)

    ############


    register(microcircuit, :f, 1, f; autodiff = true)
    register(microcircuit, :sigmoid, 1, sigmoid; autodiff = true)
    register(microcircuit, :g, 2, g; autodiff = true)
    ##############
    # Constraints for solving





    #s Constraints
    @NLconstraint(microcircuit, cssst, ssst == τ_GABA*γgaba*Rsst)
    @NLconstraint(microcircuit, csvip, svip == τ_GABA*γgaba*Rvip)
    @NLconstraint(microcircuit, cspv, spv == τ_GABA*γgaba*Rpv)
    @NLconstraint(microcircuit, cssstdend, sssttodend == τ_GABAdend*γgaba*Rsst)

    @NLconstraint(microcircuit, cspc, 0.0 == -spc/τ_NMDA + γnmda*(1.0 - spc)*Rpc)
    @NLconstraint(microcircuit, csint, 0.0 == -sint/τ_NMDA + γnmda*(1.0 - sint)*Rint)

    # firing rate constraints
    @NLconstraint(microcircuit, crsst, -Rsst + f(cSST*Isst + r0SST)==0.0)
    @NLconstraint(microcircuit, crpv, -Rpv + f(cPV*Ipv + r0PV)==0.0)
    @NLconstraint(microcircuit, crvip, -Rvip + f(cVIP*Ivip + r0VIP) ==0.0)
    @NLconstraint(microcircuit, csoma, -Rpc + (a*Isoma-b)/(1.0 -exp(-d*(a*Isoma -b))) == 0.0)
    @NLconstraint(microcircuit,cint, α*Rint + (a*Iint-b)/(1.0 -exp(-d*(a*Iint -b))) == 0.0)

    #current Constraints


    @NLconstraint(microcircuit, cIdend, Idend == g(Iexc,Iinh))
    @NLconstraint(microcircuit, cIexcdend, Iexc == Ibg_dend)
    @NLconstraint(microcircuit, cIinhdend,  Iinh == g_ssttodend*sssttodend)
    @NLconstraint(microcircuit, cIsst, Isst == Ibg_sst + g_viptosst*svip + g_inttosst*sint + g_pctosst*spc)
    @NLconstraint(microcircuit, cIpv, Ipv == Ibg_pv + g_pctopv*spc + g_ssttopv*ssst +g_pvtopv*spv)
    @NLconstraint(microcircuit, cIvip, Ivip == I_bg_vip + g_ssttovip*ssst +g_pctovip*spc)
    @NLconstraint(microcircuit, cIsoma, Isoma == Ibg_soma + Istim + Idend + g_pvtopc*spv +g_pctopc*spc + gA*sA)
    @NLconstraint(microcircuit, cIint, Iint == g_pctoint*spc + Ibg_int)


    #adaptation
    @NLconstraint(microcircuit, csA, sA == Rpc*τA)
    ###########


    optimize!(microcircuit)

return value(Rpc),value(Rint),value(Rvip),value(Rsst),value(Rpv)

end 


function facilitation_circuit_solution(Istim::Float64, param)
    ### constant values
    τ_AMPA = param["τ_AMPA"]
    τ_GABA = param["τ_GABA"]
    τ_int = param["τ_int"]
    τ_GABAdend = param["τ_GABAdend"]
    τ_NMDA = param["τ_NMDA"]

    r0PV = param["r0PV"]
    r0SST = param["r0SST"]
    r0VIP = param["r0VIP"]

    cPV = param["cPV"]
    cSST = param["cSST"]
    cVIP = param["cVIP"]

    c1 = param["c1"]
    c2 =  param["c2"]
    c3 = param["c3"]
    c4 = param["c4"]
    c5 = param["c5"]
    c6 = param["c6"]

    τadapt = param["τadapt"]

    Ibg_soma =  param["Ibg_soma"]
    Ibg_dend = param["Ibg_dend"]
    Ibg_pv = param["Ibg_pv"]
    Ibg_sst = param["Ibg_sst"]
    I_bg_vip = param["I_bg_vip"]
    Ibg_int= param["Ibg_int"]

    a = param["a"]
    b = param["b"]
    d = param["d"]

    γnmda = param["γnmda"]
    γgaba = param["γgaba"]


    g_ssttodend = param["g_ssttodend"]
    g_pvtopc = param["g_pvtopc"]
    g_pctopc = param["g_pctopc"]
    g_pctovip = param["g_pctovip"]
    g_ssttovip = param["g_ssttovip"]
    g_pctosst = param["g_pctosst"]
    g_viptosst = param["g_viptosst"]
    g_pctopv = param["g_pctopv"]
    g_ssttopv = param["g_ssttopv"]
    g_pvtopv = param["g_pvtopv"]
    g_inttosst = param["g_inttosst"]
    g_pctoint = param["g_pctoint"]



    α = param["α"]

    τA = param["τA"]
    gA = param["gA"]

    U_inttosst = param["U_inttosst"]
    τ_facilitation = param["τ_facilitation"]



    ######################


    ####
    # variables declarations

    microcircuit = Model(Ipopt.Optimizer)

    @variable(microcircuit, Rsst >= 0.0)
    @variable(microcircuit, Rpv >= 0.0)
    @variable(microcircuit, Rvip >= 0.0)
    @variable(microcircuit, Rpc >= 0.0)
    @variable(microcircuit, Rint >= 0.0)

    @variable(microcircuit, sssttodend >= 0.0)
    @variable(microcircuit, spc >= 0.0)
    @variable(microcircuit, ssst >= 0.0)
    @variable(microcircuit, svip >= 0.0)
    @variable(microcircuit, spv >= 0.0)
    @variable(microcircuit, sint >= 0.0)

    @variable(microcircuit, Isoma)
    @variable(microcircuit, Idend )
    @variable(microcircuit, Iexc >= 0.0)
    @variable(microcircuit, Iinh <= 0.0)
    @variable(microcircuit, Ivip)
    @variable(microcircuit, Isst)
    @variable(microcircuit, Ipv)
    @variable(microcircuit,Iint)

    #@variable(microcircuit,sA)
    @variable(microcircuit, u_inttosst)
    ############


    register(microcircuit, :f, 1, f; autodiff = true)
    register(microcircuit, :sigmoid, 1, sigmoid; autodiff = true)
    register(microcircuit, :g, 2, g; autodiff = true)
    ##############
    # Constraints for solving





    #s Constraints
    @NLconstraint(microcircuit, cssst, ssst == τ_GABA*γgaba*Rsst)
    @NLconstraint(microcircuit, csvip, svip == τ_GABA*γgaba*Rvip)
    @NLconstraint(microcircuit, cspv, spv == τ_GABA*γgaba*Rpv)
    @NLconstraint(microcircuit, cssstdend, sssttodend == τ_GABAdend*γgaba*Rsst)

    @NLconstraint(microcircuit, cspc, 0.0 == -spc/τ_NMDA + γnmda*(1.0 - spc)*Rpc)
    @NLconstraint(microcircuit, csint, 0.0 == -sint/τ_NMDA + γnmda*(1.0 - sint)*u_inttosst*2.5*Rint)

    # firing rate constraints
    @NLconstraint(microcircuit, crsst, -Rsst + f(cSST*Isst + r0SST)==0.0)
    @NLconstraint(microcircuit, crpv, -Rpv + f(cPV*Ipv + r0PV)==0.0)
    @NLconstraint(microcircuit, crvip, -Rvip + f(cVIP*Ivip + r0VIP) ==0.0)
    @NLconstraint(microcircuit, csoma, -Rpc + (a*Isoma-b)/(1.0 -exp(-d*(a*Isoma -b))) == 0.0)
    @NLconstraint(microcircuit,cint, α*Rint + (a*Iint-b)/(1.0 -exp(-d*(a*Iint -b))) == 0.0)

    #current Constraints


    @NLconstraint(microcircuit, cIdend, Idend == g(Iexc,Iinh))
    @NLconstraint(microcircuit, cIexcdend, Iexc == Ibg_dend)
    @NLconstraint(microcircuit, cIinhdend,  Iinh == g_ssttodend*sssttodend)
    @NLconstraint(microcircuit, cIsst, Isst == Ibg_sst + g_viptosst*svip + g_inttosst*sint + g_pctosst*spc)
    @NLconstraint(microcircuit, cIpv, Ipv == Ibg_pv + g_pctopv*spc + g_ssttopv*ssst +g_pvtopv*spv)
    @NLconstraint(microcircuit, cIvip, Ivip == I_bg_vip + g_ssttovip*ssst +g_pctovip*spc)
    @NLconstraint(microcircuit, cIsoma, Isoma == Ibg_soma + Istim + Idend + g_pvtopc*spv +g_pctopc*spc )
    @NLconstraint(microcircuit, cIint, Iint == g_pctoint*spc + Ibg_int)


    #adaptation
   # @NLconstraint(microcircuit, csA, sA == Rpc*τA)
    
    # facilitation
    @NLconstraint(microcircuit, cfint, 0 == (U_inttosst - u_inttosst)/τ_facilitation + U_inttosst*(1-u_inttosst)*Rint)
    

    ###########


    optimize!(microcircuit)

return value(Rpc),value(Rint),value(Rvip),value(Rsst),value(Rpv)
end 



function triple_solution()



end




Rpc,Rint,Rsst,Rvip,Rpv = facilitation_circuit_solution(0.1,param)


#TODo some figures now 

stimlist = 0.0:0.05:1.0

using Plots
plotlyjs()

Rpc_list_facilitation = []
Rsst_list_facilitation = []
Rpc_list = []
Rsst_list = []
Rpc_list_adaptation = []
Rsst_list_adaptation = []
param["U_inttosst"] = 0.2

for stim in stimlist
    Rpc,Rint,Rvip,Rsst,Rpv = facilitation_circuit_solution(stim,param)
    push!(Rpc_list_facilitation, Rpc)
    push!(Rsst_list_facilitation, Rsst)

    Rpc,Rint,Rvip,Rsst,Rpv = standard_circuit_solution(stim,param)
    push!(Rpc_list, Rpc)
    push!(Rsst_list, Rsst)

    Rpc,Rint,Rvip,Rsst,Rpv = adaptation_circuit_solution(stim,param)
    push!(Rpc_list_adaptation, Rpc)
    push!(Rsst_list_adaptation, Rsst)

end


scatter(Rpc_list_facilitation,Rsst_list_facilitation,zcolor=stimlist,color=:blues,label="Facilitation")
scatter!(Rpc_list_adaptation,Rsst_list_adaptation,zcolor=stimlist,color=:reds,label="Adaptation")
scatter!(Rpc_list,Rsst_list,zcolor=stimlist,color=:greens,label="Standard")

xlabel!("Soma FR")
ylabel!("SST FR")

scatter(Rpc_list_facilitation - Rpc_list,Rsst_list_facilitation - Rsst_list,zcolor=stimlist,color=:blues,label="Facilitaiton")
scatter!(Rpc_list_adaptation - Rpc_list,Rsst_list_adaptation - Rsst_list,zcolor=stimlist,color=:reds,label="Adaptation")
plot!(legend=false)
xlabel!("Δ FR soma")
ylabel!("Δ FR sst")


function saving_FR_stimlist(stimlist,param;temp = 0.2)

    Rpc_list_facilitation = []
    Rsst_list_facilitation = []
    Rpc_list = []
    Rsst_list = []
    Rpc_list_adaptation = []
    Rsst_list_adaptation = []

    Rvip_list = []
    Rint_list = []
    Rpv_list = []


    Rvip_list_facilitation = []
    Rint_list_facilitation = []
    Rpv_list_facilitation = []

    
    Rvip_list_adaptation = []
    Rint_list_adaptation = []
    Rpv_list_adaptation = []


    param["U_inttosst"] = temp
    
    for stim in stimlist
        Rpc,Rint,Rvip,Rsst,Rpv = facilitation_circuit_solution(stim,param)
        push!(Rpc_list_facilitation, Rpc)
        push!(Rsst_list_facilitation, Rsst)
        push!(Rvip_list_facilitation, Rvip)
        push!(Rpv_list_facilitation, Rpv)
        push!(Rint_list_facilitation, Rint)



        Rpc,Rint,Rvip,Rsst,Rpv = standard_circuit_solution(stim,param)
        push!(Rpc_list, Rpc)
        push!(Rsst_list, Rsst)
        push!(Rvip_list, Rvip)
        push!(Rpv_list, Rpv)
        push!(Rint_list, Rint)

    
        Rpc,Rint,Rvip,Rsst,Rpv = adaptation_circuit_solution(stim,param)
        push!(Rpc_list_adaptation, Rpc)
        push!(Rsst_list_adaptation, Rsst)
        push!(Rvip_list_adaptation, Rvip)
        push!(Rpv_list_adaptation, Rpv)
        push!(Rint_list_adaptation, Rint)



    end

    list_FR = [Rpc_list, Rsst_list,Rvip_list, Rint_list, Rpv_list] 
    list_FR_adaptation = [Rpc_list_adaptation, Rsst_list_adaptation,Rvip_list_adaptation, Rint_list_adaptation, Rpv_list_adaptation] 
    list_FR_facilitation = [Rpc_list_facilitation, Rsst_list_facilitation,Rvip_list_facilitation, Rint_list_facilitation, Rpv_list_facilitation] 

    return list_FR, list_FR_adaptation, list_FR_facilitation

end


list_FR, list_FR_adaptation, list_FR_facilitation = saving_FR_stimlist(stimlist,param)


scatter(list_FR[3], zcolor = stimlist, color=:blues)
scatter!(list_FR_facilitation[3], zcolor = stimlist, color=:greens)
scatter!(list_FR_adaptation[3], zcolor = stimlist, color=:reds)

scatter(list_FR[4], zcolor = stimlist, color=:blues)
scatter!(list_FR_facilitation[4], zcolor = stimlist, color=:greens)
scatter!(list_FR_adaptation[4], zcolor = stimlist, color=:reds)

scatter(list_FR[5], zcolor = stimlist, color=:blues)
scatter!(list_FR_facilitation[5], zcolor = stimlist, color=:greens)
scatter!(list_FR_adaptation[5], zcolor = stimlist, color=:reds)

scatter(list_FR[2], zcolor = stimlist, color=:blues)
scatter!(list_FR_facilitation[2], zcolor = stimlist, color=:greens)
scatter!(list_FR_adaptation[2], zcolor = stimlist, color=:reds)

scatter(list_FR[1], zcolor = stimlist, color=:blues)
scatter!(list_FR_facilitation[1], zcolor = stimlist, color=:greens)
scatter!(list_FR_adaptation[1], zcolor = stimlist, color=:reds)


function plot_Δ_stim(list_FR,list_FR_adaptation,list_FR_facilitation,stimlist)

    fig = plot(layout = (2,3))

    scatter!(fig[1],stimlist,list_FR_facilitation[1] - list_FR[1], zcolor = stimlist, color=:greens)
    scatter!(fig[1],stimlist,list_FR_adaptation[1] - list_FR[1], zcolor = stimlist, color=:reds)

    scatter!(fig[2],stimlist,list_FR_facilitation[2] - list_FR[2], zcolor = stimlist, color=:greens)
    scatter!(fig[2],stimlist,list_FR_adaptation[2] - list_FR[2], zcolor = stimlist, color=:reds)


    scatter!(fig[3],stimlist,list_FR_facilitation[3] - list_FR[3], zcolor = stimlist, color=:greens)
    scatter!(fig[3],stimlist,list_FR_adaptation[3] - list_FR[3], zcolor = stimlist, color=:reds)


    scatter!(fig[4],stimlist,list_FR_facilitation[4] - list_FR[4], zcolor = stimlist, color=:greens)
    scatter!(fig[4],stimlist,list_FR_adaptation[4] - list_FR[4], zcolor = stimlist, color=:reds)
    plot!(fig, legend = false)
    title!(fig[1], "Δ PC")
    title!(fig[2], "Δ SST")
    title!(fig[3], "Δ VIP")

    title!(fig[4], "Δ Integrator")
    scatter!(fig[5],stimlist,list_FR_facilitation[5] - list_FR[5], zcolor = stimlist, color=:greens)
    scatter!(fig[5],stimlist,list_FR_adaptation[5] - list_FR[5], zcolor = stimlist, color=:reds)
    ylabel!("Δ firing rate")
    xlabel!("Stimulus to soma")
    title!(fig[5], "Δ PV")
   # plot!(colorbar=false)
    return fig
end



list_FR, list_FR_adaptation, list_FR_facilitation = saving_FR_stimlist(stimlist,param;temp = 0.2)


figΔ = plot_Δ_stim(list_FR,list_FR_adaptation,list_FR_facilitation,stimlist)
savefig(figΔ, plotsdir("4-Theoretical-adaptation","figΔ-firing-rate-stimlist.png"))




fig
figΔ


###############
# Variation of time constants
##################


time_constant_list = 0.1:0.1:2.0

Rpc_list_facilitation = []
Rsst_list_facilitation = []

Rpc_list_adaptation = []
Rsst_list_adaptation = []
stim = 0.4
for tc in time_constant_list
    param["τA"] = tc
    param["τ_facilitation"] = tc
    Rpc,Rint,Rvip,Rsst,Rpv = facilitation_circuit_solution(stim,param)
    push!(Rpc_list_facilitation, Rpc)
    push!(Rsst_list_facilitation, Rsst)

    Rpc,Rint,Rvip,Rsst,Rpv = adaptation_circuit_solution(stim,param)
    push!(Rpc_list_adaptation, Rpc)
    push!(Rsst_list_adaptation, Rsst)

end


function saving_FR_timeconstant(timelist,param;temp = 0.2,stim = 0.4)

    Rpc_list_facilitation = []
    Rsst_list_facilitation = []
    Rpc_list = []
    Rsst_list = []
    Rpc_list_adaptation = []
    Rsst_list_adaptation = []

    Rvip_list = []
    Rint_list = []
    Rpv_list = []


    Rvip_list_facilitation = []
    Rint_list_facilitation = []
    Rpv_list_facilitation = []

    
    Rvip_list_adaptation = []
    Rint_list_adaptation = []
    Rpv_list_adaptation = []


    param["U_inttosst"] = temp
    
    for tc in time_constant_list
        param["τA"] = tc
        param["τ_facilitation"] = tc
        Rpc,Rint,Rvip,Rsst,Rpv = facilitation_circuit_solution(stim,param)
        push!(Rpc_list_facilitation, Rpc)
        push!(Rsst_list_facilitation, Rsst)
        push!(Rvip_list_facilitation, Rvip)
        push!(Rpv_list_facilitation, Rpv)
        push!(Rint_list_facilitation, Rint)



        Rpc,Rint,Rvip,Rsst,Rpv = standard_circuit_solution(stim,param)
        push!(Rpc_list, Rpc)
        push!(Rsst_list, Rsst)
        push!(Rvip_list, Rvip)
        push!(Rpv_list, Rpv)
        push!(Rint_list, Rint)

    
        Rpc,Rint,Rvip,Rsst,Rpv = adaptation_circuit_solution(stim,param)
        push!(Rpc_list_adaptation, Rpc)
        push!(Rsst_list_adaptation, Rsst)
        push!(Rvip_list_adaptation, Rvip)
        push!(Rpv_list_adaptation, Rpv)
        push!(Rint_list_adaptation, Rint)



    end

    list_FR = [Rpc_list, Rsst_list,Rvip_list, Rint_list, Rpv_list] 
    list_FR_adaptation = [Rpc_list_adaptation, Rsst_list_adaptation,Rvip_list_adaptation, Rint_list_adaptation, Rpv_list_adaptation] 
    list_FR_facilitation = [Rpc_list_facilitation, Rsst_list_facilitation,Rvip_list_facilitation, Rint_list_facilitation, Rpv_list_facilitation] 

    return list_FR, list_FR_adaptation, list_FR_facilitation

end

list_FR_TC, list_FR_TC_A, list_FR_TC_F = saving_FR_timeconstant(time_constant_list,param)


scatter(list_FR_TC[1])
scatter!(list_FR_TC_A[1])
scatter!(list_FR_TC_F[1])


scatter(list_FR_TC[2])
scatter!(list_FR_TC_A[2])
scatter!(list_FR_TC_F[2])


scatter(list_FR_TC[3])
scatter!(list_FR_TC_A[3])
scatter!(list_FR_TC_F[3])


scatter(list_FR_TC[4])
scatter!(list_FR_TC_A[4])
scatter!(list_FR_TC_F[4])



scatter(list_FR_TC[5])
scatter!(list_FR_TC_A[5])
scatter!(list_FR_TC_F[5])


scatter(Rpc_list_facilitation,Rsst_list_facilitation,zcolor=time_constant_list,color=:blues,label="Facilitation")
scatter!(Rpc_list_adaptation,Rsst_list_adaptation,zcolor=time_constant_list,color=:reds,label="Adaptation")
xlabel!("Soma FR")
ylabel!("SST FR")
plot!(legend=false)