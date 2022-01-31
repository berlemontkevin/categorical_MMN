#This file will define the parameters of 
# the simulations

import numpy as np

PARAMS_Integrator = {
    'Ncells': np.int32(128),
    'alpha': 0.5,
    'a': 135.0,
    'b': 54.0,
    'd': 0.308, #seconds
    'Ibg': 310 * 0.001, #nA 310
    'sigma_con': 43.2,
    'Jmin': -0.5,
    'Jmax':1.2,#1.43,
    'tau':0.2, #seconds
    'tau_NMDA': 60 * 0.001,
    'gamma_NMDA': 0.641 ,
    # 'weight_to_dend': 0.4,
    # 'weight_to_pv': 0.056,
    # 'weight_to_vip': 0.045,
    # 'weight_to_ndf': 10.8,#0.025,
    # 'sigma_to_dend': 43.2,
    # 'sigma_to_ndf': 1.0,#43.2,
    # 'sigma_to_pv': 43.2,
    # 'sigma_to_vip': 43.2,
    'depression_to_vip': True,
    # 'depression_to_ndf': True,
    'depression_to_dend': True,
    'depression_to_pv': True,
    'fD': 0.2,
    'tau_depression': 0.2,
    'mult_nmda_depression': 2.0,
    'dt': 0.001
}

PARAMS_Synapses_Integrator = {
    'tau_nmda': 60 * 0.001,
    'gamma_nmda': 0.641 ,
    'weight_to_dend': 0.7,#0.5
    'weight_to_pv': 0.066,
    'weight_to_vip': 0.035,
    'weight_to_ndf': 0.155/30,
    'sigma_to_dend': 43.2,
    'sigma_to_ndf': 1.2,
    'sigma_to_pv': 43.2,
    'sigma_to_vip': 43.2,
    'dt': 0.001,
    'gamma': 0.01,
    'lambda_dec': 0.2,
    'wmax': 0.15/15
    }

PARAMS_PC= {
    'Ncells': np.int32(128),
    'a': 135.0,
    'b': 54.0,
    'd': 0.308, #seconds
    'Ibg_dend': 100 * 0.001, #nA
    'tau_dend': 0.002, #seconds
    'Ibg_soma': 0.15, #nA
    'tau_soma': 0.002, #seconds
    'weight_to_soma': 0.10,#0.18,
    'weight_to_vip': 0.158,
    'weight_to_sst': 0.8835,
    'weight_to_pv':0.18,
    'weight_to_integrator': 0.40,#0.15
    'sigma_to_sst': 360.0,#360.0,#43.2,
    'sigma_to_pv': 43.2,
    'sigma_to_integrator': 1.0,#43.2,
    'sigma_to_pc': 43.2,
    'facilitation_to_vip': True,
    'facilitation_to_sst': True,
    'depression_to_pv': True,
    'fD': 0.2,
    'tau_depression': 0.2,
    'mult_nmda_depression': 2.0,
    'U_facilitation': 0.1,
    'tau_facilitation': 1.5,
    'mult_nmda_facilitation': 2.0,
    'tau_NMDA': 60 * 0.001,
    'gamma_NMDA': 0.641 ,
    'tau_adaptation': 0.1,
    'gA': -0.01,
    'dt': 0.001
}

PARAMS_PV= {
    'Ncells': np.int32(128),
    'c_I': 330.0,
    'r0':-95.0,
    'Ibg': 290 * 0.001, #nA 250
    'tau':0.002, #seconds
    'weight_to_soma': -1.3,
    'weight_to_pv': -0.18,
    'depression_to_pv': True,
    'fD': 0.2,
    'tau_depression': 0.2,
    'mult_gaba_depression': 2.0,
    'U_facilitation': 0.1,
    'tau_facilitation': 1.5,
    'mult_gaba_facilitation': 2.0,
    'tau_GABA': 0.005,
    'gamma_GABA': 0.2,
    'dt': 0.001
}

PARAMS_SST= {
    'Ncells': np.int32(128),
    'c_I': 132.0,
    'r0':-33.0,
    'Ibg': 290 * 0.001, #nA
    'tau':0.002, #seconds
    'weight_to_vip': -0.1,
    'weight_to_dend': -0.19,
    'weight_to_pv': -0.17,
    'weight_to_ndf': -0.4,#-0.4,
    'facilitation_to_vip': True,
    'fD': 0.2,
    'tau_depression': 0.2,
    'mult_gaba_depression': 2.0,
    'U_facilitation': 0.1,
    'tau_facilitation': 1.5,
    'mult_gaba_facilitation': 2.0,
    'tau_GABA': 0.005,
    'gamma_GABA': 0.2,
    'tau_adaptation': 0.1,
    'gA': -0.01,
    'dt': 0.001
}

PARAMS_VIP= {
    'Ncells': np.int32(128),
    'c_I': 132.0,
    'r0':-33.0,
    'Ibg': 290 * 0.001, #nA
    'tau':0.002, #seconds
    'weight_to_sst': -0.17,
    'facilitation_to_sst': True,
    'fD': 0.2,
    'tau_depression': 0.2,
    'mult_gaba_depression': 2.0,
    'U_facilitation': 0.1,
    'tau_facilitation': 1.5,
    'mult_gaba_facilitation': 2.0,
    'tau_GABA': 0.005,
    'gamma_GABA': 0.2,
    'dt': 0.001
}

PARAMS_NDF= {
    'Ncells': np.int32(128),
    'c_I': 132.0,
    'r0':-33.0,
    'Ibg': 290 * 0.001, #nA
    'tau':0.002, #seconds
    'weight_to_dend': -0.5*0.5/30,#0.89,
    'weight_to_pv': -0.0,#-0.47,
    'depression_to_dend': True,
    'fD': 0.2,
    'tau_depression': 0.2,
    'mult_gaba_depression': 2.0,
    'U_facilitation': 0.1,
    'tau_facilitation': 1.5,
    'mult_gaba_facilitation': 2.0,
    'tau_GABA': 0.005,
    'gamma_GABA': 0.2,
    'dt': 0.001,
    'wmax': 0.89/15,#/20
    'gamma': 0.007,#0.005 ,
    'lambda_dec': 0.2
}

PARAMS_Synapses= {
    'cross_integrator_to_vip_depression': True,
    'cross_integrator_to_pv_depression': True,
    'cross_integrator_to_sst_facilitation': True,
    'cross_integrator_to_dend_depression': True,
    'cross_soma_to_sst_facilitation': True
}

PARAMS_Simulation = {
    'dt': 0.001,
    't_start': 0.0,
    't_total': 8.0
}

PARAMS_ALL = {
    'Synapses': PARAMS_Synapses,
    'Simulation': PARAMS_Simulation,
    'PARAMS_INT': PARAMS_Integrator,
    'PARAMS_PC': PARAMS_PC,
    'PARAMS_PV': PARAMS_PV,
    'PARAMS_SST': PARAMS_SST,
    'PARAMS_VIP': PARAMS_VIP,
    'PARAMS_NDF': PARAMS_NDF
}

PARAMS_Stimulus = {
    'strength': 0.2,
    'Tinter': 0.5,
    'Tstim': 0.2,
    'nbr_rep': 3,
    'std_id': 20,
    'dev_id': 88
}