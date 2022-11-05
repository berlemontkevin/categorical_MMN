#This file will define the parameters of 
# the simulations

import numpy as np

PARAMS_Integrator = {
    'Ncells': 128,
    'alpha': 0.6,
    'a': 135.0,
    'b': 54.0,
    'd': 0.308, #seconds
    'Ibg': 310 * 0.001, #nA 310
    'sigma_con': 43.2,
    'Jmin': -0.5,#0.5
    'Jmax':1.2,#1.43,
    'tau':0.4, #seconds
    'tau_NMDA': 60 * 0.001,
    'gamma_NMDA': 0.641 ,
    'depression_to_vip': True,
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
    'weight_to_ndf': 0.155,#/50,
    'sigma_to_dend': 43.2,
    'sigma_to_ndf': 43.2,#1.2,
    'sigma_to_pv': 43.2,
    'sigma_to_vip': 43.2,
    'dt': 0.001,
    'gamma': 0.01,
    'lambda_dec': 10.0,#0.2,
    'wmax': 2.0*0.15/20,
    'bool_plasticity': 1
    }

PARAMS_PC= {
    'Ncells': 128,
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
    'weight_to_integrator': 0.50,#0.15
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
    'Ncells': 128,
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
    'Ncells': 128,
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
    'dt': 0.001,
    'wmax': 0.19/15,#/20
    'gamma': 0.007,#0.005 ,
    'lambda_dec': 0.2
}

PARAMS_VIP= {
    'Ncells': 128,
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
    'Ncells': 128,
    'c_I': 132.0,
    'r0':-33.0,
    'Ibg': 290 * 0.001, #nA
    'tau':0.002, #seconds
    'weight_to_dend': -0.5*0.5,#/30,#0.89,
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
    'wmax': 2.0*0.89/15,#/20
    'gamma': 0.007,#0.005 ,
    'lambda_dec': 1.0,#0.2
    'bool_plasticity': 1,
    'sigma_to_dend': 43.2#1.2,
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
    't_total': 16.0
}



PARAMS_Stimulus = {
    'type': 'deterministic_MMN',
    'strength_std': 0.2,
    'strength_dev': 0.2,
    'sigma_std': 43.2,
    'sigma_dev': 43.2,
    'Tinter': 0.5,
    'Tstim': 0.2,
    'nbr_rep_std': 6,
    'nbr_rep_dev': 2,
    'std_id': 20, # the idx of the neuron receiving the stimulus
    'dev_id': 88,
    'Tresting': 4.0 #initial resting time in seconds
}

PARAMS_ALL = {
    'Synapses': PARAMS_Synapses,
    'Simulation': PARAMS_Simulation,
    'PARAMS_INT': PARAMS_Integrator,
    'PARAMS_PC': PARAMS_PC,
    'PARAMS_PV': PARAMS_PV,
    'PARAMS_SST': PARAMS_SST,
    'PARAMS_VIP': PARAMS_VIP,
    'PARAMS_NDF': PARAMS_NDF,
    'Stimulus': PARAMS_Stimulus
}