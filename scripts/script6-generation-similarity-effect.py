import os
import sys
import importlib

project_dir = '/home/kb3856/1-project-categorical-MMN'
src_dir = os.path.join(project_dir, 'src')
sys.path.append(src_dir)
import parameters as params
import classes as cl
import basic_functions as bf
import numpy as np

import pickle
import json

from itertools import product
from multiprocessing import Pool

data_dir = os.path.join(project_dir, 'datatemp')


def create_list_dict(d):
    ''' This funciton returns all the possible combinations of elements within the keys
    '''
    list_dict = []
    key_list = list(d.keys())
    length_list = [len(d[key]) for key in key_list]
    nbr_list = len(key_list)
    
    all_lists = list(product(*[d[key] for key in key_list]))
    return all_lists
   
dict_param = {
    'adaptation_list' : [-0.01,-0.03,-0.05,0.0],#[-0.01,-0.02,-0.03,-0.04,-0.05,0.0],
'top_down_feedback_list' : [0.0,0.5,1.0,1.2,1.4,1.6,1.8,2.0],
'Tinter_list' : [0.5,0.7,1.0,1.5,2.0,3.0,5.0],
'dev_id_list' : [88,78,68,58,48,38],
'delay_list': [50.0,100.0],
'ndf_plasticity': [1],
'int_plasticity': [1],
'bool_sigma_ndf': [1],
'bool_sigma_int': [1]    
}
    
all_dicts = create_list_dict(dict_param)


def run_parallel(list_param):
    params_stim = params.PARAMS_Stimulus.copy()
    params_sim = params.PARAMS_Simulation.copy()
    params_int = params.PARAMS_Integrator.copy()
    params_sst = params.PARAMS_SST.copy()
    params_ndf = params.PARAMS_NDF.copy()
    params_pc = params.PARAMS_PC.copy()
    params_vip = params.PARAMS_VIP.copy()
    params_pv = params.PARAMS_PV.copy()
    params_syn = params.PARAMS_Synapses_Integrator.copy()

    params_stim['prob_std'] = 0.8
    params_sim['t_total'] = 200.0#40.0 #in seconds, constant as the task is probabilistic
    # params.PARAMS_Simulation['N_t'] = int(params.PARAMS_Simulation['t_total']/self.PARAMS_Simulation['dt'])
    params_stim['list_std'] = [25, 50, 75, 100, 125]
    delay = list_param[4]*0.001

    for count_type,type_stim in enumerate(['multi_control_paradigm', 'probabilistic_MMN']):
        
        params_stim['type'] = type_stim
        params_pc['gA'] = list_param[0] #values are coherent with other papers
        params_sst['gA'] = list_param[0] #values are coherent with other papers

        params_stim['dev_id'] = list_param[3]

        if list_param[7] == 1:
            params_ndf['sigma_to_dend'] = 43.2
            params_ndf['weight_to_dend'] =  -2.0*0.89/15*20,#/20']

        else:
            params_ndf['sigma_to_dend'] = 1.2
            params_ndf['weight_to_dend'] = -2.0*0.89/15,#/20']

        if list_param[8] == 1:
            params_syn['wmax'] = list_param[1]*0.15/20 #it corresponds to a multiplier
            params_syn['weight_to_ndf'] = list_param[1]*0.15/20*20 #it corresponds to a multiplier
            params_syn['sigma_to_ndf'] = 43.2

        else:
            params_syn['wmax'] = list_param[1]*0.15/20 #it corresponds to a multiplier
            params_syn['weight_to_ndf'] = list_param[1]*0.15/20 #it corresponds to a multiplier
            params_syn['sigma_to_bdf'] = 1.2  
     
        params_stim['Tinter'] = list_param[2]
        
        params_syn['bool_plasticity'] = list_param[6]
        params_ndf['bool_plasticity'] = list_param[5]
        
        
        stim = cl.Stimulus(params_stim, params_sim, params_pc['Ncells'])

        my_network = cl.Network(params_int, params_syn, params_pc, params_pv, params_sst, params_vip, params_ndf, params_sim, params_stim)

        my_network.full_dynamics(stim.stimulus)


        if type_stim == 'multi_control_paradigm':
                            
                            mean_proba_control, sum_control = my_network.compute_mean_firing_rate(delay=delay)
                            saving_multi_control = mean_proba_control
                            
        elif type_stim == 'probabilistic_MMN':  
            
                            mean_std, mean_dev, sum_std, sum_dev = my_network.compute_mean_firing_rate(stim, delay = delay)
                            saving_proba_control = mean_std
                            
                            saving_proba_dev = mean_dev
        
    
    save_dict = {
        'PARAMS': my_network.PARAMS_ALL,
        'multi_control': saving_multi_control,
        'proba_control': saving_proba_control,
        'proba_dev': saving_proba_dev,
        'sum_control': sum_control,
        'sum_dev': sum_dev,
        'sum_std': sum_std        
    }
    name = 'script6_mean_firing_rates_plast_int'+str(list_param[6]) + '_plast_ndf'+ str(list_param[5])+'dev_id'+str(list_param[3])+'adapt_'+str(list_param[0])+'topdown_'+str(list_param[1])+'Tinter_'+str(list_param[2])+'delay_'+str(list_param[4])+'broadtodend_'+str(list_param[7])+'broadtondf_'+str(list_param[8])+'.json'
    file_dir = os.path.join(data_dir, name)
    
    with open(file_dir, 'w') as handle:
                json.dump(save_dict, handle, cls = cl.NumpyEncoder)

    return #print(name)




n_threads = 30
pool = Pool(processes = n_threads)
pool.map(run_parallel, all_dicts)