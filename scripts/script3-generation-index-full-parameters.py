import os
import sys
import importlib

project_dir = '/mnt/c/Users/kevin/OneDrive/2-code/1-Research_projects/1-project-categorical-MMN'
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

data_dir = os.path.join(project_dir, 'data/sims/script3')


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
    'adaptation_list' : [-0.01,-0.02,-0.03,-0.04,-0.04,0.0],
'top_down_feedback_list' : [0.0,0.5,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0],
'Tinter_list' : [0.2, 0.5,0.7,1.0,2.0],
'dev_id_list' : [88,68,58,38,28],
'delay_list': [50.0,100.0]    
}
    
all_dicts = create_list_dict(dict_param)


def run_parallel(list_param):
    params_temp = params.PARAMS_Stimulus.copy()
    params_temp['prob_std'] = 0.8
    params_temp['t_total'] = 0.01#40.0 #in seconds, constant as the task is probabilistic
    # params.PARAMS_Simulation['N_t'] = int(params.PARAMS_Simulation['t_total']/self.PARAMS_Simulation['dt'])
    params_temp['list_std'] = [25, 50, 75, 100, 125]
    delay = list_param[4]

    for count_type,type_stim in enumerate(['multi_control_paradigm', 'probabilistic_MMN']):
        
        params_temp['type'] = type_stim
        params_temp['gA'] = list_param[0] #values are coherent with other papers
        params_temp['dev_id'] = list_param[3]
        params_temp['wmax'] = list_param[1]*0.15/20 #it corresponds to a multiplier
        params_temp['Tinter'] = list_param[2]
        
        stim = cl.Stimulus(params.PARAMS_Stimulus, params.PARAMS_Simulation, params.PARAMS_PC['Ncells'])

        my_network = cl.Network(params_temp.PARAMS_Integrator, params_temp.PARAMS_Synapses_Integrator, params_temp.PARAMS_PC, params_temp.PARAMS_PV, params_temp.PARAMS_SST, params_temp.PARAMS_VIP, params_temp.PARAMS_NDF, params_temp.PARAMS_Simulation, params_temp.PARAMS_Stimulus)

        my_network.full_dynamics(stim.stimulus)


        if type_stim == 'multi_control_paradigm':
                            
                            mean_proba_control = my_network.compute_mean_firing_rate(delay=delay)
                            saving_multi_control = mean_proba_control
                            
        elif type_stim == 'probabilistic_MMN':  
            
                            mean_std, mean_dev = my_network.compute_mean_firing_rate(stim, delay = delay)
                            saving_proba_control = mean_std
                            
                            saving_proba_dev = mean_dev
        
    
    save_dict = {
        'PARAMS': my_network.PARAMS_ALL,
        'multi_control': saving_multi_control,
        'proba_control': saving_proba_control,
        'proba_dev': saving_proba_dev        
    }
    name = 'script2_mean_firing_rates_dev_id'+str(list_param[3])+'adapt_'+str(list_param[0])+'topdown_'+str(list_param[1])+'Tinter_'+str(list_param[2])+'.json'
    file_dir = os.path.join(data_dir, name)
    
    with open(file_dir, 'w') as handle:
                json.dump(save_dict, handle, cls = cl.NumpyEncoder)

    return print(name)




    pool.map(run_parallel, all_dicts)