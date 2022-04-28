import os
import sys
import importlib

#project_dir = '/mnt/c/Users/kevin/OneDrive/2-code/1-Research_projects/1-project-categorical-MMN'
project_dir = '/home/kb3856/1-project-categorical-MMN'
src_dir = os.path.join(project_dir, 'src')
sys.path.append(src_dir)
import parameters as params
import classes as cl
import basic_functions as bf
import numpy as np

import pickle
import json

# data_dir = os.path.join(project_dir, 'data/sims/script2')
data_dir = os.path.join(project_dir, 'datatemp')

# This script will generate the data to try to understand the behavior of the adaptation and prediction index
# Two parameters will be varied: adaptation strength and feedback strength
# 
#


#define the parameters being varied
# two delay in the MMN are tested 50ms and 100ms
delay = 100
adaptation_list = [-0.001, -0.005, -0.01, -0.05]
top_down_feedback_list = [1.0, 2.0, 3.0, 4.0]
Tinter_list = [0.2, 0.5,0.7,1.0,1.5,2.0]
dev_id_list = [28]#[88,78,68,58,48,38,28]


params.PARAMS_Stimulus['prob_std'] = 0.8
params.PARAMS_Simulation['t_total'] = 120.0#40.0 #in seconds, constant as the task is probabilistic
# params.PARAMS_Simulation['N_t'] = int(params.PARAMS_Simulation['t_total']/self.PARAMS_Simulation['dt'])
params.PARAMS_Stimulus['list_std'] = [25, 50, 75, 100, 125]
#saving variables
saving_multi_control = np.zeros((len(adaptation_list), len(top_down_feedback_list), len(Tinter_list)))
saving_proba_control = np.zeros((len(adaptation_list), len(top_down_feedback_list), len(Tinter_list)))
saving_proba_dev = np.zeros((len(adaptation_list), len(top_down_feedback_list), len(Tinter_list)))


for count_d,d in enumerate(dev_id_list):
    
    params.PARAMS_Stimulus['dev_id'] = d
    for count_type,type_stim in enumerate(['multi_control_paradigm', 'probabilistic_MMN']):
        
        params.PARAMS_Stimulus['type'] = type_stim
        
        for count_a,a in enumerate(adaptation_list):
            
            params.PARAMS_PC['gA'] = a #values are coherent with other papers
        
        
            for count_topdown,t in enumerate(top_down_feedback_list):
            
                params.PARAMS_Synapses_Integrator['wmax'] = t*0.15/20 #it corresponds to a multiplier
            
                for count_tinter,i in enumerate(Tinter_list):
                    params.PARAMS_Stimulus['Tinter'] = i

                

                    stim = cl.Stimulus(params.PARAMS_Stimulus, params.PARAMS_Simulation, params.PARAMS_PC['Ncells'])

                    my_network = cl.Network(params.PARAMS_Integrator, params.PARAMS_Synapses_Integrator, params.PARAMS_PC, params.PARAMS_PV, params.PARAMS_SST, params.PARAMS_VIP, params.PARAMS_NDF, params.PARAMS_Simulation, params.PARAMS_Stimulus)

                    my_network.full_dynamics(stim.stimulus)


                    if type_stim == 'multi_control_paradigm':
                            
                            mean_proba_control = my_network.compute_mean_firing_rate()
                            saving_multi_control[count_a, count_topdown, count_tinter] = mean_proba_control
                            
                    elif type_stim == 'probabilistic_MMN':  
                            mean_std, mean_dev = my_network.compute_mean_firing_rate(stim)
                            saving_proba_control[count_a, count_topdown, count_tinter] = mean_std
                            
                            saving_proba_dev[count_a, count_topdown, count_tinter] = mean_dev
                            
                            
    save_dict = {
        'PARAMS': my_network.PARAMS_ALL,
        'multi_control': saving_multi_control,
        'proba_control': saving_proba_control,
        'proba_dev': saving_proba_dev,
        'dev_id': d
        
    }
    file_dir = os.path.join(data_dir, 'script2_delay_'+str(delay)+'mean_firing_rates_dev_id'+str(d)+'.json')
    
    with open(file_dir, 'w') as handle:
                json.dump(save_dict, handle, cls = cl.NumpyEncoder)
                    

                              