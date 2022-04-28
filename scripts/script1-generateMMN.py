import os
import sys
import importlib

project_dir = '/mnt/c/Users/kevin/OneDrive/2-code/1-Research_projects/1-project-categorical-MMN'
src_dir = os.path.join(project_dir, 'src')
sys.path.append(src_dir)
import parameters as params
import classes as cl
import basic_functions as bf


data_dir = os.path.join(project_dir, 'data/sims/script1')

# define the list of parameters to be varied
Tinter_list = [0.2, 0.5,0.7,1.0,1.5,2.0]
nbr_rep_list = [1,2,3,4,5,6,7,8,9,10]
dev_id_list = [88,78,68,58,48,38,28]

for i in Tinter_list:

    for n in nbr_rep_list:
        for d in dev_id_list:
        
                 
                params.PARAMS_Stimulus['Tinter'] = i
                params.PARAMS_Stimulus['nbr_rep_std'] = n
                params.PARAMS_Stimulus['dev_id'] = d

                params.PARAMS_Simulation['t_total'] = (params.PARAMS_Stimulus['Tinter'] + params.PARAMS_Stimulus['Tstim'])*(params.PARAMS_Stimulus['nbr_rep_std']+ params.PARAMS_Stimulus['nbr_rep_dev']+1) + params.PARAMS_Stimulus['Tresting']

                stim = cl.Stimulus(params.PARAMS_Stimulus, params.PARAMS_Simulation, params.PARAMS_PC['Ncells'])

                my_network = cl.Network(params.PARAMS_Integrator, params.PARAMS_Synapses_Integrator, params.PARAMS_PC, params.PARAMS_PV, params.PARAMS_SST, params.PARAMS_VIP, params.PARAMS_NDF, params.PARAMS_Simulation, params.PARAMS_Stimulus)

                my_network.full_dynamics(stim.stimulus)

                file_dir = os.path.join(data_dir, 'script1_Tinter' +str(i) + '_nbr_rep_'+str(n)+'_dev_id_'+str(d)+'.json')
                my_network.write_firing_rates_in_file(file_dir)


