# These functions will be used to calculate the dynamics of the system

import numpy as np
import basic_functions as bf
import numba as nb
import initialization as init
import parameters as par
# defines the structures of the network

#define the rates
rate_integrator = np.ones(par.PARAMS_Integrator['Ncells'])
rate_vip = np.ones(par.PARAMS_Integrator['Ncells'])
rate_sst = np.ones(par.PARAMS_Integrator['Ncells'])
rate_pv = np.ones(par.PARAMS_Integrator['Ncells'])
rate_pc = np.ones(par.PARAMS_Integrator['Ncells'])
rate_ndf = np.ones(par.PARAMS_Integrator['Ncells'])
VARIABLES_RATES = {
    'rate_integrator': rate_integrator,
    'rate_vip': rate_vip,
    'rate_sst': rate_sst,
    'rate_pv': rate_pv,
    'rate_pc': rate_pc,
    'rate_ndf': rate_ndf
}

#define the currents
current_integrator = np.ones(par.PARAMS_Integrator['Ncells'])
current_pc = np.ones(par.PARAMS_Integrator['Ncells'])
current_sst = np.ones(par.PARAMS_Integrator['Ncells'])
current_pv = np.ones(par.PARAMS_Integrator['Ncells'])
current_vip = np.ones(par.PARAMS_Integrator['Ncells'])
current_ndf = np.ones(par.PARAMS_Integrator['Ncells'])
current_dend_Iexc = np.ones(par.PARAMS_Integrator['Ncells'])
current_dend_Iinh = np.ones(par.PARAMS_Integrator['Ncells'])
VARIABLES_CURRENTS = {
    'current_integrator': current_integrator,
    'current_pc': current_pc,
    'current_sst': current_sst,
    'current_pv': current_pv,
    'current_vip': current_vip,
    'current_ndf': current_ndf,
    'current_dend_Iexc': current_dend_Iexc,
    'current_dend_Iinh': current_dend_Iinh
}


# synaptic variables
#TODO a better definition of the variable that's modulated
VARIABLES_SYNAPTIC = {
    'W_integrator_to_dend': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_integrator_to_pv': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_integrator_to_ndf': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_integrator_to_vip': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    # weights for pc
    'W_pc_to_pc': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_pc_to_pv': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_pc_to_sst': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_pc_to_vip': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_pc_to_integrator': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    # weights for pv
    'W_pv_to_pc': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_pv_to_pv': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_pv_to_sst': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    # weights for sst
    'W_sst_to_dend': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_sst_to_vip': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_sst_to_pv': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    'W_sst_to_ndf': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    # weights for vip
    'W_vip_to_sst': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells'])),
    # weights for ndf
    'W_ndf_to_dend': 0.1*np.ones((par.PARAMS_Integrator['Ncells'],par.PARAMS_Integrator['Ncells']))
}


@njit(cache=true)
def update_dend_to_soma(
    VARIABLES_RATES,
    VARIABLES_CURRENTS,
    PARAMS_PC,
    PARAMS_Simulation
    ):

    ''' This funciton will update the current that is send by the dendrites to the pc
    '''
    for i in range(PARAMS_PC['Ncells']):
        VARIABLES_CURRENTS['current_pc'][i] = VARIABLES_CURRENTS['current_pc'][i] + \
            bf.current_from_dend(VARIABLES_CURRENTS['current_dend_Iinh'][i],VARIABLES_CURRENTS['current_dend_Iexc'][i],PARAMS_PC['c1'],PARAMS_PC['c2'],PARAMS_PC['c3'], PARAMS_PC['c4'], PARAMS_PC['c5'], PARAMS_PC['c6'])
            
    return None

#region update all the rates
@njit(cache=true)
def update_rate_pc(VARIABLES_RATES,VARIABLES_CURRENTS,PARAMS_PC,PARAMS_Simulation):
    ''' This function will update the rate of the pc
    '''
    for i in range(PARAMS_PC['Ncells']):
        
        fI = bf.abott_fI_curve(VARIABLES_CURRENTS['current_pc'][i], PARAMS_PC['a'], PARAMS_PC['b'], PARAMS_PC['d'])
        fI = -VARIABLES_RATES['rate_pc'][i] + fI
        
        VARIABLES_RATES['rate_pc'][i] = bf.euler_update(VARIABLES_RATES['rate_pc'][i],fI, PARAMS_Simulation['dt'], tau =  PARAMS_PC['tau_soma'])
    
    return None

@njit(cache=true)
def update_rate_vip(VARIABLES_RATES,VARIABLES_CURRENTS,PARAMS_VIP,PARAMS_Simulation):
    ''' This function will update the rate of the vip
    '''
    for i in range(PARAMS_VIP['Ncells']):
        
        fI = bf.relu(VARIABLES_CURRENTS['current_vip'][i]*PARAMS_VIP['cI'] + PARAMS_VIP['r0'])
        fI = -VARIABLES_RATES['rate_vip'][i] + fI
        
        VARIABLES_RATES['rate_vip'][i] = bf.euler_update(VARIABLES_RATES['rate_vip'][i],fI, PARAMS_Simulation['dt'], tau = PARAMS_VIP['tau_soma'])
    
    return None


@njit(cache=true)
def update_rate_sst(VARIABLES_RATES,VARIABLES_CURRENTS,PARAMS_SST,PARAMS_Simulation):
    ''' This function will update the rate of the sst
    '''
    for i in range(PARAMS_SST['Ncells']):
        
        fI = bf.relu(VARIABLES_CURRENTS['current_sst'][i]*PARAMS_SST['cI'] + PARAMS_SST['r0'])
        fI = -VARIABLES_RATES['rate_sst'][i] + fI
        
        VARIABLES_RATES['rate_sst'][i] = bf.euler_update(VARIABLES_RATES['rate_sst'][i],fI, PARAMS_Simulation['dt'], tau = PARAMS_SST['tau_soma'])
    
    return None

@njit(cache=true)
def update_rate_ndf(VARIABLES_RATES,VARIABLES_CURRENTS,PARAMS_NDF,PARAMS_Simulation):
    ''' This function will update the rate of the ndf
    '''
    for i in range(PARAMS_NDF['Ncells']):
        
        fI = bf.relu(VARIABLES_CURRENTS['current_ndf'][i]*PARAMS_NDF['cI'] + PARAMS_NDF['r0'])
        fI = -VARIABLES_RATES['rate_ndf'][i] + fI
        
        VARIABLES_RATES['rate_ndf'][i] = bf.euler_update(VARIABLES_RATES['rate_ndf'][i],fI, PARAMS_Simulation['dt'], tau = PARAMS_NDF['tau_soma'])
    
    return None

@njit(cache=true)
def update_rate_pv(VARIABLES_RATES,VARIABLES_CURRENTS,PARAMS_PV,PARAMS_Simulation):
    ''' This function will update the rate of the pv
    '''
    for i in range(PARAMS_PV['Ncells']):
        
        fI = bf.relu(VARIABLES_CURRENTS['current_pv'][i]*PARAMS_PV['cI'] + PARAMS_PV['r0'])
        fI = -VARIABLES_RATES['rate_pv'][i] + fI
        
        VARIABLES_RATES['rate_pv'][i] = bf.euler_update(VARIABLES_RATES['rate_pv'][i],fI, PARAMS_Simulation['dt'], tau = PARAMS_PV['tau_soma'])
    
    return None

@njit(cache=true)
def update_rate_integrator(VARIABLES_RATES,VARIABLES_CURRENTS,PARAMS_Integrator,PARAMS_Simulation):
    ''' This function will update the rate of the integrator
    '''
    for i in range(PARAMS_Integrator['Ncells']):
        
        fItemp = bf.abott_fI_curve(VARIABLES_CURRENTS['current_integrator'][i], PARAMS_Integrator['a'], PARAMS_Integrator['b'], PARAMS_Integrator['d'])
        fI = bf.deriv_integrator(VARIABLES_RATES['rate_integrator'][i],PARAMS_Integrator['alpha'] , fItemp, PARAMS_Integrator['tau_soma'])
        
        VARIABLES_RATES['rate_integrator'][i] = bf.euler_update(VARIABLES_RATES['rate_integrator'][i],fI, PARAMS_Simulation['dt'])
    
    return None

#endregion

#region update all the currents

@njit(cache=true)
def update_current_integrator(VARIABLES_CURRENTS,VARIABLES_RATES,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the current sent to the integrator
    '''
    
    VARIABLES_CURRENTS['current_integrator'] = PARAMS_Integrator['Ibg'] + WEIGHTS['W_pc_to_integrator'] @ VARIABLES_RATES['W_pc_to_integrator']
    
    
    return None

@njit(cache=true)
def update_current_vip(VARIABLES_CURRENTS,VARIABLES_RATES,PARAMS_VIP,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the current sent to the vip
    '''
    
    VARIABLES_CURRENTS['current_vip'] = PARAMS_VIP['Ibg'] + WEIGHTS['W_pc_to_vip'] @ VARIABLES_RATES['W_pc_to_vip'] + WEIGHTS['W_sst_to_vip'] @ VARIABLES_RATES['W_sst_to_vip'] + WEIGHTS['W_integrator_to_vip'] @ VARIABLES_RATES['W_integrator_to_vip']
    
    
    return None

@njit(cache=true)
def update_current_sst(VARIABLES_CURRENTS,VARIABLES_RATES,PARAMS_SST,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the current sent to the sst
    '''
    
    VARIABLES_CURRENTS['current_sst'] = PARAMS_SST['Ibg'] + WEIGHTS['W_pc_to_sst'] @ VARIABLES_RATES['W_pc_to_sst'] + WEIGHTS['W_vip_to_sst'] @ VARIABLES_RATES['W_vip_to_sst'] 
    
    
    return None

@njit(cache=true)
def update_current_ndf(VARIABLES_CURRENTS,VARIABLES_RATES,PARAMS_NDF,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the current sent to the ndf
    '''
    
    VARIABLES_CURRENTS['current_ndf'] = PARAMS_NDF['Ibg'] + WEIGHTS['W_sst_to_ndf'] @ VARIABLES_RATES['W_sst_to_ndf'] + WEIGHTS['W_integrator_to_ndf'] @ VARIABLES_RATES['W_integrator_to_ndf']
    
    
    return None

@njit(cache=true)
def update_current_pv(VARIABLES_CURRENTS,VARIABLES_RATES,PARAMS_PV,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the current sent to the pv
    '''
    
    VARIABLES_CURRENTS['current_pv'] = PARAMS_PV['Ibg'] + WEIGHTS['W_sst_to_pv'] @ VARIABLES_RATES['W_sst_to_pv'] + WEIGHTS['W_integrator_to_pv'] @ VARIABLES_RATES['W_integrator_to_pv']
    
    
    return None

@njit(cache = true)
def update_current_pc(VARIABLES_CURRENTS,VARIABLES_RATES,PARAMS_PC,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the current sent to the pc
    '''
    
    VARIABLES_CURRENTS['current_dend_Iexc'] = PARAMS_PC['Ibg_dend'] + WEIGHTS['W_integrator_to_dend'] @ VARIABLES_RATES['W_integrator_to_dend']
    
    VARIABLES_CURRENTS['current_dend_Iinh'] =  WEIGHTS['W_sst_to_dend'] @ VARIABLES_RATES['W_sst_to_dend'] + WEIGHTS['W_ndf_to_dend'] @ VARIABLES_RATES['W_ndf_to_dend']
    
    VARIABLES_CURRENTS['current_pc'] = PARAMS_PC['Ibg_soma'] +  WEIGHTS['W_pv_to_pc'] @ VARIABLES_RATES['W_pv_to_pc'] + WEIGHTS['W_pc_to_pc'] @ VARIABLES_RATES['W_pc_to_pc'] 
    
    update_dend_to_soma(VARIABLES_RATES, VARIABLES_CURRENTS, PARAMS_PC, PARAMS_Simulation)
    
    
    return None

@njit(cache=true)
def update_all_currents(VARIABLES_CURRENTS,VARIABLES_SYNAPTIC,PARAMS_Simulation,WEIGHTS, PARAMS_PC, PARAMS_Integrator, PARAMS_SST, PARAMS_VIP, PARAMS_NDF, PARAMS_PV):
    ''' This function will update all the currents in the network
    '''
    
    update_current_integrator(VARIABLES_CURRENTS,VARIABLES_SYNAPTIC,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS)
    update_current_ndf(VARIABLES_CURRENTS, VARIABLES_SYNAPTIC, PARAMS_NDF, PARAMS_Simulation, WEIGHTS)
    update_current_pv(VARIABLES_CURRENTS, VARIABLES_SYNAPTIC, PARAMS_PV, PARAMS_Simulation, WEIGHTS)
    update_current_sst(VARIABLES_CURRENTS, VARIABLES_SYNAPTIC, PARAMS_SST, PARAMS_Simulation, WEIGHTS)
    update_current_vip(VARIABLES_CURRENTS, VARIABLES_SYNAPTIC, PARAMS_VIP, PARAMS_Simulation, WEIGHTS)
    update_current_pc(VARIABLES_CURRENTS, VARIABLES_SYNAPTIC, PARAMS_PC, PARAMS_Simulation, WEIGHTS)
    
    return None


#endregion

#region update all the weights


# in this code, the matrix W has to be multiplied by the variables stored in the dictionary VARIABLES_S
# after, a simple term-term multiplication will be enough
# TODO adding facilitaiton ,depression and so on

@njit(cache=true)
def update_s_pc_to_integrator(VARIABLES_S,VARIABLES_RATES,PARAMS_Integrator,PARAMS_PC,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from s to integrator
    '''
    
    VARIABLES_S['W_pc_to_integrator'] = VARIABLES_S['W_pc_to_integrator'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_PC['tau_nmda'] + PARAMS_PC['gamma_nmda'] * (1.0 - VARIABLES_S['W_pc_to_integrator']) * VARIABLES_RATES['rate_pc'] )
    
    return None

@njit(cache=true)
def update_s_pc_to_sst(VARIABLES_S,VARIABLES_RATES,PARAMS_SST,PARAMS_PC,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from pc to sst
    '''
    
    VARIABLES_S['W_pc_to_sst'] = VARIABLES_S['W_pc_to_sst'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_PC['tau_nmda'] + PARAMS_PC['gamma_nmda']* (1.0 - VARIABLES_S['W_pc_to_sst']) * VARIABLES_RATES['rate_pc'] )
    
    return None

@njit(cache=true)
def update_s_pc_to_vip(VARIABLES_S,VARIABLES_RATES,PARAMS_VIP,PARAMS_PC,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from pc to vip
    '''
    
    VARIABLES_S['W_pc_to_vip'] = VARIABLES_S['W_pc_to_vip'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_PC['tau_nmda'] + PARAMS_PC['gamma_nmda'] * (1.0 - VARIABLES_S['W_pc_to_vip']) * VARIABLES_RATES['rate_pc'] )
    
    return None

@njit(cache=true)
def update_s_pc_to_pv(VARIABLES_S,VARIABLES_RATES,PARAMS_PV,PARAMS_PC,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from pc to pv
    '''
    
    VARIABLES_S['W_pc_to_pv'] = VARIABLES_S['W_pc_to_pv'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_PC['tau_nmda'] + PARAMS_PC['gamma_nmda'] * (1.0 - VARIABLES_S['W_pc_to_pv'])* VARIABLES_RATES['rate_pc'] )
    
    return None

@njit(cache=true)
def update_s_pc_to_pc(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from pc to pc
    '''
    
    VARIABLES_S['W_pc_to_pc'] = VARIABLES_S['W_pc_to_pc'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_PC['tau_nmda'] + PARAMS_PC['gamma_nmda'] * (1.0 - VARIABLES_S['W_pc_to_pc']) * VARIABLES_RATES['rate_pc'] )
    
    return None

# now update the weights to pv
@njit(cache=true)
def update_s_pv_to_pc(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_PV,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from pv to pc
    '''
    
    VARIABLES_S['W_pv_to_pc'] = VARIABLES_S['W_pv_to_pc'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_PV['tau_gaba'] + PARAMS_PV['gamma_nmda'] * VARIABLES_RATES['rate_pv'] )
    
    return None

@njit(cache=true)
def update_s_pv_to_pv(VARIABLES_S,VARIABLES_RATES,PARAMS_PV,PARAMS_PV,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from pv to pv
    '''
    
    VARIABLES_S['W_pv_to_pv'] = VARIABLES_S['W_pv_to_pv'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_PV['tau_gaba'] + PARAMS_PV['gamma_gaba'] * VARIABLES_RATES['rate_pv'] )
    
    return None

# update weights from sst
@njit(cache=true)
def update_s_sst_to_vip(VARIABLES_S,VARIABLES_RATES,PARAMS_VIP,PARAMS_SST,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from sst to vip
    '''
    
    VARIABLES_S['W_sst_to_vip'] = VARIABLES_S['W_sst_to_vip'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_SST['tau_gaba'] + PARAMS_SST['gamma_gaba'] * VARIABLES_RATES['rate_sst'] )
    
    return None

@njit(cache=true)
def update_s_sst_to_dend(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_SST,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from sst to dend
    '''
    
    VARIABLES_S['W_sst_to_dend'] = VARIABLES_S['W_sst_to_dend'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_SST['tau_gaba'] + PARAMS_SST['gamma_gaba'] * VARIABLES_RATES['rate_sst'] )
    
    return None

@njit(cache=true)
def update_s_sst_to_pv(VARIABLES_S,VARIABLES_RATES,PARAMS_PV,PARAMS_SST,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from sst to pv
    '''
    
    VARIABLES_S['W_sst_to_pv'] = VARIABLES_S['W_sst_to_pv'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_SST['tau_gaba'] + PARAMS_SST['gamma_gaba'] * VARIABLES_RATES['rate_sst'] )
    
    return None

@njit(cache=true)
def update_s_sst_to_ndf(VARIABLES_S,VARIABLES_RATES,PARAMS_NDF,PARAMS_SST,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from sst to ndf
    '''
    
    VARIABLES_S['W_sst_to_ndf'] = VARIABLES_S['W_sst_to_ndf'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_SST['tau_gaba'] + PARAMS_SST['gamma_gaba'] * VARIABLES_RATES['rate_sst'] )
    
    return None

#updte weights from vip

@njit(cache=true)
def update_s_vip_to_sst(VARIABLES_S,VARIABLES_RATES,PARAMS_SST,PARAMS_VIP,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from vip to sst
    '''
    
    VARIABLES_S['W_vip_to_sst'] = VARIABLES_S['W_vip_to_sst'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_VIP['tau_gaba'] + PARAMS_VIP['gamma_gaba'] * VARIABLES_RATES['rate_vip'] )
    
    return None

#updte weights from ndf
@njit(cache=true)
def update_s_ndf_to_dend(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_NDF,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from ndf to dend
    '''
    
    VARIABLES_S['W_ndf_to_dend'] = VARIABLES_S['W_ndf_to_dend'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_NDF['tau_gaba'] + PARAMS_NDF['gamma_gaba'] * VARIABLES_RATES['rate_ndf'] )
    
    return None


#update weights from integrator
@njit(cache=true)
def update_s_integrator_to_dend(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from integrator to dend
    '''
    
    VARIABLES_S['W_integrator_to_dend'] = VARIABLES_S['W_integrator_to_dend'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_Integrator['tau_nmda'] + PARAMS_Integrator['gamma_nmda'] * (1.0 - VARIABLES_S['W_integrator_to_dend'])* VARIABLES_RATES['rate_integrator'] )
    
    return None

@njit(cache=true)
def update_s_integrator_to_pv(VARIABLES_S,VARIABLES_RATES,PARAMS_PV,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from integrator to pv
    '''
    
    VARIABLES_S['W_integrator_to_pv'] = VARIABLES_S['W_integrator_to_pv'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_Integrator['tau_nmda'] + PARAMS_Integrator['gamma_nmda'] * (1.0 - VARIABLES_S['W_integrator_to_pv'])* VARIABLES_RATES['rate_integrator'] )
    
    return None

@njit(cache=true)
def update_s_integrator_to_vip(VARIABLES_S,VARIABLES_RATES,PARAMS_VIP,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from integrator to vip
    '''
    
    VARIABLES_S['W_integrator_to_vip'] = VARIABLES_S['W_integrator_to_vip'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_Integrator['tau_nmda'] + PARAMS_Integrator['gamma_nmda'] * (1.0 - VARIABLES_S['W_integrator_to_vip'])* VARIABLES_RATES['rate_integrator'] )
    
    return None

@njit(cache=true)
def update_s_integrator_to_ndf(VARIABLES_S,VARIABLES_RATES,PARAMS_NDF,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS):
    ''' This function will update the weights from integrator to ndf
    '''
    
    VARIABLES_S['W_integrator_to_ndf'] = VARIABLES_S['W_integrator_to_ndf'] + PARAMS_Simulation['dt'] * (-VARIABLES_S / PARAMS_Integrator['tau_nmda'] + PARAMS_Integrator['gamma_nmda'] * (1.0 - VARIABLES_S['W_integrator_to_ndf'])* VARIABLES_RATES['rate_integrator'] )
    
    return None

#update everything
def update_all_synapses(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_PV,PARAMS_SST,PARAMS_VIP,PARAMS_NDF,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS):
    ''' This function will update all the synapses
    '''
    
    update_s_pc_to_pv(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_PV,PARAMS_Simulation,WEIGHTS)
    update_s_pc_to_sst(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_SST,PARAMS_Simulation,WEIGHTS)
    update_s_pc_to_vip(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_VIP,PARAMS_Simulation,WEIGHTS)
    update_s_pc_to_integrator(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS)
    update_s_pc_to_pc(VARIABLES_S,VARIABLES_RATES,PARAMS_PC,PARAMS_PC,PARAMS_Simulation,WEIGHTS)

    
    update_s_pv_to_pv(VARIABLES_S,VARIABLES_RATES,PARAMS_PV,PARAMS_PV,PARAMS_Simulation,WEIGHTS)
    update_s_pv_to_pc(VARIABLES_S,VARIABLES_RATES,PARAMS_PV,PARAMS_PC,PARAMS_Simulation,WEIGHTS)
    
    update_s_sst_to_dend(VARIABLES_S, VARIABLES_RATES, PARAMS_PC, PARAMS_SST, PARAMS_Simulation, WEIGHTS)
    update_s_sst_to_vip(VARIABLES_S, VARIABLES_RATES, PARAMS_VIP, PARAMS_SST, PARAMS_Simulation, WEIGHTS)
    update_s_sst_to_pv(VARIABLES_S, VARIABLES_RATES, PARAMS_PV, PARAMS_SST, PARAMS_Simulation, WEIGHTS)
    update_s_sst_to_ndf(VARIABLES_S, VARIABLES_RATES, PARAMS_NDF, PARAMS_SST, PARAMS_Simulation, WEIGHTS)
    
    update_s_vip_to_sst(VARIABLES_S, VARIABLES_RATES, PARAMS_SST, PARAMS_VIP, PARAMS_Simulation, WEIGHTS)
    
    update_s_ndf_to_dend(VARIABLES_S, VARIABLES_RATES, PARAMS_PC, PARAMS_NDF, PARAMS_Simulation, WEIGHTS)
    
    update_s_integrator_to_dend(VARIABLES_S,VARIABLES_RATES,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS)
    update_s_integrator_to_pv(VARIABLES_S,VARIABLES_RATES,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS)
    update_s_integrator_to_vip(VARIABLES_S,VARIABLES_RATES,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS)
    update_s_integrator_to_ndf(VARIABLES_S,VARIABLES_RATES,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS)
    
    
    return None

#endregion

#region update everything


@njit(cache=true)
def time_step(VARIABLES_RATES, VARIABLES_SYNAPTIC, VARIABLES_CURRENTS, PARAMS_PC, PARAMS_PV, PARAMS_SST, PARAMS_VIP, PARAMS_NDF, PARAMS_Integrator, PARAMS_Simulation, WEIGHTS):
    ''' Compute one time step of the dynamical system: synapses then currents then rates '''
    
  
    update_all_synapses(VARIABLES_SYNAPTIC,VARIABLES_RATES,PARAMS_PC,PARAMS_PV,PARAMS_SST,PARAMS_VIP,PARAMS_NDF,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS)
    
    update_all_currents(VARIABLES_CURRENTS,VARIABLES_SYNAPTIC,PARAMS_PC,PARAMS_PV,PARAMS_SST,PARAMS_VIP,PARAMS_NDF,PARAMS_Integrator,PARAMS_Simulation,WEIGHTS)
    
    update_rate_integrator(VARIABLES_RATES, VARIABLES_CURRENTS, PARAMS_Integrator, PARAMS_Simulation)
    update_rate_ndf(VARIABLES_RATES, VARIABLES_CURRENTS, PARAMS_NDF, PARAMS_Simulation)
    update_rate_sst(VARIABLES_RATES, VARIABLES_CURRENTS, PARAMS_SST, PARAMS_Simulation)
    update_rate_vip(VARIABLES_RATES, VARIABLES_CURRENTS, PARAMS_VIP, PARAMS_Simulation)
    update_rate_pc(VARIABLES_RATES, VARIABLES_CURRENTS, PARAMS_PC, PARAMS_Simulation)
    update_rate_pv(VARIABLES_RATES, VARIABLES_CURRENTS, PARAMS_PV, PARAMS_Simulation)
    
    return None

def full_dynamics(VARIABLES_RATES, VARIABLES_SYNAPTIC, VARIABLES_CURRENTS, PARAMS_PC, PARAMS_PV, PARAMS_SST, PARAMS_VIP, PARAMS_NDF, PARAMS_Integrator, PARAMS_Simulation, WEIGHTS):
    ''' Compute the dynamics of the system for all the time steps '''
    
    nbr_steps = int(PARAMS_Simulation['t_total']/PARAMS_Simulation['dt'])
    for i in range(nbr_steps):
        time_step(VARIABLES_RATES, VARIABLES_SYNAPTIC, VARIABLES_CURRENTS, PARAMS_PC, PARAMS_PV, PARAMS_SST, PARAMS_VIP, PARAMS_NDF, PARAMS_Integrator, PARAMS_Simulation, WEIGHTS)
    
    
    return None


#endregion