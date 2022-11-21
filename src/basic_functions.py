import numpy as np
from numba import njit

import logging
logger = logging.getLogger(__name__)

@njit(cache=True)
def sigmoid(x):
    return 1 / (1 + np.exp(-x)) 

@njit(cache=True)
def relu(x):
    return np.maximum(x, 0)

@njit(cache=True)
def relu_deriv(x):
    return np.heaviside(x, 0)

@njit(cache=True)
def abott_fI_curve(x, a, b,d):
        temp=np.divide(a*x - b, 1.0 - np.exp(-d*(a*x - b)))
        return relu(temp)
    
    
@njit(cache=True)
def current_from_dend(Iinh,Iexc,c1,c2,c3,c4,c5,c6):
        return c1 * (-0.5 + sigmoid((Iexc - c2 * Iinh + c6 ) / (c3 * Iinh + c4))) + c5

@njit(cache=True)
def gaussian(x, mu, sigma):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))

@njit(cache=True)
def gaussian_connectivity(x,mu, Jmin, Jmax, sigma):
    return Jmin + Jmax*gaussian(x, mu, sigma)


def create_ring_connectivity(NCells1, NCells2, Jmax, sigma, Jmin = 0.0):
        ''' Create the ring connectivity from Layer 1 to Layer 2
        '''
        logger.info('Creating ring connectivity')
        theta_range1 = np.linspace(0.0, 360.0 - 360.0/NCells1, NCells1)
        theta_range2 = np.linspace(0.0, 360.0 - 360.0/NCells2, NCells2)

        the_rings = np.zeros((NCells1, NCells2))
        for count,theta0 in enumerate(theta_range2):
            min_dir = np.abs(theta_range1 - theta0)
            min_dir[min_dir>180.0] = min_dir[min_dir>180.0] - 360.0
            gauss_con = np.exp(-0.5*(min_dir)**2/sigma**2)
            ring  = Jmax * gauss_con + Jmin
            the_rings[count,:] = ring/sigma*np.sqrt(np.pi)#NCells1#*np.sqrt(2.0*sigma**2/np.pi)
        return the_rings
    
    
@njit(cache=True)    
def deriv_integrator(x,alpha,fI,tau):
    return (-x + alpha*x + fI )/tau   
    
    
@njit(cache=True)
def euler_update(x, f, dt, tau = 1.0):
    return relu(x + dt * f / tau)

@njit(cache=True)
def leaky_integrator(x, f, dt, tau = 1.0, alpha = 1.0):
    return euler_update(x, (-alpha*x + f), dt, tau)

def generate_ring_stim(center, sigma, N, tau):
    ''' Fucntion that generates an exponential stimulus
    '''
    logger.info('Generating ring stimulus')
    theta_range1 = np.linspace(0.0, 360.0 - 360.0/N, N)
    center = center*360.0/N
    stim = np.zeros(N)
    min_dir = np.abs(theta_range1 - center)
    min_dir[min_dir>180.0] = min_dir[min_dir>180.0] - 360.0
  
    stim = tau*np.exp(-(min_dir)**2/(2*sigma**2))
    return stim 

@njit(cache=True)
def soft_bound(x, xmax,beta = 1.0):
    return np.power(np.sign(xmax -x)*np.abs((xmax -x)),beta)

@njit(cache=True)
def Hebb_with_decay(w, rpre, rpost, gamma, lambda_dec, wmax = 1.0, beta = 1.0):
    ''' Return de Hebbian update rule with decay
    It takes am atrix into argument
    '''
    rtot = np.outer(rpre,rpost)
    # let's use the fact that for now, connections are sparse from ndf to dend
    
    return gamma * rtot * soft_bound(w, wmax, beta) - lambda_dec*w

@njit(cache=True)
def anti_Hebb_with_decay(w, rpre, rpost, gamma, lambda_dec, wmax = 1.0, beta = 1.0):
    ''' Return de anti Hebbian update rule with decay
    It takes am atrix into argument
    '''
    rtot = np.outer(rpre,rpost)
    # let's use the fact that for now, connections are sparse from ndf to dend
    
    return -gamma * rtot * soft_bound(w, wmax, beta) - lambda_dec*w

@njit(cache=True)
def BCM_rule(w, rpre, rpost, theta, epsilon = 0.0):
    ''' Return the BCM rule with theta as the average of the firing rates
        To compute the average of the firing raqte, the trick is to multiply by exp(-deltat/tau) at each time step
        Then the average is the weighted sum? (not technically correct but good enough?)
    '''
    router = np.outer(rpre,rpost)
    theta = (theta*np.exp(-0.001/0.1) + rpost*rpost*(1-np.exp(-0.001/0.1)))
    temp = 0.0001*(rpost - theta)*router# - epsilon*w
    
    
    return temp, theta


@njit(cache=True)
def threshold_matrix(m,threshold):
    ''' Return the matrix with saturated values above threshold
    '''
    return np.minimum(m,threshold)


# @njit(cache=True)
def compute_list_time_input(dt, ttotal, tinter, tresting, tstim, delay = 0.1):
    # temp = [t for t in np.arange(tresting, ttotal, dt)]
    
    max_iter = int(np.floor((ttotal - tresting)/(tinter+tstim)))
    temp = np.arange(tresting + tinter + delay,ttotal,tinter + tstim)/dt
    temp = temp.astype(int)
    return temp

def compute_MMN_effect(pc_layer,PARAMS_Stimulus):
    ''' Compute the MMN effect of the network
    '''
    # PARAMS_Stimulus = {'tau': 0.1, 't_stim': 0.1, 't_rest': 0.1, 't_inter': 0.1, 'N_stim': 1, 'N_rest': 1, 'N_inter': 1}
    tau = PARAMS_Stimulus['tau']
    t_stim = PARAMS_Stimulus['t_stim']
    t_rest = PARAMS_Stimulus['t_rest']
    t_inter = PARAMS_Stimulus['t_inter']
    N_stim = PARAMS_Stimulus['N_stim']
    N_rest = PARAMS_Stimulus['N_rest']
    N_inter = PARAMS_Stimulus['N_inter']
    dt = PARAMS_Stimulus['dt']
    ttotal = PARAMS_Stimulus['ttotal']
    tresting = PARAMS_Stimulus['tresting']
    tstim = PARAMS_Stimulus['tstim']
    tinter = PARAMS_Stimulus['tinter']
    
    return pc_layer[:,nbr_rep-1] - pc_layer[:,nbr_rep]
