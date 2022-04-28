import numpy as np
from numba.experimental import jitclass
from numba import float64, int32
import pandas as pd
import basic_functions as bf
import pickle
import json

#region General Neuron


class NeuronLayer:
    '''this class will define the general methods of a neuron
    '''
    def __init__(self):
        self.rate_now = None
        self.input_current = None
        
    
    def get_rate(self):
        return self.rate_now
    
    
    def get_input_current(self):
        return self.input_current
    
    
    def set_input_current(self, external_current):
        self.input_current = external_current
        
    
    def set_rate(self, new_rate):
        self.rate_now = new_rate
    
    
#endregion

#region IntegratorLayer class
#region jitclass Integrator class
spec = [
    ('Ncelss', int32),
    ('alpha', float64),
    ('a', float64),
    ('b', float64),
    ('d', float64),
    ('Ibg', float64),
    ('tau',float64),
    ('sigma_con', float64),
    ('Jmin', float64),
    ('Jmax', float64),
    ('rate_now',float64[:]),
    ('input_current', float64[:]),
    ('S_NMDA', float64[:]),
    ('weights', float64[:]),
    ('weights_NMDA', float64[:])
    ]
#endregion
# @jitclass(spec)
class Integrator_Layer(NeuronLayer):
    def __init__(self,  Ncells, alpha, a, b, d, Ibg,tau, sigma_con, Jmin, Jmax,**kwargs):
                     self.Ncells = Ncells
                     self.alpha = alpha
                     self.a = a
                     self.b = b
                     self.d = d
                     self.Ibg = Ibg
                     self.tau = tau
                     
                     # synapse parameters
                    #  self.tau_nmda = tau_nmda
                    #  self.gamma_nmda = gamma_nmda
                     
                    # connectivity parameters
                     self.sigma_con = sigma_con
                     self.Jmin = Jmin
                     self.Jmax = Jmax
                     
                     self.Iext = np.zeros(self.Ncells)
                     
                     self.input_current = np.zeros(self.Ncells)
                     self.rate_now = 0.0*np.ones(self.Ncells)
                     self.S_NMDA = 0.1*np.ones(self.Ncells)

                     self.weights = np.zeros((self.Ncells, self.Ncells))
                     self.weights_NMDA = np.zeros( self.Ncells)
                     
                    #  self.create_ring_connectivity() #create the gaussian matrix of connections
                     
                     
    

        
    def create_ring_connectivity(self):
            
        theta_range = np.linspace(0.0, 360.0 - 360.0/self.Ncells, self.Ncells)
        the_rings = np.zeros((self.Ncells, self.Ncells))
        for count,theta0 in enumerate(theta_range):
            min_dir = np.abs(theta_range - theta0)
            min_dir[min_dir>180.0] = min_dir[min_dir>180.0] - 360.0
            gauss_con = np.exp(-0.5*(min_dir)**2/self.sigma_con**2)
            ring  = self.Jmax * gauss_con + self.Jmin
            the_rings[count,:] = ring/self.Ncells
        self.weights = the_rings
     
                        
    def fI_curve(self):
        self.set_input_current(self.get_input_current() + self.Ibg + self.Iext) 
        self.set_input_current(bf.relu(self.input_current))
        
        return bf.abott_fI_curve(self.input_current, self.a, self.b, self.d)
                     
    # def dt_NMDA(self):
    #     return -self.S_NMDA/self.tau_nmda + self.gamma_nmda*(1 - self.S_NMDA)*self.rate_now 
    
    
    # def dt_rate(self):
    #     return self.alpha*self.rate_now + self.fI_curve()
    
    
    # def update_r(self,dt):
    #     temp = (-self.get_rate() + self.dt_rate())*dt/self.tau 
    #     self.set_rate(self.get_rate() +  temp)
    
    
    def update_input_current(self, external_current):
        self.set_input_current(external_current + self.Ibg)
        
    # def dr_dt(self,dt,external_stimulus):
    #     self.S_NMDA += self.dt_NMDA()*dt
    #     self.S_NMDA = np.minimum(self.S_NMDA,1.0)
    #     self.weights_NMDA = self.weights @ self.S_NMDA
    #     temp_current =  external_stimulus  + self.weights_NMDA  
    #     self.update_input_current(temp_current)
    #     self.update_r(dt)
     
      
    def dr_dt(self,dt):
        # self.update_r(dt)
        temp =  bf.deriv_integrator(self.rate_now, self.alpha, self.fI_curve(),self.tau)
        temp_rate = bf.euler_update(self.rate_now,temp, dt)
        self.set_rate(temp_rate)
    
    
    def get_weights(self):
        return self.weights
    
#endregion
    

#region PC_Layer class
class PC_Layer(NeuronLayer):
    def __init__(self, Ncells, tau_soma, tau_dend, tau_NMDA, gamma_NMDA,Ibg_soma, a, b, d,tau_adaptation, gA,**kwargs):
        self.Ncells = Ncells
        self.tau_soma = tau_soma
        self.tau_dend = tau_dend
        self.tau_nmda = tau_NMDA
        self.gamma_nmda = gamma_NMDA
       
        self.Iext = np.zeros(self.Ncells)
       
        self.rate_now = 5.0*np.ones(self.Ncells)
        self.S_NMDA = 0.1*np.ones(self.Ncells)
        self.weights = np.zeros((self.Ncells, self.Ncells))
        self.weights_NMDA = np.zeros( self.Ncells)
        
        self.input_current = np.zeros(self.Ncells)
        self.Iexc_to_dend = np.zeros(self.Ncells)
        self.Iinh_to_dend = np.zeros(self.Ncells)
        self.Idend_output = np.zeros(self.Ncells)
        
        self.param_dend = { 'c1' : 0.12,
                            'c2' : -7.0,
                            'c3': -0.482,
                            'c4': 0.00964,
                            'c5': 0.19624,
                            'c6': 0.0}
        
        self.a = a
        self.b = b
        self.d = d
        
        self.Ibg_soma = Ibg_soma
        
        self.sA = np.zeros(self.Ncells)
        self.tau_adaptation = tau_adaptation
        self.gA = gA
        
        
    def create_ring_connectivity(self):
        theta_range = np.linspace(0.0, 360.0 - 360.0/self.Ncells, self.Ncells)
        the_rings = np.zeros((self.Ncells, self.Ncells))
        for count,theta0 in enumerate(theta_range):
            min_dir = np.abs(theta_range - theta0)
            min_dir[min_dir>180.0] = min_dir[min_dir>180.0] - 360.0
            gauss_con = np.exp(-0.5*(min_dir)**2/self.sigma_con**2)
            ring  = self.Jmax * gauss_con + self.Jmin
            the_rings[count,:] = ring/self.Ncells
        self.weights = the_rings
    
    
    def get_exc_input_current(self):
        return self.Iexc
    
    
    def get_inh_input_current(self):
        return self.Iinh
    
    
    def set_exc_input_current(self, external_current):
        self.Iexc = external_current
        
    
    def set_inh_input_current(self, external_current):
        self.Iinh = external_current
        
        
    def current_output_dendrite(self,Iinh,Iexc,c1,c2,c3,c4,c5,c6):
        self.Idend_output = bf.current_from_dend(Iinh,Iexc,c1,c2,c3,c4,c5,c6)
        return self.Idend_output
    
    
    def fI_curve(self,dt):
        self.sA = bf.euler_update(self.sA, -self.sA/self.tau_adaptation + self.get_rate(), dt)
        self.Idend_output = self.current_output_dendrite(self.Iinh_to_dend,self.Iexc_to_dend,**self.param_dend)
        self.set_input_current(  self.get_input_current() + self.Ibg_soma + self.Idend_output  + self.Iext + self.gA*self.sA)
        self.set_input_current(bf.relu(self.get_input_current()))
        
      
        return bf.abott_fI_curve(self.get_input_current(), self.a, self.b, self.d)
    
        
    def dr_dt(self,dt):
        self.set_rate(np.maximum(self.rate_now,0.0))
              
        # temp = (-self.get_rate() + self.fI_curve(dt))*dt/self.tau_soma 
        self.set_rate(bf.leaky_integrator(self.rate_now, self.fI_curve(dt),dt, self.tau_soma)) 
     
#endregion   

#region PV_Layer class
class PV_Network(NeuronLayer):
    '''
    This class defines the arguments and methods of the PV method
    '''
    
    def __init__(self, Ncells, c_I, r0, Ibg, tau, **kwargs):
        self.Ncells = Ncells
        self.c_I = c_I
        self.r0 = r0
        self.Ibg = Ibg
        self.tau = tau
        
        self.input_current = np.zeros(self.Ncells)
        self.rate_now = np.zeros(self.Ncells)
        self.weights = np.zeros((self.Ncells, self.Ncells))
        
    
    def fI_curve(self):
        self.input_current += self.Ibg
        return bf.relu(self.c_I * self.get_input_current() + self.r0)
    
    
    def dr(self):
        return (-self.get_rate() + self.fI_curve())/self.tau
    
    
    def dr_dt(self,dt):
        self.set_rate(bf.leaky_integrator(self.rate_now, self.fI_curve(), dt, self.tau))
        # self.set_rate(bf.euler_update(self.rate_now, , dt))
        
#endregion       
        
#region VIP_Network        
class VIP_Network(NeuronLayer):
    '''
    This class defines the arguments and methods of the VIP method
    '''
    
    
    def __init__(self, Ncells, c_I, r0, Ibg, tau, **kwargs):
        self.Ncells = Ncells
        self.c_I = c_I
        self.r0 = r0
        self.Ibg = Ibg
        self.tau = tau
        
        self.input_current = np.zeros(self.Ncells)
        self.rate_now = np.zeros(self.Ncells)
        self.weights = np.zeros((self.Ncells, self.Ncells))
    
    
    
    def fI_curve(self):
        self.input_current += self.Ibg
        return bf.relu(self.c_I * self.get_input_current() + self.r0)
    
    
    def dr(self):
        return (-self.get_rate() + self.fI_curve())/self.tau

    
    
    def dr_dt(self,dt):
        self.set_rate(bf.leaky_integrator(self.rate_now, self.fI_curve(), dt, self.tau))

        # self.set_rate(bf.euler_update(self.rate_now, self.dr(), dt))
#endregion         

#region SST_Network
class SST_Network(NeuronLayer):
    '''
    This class defines the arguments and methods of the SST method
    '''
    
    
    def __init__(self, Ncells, c_I, r0, Ibg, tau, **kwargs):
        self.Ncells = Ncells
        self.c_I = c_I
        self.r0 = r0
        self.Ibg = Ibg
        self.tau = tau
        
        self.input_current = np.zeros(self.Ncells)
        self.rate_now = np.zeros(self.Ncells)
        self.weights = np.zeros((self.Ncells, self.Ncells))

        self.sA = np.zeros(self.Ncells)
        self.tau_adaptation = kwargs['tau_adaptation']
        self.gA = kwargs['gA']
    
    
    def fI_curve(self,dt):
        self.sA = bf.euler_update(self.sA, -self.sA/self.tau_adaptation + self.get_rate(), dt)
        self.input_current += self.Ibg + self.gA * self.sA
        return bf.relu(self.c_I * self.get_input_current() + self.r0)
    
    
    def dr(self,dt):
        return (-self.get_rate() + self.fI_curve(dt))/self.tau

    
    
    def dr_dt(self,dt):
        self.set_rate(bf.leaky_integrator(self.rate_now, self.fI_curve(dt), dt, self.tau))

        # self.set_rate(bf.euler_update(self.rate_now, self.dr(dt), dt))
#endregion            


#region NDF_Network
class NDF_Network(NeuronLayer):
    '''
    This class defines the arguments and methods of the NDF method
    '''
    
    
    def __init__(self, Ncells, c_I, r0, Ibg, tau, **kwargs):
        self.Ncells = Ncells
        self.c_I = c_I
        self.r0 = r0
        self.Ibg = Ibg
        self.tau = tau
        
        self.input_current = np.zeros(self.Ncells)
        self.rate_now = np.zeros(self.Ncells)
        self.weights = np.zeros((self.Ncells, self.Ncells))
        self.theta = np.zeros(self.Ncells)
     
    
    def fI_curve(self):
        self.input_current += self.Ibg
        return bf.relu(self.c_I * self.get_input_current() + self.r0)
    
    
    def dr(self):
        return (-self.get_rate() + self.fI_curve())/self.tau

    
    
    def dr_dt(self,dt):
        
        self.set_rate(bf.leaky_integrator(self.rate_now, self.fI_curve(), dt, self.tau))

        # self.set_rate(bf.euler_update(self.rate_now, self.dr(), dt))
#endregion          

        
        
class Network:
    ''' This class will define the network will all the different layers, interneurons and connections.
    It will define the dynamics methods too.
    '''
    
    def __init__(self, PARAMS_Integrator,PARAMS_Synapses_Integrator, PARAMS_PC, PARAMS_PV, PARAMS_SST, PARAMS_VIP, PARAMS_NDF, PARAMS_Simulation, PARAMS_Stimulus, **kwargs):
        
        self.PARAMS_Integrator = {**PARAMS_Integrator, **PARAMS_Synapses_Integrator}
        self.PARAMS_PC = PARAMS_PC
        self.PARAMS_PV = PARAMS_PV
        self.PARAMS_SST = PARAMS_SST
        self.PARAMS_VIP = PARAMS_VIP
        self.PARAMS_NDF = PARAMS_NDF
        
        self.integrator_layer = Integrator_Layer(**PARAMS_Integrator)
        self.pc_layer = PC_Layer(**PARAMS_PC)
        self.pv_layer = PV_Network(**PARAMS_PV)
        self.ndf_layer = NDF_Network(**PARAMS_NDF)
        self.sst_layer = SST_Network(**PARAMS_SST)
        self.vip_layer = VIP_Network(**PARAMS_VIP)
        
        self.PARAMS_Simulation = PARAMS_Simulation
        self.PARAMS_Stimulus = PARAMS_Stimulus
        
        #region Define the connections between the different network
        self.W_sst_to_vip = np.zeros((self.sst_layer.Ncells, self.vip_layer.Ncells))
        self.W_sst_to_dend = np.zeros((self.sst_layer.Ncells, self.pc_layer.Ncells))
        self.W_sst_to_pv = np.zeros((self.sst_layer.Ncells, self.pv_layer.Ncells))
        self.W_sst_to_ndf = np.zeros((self.sst_layer.Ncells, self.ndf_layer.Ncells))
        
        self.W_vip_to_sst = np.zeros((self.vip_layer.Ncells, self.sst_layer.Ncells))
        
        self.W_ndf_to_dend = np.zeros((self.ndf_layer.Ncells, self.pc_layer.Ncells))
        self.W_ndf_to_pv = np.zeros((self.ndf_layer.Ncells, self.sst_layer.Ncells))
        
        self.W_pv_to_pc = np.zeros((self.pv_layer.Ncells, self.pc_layer.Ncells))
        self.W_pv_to_pv = np.zeros((self.pv_layer.Ncells, self.pv_layer.Ncells))
        
        self.W_pc_to_pc = np.zeros((self.pc_layer.Ncells, self.pc_layer.Ncells))
        self.W_pc_to_vip = np.zeros((self.pc_layer.Ncells, self.vip_layer.Ncells))
        self.W_pc_to_sst = np.zeros((self.pc_layer.Ncells, self.sst_layer.Ncells))
        self.W_pc_to_pv = np.zeros((self.pc_layer.Ncells, self.pv_layer.Ncells))
        self.W_pc_to_integrator = np.zeros((self.pc_layer.Ncells, self.integrator_layer.Ncells))
        
        self.W_integrator_to_dend = np.zeros((self.integrator_layer.Ncells, self.pc_layer.Ncells))
        self.W_integrator_to_pv = np.zeros((self.integrator_layer.Ncells, self.pv_layer.Ncells))
        self.W_integrator_to_vip = np.zeros((self.integrator_layer.Ncells, self.vip_layer.Ncells))
        self.W_integrator_to_ndf = np.zeros((self.integrator_layer.Ncells, self.ndf_layer.Ncells))
        self.W_integrator_to_integrator = np.zeros((self.integrator_layer.Ncells, self.integrator_layer.Ncells))
        
        #now definition of the S variables for all the connections
        self.S_sst_to_vip = np.zeros((self.sst_layer.Ncells, self.vip_layer.Ncells))
        self.S_sst_to_dend = np.zeros((self.sst_layer.Ncells, self.pc_layer.Ncells))
        self.S_sst_to_pv = np.zeros((self.sst_layer.Ncells, self.pv_layer.Ncells))
        self.S_sst_to_ndf = np.zeros((self.sst_layer.Ncells, self.ndf_layer.Ncells))
        self.S_ndf_to_pv = np.zeros((self.ndf_layer.Ncells, self.sst_layer.Ncells))

        self.S_vip_to_sst = np.zeros((self.vip_layer.Ncells, self.sst_layer.Ncells))
        
        self.S_ndf_to_dend = np.zeros((self.ndf_layer.Ncells, self.pc_layer.Ncells))
        
        self.S_pv_to_pc = np.zeros((self.pv_layer.Ncells, self.pc_layer.Ncells))
        self.S_pv_to_pv = np.zeros((self.pv_layer.Ncells, self.pv_layer.Ncells))
        
        self.S_pc_to_pc = np.zeros((self.pc_layer.Ncells, self.pc_layer.Ncells))
        self.S_pc_to_vip = np.zeros((self.pc_layer.Ncells, self.vip_layer.Ncells))
        self.S_pc_to_sst = np.zeros((self.pc_layer.Ncells, self.sst_layer.Ncells))
        self.S_pc_to_pv = np.zeros((self.pc_layer.Ncells, self.pv_layer.Ncells))
        self.S_pc_to_integrator = np.zeros((self.pc_layer.Ncells, self.integrator_layer.Ncells))
        
        self.S_integrator_to_dend = np.zeros((self.integrator_layer.Ncells, self.pc_layer.Ncells))
        self.S_integrator_to_pv = 0.1*np.ones((self.integrator_layer.Ncells, self.pv_layer.Ncells))
        self.S_integrator_to_vip = np.zeros((self.integrator_layer.Ncells, self.vip_layer.Ncells))
        self.S_integrator_to_ndf = np.zeros((self.integrator_layer.Ncells, self.ndf_layer.Ncells))
        self.S_integrator_to_integrator = np.zeros((self.integrator_layer.Ncells, self.integrator_layer.Ncells))
        
        #now definition of the u or d variables for all the connections
        self.U_sst_to_vip = np.zeros((self.sst_layer.Ncells, self.vip_layer.Ncells))
        self.U_sst_to_dend = np.zeros((self.sst_layer.Ncells, self.pc_layer.Ncells))
        self.U_sst_to_pv = np.zeros((self.sst_layer.Ncells, self.pv_layer.Ncells))
        self.U_sst_to_ndf = np.zeros((self.sst_layer.Ncells, self.ndf_layer.Ncells))
        self.U_ndf_to_pv = np.zeros((self.ndf_layer.Ncells, self.sst_layer.Ncells))

        
        self.U_vip_to_sst = np.zeros((self.vip_layer.Ncells, self.sst_layer.Ncells))
        
        self.U_ndf_to_dend = np.zeros((self.ndf_layer.Ncells, self.pc_layer.Ncells))
        
        self.U_pv_to_pc = np.zeros((self.pv_layer.Ncells, self.pc_layer.Ncells))
        self.U_pv_to_pv = np.zeros((self.pv_layer.Ncells, self.pv_layer.Ncells))
        
        self.U_pc_to_pc = np.zeros((self.pc_layer.Ncells, self.pc_layer.Ncells))
        self.U_pc_to_vip = np.zeros((self.pc_layer.Ncells, self.vip_layer.Ncells))
        self.U_pc_to_sst = np.zeros((self.pc_layer.Ncells, self.sst_layer.Ncells))
        self.U_pc_to_pv = np.zeros((self.pc_layer.Ncells, self.pv_layer.Ncells))
        self.U_pc_to_integrator = np.zeros((self.pc_layer.Ncells, self.integrator_layer.Ncells))
        
        self.U_integrator_to_dend = np.zeros((self.integrator_layer.Ncells, self.pc_layer.Ncells))
        self.U_integrator_to_pv = 0.1*np.ones((self.integrator_layer.Ncells, self.pv_layer.Ncells))
        self.U_integrator_to_vip = np.zeros((self.integrator_layer.Ncells, self.vip_layer.Ncells))
        self.U_integrator_to_ndf = np.zeros((self.integrator_layer.Ncells, self.ndf_layer.Ncells))
        self.U_integrator_to_integrator = np.zeros((self.integrator_layer.Ncells, self.integrator_layer.Ncells))
        
        
        self.W_to_dend = None #to be updated later (check Pierre's code to find a nice initialization)
        self.W_to_pv = None
        self.W_to_vip = None
        self.W_to_sst = None
        self.W_to_ndf = None
        self.W_to_pc = None
        self.W_to_integrator = None
        
        self.create_internal_connectivity()
        
        #endregion
        self.PARAMS_Simulation ['N_t'] = int(self.PARAMS_Simulation['t_total']/self.PARAMS_Simulation['dt'])
        self.save_pc_layer = np.zeros((self.pc_layer.Ncells, self.PARAMS_Simulation['N_t']))
        self.save_pv_layer = np.zeros((self.pv_layer.Ncells, self.PARAMS_Simulation['N_t']))
        self.save_sst_layer = np.zeros((self.sst_layer.Ncells, self.PARAMS_Simulation['N_t']))
        self.save_vip_layer = np.zeros((self.vip_layer.Ncells, self.PARAMS_Simulation['N_t']))
        self.save_ndf_layer = np.zeros((self.ndf_layer.Ncells, self.PARAMS_Simulation['N_t']))
        self.save_integrator_layer = np.zeros((self.integrator_layer.Ncells, self.PARAMS_Simulation['N_t']))
        
        # saving some weights
        self.save_W_standard_integrator_to_ndf = np.zeros((self.ndf_layer.Ncells, self.PARAMS_Simulation['N_t']))
        self.save_W_dev_integrator_to_ndf = np.zeros((self.ndf_layer.Ncells, self.PARAMS_Simulation['N_t']))
        
        self.save_W_standard_ndf_to_dend = np.zeros((self.pc_layer.Ncells, self.PARAMS_Simulation['N_t']))
        self.save_W_dev_ndf_to_dend = np.zeros((self.pc_layer.Ncells, self.PARAMS_Simulation['N_t']))
        
        self.save_W_standard_integrator_to_dend = np.zeros((self.pc_layer.Ncells, self.PARAMS_Simulation['N_t']))
        self.save_W_dev_integrator_to_dend = np.zeros((self.pc_layer.Ncells, self.PARAMS_Simulation['N_t']))
        
        
        
        self.PARAMS_ALL = { 
    'Simulation': self.PARAMS_Simulation,
    'PARAMS_INT': self.PARAMS_Integrator,
    'PARAMS_PC': self.PARAMS_PC,
    'PARAMS_PV': self.PARAMS_PV,
    'PARAMS_SST': self.PARAMS_SST,
    'PARAMS_VIP': self.PARAMS_VIP,
    'PARAMS_NDF': self.PARAMS_NDF,
    'Stimulus': self.PARAMS_Stimulus}
    
    #region Define connectivity matrices for the network
    def create_connectivity_from_SST(self):
        '''create the weights matrices going from SST neurons
        '''
        np.fill_diagonal(self.W_sst_to_vip, self.PARAMS_SST['weight_to_vip'])
        np.fill_diagonal(self.W_sst_to_dend, self.PARAMS_SST['weight_to_dend'])
        np.fill_diagonal(self.W_sst_to_pv, self.PARAMS_SST['weight_to_pv'])
        np.fill_diagonal(self.W_sst_to_ndf, self.PARAMS_SST['weight_to_ndf'])
        
    def create_connectivity_from_VIP(self):
        '''create the weights matrices going from VIP neurons
        '''
        np.fill_diagonal(self.W_vip_to_sst, self.PARAMS_VIP['weight_to_sst'])
        
    def create_connectivity_from_NDF(self):
        '''create the weights matrices going from NDF neurons
        '''
        # np.fill_diagonal(self.W_ndf_to_dend, self.PARAMS_NDF['weight_to_dend'])
        self.W_ndf_to_dend[:, :] = bf.create_ring_connectivity(self.ndf_layer.Ncells, self.pc_layer.Ncells, self.PARAMS_NDF['weight_to_dend'], self.PARAMS_NDF['sigma_to_dend'])
        
        np.fill_diagonal(self.W_ndf_to_pv, self.PARAMS_NDF['weight_to_pv'])

    def create_connectivity_from_PV(self):
        '''create the weights matrices going from PV neurons
        '''
        np.fill_diagonal(self.W_pv_to_pc, self.PARAMS_PV['weight_to_soma'])
        np.fill_diagonal(self.W_pv_to_pv, self.PARAMS_PV['weight_to_pv'])
        
        
    
    def create_connectivity_from_PC(self):
        ''' create the weights matrices going from PC neurons
        '''
        
        np.fill_diagonal(self.W_pc_to_vip, self.PARAMS_PC['weight_to_vip'])

        self.W_pc_to_sst[:, :] = bf.create_ring_connectivity(self.pc_layer.Ncells, self.sst_layer.Ncells, self.PARAMS_PC['weight_to_sst'], self.PARAMS_PC['sigma_to_sst'])

        self.W_pc_to_pv[:, :] = bf.create_ring_connectivity(self.pc_layer.Ncells, self.pv_layer.Ncells, self.PARAMS_PC['weight_to_pv'], self.PARAMS_PC['sigma_to_pv'])

        self.W_pc_to_integrator[:, :] = bf.create_ring_connectivity(self.pc_layer.Ncells, self.integrator_layer.Ncells, self.PARAMS_PC['weight_to_integrator'], self.PARAMS_PC['sigma_to_integrator'])
        
        self.W_pc_to_pc[:, :] = bf.create_ring_connectivity(self.pc_layer.Ncells, self.pc_layer.Ncells, self.PARAMS_PC['weight_to_soma'], self.PARAMS_PC['sigma_to_pc'])

    def create_connectivity_from_Integrator(self):
        ''' create the weights matrices going from Integrator neurons
        '''
      
        self.W_integrator_to_dend[:, :] = bf.create_ring_connectivity(self.integrator_layer.Ncells, self.pc_layer.Ncells, self.PARAMS_Integrator['weight_to_dend'], self.PARAMS_Integrator['sigma_to_dend'])
        
        self.W_integrator_to_ndf[:, :] = bf.create_ring_connectivity(self.integrator_layer.Ncells, self.ndf_layer.Ncells, self.PARAMS_Integrator['weight_to_ndf'], self.PARAMS_Integrator['sigma_to_ndf'])
        
        self.W_integrator_to_pv[:, :] = bf.create_ring_connectivity(self.integrator_layer.Ncells, self.pv_layer.Ncells, self.PARAMS_Integrator['weight_to_pv'], self.PARAMS_Integrator['sigma_to_pv'])
        
        self.W_integrator_to_vip[:, :] = bf.create_ring_connectivity(self.integrator_layer.Ncells, self.vip_layer.Ncells, self.PARAMS_Integrator['weight_to_vip'], self.PARAMS_Integrator['sigma_to_vip'])
        
        self.W_integrator_to_integrator[:, :] = bf.create_ring_connectivity(self.integrator_layer.Ncells, self.integrator_layer.Ncells, self.PARAMS_Integrator['Jmax'], self.PARAMS_Integrator['sigma_con'], Jmin = self.PARAMS_Integrator['Jmin']) / self.integrator_layer.Ncells
        for i in range(self.integrator_layer.Ncells):
            self.W_integrator_to_integrator[i, i] = 0
    
    
    def create_internal_connectivity(self):
        ''' for now the only broad connections are E to E,PV and SST, as well as Int to Int, N, PV, VIP, E
        '''
        
        self.create_connectivity_from_PC()
        self.create_connectivity_from_Integrator()
        self.create_connectivity_from_SST()
        self.create_connectivity_from_PV()
        self.create_connectivity_from_VIP()
        self.create_connectivity_from_NDF()
    #endregion
       
       
    #region get the connections towards all types of neurons
    
    # 
    def get_connections_to_Integrator(self):
        ''' Get the connections to the Integrator as an arrays of tuple (pre, post, weight)
        '''
        temp_pre = []
        temp_s = []
        temp_weights = []
        
        temp_pre.append(self.pc_layer.rate_now)
        temp_s.append(self.S_pc_to_integrator)
        temp_weights.append(self.W_pc_to_integrator)
        
        temp_pre.append(self.integrator_layer.rate_now)
        temp_s.append(self.S_integrator_to_integrator)
        temp_weights.append(self.W_integrator_to_integrator)
        
        # temp.append((self.pc_layer.rate_now, self.integrator_layer.rate_now, self.W_pc_to_integrator))
        return temp_pre, temp_s, temp_weights
    
    # 
    def get_connections_to_PC(self):
        ''' Get the connections to the PC as an arrays of tuple (pre, post, weight)
        '''
        temp_pre = []
        temp_s = []
        temp_weights = []
        
        temp_pre.append(self.pv_layer.rate_now)
        temp_s.append(self.S_pv_to_pc)
        temp_weights.append(self.W_pv_to_pc)
        
        temp_pre.append(self.pc_layer.rate_now)
        temp_s.append(self.S_pc_to_pc)
        temp_weights.append(self.W_pc_to_pc)
        
        # temp.append((self.pv_layer.rate_now, self.pc_layer.rate_now, self.W_pv_to_pc))
        # temp.append((self.pc_layer.rate_now, self.pc_layer.rate_now, self.W_pc_to_pc))
        return temp_pre, temp_s, temp_weights
    
    # 
    def get_connections_to_PV(self):
        ''' Get the connections to the PV as an arrays of tuple (pre, post, weight)
        '''
        
        temp_pre = []
        temp_s = []
        temp_weights = []
        
        temp_pre.append(self.sst_layer.rate_now)
        temp_s.append(self.S_sst_to_pv)
        temp_weights.append(self.W_sst_to_pv)
        # temp.append((self.sst_layer.rate_now, self.pv_layer.rate_now, self.W_sst_to_pv))
        
        temp_pre.append(self.pv_layer.rate_now)
        temp_s.append(self.S_pv_to_pv)
        temp_weights.append(self.W_pv_to_pv)
        # temp.append((self.pv_layer.rate_now, self.pv_layer.rate_now, self.W_pv_to_pv))
        
        temp_pre.append(self.pc_layer.rate_now)
        temp_s.append(self.S_pc_to_pv)
        temp_weights.append(self.W_pc_to_pv)
        # temp.append((self.pc_layer.rate_now, self.pv_layer.rate_now, self.W_pc_to_pv))
        
        temp_pre.append(self.integrator_layer.rate_now)
        temp_s.append(self.S_integrator_to_pv)
        temp_weights.append(self.W_integrator_to_pv)
        # temp.append((self.integrator_layer.rate_now, self.pv_layer.rate_now, self.W_integrator_to_pv))
        
        temp_pre.append(self.ndf_layer.rate_now)
        temp_s.append(self.S_ndf_to_pv)
        temp_weights.append(self.W_ndf_to_pv)
        
        return temp_pre, temp_s, temp_weights
      
    #   
    def get_connections_to_SST(self):
        ''' Get the connections to the SST as an arrays of tuple (pre, post, weight)
        '''
        temp_pre = []
        temp_s = []
        temp_weights = []
        
        temp_pre.append(self.vip_layer.rate_now)
        temp_s.append(self.S_vip_to_sst)
        temp_weights.append(self.W_vip_to_sst)
        # temp.append((self.vip_layer.rate_now, self.sst_layer.rate_now, self.W_vip_to_sst))
        
        temp_pre.append(self.pc_layer.rate_now)
        temp_s.append(self.S_pc_to_sst)
        temp_weights.append(self.W_pc_to_sst)
        # temp.append((self.pc_layer.rate_now, self.sst_layer.rate_now, self.W_pc_to_sst))
        
        
        
        return temp_pre, temp_s, temp_weights
    
    # 
    def get_connections_to_VIP(self):
        ''' Get the connections to the VIP as an arrays of tuple (pre, post, weight)
        '''
        temp_pre = []
        temp_s = []
        temp_weights = []
        
        temp_pre.append(self.sst_layer.rate_now)
        temp_s.append(self.S_sst_to_vip)
        temp_weights.append(self.W_sst_to_vip)
        # temp.append((self.sst_layer.rate_now, self.vip_layer.rate_now, self.W_sst_to_vip))
        
        temp_pre.append(self.pc_layer.rate_now)
        temp_s.append(self.S_pc_to_vip)
        temp_weights.append(self.W_pc_to_vip)
        # temp.append((self.pc_layer.rate_now, self.vip_layer.rate_now, self.W_pc_to_vip))
        
        temp_pre.append(self.integrator_layer.rate_now)
        temp_s.append(self.S_integrator_to_vip)
        temp_weights.append(self.W_integrator_to_vip)
        # temp.append((self.integrator_layer.rate_now, self.vip_layer.rate_now, self.W_integrator_to_vip))
        return temp_pre, temp_s, temp_weights
    
    # 
    def get_connections_to_NDF(self):
        ''' Get the connections to the NDF as an arrays of tuple (pre, post, weight)
        '''
        temp_pre = []
        temp_s = []
        temp_weights = []
        
        temp_pre.append(self.sst_layer.rate_now)
        temp_s.append(self.S_sst_to_ndf)
        temp_weights.append(self.W_sst_to_ndf)
        # temp.append((self.sst_layer.rate_now, self.ndf_layer.rate_now, self.W_sst_to_ndf))
        
        temp_pre.append(self.integrator_layer.rate_now)
        temp_s.append(self.S_integrator_to_ndf)
        temp_weights.append(self.W_integrator_to_ndf)
        # temp.append((self.integrator_layer.rate_now, self.ndf_layer.rate_now, self.W_integrator_to_ndf))
        return temp_pre, temp_s, temp_weights
    
    # 
    def get_connections_to_dend(self):
        ''' Get the connections to the dendrites as an arrays of tuple (pre, post, weight)
        '''
        temp_inh = []
        temp_exc = []
        temp_inh.append((self.sst_layer.rate_now, self.S_sst_to_dend, self.W_sst_to_dend))
        temp_exc.append((self.integrator_layer.rate_now, self.S_integrator_to_dend, self.W_integrator_to_dend))
        temp_inh.append((self.ndf_layer.rate_now, self.S_ndf_to_dend, self.W_ndf_to_dend))
        return temp_inh,temp_exc
    
    #endregion
       
    #region Update of synaptic conductances
    
    
    def update_GABA_synaptic_conductances(self, PARAMS,weights, presynaptic_rates, u = 0.0, d = 0.0,facilitation = False, depression = False):
        ''' 
        Compute the derivative of the GABA synaptic conductances and update the conductances
        
        :param facilitation:  Boolean to indicate if the synapse is subject to facilitation
        :param depression: Boolean to indicate if the synapse is subject to depression
        '''
        
        #region Compute the GABA conductance derivatives
        
        
        if facilitation:
            temp = PARAMS['U_facilitation'] * (1.0 - u) * presynaptic_rates + (PARAMS['U_facilitation'] - u)/PARAMS['tau_facilitation']
            u = bf.euler_update(u, temp, PARAMS['dt'])
            
            temp_ds = -weights/PARAMS['tau_GABA'] + PARAMS['gamma_GABA'] * PARAMS['mult_gaba_facilitation'] * u * presynaptic_rates
            return bf.euler_update(weights, temp_ds, PARAMS['dt']), u
            
        elif depression:
            temp =   (1.0 - d)/PARAMS['tau_depression'] - d * presynaptic_rates * (1.0 - PARAMS['fD'])    
            d = bf.euler_update(d, temp, PARAMS['dt'])
            
            temp_ds = - weights/PARAMS['tau_GABA'] + PARAMS['gamma_GABA'] * PARAMS['mult_gaba_depression'] * d * presynaptic_rates
            return bf.euler_update(weights, temp_ds, PARAMS['dt']), d
            
        else:
            temp_ds = - weights/PARAMS['tau_GABA'] + PARAMS['gamma_GABA'] * presynaptic_rates
            return bf.euler_update(weights, temp_ds, PARAMS['dt'])
        
        #endregion
    
    
    def update_NMDA_synaptic_conductances(self, PARAMS, weights, presynaptic_rates, u = 0.0, d = 0.0, facilitation = False, depression = False):
        ''' 
        Compute the derivative of the NMDA synaptic conductances and update the conductances
        
        :param facilitation:  Boolean to indicate if the synapse is subject to facilitation
        :param depression: Boolean to indicate if the synapse is subject to depression
        '''
        
        #region Compute the NMDA conductance derivatives
        
        if facilitation:
            temp = PARAMS['U_facilitation'] * (1.0 - u) * presynaptic_rates + (PARAMS['U_facilitation'] - u)/PARAMS['tau_facilitation']
            u = bf.euler_update(u, temp, PARAMS['dt'])
            
            temp_ds = -weights/PARAMS['tau_NMDA'] + PARAMS['gamma_NMDA'] * PARAMS['mult_nmda_facilitation'] * u * presynaptic_rates * (1.0 - weights)
            return bf.euler_update(weights, temp_ds, PARAMS['dt']), u
            
        elif depression:
            temp =   (1.0 - d)/PARAMS['tau_depression'] - d * presynaptic_rates * (1.0 - PARAMS['fD'])    
            d = bf.euler_update(d, temp, PARAMS['dt'])
            
            temp_ds = - weights/PARAMS['tau_NMDA'] + PARAMS['gamma_NMDA'] * PARAMS['mult_nmda_depression'] * d * presynaptic_rates * (1.0 - weights)
            return bf.euler_update(weights, temp_ds, PARAMS['dt']), d
            
        else:
            temp_ds = - weights/PARAMS['tau_NMDA'] + PARAMS['gamma_NMDA']  * (1.0 - weights  ) * presynaptic_rates
            
            return bf.euler_update(weights, temp_ds, PARAMS['dt'])
            
        
        #endregion
    #TODO update the functions with facilitation and depression
     # TODO there is the issue of the different time constant for dendrites
   
    def update_conductances_from_PV(self):
        ''' Update all the conductances from the PV layer
        '''
        
        self.S_pv_to_pc = self.update_GABA_synaptic_conductances(self.PARAMS_PV, self.S_pv_to_pc, self.pv_layer.rate_now)
        
        self.S_pv_to_pv, self.U_pv_to_pv = self.update_GABA_synaptic_conductances(self.PARAMS_PV, self.S_pv_to_pv, self.pv_layer.rate_now, depression = self.PARAMS_PV['depression_to_pv'], d = self.U_pv_to_pv)
        
    def update_conductances_from_VIP(self):
        ''' Update all the conductances from the VIP layer
        '''
        self.S_vip_to_sst, self.U_vip_to_sst = self.update_GABA_synaptic_conductances(self.PARAMS_VIP, self.S_vip_to_sst, self.vip_layer.rate_now, facilitation = self.PARAMS_VIP['facilitation_to_sst'], u = self.U_vip_to_sst)
    
    def update_conductances_from_SST(self):
        ''' Update all the conductances from the SST layer
        '''
        self.S_sst_to_dend = self.update_GABA_synaptic_conductances(self.PARAMS_SST, self.S_sst_to_dend, self.sst_layer.rate_now)
        
        self.S_sst_to_pv = self.update_GABA_synaptic_conductances(self.PARAMS_SST, self.S_sst_to_pv, self.sst_layer.rate_now)
        
        self.S_sst_to_vip, self.U_sst_to_vip = self.update_GABA_synaptic_conductances(self.PARAMS_SST, self.S_sst_to_vip, self.sst_layer.rate_now, facilitation = self.PARAMS_SST['facilitation_to_vip'], u = self.U_sst_to_vip)
        
        self.S_sst_to_ndf = self.update_GABA_synaptic_conductances(self.PARAMS_SST, self.S_sst_to_ndf, self.sst_layer.rate_now)
    
    def update_conductances_from_Integrator(self):
        ''' Update all the conductances from the Integrator layer
        '''
        self.S_integrator_to_dend, self.U_integrator_to_dend = self.update_NMDA_synaptic_conductances(self.PARAMS_Integrator, self.S_integrator_to_dend, self.integrator_layer.rate_now, depression = self.PARAMS_Integrator['depression_to_dend'], d = self.U_integrator_to_dend)
        
        self.S_integrator_to_pv, self.U_integrator_to_pv = self.update_NMDA_synaptic_conductances(self.PARAMS_Integrator, self.S_integrator_to_pv, self.integrator_layer.rate_now, depression = self.PARAMS_Integrator['depression_to_pv'], d = self.U_integrator_to_pv)
        
        self.S_integrator_to_vip, self.U_integrator_to_vip = self.update_NMDA_synaptic_conductances(self.PARAMS_Integrator, self.S_integrator_to_vip, self.integrator_layer.rate_now, depression = self.PARAMS_Integrator['depression_to_vip'], d = self.U_integrator_to_vip)
        
        self.S_integrator_to_ndf = self.update_NMDA_synaptic_conductances(self.PARAMS_Integrator, self.S_integrator_to_ndf, self.integrator_layer.rate_now)
        
        self.S_integrator_to_integrator = self.update_NMDA_synaptic_conductances(self.PARAMS_Integrator, self.S_integrator_to_integrator, self.integrator_layer.rate_now)
        
    def update_conductances_from_NDF(self):
        ''' Update all the conductances from the NDF layer
        '''
        self.S_ndf_to_dend, self.U_ndf_to_dend = self.update_GABA_synaptic_conductances(self.PARAMS_NDF, self.S_ndf_to_dend, self.ndf_layer.rate_now, depression = self.PARAMS_NDF['depression_to_dend'], d = self.U_ndf_to_dend)
        
        self.S_ndf_to_pv = self.update_GABA_synaptic_conductances(self.PARAMS_NDF, self.S_ndf_to_pv, self.ndf_layer.rate_now)
        
    def update_conductances_from_PC(self):
        ''' Update all the conductances from the PC layer
        '''
        self.S_pc_to_pc = self.update_NMDA_synaptic_conductances(self.PARAMS_PC, self.S_pc_to_pc, self.pc_layer.rate_now)
        
        self.S_pc_to_pv, self.U_pc_to_pv = self.update_NMDA_synaptic_conductances(self.PARAMS_PC, self.S_pc_to_pv, self.pc_layer.rate_now, depression = self.PARAMS_PC['depression_to_pv'] , d = self.U_pc_to_pv)
        
        self.S_pc_to_vip, self.U_pc_to_vip = self.update_NMDA_synaptic_conductances(self.PARAMS_PC, self.S_pc_to_vip, self.pc_layer.rate_now, facilitation = self.PARAMS_PC['facilitation_to_vip'] , u = self.U_pc_to_vip)
        
        self.S_pc_to_sst, self.U_pc_to_sst = self.update_NMDA_synaptic_conductances(self.PARAMS_PC, self.S_pc_to_sst, self.pc_layer.rate_now, facilitation = self.PARAMS_PC['facilitation_to_sst'] , u = self.U_pc_to_sst)
        
        self.S_pc_to_integrator = self.update_NMDA_synaptic_conductances(self.PARAMS_PC, self.S_pc_to_integrator, self.pc_layer.rate_now)
        
    def update_conductances(self):
        ''' Update all the conductances
        '''
        self.update_conductances_from_PV()
        self.update_conductances_from_VIP()
        self.update_conductances_from_SST()
        self.update_conductances_from_Integrator()
        self.update_conductances_from_NDF()
        self.update_conductances_from_PC()
        
    #endregion        
        
        
    #region Update the weights 
    def update_weights_from_NDF(self):
        ''' Update the weights from the NDF layer using the Hebb with decay rule
        '''
    
        temp = bf.Hebb_with_decay(-self.W_ndf_to_dend, self.ndf_layer.get_rate(), self.pc_layer.get_rate(), self.PARAMS_NDF['gamma'], self.PARAMS_NDF['lambda_dec'], wmax = self.PARAMS_NDF['wmax'])
        self.W_ndf_to_dend = self.W_ndf_to_dend - self.PARAMS_NDF['dt']*temp*self.PARAMS_NDF['bool_plasticity']
    
        temp = bf.Hebb_with_decay(self.W_integrator_to_ndf, self.integrator_layer.get_rate(), self.ndf_layer.get_rate(), self.PARAMS_Integrator['gamma'], self.PARAMS_Integrator['lambda_dec'], wmax = self.PARAMS_Integrator['wmax'])
        self.W_integrator_to_ndf = self.W_integrator_to_ndf + self.PARAMS_Integrator['dt']*temp*self.PARAMS_Integrator['bool_plasticity']
    
  
    
    
    #endregion
   
    
    #region Update the currents to all layers of the network
    def update_current_to_Integrator(self):
            ''' Compute the current towards the integrator layer
            '''
            
            self.integrator_layer.input_current[:] = 0.0
            temp_pre, temp_s, temp_weights = self.get_connections_to_Integrator()
            for i in range(len(temp_pre)):
                self.integrator_layer.input_current += np.sum(temp_weights[i] * temp_s[i],axis=1)
            
    def update_current_to_PC(self):
        ''' Compute the current towards the PC layer
        '''
        
        self.pc_layer.input_current[:] = 0.0
        temp_pre, temp_s, temp_weights = self.get_connections_to_PC()
        for i in range(len(temp_pre)):
            self.pc_layer.input_current += np.sum(temp_weights[i] * temp_s[i],axis=1)
            
        self.pc_layer.Iexc_to_dend[:] = 0.0
        self.pc_layer.Iinh_to_dend[:] = 0.0
        temp_inh, temp_exc = self.get_connections_to_dend()
        for i in range(len(temp_inh)):
            self.pc_layer.Iinh_to_dend += np.sum(temp_inh[i][2] * temp_inh[i][1], axis =1)
        for i in range(len(temp_exc)):
            self.pc_layer.Iexc_to_dend += np.sum(temp_exc[i][2] * temp_exc[i][1], axis =1)
        
    def update_current_to_PV(self):
        ''' Compute the current towards the PV layer
        '''
        
        self.pv_layer.input_current[:] = 0.0
        temp_pre, temp_s, temp_weights = self.get_connections_to_PV()
        for i in range(len(temp_pre)):
            self.pv_layer.input_current += np.sum(temp_weights[i] * temp_s[i],axis=1)
            
    def update_current_to_SST(self):
        ''' Compute the current towards the SST layer
        '''
        
        self.sst_layer.input_current[:] = 0.0
        temp_pre, temp_s, temp_weights = self.get_connections_to_SST()
        for i in range(len(temp_pre)):
            self.sst_layer.input_current += np.sum(temp_weights[i] * temp_s[i],axis=1)
            
    def update_current_to_VIP(self):
        ''' Compute the current towards the VIP layer
        '''
        
        self.vip_layer.input_current[:] = 0.0
        temp_pre, temp_s, temp_weights = self.get_connections_to_VIP()
        for i in range(len(temp_pre)):
            self.vip_layer.input_current += np.sum(temp_weights[i] * temp_s[i],axis=1)
            
    def update_current_to_NDF(self):
        ''' Compute the current towards the NDF layer
        '''
        
        self.ndf_layer.input_current[:] = 0.0
        temp_pre, temp_s, temp_weights = self.get_connections_to_NDF()
        for i in range(len(temp_pre)):
            self.ndf_layer.input_current += np.sum(temp_weights[i] * temp_s[i],axis=1)
            
    
    def update_all_currents(self):
        ''' Update all currents in the network
        '''
        
        self.update_conductances()
        
        self.update_current_to_Integrator()
        self.update_current_to_PC()
        self.update_current_to_PV()
        self.update_current_to_SST()
        self.update_current_to_VIP()
        self.update_current_to_NDF()
    #endregion
    
    #region Update the rates of all layers of the network
    def time_step(self):
        ''' Time step of the network
        '''
        
        # update the currents of all layers
        self.update_all_currents()
        
        # update the rates of all layers
        self.integrator_layer.dr_dt(self.PARAMS_Simulation['dt'])
        self.pc_layer.dr_dt(self.PARAMS_Simulation['dt'])
        self.pv_layer.dr_dt(self.PARAMS_Simulation['dt'])
        self.sst_layer.dr_dt(self.PARAMS_Simulation['dt'])
        self.vip_layer.dr_dt(self.PARAMS_Simulation['dt'])
        self.ndf_layer.dr_dt(self.PARAMS_Simulation['dt'])
        self.update_weights_from_NDF() # update the weights from the NDF layer
        
    #TODO adding external inputs and so on 
    def full_dynamics(self,stim):
        ''' Full dynamics of the network. Until t_total with time_step dt
        '''
        
        for i in range(self.PARAMS_Simulation['N_t']):
            self.pc_layer.Iext[:] = stim[:,i]
         
            self.time_step()
            self.save_all_layers(i)
            
            
            
    def save_all_layers(self,i):
        ''' Save all layers of the network
        '''
        
        self.save_pc_layer[:,i] = self.pc_layer.get_rate()
        self.save_pv_layer[:,i] = self.pv_layer.get_rate()
        self.save_sst_layer[:,i] = self.sst_layer.get_rate()
        self.save_vip_layer[:,i] = self.vip_layer.get_rate()
        self.save_ndf_layer[:,i] = self.ndf_layer.get_rate()
        self.save_integrator_layer[:,i] = self.integrator_layer.get_rate()
        
        # save the weights received by the neuron selective to it
        self.save_W_dev_integrator_to_ndf[:,i] = self.W_integrator_to_ndf[self.PARAMS_Stimulus['dev_id']]
        self.save_W_standard_integrator_to_ndf[:,i] = self.W_integrator_to_ndf[self.PARAMS_Stimulus['std_id']]
        
        self.save_W_dev_integrator_to_dend[:,i] = self.W_integrator_to_dend[self.PARAMS_Stimulus['dev_id']]
        self.save_W_standard_integrator_to_dend[:,i] = self.W_integrator_to_dend[self.PARAMS_Stimulus['std_id']]
        
        self.save_W_dev_ndf_to_dend[:,i] = self.W_ndf_to_dend[self.PARAMS_Stimulus['dev_id']]
        self.save_W_standard_ndf_to_dend[:,i] = self.W_ndf_to_dend[self.PARAMS_Stimulus['std_id']]
        
        
        
        
    def write_firing_rates_in_file(self, name):
        ''' Write all the savedfiring rates in a pickle file
        '''
        
        temp = bf.compute_list_time_input(self.PARAMS_Simulation['dt'],self.PARAMS_Simulation['t_total'], self.PARAMS_Stimulus['Tinter'], self.PARAMS_Stimulus['Tresting'], self.PARAMS_Stimulus['Tstim'])
        
        save_dict = {
            'pc_firing_rates': self.save_pc_layer[:,temp],
            'sst_firing_rates': self.save_sst_layer[:,temp],
            'vip_firing_rates': self.save_vip_layer[:,temp],
            'ndf_firing_rates': self.save_ndf_layer[:,temp],
            'pv_firing_rates': self.save_pv_layer[:,temp],
            'integrator_firing_rates': self.save_integrator_layer[:,temp],
            'params_all': self.PARAMS_ALL
        }
        
        # Store data (serialize)    
        with open(name, 'w') as handle:
            json.dump(save_dict, handle, cls = NumpyEncoder)
        
    def compute_mean_firing_rate(self,stim = None, delay = 0.1):
        ''' '''
        if self.PARAMS_Stimulus['type'] not in ['multi_control_paradigm','probabilistic_MMN','deterministic_MMN']:
            raise ValueError('The function compute_mean_firing_rate is only defined for multi_control_paradigm')
        elif self.PARAMS_Stimulus['type'] == 'multi_control_paradigm':
            temp = bf.compute_list_time_input(self.PARAMS_Simulation['dt'],self.PARAMS_Simulation['t_total'], self.PARAMS_Stimulus['Tinter'], self.PARAMS_Stimulus['Tresting'], self.PARAMS_Stimulus['Tstim'], delay=delay)
            
            mean_temp = np.zeros(len(temp))
            mean_sum = np.zeros(len(temp))
            for i in range(len(temp)):
                mean_temp[i] = np.max(self.save_pc_layer[:,temp[i]])

                mean_sum[i] = np.sum(self.save_pc_layer[:,temp[i]])
            
            return np.mean(mean_temp), np.mean(mean_sum)
        
        elif self.PARAMS_Stimulus['type'] == 'probabilistic_MMN':
            std_id = self.PARAMS_Stimulus['std_id']
            dev_id = self.PARAMS_Stimulus['dev_id']
            
            temp = bf.compute_list_time_input(self.PARAMS_Simulation['dt'],self.PARAMS_Simulation['t_total'], self.PARAMS_Stimulus['Tinter'], self.PARAMS_Stimulus['Tresting'], self.PARAMS_Stimulus['Tstim'], delay = delay)
            
            is_std = [i for i in range(len(temp)) if np.argmax(self.save_pc_layer[:,temp[i]]) == std_id]
            is_dev = [i for i in range(len(temp)) if np.argmax(self.save_pc_layer[:,temp[i]]) == dev_id]
            
            mean_temp_std = np.zeros(len(is_std))
            mean_temp_dev = np.zeros(len(is_dev))
            
            mean_sum_std = np.zeros(len(is_std))
            mean_sum_dev = np.zeros(len(is_dev))
            
            for i in range(len(is_std)):
                mean_temp_std[i] = np.max(self.save_pc_layer[:,temp[is_std[i]]])
                mean_sum_std[i] = np.sum(self.save_pc_layer[:,temp[is_std[i]]])
                
            for i in range(len(is_dev)):
                mean_temp_dev[i] = np.max(self.save_pc_layer[:,temp[is_dev[i]]])
                mean_sum_dev[i] = np.sum(self.save_pc_layer[:,temp[is_dev[i]]])
                
            return np.mean(mean_temp_std), np.mean(mean_temp_dev), np.mean(mean_sum_std), np.mean(mean_sum_dev)
        
        elif self.PARAMS_Stimulus['type'] == 'deterministic_MMN':
           
            temp = bf.compute_list_time_input(self.PARAMS_Simulation['dt'],self.PARAMS_Simulation['t_total'], self.PARAMS_Stimulus['Tinter'], self.PARAMS_Stimulus['Tresting'], self.PARAMS_Stimulus['Tstim'])
        
            save_dict = {
            'pc_firing_rates': self.save_pc_layer[:,temp],
            'sst_firing_rates': self.save_sst_layer[:,temp],
            'vip_firing_rates': self.save_vip_layer[:,temp],
            'ndf_firing_rates': self.save_ndf_layer[:,temp],
            'pv_firing_rates': self.save_pv_layer[:,temp],
            'integrator_firing_rates': self.save_integrator_layer[:,temp],
            'params_all': self.PARAMS_ALL
        }
            
            
            Tinter = self.PARAMS_Stimulus['Tinter']
            Tstim = self.PARAMS_Stimulus['Tstim']
            nbr_rep_std = self.PARAMS_Stimulus['nbr_rep_std']
            nbr_rep_dev = self.PARAMS_Stimulus['nbr_rep_dev']
            std_id = self.PARAMS_Stimulus['std_id']
            dev_id = self.PARAMS_Stimulus['dev_id']
            Tresting = self.PARAMS_Stimulus['Tresting']
            strength_std = self.PARAMS_Stimulus['strength_std']
            strength_dev = self.PARAMS_Stimulus['strength_dev']
            
            sigma_std = self.PARAMS_Stimulus['sigma_std']
            sigma_dev = self.PARAMS_Stimulus['sigma_dev']
                
            fr_std = np.sum(self.save_pc_layer[:,temp[nbr_rep_std-1]])    
            fr_dev = np.sum(self.save_pc_layer[:,temp[nbr_rep_std]])    

                
            
            return fr_std, fr_dev
            
            
            
            
    #endregion
    
            
#region Stimulus class
# This region will define the stimulus class and the different methods that can apply to it

class Stimulus:
    ''' Class that defines the stimulus
    '''
    
    def __init__(self, PARAMS_Stimulus, PARAMS_Simulation, Ncells):
        ''' Initialize the stimulus
        and define all the parameters
        '''
        
        self.PARAMS_Stimulus = PARAMS_Stimulus
        self.PARAMS_Simulation = PARAMS_Simulation
        # self.PARAMS_Simulation['t_total'] = (self.PARAMS_Stimulus['Tinter'] + self.PARAMS_Stimulus['Tstim'])*(self.PARAMS_Stimulus['nbr_rep_std']+ self.PARAMS_Stimulus['nbr_rep_dev']) + self.PARAMS_Stimulus['Tresting']
        
        self.PARAMS_Stimulus['N_t'] = int(self.PARAMS_Simulation['t_total']/self.PARAMS_Simulation['dt'])
        self.stimulus = np.zeros((Ncells, self.PARAMS_Stimulus['N_t'])) #full matrix of stimuli for every neurons
        self.Ncells = Ncells
        
        if self.PARAMS_Stimulus['type'] == 'deterministic_MMN':
            self.generate_deterministic_MMN()
            
        if self.PARAMS_Stimulus['type'] == 'oscillation_pattern_MMN':
            self.generate_pattern_MMN()
            
        if self.PARAMS_Stimulus['type'] == 'multi_control_paradigm':
            self.generate_multi_control()
            
        if self.PARAMS_Stimulus['type'] == 'probabilistic_MMN':
            self.generate_probabilistic_MMN()
            
#region define the deterministic MMN stimulus        
    def generate_deterministic_MMN(self):
        ''' This function generates  a stimulus that will allow to test for MMN signal
        Let's first assume that the signal is deterministic in terms of oscillations 
        '''
        Tinter = self.PARAMS_Stimulus['Tinter']
        Tstim = self.PARAMS_Stimulus['Tstim']
        nbr_rep_std = self.PARAMS_Stimulus['nbr_rep_std']
        nbr_rep_dev = self.PARAMS_Stimulus['nbr_rep_dev']
        std_id = self.PARAMS_Stimulus['std_id']
        dev_id = self.PARAMS_Stimulus['dev_id']
        Tresting = self.PARAMS_Stimulus['Tresting']
        strength_std = self.PARAMS_Stimulus['strength_std']
        strength_dev = self.PARAMS_Stimulus['strength_dev']
        
        sigma_std = self.PARAMS_Stimulus['sigma_std']
        sigma_dev = self.PARAMS_Stimulus['sigma_dev']
        
        
        nbr_stim_std = 0
        nbr_stim_dev = -nbr_rep_std

        for t in range(self.PARAMS_Stimulus['N_t']):
            if t < Tinter/self.PARAMS_Simulation['dt'] + nbr_stim_std*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                self.stimulus[:,t] = 0.0
                
            elif t < (Tinter + Tstim)/self.PARAMS_Simulation['dt'] + nbr_stim_std*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                if nbr_stim_std < nbr_rep_std:
                    
                    self.stimulus[:,t] = bf.generate_ring_stim(std_id, sigma_std, self.Ncells, strength_std)
               
                else:  
                    if nbr_stim_dev < nbr_rep_dev:
    
                        self.stimulus[:,t] = bf.generate_ring_stim(dev_id, sigma_dev, self.Ncells, strength_dev)
                    
                    else:
                        self.stimulus[:,t] = bf.generate_ring_stim(std_id, sigma_std, self.Ncells, strength_std)

                    
            if t > (Tinter + Tstim)/self.PARAMS_Simulation['dt'] + nbr_stim_std*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                
                nbr_stim_std += 1
                nbr_stim_dev +=1
#endregion

    def generate_pattern_MMN(self):
        ''' This function generates a pattern stimulus ABABABAB that will allow to test the MMN in this context
        '''
        Tinter = self.PARAMS_Stimulus['Tinter']
        Tstim = self.PARAMS_Stimulus['Tstim']
        nbr_rep_std = self.PARAMS_Stimulus['nbr_rep_std']
        nbr_rep_dev = self.PARAMS_Stimulus['nbr_rep_dev']
        std_id = self.PARAMS_Stimulus['std_id']
        dev_id = self.PARAMS_Stimulus['dev_id']
        Tresting = self.PARAMS_Stimulus['Tresting']
        strength_std = self.PARAMS_Stimulus['strength_std']
        strength_dev = self.PARAMS_Stimulus['strength_dev']
        
        sigma_std = self.PARAMS_Stimulus['sigma_std']
        sigma_dev = self.PARAMS_Stimulus['sigma_dev']
        
        
        nbr_stim_std = 0
        nbr_stim_dev = -nbr_rep_std
        nbr_oscillation = 0
        # this time the stimulus is A_B_
        # the deviant is B_B
        # the trick is to have an oscillation varaible?
        for t in range(self.PARAMS_Stimulus['N_t']):
            
            if t < Tinter/self.PARAMS_Simulation['dt'] + nbr_stim_std*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                
                self.stimulus[:,t] = 0.0
                
            elif t < (Tinter + Tstim)/self.PARAMS_Simulation['dt'] + nbr_stim_std*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                
                if( nbr_oscillation == 0 ):
                    if nbr_stim_dev < nbr_rep_std:
                        self.stimulus[:,t] = bf.generate_ring_stim(std_id, sigma_std, self.Ncells, strength_std)
                    else:
                        self.stimulus[:,t] = bf.generate_ring_stim(dev_id, sigma_dev, self.Ncells, strength_dev)

                else:
                    self.stimulus[:,t] = bf.generate_ring_stim(dev_id, sigma_dev, self.Ncells, strength_dev)
                    
            if t > (Tinter + Tstim)/self.PARAMS_Simulation['dt'] + nbr_stim_std*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                    
                    if nbr_stim_dev > nbr_rep_std:
                        nbr_stim_dev = -nbr_rep_std
                        
                    nbr_stim_std += 1
                    nbr_stim_dev +=1
                    if( nbr_oscillation == 0):
                        
                        nbr_oscillation = 1
                    else:
                        nbr_oscillation = 0
                        
                        
    def generate_multi_control(self):
        ''' Generate the stimuli distribution in the case of the multi-control paradigm '''
        Tinter = self.PARAMS_Stimulus['Tinter']
        list_std = self.PARAMS_Stimulus['list_std']
        Tstim = self.PARAMS_Stimulus['Tstim']
        std_id = self.PARAMS_Stimulus['std_id']
        Tresting = self.PARAMS_Stimulus['Tresting']
        strength_std = self.PARAMS_Stimulus['strength_std']
        nbr_rep_std = self.PARAMS_Stimulus['nbr_rep_std']
        sigma_std = self.PARAMS_Stimulus['sigma_std']
        nbr_stim_std = 0
        
        i=0
        stim_idx = np.random.randint(0,len(list_std))
        for t in range(self.PARAMS_Stimulus['N_t']):
                if t < Tinter/self.PARAMS_Simulation['dt'] + i*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
           
                    self.stimulus[:,t] = 0.0
           
                elif t < (Tinter + Tstim)/self.PARAMS_Simulation['dt'] + i*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                    self.stimulus[:,t] = bf.generate_ring_stim(list_std[stim_idx], sigma_std, self.Ncells, strength_std)
                if t > (Tinter + Tstim)/self.PARAMS_Simulation['dt'] + i*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                    i += 1
                    # nbr_stim_std += 1
                    stim_idx = np.random.randint(0,len(list_std))
           
    def generate_probabilistic_MMN(self):
        ''' Generate the stimulus for the MMN task in a probabilistic way'''
        Tinter = self.PARAMS_Stimulus['Tinter']
        Tstim = self.PARAMS_Stimulus['Tstim']
        # nbr_rep_std = self.PARAMS_Stimulus['nbr_rep_std']
        # nbr_rep_dev = self.PARAMS_Stimulus['nbr_rep_dev']
        std_id = self.PARAMS_Stimulus['std_id']
        dev_id = self.PARAMS_Stimulus['dev_id']
        Tresting = self.PARAMS_Stimulus['Tresting']
        strength_std = self.PARAMS_Stimulus['strength_std']
        strength_dev = self.PARAMS_Stimulus['strength_dev']
        
        sigma_std = self.PARAMS_Stimulus['sigma_std']
        sigma_dev = self.PARAMS_Stimulus['sigma_dev']
        
        prob_std = self.PARAMS_Stimulus['prob_std']
        prob_dev = 1.0 - prob_std
        
        i=0
        stim_idx = np.random.uniform()
        
        
        for t in range(self.PARAMS_Stimulus['N_t']):
            
                if t < Tinter/self.PARAMS_Simulation['dt'] + i*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
           
                    self.stimulus[:,t] = 0.0
           
                elif t < (Tinter + Tstim)/self.PARAMS_Simulation['dt'] + i*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
        
                    if stim_idx < prob_std:
                    
                       self.stimulus[:,t] = bf.generate_ring_stim(std_id, sigma_std, self.Ncells, strength_std)
                    
                    else:
                        
                        self.stimulus[:,t] = bf.generate_ring_stim(dev_id, sigma_dev, self.Ncells, strength_dev)
                        
                if t > (Tinter + Tstim)/self.PARAMS_Simulation['dt'] + i*(Tstim+Tinter)/self.PARAMS_Simulation['dt'] + Tresting/self.PARAMS_Simulation['dt']:
                    
                    i += 1
                    # nbr_stim_std += 1
                    stim_idx = np.random.uniform()
        
                    
                    
#endregion
        
    
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)