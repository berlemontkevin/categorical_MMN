# This file will have all the initialization functions for the neural networks
import numpy as np
import basic_functions as bf

# dictionnary with all the weights
Dict_Weights = {
    'W_integrator_to_dend': None,
    'W_integrator_to_pv': None,
    'W_integrator_to_ndf': None,
    'W_integrator_to_vip': None,
    # weights for pc
    'W_pc_to_pc': None,
    'W_pc_to_pv': None,
    'W_pc_to_sst': None,
    'W_pc_to_vip': None,
    'W_pc_to_integrator': None,
    # weights for pv
    'W_pv_to_pc': None,
    'W_pv_to_pv': None,
    'W_pv_to_sst': None,
    # weights for sst
    'W_sst_to_dend': None,
    'W_sst_to_vip': None,
    'W_sst_to_pv': None,
    'W_sst_to_ndf': None,
    # weights for vip
    'W_vip_to_sst': None,
    # weights for ndf
    'W_ndf_to_dend': None
}

#TODO maybe keep a list of all then ame of the weights?

#region Initialization of the weights



def create_connectivity_from_SST(Dict_Weights, PARAMS_SST):
        '''create the weights matrices going from SST neurons
        '''
        
        Dict_Weights['W_sst_to_vip'] = np.zeros((PARAMS_SST['Ncells'], PARAMS_SST['Ncells']))
        np.fill_diagonal(Dict_Weights['W_sst_to_vip'], PARAMS_SST['weight_to_vip'])
        
        
        Dict_Weights['W_sst_to_dend'] = np.zeros((PARAMS_SST['Ncells'], PARAMS_SST['Ncells']))
        np.fill_diagonal(Dict_Weights['W_sst_to_dend'], PARAMS_SST['weight_to_dend'])
        
        
        Dict_Weights['W_sst_to_pv'] = np.zeros((PARAMS_SST['Ncells'], PARAMS_SST['Ncells']))
        np.fill_diagonal(Dict_Weights['W_sst_to_pv'], PARAMS_SST['weight_to_pv'])
        
        
        Dict_Weights['W_sst_to_ndf'] = np.zeros((PARAMS_SST['Ncells'], PARAMS_SST['Ncells']))
        np.fill_diagonal(Dict_Weights['W_sst_to_ndf'], PARAMS_SST['weight_to_ndf'])
        
def create_connectivity_from_VIP(Dict_Weights, PARAMS_VIP):
        '''create the weights matrices going from VIP neurons
        '''
        Dict_Weights['W_vip_to_sst'] = np.zeros((PARAMS_VIP['Ncells'], PARAMS_VIP['Ncells']))
        np.fill_diagonal(Dict_Weights['W_vip_to_sst'], PARAMS_VIP['weight_to_sst'])
        
def create_connectivity_from_NDF(Dict_Weights, PARAMS_NDF):
        '''create the weights matrices going from NDF neurons
        '''
        
        Dict_Weights['W_ndf_to_dend'] = np.zeros((PARAMS_NDF['Ncells'], PARAMS_NDF['Ncells']))
        np.fill_diagonal(Dict_Weights['W_ndf_to_dend'], PARAMS_NDF['weight_to_dend'])
        
def create_connectivity_from_PV(Dict_Weights, PARAMS_PV):
        '''create the weights matrices going from PV neurons
        '''
        
        Dict_Weights['W_pv_to_pc'] = np.zeros((PARAMS_PV['Ncells'], PARAMS_PV['Ncells']))
        np.fill_diagonal(Dict_Weights['W_pv_to_pc'], PARAMS_PV['weight_to_soma'])
        
        Dict_Weights['W_pv_to_pv'] = np.zeros((PARAMS_PV['Ncells'], PARAMS_PV['Ncells']))
        np.fill_diagonal(Dict_Weights['W_pv_to_pv'], PARAMS_PV['weight_to_pv'])
        
        
    
def create_connectivity_from_PC(Dict_Weights, PARAMS_PC):
        ''' create the weights matrices going from PC neurons
        '''
        
        Dict_Weights['W_pc_to_vip'] = np.zeros((PARAMS_PC['Ncells'], PARAMS_PC['Ncells']))
        np.fill_diagonal(Dict_Weights['W_pc_to_vip'], PARAMS_PC['weight_to_vip'])

        Dict_Weights['W_pc_to_sst'] = np.zeros((PARAMS_PC['Ncells'], PARAMS_PC['Ncells']))
        Dict_Weights['W_pc_to_sst'][:, :] = bf.create_ring_connectivity(PARAMS_PC['Ncells'],PARAMS_PC['Ncells'], PARAMS_PC['weight_to_sst'], PARAMS_PC['sigma_to_sst'])

        Dict_Weights['W_pc_to_pv'] = np.zeros((PARAMS_PC['Ncells'], PARAMS_PC['Ncells']))
        Dict_Weights['W_pc_to_pv'][:, :] = bf.create_ring_connectivity(PARAMS_PC['Ncells'], PARAMS_PC['Ncells'], PARAMS_PC['weight_to_pv'], PARAMS_PC['sigma_to_pv'])

        Dict_Weights['W_pc_to_integrator'] = np.zeros((PARAMS_PC['Ncells'], PARAMS_PC['Ncells']))
        Dict_Weights['W_pc_to_integrator'][:, :] = bf.create_ring_connectivity(PARAMS_PC['Ncells'], PARAMS_PC['Ncells'], PARAMS_PC['weight_to_integrator'], PARAMS_PC['sigma_to_integrator'])
        
        Dict_Weights['W_pc_to_pc'] = np.zeros((PARAMS_PC['Ncells'], PARAMS_PC['Ncells']))

        Dict_Weights['W_pc_to_pc'][:, :] = bf.create_ring_connectivity(PARAMS_PC['Ncells'], PARAMS_PC['Ncells'], PARAMS_PC['weight_to_soma'], PARAMS_PC['sigma_to_pc'])

def create_connectivity_from_Integrator(Dict_Weights, PARAMS_INT):
        ''' create the weights matrices going from Integrator neurons
        '''
      
        Dict_Weights['W_integrator_to_dend'] = np.zeros((PARAMS_INT['Ncells'], PARAMS_INT['Ncells']))
        Dict_Weights['W_integrator_to_dend'][:, :] = bf.create_ring_connectivity(PARAMS_INT['Ncells'], PARAMS_INT['Ncells'], PARAMS_INT['weight_to_dend'],PARAMS_INT['sigma_to_dend'])
        
        Dict_Weights['W_integrator_to_ndf'] = np.zeros((PARAMS_INT['Ncells'], PARAMS_INT['Ncells']))
        Dict_Weights['W_integrator_to_ndf'][:, :] = bf.create_ring_connectivity(PARAMS_INT['Ncells'],PARAMS_INT['Ncells'], PARAMS_INT['weight_to_ndf'], PARAMS_INT['sigma_to_ndf'])
        
        Dict_Weights['W_integrator_to_pv'] = np.zeros((PARAMS_INT['Ncells'], PARAMS_INT['Ncells']))
        Dict_Weights['W_integrator_to_pv'][:, :] = bf.create_ring_connectivity(PARAMS_INT['Ncells'],PARAMS_INT['Ncells'], PARAMS_INT['weight_to_pv'], PARAMS_INT['sigma_to_pv'])
        
        Dict_Weights['W_integrator_to_vip'] = np.zeros((PARAMS_INT['Ncells'], PARAMS_INT['Ncells']))
        Dict_Weights['W_integrator_to_vip'][:, :] = bf.create_ring_connectivity(PARAMS_INT['Ncells'],PARAMS_INT['Ncells'], PARAMS_INT['weight_to_vip'], PARAMS_INT['sigma_to_vip'])
        
    

def create_internal_connectivity(Dict_Weights, PARAMS_PC, PARAMS_INT, PARAMS_SST, PARAMS_PV, PARAMS_NDF, PARAMS_VIP, **kwargs):
        ''' for now the only broad connections are E to E,PV and SST, as well as Int to Int, N, PV, VIP, E
        '''
        
        create_connectivity_from_PC(Dict_Weights, PARAMS_PC)
        create_connectivity_from_Integrator(Dict_Weights, PARAMS_INT)
        create_connectivity_from_SST(Dict_Weights, PARAMS_SST)
        create_connectivity_from_PV(Dict_Weights, PARAMS_PV)
        create_connectivity_from_VIP(Dict_Weights, PARAMS_VIP)
        create_connectivity_from_NDF(Dict_Weights, PARAMS_NDF)
        
        
#endregion
       
       
    # #region get the connections towards all types of neurons
    
    # # 
    # def get_connections_to_Integrator(self):
    #     ''' Get the connections to the Integrator as an arrays of tuple (pre, post, weight)
    #     '''
    #     temp_pre = []
    #     temp_post = []
    #     temp_weights = []
        
    #     temp_pre.append(self.integrator_layer.rate_now)
    #     temp_post.append(self.integrator_layer.rate_now)
    #     temp_weights.append(self.W_pc_to_integrator)
    #     # temp.append((self.pc_layer.rate_now, self.integrator_layer.rate_now, self.W_pc_to_integrator))
    #     return temp_pre, temp_post, temp_weights
    
    # # 
    # def get_connections_to_PC(self):
    #     ''' Get the connections to the PC as an arrays of tuple (pre, post, weight)
    #     '''
    #     temp_pre = []
    #     temp_post = []
    #     temp_weights = []
        
    #     temp_pre.append(self.pv_layer.rate_now)
    #     temp_post.append(self.pc_layer.rate_now)
    #     temp_weights.append(self.W_pv_to_pc)
        
    #     temp_pre.append(self.pc_layer.rate_now)
    #     temp_post.append(self.pc_layer.rate_now)
    #     temp_weights.append(self.W_pc_to_pc)
        
    #     # temp.append((self.pv_layer.rate_now, self.pc_layer.rate_now, self.W_pv_to_pc))
    #     # temp.append((self.pc_layer.rate_now, self.pc_layer.rate_now, self.W_pc_to_pc))
    #     return temp_pre, temp_post, temp_weights
    
    # # 
    # def get_connections_to_PV(self):
    #     ''' Get the connections to the PV as an arrays of tuple (pre, post, weight)
    #     '''
        
    #     temp_pre = []
    #     temp_post = []
    #     temp_weights = []
        
    #     temp_pre.append(self.sst_layer.rate_now)
    #     temp_post.append(self.pv_layer.rate_now)
    #     temp_weights.append(self.W_sst_to_pv)
    #     # temp.append((self.sst_layer.rate_now, self.pv_layer.rate_now, self.W_sst_to_pv))
        
    #     temp_pre.append(self.pv_layer.rate_now)
    #     temp_post.append(self.pv_layer.rate_now)
    #     temp_weights.append(self.W_pv_to_pv)
    #     # temp.append((self.pv_layer.rate_now, self.pv_layer.rate_now, self.W_pv_to_pv))
        
    #     temp_pre.append(self.pc_layer.rate_now)
    #     temp_post.append(self.pv_layer.rate_now)
    #     temp_weights.append(self.W_pc_to_pv)
    #     # temp.append((self.pc_layer.rate_now, self.pv_layer.rate_now, self.W_pc_to_pv))
        
    #     temp_pre.append(self.integrator_layer.rate_now)
    #     temp_post.append(self.pv_layer.rate_now)
    #     temp_weights.append(self.W_integrator_to_pv)
    #     # temp.append((self.integrator_layer.rate_now, self.pv_layer.rate_now, self.W_integrator_to_pv))
    #     return temp_pre, temp_post, temp_weights
      
    # #   
    # def get_connections_to_SST(self):
    #     ''' Get the connections to the SST as an arrays of tuple (pre, post, weight)
    #     '''
    #     temp_pre = []
    #     temp_post = []
    #     temp_weights = []
        
    #     temp_pre.append(self.vip_layer.rate_now)
    #     temp_post.append(self.sst_layer.rate_now)
    #     temp_weights.append(self.W_vip_to_sst)
    #     # temp.append((self.vip_layer.rate_now, self.sst_layer.rate_now, self.W_vip_to_sst))
        
    #     temp_pre.append(self.pc_layer.rate_now)
    #     temp_post.append(self.sst_layer.rate_now)
    #     temp_weights.append(self.W_pc_to_sst)
    #     # temp.append((self.pc_layer.rate_now, self.sst_layer.rate_now, self.W_pc_to_sst))
    #     return temp_pre, temp_post, temp_weights
    
    # # 
    # def get_connections_to_VIP(self):
    #     ''' Get the connections to the VIP as an arrays of tuple (pre, post, weight)
    #     '''
    #     temp_pre = []
    #     temp_post = []
    #     temp_weights = []
        
    #     temp_pre.append(self.sst_layer.rate_now)
    #     temp_post.append(self.vip_layer.rate_now)
    #     temp_weights.append(self.W_sst_to_vip)
    #     # temp.append((self.sst_layer.rate_now, self.vip_layer.rate_now, self.W_sst_to_vip))
        
    #     temp_pre.append(self.pc_layer.rate_now)
    #     temp_post.append(self.vip_layer.rate_now)
    #     temp_weights.append(self.W_pc_to_vip)
    #     # temp.append((self.pc_layer.rate_now, self.vip_layer.rate_now, self.W_pc_to_vip))
        
    #     temp_pre.append(self.integrator_layer.rate_now)
    #     temp_post.append(self.vip_layer.rate_now)
    #     temp_weights.append(self.W_integrator_to_vip)
    #     # temp.append((self.integrator_layer.rate_now, self.vip_layer.rate_now, self.W_integrator_to_vip))
    #     return temp_pre, temp_post, temp_weights
    
    # # 
    # def get_connections_to_NDF(self):
    #     ''' Get the connections to the NDF as an arrays of tuple (pre, post, weight)
    #     '''
    #     temp_pre = []
    #     temp_post = []
    #     temp_weights = []
        
    #     temp_pre.append(self.sst_layer.rate_now)
    #     temp_post.append(self.ndf_layer.rate_now)
    #     temp_weights.append(self.W_sst_to_ndf)
    #     # temp.append((self.sst_layer.rate_now, self.ndf_layer.rate_now, self.W_sst_to_ndf))
        
    #     temp_pre.append(self.integrator_layer.rate_now)
    #     temp_post.append(self.ndf_layer.rate_now)
    #     temp_weights.append(self.W_integrator_to_ndf)
    #     # temp.append((self.integrator_layer.rate_now, self.ndf_layer.rate_now, self.W_integrator_to_ndf))
    #     return temp_pre, temp_post, temp_weights
    
    # # 
    # def get_connections_to_dend(self):
    #     ''' Get the connections to the dendrites as an arrays of tuple (pre, post, weight)
    #     '''
    #     temp_inh = []
    #     temp_exc = []
    #     temp_inh.append((self.sst_layer.rate_now, self.pc_layer.rate_now, self.W_sst_to_dend))
    #     temp_exc.append((self.integrator_layer.rate_now, self.pc_layer.rate_now, self.W_integrator_to_dend))
    #     temp_inh.append((self.ndf_layer.rate_now, self.pc_layer.rate_now, self.W_ndf_to_dend))
    #     return temp_inh,temp_exc
    
    # #endregion