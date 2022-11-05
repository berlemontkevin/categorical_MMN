import os
import torch
from torch import nn
import logging


logging.basicConfig(filename='sample.log', level=logging.INFO)


# def my_update(current_state):

#     # update the state


#     return next_state

# def fI_curve():
#     # compute the fI curve
#     return


# could one of the idea to just have a huge matrix of all the possible combinations of the parameters

class EngelModel(nn.Module):
    def __init__(self, input_size, association_size):
        super(EngelModel, self).__init__()
        self.input_size = input_size
        self.association_size = association_size
        self.decision_size = 2

        self.tauS_sensory = 0.06
        self.gammaS_sensory = 0.641
        
        
        self.a = 270 # Hz na-1
        self.b = 108.0 #Hz
        self.d = 0.154 #s

        # synaptic connections from all the layers
        self.W_sens_to_assoc = torch.nn.Parameter(torch.randn(input_size, association_size), requires_grad=False)
    
        self.W_assoc_to_decision = torch.nn.Parameter(torch.randn(association_size, self.decision_size), requires_grad=False)
        self.W_decision_to_assoc = torch.nn.Parameter(torch.randn(self.decision_size, association_size), requires_grad=False)

        #recurrent synaptic connections
        self.W_sens_to_sens = torch.randn(input_size, input_size)
        self.W_assoc_to_assoc = torch.randn(association_size, association_size)
        self.W_decision_to_decision = torch.randn(self.decision_size, self.decision_size)


    def fIcurve(self, I, a, b, d):
        # compute the fI curve

        return (a*I - b)/(1.0 + torch.exp(-d*(a*I-b)))


    def my_update_s(self,s,r):
        # update the state
        
        return s + self.dt*(-s/self.tauS_sensory + self.gammaS_sensory*r*(1.0-s))


    def ring_connectivity(self, Jminus, Jplus, sigma,  W):
        # this function will compute the weights of a guassian ring profile connectivity
        # the weights are symmetric
        num_neurons = W.shape[0]
        # transform an index into an angle
        angle = 2.0*torch.pi*torch.arange(num_neurons)/num_neurons
        # compute the distance between the angle of the neurons and neuron i modulo pi degrees
        for i in range(num_neurons):
            distance = torch.minimum(torch.abs(angle - angle[i]) , torch.abs(angle - angle[i] - 2.0*torch.pi))
            # compute the gaussian profile
            gaussian = torch.exp(-distance**2/(2.0*sigma**2))
            # compute the weights
            W[i,:] = Jplus*gaussian + Jminus



    def reset_weights(self):
        # function that should initialize the wiehgts the following
        # assoc to dec and dec to assoc are uniform between 0.25 and 0.75
        # sens to assoc is gaussian
        #decision circuit is a standard WW model
        # senspry circuit is a Gaussian profile
        # assocpry circuit is a Gaussian profile

        self.W_assoc_to_decision = torch.nn.Parameter(0.5*torch.rand(self.association_size, self.decision_size) + 0.25)
        self.W_decision_to_assoc = torch.nn.Parameter(0.5*torch.rand(self.decision_size, self.association_size) + 0.25)

        self.ring_connectivity(0.0, 1.0, 43.2/380.0*2.0*torch.pi, self.W_sens_to_assoc)
        self.ring_connectivity(-0.5,1.43,43.2/380.0*2.0*torch.pi, self.W_sens_to_sens)
        self.ring_connectivity(-10.0,-0.4,43.2/380.0*2.0*torch.pi, self.W_assoc_to_assoc)

        # reset the weights of the attractor/decision layer according to the Wong and Wang parameters
        self.W_decision_to_decision[0,0] = self.W_decision_to_decision[1,1] = 0.3725
        self.W_decision_to_decision[0,1] = self.W_decision_to_decision[1,0] = -0.1137





    def forward(self, I_input):
# TODO add the input as a temporal input (to generate noise in an efficient way)

        # init the state
        s_sensory = torch.zeros(self.input_size)
        s_association = torch.zeros(self.association_size)
        s_decision = torch.zeros(self.decision_size)

        r_sensory = torch.zeros(self.input_size)
        r_association = torch.zeros(self.association_size)
        r_decision = torch.zeros(self.decision_size)

        # compute temporary variables
        I_sens_to_assoc = torch.matmul(s_sensory, self.W_sens_to_assoc)/self.input_size
        I_assoc_to_decision = torch.matmul(s_association, self.W_assoc_to_decision)/self.association_size

        I_rec_sens = torch.matmul(s_sensory, self.W_sens_to_sens)/self.input_size
        I_rec_assoc = torch.matmul(s_association, self.W_assoc_to_assoc)/self.association_size
        I_rec_decision = torch.matmul(s_decision, self.W_decision_to_decision)/self.decision_size

        # compute the next state

        I_to_sensory = self.I_fixed_sensory + I_input + I_rec_sens
        r_sensory = self.fIcurve(I_to_sensory, self.a, self.b, self.d)
        s_sensory = self.my_update_s(s_sensory, r_sensory)

        I_to_association = self.I_fixed_association + I_sens_to_assoc + I_rec_assoc
        r_association = self.fIcurve(I_to_association, self.a, self.b, self.d)
        s_association = self.my_update_s(s_association, r_association)

        I_to_decision = self.I_fixed_association + I_assoc_to_decision + I_rec_decision
        r_decision = self.fIcurve(I_to_decision, self.a, self.b, self.d)
        s_decision = self.my_update_s(s_decision, r_decision)

        I_sens_to_assoc = torch.matmul(s_sensory, self.W_sens_to_assoc)
        I_assoc_to_decision = torch.matmul(s_association, self.W_assoc_to_decision)

        I_rec_sens = torch.matmul(s_sensory, self.W_sens_to_sens)/self.input_size
        I_rec_assoc = torch.matmul(s_association, self.W_assoc_to_assoc)/self.association_size
        I_rec_decision = torch.matmul(s_decision, self.W_decision_to_decision)/self.decision_size
