# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 16:06:40 2024

@author: young
"""
#------------------EC-----------------------------------------------------------
import numpy as np
import scipy 
import torch 
import torch.nn as nn
from torch.autograd import Variable

#==============================================================================
class EigenVectorContinuation():
    """
    H_constructor : a function which returns H(c) matrix 
    train_c       : list of training points 
    level_index   : n-th level to eigen-vector continuation 
    
                    
    For example:
        
    ec_a = EigenVectorContinuation(get_Htot,[0.5,1.5],0)
    ec_a.prepare()
    ec_a.predict(1.0)                 
          
    works to do:  
        (1) use level_index list 
        (2) use list of c for prediction                 
    """
    def __init__(self,H_constructor,train_c,level_index=0):
        self.H_constructor = H_constructor # H(c) matrix constructor 
        self.train_c = np.array(train_c) # training input parameters 
        self.dim_ec = len(train_c) 
        self.train_vectors =None
        self.level_index = 0 # level to use for EC          
        
    def prepare(self,):
        eigvals = np.zeros_like(self.train_c)
        eigvecs = [ ] 
        for i,c in enumerate(self.train_c):
            Htot = self.H_constructor(c)
            ee,vv = scipy.linalg.eigh(Htot) #sorted eigenvalue and eigenvector   
            eigvals[i] = ee[self.level_index] # n-th eigenvalue
            eigvecs.append( vv[:,self.level_index] ) # n-th eigenvector 
        eigvecs = np.array(eigvecs)     
        self.train_vectors = eigvecs 
        return  
    
    def get_reduced_Hamiltonian(self,c):
        #--EC Hamiltonian and Norm matrix
        Htot = self.H_constructor(c) 
        H_EC = np.zeros( (self.dim_ec,self.dim_ec) )
        N_EC = np.zeros( (self.dim_ec,self.dim_ec) )
        for i in range(self.dim_ec):
            vL = self.train_vectors[i]
            for j in range(self.dim_ec):
                vR = self.train_vectors[j]
                H_EC[i,j] = np.matmul(vL, np.matmul(Htot,vR))
                N_EC[i,j] = np.matmul(vL, vR)
        return H_EC, N_EC 
    
    def predict(self,c):
        H_EC, N_EC = self.get_reduced_Hamiltonian(c)
        ee, vv = scipy.linalg.eigh(H_EC,N_EC)
        return ee[self.level_index], vv[:,self.level_index]

#============================================================================== 
def random_Hermitian_variable(ndim):   
    #----construct Hermitian matrix from input ndim**2 variable 
    var = Variable(torch.randn(ndim),requires_grad=True)
    var2 = Variable(torch.randn(int(ndim*(ndim-1)/2),dtype=torch.complex64),requires_grad=True)
    mm = torch.zeros((ndim,ndim))
    indices = np.diag_indices(ndim)
    mm[indices] = var 
    mm = mm.type(torch.complex64)
    indices = np.triu_indices(ndim,1) # indices of upper traiangle without diagonal 
    mm[indices] = var2
    mm = mm + torch.transpose(mm.conj(),0,1) 
    return mm  , (var,var2)
   
class ParametricMatrixModel():
    """
      ndim : dimension of PMM 
      train_data : set of c_list, eigenvalues 
               
      1. create Variables for PMM 
      2. create Hermitian Matrix 
    """
    def __init__(self,ndim,ansatz_constructor,n_mats,optimizer=None,loss_type=None):
        self.ansatz = ansatz_constructor  
        self.optim = None 
        self.get_loss = None 
        self.n_mats = n_mats # number of Hermitian matrices for ansatz 
        return 
    def prepare(self,):
        return 
    def train(self,train_data,epochs=5000):
        return 
    def predict(self,c):
        return 