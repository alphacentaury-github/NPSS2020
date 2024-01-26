# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 15:12:00 2024

@author: young
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt

#-----------setup lattice----------------------------------
hbarc = 197.3
L = 40
cutoff = 200.0 # MeV 
a_fm = hbarc/cutoff ;print('a = {} fm L={} fm'.format(a_fm, L*a_fm))
mass = 938.92/cutoff #l.u.
iKin = 0             # kinectic term option 

r = np.arange(L)
nx = r % L
dx = (nx +L/2.0) % L - L/2.0
dr = np.sqrt(dx**2)
qx = dx*(2*np.pi)/L  # momentum qx in range(-2pi/L,2pi/L) 
q2 = qx**2

#--------setup Hfree
def get_HKin(iKin,test=False):
    """
    construct K(i,j) matrix with approximation
    iKin=0 case, one have to use FFT. 
    
    Here we are not using sparse matrix.
    """
    if iKin==0:
        diffs = np.fft.ifft(qx**2)
        temp = np.eye(L)*0j
        for i in range(L):
            for j in range(L):
                temp[i,j] = diffs[(i-j)%L]
        Hfree = np.real(temp)
    else:
        if iKin==1: w0, w1, w2, w3 = [1.0,1.0,0.0,0.0]
        if iKin==2: w0, w1, w2, w3 = [49./36.,3./2.,3./20.,1./90.]
         
        Hfree = np.eye(L)*w0*2 # w0
        Hfree[r, (r+1)%L] = w1*(-1)
        Hfree[r, (r-1)%L] = w1*(-1)
        Hfree[r, (r+2)%L] = w2
        Hfree[r, (r-2)%L] = w2
        Hfree[r, (r+3)%L] = w3*(-1)
        Hfree[r, (r-3)%L] = w3*(-1)
        # in case of 1D, no overlap between hopping. 
        # However, in other dimension, one have to use summation. 
    
    if test:
        ee,vv = scipy.linalg.eig(Hfree)
        print('p^2  Hkin:{} '.format(np.sort(ee)))
        if iKin==0: print('p^2  exa :{} '.format(np.sort(qx**2) )) 
        if iKin>0:
            print('p^2  exa :{} '.format(np.sort(2*(w0-w1*np.cos(qx)+w2*np.cos(2*qx)-w3*np.cos(3*qx))) ))
        
    return Hfree/(2.0*mass)

def get_sqwell_pot(c=1,R0=3.0/a_fm,V0 = -20./cutoff):    
    # square-well-potential case 
    print('square well V0={} MeV, R0={} fm'.format(c*V0*cutoff,R0*a_fm))
    H_V = np.zeros((L,L))
    H_V[ abs(dx)< R0,abs(dx)< R0 ] = c*V0
    return H_V 

def get_Htot(potential_constructor,*args,**kwargs):
    """
    change this routine for general problem
    """
    H_V = potential_constructor(*args,**kwargs)
    return Hfree+H_V 

def test_sqwell_sol():
    #--test calculation with known answer
    #----set potential
    R0=3.0/a_fm
    V0 = -20./mass/R0**2
    #---Wiki sol-------
    vv = np.array([1.28, 2.54, 3.73])
    print('E_continuum={} MeV'.format( (2*vv**2/mass/(2*R0)**2 +V0)*cutoff))

    Htot = get_Htot(get_sqwell_pot,1,R0,V0)
    ee = scipy.linalg.eigvalsh(Htot) 
    print('E_lattice={} MeV'.format(ee[:4]*cutoff))
    return 
#---------------------test PMM--------------------------------------------

#===============================================================================
if __name__ == "__main__":
    """
    1d problem test with for EC and PMM
    """    
    #----initialization--------------
    Hfree = get_HKin(iKin)

    #-----test
    print('test problem')
    test_sqwell_sol()
    
    #----- PMM test Problem-------------------------- 
    #-----define Hamiltonian constructor H(c) 
    R0 = 3.0/a_fm # fm -> l.u 
    V0 = -150.0/cutoff # MeV -> l.u  
    getHc = lambda c: get_Htot(get_sqwell_pot,c,R0,V0)
    
    #---list of train data 
    clist = [0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9]
    data=[]
    for c in clist:
        Htot = getHc(c)
        data.append(scipy.linalg.eigvalsh(Htot)[:2])
    data = np.array(data)
    print(data)    
    
    #-----construct torch model
    import torch 
    import torch.nn as nn 
    import torch.optim as optim
    
    class MyModule(nn.Module):
        """
        model of 2d PMM 
        """
        def __init__(self,):
            super(MyModule,self).__init__()
            self.var1 = torch.Variable([1.0,1.0,1.0,1.0],require_grad=True)
            self.var2 = torch.Variable([1.0,1.0,1.0,1.0],require_grad=True)
        def forward(self,input):
            def mm_var(var):
                s0 = torch.tensor( [ [1.0,0.0],[0.0,1.0]])
                sx = torch.tensor( [ [0.0,1.0],[1.0,0.0]])
                sy = torch.tensor( [ [0.0+0j,-1j],[1j,0.0+0j]],dtype=torch.complex64)
                sz = torch.tensor( [ [1.0,0.0],[0.0,-1.0]])
                return s0*var[0]+sx*var[1]+sy*var[2]+sz*var[3]
            M1 = mm_var(self.var1)
            M2 = mm_var(self.var2)
            Mtot = M1 + input *M2
            dd = torch.linalg.eigvalsh(Mtot)
            return dd 
    
    
    