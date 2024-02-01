# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 15:12:00 2024

@author: young
"""
import numpy as np
import torch
import torch.nn as nn
from torch import FloatTensor
from torch.autograd import Variable
import scipy
import matplotlib.pyplot as plt

from EigenVectorContinuation import * 

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
        ee,vv = scipy.linalg.eigh(Hfree)
        print('p^2  Hkin:{} '.format(np.sort(ee)))
        if iKin==0: print('p^2  exa :{} '.format(np.sort(qx**2) )) 
        if iKin>0:
            print('p^2  exa :{} '.format(np.sort(2*(w0-w1*np.cos(qx)+w2*np.cos(2*qx)-w3*np.cos(3*qx))) ))
        
    return Hfree/(2.0*mass)

## Test calculation with the answer 

def WS_pot(x,V0,R0,a=1.0):
    # evrything  in l.u.
    # x must be a distance, not coordinate  
    return V0/(1.+np.exp((abs(x)-R0)/a))
def sqwell_pot(x,V0,R0):
    # evrything  in l.u.
    # x must be a distance, not coordinate  
    return V0*(abs(x) <= R0)

#----set potential
def test_case_1d():
    R0 = 3.0/a_fm # 2 fm -> l.u 
    #V0 = -10.0/cutoff # 10 MeV -> l.u 
    V0 = -20./mass/R0**2
    print('V0={} MeV, R0={} fm'.format(V0*cutoff,R0*a_fm))

    #---Wiki sol-------
    vv = np.array([1.28, 2.54, 3.73])
    print('E_continuum={} MeV'.format( (2*vv**2/mass/(2*R0)**2 +V0)*cutoff))

    #---solve S.E. 
    Hfree = get_HKin(iKin)
    H_V = np.zeros((L,L))
    #H_V[ abs(dx)< R0,abs(dx)< R0 ] = V0    
    VV = sqwell_pot(dx,V0,R0)
    H_V[r,r] = VV 
    
    Htot = Hfree+H_V

    #----get first 20 eigenvalues
    #ee = scipy.linalg.eigvals(Htot)
    #print(np.sort(ee)[:20]*cutoff)
    #--or use torch
    print('E_lattice={} MeV'.format(scipy.linalg.eigvalsh(Htot)[:4]*cutoff))
    #print('E_lattice={} MeV'.format(torch.linalg.eigvalsh(torch.Tensor(Htot))[:4]*cutoff))
    return 

#---prepare data
#clist = [0.5,0.6,0.7,0.8,0.9,1.5]
#clist = [0.1,0.3,0.5,0.8]
clist = [0.5,0.52,0.55,0.6,0.62,0.65,0.7,0.72,0.75,0.8,0.82,0.85,0.9,1.5] # a training points for PMM 
c_p = np.arange(0.1,2.5,0.1) # full test range 

#---construct Hamiltonian 
R0 = 3.0/a_fm # fm -> l.u 
V0 = -100.0/cutoff # MeV -> l.u 
a0 = 1.0 # in l.u. 

Hfree = get_HKin(iKin) 

def get_Htot(c):
    H_V = np.zeros((L,L))
    VV = WS_pot(dx,V0,c*R0,a0)
    H_V[r,r] = VV 
    Htot = Hfree+H_V
    return Htot 

#------training data-------------------------------------------------- 
data=[]
for c in clist:
    Htot = get_Htot(c) 
    data.append(torch.linalg.eigvalsh(torch.Tensor(Htot)).numpy()[:2])
data = np.array(data)
 
data_tensor = torch.tensor(data)

#-----eigenvector continuation----------------------------------------
## Eigenvector Continuation test
#--prepare EC training data 

clist_ec = [0.5,1.5]

ec_test0 = EigenVectorContinuation(get_Htot,clist_ec,0)
ec_test0.prepare()  
ec_test1 = EigenVectorContinuation(get_Htot,clist_ec,1)
ec_test1.prepare()
ec_evals0=[]
ec_evals1=[]
ec_vects0=[]
ec_vects1=[]
ec_Htot =[]

for c in clist_ec:
    Htot = get_Htot(c) 
    ee,vv = ec_test0.predict(c)
    ec_evals0.append(ee)
    ec_vects0.append(vv)
    ee,vv = ec_test1.predict(c)
    ec_evals1.append(ee)
    ec_vects1.append(vv)

ec_evals0 = np.array(ec_evals0)
ec_vects0 = np.array(ec_vects0) # gs eigen vectors at c1,c2 
ec_evals1 = np.array(ec_evals1)
ec_vects1 = np.array(ec_vects1) # ex eigen vectors at c1,c2 
#-------------------------------------------------------------------------------------------
# test of ground and excited states together 
def EC_both(c):
    #----prediction by EC
    H_V = np.zeros((L,L))
    #H_V[ abs(dx)< R0,abs(dx)< R0 ] = c*V0
    VV = WS_pot(dx,V0,c*R0,a0)
    H_V[r,r] = VV 
    Htot = torch.Tensor(Hfree+H_V) # H(c) 
    #---construct reduced matrix and norm matrix 
    Mc = np.zeros((4,4))
    Nc = np.zeros((4,4))
    for i in range(4):
        if i in [0,1]:
            vL = ec_vects0[i,:]
        else:
            vL = ec_vects1[i-2,:]
        for j in range(4):
            if j in [0,1]:
                vR = ec_vects0[j,:]
            else:
                vR = ec_vects1[j-2,:]
            Mc[i,j] = np.matmul(vL,np.matmul(Htot,vR))
            Nc[i,j] = np.matmul(vL,vR)
    ee,vv = scipy.linalg.eigh(Mc,Nc) 
    return ee, vv, Mc, Nc 
    
### 2d PMM================================================================
# This 2d PMM is trained for gs and 1st excited states
# And Hermitian Matrices are prepared by using Pauli matrices...
# However, that might be a problem. 
# One may need to increase the dimension of PMM 
# or use separate PMM for g.s. and excited state. 
s0 = torch.tensor( [ [1.0,0.0],[0.0,1.0]])
sx = torch.tensor( [ [0.0,1.0],[1.0,0.0]])
sy = torch.tensor( [ [0.0+0j,-1j],[1j,0.0+0j]],dtype=torch.complex64)
sz = torch.tensor( [ [1.0,0.0],[0.0,-1.0]])

var1 = Variable(torch.tensor([0.5,0.5,0.5,0.5]),requires_grad=True)
var2 = Variable(torch.tensor([0.5,0.5,0.5,0.5]),requires_grad=True)
var3 = Variable(torch.tensor([0.5,0.5,0.5,0.5]),requires_grad=True)
var4 = Variable(torch.tensor([0.5,0.5,0.5,0.5]),requires_grad=True)
var5 = Variable(torch.tensor([0.5,0.5,0.5,0.5]),requires_grad=True)

optim = torch.optim.Adam([var1,var2],lr=2.e-3)
get_loss = nn.MSELoss(reduction='sum')

def mm_var(var):
    return s0*var[0]+sx*var[1]+sy*var[2]+sz*var[3]
    
def get_Mtot(c): # ansatz of PMM 
    M1 = mm_var(var1)
    M2 = mm_var(var2)
    Mtot = M1 + M2/c 
    return Mtot 

#=====repeat as many to reduce the loss 
for ii in range(6000):
  dd = torch.zeros_like(data_tensor[:,0]) # use data_tensor to use both energies 
  for i,c in enumerate(clist):
      Mtot = get_Mtot(c)
      dd[i] = torch.linalg.eigvalsh(Mtot)[:1] # [:2] to use both energies  
    
  loss = get_loss(dd,data_tensor[:,0]) # use data_tensor to use both energies 
  optim.zero_grad()
  loss.backward()
  optim.step()
  if ii %1000 ==0:  print('epoch={} loss={}'.format(ii,loss))
print(loss)  

out_t=[]
out_p=[]
out_ec0=[]
out_ec1=[]

ec_mats=[] 
for c in c_p:
  H_V = np.zeros((L,L))
  #H_V[ abs(dx)< R0,abs(dx)< R0 ] = c*V0
  VV = WS_pot(dx,V0,c*R0,a0)
  H_V[r,r] = VV 
  Htot = Hfree + H_V
  out_t.append(torch.linalg.eigvalsh(torch.Tensor(Htot)).numpy()[:2])
    
  with torch.inference_mode():
    Mtot = get_Mtot(c)
    out_p.append(torch.linalg.eigvalsh(Mtot).numpy())
  #--EC results 
  temp1, temp2 = ec_test0.get_reduced_Hamiltonian(c)
  ee,vv = ec_test0.predict(c)
  out_ec0.append(np.real(ee))
  isqnn = scipy.linalg.fractional_matrix_power(temp2,-0.5) #sqrt(N^-1)
  ec_mats.append( isqnn @ temp1 @ isqnn ) # Hermitian EC Hamiltonian 
  
  ee,vv = ec_test1.predict(c)
  out_ec1.append(ee)
  
out_t = np.array(out_t)
out_p = np.array(out_p)
out_ec0 = np.array(out_ec0)
out_ec1 = np.array(out_ec1)
#---------------------------------------
#tt =[]
#for i,c in enumerate(c_p):
#    tt.append(ec_mats[i][0,0])
#tt=np.array(tt)

##f = lambda x, a, b,c,d,e : a*x**2+b*x+c+d/x+e/x**2
#f = lambda x, a, b,c,d,e : a*x**2+b*x+c+d/x
##f = lambda x, b,c,d : b*x+c+d/x
##f = lambda x, a, b,c : a*np.exp(-b*x)+c
#popt,pcov =scipy.optimize.curve_fit(f,c_p,tt)
#plt.plot(c_p,tt,'.',label='EC')
#plt.plot(c_p,f(c_p,*popt),label='fit')
#plt.xlabel('c');plt.ylabel(r'$M_{00}$')
#plt.legend()     
#---------------------------------------
#out_ec_test=[]
#for c in c_p:
#    ee,vv,Mc,Nc = EC_test(c)
#    #ee,vv,Mc,Nc = EC_both(c)
#    out_ec_test.append(np.sort(np.real(ee)))
#out_ec_test= np.array(out_ec_test)    

#-----plot------------------------------
plt.plot(clist,data[:,0],'g*')
plt.plot(clist,data[:,1],'g*')
plt.plot(c_p,out_t[:,0],'r--',label='true e0' )
plt.plot(c_p,out_t[:,1],'b--',label='true e1' )
plt.plot(c_p,out_p[:,0],'r',label='PMM e0' )
#plt.plot(c_p,out_p[:,1],'b',label='PMM e1' )
plt.plot(c_p,out_ec0,'y-.',label='EC e0')
plt.plot(c_p,out_ec1,'y-.',label='EC e1')
#plt.plot(c_p,out_ec_test[:,0],'g-.',label='EC_test e0')
#plt.plot(c_p,out_ec_test[:,1],'g-.',label='EC_test e1')

plt.xlabel('c');plt.ylabel('E [l.u.]')
plt.ylim([-0.5,0.1])
plt.legend()