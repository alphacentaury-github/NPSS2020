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

## Test calculation with the answer 

#----set potential
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
H_V[ abs(dx)< R0,abs(dx)< R0 ] = V0
Htot = Hfree+H_V

#----get first 20 eigenvalues
#ee = scipy.linalg.eigvals(Htot)
#print(np.sort(ee)[:20]*cutoff)
#--or use torch
print('E_lattice={} MeV'.format(torch.linalg.eigvalsh(torch.Tensor(Htot))[:4]*cutoff))

## Test PMM with square-well potential problem 

#---prepare data
clist = [0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9]
#clist = [0.1,0.3,0.5,0.8]
Hfree = get_HKin(iKin)
R0 = 3.0/a_fm # fm -> l.u 
V0 = -150.0/cutoff # MeV -> l.u 
data=[]
for c in clist:
    H_V = np.zeros((L,L))
    H_V[ abs(dx)< R0,abs(dx)< R0 ] = c*V0
    Htot = Hfree+H_V
    data.append(torch.linalg.eigvalsh(torch.Tensor(Htot)).numpy()[:2])
data = np.array(data)
print(data)    
data_tensor = torch.tensor(data)

## Eigenvector Continuation test
#--prepare EC training data 
ec_evals0=[]
ec_evals1=[]
ec_vects0=[]
ec_vects1=[]
ec_mats =[]
for c in [0.5,1.5]:
    H_V = np.zeros((L,L))
    H_V[ abs(dx)< R0,abs(dx)< R0 ] = c*V0
    Htot = Hfree+H_V
    #ee,vv = np.linalg.eigh(Htot)
    ee,vv = scipy.linalg.eigh(Htot)
    ec_mats.append(Htot)
    ec_evals0.append(ee[0]) # gs(?) it is not clear whether the ordering is correct...
    ec_vects0.append(vv[:,0])
    ec_evals1.append(ee[1]) # gs(?) it is not clear whether the ordering is correct...
    ec_vects1.append(vv[:,1])
    # print('test eigen:{} <H>:{}'.format(ee[0], vv[:,0]@ Htot@vv[:,0]  ) )
ec_mats = np.array(ec_mats)
ec_evals0 = np.array(ec_evals0)
ec_vects0 = np.array(ec_vects0)
ec_evals1 = np.array(ec_evals1)
ec_vects1 = np.array(ec_vects1)
#----consistency test of eigenvector, <H>
for i in range(2):
    print('{}: eigv={} <H>={}'.format(i,ec_evals0[i], ec_vects0[i,:]@ ec_mats[i] @ ec_vects0[i,:]))

def EC_test0(c):
    #----prediction by EC
    H_V = np.zeros((L,L))
    H_V[ abs(dx)< R0,abs(dx)< R0 ] = c*V0
    Htot = torch.Tensor(Hfree+H_V) # H(c) 
    #---construct reduced matrix and norm matrix 
    Mc = np.zeros((2,2))
    Nc = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            Mc[i,j] = np.matmul(ec_vects0[i,:],np.matmul(Htot,ec_vects0[j,:]))
            Nc[i,j] = np.matmul(ec_vects0[i,:],ec_vects0[j,:])
    ee,vv = scipy.linalg.eig(Mc,Nc) 
    return ee, vv

def EC_test1(c):
    #----prediction by EC
    H_V = np.zeros((L,L))
    H_V[ abs(dx)< R0,abs(dx)< R0 ] = c*V0
    Htot = torch.Tensor(Hfree+H_V) # H(c) 
    #---construct reduced matrix and norm matrix 
    Mc = np.zeros((2,2))
    Nc = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            Mc[i,j] = np.matmul(ec_vects1[i,:],np.matmul(Htot,ec_vects1[j,:]))
            Nc[i,j] = np.matmul(ec_vects1[i,:],ec_vects1[j,:])
    ee,vv = scipy.linalg.eig(Mc,Nc) 
    return ee, vv

ee,vv=EC_test0(1.5)
print('EC e={}'.format(ee[1]))

### 2d PMM

s0 = torch.tensor( [ [1.0,0.0],[0.0,1.0]])
sx = torch.tensor( [ [0.0,1.0],[1.0,0.0]])
sy = torch.tensor( [ [0.0+0j,-1j],[1j,0.0+0j]],dtype=torch.complex64)
sz = torch.tensor( [ [1.0,0.0],[0.0,-1.0]])

var1 = Variable(torch.tensor([0.5,0.5,0.5,0.5]),requires_grad=True)
var2 = Variable(torch.tensor([0.5,0.5,0.5,0.5]),requires_grad=True)
optim = torch.optim.Adam([var1,var2],lr=5.e-3)
#optim = torch.optim.SGD([var1,var2])
get_loss = nn.MSELoss(reduction='sum')

def mm_var(var):
    return s0*var[0]+sx*var[1]+sy*var[2]+sz*var[3]

for ii in range(5000):
  dd = torch.zeros_like(data_tensor)
  for i,c in enumerate(clist):
      M1 = mm_var(var1)
      M2 = mm_var(var2)
      Mtot = M1 + c *M2
      dd[i] = torch.linalg.eigvalsh(Mtot)[:2]
    
  loss = get_loss(dd,data_tensor)
  optim.zero_grad()
  loss.backward()
  optim.step()
  if ii %1000 ==0:  print('epoch={} loss={}'.format(ii,loss))
print(loss)  

out_t=[]
out_p=[]
out_ec0=[]
out_ec1=[]
c_p = np.arange(0.05,2.5,0.1)
for c in c_p:
  H_V = np.zeros((L,L))
  H_V[ abs(dx)< R0,abs(dx)< R0 ] = c*V0
  Htot = Hfree + H_V
  out_t.append(torch.linalg.eigvalsh(torch.Tensor(Htot)).numpy()[:2])
    
  with torch.inference_mode():
    M1 = mm_var(var1)
    M2 = mm_var(var2)
    Mtot = M1 + c *M2
    out_p.append(torch.linalg.eigvalsh(Mtot).numpy())
  #--EC results 
  ee,vv = EC_test0(c)
  out_ec0.append(np.real( np.min(ee) ))
  ee,vv = EC_test1(c)
  out_ec1.append(np.real( np.min(ee) ))
out_t = np.array(out_t)
out_p = np.array(out_p)
out_ec0 = np.array(out_ec0)
out_ec1 = np.array(out_ec1)

#-----plot------------------------------
plt.plot(clist,data[:,0],'g*')
plt.plot(clist,data[:,1],'g*')
plt.plot(c_p,out_t[:,0],'r--',label='true e0' )
plt.plot(c_p,out_t[:,1],'b--',label='true e1' )
plt.plot(c_p,out_p[:,0],'r',label='PMM e0' )
plt.plot(c_p,out_p[:,1],'b',label='PMM e1' )
plt.plot(c_p,out_ec0,'y-.',label='EC e0')
plt.plot(c_p,out_ec1,'y-.',label='EC e1')
plt.xlabel('c');plt.ylabel('E [l.u.]')
plt.legend()