# -*- coding: utf-8 -*-
"""
 Try to solve deuteron as C.C. problem 
 in relative momentum coordinates. 
 
 In relative coordinate C.M motion is no concern and 
 problem is a one-body problem. 
 particle in J=1, M=0 state. 

 Hamiltonian is constructed as  
 <p_i,a| H |p_j,b>, i is momentum index, a= S or D for spin index
                                           
 # However strangely there are some factors for potential  
 <p_i,a| H |p_j,b> = p_i^2 *delta_{ij}
                    + sqrt(w_i w_j) p_i p_j V(p_i,a,p_j,b)
 # engenvalue of this matrix gives deuteron energy                   
                    
 # because it is only a one-particle system,
   reference state can be |p=0,a=S> and 
   other states can be thought as 1p1h excited state
   only single amplitude t^a_0 is available. 
   
   Thus, e^T|p=0,S> = |p=0,S> +\sum_a t^a_0 |p_a,s_a>    
   This seems to be a simple linear combination of momentum basis.  
   Thus, t^a_0 is expected to be a eigen-vector of Hamiltonian. 
   
"""
from srg_nn import * 
#-------------------
# construct Hamiltonian
#-------------------
mom_tmp = read_mesh("n3lo500_3s1.meq")
momenta = np.concatenate([mom_tmp,mom_tmp])
weights = uniform_weights(momenta)
dim     = len(momenta)

# set up p^2 (kinetic energy in units where h^2/2\mu = 1)
# p itself is in fm^-1 units. 
# to return MeV units, one have to multiply hbarm 
T = diag(momenta*momenta) # fm^-2 units. 

# set up interaction matrix in coupled channels:
#
#    /   V_{3S1}            V_{3S1-3D1}  \
#    \   V_{3S1-3D1}^\dag   V_{3D1}      /

# read individual partial waves
partial_waves=[]
for filename in ["n3lo500_3s1.meq", "n3lo500_3d1.meq", "n3lo500_3sd1.meq"]:
  partial_waves.append(read_interaction(filename))
  # print(partial_waves[-1].shape

# assemble coupled channel matrix
V = np.vstack((np.hstack((partial_waves[0], partial_waves[2])), 
               np.hstack((np.transpose(partial_waves[2]), partial_waves[1]))
              ))

# switch to scattering units
# To change units where h^2/2\mu = 1. 
V = V/hbarm

# set up conversion matrix for V: this is used to absorb momentum^2 and
# weight factors into V, so that we can use the commutator routines for
# eta and the derivative as is
conversion_matrix = np.zeros_like(T)
for i in range(dim):
  for j in range(dim):
    # Regularize the conversion matrix at zero momentum - set elements
    # to machine precision so we can invert the matrix for plots etc.
    # Note that momentum values are positive, by construction.
    qiqj = max(np.finfo(float).eps, momenta[i]*momenta[j])
    conversion_matrix[i,j] = qiqj*sqrt(weights[i]*weights[j])

V *= conversion_matrix  # eq. (10.68) of book 

# Hamiltonian matrix 
H = (T+V)*hbarm
#----------------------
# C.C. 
#---------------------
# single amplitude and normal ordered F_N 
ts = np.zeros((dim,dim))
Eref = H[0,0]
fs = H #- np.diag(np.ones(dim)*Eref)
# initialize ts 
for a in range(1,dim):
    ts[a,0] = 1.0
    
def ccsdenergy(fs,ts):
  ECCSD = 0.0
  i = 0 
  for a in range(1,dim):
      ECCSD += fs[i,a]*ts[a,i]
  return ECCSD

def makeT1(fs,ts):
    tsnew = np.zeros((dim,dim))
    for a in range(1,dim):
        hh = fs[a,0]
        for c in range(1,dim):
            hh += fs[a,c]*ts[c,0]
            hh -= fs[0,0]*ts[a,0]
            hh -= fs[0,c]*ts[c,0]*ts[a,0]
        #---update ts 
        tsnew[a,0] = ts[a,0]+hh/(fs[a,a]-fs[0,0])
    return tsnew     
#================
# MAIN LOOP
# CCSD iteration
#================
print('---init----')
print("FCI = %14.8f"%(eigvalsh(H)[0]))
print("Eref = %14.8f Ecorr = %14.8f"%(Eref,ccsdenergy(fs,ts)))


ECCSD = 1.0
DECC = 1.0
while DECC > 1.0e-12: # arbitrary convergence criteria
  OLDCC = ECCSD
  tsnew = makeT1(fs,ts)
  ts = tsnew
  ECCSD = ccsdenergy(fs,ts)
  DECC = abs(ECCSD - OLDCC)
  print('%14.8f %14.8f'%(ECCSD,DECC)  )
  
print(' ??? Why it does not work? ??')  
  