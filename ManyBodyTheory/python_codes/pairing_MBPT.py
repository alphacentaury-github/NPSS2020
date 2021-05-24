#!/usr/bin/python
#===================
# MBPT2 
#=====================
from sympy import *
from pylab import *
import matplotlib.pyplot as plt

#---Hamiltonian-------------
# E(p,sigma) = delta*(p-1)
def h0(p,q):
    if p == q:
        p1, s1 = states[p]
        return (p1 - 1)
    else:
        return 0
    
# Fock term f
def f(p,q):
    if p == q:
        return 0
    s = h0(p,q)
    for i in below_fermi:
        s += assym(p,i,q,i)
        return s

#---- v(p,q,r,s) = -0.5*g*( p+,p-,q-,q+) 
def assym(p,q,r,s):
    p1, s1 = states[p]
    p2, s2 = states[q]
    p3, s3 = states[r]
    p4, s4 = states[s]

    if p1 != p2 or p3 != p4:
        return 0
    if s1 == s2 or s3 == s4:
        return 0
    if s1 == s3 and s2 == s4:
        return -g/2.
    if s1 == s4 and s2 == s3:
        return g/2.

#---sum of s.p. energies 
# eps = e_i+e_j-e_a-e_b 
def eps(holes, particles):
    E = 0
    for h in holes:
        p, s = states[h]
        E += (p-1)
    for p in particles:
        p, s = states[p]
        E -= (p-1)
    return E


def MBPT2(assym=assym,eps=eps,above_fermi=above_fermi,below_fermi=below_fermi):
  """
   compute MBPT2 diagrams 
   for given particles and holes list and 
      one-body energy difference eps (which is from one-body operator f )
      two-body operator assym 
  """
  # Diagram 1 in MBPT2 
  s1 = 0
  for a in above_fermi:
    for b in above_fermi:
        for i in below_fermi:
            for j in below_fermi:
                s1 += 0.25*assym(a,b,i,j)*assym(i,j,a,b)/eps((i,j),(a,b))


  # Diagram 2 is zero for pairing model 
  s2 = 0
  # Diagram 3 in MBPT3 
  s3 = 0
  for a in above_fermi:
     for b in above_fermi:
       for c in above_fermi:
           for i in below_fermi:
               for j in below_fermi:
                   for k in below_fermi:
                       s3 += assym(i,j,a,b)*assym(a,c,j,k)*assym(b,k,c,i)/eps((i,j),(a,b))/eps((k,j),(a,c))

  # Diagram 4 in MBPT3 
  s4 = 0
  for a in above_fermi:
    for b in above_fermi:
        for c in above_fermi:
            for d in above_fermi:
                for i in below_fermi:
                    for j in below_fermi:
                        s4 += 0.125*assym(i,j,a,b)*assym(a,b,c,d)*assym(c,d,i,j)/eps((i,j),(a,b))/eps((i,j),(c,d))

  # Diagram 5 MBPT3
  s5 = 0
  for a in above_fermi:
    for b in above_fermi:
        for i in below_fermi:
            for j in below_fermi:
                for k in below_fermi:
                    for l in below_fermi:
                        s5 += 0.125*assym(i,j,a,b)*assym(k,l,i,j)*assym(a,b,k,l)/eps((i,j),(a,b))/eps((k,l),(a,b))

  # Diagram 6,7 are zero for pairing model                         
  s6 = 0 
  s7 = 0
  # Diagram 8 are zero for pairing model
  s8 = 0
  for a in above_fermi:
    for b in above_fermi:
        for i in below_fermi:
            for j in below_fermi:
                for k in below_fermi:
                    s8 -= 0.5*assym(i,j,a,b)*assym(a,b,i,k)*f(k,j)/eps((i,j),(a,b))/eps((i,k),(a,b))

  # Diagram 9 are zero for pairing model  
  s9 = 0
  for a in above_fermi:
    for b in above_fermi:
        for c in above_fermi:
            for i in below_fermi:
                for j in below_fermi:
                    s9 += 0.5*assym(i,j,a,b)*assym(a,c,i,j)*f(b,c)/eps((i,j),(a,b))/eps((i,j),(a,c))
  return s1,s2,s3,s4,s5,s6,s7,s8,s9                   

if __name__ == "__main__":
  #---s.p. levels/states---------
  below_fermi = (0,1,2,3)
  above_fermi = (4,5,6,7)
  states = [(1,1),(1,-1),(2,1),(2,-1),(3,1),(3,-1),(4,1),(4,-1)] # (p,sigma)
  g = Symbol('g')
  
  ga = linspace(-1,1,20)
  e1 = []
  corr2 = []
  corr3 = []
  
  (s1,s2,s3,s4,s5,s6,s7,s8,s9) = MBPT2()
        
  for g_val in ga:
    H1 = matrix([[2-g_val , -g_val/2.,  -g_val/2., -g_val/2., -g_val/2.,     0],
                 [-g_val/2.,   4-g_val,  -g_val/2., -g_val/2.,    0., -g_val/2.],
                 [-g_val/2., -g_val/2.,    6-g_val,     0, -g_val/2., -g_val/2.],
                 [-g_val/2., -g_val/2.,      0,   6-g_val, -g_val/2., -g_val/2.],
                 [-g_val/2.,     0,  -g_val/2., -g_val/2.,   8-g_val, -g_val/2.],
                 [0    , -g_val/2.,  -g_val/2., -g_val/2., -g_val/2.,  10-g_val]])

    u1, v1 = linalg.eig(H1)
    e1.append(min(u1)) # exact sol 
    corr2.append((s1).subs(g,g_val)) # MBPT2 # subs because 'g' is symbol 
    corr3.append((s1+s3+s4+s5).subs(g,g_val)) # MBPT3 # subs because 'g' is symbol

  exact = e1 - (2-ga) # correlation energy = E - E_0 

  plt.axis([-1,1,-0.5,0.05])
  plt.xlabel(r'Interaction strength, $g$', fontsize=16)
  plt.ylabel(r'Correlation energy', fontsize=16)
  fci = plt.plot(ga, exact,'b-*',linewidth = 2.0, label = 'Exact')
  mbpt2 = plt.plot(ga, corr2,'r:.', linewidth = 2.0, label = 'MBPT2')
  mbpt3 = plt.plot(ga, corr3, 'm:v',linewidth = 2.0, label = 'MBPT3')
  plt.legend()
  #plt.savefig('perturbationtheory.pdf', format='pdf')
  plt.show()