# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:36:17 2018

@author: Y.-H. Song

Program to solve Brueckner equation in symmetric nuclear matter

Following examples given in Haftel and Tabakin paper
  NPA158(1970)1 

To do : instead of effective mass approximation,
        I may use the full s.p. potential U(k). 
"""
import numpy as np
import matplotlib.pyplot as plt 
import scipy.integrate
import scipy.linalg 
import scipy.optimize 

hbarc =197.3269788 # MeV.fm
alpha_em =1./137.035999 # fine structure constrant 
amu = 931.4940954 # MeV/c^2 
mpi = 138.03898   # MeV/c^2 isoscalar pion mass 
mN  = 938.9187       # MeV/c^2 nucleon mass

# gauss lengendre quadrature (-1,1)
# leggau(n) function can be used
def leggau(n,xmin=-1,xmax=+1):
    """ Gaussian quadrature in range (xmin,xmax)
        from quadrature (-1,1)
    """ 
    x,w=np.polynomial.legendre.leggauss(n)
    xp=(xmax-xmin)*0.5*(x+1.0)+xmin
    wp=(xmax-xmin)*0.5*w 
    return xp,wp

def NN_interaction(k,kp,L=0,Lp=0,S=0,J=0,opt=0):
    """
    potential in momentum space, 
    < Lp,pk|V^{S,J}|L, k> 
    (spin conserving interaction at the moment)
    k and kp are momentum in fm^-1 units. 
    <V> is returned in fm power
    """
    if opt==0: # simple separable potential eq.(3.15) in Haftel-Tabakin  
        alph=3.6; a=1.2;
        g0=alph/(k**2+a**2)
        g0p=alph/(kp**2+a**2)
        return -g0*g0p 
    else :
        print('ERROR: interaction not available ')
        return 'ERR'
        
def energy_2N_U(K,q,kF,interp_U,nleg=50):
    """
    angle averaged two nucleon energy under mean potential 
    corresponds to eq. (3.8) in Haftel and Tabakin.
    E(k,K) = (k^2/m+K^2/m+U(|k+K|)+U(|k-K|) )
    eq.(12) of NPA414(1984)_Rikus
    
    E(k,K) = k^2/m + K^2/m + fac1/fac2 
    
    fac1 = \int_{-1}^1 dx Q(K,k,k_F)*(U(|K+k|+U(|K-k|))
    fac2 = \int_{-1}^1 dx Q(K,k,k_F)
    where Q(K,q,k_F) is angle dependence factor.
     
    K : cm momentum of two nucleon in fm^-1 unit
    q : array, rel mom of 2N in fm^-1 unit
    kF : Fermi momentum of nuclear matter in fm^-1 unit 
    interp_U: interpolation function, U(p) in MeV units
    
    return Energy in fm^-1 unit. 
    """
    def Q(K,q,kF,x):
        case0 = np.array(K**2+q**2+2*K*q*x > kF**2)*1.0 # |K+q|> kF
        case1 = np.array(K**2+q**2-2*K*q*x > kF**2)*1.0 # |K-q|> kF
        return case0*case1 
    # fac1  
    xx, wxx = leggau(nleg,-1.0,1.0)
    Kplusq = np.sqrt(K**2+q**2+2*K*q*xx)
    Kmnusq = np.sqrt(K**2+q**2-2*K*q*xx)
    sums1 = 0.0 
    sums2 = 0.0
    for i in range(nleg):
        sums1 = sums1 + wxx[i]*Q(K,q,kF,xx[i])*(interp_U(Kplusq)+interp_U(Kmnusq))
        sums2 = sums2 + wxx[i]*Q(K,q,kF,xx[i]) 
    out = q**2+K**2 + sums1/sums2*mN/hbarc**2    
    return out 



def av_ang_Q(K,kp,kF):
    """
    eq.(3.9) in REF
    
    input K,kp,kF are in fm**-1 unit
    
    vectorized for kp
    """
    if (K > kF) :
        raise ValueError('K or k is too larger than kF')
    if  (K<0)or np.array(kp < 0).any() :    
        raise ValueError('negative k or K')
        
    k1=np.sqrt(kF**2-K**2)
    k2=kF+K
    # vectorized comparison 
    case1 = np.array(kp <= k1)*0.0
    case2 = np.array(kp >= k2)*1.0
    case3 = np.array(k1 < kp)*np.array(kp <= k2)*(K**2+kp**2-kF**2)/(2.0*K*kp)
    return ( case1+case2+case3 )

def example_G00_exact(K,k,k0,kF,effective_ratio=0.8,effective_U0=80.0
                      ,method='leggau',nleg=200):
    """
    eq. (3.16) in REF
    
    Analytic solution of G-matrix for 1S0 seperable potential
    
    input K,k,k0 in fm**(-1) units
    return G00(K,k,k0) in fm unit
    
    effective_ratio=M*/M 
    effective_U0 in MeV unit
    
    """
    fac1=example_pot(k,k0)
    #--integration over k in [0,infty]
    if (method=='leggau'):
      # use change of variable k=C tan x with x in [0,pi/2] 
      xx, wxx = leggau(nleg,0.,np.pi/2)    
      kk  = 0.75*np.tan(xx)
      wkk = wxx*0.75/np.cos(xx)**2
      integrand=wkk*2.0/np.pi*kk**2*example_pot(kk,kk)*av_ang_Q(K,kk,kF)
      denom=(effective_2e(K,kk,'upper',effective_ratio,effective_U0)
          -effective_2e(K,k0,'lower',effective_ratio,effective_U0) )
    elif (method=='laggau'):
      # in case of laggau, weight function is e^-x, thus
      #  require modification of integrand 
      kk, wkk = np.polynomial.laguerre.laggauss(nleg)
      integrand= wkk*np.exp(kk)*( 
            2.0/np.pi*kk**2*example_pot(kk,kk)*av_ang_Q(K,kk,kF) )
      denom=(effective_2e(K,kk,'upper',effective_ratio,effective_U0)
          -effective_2e(K,k0,'lower',effective_ratio,effective_U0) )
            
    integr=np.sum(integrand/denom) 
    #---second part
    fac2=1.0+integr   # since v(k,k)=- g0(k)^2
    return fac1/fac2

def example_G00_invmat(K,k0,kF,effective_ratio=0.8,effective_U0=80.0
                      ,method='leggau',nleg=50):
    """
    eq. (3.14) in REF
    
    Analytic solution of G-matrix for 1S0 seperable potential
    
    input K,k,k0 in fm**(-1) units
    return G00(k,k0;K) in fm unit as an array 
    
    effective_ratio=M*/M 
    effective_U0 in MeV unit
    
    """   
    ##define quadrature
    if (method=='leggau'):
      xx, wxx = leggau(nleg,0.,np.pi/2)    
      kk  = 0.75*np.tan(xx)
      wkk = wxx*0.75/np.cos(xx)**2
      
    if (method=='laggau'):
      # in case of laggau, weight function is e^-x, thus
      #  require modification of integrand 
      kk, wkk = np.polynomial.laguerre.laggauss(nleg)
      wkk = wkk*np.exp(kk) # include negating factor for weight 
      
    # define n+1 array k_i and u_j 
    ki=np.ones(nleg+1)
    ki[0:nleg]=kk
    ki[nleg]=k0
    ui=np.ones(nleg+1)
    ui[0:nleg]=(2.0/np.pi*wkk*kk**2*av_ang_Q(K,kk,kF)
            /(effective_2e(K,kk,'upper',effective_ratio,effective_U0)
             -effective_2e(K,k0,'lower',effective_ratio,effective_U0))  )
    ui[nleg]=0.0
    vv=np.zeros((nleg+1,nleg+1))
      
    for i in np.arange(nleg+1):
        for j in np.arange(nleg+1):
          vv[i,j]=example_pot(ki[i],ki[j])
    mat=np.eye(nleg+1)+np.array(
      [[ui[j]*vv[i,j] for j in np.arange(nleg+1)] for i in np.arange(nleg+1)]  )
    sol=scipy.linalg.solve(mat,vv[:,nleg])
    return sol
      
def test_example_G00(kF=1.45,k0_vals=[1.015,0.725,0.435]
                     ,effective_ratio=0.8,effective_U0=80):
  # comparison between leggau method and laggau method in eq. (3.16)
  K=np.arange(0.01,kF-0.01,0.1)
  plt.figure() 
  plt.ylim(-10,0)
  plt.xlim(0,1.5)
  for k0 in k0_vals:
    G00_0=np.array([example_G00_exact(i,k0, k0,kF
            ,effective_ratio,effective_U0
            ,method='leggau',nleg=50) for i in K])
    plt.plot(K,G00_0,'*')  
    # comparison with  (3.16) and G00 from matrix inversion
    gg=np.array([ example_G00_invmat(i,k0,kF
                      ,effective_ratio=0.8,effective_U0=80.0
                      ,method='leggau',nleg=50)[-1] for i in K])  
    plt.plot(K,gg)  
  return 

def example_U(k, kF, nleg=50):
    """
      compute mean s.p. potential  U(k) eq.(3.20)
      using the G00 of eq.(3.16)
      
      k is not an array 
    """
    if k > kF :
        raise ValueError('ERROR:s.p. potential above Fermi level is attempted')
    #--define k0 quadrature    
    k0_1, wk_1 = leggau(nleg,0.,0.5*(kF-k) ) # 2k_0<=k_F-k_mu
    Kav_1 = np.sqrt(k**2+k0_1**2)
    k0_2, wk_2 = leggau(nleg,0.5*(kF-k),0.5*(kF+k)) # k_F-k_mu<=2k_0<=k_F+k_mu
    Kav_2 = np.sqrt(k**2+k0_2**2-0.25*(2*k0_2+k-kF)*(2*k0_2+k+kF) )
    
    fac1=8.0/np.pi*hbarc**2/mN # np.pi may be missed 
    fac2=3.0 # (2J+1)(2T+1), 3.0 if only 1S0, 6.0 if 3S1 is also included
    # k0 integration          
    integrand1=np.array([ wk_1[i]*k0_1[i]**2
                         *example_G00_exact(Kav_1[i],k0_1[i],k0_1[i],kF
                         ,effective_ratio=0.8,effective_U0=80.0
                         ,method='leggau',nleg=200)  
                         for i in np.arange(nleg)])

    # HT expression seems to contain typos.
    integrand2=np.array([ wk_2[i]*k0_2[i]/(2.*k)
                 *(0.25*(kF**2-k**2)-k0_2[i]*(k0_2[i]-k))
                 *example_G00_exact(Kav_2[i],k0_2[i],k0_2[i],kF
                         ,effective_ratio=0.8,effective_U0=80.0
                         ,method='leggau',nleg=200)
                 for i in np.arange(nleg)] )
    fac3=np.sum(integrand1)+ np.sum(integrand2)
    
    return fac1*fac2*fac3

def example_U_2(k, kF,effective_ratio=0.8,effective_U0=80, nleg=50):
    """
      compute mean s.p. potential  U(k) eq.(3.20)
      using the G00 of eq.(3.16)
      
      k is not an array 
    """
    if k > kF :
        raise ValueError('ERROR:s.p. potential above Fermi level is attempted')
    #--define k0 quadrature    
    k0_1, wk_1 = leggau(nleg,0.,0.5*(kF-k) ) # 2k_0<=k_F-k_mu
    Kav_1 = np.sqrt(k**2+k0_1**2)
    k0_2, wk_2 = leggau(nleg,0.5*(kF-k),0.5*(kF+k)) # k_F-k_mu<=2k_0<=k_F+k_mu
    Kav_2 = np.sqrt(k**2+k0_2**2-0.25*(2*k0_2+k-kF)*(2*k0_2+k+kF) )
    
    fac1=8*hbarc**2/mN/np.pi 
    fac2=3.0 # (2J+1)(2T+1), 3.0 if only 1S0, 6.0 if 3S1 is also included
    # k0 integration          
    integrand1=np.array([ wk_1[i]*k0_1[i]**2
                         *example_G00_exact(Kav_1[i],k0_1[i],k0_1[i],kF
                         ,effective_ratio,effective_U0
                         ,method='leggau',nleg=40)  
                         for i in np.arange(nleg)])

    integrand2=1./(2.0*k)*np.array([ wk_2[i]*k0_2[i]
                 *( (kF**2-k**2)/4.-k0_2[i]*(k0_2[i]-k))
                 *example_G00_exact(Kav_2[i],k0_2[i],k0_2[i],kF
                         ,effective_ratio,effective_U0
                         ,method='leggau',nleg=40)
                 for i in np.arange(nleg)] )
    fac3=np.sum(integrand1)+ np.sum(integrand2)  
    # How to update effective ratio and U0 from U(k)?? 
    return fac1*fac2*fac3 #, fac1*fac2*np.sum(integrand1), fac1*fac2*np.sum(integrand2)

def test_example_U(kF,nleg=50):
    """
    test plot U(k)
    """
    k  = np.arange(0.01,kF,0.05)
    #Uk = np.array([example_U(i,kF,nleg=nleg) for i in k])
    Uk = np.array([example_U_2(i,kF,nleg=nleg) for i in k])
    plt.figure()
    plt.plot(k,Uk)
    plt.ylim(-250,-50)
    return k, Uk
#============================================================
if __name__ == '__main__':
    test_example_G00(1.45) # Fig.1
    test_example_U(1.45,nleg=50) # Fig.2 
    