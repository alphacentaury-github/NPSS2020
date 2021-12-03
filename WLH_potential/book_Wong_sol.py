# -*- coding: utf-8 -*-
"""
Solutions to the Introductory nuclear physics by S. Wong
"""


import scipy 
import scipy.special  
import numpy as np 

def example_5_1():
    """
    5-1. 20Ne(2+,1.634 MeV) lifetime = 0.655 ps. 
     What is the B(E2) and decay probability W ? 
     
     One can obtain W , assume 100% branching ratio,
     W = ln 2/ T_{1/2} 
       
     Then, B(E2) can be obtained by using eq. (5.27) 
     and below.
     
     According to NRV.
     
     lifetime is 0.73 ps. 
     B(E2) = 0.0333 e^2 b^2 ...
        not consistent with my results.
    """
    alpha=1/137.03
    hc = 197.32      # MeV.fm 
    hbar = 6.582e-22 # MeV.salpha=1/137.03
    
    E_g = 1.634
    T_half = 0.655e-12 # sec 
    W_g = np.log(2)/T_half # decay probability 
    
    lam = 2
    dfac = scipy.special.factorial2(2*lam+1)
    factor = alpha*hc*(8*np.pi*(lam+1))/(lam*dfac**2)/(hc)**(2*lam+1)/(hbar)
    B_E2 = W_g/(factor)/E_g**(2*lam+1)
    ans = "B(E2)={:e}  W(E2)={:e} ".format(B_E2, W_g) 
    return ans 