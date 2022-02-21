# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 10:54:13 2019

@author: notebook3
"""

#   format    :  a1,i3,i5,i5,i5,1x,a3,a4,1x,f13.5,f11.5,f11.3,f9.3,1x,a2,f11.3,f9.3,1x,i3,1x,f12.5,f11.5
#                cc NZ  N  Z  A    el  o     mass  unc binding unc     B  beta  unc    atomic_mass   unc
#
##a1,   i3,       i5,     i5,       i5,        1x,   a3,        a4, 
# cc    NZ         N       Z        A                EL          O
#(ll[0], ll[1:4], ll[4:9], ll[9:14], ll[14:19],ll[19],ll[20:23],ll[23:27] )
##1x,      f13.5  ,f11.5     ,f11.3      ,f9.3     ,1x   ,  a2,   
#          mass    unc        binding     unc               B      
#(ll[27],ll[28:41], ll[41:52], ll[52:63], ll[63:72], ll[72], ll[73:75] )
##f11.3,    f9.3,      1x,   i3,       1x,       f12.5,      f11.5
# beta      unc              atomic_mass                    unc
#(ll[75:86], ll[86:95],ll[95],ll[96:99], ll[99], ll[100:112],ll[112:123] ) 
# * represent uncalculable 
# # instead of . indicate non-experimental estimate 

# convert integer to array 
def integer_to_bin_array(N,dim_b=8):
    """ Convert integers into dim_b bits  binary
    """
    import numpy as np
    bin_str = bin(N)[2:] # to binary string 
    if len(bin_str) > dim_b :
        raise ValueError('{} is too large to be {} bits'.format(N,dim_b))
    bin_str=bin_str.zfill(dim_b) #padding zero 
    bin_array=[float(j) for j in bin_str] 
    return bin_array

def read_mass16(filename='mass16.txt',include_guess=False):
  """
    Read AME mass data file and return a dictionary 
    
    include_guess=False : only make a dictionary with measured data 
                 =True  : dictionary includes guessed data 
  """
  import numpy as np
  import pandas as pd
  def check_number(x):
    # if string is number, return real number 
    # if not, retrun False 
    try:
        y= float(x)
        return (y , 1)
    except :
        return (x, 0) 

  file = open(filename,'r')
  llines = file.readlines()
  file.close() 
  nl_comment=38 # end of commenting lines 
        
  AME16={'NZ':[],'N':[],'Z':[],
       'A':[],'EL':[],'O':[],
       'mass_excess':[],'mass_excess_measured':[],'mass_excess_unc':[],
       'BE/A':[],'BE/A_measured':[],'BE/A_unc':[],
       'B':[],
       'beta':[],'beta_measured':[],'beta_unc':[],
       'atomic_mass':[],'atomic_mass_measured':[],'atomic_mass_unc':[],
       'N_bits':[],'Z_bits':[] 
        }
             
  for ll in llines[nl_comment+1: ]:
    (NZ,dummy)= check_number(ll[1:4])
    (N, dummy)= check_number(ll[4:9])
    (Z, dummy)= check_number(ll[9:14])
    (A,dummy) = check_number(ll[14:19])
    EL = ll[20:23]
    O  = ll[23:27] # ?
    (mass_excess,mass_excess_measured) = check_number(ll[28:41])  # mass_excess
    (mass_excess_unc,dummy)= check_number(ll[41:52])
    (binding,binding_measured) = check_number(ll[52:63])
    (binding_unc,dummy) = check_number(ll[63:72])
    B = ll[73:75] #?
    (beta,beta_measured) = check_number(ll[75:86])
    (beta_unc,dummy) = check_number(ll[86:95])
    (amass,dummy) = check_number(ll[96:99])
    (amass2,amass_measured) = check_number(ll[100:112])
    if amass_measured:
        amass = amass+amass2/10**6 
    else: 
        amass = str(amass)+' '+amass2.strip() 
    (amass_unc,dummy) = check_number(ll[112:122])
    N_bits = integer_to_bin_array(int(N))
    Z_bits = integer_to_bin_array(int(Z))
    if (not include_guess) and (not amass_measured): # only take measured values 
      pass 
    else : 
      AME16['NZ'].append(NZ)
      AME16['N'].append(N)
      AME16['Z'].append(Z)
      AME16['A'].append(A)
      AME16['EL'].append(EL)
      AME16['O'].append(O)
      AME16['mass_excess'].append(mass_excess)
      AME16['mass_excess_measured'].append(mass_excess_measured)
      AME16['mass_excess_unc'].append(mass_excess_unc)
      AME16['BE/A'].append(binding)
      AME16['BE/A_measured'].append(binding_measured)
      AME16['BE/A_unc'].append(binding_unc)
      AME16['B'].append(B)
      AME16['beta'].append(beta)
      AME16['beta_measured'].append(beta_measured)
      AME16['beta_unc'].append(beta_unc)
      AME16['atomic_mass'].append(amass)
      AME16['atomic_mass_measured'].append(amass_measured)
      AME16['atomic_mass_unc'].append(amass_unc)  
      AME16['N_bits'].append(N_bits)  
      AME16['Z_bits'].append(Z_bits)    
  return AME16      
# make a dataframe 
# AME16=read_mass16();test_df=pd.DataFrame.from_dict(AME16)