"""
  Generate scattering data for optical potential  
  
  Variation in 
  (1) projectile and target energy 
      ZP, AP, ZT, AT, Elab 
  (2) optical potential parameters 
      (assume only real/imaginary volume terms.) 
      V, rv, av, W, rw, aw , rc      
  
  Use Fresco to compute elastic scattering cross section 
  
  It may be better to use ratio to Rutherford as an output. 
    
  To do: 
     decide how the data should be stored. 
     decide range of projectile/target/energy  
     decide range of optical potential parameters 
"""
import os
import shutil
import sys
import numpy as np
import numpy.linalg as npla

from subprocess import (call,Popen)
import myutil
from myutil import read_fresco_res ,clean_comm

import matplotlib.pyplot as plt

import run_fresco_v2 as frun 

element_names = ["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
		 "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
		 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
		 "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er",
		 "Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
		 "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
		 "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og",
		 "119","120","121","122","123","124","125","126","127","128","129","130"]
#=============================================================================
#def get_fresco_result():
#def chck_fresco_out(fname=''):
#def get_elastic_result_from_fresco_out(fname='_output.out'):
#------------------------------------------------------------------------------
#def run_fresco_from_input_txt(fresco_input_txt,
#                              fresco_path='fresco.exe',
#                              fresco_input_path='_test.in',
#                              fresco_output_path='_test.out'):

def elastic_input(AP,ZP,AT,ZT,Elab,V, rv, av, W, rw, aw, rc):
    """
    generate input text of Fresco for elastic scattering
    """
    default_input_txt=( 
    '{0:}+{1:} elastic scattering at Elab={2:}\n'+
    'NAMELIST\n'+
    '&FRESCO\n'+ 
    'HCM=0.1  RMATCH=20.0  JTMIN=0.0  JTMAX=200.0  THMIN=0.0\n' +
    'THMAX=180.0  THINC=1.0  CHANS=1  SMATS=2  XSTABL=1 \n'+
    'ELAB={2:}    / \n'+
    '&PARTITION  NAMEP={0:} MASSP={3:} ZP={4:} NAMET={1:} MASST={5:} ZT={6:} NEX=1 QVAL=0.0/\n'+ 
    '&STATES  JP = 0.0 EP = 0.0 BANDP = 1 JT = 0.0 ET = 0.0 BANDT = 1 CPOT = 1 / \n'+
    '&PARTITION / \n'+
    '&POT KP=1 TYPE=0 SHAPE=0 AP={3:} AT={5:} RC={7:} / \n'+
    '&POT KP=1 TYPE=1 SHAPE=0 P1={8:} P2={9:} P3={10:} P4={11:} P5={12:} P6={13:} P7=0.000 /\n'+ 
    '&POT / \n'+
    '&OVERLAP / \n'+
    '&COUPLING /\n')
    nameP = '{}'.format(AP)+element_names[ZP]
    nameT = '{}'.format(AT)+element_names[ZT]
    input_txt = default_input_txt.format(nameP,nameT,Elab,AP,ZP,AT,ZT,rc,V,rv,av,W,rw,aw) 
    return input_txt            

def elastic_run(AP,ZP,AT,ZT,Elab,V, rv, av, W, rw, aw, rc):
    """
    run fresco and get results 
    """
    fresco_input_txt = elastic_input(AP,ZP,AT,ZT,Elab,V, rv, av, W, rw, aw, rc)
    output = frun.run_fresco_from_input_txt(fresco_input_txt,
                                  fresco_path='fresco.exe',
                                  fresco_input_path='_test.in',
                                  fresco_output_path='_test.out')
    return output 
    
#-----------------------------------------------------------------------------
if __name__ == '__main__':    
    #main routine 
    for AP in [1]:
        for ZP in [1]:
            for AT in [12]:
                for ZT in [6]:
                    for Elab in [10.0,12.0]:
                        for V in [10.0,20.0,30.0,40.0]:
                            for rv in [1.2,1.25,1.3]:
                                for av in [0.6,0.65,0.7]:
                                    for W in [5.0,10.0,15.0,20.0]:
                                        for rw in [1.2,1.25,1.3]:
                                            for aw in [0.65]:
                                                for rc in [1.2]:
                                                    out = elastic_run(AP,ZP,AT,ZT,Elab,V, rv, av, W, rw, aw, rc)
                                                    plt.semilogy(out[0][:,0],out[0][:,1])
                                                
    
    
    