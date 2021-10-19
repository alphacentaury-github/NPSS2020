# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:21:28 2021

@author: Y-H Song

SPE reader
  : read single particle energy files and get the self energy
  For better low density region extrapolation,
  change variables from (k,rho,iso) -> (k,rho_p,rho_n)
  
  then sigma(k, rho_p=0,rho_n->0)

#list_spe = 1, 2, 2I
#list_pn =  0, 1  (proton,neutron)
#list_dens = 0.01- 0.32 fm^-3 in 0.01 steps
#list_iso = (n_n-n_p)/(n_n+n_p) = 0-0.6 in 0.1 steps
#list_momentum = 0.02-3.92 fm^-1 in 0.1 steps

#note the spe1 include kinetic energy k^2/2M
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from collections import Sequence 

dd = "./SPE_holt/" # file path
hc=197.3269 # MeV fm
amu=931.494 #MeV
mp=938.27/hc #fm^-1
mn=939.5653/hc #fm^-1
mN =0.5*(mp+mn) # average nucleon mass fm^-1 

#--Holt data grids
iso_list = np.arange(0,0.7,0.1) # 7 values
dens_list = np.arange(0.01,0.33,0.01) # 32 values
mom_list = np.arange(0.02,4.02,0.1) # 40 values

def gauleg(a,b,n):
    """
    Gauss Legendre quadrature points/weights in (a,b) of degree 2^n-1

    Parameters
    ----------
    a : float
        lower bound
    b : float
        upper bound
    n : integer
        degree
    Returns
    -------
    x : ndarray
       sample points
    w : ndarray
       weights
    """
    x, y = np.polynomial.legendre.leggauss(n)
    x = ((b-a)*x+b+a)*0.5
    w = (b-a)*0.5* y
    return (x,w)

# Shape forms
def shape_WS(r,R0=1.0,a0=1.0,V0=1.0,PWR=1.0,IDER=0):
    """
    return Woods-Saxon shape or derivative

    IDER==0 case:
        f(r)=V0/(1+exp((r-R0)/a0) )^PWR

    IDER!=0 case:
        df(r)/dr

    """
    EXPR = np.exp((r-R0)/a0)
    if IDER ==0:
        return V0/(1.0+EXPR)**PWR
    else:
        return -PWR*V0*EXPR/(a0*(1.0+ EXPR)**(PWR+1))

def normalize_density(f_rho,norm,r_range=(0.,20.),num_quad=70):
    """
    compute normalization factor 'C'
    so that the

    norm = 4pi\int dr r^2 f_rho(r)*C

    f_rho : density function of r
    norm  : normalization factor
    """
    (xr,wr) = gauleg(r_range[0],r_range[1],num_quad)
    rho = f_rho(xr)
    vol = np.sum(wr*rho*xr**2)*4.0*np.pi
    norm_factor = norm/vol
    return norm_factor

def Sao_Paulo_density(AP,ZP):
       """
       return WS potential parameter of density
       based on Sao-Paulo parametrization
       """
       NP = AP - ZP
       # rho0=0.091 # fm^{-3}
       R_p = 1.81*ZP**(1./3.)-1.12
       a_p = 0.47-0.00083*ZP
       R_n = 1.49*NP**(1./3.)-0.79
       a_n = 0.47+0.00046*NP
       #---proton density
       f_rho = lambda r : shape_WS(r,R0=R_p,a0=a_p,V0=1.0)
       C_p = normalize_density(f_rho, ZP,r_range=(0.,20.),num_quad=70)
       rho_p = lambda r : shape_WS(r,R0=R_p,a0=a_p,V0=C_p)
       #---neutron density
       f_rho = lambda r : shape_WS(r,R0=R_n,a0=a_n,V0=1.0)
       C_n = normalize_density(f_rho, NP,r_range=(0.,20.),num_quad=70)
       rho_n = lambda r : shape_WS(r,R0=R_n,a0=a_n,V0=C_n)
       return (rho_p,rho_n)

def read_spe_file(pn,dens,iso):
    # pn = 0 or 1
    # dens = 100*rho
    # iso = 10*delta

    #check input
    if not(pn in [0,1]):
        raise ValueError("Error: input pn must be 0 or 1")
    if not(dens in range(1,33)):
        raise ValueError("input dens must be within [1,32]")
    if not(iso in range(0,7)):
        raise ValueError("input iso must be within [0,6]")
    #read file
    fname1=dd+"spe1_pn{:}_n{:}_iso{:}.dat".format(pn,dens,iso)
    fname2=dd+"spe2_pn{:}_n{:}_iso{:}.dat".format(pn,dens,iso)
    fname2i=dd+"spe2I_pn{:}_n{:}_iso{:}.dat".format(pn,dens,iso)
    out1 = np.loadtxt(fname1)
    out2 = np.loadtxt(fname2)
    out3 = np.loadtxt(fname2i)

    out = np.column_stack((np.ones(len(out1[:,0]))*dens*0.01,
                           np.ones(len(out1[:,0]))*iso*0.1,
                           out1[:,0],
                           out1[:,1]-out1[:,0]**2/(2*mN)*hc,
                           out2[:,1],
                           out3[:,1]
                           ) )
    return out

def read_all_spe_file():
    # read all files and return all table as a dictionary
    # alldat2[np][dens-index][iso-index][k-index] = [dens, iso, k, spe1,spe2,spe2I]
    alldat2={}
    alldat2[0]=[]
    alldat2[1]=[]
    for pn in [0,1]:
        temp2=[]
        for dens in range(1,33):
            temp=[]
            for iso in range(0,7):
                try:
                    out = read_spe_file(pn,dens,iso)
                    # alldat2[pn].append(out)
                    temp.append(out)
                except ValueError:
                    raise ValueError("Error! reading speX_pn{:}_n{:}_iso{:}.dat".format(
                        pn,dens,iso) )
            temp2.append(temp)
        alldat2[pn]=temp2    
    return alldat2

def add_derivative_spe(alldat2):
    """
    from given input data 
    
    alldat2[np][dens-index][iso-index][k-index] = [dens, iso, k, spe1,spe2,spe2I]
    
    add columns of derivative values 
    
    out[np][dens-index][iso-index][k-index] = [dens, iso, k, spe1,spe2,spe2I,dV/dk]
    """
    step = mom_list[1] - mom_list[0]
    dVdk = np.zeros(len(mom_list))
    for pn in [0,1]:
        for dens in range(len(dens_list)):
            for iso in range(len(iso_list)):
                dVdk = np.zeros(len(mom_list))
                for k in range(len(mom_list)):
                    (rho,delta,kk,spe1,spe2,spe2i) = alldat2[pn][dens][iso][k]
                    # compute dimensionless d(spe1+spe2)/dk 
                    V = (spe1+spe2)
                    if k==0 : 
                        V_plus = (alldat2[pn][dens][iso][k+1][3]+alldat2[pn][dens][iso][k+1][4])
                        dVdk[k] = (V_plus-V)/(step*hc)  #Careful for dimension
                        #--alternative form--
                        #V_2p = (alldat2[pn][dens][iso][k+2][3]+alldat2[pn][dens][iso][k+2][4])
                        #V_p  = (alldat2[pn][dens][iso][k+1][3]+alldat2[pn][dens][iso][k+1][4])
                        #dVdk[k] = (-0.5*V_2p+2*V_p-1.5*V)/(step*hc) 
                    elif k== len(mom_list)-1: 
                        V_minus = (alldat2[pn][dens][iso][k-1][3]+alldat2[pn][dens][iso][k-1][4])
                        dVdk[k] = (V-V_minus)/(step*hc)
                        #--alternative form--
                        #V_2m = (alldat2[pn][dens][iso][k-2][3]+alldat2[pn][dens][iso][k-2][4])
                        #V_m  = (alldat2[pn][dens][iso][k-1][3]+alldat2[pn][dens][iso][k-1][4])
                        #dVdk[k] = (0.5*V_2m-2*V_m+1.5*V)/(step*hc) 
                    else : 
                        V_plus = (alldat2[pn][dens][iso][k+1][3]+alldat2[pn][dens][iso][k+1][4])
                        V_minus = (alldat2[pn][dens][iso][k-1][3]+alldat2[pn][dens][iso][k-1][4])
                        dVdk[k] = (V_plus-V_minus)/(2*step*hc)
                alldat2[pn][dens][iso] = np.column_stack((alldat2[pn][dens][iso],dVdk ))
    return alldat2  

def read_all_spe_files():
    alldat2 = read_all_spe_file() 
    data = add_derivative_spe(alldat2)
    return data 

#----2021-10-18 working here 
    
def prepare_interp_data(alldat):
    # re-arrange data points and values for multi-dimensional interpolation
    # alldat is an output from read_all_spe_files 
    p_data  = alldat[0] #pn=0
    n_data  = alldat[1] #pn=1

    newdat = []
    points = (dens_list, iso_list, mom_list)
    newdat.append(points) 
    sig_p_re = np.zeros( (len(dens_list),len(iso_list),len(mom_list)) )
    sig_p_im = np.zeros( (len(dens_list),len(iso_list),len(mom_list)) )
    sig_n_re = np.zeros( (len(dens_list),len(iso_list),len(mom_list)) )
    sig_n_im = np.zeros( (len(dens_list),len(iso_list),len(mom_list)) )
    dsig_p     = np.zeros( (len(dens_list),len(iso_list),len(mom_list)) )
    dsig_n     = np.zeros( (len(dens_list),len(iso_list),len(mom_list)) )
    for dens in range(len(dens_list)):
        for iso in range(len(iso_list)):
            sig_p_re[dens][iso][:] = p_data[dens][iso][:,3]+ p_data[dens][iso][:,4]
            sig_p_im[dens][iso][:] = p_data[dens][iso][:,5]
            dsig_p[dens][iso][:] = p_data[dens][iso][:,6]
            sig_n_re[dens][iso][:] = n_data[dens][iso][:,3]+ n_data[dens][iso][:,4]
            sig_n_im[dens][iso][:] = n_data[dens][iso][:,5]
            dsig_n[dens][iso][:] = n_data[dens][iso][:,6]
            
    return (points, sig_p_re, sig_p_im, dsig_p, sig_n_re, sig_n_im,dsig_n)

def SE_interp(rho,delta,k,mode=None,data=None):
    """
    self energy interpolation 
    mode=0 : Re(Sigma_p)
        =1 : Im(Sigma_p)
        =2 : d/dk Re(Sigma_p)
        =3 : Re(Sigma_n)
        =4 : Im(Sigma_n)
        =5 : d/dk Re(Sigma_n)
    However, the extrapolation does not work!!
    
    upto here no sign change from the data. 
    """
    points=data[0]
    values=data[mode+1]
    
    if (delta<0 and abs(delta)<=0.6 ): #negative delta case 
        if mode in [0,1]: #proton case
            values=data[mode+1+3]
            delta = -delta 
        else :  #neutron case   
            values=data[mode+1-3]
            delta = -delta
    if (rho< dens_list[0]):
        out = (interpn(points,values,(dens_list[0],delta,k))/dens_list[0]*rho )         
    else:
        out = interpn(points,values,(rho,delta,k)) 
    return out 
    

#===========================================================================
if __name__ == '__main__':
    pn = 1
    dens = 30 # 0.01
    iso = 0
    out = read_spe_file(pn,dens,iso)
    plt.title("pn={},dens={},iso={}".format(pn,dens*0.01,iso*0.1))
    plt.plot(out[:,2],out[:,3],'.',label='LO')
    plt.plot(out[:,2],out[:,4],'.',label='NLO real')
    plt.plot(out[:,2],out[:,5],'.',label='NLO imag')
    plt.plot(out[:,2], out[:,3]+out[:,4],'.',label='real LO+NLO')
    plt.xlabel('k [fm^-1]')
    plt.ylabel('MeV')
    plt.legend()

    alldat = read_all_spe_files()
    data = prepare_interp_data(alldat)

    print('==multi-dimensional linear interpolation==')
    kk = np.linspace(0.1,3.2,100)
    plt.plot(kk,SE_interp(dens*0.01,iso*0.1,kk,mode=0,data=data))
    plt.plot(kk,SE_interp(dens*0.01,iso*0.1,kk,mode=1,data=data))
    plt.plot(kk,SE_interp(dens*0.01,iso*0.1,kk,mode=2,data=data))
    plt.plot(kk,SE_interp(dens*0.01,iso*0.1,kk,mode=3,data=data))
    #plt.plot(kk,SE_interp(dens*0.01,-iso*0.1,kk,mode=3,data=interp_data))
    plt.legend()

    def get_eff_mass(rho,delta,k,charge=0,data=None):
        """
        return M*/M          
        """
        if charge==0:
            mode = 5 # d/dk Re Sig_n
        elif charge==1:
            mode = 2 # d/dk Re Sig_p        
        dVk =  SE_interp(rho,delta,k,mode = mode,data=data)   
        out = 1./(1 + mN/k*dVk)
        #if (out < 0.).any() :
        #    raise ValueError('Error: negative effective mass')
        return out 

    print('==test effective mass==')
    rho = 0.3; delta = 0;
    plt.plot(kk,SE_interp(rho,delta,kk,mode=0,data=data))
    emass = get_eff_mass(rho,delta,kk,charge=0,data=data)
    plt.plot(kk,emass*10)
    print(r'$\Sigma$={}'.format(emass)  )

    # #----test extrapolation
    # # (1) scipy.interpn cannot extrapolate except for constant value
    # # (2) gaussian process can not be used for ~9000 data points
    # # (3) neural network for extrapolation?
    # from sklearn.gaussian_process import GaussianProcessRegressor
    # from sklearn.gaussian_process.kernels import RBF, WhiteKernel
    # from sklearn.neural_network import MLPRegressor
    # X = []
    # Y = []
    # mode = 0
    # points=interp_data[0]
    # values=interp_data[mode+1]
    # for i in range(len(points[0])):
        # for j in range(len(points[1])):
            # for k in range(len(points[2])):
                # X.append([points[0][i],points[1][j],points[2][k]])
                # Y.append(values[i][j][k])
    # X_train= np.array(X)
    # Y_train= np.array(Y)
    # #kernel = 1.0 * RBF(length_scale=100.0, length_scale_bounds=(1e-2, 1e3))+ WhiteKernel(noise_level=1, noise_level_bounds=(1e-10, 1e+1))
    # #gp = GaussianProcessRegressor(kernel=kernel,alpha=0.0).fit(X, Y) # not practical
    # #regr = MLPRegressor(hidden_layer_sizes=(100,50),random_state=1, max_iter=900).fit(X, Y)

    # #----optical potential in nuclear matter
    def opt_pot_NM(E_kin,rho,delta,charge=0,data=data):
        """
         U_{np}(E,rho,delta)=V(E,rho,delta)+ i W(E,rho,delta)

        # E_kin in MeV unit
        # rho in fm^-3 unit
        
        no sign change in imaginary part up to here 
        """
        k = np.sqrt(2*mN*E_kin/hc) #fm^-1

        if charge==1:
            V = SE_interp(rho,delta,k,mode=0,data=data)
            m_ratio = get_eff_mass(rho,delta,k,charge,data=data)
            W = SE_interp(rho,delta,k,mode=1,data=data)*m_ratio
        elif charge==0:
            V = SE_interp(rho,delta,k,mode=3,data=data)
            m_ratio = get_eff_mass(rho,delta,k,charge,data=data)
            W = SE_interp(rho,delta,k,mode=4,data=data)*m_ratio
        return (V,W)

    def opt_pot_NM2(E_kin,rho,delta,charge=None,data=None):
        # for the case of one input is an array 
        if isinstance(rho,(Sequence, np.ndarray)):
            v=[];w=[];
            for i in rho:
                (a,b) = opt_pot_NM(E_kin,i,delta,charge,data)
                v.append(a[0])
                w.append(b[0])
            return (v,w) 
        elif isinstance(delta,(Sequence,np.ndarray)):    
            v=[];w=[];
            for i in delta:
                (a,b) = opt_pot_NM(E_kin,rho,i,charge,data)
                v.append(a[0])
                w.append(b[0])
            return (v,w) 
        else: 
            return opt_pot_NM(E_kin,rho,delta,charge,data)


    # #---test NM optical
    plt.figure()
    plt.title('neutron-NM optical E=85 MeV')
    plt.xlabel('rho')
    plt.ylabel('MeV')
    x = np.arange(0.01,0.3,0.01) #density 
    (v,w) = opt_pot_NM2(85, x, 0,charge=0,data=data)
    plt.plot( x, v,label='V(delta=0)')
    plt.plot( x, w,label='W(delta=0)')
    (v,w) = opt_pot_NM2(85, x, 0.6,charge=1,data=data)
    plt.plot( x, v,label='V(delta=0.6)')
    plt.plot( x, w,label='W(delta=0.6)')
    plt.legend()
        
    # #-- delta dependence plot 
    plt.figure()
    plt.title('neutron-NM optical E=85 MeV')
    plt.xlabel('delta')
    plt.ylabel('MeV')
    x = np.arange(0,0.7,0.1) #delta array 
    (v,w) = opt_pot_NM2(85, 0.01, x,charge=0,data=data)
    plt.plot( x, v,label='V(rho=0.01)')
    plt.plot( x, w,label='W(rho=0.01)')
    (v,w) = opt_pot_NM2(85, 0.32, x,charge=0,data=data)
    plt.plot( x, v,label='V(rho=0.32)')
    plt.plot( x, w,label='W(rho=0.32)')
    plt.legend() 
    
    # # energy dependence 
    plt.figure() 
    plt.title('neutron-NM optical delta=0')
    plt.xlabel('E [MeV]')
    plt.ylabel('MeV')
    x = np.arange(0.04,3.9,0.1) # fm^-1 unit 
    xe = x**2/(2*mN)*hc # MeV unit 
    (v,w) = opt_pot_NM2(xe, 0.3, 0 ,charge=0,data=data)
    plt.plot( xe, v,label='V(rho=0.3,delta=0)')
    plt.plot( xe, w,label='W(rho=0.3,delta=0)')
    (v,w) = opt_pot_NM2(xe, 0.02, 0 ,charge=0,data=data)
    plt.plot( xe, v,label='V(rho=0.02,delta=0)')
    plt.plot( xe, w,label='W(rho=0.02,delta=0)')
    plt.legend() 
    
    
    # #---optical potential for finite nuclei
    def opt_pot_LDA(E_kin,f_dens_p,f_dens_n,charge=0,r_range=[0,20,0.5],data=None):
        """
         f_dens_p(r), f_dens_n(r) is an interpolating function
        """
        rmin,rmax,rstep = r_range
        R = np.arange(rmin,rmax+rstep,rstep)
        U = []
        for r in R:
            rho_p = f_dens_p(r)
            rho_n = f_dens_n(r)
            rho = rho_p + rho_n
            rho_np = (rho_n -rho_p)
            delta = rho_np/rho
            (v,w) = opt_pot_NM2(E_kin,rho,delta,charge,data)
            U.append(v+1j*w)
        return (R,np.array(U)) 
    # #---test with Sao-Paulo density
    (rho_p,rho_n) = Sao_Paulo_density(40, 20)
    R=np.arange(0.,12.0,0.1)
    plt.figure()
    plt.title('Sao-Paulo density 40Ca')
    plt.plot(R,rho_p(R),label='p')
    plt.plot(R,rho_n(R),label='n')
    plt.plot(R, (rho_n(R)-rho_p(R))/(rho_n(R)+rho_p(R))*0.1,label='delta*0.1')
    plt.xlim(0,12)
    plt.legend()
    plt.figure()
    (R,U) = opt_pot_LDA(85,rho_p,rho_n,charge=0,r_range=[0,12.0,0.1],data=data)
    plt.plot(R,np.real(U),label='real')
    plt.plot(R,np.imag(U),label='imag')
    plt.xlim(0,12)
    plt.legend()












