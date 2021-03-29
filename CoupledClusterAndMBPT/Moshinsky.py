import numpy
import numpy as np
import scipy 
from sympy.physics.sho import E_nl
from sympy import symbols
from sympy.physics.sho import R_nl
from sympy import var
from sympy.physics.wigner import wigner_9j
from sympy.physics.wigner import wigner_6j
from sympy.physics.wigner import wigner_3j

import scipy.special
from scipy.special import factorial
from functools import lru_cache

hbarc = 197.326968
amu = 931.4940954 # MeV
mass_p  = 938.272    
mass_n  = 939.5653
mass_N = (mass_p+mass_n)/2.0
alpha = 1.0/137.035

def Radial_HO3D(n,l,r,nu=0.5):
    """
    Radial Wave function of 3D S.H.O. 
    quantum number n corresponds to 
    here  E=(2n+l+3/2)omega 
    nu = m omega/(2hbar) 
        
    wave ~ R_nl(r) Y_lm 
         ~ u_nl(r)/r Y_lm 
         
    Default nu corresponds to dimensionless case.      
    """
    import scipy.special
    k=n;
    alpha=l+1./2;
    temp = r**l*np.exp(-nu*r**2)*scipy.special.assoc_laguerre(2*nu*r**2, k, alpha);
    norm = np.sqrt(2*(2*nu)**(l+3./2)*scipy.special.factorial(k)/scipy.special.gamma(k+l+3./2))
    return norm*temp

def b_HO3D(omega):
    """
    b=sqrt(hbarc/(m.omega)) for nucleon
    
    omega is in MeV units
    b is in fm units. 
    """ 
    return np.sqrt(1./(mass_N*omega))*hbarc 
    
def RacahW(j1,j2,j1p,j2p,J,k):
    return (-1)**(j1+j2+j1p+j2p)*float(wigner_6j(j1,j2,J,j2p,j1p,k,prec=64))
    
#   Moshinsky Bracket of n1=n2=0  
def Mosh_Afac(l1,l,l2,L,x):
    """ A(l_1,l,l_2,L,x) of eq.(64) in Moshinsky paper. 
        input arguments are all integers.
        
        Validity checck on x is done here.
    """
    xmin=max(abs(l-l1),abs(L-l2))
    xmax=min(l+l1,L+l2)
    if (x > xmax) or (x < xmin):
        return 0.0
    norm1=np.sqrt( factorial(l1+l+x+1)*factorial(l1+l-x)*factorial(l1+x-l)/factorial(l+x-l1) )
    norm2=np.sqrt( factorial(l2+L+x+1)*factorial(l2+L-x)*factorial(l2+x-L)/factorial(L+x-l2) )
    qmax = min(l+l1,L+l2)
    qmin = max(x,l1-l,l2-L)
    sums = 0.0
    for q in range(qmin,qmax+1):
        if (-1)**(l+q-l1)==1:
            term1=(-1)**((l+q-l1)/2)
            term2=factorial(l+q-l1)/factorial((l+q-l1)/2)/factorial( (l+l1-q)/2)
            term3=1./factorial(q-x)/factorial(q+x+1)
            term4=factorial(L+q-l2)/factorial((L+q-l2)/2)/factorial((L+l2-q)/2)
            sums= sums+term1*term2*term3*term4  
    return norm1*norm2*sums 

def Mosh_B00(l1,l2,lam,n,l,N,L):
    """
      Moshinsky Bracket with n1=n2=0
      lam is total angular momenta.
      
      Validity check on argument is not done 
    """
    if (l1<0 or l2 <0 or lam <0 or n <0 or l<0 or N<0 or L<0):
        return 0.0
    norm = factorial(l1)*factorial(l2)/factorial(2*l1)/factorial(2*l2)
    norm = norm*(2*l+1)*(2*L+1)/(2**(l+L))
    norm = norm*factorial(n+l)/factorial(n)/factorial(2*n+2*l+1)
    norm = norm*factorial(N+L)/factorial(N)/factorial(2*N+2*L+1)
    norm = np.sqrt(norm)*(-1)**(n+l+L-lam)
    
    xmin=max(abs(l-l1),abs(L-l2))
    xmax=min(l+l1,L+l2)
    if (xmax < xmin ) :
        return 0.0 
    sums=0.0
    for x in range(xmin,xmax+1):
        sums = sums +(2*x+1)*Mosh_Afac(l1,l,l2,L,x)*RacahW(l,L,l1,l2,lam,x) 
    return norm*sums 

## Moshinsky Bracket for non-zero n using Recurrence relation
def ME_negrsqr(n,l,N,L,lam,np,lp,Np,Lp,opt):
    """
    matrix element of 
    < nl,NL,lam|-r^2_{opt}|np,lp,Np,Lp,lam> 
    
    from Table 1 of Moshinsky paper.
    
    opt= 1 or 2
    """
    sgn =(-1)**(opt+1) # opt=2-> negative sign 
    out=0.0
    if( n<0 or l<0 or N<0 or L <0 or lam<0 or np<0 or lp<0 or Np<0 or Lp<0):
        return 0.0
    if ((np==n-1) and (lp==l) and (Np==N) and (Lp==L)):
        out= 0.5* numpy.sqrt(n*(n+l+0.5))  # note np is not numpy here 
    if ((np==n) and (lp==l) and (Np==N-1) and (Lp==L)):
        out= 0.5* numpy.sqrt(N*(N+L+0.5))    
    if ((np==n-1) and (lp==l+1) and (Np==N-1) and (Lp==L+1)):
        out= sgn*numpy.sqrt(n*N*(l+1)*(L+1))*(-1)**(lam+L+l)*RacahW(l,l+1,L,L+1,1,lam)
    if ((np==n-1) and (lp==l+1) and (Np==N) and (Lp==L-1)):
        out= sgn*numpy.sqrt(n*(N+L+0.5)*(l+1)*L)*(-1)**(lam+L+l)*RacahW(l,l+1,L,L-1,1,lam)    
    if ((np==n) and (lp==l-1) and (Np==N-1) and (Lp==L+1)):
        out= sgn*numpy.sqrt((n+l+0.5)*N*l*(L+1) )*(-1)**(lam+L+l)*RacahW(l,l-1,L,L+1,1,lam)    
    if ((np==n) and (lp==l-1) and (Np==N) and (Lp==L-1)):
        out= sgn*numpy.sqrt((n+l+0.5)*(N+L+0.5)*l*L )*(-1)**(lam+L+l)*RacahW(l,l-1,L,L-1,1,lam)    
    return out    

# use decorator to store recursive calculation results
@lru_cache(maxsize=None)          
def Moshinsky_Bracket(n,l,N,L,lam,n1,l1,n2,l2):
    """
    Moshinsky Braket 
    < nl,NL;lam | n1 l1,n2,l2;lam> 
    
    lam is total angular momentum 
    
    use recurrence relation of eq.(39) of Moshinsky paper. 
    
    Use always (n1 < n2) and (l1 < l2) 
               (n < N) and ( l < L ) by symmetry 
    """
    rho1=2*n+l+2*N+L
    rho2=2*n1+l1+2*n2+l2
    if (rho1 != rho2) or ( (-1)**(l+L) != (-1)**(l1+l2) ):
        return 0.0 
    
    if ( l1 > l2 ): # Use symmetry 
        return (-1)**(L-lam)*Moshinsky_Bracket(n,l,N,L,lam,n2,l2,n1,l1)
    if ( l > L) : # Use symmetry 
        return (-1)**(l1-lam)*Moshinsky_Bracket(N,L,n,l,lam,n1,l1,n2,l2)
    
    if (n1==0 and n2==0):
        return Mosh_B00(l1,l2,lam,n,l,N,L)
    if (n1 > 0): # reduce n1 value
        norm=1./np.sqrt(n1*(n1+l1+0.5))
        sums= 0.0
        cases=[ [n-1,l,N,L],[n,l,N-1,L]
               ,[n-1,l+1,N-1,L+1],[n-1,l+1,N,L-1]
               ,[n,l-1,N-1,L+1],[n,l-1,N,L-1]] 
        for case in cases:
            (npp,lpp,Npp,Lpp)=case
            term1=ME_negrsqr(n,l,N,L,lam,npp,lpp,Npp,Lpp,1)
            term2=Moshinsky_Bracket(npp,lpp,Npp,Lpp,lam,n1-1,l1,n2,l2)
            sums = sums +term1*term2 
        return norm*sums
    if (n2 > 0): #reduce n2 value
        norm=1./np.sqrt(n2*(n2+l2+0.5))
        sums=0.0
        cases=[ [n-1,l,N,L],[n,l,N-1,L]
               ,[n-1,l+1,N-1,L+1],[n-1,l+1,N,L-1]
               ,[n,l-1,N-1,L+1],[n,l-1,N,L-1]] 
        for case in cases:
            (npp,lpp,Npp,Lpp)=case
            term1=ME_negrsqr(n,l,N,L,lam,npp,lpp,Npp,Lpp,2)
            term2=Moshinsky_Bracket(npp,lpp,Npp,Lpp,lam,n1,l1,n2-1,l2)
            sums = sums +term1*term2 
        return norm*sums 
    else :
        return 0.0 
        
@lru_cache(maxsize=None)
def BMT_coef(n,l,npp,lpp,p):
    """
    Brody-Moshinsky-Talmi coefficient B(n,l,n',l',p)
    from eq.(23) of Brody-Moshinsky paper. 
    
    nonzero for  (l+l')/2<=p <= n+npp+(l+l')/2 
    """
    if (n<0 or l<0 or npp<0 or lpp<0 or p<0):
        return 0.0 
    
    term1 = (-1)**(p-(l+lpp)/2)*factorial(2*p+1)/(2**(n+npp)*factorial(p))
    term2 = np.sqrt( factorial(n)*factorial(npp)*factorial(2*n+2*l+1)*factorial(2*npp+2*lpp+1)
                    /factorial(n+l)/factorial(npp+lpp))
    alpha = int( max(0,p-(l+lpp)/2-npp) )
    beta  = int( min(n,p-(l+lpp)/2 ) )
    sums=0.0
    for k in range(alpha,beta+1):
        temp1 = factorial(l+k)*factorial(p-(l-lpp)/2-k)
        temp2 = factorial(k)*factorial(2*l+1+2*k)*factorial(n-k)
        temp3 = factorial(2*p-l+lpp+1-2*k)*factorial(npp-p+(l+lpp)/2+k)*factorial(p-(l+lpp)/2-k)        
        sums = sums+ temp1/temp2/temp3 
    return term1*term2*sums 

@lru_cache(maxsize=None)
def Talmi_integral(potential,p,b,Rmax=1000.):
    """
    Compute Talmi integral for given 
    potential name of V(r)
    p is order of Talmi-integral 
    b=sqrt(hbar/(m omega)) is a H.O. scale parameter in unit of fm 
    """
    from scipy.integrate import quad
    from scipy.special import gamma
    
    try: # for special cases 
        if (potential=='coulomb'):
            talmi = 1.43996/(np.sqrt(2)*b)*gamma(p+1)/gamma(p+3/2)
            return talmi        
    except:
        pass
    
    scale_factor=np.sqrt(2)*b
    def integrand(R):
        return  R**(2*p+2)*np.exp(-R**2)*potential(R*scale_factor)
    
    term1 = 2/gamma(p+3/2)
    term2 = quad(integrand,0,Rmax)[0]
    #print(' radial integration is done')
    return term1*term2    

@lru_cache(maxsize=None)
def TBME(a,b,c,d,J,T,potential=None,b_HO=None):
    """
    compute 2-body matrix elements of 
    < a, b| V| c, d>^{JT}_{nas}
    in H.O. basis 
    
    a(b,c,d) : a tuple (n,l,j)  j= integer or half integer 
    J        : total angular momentum 
    T        : total isospin 
    
    potential is given as a dictionary 
    
    potential[op_type]= (radial function)(r) 
    with 
    op_type  : 'SO'  spin singlet, spatial odd 
               'TO'  spin triplet, spatial odd
               'SE'  spin singlet, spatial even
               'TE'  spin triplet spatial even
               'LSO' spin-orbit spatial odd
               'LSE' spin-orbit spatial even
               'TNO' tensor spatial odd
               'TNE' tensor spatial even 
    """
    (na,la,ja) = a
    (nb,lb,jb) = b
    (nc,lc,jc) = c
    (nd,ld,jd) = d
    
    factor=2*np.sqrt((2*ja+1)*(2*jb+1)*(2*jc+1)*(2*jd+1) )
    if (a==b) : 
        factor=factor/np.sqrt(2.0)
    if (c==d) :
        factor=factor/np.sqrt(2.0)
    
    sums=0.0
    for S in [0,1]:
        Lam_min = max(abs(la-lb),abs(J-S))
        Lam_max = min(la+lb,J+S)
        Lamp_min = max(abs(lc-ld),abs(J-S))
        Lamp_max = min(lc+ld,J+S)
        for Lam in range(Lam_min,Lam_max+1):
            fac1 = float(wigner_9j(la,0.5,ja,lb,0.5,jb,Lam,S,J,prec=64))
            if (abs(fac1) < 1.0e-24): # skip zero
                continue  
            for Lamp in range(Lamp_min,Lamp_max+1): 
                fac2 = float(wigner_9j(lc,0.5,jc,ld,0.5,jd,Lamp,S,J,prec=64))
                if (abs(fac2) < 1.0e-24): #skip zero
                    continue 
                fac_Lam_Lamp_S = ( (-1)**(Lam-Lamp)*(2*Lam+1)*(2*Lamp+1)*(2*S+1)
                            *fac1*fac2 )
                #print('S=%3i Lam=%3i Lamp=%3i '%(S,Lam,Lamp)) #for test 
                temp = TBME_LS(a,b,c,d,J,T,Lam,Lamp,S,potential,b_HO)
                sums = sums + fac_Lam_Lamp_S * temp
                #print( 'fac_Lam_Lamp_S=%f tbme_LS=%f factor=%f'%(fac_Lam_Lamp_S,temp,factor  ))
    return sums*factor 

def TBME_LS(a,b,c,d,J,T,Lam=None,Lamp=None,S=None,potential=None,b_HO=None):
    """
    Matrix elemenent in LS scheme. 
    
    < (a,b) (Lam,S) J| V |(c,d) (Lamp,S) J > 
    """
    (na,la,ja) = a
    (nb,lb,jb) = b
    (nc,lc,jc) = c
    (nd,ld,jd) = d
    rho = 2*na+la+2*nb+lb 
    rhop= 2*nc+lc+2*nd+ld
    
    sums = 0.0
    for L in range(0,rho+1):
        l_min = abs(Lam-L)
        l_max = min(L+Lam,abs(rho-L))
        lp_min = abs(Lamp-L)
        lp_max = min(L+Lamp,abs(rhop-L))
        if (-1)**(l_min+L)-(-1)**(la+lb) != 0 : # parity check 
            l_min = l_min+1
        if (-1)**(lp_min+L)-(-1)**(lc+ld) !=0 :
            lp_min = lp_min +1
        for l in range(l_min,l_max+1,2):
            if (-1)**(l+S+T) == 1 : #  l+S+T=even case
                continue    
            for lp in range(lp_min,lp_max+1,2):
                if (-1)**(lp+S+T)==1 : # lp+S+T=even
                    continue 
                Nmin = 0
                Nmax = min( int((rho-l-L)/2),int((rhop-lp-L)/2) )
                for N in range(Nmin,Nmax+1):
                    n = int((rho-l-L-2*N)/2) 
                    np = int((rhop-lp-L-2*N)/2) 
                    fac_nlNL = Moshinsky_Bracket(n,l,N,L,Lam,na,la,nb,lb)*Moshinsky_Bracket(np,lp,N,L,Lamp,nc,lc,nd,ld)
                    if (abs(fac_nlNL) < 1.0e-24): #skip zero 
                        continue
                    jmin = max(abs(l-S),abs(lp-S)) 
                    jmax = min(l+S,lp+S)           
                    for j in range(jmin,jmax+1):
                        fac_j = (2*j+1)*float(wigner_6j(L,l,Lam,S,J,j,prec=64))*float(wigner_6j(L,lp,Lamp,S,J,j,prec=64))
                        if (abs(fac_j) < 1.0e-24): #skip zero 
                            continue                         
                        fac_tbme=TBME_rel(n,l,np,lp,S,j,T,potential,b_HO) 
                        sums = sums + fac_nlNL*fac_j *fac_tbme
                        #print('n=%3i l=%3i  np=%3i   lp=%3i j=%3i N=%3i  L=%3i  fac_nlNL=%f fac_j=%f fac_tbme=%f'%(n,l,np,lp,j,N,L,fac_nlNL,fac_j,fac_tbme) )
    return sums

def TBME_rel(n,l,np,lp,S,j,T,potential=None,b_HO=None):
    """
    potential matrix element in H.O. basis in relative coordinate  
    < n,l,S; j | V | np,lp, S; j> 
    
    potential is given as a dictionary radial function is given in function of r = (r1-r2)
    
    b_HO is a H.O. scale of basis 
    """
    sums=0.0
    if potential: # potential is specified
        for op in potential: # sum over operators
            factor=0.0
            if (op=='SO' and (S,T)==(0,0) and l==lp):
                sums = sums + RadialIntegral(n,l,np,lp, potential[op],b_HO )
            if (op=='TO' and (S,T)==(1,1) and l==lp):
                sums = sums + RadialIntegral(n,l,np,lp, potential[op],b_HO )
            if (op=='SE' and (S,T)==(0,1) and l==lp):
                sums = sums + RadialIntegral(n,l,np,lp, potential[op],b_HO )
            if (op=='TE' and (S,T)==(1,0) and l==lp):
                sums = sums + RadialIntegral(n,l,np,lp, potential[op],b_HO )
            if ((op=='LSO' and (S,T)==(1,1) and l==lp) 
                   or (op=='LSE' and (S,T)==(1,0) and l==lp) ):  # LS operator 
                factor=0.5*(j*(j+1)-l*(l+1)-S*(S+1))
                sums = sums + factor*RadialIntegral(n,l,np,lp, potential[op],b_HO )
            if ((op=='TNO' and (S,T)==(1,1) ) or (op=='TNE' and (S,T)==(1,0) ) ): # tensor operator
                factor=0.0
                if ( l==j and lp==j):
                    factor=2.0
                if (l==j-1 and lp==j-1):
                    factor=-2*(j-1)/(2*j+1)
                if (l==j+1 and l==j+1):
                    factor=-2*(j+2)/(2*j+1)
                if (l==j-1 and lp==j+1) or (l==j+1 and lp==j-1):
                    factor=6*np.sqrt(j*(j+1))/(2*j+1)
                sums = sums + factor*RadialIntegral(n,l,np,lp, potential[op],b_HO )    
    else : # potential is not given, assume constant 1
        if (n,l)==(np,lp) :
            return 1.0 
    return sums 

def RadialIntegral(n,l,npp,lpp, pot=None,b_HO=None):
    """
    Assumming known potential/function name 
    compute < n,l| pot(r)| np,lp> in H.O. basis with given b_HO. 
    
    Use Talmi-integral  
    < n,l| pot(r)| np,lp> = sum_p BMT_coef(...,p) *Talmi(...,p) 
    
    pot is given as a dictionary either as 
    pot={name: radial_function,parameter: additional_parameter}
    """
    sums = 0.0 
    p_min = int( (l+lpp)/2 ) # integer because l and lp have same parity
    p_max = p_min + n + npp
    #----special case of delta potential 
    # originally the angular part 1/(4 pi) should be part of TBME_rel as a factor
    # but is included here
    try:
        if (pot['potname']=='delta'):
            factor=pot['parameter']
            if (p_min==0): # ie. l=lp=0
                temp1=1./(4*np.pi)
                temp2=np.sqrt(factorial(2*n+1)*factorial(2*npp+1))/(2**(n+npp)*factorial(n)*factorial(npp))
                temp3=np.sqrt(2/np.pi)/b_HO**3
                return temp1*temp2*temp3*factor 
            else:
                return 0.0   
        if (pot['potname']=='Hrel'):
            # this case results are in units of (2 hbar omega/A) ?
            out=0.0
            if (n==npp and l==lpp):
                out = (2*n+l+1.5)
            return out/2  # overall 1/A factor?
        if (pot['potname']=='Trel'):
            # in units of (hbar omega/A) ?
            out=0.0
            if ((n,l)==(npp,lpp)):
                out = (2*n+l+1.5)
            if (n==npp-1 and l==lpp+2):
                out = np.sqrt(2*npp*+2*lpp+3)
            if (npp==n-1 and lpp==l+2):
                out = np.sqrt(2*n*+2*l+3)    
  
            if (n==npp and l==lpp+2):
                out = -np.sqrt((npp+lpp+1.5)*(npp+lpp+2.5))    
            if (n==npp and lpp==l+2):
                out = -np.sqrt((n+l+1.5)*(n+l+2.5))    
    
            if (n==npp+1 and l==lpp):
                out = np.sqrt((npp+1)*(npp+lpp+1.5))
            if (npp==n+1 and l==lpp):
                out = np.sqrt((n+1)*(n+l+1.5))
            return out/2 # overall 1/A factor ?   
    except:
        pass
    #----end special case 
    potname=pot['potname'] # external function 
    for p in range(p_min,p_max+1):
        coef = BMT_coef(n,l,npp,lpp,p)
        talmi = Talmi_integral(potname,p,b_HO) # "coulomb" is treated in special way. 
        sums = sums +coef*talmi 
    return sums     