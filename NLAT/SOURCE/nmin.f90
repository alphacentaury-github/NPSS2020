!    CREATED BY: LUKE TITUS
!    LAST UPDATED: 10/10/2014
!    THIS PROGRAM SOLVES THE BOUND STATE EQUATION OF THE FORM y''-[K+f(r)]y(r)+S(r)=0 USING THE NUMEROV METHOD
!    OUTPUT...
!      BOUND STATE ENERGY  
!      BOUND STATE WAVE FUNCTION

     Subroutine nm_in(L,nmax,h,vpot,SchEqConst,E,logd,source,y,nmatch,const,eta)
     
     Implicit None
!!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!!
     real(8)::h,E,hp,SchEqConst,maximum,Vcent,y1,y2,y3,witt1,witt2,witt3,anc,k,const,eta
     real(8)::y(nmax),logd,source(nmax),f(nmax),vpot(nmax)
     integer::nmax,n,nmatch,calc,L

!    INITIAL CONDITIONS
     hp=h*h
     k=sqrt(-SchEqConst*E)
     f(nmax)=SchEqConst*(E-vpot(nmax))-L*(L+1)/(nmax*h)**2
     f(nmax-1)=SchEqConst*(E-vpot(nmax-1))-L*(L+1)/((nmax-1)*h)**2  

!    SOLVE FOR THE WAVE FUNCTION
     do n=(nmax-1),2,-1
       f(n-1)=SchEqConst*(E-vpot(n-1))-L*(L+1)/((h*(n-1))**2)
       y(n-1)=(y(n)*(2.0-(5.0*hp/6.0)*f(n))-y(n+1)*(1+(hp/12.0)*f(n+1))-(hp/12.0)*(source(n-1)+10.0*source(n)+source(n+1))) &
              /(1+(hp/12.0)*f(n-1))       
     end do

!    CALCULATE THE LOGARITHMIC DERIVATIVE
     LOGD=147.0*Y(nmatch)-360.0*Y(nmatch-1)+450.0*Y(nmatch-2)-400.0*Y(nmatch-3)+225.0*Y(nmatch-4)-72.*Y(nmatch-5)+10.*Y(nmatch-6)
     LOGD=LOGD/(60.0*Y(nmatch)*H)
    
     return
     end subroutine nm_in
