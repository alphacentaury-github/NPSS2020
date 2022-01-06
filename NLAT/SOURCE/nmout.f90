!    CREATED BY: LUKE TITUS
!    LAST UPDATED: 10/10/2014
!    THIS PROGRAM SOLVES THE BOUND STATE EQUATION OF THE FORM y''-[K^2+f(r)]y(r)+S(r)=0 USING THE NUMEROV METHOD
!    OUTPUT...
!      BOUND STATE ENERGY  
!      BOUND STATE WAVE FUNCTION

     Subroutine nm_out(L,nmax,h,vpot,C,E,logd,integral,y,nmatch)
 
     Implicit None
!!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!!     
     Real*8::h,E,hp,C,maximum,Vcent,vpot(nmax),y(nmax),logd,integral(nmax),f(nmax)
     integer::nmax,n,nmin,nmatch,L

!    INITIAL CONDITIONS
     hp=h*h 
     y(1)=(h)**(L+1.0)
     y(2)=(2.0*h)**(L+1.0)
     f(1)=C*(E-vpot(1))-L*(L+1)/(1.0*h)**2
     f(2)=C*(E-vpot(2))-L*(L+1)/(2.0*h)**2   
     nmin=2

!    INITIAL CONDITIONS FOR LARGE L          
     if (L>0) then
       nmin = int(2*L)
       do n=1,(nmin-2)
          y(n)=0
          f(n)=0
       end do
       y(nmin-1)=((nmin-1.0)*h)**(L+1.0)
       y(nmin)=(nmin*h)**(L+1.0)
       f(nmin-1)=C*(E-vpot(nmin-1))-L*(L+1)/((nmin-1.0)*h)**2
       f(nmin)=C*(E-vpot(nmin))-L*(L+1)/(nmin*h)**2
     end if

!    SOLVE FOR THE WAVE FUNCTION
     do n=nmin,(nmax-1)
       f(n+1)=C*(E-vpot(n+1))-L*(L+1)/((h*(n+1))**2)
       y(n+1)=(y(n)*(2.0-(5.0*hp/6.0)*f(n))-y(n-1)*(1+(hp/12.0)*f(n-1))-(hp/12.0)*(integral(n-1)+10.0*integral(n)+integral(n+1))) &
              /(1+(hp/12.0)*f(n+1))
     end do

!    CALCULATE THE LOGARITHMIC DERIVATIVE
     LOGD=147.0*Y(nmatch)-360.0*Y(nmatch-1)+450.0*Y(nmatch-2)-400.0*Y(nmatch-3)+225.0*Y(nmatch-4)-72.*Y(nmatch-5)+10.*Y(nmatch-6)
     LOGD=LOGD/(60.0*Y(nmatch)*H)
    
     return
     end subroutine nm_out
