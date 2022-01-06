! Created on 8/27/2012
! Created by Luke Titus
! This program solves differential equations of the form y''+f(r)y(r)+S(r)=0
!   using the Numerov Method


  Subroutine nm(L,nmax,h,vpot,C,E,logd,source,y)
  implicit none     
     
  Real(8)::h,E,hp,C,maximum,Vcent,R
  Complex*16::vpot(nmax),logd,source(nmax),f(nmax),y(nmax)
  integer::nmax,n,nmin,L,nmax_temp,Source_Index

! INITIAL CONDITIONS
  hp=h*h 
  y(1)=(h)**(L+1.0)
  y(2)=(2.0*h)**(L+1.0)
  f(1)=C*(E-vpot(1))-L*(L+1.0)/(1.0*h)**2
  f(2)=C*(E-vpot(2))-L*(L+1.0)/(2.0*h)**2   
  nmin=2

! INITIAL CONDITIONS FOR LARGE L          
  if (L>0) then
    nmin = int(2*L)
    do n=1,(nmin-2)
      y(n)=0
      f(n)=0
    end do
    y(nmin-1)=((nmin-1.0)*h)**(L+1.0)
    y(nmin)=(nmin*h)**(L+1.0)
    f(nmin-1)=C*(E-vpot(nmin-1))-L*(L+1.0)/((nmin-1.0)*h)**2.0
    f(nmin)=C*(E-vpot(nmin))-L*(L+1)/(nmin*h)**2.0
  end if


! SOLVE FOR THE WAVE FUNCTION
  do n=nmin,(nmax-1)
    f(n+1)=C*(E-vpot(n+1))-L*(L+1)/((h*(n+1))**2)
    y(n+1)=(y(n)*(2.0-(5.0*hp/6.0)*f(n))-y(n-1)*(1+(hp/12.0)*f(n-1))-(hp/12.0) &
           *(source(n-1)+10.0*source(n)+source(n+1)))/(1+(hp/12.0)*f(n+1))
  end do


! CALCULATE THE LOGARITHMIC DERIVATIVE
  LOGD=147.0*Y(NMAX)-360.0*Y(NMAX-1)+450.0*Y(NMAX-2)-400.0*Y(NMAX-3)+225.0*Y(NMAX-4)-72.0*Y(NMAX-5)+10.0*Y(NMAX-6)
  LOGD=LOGD/(60.0*Y(NMAX)*H)
    

return
end subroutine nm
