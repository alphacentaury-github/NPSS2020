  subroutine nonlocal_binding_potential(BoundParameters,Unl,r,rprime,UlocNL,StepSize)
  implicit none

  integer::WhatSystem,WhatPot,L
  real(8)::Vnucl,Wnucl,Vsurf,Wsurf,Vspin,Wspin,Vcoul,Z1,Z2,StepSize,hbarc,A,r,rprime,p,Unl,Ulocnl
  real(8)::Vv,Rv,av,Wv,Rwv,awv,Vd,Rvd,avd,Wd,Rwd,awd,Vso,Rso,aso,Wso,Rwso,awso,Rcoul,J,SpinFragment
  real(8)::BoundParameters(1:15,1:100)
  complex*16::i


  i=cmplx(0.0,1.0)
  hbarc=197.3269718

  Z1=BoundParameters(1,3)
  Z2=BoundParameters(1,4)
  SpinFragment=BoundParameters(1,5)
  A=BoundParameters(1,2)
  L=nint(BoundParameters(1,9))
  J=BoundParameters(1,10)
  WhatSystem=nint(BoundParameters(2,1))
  WhatPot=nint(BoundParameters(2,2))

! WhatSystem = 1 = Deuteron Bound State
! WhatSystem = 2 = Nucleon Bound State
  if (WhatSystem==2) then

    Vnucl=0.0
    Vsurf=0.0
    Vspin=0.0
    Vcoul=0.0

!   LOCAL PART OF NL: N+A BOUND STATE
    Vv=-BoundParameters(5,1)
    Rv=BoundParameters(5,2)*A**(1.0/3.0)
    av=BoundParameters(5,3)

    Vd=-BoundParameters(5,4)
    Rvd=BoundParameters(5,5)*A**(1.0/3.0)
    avd=BoundParameters(5,6)

    Vso=-BoundParameters(5,7)
    Rso=BoundParameters(5,8)*A**(1.0/3.0)
    aso=BoundParameters(5,9)

    Rcoul=BoundParameters(5,10)*A**(1.0/3.0)

!    write(*,*) 'Local Part of NL'
!    write(*,*) 'Vv, Rv, av = ', Vv, Rv, av
!    write(*,*) 'Vd, Rvd, avd = ', Vd, Rvd, avd
!    write(*,*) 'Vso, Rso, aso = ', Vso, Rso, aso
!    write(*,*) ''

    if (abs(Vv)>0.0) then
      Vnucl=Vv/(1+exp((r-Rv)/av))
    end if
    if (abs(Vd)>0.0) then
      Vsurf=(4*Vd*exp((r-Rvd)/avd))/(1+exp((r-Rvd)/avd))**(2.0) 
    end if
    if (abs(Vso)>0.0) then
      Vspin=(J*(J+1.0)-L*(L+1.0)-SpinFragment*(SpinFragment+1.0))*(Vso/(aso*r)) &
            *((hbarc/139.6)**2.0)*(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
    end if
    if(r<Rcoul) then
      Vcoul=((Z1*Z2*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
    end if
    if(r>=Rcoul)then
      Vcoul=(Z1*Z2*1.43997)/r
    end if

!    write(*,*) 'R = ', R
!    write(*,*) 'J,L,Ip = ', J,L,SpinFragment
!    write(*,*) 'Vcoul = ', Vcoul
!    write(*,*) 'Vnucl = ', Vnucl
!    write(*,*) 'Vsurf = ', Vsurf
!    write(*,*) 'Vspin = ', Vspin
!    write(*,*) ''

    UlocNL=Vcoul+Vnucl+Vsurf+Vspin

!   WhatPot = 1 = User-defined
    if (WhatPot==1) then

      Vnucl=0.0
      Vsurf=0.0
      Vspin=0.0
      Vcoul=0.0

      p=(r+rprime)/2.0

!     NONLOCAL PART OF POTENTIAL: N+A
      Vv=-BoundParameters(6,1)
      Rv=BoundParameters(6,2)*A**(1.0/3.0)
      av=BoundParameters(6,3)

      Vd=-BoundParameters(6,4)
      Rvd=BoundParameters(6,5)*A**(1.0/3.0)
      avd=BoundParameters(6,6)

      Vso=-BoundParameters(6,7)
      Rso=BoundParameters(6,8)*A**(1.0/3.0)
      aso=BoundParameters(6,9)

!      write(*,*) 'Nonlocal Part'
!      write(*,*) 'Vv, Rv, av = ', Vv, Rv, av
!      write(*,*) 'Vd, Rvd, avd = ', Vd, Rvd, avd
!      write(*,*) 'Vso, Rso, aso = ', Vso, Rso, aso

!     NONLOCAL POTETENTIAL
      if (abs(Vv)>0.0) then
        Vnucl=Vv/(1+exp((p-Rv)/av))
      end if
      if (abs(Vd)>0.0) then
        Vsurf=(4*Vd*exp((p-Rvd)/avd))/(1+exp((p-Rvd)/avd))**(2.0) 
      end if
      if (abs(Vso)>0.0) then
        Vspin=(J*(J+1.0)-L*(L+1.0)-SpinFragment*(SpinFragment+1.0))*(Vso/(aso*p)) &
              *((hbarc/139.6)**2.0)*(exp((p-Rso)/aso))/((1+exp((p-Rso)/aso))**(2.0)) 
      end if

      Unl=Vnucl+Vsurf+Vspin

!   END OF USER-DEFINED NONLOCAL POTENTIAL
    end if

! END OF N+A BOUND STATE
  end if

  return
  end subroutine
