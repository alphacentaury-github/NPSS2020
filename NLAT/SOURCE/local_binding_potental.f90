! CREATED BY: LUKE TITUS
! LAST UPDATED: 2/2/2016
! THIS PROGRAM OUTPUTS THE LOCAL BINDING POTENTIAL CHOSEN BY THE USER

  subroutine local_binding_potential(BoundParameters,vpot,nmax,StepSize)
	
  IMPLICIT NONE
!!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!!
  integer::n,nmax,L,WhatSystem,WhatPot
  real(8)::vpot(nmax),StepSize
  real(8)::BoundParameters(1:15,1:100)
  real(8)::r,ChargeCore,ChargeFragment,MassCore,SpinFragment,J,hbarc
  real(8)::Vv,Rv,av,Vd,Rd,ad,Vso,Rso,aso,Rcoul,Vnucl,Vsurf,Vspin,Vcoul

  vpot(:)=0.0
  ChargeFragment=BoundParameters(1,3)
  ChargeCore=BoundParameters(1,4)
  MassCore=BoundParameters(1,2)
  SpinFragment=BoundParameters(1,5)
  L=nint(BoundParameters(1,9))
  J=BoundParameters(1,10)
  WhatSystem=nint(BoundParameters(2,1))
  WhatPot=nint(BoundParameters(2,2))
  hbarc=197.32705

!  write(*,*) 'ChargeFragment=',ChargeFragment
!  write(*,*) 'ChargeCore=',ChargeCore
!  write(*,*) 'MassCore=',MassCore
!  write(*,*) 'SpinFragment=',SpinFragment
!  write(*,*) 'L=',L
!  write(*,*) 'J=',J
!  write(*,*) 'WhatSystem=',WhatSystem
!  write(*,*) 'WhatPot=',WhatPot

! WhatSystem = 1 = Deuteron Bound State
! WhatSystem = 2 = Nucleon Bound State
  if (WhatSystem==1) then

!   WhatPot = 1 = Central Gaussian
    if (WhatPot==1) then
      do n=1,nmax
        vpot(n)=-71.85*exp(-(n*StepSize/1.494)**2.0)
      end do      
    end if

! END DEUTERON
  end if

! WhatSystem = 1 = Deuteron Bound State
! WhatSystem = 2 = Nucleon Bound State
  if (WhatSystem==2) then
   
!   WhatPot = 1 = WOODS-SAXON
    if (WhatPot==1) then

      Vv=-BoundParameters(4,1)
      Rv=BoundParameters(4,2)*(MassCore**(1.0/3.0))
      av=BoundParameters(4,3)

      Vd=-BoundParameters(4,4)
      Rd=BoundParameters(4,5)*(MassCore**(1.0/3.0))
      ad=BoundParameters(4,6)

      Vso=-BoundParameters(4,7)
      Rso=BoundParameters(4,8)*(MassCore**(1.0/3.0))
      aso=BoundParameters(4,9)

      Rcoul=BoundParameters(4,10)*(MassCore**(1.0/3.0))

      do n=1,nmax           
        r=n*StepSize
!       NUCLEAR POTENTIAL
        if (abs(Vv)>0.0) then
          Vnucl=Vv/(1+exp((r-Rv)/av))
        end if
        if (abs(Vd)>0.0) then
          Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0)   
        end if
        if (abs(Vso)>0.0) then
          Vspin=(J*(J+1.0)-L*(L+1.0)-SpinFragment*(SpinFragment+1.0)) &
                *(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0)) 
        end if
!       COULOMB POTENTIAL
        if(r<Rcoul) then
          Vcoul=((ChargeCore*ChargeFragment*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
        end if
        if(R >=Rcoul)then
          Vcoul=(ChargeCore*ChargeFragment*1.43997)/r
        end if
!       DEFINE THE LOCAL POTENTIAL
        vpot(n)=Vcoul+Vnucl+Vsurf+Vspin
      end do
    
!   END WOODS-SAXON
    end if

! END NUCLEON BOUND STATE
  end if

  return
  end subroutine
  
