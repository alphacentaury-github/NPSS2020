  subroutine nonlocal_scattering_potential(ScatParameters,Unl1,Unl2,L,jp,r,rprime,UlocNL1,UlocNL2)

  implicit none
!!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!!   
  integer::n,nmax,L
  integer::WhatSystem,WhatPot,WhatPreDefPot,PotType
  real(8)::scattering_parameters(1:4,1:500),Vnucl,Wnucl,Vsurf,Wsurf,Vspin,Wspin,Vcoul,Zp,Zt,hbarc,A,r,rprime,P
  real(8)::Vv,Rv,av,Wv,Rwv,awv,Vd,Rvd,avd,Wd,Rwd,awd,Vso,Rso,aso,Wso,Rwso,awso,Rcoul,jp,Ip
  real(8)::ScatParameters(1:15,1:100)
  complex*16::i,Unl,Ulocnl1,UlocNL2,Unl1,Unl2


  i=cmplx(0.0,1.0)
  hbarc=197.3269718

  Zp=ScatParameters(1,3)
  Zt=ScatParameters(1,4)
  Ip=ScatParameters(1,5)
  A=ScatParameters(1,2)
  WhatSystem=nint(ScatParameters(2,1))

  UlocNL1=0
  UlocNL2=0
  Unl1=0
  Unl2=0

! WhatSystem = 1 = Deuteron scattering state
! WhatSystem = 2 = Nucleon scattering state
  if (WhatSystem == 1) then
    PotType=nint(ScatParameters(3,2))
    WhatPot=nint(ScatParameters(3,3))
    WhatPreDefPot=nint(ScatParameters(3,4))
  end if
  if (WhatSystem == 2) then
    WhatPot=nint(ScatParameters(3,3))
    WhatPreDefPot=nint(ScatParameters(3,4))
  end if

!  write(*,*) 'WhatSystem = ',WhatSystem
!  write(*,*) 'PotType = ', PotType
!  write(*,*) 'WhatPot = ', WhatPot
!  write(*,*) 'WhatPreDefPot = ', WhatPreDefPot

! WhatSystem = 1 = Deuteron Scattering State
! WhatSystem = 2 = Nucleon Scattering State
  if (WhatSystem == 1) then

!   PotType = 1 = Deuteron Optical Potential
!   Pottype = 2 = Nucleon Potentials
    if (PotType == 1) then

!     WhatPot = 1 = User defined deuteron nonlocal optical potential
!     WhatPot = 2 = Pre-defined deuteron nonlocal optical potential
!     WhatPot = 3 = Read in deuteron nonlocal optical potential
      if (WhatPot == 1) then

!       DEFINED THE LOCAL PART OF THE NONLOCAL POTENTIAL
        Vv=-ScatParameters(5,1)
        Rv=ScatParameters(5,2)*A**(1.0/3.0)
        av=ScatParameters(5,3)
        Wv=-ScatParameters(5,4)
        Rwv=ScatParameters(5,5)*A**(1.0/3.0)
        awv=ScatParameters(5,6)
        Vd=-ScatParameters(5,7)
        Rvd=ScatParameters(5,8)*A**(1.0/3.0)
        avd=ScatParameters(5,9)
        Wd=-ScatParameters(5,10)
        Rwd=ScatParameters(5,11)*A**(1.0/3.0)
        awd=ScatParameters(5,12)
        Vso=-ScatParameters(5,13)
        Rso=ScatParameters(5,14)*A**(1.0/3.0)
        aso=ScatParameters(5,15)
        Wso=-ScatParameters(5,16)
        Rwso=ScatParameters(5,17)*A**(1.0/3.0)
        awso=ScatParameters(5,18)
        Rcoul=ScatParameters(5,19)*A**(1.0/3.0)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

!        write(*,*) 'Vv,Rv,av=',Vv,Rv,av
!        write(*,*) 'Wv,Rwv,awv=',Wv,Rwv,awv
!        write(*,*) 'Vd,Rd,ad=',Vd,Rvd,avd
!        write(*,*) 'Wd,Rwd,awd=',Wd,Rwd,awd
!        write(*,*) 'Vso,Rso,aso=',Vso,Rso,aso
!        write(*,*) 'Wso,Rwso,awso=',Wso,Rwso,awso
!        write(*,*) 'Rcoul=',Rcoul

!       LOCAL PART OF NONLOCAL POTETENTIAL
        if (abs(Vv)>0.0) then
          Vnucl=Vv/(1+exp((r-Rv)/av))
        end if
        if (abs(Wv)>0.0) then
          Wnucl=Wv/(1+exp((r-Rwv)/awv))    
        end if
        if (abs(Vd)>0.0) then
          Vsurf=(4*Vd*exp((r-Rvd)/avd))/(1+exp((r-Rvd)/avd))**(2.0)  
        end if
        if (abs(Wd)>0.0) then
          Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
        end if
        if (abs(Vso)>0.0) then
          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                  /((1+exp((r-Rso)/aso))**(2.0))
        end if
        if (abs(Wso)>0.0) then
          Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*r))*((hbarc/139.6)**2.0)*(exp((r-Rwso)/awso)) &
                  /((1+exp((r-Rwso)/awso))**(2.0))
        end if 

!       LOCAL COULOMB POTENTIAL OF THE NONLOCAL POTENTIAL
        if (abs(Rcoul)>0.0) then
          if(r<Rcoul) then
            Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
          end if
          if(r>=Rcoul)then
            Vcoul=(Zp*Zt*1.43997)/r
          end if
        end if

!        UlocNL1=Vspin+i*Wspin
!        UlocNL1=0.0
        UlocNL2=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)+Vcoul

!       DEFINED THE NONLOCAL PART OF THE POTENTIAL
        p=(r+rprime)/2.0

        Vv=-ScatParameters(6,1)
        Rv=ScatParameters(6,2)*A**(1.0/3.0)
        av=ScatParameters(6,3)
        Wv=-ScatParameters(6,4)
        Rwv=ScatParameters(6,5)*A**(1.0/3.0)
        awv=ScatParameters(6,6)
        Vd=-ScatParameters(6,7)
        Rvd=ScatParameters(6,8)*A**(1.0/3.0)
        avd=ScatParameters(6,9)
        Wd=-ScatParameters(6,10)
        Rwd=ScatParameters(6,11)*A**(1.0/3.0)
        awd=ScatParameters(6,12)
        Vso=-ScatParameters(6,13)
        Rso=ScatParameters(6,14)*A**(1.0/3.0)
        aso=ScatParameters(6,15)
        Wso=-ScatParameters(6,16)
        Rwso=ScatParameters(6,17)*A**(1.0/3.0)
        awso=ScatParameters(6,18)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

!       NONLOCAL POTETENTIAL
        if (abs(Vv)>0.0) then
          Vnucl=Vv/(1+exp((p-Rv)/av)) 
        end if
        if (abs(Wv)>0.0) then
          Wnucl=Wv/(1+exp((p-Rwv)/awv))    
        end if
        if (abs(Vd)>0.0) then
          Vsurf=(4*Vd*exp((p-Rvd)/avd))/(1+exp((p-Rvd)/avd))**(2.0)  
        end if
        if (abs(Wd)>0.0) then
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  
        end if
        if (abs(Vso)>0.0) then
          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*p))*((hbarc/139.6)**2.0)*(exp((p-Rso)/aso)) &
                  /((1+exp((p-Rso)/aso))**(2.0))
        end if
        if (abs(Wso)>0.0) then
          Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*p))*((hbarc/139.6)**2.0)*(exp((p-Rwso)/awso)) &
                  /((1+exp((p-Rwso)/awso))**(2.0))
        end if 

        Unl2=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)

!     END USER DEFINED DEUTERON NONLOCAL OPTICAL POTENTIAL
      end if

!   END DEUTERON NONLOCAL OPTICAL POTENTIAL
    end if

!   PotType = 1 = Deuteron nonlocal optical potential
!   Pottype = 2 = Adiabatic potential for the deuteron
    if (PotType == 2) then

!     WhatPot = 1 = User defined deuteron nonlocal optical potential
!     WhatPot = 2 = Pre-defined deuteron nonlocal optical potential
!     WhatPot = 3 = Read in deuteron nonlocal optical potential
      if (WhatPot == 1) then
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       NEUTRON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Vv=-ScatParameters(9,1)
        Rv=ScatParameters(9,2)*A**(1.0/3.0)
        av=ScatParameters(9,3)
        Wv=-ScatParameters(9,4)
        Rwv=ScatParameters(9,5)*A**(1.0/3.0)
        awv=ScatParameters(9,6)
        Vd=-ScatParameters(9,7)
        Rvd=ScatParameters(9,8)*A**(1.0/3.0)
        avd=ScatParameters(9,9)
        Wd=-ScatParameters(9,10)
        Rwd=ScatParameters(9,11)*A**(1.0/3.0)
        awd=ScatParameters(9,12)
        Vso=-ScatParameters(9,13)
        Rso=ScatParameters(9,14)*A**(1.0/3.0)
        aso=ScatParameters(9,15)
        Wso=-ScatParameters(9,16)
        Rwso=ScatParameters(9,17)*A**(1.0/3.0)
        awso=ScatParameters(9,18)
        Rcoul=ScatParameters(9,19)*A**(1.0/3.0)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

!       LOCAL PART OF NONLOCAL POTETENTIAL
        if (abs(Vv)>0.0) then
          Vnucl=Vv/(1+exp((r-Rv)/av))
        end if
        if (abs(Wv)>0.0) then
          Wnucl=Wv/(1+exp((r-Rwv)/awv))    
        end if
        if (abs(Vd)>0.0) then
          Vsurf=(4*Vd*exp((r-Rvd)/avd))/(1+exp((r-Rvd)/avd))**(2.0)  
        end if
        if (abs(Wd)>0.0) then
          Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
        end if
        if (abs(Vso)>0.0) then
          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                  /((1+exp((r-Rso)/aso))**(2.0))
        end if
        if (abs(Wso)>0.0) then
          Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*r))*((hbarc/139.6)**2.0)*(exp((r-Rwso)/awso)) &
                  /((1+exp((r-Rwso)/awso))**(2.0))
        end if 

        UlocNL1=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PROTON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Vv=-ScatParameters(10,1)
        Rv=ScatParameters(10,2)*A**(1.0/3.0)
        av=ScatParameters(10,3)
        Wv=-ScatParameters(10,4)
        Rwv=ScatParameters(10,5)*A**(1.0/3.0)
        awv=ScatParameters(10,6)
        Vd=-ScatParameters(10,7)
        Rvd=ScatParameters(10,8)*A**(1.0/3.0)
        avd=ScatParameters(10,9)
        Wd=-ScatParameters(10,10)
        Rwd=ScatParameters(10,11)*A**(1.0/3.0)
        awd=ScatParameters(10,12)
        Vso=-ScatParameters(10,13)
        Rso=ScatParameters(10,14)*A**(1.0/3.0)
        aso=ScatParameters(10,15)
        Wso=-ScatParameters(10,16)
        Rwso=ScatParameters(10,17)*A**(1.0/3.0)
        awso=ScatParameters(10,18)
        Rcoul=ScatParameters(10,19)*A**(1.0/3.0)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

!       LOCAL PART OF NONLOCAL POTETENTIAL
        if (abs(Vv)>0.0) then
          Vnucl=Vv/(1+exp((r-Rv)/av))
        end if
        if (abs(Wv)>0.0) then
          Wnucl=Wv/(1+exp((r-Rwv)/awv))    
        end if
        if (abs(Vd)>0.0) then
          Vsurf=(4*Vd*exp((r-Rvd)/avd))/(1+exp((r-Rvd)/avd))**(2.0)  
        end if
        if (abs(Wd)>0.0) then
          Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
        end if
        if (abs(Vso)>0.0) then
          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                  /((1+exp((r-Rso)/aso))**(2.0))
        end if
        if (abs(Wso)>0.0) then
          Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*r))*((hbarc/139.6)**2.0)*(exp((r-Rwso)/awso)) &
                  /((1+exp((r-Rwso)/awso))**(2.0))
        end if 

!       LOCAL COULOMB POTENTIAL OF THE NONLOCAL POTENTIAL
        if (abs(Rcoul)>0.0) then
          if(r<Rcoul) then
            Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
          end if
          if(r>=Rcoul)then
            Vcoul=(Zp*Zt*1.43997)/r
          end if
        end if

        UlocNL2=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)+Vcoul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       NEUTRON POTENTIAL: NONLOCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       DEFINED THE NONLOCAL PART OF THE POTENTIAL
        p=(r+rprime)/2.0

        Vv=-ScatParameters(11,1)
        Rv=ScatParameters(11,2)*A**(1.0/3.0)
        av=ScatParameters(11,3)
        Wv=-ScatParameters(11,4)
        Rwv=ScatParameters(11,5)*A**(1.0/3.0)
        awv=ScatParameters(11,6)
        Vd=-ScatParameters(11,7)
        Rvd=ScatParameters(11,8)*A**(1.0/3.0)
        avd=ScatParameters(11,9)
        Wd=-ScatParameters(11,10)
        Rwd=ScatParameters(11,11)*A**(1.0/3.0)
        awd=ScatParameters(11,12)
        Vso=-ScatParameters(11,13)
        Rso=ScatParameters(11,14)*A**(1.0/3.0)
        aso=ScatParameters(11,15)
        Wso=-ScatParameters(11,16)
        Rwso=ScatParameters(11,17)*A**(1.0/3.0)
        awso=ScatParameters(11,18)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

!       NONLOCAL POTETENTIAL
        if (abs(Vv)>0.0) then
          Vnucl=Vv/(1+exp((p-Rv)/av)) 
        end if
        if (abs(Wv)>0.0) then
          Wnucl=Wv/(1+exp((p-Rwv)/awv))    
        end if
        if (abs(Vd)>0.0) then
          Vsurf=(4*Vd*exp((p-Rvd)/avd))/(1+exp((p-Rvd)/avd))**(2.0)  
        end if
        if (abs(Wd)>0.0) then
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  
        end if
        if (abs(Vso)>0.0) then
          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*p))*((hbarc/139.6)**2.0)*(exp((p-Rso)/aso)) &
                  /((1+exp((p-Rso)/aso))**(2.0))
        end if
        if (abs(Wso)>0.0) then
          Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*p))*((hbarc/139.6)**2.0)*(exp((p-Rwso)/awso)) &
                  /((1+exp((p-Rwso)/awso))**(2.0))
        end if 

        Unl1=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PROTON POTENTIAL: NONLOCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       DEFINED THE NONLOCAL PART OF THE POTENTIAL
        p=(r+rprime)/2.0

        Vv=-ScatParameters(12,1)
        Rv=ScatParameters(12,2)*A**(1.0/3.0)
        av=ScatParameters(12,3)
        Wv=-ScatParameters(12,4)
        Rwv=ScatParameters(12,5)*A**(1.0/3.0)
        awv=ScatParameters(12,6)
        Vd=-ScatParameters(12,7)
        Rvd=ScatParameters(12,8)*A**(1.0/3.0)
        avd=ScatParameters(12,9)
        Wd=-ScatParameters(12,10)
        Rwd=ScatParameters(12,11)*A**(1.0/3.0)
        awd=ScatParameters(12,12)
        Vso=-ScatParameters(12,13)
        Rso=ScatParameters(12,14)*A**(1.0/3.0)
        aso=ScatParameters(12,15)
        Wso=-ScatParameters(12,16)
        Rwso=ScatParameters(12,17)*A**(1.0/3.0)
        awso=ScatParameters(12,18)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

!       NONLOCAL POTETENTIAL
        if (abs(Vv)>0.0) then
          Vnucl=Vv/(1+exp((p-Rv)/av)) 
        end if
        if (abs(Wv)>0.0) then
          Wnucl=Wv/(1+exp((p-Rwv)/awv))    
        end if
        if (abs(Vd)>0.0) then
          Vsurf=(4*Vd*exp((p-Rvd)/avd))/(1+exp((p-Rvd)/avd))**(2.0)  
        end if
        if (abs(Wd)>0.0) then
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  
        end if
        if (abs(Vso)>0.0) then
          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*p))*((hbarc/139.6)**2.0)*(exp((p-Rso)/aso)) &
                  /((1+exp((p-Rso)/aso))**(2.0))
        end if
        if (abs(Wso)>0.0) then
          Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*p))*((hbarc/139.6)**2.0)*(exp((p-Rwso)/awso)) &
                  /((1+exp((p-Rwso)/awso))**(2.0))
        end if 

        Unl2=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)

!     END USER DEFINED NUCLEON OPTICAL POTENTIALS FOR THE DEUTERON
      end if

!     WhatPot = 1 = User defined nucleon nonlocal optical potential for the deuteron
!     WhatPot = 2 = Pre-defined nucleon nonlocal optical potential for the deuteron
!     WhatPot = 3 = Read in nucleon nonlocal optical potential for the deuteron
      if (WhatPot == 2) then

!       WhatPreDefPot = 1 = Perey-Buck
!       WhatPreDefPot = 2 = TPM
        if (WhatPreDefPot == 1) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         NEUTRON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          Vso=-7.180
          Rso=1.22*A**(1.0/3.0)
          aso=0.65

          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                    /((1+exp((r-Rso)/aso))**(2.0))

          UlocNL1=Vspin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         PROTON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          Vso=-7.180
          Rso=1.22*A**(1.0/3.0)
          aso=0.65

          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                    /((1+exp((r-Rso)/aso))**(2.0))

          Rcoul=1.22*A**(1.0/3.0)
          if(r<Rcoul) then
            Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
          end if
          if(r>=Rcoul)then
            Vcoul=(Zp*Zt*1.43997)/r
          end if

          UlocNL2=Vspin+Vcoul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         NEUTRON POTENTIAL: NONLOCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         DEFINED THE NONLOCAL PART OF THE POTENTIAL
          p=(r+rprime)/2.0

          Vv=-71.000
          Rv=1.22*A**(1.0/3.0)
          av=0.65
          Wd=-15.000
          Rwd=1.22*A**(1.0/3.0)
          awd=0.47

          Vnucl=Vv/(1+exp((p-Rv)/av)) 
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  

          Unl1=Vnucl+i*Wsurf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         PROTON POTENTIAL: NONLOCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          Vv=-71.000
          Rv=1.22*A**(1.0/3.0)
          av=0.65
          Wd=-15.000
          Rwd=1.22*A**(1.0/3.0)
          awd=0.47

          Vnucl=Vv/(1+exp((p-Rv)/av)) 
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  

          Unl2=Vnucl+i*Wsurf

!       END OF PEREY-BUCK PRE-DEFINED NONLOCAL POTENTIAL
        end if

!       WhatPreDefPot = 1 = Perey-Buck
!       WhatPreDefPot = 2 = TPM
        if (WhatPreDefPot == 2) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         NEUTRON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          Vso=-9.00
          Rso=1.10*A**(1.0/3.0)
          aso=0.59

          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                    /((1+exp((r-Rso)/aso))**(2.0))

          UlocNL1=Vspin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         PROTON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          Vso=-8.130
          Rso=1.02*A**(1.0/3.0)
          aso=0.59

          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                    /((1+exp((r-Rso)/aso))**(2.0))

          Rcoul=1.34*A**(1.0/3.0)
          if(r<Rcoul) then
            Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
          end if
          if(r>=Rcoul)then
            Vcoul=(Zp*Zt*1.43997)/r
          end if

          UlocNL2=Vspin+Vcoul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         NEUTRON POTENTIAL: NONLOCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         DEFINE THE NONLOCAL PART OF THE POTENTIAL
          p=(r+rprime)/2.0

          Vv=-70.00
          Rv=1.25*A**(1.0/3.0)
          av=0.61
          Wv=-1.39
          Rwv=1.17*A**(1.0/3.0)
          awv=0.55
          Wd=-21.11
          Rwd=1.15*A**(1.0/3.0)
          awd=0.46

          Vnucl=Vv/(1+exp((p-Rv)/av)) 
          Wnucl=Wv/(1+exp((p-Rwv)/awv))  
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  

          Unl1=Vnucl+i*(Wnucl+Wsurf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         PROTON POTENTIAL: NONLOCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          Vv=-70.95
          Rv=1.29*A**(1.0/3.0)
          av=0.58
          Wv=-9.03
          Rwv=1.24*A**(1.0/3.0)
          awv=0.50
          Wd=-15.74
          Rwd=1.20*A**(1.0/3.0)
          awd=0.45

          Vnucl=Vv/(1+exp((p-Rv)/av)) 
          Wnucl=Wv/(1+exp((p-Rwv)/awv))  
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  

          Unl2=Vnucl+i*(Wnucl+Wsurf)

!       END OF TPM PRE-DEFINED NONLOCAL NUCLEON POTENTIAL FOR THE DEUTERON
        end if

      end if

!   END NUCLEON OPTICAL POTENTIALS FOR THE DEUTERON
    end if

! END DEUTERON SCATTERING STATE
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NUCLEON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! WhatSystem = 1 = Deuteron Scattering State
! WhatSystem = 2 = Nucleon Scattering State
  if (WhatSystem == 2 ) then

!   WhatPot = 1 = User-defined nonlocal potential
!   WhatPot = 2 = Pre-defined nonlocal potential
!   WhatPot = 3 = Read in nonlocal potential
    if (WhatPot==1) then

!     DEFINED THE LOCAL PART OF THE NONLOCAL POTENTIAL
      Vv=-ScatParameters(5,1)
      Rv=ScatParameters(5,2)*A**(1.0/3.0)
      av=ScatParameters(5,3)
      Wv=-ScatParameters(5,4)
      Rwv=ScatParameters(5,5)*A**(1.0/3.0)
      awv=ScatParameters(5,6)
      Vd=-ScatParameters(5,7)
      Rvd=ScatParameters(5,8)*A**(1.0/3.0)
      avd=ScatParameters(5,9)
      Wd=-ScatParameters(5,10)
      Rwd=ScatParameters(5,11)*A**(1.0/3.0)
      awd=ScatParameters(5,12)
      Vso=-ScatParameters(5,13)
      Rso=ScatParameters(5,14)*A**(1.0/3.0)
      aso=ScatParameters(5,15)
      Wso=-ScatParameters(5,16)
      Rwso=ScatParameters(5,17)*A**(1.0/3.0)
      awso=ScatParameters(5,18)
      Rcoul=ScatParameters(5,19)*A**(1.0/3.0)

      Vnucl=0.0
      Wnucl=0.0
      Vsurf=0.0
      Wsurf=0.0
      Vspin=0.0
      Wspin=0.0

!      write(*,*) 'Vv,Rv,av=',Vv,Rv,av
!      write(*,*) 'Wv,Rwv,awv=',Wv,Rwv,awv
!      write(*,*) 'Vd,Rd,ad=',Vd,Rvd,avd
!      write(*,*) 'Wd,Rwd,awd=',Wd,Rwd,awd
!      write(*,*) 'Vso,Rso,aso=',Vso,Rso,aso
!      write(*,*) 'Wso,Rwso,awso=',Wso,Rwso,awso
!      write(*,*) 'Rcoul=',Rcoul

!     LOCAL PART OF NONLOCAL POTETENTIAL
      if (abs(Vv)>0.0) then
        Vnucl=Vv/(1+exp((r-Rv)/av))
      end if
      if (abs(Wv)>0.0) then
        Wnucl=Wv/(1+exp((r-Rwv)/awv))    
      end if
      if (abs(Vd)>0.0) then
        Vsurf=(4*Vd*exp((r-Rvd)/avd))/(1+exp((r-Rvd)/avd))**(2.0)  
      end if
      if (abs(Wd)>0.0) then
        Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
      end if
      if (abs(Vso)>0.0) then
        Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                /((1+exp((r-Rso)/aso))**(2.0))
      end if
      if (abs(Wso)>0.0) then
        Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*r))*((hbarc/139.6)**2.0)*(exp((r-Rwso)/awso)) &
                /((1+exp((r-Rwso)/awso))**(2.0))
      end if 

!     LOCAL COULOMB POTENTIAL OF THE NONLOCAL POTENTIAL
      if (abs(Rcoul)>0.0) then
        if(r<Rcoul) then
          Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
        end if
        if(r>=Rcoul)then
          Vcoul=(Zp*Zt*1.43997)/r
        end if
      end if

!      UlocNL1=Vspin+i*Wspin
      UlocNL2=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)+Vcoul

!     DEFINED THE NONLOCAL PART OF THE POTENTIAL
      p=(r+rprime)/2.0

      Vv=-ScatParameters(6,1)
      Rv=ScatParameters(6,2)*A**(1.0/3.0)
      av=ScatParameters(6,3)
      Wv=-ScatParameters(6,4)
      Rwv=ScatParameters(6,5)*A**(1.0/3.0)
      awv=ScatParameters(6,6)
      Vd=-ScatParameters(6,7)
      Rvd=ScatParameters(6,8)*A**(1.0/3.0)
      avd=ScatParameters(6,9)
      Wd=-ScatParameters(6,10)
      Rwd=ScatParameters(6,11)*A**(1.0/3.0)
      awd=ScatParameters(6,12)
      Vso=-ScatParameters(6,13)
      Rso=ScatParameters(6,14)*A**(1.0/3.0)
      aso=ScatParameters(6,15)
      Wso=-ScatParameters(6,16)
      Rwso=ScatParameters(6,17)*A**(1.0/3.0)
      awso=ScatParameters(6,18)

!      write(*,*) 'Vv,Rv,av=',Vv,Rv,av
!      write(*,*) 'Wv,Rwv,awv=',Wv,Rwv,awv
!      write(*,*) 'Vd,Rd,ad=',Vd,Rvd,avd
!      write(*,*) 'Wd,Rwd,awd=',Wd,Rwd,awd
!      write(*,*) 'Vso,Rso,aso=',Vso,Rso,aso
!      write(*,*) 'Wso,Rwso,awso=',Wso,Rwso,awso

      Vnucl=0.0
      Wnucl=0.0
      Vsurf=0.0
      Wsurf=0.0
      Vspin=0.0
      Wspin=0.0

!     NONLOCAL POTETENTIAL
      if (abs(Vv)>0.0) then
        Vnucl=Vv/(1+exp((p-Rv)/av))
      end if
      if (abs(Wv)>0.0) then
        Wnucl=Wv/(1+exp((p-Rwv)/awv))    
      end if
      if (abs(Vd)>0.0) then
        Vsurf=(4*Vd*exp((p-Rvd)/avd))/(1+exp((p-Rvd)/avd))**(2.0)  
      end if
      if (abs(Wd)>0.0) then
        Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  
      end if
      if (abs(Vso)>0.0) then
        Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*p))*((hbarc/139.6)**2.0)*(exp((p-Rso)/aso)) &
                /((1+exp((p-Rso)/aso))**(2.0))
      end if
      if (abs(Wso)>0.0) then
        Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*p))*((hbarc/139.6)**2.0)*(exp((p-Rwso)/awso)) &
                /((1+exp((p-Rwso)/awso))**(2.0))
      end if 

      Unl2=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)
!      Unl=Vnucl+i*Wsurf

!   END USER-DEFINED NONLOCAL POTENTIAL
    end if

!   WhatPot = 1 = User-defined nonlocal potential
!   WhatPot = 2 = Pre-defined nonlocal potential
    if (WhatPot==2) then

!     WhatPreDefPot = 1 = Perey-Buck
!     WhatPreDefPot = 2 = TPM
      if (WhatPreDefPot == 1) then
        p=(r+rprime)/2.0

        Vv=-71.0
        Rv=1.220*A**(1.0/3.0)
        av=0.65
        Wd=-15.0
        Rwd=1.22*A**(1.0/3.0)
        awd=0.47
        Vso=-7.180
        Rso=1.220*A**(1.0/3.0)
        aso=0.65
        Rcoul=1.22*A**(1.0/3.0)

!       NONLOCAL POTETENTIAL
        Vnucl=Vv/(1+exp((p-Rv)/av))  
        Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  
        Unl1=Vnucl+i*Wsurf
!        Unl2=Vnucl+i*Wsurf

!       LOCAL COULOMB POTENTIAL OF THE NONLOCAL POTENTIAL
        if(r<Rcoul) then
          Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
        end if
        if(r>=Rcoul)then
          Vcoul=(Zp*Zt*1.43997)/r
        end if
!       LOCAL SPIN-ORBIT OF THE NONLOCAL POTENTIAL
        Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                /((1+exp((r-Rso)/aso))**(2.0))
        UlocNL1=Vcoul+Vspin
!        UlocNL2=Vcoul+Vspin

!        write(*,*) 'Unl=',Unl
!        write(*,*) 'UlocNL1=',UlocNL1
!        write(*,*) 'UlocNL2=',UlocNL2
!        write(*,*) ''

!     END PEREY-BUCK
      end if

!     WhatPreDefPot = 1 = Perey-Buck
!     WhatPreDefPot = 2 = TPM
      if (WhatPreDefPot == 2) then

!       NEUTERON TPM NONLOCAL POTENTIAL
        if(nint(Zp)==0) then

          Vso=-9.00
          Rso=1.10*A**(1.0/3.0)
          aso=0.59

          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                    /((1+exp((r-Rso)/aso))**(2.0))

          UlocNL1=Vspin

!         DEFINE THE NONLOCAL PART OF THE POTENTIAL
          p=(r+rprime)/2.0

          Vv=-70.00
          Rv=1.25*A**(1.0/3.0)
          av=0.61
          Wv=-1.39
          Rwv=1.17*A**(1.0/3.0)
          awv=0.55
          Wd=-21.11
          Rwd=1.15*A**(1.0/3.0)
          awd=0.46

          Vnucl=Vv/(1+exp((p-Rv)/av)) 
          Wnucl=Wv/(1+exp((p-Rwv)/awv))  
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  

          Unl1=Vnucl+i*(Wnucl+Wsurf)          

        end if

!       PROTON TPM NONLOCAL POTENTIAL
        if(nint(Zp)==1) then

          Vso=-8.130
          Rso=1.02*A**(1.0/3.0)
          aso=0.59

          Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                    /((1+exp((r-Rso)/aso))**(2.0))

          Rcoul=1.34*A**(1.0/3.0)
          if(r<Rcoul) then
            Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
          end if
          if(r>=Rcoul)then
            Vcoul=(Zp*Zt*1.43997)/r
          end if

          UlocNL2=Vspin+Vcoul

!         DEFINE THE NONLOCAL PART OF THE POTENTIAL
          p=(r+rprime)/2.0

          Vv=-70.95
          Rv=1.29*A**(1.0/3.0)
          av=0.58
          Wv=-9.03
          Rwv=1.24*A**(1.0/3.0)
          awv=0.50
          Wd=-15.74
          Rwd=1.20*A**(1.0/3.0)
          awd=0.45

          Vnucl=Vv/(1+exp((p-Rv)/av)) 
          Wnucl=Wv/(1+exp((p-Rwv)/awv))  
          Wsurf=(4*Wd*exp((p-Rwd)/awd))/(1+exp((p-Rwd)/awd))**(2.0)  

          Unl2=Vnucl+i*(Wnucl+Wsurf)

        end if

!     END TPM
      end if

!   END PRE-DEFINED NONLOCAL POTENTIAL
    end if

!   WhatPot = 1 = User-defined nonlocal potential
!   WhatPot = 2 = Pre-defined nonlocal potential
!   WhatPot = 3 = Read in nonlocal potential
    if (WhatPot == 3) then

!     DEFINED THE LOCAL PART OF THE NONLOCAL POTENTIAL
      Vv=-ScatParameters(5,1)
      Rv=ScatParameters(5,2)*A**(1.0/3.0)
      av=ScatParameters(5,3)
      Wv=-ScatParameters(5,4)
      Rwv=ScatParameters(5,5)*A**(1.0/3.0)
      awv=ScatParameters(5,6)
      Vd=-ScatParameters(5,7)
      Rvd=ScatParameters(5,8)*A**(1.0/3.0)
      avd=ScatParameters(5,9)
      Wd=-ScatParameters(5,10)
      Rwd=ScatParameters(5,11)*A**(1.0/3.0)
      awd=ScatParameters(5,12)
      Vso=-ScatParameters(5,13)
      Rso=ScatParameters(5,14)*A**(1.0/3.0)
      aso=ScatParameters(5,15)
      Wso=-ScatParameters(5,16)
      Rwso=ScatParameters(5,17)*A**(1.0/3.0)
      awso=ScatParameters(5,18)
      Rcoul=ScatParameters(5,19)*A**(1.0/3.0)

      Vnucl=0.0
      Wnucl=0.0
      Vsurf=0.0
      Wsurf=0.0
      Vspin=0.0
      Wspin=0.0

!      write(*,*) 'Vv,Rv,av=',Vv,Rv,av
!      write(*,*) 'Wv,Rwv,awv=',Wv,Rwv,awv
!      write(*,*) 'Vd,Rd,ad=',Vd,Rvd,avd
!      write(*,*) 'Wd,Rwd,awd=',Wd,Rwd,awd
!      write(*,*) 'Vso,Rso,aso=',Vso,Rso,aso
!      write(*,*) 'Wso,Rwso,awso=',Wso,Rwso,awso
!      write(*,*) 'Rcoul=',Rcoul

!     LOCAL PART OF NONLOCAL POTETENTIAL
      if (abs(Vv)>0.0) then
        Vnucl=Vv/(1+exp((r-Rv)/av))
      end if
      if (abs(Wv)>0.0) then
        Wnucl=Wv/(1+exp((r-Rwv)/awv))    
      end if
      if (abs(Vd)>0.0) then
        Vsurf=(4*Vd*exp((r-Rvd)/avd))/(1+exp((r-Rvd)/avd))**(2.0)  
      end if
      if (abs(Wd)>0.0) then
        Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
      end if
      if (abs(Vso)>0.0) then
        Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0)*(exp((r-Rso)/aso)) &
                /((1+exp((r-Rso)/aso))**(2.0))
      end if
      if (abs(Wso)>0.0) then
        Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(awso*r))*((hbarc/139.6)**2.0)*(exp((r-Rwso)/awso)) &
                /((1+exp((r-Rwso)/awso))**(2.0))
      end if 

!     LOCAL COULOMB POTENTIAL OF THE NONLOCAL POTENTIAL
      if (abs(Rcoul)>0.0) then
        if(r<Rcoul) then
          Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
        end if
        if(r>=Rcoul)then
          Vcoul=(Zp*Zt*1.43997)/r
        end if
      end if

!      UlocNL1=Vspin+i*Wspin
      UlocNL2=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)+Vcoul      

    end if


! END NUCLEON SCATTERING STATE
  end if


  return
  end subroutine
