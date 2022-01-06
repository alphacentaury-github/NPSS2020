!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CREATED BY: LUKE TITUS
! LAST MODIFIED: 2/4/2016
!
! This code outputs the local nuclear potentials Upot1 and Upot2, the spin-orbit
! potentials Spin1 and Spin2, and the Coulomb potential. The potentials labeled '1' are for the 
! nucleon scattering state. If a deuteron scattering state is considered, then the '1' labels
! the neutron, and the '2' labels the proton.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine local_scattering_potential(ScatParameters,Upot1,Upot2,nmax,L,jp,Coulomb,Spin1,Spin2,StepSize,E)
  implicit none

  integer::n,nmax,L
  integer::WhatSystem,WhatPot,PotType,WhatPreDefPot
  real(8)::scattering_parameters(1:4,1:500),Vnucl,Wnucl,Vsurf,Wsurf,Vspin,Wspin,Vcoul,Zp,Zt,StepSize,r,hbarc,A,Nt,E
  real(8)::Vv,Rv,av,Wv,Rwv,awv,Vd,Rd,ad,Wd,Rwd,awd,Vso,Rso,aso,Wso,Rwso,awso,Rcoul,jp,Ip,Coulomb(nmax)
  real(8)::ScatParameters(1:15,1:100)
  real(8)::SpinProjectile,MagicNumbers(6),mu,sum
  real(8)::v1n,v2n,v3n,v4n,w1n,w2n,d1n,d2n,d3n,vso1n,vso2n,wso1n,wso2n,Efn,v1p,v2p,v3p,v4p,w1p,w2p,d1p,d2p,d3p 
  real(8)::vso1p,vso2p,wso1p,wso2p,Efp,Vcp,Rc,Ec
  complex*16::i,vpot(nmax),Upot1(nmax),Upot2(nmax),Spin1(nmax),Spin2(nmax)


  i=cmplx(0.0,1.0)
  hbarc=197.3269718

  Zp=ScatParameters(1,3)
  Zt=ScatParameters(1,4)
  Ip=ScatParameters(1,5)
  A=ScatParameters(1,2)
  WhatSystem=nint(ScatParameters(2,1))
  Nt=A-Zt  


! FIX THIS!!! IT IS CONFUSING

! WhatSystem = 1 = Deuteron Scattering State
! WhatSystem = 2 = Nucleon Scattering State
  if (WhatSystem == 1) then

    PotType=nint(ScatParameters(2,2))

!   PotType = 1 = Local Deuteron Optical Potential
!   PotType = 2 = Local Nucleon Optical Potentials
    if (PotType == 1) then
      WhatPot=nint(ScatParameters(2,3))
      WhatPreDefPot=nint(ScatParameters(2,4))
    end if
    if (PotType == 2) then
      WhatPot=nint(ScatParameters(2,5))
      WhatPreDefPot=nint(ScatParameters(2,6))
    end if
  end if

! WhatSystem = 1 = Deuteron Scattering State
! WhatSystem = 2 = Nucleon Scattering State
  if (WhatSystem == 2) then
    WhatPot=nint(ScatParameters(2,2))
    WhatPreDefPot=nint(ScatParameters(2,3))
  end if

!  write(*,*) 'StepSize = ', StepSize
!  write(*,*) 'Zp,Zt=',Zp,Zt
!  write(*,*) 'Ip = ', Ip
!  write(*,*) 'A = ', A
!  write(*,*) 'WhatSystem = ', WhatSystem
!  write(*,*) 'PotType = ', PotType
!  write(*,*) 'WhatPot = ', WhatPot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEUTERON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! WhatSystem = 1 = Deuteron Scattering State
! WhatSystem = 2 = Nucleon Scattering State
  if (WhatSystem == 1) then
    Upot2(:)=0.0
    Spin2(:)=0.0

!   PotType = 1 = Local Deuteron Optical Potential
!   PotType = 2 = Local Nucleon Optical Potentials
    if (PotType == 1) then

!     WhatPot = 1 = User-defined local potential 
!     WhatPot = 2 = Pre-defined local potential
      if (WhatPot==1) then

        Vv=-ScatParameters(4,1)
        Rv=ScatParameters(4,2)*A**(1.0/3.0)
        av=ScatParameters(4,3)
        Wv=-ScatParameters(4,4)
        Rwv=ScatParameters(4,5)*A**(1.0/3.0)
        awv=ScatParameters(4,6)
        Vd=-ScatParameters(4,7)
        Rd=ScatParameters(4,8)*A**(1.0/3.0)
        ad=ScatParameters(4,9)
        Wd=-ScatParameters(4,10)
        Rwd=ScatParameters(4,11)*A**(1.0/3.0)
        awd=ScatParameters(4,12)
        Vso=-ScatParameters(4,13)
        Rso=ScatParameters(4,14)*A**(1.0/3.0)
        aso=ScatParameters(4,15)
        Wso=-ScatParameters(4,16)
        Rwso=ScatParameters(4,17)*A**(1.0/3.0)
        awso=ScatParameters(4,18)
        Rcoul=ScatParameters(4,19)*A**(1.0/3.0)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

        do n=1,nmax
          r=n*StepSize
!         NUCLEAR POTENTIAL
          if (abs(Vv)>0.0) then
            Vnucl=Vv/(1+exp((r-Rv)/av))
          end if
          if (abs(Wv)>0.0) then
            Wnucl=Wv/(1+exp((r-Rwv)/awv))
          end if
          if (abs(Vd)>0.0) then
            Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
          end if
          if (abs(Wd)>0.0) then 
            Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
          end if
          if (abs(Vso)>0.0) then
            Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                  *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
          end if
          if (abs(Wso)>0.0) then
            Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                  *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
          end if
!         COULOMB POTENTIAL
          if(r<Rcoul) then
            Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
          end if
          if(r>=Rcoul)then
            Vcoul=(Zp*Zt*1.43997)/r
          end if

!          Upot1(n)=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)
          Upot1(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
          Coulomb(n)=Vcoul
          Spin1(n)=Vspin+i*Wspin

        end do

!     USER DEFINED LOCAL DEUTERON OPTICAL POTENTIAL
      end if

!     WhatPot = 1 = User-defined local potential 
!     WhatPot = 2 = Pre-defined local potential
      if (WhatPot==2) then

!       WhatPreDefPot = 1 = Daehnick local deuteorn optical potential
        if (WhatPreDefPot == 1) then

          MagicNumbers(1)=8
          MagicNumbers(2)=20
          MagicNumbers(3)=28
          MagicNumbers(4)=50
          MagicNumbers(5)=82
          MagicNumbers(6)=126
          sum=0
          
          do n=1,6
            mu=((MagicNumbers(n)-Nt)/2)**2.0
            sum=sum+exp(-mu)
          end do

          Vv=-(88.5-0.26*E+0.88*Zt*A**(-1.0/3.0))
          rv=1.17
          Rv=rv*A**(1.0/3.0)
          av=0.709+0.0017*E
          Wv=-(12.2+0.026*E)*(1-exp(-(E/100)**2.0))
          rwv=1.325
          Rwv=rwv*A**(1.0/3.0)
          awv=0.53+0.07*A**(1.0/3.0)-0.04*sum
          Wd=-(12.2+0.026*E)*exp(-(E/100)**2)
          Rwd=Rwv
          awd=awv
          Vso=(7.33-0.029*E)
          rso=1.07
          Rso=rso*A**(1.0/3.0)
          aso=0.66
          Rcoul=1.30*A**(1.0/3.0)
          
!          write(*,*) 'E = ',E
!          write(*,*) 'Vv = ', Vv
!          write(*,*) 'rv = ', Rv
!          write(*,*) 'av = ', av
!          write(*,*) 'Wv = ', Wv
!          write(*,*) 'Rwv = ', Rwv
!          write(*,*) 'awv = ', awv
!          write(*,*) 'Wd = ', Wd
!          write(*,*) 'Rwd = ', Rwd
!          write(*,*) 'awd = ', awd
!          write(*,*) 'Vso = ', Vso
!          write(*,*) 'Rso = ', Rso
!          write(*,*) 'aso = ', aso
!          write(*,*) 'Rcoul = ', Rcoul

          Vnucl=0.0
          Wnucl=0.0
          Vsurf=0.0
          Wsurf=0.0
          Vspin=0.0
          Wspin=0.0

          do n=1,nmax
            r=n*StepSize
            Vnucl=Vv/(1+exp((r-Rv)/av))
            Wnucl=Wv/(1+exp((r-Rwv)/awv))
            Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
!           COULOMB POTENTIAL
            if(r<Rcoul) then
              Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
            end if
            if(r>=Rcoul)then
              Vcoul=(Zp*Zt*1.43997)/r
            end if

!            Upot1(n)=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)
            Upot1(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Coulomb(n)=Vcoul
            Spin1(n)=Vspin+i*Wspin

          end do

!       END DAEHNICK LOCAL DEUTERON OPTICAL POTENTIAL
        end if

!     END PRE-DEFINED LOCAL DEUTERON OPTICAL POTENTIAL
      end if

!   END LOCAL DEUTERON OPTICAL POTENTIAL
    end if

!   PotType = 1 = Local Deuteron Optical Potential
!   PotType = 2 = Local Nucleon Optical Potentials
    if (PotType == 2) then

!     WhatPot = 1 = User-defined local potential 
!     WhatPot = 2 = Pre-defined local potential
      if (WhatPot==1) then

!       NEUTRON POTENTIAL
        Vv=-ScatParameters(7,1)
        Rv=ScatParameters(7,2)*A**(1.0/3.0)
        av=ScatParameters(7,3)
        Wv=-ScatParameters(7,4)
        Rwv=ScatParameters(7,5)*A**(1.0/3.0)
        awv=ScatParameters(7,6)
        Vd=-ScatParameters(7,7)
        Rd=ScatParameters(7,8)*A**(1.0/3.0)
        ad=ScatParameters(7,9)
        Wd=-ScatParameters(7,10)
        Rwd=ScatParameters(7,11)*A**(1.0/3.0)
        awd=ScatParameters(7,12)
        Vso=-ScatParameters(7,13)
        Rso=ScatParameters(7,14)*A**(1.0/3.0)
        aso=ScatParameters(7,15)
        Wso=-ScatParameters(7,16)
        Rwso=ScatParameters(7,17)*A**(1.0/3.0)
        awso=ScatParameters(7,18)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

        do n=1,nmax
          r=n*StepSize
!         NUCLEAR POTENTIAL
          if (abs(Vv)>0.0) then
            Vnucl=Vv/(1+exp((r-Rv)/av))
          end if
          if (abs(Wv)>0.0) then
            Wnucl=Wv/(1+exp((r-Rwv)/awv))
          end if
          if (abs(Vd)>0.0) then
            Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
          end if
          if (abs(Wd)>0.0) then 
            Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
          end if
          if (abs(Vso)>0.0) then
            Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                  *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
          end if
          if (abs(Wso)>0.0) then
            Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                  *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
          end if

          Upot1(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
          Spin1(n)=Vspin+i*Wspin

        end do

!       PROTON POTENTIAL
        Vv=-ScatParameters(8,1)
        Rv=ScatParameters(8,2)*A**(1.0/3.0)
        av=ScatParameters(8,3)
        Wv=-ScatParameters(8,4)
        Rwv=ScatParameters(8,5)*A**(1.0/3.0)
        awv=ScatParameters(8,6)
        Vd=-ScatParameters(8,7)
        Rd=ScatParameters(8,8)*A**(1.0/3.0)
        ad=ScatParameters(8,9)
        Wd=-ScatParameters(8,10)
        Rwd=ScatParameters(8,11)*A**(1.0/3.0)
        awd=ScatParameters(8,12)
        Vso=-ScatParameters(8,13)
        Rso=ScatParameters(8,14)*A**(1.0/3.0)
        aso=ScatParameters(8,15)
        Wso=-ScatParameters(8,16)
        Rwso=ScatParameters(8,17)*A**(1.0/3.0)
        awso=ScatParameters(8,18)
        Rcoul=ScatParameters(8,19)*A**(1.0/3.0)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

        do n=1,nmax
          r=n*StepSize
!         NUCLEAR POTENTIAL
          if (abs(Vv)>0.0) then
            Vnucl=Vv/(1+exp((r-Rv)/av))
          end if
          if (abs(Wv)>0.0) then
            Wnucl=Wv/(1+exp((r-Rwv)/awv))
          end if
          if (abs(Vd)>0.0) then
            Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
          end if
          if (abs(Wd)>0.0) then 
            Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
          end if
          if (abs(Vso)>0.0) then
            Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                  *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
          end if
          if (abs(Wso)>0.0) then
            Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                  *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
          end if
!         COULOMB POTENTIAL
          if(r<Rcoul) then
            Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
          end if
          if(r>=Rcoul)then
            Vcoul=(Zp*Zt*1.43997)/r
          end if

          Upot2(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
          Coulomb(n)=Vcoul
          Spin2(n)=Vspin+i*Wspin

        end do

!     END USER DEFINED NUCLEON POTENTIALS FOR THE DEUTERON
      end if

!     WhatPot = 1 = User-defined local potential 
!     WhatPot = 2 = Pre-defined local potential
      if (WhatPot==2) then

!       EVALUATE THE NUCLEON POTENTIALS AT HALF THE DEUTERON ENERGY
        E=E/2

!       WhatPreDefPot = 1 = Koning-Delaroche
!       WhatPreDefPot = 2 = Chapel Hill
        if (WhatPreDefPot == 1) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         KONING-DELAROCHE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         NEUTRON POTENTIAL
          v1n=59.3-21*(Nt-Zt)/A-0.024*A
          v2n=0.007228-(1.48E-6)*A
          v3n=(1.994E-5)-(2E-8)*A
          v4n=7E-9
          w1n=12.195+0.0167*A
          w2n=73.55+0.0795*A
          d1n=16-16*(Nt-Zt)/A
          d2n=0.018+0.003802/(1+exp((A-156)/8))
          d3n=11.5
          vso1n=5.922+0.003*A
          vso2n=0.004
          wso1n=-3.1
          wso2n=160
          Efn=-11.2814+0.02646*A
		
!         VOLUME
          Vv=-(v1n*(1-v2n*(E-Efn)+v3n*(E-Efn)**2-v4n*(E-Efn)**3))
          Wv=-(w1n*((E-Efn)**2/((E-Efn)**2+w2n**2)))
          Rv=1.3039-0.4054*A**(-0.3333333)
          Rv=Rv*A**(1.0/3.0)
          Rwv=Rv
          av=0.6778-(1.487E-4)*A
          awv=av

!         SURFACE
          Vd=0.0
          Wd=-(d1n*((E-Efn)**2/((E-Efn)**2+d3n**2))*exp(-d2n*(E-Efn)))
          Rwd=1.3424-0.01585*A**(0.3333333)
          Rwd=Rwd*A**(1.0/3.0)
          awd=0.5446-(1.656E-4)*A

!         SPIN-ORBIT
          Vso=-(vso1n*exp(-vso2n*(E-Efn)))
          Wso=-(wso1n*((E-Efn)**2/((E-Efn)**2+wso2n**2)))
          Rso=1.1854-0.647*A**(-0.3333333)
          Rso=Rso*A**(1.0/3.0)
          Rwso=Rso
          aso=0.59
          awso=aso
	
  !        print *, "Potentials given as: depth, radius, diffuseness"
  !        print *, "Real central potential:"  
  !        print *, Vv, Rv, av
  !        print *, "Imaginary central potential:"
  !        print *, Wv, Rwv, awv
  !        print *, "Imaginary surface potential:"
  !        print *, Wd, Rwd, awd
  !        print *, "Real spin-orbit potential:"
  !        print *, Vso, Rso, aso
  !        print *, "Imaginary spin-orbit potential:"
  !        print *, Wso, Rwso, awso         

          do n=1,nmax
            r=n*StepSize
!           NUCLEAR POTENTIAL
            if (abs(Vv)>0.0) then
              Vnucl=Vv/(1+exp((r-Rv)/av))
            end if
            if (abs(Wv)>0.0) then
              Wnucl=Wv/(1+exp((r-Rwv)/awv))
            end if
            if (abs(Vd)>0.0) then
              Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
            end if
            if (abs(Wd)>0.0) then 
              Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            end if
            if (abs(Vso)>0.0) then
              Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
            end if
            if (abs(Wso)>0.0) then
              Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
            end if

            Upot1(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Spin1(n)=Vspin+i*Wspin

          end do

!         PROTON POTENTIAL
          v1p=59.3+21*(Nt-Zt)/A-0.024*A
          v2p=0.007067+(4.23E-6)*A
          v3p=(1.729E-5)+(1.136E-8)*A
          v4p=7E-9
          w1p=14.667+0.009629*A
          w2p=73.55+0.0795*A
          d1p=16+16*(Nt-Zt)/A
          d2p=0.0180+0.003802/(1+exp((A-156)/8))
          d3p=11.5
          vso1p=5.922+0.003*A
          vso2p=0.004
          wso1p=-3.1
          wso2p=160
          Efp=-8.4075+0.01378*A
          Rcoul=1.198+0.697*A**(-0.66666)+12.994*A**(-1.66666)
!          Rcoul=Rcoul*A**(1.0/3.0)
          Vcp=(1.73/Rcoul)*Zt*A**(-0.3333333)		
          Rcoul=Rcoul*A**(1.0/3.0)


!         VOLUME		
          Vv=-(v1p*(1-v2p*(E-Efp)+v3p*(E-Efp)**2-v4p*(E-Efp)**3)+vcp*v1p*(v2p-2*v3p*(E-Efp)+3*v4p*(E-Efp)**2))
          Wv=-(w1p*((E-Efp)**2/((E-Efp)**2+w2p**2)))
          Rv=1.3039-0.4054*A**(-0.3333333)
          Rv=rv*A**(1.0/3.0)
          Rwv=Rv
          av=0.6778-(1.487E-4)*A
          awv=av

!         SURFACE
          Vd=0.0
          Wd=-(d1p*((E-Efp)**2/((E-Efp)**2+d3p**2))*exp(-d2p*(E-Efp)))
          Rwd=1.3424-0.01585*A**(0.3333333)
          Rwd=rwd*A**(1.0/3.0)
          awd=0.5187+(5.205E-4)*A

!         SPIN-ORBIT
          Vso=-(vso1p*exp(-vso2p*(E-Efp)))
          Wso=-(wso1p*((E-Efp)**2/((E-Efp)**2+wso2p**2)))
          Rso=1.1854-0.647*A**(-0.33333333)
          Rso=Rso*A**(1.0/3.0)
          Rwso=Rso
          aso=0.59
          awso=aso
	
 !         print *, "Potentials given as: depth, radius, diffuseness"
 !         print *, "Real central potential:"
 !         print *, Vv, Rv, av
 !         print *, "Imaginary central potential:"
 !         print *, Wv, Rwv, awv
 !         print *, "Imaginary surface potential:"
 !         print *, Wd, Rwd, awd
 !         print *, "Real spin-orbit potential:"
 !         print *, Vso, Rso, aso
 !         print *, "Imaginary spin-orbit potential:"
 !         print *, Wso, Rwso, aso
 !         print *, "Coulomb Radius:"
 !         print *, Rcoul

          do n=1,nmax
            r=n*StepSize
!           NUCLEAR POTENTIAL
            if (abs(Vv)>0.0) then
              Vnucl=Vv/(1+exp((r-Rv)/av))
            end if
            if (abs(Wv)>0.0) then
              Wnucl=Wv/(1+exp((r-Rwv)/awv))
            end if
            if (abs(Vd)>0.0) then
              Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
            end if
            if (abs(Wd)>0.0) then 
              Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            end if
            if (abs(Vso)>0.0) then
              Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
            end if
            if (abs(Wso)>0.0) then
              Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
            end if
!           COULOMB POTENTIAL
            if(r<Rcoul) then
              Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
            end if
            if(r>=Rcoul)then
              Vcoul=(Zp*Zt*1.43997)/r
            end if

!            Upot1(n)=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)
            Upot2(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Coulomb(n)=Vcoul
            Spin2(n)=Vspin+i*Wspin

          end do

!       END KONING DELAROCHE
        end if

!       WhatPreDefPot = 1 = Koning-Delaroche
!       WhatPreDefPot = 2 = Chapel Hill
        if (WhatPreDefPot == 2) then		

!         REAL VOLUME
          Vv=-(52.9-13.1*((Nt-Zt)/A)+E*(-0.299))
          Rv=1.25-0.225/(A**(1.0/3.0))
          Rv=Rv*A**(1.0/3.0)
          av=0.69

!         IMAGINARY VOLUME
          Wv=-(7.8/(1+exp((35.0-E)/16.0)))
          Rwv=1.33-0.42/(A**(1.0/3.0))
          Rwv=Rwv*A**(1.0/3.0)
          awv=0.69

!         IMAGINARY SURFACE
          Wd=-((10.0-18.0*(Nt-Zt)/A)/(1+exp((E-Ec-36.0)/37.0)))
          Rwd=1.33-0.42/(A**(1.0/3.0))
          Rwd=Rwd*A**(1.0/3.0)
          awd=0.69
          Vd=0.0

!         REAL SPIN-ORBIT
          Vso=-5.9
          Rso=1.34-1.2/(A**(1.0/3.0))
          Rso=Rso*A**(1.0/3.0)
          aso=0.63
          Wso=0.0

!          write(*,*),'CH89 Potential (n,n)'
!          write(*,'(A10,4f6.1)'),'N,Z,A,E =',Nt,Zt,A,E
!          write(*,'(A11,3f8.3)'),'Vv,rv,av =',Vv,Rv,av
!          write(*,'(A13,3f8.3)'),'Wv,rwv,awv =',Wv,rwv,awv
!          write(*,'(A13,3f8.3)'),'Wd,rwd,awd =',Wd,Rwd,awd
!          write(*,'(A14,3f8.3)'),'Vso,rso,aso =',Vso,Rso,aso
	
          do n=1,nmax
            r=n*StepSize
!           NUCLEAR POTENTIAL
            if (abs(Vv)>0.0) then
              Vnucl=Vv/(1+exp((r-Rv)/av))
            end if
            if (abs(Wv)>0.0) then
              Wnucl=Wv/(1+exp((r-Rwv)/awv))
            end if
            if (abs(Vd)>0.0) then
              Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
            end if
            if (abs(Wd)>0.0) then 
              Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            end if
            if (abs(Vso)>0.0) then
              Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
            end if
            if (abs(Wso)>0.0) then
              Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
            end if

            Upot1(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Spin1(n)=Vspin+i*Wspin

          end do

          Rc=1.238*A**(1.0/3.0)+0.116
          Ec=1.73*Zt/Rc

!         REAL VOLUME
          Vv=-(52.9+13.1*((Nt-Zt)/A)+(E-Ec)*(-0.299))
          Rv=1.25-0.225/(A**(1.0/3.0))
          Rv=Rv*A**(1.0/3.0)
          av=0.69

!         IMAGINARY VOLUME
          Wv=-(7.8/(1+exp((35.0-E+Ec)/16.0)))
          Rwv=1.33-0.42/(A**(1.0/3.0))
          Rwv=Rwv*A**(1.0/3.0)
          awv=0.69

!         IMAGINARY SURFACE
          Wd=-((10.0+18.0*(Nt-Zt)/A)/(1+exp((E-Ec-36.0)/37.0)))
          Rwd=1.33-0.42/(A**(1.0/3.0))
          Rwd=Rwd*A**(1.0/3.0)
          awd=0.69
          Vd=0.0

!         REAL SPIN-ORBIT
          Vso=-5.9
          Rso=1.34-1.2/(A**(1.0/3.0))
          Rso=Rso*A**(1.0/3.0)
          aso=0.63
          Wso=0.0

          Rcoul=1.238+0.116/(A**(1.0/3.0))
          Rcoul=Rcoul*A**(1.0/3.0)

!          write(*,*),'CH89 Potential (p,p)'
!          write(*,'(A10,4f6.1)'),'N,Z,A,E =',Nt,Zt,A,E
!          write(*,'(A11,3f8.3)'),'Vv,rv,av =',Vv,Rv,av
!          write(*,'(A13,3f8.3)'),'Wv,rwv,awv =',Wv,rwv,awv
!          write(*,'(A13,3f8.3)'),'Wd,rwd,awd =',Wd,Rwd,awd
!          write(*,'(A14,3f8.3)'),'Vso,rso,aso =',Vso,Rso,aso
!          write(*,'(A8,f8.3)'),'rcoul =',Rcoul

          do n=1,nmax
            r=n*StepSize
!           NUCLEAR POTENTIAL
            if (abs(Vv)>0.0) then
              Vnucl=Vv/(1+exp((r-Rv)/av))
            end if
            if (abs(Wv)>0.0) then
              Wnucl=Wv/(1+exp((r-Rwv)/awv))
            end if
            if (abs(Vd)>0.0) then
              Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
            end if
            if (abs(Wd)>0.0) then 
              Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            end if
            if (abs(Vso)>0.0) then
              Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
            end if
            if (abs(Wso)>0.0) then
              Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
            end if
!           COULOMB POTENTIAL
            if(r<Rcoul) then
              Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
            end if
            if(r>=Rcoul)then
              Vcoul=(Zp*Zt*1.43997)/r
            end if

!            Upot1(n)=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)
            Upot2(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Coulomb(n)=Vcoul
            Spin2(n)=Vspin+i*Wspin

          end do

!       END CHAPEL HILL
        end if

!       NUCLEON POTENTIALS WERE EVALUATED AT HALF THE DEUTERON ENERGY
!       RESTORE THE ORIGINAL VALUE FOR THE BEAM ENERGY
        E=E*2.0

!     END PRE-DEFINED LOCAL NUCLEON POTENTIALS FOR THE DEUTERON
      end if

!   END LOCAL NUCLEON OPTICAL POTENTIALS FOR THE DEUTERON
    end if

! END DEUTERON SCATTERING STATE
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NUCLEON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! WhatSystem = 1 = Deuteron Scattering State
! WhatSystem = 2 = Nucleon Scattering State
  if (WhatSystem == 2) then
    Upot2(:)=0.0
    Spin2(:)=0.0

!   WhatPot = 1 = User defined local potential  
!   WhatPot = 2 = Pre-defined local potential
    if (WhatPot==1) then

        Vv=-ScatParameters(4,1)
        Rv=ScatParameters(4,2)*A**(1.0/3.0)
        av=ScatParameters(4,3)
        Wv=-ScatParameters(4,4)
        Rwv=ScatParameters(4,5)*A**(1.0/3.0)
        awv=ScatParameters(4,6)
        Vd=-ScatParameters(4,7)
        Rd=ScatParameters(4,8)*A**(1.0/3.0)
        ad=ScatParameters(4,9)
        Wd=-ScatParameters(4,10)
        Rwd=ScatParameters(4,11)*A**(1.0/3.0)
        awd=ScatParameters(4,12)
        Vso=-ScatParameters(4,13)
        Rso=ScatParameters(4,14)*A**(1.0/3.0)
        aso=ScatParameters(4,15)
        Wso=-ScatParameters(4,16)
        Rwso=ScatParameters(4,17)*A**(1.0/3.0)
        awso=ScatParameters(4,18)
        Rcoul=ScatParameters(4,19)*A**(1.0/3.0)

        Vnucl=0.0
        Wnucl=0.0
        Vsurf=0.0
        Wsurf=0.0
        Vspin=0.0
        Wspin=0.0

        do n=1,nmax
          r=n*StepSize
!         NUCLEAR POTENTIAL
          if (abs(Vv)>0.0) then
            Vnucl=Vv/(1+exp((r-Rv)/av))
          end if
          if (abs(Wv)>0.0) then
            Wnucl=Wv/(1+exp((r-Rwv)/awv))
          end if
          if (abs(Vd)>0.0) then
            Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
          end if
          if (abs(Wd)>0.0) then 
            Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
          end if
          if (abs(Vso)>0.0) then
            Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                  *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
          end if
          if (abs(Wso)>0.0) then
            Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                  *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
          end if
!         COULOMB POTENTIAL
          if(r<Rcoul) then
            Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
          end if
          if(r>=Rcoul)then
            Vcoul=(Zp*Zt*1.43997)/r
          end if

          Upot1(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
          Coulomb(n)=Vcoul
          Spin1(n)=Vspin+i*Wspin

        end do

!   END USER-DEFINED LOCAL POTENTIAL
    end if

!   WhatPot = 1 = User defined local potential  
!   WhatPot = 2 = Pre-defined local potential
    if (WhatPot == 2) then

!     WhatPreDefPot = 1 = Koning-Delaroche
!     WhatPreDefPot = 2 = Chapel Hill
      if (WhatPreDefPot == 1) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     KONING-DELAROCHE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       IF A NEUTRON PROJECTIL
        if (Zp == 0) then
          v1n=59.3-21*(Nt-Zt)/A-0.024*A
          v2n=0.007228-(1.48E-6)*A
          v3n=(1.994E-5)-(2E-8)*A
          v4n=7E-9
          w1n=12.195+0.0167*A
          w2n=73.55+0.0795*A
          d1n=16-16*(Nt-Zt)/A
          d2n=0.018+0.003802/(1+exp((A-156)/8))
          d3n=11.5
          vso1n=5.922+0.003*A
          vso2n=0.004
          wso1n=-3.1
          wso2n=160
          Efn=-11.2814+0.02646*A
		
!         VOLUME
          Vv=-(v1n*(1-v2n*(E-Efn)+v3n*(E-Efn)**2-v4n*(E-Efn)**3))
          Wv=-(w1n*((E-Efn)**2/((E-Efn)**2+w2n**2)))
          Rv=1.3039-0.4054*A**(-0.3333333)
          Rv=Rv*A**(1.0/3.0)
          Rwv=Rv
          av=0.6778-(1.487E-4)*A
          awv=av

!         SURFACE
          Vd=0.0
          Wd=-(d1n*((E-Efn)**2/((E-Efn)**2+d3n**2))*exp(-d2n*(E-Efn)))
          Rwd=1.3424-0.01585*A**(0.3333333)
          Rwd=Rwd*A**(1.0/3.0)
          awd=0.5446-(1.656E-4)*A

!         SPIN-ORBIT
          Vso=-(vso1n*exp(-vso2n*(E-Efn)))
          Wso=-(wso1n*((E-Efn)**2/((E-Efn)**2+wso2n**2)))
          Rso=1.1854-0.647*A**(-0.3333333)
          Rso=Rso*A**(1.0/3.0)
          Rwso=Rso
          aso=0.59
          awso=aso
	
!          write(*,*) "Neutrons"
!          write(*,*) "Energy"
!	  write(*,*) E
!          write(*,*) "Potentials given as: depth, radius, diffuseness"
!          write(*,*) "Real central potential:"  
!          write(*,*) Vv, Rv, av
!          write(*,*) "Imaginary central potential:"
!          write(*,*) Wv, Rwv, awv
!          write(*,*) "Imaginary surface potential:"
!          write(*,*) Wd, Rwd, awd
!          write(*,*) "Real spin-orbit potential:"
!          write(*,*) Vso, Rso, aso
!          write(*,*) "Imaginary spin-orbit potential:"
!          write(*,*) Wso, Rwso, awso         
!	  write(*,*) "Coulomb Redius:"
!	  write(*,*) Rcoul


          do n=1,nmax
            r=n*StepSize
!           NUCLEAR POTENTIAL
            if (abs(Vv)>0.0) then
              Vnucl=Vv/(1+exp((r-Rv)/av))
            end if
            if (abs(Wv)>0.0) then
              Wnucl=Wv/(1+exp((r-Rwv)/awv))
            end if
            if (abs(Vd)>0.0) then
              Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
            end if
            if (abs(Wd)>0.0) then 
              Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            end if
            if (abs(Vso)>0.0) then
              Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
            end if
            if (abs(Wso)>0.0) then
              Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
            end if

            Vcoul=0.0

            Upot1(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Coulomb(n)=Vcoul
            Spin1(n)=Vspin+i*Wspin

          end do

!       END OF NEUTRON KONING-DELAROCHE POTENTIAL
        end if

!       IF A PROTON PROJECTILE IN KONING-DELAROCHE POTENTIAL
        if (Zp == 1) then
          v1p=59.3+21*(Nt-Zt)/A-0.024*A
          v2p=0.007067+(4.23E-6)*A
          v3p=(1.729E-5)+(1.136E-8)*A
          v4p=7E-9
          w1p=14.667+0.009629*A
          w2p=73.55+0.0795*A
          d1p=16+16*(Nt-Zt)/A
          d2p=0.0180+0.003802/(1+exp((A-156)/8))
          d3p=11.5
          vso1p=5.922+0.003*A
          vso2p=0.004
          wso1p=-3.1
          wso2p=160
          Efp=-8.4075+0.01378*A
          Rcoul=1.198+0.697*A**(-0.66666)+12.994*A**(-1.66666)
!          Rcoul=Rcoul*A**(1.0/3.0)
          Vcp=(1.73/Rcoul)*Zt*A**(-0.3333333)		
          Rcoul=Rcoul*A**(1.0/3.0)

!         VOLUME		
          Vv=-(v1p*(1-v2p*(E-Efp)+v3p*(E-Efp)**2-v4p*(E-Efp)**3)+Vcp*v1p*(v2p-2*v3p*(E-Efp)+3*v4p*(E-Efp)**2))
          Wv=-(w1p*((E-Efp)**2/((E-Efp)**2+w2p**2)))
          Rv=1.3039-0.4054*A**(-0.3333333)
          Rv=Rv*A**(1.0/3.0)
          Rwv=Rv
          av=0.6778-(1.487E-4)*A
          awv=av

!         SURFACE
          Vd=0.0
          Wd=-(d1p*((E-Efp)**2/((E-Efp)**2+d3p**2))*exp(-d2p*(E-Efp)))
          Rwd=1.3424-0.01585*A**(0.3333333)
          Rwd=Rwd*A**(1.0/3.0)
          awd=0.5187+(5.205E-4)*A

!         SPIN-ORBIT
          Vso=-(vso1p*exp(-vso2p*(E-Efp)))
          Wso=-(wso1p*((E-Efp)**2/((E-Efp)**2+wso2p**2)))
          Rso=1.1854-0.647*A**(-0.33333333)
          Rso=Rso*A**(1.0/3.0)
          Rwso=Rso
          aso=0.59
          awso=aso
	
!          write(*,*) "Protons"
!	  write(*,*) "Energy Efp"
!	  write(*,*) E, Efp
!          write(*,*) "Potentials given as: depth, radius, diffuseness"
!          write(*,*) "Real central potential:"  
!          write(*,*) Vv, Rv, av
!          write(*,*) "Imaginary central potential:"
!          write(*,*) Wv, Rwv, awv
!          write(*,*) "Imaginary surface potential:"
!          write(*,*) Wd, Rwd, awd
!          write(*,*) "Real spin-orbit potential:"
!          write(*,*) Vso, Rso, aso
!          write(*,*) "Imaginary spin-orbit potential:"
!          write(*,*) Wso, Rwso, awso         
!	  write(*,*) "Coulomb Redius:"
!	  write(*,*) Rcoul
!	  write(*,*) "N Z A"
!	  write(*,*) Nt, Zt, A

          do n=1,nmax
            r=n*StepSize
!           NUCLEAR POTENTIAL
            if (abs(Vv)>0.0) then
              Vnucl=Vv/(1+exp((r-Rv)/av))
            end if
            if (abs(Wv)>0.0) then
              Wnucl=Wv/(1+exp((r-Rwv)/awv))
            end if
            if (abs(Vd)>0.0) then
              Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
            end if
            if (abs(Wd)>0.0) then 
              Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            end if
            if (abs(Vso)>0.0) then
              Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
            end if
            if (abs(Wso)>0.0) then
              Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
            end if
!           COULOMB POTENTIAL
            if(r<Rcoul) then
              Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
            end if
            if(r>=Rcoul)then
              Vcoul=(Zp*Zt*1.43997)/r
            end if

!            Upot1(n)=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)
            Upot2(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Coulomb(n)=Vcoul
            Spin2(n)=Vspin+i*Wspin

          end do

!       END OF PROTON KONING-DELAROCHE POTENTIAL
        end if

!     END KONING-DELAROCHE POTENTIAL 
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     CHAPEL HILL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     WhatPreDefPot = 1 = Koning-Delaroche
!     WhatPreDefPot = 2 = Chapel Hill
      if (WhatPreDefPot == 2) then

!       NEUTRON CHAPEL HILL POTENTIAL
        if (Zp == 0) then		

!         REAL VOLUME
          Vv=-(52.9-13.1*((Nt-Zt)/A)+E*(-0.299))
          Rv=1.25-0.225/(A**(1.0/3.0))
          Rv=Rv*A**(1.0/3.0)
          av=0.69

!         IMAGINARY VOLUME
          Wv=-(7.8/(1+exp((35.0-E)/16.0)))
          Rwv=1.33-0.42/(A**(1.0/3.0))
          Rwv=Rwv*A**(1.0/3.0)
          awv=0.69

!         IMAGINARY SURFACE
          Wd=-((10.0-18.0*(Nt-Zt)/A)/(1+exp((E-Ec-36.0)/37.0)))
          Rwd=1.33-0.42/(A**(1.0/3.0))
          Rwd=Rwd*A**(1.0/3.0)
          awd=0.69
          Vd=0.0

!         REAL SPIN-ORBIT
          Vso=-5.9
          Rso=1.34-1.2/(A**(1.0/3.0))
          Rso=Rso*A**(1.0/3.0)
          aso=0.63
          Wso=0.0

!          write(*,*),'CH89 Potential (n,n)'
!          write(*,'(A10,4f6.1)'),'N,Z,A,E =',Nt,Zt,A,E
!          write(*,'(A11,3f8.3)'),'Vv,rv,av =',Vv,Rv,av
!          write(*,'(A13,3f8.3)'),'Wv,rwv,awv =',Wv,rwv,awv
!          write(*,'(A13,3f8.3)'),'Wd,rwd,awd =',Wd,Rwd,awd
!          write(*,'(A14,3f8.3)'),'Vso,rso,aso =',Vso,Rso,aso
	
          do n=1,nmax
            r=n*StepSize
!           NUCLEAR POTENTIAL
            if (abs(Vv)>0.0) then
              Vnucl=Vv/(1+exp((r-Rv)/av))
            end if
            if (abs(Wv)>0.0) then
              Wnucl=Wv/(1+exp((r-Rwv)/awv))
            end if
            if (abs(Vd)>0.0) then
              Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
            end if
            if (abs(Wd)>0.0) then 
              Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            end if
            if (abs(Vso)>0.0) then
              Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
            end if
            if (abs(Wso)>0.0) then
              Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
            end if

            Upot1(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Spin1(n)=Vspin+i*Wspin

          end do

!       END NEUTRON CHAPEL HILL POTENTIAL
        end if

!       PROTON CHAPEL HILL POTENTIAL
        if (Zp == 1) then

          Rc=1.238*A**(1.0/3.0)+0.116
          Ec=1.73*Zt/Rc

!         REAL VOLUME
          Vv=-(52.9+13.1*((Nt-Zt)/A)+(E-Ec)*(-0.299))
          Rv=1.25-0.225/(A**(1.0/3.0))
          Rv=Rv*A**(1.0/3.0)
          av=0.69

!         IMAGINARY VOLUME
          Wv=-(7.8/(1+exp((35.0-E+Ec)/16.0)))
          Rwv=1.33-0.42/(A**(1.0/3.0))
          Rwv=Rwv*A**(1.0/3.0)
          awv=0.69

!         IMAGINARY SURFACE
          Wd=-((10.0+18.0*(Nt-Zt)/A)/(1+exp((E-Ec-36.0)/37.0)))
          Rwd=1.33-0.42/(A**(1.0/3.0))
          Rwd=Rwd*A**(1.0/3.0)
          awd=0.69
          Vd=0.0

!         REAL SPIN-ORBIT
          Vso=-5.9
          Rso=1.34-1.2/(A**(1.0/3.0))
          Rso=Rso*A**(1.0/3.0)
          aso=0.63
          Wso=0.0

          Rcoul=1.238+0.116/(A**(1.0/3.0))
          Rcoul=Rcoul*A**(1.0/3.0)

!          write(*,*),'CH89 Potential (p,p)'
!          write(*,'(A10,4f6.1)'),'N,Z,A,E =',Nt,Zt,A,E
!          write(*,'(A11,3f8.3)'),'Vv,rv,av =',Vv,Rv,av
!          write(*,'(A13,3f8.3)'),'Wv,rwv,awv =',Wv,rwv,awv
!          write(*,'(A13,3f8.3)'),'Wd,rwd,awd =',Wd,Rwd,awd
!          write(*,'(A14,3f8.3)'),'Vso,rso,aso =',Vso,Rso,aso
!          write(*,'(A8,f8.3)'),'rcoul =',Rcoul

          do n=1,nmax
            r=n*StepSize
!           NUCLEAR POTENTIAL
            if (abs(Vv)>0.0) then
              Vnucl=Vv/(1+exp((r-Rv)/av))
            end if
            if (abs(Wv)>0.0) then
              Wnucl=Wv/(1+exp((r-Rwv)/awv))
            end if
            if (abs(Vd)>0.0) then
              Vsurf=(4*Vd*exp((r-Rd)/ad))/(1+exp((r-Rd)/ad))**(2.0) 
            end if
            if (abs(Wd)>0.0) then 
              Wsurf=(4*Wd*exp((r-Rwd)/awd))/(1+exp((r-Rwd)/awd))**(2.0)  
            end if
            if (abs(Vso)>0.0) then
              Vspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Vso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rso)/aso))/((1+exp((r-Rso)/aso))**(2.0))
            end if
            if (abs(Wso)>0.0) then
              Wspin=(jp*(jp+1.0)-L*(L+1.0)-Ip*(Ip+1.0))*(Wso/(aso*r))*((hbarc/139.6)**2.0) &
                    *(exp((r-Rwso)/awso))/((1+exp((r-Rwso)/awso))**(2.0))
            end if
!           COULOMB POTENTIAL
            if(r<Rcoul) then
              Vcoul=((Zp*Zt*1.43997)/(2.0*Rcoul))*(3.0-(r**2.0)/(Rcoul**2.0)) 
            end if
            if(r>=Rcoul)then
              Vcoul=(Zp*Zt*1.43997)/r
            end if

!            Upot1(n)=Vnucl+Vsurf+Vspin+i*(Wnucl+Wsurf+Wspin)
            Upot2(n)=Vnucl+Vsurf+i*(Wnucl+Wsurf)
            Coulomb(n)=Vcoul
            Spin2(n)=Vspin+i*Wspin

          end do

!       END PROTON CHAPEL HILL POTENTIAL
        end if

!     END CHAPEL HILL POTENTIAL
      end if

!   END PRE-DEFINED NUCLEON OPTICAL POTENTIALS
    end if

! END SINGLE NUCLEON LOCAL POTENTIAL
  end if 

  return
  end subroutine
