  subroutine calc_source(ScatParameters,vpot,nmax,L,jp,Lmax,wf,SchEqConst,source,calc, &
                         npoints,Kernel,KernelRead,StepSize,integral,Rmatch,WhatSystem)
  implicit none

  integer::n,nmax,mm1,ifail,q,npoints,v,maxL,calc,Lmax,L,WhatPot,WhatPreDefPot,WhatSystem
  real(8)::StepSize,pi,hbarc,r,SchEqConst,rmatch,eta,time,beta,rprime,exponent
  real(8)::Zp,nmm,hl,kl,z,maximum,minimum,jp,jtot
!  real(8)::g(1:nmax,1:npoints)
  real(8)::ScatParameters(1:15,1:100)
  real(8)::points(npoints),weights(npoints)
  complex*16::i,logd,diff,arg,constant,k,integral(nmax),vpot(nmax),wf(nmax),csj(0:100),source(nmax)
  complex*16::Unl,Unl1,Unl2,UlocNL,UlocNL1,UlocNL2
  complex*16::Kernel(1:nmax,1:npoints)
  complex*16::KernelRead(1:nmax,1:nmax)

!  open(2000,file='TempPot.txt')

  maxL=max(Lmax,20)
  pi=acos(-1.0)
  i=cmplx(0.0,1.0)

  WhatPot=nint(ScatParameters(3,3))
  WhatPreDefPot=nint(ScatParameters(3,4))
  Zp=ScatParameters(1,3)

! WhatSystem = 1 = Deuteron scattering state
! WhatSystem = 2 = Nucleon scattering state
  if (WhatSystem == 1) then
    beta=ScatParameters(6,19)
! END DEUTERON SCATTERING STATE
  end if
  if (WhatSystem == 2) then
!   WhatPot = 1 = User defined nonlocal potential
!   WhatPot = 2 = Pre-defined nonlocal potential
!   WhatPot = 3 = Read in nonlocal potential
    if (WhatPot==1) then
      beta=ScatParameters(6,19)
    end if
    if (WhatPot==2) then
!      WhatPreDefPot = 1 = Perey-Buck
!      WhatPreDefPot = 2 = TPM
       if (WhatPreDefPot ==1 ) then
         beta=0.85
       end if
       if (WhatPreDefPot == 2) then
         if (nint(Zp)==0) then
           beta=0.90
         end if
         if (nint(Zp)==1) then
           beta=0.88
         end if
       end if
    end if
! END NUCLEON SCATTERING STATE
  end if

!  write(*,*) 'beta=',beta

! WhatPot = 1 = User defined nonlocal potential
! WhatPot = 2 = Pre-defined nonlocal potential
! WhatPot = 3 = Read in nonlocal potential
  if (WhatPot == 1 .or. WhatPot == 2) then
    do n=1,nmax
      r=n*StepSize
      rprime=0.0
      integral(n)=0.0 
!     SET UP MESH POINTS
!     INTEGRAL IS GAUSSIAN CENTERED ABOUT r=rprime WITH THE SPREAD IN rprime
!     USING 'npoints' MESH POINTS CENTERED ABOUT r=rprime AND EXTENDING 5 BETA ABOVE AND BELOW   
      minimum=max(r-5.0*beta,0.0)             
      maximum=min(r+5.0*beta,Rmatch)
      call gauss(minimum,maximum,npoints,points,weights)        

!     CALCULATE NON-LOCAL POTENTIAL: LOOP OVER r' TO DO INTEGRAL
      do q=1,npoints
        rprime=points(q)
        call nonlocal_scattering_potential(ScatParameters,Unl1,Unl2,L,jp,R,Rprime,UlocNL1,UlocNL2)

        UlocNL=UlocNL1+UlocNL2
        Unl=Unl1+Unl2

        v=int(points(q)/StepSize)
        if (v==0) then
          v=1
        end if
        if(calc==0) then
!         CALCULATE THE NON-LOCAL INTEGRAL: IF z>700 THEN j_L(-i*z) IS TOO BIG
          z=2.0*r*rprime/(beta**2.0)
          arg=cmplx(0.0,-z)
          if(z<=700) then
            call csphjy(maxL,arg,nmm,csj)
            kl=2.0*(i**(L))*z*csj(L)
            exponent=(r**2.0+rprime**2.0)/(beta**2.0)
            hl=(1.0/(beta*pi**(1.0/2.0)))*(exp(-exponent))*kl        
          end if
          if(z>700)then
            hl=(1.0/(beta*pi**(1.0/2.0)))*(exp(-((r-rprime)/(beta))**2.0))
          end if
          Kernel(n,q)=hl    
        end if 
        integral(n)=integral(n)+weights(q)*wf(v)*Kernel(n,q)*Unl
!     END q LOOP 
      end do
         
      source(n)=-SchEqConst*(integral(n)+UlocNL*wf(n)-vpot(n)*wf(n))  

!   END n LOOP
    end do

  end if

! WhatPot = 1 = User defined nonlocal potential
! WhatPot = 2 = Pre-defined nonlocal potential
! WhatPot = 3 = Read in nonlocal potential
  if (WhatPot == 3) then
    do n=1,nmax
      r=n*StepSize
      rprime=0.0
      integral(n)=0.0 

!     CALCULATE NON-LOCAL POTENTIAL: LOOP OVER r' TO DO INTEGRAL
      do q=1,nmax
        rprime=q*StepSize
        call nonlocal_scattering_potential(ScatParameters,Unl1,Unl2,L,jp,R,Rprime,UlocNL1,UlocNL2)
        UlocNL=UlocNL1+UlocNL2
        integral(n)=integral(n)+StepSize*wf(q)*KernelRead(n,q)
      end do
      source(n)=-SchEqConst*(integral(n)+UlocNL*wf(n)-vpot(n)*wf(n)) 
    end do
  end if

  return
  end subroutine calc_source
