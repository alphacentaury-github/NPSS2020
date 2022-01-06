  subroutine calc_binding_source(BoundParameters,vpot,nmax,L,wfin,wfout,SchEqConst, &
       integralin,integralout,calc,npoints,Kernel,KernelRead,StepSize,Rmatch)
  implicit none

  integer::nmax,calc,npoints,n,q,maxL,v,L,WhatPot
  real(8)::StepSize,pi,beta,r,rprime,maximum,minimum,z,hl,kl,exponent,rmatch,nmm,SchEqConst,Unl,UlocNL
  real(8)::vpot(1:nmax),Kernel(1:nmax,1:npoints),KernelRead(1:nmax,1:nmax)
  real(8)::integralin(1:nmax),integralout(1:nmax),wfin(1:nmax),wfout(1:nmax),integral(1:nmax)
  real(8)::BoundParameters(1:15,1:100)
  real(8)::points(npoints),weights(npoints)
  complex*16::i,csj(0:30),arg

!  open(3000,file='NLpotBound.txt')

  pi=acos(-1.0)
  i=cmplx(0.0,1.0)

  maxL=max(20,L)
  nmm=0

  integralin(:)=0.0
  integralout(:)=0.0

  beta=BoundParameters(6,10)
  WhatPot=nint(BoundParameters(3,2))

! WhatPot = 1 = User defined nonlocal potential
! WhatPot = 2 = Read in nonlocal potential
  if (WhatPot == 1) then
    do n=1,nmax
      r=n*StepSize
      rprime=0.0
      integralin(n)=0.0
      integralout(n)=0.0
      minimum=max(r-5.0*beta,0.0)
      maximum=max(r+5.0*beta,rmatch)
      call gauss(minimum,maximum,npoints,points,weights)  
      do q=1,npoints
        rprime=points(q)
        call nonlocal_binding_potential(BoundParameters,Unl,r,rprime,UlocNL,StepSize)
        v=nint(points(q)/StepSize)
        if (v==0) then
          v=1
        end if
        if(calc==0) then
!         CALCULATE THE NON-LOCAL INTEGRAL
!         IF z>700 THEN j_L(-i*z) IS TOO BIG
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

        integralin(n)=integralin(n)+weights(q)*wfin(v)*Kernel(n,q)*Unl
        integralout(n)=integralout(n)+weights(q)*wfout(v)*Kernel(n,q)*Unl
     
!        if (calc==0) then
!          write(3000,'(2F6.2,F11.4)') n*StepSize,q*StepSize,real(gg(n,q)*Unl)
!        end if

      end do

      integralin(n)=-SchEqConst*(integralin(n)+UlocNL*wfin(n)-vpot(n)*wfin(n))  
      integralout(n)=-SchEqConst*(integralout(n)+UlocNL*wfout(n)-vpot(n)*wfout(n))

    end do
  end if

! WhatPot = 1 = User defined nonlocal potential
! WhatPot = 2 = Read in nonlocal potential
  if (WhatPot == 2) then
    do n=1,nmax
      r=n*StepSize
      rprime=0.0
      integralin(n)=0.0
      integralout(n)=0.0  
      do q=1,nmax
        rprime=q*StepSize
        call nonlocal_binding_potential(BoundParameters,Unl,r,rprime,UlocNL,StepSize)
        integralin(n)=integralin(n)+StepSize*wfin(q)*KernelRead(n,q)
        integralout(n)=integralout(n)+StepSize*wfout(q)*KernelRead(n,q)
      end do

      integralin(n)=-SchEqConst*(integralin(n)+UlocNL*wfin(n)-vpot(n)*wfin(n))  
      integralout(n)=-SchEqConst*(integralout(n)+UlocNL*wfout(n)-vpot(n)*wfout(n))

    end do
  end if

  return
  end subroutine
