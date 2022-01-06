  subroutine scattering_state(ScatParameters,common,accuracy,ReturnWFs,Lmax,NmaxReturn,k,eta, &
                              DeuteronBoundWF,Vnp,SH,SH_index_max,SH_Lmax,SH_step, &
                              Jp_Max,ReducedMass,Ecm,printing,Directory)
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::n,NmaxLong,nmatch,ChargeTarget,ChargeProjectile,PotType,NmaxReturn
  integer::NonLoc,ParityTarget,ParityProjectile,L,Lmax,IFAIL,MM1
  integer::calc,npoints,SH_index_max,SH_Lmax,q,WhatPot
  integer::nn,NmaxShort,Jp_Max,WhatSystem,LocPotType,printing(1:100)
  integer::print_DeuteronLocalIntegral,print_NucleonLocalIntegral
  integer::print_DeuteronNonlocalIntegral,print_NucleonNonlocalIntegral
  integer::print_DeuteronLocalSmatrix,print_NucleonLocalSmatrix
  integer::print_DeuteronNonlocalSmatrix,print_NucleonNonlocalSmatrix
  real(8)::r,time,E,diff,hbarc,Rmax,Rmatch,eta,convergence,Ecm
  real(8)::MassTarget,MassProjectile,ReducedMass,SchEqConst,SpinTarget,SpinProjectile,Elab,jpmax
  real(8)::pi,maxL,jp,StepSizeShort,RealPart,ImagPart,Rprime
  real(8)::dWF_Send(NmaxReturn),Vnp(NmaxReturn),DeuteronBoundWF(NmaxReturn)
  real(8)::SH_step,SH(0:SH_Lmax,1:SH_index_max)
  real(8)::SlopeReal,SlopeImag,OldStepSize,StepSize,StepSizeLong,ReductionFactor
  real(8)::FC(0:Lmax),FCP(0:Lmax),GC(0:Lmax),GCP(0:Lmax)
  real(8),allocatable::Coulomb(:)
  real(8)::ScatParameters(1:15,1:100),common(1:100),accuracy(1:100)
  complex*16::i,logd,wfext,constant,k,Hplus,Hplusp,Hminus,Hminusp,Rmatrix
  complex*16::Rold,csj(0:50),arg
  complex*16::ReturnWFs(0:Lmax,0:Jp_Max,1:NmaxReturn),SM(0:Lmax,0:Jp_Max)
  complex*16,allocatable::vpot(:),wf(:),source(:),Kernel(:,:),KernelRead(:,:)
  complex*16,allocatable::Upot1(:),Upot2(:),wfSend(:),TempWF(:),NLint(:),Spin1(:),Spin2(:)
  character(LEN=50)::Directory,extension
  character(LEN=100)::filename

  print_DeuteronLocalIntegral=printing(7)
  print_NucleonLocalIntegral=printing(8)
  print_DeuteronNonlocalIntegral=printing(9)
  print_NucleonNonlocalIntegral=printing(10)
  print_DeuteronLocalSmatrix=printing(11)
  print_NucleonLocalSmatrix=printing(12)
  print_DeuteronNonlocalSmatrix=printing(13)
  print_NucleonNonlocalSmatrix=printing(14)

  if (print_DeuteronLocalIntegral == 1) then
!    open(300,file='DeuteronLocalIntegral.txt')
    extension='DeuteronLocalIntegral.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(300,file=trim(filename))
  end if
  if (print_NucleonLocalIntegral == 1) then
!    open(301,file='NucleonLocalIntegral.txt')
    extension='NucleonLocalIntegral.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(301,file=trim(filename))
  end if
  if (print_DeuteronNonlocalIntegral == 1) then
!    open(302,file='DeuteronNonlocalIntegral.txt')
    extension='DeuteronNonlocalIntegral.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(302,file=trim(filename))
  end if
  if (print_NucleonNonlocalIntegral == 1) then
!    open(303,file='NucleonNonlocalIntegral.txt')
    extension='NucleonNonlocalIntegral.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(303,file=trim(filename))
  end if
  if (print_DeuteronLocalSmatrix == 1) then
!    open(304,file='DeuteronLocalSmatrix.txt')
    extension='DeuteronLocalSmatrix.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(304,file=trim(filename))
  end if
  if (print_NucleonLocalSmatrix == 1) then
!    open(305,file='NucleonLocalSmatrix.txt')
    extension='NucleonLocalSmatrix.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(305,file=trim(filename))
  end if
  if (print_DeuteronNonlocalSmatrix == 1) then
!    open(306,file='DeuteronNonlocalSmatrix.txt')
    extension='DeuteronNonlocalSmatrix.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(306,file=trim(filename))
  end if
  if (print_NucleonNonlocalSmatrix == 1) then
!    open(307,file='NucleonNonlocalSmatrix.txt')
    extension='NucleonNonlocalSmatrix.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(307,file=trim(filename))
  end if

  WhatSystem=nint(ScatParameters(2,1))
  Rmax=common(3)
  WhatPot = ScatParameters(3,3)
!  write(*,*) 'WhatPot = ', WhatPot

! WhatPot = 1 = User defined nonlocal potential
! WhatPot = 2 = Pre-defined nonlocal potential
! WhatPot = 3 = Read in nonlocal potential
  if (WhatPot == 3) then
    open(308,file='NLpotScat.txt')
  end if

! WhatSystem = 1 = Deuteron scattering state
! WhatSystem = 2 = Nucleon scattering state
  if (WhatSystem == 1) then
    write(*,*) ''
    write(*,*) 'Deuteron Scattering State'
    NonLoc=nint(ScatParameters(3,1))
    PotType=nint(ScatParameters(3,2))
    LocPotType=nint(ScatParameters(2,2))
  else
    write(*,*) ''
    write(*,*) 'Nucleon Scattering State'
    NonLoc=nint(ScatParameters(3,1))
  end if

! PotType = 1 = Two-body potential
! PotType = 2 = Adiabatic potential
  if (WhatSystem == 1 .and. PotType == 2 .and. Nonloc==1) then
    StepSizeLong=accuracy(11)
    StepSizeShort = common(2)
    NmaxLong=nint(Rmax/StepSizeLong)
    NmaxShort=nint(Rmax/StepSizeShort)
  else
    StepSizeLong=common(2)
    NmaxLong=nint(Rmax/StepSizeLong)
    StepSizeShort=StepSizeLong
    NmaxShort=NmaxLong
  end if

!  write(*,*) 'WhatSystem, Ecm = ', WhatSystem, Ecm
  write(*,*) 'Ecm = ', Ecm

  Rmatch=Rmax 
  MassProjectile=ScatParameters(1,1)
  MassTarget=ScatParameters(1,2)
  ChargeProjectile=ScatParameters(1,3)
  ChargeTarget=ScatParameters(1,4)
  SpinProjectile=ScatParameters(1,5)
  SpinTarget=ScatParameters(1,6)
  ParityProjectile=ScatParameters(1,7)
  ParityTarget=ScatParameters(1,8)
  ReducedMass=((MassProjectile*MassTarget)/(MassProjectile+MassTarget))*accuracy(1)
  convergence=accuracy(10)
  npoints=nint(accuracy(2))

  Elab=((MassTarget+MassProjectile)/(MassTarget))*Ecm

!  write(*,*) ''
!  write(*,*) 'WhatSystem = ', WhatSystem
!  write(*,*) 'PotType = ', PotType
!  write(*,*) 'Rmatch = ', Rmatch
!  write(*,*) 'MassProjectile = ', MassProjectile
!  write(*,*) 'MassTarget = ', MassTarget
!  write(*,*) 'ChargeProjectile = ', ChargeProjectile
!  write(*,*) 'ChargeTarget = ', ChargeTarget
!  write(*,*) 'SpinProjectile = ', SpinProjectile
!  write(*,*) 'SpinTarget = ', SpinTarget
!  write(*,*) 'ParityProjectile = ', ParityProjectile
!  write(*,*) 'ParityTarget = ', ParityTarget
!  write(*,*) 'Ecm = ', Ecm
!  write(*,*) 'ReducedMass = ', ReducedMass
!  write(*,*) 'convergence = ', convergence
!  write(*,*) 'NonLoc = ', NonLoc
!  write(*,*) 'npoints = ', npoints
!  write(*,*) 'Rmax = ',Rmax
!  write(*,*) 'StepSizeLong = ', StepSizeLong
!  write(*,*) 'StepSizeShort = ', StepSizeShort
!  write(*,*) 'NmaxLong = ', NmaxLong
!  write(*,*) 'NmaxShort = ', NmaxShort
!  write(*,*) 'LocPotType = ', LocPotType

  IFAIL=0
  MM1=1
  maxL=real(Lmax)

  pi=acos(-1.0)
  i=cmplx(0.0,1.0)
  hbarc=197.3269718
  SchEqConst=(2.0*ReducedMass)/(hbarc**2.0)
  k=sqrt(SchEqConst*Ecm)
  eta=(ChargeTarget*ChargeProjectile*1.43997*ReducedMass)/((hbarc**2)*k)

!  write(*,*) 'eta,k=',eta,k
  
  jpmax=real(Lmax+SpinProjectile)

!  write(*,*) 'Lmax,jpmax =',Lmax,jpmax

! DEFINE ALL OF THE ARRAYS
  allocate(vpot(1:NmaxLong))
  allocate(wf(1:int(NmaxShort)))
  allocate(wfSend(1:int(NmaxShort)))
  allocate(TempWF(1:int(NmaxShort)))
  allocate(source(1:NmaxLong))
  allocate(Kernel(1:NmaxLong,1:npoints))
  allocate(KernelRead(1:NmaxLong,1:NmaxLong))
  allocate(Upot1(1:NmaxLong))
  allocate(Upot2(1:NmaxLong))
  allocate(Coulomb(1:NmaxLong))
  allocate(NLint(1:NmaxLong))
  allocate(Spin1(1:NmaxLong))
  allocate(Spin2(1:NmaxLong))

  source(:)=0.0

! LOOP OVER ALL PARTIAL WAVES
  do L=0,Lmax

    do Jp=(L-SpinProjectile),(L+SpinProjectile)
      if (L==0 .and. Jp<0.001) then
!       DO NOTHING
      else if (Jp>-0.001) then
        write(*,*) ''
        write(*,*) 'L Jp =',L,Jp

!       CALCULATE THE COULOMB FUNCTIONS AND THEIR DERIVATIVES AT THE MATCHING RADIUS
!       THIS IS NEEDED FOR CALCULATING THE R-MATRIX ELEMENT
        CALL COULFG(k*NmaxLong*StepSizeLong,eta,0,maxL,FC,GC,FCP,GCP,1,0,IFAIL,MM1)

!       CALCULATE HANKEL FUNCTIONS AND THEIR DERIVATIVES AT THE MATCHING RADIUS
        Hplus=GC(L)+i*FC(L)
        Hminus=GC(L)-i*FC(L)
        Hplusp=GCP(L)+i*FCP(L)
        Hminusp=GCP(L)-i*FCP(L)
 
        call local_scattering_potential(ScatParameters,Upot1,Upot2,NmaxLong,L,jp, &
                                        Coulomb,Spin1,Spin2,StepSizeLong,Elab)
       
!       IF DEUTERON SCATERING STATE, CALCULATE LOCAL ADIABATIC POTENTIAL
        if(WhatSystem == 1 .and. LocPotType == 2) then
          call local_adiabatic(DeuteronBoundWF,NmaxLong,Upot1,Upot2,StepSizeLong,Rmax)
        end if

!       CALCULATE THE OPTICAL POTENTIAL
        vpot(:)=Upot1(:)+Upot2(:)+Coulomb(:)+Spin1(:)+Spin2(:)

!       SOLVE LOCAL DIFFERENTIAL EQUATION AND CALCULATE R-MATRIX AND S-MATRIX ELEMENTS     
        source(:)=0.0
        call nm(L,NmaxLong,StepSizeLong,vpot,SchEqConst,Ecm,logd,source,wf)
        Rmatrix=(1.0/(real(NmaxLong)*StepSizeLong*logd))
        SM(int(L),int(2*jp))=(Hminus-(k*NmaxLong*StepSizeLong)*Rmatrix*(Hminusp)) &
                                         /(Hplus-k*NmaxLong*StepSizeLong*Rmatrix*(Hplusp))
        write(*,*) 'SM = ', SM(int(L),int(2*jp))

        if (WhatSystem==1) then
          if (print_DeuteronLocalSmatrix == 1) then
            write(304,*), L,Jp,real(SM(int(L),int(2*jp))),aimag(SM(int(L),int(2*jp)))
          end if
        end if
        if (WhatSystem==2) then
          if (print_NucleonLocalSmatrix == 1) then
            write(305,*), L,Jp,real(SM(int(L),int(2*jp))),aimag(SM(int(L),int(2*jp)))
          end if
        end if

!       FIND THE PROPER NORMALIZATION CONSTANT FOR THE EXTERNAL LOCAL WAVE FUNCTION         
        constant=0.0
        do n=(NmaxLong-99),(NmaxLong-10)
          r=n*StepSizeLong
          CALL COULFG(k*r,eta,0,maxL,FC,GC,FCP,GCP,1,0,IFAIL,MM1)
          Hplus=GC(L)+i*FC(L)
          Hminus=GC(L)-i*FC(L)      
          wfext=(i/2.0)*(Hminus-SM(int(L),int(2*jp))*Hplus)
          constant=constant+wfext/wf(n)
        end do
        constant=constant/90.0
!        write(*,*) 'local constant =',constant

        do n=1,NmaxLong
          ReturnWFs(L,nint(2*Jp),n)=constant*wf(n)
        end do

!       PRINT LOCAL ADIABATIC POTENTIAL AND ADIABATIC SOURCE
!         'LOCAL ADIABATIC SOURCE' IS <PHI|Vnp(Up+Un)|PHI>*PSI

        if (WhatSystem==1) then
          if (print_DeuteronLocalIntegral == 1) then
            write(300,*), '#',L,Jp
            do n=1,NmaxLong
              write(300,*) n*StepSizeLong,real((Upot1(n)+Upot2(n))*constant*wf(n)),aimag((Upot1(n)+Upot2(n))*constant*wf(n))
            end do        
            write(300,*) '&'
          end if
        end if
        if (WhatSystem==2) then
          if (print_NucleonLocalIntegral == 1) then
            write(301,*), '#',L,Jp
            do n=1,NmaxLong
              write(301,*) n*StepSizeLong,real((Upot1(n)+Upot2(n))*constant*wf(n)),aimag((Upot1(n)+Upot2(n))*constant*wf(n))
            end do        
            write(301,*) '&'
          end if
        end if

!       CALCULATE THE COULOMB FUNCTIONS AND THEIR DERIVATIVES AT THE MATCHING RADIUS
!       THIS IS NEEDED FOR CALCULATING THE R-MATRIX ELEMENT
        CALL COULFG(k*NmaxLong*StepSizeLong,eta,0,maxL,FC,GC,FCP,GCP,1,0,IFAIL,MM1)

!       CALCULATE HANKEL FUNCTIONS AND THEIR DERIVATIVES AT THE MATCHING RADIUS
        Hplus=GC(L)+i*FC(L)
        Hminus=GC(L)-i*FC(L)
        Hplusp=GCP(L)+i*FCP(L)
        Hminusp=GCP(L)-i*FCP(L)
       
        if (WhatSystem == 1) then
          do n=1,NmaxShort
            dWF_Send(n)=DeuteronBoundWF(n)/(n*StepSizeShort)
          end do
        end if

        if (WhatSystem == 2) then
          do n=1,NmaxLong
            wfSend(n)=wf(n)
          end do
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       CALCULATE wfsend WHICH IS SENT TO THE NONLOCAL ADIABATIC SUBROUTINE
!         THIS PART OF CODE TRANSFORMS wf(n) INTO wfsend(n) WHERE
!         wf(n) IS IN STEPS OF h>0.01 GIVEN IN INPUT FILE, AND wfsend(n) IS IN 
!         STEPS OF 0.01. THIS IS NECESSARY SO THAT THE INTEGRALS IN THE 
!         NONLOCAL ADIABATIC POTENTIAL CAN BE CALCULATED ACCURATELY WHILE THE 
!         SOURCE CAN STILL BE CALCULATED IN STEPS OF h>0.01. THE SOURCE CALCULATED
!         IN STEPS OF h>0.01 IS THEN SENT TO THE NUMEROV METHOD SUBROUTINE
!         IN ORDER TO CALCULATE wf(n) IN STEPS OF h>0.01.
!         NMAXLONG IS FOR STEPSIZE>0.01 AND NMAXSHORT IS FOR STEPSIZE=0.01

!       CALCULATE 

        if (NonLoc==1) then
!         GENERALIZE THE ARGUMENT OF THIS IF STATEMENT
          if (abs(StepSizeLong-StepSizeShort)>0.0001) then

            SlopeReal=Real(wf(1))/StepSizeLong
            SlopeImag=aimag(wf(2))/StepSizeLong

            OldStepSize=StepSizeLong
            NmaxShort=nint(Rmax/StepSizeShort)
            ReductionFactor=nint(OldStepSize/StepSizeShort)

            do n=1,ReductionFactor
              TempWF(ReductionFactor-n+1)=SlopeReal*((ReductionFactor-n+1)*OldStepSize/ReductionFactor) &
                                          +i*SlopeImag*((ReductionFactor-n+1)*OldStepSize/ReductionFactor)
            end do

            do n=2,NmaxLong
              SlopeReal=(real(wf(n))-real(wf(n-1)))/OldStepSize
              SlopeImag=(aimag(wf(n))-aimag(wf(n-1)))/OldStepSize
              do nn=1,ReductionFactor
                TempWF(nint(ReductionFactor*n-(nn-1)))=SlopeReal*((ReductionFactor-nn+1)*OldStepSize &
                                                 /ReductionFactor)+i*SlopeImag*((ReductionFactor-nn+1) &
                                                 *OldStepSize/ReductionFactor)+wf(n-1)  
              end do
            end do
            wfSend(:)=TempWF(:)
          end if
          if ((StepSizeLong-StepSizeShort)<0.0001) then
            NmaxShort=NmaxLong
            do n=1,NmaxLong
              wfsend(n)=wf(n)
            end do
          end if

        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       SOLVE NONLOCAL EQUATION
        if (NonLoc==1) then
          diff=2.0
          calc=0

!         WhatPot = 1 = User defined nonlocal potential
!         WhatPot = 2 = Pre-defined nonlocal potential
!         WhatPot = 3 = Read in nonlocal potential
          if (WhatPot == 3) then
            do n=1,NmaxLong
              do q=1,NmaxLong
                read(308,*),R,Rprime,RealPart,ImagPart
                KernelRead(n,q)=RealPart+i*ImagPart
              end do
            end do
          end if
       
!         ITERATE UNTIL R-MATRIX AGREES TO DESIRED LEVEL OF ACCURACY
          do while (abs(diff)>convergence)

!           IF NUCLEON SCATERING STATE, CALCULATE NONLOCAL NUCLEON SOURCE             
            if (WhatSystem == 2) then
   
              do n=1,NmaxLong
                wfSend(n)=wf(n)
              end do
      
!              call calc_source(ScatParameters,Vpot,NmaxLong,L,jp,Lmax,wfSend, &
!                               SchEqConst,source,calc,npoints,g,StepSizeLong,NLint,Rmatch,WhatSystem)
              call calc_source(ScatParameters,Vpot,NmaxLong,L,jp,Lmax,wfSend, &
                               SchEqConst,source,calc,npoints,Kernel,KernelRead,StepSizeLong,NLint,Rmatch,WhatSystem)
            end if
  
!           IF DEUTERON SCATERING STATE WITH A NONLOCAL DEUTERON OPTICAL POTENTIAL
!           PotType = 1 = Two-body scattering potential
!           PotType = 2 = Adiabatic potential
            if (WhatSystem == 1 .and. PotType == 1) then
   
              do n=1,NmaxLong
                wfSend(n)=wf(n)
              end do
      
              call calc_source(ScatParameters,Vpot,NmaxLong,L,jp,Lmax,wfSend, &
                               SchEqConst,source,calc,npoints,Kernel,KernelRead,StepSizeLong,NLint,Rmatch,WhatSystem)
            end if

!           IF DEUTERON SCATERING STATE, CALCULATE NONLOCAL ADIABATIC SOURCE
            if (WhatSystem == 1 .and. PotType == 2) then
  
!             SOURCE IS THE ADIABATIC INTEGRAL PLUS THE LOCAL PART OF THE NONLOCAL EQUATION MINUS
!               THE INITIAL LOCAL POTENTIAL
!             NLint IS JUST THE INTEGRAL TERM

              call nonlocal_adiabatic(ScatParameters,accuracy,NmaxShort,dWF_Send,wfSend,NmaxLong,Source,Vpot,SchEqConst,L, &
                                      SH,SH_step,Coulomb,Vnp,StepSizeShort,StepSizeLong, &
                                      SH_index_max,SH_Lmax,NLint,jp)
            end if
  
            calc=1
            call nm(L,NmaxLong,StepSizeLong,Vpot,SchEqConst,Ecm,logd,source,wf)
            Rold=Rmatrix
            Rmatrix=(1.0/(real(NmaxLong)*StepSizeLong*logd))
            diff=abs((Rold-Rmatrix)/Rmatrix)
            call cpu_time(time)
            write(*,*) 'diff,time=',diff,time

!           CALCULATE wfsend WHICH IS SENT TO THE NONLOCAL ADIABATIC SUBROUTINE
!             THIS PART OF CODE TRANSFORMS wf(n) INTO wfsend(n) WHERE
!             wf(n) IS IN STEPS OF h>0.01 GIVEN IN INPUT FILE, AND wfsend(n) IS IN 
!             STEPS OF 0.01. THIS IS NECESSARY SO THAT THE INTEGRALS IN THE 
!             NONLOCAL ADIABATIC POTENTIAL CAN BE CALCULATED ACCURATELY WHILE THE 
!             SOURCE CAN STILL BE CALCULATED IN STEPS OF h>0.01. THE SOURCE CALCULATED
!             IN STEPS OF h>0.01 IS THEN SENT TO THE NUMEROV METHOD SUBROUTINE
!             IN ORDER TO CALCULATE wf(n) IN STEPS OF h>0.01.
!             NMAX IS FOR h>0.01 AND NMAXSEND IS FOR 0.01 STEPSIZE
            if (WhatSystem == 1) then
!             GENERALIZE THE ARGUMENT OF THIS IF STATEMENT
              if ((StepSizeLong-StepSizeShort)>0.0001) then
 
                SlopeReal=Real(wf(1))/StepSizeLong
                SlopeImag=aimag(wf(1))/StepSizeLong
 
                OldStepSize=StepSizeLong
                NmaxShort=nint(Rmax/StepSizeShort)
                ReductionFactor=nint(OldStepSize/StepSizeShort)
  
                do n=1,ReductionFactor
                  TempWF(ReductionFactor-n+1)=SlopeReal*((ReductionFactor-n+1) &
                                              *OldStepSize/ReductionFactor) &
                                              +i*SlopeImag*((ReductionFactor-n+1) &
                                              *OldStepSize/ReductionFactor)
                end do
  
                do n=2,NmaxLong
                  SlopeReal=(real(wf(n))-real(wf(n-1)))/OldStepSize
                  SlopeImag=(aimag(wf(n))-aimag(wf(n-1)))/OldStepSize
                  do nn=1,ReductionFactor
                    TempWF(ReductionFactor*n-(nn-1))=SlopeReal*((ReductionFactor-nn+1) &
                                                     *OldStepSize/ReductionFactor) &
                                                     +i*SlopeImag*((ReductionFactor-nn+1) &
                                                     *OldStepSize/ReductionFactor)+wf(n-1)  
                  end do
                end do
                wfSend(:)=TempWF(:)
              end if
!             GENERALIZE THE ARGUMENT OF THIS IF STATEMENT
              if ((StepSizeLong-StepSizeShort)<0.0001) then
                NmaxShort=NmaxLong
                do n=1,NmaxLong
                  wfsend(n)=wf(n)
                end do
              end if
            end if

!         END DIFF LOOP
          end do
          diff=2.0         
      
          SM(int(L),int(2*jp))=(Hminus-(k*NmaxLong*StepSizeLong)*Rmatrix*(Hminusp)) &
                                           /(Hplus-k*NmaxLong*StepSizeLong*Rmatrix*(Hplusp))
          write(*,*) 'SM = ', SM(L,nint(2*jp))

          if (WhatSystem==1) then
            if (print_DeuteronNonlocalSmatrix == 1) then
              write(306,*), L,Jp,real(SM(int(L),int(2*jp))),aimag(SM(int(L),int(2*jp)))
            end if
          end if
          if (WhatSystem==2) then
            if (print_NucleonNonlocalSmatrix == 1) then
              write(307,*), L,Jp,real(SM(int(L),int(2*jp))),aimag(SM(int(L),int(2*jp)))
            end if
          end if

!         FIND THE PROPER NORMALIZATION CONSTANT FOR THE EXTERNAL NONLOCAL WAVE FUNCTION         
          constant=0.0
          do n=(NmaxLong-99),(NmaxLong-10)
            r=n*StepSizeLong
            CALL COULFG(k*r,eta,0,maxL,FC,GC,FCP,GCP,1,0,IFAIL,MM1)
            Hplus=GC(L)+i*FC(L)
            Hminus=GC(L)-i*FC(L)      
            wfext=(i/2.0)*(Hminus-SM(int(L),int(2*jp))*Hplus)
            constant=constant+wfext/wf(n)
          end do
          constant=constant/90.0
!          write(*,*) 'nonlocal constant=',constant
!          write(*,*) ''

          do n=1,NmaxShort
            ReturnWFs(L,nint(2*Jp),n)=constant*wfSend(n)   
          end do

!         PRINT THE DEUTERON AND PROTON NONLOCAL SOURCE AND NONLOCAL INTEGRAL
!          'SOURCE' IS S(R) IN THE EQUATION y''(R)+f(R)y(R)+S(R)=0
!          'NONLOCAL INTEGRAL' IS  [INTEGRAL V(R,R')PSI(R')dR']
          if (WhatSystem==1) then
            if (print_DeuteronNonlocalIntegral == 1) then
              write(302,*) '#',L,Jp
              do n=1,NmaxLong
                write(302,*) n*StepSizeLong,real(constant*NLint(n)),aimag(constant*NLint(n))
              end do
              write(302,*) '&'
            end if
          end if
          if (WhatSystem==2) then
            if (print_NucleonNonlocalIntegral == 1) then
              write(303,*) '#',L,Jp
              do n=1,NmaxLong
                write(303,*) n*StepSizeLong,real(constant*NLint(n)),aimag(constant*NLint(n))
              end do
              write(303,*) '&'
            end if
          end if

!       END NONLOCAL IF
        end if

!     END Jp IF
      end if
!   END Jp LOOP
    end do
! END L LOOP
  end do

  call diffCS(SM,Lmax,SpinProjectile,SpinTarget,eta,k,Jp_Max,WhatSystem,printing,Directory)

!  NmaxLong=NmaxReturn

  return
  end subroutine scattering_state
