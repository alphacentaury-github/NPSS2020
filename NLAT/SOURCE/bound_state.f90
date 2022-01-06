! CREATED BY: LUKE TITUS
! LAST UPDATED: 2/2/2016
! THIS PROGRAM SOLVES THE BOUND STATE FOR LOCAL AND NONLOCAL POTENTIALS

  subroutine bound_state(accuracy,BoundParameters,ReturnWF,ReturnEnergy, &
                         StepSize,Rmax,nmax,printing,Directory)
	
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
  integer::n,nmax,nmatch,start,ChargeCore,ChargeFragment,nANC,q,WhatPot
  integer::NonLoc,NumberNodes,ParityCore,ParityFragment,L,calc,npoints
  integer::WhatCalc,NodesCount,sign,printing(1:100)
  integer::print_LocalBoundWF,print_NonlocalBoundWF,WhatSystem
  real(8)::r,StepSize,time,E,diff,hbarc,norm,Rmax,Rmatch,k,eta,J,LocalEnergy
  real(8)::MassCore,MassFragment,ReducedMass,SchEqConst,SpinCore,SpinFragment,const,logdout,logdin 
  real(8)::Estart,Estep,EnergyBackup,Rprime
  real(8)::Ediff,Eold,LocalWF(nmax),ReturnWF(nmax),Vnp(nmax),ReturnEnergy
  real(8)::ancnew,ancold,locANC,anc,EdiffConvergence
  real(8)::accuracy(1:100),BoundParameters(1:15,1:100),diffsave,Esave
  real(8),allocatable::whittaker(:),vpot(:),source(:),wfout(:),wfin(:),wf(:)
  real(8),allocatable::Kernel(:,:),SourceIn(:),SourceOut(:),KernelRead(:,:)
  character(LEN=50)::Directory,extension
  character(LEN=100)::filename

  print_LocalBoundWF=printing(5)
  print_NonlocalBoundWF=printing(6)

  if (print_LocalBoundWF == 1) then
!    open(200,file='LocalBoundWF.txt')
    extension='LocalBoundWF.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(200,file=trim(filename))
  end if
  if (print_NonlocalBoundWF == 1) then
!    open(201,file='NonlocalBoundWF.txt') 
    extension='NonlocalBoundWF.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(201,file=trim(filename))
  end if

  WhatSystem = nint(BoundParameters(2,1))
  WhatPot=nint(BoundParameters(3,2))

! WhatPot = 1 = User defined nonlocal potential
! WhatPot = 2 = Read in nonlocal potential
  if (WhatPot == 2) then
    open(202,file='NLpotBound.txt')
  end if

! WhatSystem = 1 = Deuteron bound state
! WhatSystem = 2 = Nucleon bound state
  if (WhatSystem==1) then
    write(*,*) ''
    write(*,*) 'Deuteron Bound State'
    Rmatch=accuracy(4)
  end if
  if (WhatSystem==2) then
    write(*,*) ''
    write(*,*) 'Nucleon Bound State'
    Rmatch=accuracy(5)
  end if

  npoints=accuracy(2)
  Estart=accuracy(6)
  Estep = accuracy(7)
  EnergyBackup=accuracy(8)
  EdiffConvergence=accuracy(9)

  MassFragment=BoundParameters(1,1)
  MassCore=BoundParameters(1,2)
  ReducedMass=((MassCore*MassFragment)/(MassCore+MassFragment))*accuracy(1)
  ChargeFragment=BoundParameters(1,3)
  ChargeCore=BoundParameters(1,4)
  SpinFragment=BoundParameters(1,5)
  SpinCore=BoundParameters(1,6)
  ParityFragment=BoundParameters(1,7)
  ParityCore=BoundParameters(1,8)
  L=nint(BoundParameters(1,9))
  J=BoundParameters(1,10)
  NumberNodes=nint(BoundParameters(1,11))
  Nonloc=nint(BoundParameters(3,1))
  calc=0
  nmatch=nint(Rmatch/StepSize)
  hbarc=197.32705
  SchEqConst=(2*ReducedMass)/(hbarc**2.0)

!  write(*,*) 'Rmax=',Rmax
!  write(*,*) 'Rmatch=',Rmatch
!  write(*,*) 'ChargeCore=',ChargeCore
!  write(*,*) 'ChargeFragment=',ChargeFragment
!  write(*,*) 'MassCore=',MassCore
!  write(*,*) 'MassFragment',MassFragment
!  write(*,*) 'Unit Mass=',accuracy(3)
!  write(*,*) 'ReducedMass=',ReducedMass
!  write(*,*) 'SpinCore=',SpinCore
!  write(*,*) 'SpinFragment=',SpinFragment
!  write(*,*) 'ParityCore=',ParityCore
!  write(*,*) 'ParityFragment=',ParityFragment
!  write(*,*) 'L=',L
!  write(*,*) 'J=',J
!  write(*,*) 'NumberNodes=',NumberNodes
!  write(*,*) 'NonLocal=',NonLoc
!  write(*,*) 'nmax=',nmax
!  write(*,*) 'npoints=',npoints
!  write(*,*) 'SchEqConst=',SchEqConst
!  write(*,*) 'nmatch=',nmatch
!  write(*,*) 'WhatCalc = ',WhatCalc
!  write(*,*) 'StepSize = ', StepSize
!  write(*,*) 'Rmax = ', Rmax
!  write(*,*) 'nmax = ', nmax

  allocate(source(1:nmax))
  allocate(SourceIn(1:nmax))
  allocate(SourceOut(1:nmax))
  allocate(vpot(1:nmax))
  allocate(whittaker(1:nmax))
  allocate(wfout(1:nmax))
  allocate(wfin(1:nmax))
  allocate(wf(1:nmax))
  allocate(Kernel(1:nmax,1:npoints))
  allocate(KernelRead(1:nmax,1:nmax))

  call local_binding_potential(BoundParameters,vpot,nmax,StepSize) 

  source(:)=0.0
  start=1
  const=1.0

! SOLVE THE LOCAL BOUND STATE PROBLEM
  E=Estart
  diff=2.0
  diffsave=2.0
  do while (E<-(2*Estep))
    E=E+Estep
    k=sqrt(-SchEqConst*E)
    eta=(ChargeCore*ChargeFragment*1.43997*ReducedMass)/((hbarc**2)*k)
    call nm_out(L,nmax,StepSize,vpot,SchEqConst,E,logdout,source,wfout,nmatch)
    call usewitt(nmax,StepSize,L,whittaker,SchEqConst,start,k,eta)
    wfin(nmax)=whittaker(nmax)
    wfin(nmax-1)=whittaker(nmax-1)
    call nm_in(L,nmax,StepSize,vpot,SchEqConst,E,logdin,source,wfin,nmatch,const,eta)
    diff=abs((logdin-logdout)/logdin)

    if (diff<diffsave) then

      const=wfout(nmatch)/wfin(nmatch)
      do n=1,nmatch
        wfin(n)=wfout(n)/const
      end do
      do n=nmatch,nmax
        wfout(n)=const*wfin(n)
      end do

!     CHECK NUMBER OF NODES
      NodesCount=1
      sign=1
      do n=1,nmax
        if (sign==1) then
          if (wfout(n)<0) then
            sign=-1
            NodesCount=NodesCount+1
          end if
        end if
        if (sign==-1) then
          if (wfout(n)>0) then
            sign=1
            NodesCount=NodesCount+1
          end if
        end if            
      end do
          
      if (NodesCount==NumberNodes) then
        Esave=E
        diffsave=diff
      end if

    end if
  end do 
  ReturnEnergy=Esave
  write(*,*) 'Local Binding Energy =',ReturnEnergy

! CALCULATE WAVE FUNCTION AT CORRECT ENERGY
  E=Esave
  k=sqrt(-SchEqConst*E)
  eta=(ChargeCore*ChargeFragment*1.43997*ReducedMass)/((hbarc**2)*k)
  call nm_out(L,nmax,StepSize,vpot,SchEqConst,E,logdout,source,wfout,nmatch)
  call usewitt(nmax,StepSize,L,whittaker,SchEqConst,start,k,eta)
  wfin(nmax)=whittaker(nmax)
  wfin(nmax-1)=whittaker(nmax-1)
  call nm_in(L,nmax,StepSize,vpot,SchEqConst,E,logdin,source,wfin,nmatch,const,eta)

! MATCH WF FOR THE OUTWARD INTEGRATION
  const=wfout(nmatch)/wfin(nmatch)    
  do n=1,nmatch
    LocalWF(n)=wfout(n)
    wfin(n)=wfout(n)/const
  end do
  do n=nmatch,nmax
    LocalWF(n)=const*wfin(n)            
    wfout(n)=const*wfin(n)
  end do

! NORMALIZE THE LOCAL WAVE FUNCTION
  norm=0.0
  do n=1,nmax
    norm=norm+StepSize*abs(LocalWF(n))**2.0
  end do
  LocalWF(:)=(1/(sqrt(norm)))*LocalWF(:)    

! PRINT LOCAL WAVE FUNCTION TO FILE 
  if (print_LocalBoundWF == 1) then
    do n=1,nmax
      write(200,*) n*StepSize,LocalWF(n)
    end do
  end if

  if (Nonloc==0) then
    ReturnWF(:)=LocalWF(:)
  end if

! FIND ANC OF LOCAL WAVE FUNCTION
  anc=0.0
  diff=2.0
  n=5*nmatch
  ancnew=LocalWF(n)/whittaker(n)
  do while (diff >=1e-4 .and. n<(nmax-1))
     n=n+1
     ancold=ancnew
     ancnew=LocalWF(n)/whittaker(n)   
     diff=abs((ancnew-ancold)/ancnew)
  end do
  nANC=n
  locANC=ancnew 
  write(*,*) 'Local ANC=',abs(locANC)
 
! SOLVE NONLOCAL BOUND STATE PROBLEM 
  if (NonLoc==1) then

!   WhatPot = 1 = User defined nonlocal potential
!   WhatPot = 2 = Read in nonlocal potential
    if (WhatPot == 2) then
      do n=1,nmax
        do q=1,nmax
          read(202,*) R,Rprime,KernelRead(n,q)
        end do
      end do
    end if

    Ediff=2.0
    do while (Ediff>EdiffConvergence)
      call calc_binding_source(BoundParameters,vpot,nmax,L,wfin,wfout,SchEqConst, &
                               SourceIn,SourceOut,calc,npoints,Kernel,KernelRead, &
                               StepSize,Rmatch)

      calc=1
      E=Esave-EnergyBackup
      diff=2.0
      diffsave=10.0
      do while(E<-(2*Estep))
        E=E+Estep
        k=sqrt(-SchEqConst*E)
        eta=(ChargeCore*ChargeFragment*1.43997*ReducedMass)/((hbarc**2)*k)
        call nm_out(L,nmax,StepSize,vpot,SchEqConst,E,logdout,SourceOut,wfout,nmatch)
        call usewitt(nmax,StepSize,L,whittaker,SchEqConst,start,k,eta)
        wfin(nmax)=whittaker(nmax)
        wfin(nmax-1)=whittaker(nmax-1)
        call nm_in(L,nmax,StepSize,vpot,SchEqConst,E,logdin,SourceIn,wfin,nmatch,const,eta)
        diff=abs((logdin-logdout)/logdin)
        if (diff<diffsave) then

          const=wfout(nmatch)/wfin(nmatch)
          do n=nmatch,nmax
            wfout(n)=const*wfin(n)
          end do    
          
!         CHECK NUMBER OF NODES
          NodesCount=1
          sign=1
          do n=1,nmax
            if (sign==1) then
              if (wfout(n)<0) then
                sign=-1
                NodesCount=NodesCount+1
              end if
            end if
            if (sign==-1) then
              if (wfout(n)>0) then
                sign=1
                NodesCount=NodesCount+1
              end if
            end if            
          end do
          
          if (NodesCount==NumberNodes) then
            Esave=E
            diffsave=diff
          end if

        end if
      end do 
      write(*,*) 'E = ',Esave
      Ediff=abs((Eold-Esave)/Esave)
      Eold=Esave

!     CALCULATE WAVE FUNCTION AT CORRECT ENERGY
      E=Esave
      k=sqrt(-SchEqConst*E)
      eta=(ChargeCore*ChargeFragment*1.43997*ReducedMass)/((hbarc**2)*k)
      call nm_out(L,nmax,StepSize,vpot,SchEqConst,E,logdout,SourceOut,wfout,nmatch)
      call usewitt(nmax,StepSize,L,whittaker,SchEqConst,start,k,eta)
      wfin(nmax)=whittaker(nmax)
      wfin(nmax-1)=whittaker(nmax-1)
      call nm_in(L,nmax,StepSize,vpot,SchEqConst,E,logdin,SourceIn,wfin,nmatch,const,eta)

!     MATCH INWARD AND OUTWARD WAVE FUNCTIONS
      const=wfout(nmatch)/wfin(nmatch)    
      do n=1,nmatch
        wfin(n)=wfout(n)/const
      end do
      do n=nmatch,nmax
        wfout(n)=const*wfin(n)
      end do

!   END Ediff LOOP
    end do

    ReturnEnergy=Esave
    write(*,*) 'Nonlocal Binding Energy = ',ReturnEnergy

!   MATCH WF FOR THE OUTWARD INTEGRATION
    const=wfout(nmatch)/wfin(nmatch)    
    do n=1,nmatch
      ReturnWF(n)=wfout(n)
    end do
    do n=nmatch,nmax
      ReturnWF(n)=const*wfin(n)            
    end do

!   NORMALIZE THE NONLOCAL WAVE FUNCTION
    norm=0.0
    do n=1,nmax
      norm=norm+StepSize*abs(ReturnWF(n))**2.0
    end do

    ReturnWF(:)=(1/(sqrt(norm)))*ReturnWF(:)

!   FIND ANC OF NONLOCAL WAVE FUNCTION
    anc=0.0
    diff=2.0
    n=5*nmatch
    ancnew=ReturnWF(n)/whittaker(n)
    do while (diff >=1e-4 .and. n<(nmax-1))
       n=n+1
       ancold=ancnew
       ancnew=ReturnWF(n)/whittaker(n)   
       diff=abs((ancnew-ancold)/ancnew)
    end do
    nANC=n
    locANC=ancnew 
    write(*,*) 'Nonlocal ANC=',abs(locANC)

    if (print_NonlocalBoundWF == 1) then
      do n=1,nmax
        write(201,*) n*StepSize,ReturnWF(n)
      end do
    end if

! END NONLOCAL CALCULATION
  end if
    
  return   	
  END subroutine bound_state
