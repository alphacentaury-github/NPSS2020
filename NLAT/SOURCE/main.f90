  program main
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::n,L,Lmax,nmax,SH_index_max,SH_index,Jp_dA_max,Jp_Nucleon_max,WhatCalc
  integer::printing(1:100),FinalBoundL
  integer::print_DeuteronBoundWF,print_NucleonBoundWF
  integer::print_DeuteronScatWFs,print_NucleonScatWFs
  real(8)::Ap_d,At_d,Ap_N,At_n,Ecm_d,Ecm_N
  real(8)::common(1:100),accuracy(1:100)
  real(8)::DeuteronBoundParameters(1:15,1:100),NucleonBoundParameters(1:15,1:100)
  real(8)::DeuteronScatParameters(1:15,1:100),NucleonScatParameters(1:15,1:100)
  real(8)::SpinNucleon,SpinDeuteron,Qvalue,Rmax,ElabEntrance
  real(8)::time,DeuteronBoundEnergy,NucleonBoundEnergy,DeuteronSpin,mass,mass_d,eta,eta_d,StepSize
  real(8)::theta,m,pi,SH_step,Jp
  real(8),allocatable::DeuteronBoundWF(:),NucleonBoundWF(:),Vnp(:),SH(:,:),yr(:),yi(:)
  complex*16,allocatable::DeuteronScatWFs(:,:,:),NucleonScatWFs(:,:,:)
  complex*16::k,k_d
  character(LEN=50)::Directory,command,extension
  character(LEN=100)::filename
  character(LEN=6)::mkdir

  call front_end(common,accuracy,DeuteronBoundParameters,NucleonBoundParameters, &
                 DeuteronScatParameters,NucleonScatParameters,printing,Directory)

  mkdir='mkdir '
  command=mkdir // trim(Directory) 
  call system(command)

  print_DeuteronBoundWF=printing(1)
  print_NucleonBoundWF=printing(2)
  print_DeuteronScatWFs=printing(3)
  print_NucleonScatWFs=printing(4)

  if (print_DeuteronBoundWF == 1) then
    extension='DeuteronBoundWF.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(100,file=trim(filename))
  end if
  if (print_NucleonBoundWF == 1) then
!    open(101,file='NucleonBoundWF.txt')
    extension='NucleonBoundWF.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(101,file=trim(filename))
  end if
  if (print_DeuteronScatWFs == 1) then
!    open(102,file='DeuteronScatWFs.txt')
    extension='DeuteronScatWFs.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(102,file=trim(filename))
  end if
  if (print_NucleonScatWFs == 1) then
!    open(103,file='NucleonScatWFs.txt')
    extension='NucleonScatWFs.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(103,file=trim(filename))
  end if

! DEFINE THE VARIABLES COMMON TO ALL SUBROUTINES
  WhatCalc=nint(common(1))
  StepSize=common(2)  
  Rmax=common(3)
  nmax=nint(Rmax/StepSize)
  ElabEntrance=common(4)
  Lmax=nint(common(5))
  SH_step=accuracy(3)
  pi=acos(-1.0)
  SH_index_max=nint(pi/SH_step)

  Qvalue=common(6)
  Ap_d=DeuteronScatParameters(1,1)
  At_d=DeuteronScatParameters(1,2)
  Ap_N=NucleonScatParameters(1,1)
  At_N=NucleonScatParameters(1,2)
  FinalBoundL=nint(NucleonBoundParameters(1,9))

! WhatCalc = 1 = Deuteron Bound State
! WhatCalc = 2 = Nucleon Bound State
! WhatCalc = 3 = Deuteron Scattering State
! WhatCalc = 4 = Nucleon Scattering State
! WhatCalc = 5 = (d,N)
! WhatCalc = 6 = (N,d)
  if (WhatCalc == 3) then
    Ecm_d=((At_d)/(Ap_d+At_d))*ElabEntrance
  end if
  if (WhatCalc == 4) then
    Ecm_N=((At_N)/(Ap_N+At_N))*ElabEntrance
  end if
  if (WhatCalc == 5) then
    Ecm_d=((At_d)/(Ap_d+At_d))*ElabEntrance
    Ecm_N=Ecm_d+Qvalue
  end if
  if (WhatCalc == 6) then
    Ecm_N=((At_N)/(Ap_N+At_N))*ElabEntrance
    Ecm_d=Ecm_N+Qvalue
  end if

!  write(*,*) 'WhatCalc = ',WhatCalc
!  write(*,*) 'StepSize = ', StepSize
!  write(*,*) 'Rmax = ', Rmax
!  write(*,*) 'nmax = ', nmax
!  write(*,*) 'ElabEntrance = ', ElabEntrance
!  write(*,*) 'Lmax = ', Lmax
!  write(*,*) 'SH_step=', SH_step  
!  write(*,*) 'SH_index_max = ', SH_index_max
!  write(*,*) 'Qvalue = ', Qvalue

  allocate(Vnp(1:nmax))
  allocate(SH(0:Lmax,0:SH_index_max))
  allocate(yr(0:Lmax))
  allocate(yi(0:Lmax))
  allocate(DeuteronBoundWF(1:nmax))
  allocate(NucleonBoundWF(1:nmax))

  if (WhatCalc==3 .or. WhatCalc == 5 .or. WhatCalc == 6) then
    SpinDeuteron=DeuteronScatParameters(1,5)
    Jp_dA_max=nint(2*(Lmax+SpinDeuteron))
    allocate(DeuteronScatWFs(0:Lmax,0:Jp_dA_max,1:nmax))
  end if

  if (WhatCalc==4 .or. WhatCalc == 5 .or. WhatCalc == 6) then
    SpinNucleon=NucleonScatParameters(1,5)
    Jp_Nucleon_max=nint(2*(Lmax+SpinNucleon))
    allocate(NucleonScatWFs(0:Lmax,0:Jp_Nucleon_max,1:nmax))
  end if

! GET FROM SUBROUTINE
  do n=1,nmax
    Vnp(n)=-71.85*exp(-(n*StepSize/1.494)**2.0)
  end do 

! CONSTRUCT SPHERICAL HARMONICS
! SH DEFINED AS REAL SINCE M=0 !!!
  m=0
  do SH_index=1,SH_index_max
    theta=SH_index*SH_step
    call spherical_harmonic(Lmax,int(m),theta,0.0,yr,yi)  
    do L=0,Lmax
      SH(L,SH_index)=yr(L) 
    end do
  end do

! DEUTERON BOUND STATE
  if (WhatCalc==1 .or. WhatCalc ==3 .or. WhatCalc==5 .or. WhatCalc==6) then
    call bound_state(accuracy,DeuteronBoundParameters,DeuteronBoundWF,DeuteronBoundEnergy, &
                     StepSize,Rmax,nmax,printing,Directory)
    write(*,*) 'Deuteron Bound Energy=',DeuteronBoundEnergy
    if (print_DeuteronBoundWF == 1) then
      do n=1,nmax
        write(100,*) n*StepSize,DeuteronBoundWF(n)
      end do
    end if
  end if

! NUCLEON BOUND STATE
  if (WhatCalc==2 .or. WhatCalc==5 .or. WhatCalc==6) then
    call bound_state(accuracy,NucleonBoundParameters,NucleonBoundWF,NucleonBoundEnergy, &
                     StepSize,Rmax,nmax,printing,Directory)
    write(*,*) 'Nucleon Bound Energy = ',NucleonBoundEnergy  
    if (print_NucleonBoundWF == 1) then
      do n=1,nmax
        write(101,*) n*StepSize,NucleonBoundWF(n)     
      end do
    end if
  end if


! DEUTERON SCATTERING STATE
  if (WhatCalc==3 .or. WhatCalc==5 .or. WhatCalc==6) then

    call scattering_state(DeuteronScatParameters,common,accuracy,DeuteronScatWFs,Lmax,nmax,k_d,eta_d, &
                          DeuteronBoundWF,Vnp,SH,SH_index_max,Lmax,SH_step, &
                          Jp_dA_max,mass_d,Ecm_d,printing,Directory)

!   PRINT THE DEUTERON SCATTERING WAVE FUNCTIONS
    if (print_DeuteronScatWFs == 1) then
      do L=0,Lmax
        do Jp=(L-SpinDeuteron),(L+SpinDeuteron)
          if (L==0 .and. Jp<0.001) then
!           DO NOTHING
          else if (Jp>-0.001) then
            write(102,*) '#',L,Jp
            do n=1,nmax-1
              write(102,*) n*StepSize,real(DeuteronScatWFs(L,nint(2*Jp),n)),aimag(DeuteronScatWFs(L,nint(2*Jp),n))
            end do
            write(102,*) '&'
          end if
        end do
      end do
    end if  

  end if

! NUCLEON SCATTERING STATE
  if (WhatCalc==4 .or. WhatCalc==5 .or. WhatCalc==6) then
    call scattering_state(NucleonScatParameters,common,accuracy,NucleonScatWFs,Lmax,nmax,k,eta, &
                          DeuteronBoundWF,Vnp,SH,SH_index_max,Lmax,SH_step, &
                          Jp_Nucleon_max,mass,Ecm_N,printing,Directory)

!   PRINT THE NUCLEON SCATTERING WAVE FUNCTIONS
    if (print_NucleonScatWFs == 1) then
      do L=0,Lmax
        do Jp=L-SpinNucleon,L+SpinNucleon
          if (Jp>-0.001) then
            write(103,*) '#',L,Jp
            do n=1,nmax
              write(103,*) n*StepSize,real(NucleonScatWFs(L,nint(2*Jp),n)),aimag(NucleonScatWFs(L,nint(2*Jp),n))
            end do
            write(103,*) '&'
          end if
        end do
      end do
    end if

  end if

  if (WhatCalc==5 .or. WhatCalc==6) then
    call tmatrix (accuracy,common,DeuteronBoundParameters,NucleonBoundParameters, &
                  DeuteronScatParameters, &
                  NucleonScatParameters,DeuteronBoundWF,NucleonBoundWF,NucleonScatWFs,DeuteronScatWFs, &
                  k,k_d,mass,mass_d,Lmax,nmax,eta,eta_d,StepSize,Vnp,SH_index_max,SH_step, &
                  Jp_dA_max,Jp_Nucleon_max,WhatCalc,printing,FinalBoundL,Directory)
  end if

  call cpu_time(time)
  write (*,*) 'Total Run Time (sec) =',time

  end program main
