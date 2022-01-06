  subroutine front_end(common,accuracy,DeuteronBoundParameters,NucleonBoundParameters, &
                       DeuteronScatParameters,NucleonScatParameters,printing,Directory)
!  program front_end
  
  implicit none
!!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!! !!!!!!!!!!    
  integer::n,WhatAccuracy,WhatCalc,WhatPot,NonLoc,PotType,WhatPreDefPot
  integer::printing(1:100)
  real(8)::common(1:100),accuracy(1:100)
  real(8)::DeuteronBoundParameters(1:15,1:100),NucleonBoundParameters(1:15,1:100)
  real(8)::DeuteronScatParameters(1:15,1:100),NucleonScatParameters(1:15,1:100)
  character(LEN=50)::Directory

! INITIALIZE ALL ARRAYS
  common(:)=0.0
  DeuteronBoundParameters(:,:)=0.0
  NucleonBoundParameters(:,:)=0.0
  DeuteronScatParameters(:,:)=0.0
  NucleonScatParameters(:,:)=0.0
  accuracy(:)=0.0
  printing(:)=0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEFINE THE COMMON ARRAY 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1 = DEUTERON BOUND STATE
! 2 = N+A BOUND STATE
! 3 = DEUTERON SCATTERING STATE
! 4 = NUCLEON SCATTERING STATE
! 5 = (d,N) TRANSFER
! 6 = (N,d) TRANSFER
  read(*,*) common(1)
  WhatCalc=nint(common(1))

! STEP SIZE
  read(*,*) common(2)

! MAXIMUM RADIUS FOR CALCULATION
  read(*,*) common(3)

! IF SCATTERING STATES
  if (WhatCalc==3 .or. WhatCalc==4 .or. WhatCalc==5 .or. WhatCalc==6) then

!   TOTAL BEAM ENERGY OF INITIAL STATE IN THE LAB FRAME [MeV]
    read(*,*) common(4)

!   L-MAX
    read(*,*) common(5)

  end if
  
! IF TRANSFER REACTIONS
  if (WhatCalc==5 .or. WhatCalc==6) then

!   Q-VALUE FOR REACTION
    read(*,*) common(6)

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEUTERON BOUND STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! IF A DEUTERON BOUND STATE IS INVOLVED, DEFINTE WHAT THE SYSTEM IS
! ALWAYS USE DEFAULT VALUES FOR THE DEUTERON
  if (WhatCalc==1 .or. WhatCalc == 3 .or. WhatCalc==5 .or. WhatCalc==6) then

!   MASS OF FRAGMENT
    DeuteronBoundParameters(1,1)=1

!   MASS OF CORE
    DeuteronBoundParameters(1,2)=1

!   CHARGE FRAGMENT
    DeuteronBoundParameters(1,3)=0

!   CHARGE CORE
    DeuteronBoundParameters(1,4)=1

!   SPIN FRAGMENT
    DeuteronBoundParameters(1,5)=0.5

!   SPIN CORE
    DeuteronBoundParameters(1,6)=0.5

!   PARITY FRAGMENT
    DeuteronBoundParameters(1,7)=1

!   PARITY CORE
    DeuteronBoundParameters(1,8)=1

!   L OF BOUND STATE
    DeuteronBoundParameters(1,9)=0

!   J OF BOUND STATE
    DeuteronBoundParameters(1,10)=0.5

!   NUMBER OF NODES OF THE BOUND STATE
    DeuteronBoundParameters(1,11)=1

!   WHAT SYSTEM (1=n+p, 2=N+A)
    DeuteronBoundParameters(2,1)=1

!   LOCAL POTENTIAL TO USE
!   1 = CENTRAL GAUSSIAN
    read(*,*) DeuteronBoundParameters(2,2)

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NUCLEON BOUND STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! N+A BOUND STATE
  if (WhatCalc==2 .or. WhatCalc==5 .or. WhatCalc==6) then

!   MASS OF FRAGMENT AND CORE
    read(*,*) NucleonBoundParameters(1,1),NucleonBoundParameters(1,2)

!   CHARGE FRAGMENT AND CORE
    read(*,*) NucleonBoundParameters(1,3),NucleonBoundParameters(1,4)

!   SPIN FRAGMENT AND CORE
    read(*,*) NucleonBoundParameters(1,5),NucleonBoundParameters(1,6)

!   PARITY FRAGMENT AND CORE
    read(*,*) NucleonBoundParameters(1,7),NucleonBoundParameters(1,8)

!   L AND J OF BOUND STATE
    read(*,*) NucleonBoundParameters(1,9),NucleonBoundParameters(1,10)

!   NUMBER OF NODES OF THE BOUND STATE
    read(*,*) NucleonBoundParameters(1,11)

!   WHAT SYSTEM (1=n+p, 2=N+A)
    NucleonBoundParameters(2,1)=2

!   LOCAL POTENTIAL TO USE
    read(*,*) NucleonBoundParameters(2,2)
    WhatPot=nint(NucleonBoundParameters(2,2))

!   1 = USER DEFINED LOCAL POTENTIAL
    if (WhatPot==1) then
      
!     REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonBoundParameters(4,1), NucleonBoundParameters(4,2), NucleonBoundParameters(4,3)

!     REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonBoundParameters(4,4), NucleonBoundParameters(4,5), NucleonBoundParameters(4,6)

!     REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonBoundParameters(4,7), NucleonBoundParameters(4,8), NucleonBoundParameters(4,9)

!     COULOMB RADIUS
      read(*,*) NucleonBoundParameters(4,10)

    end if

!   INCLUDE NONLOCALITY? (0=NO, 1=YES)
    read(*,*) NucleonBoundParameters(3,1)
    NonLoc=nint(NucleonBoundParameters(3,1))

!   IF INCLUDING NONLOCALITY
    if (NonLoc==1) then

!     WHAT NONLOCAL POTENTIAL (1=USER DEFINED, 2=READ IN)
      read(*,*) NucleonBoundParameters(3,2)
      WhatPot=nint(NucleonBoundParameters(3,2))

!     IF USER DEFINED NONLOCAL POTENTIAL
      if (WhatPot==1) then

!       LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(5,1), NucleonBoundParameters(5,2), NucleonBoundParameters(5,3)

!       LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(5,4), NucleonBoundParameters(5,5), NucleonBoundParameters(5,6)

!       LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(5,7), NucleonBoundParameters(5,8), NucleonBoundParameters(5,9)

!       LOCAL PART OF NL: COULOMB RADIUS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(5,10)

!       NONLOCAL PART: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(6,1), NucleonBoundParameters(6,2), NucleonBoundParameters(6,3)

!       NONLOCAL PART: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(6,4), NucleonBoundParameters(6,5), NucleonBoundParameters(6,6)

!       NONLOCAL PART: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(6,7), NucleonBoundParameters(6,8), NucleonBoundParameters(6,9)

!       NONLOCAL PART: BETA: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(6,10)

!     END NONLOCAL PART OF POTENTIAL
      end if

!     IF USER DEFINED NONLOCAL POTENTIAL
      if (WhatPot==2) then

!       LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(5,1), NucleonBoundParameters(5,2), NucleonBoundParameters(5,3)

!       LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(5,4), NucleonBoundParameters(5,5), NucleonBoundParameters(5,6)

!       LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(5,7), NucleonBoundParameters(5,8), NucleonBoundParameters(5,9)

!       LOCAL PART OF NL: COULOMB RADIUS: N+A BOUND STATE
        read(*,*) NucleonBoundParameters(5,10)

!     END NONLOCAL PART OF POTENTIAL
      end if

!   END NONLOCALITY
    end if 

! END NUCLEON BOUND STATE
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEUTERON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (WhatCalc==3 .or. WhatCalc==5 .or. WhatCalc==6) then

!   MASS OF PROJECTILE AND TARGET
    read(*,*) DeuteronScatParameters(1,1),DeuteronScatParameters(1,2)

!   CHARGE OF PROJECTILE AND TARGET
    read(*,*) DeuteronScatParameters(1,3),DeuteronScatParameters(1,4)

!   SPIN OF PROJECTILE AND TARGET
    read(*,*) DeuteronScatParameters(1,5),DeuteronScatParameters(1,6)

!   PARITY OF PROJECTILE AND TARGET
    read(*,*) DeuteronScatParameters(1,7),DeuteronScatParameters(1,8)

!   WHAT SYSTEM (1=d+A, 2=N+A)
    DeuteronScatParameters(2,1)=1

!   LOCAL POTENTIAL TYPE (1=DEUTERON OPTICAL POTENTIAL, 2=NUCLEON POTENTIALS)
    read(*,*) DeuteronScatParameters(2,2)
    PotType=nint(DeuteronScatParameters(2,2))

!   PotType = 1 = Deuteron Optical Potential
!   Pottype = 2 = Nucleon Potentials
    if (PotType==1) then

!     LOCAL POTENTIAL TO USE
      read(*,*) DeuteronScatParameters(2,3)
      WhatPot=nint(DeuteronScatParameters(2,3))

!     WhatPot = 1 = User-defined deuteron local optical potential for scattering
!     WhatPot = 2 = Pre-defined deuteron local optical potential
      if (WhatPot==1) then

!       LOCAL REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(4,1), DeuteronScatParameters(4,2), DeuteronScatParameters(4,3)    

!       LOCAL IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(4,4), DeuteronScatParameters(4,5), DeuteronScatParameters(4,6)   

!       LOCAL REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(4,7), DeuteronScatParameters(4,8), DeuteronScatParameters(4,9)   

!       LOCAL IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(4,10), DeuteronScatParameters(4,11), DeuteronScatParameters(4,12)   

!       LOCAL REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(4,13), DeuteronScatParameters(4,14), DeuteronScatParameters(4,15)   

!       LOCAL IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(4,16), DeuteronScatParameters(4,17), DeuteronScatParameters(4,18)   

!       LOCAL COULOMB RADIUS
        read(*,*) DeuteronScatParameters(4,19)

!     END OF LOCAL POTENTIAL
      end if

!     WhatPot = 1 = User-defined deuteron local optical potential for scattering
!     WhatPot = 2 = Pre-defined deuteron local optical potential
      if (WhatPot==2) then
        read(*,*) DeuteronScatParameters(2,4)
      end if

!   END DEUTERON LOCAL OPTICAL POTENTIAL TYPE
    end if

!   PotType = 1 = Deuteron Optical Potential
!   PotType = 2 = Adiabatic Potentials
    if (PotType == 2) then

!     LOCAL POTENTIAL TO USE
      read(*,*) DeuteronScatParameters(2,5)
      WhatPot=nint(DeuteronScatParameters(2,5))

!     WhatPot = 1 = User defined nucleon potentials for the deuteron
!     WhatPot = 2 = Pre-defined nucleon potentials for the deuteron
      if (WhatPot==1) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       NEUTRON POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       LOCAL REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(7,1), DeuteronScatParameters(7,2), DeuteronScatParameters(7,3)    

!       LOCAL IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(7,4), DeuteronScatParameters(7,5), DeuteronScatParameters(7,6)   

!       LOCAL REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(7,7), DeuteronScatParameters(7,8), DeuteronScatParameters(7,9)   

!       LOCAL IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(7,10), DeuteronScatParameters(7,11), DeuteronScatParameters(7,12)   

!       LOCAL REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(7,13), DeuteronScatParameters(7,14), DeuteronScatParameters(7,15)   

!       LOCAL IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(7,16), DeuteronScatParameters(7,17), DeuteronScatParameters(7,18) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       PROTON POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       LOCAL REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(8,1), DeuteronScatParameters(8,2), DeuteronScatParameters(8,3)   

!       LOCAL IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(8,4), DeuteronScatParameters(8,5), DeuteronScatParameters(8,6)   

!       LOCAL REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(8,7), DeuteronScatParameters(8,8), DeuteronScatParameters(8,9)   

!       LOCAL IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(8,10), DeuteronScatParameters(8,11), DeuteronScatParameters(8,12)   

!       LOCAL REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(8,13), DeuteronScatParameters(8,14), DeuteronScatParameters(8,15)   

!       LOCAL IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) DeuteronScatParameters(8,16), DeuteronScatParameters(8,17), DeuteronScatParameters(8,18)   

!       LOCAL COULOMB RADIUS
        read(*,*) DeuteronScatParameters(8,19)

!     END USER DEFINED LOCAL NUCLEON POTENTIALS FOR THE DEUTERON
      end if      
   
!     WhatPot = 1 = User defined nucleon potentials in the adiabatic potential
!     WhatPot = 2 = Pre-defined nucleon potentials in the adiabatic potential
      if (WhatPot==2) then
      
!        WHAT PRE-DEFINED LOCAL NUCLEON POTENTIALS FOR THE DEUTERON (1=KD, 2=CH89)
         read(*,*) DeuteronScatParameters(2,6)
   
      end if

!   END NUCLEON POTENTIALS FOR THE ADIABATIC POTENTIAL
    end if

!   INCLUDE NONLOCALITY? (0=NO, 1=YES)
    read(*,*) DeuteronScatParameters(3,1)
    NonLoc=nint(DeuteronScatParameters(3,1))

!   IF NONLOCALITY IS INCLUDED
    if (NonLoc == 1) then

!     POTENTIAL TYPE? (1=TWO-BODY SCATTERING, 2=ADIABATIC POTENTIAL)
      read(*,*) DeuteronScatParameters(3,2)
      PotType=nint(DeuteronScatParameters(3,2))

!     PotType = 1 = Deuteron nonlocal optical potential
!     Pottype = 2 = Nucleon potentials for the adiabatic potential
      if (PotType ==1) then

!       WHAT NONLOCAL POTENTIAL (1=User defined, 2=Pre-defined, 3=Read in)
        read(*,*) DeuteronScatParameters(3,3)
        WhatPot=nint(DeuteronScatParameters(3,3))

!       WhatPot = 1 = User defined deuteron nonlocal optical potential
!       WhatPot = 2 = Pre-defined deuteron nonlocal potential
!       WhatPot = 3 = Read in
        if (WhatPot == 1) then

!         LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(5,1), DeuteronScatParameters(5,2), DeuteronScatParameters(5,3)    

!         LOCAL PART OF NL: IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(5,4), DeuteronScatParameters(5,5), DeuteronScatParameters(5,6)   

!         LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(5,7), DeuteronScatParameters(5,8), DeuteronScatParameters(5,9)   

!         LOCAL PART OF NL: IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(5,10), DeuteronScatParameters(5,11), DeuteronScatParameters(5,12)   

!         LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(5,13), DeuteronScatParameters(5,14), DeuteronScatParameters(5,15)   

!         LOCAL PART OF NL: IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(5,16), DeuteronScatParameters(5,17), DeuteronScatParameters(5,18)   

!         LOCAL PART OF NL: COULOMB RADIUS
          read(*,*) DeuteronScatParameters(5,19)

!         NONLOCAL REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(6,1), DeuteronScatParameters(6,2), DeuteronScatParameters(6,3)    

!         NONLOCAL IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(6,4), DeuteronScatParameters(6,5), DeuteronScatParameters(6,6)   

!         NONLOCAL REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(6,7), DeuteronScatParameters(6,8), DeuteronScatParameters(6,9)   

!         NONLOCAL IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(6,10), DeuteronScatParameters(6,11), DeuteronScatParameters(6,12)   

!         NONLOCAL REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(6,13), DeuteronScatParameters(6,14), DeuteronScatParameters(6,15)   

!         NONLOCAL IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(6,16), DeuteronScatParameters(6,17), DeuteronScatParameters(6,18)   

!         BETA: RANGE OF NONLOCALITY
          read(*,*) DeuteronScatParameters(6,19)

!       END USER DEFINED DEUTERON NONLOCAL POTENTIAL
        end if

!     END DEUTERON NONLOCAL OPTICAL POTENTIAL
      end if

!     PotType = 1 = Deuteron optical potential
!     PotType = 2 = Adiabatic potential
      if (PotType == 2 ) then

!       WHAT NONLOCAL POTENTIAL
        read(*,*) DeuteronScatParameters(3,3)
        WhatPot=nint(DeuteronScatParameters(3,3))
!        write(*,*) 'WhatPot = ', WhatPot

!       WhatPot = 1 = User defined nonlocal potential
!       WhatPot = 2 = Pre-defined nonlocal potential
!       WhatPot = 3 = Read in
        if (WhatPot == 1 ) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         NEUTRON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(9,1), DeuteronScatParameters(9,2), DeuteronScatParameters(8,3)    

!         LOCAL PART OF NL: IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(9,4), DeuteronScatParameters(9,5), DeuteronScatParameters(9,6)   

!         LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(9,7), DeuteronScatParameters(9,8), DeuteronScatParameters(9,9)   

!         LOCAL PART OF NL: IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(9,10), DeuteronScatParameters(9,11), DeuteronScatParameters(9,12)   

!         LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(9,13), DeuteronScatParameters(9,14), DeuteronScatParameters(9,15)   

!         LOCAL PART OF NL: IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(9,16), DeuteronScatParameters(9,17), DeuteronScatParameters(9,18) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         PROTON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(10,1), DeuteronScatParameters(10,2), DeuteronScatParameters(10,3)    

!         LOCAL PART OF NL: IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(10,4), DeuteronScatParameters(10,5), DeuteronScatParameters(10,6)   

!         LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(10,7), DeuteronScatParameters(10,8), DeuteronScatParameters(10,9)   

!         LOCAL PART OF NL: IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(10,10), DeuteronScatParameters(10,11), DeuteronScatParameters(10,12)   

!         LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(10,13), DeuteronScatParameters(10,14), DeuteronScatParameters(10,15)   

!         LOCAL PART OF NL: IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(10,16), DeuteronScatParameters(10,17), DeuteronScatParameters(10,18)   

!         LOCAL PART OF NL: COULOMB RADIUS
          read(*,*) DeuteronScatParameters(10,19)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         NEUTRON POTENTIAL: NONLOCAL POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(11,1), DeuteronScatParameters(11,2), DeuteronScatParameters(11,3)    

!         LOCAL PART OF NL: IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(11,4), DeuteronScatParameters(11,5), DeuteronScatParameters(11,6)   

!         LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(11,7), DeuteronScatParameters(11,8), DeuteronScatParameters(11,9)   

!         LOCAL PART OF NL: IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(11,10), DeuteronScatParameters(11,11), DeuteronScatParameters(11,12)   

!         LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(11,13), DeuteronScatParameters(11,14), DeuteronScatParameters(11,15)   

!         LOCAL PART OF NL: IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(11,16), DeuteronScatParameters(11,17), DeuteronScatParameters(11,18)

!         BETA: RANGE OF NONLOCALITY FOR THE NEUTRON
          read(*,*) DeuteronScatParameters(11,19)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         PROTON POTENTIAL: NONLOCAL POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(12,1), DeuteronScatParameters(12,2), DeuteronScatParameters(12,3)    

!         LOCAL PART OF NL: IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(12,4), DeuteronScatParameters(12,5), DeuteronScatParameters(12,6)   

!         LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(12,7), DeuteronScatParameters(12,8), DeuteronScatParameters(12,9)   

!         LOCAL PART OF NL: IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(12,10), DeuteronScatParameters(12,11), DeuteronScatParameters(12,12)   

!         LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(12,13), DeuteronScatParameters(12,14), DeuteronScatParameters(12,15)   

!         LOCAL PART OF NL: IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
          read(*,*) DeuteronScatParameters(12,16), DeuteronScatParameters(12,17), DeuteronScatParameters(12,18)   

!         BETA: RANGE OF NONLOCALITY FOR THE PROTON
          read(*,*) DeuteronScatParameters(12,19)

!       END USER DEFINED NUCLEON NONLOCAL POTENTIAL FOR DEUTERON
        end if

!       WhatPot = 1 = User defined nonlocal potential
!       WhatPot = 2 = Pre-defined nonlocal potential
!       WhatPot = 3 = Read in
        if (WhatPot==2) then

!         WHAT PRE-DEFINED NONLOCAL POTENTIAL (1=PB, 2=TPM)
          read(*,*) DeuteronScatParameters(3,4)

        end if

!     END NUCLEON OPTICAL POTENTIALS FOR THE DEUTERON
      end if

!   END OF NONLOCAL POTENTIAL
    end if

! END OF DEUTERON SCATTERING STATE
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NUCLEON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (WhatCalc==4 .or. WhatCalc==5 .or. WhatCalc==6) then

!   MASS OF PROJECTILE AND TARGET
    read(*,*) NucleonScatParameters(1,1),NucleonScatParameters(1,2)

!   CHARGE OF PROJECTILE AND TARGET
    read(*,*) NucleonScatParameters(1,3),NucleonScatParameters(1,4)

!   SPIN OF PROJECTILE AND TARGET
    read(*,*) NucleonScatParameters(1,5),NucleonScatParameters(1,6)

!   PARITY OF PROJECTILE AND TARGET
    read(*,*) NucleonScatParameters(1,7),NucleonScatParameters(1,8)

!   WHAT SYSTEM (1=d+A, 2=N+A)
    NucleonScatParameters(2,1)=2

!   LOCAL POTENTIAL TO USE
    read(*,*) NucleonScatParameters(2,2)
    WhatPot=nint(NucleonScatParameters(2,2))
!    write(*,*) 'WhatPot=',WhatPot

!   WhatPot = 1 = User defined local potential for N+A scattering
!   WhatPot = 2 = Pre-defined local potential for N+A scattering
    if (WhatPot==1) then

!     LOCAL REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonScatParameters(4,1), NucleonScatParameters(4,2), NucleonScatParameters(4,3)    

!     LOCAL IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonScatParameters(4,4), NucleonScatParameters(4,5), NucleonScatParameters(4,6)   

!     LOCAL REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonScatParameters(4,7), NucleonScatParameters(4,8), NucleonScatParameters(4,9)   

!     LOCAL IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonScatParameters(4,10), NucleonScatParameters(4,11), NucleonScatParameters(4,12)   

!     LOCAL REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonScatParameters(4,13), NucleonScatParameters(4,14), NucleonScatParameters(4,15)   

!     LOCAL IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
      read(*,*) NucleonScatParameters(4,16), NucleonScatParameters(4,17), NucleonScatParameters(4,18)   

!     LOCAL COULOMB RADIUS
      read(*,*) NucleonScatParameters(4,19)

!   END OF LOCAL POTENTIAL
    end if

!   WhatPot = 1 = User defined local potential for N+A scattering
!   WhatPot = 2 = Pre-defined local potential for N+A scattering
    if (WhatPot==2) then
      read(*,*) NucleonScatParameters(2,3)
    end if

!   INCLUDE NONLOCALITY? (0=NO, 1=YES)
    read(*,*) NucleonScatParameters(3,1)
    NonLoc=nint(NucleonScatParameters(3,1))
!    write(*,*) 'NonLoc = ', Nonloc

!   POTENTIAL TYPE
!   PotType = 1 = Two-body scattering state
    NucleonScatParameters(3,2)=1

!   IF INCLUDING NONLOCALITY
    if (NonLoc==1) then
      
!     WHAT NONLOCAL POTENTIAL (1=USER DEFINED, 2=PRE-DEFINED, 3=READ IN)
      read(*,*) NucleonScatParameters(3,3)
      WhatPot=nint(NucleonScatParameters(3,3))

!     WhatPot = 1 = User defined nonlocal potential
!     WhatPot = 2 = Pre-defined nonlocal potential 
!     WhatPot = 3 = Read in
      if (WhatPot == 1) then

!       LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,1), NucleonScatParameters(5,2), NucleonScatParameters(5,3)    

!       LOCAL PART OF NL: IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,4), NucleonScatParameters(5,5), NucleonScatParameters(5,6)     

!       LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,7), NucleonScatParameters(5,8), NucleonScatParameters(5,9)   

!       LOCAL PART OF NL: IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,10), NucleonScatParameters(5,11), NucleonScatParameters(5,12)   

!       LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,13), NucleonScatParameters(5,14), NucleonScatParameters(5,15)   

!       LOCAL PART OF NL: IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,16), NucleonScatParameters(5,17), NucleonScatParameters(5,18)   

!       LOCAL PART OF NL: COULOMB RADIUS
        read(*,*) NucleonScatParameters(5,19)

!       NONLOCAL REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(6,1), NucleonScatParameters(6,2), NucleonScatParameters(6,3)    

!       NONLOCAL IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(6,4), NucleonScatParameters(6,5), NucleonScatParameters(6,6)     

!       NONLOCAL REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(6,7), NucleonScatParameters(6,8), NucleonScatParameters(6,9)   

!       NONLOCAL IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(6,10), NucleonScatParameters(6,11), NucleonScatParameters(6,12)   

!       NONLOCAL REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(6,13), NucleonScatParameters(6,14), NucleonScatParameters(6,15)   
 
!       NONLOCAL IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(6,16), NucleonScatParameters(6,17), NucleonScatParameters(6,18)   
 
!       BETA: RANGE OF NONLOCALITY
        read(*,*) NucleonScatParameters(6,19)

!     END USER DEFINED NONLOCAL POTENTIAL
      end if

!     WhatPot = 1 = User defined nonlocal potential
!     WhatPot = 2 = Pre-defined nonlocal potential 
!     WhatPot = 3 = Read in
      if (WhatPot == 2) then

!       WHAT PRE-DEFINED NONLOCAL POTENTIAL (1=Perey-Buck, 2=TPM)
        read(*,*) NucleonScatParameters(3,4)     

!     END PRE-DEFINED NONLOCAL POTENTIAL
      end if

!     WhatPot = 1 = User defined nonlocal potential
!     WhatPot = 2 = Pre-defined nonlocal potential 
!     WhatPot = 3 = Read in
      if (WhatPot == 3) then

!       LOCAL PART OF NL: REAL VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,1), NucleonScatParameters(5,2), NucleonScatParameters(5,3)    

!       LOCAL PART OF NL: IMAG VOLUME DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,4), NucleonScatParameters(5,5), NucleonScatParameters(5,6)     

!       LOCAL PART OF NL: REAL SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,7), NucleonScatParameters(5,8), NucleonScatParameters(5,9)   

!       LOCAL PART OF NL: IMAG SURFACE DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,10), NucleonScatParameters(5,11), NucleonScatParameters(5,12)   

!       LOCAL PART OF NL: REAL SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,13), NucleonScatParameters(5,14), NucleonScatParameters(5,15)   

!       LOCAL PART OF NL: IMAG SPIN-ORBIT DEPTH, RADIUS, AND DIFFUSENESS
        read(*,*) NucleonScatParameters(5,16), NucleonScatParameters(5,17), NucleonScatParameters(5,18)   

!       LOCAL PART OF NL: COULOMB RADIUS
        read(*,*) NucleonScatParameters(5,19)   

      end if

!   END OF NONLOCALITY
    end if

! END OF NUCLEON SCATTERING STATE
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEFINE ACCURACY ARRAY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Mass unit
  read(*,*) accuracy(1)

! Number of Mesh Points for Gaussian NL Integration
  read(*,*) accuracy(2)

! Step in theta when calculating the spherical harmonics
  read(*,*) accuracy(3)

! Matching Radius for for n+p Bound State
  read(*,*) accuracy(4)

! Matching radius for nucleon bound state
  read(*,*) accuracy(5)

! Energy to start searching for bound state
  read(*,*) accuracy(6)

! Energy step when scanning for bound state
  read(*,*) accuracy(7)

! Energy to back up after each iteration in nonlocal bound state
  read(*,*) accuracy(8)
 
! Percent diff in energy at convergence of nonlocal bound state
  read(*,*) accuracy(9)

! % Diff of Log Deriv. at Convergence for nonlocal scattering state
  read(*,*) accuracy(10)

! Step size for the source in d+A nonlocal integration
  read(*,*) accuracy(11)

! Maximum radius to calculate d+A nonlocal source 
  read(*,*) accuracy(12)

! Maximum radius of dr integral in nonlocal adiabatic source
  read(*,*) accuracy(13)

! Maximum radius of ds integral in nonlocal adiabatic source
  read(*,*) accuracy(14)

! Number Of Mesh Points for d+A Source, dr Integral
  read(*,*) accuracy(15)

! Number Of Mesh Points for d+A Source, theta_r Integral
  read(*,*) accuracy(16)

! Number Of Mesh Points for d+A Source, ds Integral
  read(*,*) accuracy(17)

! Number Of Mesh Points for d+A Source, theta_s Integral 
  read(*,*) accuracy(18)

! Number Of Mesh Points for d+A Source, phi_s Integral
  read(*,*) accuracy(19)

! Maximum value of dR integral in calculation of T-matrix
  read(*,*) accuracy(20)

! Maximum value of dr integral in calculation of T-matrix
  read(*,*) accuracy(21)

! Number Of Mesh Points for T-Matrix Angular Integral 
  read(*,*) accuracy(22)

! Number Of Mesh Points for T-Matrix dR Radial Integral
  read(*,*) accuracy(23)

! Number Of Mesh Points for T-Matrix dr Radial Integral
  read(*,*) accuracy(24)

! Angular step when calculating transfer cross section
  read(*,*) accuracy(25)

! Angular step when calculating elastic cross section
  read(*,*) accuracy(26)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEFINE PRINTING ARRAY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Name of directory to print output files
  read(*,*) Directory
!  write(*,*) 'Directory = ', Directory

! DeuteronBoundWF.txt
  read(*,*) printing(1)

! NucleonBoundWF.txt
  read(*,*) printing(2)

! DeuteronScatWFs.txt
  read(*,*) printing(3)

! NucleonScatWFs.txt
  read(*,*) printing(4)

! LocalBoundWF.txt
  read(*,*) printing(5)

! NonlocalBoundWF.txt
  read(*,*) printing(6)

! DeuteronLocalIntegral.txt
  read(*,*) printing(7)

! NucleonLocalIntegral.txt
  read(*,*) printing(8)

! DeuteronNonlocalIntegral.txt
  read(*,*) printing(9)

! NucleonNonlocalIntegral.txt
  read(*,*) printing(10)

! DeuteronLocalSmatrix.txt
  read(*,*) printing(11)

! NucleonLocalSmatrix.txt
  read(*,*) printing(12)

! DeuteronNonlocalSmatrix.txt
  read(*,*) printing(13)

! NucleonNonlocalSmatrix.txt
  read(*,*) printing(14)

! DeuteronRatioToRuth.txt
  read(*,*) printing(15)

! NucleonRatioToRuth.txt
  read(*,*) printing(16)

! DeuteornElasticCS.txt
  read(*,*) printing(17)

! NucleonElasticCS.txt
  read(*,*) printing(18)

! TransferCS.txt
  read(*,*) printing(19)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RETURN
  end subroutine front_end
