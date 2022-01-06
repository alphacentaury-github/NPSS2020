  program make_input
  implicit none

  integer::What_Calc,What_Accuracy,What_QV,WhatPot,IOSTAT,EndOfFile,Nonloc,PotType
  integer::What_Print
  real(8)::Common(1:100), Accuracy(1:100), System(1:100), Printing(1:100)
  real(8)::Item1, Item2, Item3, Item4, Item5
  character(len=200)::Char1,Char2,Char3,Char4,Char5,Char6
  character(len=200)::Val1,Val2,Val3,Val4,Val5,Val6
  character(len=200)::Format1,Format2,Format3,Format4

  Format1='(A24,A)'
  Format2='(A6,A18,A)'
  Format3='(A9,A7,A8,A)'

  open(1,file='inputfile.in')
  open(2,file='temp.in')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHOOSE THE CALCULATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*) 'What kind of calculation?'
  write(*,*) '[1] Deuteron Bound State'
  write(*,*) '[2] Nucleon Bound State'
  write(*,*) '[3] Deuteron Scattering State'
  write(*,*) '[4] Nucleon Scattering State'
  write(*,*) '[5] (d,N) Transfer'
  write(*,*) '[6] (N,d) Transfer'

  read(2,*,IOSTAT=EndOfFile) What_Calc
  if (EndOfFile /= 0) then
    read(*,*) What_Calc
  end if
  write (val1,'(I3)') What_Calc
  char1='What Calculation? 1=n+p (B), 2=N+A (B), 3=d+A (S), 4=N+A (S), 5=(d,N), 6=(N,d)'
  write(1,Format1) adjustl(val1), adjustl(trim(char1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEFINE THE COMMON ARRAY 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(1,*) ''

  write(*,*) ''  
  write(*,*) 'Enter Step Size [fm]'    
  read(2,*,IOSTAT=EndOfFile) Item1
  if (EndOfFile /= 0) then 
    read(*,*) Item1
  end if
  Char1='Step Size'
  write (Val1,'(F7.3)') Item1
  write(1,Format1) adjustl(val1), adjustl(trim(char1))

  write(*,*) ''
  write(*,*) 'Enter Maximum Radius [fm]' 
  read(2,*,IOSTAT=EndOfFile) Item1
  if (EndOfFile /= 0) then 
    read(*,*) Item1
  end if
  write (Char1, '(A)')  'Maximum Radius'
  write (Val1,'(F7.3)') Item1
  write(1,Format1) adjustl(val1), adjustl(trim(char1))


  if (What_Calc == 3 .or. What_Calc == 4 .or. What_Calc == 5 .or. What_Calc == 6) then
    write(*,*) ''
    write(*,*) 'Total Beam Energy of Initial State in the Lab Frame [MeV]'
    read(2,*,IOSTAT=EndOfFile) Item1
    if (EndOfFile /= 0) then 
      read(*,*) Item1
    end if
    write (Char1, '(A)') 'Elab'
    write (Val1,'(F7.3)') Item1
    write(1,Format1) adjustl(val1), adjustl(trim(char1))

    write(*,*) ''
    write(*,*) 'L-Max'
    read(2,*,IOSTAT=EndOfFile) Item1
    if (EndOfFile /= 0) then 
      read(*,*) Item1
    end if
    write (Char1, '(A)')  'Lmax' 
    write (Val1,'(I4)') nint(Item1)
    write(1,Format1) adjustl(val1), adjustl(trim(char1))
  end if

  if (What_Calc == 5 .or. What_Calc == 6) then

    write(*,*) ''
    write(*,*) 'What is the Q-Value of the reaction?' 
    read(2,*,IOSTAT=EndOfFile) Item1
    if (EndOfFile /= 0) then 
      read(*,*) Item1
    end if
    write (Char1, '(A)') 'Q-Value'
    write (Val1,'(F7.3)') Item1
    write(1,Format1) adjustl(val1), adjustl(trim(char1))
    
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEUTERON BOUND STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (What_Calc==1 .or. What_Calc == 3 .or. What_Calc==5 .or. What_Calc==6) then

    write(1,*) ''

    write (*,*) ''
    write(*,*) 'What Pre-Defined Deuteron Binding Potential?'
    write(*,*) '[1]=Gaussian'
    read(2,*,IOSTAT=EndOfFile) WhatPot
    if (EndOfFile /= 0) then
      read(*,*) WhatPot
    end if
    write (Char1, '(A)')  'Pre-Defined Deuteron Binding Potential. [1]=Gaussian'
    write (Val1,'(I4)') WhatPot
    write(1,Format1) adjustl(val1), adjustl(trim(char1))        
  
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NUCLEON BOUND STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (What_Calc==2 .or. What_Calc==5 .or. What_Calc==6) then

    write(1,*) ''

    write (*,*) ''
    write(*,*) 'NUCLEON BOUND STATE'
    write(*,*) 'Mass number of fragment and core'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if      
    write (Char1, '(A)')  'Mass number of fragment and core (N+A BOUND STATE)'
    write (Val1,'(I4)') nint(Item1)
    write (Val2,'(I4)') nint(Item2)
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Charge of fragment and core'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Charge of fragment and core'
    write (Val1,'(I4)') nint(Item1)
    write (Val2,'(I4)') nint(Item2)
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Spin of fragment and core'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Spin of fragment and core'
    write (Val1,'(F7.1)') Item1
    write (Val2,'(F7.1)') Item2
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Parity of fragment and core'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Parity of fragment and core'
    write (Val1,'(F7.1)') Item1
    write (Val2,'(F7.1)') Item2
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'L and J of Bound State'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if
    write (Char1, '(A)')  'L and J of Bound State'
    write (Val1,'(I4)') nint(Item1)
    write (Val2,'(F7.1)') Item2
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Number of nodes in bound state? (Minimum of 1)'
    read(2,*,IOSTAT=EndOfFile) Item1
    if (EndOfFile /= 0) then 
      read(*,*) Item1
    end if
    write (Char1, '(A)')  'Number of nodes in bound state?'
    write (Val1,'(I4)') nint(Item1)
    write(1,Format1) adjustl(val1), adjustl(trim(char1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   LOCAL POTENTIAL: N+A BOUND STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*,*) ''
    write(*,*) 'What Local Potential To Use? (1=User defined)'
    read(2,*,IOSTAT=EndOfFile) WhatPot
    if (EndOfFile /= 0) then
      read(*,*) WhatPot
    end if
    write (Char1, '(A)')  'What Local Potential To Use? (1=User defined)'
    write (Val1,'(I4)') WhatPot
    write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!   WhatPot = 1 = User defined potential
    if (WhatPot == 1) then  

      write (*,*) ''
      write(*,*) 'Local Real Volume Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Real Volume Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))      

      write (*,*) ''
      write(*,*) 'Local Real Surface Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Real Surface Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

      write (*,*) ''
      write(*,*) 'Local Real Spin-Orbit Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Real Spin-Orbit Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

      write (*,*) ''
      write(*,*) 'Local Coulomb Radius'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Local Coulomb Radius'
      write (Val1,'(F7.3)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!   END USER DEFINED LOCAL POTENTIAL: N+A BOUND STATE
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NONLOCAL POTENTIAL: NUCLEON BOUND STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*,*) ''
    write(*,*) 'Use Nonlocality? [0]=No, [1]=Yes'
    read(2,*,IOSTAT=EndOfFile) Nonloc
    if (EndOfFile /= 0) then 
      read(*,*) Nonloc
    end if
    write (Char1, '(A)')  'Use Nonlocality? [0]=No, [1]=Yes'
    write (Val1,'(I4)') Nonloc
    write(1,Format1) adjustl(val1), adjustl(trim(char1))

    if (Nonloc==1) then

      write (*,*) ''
      write(*,*) 'What Nonlocal Potential? (1=User defined, 2=Read in)'
      read(2,*,IOSTAT=EndOfFile) WhatPot
      if (EndOfFile /= 0) then 
        read(*,*) WhatPot
      end if
      write (Char1, '(A)')  'What Nonlocal Potential? (1=User defined, 2=Read in)'
      write (Val1,'(I4)') WhatPot
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     WhatPot = 1 = User defined nonlocal binding potential
!     WhatPot = 2 = Read in
      if (WhatPot==1) then

        write (*,*) ''
        write(*,*) 'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))   

        write (*,*) ''
        write(*,*) 'Local Part of NL Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Part of NL: Coulomb Radius'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'Local Part of NL: Coulomb Radius'
        write (Val1,'(F7.3)') Item1
        write(1,Format1) adjustl(val1), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Nonlocal Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if  
        write (Char1, '(A)')  'Nonlocal Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Nonlocal Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Nonlocal Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Nonlocal Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Nonlocal Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Beta: Range of Nonlocality'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'Beta: Range of Nonlocality'
        write (Val1,'(F7.3)') Item1
        write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     END USER DEFINED NONLOCAL POTENTIAL
      end if

!     WhatPot = 1 = User defined nonlocal binding potential
!     WhatPot = 2 = Read in
      if (WhatPot==2) then

        write (*,*) ''
        write(*,*) 'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))   

        write (*,*) ''
        write(*,*) 'Local Part of NL Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Part of NL: Coulomb Radius'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'Local Part of NL: Coulomb Radius'
        write (Val1,'(F7.3)') Item1
        write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!     END READ IN NONLOCAL POTENTIAL
      end if

!   END NONLOCALITY
    end if

! END OF N+A BOUND STATE
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEUTERON SCATTERING STATE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (What_Calc==3 .or. What_Calc==5 .or. What_Calc==6) then

    write(1,*) ''

    write (*,*) ''
    write(*,*) 'DEUTERON SCATTERING STATE'
    write(*,*) 'Mass number of projectile and target'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if      
    write (Char1, '(A)')  'Mass number of projectile and target (DEUTERON SCATTERING STATE)'
    write (Val1,'(I4)') nint(Item1)
    write (Val2,'(I4)') nint(Item2)
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Charge of projectile and target'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Charge of projectile and target'
    write (Val1,'(I4)') nint(Item1)
    write (Val2,'(I4)') nint(Item2)
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Spin of projectile and target'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Spin of projectile and target'
    write (Val1,'(F7.1)') Item1
    write (Val2,'(F7.1)') Item2
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Parity of projectile and target'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Parity of projectile and target'
    write (Val1,'(F7.1)') Item1
    write (Val2,'(F7.1)') Item2
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   LOCAL POTENTIAL: DEUTERON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*,*) ''
    write(*,*) 'What Local Potential? (1=Deuteron Optical Potential, 2=Nucleon Potentials)'
    read(2,*,IOSTAT=EndOfFile) PotType
    if (EndOfFile /= 0) then
      read(*,*) PotType
    end if
    write (Char1, '(A)')  'What Local Potential? (1=Deuteron Optical Potential, 2=Nucleon Potentials)'
    write (Val1,'(I4)') PotType
    write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   LOCAL POTENTIAL: DEUTERON OPTICAL POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   PotType = 1 = Deuteorn Optical Potential
!   PotType = 2 = Nucleon Optical Potentials
    if (PotType == 1) then
   
      write (*,*) ''
      write(*,*) 'What Local Potential To Use? (1=User defined, 2=Pre-defined)'
      read(2,*,IOSTAT=EndOfFile) WhatPot
      if (EndOfFile /= 0) then
        read(*,*) WhatPOt
      end if
      write (Char1, '(A)')  'What Local Potential To Use? (1=User defined, 2=Pre-defined)'
      write (Val1,'(I4)') WhatPot
      write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!     WhatPot = 1 = User defined local deuteron optical potential
!     WhatPot = 2 = Pre-defined local deuteron optical potential
      if (WhatPot == 1) then

        write (*,*) ''
        write(*,*) 'Local Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Imag Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if 
        write (Char1, '(A)')  'Local Imag Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Local Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Imag Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Imag Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Imag Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Imag Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Coulomb Radius'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'Local Coulomb Radius'
        write (Val1,'(F7.3)') Item1
        write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!     END USER DEFINED LOCAL POTENTIAL
      end if

!     WhatPot = 1 = User defined local deuteron optical potential
!     WhatPot = 2 = Pre-defined local deuteron optical potential
      if (WhatPot==2) then
        write (*,*) ''
        write(*,*) 'What Pre-Defined Local Potential? (1=Daehnik)'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'What Pre-Defined Local Potential? (1=Daehnik)'
        write (Val1,'(I4)') nint(Item1)
        write(1,Format1) adjustl(val1), adjustl(trim(char1))   
      end if

!   END USER DEFINED LOCAL DEUTERON OPTICAL POTENTIAL
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   LOCAL POTENTIAL: NUCLEON OPTICAL POTENTIALS FOR THE DEUTERON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   PotType = 1 = Deuteorn Optical Potential
!   PotType = 2 = Nucleon Optical Potentials
    if (PotType == 2) then
   
      write (*,*) ''
      write(*,*) 'What Local Nucleon Potentials To Use? (1=User defined, 2=Pre-defined)'
      read(2,*,IOSTAT=EndOfFile) WhatPot
      if (EndOfFile /= 0) then
        read(*,*) WhatPot
      end if
      write (Char1, '(A)')  'What Local Nucleon Potentials To Use? (1=User defined, 2=Pre-defined)'
      write (Val1,'(I4)') WhatPot
      write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     LOCAL POTENTIAL: NEUTRON OPTICAL POTENTIAL FOR THE DEUTERON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     WhatPot = 1 = User defined local nucleon optical potentials for the deuteron
!     WhatPot = 2 = Pre-defined local nucleon optical potentials for the deuteron
      if (WhatPot == 1) then

        write (*,*) ''
        write(*,*) 'NEUTRON'
        write(*,*) 'Local Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  '(NEUTRON) Local Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Imag Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if 
        write (Char1, '(A)')  'Local Imag Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Local Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Imag Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Imag Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Imag Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Imag Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       LOCAL POTENTIAL: PROTON OPTICAL POTENTIAL FOR THE DEUTERON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write (*,*) ''
        write(*,*) 'PROTON'
        write(*,*) 'Local Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  '(PROTON) Local Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Imag Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if 
        write (Char1, '(A)')  'Local Imag Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Local Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Imag Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Imag Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Imag Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Imag Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Coulomb Radius'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'Local Coulomb Radius'
        write (Val1,'(F7.3)') Item1
        write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     END USER DEFINED LOCAL PROTON POTENTIAL FOR THE DEUTERON
      end if

!     WhatPot = 1 = User defined local nucleon optical potentials for the deuteron
!     WhatPot = 2 = Pre-defined local nucleon optical potentials for the deuteron
      if (WhatPot==2) then
        write (*,*) ''
        write(*,*) 'What Pre-Defined Local Nucleon Potentials for the Deuteron? (1=KD, 2=CH89)'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'What Pre-Defined Local Nucleon Potentials for the Deuteron? (1=KD, 2=CH89)'
        write (Val1,'(I4)') nint(Item1)
        write(1,Format1) adjustl(val1), adjustl(trim(char1)) 
      end if

!   END LOCAL DEUTERON POTENTIAL
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NONLOCALITY: DEUTERON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*,*) ''
    write(*,*) 'Use Nonlocality? [0]=No, [1]=Yes'
    read(2,*,IOSTAT=EndOfFile) Nonloc
    if (EndOfFile /= 0) then 
      read(*,*) Nonloc
    end if
    write (Char1, '(A)')  'Use Nonlocality? [0]=No, [1]=Yes'
    write (Val1,'(I4)') Nonloc
    write(1,Format1) adjustl(val1), adjustl(trim(char1))

    if (Nonloc==1) then

      write (*,*) ''
      write(*,*) 'What Type of Nonlocal Potential? (1=Deuteron Potential, 2=Nucleon Potentials)'
      read(2,*,IOSTAT=EndOfFile) PotType
      if (EndOfFile /= 0) then
        read(*,*) PotType
      end if
      write (Char1, '(A)')  'What Type of Nonlocal Potential? (1=Deuteron Potential, 2=Nucleon Potentials)'
      write (Val1,'(I4)') PotType
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     PotType = 1 = Deuteron Optical Potential
!     PotType = 2 = Nucleon Optical Potentials
      if (PotType == 1) then

        write (*,*) ''
        write(*,*) 'What Nonlocal Deuteron Optical Potential To Use? (1=User defined, 2=Pre-defined)'
        read(2,*,IOSTAT=EndOfFile) WhatPot
        if (EndOfFile /= 0) then
          read(*,*) WhatPot
        end if
        write (Char1, '(A)')  'What Nonlocal Deuteron Optical Potential To Use? (1=User defined, 2=Pre-defined)'
        write (Val1,'(I4)') WhatPot
        write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!       WhatPot = 1 = User defined
!       WhatPot = 2 = Pre-defined
!       WhatPot = 3 = Read in
        if (WhatPot == 1) then

          write (*,*) ''
          write(*,*) 'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Imag Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if 
          write (Char1, '(A)')  'Local Part of NL: Imag Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

          write (*,*) ''
          write(*,*) 'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Imag Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Imag Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Real Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Real Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Local Part of NL: Imag Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Imag Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Coulomb Radius'
          read(2,*,IOSTAT=EndOfFile) Item1
          if (EndOfFile /= 0) then 
            read(*,*) Item1
          end if
          write (Char1, '(A)')  'Local Part of NL: Coulomb Radius'
          write (Val1,'(F7.3)') Item1
          write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

          write (*,*) '' 
          write(*,*) 'Nonlocal Real Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if  
          write (Char1, '(A)')  'Nonlocal Real Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Imag Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if 
          write (Char1, '(A)')  'Nonlocal Imag Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

          write (*,*) ''
          write(*,*) 'Nonlocal Real Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Real Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Imag Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Imag Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Real Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Real Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Nonlocal Imag Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Imag Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3 
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Beta. Nonlocality parameter [fm]'
          read(2,*,IOSTAT=EndOfFile) Item1
          if (EndOfFile /= 0) then 
            read(*,*) Item1
          end if
          write (Char1, '(A)')  'Beta. Nonlocality parameter'
          write (Val1,'(F7.3)') Item1
          write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!       END USER DEFINED DEUTERON NONLOCAL OPTICAL POTENTIAL
        end if

!     END DEUTERON NONLOCAL OPTICAL POTENTIAL
      end if

!     PotType = 1 = Deuteron optical potential
!     PotType = 2 = Nucleon optical potentials for the deuteron
      if (PotType == 2) then

        write (*,*) ''
        write(*,*) 'What Nonlocal Nucleon Potentials To Use? (1=User defined, 2=Pre-defined)'
        read(2,*,IOSTAT=EndOfFile) WhatPot
        if (EndOfFile /= 0) then
          read(*,*) WhatPot
        end if
        write (Char1, '(A)')  'What Local Nucleon Potentials To Use? (1=User defined, 2=Pre-defined)'
        write (Val1,'(I4)') WhatPot
        write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!       WhatPot = 1 = User defined nonlocal nucleon optical potentials
!       WhatPot = 2 = Pre-defined nonlocal nucleon optical potentials
!       WhatPot = 3 = Read in
        if (WhatPot == 1) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         NEUTRON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          write (*,*) ''
          write(*,*) 'Local Part of NL: Neutron Real Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Neutron Real Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Neutron Imag Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if 
          write (Char1, '(A)')  'Local Part of NL: Neutron Imag Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

          write (*,*) ''
          write(*,*) 'Local Part of NL: Neutron Real Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Neutron Real Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Neutron Imag Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Neutron Imag Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Neutron Real Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Neutron Real Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Local Part of NL: Neutron Imag Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Neutron Imag Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         PROTON POTENTIAL: LOCAL PART OF NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          write (*,*) ''
          write(*,*) 'Local Part of NL: Proton Real Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if 
          write (Char1, '(A)')  'Local Part of NL: Proton Real Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Proton Imag Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if 
          write (Char1, '(A)')  'Local Part of NL: Proton Imag Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

          write (*,*) ''
          write(*,*) 'Local Part of NL: Proton Real Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Proton Real Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Proton Imag Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Proton Imag Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Local Part of NL: Proton Real Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Proton Real Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Local Part of NL: Proton Imag Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Local Part of NL: Proton Imag Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Local Part of NL: Proton Coulomb Radius'
          read(2,*,IOSTAT=EndOfFile) Item1
          if (EndOfFile /= 0) then 
            read(*,*) Item1
          end if
          write (Char1, '(A)')  'Local Part of NL: Proton Coulomb Radius'
          write (Val1,'(F7.3)') Item1
          write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         NEUTRON POTENTIAL: NONLOCAL POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          write (*,*) ''
          write(*,*) 'Nonlocal Neutron Real Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Neutron Real Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Neutron Imag Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if 
          write (Char1, '(A)')  'Nonlocal Neutron Imag Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

          write (*,*) ''
          write(*,*) 'Nonlocal Neutron Real Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Neutron Real Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Neutron Imag Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Neutron Imag Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Neutron Real Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Neutron Real Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Nonlocal Neutron Imag Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Neutron  Imag Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Neutron Nonlocality Range: Beta'
          read(2,*,IOSTAT=EndOfFile) Item1
          if (EndOfFile /= 0) then 
            read(*,*) Item1
          end if
          write (Char1, '(A)')  'Neutron Beta'
          write (Val1,'(F7.3)') Item1
          write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         PROTON POTENTIAL: NONLOCAL POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          write (*,*) ''
          write(*,*) 'Nonlocal Proton Real Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Proton Real Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1 
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Proton Imag Volume Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if 
          write (Char1, '(A)')  'Nonlocal Proton Imag Volume Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

          write (*,*) ''
          write(*,*) 'Nonlocal Proton Real Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Proton Real Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Proton Imag Surface Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Proton Imag Surface Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

          write (*,*) ''
          write(*,*) 'Nonlocal Proton Real Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Proton Real Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Nonlocal Proton Imag Spin-Orbit Depth, Radius, Diffuseness'
          read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
          if (EndOfFile /= 0) then 
            read(*,*) Item1,Item2,Item3
          end if
          write (Char1, '(A)')  'Nonlocal Proton Imag Spin-Orbit Depth, Radius, Diffuseness'
          write (Val1,'(F7.3)') Item1
          write (Val2,'(F7.3)') Item2
          write (Val3,'(F7.3)') Item3
          write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

          write (*,*) ''
          write(*,*) 'Proton Nonlocality Range: Beta'
          read(2,*,IOSTAT=EndOfFile) Item1
          if (EndOfFile /= 0) then 
            read(*,*) Item1
          end if
          write (Char1, '(A)')  'Proton Beta'
          write (Val1,'(F7.3)') Item1
          write(1,Format1) adjustl(val1), adjustl(trim(char1))

!       END USER DEFINED NUCLEON NONLOCAL POTENTIALS FOR DEUTERON
        end if

!       WhatPot = 1 = User defined nonlocal nucleon optical potentials
!       WhatPot = 2 = Pre-defined nonlocal nucleon optical potentials
!       WhatPot = 3 = Read in
        if (WhatPot == 2) then

          write (*,*) ''
          write(*,*) 'What Pre-Defined Local Potential? (1=Perey-Buck, 2=TPM)' 
          read(2,*,IOSTAT=EndOfFile) Item1
          if (EndOfFile /= 0) then 
            read(*,*) Item1
          end if
          write (Char1, '(A)')  'What Pre-Defined Local Potential? (1=Perey-Buck, 2=TPM)'
          write (Val1,'(I4)') nint(Item1)
          write(1,Format1) adjustl(val1), adjustl(trim(char1))  
        end if

!     END NONLOCAL NUCLEON OPTICAL POTENTIALS FOR THE DEUTERON
      end if

!   END NONLOCAL POTENTIAL
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

  if (What_Calc==4 .or. What_Calc==5 .or. What_Calc==6) then

    write(1,*) ''

    write (*,*) ''
    write(*,*) 'NUCLEON SCATTERING STATE'
    write(*,*) 'Mass number of projectile and target'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if      
    write (Char1, '(A)')  'Mass number of projectile and target (NUCLEON SCATTERING STATE)'
    write (Val1,'(I4)') nint(Item1)
    write (Val2,'(I4)') nint(Item2)
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Charge of projectile and target'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Charge of projectile and target'
    write (Val1,'(I4)') nint(Item1)
    write (Val2,'(I4)') nint(Item2)
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Spin of projectile and target'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Spin of projectile and target'
    write (Val1,'(F7.1)') Item1
    write (Val2,'(F7.1)') Item2
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

    write (*,*) ''
    write(*,*) 'Parity of projectile and target'
    read(2,*,IOSTAT=EndOfFile) Item1,Item2
    if (EndOfFile /= 0) then 
      read(*,*) Item1,Item2
    end if 
    write (Char1, '(A)')  'Parity of projectile and target'
    write (Val1,'(F7.1)') Item1
    write (Val2,'(F7.1)') Item2
    write(1,Format2) adjustl(val1), adjustl(val2), adjustl(trim(char1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   LOCAL POTENTIAL: NUCLEON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*,*) ''
    write(*,*) 'What Local Potential To Use? (1=User defined, 2=Pre defined)'
    read(2,*,IOSTAT=EndOfFile) WhatPot
    if (EndOfFile /= 0) then
      read(*,*) WhatPot
    end if
    write (Char1, '(A)')  'What Local Potential To Use? (1=User defined, 2=Pre-defined)'
    write (Val1,'(I4)') WhatPot
    write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!   WhatPot = 1 = User defined potential
!   WhatPot = 2 = Pre-defined potential
    if (WhatPot == 1) then

      write (*,*) ''
      write(*,*) 'Local Real Volume Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Real Volume Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

      write (*,*) ''
      write(*,*) 'Local Imag Volume Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Imag Volume Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

      write (*,*) ''
      write(*,*) 'Local Real Surface Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Real Surface Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

      write (*,*) ''
      write(*,*) 'Local Imag Surface Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Imag Surface Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

      write (*,*) ''
      write(*,*) 'Local Real Spin-Orbit Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Real Spin-Orbit Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

      write (*,*) ''
      write(*,*) 'Local Imag Spin-Orbit Depth, Radius, Diffuseness'
      read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
      if (EndOfFile /= 0) then 
        read(*,*) Item1,Item2,Item3
      end if
      write (Char1, '(A)')  'Local Imag Spin-Orbit Depth, Radius, Diffuseness'
      write (Val1,'(F7.3)') Item1
      write (Val2,'(F7.3)') Item2
      write (Val3,'(F7.3)') Item3
      write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

      write (*,*) ''
      write(*,*) 'Local Coulomb Radius'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Local Coulomb Radius'
      write (Val1,'(F7.3)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!   END USER DEFINED LOCAL POTENTIAL
    end if


!   WhatPot = 1 = User defined potential
!   WhatPot = 2 = Pre-defined potential
    if (WhatPot==2) then

      write (*,*) ''
      write(*,*) 'What Pre-Defined Local Potential? (1=KD, 2=CH89)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'What Pre-Defined Local Potential? (1=KD, 2=CH89)'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))      

!   END PRE-DEFINED POTENTIAL
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   LOCAL PART OF NONLOCAL POTENTIAL: NUCLEON SCATTERING STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*,*) ''
    write(*,*) 'Use Nonlocality? [0]=No, [1]=Yes'
    read(2,*,IOSTAT=EndOfFile) Nonloc
    if (EndOfFile /= 0) then 
      read(*,*) Nonloc
    end if
    write (Char1, '(A)')  'Use Nonlocality? [0]=No, [1]=Yes'
    write (Val1,'(I4)') Nonloc
    write(1,Format1) adjustl(val1), adjustl(trim(char1))    

!   LOCAL PART OF NONLOCAL POTENTIAL
    if (Nonloc==1) then

      write (*,*) ''
      write(*,*) 'What Nonlocal Potential? (1=User defined, 2 = Pre-defined 3=Read in)'
      read(2,*,IOSTAT=EndOfFile) WhatPot
      if (EndOfFile /= 0) then 
        read(*,*) WhatPot
      end if
      write (Char1, '(A)')  'What Nonlocal Potential? (1=User defined, 2=Pre-defined, 3=Read in)'
      write (Val1,'(I4)') WhatPot
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     WhatPot = 1 = User defined potential
!     WhatPot = 2 = Pre-defined potential
!     WhatPot = 3 = Read in
      if (WhatPot == 1) then

        write (*,*) ''
        write(*,*) 'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  
  
        write (*,*) ''
        write(*,*) 'Local Part of NL: Imag Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Imag Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  
  
        write (*,*) ''
        write(*,*) 'Local Part of NL: Imag Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Imag Surface Depth, Radius, Diffuseness' 
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Part of NL Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2  
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Part of NL: Imag Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Imag Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Part of NL: Coulomb Radius'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'Local Part of NL: Coulomb Radius'
        write (Val1,'(F7.3)') Item1
        write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       N+A NONLOCAL POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write (*,*) ''
        write(*,*) 'Nonlocal Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Nonlocal Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Nonlocal Imag Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Nonlocal Imag Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Nonlocal Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Nonlocal Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Nonlocal Imag Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Nonlocal Imag Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Nonlocal Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Nonlocal Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Nonlocal Imag Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Nonlocal Imag Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Beta: Range of Nonlocality'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'Beta: Range of Nonlocality'
        write (Val1,'(F7.3)') Item1
        write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     END USER DEFINED NONLOCAL POTENTIAL
      end if

!     WhatPot = 1 = User defined potential
!     WhatPot = 2 = Pre-defined potential
!     WhatPot = 3 = Read in
      if (WhatPot==2) then

        write (*,*) ''
        write(*,*) 'What Pre-Defined Nonlocal Potential? (1=Perey-Buck, 2=TPM)'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'What Pre-Defined Nonlocal Potential? (1=Perey-Buck, 2=TPM)'
        write (Val1,'(I4)') nint(Item1)
        write(1,Format1) adjustl(val1), adjustl(trim(char1))      

!     END PRE-DEFINED POTENTIAL
      end if

!     WhatPot = 1 = User defined potential
!     WhatPot = 2 = Pre-defined potential
!     WhatPot = 3 = Read in
      if (WhatPot==3) then
        write (*,*) ''
        write(*,*) 'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  
  
        write (*,*) ''
        write(*,*) 'Local Part of NL: Imag Volume Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Imag Volume Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))     

        write (*,*) ''
        write(*,*) 'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Surface Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  
  
        write (*,*) ''
        write(*,*) 'Local Part of NL: Imag Surface Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Imag Surface Depth, Radius, Diffuseness' 
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1))  

        write (*,*) ''
        write(*,*) 'Local Part of NL Real Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Real Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2  
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Part of NL: Imag Spin-Orbit Depth, Radius, Diffuseness'
        read(2,*,IOSTAT=EndOfFile) Item1,Item2,Item3
        if (EndOfFile /= 0) then 
          read(*,*) Item1,Item2,Item3
        end if
        write (Char1, '(A)')  'Local Part of NL: Imag Spin-Orbit Depth, Radius, Diffuseness'
        write (Val1,'(F7.3)') Item1
        write (Val2,'(F7.3)') Item2
        write (Val3,'(F7.3)') Item3
        write(1,Format3) adjustl(val1), adjustl(val2), adjustl(val3), adjustl(trim(char1)) 

        write (*,*) ''
        write(*,*) 'Local Part of NL: Coulomb Radius'
        read(2,*,IOSTAT=EndOfFile) Item1
        if (EndOfFile /= 0) then 
          read(*,*) Item1
        end if
        write (Char1, '(A)')  'Local Part of NL: Coulomb Radius'
        write (Val1,'(F7.3)') Item1
        write(1,Format1) adjustl(val1), adjustl(trim(char1))  

!     END READ IN POTENTIAL
      end if

!   END NONLOCAL
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

    write(1,*) ''

    write(*,*) ''
    write(*,*) 'Accuracy Parameters (0=Use Default, 1=Use Custom)'
    read(2,*,IOSTAT=EndOfFile) What_Accuracy
    if (EndOfFile /= 0) then 
      read(*,*) What_Accuracy
    end if
!    write (Char1, '(A)')  'What Accuracy Parameters? [0]=Default [1]=Custom'
!    write (Val1,'(I4)') What_Accuracy
!    write(1,Format1) adjustl(val1), adjustl(trim(char1))


!   WhatAccuracy = 0 = Default accuracy parameters
!   WhatAccuracy = 1 = Custom accuracy parameters
    if (What_Accuracy==0) then

!     accuracy(1)
      write (Char1, '(A)')  'Mass Unit'
      write (Val1,'(F9.4)') 931.494
      write(1,Format1) adjustl(val1), adjustl(trim(char1))   

!     accuracy(2)
      write (Char1, '(A)')  'Number of Mesh Points for Gaussian NL Integration'
      write (Val1,'(I4)') 20
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(3)
      write (Char1, '(A)')  'Step in theta when calculating the spherical harmonics'
      write (Val1,'(F8.5)') 0.00001
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(4)
      write (Char1, '(A)')  'Rmatch for n+p Bound State'
      write (Val1,'(F8.5)') 2.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(5)
      write (Char1, '(A)')  'Matching radius for nucleon bound state'
      write (Val1,'(F8.5)') 2.5
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(6)
      write (Char1, '(A)')  'Energy to start searching for bound state'
      write (Val1,'(F9.5)') -20.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(7)
      write (Char1, '(A)')  'Energy step when scanning for bound state'
      write (Val1,'(F8.5)') 0.001
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(8)
      write (Char1, '(A)')  'Energy to back up after each iteration in nonlocal bound state'
      write (Val1,'(F8.5)') 20.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(9)
      write (Char1, '(A)')  'Percent diff in energy at convergence of nonlocal bound state'
      write (Val1,'(F8.5)') 0.001
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(10)
      write (Char1, '(A)')  '% Diff of Log Deriv. at Convergence'
      write (Val1,'(F8.5)') 0.001
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(11)
      write (Char1, '(A)')  'Step size for the source in d+A nonlocal integration'
      write (Val1,'(F8.5)') 0.05
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(12)
      write (Char1, '(A)')  'Maximum radius to calculate nonlocal adiabatic source'
      write (Val1,'(F8.5)') 15.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(13)
      write (Char1, '(A)')  'Maximum radius of dr integral in nonlocal adiabatic source'
      write (Val1,'(F8.5)') 4.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(14)
      write (Char1, '(A)')  'Maximum radius of ds integral in nonlocal adiabatic source'
      write (Val1,'(F8.5)') 1.2
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(15)
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, dr Integral'
      write (Val1,'(I4)') 15
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(16)
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, theta_r Integral'
      write (Val1,'(I4)') 2
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(17)
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, ds Integral'
      write (Val1,'(I4)') 30
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(18)
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, theta_s Integral'
      write (Val1,'(I4)') 30
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(19)
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, phi_s Integral'
      write (Val1,'(I4)') 6
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(20)
      write (Char1, '(A)')  'Maximum value of dR integral in calculation of T-matrix'
      write (Val1,'(F8.5)') 20.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(21)
      write (Char1, '(A)')  'Maximum value of dr integral in calculation of T-matrix'
      write (Val1,'(F8.5)') 20.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(22)
      write (Char1, '(A)')  'Number Of Mesh Points for T-Matrix Angular Integral'
      write (Val1,'(I4)') 30
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(23)
      write (Char1, '(A)')  'Number Of Mesh Points for T-Matrix dR Radial Integral'
      write (Val1,'(I4)') 30
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(24)
      write (Char1, '(A)')  'Number Of Mesh Points for T-Matrix dr Radial Integral'
      write (Val1,'(I4)') 30
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(25)
      write (Char1, '(A)')  'Angular step when calculating transfer cross section'
      write (Val1,'(F8.5)') 1.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(26)
      write (Char1, '(A)')  'Angular step when calculating elastic cross section'
      write (Val1,'(F8.5)') 1.0
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

    end if


!   WhatAccuracy = 0 = Default accuracy parameters
!   WhatAccuracy = 1 = Custom accuracy parameters
    if (What_Accuracy==1) then

!     accuracy(1)
      write (*,*) ''
      write(*,*) 'Mass Unit [MeV/c^2]? (Default = 931.494)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Mass Unit'
      write (Val1,'(F9.4)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))   

!     accuracy(2)
      write (*,*) ''
      write(*,*) 'Number of Mesh Points for Gaussian NL Integration (Default = 20)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number of Mesh Points for Gaussian NL Integration'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(3)
      write (*,*) ''
      write(*,*) 'Step in theta when calculating the spherical harmonics (Default = 0.00001)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Step in theta when calculating the spherical harmonics'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(4)
      write (*,*) ''
      write(*,*) 'Matching Radius for for n+p Bound State [fm]? (Default = 2.0)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Rmatch for n+p Bound State'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(5)
      write (*,*) ''
      write(*,*) 'Matching radius for nucleon bound state [fm]? (Default = 2.5)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Matching radius for nucleon bound state'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(6)
      write (*,*) ''
      write(*,*) 'Energy to start searching for bound state (Default = -20)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Energy to start searching for bound state'
      write (Val1,'(F9.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(7)
      write (*,*) ''
      write(*,*) 'Energy step when scanning for bound state (Default = 0.001)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Energy step when scanning for bound state'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(8)
      write (*,*) ''
      write(*,*) 'Energy to back up after each iteration in nonlocal bound state (Default = 20)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Energy to back up after each iteration in nonlocal bound state'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(9)
      write (*,*) ''
      write(*,*) 'Percent diff in energy at convergence of nonlocal bound state (Default = 0.001)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Percent diff in energy at convergence of nonlocal bound state'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(10)
      write (*,*) ''
      write(*,*) '% Diff of Log Deriv. at Convergence for nonlocal scattering state (Default = 0.001)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  '% Diff of Log Deriv. at Convergence'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(11)
      write (*,*) ''
      write(*,*) 'Step size for the source in d+A nonlocal integration (Default=0.05)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Step size for the source in d+A nonlocal integration'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(12)
      write (*,*) ''
      write(*,*) 'Maximum radius to calculate d+A nonlocal source (Default=12)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Maximum radius to calculate nonlocal adiabatic source'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(13)
      write (*,*) ''
      write(*,*) 'Maximum radius of dr integral in nonlocal adiabatic source'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Maximum radius of dr integral in nonlocal adiabatic source'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(14)
      write (*,*) ''
      write(*,*) 'Maximum radius of ds integral in nonlocal adiabatic source'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Maximum radius of ds integral in nonlocal adiabatic source'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(15)
      write (*,*) ''
      write(*,*) 'Number Of Mesh Points for d+A Source, dr Integral'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, dr Integral'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(16)
      write (*,*) ''
      write(*,*) 'Number Of Mesh Points for d+A Source, theta_r Integral'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, theta_r Integral'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(17)
      write (*,*) ''
      write(*,*) 'Number Of Mesh Points for d+A Source, ds Integral'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, ds Integral'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(18)
      write (*,*) ''
      write(*,*) 'Number Of Mesh Points for d+A Source, theta_s Integral'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, theta_s Integral'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(19)
      write (*,*) ''
      write(*,*) 'Number Of Mesh Points for d+A Source, phi_s Integral'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number Of Mesh Points for d+A Source, phi_s Integral'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(20)
      write (*,*) ''
      write(*,*) 'Maximum value of dR integral in calculation of T-matrix'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Maximum value of dR integral in calculation of T-matrix'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(21)
      write (*,*) ''
      write(*,*) 'Maximum value of dr integral in calculation of T-matrix'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Maximum value of dr integral in calculation of T-matrix'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(22)
      write (*,*) ''
      write(*,*) 'Number Of Mesh Points for T-Matrix Angular Integral'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number Of Mesh Points for T-Matrix Angular Integral'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(23)
      write (*,*) ''
      write(*,*) 'Number Of Mesh Points for T-Matrix dR Radial Integral'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number Of Mesh Points for T-Matrix dR Radial Integral'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(24)
      write (*,*) ''
      write(*,*) 'Number Of Mesh Points for T-Matrix dr Radial Integral'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Number Of Mesh Points for T-Matrix dr Radial Integral'
      write (Val1,'(I4)') nint(Item1)
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(25)
      write (*,*) ''
      write(*,*) 'Angular step when calculating transfer cross section (Default = 1)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Angular step when calculating transfer cross section'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

!     accuracy(26)
      write (*,*) ''
      write(*,*) 'Angular step when calculating elastic cross section (Default = 1)'
      read(2,*,IOSTAT=EndOfFile) Item1
      if (EndOfFile /= 0) then 
        read(*,*) Item1
      end if
      write (Char1, '(A)')  'Angular step when calculating elastic cross section'
      write (Val1,'(F8.5)') Item1
      write(1,Format1) adjustl(val1), adjustl(trim(char1))

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PRINTING OPTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(1,*) ''

  write(*,*) ''
  write(*,*) 'Name of directory to print output files'
  read(2,*,IOSTAT=EndOfFile) Char1
  if (EndOfFile /= 0) then 
    read(*,*) Char1
  end if
  write(*,*) Char1
  write(1,'(A)') trim(adjustl(char1))

  write(*,*) ''
  write(*,*) 'Printing Options (0=Print Nothing, 1=Print Everything)'
  read(2,*,IOSTAT=EndOfFile) What_Print
  if (EndOfFile /= 0) then 
    read(*,*) What_Print
  end if

! What_Print = 0 = Print Nothing
! What_Print = 1 = Print Everything
  if (What_Print==0) then

!     printing(1)
      write (Char1, '(A)')  'DeuteronBoundWF.txt (0=DO NOT PRINT, 1=PRINT)'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(2)
      write (Char1, '(A)')  'NucleonBoundWF.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(3)
      write (Char1, '(A)')  'DeuteronScatWFs.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(4)
      write (Char1, '(A)')  'NucleonScatWFs.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(5)
      write (Char1, '(A)')  'LocalBoundWF.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(6)
      write (Char1, '(A)')  'NonlocalBoundWF.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(7)
      write (Char1, '(A)')  'DeuteronLocalIntegral.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(8)
      write (Char1, '(A)')  'NucleonLocalIntegral.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(9)
      write (Char1, '(A)')  'DeuteronNonlocalIntegral.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(10)
      write (Char1, '(A)')  'NucleonNonlocalIntegral.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(11)
      write (Char1, '(A)')  'DeuteronLocalSmatrix.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(12)
      write (Char1, '(A)')  'NucleonLocalSmatrix.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(13)
      write (Char1, '(A)')  'DeuteronNonlocalSmatrix.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(14)
      write (Char1, '(A)')  'NucleonNonlocalSmatrix.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(15)
      write (Char1, '(A)')  'DeuteronRatioToRuth.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(16)
      write (Char1, '(A)')  'NucleonRatioToRuth.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(17)
      write (Char1, '(A)')  'DeuteronElasticCS.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(18)
      write (Char1, '(A)')  'NucleonElasticCS.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(19)
      write (Char1, '(A)')  'TransferCS.txt'
      write (Val1,'(I4)') 0
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

  end if

! What_Print = 0 = Print Nothing
! What_Print = 1 = Print Everything
  if (What_Print==1) then

!     printing(1)
      write (Char1, '(A)')  'DeuteronBoundWF.txt (0=DO NOT PRINT, 1=PRINT)'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(2)
      write (Char1, '(A)')  'NucleonBoundWF.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(3)
      write (Char1, '(A)')  'DeuteronScatWFs.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(4)
      write (Char1, '(A)')  'NucleonScatWFs.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(5)
      write (Char1, '(A)')  'LocalBoundWF.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(6)
      write (Char1, '(A)')  'NonlocalBoundWF.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(7)
      write (Char1, '(A)')  'DeuteronLocalIntegral.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(8)
      write (Char1, '(A)')  'NucleonLocalIntegral.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(9)
      write (Char1, '(A)')  'DeuteronNonlocalIntegral.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(10)
      write (Char1, '(A)')  'NucleonNonlocalIntegral.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(11)
      write (Char1, '(A)')  'DeuteronLocalSmatrix.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(12)
      write (Char1, '(A)')  'NucleonLocalSmatrix.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(13)
      write (Char1, '(A)')  'DeuteronNonlocalSmatrix.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(14)
      write (Char1, '(A)')  'NucleonNonlocalSmatrix.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(15)
      write (Char1, '(A)')  'DeuteronRatioToRuth.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(16)
      write (Char1, '(A)')  'NucleonRatioToRuth.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(17)
      write (Char1, '(A)')  'DeuteronElasticCS.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(18)
      write (Char1, '(A)')  'NucleonElasticCS.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

!     printing(19)
      write (Char1, '(A)')  'TransferCS.txt'
      write (Val1,'(I4)') 1
      write(1,Format1) adjustl(val1), adjustl(trim(char1)) 

      write(1,*) ''

  end if

  write(*,*) ''
  write(*,*) 'Thanks'

  end program
