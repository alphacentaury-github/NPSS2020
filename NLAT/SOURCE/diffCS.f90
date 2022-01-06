  subroutine diffCS(SM,Lmax,SpinProjectile,SpinTarget,eta,k,Jp_Max,channel,printing,Directory)
  implicit none

  integer::n,L,Lmax,Jp_Max,channel,printing(1:100)
  integer::print_DeuteronRatioToRuth,print_DeuteronElasticCS
  integer::print_NucleonRatioToRuth,print_NucleonElasticCS
  real(8)::yr(0:Lmax),yi(0:Lmax),SpinProjectile,SpinTarget,upi,up,uti,ut,pi,mi,thetastep,eta,Cphase,ruth,phase(0:Lmax)
  real(8)::mpi,mtot,mp,m,jp,jtot,j1,m1,j2,m2,j3,m3,cgc1,cgc2,cgc3,cgc4,x,theta,phi,CrossSection(1:180),error
  complex*16::i,SM(0:Lmax,0:Jp_Max),SH(0:Lmax),k,terms(0:Lmax,0:Jp_Max),fn,fn1,fc
  character(LEN=50)::Directory,extension
  character(LEN=100)::filename

  print_DeuteronRatioToRuth=printing(15)
  print_DeuteronElasticCS=printing(16)
  print_NucleonRatioToRuth=printing(17)
  print_NucleonElasticCS=printing(18)

!  write(*,*) 'print_NucleonElasticCS =',print_NucleonElasticCS

  if (print_DeuteronRatioToRuth == 1) then
!    open(400,file='DeuteronRatioToRuth.txt')
    extension='DeuteronRatioToRuth.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(400,file=trim(filename))
  end if
  if (print_DeuteronElasticCS == 1) then
!    open(401,file='DeuteronElasticCS.txt')
    extension='DeuteronElasticCS.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(401,file=trim(filename))
  end if
  if (print_NucleonRatioToRuth == 1) then
!    open(402,file='NucleonRatioToRuth.txt')
    extension='NucleonRatioToRuth.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(402,file=trim(filename))
  end if
  if (print_NucleonElasticCS == 1) then
!   open(403,file='NucleonElasticCS.txt')
    extension='NucleonElasticCS.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(403,file=trim(filename))
  end if

  pi=acos(-1.0)
  i=cmplx(0.0,1.0)
  mi=0.0
  thetastep=1.0
  CrossSection(:)=0.0  
  SpinTarget=0.0

  do L=0,Lmax
    call CoulombPhase(eta,L,Cphase)
    phase(L)=Cphase
!    write(*,*) 'L,Cphase', L, Cphase
  end do

  do upi=-SpinProjectile,SpinProjectile
    do up=-SpinProjectile,SpinProjectile
      do uti=-SpinTarget,SpinTarget
        do ut=-SpinTarget,SpinTarget
          mpi=upi
          mtot=upi+uti
          mp=mtot-ut
          m=mp-up
          
          do L=0,Lmax
            do Jp=(L-SpinProjectile),(L+SpinProjectile)
              if (Jp>-0.001) then
                jtot=Jp
            
                if (abs(upi)<=jp .and. abs(upi+uti)<=jtot .and. abs(upi+uti-ut-up)<=L .and. abs(upi+uti-ut)<=jp) then
                  j1=L
                  m1=mi
                  j2=SpinProjectile
                  m2=upi
                  j3=jp
                  m3=mpi
                  call cg(j1,m1,j2,m2,j3,m3,cgc1)

                  j1=jp
                  m1=mpi
                  j2=SpinTarget
                  m2=uti
                  j3=jtot
                  m3=mtot
                  call cg(j1,m1,j2,m2,j3,m3,cgc2)

                  j1=L
                  m1=m
                  j2=SpinProjectile
                  m2=up
                  j3=jp
                  m3=mp
                  call cg(j1,m1,j2,m2,j3,m3,cgc3)

                  j1=jp
                  m1=mp
                  j2=SpinTarget
                  m2=ut
                  j3=jtot
                  m3=mtot
                  call cg(j1,m1,j2,m2,j3,m3,cgc4)             
              
                  terms(L,nint(2*Jp))=((i*sqrt(pi))/k)*cgc1*cgc2*cgc3*cgc4 &
                    *exp(2.0*i*phase(L))*(1-SM(L,nint(2*Jp)))*sqrt(2.0*L+1.0)

                end if

                if (abs(upi)>jp .or. abs(upi+uti)>jtot .or. abs(upi+uti-ut-up)>L .or. abs(upi+uti-ut)>jp) then
                  terms(int(L),int(2*jp))=0.0
                end if

!             END Jp IF
              end if
!           END Jp LOOP
            end do
!         END L LOOP
          end do
          
!         CALCULATES DIFFERENTIAL CROSS SECTION              
!         LOOP OVER ALL ANGLES
          do n=1,180
            x=n*thetastep
            theta=(x/180.0)*pi  
            fn=0
            phi=0.0 
            mi=0

!           CALCULATES ALL SPHERICAL HARMONICS AT A GIVEN ANGLE 
            if (Lmax-abs(nint(m))>0.0001) then
              call spherical_harmonic(Lmax,nint(m),theta,phi,yr,yi)
            end if
            do L=0,Lmax
              SH(L)=yr(L)+i*yi(L) 
            end do

            fn=0
            fn1=0        
!           LOOP OVER ALL L,JP,JTOT VALUES
            do L=0,Lmax
              do Jp=(L-SpinProjectile),(L+SpinProjectile)
                if (Jp>-0.001) then
                  jtot=Jp
                  fn1=terms(L,nint(2*Jp))*SH(L)
                  fn=fn+(fn1)
!               END Jp IF
                end if
!             END Jp LOOP
              end do
!           END OF L LOOP
            end do

            if(upi==up .and. uti==ut) then
              fc=((-eta)/(2*k*sin(theta/2.0)**2))*exp(-i*eta*log(sin(theta/2.0)**2)+2*i*phase(0))
              CrossSection(n)=CrossSection(n)+abs(fn+fc)**2
            end if
         
            if(upi.ne.up .or. uti.ne.ut) then
              CrossSection(n)=CrossSection(n)+abs(fn)**2
            end if

!         END OF THETA LOOP
          end do

!       END ut LOOP
        end do

!     END uti LOOP
      end do

!   END up LOOP
    end do

! END upi LOOP
  end do

  do n=1,180        
    x=n*thetastep
    theta=(x/180.0)*pi   
    CrossSection(n)=10.0*CrossSection(n)/((2*SpinProjectile+1)*(2*SpinTarget+1))

!   CALCULATE RUTHERFORD CROSS SECTION
    ruth=(eta**2)/(4*(k**2)*(sin(theta/2.0))**4)

    if (channel==1) then
      error=(CrossSection(n)/(10.0*ruth))*0.001
      if (print_DeuteronRatioToRuth == 1) then
        write(400,*),x,(CrossSection(n)/(10.0*ruth)),error
      end if
      error=CrossSection(n)*0.001
      if (print_DeuteronElasticCS == 1) then
        write(401,*),x,CrossSection(n),error
      end if
    end if

    if (channel==2) then
      error=(CrossSection(n)/(10.0*ruth))*0.001
      if (print_NucleonRatioToRuth == 1) then
        write(402,*),x,(CrossSection(n)/(10.0*ruth)),error
      end if
      error=CrossSection(n)*0.001
      if (print_NucleonElasticCS == 1) then
        write(403,*),x,CrossSection(n),error
      end if
    end if

  end do

  return
  end subroutine
