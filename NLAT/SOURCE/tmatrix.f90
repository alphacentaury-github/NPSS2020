

  subroutine tmatrix (accuracy,common,DeuteronBoundParameters,NucleonBoundParameters,DeuteronScatParameters, &
                      NucleonScatParameters,DeuteronBoundWF,NucleonBoundWF,NucleonScatWFs,DeuteronScatWFs, &
                      k,k_d,mass,mass_d,Lmax,nmax,eta,eta_d,StepSize,Vnp,SH_index_max,SH_step,Jp_dA_max,Jp_pB_max, &
                      WhatCalc,printing,FinalBoundL,Directory)
  IMPLICIT NONE
  integer::n,L,npoints,rnp_index,theta_cm,Lmax
  integer::RpB_index,rnA_index,u_index,rnA_index1,rnA_npoints,RdA_index
  integer::theta_index,theta_RdA_index,RpB_npoints,RpB_index1,s,nmax
  integer::SH_index_max,SH_index,FinalBoundL,theta_rnA_index
  integer::Li,Lf,InitialBoundL,Mk
  integer::Jp_dA_max,Jp_pB_max,gg,calc,WhatCalc,SH_Lmax,printing(1:100)
  integer::print_TransferCS
  real(8)::StepSize,TransferCS,CS_step
  real(8)::accuracy(1:100),common(1:100)
  real(8)::DeuteronBoundParameters(1:15,1:100),NucleonBoundParameters(1:15,1:100)
  real(8)::DeuteronScatParameters(1:15,1:100),NucleonScatParameters(1:15,1:100)
  real(8)::R,pi,time,minimum,maximum,eta,eta_d
  real(8)::RpB_max,rnA_max,RpB,h_RpB,u,theta,rnA,RdA,rnp,MA,theta_RdA,jp,mass,mass_d
  real(8)::Ip,It,Ip_d,It_d,hbarc,MB,Md,Mn,mass_factor,SH_step
  real(8)::u_points(200),u_weights(200),rnA_weights(200),rnA_points(200),yr(0:50),yi(0:50)
  real(8)::RpB_weights(200),RpB_points(200),DeuteronBoundWF(nmax)
  real(8)::Vnp(nmax),cgc4,cgc5
  real(8)::phase(0:Lmax),phase_d(0:Lmax)
  real(8)::DeuteronWF(nmax),NucleonBoundWF(nmax),m
  real(8)::SHm(0:50,-FinalBoundL:FinalBoundL,1:SH_index_max),sum_integral,cgc,cgc3
  real(8)::j1,m1,j2,m2,j3,m3,cgc1,cgc2,ProtonSpin,NeutronSpin
  real(8)::a,b,c,d,e,f,g,h,j,ninej,mi,mf,phi
  real(8)::Cphase,Qmin,Qmax,wig9j,cleb6
  real(8)::Jpi,Jpf,Q,MQ,mg,ninej1,ninej2,jf,ji,eps,constant
  complex*16::NucleonScatWFs(0:Lmax,0:Jp_pB_Max,1:nmax)
  complex*16::DeuteronScatWFs(0:Lmax,0:Jp_dA_Max,1:nmax)
  complex*16::i,integral,k,k_d,amplitude,term,ratio,Tint,g_sum,LJ_term
  complex*16::In(0:Lmax,Jp_dA_Max,0:Lmax,Jp_pB_Max),SH_sum
!  complex*16::T(1:3,-3:3,-3:3)
  complex*16,allocatable::T(:,:,:)
  character(LEN=50)::Directory,extension
  character(LEN=100)::filename

  print_TransferCS=printing(19)

  if (print_TransferCS == 1) then
!    open(500,file='TransferCS.txt')
    extension='TransferCS.txt'
    filename=trim(Directory) // "/" // trim(extension)
    open(500,file=trim(filename))
  end if

  pi=acos(-1.0)
  i=cmplx(0.0,1.0)
  hbarc=197.325
  eps=1e-6
  SH_Lmax=50

  SH_step=accuracy(3)
  RpB_max=accuracy(20)
  rnA_max=accuracy(21)
  npoints=accuracy(22)
  RpB_npoints=accuracy(23)
  rnA_npoints=accuracy(24)
  CS_step=accuracy(26)

  SH_index_max=nint(pi/SH_step)
  InitialBoundL=nint(DeuteronBoundParameters(1,9))
  ji=DeuteronBoundParameters(1,10)
!  FinalBoundL=nint(NucleonBoundParameters(1,9))
!  write(*,*) 'FinalBoundL = ',FinalBoundL
  jf=NucleonBoundParameters(1,10)
  ProtonSpin=NucleonScatParameters(1,5)
  NeutronSpin=NucleonBoundParameters(1,5)
  MA=DeuteronScatParameters(1,2)
  MB=NucleonScatParameters(1,2)
  Md=DeuteronScatParameters(1,1)
  Mn=NucleonScatParameters(1,1)
 
!  Ip=ProtonSpin
  Ip=NucleonScatParameters(1,5)
  It=NucleonScatParameters(1,6)
  Ip_d=DeuteronScatParameters(1,5)
  It_d=DeuteronScatParameters(1,6)

  mass_factor=(1-(Mn*MA)/(Md*MB))

  write(*,*) ''
  write(*,*) 'Transfer Cross Section'
!  write(*,*) 'SH_step=',SH_step
!  write(*,*) 'SH_index_max = ', SH_index_max
!  write(*,*) 'k,k_d=',k,k_d
!  write(*,*) 'mass, mass_d = ', mass, mass_d
!  write(*,*) 'Lmax = ', Lmax
!  write(*,*) 'nmax = ', nmax
!  write(*,*) 'eta = ', eta
!  write(*,*) 'eta_d = ', eta_d
!  write(*,*) 'StepSize = ', StepSize
!  write(*,*) 'Jp_dA_max = ', Jp_dA_max
!  write(*,*) 'Jp_pB_max = ', Jp_pB_max
!  write(*,*) 'InitialBoundL=',InitialBoundL
!  write(*,*) 'ji=',ji
!  write(*,*) 'FinalBoundL=',FinalBoundL
!  write(*,*) 'jf=',jf
!  write(*,*) 'ProtonSpin=',ProtonSpin
!  write(*,*) 'NeutronSpin=',NeutronSpin
!  write(*,*) 'MA=',MA
!  write(*,*) 'MB=',MB
!  write(*,*) 'Md=',Md
!  write(*,*) 'Mn=',Mn  
!  write(*,*) 'Ip = ', Ip
!  write(*,*) 'It=',It
!  write(*,*) 'Ip_d=',Ip_d
!  write(*,*) 'It_d=',It_d
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENRALIZE BOUNDS ON SPHERICAL HARMONICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONSTRUCT SPHERICAL HARMONICS FOR EACH LM
!  do M=0,5
  do M=0,FinalBoundL
    do SH_index=1,SH_index_max
      theta=SH_index*SH_step
      call spherical_harmonic(SH_Lmax,nint(m),theta,0.0,yr,yi)  
      do L=0,SH_Lmax 
        SHm(L,M,SH_index)=yr(L)
      end do
    end do
  end do
  
!! GET RID OF THIS
!! NORMALIZE SO THAT <PHI|PHI> = 1
!  integral=0.0
!  do n=1,nmax
!    r=n*StepSize
!    integral=integral+StepSize*(DeuteronBoundWF(n)**2.0)
!  end do
!  print *,'integral=',integral
!  DeuteronBoundWF(:)=DeuteronBoundWF(:)/(sqrt(abs(integral)))

! CHANGE PHI_d -> PHI_d/R
  do n=1,nmax
    DeuteronBoundWF(n)=DeuteronBoundWF(n)/(n*StepSize)
  end do

! CHANGE PHI_nA -> PHI_nA/R
  do n=1,nmax
    NucleonBoundWF(n)=NucleonBoundWF(n)/(n*StepSize)
  end do

! CALCULATE COULOMB PHASE FOR PROTON CHANNEL  
  do L=0,Lmax
    call CoulombPhase(eta,L,Cphase)
    phase(L)=Cphase
  end do
! CALCULATE COULOMB PHASE FOR DEUTERON CHANNEL  
  do L=0,Lmax
    call CoulombPhase(eta_d,L,Cphase)
    phase_d(L)=Cphase
  end do

  do L=Lmax,0,-1
    phase(L)=phase(L)-phase(0)
    phase_d(L)=phase_d(L)-phase_d(0)
!    write(*,*) 'phase,phase_d=',phase(L),phase_d(L)
  end do

  call logfac(100)

! CONSTRUCT GRID FOR INTEGRATION
  minimum=-1.0
  maximum=1.0
  call gauss(minimum,maximum,npoints,u_points,u_weights)

  minimum=StepSize
  maximum=rnA_max
  call gauss(minimum,maximum,rnA_npoints,rnA_points,rnA_weights)

  minimum=StepSize
  maximum=RpB_max
  call gauss(minimum,maximum,RpB_npoints,RpB_points,RpB_weights)

!  Qmin=abs(Ip_d-ProtonSpin)
!  Qmax=Ip_d+ProtonSpin
  Qmin=abs(Ip_d-Ip)
  Qmax=Ip_d+Ip
  allocate(T(1:nint(2*Qmax),-nint(2*Qmax):nint(2*Qmax),-nint(2*jf):nint(2*jf)))

  constant=32*pi**(3.0)/(sqrt(4.0*pi)*k*k_d)*(-1.0)**(3.0*ProtonSpin+3.0*NeutronSpin+ji+Ip_d+3*jf) &
           *sqrt(2.0*jf+1.0)*sqrt(2.0*Ip_d+1.0)/(sqrt(2.0*FinalBoundL+1.0))

! CALCULATE THE TRANSFER CROSS SECTION
  do theta_cm=1,179,CS_step

!    call cpu_time(time)
!    write(*,*) 'time=',time

    theta_index=nint(theta_cm*pi/(SH_step*180.0))

    sum_integral=0.0

!    calc=1

!   ONLY NEED Q=0.5,0.5 FOR NO SPIN-ORBIT CASE
!   Q=1.5 NEEDED WHEN SPIN-ORBIT IS USED
    do Q=Qmin,Qmax
      do MQ=-Q,Q
        do mf=-jf,jf
           
          if (Q==Qmin .and. MQ==-Q .and. mf==-jf) then
            calc=1
          else
            calc=0
          end if

          integral=0.0
!         Jpi=JTOT WHEN SPIN OF TARGET IN DEUTERON CHANNEL IS 0
          do Jpi=0,(Lmax+Ip_d)
            do Li=(Jpi-Ip_d), (Jpi+Ip_d)
              if (Li > -eps .and. Li <= Lmax) then
                do Lf=(Li-FinalBoundL),(Li+FinalBoundL)
                  if (Lf >= 0 .and. Lf <= Lmax) then
                    do Jpf=(Lf-Ip),(Lf+Ip)                        
                      if (Jpf > -eps) then

                        LJ_term=0.0

                        a=ProtonSpin
                        b=Ip_d
                        c=ji
                        d=Lf
                        e=Li
                        f=FinalBoundL
                        g=Jpf
                        h=Jpi
                        j=jf
                        ninej1=wig9j(a,b,c,d,e,f,g,h,j)

                        LJ_term=exp(i*(phase(Lf)+phase_d(Li)))*(i**(Lf-Li)) &
                                *sqrt(2.0*Li+1)*sqrt(2.0*Lf+1)*(2.0*Jpi+1.0)*(2.0*Jpf+1.0)*ninej1


                        if (abs(ninej1)> eps) then
                          
                          if (calc==1) then
                            Tint=0.0
 
!                           LOOP OVER RpB, RELATIVE VECTOR BETWEEN PROTON AND TARGET IN FINAL STATE
                            do RpB_index=1,RpB_npoints
                              RpB=RpB_points(RpB_index)
                              RpB_index1=nint(RpB/StepSize)
  
!                             LOOP OVER r, RELATIVE VECTOR BETWEEN PROTON AND NEUTERON          
                              do rnA_index=1,rnA_npoints
                                rnA=rnA_points(rnA_index)
                                rnA_index1=nint(rnA/StepSize)
   
!                               LOOP THROUGH u=COS(THETA)
                                do u_index=1,npoints
                                  u=u_points(u_index)
                                  theta=acos(u)
                                  theta_rnA_index=nint(theta/SH_step+1)
      
!                                 CALCULATE MAGNITUDE OF THE VECTORS, AND ANGLE OF RdA
                                  rnp=sqrt(  ((MA/MB)*rnA*sin(theta))**2.0   +  ((MA/MB)*rnA*cos(theta)-RpB)**2.0)
                                  RdA=sqrt((mass_factor*rnA*sin(theta))**2.0 &
                                      +(rnA*cos(theta)*mass_factor+(Mn/Md)*RpB)**2.0)
                                  theta_RdA=acos( (rnA*cos(theta)*mass_factor+(Mn/Md)*RpB)/RdA)        
                                           
!                                 CALCULATE INDEX OF THE VARIOUS ARRAYS
                                  theta_RdA_index=nint(theta_RdA/SH_step+1)           
                                  rnp_index=nint(rnp/StepSize)
                                  RdA_index=nint(RdA/StepSize)
   
!                                 MAKE SURE THE INDICES ARE NOT 0 
                                  if (rnp_index <1 ) then
                                    rnp_index=1
                                  end if
                                  if (rnp_index > nmax) then             
                                    rnp_index=nmax
                                  end if
                                  if (RdA_index < 1) then
                                    RdA_index=1
                                  end if
                                    
                                  SH_sum=0.0
                                  do Mk=0,FinalBoundL
                                     
                                    if (Mk==0) then
                                      j1=Lf
                                      m1=0.0
                                      j2=Li
                                      m2=-Mk
                                      j3=FinalBoundL
                                      m3=-Mk
                                      cgc3=cleb6(j1,m1,j2,m2,j3,m3)
                                      SH_sum=SH_sum+cgc3*SHm(Li,Mk,theta_RdA_index)*SHm(FinalBoundL,Mk,theta_rnA_index)
                                    end if
                                      
                                    if (Mk .ne. 0 .and. mod(Lf+Li-FinalBoundL,2)==0) then
                                      j1=Lf
                                      m1=0.0
                                      j2=Li
                                      m2=-Mk
                                      j3=FinalBoundL
                                      m3=-Mk
                                      cgc3=cleb6(j1,m1,j2,m2,j3,m3)
                                      SH_sum=SH_sum+(1.0+(-1.0)**(Lf+Li-FinalBoundL))*cgc3 &
                                             *SHm(Li,Mk,theta_RdA_index) &
                                             *SHm(FinalBoundL,Mk,theta_rnA_index)
                                    end if

                                  end do

                                  term=NucleonBoundWF(rnA_index1)*DeuteronBoundWF(rnp_index)*Vnp(rnp_index) &
                                       *NucleonScatWFs(Lf,nint(2*Jpf),RpB_index1) &
                                       *DeuteronScatWFs(Li,nint(2*Jpi),RdA_index) &
                                       *((rnA**2.0)*(RpB)/(RdA))*SH_sum &
                                       *u_weights(u_index)*rnA_weights(rnA_index)*RpB_weights(RpB_index) 
    
                                  Tint = Tint+term
   
!                               END u LOOP
                                end do
    
!                             END rnA LOOP
                              end do
   
!                           END RpB LOOP
                            end do
                          
                            In(Li,nint(2*Jpi),Lf,nint(2*Jpf))=Tint
!                            write(*,*) 'Tint=', Tint
                          
!                         END CALC IF
                          end if
   
                          g_sum=0.0
                          do gg=0,abs(Qmax+jf)

                             mg=MQ-mf

                             if ((abs(mg)-gg)<eps) then 
                           
                              a=Lf
                              b=Li
                              c=gg
                              d=Jpf
                              e=Jpi
                              f=jf
                              g=ProtonSpin
                              h=Ip_d
                              j=Q
                              ninej2=wig9j(a,b,c,d,e,f,g,h,j)
 
                              if (abs(ninej2)> eps) then
   
                                j1=gg
                                m1=mg
                                j2=jf 
                                m2=mf
                                j3=Q
                                m3=MQ
                                cgc1=cleb6(j1,m1,j2,m2,j3,m3)
                              
                                j1=Lf
                                m1=mg
                                j2=Li
                                m2=0.0
                                j3=gg
                                m3=mg
                                cgc2=cleb6(j1,m1,j2,m2,j3,m3)

                                if (mg<0) then
                                  g_sum=g_sum+sqrt(2.0*gg+1.0)*ninej2*cgc1*cgc2*(-1.0)**(abs(mg))*SHm(Lf,abs(mg),theta_index)
                                else
                                  g_sum=g_sum+sqrt(2.0*gg+1.0)*ninej2*cgc1*cgc2*SHm(Lf,mg,theta_index)
                                end if

!                             END ninej2 IF
                              end if

!                           END mg IF
                            end if

!                         END gg LOOP
                          end do

                          integral=integral+constant*LJ_term*g_sum*In(Li,nint(2*Jpi),Lf,nint(2*Jpf))
                         
!                       END ninej1 IF
                        end if

!                     END Jpf IF
                      end if

!                   END Jpf LOOP
                    end do

!                 END Lf IF
                  end if

!               END Lf LOOP
                end do

!             END Li IF
              end if

!           END Li LOOP
            end do

!         END Jpi LOOP
          end do

!          I(Li,nint(2*Jpi),Lf,2*Jpf)=Tint

          T(nint(2*Q),nint(2*MQ),nint(2*mf))=integral
          sum_integral=sum_integral+abs(integral)**2.0

!       END mf LOOP
        end do

!     END MQ LOOP
      end do

!   END Q LOOP
    end do

    TransferCS=10.0*abs((k/k_d)*(mass*mass_d/(4.0*pi**2.0*hbarc**4.0)))*sum_integral &
        *((2.0*It+1.0)/((2.0*Ip_d+1.0)*(2.0*jf+1.0)*(2.0*It_d+1.0)))

!   WhatCalc = 5 = (d,N)
!   WhatCalc = 6 = (N,d)
    if (WhatCalc == 5) then
      if (print_TransferCS == 1) then
        write(500,*),theta_cm,TransferCS
      end if
      write(*,*),theta_cm,TransferCS
    end if
    if (WhatCalc == 6) then
      if (print_TransferCS == 1) then
        write(500,*),theta_cm,abs((k_d**2.0*(2.0*Ip_d+1)*(2.0*It_d+1))/(k**2.0*(2.0*Ip+1)*(2.0*It+1))*TransferCS)
      end if
      write(*,*),theta_cm,abs((k_d**2.0*(2.0*Ip_d+1)*(2.0*It_d+1))/(k**2.0*(2.0*Ip+1)*(2.0*It+1))*TransferCS)
    end if


! END THETA_CM LOOP
  end do

  call cpu_time(time)  
!  print *,'time in tmatrix=',time

  end subroutine tmatrix
