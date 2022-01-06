  subroutine nonlocal_adiabatic(ScatParameters,accuracy,NmaxShort,DeuteronBoundWF,ScatWF,NmaxLong,Source,Vpot,SchEqConst,L, &
                                SH,SH_step,Coulomb,Vnp,StepSizeShort,StepSizeLong, &
                                SH_index_max,SH_Lmax,NLint,jp)
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::n,L,R_index_max,lr_index_max,WhatPot,WhatPreDefPot,Lmax
  integer::NmaxShort,NmaxLong,q,R_index,SH_index_max,SH_Lmax
  integer::phi_npoints,lr_npoints,s_npoints,lr_index,lr_index1,s_index,ur_index,us_index
  integer::phi_index,Pot_Index1_Rn,Pot_Index1_Rp,Pot_Index2_Rn,Pot_Index2_Rp,Bound_Index_p
  integer::Bound_Index_n,WF_Index,SH_index,ur_npoints,us_npoints
  integer::R_WF_Index
  real(8)::ur_points(200),ur_weights(200),phi_points(200),phi_weights(200),angle_factor,StepSizeSend
  real(8)::lr_points(200),lr_weights(200),s_points(200),s_weights(200),us_points(200),us_weights(200)
  real(8)::pi,StepSizeLong,R_max,lr_max,s_max,beta,Vnp(NmaxShort),StepSizeShort
  real(8)::DeuteronBoundWF(NmaxShort),integral,lr,R,Rprime,p,minimum,maximum,m,time,s,ur,theta_r
  real(8)::us,theta_s,phi_s,theta_Rs,Y_Rs,lr_dot_R,Rn_dot_s,Rp_dot_s,lr_dot_s,R_dot_s
  real(8)::Rn,Rp,Rn2s,Rp2s,exponent,cos_alpha,factor,SchEqConst,SH_step
  real(8)::SH(0:SH_Lmax,1:SH_index_max),Coulomb(NmaxLong)
  real(8)::StepSize,jp,NeutronBeta,ProtonBeta
  real(8)::ScatParameters(1:15,1:100),accuracy(1:100)
  complex*16::i,integral_n,integral_p,potential_sum
  complex*16::Source(NmaxLong),Vpot(NmaxLong),Unl,UlocNL,LocalPartOfNL(NmaxShort)
  complex*16::TempScatWF(NmaxShort),ScatWF(NmaxShort),NLint(NmaxLong),UlocNL1,UlocNL2
  complex*16::LocalPartOfNL1(NmaxShort),LocalPartOfNL2(NmaxShort),Unl1,Unl2
  complex*16::NeutronUnl(1:NmaxShort,1:NmaxShort),ProtonUnl(1:NmaxShort,1:NmaxShort)

!  open(1100,file='TempSource.txt')

!  do n=1,NmaxShort
!    write(1101,*) n*StepSizeShort,DeuteronBoundWF(n)
!  end do  

  pi=acos(-1.0)
  i=cmplx(0.0,1.0)
  NLint(:)=0.0
  Source(:)=0.0

!  write(*,*) ''
!  write(*,*) 'Inside nonlocal adiabatic'
!  write(*,*) 'StepSizeShort = ', StepSizeShort
!  write(*,*) 'NmaxShort = ', NmaxShort
!  write(*,*) 'StepSizeLong = ', StepSizeLong

! MAXIMUM RADIUS TO CALCULATE SOURCE
  R_max=accuracy(12)

! MAXIMUM RADIUS OF r AND s INTEGRALS
  lr_max=accuracy(13)
  s_max=accuracy(14)

  lr_npoints=accuracy(15)
  ur_npoints=accuracy(16)
  s_npoints=accuracy(17)
  us_npoints=accuracy(18)
  phi_npoints=accuracy(19)


!  write(*,*) 'R_max = ', R_max
!  write(*,*) 'lr_max = ', lr_max
!  write(*,*) 's_max = ', s_max
!  write(*,*) 'lr_npoints = ', lr_npoints
!  write(*,*) 'ur_npoints = ', ur_npoints
!  write(*,*) 's_npoints = ', s_npoints
!  write(*,*) 'us_npoints = ', us_npoints
!  write(*,*) 'phi_npoints = ', phi_npoints

!! CONSTRUCT SPHERICAL HARMONICS
!! SH DEFINED AS REAL SINCE M=0 !!!
!  m=0
!  do SH_index=1,SH_index_max
!    theta=SH_index*SH_step
!!    call spherical_harmonic(Lmax,nint(m),theta,0.0,yr,yi)  
!    call spherical_harmonic(Lmax,nint(m),theta,0.0,yr,yi)  
!    do L=0,Lmax
!      SH(L,SH_index)=yr(L) 
!    end do
!  end do

  WhatPot=nint(ScatParameters(3,3))
  WhatPreDefPot=nint(ScatParameters(3,4))

! WhatPot = 1 = User defined nonlocal nucleon potentials for the deuteron
! WhatPot = 2 = Pre-defined nonlocal nucleon potentials for the deuteron
  if (WhatPot==1) then
    NeutronBeta=ScatParameters(11,19)
    ProtonBeta=ScatParameters(12,19)
  end if
  if (WhatPot==2) then
!   WhatPreDefPot = 1 = Perey-Buck
!   WhatPreDefPot = 2 = TPM
    if (WhatPreDefPot==1) then
      NeutronBeta=0.85
      ProtonBeta=0.85
    end if
    if (WhatPreDefPot==2) then
      NeutronBeta=0.90
      ProtonBeta=0.88
    end if
  end if

  R_index_max=nint(R_max/StepSizeLong)
  lr_index_max=nint(lr_max/StepSizeShort)

! GET RID OF THIS. PUT THE INTEGRAL TERM AT THE END LIKE IN THE LOCAL CASE
! NORMALIZE DEUTERON WAVE FUNCTION  ->   <PHI|Vnp|PHI>=-1
  integral=0.0
  do n=1,NmaxShort
    lr=n*StepSizeShort
    integral=integral+StepSizeShort*(lr**2.0)*Vnp(n)*(DeuteronBoundWF(n)**2.0)
  end do
!  write(*,*) 'integral = ', integral
!  DeuteronBoundWF(:)=DeuteronBoundWF(:)/(sqrt(abs(integral)))

! CHANGE PSI -> PSI/R TO BE USED IN CALCULATIONS BELOW
  do n=1,NmaxShort
    TempScatWF(n)=ScatWF(n)/(n*StepSizeShort) 
  end do

! CONSTRUCT THE WOODS-SAXON PART OF THE NONLOCAL POTENTIAL
  do n=1,NmaxShort
    R=n*StepSizeShort
    do q=1,NmaxShort          
      Rprime=q*StepSizeShort
      call nonlocal_scattering_potential(ScatParameters,Unl1,Unl2,L,jp,R,Rprime,UlocNL1,UlocNL2)
      NeutronUnl(n,q)=Unl1
      ProtonUnl(n,q)=Unl2
    end do
    LocalPartOfNL1(n)=UlocNL1
    LocalPartOfNL2(n)=UlocNL2
  end do

  LocalPartOfNL(:)=LocalPartOfNL1(:)+LocalPartOfNL2(:)

! ALMOST NO SENSITIVITY FOR 2 -> 10
! VERY DIFFERENT FOR 1 VS. 2
! 2 -> 4 NEARLY NO EFFECT ON TRANSFER CS
  minimum=-1.0
  maximum=1.0
  call gauss(minimum,maximum,ur_npoints,ur_points,ur_weights)

! 25 -> 40: SOME EFFECT ON TRANSFER CS. SIMILAR TO S
  minimum=-1.0
  maximum=1.0
  call gauss(minimum,maximum,us_npoints,us_points,us_weights)

! 10 -> 30: NO IMPROVEMENT
  minimum=0.0
  maximum=2.0*pi
  call gauss(minimum,maximum,phi_npoints,phi_points,phi_weights)

! 15 -> 25 DOES NEARLY NOTHING TO TRANSFER CS AT 0.05
  minimum=StepSizeShort
  maximum=lr_max
  call gauss(minimum,maximum,lr_npoints,lr_points,lr_weights)

! 20 -> 30 HELPS WIGGLES AT SMALL R FOR LARGE L
! PUTTING THIS AT 200 COMPLETELY ALMOST REMOVES WIGGLES AT SMALL R FOR LARGE L
! 20 -> 30 DOES A LITTLE BIT TO THE TRANSFER CS AT 0.05
! ONLY S_NPOINTS >= 25 IS SUFFICIENT IN LOGARITHMIC VIEW
  minimum=StepSizeShort
  maximum=s_max
  call gauss(minimum,maximum,s_npoints,s_points,s_weights)  

! LOOP OVER R, RELATIVE VECTOR BETWEEN DEUTERON AND TARGET
  do R_index=1,R_index_max
    R=R_index*StepSizeLong
    R_WF_Index=nint(R/StepsizeShort)

    integral_n=0.0
    integral_p=0.0

!   LOOP OVER r, RELATIVE VECTOR BETWEEN PROTON AND NEUTERON          
    do lr_index=1,lr_npoints
      lr=lr_points(lr_index)
      lr_index1=nint(lr/StepSizeShort)
        
!     LOOP THROUGH ur=COS(THETA_r)
      do ur_index=1,ur_npoints
        ur=ur_points(ur_index)
        theta_r=acos(ur)
        lr_dot_R=lr*R*ur
        Rn=sqrt(lr**2.0/(4.0)+R**2.0-lr_dot_R)
        Rp=sqrt(lr**2.0/(4.0)+R**2.0+lr_dot_R)
        Pot_Index1_Rn=nint(Rn/StepSizeShort)
        Pot_Index1_Rp=nint(Rp/StepSizeShort)

!        write(*,*) '1Rn,1Rp=',Pot_Index1_Rn,Pot_Index1_Rp

        if (Pot_Index1_Rn <1) then
          Pot_Index1_Rn=1
        end if
        if (Pot_Index1_Rp <1) then
          Pot_Index1_Rp=1
        end if

!       LOOP THROUGH s=R'-R
        do s_index=1,s_npoints
          s=s_points(s_index)

!         LOOP THROUGH us=COS(THETA_s)
          do us_index=1,us_npoints
            us=us_points(us_index)
            theta_s=acos(us)
            R_dot_s=R*s*us
            WF_Index=nint(sqrt(R**2.0+s**2.0+2.0*R_dot_s)/StepSizeShort)
            if (WF_Index<1) then
              WF_Index=1
            end if

!           FIND ANGLE OF R+s
            theta_Rs=acos((R+s*us)/(sqrt(s**2.0+R**2.0+2.0*s*R*us)))
            SH_index=nint(theta_Rs/SH_step)
            if (SH_index < 1) then
              SH_index=1
            end if          
            Y_Rs=SH(L,SH_index)

!           LOOP THROUGH PHI_s
            do phi_index=1,phi_npoints

              phi_s=phi_points(phi_index)
              angle_factor=sin(theta_s)*sin(theta_r)*cos(phi_s)+cos(theta_s)*cos(theta_r)

!             CALCULATE DOT PRODUCTS
              Rn_dot_s=-((lr*s)/2.0)*angle_factor+s*R*us
              Rp_dot_s=((lr*s)/2.0)*angle_factor+s*R*us
              lr_dot_s=lr*s*angle_factor

!             CALCULATE MAGNITUDE OF THE VECTORS
              Rn2s=sqrt(Rn**2.0+4.0*s**2.0+4.0*Rn_dot_s)
              Rp2s=sqrt(Rp**2.0+4.0*s**2.0+4.0*Rp_dot_s)

!             CALCULATE INDEX OF THE ARRAY FOR THE POTENTIAL, BOUND STATE, AND SCATTERING STATE
              Pot_index2_Rn=nint(Rn2s/StepSizeShort)
              Pot_index2_Rp=nint(Rp2s/StepSizeShort)
              Bound_Index_n=nint(sqrt(lr**2.0+4.0*s**2.0-4.0*lr_dot_s)/StepSizeShort)
              Bound_Index_p=nint(sqrt(lr**2.0+4.0*s**2.0+4.0*lr_dot_s)/StepSizeShort)

!             MAKE SURE THE INDICES ARE NOT LESS THAN 1
              if (Pot_Index2_Rn < 1) then
                Pot_Index2_Rn=1
              end if
              if (Pot_Index2_Rp < 1) then
                Pot_Index2_Rp=1
              end if
              if (Bound_Index_n < 1) then
                Bound_Index_n=1
              end if
              if (Bound_Index_p < 1) then
                Bound_Index_p=1
              end if          


!             CALCULATE THE CONTRIBUTION TO THE NEUTRON INTEGRAL
              cos_alpha=(Rn**2.0+2.0*Rn_dot_s)/(Rn*Rn2s)
              exponent=exp(-(Rn**2.0+Rn2s**2.0)/(NeutronBeta**2.0)+(2*Rn*Rn2s/NeutronBeta**2.0)*cos_alpha)
              potential_sum=(NeutronUnl(Pot_Index1_Rn,Pot_Index2_Rn)/(pi**(3.0/2.0)*NeutronBeta**3.0))*exponent 
              integral_n=integral_n+ur_weights(ur_index)*us_weights(us_index)*lr_weights(lr_index) &
                         *s_weights(s_index)*(lr**2.0)*(s**2.0)*Vnp(lr_index1) &
                         *DeuteronBoundWF(lr_index1)*DeuteronBoundWF(Bound_Index_n) &
                         *TempScatWF(WF_Index)*phi_weights(phi_index)*Y_Rs*potential_sum

!             CALCULATE THE CONTRIBUTION TO THE PROTON INTEGRAL
              cos_alpha=(Rp**2.0+2.0*Rp_dot_s)/(Rp*Rp2s)
              exponent=exp(-(Rp**2.0+Rp2s**2.0)/(ProtonBeta**2.0)+(2*Rp*Rp2s/ProtonBeta**2.0)*cos_alpha)
              potential_sum=(ProtonUnl(Pot_Index1_Rp,Pot_Index2_Rp)/(pi**(3.0/2.0)*ProtonBeta**3.0))*exponent 
              integral_p=integral_p+ur_weights(ur_index)*us_weights(us_index)*lr_weights(lr_index) &
                         *s_weights(s_index)*(lr**2.0)*(s**2.0)*Vnp(lr_index1) &
                         *DeuteronBoundWF(lr_index1)*DeuteronBoundWF(Bound_Index_p) &
                         *TempScatWF(WF_Index)*phi_weights(phi_index)*Y_Rs*potential_sum



!              write(*,*) 'DeuteronBoundWF(lr_index1) = ', DeuteronBoundWF(lr_index1)
!              write(*,*) 'DeuteronBoundWF(Bound_Index_p) = ', DeuteronBoundWF(Bound_Index_p)
!              write(*,*) 'TempScatWF = ', TempScatWF(WF_Index)
!              write(*,*) 'Y_Rs = ', Y_Rs
!              write(*,*) 'potential_sum =', potential_sum

!!             CALCULATE THE CONTRIBUTION TO THE NEUTRON INTEGRAL
!              cos_alpha=(Rn**2.0+2.0*Rn_dot_s)/(Rn*Rn2s)
!              exponent=exp(-(Rn**2.0+Rn2s**2.0)/(beta**2.0)+(2*Rn*Rn2s/beta**2.0)*cos_alpha)
!              potential_sum=(vnl(Pot_Index1_Rn,Pot_Index2_Rn)/(pi**(3.0/2.0)*beta**3.0))*exponent 
!              integral_n=integral_n+ur_weights(ur_index)*us_weights(us_index)*lr_weights(lr_index) &
!                         *s_weights(s_index)*(lr**2.0)*(s**2.0)*Vnp(lr_index1) &
!                         *DeuteronBoundWF(lr_index1)*DeuteronBoundWF(Bound_Index_n) &
!                         *TempScatWF(WF_Index)*phi_weights(phi_index)*Y_Rs*potential_sum
!  
!!             CALCULATE THE CONTRIBUTION TO THE PROTON INTEGRAL
!              cos_alpha=(Rp**2.0+2.0*Rp_dot_s)/(Rp*Rp2s)
!              exponent=exp(-(Rp**2.0+Rp2s**2.0)/(beta**2.0)+(2*Rp*Rp2s/beta**2.0)*cos_alpha)
!              potential_sum=(vnl(Pot_Index1_Rp,Pot_Index2_Rp)/(pi**(3.0/2.0)*beta**3.0))*exponent 
!              integral_p=integral_p+ur_weights(ur_index)*us_weights(us_index)*lr_weights(lr_index) &
!                         *s_weights(s_index)*(lr**2.0)*(s**2.0)*Vnp(lr_index1) &
!                         *DeuteronBoundWF(lr_index1)*DeuteronBoundWF(Bound_Index_p) &
!                         *TempScatWF(WF_Index)*phi_weights(phi_index)*Y_Rs*potential_sum




!           PHI_s LOOP
            end do

!         us=COS(THETA_s) LOOP
          end do

!       ur=COS(THETA_r) LOOP
        end do

!     s=R'-R LOOP
      end do
     
!   lr LOOP
    end do

    factor=8.0/sqrt(2.0*L+1.0)
   
    NLint(R_index)=(integral_n*(-R*factor)*pi**(0.5)+integral_p*(-R*factor)*pi**(0.5))/abs(integral)

    Source(R_index)=-SchEqConst*(NLint(R_index) &
                    +LocalPartOfNL(R_WF_index)*R*TempScatWF(R_WF_index)-Vpot(R_index)*R*TempScatWF(R_WF_index))

! R LOOP
  end do

  return
  end subroutine
