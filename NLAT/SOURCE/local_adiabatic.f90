  subroutine local_adiabatic(DeuteronBoundWF,nmax,Upot1,Upot2,StepSize,Rmax)
  implicit none

  integer::n,nmax,lr_npoints,u_npoints,R_index_max,lr_index_max,R_index,lr_index,lr_index1
  integer::ur_index,Pot_Index1_Rn,Pot_Index1_Rp,WF_Index
  real(8)::DeuteronBoundWF(nmax),Vnp(nmax),pi,R,ur,theta_r
  real(8)::Rmax,lr_max,integral,lr,maximum,minimum
  real(8)::u_points(200),u_weights(200),lr_points(200),lr_weights(200)
  real(8)::lr_dot_R,Rn,Rp,StepSize
  complex*16::i,vpot(nmax),potential(nmax),Upot1(nmax),Upot2(nmax)
  complex*16::nAdiabatic(nmax),pAdiabatic(nmax)

  pi=acos(-1.0)
  i=cmplx(0.0,1.0)

  R_index_max=nint(Rmax/StepSize)

! GET FROM INPUTFILE
  lr_max=5.0
  lr_npoints=20
  u_npoints=20

! GET FROM SUBROUTINE
  do n=1,nmax
    Vnp(n)=-71.85*exp(-(n*StepSize/1.494)**2.0)
  end do   

! NORMALIZE DEUTERON WAVE FUNCTION WITH Vnp
! <PHI|Vnp|PHI>=-1
  integral=0.0
  do n=1,nmax
    lr=n*StepSize
    integral=integral+StepSize*Vnp(n)*(DeuteronBoundWF(n)**2.0)
  end do
!  write(*,*) 'integral in local adiabatic = ', integral
!  DeuteronBoundWF(:)=DeuteronBoundWF(:)/(sqrt(abs(integral)))

  write(*,*) 'integral=',integral

! CONSTRUCT GRID FOR ANGULAR INTEGRATION
  minimum=-1.0
  maximum=1.0
  call gauss(minimum,maximum,u_npoints,u_points,u_weights)

  minimum=StepSize
  maximum=lr_max
  call gauss(minimum,maximum,lr_npoints,lr_points,lr_weights)

! LOOP THROUGH R, DISTANCE BETWEEN TARGET AND DEUTERON
  do R_index=1,R_index_max
    R=R_index*StepSize

    nAdiabatic(R_index)=0.0
    pAdiabatic(R_index)=0.0

!   LOOP THROUGH lr: DISTANCE BETWEEN PROTON AND NEUTRON
    do lr_index=1,lr_npoints
      lr=lr_points(lr_index)
      lr_index1=nint(lr/StepSize)

!     LOOP THROUGH cos(theta_r)=ur: ANGLE BETWEEN lr AND R
      do ur_index=1,u_npoints
        ur=u_points(ur_index)
        theta_r=acos(ur)

!       CALCULATE DOT PRODUCTS
        lr_dot_R=lr*R*ur

!       CALCULATE MAGNITUDE OF THE VECTORS
        Rn=sqrt(lr**2.0/(4.0)+R**2.0-lr_dot_R)
        Rp=sqrt(lr**2.0/(4.0)+R**2.0+lr_dot_R)

!       CALCULATE INDEX OF THE ARRAY FOR THE POTENTIAL, BOUND STATE, AND SCATTERING STATE
        Pot_Index1_Rn=nint(Rn/StepSize)
        Pot_Index1_Rp=nint(Rp/StepSize)
        WF_Index=nint(R/StepSize)

!       MAKE SURE THE INDICES ARE NOT LESS THAN 1
        if (Pot_Index1_Rn <1) then
          Pot_Index1_Rn=1
        end if
        if (Pot_Index1_Rp <1) then
          Pot_Index1_Rp=1
        end if        
        if (WF_Index<1) then
          WF_Index=1
        end if
        if (Pot_Index1_Rn > R_index_max) then
          Pot_Index1_Rn = R_index_max
        end if
        if (Pot_Index1_Rp > R_index_max) then
          Pot_Index1_Rp = R_index_max
        end if


!       CALCULATE THE CONTRIBUTION TO THE NEUTRON AND PROTON INTEGRAL
        nAdiabatic(R_index)=nAdiabatic(R_index)+(u_weights(ur_index)*lr_weights(lr_index) &
                            *Vnp(lr_index1)*(DeuteronBoundWF(lr_index1)**2.0)*Upot1(Pot_Index1_Rn))

        pAdiabatic(R_index)=pAdiabatic(R_index)+(u_weights(ur_index)*lr_weights(lr_index) &
                            *Vnp(lr_index1)*(DeuteronBoundWF(lr_index1)**2.0)*Upot2(Pot_Index1_Rp))

!     cos(theta_r)=ur LOOP
      end do

!   lr LOOP
    end do

! R LOOP
  end do

  Upot1(:)=(nAdiabatic(:)/(2*integral))
  Upot2(:)=(pAdiabatic(:)/(2*integral))


  return
  end subroutine local_adiabatic
