
  subroutine CoulombPhase(eta,L,phase)
  implicit none
    
  integer::L,n
  real(8)::C,delC(0:5000),diff,eta,euler,phase

!  eta=0.662663009801078
  euler=0.5772156649
!  L=3

!  write(*,*) 'L,eta=',L,eta

!  do L=0,20

    C=0.0
    if (L>0) then
      do n=1,L
        C=C+1.0/real(n)
      end do
    end if
    if (L==0) then
      C=0.0
    end if

    delc(0)=eta/(1.0+real(L))-atan(eta/(1.0+real(L)))
    do n=1,5000
      delc(n)=delc(n-1)+eta/(1.0+real(L)+real(n))-atan(eta/(1.0+real(L)+real(n)))
!      write(*,*) n,delc(n)
    end do
!    GammaFunction=delc(5000)+eta*(C-euler)
    phase=delc(5000)+eta*(C-euler)
    
!    write(*,*) 'L,phase=',L,phase
!  end do

  return
  end subroutine
