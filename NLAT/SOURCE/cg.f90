!***********************************************************************************************************************************
!
!                                                               C G
!
!  Program:      CG
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         April 20, 2005
!
!  Language:     Fortran-90
!
!  Version:      1.00a
!
!  Description:  Calculates Clebsch-Gordan coefficients.
!
!  Files:        Source files:
!
!                   cg.f90                   Main program
!
!  Notes:
!
!***********************************************************************************************************************************

      Subroutine CG(j1,m1,j2,m2,j3,m3,c)

      IMPLICIT NONE

      INTEGER :: I, K
      DOUBLE PRECISION :: J1, J2, J3, M1, M2, M3, C, SUMK, TERM
!      DOUBLE PRECISION:: SUMK,TERM
      DOUBLE PRECISION, DIMENSION(0:99) :: FACT

      LOGICAL :: ISFRAC

!-----------------------------------------------------------------------------------------------------------------------------------

!
!     Compute table of factorials.
!

      FACT(0) = 1.0D0

      DO I = 1, 99
         FACT(I) = I * FACT(I-1)
      END DO

!
!     Get user input.
!

!      WRITE (UNIT=*, FMT='(/A)', ADVANCE='NO') '  Enter j1:  '
!      READ (UNIT=*, FMT=*) J1

!      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter j2:  '
!      READ (UNIT=*, FMT=*) J2

!      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter j3:  '
!      READ (UNIT=*, FMT=*) J3

!      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter m1:  '
!      READ (UNIT=*, FMT=*) M1

!      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter m2:  '
!      READ (UNIT=*, FMT=*) M2

!      M3 = M1 + M2

!
!     Check for invalid input.
!

      IF (ISFRAC(J1+J2+J3) .OR. ISFRAC(J1+M1)     .OR. ISFRAC(J2+M2) .OR.  &
          ISFRAC(J3+M3)    .OR. ISFRAC(-J1+J3-M2) .OR. ISFRAC(-J2+J3+M1)) THEN
         WRITE (UNIT=*, FMT='(/A)') ' Invalid input.'
         STOP
      END IF

!
!     Check for conditions that give C = 0.
!

      IF ( (J3 .LT. ABS(J1-J2)) .OR.  &
           (J3 .GT. (J1+J2))    .OR.  &
           (ABS(M1) .GT. J1)    .OR.  &
           (ABS(M2) .GT. J2)    .OR.  &
           (ABS(M3) .GT. J3)) THEN
         C = 0.0D0
      ELSE

!
!     Compute Clebsch-Gordan coefficient.
!

         C = SQRT((J3+J3+1)/FACT(NINT(J1+J2+J3+1)))
         C = C * SQRT(FACT(NINT(J1+J2-J3))*FACT(NINT(J2+J3-J1))*FACT(NINT(J3+J1-J2)))
         C = C * SQRT(FACT(NINT(J1+M1))*FACT(NINT(J1-M1))*FACT(NINT(J2+M2))*FACT(NINT(J2-M2))*FACT(NINT(J3+M3))*FACT(NINT(J3-M3)))
         SUMK = 0.0D0
         DO K = 0, 99
            IF (J1+J2-J3-K .LT. 0.0D0) CYCLE
            IF (J3-J1-M2+K .LT. 0.0D0) CYCLE
            IF (J3-J2+M1+K .LT. 0.0D0) CYCLE
            IF (J1-M1-K    .LT. 0.0D0) CYCLE
            IF (J2+M2-K    .LT. 0.0D0) CYCLE
            TERM = FACT(NINT(J1+J2-J3-K))*FACT(NINT(J3-J1-M2+K))*FACT(NINT(J3-J2+M1+K))*FACT(NINT(J1-M1-K))*  &
               FACT(NINT(J2+M2-K))*FACT(K)
            IF (MOD(K,2) .EQ. 1) TERM = -TERM
            SUMK = SUMK + 1.0D0/TERM
         END DO
         C = C * SUMK
      END IF

!
!     Print result.
!

!      WRITE (UNIT=*, FMT='(/1X,F15.10)') C

!      STOP

      RETURN
 
      END subroutine CG





!***********************************************************************************************************************************
!  ISFRAC
!
!  Return .TRUE. if the argument has a fractional part, and .FALSE. if the argument is an integer.
!***********************************************************************************************************************************

      FUNCTION ISFRAC (X) RESULT (Y)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: X
      LOGICAL :: Y
      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-8


      IF ((ABS(X)-INT(ABS(X))) .GT. EPS) THEN
         Y = .TRUE.
      ELSE
         Y = .FALSE.
      END IF

      RETURN

      END FUNCTION ISFRAC
