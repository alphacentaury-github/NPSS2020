$LARGE
C==============================================================C
C   CALCULATION OF THE RGM-KERNEL FOR THE N - ALPHA - SYSTEM   C
C==============================================================C
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,S,M,N,PO,JL,FRAG,NMAX,STSANZ
      PARAMETER ( NMAX = 140 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      H2M   = 20.735D0
      PI    = 3.141592653589793D0
      ELEML = DSQRT(1.4409D0)
C---------------------------------------C
C READ IN THE VARIOUS PARAMETERS        C
C---------------------------------------C
      OPEN(5,FILE='NALPHA.INP')
C
      READ(5,*)    FRAG
      READ(5,100)
      READ(5,*)    PO
      READ(5,100)
      READ(5,*)    J,L,S
                   JL=J-2*L
                   IF((JL.NE.1).AND.(JL.NE.-1)) THEN
                            WRITE(*,*) 'VALUE FOR J IS NOT CORRECT'
                                               STOP
                   ENDIF
      READ(5,*)
      READ(5,*)    MM
      READ(5,*)    VO,LAM
      READ(5,*)    HSTEP,STSANZ
                    IF(STSANZ.GT.NMAX) THEN
                        WRITE(*,*) 'NUMBER OF MESHPOINTS IS TOO LARGE'
                        STOP
                    ENDIF
C
      CLOSE(5)
C
  100 FORMAT(1X,//)
C
      IF(PO.EQ.1) THEN
              V0(1) = 29.050D0
              V0(2) = 66.920D0
              V0(3) =  0.000D0
              BK(1) =  0.292D0
              BK(2) =  0.415D0
              BK(3) =  0.500D0
C
              PW(1) = 0.25D0*MM
              PM(1) = 0.25D0*(2.D0-MM)
              PB(1) =-0.25D0*MM
              PH(1) =-0.25D0*(2.D0-MM)
C
              PW(2) = 0.25D0*MM
              PM(2) = 0.25D0*(2.D0-MM)
              PB(2) = 0.25D0*MM
              PH(2) = 0.25D0*(2.D0-MM)
C
              PW(3) = 0.D0
              PM(3) = 0.D0
              PB(3) = 0.D0
              PH(3) = 0.D0
C
              A       =   0.514D0
C
                  ELSE
         IF(PO.EQ.2) THEN
                YY  =  0.94D0
              V0(1) = 72.98D0
              V0(2) = 72.98D0
              V0(3) =  0.D0
              BK(1) =  0.46D0
              BK(2) =  0.46D0
              BK(3) =  0.46D0
C
              PW(1) =             YY*(1.D0+MM)/4.D0
              PM(1) =             YY*(1.D0+MM)/4.D0
              PB(1) =             YY*(1.D0-MM)/4.D0
              PH(1) =             YY*(1.D0-MM)/4.D0
C
              PW(2) = (1.D0-YY)*(1.D0-3.D0*MM)/6.D0
              PM(2) = (1.D0-YY)*(1.D0+3.D0*MM)/3.D0
              PB(2) = (1.D0-YY)*(1.D0+3.D0*MM)/6.D0
              PH(2) = (1.D0-YY)*(1.D0-3.D0*MM)/3.D0
C
              PW(3) = 0.D0
              PM(3) = 0.D0
              PB(3) = 0.D0
              PH(3) = 0.D0
C
              A       = 0.540D0
C
                     ELSE
              V0(1) =    6.000D0
              V0(2) =  546.000D0
              V0(3) =-1655.000D0
              BK(1) =    0.160D0
              BK(2) =    1.127D0
              BK(3) =    3.400D0
C
              PW(1) =-0.2361D0
              PM(1) = 1.1528D0
              PB(1) = 0.5972D0
              PH(1) =-0.5139D0
C
              PW(2) = 0.8297D0-MM
              PM(2) =      MM
              PB(2) = 0.1401D0
              PH(2) = 0.0302D0
C
              PW(3) = 0.4474D0
              PM(3) = 0.3985D0
              PB(3) = 0.1015D0
              PH(3) = 0.0526D0
C
              A       =     0.580D0
C
         ENDIF
      ENDIF
C------------------------------------ C
C  CALCULATION OF THE MATRIX ELEMENTS C
C------------------------------------ C
      CALL INITIA
      CALL BIND(EALPHA)
C
      WRITE(*,*) 'CALCULATION OF THE KINETIC ENERGY KERNEL'
      CALL   TKIN(L)

      WRITE(*,*) 'CALCULATION OF THE KERNEL FOR CENTRAL FORCE'
        CALL   ZENTRA(L)

      WRITE(*,*) 'CALCULATION OF THE COULOMB KERNEL'
          CALL   COUPOT(L,FRAG)

      WRITE(*,*) 'CALCULATION OF THE LS-KERNEL'
            CALL   SPIORB(JL,L)

      WRITE(*,*) 'CALCULATION OF THE NORM KERNEL'
              CALL   NORM(L)
C
      CALL WRITIN(EALPHA)
C
      STOP
      END
C*************************************C
C  INITIALIZE THE MATRIX MAT          C
C*************************************C
      SUBROUTINE INITIA
      INTEGER*4 M,N,NMAX
      PARAMETER ( NMAX = 140 )
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /MATRIX/ MAT
C
          DO 10 N=1,NMAX
           DO 10 M=1,NMAX
                   MAT(M,N) = 0.D0
   10     CONTINUE
      RETURN
      END
C***************************************C
C  BINDING ENERGY OF THE ALPHA-PARTICLE C
C***************************************C
      SUBROUTINE BIND(EALPHA)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
C
        KONST       = 0.D0
      DO 10 K=1,3
        KONSTS      = DSQRT(A/(A+2.D0*BK(K)))
        KONSTS      = KONSTS*KONSTS*KONSTS
        KONST       = KONST+6.D0*(PW(K)+PM(K))*KONSTS*V0(K)
   10 CONTINUE
        KONSTS      = ELEML*ELEML*DSQRT((2.D0*A)/PI)
        EALPHA      = (9.D0/2.D0)*H2M*A-KONST+KONSTS
C
      RETURN
      END
C*****************C
C  NORM KERNEL    C
C*****************C
      SUBROUTINE NORM(L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,NMAX,STSANZ
      PARAMETER ( NMAX = 140 )
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      DO 10 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 10 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        KONST       = (4.D0/5.D0)*DSQRT((4.D0*A)/(3.D0*PI))
        EG1         = KONST*KONST*KONST
        AG1         = (34.D0/75.D0)*A
        BG1         = AG1
        CG1         = (32.D0/75.D0)*A
C
        MAT(N+1,M)  =-EG1*S(AG1,BG1,CG1,L,X,Y)
   10 CONTINUE
C
      RETURN
      END
C*********************************C
C  KINETIC ENERGY KERNEL          C
C*********************************C
      SUBROUTINE TKIN(L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,NMAX,STSANZ
      PARAMETER ( NMAX = 140 )
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      DO 10 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 10 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        KONST       = (4.D0/5.D0)*DSQRT((4.D0*A)/(3.D0*PI))
        EG1         = KONST*KONST*KONST
        AG1         = (34.D0/75.D0)*A
        BG1         = AG1
        CG1         = (32.D0/75.D0)*A
        DG1         = (47.D0/5.D0)*A
        FG1         =-(1216.D0/1125.D0)*A*A
        FG1S        = FG1
        GG1         =-(1568.D0/1125.D0)*A*A
        TERM        =-H2M*((DG1+FG1*X*X+FG1S*Y*Y)*S(AG1,BG1,CG1,L,X,Y)+
     1                GG1*X*Y*T(AG1,BG1,CG1,L,X,Y))
C--------------------------------------------------C
C  NON-LOCAL PART OF KINETIC ENERGY KERNELS        C
C--------------------------------------------------C
        MAT(M,N+1)  = MAT(M,N+1)+EG1*TERM
   10 CONTINUE
C
      RETURN
      END
C***************************C
C  KERNEL FOR CENTRAL FORCE C
C***************************C
      SUBROUTINE ZENTRA(L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,NMAX,STSANZ
      PARAMETER ( NMAX = 140 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C------------------------------------------------C
C  DIRECT NUCLEAR POTENTIAL                      C
C------------------------------------------------C
        IN          = STSANZ+1
      DO 10 M=1,IN
        X           = DBLE(M)*HSTEP
       DO 10 K=1,3
        KONST       = DSQRT(4.D0*A/(4.D0*A+3.D0*BK(K)))
        KONST       = KONST*KONST*KONST
        W15G0       = 4.D0*PW(K)-PM(K)+2.D0*PB(K)-2.D0*PH(K)
        MAT(M,M)    = MAT(M,M)-V0(K)*W15G0*KONST*
     1                   DEXP(-((4.D0*A*BK(K))/(4.D0*A+3.D0*BK(K)))*X*X)
   10 CONTINUE
C
      DO 11 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 11 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        KONST       = (4.D0/5.D0)*DSQRT((4.D0*A)/(3.D0*PI))
        EG1         = KONST*KONST*KONST
        AG1         = (34.D0/75.D0)*A
        BG1         = AG1
        CG1         = (32.D0/75.D0)*A
        TERM        = 0.D0
C
         DO 12  K=1,3
          W12G1       =-3.D0*(PW(K)+PM(K))
          KONST       = DSQRT((3.D0*A)/(3.D0*A+2.D0*BK(K)))
          E12G1       = KONST*KONST*KONST
          A12G1       = (34.D0*A*A+ 28.D0*A*BK(K))/(75.D0*A+50.D0*BK(K))
          B12G1       = (34.D0*A*A+108.D0*A*BK(K))/(75.D0*A+50.D0*BK(K))
          C12G1       = (32.D0*A*A+ 64.D0*A*BK(K))/(75.D0*A+50.D0*BK(K))
C
          W25G1       = W12G1
          E25G1       = E12G1
          A25G1       = B12G1
          B25G1       = A12G1
          C25G1       = C12G1
C
          W15G1       =-PW(K)+4.D0*PM(K)-2.D0*PB(K)+2.D0*PH(K)
          E15G1       = 1.D0
          A15G1       = (34.D0/75.D0)*A+(16.D0/25.D0)*BK(K)
          B15G1       = A15G1
          C15G1       = (32.D0/75.D0)*A-(32.D0/25.D0)*BK(K)
C
          W23G1       =-3.D0*(PW(K)+PM(K))
          KONST       = DSQRT(A/(A+2.D0*BK(K)))
          E23G1       = KONST*KONST*KONST
          A23G1       = AG1
          B23G1       = BG1
          C23G1       = CG1
C
          KONST1      = W12G1*E12G1*S(A12G1,B12G1,C12G1,L,X,Y)
          KONST2      = W25G1*E25G1*S(A25G1,B25G1,C25G1,L,X,Y)
          KONST3      = W15G1*E15G1*S(A15G1,B15G1,C15G1,L,X,Y)
          KONST4      = W23G1*E23G1*S(A23G1,B23G1,C23G1,L,X,Y)
          TERM        = TERM-V0(K)*(KONST1+KONST2+KONST3+KONST4)
   12    CONTINUE
C---------------------------------------------C
C  NON-LOCAL PART OF CENTRAL POTENTIAL KERNEL C
C---------------------------------------------C
          MAT(M,N+1)  = MAT(M,N+1)+EG1*TERM
   11 CONTINUE
C
      RETURN
      END
C****************C
C  COULOMBKERNEL C
C****************C
      SUBROUTINE COUPOT(L,FRAG)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,NMAX,FRAG,SKIP,STSANZ,IFAIL
      PARAMETER ( NMAX = 140 , SKIP = 1)
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C----------------------------------------------C
C  DIRECT COULOMB POTENTIAL IN P+ALPHA CHANNEL C
C----------------------------------------------C
      IF(FRAG.EQ.1) THEN
       IN          = STSANZ+1
       DO 10 M=1,IN
        X           = DBLE(M)*HSTEP
        KONST       = DSQRT((4.D0*A)/3.D0)
        MAT(M,M)    = MAT(M,M) +
     1                      ((2.D0*ELEML*ELEML)/X)*S15AEF(KONST*X,IFAIL)
   10  CONTINUE
      ENDIF
C---------------------------------------C
C FOR SKIP=1, CALCULATION OF THE        C
C NON-LOCAL COULOMB POTENTIAL           C
C---------------------------------------C
      IF(SKIP.EQ.0) RETURN
C
      DO 11 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 11 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        KONST       = (4.D0/5.D0)*DSQRT((4.D0*A)/(3.D0*PI))
        EG1         = KONST*KONST*KONST
        AG1         = (34.D0/75.D0)*A
        BG1         = AG1
        CG1         = (32.D0/75.D0)*A
C--------------------------------------------------C
C  NON-LOCAL COULOMB KERNEL IN N+ALPHA CHANNEL     C
C--------------------------------------------------C
        IF(FRAG.EQ.1) THEN
            KONST       = 2.D0*PI*ELEML*ELEML
            LI          = DSQRT((3.D0*A)/2.D0)*(16.D0/15.D0)
            MI          = DSQRT((3.D0*A)/2.D0)*(4.D0/15.D0)
            VOR         = DSQRT((3.D0*A)/2.D0)
            TERM        = (5.D0/4.D0)*COULI(AG1,BG1,1.D0,-1.D0,CG1,L,X,
     1                    Y,2)+VOR*COULI(AG1,BG1,LI,MI,CG1,L,X,Y,1)+VOR
     2                    *COULI(AG1,BG1,MI,LI,CG1,L,X,Y,1)
            TERM        =-EG1*KONST*TERM
                       ELSE
            KONST       = ELEML*ELEML*DSQRT((2.D0*A)/PI)
            TERM        =-EG1*KONST*S(AG1,BG1,CG1,L,X,Y)
        ENDIF
C
          MAT(M,N+1)  = MAT(M,N+1)+TERM
   11 CONTINUE
C
      RETURN
      END
C**************C
C  LS - KERNEL C
C**************C
      SUBROUTINE SPIORB(JL,L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,JL,NMAX,STSANZ
      PARAMETER ( NMAX = 140 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      IF((VO.EQ.0.D0).OR.(L.EQ.0)) RETURN
C
C--------------------------------------------------C
C  DIRECT SPIN-_ORBIT POTENTIAL IN N+ALPHA CHANNEL  C
C--------------------------------------------------C
      IN          = STSANZ+1
      DO 10 M=1,IN
         X        = DBLE(M)*HSTEP
         KONST1   = (4.D0*A)/(4.D0*A+3.D0*LAM)
         KONST    = DSQRT(KONST1)
         KONST    = KONST*KONST*KONST*KONST*KONST
         MAT(M,M) = MAT(M,M)-(5.D0/2.D0)*VO*JLS(JL,L)*KONST*
     1                                           DEXP(-(KONST1*LAM)*X*X)
   10 CONTINUE
C
      DO 11 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 11 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
            KONST = (4.D0/5.D0)*DSQRT((4.D0*A)/(3.D0*PI))
            EG1   =  KONST*KONST*KONST
            SA    = (34.D0/75.D0)*A+(48.D0/75.D0)*LAM
            SB    =  SA
            SC    = (32.D0/75.D0)*A-(32.D0/25.D0)*LAM
            TERM  =  VO*EG1*((5.D0*A)/(A-3.D0*LAM))*S(SA,SB,SC,L,X,Y)
C-------------------------------------------------C
C  NON-LOCAL SPIN-ORBIT POTENTIAL                 C
C-------------------------------------------------C
          MAT(M,N+1)  = MAT(M,N+1)+JLS(JL,L)*TERM
   11 CONTINUE
C
      RETURN
      END
C**************************************C
C  WRITE MATRIX MAT ON A FILE          C
C**************************************C
      SUBROUTINE WRITIN(EALPHA)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,NMAX,STSANZ
      PARAMETER ( NMAX = 140 )
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      OPEN(20,FILE='NALPHA.OUT',STATUS='NEW',FORM='UNFORMATTED')
C
      WRITE(20) EALPHA
C
      IN = STSANZ+1
      DO 10 M=1,IN
         WRITE(20) MAT(M,M)
   10 CONTINUE
C
      DO 11 M=1,STSANZ
       DO 11 N=M,STSANZ
         WRITE(20) MAT(M,N+1),MAT(N+1,M)
   11 CONTINUE
      CLOSE(20)
C
      RETURN
      END
C---------------------------------------------------------C
C     J=1(J VON P-A IST L+1/2) J=-1(J VON P-A IST L-1/2)  C
C---------------------------------------------------------C
      DOUBLE PRECISION FUNCTION JLS(JL,L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 L,I,JL
C
         TERM = DBLE(L)
         I = JL*L
         IF(I.LT.0) TERM = -DBLE(L+1)
         JLS = TERM
C
      RETURN
      END
