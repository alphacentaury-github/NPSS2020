$LARGE
C================================================C
C   THE PROGRAM CALCULATES THE RGM KERNELS       C
C   FOR THE ALPHA + ALPHA SYSTEM                 C
C================================================C
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,S,M,N,PO,NMAX,STSANZ
      PARAMETER ( NMAX = 240 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      H2M   = 20.735D0
      PI    = 3.141592653589793D0
      ELEML = DSQRT(1.4409D0)
C---------------------------------------C
C READ IN THE VARIOUS PARAMETERS        C
C---------------------------------------C
      OPEN(5,FILE='ALPHALPH.INP')
C
      READ(5,100)
      READ(5,*)    PO
      READ(5,200)
      READ(5,*)    J,L,S
      READ(5,*)
      READ(5,*)    MM
      READ(5,*)    HSTEP,STSANZ
                    IF(STSANZ.GT.NMAX) THEN
                        WRITE(*,*) 'NUMBER OF MESHPOINTS IS TOO LARGE'
                        STOP
                    ENDIF
C
      CLOSE(5)
C
  100 FORMAT(1X,/)
  200 FORMAT(1X,//)
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
C--------------------------------------C
C  CALCULATION OF THE MATRIXELEMENTS   C
C--------------------------------------C
      CALL INITIA
      CALL BIND(EALPHA)
C
      IF(MOD(L,2).EQ.0) THEN
      WRITE(*,*) 'CALCULATION OF THE KINETIC ENERGY KERNEL'
                      CALL   TKIN(L)

      WRITE(*,*) 'CALCULATION OF THE KERNEL FOR CENTRAL FORCE'
                        CALL   ZENTRA(L)

      WRITE(*,*) 'CALCULATION OF THE COULOMB KERNEL'
                          CALL   COUPOT(L)

      WRITE(*,*) 'CALCULATION OF THE NORM KERNEL'
                            CALL   NORM(L)
      ENDIF
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
      PARAMETER ( NMAX = 240 )
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
C****************************************C
C  BINDING ENERGY OF THE ALPHA PARTICLE  C
C****************************************C
      SUBROUTINE BIND(EALPHA)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK
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
      PARAMETER ( NMAX = 240 )
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
        KONST       = DSQRT((8.D0*A)/(3.D0*PI))
        EG1         = KONST*KONST*KONST
        KONST       = DSQRT((2.D0*A)/PI)
        EG2         = 12.D0*PI*KONST*KONST*KONST
        AG1         = (5.D0/3.D0)*A
        CG1         =-(8.D0/3.D0)*A
C
        MAT(N+1,M)  =-4.D0*EG1*S(AG1,AG1,CG1,L,X,Y)
C
        IF(L.EQ.0) MAT(N+1,M) = MAT(N+1,M)+EG2*X*Y*DEXP(-A*(X*X+Y*Y))
   10 CONTINUE
C
      RETURN
      END
C*********************************C
C  KINETIC ENERGY KERNELS         C
C*********************************C
      SUBROUTINE TKIN(L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,NMAX,STSANZ
      PARAMETER ( NMAX = 240 )
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
        KONST       = DSQRT((8.D0*A)/(3.D0*PI))
        EG1         = KONST*KONST*KONST
        AG1         = (5.D0/3.D0)*A
        CG1         =-(8.D0/3.D0)*A
        DG1         = 62.D0*A
        FG1         =-(224.D0/9.D0)*A*A
        GG1         = (416.D0/9.D0)*A*A
C
        KONST       = DSQRT((2.D0*A)/PI)
        EG2         = 4.D0*PI*KONST*KONST*KONST
        DG2         =-(81.D0/2.D0)*A
        FG2         = 6.D0*A*A
C--------------------------------------------------C
C  NON-LOCAL KINETIC ENERGY KERNELS                C
C--------------------------------------------------C
        TERM        =-H2M*((DG1+FG1*X*X+FG1*Y*Y)*S(AG1,AG1,CG1,L,X,Y)+
     1                GG1*X*Y*T(AG1,AG1,CG1,L,X,Y))
        MAT(M,N+1)  = MAT(M,N+1)+EG1*TERM
C
        IF(L.EQ.0) THEN
        TERMS       =-H2M*(DG2+FG2*X*X+FG2*Y*Y)*X*Y*DEXP(-A*(X*X+Y*Y))
        MAT(M,N+1)  = MAT(M,N+1)+EG2*TERMS
        ENDIF
C
   10 CONTINUE
C
      RETURN
      END
C***********************************************C
C  CONTRIBUTIONS FROM THE NN CENTRAL POTENTIAL  C
C***********************************************C
      SUBROUTINE ZENTRA(L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,NMAX,STSANZ
      PARAMETER ( NMAX = 240 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C---------------C
C  LOCAL PART   C
C---------------C
        IN          = STSANZ+1
      DO 10 M=1,IN
        X           = DBLE(M)*HSTEP
       DO 10 K=1,3
        KONST       = DSQRT(2.D0*A/(2.D0*A+3.D0*BK(K)))
        KONST       = KONST*KONST*KONST
        W15G0       = 4.D0*(4.D0*PW(K)-PM(K)+2.D0*PB(K)-2.D0*PH(K))
        MAT(M,M)    = MAT(M,M)-V0(K)*W15G0*KONST*
     1                   DEXP(-((2.D0*A*BK(K))/(2.D0*A+3.D0*BK(K)))*X*X)
   10 CONTINUE
C
      E1 = DSQRT(A*A/PI)
      E1 = DSQRT(2.D0)*E1*E1*E1

      DO 11 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 11 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        TERM = 0.D0
        DO 12 K=1,3
C
         EG1 =-2.D0**7*E1/DSQRT(3.D0*(A+2.D0*BK(K))**3)
         EG2 = 2.D0**6*E1/DSQRT((3.D0*A)**3)
         EG3 = 3.D0*2.D0**6*E1/DSQRT((3.D0*A+4.D0*BK(K))**3)
         EG4 =-3.D0*2.D0**7*E1/DSQRT((3.D0*A+4.D0*BK(K))**3)
         EG5 = 3.D0*2.D0**3*E1/DSQRT((A+BK(K))**3)
C
         AG1 = (5.D0/3.D0)*A
         CG1 = (8.D0/3.D0)*A
         AG2 = (5.D0*A+12.D0*BK(K))/3.D0
         CG2 = (8.D0*A+24.D0*BK(K))/3.D0
         AG3 = (5.D0*A*A+8.D0*A*BK(K))/(3.D0*A+4.D0*BK(K))
         CG3 =  8.D0*A*(A+BK(K))/(3.D0*A+4.D0*BK(K))
         AG4 =  AG3
         BG4 = (5.D0*A*A+12.D0*A*BK(K))/(3.D0*A+4.D0*BK(K))
         CG4 =  8.D0*A*(A+2.D0*BK(K))/(3.D0*A+4.D0*BK(K))
         AG5 =  A*(A+2.D0*BK(K))/(A+BK(K))
         CG5 =  2.D0*A*BK(K)/(A+BK(K))
C
        WAA1 = 1.D0*PW(K)+1.D0*PM(K)
        WAA2 =-1.D0*PW(K)+4.D0*PM(K)-2.D0*PB(K)+2.D0*PH(K)
        WAA3 =-3.D0*PW(K)+2.D0*PM(K)-2.D0*PB(K)+2.D0*PH(K)
        WAA4 = 1.D0*PW(K)+1.D0*PM(K)
        WAA5 = 2.D0*PW(K)-3.D0*PM(K)+2.D0*PB(K)-2.D0*PH(K)
C
          KONST1      = WAA1*EG1*S(AG1,AG1,CG1,L,X,Y)
          KONST2      = WAA2*EG2*S(AG2,AG2,CG2,L,X,Y)
          KONST3      = WAA3*EG3*S(AG3,AG3,CG3,L,X,Y)
          KONST4      = WAA4*EG4*S(AG4,BG4,CG4,L,X,Y)
          KONST5      = WAA4*EG4*S(BG4,AG4,CG4,L,X,Y)
          KONST6      = WAA5*EG5*S(AG5,AG5,CG5,L,X,Y)
          TERM        = TERM-V0(K)*(KONST1+KONST2+KONST3+KONST4+KONST5+
     1                              KONST6)
C
        IF(L.EQ.0) THEN
C
         EG6 = 12.D0*2.D0**6*PI*E1/DSQRT(2.D0*(2.D0*A+3.D0*BK(K))**3)
         EG7 = 12.D0*2.D0**3*PI*E1/DSQRT((A+2.D0*BK(K))**3)
C
         AG6 = A*(2.D0*A+5.D0*BK(K))/(2.D0*A+3.D0*BK(K))
         BG6 = A
C
          KONST7      = EG6*WAA1*X*Y*DEXP(-(AG6*X*X+BG6*Y*Y))
          KONST8      = EG6*WAA1*X*Y*DEXP(-(BG6*X*X+AG6*Y*Y))
          KONST9      = EG7*WAA1*X*Y*DEXP(-A*(X*X+Y*Y))
C
          TERM        = TERM-V0(K)*(KONST7+KONST8+KONST9)
        ENDIF
C
   12    CONTINUE
C-------------------------------------------------C
C  NN CENTRAL POTENTIAL: NON-LOCAL CONTRIBUTIONS  C
C-------------------------------------------------C
          MAT(M,N+1)  = MAT(M,N+1)+TERM
   11 CONTINUE
C
      RETURN
      END
C*******************C
C  COULOMB KERNELS  C
C*******************C
      SUBROUTINE COUPOT(L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,NMAX,SKIP,STSANZ,IFAIL
      PARAMETER ( NMAX = 240 , SKIP = 1)
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C----------------------------C
C  DIRECT COULOMB-POTENTIAL  C
C----------------------------C
      IN          = STSANZ+1
      DO 10 M=1,IN
        X           = DBLE(M)*HSTEP
        KONST       = DSQRT((2.D0*A)/3.D0)
        MAT(M,M)    = MAT(M,M) +
     1                      ((4.D0*ELEML*ELEML)/X)*S15AEF(KONST*X,IFAIL)
   10 CONTINUE
C---------------------------------------C
C FOR SKIP=1: CALCULATION OF THE        C
C NON-LOCAL COULOMB POTENTIAL           C
C---------------------------------------C
      IF(SKIP.EQ.0) RETURN
C
        E1  = 4.D0*PI*ELEML*ELEML
        EC1 =-16.D0*E1*A*A*DSQRT(2.D0*PI)/(3.D0*PI*PI)
        EC3 = 2.5D0*EC1
        EC4 =-4.D0*E1*A*A*DSQRT(PI*A)/(3*PI*PI*DSQRT(6))
        EC5 = 8.D0*E1*A*A*DSQRT(2.D0*PI)/(PI*PI)
        EC6 =-32.D0*ELEML*ELEML*A*A/(3.D0*DSQRT(3.D0)*PI*PI)
C
        AC1 = (5.D0/3.D0)*A
        LI1 = 2.D0*DSQRT(A/3.D0)
        MI1 =-DSQRT(A/3.D0)
        CI1 =-(8.D0/3.D0)*A
        LI2 = 1.D0
        MI2 =-1.D0
        AC5 = A
        LI5 = DSQRT(A)
        CI5 = 0.D0
C
      DO 11 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 11 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        KONST1 = EC1*COULI(AC1,AC1, LI1, MI1,CI1,L,X,Y,1)
        KONST2 = EC1*COULI(AC1,AC1,-MI1,-LI1,CI1,L,X,Y,1)
        KONST3 = EC3*COULI(AC1,AC1,-MI1,-MI1,CI1,L,X,Y,1)
        KONST4 = EC4*COULI(AC1,AC1, LI2, MI2,CI1,L,X,Y,2)
        KONST5 = EC5*COULI(AC5,AC5, LI5,-LI5,CI5,L,X,Y,1)
C
        KONST6 = EC6*S(AC1,AC1,CI1,L,X,Y)
C------------------------------------------------------C
C  NON-LOCAL COULOMB KERNELS IN ALPHA-ALPHA CHANNEL    C
C------------------------------------------------------C
        MAT(M,N+1)  = MAT(M,N+1)+(KONST1+KONST2+KONST3+KONST4+KONST5+
     1                            KONST6)

       IF(L.EQ.0) THEN
C
        EC7 = 32.D0*ELEML*ELEML*A*A/PI
        EC8 = 16.D0*E1*A*A*DSQRT(PI)/(DSQRT(3)*PI*PI)
        LI8 = DSQRT(2.D0*A/3.D0)
C
        KONST7 = EC7*X*Y*DEXP(-A*(X*X+Y*Y))
        KONST8 = EC8*DEXP(-A*(X*X+Y*Y))*(Y*S15AEF(LI8*X,IFAIL)/LI8+
     1                                        X*S15AEF(LI8*Y,IFAIL)/LI8)
C
        MAT(M,N+1)  = MAT(M,N+1)+KONST7+KONST8
       ENDIF
C
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
      PARAMETER ( NMAX = 240 )
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      OPEN(20,FILE='ALPHALPH.OUT',STATUS='NEW',FORM='UNFORMATTED')
C
      WRITE(20) 2.D0*EALPHA
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
