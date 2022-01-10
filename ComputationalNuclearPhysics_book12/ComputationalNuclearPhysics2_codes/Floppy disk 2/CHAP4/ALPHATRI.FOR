$LARGE
C===============================================================C
C   THIS PROGRAM CALCULATES THE RGM KERNELS FOR THE             C
C   ALPHA + HELIUM3 AND ALPHA + TRITON SYSTEMS                  C
C===============================================================C
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,S,M,N,PO,JL,FRAG,NMAX,STSANZ
      PARAMETER ( NMAX = 240 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,VT,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      H2M   = 20.735D0
      PI    = 3.141592653589793D0
      ELEML = DSQRT(1.4409D0)
C---------------------------------------C
C READ IN THE VARIOUS KERNELS           C
C---------------------------------------C
      OPEN(5,FILE='ALPHATRI.INP')
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
      READ(5,*)    VO,VT,LAM
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
                                 A  = 0.514D0
                                 G  = 0.378D0
                   IF(FRAG.EQ.2) G  = 0.367D0
C
                  ELSE
        IF(PO.EQ.2) THEN
                YY  =  0.94D0
              V0(1) = 72.98D0
              V0(2) = 72.98D0
              V0(3) =  0.00D0
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
                                 A  = 0.540D0
                                 G  = 0.690D0
                   IF(FRAG.EQ.2) G  = 0.690D0
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
                                 A  = 0.574D0
                                 G  = 0.454D0
                   IF(FRAG.EQ.2) G  = 0.454D0
        ENDIF
      ENDIF
C-------------------------------------C
C  CALCULATION OF THE MATRIXELEMENTS  C
C-------------------------------------C
      CALL INITIA
      CALL BIND(EALPHA,ETHREE)
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
      CALL WRITIN(EALPHA,ETHREE)
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
C************************************C
C  BINDING ENERGIES OF THE CLUSTERS  C
C************************************C
      SUBROUTINE BIND(EALPHA,ETHREE)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,VT,LAM
C---------------------------------------C
C INTERNAL ENERGY OF THE ALPHA-PARTICLE C
C---------------------------------------C
        KONST       = 0.D0
      DO 10 K=1,3
        KONSTS      = DSQRT(A/(A+2.D0*BK(K)))
        KONSTS      = KONSTS*KONSTS*KONSTS
        KONST       = KONST+6.D0*(PW(K)+PM(K))*KONSTS*V0(K)
   10 CONTINUE
        KONSTS      = ELEML*ELEML*DSQRT((2.D0*A)/PI)
        EALPHA      = (9.D0/2.D0)*H2M*A-KONST+KONSTS
C----------------------------------------------C
C INTERNAL ENERGY OF HELIUM3, TRITON           C
C----------------------------------------------C
         KONST1 = 0.D0
      DO 12 K=1,3
         KONST2 = DSQRT(G/(G+2.D0*BK(K)))
         KONST2 = KONST2*KONST2*KONST2
         KONST1 = KONST1+V0(K)*(3.D0*PW(K)+3.D0*PM(K))*KONST2
   12 CONTINUE
         ETHREE = 3.D0*H2M*G-KONST1+DKRDEL(FRAG,2)*(ELEML*ELEML)*
     1            DSQRT((2.D0*G)/PI)
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
      COMMON /BREITE/ A,G
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      ZAHL = 6.D0/7.D0
C
      U1  =  36.D0*A*A+97.D0*A*G+36.D0*G*G
      V1  =  72.D0*A*A+96.D0*A*G+72.D0*G*G
      L1  =  56.D0*A+63.D0*G
      U2  =   3.D0*A*A+31.D0*A*G+3.D0*G*G
      V2  =   6.D0*A*A-36.D0*A*G+6.D0*G*G
      L2  =  14.D0*A+21.D0*G
      GAM =  DSQRT(12.D0*G*G/(PI*A))
      GAM = (12.D0/7.D0)**3*GAM*GAM*GAM
C
      E1  = DSQRT(A*A/(G*(8.D0*A+9.D0*G)))
      E1  =-3.D0*GAM*E1*E1*E1
      A1  = U1/L1
      C1  = V1/L1
C
      E2  = DSQRT(A*A*A/((A+G)**2*(2.D0*A+3.D0*G)))
      E2  = 3.D0*GAM*E2*E2*E2
      A2  = U2/L2
      C2  = V2/L2
C
      E3  =-GAM*(4.D0/3.D0)**3*(A/(A+G))**6
      A3  = 25.D0*A/7.D0
      C3  =-48.D0*A/7.D0
C
      DO 10 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 10 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        TERM1 = E1*S(ZAHL*A1,ZAHL*A1,-ZAHL*C1,L,X,Y)
        TERM2 = E2*S(ZAHL*A2,ZAHL*A2,-ZAHL*C2,L,X,Y)
        TERM3 = E3*S(ZAHL*A3,ZAHL*A3,-ZAHL*C3,L,X,Y)
C
        MAT(N+1,M)  =  TERM1+TERM2+TERM3
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
      PARAMETER ( NMAX = 240 )
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      ZAHL =  6.D0/7.D0
C
      U1  =  36.D0*A*A+97.D0*A*G+36.D0*G*G
      V1  =  72.D0*A*A+96.D0*A*G+72.D0*G*G
      L1  =  56.D0*A+63.D0*G
      U2  =   3.D0*A*A+31.D0*A*G+3.D0*G*G
      V2  =   6.D0*A*A-36.D0*A*G+6.D0*G*G
      L2  =  14.D0*A+21.D0*G
      GAM =  DSQRT(12.D0*G*G/(PI*A))
      GAM = (12.D0/7.D0)**3*GAM*GAM*GAM
C
      E1  =  DSQRT(A*A/(G*(8.D0*A+9.D0*G)))
      E1  =  3.D0*GAM*E1*E1*E1
      D1  = (360.D0*A*A+921.D0*A*G+297.D0*G*G)/L1
      F1  =- ZAHL*(4608.D0*A**4+14976.D0*A**3*G+24156.D0*A**2*G**2+
     1             14472.D0*A*G**3+4860.D0*G**4)/L1**2
      G1  =  ZAHL*(9216.D0*A**4+29952.D0*A**3*G+31848.D0*A**2*G**2+
     1             28944.D0*A*G**3+9720.D0*G**4)/L1**2
      A1  =  U1/L1
      C1  =  V1/L1
C
      E2  =  DSQRT(A*A*A/((A+G)**2*(2.D0*A+3.D0*G)))
      E2  =- 3.D0*GAM*E2*E2*E2
      D2  = (102.D0*A**3+729.D0*A**2*G+792.D0*A*G**2+81.D0*G**3)/
     1      (2.D0*L2*(A+G))
      F2  =- ZAHL*(60.D0*A**4+204.D0*A**3*G+2301.D0*A**2*G**2+162.D0*
     1             A*G**3+81.D0*G**4)/L2**2
      G2  =  ZAHL*(120.D0*A**4+408.D0*A**3*G-3630.D0*A**2*G**2+324.D0*
     1             A*G**3+162.D0*G**4)/L2**2
      A2  =  U2/L2
      C2  =  V2/L2
C
      E3  =  GAM*(4.D0/3.D0)**3*(A/(A+G))**6
      D3  = (96.D0*A*A+180.D0*A*G)/(7.D0*(A+G))
      F3  =- 1376.D0*ZAHL*A*A/49.D0
      G3  =- 2736.D0*ZAHL*A*A/49.D0
      A3  =  25.D0*A/7.D0
      C3  =- 48.D0*A/7.D0
C
      DO 10 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 10 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        TERM1      =(D1+F1*(X*X+Y*Y))*S(ZAHL*A1,ZAHL*A1,-ZAHL*C1,L,X,Y)+
     1               G1*X*Y*T(ZAHL*A1,ZAHL*A1,-ZAHL*C1,L,X,Y)
        TERM2      =(D2+F2*(X*X+Y*Y))*S(ZAHL*A2,ZAHL*A2,-ZAHL*C2,L,X,Y)+
     1               G2*X*Y*T(ZAHL*A2,ZAHL*A2,-ZAHL*C2,L,X,Y)
        TERM3      =(D3+F3*(X*X+Y*Y))*S(ZAHL*A3,ZAHL*A3,-ZAHL*C3,L,X,Y)+
     1               G3*X*Y*T(ZAHL*A3,ZAHL*A3,-ZAHL*C3,L,X,Y)
C--------------------------------------------------C
C  NON-LOCAL KINETIC ENERGY KERNEL                 C
C--------------------------------------------------C
        MAT(M,N+1)  = MAT(M,N+1)-H2M*(E1*TERM1+E2*TERM2+E3*TERM3)
   10 CONTINUE
C
      RETURN
      END
C********************************************C
C  CENTRAL NN-POTENTIAL: LOCAL CONTRIBUTION  C
C********************************************C
      SUBROUTINE ZENTRA(L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,NMAX,STSANZ
      PARAMETER ( NMAX = 240 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,VT,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C-----------------------------------------------------C
C  LOCAL POTENTIAL IN THE ALPHA+HE3,TRITON CHANNEL    C
C-----------------------------------------------------C
        IN          = STSANZ+1
      DO 10 M=1,IN
        X           = DBLE(M)*HSTEP
       DO 10 K=1,3
        KONST    = DSQRT(12.D0*A*G/(12.D0*A*G+BK(K)*(8.D0*A+9.D0*G)))
        KONST    = KONST*KONST*KONST
       KONSTS    = 12.D0*A*G*BK(K)/(12.D0*A*G+BK(K)*(8.D0*A+9.D0*G))
        W15G0    = 4.D0*PW(K)-PM(K)+2.D0*PB(K)-2.D0*PH(K)
        MAT(M,M) = MAT(M,M)-3.D0*V0(K)*W15G0*KONST*DEXP(-KONSTS*X*X)
   10 CONTINUE
C
      DO 11 K=1,3
C
      ZAHL = 6.D0/7.D0
C
      U1  =  36.D0*A*A+97.D0*A*G+36.D0*G*G
      V1  =  72.D0*A*A+96.D0*A*G+72.D0*G*G
      L1  =  56.D0*A+63.D0*G
      U2  =   3.D0*A*A+31.D0*A*G+3.D0*G*G
      V2  =   6.D0*A*A-36.D0*A*G+6.D0*G*G
      L2  =  14.D0*A+21.D0*G
      GAM =  DSQRT(12.D0*G*G/(PI*A))
      GAM = (12.D0/7.D0)**3*GAM*GAM*GAM
      MU5 =  A*G*L1+14.D0*BK(K)*(A+G)*(2.D0*A+3.D0*G)
      MU6 =  G*L1+14.D0*BK(K)*(2.D0*A+9.D0*G)
      MU7 =  A*L1+14.D0*BK(K)*(8.D0*A+3.D0*G)
      MU8 =  2.D0*(A+G)*L2+7.D0*BK(K)*(13.D0*A+15.D0*G)
      MU9 =  4.D0*A*L2+42.D0*BK(K)*(A+G)
      M10 =  L2*(A+G+2.D0*BK(K))
      M11 =  4.D0*A*(A+G)*L2+14.D0*BK(K)*(10.D0*A*A+15.D0*A*G+3.D0*G*G)
      M12 =  21.D0*(A+G)+56.D0*BK(K)
      M13 =  21.D0*(A+G)+28.D0*BK(K)
C
         AA = 4.D0*PW(K)-1.D0*PM(K)+2.D0*PB(K)-2.D0*PH(K)
         BB =-1.D0*PW(K)+4.D0*PM(K)-2.D0*PB(K)+2.D0*PH(K)
C
      WAG1  =-(AA+BB)*((A/(A+2.D0*BK(K)))**(3.D0/2.D0)+(1.D0/3.D0)*
     1                 (G/(G+2.D0*BK(K)))**(3.D0/2.D0))
      WAG2  =-(AA+BB)*((2.D0/3.D0)*((A+G)/(A+G+4.D0*BK(K)))**(3.D0/2.D0)
     1                +(1.D0/3.D0)*(A/(A+2.D0*BK(K)))**(3.D0/2.D0))
      WAG3  =- 2.D0*(AA+BB)*((A+G)/(A+G+4.D0*BK(K)))**(3.D0/2.D0)
      WAG4  =  3.D0*BB*GAM
      WAG5  = (2.D0*BB-4.D0*AA)*GAM
      WAG6  =- 4.D0*(AA+BB)*GAM
      WAG7  =- 6.D0*(AA+BB)*GAM
      WAG8  =  4.D0*(AA+BB)*GAM
      WAG9  = (AA-2.D0*BB)*GAM
      WAG10 = (2.D0*AA-4.D0*BB)*GAM
      WAG11 =  8.D0*(AA+BB)*GAM
      WAG12 = (BB-2.D0*AA)*GAM
      WAG13 =- 2.D0*(AA+BB)*GAM
C
      E1 = DSQRT(A*A/(G*(8.D0*A+9.D0*G)))
      E1 = 3.D0*GAM*E1*E1*E1
      A1 = U1/L1
      C1 = V1/L1
C
      E2 = DSQRT(A*A*A/((A+G)**2*(2.D0*A+3.D0*G)))
      E2 =-3.D0*GAM*E2*E2*E2
      A2 = U2/L2
      C2 = V2/L2
C
      E3 = GAM*(4.D0/3.D0)**3*(A/(A+G))**6
      A3 = 25.D0*A/7.D0
      C3 =-48.D0*A/7.D0
C
      E4 = DSQRT(7.D0*A*A/(G*L1))
      E4 = E4*E4*E4
      A4 = U1/L1+(24.D0/7.D0)*BK(K)
      C4 = V1/L1+(48.D0/7.D0)*BK(K)
C
      E5 = DSQRT(7.D0*A**3/MU5)
      E5 = E5*E5*E5
      A5 = (A*G*U1+2.D0*BK(K)*(A+G)*(9.D0*A*A+46.D0*A*G+12.D0*G*G))/MU5
      C5 = (A*G*V1+2.D0*BK(K)*(A+G)*(18.D0*A*A-6.D0*A*G+24.D0*G*G))/MU5
C
      E6 = DSQRT(7.D0*A**2/MU6)
      E6 = E6*E6*E6
      A6 = (G*U1+BK(K)*(18.D0*A*A+278.D0*A*G+126.D0*G*G))/MU6
      B6 = (G*U1+BK(K)*(18.D0*A*A+110.D0*A*G+126.D0*G*G))/MU6
      C6 = (G*V1+BK(K)*(36.D0*A*A+192.D0*A*G+252.D0*G*G))/MU6
C
      E7 = DSQRT(7.D0*A**3/(G*MU7))
      E7 = E7*E7*E7
      A7 = (A*U1+BK(K)*(120.D0*A*A+278.D0*A*G+24.D0*G*G))/MU7
      B7 = (A*U1+BK(K)*(120.D0*A*A+110.D0*A*G+24.D0*G*G))/MU7
      C7 = (A*V1+BK(K)*(240.D0*A*A+192.D0*A*G+48.D0*G*G))/MU7
C
      E8 = DSQRT(14.D0*A**3/((A+G)*MU8))
      E8 = E8*E8*E8
      A8 = (2.D0*(A+G)*U2+BK(K)*( 57.D0*A*A+131.D0*A*G+24.D0*G*G))/MU8
      B8 = (2.D0*(A+G)*U2+BK(K)*(141.D0*A*A+215.D0*A*G+24.D0*G*G))/MU8
      C8 = (2.D0*(A+G)*V2-BK(K)*( 96.D0*A*A+144.D0*A*G-48.D0*G*G))/MU8
C
      E9 = DSQRT(28.D0*A**4/((A+G)**2*MU9))
      E9 = E9*E9*E9
      A9 = (4.D0*A*U2+BK(K)*(156.D0*A*A+162.D0*A*G+ 6.D0*G*G))/MU9
      C9 = (4.D0*A*V2-BK(K)*(276.D0*A*A+264.D0*A*G-12.D0*G*G))/MU9
C
      E10 = DSQRT(7.D0*A**3/((A+G)*M10))
      E10 = E10*E10*E10
      A10 = ((A+G)*U2+BK(K)*(18.D0*A*A+92.D0*A*G+24.D0*G*G))/M10
      C10 = ((A+G)*V2+BK(K)*(36.D0*A*A-12.D0*A*G+48.D0*G*G))/M10
C
      E11 = DSQRT(28.D0*A**4/((A+G)*M11))
      E11 = E11*E11*E11
      A11 = (4.D0*A*(A+G)*U2+BK(K)*(42.D0*A**3+424.D0*A**2*G+288.D0*A*
     1                              G**2+6.D0*G**3))/M11
      B11 = (4.D0*A*(A+G)*U2+BK(K)*(42.D0*A**3+256.D0*A**2*G+120.D0*A*
     1                              G**2+6.D0*G**3))/M11
      C11 = (4.D0*A*(A+G)*V2+BK(K)*(84.D0*A**3-300.D0*A**2*G-180.D0*A*
     1                              G**2+12.D0*G**3))/M11
C
      E12 = DSQRT(112.D0*A**4/(3.D0*(A+G)**3*M12))
      E12 = E12*E12*E12
      A12 = (75.D0*A*(A+G)+BK(K)*(208.D0*A+8.D0*G))/M12
      C12 =-(144.D0*A*(A+G)+BK(K)*(368.D0*A-16.D0*G))/M12
C
      E13 = DSQRT(112.D0*A**4/(3.D0*(A+G)**3*M13))
      E13 = E13*E13*E13
      A13 =  (75.D0*A*(A+G)+BK(K)*(228.D0*A+128.D0*G))/M13
      B13 =  (75.D0*A*(A+G)+BK(K)*(172.D0*A+ 72.D0*G))/M13
      C13 =-(144.D0*A*(A+G)+BK(K)*(384.D0*A+192.D0*G))/M13
C
       DO 11 M=1,STSANZ
         X = DBLE(M)*HSTEP
        DO 11 N=M,STSANZ
         Y = DBLE(N)*HSTEP
C
        KON1 =       E1 *WAG1 *S(ZAHL*A1 ,ZAHL*A1 ,-ZAHL*C1 ,L,X,Y)+
     1               E2 *WAG2 *S(ZAHL*A2 ,ZAHL*A2 ,-ZAHL*C2 ,L,X,Y)+
     2               E3 *WAG3 *S(ZAHL*A3 ,ZAHL*A3 ,-ZAHL*C3 ,L,X,Y)+
     3               E4 *WAG4 *S(ZAHL*A4 ,ZAHL*A4 ,-ZAHL*C4 ,L,X,Y)+
     4               E5 *WAG5 *S(ZAHL*A5 ,ZAHL*A5 ,-ZAHL*C5 ,L,X,Y)+
     5         0.5D0*E6 *WAG6 *S(ZAHL*A6 ,ZAHL*B6 ,-ZAHL*C6 ,L,X,Y)+
     6         0.5D0*E6 *WAG6 *S(ZAHL*B6 ,ZAHL*A6 ,-ZAHL*C6 ,L,X,Y)
C
        KON2 = 0.5D0*E7 *WAG7 *S(ZAHL*A7 ,ZAHL*B7 ,-ZAHL*C7 ,L,X,Y)+
     1         0.5D0*E7 *WAG7 *S(ZAHL*B7 ,ZAHL*A7 ,-ZAHL*C7 ,L,X,Y)+
     2         0.5D0*E8 *WAG8 *S(ZAHL*A8 ,ZAHL*B8 ,-ZAHL*C8 ,L,X,Y)+
     3         0.5D0*E8 *WAG8 *S(ZAHL*B8 ,ZAHL*A8 ,-ZAHL*C8 ,L,X,Y)+
     4               E9 *WAG9 *S(ZAHL*A9 ,ZAHL*A9 ,-ZAHL*C9 ,L,X,Y)+
     5               E10*WAG10*S(ZAHL*A10,ZAHL*A10,-ZAHL*C10,L,X,Y)+
     6         0.5D0*E11*WAG11*S(ZAHL*A11,ZAHL*B11,-ZAHL*C11,L,X,Y)
C
        KON3 = 0.5D0*E11*WAG11*S(ZAHL*B11,ZAHL*A11,-ZAHL*C11,L,X,Y)+
     1               E12*WAG12*S(ZAHL*A12,ZAHL*A12,-ZAHL*C12,L,X,Y)+
     2         0.5D0*E13*WAG13*S(ZAHL*A13,ZAHL*B13,-ZAHL*C13,L,X,Y)+
     3         0.5D0*E13*WAG13*S(ZAHL*B13,ZAHL*A13,-ZAHL*C13,L,X,Y)
C------------------------------------------------C
C  CENTRAL NN-POTENTIAL: NON-LOCAL CONTRIBUTION  C
C------------------------------------------------C
          MAT(M,N+1)  = MAT(M,N+1)-V0(K)*(KON1+KON2+KON3)
   11 CONTINUE
C
      RETURN
      END
C*******************C
C  COULOMB KERNEls  C
C*******************C
      SUBROUTINE COUPOT(L,FRAG)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,NMAX,FRAG,SKIP,STSANZ
      PARAMETER ( NMAX = 240 , SKIP = 1)
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C-----------------------------------------------------------C
C  direct COULOMB POTENTIAL IN THE ALPHA+HE,TRITON CHANNEL  C
C-----------------------------------------------------------C
       IN          = STSANZ+1
       DO 10 M=1,IN
        X           = DBLE(M)*HSTEP
        KONST       = DSQRT(12.D0*A*G/(8.D0*A+9.D0*G))
        MAT(M,M)    = MAT(M,M) + ((2.D0*DFLOAT(FRAG)*ELEML*
     1                                    ELEML)/X)*S15AEF(KONST*X,IFAL)
   10  CONTINUE
C---------------------------------------C
C FOR SKIP=1:  CALCULATION OF THE       C
C NON-LOCAL COULOMB POTENTIAL             C
C---------------------------------------C
      IF(SKIP.EQ.0) RETURN
C
       ZAHL = 6.D0/7.D0
C
        U1  =  36.D0*A*A+97.D0*A*G+36.D0*G*G
        V1  =  72.D0*A*A+96.D0*A*G+72.D0*G*G
        L1  =  56.D0*A+63.D0*G
        U2  =   3.D0*A*A+31.D0*A*G+3.D0*G*G
        V2  =   6.D0*A*A-36.D0*A*G+6.D0*G*G
        L2  =  14.D0*A+21.D0*G
        GAM =  DSQRT(12.D0*G*G/(PI*A))
        GAM = (12.D0/7.D0)**3*GAM*GAM*GAM
C
        E1  = DSQRT(A*A/(G*(8.D0*A+9.D0*G)))
        E1  =-3.D0*GAM*E1*E1*E1
        A1  = U1/L1
        C1  = V1/L1
        E2  = DSQRT(A*A*A/((A+G)**2*(2.D0*A+3.D0*G)))
        E2  = 3.D0*GAM*E2*E2*E2
        A2  = U2/L2
        C2  = V2/L2
        E3  =-GAM*(4.D0/3.D0)**3*(A/(A+G))**6
        A3  = 25.D0*A/7.D0
        C3  =-48.D0*A/7.D0
C
      EC121 = 2.D0*DSQRT(A*(8.D0*A+9.D0*G)/(2.D0*(8.D0*A+3.D0*G)))/3.D0
      EC261 = 2.D0*DSQRT(A*G*(8.D0*A+9.D0*G)/((A+G)*(4.D0*A+6.D0*G)))
      EC151 = DSQRT(2.D0*(A+G))/6.D0
      EC231 = DSQRT(2.D0*A/PI)/(6.D0*PI)
      EC152 = DSQRT(2.D0*(A+G))
      EC132 = (2.D0/3.D0)*DSQRT(2.D0*A*(A+G)*(2.D0*A+3.D0*G)/(10.D0*A*A+
     1         15.D0*A*G+3.D0*G*G))
      EC372 = (2.D0/3.D0)*DSQRT(2.D0*A*(2.D0*A+3.D0*G)/(3.D0*(A+G)))
      EC342 = EC231
      EC153 = DSQRT(6.D0*(A+G))
      EC143 = DSQRT(3.D0*(A+G))/2.D0
      EC161 = 2.D0*DSQRT(G*(8.D0*A+9.D0*G)/(2.D0*(2.D0*A+9.D0*G)))/3.D0
      EC671 = DSQRT(2.D0*G/PI)/(6.D0*PI)
      EC172 = (2.D0/3.D0)*DSQRT(2.D0*(A+G)*(2.D0*A+3.D0*G)/(13.D0*A+
     1         15.D0*G))
      EC122 = DSQRT((A+G)/PI)/(6.D0*PI)
      EC123 = DSQRT((A+G)/PI)/(2.D0*PI)
C
        LI1 = (6.D0/7.D0)*(4.D0*A+8.D0*G)*DSQRT(2.D0*A/((8.D0*A+9.D0*G)*
     1        (8.D0*A+3.D0*G)))
        MI1 =-LI1*(4.D0*A+G)/(4.D0*A+8.D0*G)
        LI2 = (6.D0/7.D0)*(15.D0*A+9.D0*G)*DSQRT(G/((4.D0*A+18.D0*G)*
     1        (8.D0*A+9.D0*G)))
        MI2 =-LI2*(A+9.D0*G)/(15.D0*A+9.D0*G)
        LI3 = 6.D0*DSQRT(A*G*(A+G)/((4.D0*A+6.D0*G)*(8.D0*A+9.D0*G)))
        LI4 = (3.D0/7.D0)*DSQRT(2.D0*(A+G))
        LI5 = (6.D0/7.D0)*(A+5.D0*G)*DSQRT(2.D0*A*(A+G)/((2.D0*A+3.D0*G)
     1        *(10.D0*A*A+15.D0*A*G+3.D0*G*G)))
        MI5 =-LI5*(A-2.D0*G)/(A+5.D0*G)
        LI6 = (3.D0/7.D0)*(9.D0*A+3.D0*G)*DSQRT(2.D0*(A+G)/((2.D0*A+3.D0
     1        *G)*(13.D0*A+15.D0*G)))
        MI6 = LI6*(5.D0*A-3.D0*G)/(9.D0*A+3.D0*G)
        LI7 = DSQRT(6.D0*A*(A+G)/(2.D0*A+3.D0*G))
        LI8 = DSQRT(6.D0*(A+G))/7.D0
        LI9 = (8.D0/7.D0)*DSQRT(3.D0*(A+G))
        MI9 = (6.D0/7.D0)*DSQRT(3.D0*(A+G))
C
         F1 = 1.D0/2.D0
         F2 = 2.D0/3.D0
         F3 = 1.D0/3.D0
         F4 = 1.D0/4.D0
C
      DO 11 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 11 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
        IF(FRAG.EQ.1) THEN
C
            KON1 = F1*EC121*COULI(A1,A1,LI1, MI1,ZAHL*C1,L,X,Y,1)+
     1             F1*EC121*COULI(A1,A1,MI1, LI1,ZAHL*C1,L,X,Y,1)+
     2             F2*EC261*COULI(A1,A1,LI3, LI3,ZAHL*C1,L,X,Y,1)+
     3             F1*EC151*COULI(A1,A1,LI4,-LI4,ZAHL*C1,L,X,Y,2)
C
            KON2 = 2.D0*EC231*S(A1,A1,-ZAHL*C1,L,X,Y)
C
            KON3 = F3*EC152*COULI(A2,A2,LI4,-LI4,ZAHL*C2,L,X,Y,1)+
     1                EC132*COULI(A2,A2,LI5, MI5,ZAHL*C2,L,X,Y,1)+
     2                EC132*COULI(A2,A2,MI5, LI5,ZAHL*C2,L,X,Y,1)+
     3                EC372*COULI(A2,A2,LI7, LI7,ZAHL*C2,L,X,Y,1)
C
            KON4 = EC342*S(A2,A2,-ZAHL*C2,L,X,Y)
C
            KON5 = F4*EC153*COULI(A3,A3,LI8,-LI8,ZAHL*C3,L,X,Y,1)+
     1                EC143*COULI(A3,A3,LI9, MI9,ZAHL*C3,L,X,Y,1)+
     3                EC143*COULI(A3,A3,MI9, LI9,ZAHL*C3,L,X,Y,1)
C
            TERM =-2.D0*PI*ELEML*ELEML*(E1*(KON1+KON2)+E2*(KON3+KON4)+
     1                                  E3*KON5)
                      ELSE
C
            KON1 = EC121*COULI(A1,A1,LI1, MI1,ZAHL*C1,L,X,Y,1)+
     1             EC121*COULI(A1,A1,MI1, LI1,ZAHL*C1,L,X,Y,1)+
     2             EC161*COULI(A1,A1,LI2, MI2,ZAHL*C1,L,X,Y,1)+
     3             EC161*COULI(A1,A1,MI2, LI2,ZAHL*C1,L,X,Y,1)+
     4             EC261*COULI(A1,A1,LI3, LI3,ZAHL*C1,L,X,Y,1)+
     5             EC151*COULI(A1,A1,LI4,-LI4,ZAHL*C1,L,X,Y,2)
C
            KON2 = (EC231+EC671)*S(A1,A1,-ZAHL*C1,L,X,Y)
C
            KON3 = EC152*COULI(A2,A2,LI4,-LI4,ZAHL*C2,L,X,Y,1)+
     1             EC132*COULI(A2,A2,LI5, MI5,ZAHL*C2,L,X,Y,1)+
     2             EC132*COULI(A2,A2,MI5, LI5,ZAHL*C2,L,X,Y,1)+
     3             EC172*COULI(A2,A2,LI6, MI6,ZAHL*C2,L,X,Y,1)+
     4             EC172*COULI(A2,A2,MI6, LI6,ZAHL*C2,L,X,Y,1)+
     5             EC372*COULI(A2,A2,LI7, LI7,ZAHL*C2,L,X,Y,1)
C
            KON4 = 2.D0*EC122*S(A2,A2,-ZAHL*C2,L,X,Y)
C
            KON5 = EC153*COULI(A3,A3,LI8,-LI8,ZAHL*C3,L,X,Y,1)
C
            KON6 = 2.D0*EC123*S(A3,A3,-ZAHL*C3,L,X,Y)
C
            TERM =-2.D0*PI*ELEML*ELEML*(E1*(KON1+KON2)+E2*(KON3+KON4)+
     1                                  E3*(KON5+KON6))

        ENDIF
C
        MAT(M,N+1)  = MAT(M,N+1)+TERM
   11 CONTINUE
C
      RETURN
      END
C****************C
C  LS - KERNELS  C
C****************C
      SUBROUTINE SPIORB(JL,L)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,JL,NMAX,STSANZ
      PARAMETER ( NMAX = 240 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),MAT(NMAX+1,NMAX+1)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,VT,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      IF((VO+VT.EQ.0.D0).OR.(L.EQ.0)) RETURN
C----------------------------------------------------------C
C  LOCAL SPIN-ORBIT POTENTIAL IN ALPHA+HE, TRITON CHANNEL  C
C----------------------------------------------------------C
      IN          = STSANZ+1
      DO 10 M=1,IN
         X        = DBLE(M)*HSTEP
         KONST1   = 12.D0*A*G/(12.D0*A*G+LAM*(8.D0*A+9.D0*G))
         KONST    = DSQRT(KONST1)
         KONST    = KONST*KONST*KONST*KONST*KONST
         MAT(M,M) = MAT(M,M)-(7.D0/6.D0)*VO*JLS(JL,L)*KONST*
     1                                           DEXP(-(KONST1*LAM)*X*X)
   10 CONTINUE
C
      ZAHL  = 6.D0/7.D0
C
       VOR  = 12.D0/7.D0
       VOR  = VOR*VOR*VOR
C
       T151 = 8.D0*A+9.D0*G
       E151 = 48.D0*LAM/7.D0+(72.D0*A*A+96.D0*A*G+72.D0*G*G)/(7.D0*T151)
       S151 =-48.D0*A*G/(T151*E151)
       G151 = DSQRT(12.D0*A*G/(PI*T151))
       G151 = G151*G151*G151
       D151 = 24.D0*LAM/7.D0+(36.D0*A*A+97.D0*A*G+36.D0*G*G)/(7.D0*T151)
C
       T261 = A*G*(8.D0*A+9.D0*G)+2.D0*LAM*(A+G)*(2.D0*A+3.D0*G)
       E261 = (A*G*(72.D0*A*A+96.D0*A*G+72.D0*G*G)+2.D0*LAM*(A+G)*(18.D0
     1        *A*A-6.D0*A*G+24.D0*G*G))/(7.D0*T261)
       S261 = 12.D0*A*G*(A+G)**2/(T261*E261)
       G261 = DSQRT(12.D0*A*A*G*G/(PI*T261))
       G261 = G261*G261*G261
       D261 = (A*G*(36.D0*A*A+97.D0*A*G+36.D0*G*G)+2.D0*LAM*(A+G)*(9.D0
     1        *A*A+46.D0*A*G+12.D0*G*G))/(7.D0*T261)
C
       T152 = (2.D0*A+3.D0*G)*(A+G+2.D0*LAM)
       E152 = ((A+G)*(6.D0*A*A-36.D0*A*G+6.D0*G*G)+LAM*(36.D0*A*A-12.D0
     1        *A*G+48.D0*G*G))/(7.D0*T152)
       S152 =-12.D0*A*G*(A+G)/(T152*E152)
       G152 = DSQRT(12.D0*A*A*G*G/(PI*(A+G)*T152))
       G152 = G152*G152*G152
       D152 = ((A+G)*(3.D0*A*A+31.D0*A*G+3.D0*G*G)+LAM*(18.D0*A*A+92.D0
     1        *A*G+24.D0*G*G))/(7.D0*T152)
C
       T372 = 4.D0*A*(2.D0*A+3.D0*G)+6.D0*LAM*(A+G)
       E372 = (4.D0*A*(6.D0*A*A-36.D0*A*G+6.D0*G*G)-LAM*(276.D0*A*A+
     1        264.D0*A*G-12.D0*G*G))/(7.D0*T372)
       S372 = 12.D0*A*(A+G)**2/(T372*E372)
       G372 = DSQRT(48.D0*A*A*A*G*G/(PI*(A+G)**2*T372))
       G372 = G372*G372*G372
       D372 = (4.D0*A*(3.D0*A*A+31.D0*A*G+3.D0*G*G)+LAM*(156.D0*A*A+
     1        162.D0*A*G+6.D0*G*G))/(7.D0*T372)
C
       T153 = 3.D0*(A+G)+8.D0*LAM
       E153 =-(144.D0*A*(A+G)+LAM*(368.D0*A-16.D0*G))/(7.D0*T153)
       S153 =-16.D0*A*(A+G)/(T153*E153)
       G153 = DSQRT(64.D0*A*A*A*G*G/(PI*(A+G)**3*T153))
       G153 = G153*G153*G153
       D153 = (75.D0*A*(A+G)+LAM*(208.D0*A+8.D0*G))/(7.D0*T153)
C
      DO 11 M=1,STSANZ
        X = DBLE(M)*HSTEP
       DO 11 N=M,STSANZ
        Y = DBLE(N)*HSTEP
C
         KON1 = S151*G151*S(ZAHL*D151,ZAHL*D151,-ZAHL*E151,L,X,Y)
         KON2 = S261*G261*S(ZAHL*D261,ZAHL*D261,-ZAHL*E261,L,X,Y)
         KON3 = S152*G152*S(ZAHL*D152,ZAHL*D152,-ZAHL*E152,L,X,Y)
         KON4 = S372*G372*S(ZAHL*D372,ZAHL*D372,-ZAHL*E372,L,X,Y)
         KON5 = S153*G153*S(ZAHL*D153,ZAHL*D153,-ZAHL*E153,L,X,Y)
C-------------------------------------------------C
C  SPIN-ORBIT POTENTIAL: NON-LOCAL PART           C
C-------------------------------------------------C
          KON = (VO+3.D0*VT)*KON1+2.D0*(VO+3.D0*VT)*KON2-4.D0*VO*KON3-
     1           6.D0*VT*KON4+3.D0*(VO-VT)*KON5
C
          MAT(M,N+1)  = MAT(M,N+1)+JLS(JL,L)*VOR*KON
   11 CONTINUE
C
      RETURN
      END
C**************************************C
C  WRITE THE MATRIX MAT ON A FILE      C
C**************************************C
      SUBROUTINE WRITIN(EALPHA,ETHREE)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,IN,NMAX,STSANZ
      PARAMETER ( NMAX = 240 )
      REAL*8    MAT(NMAX+1,NMAX+1)
C
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      OPEN(20,FILE='ALPHATRI.OUT',STATUS='NEW',FORM='UNFORMATTED')
C
      WRITE(20) EALPHA+ETHREE
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
C-----------------------------------------------------------------C
C     J=1(J VON A-HE3,T IST L+1/2) J=-1(J VON A-HE3,T IST L-1/2)  C
C-----------------------------------------------------------------C
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
