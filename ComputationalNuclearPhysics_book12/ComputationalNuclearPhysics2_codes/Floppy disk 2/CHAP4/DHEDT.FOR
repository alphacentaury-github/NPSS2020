C=============================================================C
C THE PROGRAM CALCULATES RGM-KERNELS FOR THE                  C
C DEUTERON + HELIUM3, DEUTERON + TRITON SYSTEMS               C
C=============================================================C
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,PO,JL,DIS,FRAG,SBRA,SKET,NMAX,KMAX,
     1          STSANZ,sbra1,sket1,ll1
C-------------------------------------------C
C  NMAX: MAXIMUM NUMBER OF MESH POINTS      C
C  KMAX: MAXIMUM NUMBER OF CHANNELS         C
C-------------------------------------------C
      PARAMETER ( KMAX =  3 , NMAX = 70 )
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),B(3),AVOR(3,3),
     1          EDEUTE(3),MAT(NMAX+1,NMAX+1,KMAX,KMAX)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G,B,AVOR
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
      OPEN(5,FILE='DHEDT.INP')
C
      READ(5,*)     FRAG
      READ(5,100)
      READ(5,*)     PO
      READ(5,100)
      READ(5,*)     DIS
      READ(5,100)
      READ(5,*)     J,L,SbRA1,SKET1

          sbra=2
         if (sbra1.eq.1) sbra=1
          sket=2
         if (sket1.eq.1) sket=1

          ll1=2*l
         if (j.eq.(ll1+3)) jl=2
         if (j.eq.(ll1+1)) jl=1
         if (j.eq.(ll1-1)) jl=-1
         if (j.eq.(ll1-3)) jl=-2




                    IF(((SBRA*SKET.EQ.1).AND.(IABS(JL).NE.1)).OR.
     1                      ((SBRA*SKET.EQ.2).AND.(IABS(JL).NE.1)).OR.
     2                           ((SBRA*SKET.EQ.4).AND.(IABS(JL).GT.2)))
     3                                             THEN
                            WRITE(*,*) 'VALUE FOR J IS NOT CORRECT'
                                                   STOP
                    ENDIF
      READ(5,200)
      READ(5,*)     MM
      READ(5,*)     VO,LAM
      READ(5,*)     HSTEP,STSANZ
                    IF(STSANZ.GT.NMAX) THEN
                     WRITE(*,*) 'NUMBER OF MESHPOINTS IS TOO LARGE'
                     STOP
                    ENDIF
      CLOSE(5)
C
  100 FORMAT(1X,//)
  200 FORMAT(1X,///)
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
                                 G  = 0.367D0
                   IF(FRAG.NE.2) G  = 0.378D0
                 B(1) = 0.07284D0
                 B(2) = 0.36570D0
                 B(3) = 1.46960D0
C
              AVOR(1,1) = 0.0438511974D0
              AVOR(2,1) = 0.1592095494D0
              AVOR(3,1) = 0.2519553542D0
C
              AVOR(1,2) =-0.1214143990D0
              AVOR(2,2) = 0.3193592544D0
              AVOR(3,2) = 0.1198678585D0
C
              AVOR(1,3) = 0.0504353095D0
              AVOR(2,3) =-0.5290197661D0
              AVOR(3,3) = 1.4174135447D0
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
                   IF(FRAG.NE.2) G  = 0.690D0
                 B(1) = 0.320D0
                 B(2) = 0.440D0
                 B(3) = 0.600D0
C
              AVOR(1,1) = 1.1130824469D0
              AVOR(2,1) =-2.2102402241D0
              AVOR(3,1) = 1.5368531890D0
C
              AVOR(1,2) =-2.4016415961D0
              AVOR(2,2) = 4.2766276372D0
              AVOR(3,2) =-1.5562673117D0
C
              AVOR(1,3) = 3.8234300316D0
              AVOR(2,3)=-10.5292546988D0
              AVOR(3,3) = 7.3266860773D0
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
                                 A  = 0.580D0
                                 G  = 0.460D0
                   IF(FRAG.NE.2) G  = 0.460D0
                 B(1) = 0.120D0
                 B(2) = 1.200D0
                 B(3) = 1.400D0
C
              AVOR(1,1) = 0.0924494391D0
              AVOR(2,1) = 1.7517299019D0
              AVOR(3,1) =-1.5374092840D0
C
              AVOR(1,2) = 0.1468048247D0
              AVOR(2,2) =-3.9229456793D0
              AVOR(3,2) = 3.4827875051D0
C
              AVOR(1,3) = 0.0502491686D0
              AVOR(2,3)= -8.9267477273D0
              AVOR(3,3) =10.2133481498D0
        ENDIF
      ENDIF
C-------------------------------------C
C  CALCULATION OF THE MATRIXELEMENTS  C
C-------------------------------------C
      CALL INITIA
      CALL BIND(FRAG,NF,ETHREE,EDEUTE)
C
      IF(SBRA.EQ.SKET) THEN
      WRITE(*,*) 'CALCULATION OF THE KINETIC ENERGY KERNEL'
              CALL   TKIN(DIS,L,SBRA,NF)

      WRITE(*,*) 'CALCULATION OF THE KERNEL FOR CENTRAL FORCE'
                CALL   ZENTRA(DIS,L,SBRA,NF)

      WRITE(*,*) 'CALCULATION OF THE COULOMB KERNEL'
                  CALL   COUPOT(DIS,FRAG,L,SBRA,NF)

      WRITE(*,*) 'CALCULATION OF THE LS-KERNEL'
                    CALL   SPIORB(DIS,JL,L,SBRA,SKET,NF)

      WRITE(*,*) 'CALCULATION OF THE NORM KERNEL'
                      CALL   NORM(DIS,L,SBRA,NF)

                        ELSE
      WRITE(*,*) 'CALCULATION OF THE LS-KERNEL'
                          CALL  SPIORB(DIS,JL,L,SBRA,SKET,NF)
      ENDIF
C
      CALL WRITIN(DIS,SBRA,SKET,ETHREE,EDEUTE)
C
      STOP
      END
C*************************************C
C  INITIALIZE THE MATRIX MAT          C
C*************************************C
      SUBROUTINE INITIA
      INTEGER*4 I,J,M,N,KMAX,NMAX
      PARAMETER ( KMAX =  3 , NMAX = 70 )
      REAL*8    MAT(NMAX+1,NMAX+1,KMAX,KMAX)
C
      COMMON /MATRIX/ MAT
C
               DO 10 J=1,KMAX
                 DO 10 I=1,KMAX
                   DO 10 N=1,NMAX+1
                     DO 10 M=1,NMAX+1
                         MAT(M,N,I,J) = 0.D0
   10          CONTINUE
      RETURN
      END
C******************************************************************C
C  THIS SUBROUTINE CALCULATES THE BINDING ENERGIES AND THE         C
C  NORMALIZATION CONSTANTS FOR THE VARIOUS CLUSTERS                C
C******************************************************************C
      SUBROUTINE BIND(FRAG,NF,ETHREE,EDEUTE)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,FRAG
      REAL*8    V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),B(3),AVOR(3,3),
     1          EDEUTE(3)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
C---------------------------C
C  NORMALIZATION CONSTANTS  C
C---------------------------C
          NF  = (6.D0*PI)/(DSQRT(3.D0)*G)
          NF  =     NF*NF*NF
C---------------------------------C
C INTERNAL ENERGY OF THE DEUTERON C
C---------------------------------C
      DO 10 K1=1,3
       EDEUTE(K1) = 0.D0
        DO 10 I=1,3
         DO 10 J=1,3
            KONST1 = DSQRT(PI/(B(I)+B(J)))
            KONST1 = AVOR(I,K1)*AVOR(J,K1)*KONST1*KONST1*KONST1
            KONST2 = H2M*((3.D0*B(I)*B(J))/(B(I)+B(J)))
            KONST3 = 0.D0
C
           DO 11 K=1,3
            KONST4 = DSQRT((B(I)+B(J))/(B(I)+B(J)+4.D0*BK(K)))
            KONST4 = KONST4*KONST4*KONST4
            KONST3 = KONST3+V0(K)*(PW(K)+PM(K)+PB(K)+PH(K))*KONST4
   11      CONTINUE
C
            EDEUTE(K1) = EDEUTE(K1)+KONST1*(KONST2-KONST3)
   10  CONTINUE
C-------------------------------------------------C
C INTERNAL ENERGY OF THE HELIUM3, TRITON CLUSTERS C
C-------------------------------------------------C
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
      SUBROUTINE NORM(DIS,L,SBRA,NF)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,K1,K2,DIS,SBRA,STSANZ,KMAX,NMAX
      PARAMETER ( KMAX =  3 , NMAX = 70 )
      REAL*8    AVOR(3,3),B(3),T1(2),T2(2),MAT(NMAX+1,NMAX+1,KMAX,KMAX)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
            ZAHL      = 3.D0/5.D0
C----------------------------C
C  NON-LOCAL NORM KERNEL     C
C----------------------------C
            T1(1)     = 0.5D0
            T2(1)     = 0.5D0
            T1(2)     = 2.D0
            T2(2)     =-1.D0
C
      DO 10 I=1,3
       DO 10 J=1,3
C
      E1    = (36.D0/5.D0)*DSQRT(PI/(2.D0*G))*DSQRT((2.D0*PI)/
     1        (3.D0*G+2.D0*B(J)+2.D0*B(I)))
      E1    =  E1*E1*E1
      E2    = (36.D0/5.D0)*DSQRT(PI/(G+B(J)))*DSQRT(PI/(G+B(I)))
      E2    = E2*E2*E2
      X1    = 6.D0*G*G+20.D0*G*B(J)+5.D0*G*B(I)+6.D0*B(J)*B(I)
      X2    = 6.D0*G*G+5.D0*G*B(J)+20.D0*G*B(I)+6.D0*B(J)*B(I)
      X3    = 12.D0*G*G+12.D0*B(J)*B(I)
      YY    = 5.D0*(3.D0*G+2.D0*B(J)+2.D0*B(I))
      A1    = X1/YY
      A2    = (13.D0/5.D0)*G
      B1    = X2/YY
      B2    = (13.D0/5.D0)*G
      C1    = X3/YY
      C2    =-(24.D0/5.D0)*G
C
        DO 11 M=1,STSANZ
           X = DBLE(M)*HSTEP
         DO 11 N=(1/DIS)*(M-1)+1,STSANZ
           Y = DBLE(N)*HSTEP
C
          WERT =-(T1(SBRA)*E1*S(ZAHL*A1,ZAHL*B1,-ZAHL*C1,L,X,Y)+
     1            T2(SBRA)*E2*S(ZAHL*A2,ZAHL*B2,-ZAHL*C2,L,X,Y))/NF
C
          DO 12 K1=1,DIS
            IF(M.LE.N) MAT(N+1,M,K1,K1) = MAT(N+1,M,K1,K1)+
     1                                        AVOR(I,K1)*AVOR(J,K1)*WERT
           DO 12 K2=K1+1,DIS
            IF(M.LE.N) THEN
               MAT(N+1,M,K2,K1) = MAT(N+1,M,K2,K1)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
                       ELSE
               MAT(N,M+1,K2,K1) = MAT(N,M+1,K2,K1)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
            ENDIF
   12     CONTINUE
C
   11   CONTINUE
C
   10 CONTINUE
C
      RETURN
      END
C*********************************C
C  KINETIC ENERGY KERNELS         C
C*********************************C
      SUBROUTINE TKIN(DIS,L,SBRA,NF)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,K1,K2,DIS,SBRA,NMAX,KMAX,STSANZ
      PARAMETER ( KMAX =  3 , NMAX = 70 )
      REAL*8    AVOR(3,3),B(3),T1(2),T2(2),MAT(NMAX+1,NMAX+1,KMAX,KMAX)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
            ZAHL      = 3.D0/5.D0
            T1(1)     = 0.5D0
            T2(1)     = 0.5D0
            T1(2)     = 2.D0
            T2(2)     =-1.D0
C
      DO 10 I=1,3
       DO 10 J=1,3
C
      E1    = (36.D0/5.D0)*DSQRT(PI/(2.D0*G))*DSQRT((2.D0*PI)/
     1        (3.D0*G+2.D0*B(J)+2.D0*B(I)))
      E1    =  E1*E1*E1
      E2    = (36.D0/5.D0)*DSQRT(PI/(G+B(J)))*DSQRT(PI/(G+B(I)))
      E2    = E2*E2*E2
      X1    = 6.D0*G*G+20.D0*G*B(J)+5.D0*G*B(I)+6.D0*B(J)*B(I)
      X2    = 6.D0*G*G+5.D0*G*B(J)+20.D0*G*B(I)+6.D0*B(J)*B(I)
      X3    = 12.D0*G*G+12.D0*B(J)*B(I)
      YY    = 5.D0*(3.D0*G+2.D0*B(J)+2.D0*B(I))
      A1    = X1/YY
      A2    = (13.D0/5.D0)*G
      B1    = X2/YY
      B2    = (13.D0/5.D0)*G
      C1    = X3/YY
      C2    =-(24.D0/5.D0)*G
C
      TA    = (G*G*(12.D0*G+4.D0*B(J)+B(I))+B(I)*B(I)*(10.D0*G+
     1        12.D0*B(J)))/YY-(5.D0*X1*(3.D0*G*G+4.D0*B(I)*
     2        B(I)))/(YY*YY)
      TB    = (G*G*(12.D0*G+9.D0*B(J)+36.D0*B(I))+B(I)*B(I)*(40.D0*G+
     1        12.D0*B(J)))/YY-(5.D0*X2*(3.D0*G*G+4.D0*B(I)*B(I)))/
     2        (YY*YY)
      TC    = (G*G*(24.D0*G-12.D0*B(J)+12.D0*B(I))+B(I)*B(I)*
     1        (24.D0*B(J)))/YY-(5.D0*X3*(3.D0*G*G+4.D0*B(I)*
     2        B(I)))/(YY*YY)
      D1    = 4.5D0*G+3.D0*B(I)-(9.D0*G*G+12.D0*B(I)*B(I))/
     1        (6.D0*G+4.D0*B(J)+4.D0*B(I))+3.D0*A1
C
      F1    =-(3.D0/5.D0)*TA-(6.D0/5.D0)*A1*A1
      F1S   =-(3.D0/5.D0)*TB-0.3D0*C1*C1
      G1    = (3.D0/5.D0)*TC+(6.D0/5.D0)*A1*C1
      D2    = 6.D0*G+3.D0*B(I)-3.D0*G*G/(G+B(J))-3.D0*B(I)*B(I)
     1        /(G+B(I))+3.D0*A2
      F2    =-(24.D0/25.D0)*G*G-(6.D0/5.D0)*A2*A2
      F2S   =-(54.D0/25.D0)*G*G-0.3D0*C2*C2
      G2    =-(72.D0/25.D0)*G*G+(6.D0/5.D0)*A2*C2
C---------------------------------------------C
C   NON-LOCAL KINTEIC ENERGY KERNELS          C
C---------------------------------------------C
        DO 11 M=1,STSANZ
            X = DBLE(M)*HSTEP
         DO 11 N=(1/DIS)*(M-1)+1,STSANZ
            Y = DBLE(N)*HSTEP
C
       TERM1 = (T1(SBRA)*E1)*((D1+F1*X*X+F1S*Y*Y)*S(ZAHL*A1,ZAHL*B1,
     1         -ZAHL*C1,L,X,Y)+G1*X*Y*T(ZAHL*A1,ZAHL*B1,-ZAHL*C1,L,X,Y))
       TERM2 = (T2(SBRA)*E2)*((D2+F2*X*X+F2S*Y*Y)*S(ZAHL*A2,ZAHL*B2,
     1         -ZAHL*C2,L,X,Y)+G2*X*Y*T(ZAHL*A2,ZAHL*B2,-ZAHL*C2,L,X,Y))
C
          WERT  = -H2M*(TERM1+TERM2)/NF
C
          DO 12 K1=1,DIS
            IF(M.LE.N) MAT(M,N+1,K1,K1) = MAT(M,N+1,K1,K1)+
     1                                       AVOR(I,K1)*AVOR(J,K1)*WERT
           DO 12 K2=K1+1,DIS
            IF(M.LE.N) THEN
             MAT(M,N+1,K1,K2) = MAT(M,N+1,K1,K2)+
     1                                       AVOR(I,K1)*AVOR(J,K2)*WERT
                       ELSE
             MAT(M+1,N,K1,K2) = MAT(M+1,N,K1,K2)+
     1                                       AVOR(I,K1)*AVOR(J,K2)*WERT
            ENDIF
   12     CONTINUE
C
   11   CONTINUE
C
   10 CONTINUE
C
      RETURN
      END
C**********************************C
C  KERNELS FOR CENTRAL POTENTIAL   C
C**********************************C
      SUBROUTINE ZENTRA(DIS,L,SBRA,NF)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,K1,K2,IN,DIS,SBRA,STSANZ,KMAX,NMAX
      PARAMETER ( KMAX =  3 , NMAX = 70 )
      REAL*8    W121(2),W141(2),W151(2),W231(2),W241(2),W251(2),
     1          W122(2),W132(2),W142(2),W342(2),W452(2),W140(2),
     2          W451(2),V0(3),PM(3),PB(3),PW(3),PH(3),AVOR(3,3),
     3          B(3),BK(3),MAT(NMAX+1,NMAX+1,KMAX,KMAX)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C-------------------------------------------C
C CENTRAL NN-POTENTIAL: LOCAL CONTRIBUTION  C
C S IS THE SPIN IN THE D+HE3, D+T CHANNEL   C
C-------------------------------------------C
      IN = STSANZ+1
      DO 10 K1=1,DIS
        DO 10 K2=K1,DIS
          DO 11 M=1,IN
                X = DBLE(M)*HSTEP
C
           DO 14 K=1,3
            KONST4 = 0.D0
             W140(1)  = 6.D0*PW(K)-0.5D0*PM(K)+1.D0*PB(K)-3.D0*PH(K)
             W140(2)  = 6.D0*PW(K)-2.0D0*PM(K)+4.D0*PB(K)-3.D0*PH(K)
C
            DO 15 I=1,3
             DO 15 J=1,3
              KONST1= DSQRT(PI/(B(I)+B(J)))
              KONST1= AVOR(I,K1)*AVOR(J,K2)*KONST1*KONST1*KONST1
              KONST2= DSQRT((3.D0*G*(B(I)+B(J)))/(3.D0*G*(B(I)+B(J))+
     1                BK(K)*(3.D0*G+2.D0*B(I)+2.D0*B(J))))
              KONST2= KONST2*KONST2*KONST2
              KONST3= DEXP(-((3.D0*G*BK(K)*(B(I)+B(J)))/(3.D0*G*(B(I)+
     1                B(J))+BK(K)*(3.D0*G+2.D0*B(I)+2.D0*B(J))))*X*X)
              KONST4= KONST4+KONST1*KONST2*KONST3
   15       CONTINUE
                 MAT(M,M,K1,K2) = MAT(M,M,K1,K2)-V0(K)*W140(SBRA)*KONST4
   14      CONTINUE
   11     CONTINUE
   10 CONTINUE
C
              ZAHL = 3.D0/5.D0
C-------------------------------------------------C
C   CENTRAL NN-POTENTIAL: NON-LOCAL CONTRIBUTIONS C
C-------------------------------------------------C
      DO 16 I=1,3
       DO 16 J=1,3

       E1    = (36.D0/5.D0)*DSQRT(PI/(2.D0*G))*DSQRT((2.D0*PI)/
     1         (3.D0*G+2.D0*B(J)+2.D0*B(I)))
       E1    =  E1*E1*E1
       E2    = (36.D0/5.D0)*DSQRT(PI/(G+B(J)))*DSQRT(PI/(G+B(I)))
       E2    = E2*E2*E2
       X1    = 6.D0*G*G+20.D0*G*B(J)+5.D0*G*B(I)+6.D0*B(J)*B(I)
       X2    = 6.D0*G*G+5.D0*G*B(J)+20.D0*G*B(I)+6.D0*B(J)*B(I)
       X3    = 12.D0*G*G+12.D0*B(J)*B(I)
       YY    = 5.D0*(3.D0*G+2.D0*B(J)+2.D0*B(I))
       A1    = X1/YY
       A2    = (13.D0/5.D0)*G
       B1    = X2/YY
       B2    = (13.D0/5.D0)*G
       C1    = X3/YY
       C2    =-(24.D0/5.D0)*G
C
        DO 16 K=1,3
C
        W121(1)=- 1.D0*PW(K)- 1.D0*PM(K)+ 2.D0*PB(K)+ 2.D0*PH(K)
        W141(1)=-0.5D0*PW(K)-0.5D0*PM(K)-0.5D0*PB(K)-0.5D0*PH(K)
        W151(1)=-0.5D0*PW(K)+ 6.D0*PM(K)- 3.D0*PB(K)+ 1.D0*PH(K)
        W231(1)=-0.5D0*PW(K)-0.5D0*PM(K)- 2.D0*PB(K)- 2.D0*PH(K)
        W241(1)=- 1.D0*PW(K)- 2.D0*PM(K)+ 1.D0*PB(K)- 1.D0*PH(K)
        W251(1)=- 1.D0*PW(K)- 1.D0*PM(K)+ 2.D0*PB(K)+ 2.D0*PH(K)
        W451(1)=-0.5D0*PW(K)-0.5D0*PM(K)-0.5D0*PB(K)-0.5D0*PH(K)
        W122(1)=-0.5D0*PW(K)-0.5D0*PM(K)-0.5D0*PB(K)-0.5D0*PH(K)
        W132(1)=- 1.D0*PW(K)- 1.D0*PM(K)+0.5D0*PB(K)+0.5D0*PH(K)
        W142(1)=- 2.D0*PW(K)- 1.D0*PM(K)- 1.D0*PB(K)+ 1.D0*PH(K)
        W342(1)=- 1.D0*PW(K)- 1.D0*PM(K)+0.5D0*PB(K)+0.5D0*PH(K)
        W452(1)=-0.5D0*PW(K)-0.5D0*PM(K)-0.5D0*PB(K)-0.5D0*PH(K)
C
        W121(2)=-4.D0*PW(K)-4.D0*PM(K)-1.D0*PB(K)-1.D0*PH(K)
        W141(2)=-2.D0*PW(K)-2.D0*PM(K)-2.D0*PB(K)-2.D0*PH(K)
        W151(2)=-2.D0*PW(K)+6.D0*PM(K)-3.D0*PB(K)+4.D0*PH(K)
        W231(2)=-2.D0*PW(K)-2.D0*PM(K)+1.D0*PB(K)+1.D0*PH(K)
        W241(2)=-4.D0*PW(K)+4.D0*PM(K)-2.D0*PB(K)+5.D0*PH(K)
        W251(2)=-4.D0*PW(K)-4.D0*PM(K)-1.D0*PB(K)-1.D0*PH(K)
        W451(2)=-2.D0*PW(K)-2.D0*PM(K)-2.D0*PB(K)-2.D0*PH(K)
        W122(2)= 1.D0*PW(K)+1.D0*PM(K)+1.D0*PB(K)+1.D0*PH(K)
        W132(2)= 2.D0*PW(K)+2.D0*PM(K)-1.D0*PB(K)-1.D0*PH(K)
        W142(2)= 4.D0*PW(K)-4.D0*PM(K)+5.D0*PB(K)-2.D0*PH(K)
        W342(2)= 2.D0*PW(K)+2.D0*PM(K)-1.D0*PB(K)-1.D0*PH(K)
        W452(2)= 1.D0*PW(K)+1.D0*PM(K)+1.D0*PB(K)+1.D0*PH(K)
C
      E121  = DSQRT(G*YY/(G*YY+BK(K)*(30.D0*G+5.D0*B(J)+5.D0*B(I))))
      E121  = E121*E121*E121
      A121  = (G*X1+BK(K)*(21.D0*G*G+16.D0*G*B(J)+4.D0*G*B(I)+3.D0*
     1        B(J)*B(I)))/(G*YY+BK(K)*(30.D0*G+5.D0*B(J)+5.D0*B(I)))
      B121  = (G*X2+BK(K)*(21.D0*G*G+16.D0*G*B(J)+64.D0*G*B(I)+3.D0*
     1        B(J)*B(I)))/(G*YY+BK(K)*(30.D0*G+5.D0*B(J)+5.D0*B(I)))
      C121  = (G*X3+BK(K)*(42.D0*G*G-18.D0*G*B(J)+18.D0*G*B(I)+6.D0*
     1        B(J)*B(I)))/(G*YY+BK(K)*(30.D0*G+5.D0*B(J)+5.D0*B(I)))
C
      E251  = E121
      B251  = (G*X2+BK(K)*(21.D0*G*G+16.D0*G*B(I)+4.D0*G*B(J)+
     1        3.D0*B(I)*B(J)))/(G*YY+BK(K)*(30.D0*G+5.D0*B(I)+5.D0*
     2        B(J)))
      A251  = (G*X1+BK(K)*(21.D0*G*G+16.D0*G*B(I)+64.D0*G*B(J)+
     1        B(I)*3.D0*B(J)))/(G*YY+BK(K)*(30.D0*G+5.D0*B(I)+5.D0*
     2        B(J)))
      C251  = (G*X3+BK(K)*(42.D0*G*G-18.D0*G*B(I)+18.D0*G*B(J)+
     1        B(I)*6.D0*B(J)))/(G*YY+BK(K)*(30.D0*G+5.D0*B(I)+5.D0*
     2        B(J)))
C
      E151  = 1.D0
      A151  = (X1+BK(K)*(36.D0*G+24.D0*B(J)+24.D0*B(I)))/YY
      B151  = (X2+BK(K)*(36.D0*G+24.D0*B(J)+24.D0*B(I)))/YY
      C151  = (X3+BK(K)*(72.D0*G+48.D0*B(J)+48.D0*B(I)))/YY
C
      E141  = DSQRT(YY/(YY+40.D0*BK(K)))
      E141  = E141*E141*E141
      A141  = (X1+BK(K)*(80.D0*G+24.D0*B(I)))/(YY+40.D0*BK(K))
      B141  = (X2+BK(K)*(20.D0*G+24.D0*B(I)))/(YY+40.D0*BK(K))
      C141  = (X3+BK(K)*48.D0*B(I))/(YY+40.D0*BK(K))
C
      E451  = E141
      A451  = (X1+BK(K)*(20.D0*G+24.D0*B(J)))/(YY+40.D0*BK(K))
      B451  = (X2+BK(K)*(80.D0*G+24.D0*B(J)))/(YY+40.D0*BK(K))
      C451  = (X3+BK(K)*48.D0*B(J))/(YY+40.D0*BK(K))
C
      E231  = DSQRT(G/(G+2.D0*BK(K)))
      E231  = E231*E231*E231
      A231  = A1
      B231  = B1
      C231  = C1
C
      E241  = DSQRT(G*YY/(G*YY+BK(K)*(10.D0*G+5.D0*B(J)+5.D0*B(I))))
      E241  = E241*E241*E241
      A241  = (G*X1+BK(K)*(29.D0*G*G+16.D0*G*B(J)+16.D0*G*B(I)+3.D0*
     1        B(J)*B(I)))/(G*YY+BK(K)*(10.D0*G+5.D0*B(J)+5.D0*B(I)))
      B241  = (G*X2+BK(K)*(29.D0*G*G+16.D0*G*B(J)+16.D0*G*B(I)+3.D0*
     1        B(J)*B(I)))/(G*YY+BK(K)*(10.D0*G+5.D0*B(J)+5.D0*B(I)))
      C241  = (G*X3+BK(K)*(-42.D0*G*G-18.D0*G*B(J)-18.D0*G*B(I)+6.D0*
     1        B(J)*B(I)))/(G*YY+BK(K)*(10.D0*G+5.D0*B(J)+5.D0*B(I)))
C
      E122  = DSQRT((G+B(J))/(G+B(J)+4.D0*BK(K)))
      E122  = E122*E122*E122
      A122  = (13.D0/5.D0)*G
      B122  = A122
      C122  =-(24.D0/5.D0)*G
C
      E452  = DSQRT((G+B(I))/(G+B(I)+4.D0*BK(K)))
      E452  = E452*E452*E452
      A452  = (13.D0/5.D0)*G
      B452  = A452
      C452  =-(24.D0/5.D0)*G
C
      E132  = DSQRT((G+B(J))/(G+B(J)+BK(K)))
      E132  = E132*E132*E132
      A132  = 0.2D0*(13.D0*G*(G+B(J))+BK(K)*(25.D0*G+12.D0*B(J)))/
     1        (G+B(J)+BK(K))
C
      B132  = 0.2D0*(13.D0*G*(G+B(J))+BK(K)*(40.D0*G+27.D0*B(J)))/
     1        (G+B(J)+BK(K))
      C132  =-0.2D0*(24.D0*G*(G+B(J))+BK(K)*(60.D0*G+36.D0*B(J)))/
     1        (G+B(J)+BK(K))
C
      E342  = DSQRT((G+B(I))/(G+B(I)+BK(K)))
      E342  = E342*E342*E342
      A342  = 0.2D0*(13.D0*G*(G+B(I))+BK(K)*(40.D0*G+27.D0*B(I)))/
     1        (G+B(I)+BK(K))
      B342  = 0.2D0*(13.D0*G*(G+B(I))+BK(K)*(25.D0*G+12.D0*B(I)))/
     1        (G+B(I)+BK(K))
      C342  =-0.2D0*(24.D0*G*(G+B(I))+BK(K)*(60.D0*G+36.D0*B(I)))/
     1        (G+B(I)+BK(K))
C
      E142  = DSQRT(((G+B(J))*(G+B(I)))/((G+B(J))*(G+B(I))+BK(K)*
     1        (2.D0*G+B(J)+B(I))))
      E142  = E142*E142*E142
      A142  = 0.2D0*(13.D0*G*(G+B(J))*(G+B(I))+BK(K)*(29.D0*G*G+
     1        16.D0*G*B(J)+16.D0*G*B(I)+3.D0*B(J)*B(I)))/((G+B(J))*
     2        (G+B(I))+BK(K)*(2.D0*G+B(J)+B(I)))
      B142  = A142
      C142  =-0.2D0*(24.D0*G*(G+B(J))*(G+B(I))+BK(K)*(42.D0*G*G+
     1        18.D0*G*B(J)+18.D0*G*B(I)-6.D0*B(J)*B(I)))/((G+B(J))*
     2        (G+B(I))+BK(K)*(2.D0*G+B(J)+B(I)))
C
         DO 17 M=1,STSANZ
              X = DBLE(M)*HSTEP
          DO 17 N=(1/DIS)*(M-1)+1,STSANZ
              Y = DBLE(N)*HSTEP
C
      KONST1= W121(SBRA)*E121*S(ZAHL*A121,ZAHL*B121,-ZAHL*C121,L,X,Y)+
     1        W251(SBRA)*E251*S(ZAHL*A251,ZAHL*B251,-ZAHL*C251,L,X,Y)+
     2        W151(SBRA)*E151*S(ZAHL*A151,ZAHL*B151,-ZAHL*C151,L,X,Y)
C
      KONST2= W141(SBRA)*E141*S(ZAHL*A141,ZAHL*B141,-ZAHL*C141,L,X,Y)+
     1        W451(SBRA)*E451*S(ZAHL*A451,ZAHL*B451,-ZAHL*C451,L,X,Y)+
     2        W231(SBRA)*E231*S(ZAHL*A231,ZAHL*B231,-ZAHL*C231,L,X,Y)+
     3        W241(SBRA)*E241*S(ZAHL*A241,ZAHL*B241,-ZAHL*C241,L,X,Y)
C
      KONST3= W122(SBRA)*E122*S(ZAHL*A122,ZAHL*B122,-ZAHL*C122,L,X,Y)+
     1        W452(SBRA)*E452*S(ZAHL*A452,ZAHL*B452,-ZAHL*C452,L,X,Y)+
     2        W132(SBRA)*E132*S(ZAHL*A132,ZAHL*B132,-ZAHL*C132,L,X,Y)+
     3        W342(SBRA)*E342*S(ZAHL*A342,ZAHL*B342,-ZAHL*C342,L,X,Y)+
     4        W142(SBRA)*E142*S(ZAHL*A142,ZAHL*B142,-ZAHL*C142,L,X,Y)
C
           WERT =-V0(K)*(E1*(KONST1+KONST2)+E2*KONST3)/NF
C
           DO 18 K1=1,DIS
             IF(M.LE.N) MAT(M,N+1,K1,K1)= MAT(M,N+1,K1,K1)+
     1                                        AVOR(I,K1)*AVOR(J,K1)*WERT
            DO 18 K2=K1+1,DIS
             IF(M.LE.N) THEN
              MAT(M,N+1,K1,K2)= MAT(M,N+1,K1,K2)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
                        ELSE
              MAT(M+1,N,K1,K2)= MAT(M+1,N,K1,K2)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
             ENDIF
   18      CONTINUE
C
   17    CONTINUE
C
   16 CONTINUE
C
      RETURN
      END
C**************************************C
C  KERNELS FOR THE COULOMB POTENTIAL   C
C**************************************C
      SUBROUTINE COUPOT(DIS,FRAG,L,SBRA,NF)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,K1,K2,IN,DIS,SKIP,FRAG,SBRA,STSANZ,KMAX,
     1          NMAX,IFAIL
      PARAMETER ( KMAX =  3 , NMAX = 70 ,  SKIP = 1 )
      REAL*8    AVOR(3,3),B(3),WC(2),MAT(NMAX+1,NMAX+1,KMAX,KMAX)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C--------------------------------------C
C     COULOMB POTENTIAL: LOCAL PART    C
C--------------------------------------C
      IN = STSANZ+1
      DO 10 K1=1,DIS
       DO 10 K2= K1,DIS
        DO 11 M =1,IN
              X = DBLE(M)*HSTEP
         KONST4 = 0.D0
C
         DO 12 I=1,3
          DO 12 J=1,3
           KONST1 = DSQRT(PI/(B(I)+B(J)))
           KONST1 = AVOR(I,K1)*AVOR(J,K2)*KONST1*KONST1*KONST1
           KONST2 = DSQRT((3.D0*G*(B(I)+B(J)))/(3.D0*G+2.D0*B(I)+2.D0*
     1              B(J)))
           KONST3 = S15AEF(KONST2*X,IFAIL)
           KONST4 = KONST4+KONST1*KONST3
   12    CONTINUE
                MAT(M,M,K1,K2) = MAT(M,M,K1,K2)
     1                                  +DBLE(FRAG)*ELEML*ELEML/X*KONST4
   11   CONTINUE
   10 CONTINUE
C---------------------------------------C
C FOR SKIP=1: CALCULATION OF THE        C
C NON-LOCAL COULOMB POTENTIAL           C
C---------------------------------------C
      IF(SKIP.EQ.0) RETURN
C
            ZAHL = 3.D0/5.D0
C-------------------------------------------------C
C   COULOMB POTENTIAL: NON-LOCAL PART             C
C-------------------------------------------------C
      DO 13 I=1,3
       DO 13 J=1,3
C
      E1    = (36.D0/5.D0)*DSQRT(PI/(2.D0*G))*DSQRT((2.D0*PI)/
     1        (3.D0*G+2.D0*B(J)+2.D0*B(I)))
      E1    =  E1*E1*E1
      E2    = (36.D0/5.D0)*DSQRT(PI/(G+B(J)))*DSQRT(PI/(G+B(I)))
      E2    = E2*E2*E2
      X1    = 6.D0*G*G+20.D0*G*B(J)+5.D0*G*B(I)+6.D0*B(J)*B(I)
      X2    = 6.D0*G*G+5.D0*G*B(J)+20.D0*G*B(I)+6.D0*B(J)*B(I)
      X3    = 12.D0*G*G+12.D0*B(J)*B(I)
      YY    = 5.D0*(3.D0*G+2.D0*B(J)+2.D0*B(I))
      A1    = X1/YY
      A2    = (13.D0/5.D0)*G
      B1    = X2/YY
      B2    = (13.D0/5.D0)*G
      C1    = X3/YY
      C2    =-(24.D0/5.D0)*G
C
       IF(FRAG.EQ.2) THEN
           EC    = PI*(36.D0/5.D0)
           EC    =-EC*EC*EC*ELEML*ELEML
           WC(1) =-0.5D0
           WC(2) = 1.D0
C
           KONSTS= (1.D0/(3.D0*G*G+2.D0*G*B(J)+2.D0*G*B(I)))
           P121  = KONSTS*DSQRT(1.D0/(6.D0*G+B(J)+B(I)))
C
           KONST = DSQRT(1.D0/(3.D0*G+2.D0*B(J)+2.D0*B(I)))
           P151  = (5.D0/(6.D0*G))*KONST*KONST*KONST
C
           P231  = (1.D0/(G*DSQRT(2.D0*PI*PI*PI)))*KONST*KONST*KONST
           P241  = 2.D0*KONSTS*DSQRT(1.D0/(2.D0*G+B(J)+B(I)))
C
           KONST = DSQRT(1.D0/(G+B(I)))
           P132  =-(1.D0/(G+B(J)))*KONST*KONST*KONST
C
           KONST = DSQRT(1.D0/(G+B(J)))
           KONSTS= DSQRT(1.D0/(2.D0*G+B(J)+B(I)))
           P152  =-(1.D0/((G+B(J))*(G+B(I))))*KONSTS
           P352  =-(1.D0/(G+B(I)))*KONST*KONST*KONST
C
           KONST = DSQRT(G/((6.D0*G+B(J)+B(I))*(3.D0*G+2.D0*B(J)+
     1             2.D0*B(I))))
           L1    = (3.D0/5.D0)*(-3.D0*G+2.D0*B(J)-B(I))*KONST
           M1    = (3.D0/5.D0)*(3.D0*G+3.D0*B(J)+6.D0*B(I))*KONST
C
           L2    = DSQRT(G)
           M2    =-L2
C
           KONSTS= DSQRT(G/((2.D0*G+B(J)+B(I))*(3.D0*G+2.D0*B(J)+
     1             2.D0*B(I))))
           L3    = 0.3D0*(10.D0*G+4.D0*B(J)+6.D0*B(I))*KONSTS
           M3    = 0.3D0*(10.D0*G+6.D0*B(J)+4.D0*B(I))*KONSTS
C
           L4    = (9.D0/5.D0)*DSQRT(G+B(I))
           M4    = (6.D0/5.D0)*DSQRT(G+B(I))
C
           L5    = 0.6D0*DSQRT(((G+B(J))*(G+B(I)))/(2.D0*G+B(J)+B(I)))
           M5    =-L5
C
           L6    = (3.D0/5.D0)*(3.D0*G+3.D0*B(I)+6.D0*B(J))*KONST
           M6    = (3.D0/5.D0)*(-3.D0*G+2.D0*B(I)-B(J))*KONST
C
           L7    = (6.D0/5.D0)*DSQRT(G+B(J))
           M7    = (9.D0/5.D0)*DSQRT(G+B(J))
                    ELSE
           EC    = PI*(36.D0/5.D0)
           EC    =-EC*EC*EC*ELEML*ELEML
           WC(1) =-0.5D0
           WC(2) = 1.D0
C
           P151  = DSQRT(1.D0/(3.D0*G+2.D0*B(J)+2.D0*B(I)))
           P151  = (5.D0/(6.D0*G))*P151*P151*P151
           L2    = DSQRT(G)
           M2    =-L2
C
           P152  =-(1.D0/((G+B(J))*(G+B(I))))*
     1             DSQRT(1.D0/(2.D0*G+B(J)+B(I)))
           KONST = ((G+B(J))*(G+B(I)))/(2.D0*G+B(J)+B(I))
           L5    = (3.D0/5.D0)*DSQRT(KONST)
           M5    =-L5
C
           P241  = (1.D0/G)*(1.D0/(3.D0*G+2.D0*B(J)+2.D0*B(I)))*
     1             DSQRT(1.D0/(2.D0*G+B(J)+B(I)))
           KONSTS= G*(10.D0*G+4.D0*B(J)+6.D0*B(I))*(10.D0*G+4.D0*B(J)+
     1             B(I)*6.D0)
           KONST = (2.D0*G+B(J)+B(I))*(3.D0*G+2.D0*B(J)+2.D0*B(I))
           L3    = 0.3D0*DSQRT(KONSTS/KONST)
C
           KONSTS= G*(10.D0*G+6.D0*B(J)+4.D0*B(I))*(10.D0*G+6.D0*B(J)+
     1             B(I)*4.D0)
           M3    = 0.3D0*DSQRT(KONSTS/KONST)
       ENDIF
C
        DO 14 M=1,STSANZ
            X = DBLE(M)*HSTEP
         DO 14 N=(1/DIS)*(M-1)+1,STSANZ
            Y = DBLE(N)*HSTEP
C
         IF(FRAG.EQ.2) THEN
          TERM1    = P121*COULI(ZAHL*A1,ZAHL*B1,L1,M1,-ZAHL*C1,L,X,Y,1)+
     1               P121*COULI(ZAHL*A1,ZAHL*B1,L6,M6,-ZAHL*C1,L,X,Y,1)+
     2               P151*COULI(ZAHL*A1,ZAHL*B1,L2,M2,-ZAHL*C1,L,X,Y,2)
C
          TERM2    = P241*COULI(ZAHL*A1,ZAHL*B1,L3,M3,-ZAHL*C1,L,X,Y,1)+
     1               P352*COULI(ZAHL*A2,ZAHL*B2,L4,M4,-ZAHL*C2,L,X,Y,1)+
     2               P132*COULI(ZAHL*A2,ZAHL*B2,L7,M7,-ZAHL*C2,L,X,Y,1)+
     3               P152*COULI(ZAHL*A2,ZAHL*B2,L5,M5,-ZAHL*C2,L,X,Y,1)
C
          TERM3    = P231*S(ZAHL*A1,ZAHL*B1,-ZAHL*C1,L,X,Y)
C
          TERM     = EC*(WC(SBRA)*TERM3+2.D0*PI*(TERM1+WC(SBRA)*TERM2))
                    ELSE
          TERM     = 2.D0*PI*EC*(P241*COULI(ZAHL*A1,ZAHL*B1,L3,M3,-ZAHL*
     1               C1,L,X,Y,1)+WC(SBRA)*(P151*COULI(ZAHL*A1,ZAHL*B1,L2
     2              ,M2,-ZAHL*C1,L,X,Y,2)+P152*COULI(ZAHL*A2,ZAHL*B2,L5,
     3               M5,-ZAHL*C2,L,X,Y,1)))
         ENDIF
C
          WERT = TERM/NF
C
          DO 15 K1=1,DIS
             IF(M.LE.N) MAT(M,N+1,K1,K1) = MAT(M,N+1,K1,K1)+
     1                                        AVOR(I,K1)*AVOR(J,K1)*WERT
           DO 15 K2=K1+1,DIS
            IF(M.LE.N) THEN
             MAT(M,N+1,K1,K2) = MAT(M,N+1,K1,K2)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
                                       ELSE
             MAT(M+1,N,K1,K2) = MAT(M+1,N,K1,K2)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
            ENDIF
   15     CONTINUE
C
   14   CONTINUE
C
   13 CONTINUE
C
      RETURN
      END
C***************C
C  LS-KERNELS   C
C***************C
      SUBROUTINE SPIORB(DIS,JL,L,SBRA,SKET,NF)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,K1,K2,JL,IN,DIS,SBRA,SKET,STSANZ,KMAX,NMAX
      PARAMETER ( KMAX =  3 , NMAX = 70 )
      REAL*8    AVOR(3,3),B(3),V0(3),PW(3),PM(3),PH(3),PB(3),BK(3),
     1          MAT(NMAX+1,NMAX+1,KMAX,KMAX)
C
      COMMON /KONSTA/ PI,ELEML,H2M
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      IF((VO.EQ.0.D0).OR.(L.EQ.0)) RETURN
C
C-------------------------------------C
C SPIN-ORBIT POTENTIAL: LOCAL PART    C
C-------------------------------------C
      IN = STSANZ+1
      DO 10 K1= 1,DIS
       DO 10 K2= K1,DIS
        DO 10 M= 1,IN
              X = DBLE(M)*HSTEP
         KONST4 = 0.D0
C
         DO 11 I=1,3
          DO 11 J=1,3
           KONST1= DSQRT(PI/(B(I)+B(J)))
           KONST1= AVOR(I,K1)*AVOR(J,K2)*KONST1*KONST1*KONST1
           KONST2= DSQRT((3.D0*G*(B(I)+B(J)))/(3.D0*G*(B(I)+B(J))+(3.D0
     1             *G+2.D0*B(I)+2.D0*B(J))*LAM))
           KONST2= KONST2*KONST2*KONST2*KONST2*KONST2
           KONST3= DEXP(-((3.D0*G*LAM*(B(I)+B(J)))/(3.D0*G*(B(I)+B(J))+
     1             LAM*(3.D0*G+2.D0*B(I)+2.D0*B(J))))*X*X)
           KONST4= KONST4+KONST1*KONST2*KONST3
   11    CONTINUE
           MAT(M,M,K1,K2) = MAT(M,M,K1,K2)-VO*JLS(JL,L,SBRA,SKET)*KONST4
   10 CONTINUE
C
      DO 12 I=1,3
       DO 12 J=1,3
C
        X1  = 6.D0*G*G+20.D0*G*B(J)+5.D0*G*B(I)+6.D0*B(I)*B(J)
        X2  = 6.D0*G*G+20.D0*G*B(I)+5.D0*G*B(J)+6.D0*B(I)*B(J)
        YY  = 3.D0*G+2.D0*B(I)+2.D0*B(J)
C
        DLS1  =-6.D0**3*96.D0*(3.D0/5.D0)**5*(PI/(2.D0*G))**(3.D0/2.D0)
     1          *(PI/(B(I)+B(J)))**(3.D0/2.D0)*G*(2.D0*(B(I)+B(J))/YY)**
     2          (5.D0/2.D0)
        ALS1  = (3.D0/25.D0)*(X1+LAM*(36.D0*G+24.D0*B(I)+24.D0*B(J)))/YY
        BLS1  = (3.D0/25.D0)*(X2+LAM*(36.D0*G+24.D0*B(I)+24.D0*B(J)))/YY
        CLS1  =-(36.D0/25.D0)*(G*G+B(I)*B(J)+LAM*(6.D0*G+4.D0*B(I)+4.D0*
     1           B(J)))/YY
C
        DLS2  = 6.D0**3*96.D0*(3.D0/5.D0)**5*(PI/(G+B(I)))**(3.D0/2.D0)
     1          *(PI/(G+B(J)))**(3.D0/2.D0)*G*((G+B(I))*(G+B(J))/(G*YY+
     2          LAM*(2.D0*G+B(I)+B(J))))**(5.D0/2.D0)
        ALS2  = (3.D0/25.D0)*(G*X1+LAM*(29.D0*G*G+16.D0*G*B(I)+16.D0*G*
     1          B(J)+3.D0*B(I)*B(J)))/(G*YY+LAM*(2.D0*G+B(I)+B(J)))
        BLS2  = (3.D0/25.D0)*(G*X2+LAM*(29.D0*G*G+16.D0*G*B(I)+16.D0*G*
     1          B(J)+3.D0*B(I)*B(J)))/(G*YY+LAM*(2.D0*G+B(I)+B(J)))
        CLS2  =-(18.D0/25.D0)*(G*(2.D0*G*G+2.D0*B(I)*B(J))-LAM*(7.D0*G*
     1          G+3.D0*G*B(I)+3.D0*G*B(J)-B(I)*B(J)))/(G*YY+LAM*(2.D0*G+
     2          B(I)+B(J)))
C
        DLS3  =-6.D0**3*96.D0*(3.D0/5.D0)**5*(PI/(G+B(I)))**(3.D0/2.D0)
     1          *(PI/(G+B(J)))**(3.D0/2.D0)*G*((G+B(I))*(G+B(J))/((G+
     2          B(I))*(G+B(J))+LAM*(2.D0*G+B(I)+B(J))))**(5.D0/2.D0)
        ALS3  = (3.D0/25.D0)*(13.D0*G*(G+B(I))*(G+B(J))+LAM*(29.D0*G*G+
     1          16.D0*G*B(I)+16.D0*G*B(J)+3.D0*B(I)*B(J)))/((G+B(I))*(G
     2         +B(J))+LAM*(2.D0*G+B(I)+B(J)))
        BLS3  = ALS3
        CLS3  = (18.D0/25.D0)*(4.D0*G*(G+B(I))*(G+B(J))+LAM*(7.D0*G*G+
     1          3.D0*G*B(I)+3.D0*G*B(J)-B(I)*B(J)))/((G+B(I))*(G+B(J))+
     2          LAM*(2.D0*G+B(I)+B(J)))
C-----------------------------------------C
C SPIN-ORBIT POTENTIAL: NON-LOCAL PART    C
C-----------------------------------------C
        DO 13 M=1,STSANZ
           X = DBLE(M)*HSTEP
         DO 13 N=1,STSANZ
           Y = DBLE(N)*HSTEP
C
             TERM1    = DLS1/CLS1*S(ALS1,BLS1,CLS1,L,X,Y)
             TERM2    = DLS2/CLS2*S(ALS2,BLS2,CLS2,L,X,Y)
             TERM3    = DLS3/CLS3*S(ALS3,BLS3,CLS3,L,X,Y)
C
                             TERM = 0.5D0*TERM1-0.4D0*TERM2+0.8D0*TERM3
          IF(SBRA*SKET.NE.1) TERM = 0.5D0*TERM1+0.5D0*TERM2-TERM3
C
          WERT =-VO*JLS(JL,L,SBRA,SKET)*TERM/NF
C
          IF(SBRA.EQ.SKET) THEN
C
           DO 14 K1=1,DIS
            IF(M.LE.N) MAT(M,N+1,K1,K1) = MAT(M,N+1,K1,K1)+
     1                                        AVOR(I,K1)*AVOR(J,K1)*WERT
            DO 14 K2=K1+1,DIS
             IF(M.LE.N) THEN
              MAT(M,N+1,K1,K2) = MAT(M,N+1,K1,K2)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
                        ELSE
              MAT(M+1,N,K1,K2) = MAT(M+1,N,K1,K2)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
             ENDIF
   14      CONTINUE
C
                          ELSE
C
           DO 15 K1=1,DIS
            DO 15 K2=K1,DIS
             IF(M.LE.N) THEN
              MAT(M,N+1,K1,K2) = MAT(M,N+1,K1,K2)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
                        ELSE
              MAT(M+1,N,K1,K2) = MAT(M+1,N,K1,K2)+
     1                                        AVOR(I,K1)*AVOR(J,K2)*WERT
             ENDIF
   15      CONTINUE
C
          ENDIF
C
   13   CONTINUE
C
   12 CONTINUE
C
      RETURN
      END
C**************************************C
C  WRITE MATRIX MAT ON A FILE          C
C**************************************C
      SUBROUTINE WRITIN(DIS,SBRA,SKET,ETHREE,EDEUTE)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4    I,J,K,L,M,N,K1,K2,KK,IN,NIN,DIS,SBRA,SKET,KMAX,NMAX,
     1             STSANZ
      PARAMETER    ( KMAX =  3 , NMAX = 70 )
      REAL*8       EDEUTE(3),MAT(NMAX+1,NMAX+1,KMAX,KMAX)
      CHARACTER*13 OUTNAME(6)
C
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      OUTNAME(1) = 'DHEDT.OU1'
      OUTNAME(2) = 'DHEDT.OU2'
      OUTNAME(3) = 'DHEDT.OU3'
      OUTNAME(4) = 'DHEDT.OU4'
      OUTNAME(5) = 'DHEDT.OU5'
      OUTNAME(6) = 'DHEDT.OU6'
C
      KK = 1
      IN = STSANZ+1
C
      IF(SBRA.EQ.SKET) THEN
c
       DO 10 K1=1,DIS
        DO 10 K2=K1,DIS
C
          NIN = 20+KK
          OPEN(NIN,FILE=OUTNAME(KK),STATUS='NEW',FORM='UNFORMATTED')
C
          IF(K1.EQ.K2) WRITE(NIN) EDEUTE(K1)+ETHREE
C
          DO 11 M=1,IN
             WRITE(NIN) MAT(M,M,K1,K2)
   11     CONTINUE
C
          IF(K1.EQ.K2) THEN
                   DO 12 M=1,STSANZ
                    DO 12 N=M,STSANZ
                            WRITE(NIN) MAT(M,N+1,K1,K2),MAT(N+1,M,K1,K2)
   12              CONTINUE
                       ELSE
                   DO 13 M=1,STSANZ
                    DO 13 N=1,STSANZ
                     IF(M.LE.N) THEN
                            WRITE(NIN) MAT(M,N+1,K1,K2),MAT(N+1,M,K2,K1)
                                ELSE
                            WRITE(NIN) MAT(M+1,N,K1,K2),MAT(N,M+1,K2,K1)
                     ENDIF
   13              CONTINUE
          ENDIF
          CLOSE(NIN)
C
       KK = KK+1
   10  CONTINUE
                       ELSE
       DO 14 K1=1,DIS
        DO 14 K2=K1,DIS
C
          NIN = 20+KK
          OPEN(NIN,FILE=OUTNAME(KK),STATUS='NEW',FORM='UNFORMATTED')
C
          DO 15 M=1,IN
             WRITE(NIN) MAT(M,M,K1,K2)
   15     CONTINUE
C
                   DO 16 M=1,STSANZ
                    DO 16 N=1,STSANZ
                     IF(M.LE.N) THEN
                            WRITE(NIN) MAT(M,N+1,K1,K2),0.D0
                                ELSE
                            WRITE(NIN) MAT(M+1,N,K1,K2),0.D0
                     ENDIF
   16              CONTINUE
C
          CLOSE(NIN)
C
       KK = KK+1
   14  CONTINUE
c
      ENDIF
C
      RETURN
      END
C---------------------------------------------------------------C
C  J=2 (J=L+3/2)  J=1 (J=L+1/2)  J=-1 (J=L-1/2)  J=-2 (J=L-3/2) C
C---------------------------------------------------------------C
      DOUBLE PRECISION FUNCTION JLS(JL,L,SBRA,SKET)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,JL,SBRA,SKET
C
      JLS  = 0.D0
      DREH = DBLE(L)
C
       IF((SBRA.EQ.2).AND.(SKET.EQ.2)) THEN
                        IF(JL.EQ. 2) JLS= 10.D0/3.D0*DREH
                        IF(JL.EQ. 1) JLS= 10.D0/9.D0*(DREH-3.D0)
                        IF(JL.EQ.-1) JLS=-10.D0/9.D0*(DREH+4.D0)
                        IF(JL.EQ.-2) JLS=-10.D0/3.D0*(DREH+1.D0)
                                       ELSE
           IF((SBRA.EQ.1).AND.(SKET.EQ.1)) THEN
                        IF(JL.EQ. 1) JLS= 25.D0/18.D0*DREH
                        IF(JL.EQ.-1) JLS=-25.D0/18.D0*(DREH+1.D0)
                                           ELSE
           IF(JL.EQ. 1) JLS=-5.D0/18.D0*DSQRT(DREH*(2.D0*DREH+3.D0))
           IF(JL.EQ.-1) JLS=-5.D0/18.D0*
     1                              DSQRT((2.D0*DREH-1.D0)*(DREH+1.D0))
           ENDIF
       ENDIF
C
      RETURN
      END
