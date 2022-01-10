$LARGE
C===================================================C
C  THE PROGRAM CALCULATES THE RGM-KERNELS FOR THE   C
C  COUPLING BETWEEN THE (D+HE3,T)-(P,N+A)-CHANNELS  C
C===================================================C
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4  I,J,K,L,M,N,PO,JL,DIS,SBRA,FRAG,NMAX,KMAX,STSANZ
C-------------------------------------------C
C  NMAX: MAXIMUM NUMBER OF MESH POINTS      C
C  KMAX: MAXIMUM NUMBER OF CHANNELS         C
C-------------------------------------------C
      PARAMETER  ( KMAX = 3 , NMAX = 70 )
      REAL*8     V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),B(3),AVOR(3,3),
     1           MAT(NMAX,NMAX,KMAX,2)
C
      COMMON /KONSTA/ PI,H2M,ELEML
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
      OPEN(5,FILE='COUPLE.INP')
C
      READ(5,*)     FRAG
      READ(5,100)
      READ(5,*)     PO
      READ(5,100)
      READ(5,*)     DIS
      READ(5,100)
      READ(5,*)     J,L,SBRA
                    JL=J-2*L
                    IF((JL.NE.1).AND.(JL.NE.-1)) THEN
                            WRITE(*,*) 'VALUE FOR J IS NOT CORRECT'
                                                  STOP
                    ENDIF
      READ(5,100)
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
C
      CLOSE(5)
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
      CALL NFAK(NO)
C
      IF(SBRA.EQ.1) THEN
      WRITE(*,*) 'CALCULATION OF THE KINETIC ENERGY KERNEL'
              CALL   TKIN(DIS,L,NO)

      WRITE(*,*) 'CALCULATION OF THE KERNEL FOR CENTRAL FORCE'
                CALL   ZENTRA(DIS,L,NO)

      WRITE(*,*) 'CALCULATION OF THE COULOMB KERNEL'
                  CALL   COUPOT(DIS,FRAG,L,NO)

      WRITE(*,*) 'CALCULATION OF THE LS-KERNEL'
                    CALL   SPIORB(DIS,JL,L,SBRA,NO)

      WRITE(*,*) 'CALCULATION OF THE NORM KERNEL'
                      CALL   NORM(DIS,L,NO)

                        ELSE
      WRITE(*,*) 'CALCULATION OF THE LS-KERNEL'
                          CALL  SPIORB(DIS,JL,L,SBRA,NO)
      ENDIF
C
      CALL WRITIN(DIS)
C
      STOP
      END
C*************************************C
C  INITIALIZE THE MATRIX MAT          C
C*************************************C
      SUBROUTINE INITIA
      INTEGER*4 I,J,M,N,KMAX,NMAX
      PARAMETER ( KMAX =  3 , NMAX = 70 )
      REAL*8    MAT(NMAX,NMAX,KMAX,2)
C
      COMMON /MATRIX/ MAT
C
               DO 10 J=1,2
                 DO 10 I=1,KMAX
                   DO 10 N=1,NMAX
                     DO 10 M=1,NMAX
                         MAT(M,N,I,J) = 0.D0
   10          CONTINUE
      RETURN
      END
C********************************C
C  NORMALIZATION CONSTANTS       C
C********************************C
      SUBROUTINE NFAK(NO)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N
      REAL*8    AVOR(3,3),B(3)
C
      COMMON /KONSTA/ PI,H2M,ELEML
      COMMON /BREITE/ A,G,B,AVOR
C
      KONST1 = PI/(4.D0*A)
      KONST1 = 16.D0*DSQRT(KONST1*KONST1*KONST1)
      NG     = KONST1*KONST1*KONST1
C
      KONST2 = (6.D0*PI)/(DSQRT(3.D0)*G)
      NF     = KONST2*KONST2*KONST2
C
      NO     = DSQRT(NF*NG)
C
      RETURN
      END
C*****************C
C  NORM KERNEL    C
C*****************C
      SUBROUTINE NORM(DIS,L,NO)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4  I,J,K,L,M,N,K1,DIS,NMAX,KMAX,STSANZ
      PARAMETER  ( KMAX = 3 , NMAX = 70 )
      REAL*8     AVOR(3,3),B(3),MAT(NMAX,NMAX,KMAX,2)
C
      COMMON /KONSTA/ PI,H2M,ELEML
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      DO 10 I=1,3
C
       KONST = (24.D0/5.D0)*DSQRT(4.D0/3.D0)*(PI/(A+G))
       EC0   = DSQRT(3.D0/2.D0)*KONST*KONST*KONST
       AC0   = (24.D0/25.D0)*A+(9.D0/25.D0)*B(I)
       BC0   = (6.D0/25.D0)*A+(16.D0/25.D0)*B(I)
       CC0   = (24.D0/25.D0)*(A+B(I))
C
       KONST = (24.D0/5.D0)*DSQRT(PI/(A+G))*DSQRT(PI/(A+B(I)))
       EC1   = DSQRT(3.D0/2.D0)*KONST*KONST*KONST
       AC1   = (18.D0/25.D0)*A+(3.D0/25.D0)*G
       BC1   = (2.D0/25.D0)*A+(12.D0/25.D0)*G
       CC1   =-(12.D0/25.D0)*(A+G)
C
        DO 11 M=1,STSANZ
           X = DBLE(M)*HSTEP
         DO 11 N=1,STSANZ
           Y = DBLE(N)*HSTEP
C
         WERT = (EC0*S(AC0,BC0,-CC0,L,X,Y)-EC1*S(AC1,BC1,-CC1,L,X,Y))/NO
C
          DO 12 K1=1,DIS
             MAT(M,N,K1,2) = MAT(M,N,K1,2)+AVOR(I,K1)*WERT
   12     CONTINUE
C
   11   CONTINUE
   10 CONTINUE
C
      RETURN
      END
C*********************************C
C  KINTETIC ENERGY KERNELS        C
C*********************************C
      SUBROUTINE TKIN(DIS,L,NO)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4  I,J,K,L,M,N,K1,DIS,NMAX,KMAX,STSANZ
      PARAMETER  ( KMAX = 3 , NMAX = 70 )
      REAL*8     AVOR(3,3),B(3),MAT(NMAX,NMAX,KMAX,2)
C
      COMMON /KONSTA/ PI,H2M,ELEML
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      DO 10 I=1,3
C
       KONST = (24.D0/5.D0)*DSQRT(4.D0/3.D0)*(PI/(A+G))
       EC0   = DSQRT(3.D0/2.D0)*KONST*KONST*KONST
       AC0   = (24.D0/25.D0)*A+(9.D0/25.D0)*B(I)
       BC0   = (6.D0/25.D0)*A+(16.D0/25.D0)*B(I)
       CC0   = (24.D0/25.D0)*(A+B(I))
       DC0   =-6.D0*G-3.D0*B(I)+(6.D0*G*G)/(A+G)-5.D0*AC0
       FC0   = (18.D0/25.D0)*B(I)*B(I)+(10.D0/3.D0)*AC0*AC0
       FC0S  = (32.D0/25.D0)*B(I)*B(I)+(5.D0/6.D0)*CC0*CC0
       GC0   =-(48.D0/25.D0)*B(I)*B(I)-(10.D0/3.D0)*AC0*CC0
C
       KONST = (24.D0/5.D0)*DSQRT(PI/(A+G))*DSQRT(PI/(A+B(I)))
       EC1   = DSQRT(3.D0/2.D0)*KONST*KONST*KONST
       AC1   = (18.D0/25.D0)*A+(3.D0/25.D0)*G
       BC1   = (2.D0/25.D0)*A+(12.D0/25.D0)*G
       CC1   =-(12.D0/25.D0)*(A+G)
       DC1   =-6.D0*G-3.D0*B(I)+(3.D0*G*G)/(A+G)+(3.D0*B(I)*B(I))
     1         /(A+B(I))-5.D0*AC1
       FC1   = (6.D0/25.D0)*G*G+(10.D0/3.D0)*AC1*AC1
       FC1S  = (24.D0/25.D0)*G*G+(5.D0/6.D0)*CC1*CC1
       GC1   = (24.D0/25.D0)*G*G-(10.D0/3.D0)*AC1*CC1
C----------------------------------------C
C NON-LOCAL KINETIC ENERGY KERNELS       C
C----------------------------------------C
       DO 11 M=1,STSANZ
           X = DBLE(M)*HSTEP
        DO 11 N=1,STSANZ
           Y = DBLE(N)*HSTEP
C
        KONST1    = (DC0+FC0*X*X+FC0S*Y*Y)*S(AC0,BC0,-CC0,L,X,Y)+
     1               GC0*X*Y*T(AC0,BC0,-CC0,L,X,Y)
        KONST2    = (DC1+FC1*X*X+FC1S*Y*Y)*S(AC1,BC1,-CC1,L,X,Y)+
     1               GC1*X*Y*T(AC1,BC1,-CC1,L,X,Y)
C
          WERT    = -H2M*(EC0*KONST1-EC1*KONST2)/NO
C
         DO 12 K1=1,DIS
           MAT(M,N,K1,1) = MAT(M,N,K1,1)+AVOR(I,K1)*WERT
   12    CONTINUE
C
   11  CONTINUE
   10 CONTINUE
C
      RETURN
      END
C************************C
C  CENTRAL NN-POTENTIAL  C
C************************C
      SUBROUTINE ZENTRA(DIS,L,NO)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4  I,J,K,L,M,N,K1,DIS,NMAX,KMAX,STSANZ
      PARAMETER  ( KMAX = 3 , NMAX = 70 )
      REAL*8     V0(3),PW(3),PM(3),PB(3),PH(3),BK(3),B(3),AVOR(3,3),
     1           MAT(NMAX,NMAX,KMAX,2)
C
      COMMON /KONSTA/ PI,H2M,ELEML
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      DO 10 I=1,3
C
      KONST = (24.D0/5.D0)*DSQRT(4.D0/3.D0)*(PI/(A+G))
      EC0   = DSQRT(3.D0/2.D0)*KONST*KONST*KONST
      AC0   = (24.D0/25.D0)*A+(9.D0/25.D0)*B(I)
      BC0   = (6.D0/25.D0)*A+(16.D0/25.D0)*B(I)
      CC0   = (24.D0/25.D0)*(A+B(I))
C
      KONST = (24.D0/5.D0)*DSQRT(PI/(A+G))*DSQRT(PI/(A+B(I)))
      EC1   = DSQRT(3.D0/2.D0)*KONST*KONST*KONST
      AC1   = (18.D0/25.D0)*A+(3.D0/25.D0)*G
      BC1   = (2.D0/25.D0)*A+(12.D0/25.D0)*G
      CC1   =-(12.D0/25.D0)*(A+G)
C
       DO 10 K=1,3
C
       W12C0 = 3.D0*PW(K)+3.D0*PM(K)
       E12C0 = DSQRT((A+G)/(A+G+4.D0*BK(K)))
       E12C0 = E12C0*E12C0*E12C0
       A12C0 = AC0
       B12C0 = BC0
       C12C0 = CC0
C
       W14C0 = 3.D0*PW(K)+3.D0*PM(K)
       E14C0 = DSQRT((A+G)/(A+G+(4.D0/3.D0)*BK(K)))
       E14C0 = E14C0*E14C0*E14C0
       A14C0 = (1.D0/25.D0)*((A+G)*(24.D0*A+9.D0*B(I))+BK(K)*(96.D0*A
     1        +64.D0*G+12.D0*B(I)))/((A+G)+(4.D0/3.D0)*BK(K))
       B14C0 = (1.D0/25.D0)*((A+G)*(6.D0*A+16.D0*B(I))+BK(K)*(24.D0*A
     1        +16.D0*G+(64.D0/3.D0)*B(I)))/((A+G)+(4.D0/3.D0)*BK(K))
       C14C0 = (1.D0/25.D0)*(24.D0*(A+G)*(A+B(I))+BK(K)*(96.D0*A
     1        +64.D0*G+32.D0*B(I)))/((A+G)+(4.D0/3.D0)*BK(K))
C
       W15C0 = 3.D0*PW(K)-2.D0*PM(K)+PB(K)-3.D0*PH(K)
       E15C0 = DSQRT((A+G)/(A+G+(4.D0/3.D0)*BK(K)))
       E15C0 = E15C0*E15C0*E15C0
       A15C0 = (1.D0/25.D0)*((A+G)*(24.D0*A+9.D0*B(I))+BK(K)*(36.D0*A
     1        +4.D0*G+12.D0*B(I)))/((A+G)+(4.D0/3.D0)*BK(K))
       B15C0 = (1.D0/25.D0)*((A+G)*(6.D0*A+16.D0*B(I))+BK(K)*(24.D0*A
     1        +16.D0*G+(64.D0/3.D0)*B(I)))/((A+G)+(4.D0/3.D0)*BK(K))
       C15C0 = (1.D0/25.D0)*(24.D0*(A+G)*(A+B(I))+BK(K)*(16.D0*A
     1        -16.D0*G+32.D0*B(I)))/((A+G)+(4.D0/3.D0)*BK(K))
C
       W45C0 = PW(K)+PM(K)+PB(K)+PH(K)
       E45C0 = 1.D0
       A45C0 = (24.D0/25.D0)*A+(9.D0/25.D0)*B(I)+(36.D0/25.D0)*BK(K)
       B45C0 = (6.D0/25.D0)*A+(16.D0/25.D0)*B(I)+(64.D0/25.D0)*BK(K)
       C45C0 = (24.D0/25.D0)*A+(24.D0/25.D0)*B(I)+(96.D0/25.D0)*BK(K)
C
       W12C1 =-2.D0*PW(K)-2.D0*PM(K)+PB(K)+PH(K)
       E12C1 = DSQRT((A+G)/(A+G+BK(K)))
       E12C1 = E12C1*E12C1*E12C1
       A12C1 = (1.D0/25.D0)*((A+G)*(18.D0*A+3.D0*G)+BK(K)*(27.D0*A
     1        +12.D0*G))/(A+G+BK(K))
       B12C1 = (1.D0/25.D0)*((A+G)*(2.D0*A+12.D0*G)+BK(K)*(38.D0*A
     1        +48.D0*G))/(A+G+BK(K))
       C12C1 =-(1.D0/25.D0)*(12.D0*(A+G)*(A+G)+BK(K)*(48.D0*A
     1        +48.D0*G))/(A+G+BK(K))
C
       W14C1 =-2.D0*PW(K)+3.D0*PM(K)-3.D0*PB(K)+PH(K)
       E14C1 = DSQRT((A+B(I))/(A+B(I)+BK(K)))
       E14C1 = E14C1*E14C1*E14C1
       A14C1 = (1.D0/25.D0)*((A+B(I))*(18.D0*A+3.D0*G)+BK(K)*(27.D0*A
     1        +3.D0*G+9.D0*B(I)))/(A+B(I)+BK(K))
       B14C1 = (1.D0/25.D0)*((A+B(I))*(2.D0*A+12.D0*G)+BK(K)*(18.D0*A
     1        +12.D0*G+16.D0*B(I)))/(A+B(I)+BK(K))
       C14C1 =-(1.D0/25.D0)*(12.D0*(A+B(I))*(A+G)+BK(K)*(-12.D0*A
     1        +12.D0*G-24.D0*B(I)))/(A+B(I)+BK(K))
C
       W23C1 =-PW(K)-PM(K)-PB(K)-PH(K)
       E23C1 = DSQRT((A+G)/(A+G+4.D0*BK(K)))
       E23C1 = E23C1*E23C1*E23C1
       A23C1 = AC1
       B23C1 = BC1
       C23C1 = CC1
C
       W24C1 =-4.D0*PW(K)-4.D0*PM(K)+2.D0*PB(K)+2.D0*PH(K)
       KONST = (A+G)*(A+B(I))
       KONSTS= (A+G)*(A+B(I))+BK(K)*(2.D0*A+G+B(I))
       E24C1 = DSQRT(KONST/KONSTS)
       E24C1 = E24C1*E24C1*E24C1
       KONST = (A+G)*(A+B(I))*(18.D0*A+3.D0*G)+BK(K)*(72.D0*A*A+
     1         3.D0*G*G+60.D0*A*G+54.D0*A*B(I)+39.D0*G*B(I))
       KONSTS= (A+G)*(A+B(I))+BK(K)*(2.D0*A+G+B(I))
       A24C1 =(1.D0/25.D0)*(KONST/KONSTS)
       KONST = (A+G)*(A+B(I))*(2.D0*A+12.D0*G)+BK(K)*(8.D0*A*A+
     1         12.D0*G*G+30.D0*A*G+6.D0*A*B(I)+16.D0*G*B(I))
       B24C1 = (1.D0/25.D0)*(KONST/KONSTS)
       KONST = 12.D0*(A+G)*(A+G)*(A+B(I))+BK(K)*(48.D0*A*A+
     1         12.D0*G*G+60.D0*A*G+36.D0*A*B(I)+36.D0*G*B(I))
       C24C1 =-(1.D0/25.D0)*(KONST/KONSTS)
C
       W45C1 =-PW(K)-PM(K)-PB(K)-PH(K)
       E45C1 = DSQRT((A+B(I))/(A+B(I)+4.D0*BK(K)))
       E45C1 = E45C1*E45C1*E45C1
       A45C1 = AC1
       B45C1 = BC1
       C45C1 = CC1
C-------------------------------------------------C
C CENTRAL NN-POTENTIAL: NON-LOCAL CONTRIBUTIONS   C
C-------------------------------------------------C
        DO 11 M=1,STSANZ
            X = DBLE(M)*HSTEP
         DO 11 N=1,STSANZ
            Y = DBLE(N)*HSTEP
C
          KONST1 = W12C0*E12C0*S(A12C0,B12C0,-C12C0,L,X,Y)+
     1             W14C0*E14C0*S(A14C0,B14C0,-C14C0,L,X,Y)+
     2             W15C0*E15C0*S(A15C0,B15C0,-C15C0,L,X,Y)+
     3             W45C0*E45C0*S(A45C0,B45C0,-C45C0,L,X,Y)
C
          KONST2 = W12C1*E12C1*S(A12C1,B12C1,-C12C1,L,X,Y)+
     1             W14C1*E14C1*S(A14C1,B14C1,-C14C1,L,X,Y)+
     2             W23C1*E23C1*S(A23C1,B23C1,-C23C1,L,X,Y)+
     3             W24C1*E24C1*S(A24C1,B24C1,-C24C1,L,X,Y)+
     4             W45C1*E45C1*S(A45C1,B45C1,-C45C1,L,X,Y)
C
          WERT  = -V0(K)*(EC0*KONST1+EC1*KONST2)/NO
C
          DO 12 K1=1,DIS
            MAT(M,N,K1,1) = MAT(M,N,K1,1)+AVOR(I,K1)*WERT
   12     CONTINUE
C
   11   CONTINUE
   10 CONTINUE
C
      RETURN
      END
C*****************************C
C  COULOMB POTENTIAL KERNELS  C
C*****************************C
      SUBROUTINE COUPOT(DIS,FRAG,L,NO)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,K1,DIS,SKIP,FRAG,NMAX,KMAX,STSANZ
      PARAMETER ( KMAX =  3 , NMAX = 70 ,  SKIP = 1 )
      REAL*8    AVOR(3,3),B(3),MAT(NMAX,NMAX,KMAX,2)
C
      COMMON /KONSTA/ PI,H2M,ELEML
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C---------------------------------------C
C FOR SKIP=1: CALCULATION OF THE        C
C NON-LOCAL COULOMB-KERNELS             C
C---------------------------------------C
      IF(SKIP.EQ.0) RETURN
C
      DO 10 I=1,3
C
       KONST = (24.D0/5.D0)*DSQRT(4.D0/3.D0)*(PI/(A+G))
       EC0   = DSQRT(3.D0/2.D0)*KONST*KONST*KONST
       AC0   = (24.D0/25.D0)*A+(9.D0/25.D0)*B(I)
       BC0   = (6.D0/25.D0)*A+(16.D0/25.D0)*B(I)
       CC0   = (24.D0/25.D0)*(A+B(I))
C
       KONST = (24.D0/5.D0)*DSQRT(PI/(A+G))*DSQRT(PI/(A+B(I)))
       EC1   = DSQRT(3.D0/2.D0)*KONST*KONST*KONST
       AC1   = (18.D0/25.D0)*A+(3.D0/25.D0)*G
       BC1   = (2.D0/25.D0)*A+(12.D0/25.D0)*G
       CC1   =-(12.D0/25.D0)*(A+G)
C------------------------------C
C  NON-LOCAL COULOMB POTENTIAL C
C------------------------------C
       IF(FRAG.EQ.2) THEN
           P12C0 = DSQRT((A+G)/PI)
           P15C0 = DSQRT((3.D0/4.D0)*(A+G))
           L3C   = (2.D0/5.D0)*DSQRT((3.D0/4.D0)*(A+G))
           M3C   = 2.D0*L3C
C
           P12C1 = DSQRT(A+G)
           L4C   = (3.D0/5.D0)*DSQRT(A+G)
           M4C   = 2.D0*L4C
C
           P15C1 = DSQRT(A+B(I))
           L5C   = (3.D0/5.D0)*DSQRT(A+B(I))
           M5C   =-(4.D0/3.D0)*L5C
C
           P24C1 = DSQRT(((A+G)*(A+B(I)))/(2.D0*A+G+B(I)))
           L2C   = (6.D0/5.D0)*DSQRT(((A+G)*(A+B(I)))/(2*A+G+B(I)))
           M2C   = (1.D0/3.D0)*L2C
                     ELSE
           P14C0 = DSQRT((3.D0/4.D0)*(A+G))
           L1C   = (8.D0/5.D0)*DSQRT((3.D0/4.D0)*(A+G))
           M1C   =-(1.D0/2.D0)*L1C
C
           P24C1 = DSQRT(((A+G)*(A+B(I)))/(2.D0*A+G+B(I)))
           L2C   = (6.D0/5.D0)*DSQRT(((A+G)*(A+B(I)))/(2*A+G+B(I)))
           M2C   = (1.D0/3.D0)*L2C
       ENDIF
C
       DO 11 M=1,STSANZ
           X = DBLE(M)*HSTEP
        DO 11 N=1,STSANZ
           Y = DBLE(N)*HSTEP
C
        IF(FRAG.EQ.2) THEN
           TERM1  = 2.D0*P15C0*COULI(AC0,BC0,L3C,M3C,-CC0,L,X,Y,1)
C
           TERM2  =      P12C1*COULI(AC1,BC1,L4C,M4C,-CC1,L,X,Y,1)+
     1                   P15C1*COULI(AC1,BC1,L5C,M5C,-CC1,L,X,Y,1)+
     2                   P24C1*COULI(AC1,BC1,L2C,M2C,-CC1,L,X,Y,1)
C
           TERM3  =      P12C0*S(AC0,BC0,-CC0,L,X,Y)
C
           TERM   = 2.D0*PI*(EC0*TERM1-EC1*TERM2)+EC0*TERM3
                      ELSE
C
           TERM   = 2.D0*PI*(EC0*P14C0*COULI(AC0,BC0,L1C,M1C,-CC0,L,X,
     1              Y,1)-EC1*P24C1*COULI(AC1,BC1,L2C,M2C,-CC1,L,X,Y,1))
        ENDIF
C
         WERT = ELEML*ELEML*TERM/NO
C
         DO 12 K1=1,DIS
           MAT(M,N,K1,1) = MAT(M,N,K1,1)+AVOR(I,K1)*WERT
   12    CONTINUE
C
   11  CONTINUE
   10 CONTINUE
C
      RETURN
      END
C***************C
C LS-KERNELS    C
C***************C
      SUBROUTINE SPIORB(DIS,JL,L,SBRA,NO)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N,K1,JL,DIS,SBRA,KMAX,NMAX,STSANZ
      PARAMETER ( KMAX =  3 , NMAX = 70 )
      REAL*8    AVOR(3,3),B(3),V0(3),PW(3),PM(3),PH(3),PB(3),BK(3),
     1          MAT(NMAX,NMAX,KMAX,2)
C
      COMMON /KONSTA/ PI,H2M,ELEML
      COMMON /BREITE/ A,G,B,AVOR
      COMMON /POTPAR/ V0,PW,PM,PB,PH,BK,VO,LAM
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      IF((VO.EQ.0.D0).OR.(L.EQ.0)) RETURN
C
      DO 10 I=1,3
C
       DLS1=-2.D0*(12.D0/5.D0)**4*(PI/(A+G))**(3.D0/2.D0)*(PI/(A
     1      +B(I)))**(3.D0/2.D0)*(A+G)*((A+B(I))/(A+B(I)+LAM))**
     2       (5.D0/2.D0)
       DLS2= 2.D0**6*(12.D0/5.D0)**4*(PI/(A+G))**3*(A+B(I))*((A+
     1       G)/(4.D0*LAM+3.D0*A+3.D0*G))**(5.D0/2.D0)
       ALS1= 3.D0*(3.D0*B(I)*LAM+6.D0*B(I)*A+B(I)*G+9.D0*LAM*A+LAM*G
     1      +6.D0*A**2+A*G)/(25.D0*(B(I)+LAM+A))
       BLS1= 2.D0*(8.D0*B(I)*LAM+B(I)*A+6.D0*B(I)*G+9.D0*LAM*A+6.D0*
     1       LAM*G+A**2+6.D0*A*G)/(25.D0*(B(I)+LAM+A))
       CLS1=-12.D0*(2.D0*B(I)*LAM-B(I)*A-B(I)*G+LAM*A-LAM*G-A**2-A*G)
     1       /(25.D0*(B(I)+LAM+A))
       ALS2= 3.D0*(12.D0*B(I)*LAM+9.D0*B(I)*A+9.D0*B(I)*G+36.D0*LAM*A
     1      +4.D0*LAM*G+24.D0*A**2+24.D0*A*G)/(25.D0*(4.D0*LAM+3.D0*A
     2      +3.D0*G))
       BLS2= 2.D0*(32.D0*B(I)*LAM+24.D0*B(I)*A+24.D0*B(I)*G+36.D0*LAM
     1      *A+24.D0*LAM*G+9.D0*A**2+9.D0*A*G)/(25.D0*(4.D0*LAM+3.D0*
     2       A+3.D0*G))
       CLS2=-24.D0*(4.D0*B(I)*LAM+3.D0*B(I)*A+3.D0*B(I)*G+2.D0*LAM*A-
     1       2.D0*LAM*G+3.D0*A**2+3.D0*A*G)/(25.D0*(4.D0*LAM+3.D0*A+G
     2      *3.D0))
C---------------------------------------C
C SPIN-ORBIT POTENTIAL: NON-LOCAL PART  C
C---------------------------------------C
       DO 11 M=1,STSANZ
            X = DBLE(M)*HSTEP
        DO 11 N=1,STSANZ
            Y = DBLE(N)*HSTEP
C
        IF(SBRA.EQ.1) THEN
         KON = 10.D0/DSQRT(6.D0)*DLS1/CLS1*S(ALS1,BLS1,CLS1,L,X,Y)
     1         -8.D0/DSQRT(6.D0)*DLS2/CLS2*S(ALS2,BLS2,CLS2,L,X,Y)
                      ELSE
         KON =-DSQRT(2.D0/3.D0)*DLS1/CLS1*S(ALS1,BLS1,CLS1,L,X,Y)
     1        -DSQRT(2.D0/3.D0)*DLS2/CLS2*S(ALS2,BLS2,CLS2,L,X,Y)
        ENDIF
C
         WERT =-VO*JLS(JL,L,SBRA)*KON/NO
C
         DO 12 K1=1,DIS
           MAT(M,N,K1,1) = MAT(M,N,K1,1)+AVOR(I,K1)*WERT
   12    CONTINUE
C
   11  CONTINUE
   10 CONTINUE
C
      RETURN
      END
C**************************************C
C  WRITE MATRIX MAT ON A FILE          C
C**************************************C
      SUBROUTINE WRITIN(DIS)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4    I,J,K,L,M,N,K1,IN,NIN,DIS,KMAX,NMAX,STSANZ
      PARAMETER    ( KMAX =  3 , NMAX = 70 )
      REAL*8       MAT(NMAX,NMAX,KMAX,2)
      CHARACTER*13 OUTNAME(3)
C
      COMMON /DISKPA/ HSTEP,STSANZ
      COMMON /MATRIX/ MAT
C
      IN = STSANZ+1
C
      OUTNAME(1) = 'COUPLE.OU1'
      OUTNAME(2) = 'COUPLE.OU2'
      OUTNAME(3) = 'COUPLE.OU3'
C
       DO 10 K1=1,DIS
          NIN = 20+K1
          OPEN(NIN,FILE=OUTNAME(K1),STATUS='NEW',FORM='UNFORMATTED')
C
        DO 11 M=1,IN
           WRITE(NIN)  0.D0
   11   CONTINUE

        DO 12 M=1,STSANZ
         DO 12 N=1,STSANZ
           WRITE(NIN)   MAT(M,N,K1,1),MAT(M,N,K1,2)
   12   CONTINUE
C
   10  CONTINUE
       CLOSE(NIN)
C
      RETURN
      END
C------------------------------------C
C  J=1 (J = L+1/2) J=-1 (J = L-1/2)  C
C------------------------------------C
      DOUBLE PRECISION FUNCTION JLS(JL,L,SBRA)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 L,I,JL,SBRA
C
      JLS  = 0.D0
      DREH = DBLE(L)
C
       IF(SBRA.EQ.2) THEN
                               JLS= DSQRT(DREH*(2.D0*DREH+3.D0))
                  IF(JL.EQ.-1) JLS= DSQRT((2.D0*DREH-1.D0)*(DREH+1.D0))
                     ELSE
                               JLS= DREH
                  IF(JL.EQ.-1) JLS=-(DREH+1.D0)
       ENDIF
C
      RETURN
      END
