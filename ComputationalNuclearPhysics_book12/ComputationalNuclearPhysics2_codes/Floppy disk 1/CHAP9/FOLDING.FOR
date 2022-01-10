C     PROGRAM FOLDING
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /XYW/ XZ(64),WWZ(64),XX(96),WWX(96),NZ,NXX
c     COMMON /XYW/ X1(64),W1(64),X2(96),W2(96),N1,N2
      DIMENSION YZ(64),WZ(64),YX(200),WX(200)
      DIMENSION UR(50)
      DATA PI /3.141592653589793D0/
c
C     >Standard Parameter<
      AP0=0.65D0
      AD=0.54D0
      VP0=-50.D0
      AP=AP0
      VP=VP0
C     >p-shell nuclei:<
      AP1=0.35D0
      VP1=-30.D0
C     >sd-shell nuclei<
      AP2=0.45D0
      VP2=-40.D0
C     >open output and input files<
      OPEN(UNIT=5,FILE='INPUT.DAT')
      OPEN(UNIT=25,FILE='FOLDING.DAT')
C     >input: gauss-legendre reference points<
      CALL DATIN
      NX=NXX
C
   20 READ (5,*,END=500) APOT,ADICHT,ZPOT,ZDICHT,ENERLS,LANF,LEND
      WRITE (6,*) 'e?',APOT,ADICHT,ZPOT,ZDICHT,ENERLS,LANF,LEND
      VP=VP0
      IF(APOT.LE.16) THEN
          AP=AP1
          VP=VP1
      ELSE
          IF(APOT.LE.40) THEN
              AP=AP2
              VP=VP2
          ENDIF
      ENDIF
      WRITE (25,31) AP,AD,VP
   31 FORMAT(3F10.3)
   30 FORMAT (5F10.5,2I10)
C     >the following again: target and projectile interchanged<
      DO 170 IRP=1,2
      IF (IRP .EQ. 2) GOTO 35
      WRITE (25,30) APOT,ADICHT,ZPOT,ZDICHT,ENERLS,LANF,LEND
   35 RPOT=1.25D0*APOT**(1.D0/3.D0)
      RDICHT=1.12D0*ADICHT**(1.D0/3.D0)-0.86D0*ADICHT**(-1.D0/3.D0)
      IR=APOT**(1.D0/3.D0)+ADICHT**(1.D0/3.D0)+0.5D0
      RA=DBLE(IR)-5.D0
      IF (RA .LE. 0.D0) RA=1.D0
      RE=RA+13.D0
      DR=0.5D0
C
C     >calculation of the mass-distribution<
      SRO=0.D0
      DO 11 I=1,NZ
      Y=(RDICHT+8.0)*0.5D0*(XZ(I)+1.0)
      W=(RDICHT+8.0)*0.5D0*WWZ(I)
   11 SRO=SRO+W*Y*Y/(1.D0+DEXP((Y-RDICHT)/AD))
      RONULL=ADICHT/(4*PI*SRO)
C
      WRITE (6,40) AP,AD,VP,RONULL,ADICHT,APOT
   40 FORMAT ('1AP=',F5.2,3X,'AD=',F5.2,3X,'VP=',F5.1,3X,
     1 'RONULL=',F7.3,3X,'ADICHT=',F6.2,3X,'APOT=',F6.2//)
C
      R=RA
      J=1
C
   50 DO 60 I=1,NZ
      YZ(I)=R*0.5D0*XZ(I)
   60 WZ(I)=R*0.5D0*WWZ(I)
      GA=R*0.5D0
      GB=GA+(RPOT+RDICHT)*0.5D0
      XMAX=R*0.5D0+2.*(RPOT+RDICHT)
      DO 70 I=1,NX
      II=NX+I
      YX(I)=(GB-GA)*0.5D0*XX(I)+(GA+GB)*0.5D0
      WX(I)=(GB-GA)*0.5D0*WWX(I)
      YX(II)=(XMAX-GB)*0.5D0*XX(I)+(GB+XMAX)*0.5D0
   70 WX(II)=(XMAX-GB)*0.5D0*WWX(I)
      NX=2*NX
      M=1
C
C     >folding-integration<
   80 SX=0.D0
      DO 100 IX=1,NX
      SZ=0.D0
      DO 90 IZ=1,NZ
      A1=(YX(IX)+YZ(IZ)-RPOT)/AP
      A2=(YX(IX)-YZ(IZ)-RDICHT)/AD
   90 SZ=SZ+WZ(IZ)*(YX(IX)**2-YZ(IZ)**2)/(1.+DEXP(A2))*(VP/(1.+DEXP(A1)
     1))
  100 SX=SX+WX(IX)*SZ
      U=4.*PI*RONULL/R*SX
      IF (M .EQ. 2) GOTO 120
      UA=U
      UR(J)=UA
      M=2
      NX=NX/2
      GB=XMAX+XMAX
C     >rescale the reference points:<
      DO 110 I=1,NX
      YX(I)=(GB-XMAX)*0.5D0*XX(I)+(XMAX+GB)*0.5D0
  110 WX(I)=(GB-XMAX)*0.5D0*WWX(I)
      GOTO 80
C
  120 IF (DABS(U) .GT. DABS(UA)/100.) WRITE (6,130)
  130 FORMAT (' INSUFFICIENT PRECISION')
      WRITE (6,140) R,GA,XMAX,UA,XMAX,GB,U
  140 FORMAT (' R=',F4.1,6X,'INTEGRAL FROM',F5.2,' TO',F6.2,' =',D15.7,
     1 6X,'INTEGRAL FROM',F6.2,' TO',F6.2,  ' =',D15.7)
C
      IF (DABS(RE-R) .LT. 1.D-3) GOTO 160
      R=R+DR
      J=J+1
      GOTO 50
C
C     >results u(r) into auxiliary file 25<
  160 IR=J
      WRITE (25,165) (UR(I),I=1,IR),RA
  165 FORMAT (5D15.7)
C
      HELP=APOT
      ZHELP=ZPOT
      APOT=ADICHT
      ZPOT=ZDICHT
      ADICHT=HELP
      ZDICHT=ZHELP
      AP=AP0
      VP=VP0
      IF(APOT.LE.16) THEN
          AP=AP1
          VP=VP1
      ELSE
          IF(APOT.LE.40) THEN
              AP=AP2
              VP=VP2
          ENDIF
      ENDIF
      WRITE (6,31) AP,AD,VP
170   CONTINUE
      GOTO 20
C
  500 CLOSE(UNIT=25)
      CLOSE(UNIT=5)
      STOP
      END
C
C
      SUBROUTINE DATIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /XYW/ X1(64),W1(64),X2(96),W2(96),N1,N2
c     >Gaus-Legendre reference points<
Ca)   (N1=40)
      READ (5,10) N1
      WRITE (6,*)'a?', N1
   10 FORMAT (I3)
      I1=N1/2+1
      READ (5,20) (X1(I),I=I1,N1),(W1(I),I=I1,N1)
      WRITE (6,*)'b?', (X1(I),I=I1,N1),(W1(I),I=I1,N1)
   20 FORMAT (4F18.15)
      I2=N1/2
      DO 30 I=1,I2
      IN=N1-I+1
      W1(I)=W1(IN)
   30 X1(I)=-X1(IN)
Cb)   (N2=20)
      READ (5,10) N2
      WRITE (6,*)'c?', N2
      I1=N2/2+1
      READ (5,20) (X2(I),I=I1,N2),(W2(I),I=I1,N2)
      WRITE (6,*)'d?', (X2(I),I=I1,N2),(W2(I),I=I1,N2)
      I2=N2/2
      DO 40 I=1,I2
      IN=N2-I+1
      W2(I)=W2(IN)
   40 X2(I)=-X2(IN)
C
      RETURN
      END
