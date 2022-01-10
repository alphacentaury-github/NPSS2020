C     PROGRAM FRICT
C     >standard with sticking l  and rkutdf; trajectories<
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION MASS,KRR,KFIFI
      DIMENSION TIME(600),Y(4),FUNKA(600,4),ENVER(600),TTETA(600),
     1DREVER(600),DENER(600),HF(12)
      EXTERNAL FFRICT
      COMMON MASS,ALFA,ANGM,RFRICT,RWOOD,RCOUL,A1,A2,A3,A4,A5,AG,
C    X VRAD,VTAN,ARAD,DIFPOT
     X VRAD,VTAN,ARAD,KRR,KFIFI,DIFPOT
      PI = 3.14159265358979D0
c     >friction parameter: VRAD=strength of radial friction, VTAN=tangential:<
      VRAD=4.0
      VTAN=0.01
      ARAD=0.3
C     >open input file<
      OPEN(UNIT=5,FILE='FIT.DAT')
C***  >read input data<
1     READ(5,2,END = 29) ATAR,APRO,PTAR,PPRO,ENERLS,LANF,LEND
2     FORMAT (5(F10.5),2(I10))
c     >potential parameters:<
      READ (5, 4) AG,A1,A2,A3,A4,A5,RNULL
4     FORMAT(7(F10.4))
      WRITE(6,3)   ATAR,APRO,PTAR,PPRO,ENERLS,LANF,LEND
3     FORMAT ('1',5(F10.4),2(I10))
      LDELTA = 10
      LANGM = LEND
      WRITE(6, 4) AG,RNULL,A1,A2,A3,A4,A5
      WRITE (6,1001) VRAD,VTAN,ARAD
1001  FORMAT (5X,'VRAD = ',F10.4,5X,'VTAN = ',F10.4,5X,'ARAD = ',F10.4)
c     >MASS=reduced mass<
      MASS= APRO*ATAR*931.4812D0/((APRO+ATAR)*2.997925D0*2.997925D0)
      ALFA = PPRO*PTAR*1.9732891D2/1.3703602D2
C***  >change of step size of angular momentum<
c     >L=angular momentum<
      LMARK = 0
c     >energy in cm system<
      ENER = MASS*2.997925D0*2.997925D0*ENERLS/(APRO*931.4812D0)
c     >radius of potential<
      RWOOD = RNULL*(APRO**(1./3.D0)+ATAR**(1./3.D0))
c     >radius of density<
      RCOUL = 1.12*ATAR**(1./3.D0)-0.86*ATAR**(-1./3.D0)
     + +1.12*APRO**(1./3.D0)-0.86*APRO**(-1./3.D0)
c     >friction radius<
      RFRICT = RCOUL + 2.7
      A11=APRO**(5./3.D0)
      A22=ATAR**(5./3.D0)
c     >rigid body moment of inertia<
      CONST=(2./5.)*931.4812*(A11+A22)*1.12*1.12/(2.9979**2)
      AH = ALFA/(2.D0*ENER)
      DELTAE=1.D0
 6    ANGM = 6.582183D1*DBLE(LANGM)
      EPSI = DSQRT(1.D0+2.D0*ENER*ANGM*ANGM/(MASS*ALFA*ALFA))
      SPER = ANGM*ANGM/(MASS*ALFA)
      STEPP = 1.0
      R2R = 0.D0
c     >starting distance for trajectory<
      Y(1) = 2.5D0*RFRICT
      YY2Y=(2.D0* (ENER-ALFA/Y(1))-ANGM*ANGM/(MASS*Y(1)*Y(1)))
     X /MASS
      IF(YY2Y.LE.0.0) THEN
        WRITE(*,*) 'INGOING ENERGY BELOW BARRIER AT 2.5*(RFRICT=',
     X    RFRICT,')','    ECM=',ENER,'PROJECTILE-TARGET INTERCHANGED ?'
            GOTO 1
      ENDIF
c     >initial radial velocity:<
      Y(2)=-DSQRT(YY2Y)
      Y(3) =-DATAN(DSQRT((EPSI*Y(1)/(SPER+Y(1)))**2.D0-1.D0))
c     >initial tangential velocity:<
      Y(4) = ANGM/(MASS*Y(1)*Y(1))
      AW = -DSQRT(MASS/(2.D0*ENER))*(DSQRT((Y(1)-AH)**2.D0-(AH*EPSI)
     X **2.D0)+AH*DLOG((Y(1)-AH)/(AH*EPSI)+DSQRT(((Y(1)-AH)/(AH*EPSI))
     X **2.D0-1.D0)))
      WRITE(6,10)
10    FORMAT (/7X,'MASS',16X,'ALFA',17X,'ANGM',17X,'ENER')
C10    FORMAT (/7X,'MASS',16X,'ALFA',17X,'ANGM',17X,'ENER',15X,
C     X ' AW ',16X,'STEP')
      WRITE(6,11)  MASS,ALFA,ANGM,ENER,Y(1),Y(2),Y(3),Y(4)
C      WRITE(6,11)  MASS,ALFA,ANGM,ENER,AW,STEPP,Y(1),Y(2),Y(3),Y(4),
C     X EPSI,SPER
11    FORMAT (4(D20.12)/4(D20.12))
C11    FORMAT (6(D20.12)/6(D20.12))
      WRITE(6,12)
12    FORMAT (6X,' Y(1) ',15X,' Y(2) ',15X,' Y(3) ',15X,' Y(4) '/)
C12    FORMAT (6X,' Y(1) ',15X,' Y(2) ',15X,' Y(3) ',15X,' Y(4) ',15X,
C     X 'EPSI',15X,'SPER'/)
      RMIN = 100.0
      DO 13 I = 1,600
      TIME(I) = AW
      DO 14 J = 1,4
14    FUNKA(I,J) = Y(J)
      M = I
C***  >M gives the number of steps of integration along the trajectory<
C
C***  >calculates the loss of energy and angular momentum<
      ANGLOC = MASS*Y(1)*Y(1)*Y(4)
      DREVER(I) = (ANGM-ANGLOC)/6.582183D1
      ANGSTM=(MASS*Y(1)*Y(1)*ANGM)/(CONST+MASS*Y(1)*Y(1))
Cd      IF (ANGLOC.LT.ANGSTM) GOTO 111
cd 111  VTAN=0.D0
      YYY = Y(1)-RWOOD
      YYE = 0.D0
      IF(YYY/AG.LE.30.D0) YYE = DEXP(-YYY/AG)
      ENERGY = MASS*(Y(2)*Y(2)+Y(1)*Y(1)*Y(4)*Y(4))/2.D0
     X -DLOG(1.D0+YYE)*(A1+A2*YYY+A3*YYY*YYY+A4*YYY**3+A5*YYY**4)
      IF(Y(1).GE.RCOUL)GOTO15
      ENERGY    = ENERGY   +3.D0*ALFA*(1.D0-Y(1)*Y(1)/(3.D0*RCOUL
     X *RCOUL))/(2.D0*RCOUL)
      GOTO16
15    ENERGY    = ENERGY   +ALFA/Y(1)
16    ENVER(I) = ENER-ENERGY
      IF(Y(1).GT.RMIN)GOTO17
      ECLOSE = ENERGY
      DRCLOS  = ANGLOC
      RMIN = Y(1)
      GOTO18
17    CONTINUE
18    IF(Y(1).GE.2.6*RFRICT)GOTO19
      IF(ENERGY.LE.0.D0)GOTO19
      ANGML=ANGLOC/65.82183
      ANGSTH=ANGSTM/65.82183
      WRITE (6,374) Y(1),Y(2),AW,ENERGY,ENVER(I),ANGML,
     1DREVER(I),KRR,KFIFI
 374  FORMAT (9D15.6)
C 374  FORMAT (8D15.6)

C     >integration of equs: of motion<
      CALL DIF1DF(4,FFRICT,1.D-8,STEPP,AW,Y)
C     CALL RKUTDF(4,STEPP,AW,Y,HF,FFRICT)
C     WRITE (6,374) Y(1),AB(1),AB(2),AB(3)
C374  FORMAT (4F15.6)
13    CONTINUE
19    ENER2R = ENERGY
      IM=I-1
      DO 300 KI=2,IM
      IF (ENVER(KI).GT.DELTAE)  GOTO 333
 300  CONTINUE
      GOTO 311
C     >position where the energy-loss is > DELTAE for the first time:<
 333  DE1=(DELTAE-ENVER(KI-1))/(ENVER(KI)-ENVER(KI-1))
      TNULL= DE1*(TIME(KI)-TIME(KI-1))+TIME(KI-1)
      RN= DE1*(FUNKA(KI,1)-FUNKA(KI-1,1))+FUNKA(KI-1,1)
      KI1=KI+1
      IF (FUNKA(KI,2).GT.0.D0) GOTO 311
      DO 301 II=KI1,I
      IF(FUNKA(II,1).GT.RN) GOTO 302
  301 CONTINUE
  302 TNULL1=(RN-FUNKA(II-1,1))/(FUNKA(II,1)-FUNKA(II-1,1))
      TNULLQ=TNULL1*(TIME(II)-TIME(II-1))+TIME(II-1)
      DELTAT=TNULLQ-TNULL
      WRITE (6,303)  FUNKA(KI-1,1),TIME(KI-1),
     1  FUNKA(KI,1),TIME(KI),
C    2  FUNKA(KI+1,1),TIME(KI+1),
     3  FUNKA(II-1,1),TIME(II-1),
     4  FUNKA(II,1),TIME(II),
C    5  FUNKA(II+1,1),TIME(II+1),
     6RN,TNULL,TNULLQ,DELTAT,DE1
 303  FORMAT (5X,'R1=',F10.5,5X,'T1=',F10.5/
     15X,'R0=',F10.5,5X,'T0=',F10.5/
C    25X,'R2=',F10.5,5X,'T2=',F10.5/
     35X,'R3=',F10.5,5X,'T3=',F10.5/
     45X,'R0=',F10.5,5X,'T0=',F10.5/
C    55X,'R4=',F10.5,5X,'T4=',F10.5/
     65X,'RM=',F10.5,5X,'TNULL=',F10.5,5X,'TNULLQ=',F10.5
     7,5X,'DELTAT=',F10.5,5X,F10.7)
 311  CONTINUE
      R2R = Y(1)
      FI2R = Y(3)
      DRE2R = ANGLOC
      DRE2RH = DRE2R/6.582183D1
      SPER2R =  DRE2R*DRE2R/(MASS*ALFA)
C
C***  >start end stop value of trajectory<
      DENER(LANGM) = ENVER(M)
      LCLOSE = DRCLOS /6.582183D1
      WRITE(6,20)LANGM,FUNKA(1,1),FUNKA(1,3),ENER,DRE2RH,FUNKA(M,1),
     X FUNKA(M,3),ENER2R,DREVER(M),ENVER(M),RMIN ,LCLOSE,ECLOSE
20    FORMAT (5X,'STARTVALUES  = ',5X,'L/H = ',I6,5X,' R = ',D15.6,5X,
     X 'FI = ',D15.6,5X,'ENERGY  = ',D15.6/5X,'FINALVALUE:',9X,'L/H = ',
     X D10.3,2X,'R = ',D15.6,5X,'FI = ',D15.6,5X,'ENERGY  = ',D15.6/
     X 5X,'LOSS OF ANGULAR MOMENTUM/H = ',D15.6,5X,'LOSS OF ENERGY= ',
     X D15.6/5X,'RMIN = ', D15.5,5X,'LCLOSE = ',I6,5X,'ECLOSE = ',
     X D15.6/)
C20   FORMAT (5X,'START VALUE = ',5X,'L/H = ',I6,5X,' R = ',D15.6,5X,
C
C***  >calculation of scattering angle<
      IF(R2R.LT.2.D0*RFRICT)GOTO23
      EPSIL = DSQRT(1.D0+2.D0*ENER2R*DRE2R*DRE2R/(MASS*ALFA*ALFA))
      FIVOR = DATAN(DSQRT(EPSI*EPSI-1.D0))
      FINACH = DATAN(DSQRT(EPSIL*EPSIL-1.D0))
      FI = +FI2R-DATAN(DSQRT(EPSIL*EPSIL*R2R*R2R/((SPER2R+R2R)*
     X (SPER2R+R2R))-1.D0))
      TETA = PI - FINACH-FIVOR-FI
      TFINAC  = FINACH*1.8D2/PI
      TFIVOR = FIVOR*1.8D2/PI
      TFI = FI*1.8D2/PI
C     WRITE (6,22)  TFIVOR,TFINAC ,TFI
c22   FORMAT (5X,3(D15.6))
      TTETA(LANGM) = TETA*1.8D2/PI
      GOTO24
23    TTETA(LANGM) = 300.D0
24    WRITE(6,25) TTETA(LANGM),M
25    FORMAT('+',54X,'TETA = ',D15.6,10X,'M = ',I5//1X,120('*'))
C
      IF(TTETA(LANGM).EQ.300.D0)GOTO32
      LANGM = LANGM - LDELTA
      IF(LANGM.LT.LANF) GOTO32
      GOTO 6
 32   IF (LMARK.EQ.1) GOTO 27
      LMARK = 1
      LDELTA = 2
      LANGM = LANGM+30
C     LANGM=LANGM+LDELTA
      GOTO6
C
C***  >writing of scattering angle<
27    LANGM=LANGM+LDELTA/2+1
      LANF = MAX0(LANF,LANGM)
      WRITE (6,28)  (J,TTETA(J),DENER(J),J = LANF,LEND,LDELTA)
28    FORMAT(//19X,'L',19X,'TETA(L)',21X,'ENERGYLOSS(L)    '
     X//15(14X,I6,10X,F15.6,17X,F15.6/))
      RLAMD2 = 1.9732891D2**2/(2.D0*MASS*ENER*2.997925D0**2)
C28   FORMAT(//19X,'L',19X,'TETA(L)',21X,'ENERGYLOSS(L)'
C     WRITE (99,33) ATAR,APRO,ENERLS,RLAMD2
c99   FORMAT (3F10.5,F10.8)
      SIGCF = PI*RLAMD2*LANF*LANF
      WRITE (6,30) RLAMD2,SIGCF
   30 FORMAT (/2X,'LAMDA**2 = ',F10.8,5X,'SIGCF (FM**2) = ',F12.2)
      DO 285 K=LANF,LEND,LDELTA
      TTETA(K) = 0.D0
 285  CONTINUE
      GOTO1
 29   CLOSE(UNIT=5)
      STOP
      END
C
      SUBROUTINE FFRICT(AW,A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MASS,KRR,KFIFI
      DIMENSION A(4),B(4)
      COMMON MASS,ALFA,ANGM,RFRICT,RWOOD,RCOUL,A1,A2,A3,A4,A5,AG,
C    X VRAD,VTAN,ARAD,DIFPOT
     X VRAD,VTAN,ARAD,KRR,KFIFI,DIFPOT
      IF(A(1).LE.3.5*RFRICT.AND.A(1).GE.0.0)GOTO3
      GOTO4
3     KRR = 0.D0
      AA = A(1)-RWOOD
      AAA = DEXP(-AA/AG)
      DIFPOT = -DLOG(1.0+AAA)*(A2+2.*A3*AA+3.*A4*AA*AA+4.*A5*AA*AA*AA)
     X +(A1+A2*AA+A3*AA*AA+A4*AA**3+A5*AA**4)*(AAA/(AG*(1.+AAA)))
      AA=DIFPOT**2
      KRR = VRAD*AA
      KFIFI = VTAN*AA
      IF(A(1).GE.RCOUL)GOTO1
      DIFPOT = -ALFA*A(1)/(RCOUL*RCOUL*RCOUL)+DIFPOT
      GOTO2
1     DIFPOT = -ALFA/(A(1)*A(1))+DIFPOT
2     CONTINUE
      B(1) = A(2)
      B(2) = A(1)*A(4)*A(4)-(DIFPOT+A(2)*KRR)/MASS
      B(3) = A(4)
      B(4) = -2.D0*A(2)*A(4)/A(1)-A(4)*KFIFI/MASS
4     RETURN
      END
C
      SUBROUTINE DIF1DF (N,F,EPS,H,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Y(N),YA(22),YL(22),YM(22),DY(22),DZ(22),DT(22,7),
     1D(7),S(22),X,XN,H,G,B,B1,U,V,C,TA,W,EPS
      DIMENSION EP(4)
      LOGICAL KONV,BO,KL,GR
      DATA EP/0.4D-1,0.16D-2,0.64D-4,0.256D-5/
      JTJ=0
      JTI=0
      FY=1D0
      ETA=ABS(EPS)
      IF(ETA.LT.1.D-11) ETA=1.D-11
      DO 100 I=1,N
  100 YA(I)=Y(I)
      CALL F (X,Y,DZ)
   10 XN=X+H
      BO=.FALSE.
      DO 110 I=1,N
  110 S(I)=0.D0
      M=1
      JR=2
      JS=3
      DO 260 J=1,10
      IF(.NOT.BO) GOTO 200
      D(2)=1.777777777777778D0
      D(4)=7.111111111111111D0
      D(6)=2.844444444444444D1
      GOTO 201
  200 D(2)=2.25D0
      D(4)=9.D0
      D(6)=3.6D1
  201 IF(J.LE.7) GOTO 202
      L=7
      D(7)=6.4D1
      GOTO 203
  202 L=J
      D(L)=M*M
  203 KONV=L.GT.3
      M=M+M
      G=H/DBLE(M)
      B=G+G
      M=M-1
      DO 210 I=1,N
      YL(I)=YA(I)
  210 YM(I)=YA(I)+G*DZ(I)
      DO 220 K=1,M
      CALL F (X+DBLE(K)*G,YM,DY)
      DO 220 I=1,N
      U=YL(I)+B*DY(I)
      YL(I)=YM(I)
      YM(I)=U
      U=ABS(U)
      IF(U.GT.S(I)) S(I)=U
  220 CONTINUE
      CALL F (XN,YM,DY)
      KL=L.LT.2
      GR=L.GT.5
      FS=0.
      DO 233 I=1,N
      V=DT(I,1)
      C=(YM(I)+YL(I)+G*DY(I))*0.5D0
      DT(I,1)=C
      TA=C
      IF(KL) GOTO 233
      DO 231 K=2,L
      B1=D(K)*V
      B=B1-C
      W=C-V
      U=V
      IF(B.EQ.0.D0) GOTO 230
      B=W/B
      U=C*B
      C=B1*B
  230 V=DT(I,K)
      DT(I,K)=U
  231 TA=U+TA
      IF(.NOT.KONV) GOTO 232
      IF(ABS(Y(I)-TA).GT.S(I)*ETA) KONV=.FALSE.
  232 IF(GR.OR.S(I).EQ.0.D0) GOTO 233
      FV=ABS(W)/S(I)
      IF(FS.LT.FV) FS=FV
  233 Y(I)=TA
      IF(FS.EQ.0.D0) GOTO 250
      FA=FY
      K=L-1
      FY=(EP(K)/FS)**(1D0/DBLE(L+K))
      IF(L.EQ.2) GOTO 240
      IF(FY.LT.0.7D0*FA) GOTO 250
  240 IF(FY.GT.0.7D0) GOTO 250
      JTI=JTI+1
      IF(JTI.GT.5) GOTO 250
      H=H*FY
      GOTO 10
  250 IF(KONV) GOTO 20
      D(3)=4.D0
      D(5)=1.6D1
      BO=.NOT.BO
      M=JR
      JR=JS
  260 JS=M+M
      IF (JTI.GT.5 .OR. JTJ.GT.0) GOTO 30
      JTJ=1
      H=H*0.1D0
      GOTO 10
   20 X=XN
      H=H*FY
      RETURN
   30 H=0.D0
      DO 300 I=1,N
  300 Y(I)=YA(I)
      RETURN
      END
