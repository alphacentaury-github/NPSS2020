$LARGE
      DOUBLE PRECISION FUNCTION S(A,B,C,L,X,Y)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 L
      PI  = 3.141592653589793D0
      X1  = 4.D0*PI*X*Y
      IF(MOD(L,2).EQ.1) X1=-X1
      S   = X1*EXPBES(A,B,C,L,X,Y)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION T(A,B,C,L,X,Y)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 L
      PI = 3.141592653589793D0
      X1 = 4.D0*PI
      IF(MOD(L,2).EQ.0) X1=-X1
      X2 = X1*X*Y
      X3 = DBLE(L)
      T  = ((X1*X3)/C)*EXPBES(A,B,C,L,X,Y)+X2*EXPBES(A,B,C,L+1,X,Y)
      RETURN
      END
C----------------------------------------------------------C
C  EXPBES(A,B,C,L,X,Y) = EXP(-(A*X*X+B*Y*Y)) J     (C*X*Y) C
C                                             L+1/2        C
C----------------------------------------------------------C
      DOUBLE PRECISION FUNCTION EXPBES(A,B,C,L,X,Y)
      IMPLICIT REAL*8 (A-Z)
      REAL*8    POLY(2)
      INTEGER*4 I,J,K,L,M,N
C
                      DIA = A+B-DABS(C)
      IF(DIA.LE.0.D0 ) THEN
         WRITE(*,*) 'ERROR IN SET OF PARAMETERS'
                      RETURN
      ENDIF
C
      Z      = C*X*Y
      EPSILO = Z/DABS(Z)
      U      = DABS(Z)
      XYZ    = A*X*X+B*Y*Y-U
C
      IF(XYZ.GT.100) THEN
                   EXPBES = 0.D0
                    RETURN
                     ELSE
        IF((DABS(Z).LE.DBLE(L)).OR.(DABS(Z).LE.1.D0)) THEN
               TERM1= 1.D0
                 DO 1000 K = 1,L
                     F     = DBLE(2*K+1)
                     TERM1 = TERM1*(Z/F)
 1000            CONTINUE
               IF(L.EQ.0) TERM1 = 1.D0
               VORFA = DEXP(-(A*X*X+B*Y*Y))
                     TERM2 = 0.D0
                         K = 1
   10                WERT  = 1.D0
                 DO 1001 M = 1,K
                     F     = DBLE(2*M*(2*(L+M)+1))
                     WERT  = WERT*((Z*Z)/F)
 1001            CONTINUE
               TERM2 = TERM2+VORFA*WERT
               IF(WERT.LT.1.D-10) GOTO 11
                        K = K+1
                        GOTO 10
   11          EXPBES = TERM1*(VORFA+TERM2)
                     RETURN
C
                            ELSE
          IF(DABS(Z).LE.20.D0) THEN
            POLY(1) = (DEXP(Z)-DEXP(-Z))/(2.D0*Z)
            POLY(2) =- POLY(1)/Z+(DEXP(Z)+DEXP(-Z))/(2.D0*Z)
C
               IF(L.LE.1) THEN
                      EXPBES = DEXP(-(A*X*X+B*Y*Y))*POLY(L+1)
                          ELSE
                        DO 100 K=1,L-1
                            POL     = POLY(1)-DBLE(2*K+1)*POLY(2)/Z
                            POLY(1) = POLY(2)
                            POLY(2) = POL
  100                   CONTINUE
                      EXPBES = DEXP(-(A*X*X+B*Y*Y))*POL
                            RETURN
               ENDIF
                               ELSE
               POLY(1) =             1.D0/(2.D0*U)
               POLY(2) = EPSILO*((U-1.D0)/(2.D0*U*U))
C
               IF(L.LE.1) THEN
                      EXPBES = DEXP(-XYZ)*POLY(L+1)
                          ELSE
                      DO 200 K=1,L-1
                          POL     = POLY(1)-DBLE(2*K+1)*
     1                              EPSILO*POLY(2)/U
                          POLY(1) = POLY(2)
                          POLY(2) = POL
  200                 CONTINUE
                   EXPBES = DEXP(-XYZ)*POL
                        RETURN
               ENDIF
          ENDIF
        ENDIF
      ENDIF
      END
C---------------------------------------------------------------------C
C      PROGRAM CALCULATES THE COULOMB-INTEGRALS                       C
C                                                                     C
C 1. FUN = R*R'*ERF(|LI*R+MI*R'|)*EXP(-A*R*R-B*R'*R')*EXP(-CI*R*R'*Z) C
C          *1/(|LI*R+MI*R'|)*LEGEND(L,Z)                              C
C    OR                                                               C
C                                                                     C
C 2. FUN = R*R'*1/(|LI*R+MI*R'|)*EXP(-A*R*R-B*R'*R')*EXP(-CI*R*R'*Z)  C
C          *LEGEND(L,Z)                                               C
C---------------------------------------------------------------------C
      DOUBLE PRECISION FUNCTION COULI(A,B,LI,MI,CI,L,X,Y,IND)
      IMPLICIT REAL*8 (A-Z)
      INTEGER *4      LX,INDX,L,N,IND,IFAIL
      COMMON /COULOM/ LIX,MIX,XC,U,V,TERM,LX,INDX
      EXTERNAL        D01BAZ,FUN
C
      IFAIL = 0
      N= 64
C
      LIX = LI
      MIX = MI
      CIX = CI
      U   = X
      V   = Y
      LX  = L
      INDX= IND
C---------------------------------------------C
C N   = NUMBER OF MESHPOINTS                  C
C FUN = FUNCTION TO BE INTEGRATED             C
C UG  = LOWER INTEGRAL BOUNDARY               C
C OG  = UPPER INTEGRAL BOUNDARY               C
C---------------------------------------------C
       XA   = (2.D0*A*MI-CI*LI)/(2.D0*MI)
       XB   = (2.D0*B*LI-CI*MI)/(2.D0*LI)
       XC   = CI/(2.D0*MI*LI)
       UG   = ABS(LI*X-MI*Y)
       OG   = ABS(LI*X+MI*Y)
C
       COULI = 0.D0
       TERM = XA*X*X+XB*Y*Y
       IF(TERM.GE.5.D1)  GOTO 200
C
       WERT  = D01BAF(D01BAZ,UG,OG,N,FUN,IFAIL)
C
       COULI = WERT/(LI*MI)
C
  200 RETURN
      END
C
      DOUBLE PRECISION FUNCTION FUN(Z)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4       LX,INDX,IFAIL
      COMMON /COULOM/ LIX,MIX,XC,U,V,TERM,LX,INDX
C
         ARGU  = (Z*Z-(LIX*LIX*U*U+MIX*MIX*V*V))/(2.D0*LIX*MIX*U*V)
C
      IF(INDX.EQ.2) GOTO 1000
C
      FUN = S15AEF(Z,IFAIL)*DEXP(-(XC*Z*Z+TERM))*LEGEND(LX,ARGU)
      GOTO 200
 1000 FUN = DEXP(-(XC*Z*Z+TERM))*LEGEND(LX,ARGU)
C
  200 RETURN
      END
C--------------------------------------------C
C     NORMALIZED  LEGENDRE-POLYNOMIALS       C
C--------------------------------------------C
      DOUBLE PRECISION FUNCTION LEGEND(L,X)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 N,L,K
      REAL*8    F(20),X
      N      = L+1
      F(1)   = 1.D0
      IF(N.EQ.1) GOTO 1000
      F(2)   = X
      IF(N.EQ.2) GOTO 1000
      DO 1001 K=3,N
      A      = DBLE(2*K-3)/DBLE(K-1)
      B      = DBLE(K-2)  /DBLE(K-1)
      F(K)   = A*X*F(K-1)-B*F(K-2)
 1001 CONTINUE
 1000 LEGEND = F(N)
      RETURN
      END
C-----------------C
C KRONECKER-DELTA C
C-----------------C
      DOUBLE PRECISION FUNCTION DKRDEL(M,N)
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 I,J,K,L,M,N
         DKRDEL = 0.D0
         IF(M.EQ.N) DKRDEL = 1.D0
      RETURN
      END
C-------------------------------------------------------------C
C CALCULATION OF THE COULOMB-FUNCTIONS AND THEIR DERIVATIVES  C
C-------------------------------------------------------------C
      SUBROUTINE DFCOUL(L,ETA,RO,F,FP,G,GP,IEXP,SIGMA)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(1),FP(1),G(1),GP(1),IEXP(1),SIGMA(1)
      DATA DEPI/6.283185307179586D0/
      ETAC=ETA*ETA
      L1=L+1
      CALL DFCZ0(ETA,RO,F0,FP0,G0,GP0,I,S)
      F(1)=F0
      FP(1)=FP0
      G(1)=G0
      GP(1)=GP0
      IEXP(1)=I
      SIGMA(1)=S
      IF(L)1,1,2
    1 RETURN
    2 LINF=0
      IND=0
      IF((ETA.GT.0).AND.(RO.LT.(ETA+ETA)))GO TO 21
      Z=ETA+DSQRT(ETAC+DBLE(L*(L+1)))
      IF(RO.GE.Z)GO TO 20
    7 ROINF=ETA+DSQRT(ETAC+DBLE(LINF*(LINF+1)))
      IF(RO-ROINF)3,4,4
    4 IF(LINF-L)5,6,6
    5 LINF=LINF+1
      GO TO 7
    3 IND=1
    6 LIN=LINF+1
   20 XM=1.D0
      IF(IND.EQ.0)LIN=L1
      DO 8 J=2,LIN
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      F(J)=(ZAG*F(J-1)-FP(J-1))/ZIG
      FP(J)=ZIG*F(J-1)-ZAG*F(J)
      G(J)=(ZAG*G(J-1)-GP(J-1))/ZIG
      GP(J)=ZIG*G(J-1)-ZAG*G(J)
      IEXP(J)=I
      SIG=SIGMA(J-1)+DATAN(ETA/(J-1))
      IPI=SIG/DEPI
      SIG=SIG-IPI*DEPI
      IF(SIG)60,50,70
   60 IF(SIG.LT.(-DEPI/2.D0))SIG=SIG+DEPI
      GO TO 50
   70 IF(SIG.GT.(DEPI/2.D0))SIG=SIG-DEPI
   50 SIGMA(J)=SIG
    8 XM=XM+1.D0
      IF(IND.EQ.0)RETURN
      GO TO 22
   21 LIN=1
   22 FTEST=F(LIN)
      FPTEST=FP(LIN)
      LMAX=LINF+25+IDINT(5.D0*DABS(ETA))
      IF(LMAX-L)9,10,10
    9 LMAX=L
   10 FI=1.D0
      FPI=1.D0
   13 XM=LMAX+1
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      FL=(ZAG*FI+FPI)/ZIG
      FPL=ZAG*FL-ZIG*FI
      IF(DABS(FL)-1.D15)26,27,27
   26 IF(DABS(FPL)-1.D15)28,27,27
   27 FL=FL*1.D-15
      FPL=FPL*1.D-15
   28 FI=FL
      FPI=FPL
      IF(LMAX-L)11,11,12
   12 LMAX=LMAX-1
      GO TO 13
   11 F(LMAX+1)=FL
      FP(LMAX+1)=FPL
      IF(LMAX-LINF)15,15,14
   14 GO TO 12
   15 FACT=FTEST/F(LIN)
      FACTP=FPTEST/FP(LIN)
      INDICE=I/60
      XM=LINF
      DO 16 J=LIN,L1
      F(J)=F(J)*FACT
      FP(J)=FP(J)*FACTP
   25 IF(J.EQ.1)GO TO 16
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      G(J)=(ZAG*G(J-1)-GP(J-1))/ZIG
      GP(J)=ZIG*G(J-1)-ZAG*G(J)
      IF(DABS(G(J))-1.D60)17,18,18
   17 IF(DABS(GP(J))-1.D60)19,18,18
   18 G(J)=G(J)/1.D60
      GP(J)=GP(J)/1.D60
      INDICE=INDICE+1
   19 IEXP(J)=INDICE*60
      A=FP(J)*G(J)
      B=-F(J)*GP(J)
      IF(A-1.D0)29,30,30
   29 I1=IDINT(DLOG10(A))
      I2=IDINT(DLOG10(B))
      GO TO 31
   30 I1=IDINT(DLOG10(A))+1
      I2=IDINT(DLOG10(B))+1
   31 F(J)=F(J)*1.D1**(-I2)
      FP(J)=FP(J)*1.D1**(-I1)
      SIG=SIGMA(J-1)+DATAN(ETA/(J-1))
      IPI=SIG/DEPI
      SIG=SIG-IPI*DEPI
      IF(SIG)61,51,71
   61 IF(SIG.LT.(-DEPI/2.D0))SIG=SIG+DEPI
      GO TO 51
   71 IF(SIG.GT.(DEPI/2.D0))SIG=SIG-DEPI
   51 SIGMA(J)=SIG
   16 XM=XM+1.D0
      RETURN
      END
      SUBROUTINE DFCZ0(ETA,RO,F0,FP0,G0,GP0,IEXP,SIGMA0)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A1(110),A2(110),B1(110),B2(110)
      IF(RO)2,2,1
    2 WRITE (6,1000)
 1000 FORMAT (1X,'RO IS NEGATIVE OR ZERO')
      RETURN
    1 IF(ETA-30.D0)3,3,4
    4 IF(DABS(ETA)-5.D2)28,28,29
   28 CALL YFRICA(ETA,RO,F0,FP0,G0,GP0,IEXP,SIGMA0)
   29 WRITE (6,1001)
      RETURN
 1001 FORMAT (1X,'ABSOLUTE VALUE OF ETA IS LARGER THAN 500')
      RETURN
    3 IF(ETA+8.D0)4,5,5
    5 IF(ETA)6,7,6
    7 F0=DSIN(RO)
      G0=DCOS(RO)
      FP0=G0
      GP0=-F0
      IEXP=0
      SIGMA0=0.D0
      RETURN
    6 BORNE=1.666666666666667D0*DABS(ETA)+7.5D0
      IF(RO-BORNE)8,9,9
    9 CALL YFASYM(ETA,RO,F0,FP0,G0,GP0,IEXP,SIGMA0)
      RETURN
    8 IF(ETA-10.D0)10,11,11
   10 IF(ETA)12,12,13
   13 IF(RO-2.D0)14,12,12
   11 IF(ETA-(5.D0*RO+6.D1)/7.D0)12,12,14
   12 CALL YFASYM(ETA,BORNE,F0,FP0,G0,GP0,IEXP,SIGMA0)
      H=BORNE
      DH=F0/H
      IF(ETA)20,21,21
   20 N=-0.5D0*ETA+5.D0
      GO TO 22
   21 N=ETA/5.D0+5.D0
   22 N=5*(N+1)
      Z=4.D0/H
      Y=1.D0-(ETA+ETA)*Z
      A1(N+2)=1.D-55
      A1(N+3)=0.D0
      A1(N+4)=1.D-64
      B1(N+3)=1.D-50
      B1(N+4)=1.D-68
      A2(N+2)=0.D0
      A2(N+3)=1.D-74
      A2(N+4)=1.D-53
      B2(N+3)=0.D0
      B2(N+4)=1.D-66
      M=N+2
      DI=N
      DO 23 II=2,M
      I=M-II+2
      B1(I)=B1(I+2)+ Z*(DI+1.D0)*A1(I+1)
      S=A1(I+2)+Y*(A1(I+1)-A1(I))
      Q=(DI+2.D0)*B1(I)+(DI-1.D0)*B1(I+1)
      A1(I-1)=S-Z*Q
      B2(I)=B2(I+2)+ Z*(DI+1.D0)*A2(I+1)
      S=A2(I+2)+Y*(A2(I+1)-A2(I))
      Q=(DI+2.D0)*B2(I)+(DI-1.D0)*B2(I+1)
      A2(I-1)=S-Z*Q
      IF(I.GE.N)GO TO 23
      D=-(B2(I+2)+B2(I))/(B1(I+2)+B1(I))
      DO 24 J=I,M
      A2(J)=A2(J)+D*A1(J)
      B2(J)=B2(J)+D*B1(J)
   24 CONTINUE
      A2(I-1)=A2(I-1)+D*A1(I-1)
   23 DI=DI-1.D0
      Q=A1(3)-A1(1)
      C=A2(3)-A2(1)
      C=Q/C
      X1=0.5D0*(A1(2)-C*A2(2))
      DO 25 I=3,M
      X1=X1+A1(I)-C*A2(I)
   25 CONTINUE
      X1=DH/X1
      X2=-C*X1
      DO 26 I=2,M
      B1(I)=X1*B1(I)+X2*B2(I)
      A1(I)=X1*A1(I)+X2*A2(I)
   26 CONTINUE
      A1(1)=X1*A1(1)+X2*A2(1)
      B1(1)=0.D0
      X=RO/H
      Y=2.D0*X-1.D0
      T1=1.D0
      T2=Y
      RESULT=0.5D0*A1(2)+Y*A1(3)
      DERIVE=0.5D0*B1(2)+Y*B1(3)
      DO 27 I=2,N
      TI=2.D0*Y*T2-T1
      T1=T2
      T2=TI
      RESULT=RESULT+TI*A1(I+2)
      DERIVE=DERIVE+TI*B1(I+2)
   27 CONTINUE
      F0=RESULT*RO
      FP0=DERIVE*RO+RESULT
      GO TO 30
C   EXPANSION AT THE ORIGIN
   14 PI=3.141592653589793D0
      CALL JFLGAM(1.D0,ETA,TRA,SIGMA0,NTRUC)
      IEXP=0
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      PIETA2=0.5D0*PIETA
      B=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      U0=0.D0
      U1=RO
      U=U0+U1
      UP=1.D0
      XN=2.D0
      DO 15 N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      UP=UP+XN*U2/RO
   17 IF(DABS(U2/U).LT.1.D-10)GOTO19
   18 U0=U1
      U1=U2
      XN=XN+1.D0
   15 CONTINUE
   19 F0=U/B
      FP0=UP/B
   30 CALL YFIREG(ETA,RO,G0,GP0)
      RETURN
      END
      SUBROUTINE JFDELG (XD,YD,PAR,PAI,NBCHIF)
      DOUBLE PRECISION XD,YD,PAR,PAI,TEST,C,PI
      DOUBLE PRECISION X,Y,U,V,TRA,TRA1,COSI,SINI
      DOUBLE PRECISION COS2I,SIN2I,ZMOD,DEPI
      DOUBLE PRECISION TRB,XX
      DOUBLE PRECISION RAC2,PIS4
      DOUBLE PRECISION SUPINT
      DIMENSION TEST(7),C(6)
      DATA TEST/2.9152D+7,2.2958D+3,1.4124D+2,3.9522D+1,19.6611D0,12.791
     1D0,-10.D0/
      DATA RAC2/0.3465735902799726/,PIS4/0.785398163397448/
      DATA C/8.333333333333333D-2,-8.33333333333333D-3,
     1 3.968253968253968D-3,-4.166666666666667D-3,
     2 7.575757575757576D-3, -2.109279609279609D-2/
      DATA PI/3.141592653589793/
      DATA DEPI/6.283185307179586/
      DATA SUPINT/2147483647.D0/
      X=DABS(XD)
      XX=X
      NBCHIF=15
      IF(YD)1,2,1
    1 Y=DABS(YD)
      KR=1
      I=DMOD(10.99D0-X,SUPINT)
C     TRANSLATION
      IF(I)3,3,4
    4 TRA=I
      X=X+TRA
C     LOGARITHM (X+IY) (X,Y HAVE TO BE POSITIVE)
    3 IF(X-Y)5,6,7
    5 TRA1=X/Y
      TRB=1.D0+TRA1*TRA1
      TRA=Y*DSQRT(TRB)
      SINI=1./(TRB*Y)
      COSI=SINI*TRA1
      TRA1=Y/X
      GO TO 11
    6 U=RAC2+DLOG(X)
      V=PIS4
      SINI=0.5D0/X
      COSI=SINI
      GO TO 10
    7 TRA1=Y/X
      TRB=1.D0+TRA1*TRA1
      TRA=X*DSQRT(TRB)
      COSI=1./(TRB*X)
      SINI=COSI*TRA1
   11 U=DLOG(TRA)
      V=DATAN(TRA1)
C ***********************************************C
C     ASYMPTOTIC EXPANSION FOR X LARGER THAN 10  C
C ***********************************************C
   10 PAR=U-0.5*COSI
      PAI=V+0.5*SINI
      ZMOD=X+Y
      IF(ZMOD-TEST(1))13,13,14
   13 SIN2I=(SINI*COSI)+(SINI*COSI)
      COS2I=(COSI+SINI)*(COSI-SINI)
      SINI=SIN2I
      COSI=COS2I
      K=1
      GO TO 15
   16 TRA=COSI*COS2I-SINI*SIN2I
      SINI=SINI*COS2I+COSI*SIN2I
      COSI=TRA
   15 PAR=PAR-C(K)*COSI
      PAI=PAI+C(K)*SINI
      K=K+1
      IF(ZMOD-TEST(K))16,16,14
C ************************C
C     INVERSE TRANSLATION C
C ************************C
   17 I=I-1
      X=I
      X=XX+X
   56 IF(X-Y)55,55,57
   55 TRA1=X/Y
      TRB=X*TRA1+Y
      SINI=1.D0/TRB
      COSI=TRA1/TRB
      GO TO 19
   57 TRA1=Y/X
      TRB=X+Y*TRA1
      COSI=1.D0/TRB
      SINI=TRA1/TRB
   19 PAR=PAR-COSI
      PAI=PAI+SINI
   14 IF(I)18,18,17
C *****************************C
C     CONTROL OF THE QUADRANT  C
C *****************************C
   18 IF(XD)20,61,21
   61 TRA=PI*Y
      IF(TRA-1.D-2)300,300,301
  300 TRB= TRA*(2.D0+TRA*(-2.D0+TRA*(1.333333333333333D0+
     1 TRA*(  -0.6666666666666666D0+TRA*(   0.2666666666666666D0+
     2 TRA*( -0.08888888888888888D0+TRA*   0.02539682539682540D0 ))))))
      TRB=(2.D0-TRB)/TRB
      GO TO 302
  301 TRB= DEXP(-TRA-TRA)
      TRB=(1.D0+TRB)/(1.D0-TRB)
  302 PAI=0.5D0*(1.D0/Y+PI*TRB)
   21 IF(YD)28,100,100
C ****************************C
C     CHANGE X+IY INTO -X-IY  C
C ****************************C
   20 TRA=DEXP(-DEPI*Y)
      TRB=TRA*TRA
      COS2I=DEPI*DMOD(X,1.D0)
      SIN2I=-2.D0*TRA*DCOS(COS2I)+1.D0+TRB
      PAR=PAR+COSI+DEPI*TRA*DSIN(COS2I)/SIN2I
      PAI=PAI-SINI+PI*(TRB-1.D0)/SIN2I
      IF(YD)100,100,28
   28 PAI=-PAI
C     ARGUMENT DANS -PI,PI
  100 TRA=DABS(PAI/DEPI)
      IF(TRA-1.D+15)203,204,204
  204 NBCHIF=0
      PAI=0.D0
      GO TO 29
  203 IF(TRA-1.D0)205,205,206
  206 NBCHIF=DLOG10(TRA)
      NBCHIF=14-NBCHIF
      TRA=DMOD(TRA,SUPINT)
      PAI=DMOD(TRA,1.D0)*DSIGN(DEPI,PAI)
  205 IF(DABS(PAI)-PI)29,29,207
  207 PAI=PAI-DSIGN(DEPI,PAI)
      GO TO 29
C
    2 PAI=0.D0
      IF(XD)31,32,33
C     CONDITION OF EXISTENCE
   32 WRITE (6,1000)
 1000 FORMAT (1X,'JFDELG(0) IS INFINITE ')
      GO TO 50
   31 IF(X-4503599627370496.D0)34,35,35
   35 WRITE (6,1001)
 1001 FORMAT (1X,'ARGUMENT OF JFDELG TOO LARGE ')
      GO TO 50
   34 Y=DMOD(X,SUPINT)
      IF(Y) 400,36,400
  400 IF(Y-0.99D0) 33,33,405
  405 TRA= IDINT(Y+0.1D0)
      IF(DABS(Y-TRA)-5.D-15)36,36,33
   36 WRITE (6,1002)
 1002 FORMAT (1X,'JFDELG IS INFINITE FOR NEGATIVE ARGUMENTS')
   50 PAR=1.D+74
      NBCHIF=0
      GO TO 29
C *****************C
C     TRANSLATION  C
C *****************C
   33 I=DMOD(10.99D0-X,SUPINT)
      IF(I)37,37,38
   38 TRA=I
      X=X+TRA
C ***************************C
C     ASYMPTOTIC EXPANSION   C
C ***************************C
   37 Y=DLOG(X)
      PAR=Y-0.5D0/X
      IF(X-TEST(1))39,39,43
   39 COS2I=1.D0/(X*X)
      COSI=COS2I
      K=1
      GO TO 41
   42 COSI=COSI*COS2I
   41 PAR=PAR-C(K)*COSI
      K=K+1
      IF(X-TEST(K))42,42,40
C ************************C
C     INVERSE TRANSLATION C
C ************************C
   44 I=I-1
      X=I
      X=XX+X
      PAR=PAR-1.D0/X
   40 IF(I)43,43,44
C ********************C
C     X IS NEGATIVE   C
C ********************C
   43 IF(XD)45,29,29
   45 PAR=PAR+1.D0/X
      Y=PI*DMOD(X,2.D0)
      PAR=PAR+PI*DCOS(Y)/DSIN(Y)
   29 RETURN
      END
      SUBROUTINE YFASYM(ETA,RAU,FO,FPO,GO,GPO,IEXP,SIGO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IEXP=0
      TRB=0.D0
      RAU2=RAU+RAU
      ETAC=ETA*ETA
      CALL JFLGAM(1.D0,ETA,TRA,SIGO,NTRUC)
   40 N=0
      PS=1.D0
      GS=0.D0
      PT=0.D0
      GT=1.D0-ETA/RAU
      SF=PS
      SG=GS
      SPF=PT
      SPG=GT
   45 DENOM=  DBLE(N+1)*RAU2
      AN= DBLE(N+N+1)*ETA/DENOM
      BN= (ETAC- DBLE( N*(N+1)))/DENOM
      PS1=AN*PS-BN*PT
      GS1=AN*GS-BN*GT-PS1/RAU
      PT1=AN*PT+BN*PS
      GT1=AN*GT+BN*GS-PT1/RAU
   42 SF=SF+PS1
      SG=SG+GS1
      SPF=SPF+PT1
      SPG=SPG+GT1
      N=N+1
      IF(N-17)46,48,44
   48 TRA=PS*PS+PT*PT
   44 TRB=PS1*PS1+PT1*PT1
      TEST=TRA-TRB
      IF(TEST)47,47,46
   46 PS=PS1
      GS=GS1
      PT=PT1
      GT=GT1
      TRA=TRB
      GOTO  45
   47 TETAO= RAU-ETA*DLOG (RAU2)+SIGO
      TRA= DSIN(TETAO)
      TRB=DCOS(TETAO)
      GO=SF*TRB-SPF*TRA
      GPO=SG*TRB-SPG*TRA
      FO=SPF*TRB+SF*TRA
      FPO=SPG*TRB+SG*TRA
      RETURN
      END
      SUBROUTINE YFIREG(ETA,RO,G0,GP0)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(ETA.LE.0.D0)GOTO250
      IF(ETA.LE.3.D0)GOTO251
      IF(ETA.LE.1.D1)GOTO252
      IF(ETA.LE.18.D0)GOTO253
      IF(ETA.LE.22.D0)GOTO254
      IF(RO.LE.0.3D0+(3.D1-ETA)/8.D1)GOTO200
C *********************************C
C   TAYLOR-EXPANSION AROUND RAU0   C
C *********************************C
  300 CONTINUE
      RAU0=1.666666666666667D0*DABS(ETA)+7.5D0
      CALL YFASYM(ETA,RAU0,F0,FP0,G0,GP0,IEXP,SIGMA0)
      X=RAU0-RO
      X2=X*X
      X3=X*X2
      UNR=1.D0/RAU0
      ETR0=1.D0-2.D0*ETA*UNR
      U0=G0
      U1=-X*GP0
      U2=-0.5D0*ETR0*X2*U0
      S=U0+U1+U2
      V1=U1/X
      V2=2.D0*U2/X
      T=V1+V2
      XN=3.D0
      DO10N=3,10000
C     N=N
      XN1=XN-1.D0
      XN1=XN*XN1
      U3=X*U2*UNR*(1.D0-2.D0/XN)-ETR0*U1*X2/XN1+X3*U0*UNR/XN1
      S=S+U3
      V3=XN*U3/X
      T=T+V3
   16 IF(DABS(U3/S).GT.1.D-11)GO TO 11
      IF(DABS(V3/T).LE.1.D-11)GO TO 20
   11 U0=U1
      U1=U2
      U2=U3
      XN=XN+1.D0
   10 CONTINUE
   20 G0=S
      GP0=-T
      RETURN
C ***************************C
C   EXPANSION AT THE ORIGIN  C
C ***************************C
  200 CONTINUE
      PI=3.141592653589793D0
      GA=0.577215664901533D0
      ETA2=ETA*ETA
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      PIETA2=0.5D0*PIETA
      B=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      CALL JFDELG(1.D0,ETA,PAR,PAI,NB)
      C1=ETAP*(GA+GA+DLOG(2.D0)-1.D0+PAR)
      U0=0.D0
      U1=RO
      V0=1.D0
      V1=C1*RO
      U=U0+U1
      V=V0+V1
      UP=1.D0
      VP=C1
      XN=2.D0
      DO 104N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      V2=(ETAP*RO*V1-RO2*V0-ETAP*(XN+XN-1.D0)*U2)/XN1
      V=V+V2
      UP=UP+XN*U2/RO
      VP=VP+XN*V2/RO
   17 IF(DABS(U2/U).GT.1.D-14)GOTO18
      IF(DABS(V2/V).LE.1.D-14)GOTO19
   18 U0=U1
      U1=U2
      V0=V1
      V1=V2
      XN=XN+1.D0
  104 CONTINUE
   19 GP=V+ETAP*U*DLOG(RO)
      G0=B*GP
      GP0=B*(VP+ETAP*(UP*DLOG(RO)+U/RO))
      RETURN
  250 IF(RO.LE.0.5D0*ETA+9.D0)GOTO200
      GOTO 300
  251 IF(RO.LE.2.25D0+7.35D0*(3.D0-ETA))GOTO200
      GOTO 300
  252 IF(RO.LE.1.2D0+1.5D-1*(1.D1-ETA))GOTO200
      GOTO 300
  253 IF(RO.LE.0.6D0+0.75D-1*(18.D0-ETA))GOTO200
      GOTO 300
  254 IF(RO.LE.0.4D0+0.5D-1*(22.D0-ETA))GOTO200
      GOTO 300
      END
      SUBROUTINE YFCLEN(ETA,RO,U,UP,V,VP,SIGMA0,IDIV,NN)
      IMPLICIT COMPLEX*16(A-D,F-H), REAL*8(E,P-Z)
C
      IF(NN.EQ.1)GO TO 20
C
      ETA2=ETA*ETA
      FA=DCMPLX(1.D0,ETA)
      M=0.25*ETA+4.D1
C ****************************************C
C     CHEBYSHEV POLYNOMIALS TO ORDER M    C
C ****************************************C
      K=M+2
      X=(ETA+ETA)/RO
      XX=X+X-1.D0
      T0=1.D0
      T1=XX
      XX=XX+XX
      DO 6 J=2,K
      TJ=XX*T1-T0
      T0=T1
      T1=TJ
    6 CONTINUE
      TM=T1
      TL=T0
C  **************************C
C     INITIALIZATION         C
C  **************************C
      AM=(0.D0,0.D0)
      AL=(1.D-40,1.D-40)
      BN=(0.D0,0.D0)
      BM=(1.D-40,1.D-40)
      BL=(0.D0,0.D0)
      BK=DCMPLX(4.D0*DBLE(M+1),0.D0)*AL+BM
      F=(0.D0,0.D0)
      G =(0.D0,0.D0)
      GP=(0.D0,0.D0)
C *************************C
C     BACKWARD RECURSION   C
C *************************C
      K=M
   10 R=K
      TK=XX*TL-TM
      TM=TL
      TL=TK
      HK=DCMPLX(TK,0.D0)
      C1=DCMPLX(R*(R+1.D0)-ETA2,ETA*(R+R+1.D0))
      C2=(4.D0,0.D0)*DCMPLX(R+1.D0,0.D0)
      C2=C2*DCMPLX(-R-1.D0,ETA*3.D0)
      C3=FA*DCMPLX(-R-R-4.D0,ETA)
      C4=DCMPLX((7.D0*R+5.D0)/4.D0,0.D0)
      C5=DCMPLX(R+R+2.D0,0.D0)
      C6=DCMPLX((R+3.D0)/4.D0,0.D0)
      AK=(C2*AL+C3*AM-C4*BL-C5*BM-C6*BN)/C1
      J=K/2
      J=K-J-J
      IF(J)1,2,1
    1 F=F-AK
      GO TO 3
    2 F=F+AK
    3 Z=CDABS(AK)
      G=G+HK*AK
      GP=GP+HK*BK
C
C          F=A0/2-A1+A2-A3+A4-A5+...
C
C          CONGRUENCE MODULO 10**60
C
      IF(Z-1.D60)4,5,5
    5 D=(1.D60,0.D0)
      AK=AK/D
      AL=AL/D
      AM=AM/D
      BK=BK/D
      BL=BL/D
      BM=BM/D
      BN=BN/D
      F=F/D
      G=G/D
      GP=GP/D
    4 IF(K)8,8,9
    9 D=DCMPLX(4.D0*R,0.D0)
      BJ=D*AK+BL
      AM=AL
      AL=AK
      BN=BM
      BM=BL
      BL=BK
      BK=BJ
      K=K-1
      GO TO 10
C *********************************************C
C     NORMALIZATION AND CALCULATION OF Z(RO)   C
C *********************************************C
    8 D=(0.5D0,0.D0)*AK
      F=F-D
      G=G-D
      GP=GP-(0.5D0,0.D0)*BK
      D=DCMPLX(0.D0,-ETA*DLOG(2.D0)+SIGMA0)
      AXPO=CDEXP(D)
      F=F/AXPO
      G=G/F
      GP=GP/F
C ******************************************************************C
C     CALCULATION OF REGULAR AND IRREGULAR COULOMB FUNCTIONS F, G   C
C ******************************************************************C
      D=DCMPLX(0.D0,RO-ETA*DLOG(RO))
      AXPO=CDEXP(D)
      D=G*AXPO
      GP=AXPO*(DCMPLX(0.D0,1.D0-ETA/RO)*G-DCMPLX(X/RO,0.D0)*GP)
      V=D
      D=(0.D0,-1.D0)*D
      U=D
      VP=GP
      GP=(0.D0,-1.D0)*GP
      UP=GP
      IDIV=0
      RETURN
C  ******************************C
C     EXPANSION AT THE ORIGIN    C
C  ******************************C
   20 PI=3.141592653589793D0
      XA=0.577215664901533D0
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      Z=138.15510557964276D0
      IDIV=0
      IF(DABS(PIETA)-Z)21,21,22
   22 INDG=IDINT(PIETA/Z)
      IDIV=60*INDG
      IF(ETA.LT.0) IDIV=0
      PIETA=PIETA-Z*DBLE(INDG)
   21 PIETA2=0.5D0*PIETA
      P=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      CALL JFDELG(1.D0,ETA,PAR,PAI,NB)
      Z1=ETAP*(XA+XA+DLOG(2.D0)-1.D0+PAR)
      U0=0.D0
      U1=RO
      V0=1.D0
      V1=Z1*RO
      U=U0+U1
      V=V0+V1
      UP=1.D0
      VP=Z1
      XN=2.D0
      DO 104N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      V2=(ETAP*RO*V1-RO2*V0-ETAP*(XN+XN-1.D0)*U2)/XN1
      V=V+V2
      UP=UP+XN*U2/RO
      VP=VP+XN*V2/RO
      IF(DABS(U2/U).GT.1.D-14)GOTO18
      IF(DABS(V2/V).LE.1.D-14)GOTO19
   18 U0=U1
      U1=U2
      V0=V1
      V1=V2
      XN=XN+1.D0
  104 CONTINUE
   19 PP=V+ETAP*U*DLOG(RO)
      W=U/P
      WP=UP/P
      V=P*PP
      VP=P*(VP+ETAP*(UP*DLOG(RO)+U/RO))
      U=W
      UP=WP
      RETURN
      END
      SUBROUTINE YFRICA(ETA,RAU,FO,FPO,GO,GPO,IDIV,SIGMA0)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(5),QP(5)
C  *********************************C
C         RICCATI COEFFICIENTS      C
C  *********************************C
      DATA G61,G62,G63,G64,G65,G66,G67,G68,G69,
     1     G610,G611/
     2  0.1159057617187498D-01,  0.3863525390624998D-01,
     3  0.4660034179687498D-01,  0.4858398437499998D-01,
     4  0.1156514485677080D 01,  0.5687475585937496D 01,
     5  0.1323888288225445D 02,  0.1713083224826384D 02,
     6  0.1269003295898436D 02,  0.5055236816406248D 01,
     7  0.8425394694010415D 00/
      DATA G81,G82,G83,G84,G85,G86,G87,G88,G89,
     1     G810,G811,G812,G813,G814,G815/
     2  0.1851092066083633D-01,  0.8638429641723630D-01,
     3  0.1564757823944092D 00,  0.1430139541625977D 00,
     4  0.1924622058868408D 00,  0.8500803152720129D 01,
     5  0.7265429720878595D 02,  0.3057942376817972D 03,
     6  0.7699689544836672D 03,  0.1254157054424285D 04,
     7  0.1361719536066055D 04,  0.9831831171035763D 03,
     8  0.4547869927883148D 03,  0.1222640538215636D 03,
     9  0.1455524450256709D 02/
      DATA GP61,GP62,GP63,GP64,GP65,GP66/
     2  0.2897644042968748D-01,  0.2318115234375000D 00,
     3  0.8056640625000000D 00,  0.1601562499999998D 01,
     4  0.3046875000000000D 00,  0.5624999999999998D 01/
      DATA GP81,GP82,GP83,GP84,GP85,GP86,GP87,GP88/
     1  0.6478822231292720D-01,  0.6910743713378906D 00,
     2  0.3322952270507811D 01,  0.9483032226562498D 01,
     3  0.1769653320312499D 02,  0.3478710937499998D 02,
     4  0.5020312499999999D 02,  0.7874999999999999D 02/
      DATA Q /0.4959570165D-1,0.8888888889D-2,0.2455199181D-2,
     $ 0.9108958061D-3,0.2534684115D-3/
      DATA QP /0.1728260369D0,0.3174603174D-3,0.3581214850D-2,
     1 0.3117824680D-3,0.9073966427D-3/
      CALL JFLGAM(1.D0,ETA,TRA,SIGMA0,IND)
      RAU2=RAU+RAU
      RAUC=RAU*RAU
      ETAC=ETA*ETA
      ETA2=ETA+ETA
      ETARO=ETA*RAU
      ETARO2=ETARO+ETARO
      PIETA=3.141592653589793*ETA
      IND=0
      JND=0
      IG=0
      IF(ETA)20,20,21
   20 IF(-ETARO-14.0625D0)32,22,22
   22 INDICE=1
C             RICCATI 3
      IDIV=0
      GO TO 2
   21 IF(DABS(RAU-ETA2).LE.1.D-9)GO TO 18
      IF(RAU-ETA2)30,18,31
   31 IF(RAU-ETA2-2.D1*(ETA**0.25D0))34,33,33
   33 INDICE=0
C             RICCATI  2
      IDIV=0
      GO TO 2
   32 NN=1
      GO TO 35
   34 NN=0
   35 CALL YFCLEN(ETA,RAU,FO,FPO,GO,GPO,SIGMA0,IDIV,NN)
      RETURN
   30 IF(ETARO-12.D0)32,32,23
   23 TRA=ETA2-6.75D0*(ETA**0.4D0)
      IF(RAU-TRA)6,6,24
   24 IND=1
      JND=1
      RO=RAU
      RAU=TRA
      RAU0=TRA
C             RICCATI  1
    6 X=RAU/ETA2
      U= (1.D0-X)/X
      X2=X*X
      RU= DSQRT(U)
      RX= DSQRT(X)
      TRE= 1.D0/(U*RU*ETA2)
      TRB=TRE*TRE
      FI= (DSQRT((1.D0-X)*X)+ DASIN(RX)-1.570796326794897D0)*ETA2
      TR1= -0.25D0*DLOG(U)
  602 TR2= -((9.D0*U+6.D0)*U+5.D0)/48.D0
  603 TR3= ((((-3.D0*U-4.D0)*U+6.D0)*U+12.D0)*U+5.D0)/64.D0
  604 TR4=- ((((((U+2.D0)*945.D0*U+1395.D0)*U+12300.D0)*U+25191.D0)*U
     1   +19890.D0)*U+5525.D0)/46080.D0
  605 TR5= (((((((( -27.D0*U-72.D0)*U-68.D0)*U+360.D0)*U+2190.D0)*U
     1  +4808.D0)*U+5148.D0)*U+2712.D0)*U+565.D0)/2048.D0
  606 TR6=-  (((((((((G61*U+G62)*U+G63)*U+G64)*U+G65)*U+G66)*U+G67)*U
     1                +G68)*U+G69)*U+G610)*U+G611
  607 TR7= ((((((((((((-81.*U-324.)*U-486.)*U-404.)*U+4509.)*U+52344.)
     1                 *U+233436.)*U+567864.)*U+838521.)*U+775884.)*U
     2                 +441450.)*U+141660.)*U+19675.) /6144.
  608 TR8=   ((((((((((((( G81*U+G82)*U+G83)*U+G84)*U+G85)*U+G86)*U
     1                    +G87)*U+G88)*U+G89)*U+G810)*U+G811)*U+G812)*U
     2                    +G813)*U+G814)*U+G815
      FI= FI+TRE*(TR2+TRB*(TR4+TRB*(TR6+TRB*TR8)))
      PSI=-FI
      TRA= TR1+TRB*(TR3+TRB*(TR5+TRB*TR7))
      FI=FI+TRA
      PSI=PSI+TRA
      FIP=RU*ETA2
      TRA=1.D0/(X2*U)
      TR1=0.25D0
      TRE= TRE/(X2*X2*U)
      TRB=TRB/(X2*X2)
  702 TR2= -(8.D0*X-3.D0)/32.D0
  703 TR3= (( 24.D0*X-12.D0)*X+3.D0)/64.D0
  704 TR4= ((( -1536.D0*X+704.D0)*X-336.D0)*X+63.D0)/2048.D0
  705 TR5= (((( 1920.D0*X-576.D0)*X+504.D0)*X-180.D0)*X+27.D0)/1024.D0
  706 TR6 =  (((( -GP66*X+GP65)*X-GP64)*X+GP63)*X-GP62)*X+GP61
  707 TR7= - (((((( -40320.D0*X-10560.D0)*X-13248.D0)*X+7560.D0)*X
     1   -3132.D0)*X+756.D0)*X-81.D0) / 2048.D0
  708 TR8 =- (((((( GP88*X+GP87)*X+GP86)*X-GP85)*X+GP84)*X-GP83)*X+GP82
     1             )*X-GP81
      FIP=FIP +TRE*(TR2+TRB*(TR4+TRB*(TR6+TRB*TR8)))
      PSIP=-FIP
      TRA= TRA*(TR1+TRB*(TR3+TRB*(TR5+TRB*TR7)))
      FIP=FIP+TRA
      PSIP=PSIP+TRA
      XXX=138.1551055796428D0
      INDG=IDINT(PSI/XXX)
      IDIV=60*INDG
      IF(INDG.EQ.0)GO TO 8
      PSI=PSI-XXX*DBLE(INDG)
      FI =FI +XXX*DBLE(INDG)
    8 FO=0.5D0*DEXP(FI)
      GO= DEXP( PSI)
      FPO= FO*FIP/ETA2
      GPO=GO*PSIP/ETA2
      IF(JND.EQ.0)RETURN
      RAU=RO
      GO=FO
      GPO=FPO
   27 X=RAU0-RO
      X2=X*X
      X3=X*X2
      UNR=1.D0/RAU0
      ETR0=1.D0-2.D0*ETA*UNR
      U0=GO
      U1=-X*GPO
      U2=-0.5D0*ETR0*X2*U0
      S=U0+U1+U2
      V1=U1/X
      V2=2.D0*U2/X
      T=V1+V2
      XN=3.D0
      DO10N=3,10000
C     N=N
      XN1=XN-1.D0
      XN1=XN*XN1
      U3=X*U2*UNR*(1.D0-2.D0/XN)-ETR0*U1*X2/XN1+X3*U0*UNR/XN1
      S=S+U3
      V3=XN*U3/X
      T=T+V3
   16 IF(DABS(U3/S).GT.1.D-10)GO TO 11
      IF(DABS(V3/T).LE.1.D-10)GO TO 17
   11 U0=U1
      U1=U2
      U2=U3
      XN=XN+1.D0
   10 CONTINUE
   17 IF(IG)25,26,25
   25 GO=S
      GPO=-T
      FO=HO
      FPO=HPO
      RETURN
   26 HO=S
      HPO=-T
   18 ET0=ETA**(0.166666666666667)
      ETAD=ETAC*ETAC
      ET=ETA**(0.6666666666666667)
      ET1=ET*ET
      ET2=ET1*ET1
      ET3=ET2*ET
      ET4=ETAD*ET
      ET5=ET4*ET
      FO=1.D0-Q(1)/ET1-Q(2)/ETAC-Q(3)/ET3-Q(4)/ETAD-Q(5)/ET5
      GO=1.D0+Q(1)/ET1-Q(2)/ETAC+Q(3)/ET3-Q(4)/ETAD+Q(5)/ET5
      FPO=1.D0+QP(1)/ET+QP(2)/ETAC+QP(3)/ET2+QP(4)/ETAD+QP(5)/ET4
      GPO=1.D0-QP(1)/ET+QP(2)/ETAC-QP(3)/ET2+QP(4)/ETAD-QP(5)/ET4
      FO=0.7063326373D0*ET0*FO
      GO=1.223404016D0*ET0*GO
      FPO=0.4086957323D0*FPO/ET0
      GPO=-0.7078817734D0*GPO/ET0
      IDIV=0
      IF(IND.EQ.0)RETURN
      IG=1
      RAU0=ETA2
      GO TO 27
    2 X=ETA2/RAU
      X2=X*X
      U=1.D0-X
      RU= DSQRT(U)
      U3=U*U*U
      TRD= 1.D0/(U3*ETA2*ETA2)
      TRC=X2*TRD
      TRE=1.D0/(U*RU*ETA2)
      FI= -0.25D0*DLOG(U)
      TRB=TRD/64.D0
      TR3=  (((3.D0*U-4.D0)*U-6.D0)*U+12.D0)*U-5.D0
  501 TR5= (((((((( -27.D0*U+72.D0)*U-68.D0)*U-360.D0)*U+2190.D0)*U
     1      -4808.D0)*U+5148.D0)*U-2712.D0)*U+565.D0)/32.D0
  502 TR7= (((((((((((( 81.D0*U-324.D0)*U+486.D0)*U-404.D0)*U-4509.D0)*U
     1    +52344.D0)*U-233436.D0)*U+567864.D0)*U-838521.D0)*U+775884.D0
     2       )*U-441450.D0)*U+141660.D0)*U-19675.D0)/96.D0
      FI= FI+TRB*(TR3+TRD*(TR5+TRD*TR7))
      FIP=0.25D0/U
      TRB=3.D0*TRC/(64.D0*U)
      TR3= (X-4.D0)*X+8.D0
  511 TR5= ((((9.D0*X-60.D0)*X+168.D0)*X-192.D0)*X+640.D0)/16.D0
  512 TR7= ((((((-27.D0*X+252.D0)*X-1044.D0)*X+2520.D0)*X-4416.D0)*X
     1   -3520.D0)*X-13440.D0)/32.D0
      FIP =FIP+TRB*(TR3+TRC*(TR5+TRC*TR7))
      TRA= DABS((RU-1.D0)/(RU+1.D0))
      PSI= (0.5D0*DLOG(TRA)+RU/X)*ETA2+0.785398163397448D0
      TR2= -((9.D0*U-6.D0)*U+5.D0)/48.D0
  521 TR4= ((((((U-2.D0)*945.D0*U+1395.D0)*U-12300.D0)*U+25191.D0)*U
     1   -19890.D0)*U+5525.D0)/46080.D0
  522 TR6 = (((((((((-G61*U+G62)*U-G63)*U+G64)*U-G65)*U+G66)*U
     1                -G67)*U+G68)*U-G69)*U+G610)*U-G611
  523 TR8=  ((((((((((((( G81*U-G82)*U+G83)*U-G84)*U+G85)*U-G86)*U
     1                     + G87)*U-G88)*U+G89)*U-G810)*U+G811)*U
     2                     -G812)*U+G813)*U-G814)*U+G815
      PSI= PSI+TRE*(TR2+TRD*(TR4+TRD*(TR6+TRD*TR8)))
      PSIP = -RU*ETA2/X2
      TRB=TRE*X/U
      TR2= (3.D0*X-8.D0)/32.D0
  531 TR4= - (((63.D0*X-336.D0)*X+704.D0)*X-1536.D0)/2048.D0
  532 TR6 =  (((( GP61*X-GP62)*X+GP63)*X-GP64)*X+GP65)*X-GP66
  533 TR8 =  (((((( -GP81*X+GP82)*X-GP83)*X+GP84)*X-GP85)*X+GP86)*X
     1              +GP87)*X+GP88
      PSIP =PSIP+ TRB*(TR2+TRC*(TR4+TRC*(TR6+TRC*TR8)))
      TRA= DEXP(FI)
      FO= TRA*DSIN(PSI)
      GO= TRA*DCOS(PSI)
      IF(INDICE)535,536,535
  535 TRA=FO
      FO=GO
      GO=-TRA
  536 TRA=-ETA2/RAUC
      FPO=(FIP*FO+PSIP*GO)*TRA
      GPO=(FIP*GO-PSIP*FO)*TRA
      RETURN
      END
      SUBROUTINE JFLGAM(XD,YD,PAR,PAI,NBCHIF)
      DOUBLE PRECISION XD,YD,PAR,PAI,TEST,C,HLO2PI,PI,PIS2,PIS4
      DOUBLE PRECISION X,Y,U,V,TRA,TRA1,ALO2PI,RAC2,COSI,SINI
      DOUBLE PRECISION COS2I,SIN2I,ZMOD,DEPI
      DOUBLE PRECISION XX
      DOUBLE PRECISION ALOPI
      DOUBLE PRECISION SUPINT
      DIMENSION TEST(7),C(6)
      DATA TEST/2.9152D+7,2.2958D+3,1.4124D+2,3.9522D+1,19.6611D0,12.791
     1D0,-10.D0/
      DATA C/8.333333333333333D-2,-2.777777777777777D-3,7.936507936507
     1937D-4,-5.952380952380952D-4,8.417508417508418D-4,-1.917526917526
     2918D-3/
      DATA HLO2PI/0.918938533204672/,PI/3.141592653589793/
      DATA PIS2/1.570796326794897/,PIS4/0.785398163397448/
      DATA ALO2PI/1.837877066409345/,RAC2/0.3465735902799726/
      DATA DEPI/6.283185307179586/,ALOPI/1.1447298858494002/
      DATA SUPINT/2147483647.D0/
      NBCHIF=15
      X=DABS(XD)
      XX=X
      IF(YD)1,2,1
    1 Y=DABS(YD)
      KR=1
      I=DMOD(10.99D0-X,SUPINT)
C     TRANSLATION
      IF(I)3,3,4
    4 TRA=I
      X=X+TRA
C ****************************************************C
C     LOGARITHM (X+IY) (X,Y HAVE TO BE POSITIVE)      C
C ****************************************************C
    3 IF(X-Y)5,6,7
    5 TRA1=X/Y
      IF(TRA1)8,8,9
    8 U=DLOG(Y)
      V=PIS2
      GO TO 10
    6 U=RAC2+DLOG(X)
      V=PIS4
      GO TO 10
    9 TRA=Y*DSQRT(1.D0+TRA1*TRA1)
      TRA1=Y/X
   11 U=DLOG(TRA)
      V=DATAN(TRA1)
   10 GO TO (12,19,23),KR
    7 TRA1=Y/X
      TRA=X*DSQRT(1.D0+TRA1*TRA1)
      GO TO 11
   12 KR=2
C *************************************************C
C     ASYMPTOTIC EXPANSION FOR X LARGER THAN 10    C
C *************************************************C
      TRA=X-0.5D0
      PAI=V*TRA+Y*(U-1.D0)
      PAR=-X+HLO2PI+U*TRA-V*Y
      ZMOD=X+Y
      IF(ZMOD-TEST(1))13,13,14
   13 TRA=X*X+Y*Y
      COSI=X/TRA
      SINI=Y/TRA
      SIN2I=(SINI*COSI)+(SINI*COSI)
      COS2I=(COSI+SINI)*(COSI-SINI)
      K=1
      GO TO 15
   16 TRA=COSI*COS2I-SINI*SIN2I
      SINI=SINI*COS2I+COSI*SIN2I
      COSI=TRA
   15 PAR=PAR+C(K)*COSI
      PAI=PAI-C(K)*SINI
      K=K+1
      IF(ZMOD-TEST(K))16,16,14
C ****************************C
C     INVERSE TRANSLATION     C
C ****************************C
   17 I=I-1
      X=I
      X=XX+X
      GO TO 3
   19 PAR=PAR-U
      PAI=PAI-V
   14 IF(I-1)18,60,17
   60 IF(XD)17,61,17
C  **************************C
C     CONTROL OF QUADRANT    C
C  **************************C
   18 IF(XD)20,61,21
   61 TRA=PI*Y
      IF(TRA-1.D-2)300,300,301
  300 PAR= TRA*(2.D0+TRA*(-2.D0+TRA*(1.333333333333333D0+
     1 TRA*(  -0.6666666666666666D0+TRA*(   0.2666666666666666D0+
     2 TRA*( -0.08888888888888888D0+TRA*   0.02539682539682540D0 ))))))
      TRA1=-DLOG(Y)-DLOG(PAR)
      GO TO 302
  301 PAR=1.D0-DEXP(-TRA-TRA)
      TRA1=-DLOG(Y*PAR)
  302 PAR=0.5D0*(ALO2PI-TRA+TRA1)
      PAI=PAI-PIS2
   21 IF(YD)28,100,100
C     X+IY CHANGE EN -X-IY
   20 TRA=PI*Y
      PAR=ALO2PI-U-PAR-TRA
      PAI=PI-V-PAI
      TRA=DEXP(-TRA-TRA)
      X=PI*DMOD(X,2.D0)
      SINI=(1.D0-TRA)*DCOS(X)
      COSI=(1.D0+TRA)*DSIN(X)
      KR=3
      X=DABS(COSI)
      Y=DABS(SINI)
      GO TO 3
   23 IF(COSI)24,25,25
   24 V=PI-DSIGN(V,SINI)
      GO TO 26
   25 IF(SINI)27,26,26
   27 V=-V
   26 PAR=PAR-U
      PAI=PAI-V
      IF(YD)100,100,28
   28 PAI=-PAI
C     ARGUMENT DANS -PI,PI
  100 TRA=DABS(PAI/DEPI)
      IF(TRA-1.D+15)203,204,204
  204 NBCHIF=0
      PAI=0.D0
      GO TO 29
  203 IF(TRA-1.D0)205,205,206
  206 NBCHIF=DLOG10(TRA)
      NBCHIF=14-NBCHIF
      TRA=DMOD(TRA,SUPINT)
      PAI=DMOD(TRA,1.D0)*DSIGN(DEPI,PAI)
  205 IF(DABS(PAI)-PI)29,29,207
  207 PAI=PAI-DSIGN(DEPI,PAI)
      GO TO 29
C  *******************C
C     REAL JFLGAM     C
C  *******************C
    2 PAI=0.D0
      IF(XD)31,32,33
C ******************************C
C     CONDITIONS OF EXISTENCE   C
C ******************************C
   32 WRITE (6,1000)
 1000 FORMAT (1X,'JFLGAM(0) IS INFINITE')
      GO TO 50
   31 IF(X-4503599627370496.D0)34,35,35
   35 WRITE (6,1001)
 1001 FORMAT (1X,'ARGUMENT OF JFLGAM TOO LARGE')
      GO TO 50
   34 Y=DMOD(X,SUPINT)
      IF(Y)400,36,400
  400 IF(Y-0.99D0)33,33,405
  405 TRA=IDINT(Y+0.1D0)
      IF(DABS(Y-TRA)-5.D-15)36,36,33
   36 WRITE (6,1002)
 1002 FORMAT (1X,'JFLGAM IS INFINITE FOR NEGATIVE ARGUMENTS')
   50 PAR=1.D+74
      NBCHIF=0
      GO TO 29
C *******************C
C     TRANSLATION    C
C *******************C
   33 I=DMOD(10.99D0-X,SUPINT)
      IF(I)37,37,38
   38 TRA=I
      X=X+TRA
C  **************************C
C     ASYMPTOTIC EXPANSION   C
C  **************************C
   37 Y=DLOG(X)
      PAR=-X+HLO2PI+Y*(X-0.5D0)
      IF(X-TEST(1))39,39,43
   39 COSI=1.D0/X
      COS2I=COSI*COSI
      K=1
      GO TO 41
   42 COSI=COSI*COS2I
   41 PAR=PAR+C(K)*COSI
      K=K+1
      IF(X-TEST(K))42,42,40
C ***************************C
C     INVERSE TRANSLATION    C
C ***************************C
   44 X=X-1.D0
   48 Y=DLOG(X)
      PAR=PAR-Y
      I=I-1
   40 IF(I-1)43,49,44
   49 X=DABS(XD)
      GO TO 48
C *********************C
C     X IS NEGATIVE    C
C *********************C
   43 IF(XD)45,29,29
   45 PAR=ALOPI-PAR-Y
      Y=PI*DMOD(X,2.D0)
      Y=-DSIN(Y)
      IF(Y)46,36,47
   46 Y=-Y
      PAI=PI
   47 PAR=PAR-DLOG(Y)
   29 RETURN
      END
