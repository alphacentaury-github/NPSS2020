CX    ABQPDFPOT.  DFPOT - A PROGRAM FOR THE CALCULATION OF DOUBLE FOLDED      ABQP0000
CX    1   POTENTIALS.  J. COOK.                                               ABQP0000
CX    REF. IN COMP. PHYS. COMMUN. 25 (1982) 125                               ABQP0000
CX    /*PRIORITY 10                                                           ABQP0001
CX    //LMDFPOT JOB (NK00,LM,0-30,6),COOK                                     ABQP0002
CX    /*ROUTE  PRINT REMOTE23                                                 ABQP0003
CX    //  EXEC FHCLG,CPRINT=YES,REGION.G=100K                                ABQP0004
CX    //C.SYSIN DD *                                                          ABQP0005
C                                                                       ABQP0006
C     PROGRAMME TO CALCULATE A MICROSCOPIC POTENTIAL BY                 ABQP0007
C     DOUBLE FOLDING OF AN EFFECTIVE INTERACTION WITH BOTH              ABQP0008
C     TARGET AND PROJECTILE DENSITIES BY FOURIER TRANSFORMS.            ABQP0009
C                                                                       ABQP0010
C     WRITTEN BY JULIAN COOK, WHEATSTONE PHYSICS LABORATORY,            ABQP0011
C     KING S COLLEGE LONDON.      DECEMBER 1979                         ABQP0012
C                                                                       ABQP0013
C     VERSION 2    AUGUST 1980                                          ABQP0014
C     VERSION 3    MAY    1981                                          ABQP0015
C                                                                       ABQP0016
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0017
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0018
C                                                                       ABQP0019
C     DR,IMAX RADIAL INTERVAL AND NUMBER INTERVALS FOR FOLDED POTENTIAL ABQP0020
C     DQ,KMAX MOMENTUM INTERVAL AND NUMBER OF INTERVALS FOR             ABQP0021
C             FOURIER TRANSFORMS                                        ABQP0022
C     KPRINT  CONTROLS PRINT OUT                                        ABQP0023
C                                                                       ABQP0024
      COMMON/TEMPA/ARRAY1(1000),ARRAY2(1000)                            ABQP0025
CX      INTEGER PUNFMT(5),ITITLE(80),UNDER(80),BLANK                    ABQP0026
      INTEGER IADD(3)                                                   ABQP0027
      INTEGER KTRL
      REAL*8 MSR1,MSR2,MSRV,MSRU,MSRCHK                                 ABQP0028
      REAL*8 L1HAT,L2HAT,LUHAT                                          ABQP0029
      DIMENSION RHOFT1(1000),RHOFT2(1000),VFT(1000),UFT(1000),U(1000)   ABQP0030
C                                                                       ABQP0031
C     RHOFT1   FOURIER TRANSFORM OF DENSITY 1                           ABQP0032
C     RHOFT2   FOURIER TRANSFORM OF DENSITY 2                           ABQP0033
C     VFT      FOURIER TRANSFORM OF INTERACTION                         ABQP0034
C     UFT      FOURIER TRANSFORM OF FOLDED POTENTIAL                    ABQP0035
C     U        FOLDED POTENTIAL                                         ABQP0036
C                                                                       ABQP0037
      EQUIVALENCE (UFT(1000),ARRAY1(1000))                              ABQP0038
      EQUIVALENCE (U(1000),ARRAY2(1000))                                ABQP0039
      DATA LINE,BLANK/1H=,1H /                                          ABQP0040
      DATA IADD(1),IADD(2),IADD(3)/1HA,1HD,1HD/                         ABQP0041
C                                                                       ABQP0042
C     FORTRAN STREAM IP FOR PUNCHED OUTPUT                              ABQP0043
C     FORTRAN STREAM ID FOR DISC OUTPUT                                 ABQP0044
C                                                                       ABQP0045
      DATA IP,ID/7,8/                                                   ABQP0046

      OPEN(5,FILE='input.in',STATUS='OLD')             
      OPEN(6,FILE='out.dat',STATUS='UNKNOWN')             

      PI=3.14159265358979D0                                             ABQP0047
      SQ4PI=DSQRT(4.0D0*PI)                                             ABQP0048
      F1=SQ4PI                                                          ABQP0049
      F2=SQ4PI                                                          ABQP0050
      FU=SQ4PI                                                          ABQP0051
C                                                                       ABQP0052
C     CHANGE THIS CARD IF INAPPROPRIATE.  DIVERTS EXECUTION TO 999 IF   ABQP0053
C     END OF FILE IS FOUND                                              ABQP0054
C                                                                       ABQP0055
CX    1 READ(5,10,END=999)  ITITLE                                      ABQP0056
    1 READ(5,7,END=999)  KTRL                                           ABQP0056
CX      READ(5,7) KTRL
    7 FORMAT(1I10)
C                                                                       ABQP0063
C     PRINTS HEADING                                                    ABQP0064
C     
      WRITE(6,8) KTRL                                                   ABQP0065
    8 FORMAT(//,'KTRL     =     ',1I5)
      
	IF(KTRL.NE.0) GO TO 300

      WRITE(6,9)                                                        ABQP0066
    9 FORMAT(1H1,47(1H=)//                                              ABQP0067
     *10X,28HDDD   FFF  PPP    OO   TTTTT/                              ABQP0068
     *10X,28HD  D  F    P  P  O  O    T  /                              ABQP0069
     *10X,28HD  D  FFF  PPP   O  O    T  /                              ABQP0070
     *10X,28HD  D  F    P     O  O    T  /                              ABQP0071
     *10X,28HDDD   F    P      OO     T  //                             ABQP0072
     *7X,34HDOUBLE FOLDING POTENTIAL PROGRAMME//                        ABQP0073
     *7X,35HJULIAN COOK, KING'S COLLEGE, LONDON//                       ABQP0074
     *1H ,47(1H=))                                                      ABQP0075
C                                                                       ABQP0076
      READ(5,50)DROUT,RMAX,DQ,QMAX,KPRINT,KPUNCH                        ABQP0092
   50 FORMAT(4F10.0,2I5    )                                            ABQP0093
      LMAX=IDINT(RMAX/DROUT+0.5D0)                                      ABQP0094
      IF(LMAX.GT.1000) LMAX=1000                                        ABQP0095
      RMAX=DROUT*DFLOAT(LMAX)                                           ABQP0096
      IF(DQ.EQ.0.0D0) DQ=0.025D0                                        ABQP0097
      IF(QMAX.EQ.0.0D0) QMAX=3.0D0                                      ABQP0098
      KMAX=IDINT(QMAX/DQ+0.5D0)                                         ABQP0099
      IF(KMAX.GT.1000) KMAX=1000                                        ABQP0100
      QMAX=DQ*DFLOAT(KMAX)                                              ABQP0101
      WRITE(6,51)DROUT,RMAX                                             ABQP0102
   51 FORMAT(11H FOR OUTPUT,10X,5HDR = ,F10.3,13X,7HRMAX = ,F10.3)      ABQP0103
      WRITE(6,52)DQ,QMAX                                                ABQP0104
   52 FORMAT(15H FOR TRANSFORMS,6X,5HDQ = ,F10.3,13X,7HQMAX = ,F10.3)   ABQP0105
C                                                                       ABQP0108
C     READS MULTIPOLARITIES                                             ABQP0109
C                                                                       ABQP0110
      READ(5,54) L1,L2,LU                                               ABQP0111
   54 FORMAT(3I10)                                                      ABQP0112
      WRITE(6,55) L1,L2,LU                                              ABQP0113
   55 FORMAT(////16H MULTIPOLARITIES,5X,5HL1 = ,I2,5X,5HL2 = ,I2,       ABQP0114
     *  5X,5HLU = ,I2)                                                  ABQP0115
      IF(L1.GT.0) F1=1.0D0                                              ABQP0116
      IF(L2.GT.0) F2=1.0D0                                              ABQP0117
      IF(LU.GT.0) FU=1.0D0                                              ABQP0118
      L1HAT=DSQRT(DFLOAT(2*L1+1))                                       ABQP0119
      L2HAT=DSQRT(DFLOAT(2*L2+1))                                       ABQP0120
      LUHAT=DSQRT(DFLOAT(2*LU+1))                                       ABQP0121
      C=CLEB(L1,L2,LU)                                                  ABQP0122
      WRITE(6,56) C                                                     ABQP0123
   56 FORMAT(31H0CLEBSCH-GORDAN COEFFICIENT IS ,F10.3)                  ABQP0124
      C=C*L1HAT*L2HAT*F1*F2/(LUHAT*FU*SQ4PI)                            ABQP0125
      LL=(L1-L2-LU)/2                                                   ABQP0126
      IF((LL/2)*2.NE.LL) C=-C                                           ABQP0127
C                                                                       ABQP0128
C     CALCULATES FOURIER TRANSFORM OF DENSITY DISTRIBUTION 1            ABQP0129
C                                                                       ABQP0130
      WRITE(6,60)                                                       ABQP0131
   60 FORMAT(////23H DENSITY DISTRIBUTION 1/1H ,22(1H=))                ABQP0132
      CALL FUN (L1,RHOFT1,VOL1,MSR1)                                    ABQP0133

C                                                                       ABQP0134
C     CALCULATES FOURIER TRANSFORM OF DENSITY DISTRIBUTION 2            ABQP0135
C                                                                       ABQP0136
      WRITE(6,70)                                                       ABQP0137
   70 FORMAT(////23H DENSITY DISTRIBUTION 2/1H ,22(1H=))                ABQP0138
      CALL FUN (L2,RHOFT2,VOL2,MSR2)                                    ABQP0139
CX      CALL TRANSF1 (L2,1,RHOFT2,DENST,DR,IMAX,DQ,KMAX)
C                                                                       ABQP0140
C     CALCULATES FOURIER TRANSFORM OF EFFECTIVE INTERACTION             ABQP0141
C                                                                       ABQP0142

      WRITE(6,80)                                                       ABQP0143
   80 FORMAT(////22H EFFECTIVE INTERACTION/1H ,21(1H=))                 ABQP0144
      CALL FUN (0,VFT,VOLV,MSRV)                                        ABQP0145
C                                                                       ABQP0146
C     CALCULATE FOURIER TRANSFORM OF FOLDED POTENTIAL                   ABQP0147
C                                                                       ABQP0148
      DO 100 K=1,KMAX                                                   ABQP0149
  100 UFT(K)=C*RHOFT1(K)*RHOFT2(K)*VFT(K)                               ABQP0150
C                                                                       ABQP0151
C     PRINTS FOURIER TRANSFORMS                                         ABQP0152
C                                                                       ABQP0153
      IF(KPRINT.LT.2) GO TO 110                                         ABQP0154
       WRITE(6,101)                                                     ABQP0155
  101 FORMAT(////19H FOURIER TRANSFORMS/1H ,18(1H-)//7X,1HQ,12X,        ABQP0156
     *9HRHOFT1(Q),11X,9HRHOFT2(Q),12X,6HVFT(Q),14X,6HUFT(Q)/)           ABQP0157
      DO 102 K=1,KMAX                                                   ABQP0158
      Q=DQ*FLOAT(K)                                                     ABQP0159
  102 WRITE(6,103)Q,RHOFT1(K),RHOFT2(K),VFT(K),UFT(K)                   ABQP0160
  103 FORMAT(F10.2,4(5X,1PE15.5))                                       ABQP0161
      WRITE(6,104)                                                      ABQP0162
  104 FORMAT(////)                                                      ABQP0163
  110 WRITE(6,90)                                                       ABQP0164
   90 FORMAT(////17H FOLDED POTENTIAL/1H ,16(1H=))                      ABQP0165
      DR=DROUT                                                          ABQP0166
      IMAX=LMAX                                                         ABQP0167
C                                                                       ABQP0168
C     TAKE INVERSE TRANSFORM TO GET FOLDED POTENTIAL                    ABQP0169
C                                                                       ABQP0170
      CALL TRANSF (LU,1)                                                ABQP0171
      DO 120 I=1,IMAX                                                   ABQP0172
  120 ARRAY1(I)=ARRAY2(I)                                               ABQP0173
C                                                                       ABQP0174
C     FIND VOLUME INTEGRALS OF FOLDED POTENTIAL                         ABQP0175
C                                                                       ABQP0176
      CALL INTEGR (LU,U,VOLU,VOL2U)                                     ABQP0177
      MSRU=VOL2U/VOLU                                                   ABQP0178
      WRITE(6,91) VOLU,MSRU                                             ABQP0179
   91 FORMAT(7H VOL(L),14X,1PD11.4/16H VOL(L+2)/VOL(L),5X,0PF11.3)      ABQP0180
C                                                                       ABQP0181
C     CHECK CONSISTENCY OF INTEGRAL PROPERTIES                          ABQP0182
C                                                                       ABQP0183
      C=C*FACT2(LU)/(FACT2(L1)*FACT2(L2))                               ABQP0184
      VOLCHK=C*VOL1*VOL2*VOLV                                           ABQP0185
      DVOL=100.0D0*(VOLU-VOLCHK)/VOLCHK                                 ABQP0186
      MSRCHK=DFLOAT(2*LU+3)*(MSR1/DFLOAT(2*L1+3)+MSR2/DFLOAT(2*L2+3)    ABQP0187
     *    +MSRV/3.0D0)                                                  ABQP0188
      DMSR=100.0D0*(MSRU-MSRCHK)/MSRCHK                                 ABQP0189
      WRITE(6,92)                                                       ABQP0190
   92 FORMAT(////19H CONSISTENCY CHECKS/1H ,18(1H=)/                    ABQP0191
     *27X,9HPOTENTIAL,8X,5HCHECK,9X,8HPCT DIFF)                         ABQP0192
      WRITE(6,93)VOLU,VOLCHK,DVOL                                       ABQP0193
   93 FORMAT(7H VOL(L),18X,1PD11.4,4X,D11.4,0PF15.3)                    ABQP0194
      WRITE(6,94)MSRU,MSRCHK,DMSR                                       ABQP0195
   94 FORMAT(16H VOL(L+2)/VOL(L),5X,3F15.3)                             ABQP0196
      WRITE(6,95)                                                       ABQP0197
   95 FORMAT(33H0NOTE ONLY VALID FOR L1 + L2 = LU)                      ABQP0198
C                                                                       ABQP0199
C     PRINT POTENTIAL.  ALSO WRITE ON DISC AND PUNCH IF NEEDED          ABQP0200
C                                                                       ABQP0201
      WRITE(6,125)                                                      ABQP0202
  125 FORMAT(/9X,4(1HR,12X,4HV(R),11X)/)                                ABQP0203
      CALL PRINT                                                        ABQP0204
      IF(KPUNCH.EQ.0) GO TO 210                                         ABQP0205
      IF(KPUNCH.GE.1)WRITE(ID)(U(I),I=1,LMAX)                           ABQP0206
      IF(KPUNCH.EQ.2)WRITE(IP,201)(U(I),I=1,LMAX)                       ABQP0207
  201 FORMAT(5E15.6)
  210 CONTINUE                                                          ABQP0208

      DO N=1,IMAX
	RR=DR*DFLOAT(N)
      WRITE(9,211)RR,U(N)
  211 FORMAT(4E15.6)
      END DO
      
C                                                                       ABQP0209
C     RETURN TO START NEW CALCULATION                                   ABQP0210
C                                                                       ABQP0211
      GO TO 1                                                           ABQP0212
C                                                                       ABQP0213
C     ADDS LAST CALCULATED POTENTIAL TO POTENTIAL STORED ON DISC        ABQP0214
C                                                                       ABQP0215
  300 WRITE(6,310)                                                      ABQP0216
  310 FORMAT(16H1TOTAL POTENTIAL/1H ,15(1H=)//9X,4(1HR,12X,4HV(R),      ABQP0217
     *  11X)/)                                                          ABQP0218
      READ(5,315)C1,C2                                                  ABQP0219
  315 FORMAT(2F10.0)                                                    ABQP0220
      REWIND ID                                                         ABQP0221
      READ(ID)(ARRAY1(I),I=1,LMAX)                                      ABQP0222
      DO 320 I=1,LMAX                                                   ABQP0223
      U(I)=C1*ARRAY1(I)+C2*U(I)                                         ABQP0224
  320 ARRAY1(I)=U(I)                                                    ABQP0225
      REWIND ID                                                         ABQP0226
      WRITE(ID)(U(I),I=1,LMAX)                                          ABQP0227
      CALL PRINT                                                        ABQP0228
C                                                                       ABQP0229
C     FIND VOLUME INTEGRALS OF TOTAL POTENTIAL                          ABQP0230
C                                                                       ABQP0231
      CALL INTEGR (LU,U,VOLU,VOL2U)                                     ABQP0232
      MSRU=VOL2U/VOLU                                                   ABQP0233
      WRITE(6,104)                                                      ABQP0234
      WRITE(6,91)VOLU,MSRU                                              ABQP0235
      IF(KPUNCH.LT.2)GO TO 1                                            ABQP0236
C                                                                       ABQP0239
C     RETURN TO START NEW CALCULATION                                   ABQP0240
C                                                                       ABQP0241
      GO TO 1                                                           ABQP0242
 999  STOP                                                              ABQP0243
      END                                                               ABQP0244

      SUBROUTINE FUN (L,FT,VOL,MSR)                                     ABQP0245
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0246
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0247
      COMMON/FFF/FF(1000)                                               ABQP0248
      COMMON/TEMPA/F(1000),ARRAY2(1000)                                 ABQP0249
      COMMON/DEL/V                                                      ABQP0250
      EQUIVALENCE (ARRAY2(1000),FUNC(1000))                             ABQP0251
      DIMENSION FT(1000),FUNC(1000)                                     ABQP0252
C                                                                       ABQP0253
C     FF    INDIVIDUAL CONTRIBUTION TO FUNCTION                         ABQP0254
C     F     TOTAL FUNCTION                                              ABQP0255
C     FT    FOURIER TRANSFORM OF FUNCTION                               ABQP0256
C     V     STRENGTH OF DELTA TERM                                      ABQP0257
C     ARRAY2  TEMPORARY STORAGE                                         ABQP0258
C                                                                       ABQP0259
      REAL*8 MSR,MSRF                                                   ABQP0260
C                                                                       ABQP0261
C     INITIALISE QUANTITIES                                             ABQP0262
C                                                                       ABQP0263
      V=0.0D0                                                           ABQP0264
      IDEL=0                                                            ABQP0265
      FOURPI=4.0D0*PI                                                   ABQP0266
      DO 10 I=1,1000                                                    ABQP0267
      ARRAY2(I)=0.0D0                                                   ABQP0268
      F(I)=0.0D0                                                        ABQP0269
   10 FT(I)=0.0D0                                                       ABQP0270
      VOL=0.0D0                                                         ABQP0271
      VOL2=0.0D0                                                        ABQP0272
C                                                                       ABQP0273
C     READ INTERVALS,MAXIMUM RADIUS,CORRECTION FOR PROTON SIZE,         ABQP0274
C     DENSITY DEPENDENCE FACTOR                                         ABQP0275
C                                                                       ABQP0276
      READ(5,50)DR,RMAX,ALPHA,BETA                                      ABQP0277
   50 FORMAT(4F10.0)                                                    ABQP0278
      IMAX=IDINT(RMAX/DR+0.5D0)                                         ABQP0279
      IF(IMAX.GT.1000) IMAX=1000                                        ABQP0280
      RMAX=DR*DFLOAT(IMAX)                                              ABQP0281
      WRITE(6,60)DR,RMAX,ALPHA,BETA                                     ABQP0282
   60 FORMAT(6H DR = ,F10.3,5X,7HRMAX = ,F10.3,5X,8HALPHA = ,F10.3,     ABQP0283
     *5X,7HBETA = ,F10.3)                                               ABQP0284
C                                                                       ABQP0285
C     CALCULATE TRANSITION DENSITY OR INTERACTION                       ABQP0286
C                                                                       ABQP0287
   70 READ(5,80)IOPT,A,SIGN                                             ABQP0288
   80 FORMAT(I10,2F10.0)                                                ABQP0289
      IF(IOPT.LT.0) IADD=0                                              ABQP0290
      IF(IOPT.GT.0) IADD=1                                              ABQP0291
      IOPT=IABS(IOPT)                                                   ABQP0292
      WRITE(6,90)IOPT,A,SIGN                                            ABQP0293
   90 FORMAT(/8H IOPT = ,I2,10X,4HA = ,F10.3,5X,7HSIGN = ,F10.3/        ABQP0294
     *1H ,9(1H-))                                                       ABQP0295
      IF(L.GT.0) A=A*FOURPI                                             ABQP0296
C                                                                       ABQP0297
C     CALL APPROPRIATE SUBROUTINE TO CALCULATE CONTRIBUTION             ABQP0298
C     TO THE FUNCTION                                                   ABQP0299
C                                                                       ABQP0300

      IF(IOPT.EQ.1) CALL GAUSS                                          ABQP0301
      IF(IOPT.EQ.2) CALL YUKAWA                                         ABQP0302
      IF(IOPT.EQ.3) CALL FERMI                                          ABQP0303
      IF(IOPT.EQ.4) CALL SAXON                                          ABQP0304
      IF(IOPT.EQ.5) CALL CARDIN                                         ABQP0305
      IF(IOPT.EQ.6) CALL USER                                           ABQP0306
      IF(IOPT.EQ.7) CALL DELTA  (L)                                     ABQP0307
      IF(IOPT.EQ.8) CALL TASSIE (L,VOL,VOL2)                            ABQP0308
      IF(IOPT.EQ.9) CALL HARMONIC                                       ABQP0308-1
      IF(IOPT.EQ.10)CALL USER1                                          ABQP0308-2
      IF(IOPT.EQ.11)CALL USER2                                          ABQP0308-3
      IF(IOPT.EQ.14)CALL USER3                                          ABQP0308-3
      IF(IOPT.EQ.12)CALL FOURBESS                                       ABQP0308-4
      IF(IOPT.EQ.13)CALL GAUSSIAN                                       ABQP0308-5
      IF(IOPT.NE.7) IDEL=1                                              ABQP0309
      IF(IOPT.EQ.7) GOTO 91                                             ABQP0310

C                                                                       ABQP0311
C     FIND VOLUME INTEGRALS OF CONTRIBUTION                             ABQP0312
C                                                                       ABQP0313
      CALL INTEGR (L,FF,VOLF,VOL2F)                                     ABQP0314
      GOTO 92                                                           ABQP0315
   91 VOLF=V                                                            ABQP0316
      VOL2F=0.0D0                                                       ABQP0317
   92 MSRF=VOL2F/VOLF                                                   ABQP0318
      FACT=1.0D0                                                        ABQP0319
C                                                                       ABQP0320
C     NORMALISE CONTRIBUTION                                            ABQP0321
C                                                                       ABQP0322
      IF(A.EQ.0.0D0)GO TO 97                                            ABQP0323
      FACT=A/VOLF                                                       ABQP0324
      VOLF=A                                                            ABQP0325
      VOL2F=FACT*VOL2F                                                  ABQP0326
      DO 96 I=1,IMAX                                                    ABQP0327
   96 FF(I)=FACT*FF(I)                                                  ABQP0328
C                                                                       ABQP0329
C     CHANGE SIGN OF CONTRIBUTION                                       ABQP0330
C                                                                       ABQP0331
   97 IF(SIGN.GE.0.0D0) GO TO 99                                        ABQP0332
      VOLF=-VOLF                                                        ABQP0333
      FACT=-FACT                                                        ABQP0334
      DO 98 I=1,IMAX                                                    ABQP0335
   98 FF(I)=-FF(I)                                                      ABQP0336
C                                                                       ABQP0337
C     ADD THIS CONTRIBUTION TO THE FUNCTION                             ABQP0338
C                                                                       ABQP0339
   99 DO 100 I=1,IMAX                                                   ABQP0340
  100 F(I)=F(I)+FF(I)                                                   ABQP0341
      IF(IOPT.EQ.7) V=FACT*V                                            ABQP0342
      VOL=VOL+VOLF                                                      ABQP0343
      VOL2=VOL2+VOL2F                                                   ABQP0344
      WRITE(6,105)FACT,VOLF,MSRF                                        ABQP0345
  105 FORMAT(8H FACT = ,1PD10.3,5X,9HVOL(L) = ,0PF10.3,5X,              ABQP0346
     *18HVOL(L+2)/VOL(L) = ,F10.3)                                      ABQP0347
C                                                                       ABQP0348
C     GO BACK TO CALCULATE MORE CONTRIBUTIONS IF REQUIRED               ABQP0349
C                                                                       ABQP0350
      IF(IADD.EQ.0) GO TO 70                                            ABQP0351
      IF(KPRINT.EQ.0) GOTO 109                                          ABQP0352
      WRITE(6,106)                                                      ABQP0353
  106 FORMAT(//9H FUNCTION/1H ,8(1H-)//9X,4(1HR,12X,4HF(R),11X)/)       ABQP0354
      CALL PRINT                                                        ABQP0355
      WRITE(6,107)                                                      ABQP0356
  107 FORMAT(////)                                                      ABQP0357
  109 CONTINUE                                                          ABQP0358
      IF(BETA.EQ.0.0D0) GOTO 120                                        ABQP0359
C                                                                       ABQP0360
C     CALCULATE FUNCTION TO BE USED IN EXPONENT FOR DENSITY             ABQP0361
C     DEPENDENCE.   SAME PROCEDURE AS ABOVE                             ABQP0362
C                                                                       ABQP0363
      WRITE(6,190)                                                      ABQP0364
  190 FORMAT(//36H FOLLOWING FUNCTION USED IN EXPONENT,                 ABQP0365
     *23H FOR DENSITY DEPENDENCE/1X,58(1H-)/)                           ABQP0366
      DO 200 I=1,IMAX                                                   ABQP0367
      FUNC(I)=F(I)                                                      ABQP0368
  200 F(I)=0.0D0                                                        ABQP0369
      VOL=0.0D0                                                         ABQP0370
      VOL2=0.0D0                                                        ABQP0371
  210 READ(5,80)IOPT,A,SIGN                                             ABQP0372
      IF(IOPT.LT.0) IADD=0                                              ABQP0373
      IF(IOPT.GT.0) IADD=1                                              ABQP0374
      IOPT=IABS(IOPT)                                                   ABQP0375
      WRITE(6,90)IOPT,A,SIGN                                            ABQP0376

      IF(IOPT.EQ.1) CALL GAUSS                                          ABQP0377
      IF(IOPT.EQ.2) CALL YUKAWA                                         ABQP0378
      IF(IOPT.EQ.3) CALL FERMI                                          ABQP0379
      IF(IOPT.EQ.4) CALL SAXON                                          ABQP0380
      IF(IOPT.EQ.5) CALL CARDIN                                         ABQP0381
      IF(IOPT.EQ.6) CALL USER                                           ABQP0382
      IF(IOPT.EQ.8) CALL TASSIE (0,VOL,VOL2)                            ABQP0383
      IF(IOPT.EQ.9) CALL HARMONIC                                       ABQP0384-1
      IF(IOPT.EQ.10) GOTO 265                                           ABQP0384

      CALL INTEGR(0,FF,VOLF,VOL2F)                                      ABQP0385
      MSRF=VOL2F/VOLF                                                   ABQP0386
      FACT=1.0D0                                                        ABQP0387
      IF(A.EQ.0.0D0) GOTO 230                                           ABQP0388
      FACT=A/VOLF                                                       ABQP0389
      VOLF=A                                                            ABQP0390
      VOL2F=FACT*VOL2F                                                  ABQP0391
      DO 220 I=1,IMAX                                                   ABQP0392
  220 FF(I)=FACT*FF(I)                                                  ABQP0393
  230 IF(SIGN.GE.0.0D0) GOTO 250                                        ABQP0394
      VOLF=-VOLF                                                        ABQP0395
      FACT=-FACT                                                        ABQP0396
      DO 240 I=1,IMAX                                                   ABQP0397
  240 FF(I)=-FF(I)                                                      ABQP0398
  250 DO 260 I=1,IMAX                                                   ABQP0399
  260 F(I)=F(I)+FF(I)                                                   ABQP0400
      VOL=VOL+VOLF                                                      ABQP0401
      VOL2=VOL2+VOL2F                                                   ABQP0402
      WRITE(6,105)FACT,VOLF,MSRF                                        ABQP0403
      IF(IADD.EQ.0) GOTO 210                                            ABQP0404
  265 IF(IOPT.NE.9) GOTO 269                                            ABQP0405
      DO 266 I=1,IMAX                                                   ABQP0406
  266 F(I)=FUNC(I)                                                      ABQP0407
      WRITE(6,267)                                                      ABQP0408
  267 FORMAT(49H FUNCTION IN EXPONENT IDENTICAL TO FUNCTION ABOVE)      ABQP0409
  269 DO 270 I=1,IMAX                                                   ABQP0410
  270 F(I)=FUNC(I)*DEXP(-BETA*F(I))                                     ABQP0411
C                                                                       ABQP0412
C     CALCULATION OF COMPLETE FUNCTION NOW FINISHED                     ABQP0413
C                                                                       ABQP0414
      IF(KPRINT.EQ.0) GOTO 115                                          ABQP0415
      WRITE(6,111)                                                      ABQP0416
  111 FORMAT(//18H MODIFIED FUNCTION/1H ,17(1H-)//9X,                   ABQP0417
     *4(1HR,12X,4HF(R),11X)/)                                           ABQP0418
      CALL PRINT                                                        ABQP0419
      WRITE(6,107)                                                      ABQP0420
  115 CONTINUE                                                          ABQP0421
      CALL INTEGR (L,F,VOL,VOL2)                                        ABQP0422
C                                                                       ABQP0423
C     CALCULATE FOURIER TRANSFORM OF FUNCTION                           ABQP0424
C                                                                       ABQP0425
      DO 119 K=1,KMAX                                                   ABQP0426
  119 ARRAY2(K)=0.0D0                                                   ABQP0427
  120 IF(IDEL.EQ.0) GOTO 121                                            ABQP0428
      CALL TRANSF (L,0)                                                 ABQP0429
  121 DO 125 K=1,KMAX                                                   ABQP0430
  125 FT(K)=ARRAY2(K)+V                                                 ABQP0431
      IF(ALPHA.EQ.0.0D0) GOTO 140                                       ABQP0432
C                                                                       ABQP0433
C     CORRECT FOR FINITE PROTON SIZE                                    ABQP0434
C                                                                       ABQP0435
      WRITE(6,126)                                                      ABQP0436
  126 FORMAT(41H0DENSITY CORRECTED FOR FINITE PROTON SIZE/)             ABQP0437
      ALPHA2=ALPHA**2                                                   ABQP0438
      ALPHA4=ALPHA**4                                                   ABQP0439
      DO 130 K=1,KMAX                                                   ABQP0440
      Q=DQ*DFLOAT(K)                                                    ABQP0441
  130 FT(K)=FT(K)*(ALPHA2+Q*Q)**2/ALPHA4                                ABQP0442
      VOL2=DABS(VOL2)-DABS(12.D0*VOL/ALPHA2)                            ABQP0443
  140 MSR=DABS(VOL2/VOL)                                                ABQP0444
      WRITE(6,150) VOL,MSR                                              ABQP0445
  150 FORMAT(13H0TOTAL VOL(L),14X,F10.3/                                ABQP0446
     *22H TOTAL VOL(L+2)/VOL(L),5X,F10.3)                               ABQP0447
      RETURN                                                            ABQP0448
      END                                                               ABQP0449


      SUBROUTINE GAUSS                                                  ABQP0450
C                                                                       ABQP0451
C     GAUSSIAN*R*PWR FUNCTION                                           ABQP0452
C                                                                       ABQP0453
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0454
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0455
      COMMON/FFF/FF(1000)                                               ABQP0456
      WRITE(6,30)                                                       ABQP0457
   30 FORMAT(18H GAUSSIAN FUNCTION)                                     ABQP0458
      READ(5,50)C,ALPHA,PWR                                             ABQP0459
   50 FORMAT(3F10.0)                                                    ABQP0460
      IF(C.EQ.0.0D0) C=1.0D0                                            ABQP0461
      WRITE(6,60)C,ALPHA,PWR                                            ABQP0462
   60 FORMAT(5H C = ,F10.3,5X,8HALPHA = ,F10.3,5X,6HPWR = ,F10.3)       ABQP0463
      DO 100 I=1,IMAX                                                   ABQP0464
      R=DR*DFLOAT(I)                                                    ABQP0465
  100 FF(I)=C*R**PWR*DEXP(-(R/ALPHA)**2)                                ABQP0466   
      RETURN                                                            ABQP0467
      END                                                               ABQP0468
      SUBROUTINE YUKAWA                                                 ABQP0469
C                                                                       ABQP0470
C     YUKAWA FUNCTION                                                   ABQP0471
C                                                                       ABQP0472
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0473
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0474
      COMMON/FFF/FF(1000)                                               ABQP0475
      WRITE(6,30)                                                       ABQP0476
   30 FORMAT(16H YUKAWA FUNCTION)                                       ABQP0477
      READ(5,50)V,BETA                                                  ABQP0478
   50 FORMAT(2F10.0)                                                    ABQP0479
      IF(V.EQ.0.0D0) V=1.0D0                                            ABQP0480
      WRITE(6,60)V,BETA                                                 ABQP0481
   60 FORMAT(5H V = ,F10.3,5X,7HBETA = ,F10.3)                          ABQP0482
      DO 100 I=1,IMAX                                                   ABQP0483
      R=DR*DFLOAT(I)                                                    ABQP0484
      BR=BETA*R                                                         ABQP0485
  100 FF(I)=V*DEXP(-BR)/BR                                              ABQP0486
      RETURN                                                            ABQP0487
      END                                                               ABQP0488

      SUBROUTINE FERMI                                                  ABQP0489
C                                                                       ABQP0490
C     3 PARAMETER FERMI FUNCTION                                        ABQP0491
C                                                                       ABQP0492
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0493
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0494
      COMMON/FFF/FF(1000)                                               ABQP0495
      WRITE(6,30)                                                       ABQP0496
   30 FORMAT(27H 3 PARAMETER FERMI FUNCTION)                            ABQP0497
      READ(5,50)C,RR,AR,W                                               ABQP0498
   50 FORMAT(4F10.0)                                                    ABQP0499
      IF(C.EQ.0.0D0) C=1.0D0                                            ABQP0500
      WRITE(6,60)C,RR,AR,W                                              ABQP0501
   60 FORMAT(5H C = ,F10.3,5X,5HRR = ,F10.3,5X,5HAR = ,F10.3,5X,        ABQP0502
     *4HW = ,F10.3)                                                     ABQP0503

      DO 100 I=1,IMAX                                                   ABQP0504
      R=DR*DFLOAT(I)                                                    ABQP0505
      FF(I)=C*(1.0D0+W*R*R/RR**2)/(1.0D0+DEXP((R-RR)/AR))               ABQP0506
  100	CONTINUE

      RETURN                                                            ABQP0507
      END                                                               ABQP0508

      SUBROUTINE SAXON                                                  ABQP0509
C                                                                       ABQP0510
C     SAXON-WOODS AND DERIVATIVE SAXON-WOODS FUNCTIONS RAISED TO POWER  ABQP0511
C                                                                       ABQP0512
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0513
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0514
      COMMON/FFF/FF(1000)                                               ABQP0515
      READ(5,50)C,RR,AR,PWR,IDER                                        ABQP0516
   50 FORMAT(4F10.0,I10)                                                ABQP0517
      IF(C.EQ.0.0D0) C=1.0D0                                            ABQP0518
      IF(PWR.EQ.0.0D0) PWR=1.0D0                                        ABQP0519
      IF(IDER.NE.0) GOTO 90                                             ABQP0520
      WRITE(6,60)                                                       ABQP0521
   60 FORMAT(21H SAXON-WOODS FUNCTION)                                  ABQP0522
      WRITE(6,70)C,RR,AR,PWR                                            ABQP0523
   70 FORMAT(5H C = ,F10.3,5X,5HRR = ,F10.3,5X,5HAR = ,F10.3,5X,        ABQP0524
     *6HPWR = ,F10.3)                                                   ABQP0525
      DO 80 I=1,IMAX                                                    ABQP0526
      R=DR*DFLOAT(I)                                                    ABQP0527
   80 FF(I)=C/(1.0D0+DEXP((R-RR)/AR))**PWR                              ABQP0528
      RETURN                                                            ABQP0529
   90 WRITE(6,100)                                                      ABQP0530
  100 FORMAT(32H DERIVATIVE SAXON-WOODS FUNCTION)                       ABQP0531
      WRITE(6,70)C,RR,AR,PWR                                            ABQP0532
      PWRC=PWR*C/AR                                                     ABQP0533
      PWR1=PWR+1.0D0                                                    ABQP0534
      DO 110 I=1,IMAX                                                   ABQP0535
      R=DR*DFLOAT(I)                                                    ABQP0536
      EXPR=DEXP((R-RR)/AR)                                              ABQP0537
  110 FF(I)=-PWRC*EXPR/(1.D0+EXPR)**PWR1                                ABQP0538
      RETURN                                                            ABQP0539
      END                                                               ABQP0540
      SUBROUTINE CARDIN                                                 ABQP0541
C                                                                       ABQP0542
C     SUBROUTINE TO READ FUNCTION FROM CARDS                            ABQP0543
C     INCLUDES INTERPOLATION PROCEDURE                                  ABQP0544
C                                                                       ABQP0545
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0546
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0547
      COMMON/TEMPA/F(1000),RHO(1000)                                    ABQP0548
      COMMON/FFF/RHOO(1000)                                             ABQP0549
C                                                                       ABQP0550
C     RHO   FUNCTION READ IN FROM CARDS                                 ABQP0551
C     RHOO  INTERPOLATED FUNCTION RETURNED BY SUBROUTINE                ABQP0552
C                                                                       ABQP0553
      DIMENSION INFMT(5)                                                ABQP0554
      READ(5,10)INFMT,DRIN,NVAL,ISTART                                  ABQP0555
   10 FORMAT(5A4,F10.0,2I10)                                            ABQP0556
      WRITE(6,20)INFMT,DRIN                                             ABQP0557
   20 FORMAT(11H CARD INPUT/                                            ABQP0558
     *26H FUNCTION READ IN FORMAT  ,5A4,5X,12HAT INTERVALS,F10.3)       ABQP0559
      READ(5,INFMT)(RHO(I),I=1,NVAL)                                    ABQP0560
      IF(ISTART.EQ.0) GO TO 23                                          ABQP0561
      N1=NVAL                                                           ABQP0562
      N2=NVAL+1                                                         ABQP0563
      DO 22 N=1,NVAL                                                    ABQP0564
      RHO(N2)=RHO(N1)                                                   ABQP0565
      N1=N1-1                                                           ABQP0566
   22 N2=N2-1                                                           ABQP0567
      RHO(1)=RHO(2)                                                     ABQP0568
      NVAL=NVAL+1                                                       ABQP0569
   23 CONTINUE                                                          ABQP0570
      IF(DR.EQ.DRIN) GOTO 165                                           ABQP0571
      WRITE(6,30) DR                                                    ABQP0572
   30 FORMAT(37H FUNCTION INTERPOLATED AT INTERVALS  ,F10.3)            ABQP0573
      H2=2.0*DRIN*DRIN                                                  ABQP0574
      H3=H2*DRIN                                                        ABQP0575
      N7=NVAL                                                           ABQP0576
      N6=N7-1                                                           ABQP0577
      N5=N6-1                                                           ABQP0578
      R7=DRIN*DFLOAT(N6)                                                ABQP0579
      R6=R7-DRIN                                                        ABQP0580
      R5=R6-DRIN                                                        ABQP0581
      I2=2                                                              ABQP0582
      DO 100 J=1,IMAX                                                   ABQP0583
      RX=DR*DFLOAT(J)                                                   ABQP0584
      IF(RX.LT.DRIN) GO TO 35                                           ABQP0585
      IF(RX.LT.R6) GO TO 40                                             ABQP0586
      IF(RX.GT.R7) GO TO 60                                             ABQP0587
      RHOO(J)=(RHO(N5)*(RX-R6)*(RX-R7)-RHO(N6)*(RX-R5)*(RX-R7)*2.0      ABQP0588
     *        +RHO(N7)*(RX-R5)*(RX-R6))/H2                              ABQP0589
      GO TO 100                                                         ABQP0590
   35 RHOO(J)=(RHO(2)*(RX-2.0D0*DRIN)*(RX-3.0D0*DRIN)-RHO(3)*(RX-DRIN)  ABQP0591
     *       *(RX-3.0D0*DRIN)*2.D0+RHO(4)*(RX-DRIN)*(RX-2.D0*DRIN))/H2  ABQP0592
      GOTO 100                                                          ABQP0593
   40 X2=DFLOAT(I2-2)*DRIN                                              ABQP0594
      I1=I2                                                             ABQP0595
      DO 50 I=I1,N5                                                     ABQP0596
      I2=I                                                              ABQP0597
      X1=X2                                                             ABQP0598
      X2=X1+DRIN                                                        ABQP0599
      X3=X2+DRIN                                                        ABQP0600
      X4=X3+DRIN                                                        ABQP0601
      IF(RX.LT.X2) GO TO 50                                             ABQP0602
      IF(RX.GT.X3) GO TO 50                                             ABQP0603
      RHOO(J)=(-RHO(I-1)*(RX-X2)*(RX-X3)*(RX-X4)/3.0                    ABQP0604
     *        +RHO(I)*(RX-X1)*(RX-X3)*(RX-X4)                           ABQP0605
     *        -RHO(I+1)*(RX-X1)*(RX-X2)*(RX-X4)                         ABQP0606
     *        +RHO(I+2)*(RX-X1)*(RX-X2)*(RX-X3)/3.0)/H3                 ABQP0607
      GO TO 100                                                         ABQP0608
   50 CONTINUE                                                          ABQP0609
   60 J1=J                                                              ABQP0610
      GO TO 150                                                         ABQP0611
  100 CONTINUE                                                          ABQP0612
      GO TO 200                                                         ABQP0613
  150 DO 160 J=J1,IMAX                                                  ABQP0614
  160 RHOO(J)=0.0                                                       ABQP0615
      GO TO 200                                                         ABQP0616
  165 JMAX=NVAL-1                                                       ABQP0617
      DO 170 J=1,JMAX                                                   ABQP0618
  170 RHOO(J)=RHO(J+1)                                                  ABQP0619
  200 RETURN                                                            ABQP0620
      END                                                               ABQP0621
      SUBROUTINE USER                                                   ABQP0622
C                                                                       ABQP0623
C     USER DEFINED FUNCTION                                             ABQP0624
C                                                                       ABQP0625
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0626
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0627
      COMMON/FFF/FF(1000)                                               ABQP0628
      WRITE(6,30)                                                       ABQP0629
   30 FORMAT(22H USER DEFINED FUNCTION)                                 ABQP0630

      READ(5,50)AA,BB,CC                                         
   50 FORMAT(4F10.0)                                                 
      WRITE(6,60)                                                        
   60 FORMAT(5H USER)                                   
      WRITE(6,70)AA,BB,CC                                             
   70 FORMAT(6H AA = ,F10.3,5X,4HBB= ,F10.3,5X,5HCC = ,F10.3)                                                    

C                                                                       ABQP0631
C     READ IN PARAMETERS FOR FUNCTION AND THEN DEFINE FUNCTION          ABQP0632
C     AT EACH RADIAL POINT                                              ABQP0633
C                                                                       ABQP0634
      DO 100 I=1,IMAX                                                   ABQP0635
      R=DR*DFLOAT(I)                                                    ABQP0636
C                                                                       ABQP0637
C     NEXT LINE GIVES VALUE OF FUNCTION AT RADIUS R                     ABQP0638
C                                                                       ABQP0639
  100 FF(I)=0.0D0                                                       ABQP0640
      
      DO I=1,IMAX
      R=DR*DFLOAT(I)                                                    
	A1=3.0D00/(8.00*PI**(3.0D00/2.0D00))
	FF(I)=A1*((1.0D00/AA**3)*DEXP(-R**2/(4.0D00*AA**2))-
     1      CC**2*(6.0D00*BB**2-R**2)/(4.0D00*BB**7)*
     2      DEXP(-R**2/(4.0D00*BB**2)))                                 !6Li
      END DO              
	                                                        
      RETURN                                                            ABQP0641
      END                                                               ABQP0642
      SUBROUTINE DELTA (L)                                              ABQP0643
C                                                                       ABQP0644
C     DELTA FUNCTION                                                    ABQP0645
C                                                                       ABQP0646
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0647
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0648
      COMMON/FFF/FF(1000)                                               ABQP0649
      COMMON/DEL/V                                                      ABQP0650
      WRITE(6,30)                                                       ABQP0651
   30 FORMAT(15H DELTA FUNCTION)                                        ABQP0652
      IF(L.NE.0) WRITE(6,40)                                            ABQP0653
   40 FORMAT(44H NOTE FOURIER TRANSFORM UNDEFINED FOR L.NE.0)           ABQP0654
      READ(5,50)V                                                       ABQP0655
   50 FORMAT(F10.0)                                                     ABQP0656
      IF(V.EQ.0.0D0) V=1.0D0                                            ABQP0657
      WRITE(6,60)V                                                      ABQP0658
   60 FORMAT(5H V = ,F10.3)                                             ABQP0659
      DO 100 I=1,IMAX                                                   ABQP0660
  100 FF(I)=0.0D0                                                       ABQP0661
      RETURN                                                            ABQP0662
      END                                                               ABQP0663
      SUBROUTINE TASSIE (L,VOL,VOL2)                                    ABQP0664
C                                                                       ABQP0665
C     CALCULATES TASSIE TYPE TRANSITION DENSITIES                       ABQP0666
C                                                                       ABQP0667
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0668
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0669
      COMMON/FFF/FF(1000)                                               ABQP0670
      COMMON/TEMPA/F(1000),ARRAY2(1000)                                 ABQP0671
      READ(5,10)ITYPE                                                   ABQP0672
   10 FORMAT(I10)                                                       ABQP0673
      IF(ITYPE.EQ.0.AND.L.EQ.0) WRITE(6,20)                             ABQP0674
   20 FORMAT(13H ERROR -- L=0)                                          ABQP0675
      WRITE(6,30)ITYPE                                                  ABQP0676
   30 FORMAT(34H TASSIE TRANSITION DENSITY ITYPE =,I1)                  ABQP0677
C                                                                       ABQP0678
C     CALCULATE DERIVATIVE OF FUNCTION IN F(I) AND STORE IN FF(I)       ABQP0679
C                                                                       ABQP0680
      HH=1.0D0/(12.D0*DR)                                               ABQP0681
      YY=F(IMAX-4)                                                      ABQP0682
      B=HH*(-25.D0*F(1)+48.D0*F(2)-36.D0*F(3)+16.D0*F(4)-3.D0*F(5))     ABQP0683
      C=HH*(-3.D0*F(1)-10.D0*F(2)+18.D0*F(3)-6.D0*F(4)+F(5))            ABQP0684
      DO 40 I=5,IMAX                                                    ABQP0685
      A=B                                                               ABQP0686
      B=C                                                               ABQP0687
      C=HH*(F(I-4)-F(I)+8.D0*(F(I-1)-F(I-3)))                           ABQP0688
   40 FF(I-4)=A                                                         ABQP0689
      A=HH*(-YY+6.D0*F(IMAX-3)-18.D0*F(IMAX-2)+10.D0*F(IMAX-1)          ABQP0690
     *  +3.D0*F(IMAX))                                                  ABQP0691
      FF(IMAX)=HH*(3.D0*YY-16.D0*F(IMAX-3)+36.D0*F(IMAX-2)              ABQP0692
     *  -48.D0*F(IMAX-1)+25.D0*F(IMAX))                                 ABQP0693
      FF(IMAX-1)=A                                                      ABQP0694
      FF(IMAX-2)=C                                                      ABQP0695
      FF(IMAX-3)=B                                                      ABQP0696
      DO 50 I=1,IMAX                                                    ABQP0697
   50 F(I)=0.0D0                                                        ABQP0698
      VOL=0.0D0                                                         ABQP0699
      VOL2=0.0D0                                                        ABQP0700
      IF(ITYPE.EQ.1.OR.L.LE.1) RETURN                                   ABQP0701
      EXPL=DFLOAT(L-1)                                                  ABQP0702
      DO 60 I=1,IMAX                                                    ABQP0703
      R=DR*DFLOAT(I)                                                    ABQP0704
   60 FF(I)=R**EXPL*FF(I)                                               ABQP0705
      RETURN                                                            ABQP0706
      END                                                               ABQP0707

      SUBROUTINE HARMONIC                                                    
C                                                                        
C     HARMONIC OSCILLATION FUNCTION                                              
C                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                          
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                               
      COMMON/FFF/FF(1000)                                                
      WRITE(6,30)                                                       
   30 FORMAT(23H USER1 DEFINED FUNCTION)                                  

COMMENT   AA=ALPHA
COMMENT   BB=A

      READ(5,50)AA,BB                                         
   50 FORMAT(6F10.0)                                                 
      WRITE(6,60)                                                        
   60 FORMAT(6H USER1)                                   
      WRITE(6,70)AA,BB                                             
   70 FORMAT(5HAA = ,F10.3,5X,5HBB = ,F10.3)                                                    

C                                                                        
C     READ IN PARAMETERS FOR FUNCTION AND THEN DEFINE FUNCTION           
C     AT EACH RADIAL POINT                                               
C                                                                       
      DO 100 I=1,IMAX                                                    
      R=DR*DFLOAT(I)                                                     
C                                                                       
C     NEXT LINE GIVES VALUE OF FUNCTION AT RADIUS R                      
C                                                                        
  100 FF(I)=0.0D0                                                       
      
      DO I=1,IMAX
      R=DR*DFLOAT(I)                                                    
	FF(I)=(1.0D00+AA*(R/BB)**2)*DEXP(-(R/BB)**2)                 ! 7Li CASE
      END DO              
	                                                        
      RETURN                                                            
      END                                                                


      SUBROUTINE USER1                                                    
C                                                                        
C     USER DEFINED FUNCTION                                              
C                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                          
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                               
      COMMON/FFF/FF(1000)                                                
      WRITE(6,30)                                                       
   30 FORMAT(23H USER1 DEFINED FUNCTION)                                  

      READ(5,50)AA,BB,CC,DD,EE,F1                                         
   50 FORMAT(6F10.0)                                                 
      WRITE(6,60)                                                        
   60 FORMAT(6H USER1)                                   
      WRITE(6,70)AA,BB,CC,DD,EE,F1                                             
   70 FORMAT(5HAA = ,F10.3,5X,5HBB = ,F10.3,5X,5HCC = ,F10.3,/
     1      ,5HDD = ,F10.3,5X,5HEE = ,F10.3,5X,5HF1 = ,F10.3)                                                    

C                                                                        
C     READ IN PARAMETERS FOR FUNCTION AND THEN DEFINE FUNCTION           
C     AT EACH RADIAL POINT                                               
C                                                                       
      DO 100 I=1,IMAX                                                    
      R=DR*DFLOAT(I)                                                     
C                                                                       
C     NEXT LINE GIVES VALUE OF FUNCTION AT RADIUS R                      
C                                                                        
  100 FF(I)=0.0D0                                                       
      
      DO I=1,IMAX
      R=DR*DFLOAT(I)                                                    
	FF(I)=(AA+BB*CC**2*R**2)*DEXP(-CC**2*R**2)
     1     +(DD+EE*F1**2*R**2)*DEXP(-F1**2*R**2)                     ! 9Be CASE
      END DO              
	                                                        
      RETURN                                                            
      END                                                                

      SUBROUTINE USER2                                                    
C                                                                        
C     USER DEFINED FUNCTION   (FOR  7LI)                                           
C                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                          
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                               
      COMMON/FFF/FF(1000)                                                
      WRITE(6,30)                                                       
   30 FORMAT(23H USER2 DEFINED FUNCTION)                                  

      READ(5,50)AA,BB,CC                                         
   50 FORMAT(6F10.0)                                                 
      WRITE(6,60)                                                        
   60 FORMAT(6H USER2)                                   
      WRITE(6,70)AA,BB,CC                                             
   70 FORMAT(5HAA = ,F10.5,5X,5HBB = ,F10.5,5X,5HCC = ,F10.5)                                                    

C                                                                        
C     READ IN PARAMETERS FOR FUNCTION AND THEN DEFINE FUNCTION           
C     AT EACH RADIAL POINT                                               
C                                                                       
      DO 100 I=1,IMAX                                                    
      R=DR*DFLOAT(I)                                                     
C                                                                       
C     NEXT LINE GIVES VALUE OF FUNCTION AT RADIUS R                      
C                                                                        
  100 FF(I)=0.0D0                                                       
      
      DO I=1,IMAX
      R=DR*DFLOAT(I)                                                    
	FF(I)=(AA+BB*R**2)*DEXP(-CC**2*R**2)                           ! 7Li CASE
      END DO              
	                                                        
      RETURN                                                            
      END                                                                

      SUBROUTINE USER3                                                    
C                                                                        
C     USER DEFINED FUNCTION   (FOR  8LI)                                           
C                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                          
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                               
      COMMON/FFF/FF(1000)                                                
      WRITE(6,30)                                                       
   30 FORMAT(23H USER3 DEFINED FUNCTION)                                  

      READ(5,50)AA,BB,CC,DD                                         
   50 FORMAT(6F10.0)                                                 
      WRITE(6,60)                                                        
   60 FORMAT(6H USER3)                                   
      WRITE(6,70)AA,BB,CC,DD                                             
   70 FORMAT(5HAA = ,F10.5,5X,5HBB = ,F10.5,5X,5HCC = ,F10.5)                                                    

C                                                                        
C     READ IN PARAMETERS FOR FUNCTION AND THEN DEFINE FUNCTION           
C     AT EACH RADIAL POINT                                               
C                                                                       
      DO 100 I=1,IMAX                                                    
      R=DR*DFLOAT(I)                                                     
C                                                                       
C     NEXT LINE GIVES VALUE OF FUNCTION AT RADIUS R                      
C                                                                        
  100 FF(I)=0.0D0                                                       
      
      DO I=1,IMAX
      R=DR*DFLOAT(I)                                                    
	FF(I)=AA*DEXP(-BB*R**2)+CC*DEXP(-DD*R**2)                ! 8Li CASE
      END DO              
	                                                        
      RETURN                                                            
      END                                                                


      SUBROUTINE FOURBESS                                                    
C                                                                        
C     FOURIER-BESSEL EXPANSION                                              
C                                                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                          
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                               
      COMMON/FFF/FF(1000)                                                
      WRITE(6,30)                                                       
   30 FORMAT(23H USER2 DEFINED FUNCTION)                                  

      READ(5,50) C1,C2,C3,C4,C5,C6,C7,C8,C9,C10
	READ(5,50) C11,C12,C13,C14,C15,C16,C17                                         
   50 FORMAT(10F8.5)                                                 
      WRITE(6,60)                                                        
   60 FORMAT(25H FOURIER-BESSEL EXPANSION)                                   
      WRITE(6,70)C1,C2,C3,C4,C5,C6,C7,C8,C9,C10
      WRITE(6,70)C11,C12,C13,C14,C15,C16,C17                                              
   70 FORMAT(10F10.5)                                                    

C                                                                        
C     READ IN PARAMETERS FOR FUNCTION AND THEN DEFINE FUNCTION           
C     AT EACH RADIAL POINT                                               
C                                                                       
      DO 100 I=1,IMAX                                                    
      R=DR*DFLOAT(I)
C                                                                       
C     NEXT LINE GIVES VALUE OF FUNCTION AT RADIUS R                      
C                                                                        
  100 FF(I)=0.0D0                                                       
      
	RRMAX=DR*IMAX

      DO I=1,IMAX
      R=DR*DFLOAT(I)                                                    
	X1 = 1.0D00*PI/RRMAX*R                                                     
	X2 = 2.0D00*PI/RRMAX*R                                                     
	X3 = 3.0D00*PI/RRMAX*R                                                     
	X4 = 4.0D00*PI/RRMAX*R                                                     
	X5 = 5.0D00*PI/RRMAX*R                                                     
	X6 = 6.0D00*PI/RRMAX*R                                                     
	X7 = 7.0D00*PI/RRMAX*R                                                     
	X8 = 8.0D00*PI/RRMAX*R                                                     
	X9 = 9.0D00*PI/RRMAX*R                                                     
	X10=10.0D00*PI/RRMAX*R                                                     
	X11=11.0D00*PI/RRMAX*R                                                     
	X12=12.0D00*PI/RRMAX*R                                                     

	FF(I)=C1*DSBESJ0(X1)+C2*DSBESJ0(X2)+C3*DSBESJ0(X3)+C4*DSBESJ0(X4)
     1+C5*DSBESJ0(X5)+C6*DSBESJ0(X6)+C7*DSBESJ0(X7)+C8*DSBESJ0(X8)
     2+C9*DSBESJ0(X9)+C10*DSBESJ0(X10)+C11*DSBESJ0(X11)+C12*DSBESJ0(X12)  !144SM CASE
      END DO              
	                                                        
      RETURN                                                            
      END

      SUBROUTINE GAUSSIAN                                               
C                                                                       
C     3 PARAMETER GAUSSIAN FUNCTION                                        
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              
      COMMON/FFF/FF(1000)                                               
      WRITE(6,30)                                                       
   30 FORMAT(27H 3 PARAMETER FERMI FUNCTION)                            
      READ(5,50)C,RR,AR,W                                               
   50 FORMAT(4F10.0)                                                    
      IF(C.EQ.0.0D0) C=1.0D0                                            
      WRITE(6,60)C,RR,AR,W                                              
   60 FORMAT(5H C = ,F10.3,5X,5HRR = ,F10.3,5X,5HAR = ,F10.3,5X,        
     *4HW = ,F10.3)                                                     

      DO 100 I=1,IMAX                                                   
      R=DR*DFLOAT(I)                                                    
      FF(I)=C*(1.0D0+W*R*R/RR**2)/(1.0D0+DEXP((R**2-RR**2)/AR**2))      
  100	CONTINUE

      RETURN                                                            
      END                                                               

      SUBROUTINE TRANSF (L,ITYPE)                                       ABQP0708
C                                                                       ABQP0709
C     L IS THE MULTIPOLARITY                                            ABQP0710
C                                                                       ABQP0711
C     ITYPE = 0    CALCULATES FOURIER TRANSFORM                         ABQP0712
C     ITYPE = 1    CALCULATES INVERSE TRANSFORM                         ABQP0713
C                                                                       ABQP0714
C     BASED ON FILON'S FORMULA IN M ABRAMOWITZ AND I A STEGUN (1965)    ABQP0715
C     HANDBOOK OF MATHEMATICAL FUNCTIONS, DOVER, NEW YORK.              ABQP0716
C                                                                       ABQP0717
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0718
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0719
      COMMON/TEMPA/AIN(1000),AOUT(1000)                                 ABQP0720
C                                                                       ABQP0721
C     AIN   CONTAINS THE FUNCTION                                       ABQP0722
C     AOUT  CONTAINS THE FOURIER TRANSFORM                              ABQP0723
C                                                                       ABQP0724
      REAL*8 K,KH,KH2,KH3                                               ABQP0725
      PIPI=2.0D0*PI                                                     ABQP0726
      IF(ITYPE.EQ.1) GOTO 5                                             ABQP0727
C                                                                       ABQP0728
C     SETS INTERVALS AND LIMITS FOR FOURIER TRANSFORM                   ABQP0729
C                                                                       ABQP0730
      HA=DR                                                             ABQP0731
      HB=DQ                                                             ABQP0732
      NA=IMAX                                                           ABQP0733
      NB=KMAX                                                           ABQP0734
      CONST=4.0D0*PI                                                    ABQP0735
      GOTO 15                                                           ABQP0736
C                                                                       ABQP0737
C     SETS INTERVALS AND LIMITS FOR INVERSE TRANSFORM                   ABQP0738
C                                                                       ABQP0739
    5 HA=DQ                                                             ABQP0740
      HB=DR                                                             ABQP0741
      NA=KMAX                                                           ABQP0742
      NB=IMAX                                                           ABQP0743
      CONST=1.0D0/(2.0D0*PI**2)                                         ABQP0744
   15 IF(L.GT.0) GOTO 110                                               ABQP0745
C                                                                       ABQP0746
C     INTEGRATION FOR L.EQ.0                                            ABQP0747
C                                                                       ABQP0748
      DO 100 J=1,NB                                                     ABQP0749
      K=HB*DFLOAT(J)                                                    ABQP0750
      KH=K*HA                                                           ABQP0751
      KH2=KH*KH                                                         ABQP0752
      KH3=KH2*KH                                                        ABQP0753
      SINKH=DSIN(DMOD(KH,PIPI))                                         ABQP0754
      COSKH=DCOS(DMOD(KH,PIPI))                                         ABQP0755
      A=1.0D0/KH+SINKH*COSKH/KH2-2.0D0*SINKH*SINKH/KH3                  ABQP0756
      B=(1.0D0+COSKH*COSKH)/KH2-2.0D0*SINKH*COSKH/KH3                   ABQP0757
      C=4.0D0*(SINKH/KH3-COSKH/KH2)                                     ABQP0758
      EVEN=0.0D0                                                        ABQP0759
      ODD=0.0                                                           ABQP0760
      DO 20 I=2,NA,2                                                    ABQP0761
      F=DFLOAT(I)                                                       ABQP0762
      X=F*KH                                                            ABQP0763
      R=F*HA                                                            ABQP0764
   20 EVEN=EVEN+AIN(I)*R*DSIN(DMOD(X,PIPI))                             ABQP0765
      DO 25 I=1,NA,2                                                    ABQP0766
      F=DFLOAT(I)                                                       ABQP0767
      X=F*KH                                                            ABQP0768
      R=F*HA                                                            ABQP0769
   25 ODD=ODD+AIN(I)*R*DSIN(DMOD(X,PIPI))                               ABQP0770
      R=HA*DFLOAT(NA)                                                   ABQP0771
      X=R*K                                                             ABQP0772
      AA=AIN(NA)*R                                                      ABQP0773
      S=-AA*DSIN(DMOD(X,PIPI))+2.0D0*EVEN                               ABQP0774
  100 AOUT(J)=CONST*HA*(-A*AA*DCOS(DMOD(X,PIPI))+B*S+C*ODD)/K           ABQP0775
      RETURN                                                            ABQP0776
C                                                                       ABQP0777
C     INTEGRATION FOR L.GT.0                                            ABQP0778
C                                                                       ABQP0779
  110 DO 200 J=1,NB                                                     ABQP0780
      K=HB*DFLOAT(J)                                                    ABQP0781
      KH=K*HA                                                           ABQP0782
      KH2=KH*KH                                                         ABQP0783
      KH3=KH2*KH                                                        ABQP0784
      SINKH=DSIN(DMOD(KH,PIPI))                                         ABQP0785
      COSKH=DCOS(DMOD(KH,PIPI))                                         ABQP0786
      A=1.0D0/KH+SINKH*COSKH/KH2-2.0D0*SINKH*SINKH/KH3                  ABQP0787
      B=(1.0D0+COSKH*COSKH)/KH2-2.0D0*SINKH*COSKH/KH3                   ABQP0788
      C=4.0D0*(SINKH/KH3-COSKH/KH2)                                     ABQP0789
      EVEN=0.0D0                                                        ABQP0790
      ODD=0.0D0                                                         ABQP0791
      DO 240 I=2,NA,2                                                   ABQP0792
      F=DFLOAT(I)                                                       ABQP0793
      X=F*KH                                                            ABQP0794
      R=F*HA                                                            ABQP0795
      AA=AIN(I)*R*R                                                     ABQP0796
      U=AA*SBES(L,X)                                                    ABQP0797
      V=AA*CBES(L,X)                                                    ABQP0798
  240 EVEN=EVEN+U*DSIN(DMOD(X,PIPI))+V*DCOS(DMOD(X,PIPI))               ABQP0799
      DO 250 I=1,NA,2                                                   ABQP0800
      F=DFLOAT(I)                                                       ABQP0801
      X=F*KH                                                            ABQP0802
      R=F*HA                                                            ABQP0803
      AA=AIN(I)*R*R                                                     ABQP0804
      U=AA*SBES(L,X)                                                    ABQP0805
      V=AA*CBES(L,X)                                                    ABQP0806
  250 ODD=ODD+U*DSIN(DMOD(X,PIPI))+V*DCOS(DMOD(X,PIPI))                 ABQP0807
      R=HA*DFLOAT(NA)                                                   ABQP0808
      X=R*K                                                             ABQP0809
      AA=AIN(NA)*R*R                                                    ABQP0810
      UA=AA*SBES(L,X)                                                   ABQP0811
      VA=AA*CBES(L,X)                                                   ABQP0812
      S=-UA*DSIN(DMOD(X,PIPI))-VA*DCOS(DMOD(X,PIPI))+2.0D0*EVEN         ABQP0813
  200 AOUT(J)=CONST*HA*(-A*UA*DCOS(DMOD(X,PIPI))+A*VA*DSIN(DMOD(X,PIPI))ABQP0814
     *       +B*S+C*ODD)                                                ABQP0815
      RETURN                                                            ABQP0816
      END                                                               ABQP0817
      DOUBLE PRECISION FUNCTION SBES(L,X)                               ABQP0818
C                                                                       ABQP0819
C     USED FOR DECOMPOSITION OF SPHERICAL BESSEL FUNCTION FOR L.GT.0    ABQP0820
C                                                                       ABQP0821
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0822
      A=1.0D0/X                                                         ABQP0823
      B=A/X                                                             ABQP0824
      SBES=B                                                            ABQP0825
      IF(L.EQ.1) RETURN                                                 ABQP0826
      DO 10 N=2,L                                                       ABQP0827
      C=(2.0D0*DFLOAT(N)-1.0D0)*B/X-A                                   ABQP0828
      A=B                                                               ABQP0829
   10 B=C                                                               ABQP0830
      SBES=C                                                            ABQP0831
      RETURN                                                            ABQP0832
      END                                                               ABQP0833
      DOUBLE PRECISION FUNCTION CBES(L,X)                               ABQP0834
C                                                                       ABQP0835
C     USED FOR DECOMPOSITION OF SPHERICAL BESSEL FUNCTION FOR L.GT.0    ABQP0836
C                                                                       ABQP0837
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0838
      B=1.0D0/X                                                         ABQP0839
      A=B/X                                                             ABQP0840
      LPLUS1=L+1                                                        ABQP0841
      DO 10 N=1,LPLUS1                                                  ABQP0842
      C=(3.0D0-2.0D0*DFLOAT(N))*B/X-A                                   ABQP0843
      A=B                                                               ABQP0844
   10 B=C                                                               ABQP0845
      CBES=C*(-1.0D0)**LPLUS1                                           ABQP0846
      RETURN                                                            ABQP0847
      END                                                               ABQP0848
      SUBROUTINE INTEGR (L,X,VOL,VOL2)                                  ABQP0849
C                                                                       ABQP0850
C     CALCULATES R**2 AND R**4 WEIGHTED VOLUME INTEGRALS                ABQP0851
C                                                                       ABQP0852
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0853
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0854
      DIMENSION X(1000)                                                 ABQP0855
      SUM1=0.0D0                                                        ABQP0856
      SUM2=0.0D0                                                        ABQP0857
      EL2=DFLOAT(L+2)                                                   ABQP0858
      EL4=DFLOAT(L+4)                                                   ABQP0859
      JMAX=IMAX-MOD((IMAX+1),2)                                         ABQP0860
      DO 100 J=1,JMAX                                                   ABQP0861
      R=DR*DFLOAT(J)                                                    ABQP0862
      R2=R**EL2                                                         ABQP0863
      R4=R**EL4                                                         ABQP0864
      FACTOR=2.0D0+2.0D0*DFLOAT(MOD(J,2))                               ABQP0865
      IF(J.EQ.JMAX) FACTOR=1.0D0                                        ABQP0866
      SUM1=SUM1+FACTOR*R2*X(J)                                          ABQP0867
  100 SUM2=SUM2+FACTOR*R4*X(J)                                          ABQP0868
      VOL=4.0D0*PI*SUM1*DR/3.0D0                                        ABQP0869
      VOL2=4.0D0*PI*SUM2*DR/3.0D0                                       ABQP0870
      RETURN                                                            ABQP0871
      END                                                               ABQP0872
      DOUBLE PRECISION FUNCTION CLEB(L1,L2,L)                           ABQP0873
C                                                                       ABQP0874
C     CLEBSCH-GORDAN COEFFICIENT FOR COUPLING OF L1 AND L2 TO           ABQP0875
C     GIVE RESULTANT L WITH ALL Z COMPONENTS EQUAL TO ZERO              ABQP0876
C                                                                       ABQP0877
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0878
      DIMENSION FACTL(62)                                               ABQP0879
C                                                                       ABQP0880
C     FACTL(I+1) FOR FACTORIAL I                                        ABQP0881
C                                                                       ABQP0882
      MAXL=L1+L2+L                                                      ABQP0883
      CLEB=1.0D0                                                        ABQP0884
      IF(MAXL.EQ.0) RETURN                                              ABQP0885
      CLEB=0.0D0                                                        ABQP0886
      IF((MAXL/2)*2.NE.MAXL) RETURN                                     ABQP0887
      MAXL2=MAXL+2                                                      ABQP0888
      IF(MAXL2.GT.62) RETURN                                            ABQP0889
      FACTL(1)=1.0D0                                                    ABQP0890
      DO 10 I=2,MAXL2                                                   ABQP0891
   10 FACTL(I)=DFLOAT(I-1)*FACTL(I-1)                                   ABQP0892
      LG=MAXL/2                                                         ABQP0893
      DEL=DSQRT(FACTL(L1+L2-L+1)*FACTL(L1+L-L2+1)*FACTL(L2+L-L1+1)/     ABQP0894
     *    FACTL(L1+L2+L+2))                                             ABQP0895
      CLEB=(-1.0D0)**(LG+L)*DSQRT(DFLOAT(2*L+1))*DEL*FACTL(LG+1)/       ABQP0896
     *    (FACTL(LG-L1+1)*FACTL(LG-L2+1)*FACTL(LG-L+1))                 ABQP0897
      RETURN                                                            ABQP0898
      END                                                               ABQP0899
      DOUBLE PRECISION FUNCTION FACT2(L)                                ABQP0900
C                                                                       ABQP0901
C     CALCULATES DOUBLE FACTORIAL OF ARGUMENT 2*L + 1                   ABQP0902
C                                                                       ABQP0903
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0904
      FACT2=1.0D0                                                       ABQP0905
      IF(L.EQ.0) RETURN                                                 ABQP0906
      X=1.0D0                                                           ABQP0907
      DO 10 I=1,L                                                       ABQP0908
      X=X+2.0D0                                                         ABQP0909
   10 FACT2=FACT2*X                                                     ABQP0910
      RETURN                                                            ABQP0911
      END                                                               ABQP0912

      SUBROUTINE PRINT                                                  ABQP0913
C                                                                       ABQP0914
C     PRINTS DENSITIES, INTERACTION AND POTENTIAL                       ABQP0915
C                                                                       ABQP0916
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0917
      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0918
      COMMON/TEMPA/X(1000),ARRAY2(1000)                                 ABQP0919
      N8=IMAX/4                                                         ABQP0920
      N9=N8+N8                                                          ABQP0921
      N10=N8+N9                                                         ABQP0922
      DO 10 I=1,N8                                                      ABQP0923
      I1=I+N8                                                           ABQP0924
      I2=I+N9                                                           ABQP0925
      I3=I+N10                                                          ABQP0926
      R=DR*DFLOAT(I)                                                    ABQP0927
      R1=DR*DFLOAT(I1)                                                  ABQP0928
      R2=DR*DFLOAT(I2)                                                  ABQP0929
      R3=DR*DFLOAT(I3)                                                  ABQP0930
   10 WRITE(6,20)R,X(I),R1,X(I1),R2,X(I2),R3,X(I3)                      ABQP0931
   20 FORMAT(3X,4(0PF9.2,1PE19.8))                                      ABQP0932
      N9=4*N8                                                           ABQP0933
      N10=IMAX-N9                                                       ABQP0934

	DO N=1,400
	R=DR*DFLOAT(N)
	WRITE(10,*) R, X(N)
	END DO

      IF(N10.EQ.0)RETURN                                                ABQP0935
      DO 30 I=1,N10                                                     ABQP0936
      N8=N9+I                                                           ABQP0937
      R=DR*DFLOAT(N8)                                                   ABQP0938
   30 WRITE(6,40)R,X(N8)                                                ABQP0939
   40 FORMAT(87X,0PF9.2,1PE19.8)                                        ABQP0940
      RETURN                                                            ABQP0941

      END                                                               ABQP0942
CX    //L.SYSIN DD *                                                          ABQP0943
CX    //G.FT08F001  DD UNIT=WORK,DISP=(NEW,DELETE),SPACE=(TRK,2),             ABQP0944
CX    //   DCB=(RECFM=VBS,BLKSIZE=1000),DSN=&&POT                             ABQP0945
CX    //G.SYSIN DD *                                                          ABQP0946
