      PROGRAM FDU0
C
C========================================== fitting  version  ========C
C                                                                     C
C    FDU0 Package                                                     C
C    ---- -------                                                     C
C                                                                     C
C    Shell-Model code in S and D fermion pairs truncated model space  C
C    based on the Fermion Dynamical Symmetry prescription             C
C                                                                     C
C    Codes : FDU0     - Shell-Model code                              C
C            FDSMCFP  - Pair CFP code                                 C
C            PD       - Pair and Multipole Operators Matrix Elements  C
C            FDTR     - Transitions and Moments                       C
C            LIB      - Utility programs                              C
C                                                                     C
C    Ref.  : Hua Wu and M. Vallieres                                  C
C            Phys. Rev. C39 (1989) 1066-1075                          C
C                                                                     C
C    Code History:                                                    C
C                                                                     C
C    -  Original code writen by Hua Wu, Sept. 1988                    C
C                                                                     C
C    -  Checks:                                                       C
C       sp6xsp6 program(su2,su3 limits)                               C
C       single so8 three limits                                       C
C       couple so8 so6 limit            9/20/88                       C
C                                                                     C
C    -  Fitting Procedure Added         Jan. 1989                     C
C                                                                     C
C    Computers:                                                       C
C                                                                     C
C    Works as is on VAXes and CONVEX                                  C
C    Works on IBM mainframe after changing the file names             C
C              in all open statements to: '/ fname dat'               C
C                                                                     C
C=====================================================================C
C       H =                                                           C
C  EN*JN(JN+1) + B2N*P2N.P2N + B3N*P3N.P3N +                          C
C  G0N*SN(+).SN + G2N*DN(+).DN  +                                     C
C  EP*JP(JP+1) + B2P*P2P.P2P + B3P*P3P.P3P +                          C
C  G0P*SP(+).SP + G2P*DP(+).DP +                                      C
C  ENP*JN.JP   + B2NP*P2N.P2P + B3NP*P3N.P3P                          C
C  B1N*P1N.P1N + B1P*P1P.P1P  + B1NP*P1N.P1P                          C
C  EJ*J(J+1)                                                          C
C=====================================================================C
C
      PARAMETER (MAXP=11221,MAXMULT=1752,MAXPAIR=1168,MAXN=21)
      PARAMETER (MAXNS=4200,MAXNE=400000,MAXLAN=100,MAXDNS=220)
      PARAMETER (MAXEE=50)
      NAMELIST /INP/
     &         NEUTRON_SYMMETRY,PROTON_SYMMETRY,  !symmetry information
     &         OMN,OMP,N1P,N1N,                   !shell information
     &         EP,B2P,B3P,G0P,G2P,                !hamiltonian
     &         EN,B2N,B3N,G0N,G2N,                !information
     &         ENP,B2NP,B3NP,
     &         B1P,B1N,B1NP,EJ,                   !auxiliary input param
     &         WAVE,NJ,LAN,MAXJ,minj,             !control infromation
     &         FIT,                               !fitting control
     &         EPS,                               !small number for
     &                                            !escape from fitting
     &         STEP,                              !fitting parameters
     &                                            !changing speed
     &         FIT_MODE,                          !first fit_mode,'L'
     &                                            ! or 'G'
     &         LEVELS,                            !input fitting levels i.d.
     &         EE,                                !experimental values
     &         WEIGHT                             !fitting weights
C
      CHARACTER NEUTRON_SYMMETRY*3,PROTON_SYMMETRY*3,WAVE*3,FIT_MODE*1
      INTEGER   OMN,OMP,PAIRENDN,PAIRENDP,PENDN,PENDP,MULTENDP,MULTENDN
      INTEGER   JSIZEN(0:2*MAXN),JSIZEP(0:2*MAXN)
      INTEGER   PINDEXN(MAXP),PINDEXP(MAXP)
      INTEGER   PAIRINDEXN(MAXPAIR),PAIRINDEXP(MAXPAIR)
      INTEGER   MULTINDEXN(MAXMULT),MULTINDEXP(MAXMULT)
      INTEGER   STATEJN(MAXNS),STATEAN(MAXNS)
      INTEGER   STATEJP(MAXNS),STATEAP(MAXNS)
      INTEGER   AP_LEFT,AP_RIGHT,AN_LEFT,AN_RIGHT,AN,AP,R
      REAL*8    PTABN(MAXP),PTABP(MAXP),DTERM,SIXJ,DSQRT,PHASE
      REAL      PAIRTABN(MAXPAIR),PAIRTABP(MAXPAIR)
      REAL      MULTTABN(MAXMULT),MULTTABP(MAXMULT)
      REAL      E(MAXNS),Z(MAXNS*MAXLAN)
      REAL      X(13),X0(13),B(13),A(13,13),XX(13)
      INTEGER   ICHNG(2*13)
      REAL      FITTING(MAXEE,13),GROUND_EXPECT(13),EXPECT(13),STEP
      EQUIVALENCE (X(1),EN),(X(2),B2N),(X(3),B3N),(X(4),G0N),(X(5),G2N)
      EQUIVALENCE (X(6),EP),(X(7),B2P),(X(8),B3P),(X(9),G0P),(X(10),G2P)
      EQUIVALENCE (X(11),ENP),(X(12),B2NP),(X(13),B3NP)
      REAL      LEVELS(MAXEE),EE(MAXEE),WEIGHT(MAXEE),ETH(MAXEE),PFI(13)
      INTEGER   ADJUST(13),ORDER,FIT
C
      COMMON /TAB/ PTABN,PTABP,PAIRTABN,PAIRTABP,MULTTABN,MULTTABP
      COMMON /TABEND/ PAIRENDN,PAIRENDP,PENDN,PENDP,MULTENDP,MULTENDN
      COMMON /TABINDEX/ PINDEXN,PINDEXP,PAIRINDEXN,PAIRINDEXP,
     &                  MULTINDEXN,MULTINDEXP
      COMMON /STATE/ STATEJN,STATEAN,STATEJP,STATEAP
      COMMON /OTHER/ JSIZEN,JSIZEP,N1N,N1P,LAN
C
      INDP(IAP,JP,R,IA,J)=20000000*J+200000*JP+10000*R+100*IA+IAP
      INDPAIR(IAP,J,R,IA)=1000000*R+10000*J+100*IAP+IA
      PHASE(I)=1-2*MOD(ABS(I),2)
      LTRIANGLE(J1,J2,J3)=(ABS(J1-J2).LE.J3).AND.(J1+J2.GE.J3)
C---------------------------------------------
C     CALL XUFLOW(0)       !  to remove the underflow for IBM mainframe
C                          !  can be replaced by noxuflow run time option
C-------------------------------------------
C                                                 set defaults
      NEUTRON_SYMMETRY='SP6'
      PROTON_SYMMETRY='SO8'
      OMP=10
      OMN=15
      N1P=0
      N1N=0
      EP=0.
      B2P=0.
      B3P=0.
      G0P=0.
      G2P=0.
      EN=0.
      B2N=0.
      B3N=0.
      G0N=0.
      G2N=0.
      B2NP=0.
      B3NP=0.
      ENP=0.
      WAVE='OFF'
      MINJ=0.
      MAXJ=10
      NJ=8
      LAN=100
      B1N=0.
      B1P=0.
      JOBN=0
      B1NP=0.
      FIT=0
      STEP=-1000
      FIT_MODE='L'
      DO I=1,MAXEE
       EE(I)=0.
       LEVELS(I)=0.
       WEIGHT(I)=1.
      END DO
      OPEN(1,FILE='FDSM.INP',STATUS='OLD')
      READ(1,INP)
      REWIND(1)
      IF(LAN.GT.MAXLAN) CALL DIMERR('maxlan',LAN )
      IF(PROTON_SYMMETRY.EQ.'sp6') PROTON_SYMMETRY='SP6'
      IF(PROTON_SYMMETRY.EQ.'so8') PROTON_SYMMETRY='SO8'
      IF(NEUTRON_SYMMETRY.EQ.'sp6') NEUTRON_SYMMETRY='SP6'
      IF(NEUTRON_SYMMETRY.EQ.'so8') NEUTRON_SYMMETRY='SO8'
      IF(WAVE.EQ.'on') WAVE='ON'
      IF((N1N.GT.MAXN).OR.(N1P.GT.MAXN).OR.
     &   (OMN.GT.MAXN).OR.(OMP.GT.MAXN))
     &   STOP 'N1 OR OMI TOO LARGE'
      IF((PROTON_SYMMETRY.NE.'SP6').AND.(PROTON_SYMMETRY.NE.'SO8'))
     &  STOP 'PROTON SYMMETRY TYPE CAN ONLY BE SP6 OR SO8'
      IF((NEUTRON_SYMMETRY.NE.'SP6').AND.(NEUTRON_SYMMETRY.NE.'SO8'))
     &  STOP 'PROTON SYMMETRY TYPE CAN ONLY BE SP6 OR SO8'
C                                                 handle auxiliary parameters
      IF(PROTON_SYMMETRY.EQ.'SP6') THEN
       EP=EP+3/8.0*B1P
       FBP=SQRT(3/8.0)
      ELSE
       EP=EP+0.2*B1P
       FBP=SQRT(0.2)
      END IF
      IF(NEUTRON_SYMMETRY.EQ.'SP6') THEN
       EN=EN+3/8.0*B1N
       FBN=SQRT(3/8.0)
      ELSE
       EN=EN+0.2*B1N
       FBN=SQRT(0.2)
      END IF
      ENP=ENP+FBN*FBP*B1NP
      EN=EN+EJ
      EP=EP+EJ
      ENP=ENP+2*EJ
      B1P=0.
      B1N=0.
      B1NP=0.
      EJ=0.
      IF((PROTON_SYMMETRY.EQ.'SP6').AND.(B3P.NE.0.)) THEN
       WRITE(2,*) 'B3P can only be zero when proton with sp6 symmetry'
       B3P=0.
      END IF
      IF((NEUTRON_SYMMETRY.EQ.'SP6').AND.(B3N.NE.0.)) THEN
       WRITE(2,*) 'B3N can only be zero when neutron with sp6 symmetry'
       B3N=0.
      END IF
      IF(((PROTON_SYMMETRY.NE.'SO8').OR.(NEUTRON_SYMMETRY.NE.'SO8'))
     &   .AND.(B3NP.NE.0.)) THEN
       WRITE(2,*)
     & 'B3NP can only be zero when proton or neutron with sp6 symmetry'
       B3NP=0.
      END IF
      IF(N1N.EQ.0) THEN
       ENP=0.
       B2NP=0.
       B3NP=0.
       EN=0.
       B2N=0.
       B3N=0.
       G0N=0.
       G2N=0.
      END IF
      IF(N1P.EQ.0) THEN
       ENP=0.
       B2NP=0.
       B3NP=0.
       EP=0.
       B2P=0.
       B3P=0.
       G0P=0.
       G2P=0.
      END IF
      SMALL=ABS(EP)+ABS(B2P)+ABS(B3P)+ABS(G0P)+ABS(G2P)+
     &   ABS(EN)+ABS(B2N)+ABS(B3N)+ABS(G0N)+ABS(G2N)+
     &   ABS(ENP)+ABS(B2NP)+ABS(B3NP)
      SMALL=SMALL*1E-6
      NEE=0
      DO WHILE(LEVELS(NEE+1).NE.0.)
       NEE=NEE+1
      END DO
      IF(FIT.GT.0) THEN
       JOBN=1
       DO I=1,NEE
        DO II=1,13
         FITTING(I,II)=0.
        END DO
       END DO
       DO I=1,13
        IF(X(I).NE.0.) THEN
         ADJUST(I)=1
        ELSE
         ADJUST(I)=0
        END IF
       END DO
      END IF                    !fit
      IF(WAVE.EQ.'ON') JOBN=1
      OPEN(2,FILE='FDSM.OUT',STATUS='NEW')
      WRITE(2,*) '****************INPUT IMAGE**********************'
      CALL COPY(1,2)
      WRITE(2,*) '*************************************************'
      CLOSE(1)
      IF(FIT.EQ.0) THEN
       WRITE(2,*) 'Internal parameters used by program:'
       WRITE(2,INP)
      END IF
      IF(WAVE.EQ.'ON') THEN
       OPEN(14,FILE='FDSM.WAV',FORM='UNFORMATTED',STATUS='NEW')
      END IF
C                                                 read files
      IF(N1N.NE.0) THEN
       CALL GETFILE(NEUTRON_SYMMETRY,'PAIR',OMN,N1N,PAIRENDN,MAXPAIR,
     &                           PAIRINDEXN,PAIRTABN,ANY)
       CALL GETFILE(NEUTRON_SYMMETRY,'P',OMN,N1N,PENDN,MAXP,
     &                             PINDEXN,ANY,PTABN)
       CALL GETFILE(NEUTRON_SYMMETRY,'MULT',OMN,N1N,MULTENDN,MAXMULT,
     &                                MULTINDEXN,MULTTABN,ANY)
       CALL GETFILE(NEUTRON_SYMMETRY,'JSIZE',OMN,N1N,
     &                                ANY,ANY,JSIZEN(0),ANY,ANY)
      ELSE
       JSIZEN(0)=1
      END IF
      IF(N1P.NE.0) THEN
       CALL GETFILE(PROTON_SYMMETRY,'PAIR',OMP,N1P,PAIRENDP,MAXPAIR,
     &                        PAIRINDEXP,PAIRTABP,ANY)
       CALL GETFILE(PROTON_SYMMETRY,'P',OMP,N1P,PENDP,MAXP,
     &                                  PINDEXP,ANY,PTABP)
       CALL GETFILE(PROTON_SYMMETRY,'MULT',OMP,N1P,MULTENDP,MAXMULT,
     &                             MULTINDEXP,MULTTABP,ANY)
       CALL GETFILE(PROTON_SYMMETRY,'JSIZE',OMP,N1P,
     &                       ANY,ANY,JSIZEP(0),ANY,ANY)
      ELSE
       JSIZEP(0)=1
      END IF
C
      FIP=1E20
 888  CONTINUE
      IF(WAVE.EQ.'ON') THEN
       REWIND(14)
       WRITE(14) NEUTRON_SYMMETRY,PROTON_SYMMETRY,OMN,OMP,N1N,N1P,
     &           EN,B2N,B3N,G0N,G2N,
     &           EP,B2P,B3P,G0P,G2P,
     &           ENP,B2NP,B3NP
      END IF
      FI=0.
      EMIN=0.
      DO 100 J=MINJ,MAXJ
       CALL HART(X,J,NJ,JOBN,E,Z,NS)
       IF(NS.LE.0) GOTO 100
       IF(J.EQ.0) EMIN=E(1)
       M=MIN(NJ,NS,LAN)
       DO I=1,M
        E(I)=E(I)-EMIN
       END DO
       PRINT*,'Energies:'
       PRINT '(8F10.4)',(E(I),I=1,M)
       WRITE(2,'(1X,A,I2)') 'J=',J
       WRITE(2,'(8F10.4)') (E(I),I=1,M)
       DO 501 I=1,NEE
        IF(INT(LEVELS(I)).NE.J) GOTO 501
        ORDER=NINT(10*(LEVELS(I)-J))
        ETH(I)=E(ORDER)
        FI=FI+WEIGHT(I)*(ETH(I)-EE(I))**2
501    END DO
       DO 50 I=1,NEE
        IF(INT(LEVELS(I)).NE.J) GOTO 50
        ORDER=NINT(10*(LEVELS(I)-J))
        IF(FIT.GT.0) THEN
         CALL EXPECTATION(J,Z(NS*(ORDER-1)+1),NS,ADJUST,EXPECT)
         IF(LEVELS(I).EQ.0.1) THEN
          DO II=1,13
           GROUND_EXPECT(II)=EXPECT(II)
          END DO
         END IF
         DO II=1,13
          FITTING(I,II)=EXPECT(II)
         END DO
        END IF
  50   END DO
C                                                 output wave functions
       IF(WAVE.EQ.'OFF') GOTO 100
       WRITE(14) J,NS,M
       WRITE(14) (STATEJN(I),STATEAN(I),STATEJP(I),STATEAP(I),I=1,NS)
       WRITE(14) ( E(I),(Z(II+(I-1)*NS),II=1,NS) ,I=1,M)
 100  END DO               !J
      WRITE(2,'(A,F9.3)') ' FIRST 0+ STATE ENERGY=',EMIN
      WRITE(2,*) 'FI=',FI
      PRINT*,'FI=',FI
      IF(FIT.LE.0) GOTO 999
      IF(FI.GE.FIP) THEN
       DO II=1,13
        X(II)=X0(II)
       END DO
       FIT_MODE='G'
c      point one at Õx->å, point two at Õx-> - step/td*grad(fi)å
c      in direction grad(fi), define xl,fi=fi(xl)
c      fip=fi(xl=0)=c
c      fi=fi(xl=-step)=a*xl**2+b*xl+c
c      d(fi)/d(xl)=b=td
c      so fi-fip=a*xl**2+td*xl
c      a=(fi-fip-td*xl)/(xl**2)
c      minimum point at -b/(2a)=
c      -td*xl**2/2.0/(fi-fip-td*xl)=
c      -td*step**2/2.0/(fi-fip+td*step)
c      so new will be
       STEP=STEP**2*TD/2.0/(FI-FIP+TD*STEP)   !parabola fit
       FI=FIP
      ELSE
       DO I=1,NEE
        IF(INT(LEVELS(I)).LE.MAXJ) THEN
         DO II=1,13
          FITTING(I,II)=FITTING(I,II)-GROUND_EXPECT(II)
         END DO
        ELSE
         DO II=1,13
          FITTING(I,II)=0.
         END DO
        END IF
       END DO
       TD=0.
       DO II=1,13
        IF(ADJUST(II).EQ.1) THEN
         SUM=0.
         DO I=1,NEE
          SUM=SUM+WEIGHT(I)*(ETH(I)-EE(I))*FITTING(I,II)
         END DO
         PFI(II)=2*SUM
        ELSE
         PFI(II)=0.
        END IF
        TD=TD+PFI(II)**2
       END DO
       TD=SQRT(TD)
       IF(TD.LT.EPS) THEN
        WRITE(2,*) 'minimum reached with total derivative=',TD
        PRINT*,    'minimum reached with total derivative=',TD
        GOTO 999
       END IF
      END IF                    !fit
 56   IF(FIT_MODE.EQ.'G') THEN
       IF(STEP.LT.0.) STEP=FI/TD
       WRITE(2,*) ' Derivatives used to generator new parameters:'
       WRITE(2,'(2(1X,5E14.4,/),1X,3E14.4)') (PFI(II),II=1,13)
       WRITE(2,*) 'Total derivative=',TD,'    Step=',step
       WRITE(2,*) 'Expected reduction for FI=',step*td
       WRITE(2,*) 'Expected new FI=',fi-Step*td
       PRINT*, ' Derivatives used to generator new parameters:'
       PRINT '(2(1X,5E14.4,/),1X,3E14.4)',(PFI(II),II=1,13)
       PRINT*,'Total derivative=',TD,'    Step=',step
       PRINT*,'Expected reduction for FI=',step*td
       PRINT*,'Expected new FI=',fi-Step*td
       DO II=1,13
        X0(II)=X(II)
        X(II)=X(II)-PFI(II)*STEP/TD
       END DO
       GOTO 777
      END IF
C                                                 linear fitting
      I1=0
      I2=0
      DO WHILE(I1.LT.13)
       I1=I1+1
       IF(ADJUST(I1).EQ.1) THEN
        I2=I2+1
        I3=0
        I4=0
        DO WHILE((I3.LT.13).AND.(I4.LT.I2))
         I3=I3+1
         IF(ADJUST(I3).EQ.1) THEN
          I4=I4+1
          SUM=0.
          DO I=1,NEE
           SUM=SUM+WEIGHT(I)*FITTING(I,I1)*FITTING(I,I3)
          END DO
          A(I2,i4)=SUM
          A(I4,I2)=SUM
         END IF
        END DO
       END IF
      END DO
      I1=0
      I2=0
      DO WHILE(I1.LT.13)
       I1=I1+1
       IF(ADJUST(I1).EQ.1) THEN
        I2=I2+1
        SUM=0.
        DO I=1,NEE
         SUM=SUM+EE(I)*WEIGHT(I)*FITTING(I,I1)
        END DO
        B(I2)=SUM
       END IF
      END DO
      NA=I2
C                                              solve linear equation   A Y = B
C                                              solution put in B
C     CALL LEQ2S(A,NA,B,1,13,0,ICHNG,DET,IER)
      CALL LSARG(NA,A,13,B,1,XX)
      I1=0
      I2=0
      DO WHILE(I1.LT.13)
       I1=I1+1
       IF(ADJUST(I1).EQ.1) THEN
        I2=I2+1
        X0(I1)=X(I1)
        X(I1)=XX(i2)
       END IF
      END DO
 777  CONTINUE
      FIP=FI
      FIT=FIT-1
      WRITE(2,*)
      WRITE(2,*) ' New set of parameters:'
      WRITE(2,    INP)
      PRINT*
      PRINT*
      PRINT*, '  New set of parameters:'
      PRINT INP
      GOTO 888
 999  END

      SUBROUTINE EXPECTATION(J,Z,NS,ADJUST,EXPECT)
C=====================================================================C
C                                                                     C
C EXPECT(1)=<JN(JN+1)>  EXPECT(2)=<P2N.P2N>  EXPECT(3)=<P3N.P3N>      C
C EXPECT(4)=<SN(+).SN>  EXPECT(5)=<DN(+).DN>                          C
C EXPECT(6)=<JP(JP+1)>  EXPECT(7)=<P2P.P2P>  EXPECT(8)=<P3P.P3P>      C
C EXPECT(9)=<SP(+).SP>  EXPECT(10)=<DP(+).DP>                         C
C EXPECT(11)=<JN.JP>    EXPECT(12)=<P2N.P2P> EXPECT(13)=<P3N.P3P>     C
C                                                                     C
C=====================================================================C
      REAL      EXPECT(13),Z(NS)
      INTEGER   ADJUST(13)
      PARAMETER (MAXP=11221,MAXMULT=1752,MAXPAIR=1168,MAXN=21)
      PARAMETER (MAXNS=4200,MAXNE=400000,MAXLAN=100,MAXDNS=220)
      INTEGER   PAIRENDN,PAIRENDP,PENDN,PENDP,MULTENDP,MULTENDN
      INTEGER   PINDEXN(MAXP),PINDEXP(MAXP)
      INTEGER   PAIRINDEXN(MAXPAIR),PAIRINDEXP(MAXPAIR)
      INTEGER   MULTINDEXN(MAXMULT),MULTINDEXP(MAXMULT)
      INTEGER   STATEJN(MAXNS),STATEAN(MAXNS)
      INTEGER   STATEJP(MAXNS),STATEAP(MAXNS)
      INTEGER   AP_LEFT,AP_RIGHT,AN_LEFT,AN_RIGHT,AN,AP,R
      REAL*8    PTABN(MAXP),PTABP(MAXP),DTERM,SIXJ,DSQRT,PHASE
      REAL      PAIRTABN(MAXPAIR),PAIRTABP(MAXPAIR)
      REAL      MULTTABN(MAXMULT),MULTTABP(MAXMULT)
      LOGICAL   NON_ZERO,LTRIANGLE
C
      COMMON /TAB/ PTABN,PTABP,PAIRTABN,PAIRTABP,MULTTABN,MULTTABP
      COMMON /TABEND/ PAIRENDN,PAIRENDP,PENDN,PENDP,MULTENDP,MULTENDN
      COMMON /TABINDEX/ PINDEXN,PINDEXP,PAIRINDEXN,PAIRINDEXP,
     &                  MULTINDEXN,MULTINDEXP
      COMMON /STATE/ STATEJN,STATEAN,STATEJP,STATEAP
C
      INDP(IAP,JP,R,IA,J)=20000000*J+200000*JP+10000*R+100*IA+IAP
      INDPAIR(IAP,J,R,IA)=1000000*R+10000*J+100*IAP+IA
      PHASE(I)=1-2*MOD(ABS(I),2)
      LTRIANGLE(J1,J2,J3)=(ABS(J1-J2).LE.J3).AND.(J1+J2.GE.J3)
C
      DO I=1,13
       EXPECT(I)=0.
      END DO
      DO 1 I1=1,NS
        IF(ABS(Z(I1)).LT.1E-10) GOTO 1
        JN_LEFT=STATEJN(I1)
        JP_LEFT=STATEJP(I1)
        AP_LEFT=STATEAP(I1)
        AN_LEFT=STATEAN(I1)
        DO 2 I2=1,I1
         IF(ABS(Z(I2)).LT.1E-10) GOTO 2
         ZZ=Z(I1)*Z(I2)
         IF(I1.NE.I2) ZZ=ZZ+ZZ
         JN_RIGHT=STATEJN(I2)
         JP_RIGHT=STATEJP(I2)
         AP_RIGHT=STATEAP(I2)
         AN_RIGHT=STATEAN(I2)
C                                                 neutron interaction
         NON_ZERO=(JP_LEFT.EQ.JP_RIGHT).AND.(AP_LEFT.EQ.AP_RIGHT)
     &       .AND.(JN_LEFT.EQ.JN_RIGHT)
         IF(NON_ZERO) THEN
          IF((I1.EQ.I2).AND.(ADJUST(1).EQ.1)) THEN
           EXPECT(1)=EXPECT(1)+ZZ*JN_LEFT*(JN_LEFT+1)
          END IF
          IF(ADJUST(2).EQ.1) THEN
           INDEX=INDPAIR(AN_LEFT,JN_LEFT,2,AN_RIGHT)
           TERM=MULTTABN(ISECH(MULTINDEXN,MULTENDN,INDEX,1))
           EXPECT(2)=EXPECT(2)+ZZ*TERM
          END IF
          IF(ADJUST(3).EQ.1) THEN
           INDEX=INDPAIR(AN_LEFT,JN_LEFT,3,AN_RIGHT)
           TERM=MULTTABN(ISECH(MULTINDEXN,MULTENDN,INDEX,1))
           EXPECT(3)=EXPECT(3)+ZZ*TERM
          END IF
          IF(ADJUST(4).EQ.1) THEN
           INDEX=INDPAIR(AN_LEFT,JN_LEFT,0,AN_RIGHT)
           TERM=PAIRTABN(ISECH(PAIRINDEXN,PAIRENDN,INDEX,1))
           EXPECT(4)=EXPECT(4)+ZZ*TERM
          END IF
          IF(ADJUST(5).EQ.1) THEN
           INDEX=INDPAIR(AN_LEFT,JN_LEFT,2,AN_RIGHT)
           TERM=PAIRTABN(ISECH(PAIRINDEXN,PAIRENDN,INDEX,1))
           EXPECT(5)=EXPECT(5)+ZZ*TERM
          END IF
         END IF
C                                                 proton interaction
         NON_ZERO=(JN_LEFT.EQ.JN_RIGHT).AND.(AN_LEFT.EQ.AN_RIGHT)
     &       .AND.(JP_LEFT.EQ.JP_RIGHT)
         IF(NON_ZERO) THEN
          IF((I1.EQ.I2).AND.(ADJUST(6).EQ.1)) THEN
           EXPECT(6)=EXPECT(6)+ZZ*JP_LEFT*(JP_LEFT+1)
          END IF
          IF(ADJUST(7).EQ.1) THEN
           INDEX=INDPAIR(AP_LEFT,JP_LEFT,2,AP_RIGHT)
           TERM=MULTTABP(ISECH(MULTINDEXP,MULTENDP,INDEX,1))
           EXPECT(7)=EXPECT(7)+ZZ*TERM
          END IF
          IF(ADJUST(8).EQ.1) THEN
           INDEX=INDPAIR(AP_LEFT,JP_LEFT,3,AP_RIGHT)
           TERM=MULTTABP(ISECH(MULTINDEXP,MULTENDP,INDEX,1))
           EXPECT(8)=EXPECT(8)+ZZ*TERM
          END IF
          IF(ADJUST(9).EQ.1) THEN
           INDEX=INDPAIR(AP_LEFT,JP_LEFT,0,AP_RIGHT)
           TERM=PAIRTABP(ISECH(PAIRINDEXP,PAIRENDP,INDEX,1))
           EXPECT(9)=EXPECT(9)+ZZ*TERM
          END IF
          IF(ADJUST(10).EQ.1) THEN
           INDEX=INDPAIR(AP_LEFT,JP_LEFT,2,AP_RIGHT)
           TERM=PAIRTABP(ISECH(PAIRINDEXP,PAIRENDP,INDEX,1))
           EXPECT(10)=EXPECT(10)+ZZ*TERM
          END IF
         END IF
C                                                 n-p interaction
         IF((I1.EQ.I2).AND.(ADJUST(11).EQ.1)) THEN
          EXPECT(11)=EXPECT(11)+
     &     ZZ*0.5*(J*(J+1)-JP_LEFT*(JP_LEFT+1)-JN_LEFT*(JN_LEFT+1))
         END IF
         DO R=2,3
          NON_ZERO=(ADJUST(10+R).EQ.1).AND.LTRIANGLE(JN_LEFT,R,JN_RIGHT)
     &                        .AND.LTRIANGLE(JP_LEFT,R,JP_RIGHT)
          IF(NON_ZERO) THEN
           DTERM=SIXJ(2*JN_RIGHT+1,2*R+1,2*JN_LEFT+1,
     &                2*JP_LEFT+1,2*J+1,2*JP_RIGHT+1)
           DTERM=DTERM*DSQRT((2*JP_LEFT+1)*(2*JN_LEFT+1.0D0))
     &                *PHASE(JN_RIGHT+JP_LEFT+J)
           INDEX=INDP(AN_LEFT,JN_LEFT,R,AN_RIGHT,JN_RIGHT)
           DTERM=DTERM*PTABN(ISECH(PINDEXN,PENDN,INDEX,1))
           INDEX=INDP(AP_LEFT,JP_LEFT,R,AP_RIGHT,JP_RIGHT)
           DTERM=DTERM*PTABP(ISECH(PINDEXP,PENDP,INDEX,1))
           EXPECT(10+R)=EXPECT(10+R)+ZZ*DTERM
          END IF
         END DO            !R
 2      END DO             !I2
 1    END DO               !I1
      END

      SUBROUTINE HART(X,J,NJ,JOBN,E,Z,NS)
C=====================================================================C
C     Hamiltonian Matrix Setup and Diagonalization                    C
C=====================================================================C
      PARAMETER (MAXP=11221,MAXMULT=1752,MAXPAIR=1168,MAXN=21)
      PARAMETER (MAXNS=4200,MAXNE=400000,MAXLAN=100,MAXDNS=220)
      REAL      X(13),E(MAXNS),Z(1)
      INTEGER   PAIRENDN,PAIRENDP,PENDN,PENDP,MULTENDP,MULTENDN
      INTEGER   PINDEXN(MAXP),PINDEXP(MAXP)
      INTEGER   PAIRINDEXN(MAXPAIR),PAIRINDEXP(MAXPAIR)
      INTEGER   MULTINDEXN(MAXMULT),MULTINDEXP(MAXMULT)
      INTEGER   STATEJN(MAXNS),STATEAN(MAXNS)
      INTEGER   STATEJP(MAXNS),STATEAP(MAXNS)
      INTEGER   AP_LEFT,AP_RIGHT,AN_LEFT,AN_RIGHT,AN,AP,R
      REAL*8    PTABN(MAXP),PTABP(MAXP),DTERM,SIXJ,DSQRT,PHASE
      REAL      PAIRTABN(MAXPAIR),PAIRTABP(MAXPAIR)
      REAL      MULTTABN(MAXMULT),MULTTABP(MAXMULT)
      LOGICAL   NON_ZERO,LTRIANGLE
      INTEGER   JSIZEN(0:2*MAXN),JSIZEP(0:2*MAXN),N1N,N1P
      REAL      AM(MAXNE)                               !should > maxnds**2
      REAL      WK(MAXNS+MAXLAN*(MAXLAN+1))
      INTEGER   IAM(MAXNE)
C
      COMMON /TAB/ PTABN,PTABP,PAIRTABN,PAIRTABP,MULTTABN,MULTTABP
      COMMON /TABEND/ PAIRENDN,PAIRENDP,PENDN,PENDP,MULTENDP,MULTENDN
      COMMON /TABINDEX/ PINDEXN,PINDEXP,PAIRINDEXN,PAIRINDEXP,
     &                  MULTINDEXN,MULTINDEXP
      COMMON /STATE/ STATEJN,STATEAN,STATEJP,STATEAP
      COMMON /OTHER/ JSIZEN,JSIZEP,N1N,N1P,LAN
C
      INDP(IAP,JP,R,IA,J)=20000000*J+200000*JP+10000*R+100*IA+IAP
      INDPAIR(IAP,J,R,IA)=1000000*R+10000*J+100*IAP+IA
      PHASE(I)=1-2*MOD(ABS(I),2)
      LTRIANGLE(J1,J2,J3)=(ABS(J1-J2).LE.J3).AND.(J1+J2.GE.J3)
C
      EN=X(1)
      B2N=X(2)
      B3N=X(3)
      G0N=X(4)
      G2N=X(5)
      EP=X(6)
      B2P=X(7)
      B3P=X(8)
      G0P=X(9)
      G2P=X(10)
      ENP=X(11)
      B2NP=X(12)
      B3NP=X(13)
C                                                 generate states
      NS=0
      DO 1 JN=0,2*N1N
       IF(JSIZEN(JN).LE.0) GOTO 1
       DO 2 JP=0,2*N1P
        IF(JSIZEP(JP).LE.0) GOTO 2
        IF(.NOT.LTRIANGLE(J,JP,JN)) GOTO 2
        DO AN=1,JSIZEN(JN)
         DO AP=1,JSIZEP(JP)
          NS=NS+1
          IF(NS.GT.MAXNS) CALL DIMERR('MAXNS',NS)
          STATEJN(NS)=JN
          STATEJP(NS)=JP
          STATEAP(NS)=AP
          STATEAN(NS)=AN
         END DO
        END DO
 2     END DO
 1    END DO
      IF(NS.LE.0) RETURN
      PRINT*
      PRINT*,'Now J=',J
      PRINT*,'NS=',NS
C                                                 matrix elements
      NE=0.
      DO I1=1,NS
       JN_LEFT=STATEJN(I1)
       JP_LEFT=STATEJP(I1)
       AP_LEFT=STATEAP(I1)
       AN_LEFT=STATEAN(I1)
       DO I2=1,I1
        JN_RIGHT=STATEJN(I2)
        JP_RIGHT=STATEJP(I2)
        AP_RIGHT=STATEAP(I2)
        AN_RIGHT=STATEAN(I2)
        ELEMENT=0.
C                                                 n-p interaction
        DO R=2,3
         IF(R.EQ.2) THEN
          BNP=B2NP
         ELSE
          BNP=B3NP
         END IF
         NON_ZERO=(BNP.NE.0.).AND.LTRIANGLE(JN_LEFT,R,JN_RIGHT)
     &                        .AND.LTRIANGLE(JP_LEFT,R,JP_RIGHT)
         IF(NON_ZERO) THEN
          DTERM=SIXJ(2*JN_RIGHT+1,2*R+1,2*JN_LEFT+1,
     &                2*JP_LEFT+1,2*J+1,2*JP_RIGHT+1)
          DTERM=DTERM*DSQRT((2*JP_LEFT+1)*(2*JN_LEFT+1.0D0))
     &                *PHASE(JN_RIGHT+JP_LEFT+J)
          INDEX=INDP(AN_LEFT,JN_LEFT,R,AN_RIGHT,JN_RIGHT)
          DTERM=DTERM*PTABN(ISECH(PINDEXN,PENDN,INDEX,1))
          INDEX=INDP(AP_LEFT,JP_LEFT,R,AP_RIGHT,JP_RIGHT)
          DTERM=DTERM*PTABP(ISECH(PINDEXP,PENDP,INDEX,1))
          ELEMENT=ELEMENT+BNP*DTERM
         END IF
        END DO                     !R
        IF(I1.EQ.I2) THEN
         ELEMENT=ELEMENT+
     &     ENP*0.5*(J*(J+1)-JP_LEFT*(JP_LEFT+1)-JN_LEFT*(JN_LEFT+1))
        END IF
C                                                 proton interaction
        NON_ZERO=(JN_LEFT.EQ.JN_RIGHT).AND.(AN_LEFT.EQ.AN_RIGHT)
     &       .AND.(JP_LEFT.EQ.JP_RIGHT)
        IF(NON_ZERO) THEN
         IF(I1.EQ.I2) THEN
          ELEMENT=ELEMENT+EP*JP_LEFT*(JP_LEFT+1)
         END IF
         IF(B2P.NE.0.) THEN
          INDEX=INDPAIR(AP_LEFT,JP_LEFT,2,AP_RIGHT)
          TERM=MULTTABP(ISECH(MULTINDEXP,MULTENDP,INDEX,1))
         ELEMENT=ELEMENT+B2P*TERM
         END IF
         IF(B3P.NE.0.) THEN
          INDEX=INDPAIR(AP_LEFT,JP_LEFT,3,AP_RIGHT)
          TERM=MULTTABP(ISECH(MULTINDEXP,MULTENDP,INDEX,1))
          ELEMENT=ELEMENT+B3P*TERM
         END IF
         IF(G0P.NE.0.) THEN
          INDEX=INDPAIR(AP_LEFT,JP_LEFT,0,AP_RIGHT)
          TERM=PAIRTABP(ISECH(PAIRINDEXP,PAIRENDP,INDEX,1))
          ELEMENT=ELEMENT+G0P*TERM
         END IF
         IF(G2P.NE.0.) THEN
          INDEX=INDPAIR(AP_LEFT,JP_LEFT,2,AP_RIGHT)
          TERM=PAIRTABP(ISECH(PAIRINDEXP,PAIRENDP,INDEX,1))
          ELEMENT=ELEMENT+G2P*TERM
         END IF
        END IF
C                                                 neutron interaction
        NON_ZERO=(JP_LEFT.EQ.JP_RIGHT).AND.(AP_LEFT.EQ.AP_RIGHT)
     &       .AND.(JN_LEFT.EQ.JN_RIGHT)
        IF(NON_ZERO) THEN
         IF(I1.EQ.I2) THEN
          ELEMENT=ELEMENT+EN*JN_LEFT*(JN_LEFT+1)
         END IF
         IF(B2N.NE.0.) THEN
          INDEX=INDPAIR(AN_LEFT,JN_LEFT,2,AN_RIGHT)
          TERM=MULTTABN(ISECH(MULTINDEXN,MULTENDN,INDEX,1))
         ELEMENT=ELEMENT+B2N*TERM
         END IF
         IF(B3N.NE.0.) THEN
          INDEX=INDPAIR(AN_LEFT,JN_LEFT,3,AN_RIGHT)
          TERM=MULTTABN(ISECH(MULTINDEXN,MULTENDN,INDEX,1))
          ELEMENT=ELEMENT+B3N*TERM
         END IF
         IF(G0N.NE.0.) THEN
          INDEX=INDPAIR(AN_LEFT,JN_LEFT,0,AN_RIGHT)
          TERM=PAIRTABN(ISECH(PAIRINDEXN,PAIRENDN,INDEX,1))
          ELEMENT=ELEMENT+G0N*TERM
         END IF
         IF(G2N.NE.0.) THEN
          INDEX=INDPAIR(AN_LEFT,JN_LEFT,2,AN_RIGHT)
          TERM=PAIRTABN(ISECH(PAIRINDEXN,PAIRENDN,INDEX,1))
          ELEMENT=ELEMENT+G2N*TERM
         END IF
        END IF
        IF(NS.LE.MAXDNS) THEN
         AM(I1+NS*(I2-1))=ELEMENT
         AM(I2+NS*(I1-1))=ELEMENT
        ELSE IF(ABS(ELEMENT).GT.SMALL) THEN
         NE=NE+1
         IF(NE.GT.MAXNE) CALL DIMERR('MAXNE',NE)
         AM(NE)=ELEMENT
         IAM(NE)=(NS+1)*I1+I2
        END IF
       END DO          !I2
      END DO           !I1
      PRINT*,'number of non-zero matrix elements=',ne
      PRINT*,'DIAGNOLIZATION........'
      IF(NS.LE.MAXDNS) THEN
       M=MIN(NS,NJ)
C      CALL EIGRS(AM,NS,JOBN,E,Z,NS,WK,IER)
       IF(JOBN.EQ.0) THEN
        CALL E2LSF(NS,AM,NS,E,AM,WK)
       ELSE
        CALL E2CSF(NS,AM,NS,E,Z,NS,WK)
        DO II=1,NS
         SUM=0.
         DO I=1,NS
          SUM=SUM+Z(I+(II-1)*NS)**2
         END DO
         SUM=SQRT(SUM)
         DO I=1,NS
          Z(I+(II-1)*NS)=Z(I+(II-1)*NS)/SUM
         END DO
        END DO
       END IF
      ELSE
       M=MIN(NS,NJ,LAN)
       LANT=MIN(NS,LAN)
       CALL EIGRL(AM,IAM,NE,NS,LANT,M,JOBN,E,Z,NS,WK,IER)
      END IF
      RETURN
      END
