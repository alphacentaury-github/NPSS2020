      PROGRAM QMODL2
**********************************************************************
*    THIS PROGRAM WAS WRITTEN BY:                DECEMBER, 1990      *
*                                                                    *
*       U. STRAUB AND A. FAESSLER,                                   *
*       INSTITUTE FOR THEORETICAL PHYSICS,                           *
*       UNIVERSITY OF TUEBINGEN,                                     *
*     D-7400 TUEBINGEN, GERMANY                                      *
*                                                                    *
*====================================================================*
*                                                                    *
*    PROGRAM QMODL2  COMPUTES THE NUCLEON-NUCLEON SCATTERING         *
*    IN THE NON-RELATIVISTIC QUARK MODEL USING THE RESONATING        *
*    GROUP METHOD. THIS PROGRAM CONSISTS OF TWO PARTS:               *
*                                                                    *
*    PART1: HERE ALL THE KERNELS ARE COMPUTED.                       *
*                                                                    *
*    PART2: HERE THE RESONATING GROUP EQUATIONS ARE SOLVED           *
*           USING THE KERNELS OF PART1.                              *
*                                                                    *
*   IMPORTANT VARIABLES:                                             *
*                                                                    *
*   B:       WIDTH OF GAUSSIAN QUARK CORE OF SINGLE NUCLEON (FM)     *
*   QM:      CONSITUENT QUARK MASS   (MEV)                           *
*   ALFS:    GLUON-QUARK COUPLING STRENGTH                           *
*   ALIN:    STRENGTH OF LINEAR CONFINEMENT (MEV/FM)                 *
*   AQUAD:   STRENGTH OF QUADR. CONFINEMENT (MEV/FM**2)              *
*   IOBEP:   TYPE OF MESON EXCHANGE, ALLOWED VALUES: 0, 1, 2.        *
*            0: NO MESON EXCHANGE.                                   *
*            1: MESON EXCHANGE ON BARYON LEVEL.                      *
*            2: SIGMA-MESON EXCHANGE ON BARYON LEVEL,                *
*               PION EXCHANGE ON QUARK LEVEL.                        *
*   PARP:    SCALING FACTORS FOR PION EXCHANGE ON BARYON LEVEL       *
*   PARS:       -"-              SIGMA         -"-                   *
*   PARW:       -"-              OMEGA         -"-                   *
*   PARR:       -"-              RHO           -"-                   *
*   NOCHNN:  NUMBER OF OPEN+CLOSED CHANNELS  (ALLOWED: 1, 2, 3)      *
*   SPIN:    SPIN OF NN SYSTEM (0.0  OR  1.0)                        *
*   LPW:     PARTIAL WAVE QUANTUM NUMBER  (HERE LPW=0 ONLY)          *
*   NS:      NUMBER OF MESH POINTS (MUST BE .LE. 10)                 *
*   SMIN:    FIRST MESH POINT (FM)                                   *
*   DS:      SPLITTING BETWEEN MESH POINTS (FM)                      *
*                                                                    *
*   RC:      MATCHING RADIUS (FM)                                    *
*   NE:      NUMBER OF ENERGIES (MUST BE .LE. 24)                    *
*   EREL:    LIST OF ENERGIES (IN CM-FRAME), FOR WHICH               *
*            NN SCATTERING SHOULD BE COMPUTED.                       *
*                                                                    *
*                                                                    *
*   ORDER OF CHANNELS IS AS FOLLOWS:                                 *
*                                                                    *
*    NOCHNN =  1:   ONE-CHANNEL CALCULATION                          *
*          ICHANN=1 :     <NUCNUC/NUCNUC>                            *
*                                                                    *
*    NOCHNN =  2:   TWO-CHANNEL CALCULATION                          *
*          ICHANN=1 :     <NUCNUC/NUCNUC>                            *
*                 2 :     <NUCNUC/DELDEL>                            *
*                 3 :     <DELDEL/DELDEL>                            *
*                                                                    *
*    NOCHNN =  3:   THREE-CHANNEL CALCULATION                        *
*          ICHANN=1 :     <NUCNUC/NUCNUC>                            *
*                 2 :     <NUCNUC/DELDEL>                            *
*                 3 :     <DELDEL/DELDEL>                            *
*                 4 :     <NUCNUC/HIDHID>                            *
*                 5 :     <DELDEL/HIDHID>                            *
*                 6 :     <HIDHID/HIDHID>                            *
*                                                                    *
*   PLEASE NOTE: WE ARE USING THIS ORDERING OF CHANNELS, AND NOT     *
*                THE ORDER DEFINED IN SUBROUTINE DCHAN ||||||        *
*====================================================================*
*                                                                    *
*   I M P O R T A N T:                                               *
*                                                                    *
*   THIS PROGRAM NEEDS ADDITIONAL ROUTINES FROM THE NAG LIBRARY:     *
*                                                                    *
*   FROM NAG (VERSION 13)                                            *
*   --------                                                         *
*             E01BAF (CUBIC SPLINE INTERPOLATION, DOUBLE PREC.)      *
*             E02BBF (CUBIC SPLINE EVALUATION,    DOUBLE PREC.)      *
*             S15AEF (ERROR FUNCTION, DOUBLE PRECISION)              *
*             F04ADF (SOLVES SET OF LINEAR EQUATIONS  A X = B,       *
*                     A = GIVEN N BY N MATRIX (COMPLEX)              *
*                     X = UNKNOWN N-DIM. VECTOR.                     *
*                     B = GIVEN N-DIM. VECTOR (COMPLEX)              *
*                     DOUBLE PRECISION VERSION)                      *
*                                                                    *
**********************************************************************
*
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      IMPLICIT COMPLEX*16(C)
      CHARACTER*32 TFRMT
      PARAMETER (MESH=10,MAXCHN=3,NENRGY=24)
*****************  DIMENSIONS FOR PART1 *****************
      DIMENSION SMIN(3),DS(3),RCSTD(2),RCSTE(5),PST(2)
      DIMENSION NTD(2),NTE(5),STCD(6,5),STCE(6,11),DTT(2),ATT(5)
      DIMENSION TX(100),SY(100)
*****************  DIMENSIONS FOR PART2 *****************
      DIMENSION RK(MAXCHN,MAXCHN,MESH,MESH)
      DIMENSION W(MAXCHN,MAXCHN,MESH,MESH)
      DIMENSION BIGNOR(MAXCHN,MAXCHN,MESH,MESH)
      DIMENSION BIGHAM(MAXCHN,MAXCHN,MESH,MESH)
      DIMENSION CL(MAXCHN,MAXCHN,MESH,MESH)
      DIMENSION ET(MAXCHN),CPX(MESH),CSX(MESH),CP(MAXCHN,MESH)
      DIMENSION CS(MAXCHN,MESH),ZASAV(MAXCHN*MESH,MAXCHN*MESH)
      DIMENSION ZAA(MAXCHN*MESH,MAXCHN*MESH)
      DIMENSION ZBB(MAXCHN*MESH),ZM(MAXCHN,MESH,MAXCHN)
      DIMENSION ZPCSAV(MAXCHN*MESH),ZPC(MAXCHN*MESH,MAXCHN)
      DIMENSION CC(MAXCHN,MESH,MAXCHN)
      DIMENSION KEXC(MAXCHN)
      DIMENSION CSMTR(MAXCHN,MAXCHN),CSSADJ(MAXCHN,MAXCHN)
      DIMENSION CSMAT(MAXCHN,MAXCHN,NENRGY)
      DIMENSION EINT(MAXCHN,MAXCHN),RKEY(MAXCHN)
      DIMENSION REDMAS(MAXCHN),BPRIME(MAXCHN)
      DIMENSION EREL(NENRGY),DELTA(NENRGY,MAXCHN),RKSAVE(NENRGY,MAXCHN)
      DIMENSION DEVIAT(NENRGY),NOPARR(NENRGY)
      DIMENSION ZCC(20)
      DIMENSION WKSPC(MAXCHN*MESH)
*********************
      COMMON /A/ RRLM,SWW(MESH),RC,RKE
      COMMON/KEXCO/ KEXC
      COMMON /A1/ RLM(MESH,MESH,6)
      COMMON /A2/ RGM(MESH,MESH,6)
      COMMON /A3/ OMG(12),ETA(12),NDIM
      COMMON /A4/ XCOF(13),XEXP(13),LDIM
      COMMON /B/ A(3),BS(5,3,6)
      COMMON /L/ LPW
      COMMON /WIDTH/ B
      COMMON /P/ PARP,PARS,PARW,PARR
      COMMON /KERNEL/ SW(MESH,3),RN(6),RKIN(6),RCF(6),RCOU(6),RD(6),
     &                RCMI(6),RPION,OBEP,RHAM(6)
      DATA NTD,NTE/6,9,2,1,4,4,4/
      DATA RMU/0.7D0/, RCUTOF/4.2D0/
      DATA IANTS /1/
      DATA KOBEP,HOBEP/100,0.05D0/
      DATA HBC /197.32891D0/
      DATA RCONST/0.0/
      DATA ICOUPL /1/
      DATA ITEST  /0/
*/////////////////////////////////////////////
*       PART 1 : CALCULATION OF KERNELS     //
*/////////////////////////////////////////////
*
******************************************************
*   SUPPRESS ERROR MESSAGE 208 = UNDERFLOW EXCEPTION
*     CALL ERRSET(208,0,-1,1,0,208)
******************************************************
      PHI=4.*DATAN(1D0)
*
*     INPUT
*
      READ(5,*) B
      READ(5,*) QM
      READ(5,*) ALFS
      READ(5,*) ALIN
      READ(5,*) AQUAD
      READ(5,*) IOBEP
      READ(5,*) PARP
      READ(5,*) PARS
      READ(5,*) PARW
      READ(5,*) PARR
      READ(5,*) NOCHNN
      READ(5,*) SPIN
      READ(5,*) NS
      READ(5,*) SMIN(1)
      READ(5,*) DS(1)
*
      IF (IOBEP.NE.0  .AND.  NOCHNN.GT.1) THEN
         WRITE(6,*)  ' IOBEP AND NOCHNN ARE NOT CONSISTENT ||||'
         WRITE(6,*)  ' IOBEP  .NE. 0  ===> NOCHNN .EQ. 1 '
         WRITE(6,*)  ' NOCHNN .GT. 1  ===> IOBEP  .EQ. 0 '
         WRITE(6,*)  ' PROGRAM IS STOPPED |'
         STOP
      END IF
      SMIN(2)=SMIN(1)
      SMIN(3)=SMIN(1)
      DS(2)=DS(1)
      DS(3)=DS(1)
      LPW=0
      IF (LPW-(LPW/2)*2.EQ.0) THEN
*************    LPW  =  EVEN   ***********
         TISOSP=1.0-SPIN
      ELSE
*************    LPW  =  ODD    ***********
         TISOSP=SPIN
      END IF
      IST=SPIN+TISOSP+0.001
      NC=NOCHNN*(NOCHNN+1)/2
*
      IF (ABS(ALIN).LT.1D-2) THEN
         ICONF=2
         ACONF=AQUAD
      ELSE IF (ABS(AQUAD).LT.1D-2) THEN
         ICONF=1
         ACONF=ALIN
      END IF
*
*  DIRECT AND EXCHANGE COEFFICIENTS FOR SPIN-ISOSPIN-COLOR
*
      CALL SFCME(SPIN,STCD,STCE)
*
* COEFFICIENTS FOR PION EXCHANGE
*
      CALL COFP(IST,RCSTD,RCSTE,STCD,STCE)
*
      WRITE(6,'(/,20X,A)')   '*****************************'
      WRITE(6,'(20X,A)')     '*                           *'
      IF (ABS(SPIN).LT.1D-4)
     &WRITE(6,'(20X,A)')     '*    SINGLET CALCULATION    *'
      IF (ABS(SPIN-1.0).LT.1D-4)
     &WRITE(6,'(20X,A)')     '*    TRIPLET CALCULATION    *'
      WRITE(6,'(20X,A)')     '*                           *'
      WRITE(6,'(20X,A)')     '*****************************'
      WRITE(6,'(//,10X,A,F7.3,A)')
     &   'SPIN:                              SPIN = ',SPIN
      WRITE(6,'(10X,A,I3)')
     &   'PARTIAL WAVE QUANTUM NUMBER:        LPW = ',LPW
      WRITE(6,'(10X,A,F7.3,A)')
     &   'WIDTH OF GAUSSIAN QUARK CORE:         B = ',B,' (FM)'
      WRITE(6,'(10X,A,F7.3,A)')
     &   'CONSTITUENT QUARK MASS:             QM  = ',QM,' (MEV)'
      WRITE(6,'(10X,A,F7.3,A)')
     &   'GLUON COUPLING CONSTANT:         ALFA_S = ',ALFS
      WRITE(6,'(10X,A,F7.3,A)')
     &   'LINEAR CONFINEMENT STRENGTH:      ALIN  = ',ALIN,' (MEV/FM)'
      WRITE(6,'(10X,A,F7.3,A)')
     &   'QUADR. CONFINEMENT STRENGTH:      AQUAD = ',AQUAD,
     &   ' (MEV/FM**2)'
      WRITE(6,'(12X,A,I3,A)')
     &   '===>   CONFINEMENT TYPE:        ICONF = ',ICONF
      WRITE(6,'(19X,A)')    '(1 = LIN., 2 = QUADR.)'
      WRITE(6,'(10X,A,F7.3,A)')
     &   'FIRST MESH POINT:           SMIN = S(1) = ',SMIN(1),' (FM)'
      WRITE(6,'(10X,A,F7.3,A)')
     &   'SPLITTING OF MESH POINTS: S(I+1) - S(I) = ',DS(1),' (FM)'
      WRITE(6,'(10X,A,I3)')
     &   'NUMBER OF MESH POINTS:               NS = ',NS
      WRITE(6,'(10X,A,I3)')
     &   'NUMBER OF CHANNELS:              NOCHNN = ',NOCHNN
      IF (NOCHNN.EQ.1)  WRITE(6,'(12X,A)')
     &  '===>   NUCLEON-NUCLEON CHANNEL ONLY'
      IF (NOCHNN.EQ.2)  WRITE(6,'(12X,A)')
     &  '===>   NUCLEON-NUCLEON AND DELTA-DELTA CHANNEL '
      IF (NOCHNN.EQ.3)  WRITE(6,'(12X,A)')
     &  '===>   NUCLEON-NUCLEON, DELTA-DELTA AND HIDDEN COLOR CHANNEL '
      WRITE(6,'(10X,A,I3)')
     &   'PARAMETER FOR MESON EXCHANGE:     IOBEP = ',IOBEP
      IF (IOBEP.EQ.0)  WRITE(6,'(12X,A)')
     &  '===>   NO MESON EXCHANGE AT ALL'
      IF (IOBEP.EQ.1)  WRITE(6,'(12X,A)')
     &  '===>   MESON EXCHANGE ON BARYON LEVEL '
      IF (IOBEP.EQ.2)  WRITE(6,'(12X,A,/)')
     &  '===>   SIGMA EXCH. ON BARYON LEVEL, PION EXCH. ON QUARK LEVEL'
      IF (IANTS.EQ.0 .AND. IOBEP.EQ.2)  WRITE(6,'(12X,A,/)')
     &  'NO ANTISYMMETRISATION FOR SIMULTANEOUS PION-QUARK EXCHANGE ||'
      WRITE(6,'(/,20X,A)')   '========================================'
      WRITE(6,'(20X,A)')     'I                                      I'
      WRITE(6,'(20X,A)')     'I   PART 1:  CALCULATION  OF KERNELS   I'
      WRITE(6,'(20X,A)')     'I                                      I'
      WRITE(6,'(20X,A,/)')   '========================================'
*
      WRITE(6,'(/,5X,A,/)') 'DIRECT COEFF. FOR SPIN-ISOSPIN-COLOR'
      DO 222  K=1,5
222   WRITE(6,'(1X,I5,6D13.5)') K,(STCD(J,K),J=1,6)
*
      WRITE(6,'(/,5X,A,/)') 'EXCHANGE COEFF. FOR SPIN-ISOSPIN-COLOR'
      DO 333  K=1,11
333   WRITE(6,'(1X,I5,6D13.5)') K+5,(STCE(J,K),J=1,6)
*
      WRITE(6,'(/,5X,A,/)') 'COEFF. FOR GAUSS-LAGUERRE INTEGRATION'
      DO 444  K=1,12
444   WRITE(6,'(1X,I5,6D13.5)') K,OMG(K),ETA(K)
*
*  COEFFICIENTS FOR GAUSS-LAGUERRE-INTEGRATION
*
      DO 19 IK=1,NDIM
      OMG(IK)=2*OMG(IK)/DSQRT(PHI)
      XCOF(IK+1)=OMG(IK)
      ETA(IK)=ETA(IK)**2
      XEXP(IK+1)=ETA(IK)
19    CONTINUE
      XEXP(1)=1D0
      XCOF(1)=2./DSQRT(PHI)
*  SETTING-UP MESHPOINTS IN GENERATOR-COORDINATE SPACE
      DO 1 K=1,NS
      DO 1 KC=1,3
      SW(K,KC)=SMIN(KC)+DS(KC)*(K-1)
      SWW(K)=SMIN(1)+DS(1)*(K-1)
1     CONTINUE
      B2=B*B
      DO 2 I=1,NS
      DO 2 J=1,NS
      DO 2 K1=1,3
      DO 2 K2=K1,3
      K=K1+K2
      IF(K1.EQ.1) K=K-1
      RLM(I,J,K)=(SW(I,K1)**2+SW(J,K2)**2)/B2
2     RGM(I,J,K)=SW(I,K1)*SW(J,K2)/B2
      RK0=3/(4*QM*B2)* HBC**2
      DEN=(2.*PHI*B2)**(-1.5)
      ALFD=-0.25*ALFS* HBC**3 *DEN*PHI/(QM*QM)
      VD0=2.*ALFD/3.
      ALFS=ALFS*HBC/4
      RKX= HBC**2 /(3*QM)
      FQ=DEXP(-RMU*RMU*B*B/3.)
      VYP=(3*0.08*FQ/25.D0)*HBC/(DSQRT(2.D0)*B)
      VPDEL=(3*0.08*FQ/25.D0) *RMU*HBC *4.0*PHI /RMU**3 *DEN
************************************
*  LOOP OVER MESH POINTS SI, SJ    *
************************************
      DO 3 IS=1,NS
      DO 3 ISP=1,NS
      DO 4 I=1,NC
      A(1)=RGM(IS,ISP,I)*.75
      A(2)=RGM(IS,ISP,I)*.25
      A(3)=RGM(IS,ISP,I)*.5
4     CALL BSL(LPW+1,I)
*
*     NORM
*
      DO 5 NCH=1,NC
      DT=RNORM(IS,ISP,0,NCH)
      AT=RNORM(IS,ISP,1,NCH)
      IF(NCH.EQ.1) DOT=DT
      IF(NCH.EQ.1) AOT=AT
      DTT(1)=STCD(NCH,1)*DT
      ATT(1)=STCE(NCH,1)*AT
5     RN(NCH)=DTT(1)-9*ATT(1)
*
*     KINETIC ENERGY
*
      DO 6 NCH=1,NC
      DT=TKIN(IS,ISP,0,RK0,NCH)
      AT=TKIN(IS,ISP,1,RK0,NCH)
      DTT(1)=STCD(NCH,1)*DT
      ATT(1)=STCE(NCH,1)*AT
*     WRITE(6,*) SW(IS,1),SW(ISP,1),DTT(1),ATT(1)
6     RKIN(NCH)=DTT(1)-9*ATT(1)
*
*     CONFINEMENT
*
      DO 7 NCH=1,NC
      DT=0.
      AT=0.
      IF(ACONF.EQ.0D0) GO TO 7
      DO 8 I=1,2
      DTT(I)=STCD(NCH,I+1)*VCONF(IS,ISP,0,I,B2,ACONF,NCH,ICONF)*NTD(I)
8     DT=DT+DTT(I)
      DO 9 I=1,5
      ATT(I)=STCE(NCH,I+1)*VCONF(IS,ISP,1,I,B2,ACONF,NCH,ICONF)*NTE(I)
9     AT=AT+ATT(I)
7     RCF(NCH)=DT-9*AT
*
*  FOLLOWING ONLY IF A CONSTANT TERM IS TO BE ADDED TO THE CONFINEMENT
*
      IF(RCONST.EQ.0D0) GO TO 72
      DO 71 NCH=1,NC
      DT=0.
      AT=0.
      DO 81 I=1,2
      DTT(I)=STCD(NCH,I+1)*RCONST*RNORM(IS,ISP,0,NCH)*NTD(I)
81    DT=DT+DTT(I)
      DO 91 I=1,5
      ATT(I)=STCE(NCH,I+1)*RCONST*RNORM(IS,ISP,1,NCH)*NTE(I)
91    AT=AT+ATT(I)
71    RCF(NCH)=RCF(NCH)+DT-9*AT
72    CONTINUE
*
*   SPIN INDEPENDENT 1/R -TERM (COULOMB) IN ONE-GLUON EXCHANGE
*
      DO 10 NCH=1,NC
      DT=0.
      AT=0.
      IF(ALFS.EQ.0D0) GO TO 10
      DO 11 I=1,2
      DTT(I)=STCD(NCH,I+1)*VC(IS,ISP,0,I,ALFS,B2,NCH)*NTD(I)
11    DT=DT+DTT(I)
      DO 12 I=1,5
      ATT(I)=STCE(NCH,I+1)*VC(IS,ISP,1,I,ALFS,B2,NCH)*NTE(I)
12    AT=AT+ATT(I)
10    RCOU(NCH)=DT-9*AT
*
*   SPIN INDEPENDENT DELTA-TERM IN ONE-GLUON EXCHANGE
*
      DO 61 NCH=1,NC
      DT=0.
      AT=0.
      DO 51 I=1,2
51    DT=DT+STCD(NCH,I+1)*VD(IS,ISP,0,I,ALFD,NCH)*NTD(I)
      DO 41 I=1,5
41    AT=AT+STCE(NCH,I+1)*VD(IS,ISP,1,I,ALFD,NCH)*NTE(I)
61    RD(NCH)=DT-9*AT
*
*   COLOR MAGNETIC PART OF ONE-GLUON EXCHANGE
*
      DO 13 NCH=1,NC
      DT=0.
      AT=0.
      DO 14 I=1,2
      DTT(I)=STCD(NCH,I+3)*VD(IS,ISP,0,I,VD0,NCH)*NTD(I)
14    DT=DT+DTT(I)
      DO 15 I=1,5
      ATT(I)=STCE(NCH,I+6)*VD(IS,ISP,1,I,VD0,NCH)*NTE(I)
15    AT=AT+ATT(I)
13    RCMI(NCH)=DT-9*AT
*
*     PION EXCHANGE INTERACTION
*
      DATA PST/0.777778D0,0.382716D0/
      NPST=(2+IST)/2
      RN1=RN(1)
      IF(LPW.EQ.1) RN1=RNORM(IS,ISP,0,1)-3*PST(NPST)*
     &RNORM(IS,ISP,1,1)
*
      DTY=0.
      ATY=0.
      FACTOR=1./3.
      IDELTA=1
*  IDELTA DETERMINES THE FORM OF THE CUTOFF PART OF THE PION POTENTIAL
*  IDELTA=1:  CUTOFF = DELTA FUNCTION
*  IDELTA=0:  CUTOFF = YUKAWA WITH CUTOFF MASS RCUTOF = 4.2 (1/FM)
      DO 32 I=1,2
      IF (IDELTA.EQ.1) THEN
         ORBITL = VY2(IS,ISP,0,I,VYP,1,1)
     &          - VD(IS,ISP,0,I,VPDEL,1)
      ELSE IF (IDELTA.EQ.0) THEN
         ORBITL = VY2(IS,ISP,0,I,VYP,1,1)
     &          + VY2(IS,ISP,0,I,VYP,1,2) *((RCUTOF/RMU)**2)
      END IF
 32   DTY=DTY+NTD(I)*RCSTD(I)*ORBITL
      DO 33 I=1,5
      IF (IDELTA.EQ.1) THEN
         ORBITL = VY2(IS,ISP,1,I,VYP,1,1)
     &          - VD(IS,ISP,1,I,VPDEL,1)
      ELSE IF (IDELTA.EQ.0) THEN
         ORBITL = VY2(IS,ISP,1,I,VYP,1,1)
     &          + VY2(IS,ISP,1,I,VYP,1,2) *((RCUTOF/RMU)**2)
      END IF
 33   ATY=ATY+NTE(I)*RCSTE(I)*FACTOR*ORBITL
      RPION=DTY-9*ATY
      IF(IANTS.EQ.0) RPION=DTY
      IF(IOBEP.NE.2) RPION=0.0
*
*     ONE BOSON EXCHANGE
*
      IF(IS.NE.1.OR.ISP.NE.1) GO TO 2222
      RLAM=3./(4.*B2)
* CUBIC SPLINE INTERPOLATION BY MEANS OF NAG-LIB
      CALL VOBEP2(IOBEP,KOBEP,HOBEP,TX,SY,SPIN)
*     WRITE(6,*)  'TX(1)   ,SY(1)   = ',TX(1),SY(1)
*     WRITE(6,*)  'TX(100) ,SY(100) = ',TX(100),SY(100)
2222  OBEP=0D0
      OBF=4*PHI*((3/(2.*PHI*B2))**1.5)
      SOI=SW(IS,1)
      SOJ=SW(ISP,1)
      OBF=OBF*DEXP(-RLAM*(SOI**2+SOJ**2))
      DO 23 IB=1,KOBEP
      FAKB=1.
***  WEIGHT OF FIRST POINT CHANGED FROM 0.5 TO 1.0  SINCE   ***
***  FIRST POINT OF NUMERICAL INTEGRATION (R=0) IS OMITTED. ***
***  DUE TO FACTOR R*R THE INTEGRAND IS ZERO AT THIS POINT  ***
      IF(IB.EQ.KOBEP) FAKB=.5
      X1=2*RLAM*TX(IB)*SOI
      X2=2*RLAM*TX(IB)*SOJ
      FEXP=-2.0*RLAM*TX(IB)*TX(IB)
 23   OBEP=OBEP+FAKB*SY(IB)*TX(IB)*TX(IB)*DEXP(FEXP)*
     &BSS(X1,LPW)*BSS(X2,LPW)
      OBEP=OBEP*OBF*HOBEP
      IF (IOBEP.EQ.0) OBEP=0.0
*
*       OUTPUT
*
      IF(IS.NE.0.OR.ISP.NE.1) GO TO 27
      DO 26 KO=1,KOBEP
  26    PRINT 107,TX(KO),SY(KO)
107   FORMAT(10X,2F17.5)
27    CONTINUE
      DO 67 K=1,NC
  67  RHAM(K)=RKIN(K)+RCF(K)+RCOU(K)+RD(K)+RCMI(K)
     &       +RPION+OBEP
*
      DO 205 NGAMM=1,NOCHNN
      DO 205 NDELT=NGAMM,NOCHNN
      IF ((NGAMM.EQ.1) .AND. (NDELT.EQ.1))  THEN
         INDEX=1
      ELSE IF ((NGAMM.EQ.1) .AND. (NDELT.EQ.2))  THEN
         INDEX=2
      ELSE IF ((NGAMM.EQ.2) .AND. (NDELT.EQ.2))  THEN
         INDEX=3
      ELSE IF ((NGAMM.EQ.1) .AND. (NDELT.EQ.3))  THEN
         INDEX=4
      ELSE IF ((NGAMM.EQ.2) .AND. (NDELT.EQ.3))  THEN
         INDEX=5
      ELSE IF ((NGAMM.EQ.3) .AND. (NDELT.EQ.3))  THEN
         INDEX=6
      END IF
      BIGNOR(NGAMM,NDELT,IS,ISP)=RN(INDEX)
      BIGNOR(NDELT,NGAMM,ISP,IS)=RN(INDEX)
      BIGHAM(NGAMM,NDELT,IS,ISP)=RHAM(INDEX)
      BIGHAM(NDELT,NGAMM,ISP,IS)=RHAM(INDEX)
205   CONTINUE
      IF (IS.EQ.1 .OR. ISP.EQ.1 .OR. IS.EQ.ISP)
     &          CALL KEROUT(IS,ISP,NC)
******************************************
*  END OF LOOP OVER MESH POINTS SI,SJ *
******************************************
3     CONTINUE
      RJ=SW(NS,1)
      RKEYY=RK0*(1./3.+LPW*(8*B2/(9*RJ*RJ)+32*B2*B2/(27*RJ**4)))
      EINT1=RHAM(1)/RN(1)
      IF (NC.EQ.1) THEN
         WRITE(6,'(5X,A,3D12.5)')   'EINT1: ',EINT1
         WRITE(6,'(5X,A,3D12.5)')   'RKEYY: ',RKEYY
      ELSE IF (NC.EQ.3) THEN
         EINT2=RHAM(3)/RN(3)
         WRITE(6,'(5X,A,3D12.5)')   'EINT1, EINT2: ',EINT1,EINT2
         WRITE(6,'(5X,A,3D12.5)')   'RKEYY: ',RKEYY
      ELSE IF (NC.EQ.6) THEN
         EINT2=RHAM(3)/RN(3)
         EINT3=RHAM(6)/RN(6)
         WRITE(6,'(5X,A,3D12.5)')   'EINT1, EINT2, EINT3: ',
     &                               EINT1, EINT2, EINT3
         WRITE(6,'(5X,A,3D12.5)')   'RKEYY: ',RKEYY
      END IF
*/////////////////////////////////////////////
*       PART 2 : SOLVE THE RGM EQUATION     //
*/////////////////////////////////////////////
*
**************************
*       INPUT            *
**************************
*
*     READ(5,*) ITEST
      READ(5,*) RC
      READ(5,*) NE
      IF (NE.LE.8) THEN
         READ(5,*) (EREL(I),I=1,NE)
      ELSE IF (NE.GT.8  .AND. NE.LE.16) THEN
         READ(5,*) (EREL(I),I=1,8)
         READ(5,*) (EREL(I),I=9,NE)
      ELSE IF (NE.GT.16 .AND. NE.LE.24) THEN
         READ(5,*) (EREL(I),I=1,8)
         READ(5,*) (EREL(I),I=9,16)
         READ(5,*) (EREL(I),I=17,NE)
      END IF
*
      WRITE(6,'(//,20X,A)')  '========================================'
      WRITE(6,'(20X,A)')     'I                                      I'
      WRITE(6,'(20X,A)')     'I   PART 2:  SOLUTION OF RGM EQUATION  I'
      WRITE(6,'(20X,A)')     'I                                      I'
      WRITE(6,'(20X,A,/)')   '========================================'
      IF (ITEST.NE.0)
     &WRITE(6,'(/,10X,A,I3,/)') 'TEST PARAMETER:    ITEST = ',ITEST
      WRITE(6,'(10X,A,F7.3,A)')
     &          'MATCHING RADIUS:      RC = ',RC,' (FM)'
      WRITE(6,'(10X,A,I3)')     'NUMBER OF ENERGIES:   NE = ',NE
      WRITE(6,'(10X,A,8F8.3,:,/,22X,8F8.3,:,/,22X,8F8.3)')
     &    'ENERGIES EREL:',(EREL(I),I=1,NE)
*
      DO 113 N=1,NOCHNN
      KEXC(N)=0
113   CONTINUE
*
************************************
*       LOOP  OVER  ENERGIES       *
************************************
*
      DO 733 I33=1,NE
*
*****    SET ALL MATRICES TO ZERO.
*
      DO 133 LGAMM=1,MAXCHN
      DO 133 LDELT=1,MAXCHN
      DO 133 L1=1,MESH
      DO 133 L2=1,MESH
      RK(LGAMM,LDELT,L1,L2)=0.D0
      CL(LGAMM,LDELT,L1,L2)=0.D0
      ZASAV((LGAMM-1)*MESH+L1,(LDELT-1)*MESH+L2)=0.D0
      ZAA((LGAMM-1)*MESH+L1,(LDELT-1)*MESH+L2)=0.D0
      ZBB((LGAMM-1)*MESH+L1)=0.D0
      ZM(LGAMM,L1,LDELT)=0.D0
133   CONTINUE
*
*****   DEFINE  REDUCED MASSES, RELATIVE KIN. ENERGIES AND
*****   RELATIVE MOMENTA FOR ALL  DIAGONAL  CHANNELS.
*
*********    COMPUTE  REDUCED MASSES.
*
      DO 35 NGAMM=1,NOCHNN
         REDMAS(NGAMM)=3.*QM/2.D0
*
********    DEFINE OSZILLATOR PARAMETER FOR REL. WAVEFUNCTION
*
      BPRIME(NGAMM)=DSQRT(1.5*QM*B**2/REDMAS(NGAMM))
*     BPRIME(NGAMM)=B
*
******   COMPUTE  REL. KIN. ENERGY, REL. MOMENTUM
*
         H2M1=(HBC**2)/(2*REDMAS(1))
         ER1=EREL(I33)
         RKE1=DSQRT(ER1/H2M1)
         RKSAVE(I33,1)=RKE1
         WRITE(6,'(////,10X,2A,3X,F7.3,A)')   ' FIRST CHANNEL : ',
     &      'RELATIVE KINETIC ENERGY IN C.M.SYSTEM =',ER1,' MEV'
         NOPEN=1
35    CONTINUE
****
****   STORE NUMBER OF OPEN CHANNELS FOR THAT ENERGY.
****
      NOPARR(I33)=NOPEN
*****   DEFINE  CONTROL PARAMETER  KEXC  FOR OPEN CHANNELS.
      DO 114 N=1,NOPEN
      KEXC(N)=KEXC(N)+1
114   CONTINUE
*****   DEFINE  TOTAL ENERGY
      RKEY(1)=RKEYY
      EINT(1,1)=EINT1
      ET(1)=ER1+EINT(1,1)-RKEY(1)
      IF (NOCHNN.EQ.2) THEN
         EINT(2,2)=EINT2
         ET(2)=ER1+EINT(2,2)-RKEY(1)
      ELSE IF (NOCHNN.EQ.3) THEN
         EINT(2,2)=EINT2
         EINT(3,3)=EINT3
         ET(2)=ER1+EINT(2,2)-RKEY(1)
         ET(3)=ER1+EINT(3,3)-RKEY(1)
      END IF
      WRITE(6,'(//,4X,A,D10.5)')   ' E_TOTAL IN  1. CHANNEL :',ET(1)
*******************************
*       COMPUTE  RK(...)      *
*******************************
      DO 703 NGAMM=1,NOCHNN
      DO 703 NDELT=1,NOCHNN
         ERN=ER1
         H2MN=H2M1
         BPRM2=BPRIME(NGAMM)**2
         RRLM=3./(4.*BPRM2)
      DO 703 I=1,NS
      DO 703 J=1,NS
      IF (NGAMM.EQ.NDELT.AND.NGAMM.LE.NOPEN) THEN
         RKEXEY=RKEX(I,J,LPW,ERN,BPRM2,H2MN,NGAMM)
         RK(NGAMM,NDELT,I,J)=BIGHAM(NGAMM,NDELT,I,J)
     &    -ET(NGAMM)*BIGNOR(NGAMM,NDELT,I,J) - RKEXEY
      ELSE
         RK(NGAMM,NDELT,I,J)=BIGHAM(NGAMM,NDELT,I,J) - ET(1)
     &     *BIGNOR(NGAMM,NDELT,I,J)
      END IF
*****  DECOUPLING.
      IF (ICOUPL.EQ.1) THEN
         CONTINUE
      ELSE IF (ICOUPL.EQ.0) THEN
         IF (NGAMM.NE.NDELT) THEN
            RK(NGAMM,NDELT,I,J)=0.D0
         END IF
      ELSE
         WRITE(6,*) 'WARNING: ICOUPL MUST BE 0 OR 1 !!!'
         WRITE(6,*) 'ICOUPL = ',ICOUPL
         WRITE(6,*) 'PROGRAM IS STOPPED !!!'
         STOP
      END IF
*****   END OF DECOUPLING.
703   CONTINUE
***********************************
*        P AND S CALCULATION      *
***********************************
      RRLM=3./(4.*BPRIME(1)**2)
      CALL MATCH3(NS,LPW,CPX,CSX,RKE1)
      DO 715 I=1,NS
      CP(1,I)=CPX(I)
      CS(1,I)=CSX(I)
715   CONTINUE
*
********************************************************
*       L-MATRIX CALCULATION  FOLLOWING                *
*   KAMIMURA, PROG. THEOR. PHYS. SUPP. 62 (1977) 236,  *
*              FORMULAE (4.123)                        *
********************************************************
      NS1=NS-1
      DO 708 NGAMM=1,NOCHNN
      DO 708 NDELT=1,NOCHNN
      IF (NGAMM.LE.NOPEN.AND.NDELT.LE.NOPEN) THEN
         DO 781 I=1,NS1
         DO 781 J=1,NS1
         CL(NGAMM,NDELT,I,J)=RK(NGAMM,NDELT,I,J)
     &      -(CP(NGAMM,NS)/CP(NGAMM,I))*RK(NGAMM,NDELT,NS,J)
     &      -(CP(NDELT,NS)/CP(NDELT,J))*RK(NGAMM,NDELT,I,NS)
     &      +(CP(NGAMM,NS)/CP(NGAMM,I))*(CP(NDELT,NS)/CP(NDELT,J))
     &      *RK(NGAMM,NDELT,NS,NS)
781      CONTINUE
      ELSE IF (NGAMM.LE.NOPEN.AND.NDELT.GT.NOPEN) THEN
         DO 782 I=1,NS1
         DO 782 J=1,NS
         CL(NGAMM,NDELT,I,J)=RK(NGAMM,NDELT,I,J)
     &         -(CP(NGAMM,NS)/CP(NGAMM,I))*RK(NGAMM,NDELT,NS,J)
         CL(NDELT,NGAMM,J,I)=RK(NDELT,NGAMM,J,I)
     &         -(CP(NGAMM,NS)/CP(NGAMM,I))*RK(NDELT,NGAMM,J,NS)
782      CONTINUE
      ELSE IF (NGAMM.GT.NOPEN.AND.NDELT.GT.NOPEN) THEN
         DO 783 I=1,NS
         DO 783 J=1,NS
         CL(NGAMM,NDELT,I,J)=RK(NGAMM,NDELT,I,J)
783      CONTINUE
      END IF
708   CONTINUE
**********************************
*       A MATRIX CALCULATION     *
**********************************
      CALL LOADIS(CL,ZASAV,NS1,NOCHNN,NOPEN)
********************************************************
*              M-MATRIX CALCULATION                    *
*         FOLLOWING  KAMIMURA  FORMULAE (4.123)        *
********************************************************
      DO 709 NGAMM=1,NOCHNN
      DO 709 INCID=1,NOPEN
      IF (NGAMM.LE.NOPEN) THEN
         DO 791 I=1,NS1
         ZM(NGAMM,I,INCID)= -CP(INCID,NS)*RK(NGAMM,INCID,I,NS)
     &    +(CP(NGAMM,NS)/CP(NGAMM,I))*CP(INCID,NS)*RK(NGAMM,INCID,NS,NS)
791      CONTINUE
      ELSE IF (NGAMM.GT.NOPEN) THEN
         DO 792 I=1,NS
         ZM(NGAMM,I,INCID)= -CP(INCID,NS)*RK(NGAMM,INCID,I,NS)
792      CONTINUE
      END IF
709   CONTINUE
********************************************
      IF (ITEST.EQ.I33) THEN
         CALL LTEST(RK,CL,ZM,ZASAV,CP,CS,NS,NOCHNN)
      END IF
******************************************
*      LOOP OVER ALL INCOMING WAVES      *
*        SOLVE  KAMIMURA (4.123)         *
******************************************
      DO 712 INCID=1,NOPEN
******************************************
      CALL SETZAA(ZASAV,ZAA,NS1,NOCHNN,NOPEN)
*
*     B MATRIX CALCULATION
*
      DO 711 NGAMM=1,NOCHNN
      IF (NGAMM.LE.NOPEN) THEN
         DO 111 I=1,NS1
         ZBB((NGAMM-1)*NS1+I)=ZM(NGAMM,I,INCID)
111      CONTINUE
      ELSE IF (NGAMM.GT.NOPEN) THEN
         DO 112 I=1,NS
         ZBB(NOPEN*NS1+(NGAMM-NOPEN-1)*NS+I)=ZM(NGAMM,I,INCID)
112      CONTINUE
      END IF
711    CONTINUE
*******************************************************************
*     SOLUTION OF THE LINEAR EQUATION SYSTEM  ZAA*ZPCSAV=ZBB      *
*******************************************************************
      MDIM=NOPEN*NS1+(NOCHNN-NOPEN)*NS
      MDIMAX=MAXCHN*MESH
      CALL F04ADF(ZAA,MDIMAX,ZBB,MDIMAX,MDIM,1,ZPCSAV,MDIMAX,WKSPC,IER)
*
      DO 713 N=1,NOPEN*NS1
      ZPC(N,INCID)=ZPCSAV(N)
713   CONTINUE
******************************************
712   CONTINUE
******************************************
      CALL DIVID(ZPC,CP,CC,NS1,NOPEN)
**************************
*       PHASE SHIFTS     *
*      KAMIMURA (4.86)   *
**************************
      CALL SMTRX(CC,CS,CSMTR,NS,NOPEN)
*
      WRITE(6,'(//,5X,A)')   ' S-MATRIX ELEMENTS    CSMTR(OUT,IN) '
      DO 742 M=1,NOPEN
      WRITE(6,'(/,3X,4(A,D11.5,3X,D11.5,A,5X))')
     &    ('(',CSMTR(M,N),')',N=1,NOPEN)
      DELTA(I33,M)=PHASE(CSMTR(M,M))/2.
742   CONTINUE
*
*   COMPUTE:    S-MATRIX  *  S-MATRIX ADJOINT
*
      WRITE(6,'(///,5X,A)')   ' UNITARITY OF S-MATRIX '
      DO 743 M=1,NOPEN
      DO 743 N=1,NOPEN
      CSSADJ(M,N)=0.
      DO 743 I=1,NOPEN
      CSSADJ(M,N)=CSMTR(I,M)*CONJG(CSMTR(I,N))+CSSADJ(M,N)
743   CONTINUE
      DO 744 M=1,NOPEN
      WRITE(6,'(/,3X,4(A,D11.5,3X,D11.5,A,5X))')
     &    ('(',CSSADJ(M,N),')',N=1,NOPEN)
744   CONTINUE
*
*****   COMPUTE DEVIATION OF S-MATRIX FROM UNITARITY   ****
*
      DEVIAT(I33)=0.
      DO 746 M=1,NOPEN
      DO 746 N=1,NOPEN
      IF (M.EQ.N) THEN
         ERR=ABS(1 - ABS(CSSADJ(M,N)))
      ELSE IF (M.NE.N) THEN
         ERR=ABS(CSSADJ(M,N))
      END IF
      DEVIAT(I33)=DEVIAT(I33)+ERR
746    CONTINUE
***
733   CONTINUE
************************************************
*         END OF LOOP OVER ENERGIES            *
************************************************
*
**************************************
*         OUTPUT  OF  SUMMARY        *
**************************************
*
      WRITE(6,'(////)')
      WRITE(6,'(40X,A,/,49X,A,I2,///)')
     &   'BARYON-BARYON PHASE SHIFT','LPW = ',LPW
      WRITE(6,'(52X,A)')  'DEVIATION'
      WRITE(6,'(7X,A,10X,A,10X,A,/)')  'ENERGY','PHASE-SHIFTS  DELTA'
     &   ,'FROM UNITY (PERCENT)'
      DO 724 I24=1,NE
         IF (NOPARR(I24).EQ.1)  TFRMT='(5X,D10.4,4X,D11.4,28X,F8.4)'
         IF (NOPARR(I24).EQ.2)  TFRMT='(5X,D10.4,4X,2D11.4,17X,F8.4)'
         IF (NOPARR(I24).EQ.3)  TFRMT='(5X,D10.4,4X,3D11.4,6X,F8.4)'
      WRITE(6,TFRMT)  EREL(I24),(DELTA(I24,J),J=1,NOPARR(I24))
     &                ,DEVIAT(I24)*100.
724   CONTINUE
*
999   CONTINUE
      STOP
      END
***************************************************
***************************************************
*                                                 *
*     SUBROUTINE SECTION: PART1 AND PART2         *
*                                                 *
***************************************************
***************************************************
      BLOCK DATA INTIAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      COMMON /Y1/ YALF(5),YC(5),NGAUSS
      COMMON /Y2/ RLA(5),RLC(5),RCUTOF
      COMMON /A3/ OMG(12),ETA(12),NDIM
      COMMON /A4/ XCOF(13),XEXP(13),LDIM
      DATA YALF/8.D0,2.D0,1.5D0,0.2D0,0.044D0/
      DATA YC/.13545D0,-.30767D0,.51367D0,.36939D0,.24647D0/
      DATA NGAUSS/5/
      DATA RLA/70.0D0,20.0D0,8.0D0,2.5D0,0.800D0/
      DATA RLC/0.18904D0,0.12731D0,0.25348D0,0.25526D0,0.0667947D0/
      DATA RCUTOF/4.2D0/
      DATA  OMG/.12793D0,.12583D0,.12167D0,.11550D0,.10744D0,
     &          .97610D-1,.86190D-1,.73340D-1,.59290D-1,.44270D-1,
     &          .28530D-1,.12340D-1/
      DATA  ETA/.64050D-1,.19111D0,.31504D0,.43379D0,.54542D0,
     &          .64809D0,.74012D0,.82000D0,.88641D0,.93827D0,
     &          .97472D0,.99518D0/
      DATA NDIM/12/
      DATA LDIM/13/
      END
*
*
      SUBROUTINE SFCME(SPIN,STCD,STCE)
********  DEFINE COEFFICIENTS IN SFC-SPACE
********  NO PION EXCHANGE IS INCLUDED.
********  LPW=0 IS ASSUMED.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISP0(6,16), ISP1(6,16)
      DIMENSION STCD(6,5),STCE(6,11)
*
      DATA (ISP0(K, 1),K=1,6)/ 243,   0, 243,   0,   0, 243/
      DATA (ISP0(K, 2),K=1,6)/-648,   0,-648,   0,   0,-162/
      DATA (ISP0(K, 3),K=1,6)/   0,   0,   0,   0,   0,-324/
      DATA (ISP0(K, 4),K=1,6)/ 648,   0,-648,   0,   0,-180/
      DATA (ISP0(K, 5),K=1,6)/   0,   0,   0,-108, -72, -30/
      DATA (ISP0(K, 6),K=1,6)/  -3,  12,   3, -36,  72,-189/
      DATA (ISP0(K, 7),K=1,6)/   8, -32,  -8,  96,-192,  18/
      DATA (ISP0(K, 8),K=1,6)/ -16,  64,  16,  24, -48, 288/
      DATA (ISP0(K, 9),K=1,6)/   8, -32,  -8,  96,-192, 180/
      DATA (ISP0(K,10),K=1,6)/   8, -32,  -8, -12,  24, 180/
      DATA (ISP0(K,11),K=1,6)/  -4,  16,   4,   6, -12, 315/
      DATA (ISP0(K,12),K=1,6)/-136, -32,  -8, -96,-192, 196/
      DATA (ISP0(K,13),K=1,6)/ 496, -64,-112,   4,  88, 122/
      DATA (ISP0(K,14),K=1,6)/  56,  64,  -8, -96,-192, 112/
      DATA (ISP0(K,15),K=1,6)/  56, -32,  -8,  88,  16, 112/
      DATA (ISP0(K,16),K=1,6)/   0, -24, -12, 100,  88,  42/
*
      DATA (ISP1(K, 1),K=1,6)/ 243,   0, 243,   0,   0, 243/
      DATA (ISP1(K, 2),K=1,6)/-648,   0,-648,   0,   0,-162/
      DATA (ISP1(K, 3),K=1,6)/   0,   0,   0,   0,   0,-324/
      DATA (ISP1(K, 4),K=1,6)/ 648,   0,-648,   0,   0, -72/
      DATA (ISP1(K, 5),K=1,6)/   0,   0,   0, -84,-120,   6/
      DATA (ISP1(K, 6),K=1,6)/  -3,  12,   3, -36,  72,-189/
      DATA (ISP1(K, 7),K=1,6)/   8, -32,  -8,  96,-192,  18/
      DATA (ISP1(K, 8),K=1,6)/ -16,  64,  16,  24, -48, 288/
      DATA (ISP1(K, 9),K=1,6)/   8, -32,  -8,  96,-192, 180/
      DATA (ISP1(K,10),K=1,6)/   8, -32,  -8, -12,  24, 180/
      DATA (ISP1(K,11),K=1,6)/  -4,  16,   4,   6, -12, 315/
      DATA (ISP1(K,12),K=1,6)/-136, -32,  -8, -96,-192,  88/
      DATA (ISP1(K,13),K=1,6)/ 304, -64,-304,  12,  72,  38/
      DATA (ISP1(K,14),K=1,6)/  56,  64,  -8, -96,-192,  40/
      DATA (ISP1(K,15),K=1,6)/  56, -32,  -8,  64,  64,  40/
      DATA (ISP1(K,16),K=1,6)/   8,  -8,   4,  80, 128,   0/
*
      IF (ABS(SPIN-0).LE.1D-4) THEN
      DO 100 K=1,6
      IF ((K.EQ.1) .OR. (K.EQ.3) .OR. (K.EQ.5) .OR. (K.EQ.6)) THEN
         FACTOR=1.D0 / 243.D0
      ELSE IF ((K.EQ.2) .OR. (K.EQ.4)) THEN
         FACTOR=SQRT(5.D0) / 243.D0
      END IF
      DO 110 I=1,16
         IF (I.LE.5) THEN
            STCD(K,I)=ISP0(K,I)*FACTOR
         ELSE IF(I.GT.5) THEN
            STCE(K,I-5)=ISP0(K,I)*FACTOR
         END IF
110   CONTINUE
100   CONTINUE
      END IF
*
      IF (ABS(SPIN-1).LE.1D-4) THEN
      DO 120 K=1,6
      IF ((K.EQ.1) .OR. (K.EQ.3) .OR. (K.EQ.5) .OR. (K.EQ.6)) THEN
         FACTOR=1.D0 / 243.D0
      ELSE IF ((K.EQ.2) .OR. (K.EQ.4)) THEN
         FACTOR=SQRT(5.D0) / 243.D0
      END IF
      DO 130 I=1,16
         IF (I.LE.5) THEN
            STCD(K,I)=ISP1(K,I)*FACTOR
         ELSE IF(I.GT.5) THEN
            STCE(K,I-5)=ISP1(K,I)*FACTOR
         END IF
130   CONTINUE
120   CONTINUE
      END IF
      RETURN
      END
*
*
      SUBROUTINE COFP(IST,RCSTD,RCSTE,STCD,STCE)
********  DEFINE COEFFICIENTS IN SFC-SPACE FOR ALL INTERACTIONS.
********  NOCHNN=1 IS ASSUMED.
********  NO RESTRICTION FOR LPW.
********  RCSTD, RCSTE:  M.E. FOR PION  EXCHANGE PART
********  STCD,  STCE:   M.E. FOR GLUON EXCHANGE PART
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RCSTD(1),RCSTE(1),C10(7),C00(7),C11(7)
      DIMENSION STCD(6,5),STCE(6,11),CG00(16),CG11(16)
      DATA CG00/729.0D0,-1944.0D0,0.0D0,1944.0D0,0.0D0,
     &189.0D0,-504.0D0,1008.0D0,-504.0D0,-504.0D0,252.0D0,
     &792.0D0,-144.0D0,360.0D0,360.0D0,18.0D0/
      DATA CG11/729.0D0,-1944.0D0,0.0D0,1944.0D0,0.0D0,
     &93.0D0,-248.0D0,496.0D0,-248.0D0,-248.0D0,124.0D0
     &,472.0D0,944.0D0,136.0D0,136.0D0,42.0D0/
      DATA C10/5.0D0,-0.925926D0,-1.37037D0,2.185185D0,
     &0.407407D0,0.407407D0,0.111111D0/
      DATA C00/5.0D0,2.777778D0,4.777778D0,3.444444D0,
     &3.444444D0,3.444444D0,2.333333D0/
      DATA C11/5.0D0,0.308642D0,2.604938D0,2.160500D0,
     &1.567901D0,1.567901D0,0.901234D0/
      IST1=IST+1
      GO TO (1,2,3),IST1
 1    DO 10 I=1,7
      IF(I.LE.2) RCSTD(I)=C00(I)
      IF(I.GT.2) RCSTE(I-2)=C00(I)
 10   CONTINUE
      DO 11 K=1,16
      IF(K.LE.5) STCD(1,K)=CG00(K)/729.
      IF(K.GT.5) STCE(1,K-5)=CG00(K)/729.
 11   CONTINUE
      GO TO 4
 2    DO 20 I=1,7
      IF(I.LE.2) RCSTD(I)=C10(I)
      IF(I.GT.2) RCSTE(I-2)=C10(I)
 20   CONTINUE
      GO TO 4
 3    DO 30 I=1,7
      IF(I.LE.2) RCSTD(I)=C11(I)
      IF(I.GT.2) RCSTE(I-2)=C11(I)
 30   CONTINUE
      DO 31 K=1,16
      IF(K.LE.5) STCD(1,K)=CG11(K)/729.
      IF(K.GT.5) STCE(1,K-5)=CG11(K)/729.
 31   CONTINUE
 4    RETURN
      END
*
*
      SUBROUTINE BSL(LMAX,NCH)
****************************************************************
*     COMPUTING BESSEL FUNCTION BY RECURSION RELATION          *
*     FOR THREE DIFFERENT ARGUMENTS                            *
****************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /B/ A(3),BS(5,3,6)
      DO 2 I=1,3
      X=A(I)
      BS(1,I,NCH)=DSINH(X)/X
      BS(2,I,NCH)=-(BS(1,I,NCH)/X)+DCOSH(X)/X
      IF(LMAX.EQ.1) GO TO 2
      DO 1 LK=2,LMAX
   1  BS(LK+1,I,NCH)=-((2*LK-1)*BS(LK,I,NCH)/X)+BS(LK-1,I,NCH)
  2   CONTINUE
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION RNORM(I,J,K,NCH)
*  NORM KERNEL
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (Mesh=10)
      COMMON /A1/ RLM(MESH,MESH,6)
      COMMON /B/ A(3),BS(5,3,6)
      COMMON /L/ LPW
      COMMON /WIDTH/ B
      L1=LPW+1
*   K=0 DIREKT, K=1 EXCH
      EX=-3.*RLM(I,J,NCH)/8.
      RNORM=BS(L1,K+1,NCH)*DEXP(EX)
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION TKIN(I,J,K,RK0,NCH)
*  KERNEL OF THE KINETIC ENERGY
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (Mesh=10)
      COMMON /A1/ RLM(MESH,MESH,6)
      COMMON /A2/ RGM(MESH,MESH,6)
      COMMON /B/ A(3),BS(5,3,6)
      COMMON /L/ LPW
      COMMON /WIDTH/ B
      N=2+4*K
      EX=-3*RLM(I,J,NCH)/8
      BSX=0D0
      IF(LPW.GT.0) BSX=BS(LPW,K+1,NCH)
      TKIN=(5-RLM(I,J,NCH)/4.)*BS(LPW+1,K+1,NCH)
      TK=(LPW+1)*BS(LPW+2,K+1,NCH)/(2*LPW+1)
      TK=RGM(I,J,NCH)*(TK+LPW*BSX/(2*LPW+1))/N
      TKIN=RK0*(TKIN+TK)*DEXP(EX)
      RETURN
      END
*
*
      SUBROUTINE DCHAN(NCH,N1,N2)
      N1=3
      IF(NCH.LE.3) N1=1
      IF(NCH.GT.3.AND.NCH.LE.5) N1=2
      N2=NCH-N1
      IF(NCH.LE.3) N2=N2+1
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION VCONF(I,J,K,L,C1,C2,NCH,ICONF)
*  CALLING CONFINEMENT POTENTIALS
      IMPLICIT REAL*8(A-H,O-Z)
      GO TO (1,2,1),ICONF
1      VCONF=VCFL(I,J,K,L,C1,C2,NCH)
      RETURN
2      VCONF=VCFQ(I,J,K,L,C1,C2,NCH)
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION VCFL(I,J,K,NL,X1,X2,NCH)
*   LINEAR CONFINEMENT POTENTIAL
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT REAL*8(L)
      PARAMETER (Mesh=10)
      COMMON /A1/ RLM(MESH,MESH,6)
      COMMON /A2/ RGM(MESH,MESH,6)
      COMMON /A4/ LCOF(13),LEXP(13),NLDIM
      COMMON /L/ NLPW
      COMMON /WIDTH/ B
      CALL DCHAN(NCH,N1,N2)
      IF(N1.NE.1) N1=2*N1
      IF(N2.NE.1) N2=2*N2
      PHI=4.*DATAN(1D0)
      CVL=0D0
      ARG=3*RGM(I,J,NCH)/4
      IF(K.EQ.1) ARG=RGM(I,J,NCH)/4
      EX0=-3*RLM(I,J,NCH)/8
      F1=RLM(I,J,NCH)/4
      F2=RGM(I,J,NCH)/2
      GO TO (1,2,3,3,4),NL
1     CVL=4*DEXP(EX0)*BSS(ARG,NLPW)/DSQRT(PHI)
      GO TO 8
2     DO 5 K5=1,NLDIM
      EX=EX0-LEXP(K5)*RLM(I,J,NCH)/8
      ARG1=ARG-((-1)**K)*LEXP(K5)*RGM(I,J,NCH)/4
      CVL=CVL+LCOF(K5)*DEXP(EX)*BSS(ARG1,NLPW)
      IF(K5.EQ.1) GO TO 5
      CVL=CVL+F1*DEXP(EX)*BSS(ARG1,NLPW)*LCOF(K5)
      CVL=CVL+ (-1)**K *F2*DEXP(EX)*BSS(ARG1,NLPW+1)*LCOF(K5)
5     CONTINUE
      GO TO 8
3     DO 6 K6=1,NLDIM
      CEX=EX0-LEXP(K6)*RGM(J,J,N2)/8
      IF(NL.EQ.4) CEX=EX0-LEXP(K6)*RGM(I,I,N1)/8
      F3=RGM(J,J,N2)/4
      IF(NL.EQ.4) F3=RGM(I,I,N1)/4
      ARG1=ARG
      CVL=CVL+LCOF(K6)*DEXP(CEX)*BSS(ARG1,NLPW)
      IF(K6.EQ.1) GO TO 6
      CVL=CVL+F3*LCOF(K6)*DEXP(CEX)*BSS(ARG1,NLPW)
6     CONTINUE
      GO TO 8
4     DO 7 K7=1,NLDIM
      EX=EX0-LEXP(K7)*RLM(I,J,NCH)/8
      ARG1=ARG-LEXP(K7)*RGM(I,J,NCH)/4
      CVL=CVL+LCOF(K7)*DEXP(EX)*BSS(ARG1,NLPW)
      IF(K7.EQ.1) GO TO 7
      CVL=CVL+F1*LCOF(K7)*DEXP(EX)*BSS(ARG1,NLPW)
      CVL=CVL+F2*LCOF(K7)*DEXP(EX)*BSS(ARG1,NLPW+1)
7     CONTINUE
8     VCFL=-CVL*DSQRT(X1/2)*X2
      RETURN
      END
*
*
      FUNCTION VCFQ(I,J,K,L,C1,C2,NCH)
*  QUADRATIC CONFINEMENT POTENTIAL
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (Mesh=10)
      COMMON /A1/ RLM(MESH,MESH,6)
      COMMON /A2/ RGM(MESH,MESH,6)
      COMMON /B/ A(3),BS(5,3,6)
      COMMON /L/ LPW
      COMMON /WIDTH/ B
      CALL DCHAN(NCH,N1,N2)
      IF(N1.NE.1) N1=2*N1
      IF(N2.NE.1) N2=2*N2
      EX=DEXP(-3.*RLM(I,J,NCH)/8.)
      F1=RGM(I,J,NCH)/6.
      F2=RLM(I,J,NCH)/12
      GO TO (3,1,2,2,4),L
1     IF(K.EQ.1) F1=-F1
      GO TO 4
2     F1=0.
      F2=RLM(J,J,N2)/24
      IF(L.EQ.4) F2=RLM(I,I,N1)/24
      GO TO 4
3     VCFQ=-3*C1*C2*BS(LPW+1,K+1,NCH)*EX
      GO TO 5
4     VCFQ=-3*C1*C2*(1+F2)*BS(LPW+1,K+1,NCH)*EX
      VCFQ=VCFQ-3*C1*C2*F1*BS(LPW+2,K+1,NCH)*EX
5     RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION VC(I,J,K,L,X1,X2,NCH)
*  COULOMB POTENTIAL
      PARAMETER (Mesh=10)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /A1/ RLM(MESH,MESH,6)
      COMMON /A2/ RGM(MESH,MESH,6)
      COMMON /A3/ OMG(12),ETA(12),NDIM
      COMMON /B/ A(3),BS(5,3,6)
      COMMON /L/ LPW
      COMMON /WIDTH/ B
      CALL DCHAN(NCH,N1,N2)
      IF(N1.NE.1) N1=2*N1
      IF(N2.NE.1) N2=2*N2
      PHI=4.*DATAN(1D0)
      M=2
      IF(K.EQ.0) M=1
      GO TO (1,4,2,2,4),L
1     EX=-3*RLM(I,J,NCH)/8
      F1=2*DEXP(EX)*BS(LPW+1,M,NCH)/DSQRT(PHI)
      GO TO 6
2     F1=0.
      EX=-3*RLM(I,J,NCH)/8.
      R=DSQRT(RLM(J,J,N2)/16.)
      IF(L.EQ.4) R=DSQRT(RLM(I,I,N1)/16.)
      F1=S15AEF(R,IFAIL)*DEXP(EX)/R
      F1=F1*BS(LPW+1,M,NCH)
      GO TO 6
4     F1=0.
      DO 5 IK=1,NDIM
      EX=-(3+ETA(IK))*RLM(I,J,NCH)/8.
      ARG=3-ETA(IK)
      IF(K.EQ.1) ARG=4-ARG
      IF(L.EQ.5) ARG=2-ARG
      ARG=ARG*RGM(I,J,NCH)/4
5     F1=F1+OMG(IK)*DEXP(EX)*BSS(ARG,LPW)
6     VC=X1*F1/DSQRT(2*X2)
*     IF(I.EQ.1) WRITE(6,100) I,J,K,VC
  100 FORMAT(1H ,3I3,F10.5)
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION VD(I,J,K,L,VD0,NCH)
*  DELTA POTENTIAL
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (Mesh=10)
      COMMON /A1/ RLM(MESH,MESH,6)
      COMMON /B/ A(3),BS(5,3,6)
      COMMON /L/ LPW
      COMMON /WIDTH/ B
      L1=LPW+1
      CALL DCHAN(NCH,N1,N2)
      IF(N1.NE.1) N1=2*N1
      IF(N2.NE.1) N2=2*N2
      M=2
      IF(K.EQ.0) M=1
      EX=0.
      GO TO (1,2,3,4,5),L
1     F=.375
      GO TO 6
2     F=.5
      M=3
      GO TO 6
3     F=.375
      EX=RLM(J,J,N2)/16.
      GO TO 6
4     F=.375
      EX=RLM(I,I,N1)/16.
      GO TO 6
5     F=.5
6     EX=EX+F*RLM(I,J,NCH)
      VD=VD0*DEXP(-EX)
      IF(L.EQ.5.AND.LPW.NE.0) VD=0.D0
      IF(L.EQ.5) GO TO 7
      VD=VD*BS(L1,M,NCH)
7     RETURN
      END
*
*
*
*
      SUBROUTINE VOBEP2(IOBEP,KOBEP,HOBEP,TX,SY,SPIN)
      PARAMETER (NDATA=11,NDATA4=15,LWRK=82)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RKNOT(NDATA4),RCOEF(NDATA4),RWORK(LWRK)
      DIMENSION TX(KOBEP),SY(KOBEP)
      DIMENSION ZR(NDATA),ZPIO(NDATA),ZSIG(NDATA)
      DIMENSION ZOMG(NDATA),ZRHO(NDATA),ZV(NDATA)
      DIMENSION ZPIO0(NDATA),ZSIG0(NDATA),ZOMG0(NDATA),ZRHO0(NDATA)
      DIMENSION ZPIO1(NDATA),ZSIG1(NDATA),ZOMG1(NDATA),ZRHO1(NDATA)
      COMMON /L/ LPW
      COMMON /WIDTH/ B
      COMMON /P/ PARP,PARS,PARW,PARR
************** OBE POTENTIALS FOR SPIN=0  **************
      DATA ZPIO0/ 164.407,  89.697,  26.100, -3.0270, -1.9938,
     &   -1.1240, -.66009, -.39871, -.24584, -.15399, -.097665/
      DATA ZSIG0/-373.121,-214.546,-119.468,-42.0311, -8.4209,
     &   -1.7996, -.40062, -.09173, -.02144, -.00509, -.001224/
      DATA ZOMG0/ 384.360,  87.177,  45.913,  21.485,  1.3547,
     &    .14963,  .01721, .002038, .000246, .000030,  .000004/
      DATA ZRHO0/  54.805, 13.0636,  5.5117, .924880, .007900,
     &   .001045, .000144, .000020, .000003, .000000,  .000000/
************** OBE POTENTIALS FOR SPIN=1  **************
      DATA ZPIO1/ 164.407,  89.697,  26.100, -3.0270, -1.9938,
     &   -1.1240, -.66009, -.39871, -.24584, -.15399, -.097665/
      DATA ZSIG1/-373.121,-214.546,-119.468,-42.0311, -8.4209,
     &   -1.7996, -.40062, -.09173, -.02144, -.00509, -.001224/
      DATA ZOMG1/ 12.4658,  80.534,  44.049,  14.525,  2.2037,
     &    .24341,  .02800, .003314, .000400, .000049,  .000006/
      DATA ZRHO1/-117.120,-16.2727,-11.1071,-2.21368,-.060261,
     &  -.007969,-.001098,-.000156,-.000022,-.000003,  .000000/
**************************************************************
      DATA GCOUPP/14.2/, GCOUPS/5.66/, GCOUPW/4.5/, GCOUPR/0.5/
      DATA ZR/0.04999,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0001/
      HBC=197.32891
      IF (IOBEP.EQ.0) THEN
         PARP=0.D0
         PARS=0.D0
         PARW=0.D0
         PARR=0.D0
      ELSE IF (IOBEP.EQ.1) THEN
         CONTINUE
      ELSE IF (IOBEP.EQ.2) THEN
         PARP=0.D0
         PARW=0.D0
         PARR=0.D0
      END IF
**************** CHOOSE POTENTIAL ACCORDING TO SPIN  ****
      IF (ABS(SPIN).LT.1D-4) THEN
         DO 6 I=1,NDATA
         ZPIO(I)=ZPIO0(I)
         ZSIG(I)=ZSIG0(I)
         ZOMG(I)=ZOMG0(I)
         ZRHO(I)=ZRHO0(I)
6        CONTINUE
      ELSE IF (ABS(SPIN-1.D0).LT.1D-4) THEN
         DO 7 I=1,NDATA
         ZPIO(I)=ZPIO1(I)
         ZSIG(I)=ZSIG1(I)
         ZOMG(I)=ZOMG1(I)
         ZRHO(I)=ZRHO1(I)
7        CONTINUE
      END IF
*
      DO 9  I=1,NDATA
      ZV(I)=PARP*ZPIO(I)+PARS*ZSIG(I)+PARW*ZOMG(I)+PARR*ZRHO(I)
*         WRITE(6,*) 'X, EXACT F(X):',ZR(I),ZV(I)
9     CONTINUE
      WRITE(6,'(/6X,A,F5.2)') 'PION-PARAMETER  IS',PARP
      WRITE(6,'(6X,A,F5.2)')  'SIGMA-PARAMETER IS',PARS
      WRITE(6,'(6X,A,F5.2)')  'OMEGA-PARAMETER IS',PARW
      WRITE(6,'(6X,A,F5.2)')  'RHO-PARAMETER   IS',PARR
      write(6,'(/6x,2a,4f5.2)')
     & 'SCALED COUPLING CONSTANTS:',
     & '  PARP*GP, PARS*GS, PARW*GW, PARR*GR :',
     &     PARP*gcoupp, PARS*gcoups, PARW*gcoupw, PARR*gcoupr
*****************  CUBIC SPLINE INTERPOLATION  ***************
*****************     N E W   V E R S I O N    ***************
      IFAIL=0
      CALL  E01BAF(NDATA,ZR,ZV,RKNOT,RCOEF,NDATA4,RWORK,LWRK,IFAIL)
*
      DO 10 I=1,KOBEP
*******  DEFINE THE NEW FINER GRID = TX(I)   **********
      TX(I)=HOBEP*I
*******  COMPUTE CUBIC POLYNOMIAL FOR FINER GRID    *****
      XVAL=TX(I)
      CALL E02BBF(NDATA4,RKNOT,RCOEF,XVAL,SVAL,IFAIL)
      SY(I)=SVAL
*        WRITE(6,*) 'X, APPROX. F(X):',TX(I),SY(I)
10    CONTINUE
 101  FORMAT(10F12.6)
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION BSS(X,L)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BIN(5)
      IF(X.NE.0D0) GO TO 3
      IF (L.EQ.0) BSS=1.D0
      IF (L.NE.0) BSS=0.D0
      GO TO 1
  3   BIN(1)=DSINH(X)/X
      BSS=BIN(1)
      IF(L.EQ.0) GO TO 1
      BIN(2)=-(BIN(1)/X)+DCOSH(X)/X
      BSS=BIN(2)
      IF(L.EQ.1) GO TO 1
      DO 2 LK=2,L
      BIN(LK+1)=-((2*LK-1)*BIN(LK)/X)+BIN(LK-1)
   2  BSS=BIN(LK+1)
   1  RETURN
      END
*
*
      SUBROUTINE KEROUT(IS,ISP,NC)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (Mesh=10)
      COMMON /KERNEL/ SW(MESH,3),RN(6),RKIN(6),RCF(6),RCOU(6),RD(6),
     &                RCMI(6),RPION,OBEP,RHAM(6)
      DATA IFIRST/0/
      IF(IFIRST.EQ.0) WRITE(6,100)
      IFIRST=1
      WRITE(6,*) IS,ISP
      WRITE(6,110) 'SIJ  ', (SW(IS,K),SW(ISP,K),K=1,3)
      WRITE(6,110) 'NORM ', (RN(K) ,K=1,NC)
      WRITE(6,110) 'KINE ', (RKIN(K), K=1,NC)
      WRITE(6,110) 'CONF ', (RCF(K),K=1,NC)
      WRITE(6,110) 'COUL ', (RCOU(K), K=1,NC)
      WRITE(6,110) 'DELTA', (RD(K), K=1,NC)
      WRITE(6,110) 'CMI  ', (RCMI(K), K=1,NC)
      WRITE(6,110) 'P+QEX', RPION
      WRITE(6,110) 'OBEP ', OBEP
      WRITE(6,110) 'HAM  ', (RHAM(K),K=1,NC)
      WRITE(6,'(//)')
  100 FORMAT(///10X,'RGM INTEGRAL KERNELS',//)
  110 FORMAT(5X,A,6D15.5)
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION VY2(I,J,K,L,X1,NCH,KL)
*  YUKAWA POTENTIAL FOR PION PART  AND/OR CUTOFF PART
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION YC(5),YALF(5)
      PARAMETER (Mesh=10)
      COMMON /A1/ RLM(MESH,MESH,6)
      COMMON /A2/ RGM(MESH,MESH,6)
      COMMON /A3/ OMG(12),ETA(12),NDIM
      COMMON /B/ A(3),BS(5,3,6)
      COMMON /Y1/ YA(5),YCC(5),NGAUSS
      COMMON /Y2/ RLA(5),RLC(5),RCUTOF
      COMMON /L/ LPW
      COMMON /WIDTH/ B
      CALL DCHAN(NCH,N1,N2)
      IF(N1.NE.1) N1=2*N1
      IF(N2.NE.1) N2=2*N2
      PHI=4.*DATAN(1.D0)
      GO TO(7,8),KL
  7   DO 9 K9=1,NGAUSS
      YALF(K9)=YA(K9) *(B/0.5)**2
  9   YC(K9)=YCC(K9)
      GO TO 10
  8   DO 11 K11=1,NGAUSS
      YALF(K11)=RLA(K11) *(B/0.5)**2
  11  YC(K11)=RLC(K11)
  10  VY=0.D0
      IL=LPW+1
      M=2
      IF(K.EQ.0) M=1
      DO 6 IG=1,NGAUSS
      YNUM=1+YALF(IG)
      YCA=YC(IG)/YNUM
      GO TO (1,4,2,2,4),L
 1    EX=-3.*RLM(I,J,NCH)/8.
      F1=2*DEXP(EX)*BS(IL,M,NCH)/DSQRT(PHI)
*     IF(K.EQ.0) WRITE(6,*) F1,X1,YCA
      GO TO 6
 2    F1=0.
      DO 3 IK=1,NDIM
      YNU=(YALF(IG)+ETA(IK))/YNUM
      EX=-3.*RLM(I,J,NCH)/8.
      R=YNU*RLM(J,J,N2)/16
      IF(L.EQ.4) R=YNU*RLM(I,I,N1)/16.
      EX=EX-R
 3    F1=F1+OMG(IK)*DEXP(EX)
      F1=F1*BS(IL,M,NCH)
      GO TO 6
 4    F1=0.
      DO 5 IK=1,NDIM
      YNU=(YALF(IG)+ETA(IK))/YNUM
      EX=-(3+YNU)*RLM(I,J,NCH)/8.
      ARG=3-YNU
      IF(K.EQ.1) ARG=4-ARG
      IF(L.EQ.5) ARG=2-ARG
      ARG=ARG*RGM(I,J,NCH)/4.
5     F1=F1+OMG(IK)*DEXP(EX)*BSS(ARG,LPW)
6     VY=VY+YCA*X1*F1
      IF(KL.EQ.2) VY2=-VY
      IF(KL.EQ.1) VY2=VY
*     IF(I.EQ.1) WRITE(6,100) I,J,K,VY
  100 FORMAT(1H ,3I3,F10.5)
      RETURN
      END
*
*
      FUNCTION RKEX(I,J,L,ER,BPRM2,H2M,NGAMM)
******    IMPROVED VERSION   ******
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      IMPLICIT COMPLEX*16(C)
      PARAMETER (MESH=10,MAXCHN=3)
      DIMENSION RK1(MESH,MESH,4,MAXCHN),RK2(MESH,MESH,4,MAXCHN)
************  SAVE STATEMENT HAS BEEN ADDED TO PROGRAM
************                      U.STRAUB   17.12.87
      SAVE RK1,RK2
      DIMENSION KEXC(MAXCHN)
      COMMON /A/ RRLM,SWW(MESH),RC,RKE
      COMMON/KEXCO/ KEXC
*
      IF(KEXC(NGAMM)-1) 5,5,10
    5 P4=16.*ATAN(1.)
      RLAM=3./(4.*BPRM2)
      FAK=((8*RLAM/P4)**1.5)*P4
      SI=SWW(I)
      SJ=SWW(J)
      ARG=2*SI*SJ-SI*SI-SJ*SJ-16.D0*BPRM2/3.D0*DLOG(1.D-8)
      IF (ARG.LT.0.D0) THEN
         RK1(I,J,L+1,NGAMM)=0.D0
         RK2(I,J,L+1,NGAMM)=0.D0
         GOTO 10
      END IF
      RMAX=0.5*(SI+SJ+DSQRT(ARG))
      IF (RMAX.LE.RC) THEN
         RK1(I,J,L+1,NGAMM)=0.D0
         RK2(I,J,L+1,NGAMM)=0.D0
         GOTO 10
      END IF
      N=20*(RMAX-RC)
      IF (MOD(N,2).EQ.1) N=N+1
      IF (N.EQ.0)        N=2
      RK1(I,J,L+1,NGAMM)=FAK*GQ(I,J,L,1,N,RC,RMAX)
      RK2(I,J,L+1,NGAMM)=FAK*GQ(I,J,L,2,N,RC,RMAX)
   10 RKEX=H2M*RK1(I,J,L+1,NGAMM)-ER*RK2(I,J,L+1,NGAMM)
      RETURN
      END
*
*
      FUNCTION PHASE(CX)
*******    CALCULATES PHASE OF COMPLEX NUMBER  CX .  *****
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      IMPLICIT COMPLEX*16(C)
      Y=DIMAG(CX)
      X=REAL(CX)
      PHASE=DATAN2(Y,X)
      RETURN
      END
*
*
      FUNCTION F2(X,I,L)
      IMPLICIT REAL*8(A-B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)
      PARAMETER (MESH=10)
      COMMON /A/ RLAM,SWW(MESH),RC,RKE
      RSAVE=RC
      RC=X
      S=SWW(I)
      Z=2*RLAM*S*X
      F2=X*REAL(CBESS(Z,L,1))
      RC=RSAVE
      RETURN
      END
*
*
      FUNCTION GQ(I,J,L,K,N,A,RMAX)
********    IMPROVED VERSION    ******
      IMPLICIT REAL*8(A-H,O-Z)
      STEP=(RMAX-A)/DFLOAT(N)
      X1=A
      X2=A-STEP
      N2=(N/2)-1
      IF(K.EQ.1) GQ=F2(X1,I,L)*DLR(X1,J,L)
      IF(K.EQ.2) GQ=F2(X1,I,L)*F2(X1,J,L)
      DO 1 M=1,N2
      X1=X1+2*STEP
      X2=X2+2*STEP
      IF(K.EQ.1) GQ=GQ+2.0*F2(X1,I,L)*DLR(X1,J,L)
     1   +4.0*F2(X2,I,L)*DLR(X2,J,L)
      IF(K.EQ.2) GQ=GQ+2.0*F2(X1,I,L)*F2(X1,J,L)
     1   +4.0*F2(X2,I,L)*F2(X2,J,L)
1     CONTINUE
      GQ=STEP*GQ/3.0D0
      RETURN
      END
*
*
      SUBROUTINE DIVID(ZPC,CP,CC,NS1,NOPEN)
*******    CALCULATES VARIATIONAL PARAMETERS CC  FOR   *****
*******    OPEN CHANNELS AS    CC = ZPC / CP .         *****
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      IMPLICIT COMPLEX*16(C)
      IMPLICIT COMPLEX*16(Z)
      PARAMETER (MESH=10,MAXCHN=3)
      DIMENSION ZPC(MAXCHN*MESH,MAXCHN),CC(MAXCHN,MESH,MAXCHN)
      DIMENSION CP(MAXCHN,MESH),CSUM(MAXCHN,MAXCHN)
      NS=NS1+1
      DO 1 NGAMM=1,NOPEN
      DO 1 INCID=1,NOPEN
      CSUM(NGAMM,INCID)=0.
      DO 2 I=1,NS1
      CC(NGAMM,I,INCID)=ZPC((NGAMM-1)*NS1+I,INCID)/CP(NGAMM,I)
      CSUM(NGAMM,INCID)=CSUM(NGAMM,INCID)+CC(NGAMM,I,INCID)
2     CONTINUE
      IF (NGAMM.EQ.INCID) THEN
         CC(NGAMM,NS,INCID)=1. - CSUM(NGAMM,INCID)
      ELSE IF (NGAMM.NE.INCID) THEN
         CC(NGAMM,NS,INCID)= -CSUM(NGAMM,INCID)
      END IF
1     CONTINUE
      RETURN
      END
*
*
      FUNCTION RCTE(K)
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      PARAMETER (MESH=10)
      COMMON /A/ RRLM,SWW(MESH),RC,RK
      PHI=4.*ATAN(1.)+0.*K
      RCTE=4*PHI*((2*RRLM/PHI)**0.75)
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION XSINH(RLM,R,S)
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      EX1=-RLM*((R-S)**2)
      VALUE1=DEXP(EX1)
      EX2=-RLM*((R+S)**2)
      IF (EX2.LE.-100) THEN
         VALUE2=0.D0
      ELSE
         VALUE2=DEXP(EX2)
      END IF
      XSINH=(VALUE1 - VALUE2)/2.
      RETURN
      END
*
*
      DOUBLE PRECISION FUNCTION XCOSH(RLM,R,S)
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      EX1=-RLM*((R-S)**2)
      VALUE1=DEXP(EX1)
      EX2=-RLM*((R+S)**2)
      IF (EX2.LE.-100) THEN
         VALUE2=0.D0
      ELSE
         VALUE2=DEXP(EX2)
      END IF
      XCOSH=(VALUE1 + VALUE2)/2.
      RETURN
      END
*
*
      COMPLEX FUNCTION CBESS*16(Z,L,ITYP)
      IMPLICIT REAL*8(A-B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)
      PARAMETER (MESH=10)
      DIMENSION CF(10)
      COMMON /A/ RRLM,SWW(MESH),RC,RK
*  INITIAL VALUES
      S=Z/(2.*RRLM*RC)
      K=1
      GO TO (1,2,2),ITYP
1     X1=XSINH(RRLM,RC,S)/Z
      X2=0D0
      CF(1)=DCMPLX(X1,X2)
      X1=(-X1/Z)+XCOSH(RRLM,RC,S)/Z
      CF(2)=DCMPLX(X1,X2)
      IF(L.LT.2) GO TO 5
      GO TO 3
2     IF(ITYP.EQ.3) K=-K
      X1=DCOS(Z)/Z
      X2=K*DSIN(Z)/Z
      CF(1)=DCMPLX(X1,X2)
      X1=(DCOS(Z)/(Z**2))+(DSIN(Z)/Z)
      X2=K*((DSIN(Z)/(Z**2))-DCOS(Z)/Z)
      CF(2)=DCMPLX(X1,X2)
      IF(L.LT.2) GO TO 5
*  ITERATION
3     DO 4 LK=2,L
      CF(LK+1)=((2*LK-1)*CF(LK)/Z)-CF(LK-1)
      IF(ITYP.EQ.1) CF(LK+1)=-CF(LK+1)
4     CONTINUE
5     CBESS=CF(L+1)
      RETURN
      END
*
*
      FUNCTION DLR(X,I,L)
      IMPLICIT REAL*8(A-B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)
      PARAMETER (MESH=10)
      COMMON /A/ RLAM,SWW(MESH),RC,RKE
      RSAVE=RC
      RC=X
      S=SWW(I)
      Z=2*RLAM*S*X
      FL=2*RLAM*(X**2+S**2)-3-2*L
      FL1=4*RLAM*S*X
      DLR=FL*REAL(CBESS(Z,L,1))-FL1*REAL(CBESS(Z,L+1,1))
      DLR=-2*RLAM*X*DLR
      RC=RSAVE
      RETURN
      END
*
*
      SUBROUTINE LOADIS(CL,ZASAV,NS1,NOCHNN,NOPEN)
****************************************************************
*   SAVE MATRIX FOR OPEN + DISTORTION CHANNEL CALCULATION      *
*      NOCHNN : TOTAL NUMBER OF CHANNELS (OPEN + CLOSED).      *
*      NOPEN  : NUMBER OF OPEN CHANNELS .                      *
****************************************************************
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      IMPLICIT COMPLEX*16(C)
      PARAMETER (MESH=10,MAXCHN=3)
      DIMENSION CL(MAXCHN,MAXCHN,MESH,MESH)
      DIMENSION ZASAV(MAXCHN*MESH,MAXCHN*MESH)
*
      NS=NS1+1
      DO 10 NCH1=1,NOCHNN
      DO 10 NCH2=1,NOCHNN
      IF (NCH1.LE.NOPEN.AND.NCH2.LE.NOPEN) THEN
         DO 1 I=1,NS1
         DO 1 J=1,NS1
         ZASAV(NS1*(NCH1-1)+I,NS1*(NCH2-1)+J) = CL(NCH1,NCH2,I,J)
1        CONTINUE
      ELSE IF (NCH1.LE.NOPEN.AND.NCH2.GT.NOPEN) THEN
         DO 2 I=1,NS1
         DO 2 J=1,NS
         ZASAV(NS1*(NCH1-1)+I,NOPEN*NS1+(NCH2-NOPEN-1)*NS+J)
     &   = CL(NCH1,NCH2,I,J)
         ZASAV(NOPEN*NS1+(NCH2-NOPEN-1)*NS+J,NS1*(NCH1-1)+I)
     &   = CL(NCH2,NCH1,J,I)
2        CONTINUE
      ELSE IF (NCH1.GT.NOPEN.AND.NCH2.GT.NOPEN) THEN
         DO 3 I=1,NS
         DO 3 J=1,NS
         ZASAV(NOPEN*NS1+(NCH1-NOPEN-1)*NS+I,NOPEN*NS1
     &                  +(NCH2-NOPEN-1)*NS+J) = CL(NCH1,NCH2,I,J)
3        CONTINUE
      END IF
10    CONTINUE
*
      RETURN
      END
*
*
      SUBROUTINE MATCH3(NS,L,CP,CS,RKMOM)
********* COMPUTES P AND S COEFFICIENTS FOR MATCHING
********* GAUSSIAN WITH BESSEL FUNCTION AT RADIUS RC.
      IMPLICIT REAL*8(A-B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)
      PARAMETER (MESH=10)
      DIMENSION CP(MESH),CS(MESH),Z(MESH),CJL(MESH),CJL1(MESH)
      COMMON /A/ RRLM,SWW(MESH),RC,RK
      L1=L+1
      Z1=RC*RKMOM
      CHP=CBESS(Z1,L,2)
      CHP1=CBESS(Z1,L1,2)
      CHM=CBESS(Z1,L,3)
      CHM1=CBESS(Z1,L1,3)
      DO 1 K=1,NS
      Z(K)=2*RRLM*SWW(K)
      ZK=Z(K)*RC
      CJL(K)=CBESS(ZK,L,1)
      CJL1(K)=CBESS(ZK,L1,1)
1     CONTINUE
      RKINV=2./RKMOM
      CWRKH=DCMPLX(0D0,RKINV)
      DO 2 K=1,NS
      CWRKP=2*RRLM*RC*CJL(K)*CHP-Z(K)*CJL1(K)*CHP-RKMOM*CJL(K)*CHP1
      CWRKM=2*RRLM*RC*CJL(K)*CHM-Z(K)*CJL1(K)*CHM-RKMOM*CJL(K)*CHM1
      CP(K)=CWRKH/CWRKP
      CP(K)=CP(K)/(RCTE(K)*(RC**2))
      CS(K)=CWRKM/CWRKP
2     CONTINUE
      RETURN
      END
*
*
      SUBROUTINE SETZAA(ZASAV,ZAA,NS1,NOCHNN,NOPEN)
*****    SET  ZAA    ********
      IMPLICIT COMPLEX*16(Z)
      PARAMETER (MESH=10,MAXCHN=3)
      DIMENSION ZASAV(MAXCHN*MESH,MAXCHN*MESH)
      DIMENSION ZAA(MAXCHN*MESH,MAXCHN*MESH)
      MAX=NOPEN*NS1+(NOCHNN-NOPEN)*(NS1+1)
      DO 1 I=1,MAX
      DO 1 J=1,MAX
*     ZAA(I,J)=0.5*(ZASAV(I,J)+ZASAV(J,I))
      ZAA(I,J)=ZASAV(I,J)
*     ZAA(J,I)=ZASAV(J,I)
1     CONTINUE
      RETURN
      END
*
*
      SUBROUTINE SMTRX(CC,CS,CSMTR,NS,NOPEN)
*********     CALCULATES  S-MATRIX  USING    ****
*********     FORMULAE (4.86)  OF KAMIMURA   ****
      IMPLICIT COMPLEX*16(C)
      PARAMETER (MESH=10,MAXCHN=3)
      DIMENSION CS(MAXCHN,MESH),CC(MAXCHN,MESH,MAXCHN)
      DIMENSION CSMTR(MAXCHN,MAXCHN)
      DO 1 NGAMM=1,NOPEN
      DO 1 INCID=1,NOPEN
      CSMTR(NGAMM,INCID)=0.
      DO 2 I=1,NS
      CSMTR(NGAMM,INCID)=CC(NGAMM,I,INCID)*CS(NGAMM,I)
     &           +CSMTR(NGAMM,INCID)
2     CONTINUE
1     CONTINUE
      RETURN
      END
*
*
      SUBROUTINE LTEST(RK,CL,ZM,ZA,CP,CS,NS,NOCHNN)
*****     THIS SUBROUTINE PRINTS  K-MATRIX, L-MATRIX,     *******
*****     M-MATRIX, A-MATRIX, B-MATRIX  AND               *******
*****     P AND S  COEFFICIENTS.                          *******
      IMPLICIT REAL*8(A-B,D-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      IMPLICIT COMPLEX*16(C)
      PARAMETER (MESH=10,MAXCHN=3)
      DIMENSION RK(MAXCHN,MAXCHN,MESH,MESH)
      DIMENSION CL(MAXCHN,MAXCHN,MESH,MESH)
      DIMENSION ZM(MAXCHN,MESH,MAXCHN)
      DIMENSION ZA(MAXCHN*MESH,MAXCHN*MESH)
      DIMENSION CP(MAXCHN,MESH),CS(MAXCHN,MESH)
      DIMENSION XC(MAXCHN*MESH)
      NS1=NS-1
      NDIM=(NS1+1)*NOCHNN
      NC=NOCHNN*(NOCHNN+1)/2
***********    OUTPUT
      WRITE(6,'(///,10X,A,/)')   'K MATRIX,REAL PART'
      DO 1 NGAMM=1,NOCHNN
      DO 1 NDELT=1,NOCHNN
      WRITE(6,'(/,2X,A,3X,2I5,/)')  'NGAMM,NDELT =',NGAMM,NDELT
      DO 1 I=1,NS
      DO 2 J=1,NS
2     XC(J)=REAL(RK(NGAMM,NDELT,I,J))
1     WRITE(6,'(/,2X,8(10D11.4,/))')  (XC(J),J=1,NS)
*
      WRITE(6,'(///,10X,A,/)')   'L MATRIX,REAL PART'
      DO 5 NGAMM=1,NOCHNN
      DO 5 NDELT=1,NOCHNN
      WRITE(6,'(/,2X,A,3X,2I5,/)')  'NGAMM,NDELT =',NGAMM,NDELT
      DO 5 I=1,NS
      DO 6 J=1,NS
6     XC(J)=REAL(CL(NGAMM,NDELT,I,J))
5     WRITE(6,'(/,2X,8(10D11.4,/))')  (XC(J),J=1,NS)
*
      WRITE(6,'(///,10X,A,/)')   'L MATRIX,IMAGINARY PART'
      DO 7 NGAMM=1,NOCHNN
      DO 7 NDELT=1,NOCHNN
      WRITE(6,'(/,2X,A,3X,2I5,/)')  'NGAMM,NDELT =',NGAMM,NDELT
      DO 7 I=1,NS
      DO 8 J=1,NS
8     XC(J)=DIMAG(CL(NGAMM,NDELT,I,J))
7     WRITE(6,'(/,2X,8(10D11.4,/))')  (XC(J),J=1,NS)
*
      WRITE(6,'(///,10X,A,/)')   'M MATRIX,REAL PART'
      DO 9 NGAMM=1,NOCHNN
      DO 9 INCID=1,NOCHNN
      WRITE(6,'(/,2X,A,3X,2I5,/)')  'NGAMM,INCID =',NGAMM,INCID
      DO 10 I=1,NS
10    XC(I)=REAL(ZM(NGAMM,I,INCID))
9     WRITE(6,'(/,2X,8(10D11.4,/))')  (XC(J),J=1,NS)
*
      WRITE(6,'(///,10X,A,/)')   'M MATRIX,IMAGINARY PART'
      DO 11 NGAMM=1,NOCHNN
      DO 11 INCID=1,NOCHNN
      WRITE(6,'(/,2X,A,3X,2I5,/)')  'NGAMM,INCID =',NGAMM,INCID
      DO 12 I=1,NS
12    XC(I)=DIMAG(ZM(NGAMM,I,INCID))
11    WRITE(6,'(/,2X,8(10D11.4,/))')  (XC(J),J=1,NS)
*
      WRITE(6,'(///,10X,A,/)')   'A MATRIX,REAL PART'
      DO 13 I=1,NDIM
      DO 14 J=1,NDIM
14    XC(J)=REAL(ZA(I,J))
13    WRITE(6,'(/,2X,8(10D11.4,/))')  (XC(J),J=1,NDIM)
*
      WRITE(6,'(///,10X,A,/)')   'A MATRIX,IMAGINARY PART'
      DO 15 I=1,NDIM
      DO 16 J=1,NDIM
16    XC(J)=DIMAG(ZA(I,J))
15    WRITE(6,'(/,2X,8(10D11.4,/))')  (XC(J),J=1,NDIM)
*
      WRITE(6,'(///,10X,A,/)')   'P AND S'
      DO 17 NGAMM=1,NOCHNN
      WRITE(6,'(/,2X,A,3X,2I5,/)')  'NGAMM =',NGAMM
      WRITE(6,'(/,2X,4(9X,A),/)')  'REAL(CP)','IMAG(CP)','REAL(CS)'
     &           ,'IMAG(CS)'
      DO 18 I=1,NS
      X1=REAL(CP(NGAMM,I))
      X2=DIMAG(CP(NGAMM,I))
      X3=REAL(CS(NGAMM,I))
      X4=DIMAG(CS(NGAMM,I))
18    WRITE(6,'(2X,4D17.5)')  X1,X2,X3,X4
17    CONTINUE
      RETURN
      END
