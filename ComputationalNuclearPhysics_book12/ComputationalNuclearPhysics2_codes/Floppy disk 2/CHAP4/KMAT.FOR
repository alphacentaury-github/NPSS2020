$LARGE
C============================================================C
C   This program solves the coupled channel RGM equations    C
C============================================================C
      IMPLICIT REAL*8 (A-Z)
      INTEGER   I,J,K,L,M,N,KK,IC,EANZ,NS,NMAX,NDIM,IN,II,PNMAX,PKK,
     &          MAXN1,MAXN2,MAXN3,MAXN4,MAXN5,MAXN6,NUMBER,ZTOTJ,IJ,
     &          SIGN,NTOTJ,PK,NIN,MAU
      PARAMETER ( PK = 2  , PKK = PK*(PK-1)/2 , PNMAX =  50,
     &          MAXN1 = PK*(PNMAX+1)    ,MAXN2 = 1+PKK*(PNMAX+1)    ,
     &          MAXN3 = PK*PNMAX*PNMAX  ,MAXN4 = 1+PKK*PNMAX*PNMAX  ,
     &          MAXN5 = PK*PK*(PNMAX+1)*(PNMAX+1),MAXN6=PK*PK*PNMAX )
      INTEGER   INDEX(PK,PK),NZ(PK),LL(PK),MBU(PK)
C
      CHARACTER*13 MBUNAME(PK),INPNAME(PK,PK)
C---------------------------C
C MAXN1 = K*(NS+1)          C
C MAXN2 = K*(K-1)/2*(NS+1)  C
C MAXN3 = K*NS*NS           C
C MAXN4 = K*(K-1)/2*NS*NS   C
C MAXN5 = K*K*(NS+1)*(NS+1) C
C MAXN6 = K*K*NS            C
C---------------------------C
      REAL*8                                     HB(MAXN5),VECT(MAXN6),
     &                  PDIA(MAXN3),VNLO(MAXN2),CHSPIN(PK),RINT(MAXN1),
     &               BB(MAXN1,1),NDIA(MAXN3),PCOUP(MAXN4),NCOUP(MAXN4),
     &            VDIA(MAXN3),VCOUPL(MAXN4),V(MAXN1),EINT(PK),AMUE(PK)
      COMPLEX*16 WAVE(MAXN6)
C-----------------------------------------------------C
C The program needs about                             C
C               2                                     C
C   8 K(3K+1) NS  + 8 K(5K+3) NS  +8 K(K+3)     BYTE  C
C memory.                                             C
C-----------------------------------------------------C

      OPEN(5,FILE='KMAT.INP')
C
          READ(5,53)   ZTOTJ,NTOTJ
          READ(5,*)         K
          READ(5,*)      H0,NMAX
          READ(5,*)     HSTEP,NS
                        IF(K*NS.GT.PK*PNMAX) THEN
                                WRITE(*,*) 'array-dimensions too small'
                                             STOP
                        ENDIF
          READ(5,*)
          READ(5,*)
          READ(5,*)
                   IJ = IDNINT(HSTEP/H0)
                   IF((NS+1)*HSTEP.GT.(NMAX+1)*H0) THEN
                              WRITE(*,*)  'matching radius is too large'
                                           STOP
                                                   ELSE
                      IF(HSTEP.LT.H0) THEN
                                  WRITE(*,*)  'stepsize is too small'
                                           STOP
                      ENDIF
                   ENDIF
C
          READ(5,54)
          READ(5,*)   (AMUE(I),I=1,K)
          READ(5,*)   (NZ(I),I=1,K)
          READ(5,55)
          READ(5,*)   (CHSPIN(I),I=1,K)
          READ(5,*)   (LL(I),I=1,K)
          READ(5,54)
          READ(5,*)   ((INPNAME(I,J),J=I,K),I=1,K)
          READ(5,56)
C


      DO 10 I=1,K
       DO 10 J=I,K
        SIGN    = I*10 + J
         IF(I.EQ.J) THEN
                  INDEX(I,J) =  SIGN
C
                  OPEN(SIGN,FILE=INPNAME(I,J),STATUS='OLD',
     &                 FORM='UNFORMATTED')
         ELSE
          IF(INPNAME(I,J)(1:1).EQ.'+') THEN
                  INDEX(I,J) =  SIGN
C
                  OPEN(SIGN,FILE=INPNAME(I,J)(2: ),STATUS='OLD',
     &                 FORM='UNFORMATTED')
          ELSE
           IF(INPNAME(I,J)(1:1).EQ.'0') THEN
                              INDEX(I,J) = 0
           ELSE
            IF(INPNAME(I,J)(1:1).EQ.'-') THEN
                               INDEX(I,J) = -SIGN
C
                   OPEN(SIGN,FILE=INPNAME(I,J)(2: ),STATUS='OLD',
     &                  FORM='UNFORMATTED')
            ELSE
             WRITE(*,*) 'FIRST CHARACTER OF INPUTFILENAME IS WRONG'
             STOP
            ENDIF
           ENDIF
          ENDIF
         ENDIF
   10 CONTINUE
C
          READ(5,*)         MAU
          READ(5,*)       NUMBER
          READ(5,*)   (MBUNAME(I),I=1,NUMBER)
                                DO 11 I=NUMBER+1,PK
                                    MBUNAME(I) = ' '
                                        MBU(I) = 0
   11                           CONTINUE
C
      DO 12 I=1,NUMBER
         MBU(I) = 50 + I
         OPEN(MBU(I),FILE=MBUNAME(I),STATUS='NEW',FORM='UNFORMATTED')
   12 CONTINUE

          READ(5,*)
          READ(5,*)  EANF
          READ(5,*)  ESTEP,EANZ
C
   53 FORMAT(1X,T4,I2,T8,I2)
   54 FORMAT(1X,//////)
   55 FORMAT(1X,///////)
   56 FORMAT(1X,/)
C
      OPEN(10,FILE='KMAT.OUT',STATUS='NEW')
                        TOTJ = DFLOAT(ZTOTJ)/DFLOAT(NTOTJ)
                          KK = K*(K-1)/2
                        NDIM = K*NS+K
                          IN = NS+1
C
                       WRITE(10,61) ZTOTJ,NTOTJ,HSTEP,NS
C
                    DO 13 I=1,K
                       WRITE(10,62) I,LL(I),CHSPIN(I),AMUE(I)
   13               CONTINUE
C
   61 FORMAT(1X,'TOTAL ANG. MOMENTUM  J : ',I2,'/',I2,7X,
     &       //,1X,'DISCR. STEPSIZE : ',F7.3,3X
     &       ,'NUMBER of MESHPOINTS : ',I4,/)
   62 FORMAT(1X,'CHANNEL : ',I2,3X,'ORB. ANG. MOMENTUM :',I2,3X,
     &       'SPINSTAT. FACTOR : ',F6.2,3X,'AMUE : ',F6.3)
C
                       WRITE(10,54)
C----------------------------------C
C    read matrix elements          C
C----------------------------------C
      DO 14 I=1,K
        NIN = INDEX(I,I)
                                      READ(NIN) EINT(I)
         DO 15 M=1,NMAX+1
          IF((M.GT.IJ*IN).OR.(MOD(M,IJ).NE.0)) THEN
                                       READ(NIN) XXXX
                                            ELSE
                                 READ(NIN) V(M/IJ+(I-1)*IN)
          ENDIF
   15    CONTINUE
C
         DO 16 M=1,NMAX
          DO 16 N=M,NMAX
           IF((M.LE.IJ*NS).AND.(N.LE.IJ*NS).AND.(MOD(M,IJ).EQ.0).AND.
     &       (MOD(N,IJ).EQ.0))  THEN
                       READ(NIN) PDIA(M/IJ+(N/IJ-1)*NS+(I-1)*NS**2),
     &                           NDIA(M/IJ+(N/IJ-1)*NS+(I-1)*NS**2)
                                 PDIA(N/IJ+(M/IJ-1)*NS+(I-1)*NS**2) =
     &                           PDIA(M/IJ+(N/IJ-1)*NS+(I-1)*NS**2)
                                 NDIA(N/IJ+(M/IJ-1)*NS+(I-1)*NS**2) =
     &                           NDIA(M/IJ+(N/IJ-1)*NS+(I-1)*NS**2)
                                      ELSE
                         READ(NIN) XXXX,YYYY
          ENDIF
   16    CONTINUE
   14 CONTINUE
C
         II = 1
      DO 17 I=1,K
       DO 17 J=I+1,K
        IF(INDEX(I,J).NE.0) THEN
         NIN = INDEX(I,J)
          IF(INDEX(I,J).LT.0) NIN = -INDEX(I,J)
C
            DO 18 M=1,NMAX+1
              IF((M.GT.IJ*NS).OR.(MOD(M,IJ).NE.0)) THEN
                                           READ(NIN) XXXX
                                                ELSE
                                   READ(NIN) VNLO(M/IJ+(II-1)*IN)
              ENDIF
   18       CONTINUE
C
            DO 19 N=1,NMAX
             DO 19 M=1,NMAX
              IF((M.LE.IJ*NS).AND.(N.LE.IJ*NS).AND.(MOD(M,IJ).EQ.0).AND.
     &           (MOD(N,IJ).EQ.0))  THEN
                  IF(INDEX(I,J).GT.0) THEN
                       READ(NIN) PCOUP(M/IJ+(N/IJ-1)*NS+(II-1)*NS**2),
     &                           NCOUP(M/IJ+(N/IJ-1)*NS+(II-1)*NS**2)
                                          ELSE
                       READ(NIN) PCOUP(N/IJ+(M/IJ-1)*NS+(II-1)*NS**2),
     &                           NCOUP(N/IJ+(M/IJ-1)*NS+(II-1)*NS**2)
                  ENDIF
                                           ELSE
C
                       READ(NIN) XXXX,YYYY
              ENDIF
   19       CONTINUE
                         ELSE
            DO 20 N=1,NS
                              VNLO(N+(II-1)*IN) = 0.D0
            DO 20 M=1,NS
                 PCOUP(M+(N-1)*NS+(II-1)*NS**2) = 0.D0
                 NCOUP(M+(N-1)*NS+(II-1)*NS**2) = 0.D0
   20       CONTINUE
       ENDIF
       II = II+1
   17 CONTINUE
C
                              WRITE(10,63)
   63 FORMAT(12X,45('*'),/,12X,'* S M A T R I X   AND   S - F A K T ',
     &       'O R S *',/,12X,45('*'),//)
C
              E = EANF
C
      DO 100 IC = 1,EANZ
      WRITE(*,*) 'LOOP :',IC
C--------------------------------------C
C  define potential matrix VDIA        C
C--------------------------------------C
          DO 21 I=1,K
           DO 21 M=1,NS
            DO 21 N=1,NS
              VDIA(M+(N-1)*NS+(I-1)*NS**2)    =
     &               PDIA(M+(N-1)*NS+(I-1)*NS**2) - (E+EINT(1)) *
     &                                   NDIA(M+(N-1)*NS+(I-1)*NS**2)
   21      CONTINUE
C-----------------------------------------C
C  define coupling matrix VCOUPL          C
C-----------------------------------------C
          DO 22 I=1,KK
           DO 22 N=1,NS
            DO 22 M=1,NS
              VCOUPL(M+(N-1)*NS+(I-1)*NS**2) =
     &               PCOUP(M+(N-1)*NS+(I-1)*NS**2) - (E+EINT(1)) *
     &                                  NCOUP(M+(N-1)*NS+(I-1)*NS**2)
   22      CONTINUE
C
            CALL SOLVCC( HB,BB,VDIA,VCOUPL,V,VNLO,RINT,CHSPIN,
     &                   E,EINT,LL,AMUE,NZ,NS,NDIM,HSTEP,K,IN,
     &                   TOTJ,KK,MAU,MBU,VECT,WAVE)
C
            E = E+ESTEP
  100 CONTINUE
C
       CLOSE(5)
       CLOSE(10)
C
      DO 23 I=1,NUMBER
         CLOSE(MBU(I))
   23 CONTINUE
C
      DO 24 I=1,K
        DO 24 J= I,K
            CLOSE(INDEX(I,J))
   24 CONTINUE
C
       STOP
       END
C=====================================================================C
C                                                                     C
C        SOLVES  MULTI CHANNEL  SCHROEDINGER EQUATION  WITH           C
C                     NON-LOCAL POTENTIALS                            C
C                                                                     C
C        CALLS  DFCOUL  AUTOMATICALLY                                 C
C                                                                     C
C                                                                     C
C        STRUCTURE of MATRIX:                                         C
C                                                                     C
C                       NS+1              2NS+2                3NS+3  C
C                         º                  º                   º    C
C        *---------------------------------------------------------*  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º    K(11)     º 0 º    K(12)     º 0 º    K(13)      º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º--------------  *  --------------  *  ---------------  * º  C
C   NS+1-º0 0 0 0 0 0  *  *  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 º  C
C        º--------------  0  --------------  0  ---------------  0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º    K(21)     º 0 º    K(22)     º 0 º    K(23)      º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º--------------  *  --------------  *  ---------------  * º  C
C  2NS+2-º0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 *  * 0 0 0 0 0 0 0 0 0 0 º  C
C        º--------------  0 ---------------  0  ---------------  0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º    K(31)     º 0 º    K(32)     º 0 º    K(33)      º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º              º 0 º              º 0 º               º 0 º  C
C        º--------------  *  --------------  *  ---------------  * º  C
C  3NS+3-º0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 *  * º  C
C        *---------------------------------------------------------*  C
C=====================================================================C
      SUBROUTINE SOLVCC (HB,BB,VDIA,VCOUPL,V,VNLO,RINT,CHSPIN,
     &                   E,EINT,LL,AMUE,NZ,NS,NDIM,HSTEP,K,IN,
     &                   TOTJ,KK,MAU,MBU,VECT,WAVE)
      IMPLICIT REAL*8(A-H,O-Z)
C
      INTEGER PK
      PARAMETER( PK =  3 )
C
      DIMENSION  FN1(PK),GN1(PK),EINT(K),HB(NDIM,NDIM),AMUE(K),V(IN,K)
      DIMENSION  HILF(101),HILFC(101),ER(PK),LL(K),CHNSWI(PK),CONV(PK)
      DIMENSION  SIGMA(PK),RINT(NDIM),CHSPIN(K),RSN(PK),RSM(PK),FN(PK)
      DIMENSION  XKMAT(PK,PK),VECT(NS,K,K),Q(PK),NZ(K)
      DIMENSION  XX(PK),VDIA(NS,NS,K),VCOUPL(NS,NS,KK),DFLL(PK),GN(PK)
      DIMENSION  BB(NDIM,1),VNLO(IN,KK),HILFD(101),MBU(K)
      COMPLEX*16 SMAT(PK,PK),SUMME,BL,QUO,SINT(PK),EMAT(PK,PK),SMADET,
     1           WAVE(NS,K,K),QMAT(PK,PK),III,FLUSS
C
      III     = (0.D0,1.D0)
      PI      = 3.141592653589793D0
      H2M     = 20.735D0
C
       LMIN = LL(1)
       HST2 = HSTEP*HSTEP
      ONE12 = 0.08333333333333333D0*HST2
      TEN12 = HST2-2.D0*ONE12
         KS = NS-1
         MS = NS+1
C
      DO   5  IC=1,K
           CONV(IC) = H2M/AMUE(IC)
           DFLL(IC) = DFLOAT(LL(IC)*LL(IC)+LL(IC))*CONV(IC)
             ER(IC) = E+EINT(1)-EINT(IC)
        IF(E.GT.0.D0)      WRITE(10,99) IC,ER(IC)
        IF(LMIN.GT.LL(IC)) LMIN = LL(IC)
    5 CONTINUE
C
   99 FORMAT(1X,'ENERGY IN CHANNEL NR.',I2,'    : ',F12.6,'*MEV')
                XNS = NS*HSTEP
C
          DO  1011  IC=1,K
            SIGMA(IC) = 0.D0
               FN(IC) = 0.D0
              FN1(IC) = 0.D0
               GN(IC) = 0.D0
              GN1(IC) = 0.D0
C
              Q(IC)=DSQRT(DABS(ER(IC))*(AMUE(IC)/H2M))
C
          IF(ER(IC).GT.0.D0) THEN
             CALL COUFCT(LL(IC),LMIN,Q(IC),AMUE(IC),NZ(IC),XNS,
     &                   HSTEP,FN(IC),FN1(IC),GN(IC),GN1(IC),SIGMA(IC))
          ENDIF
C
                        CHNSWI(IC) =  1.D0
              IF(ER(IC).LT.0.D0) THEN
                        CHNSWI(IC) =  0.D0
              ENDIF
 1011     CONTINUE
C---------------------------C
C declare matrix HB         C
C---------------------------C
           DO 123 M=1,NDIM
           DO 123 N=1,NDIM
            HB(M,N) =0.D0
  123      CONTINUE
C
          DO  2000  IC = 1,K
              IU=(IC-1)*MS
              DO  11  J=1,NS
                                  SIMP =  4.D0/3.D0*HSTEP
                IF(MOD(J,2).EQ.0) SIMP =  2.D0/3.D0*HSTEP
                IF(J.EQ.NS)       SIMP =  1.D0/3.D0*HSTEP
C
                  DO  15  I=1,NS
                      HILF(I) = VDIA(I,J,IC)*SIMP
   15             CONTINUE
                   HB(1+IU,J+IU)=-HILF(1)*TEN12-HILF(2)*ONE12
                  DO  16  I=2,KS
                      HB(I+IU,J+IU) = -( HILF(I-1)+HILF(I+1) ) * ONE12
     &                                               - HILF(I) * TEN12
   16             CONTINUE
                  HB(NS+IU,J+IU) =-HILF(NS-1)*ONE12-HILF(NS)*TEN12
                               X = DFLOAT(J)*HSTEP
                   HB(J+IU,J+IU) = HB(J+IU,J+IU) + ((ER(IC)-V(J,IC)-
     &                             DFLL(IC)/(X*X))*TEN12-2.D0*CONV(IC))
   11     CONTINUE
          DO  12  I=1,MS
                     X = DFLOAT(I)*HSTEP
               HILF(I) = ER(IC)-V(I,IC)-DFLL(IC)/(X*X)
   12     CONTINUE
          DO  13  I=2,NS
              HB(I-1+IU,MS+IU)= 0.D0
              HB(MS+IU,I-1+IU)= 0.D0
              HB(I-1+IU,I+IU) = HB(I-1+IU,I+IU)+(CONV(IC)+ONE12*
     &                          HILF(I))
              HB(I+IU,I-1+IU) = HB(I+IU,I-1+IU)+(CONV(IC)+ONE12*
     &                          HILF(I-1))
   13     CONTINUE
          HB(NS+IU,MS+IU) =  (CONV(IC)+ONE12*HILF(MS))*GN1(IC)
          HB(MS+IU,NS+IU) = - 1.D0*CHNSWI(IC)
          HB(MS+IU,MS+IU) =   GN(IC)+(1.D0-CHNSWI(IC))
            RSN(IC)       = -(CONV(IC)+ONE12*HILF(MS))*FN1(IC)
            RSM(IC)       = - FN(IC)
 2000     CONTINUE
C
       IA = 1
      DO  1010  IC=1,K-1
       IU = (IC-1)*MS
        DO  1020  IK=IC,K-1
         ID = IK*MS
          DO  1030  J=1,NS
                                  SIMP =  4.D0/3.D0*HSTEP
                IF(MOD(J,2).EQ.0) SIMP =  2.D0/3.D0*HSTEP
                IF(J.EQ.NS)       SIMP =  1.D0/3.D0*HSTEP
C
            DO  1040  I=1,NS
              HILFD(I) = VNLO(I,IA)*DKRDEL(I,J)+VCOUPL(I,J,IA)*SIMP
              HILFC(I) = VNLO(I,IA)*DKRDEL(I,J)+VCOUPL(J,I,IA)*SIMP
 1040       CONTINUE
              HB(1+IU,J+ID) =-HILFC(1)*TEN12-HILFC(2)*ONE12
              HB(1+ID,J+IU) =-HILFD(1)*TEN12-HILFD(2)*ONE12
            DO  1050  I=2,KS
              HB(I+IU,J+ID) =-(HILFC(I-1)+HILFC(I+1))*ONE12-HILFC(I)
     &                        *TEN12
              HB(I+ID,J+IU) =-(HILFD(I-1)+HILFD(I+1))*ONE12-HILFD(I)
     &                        *TEN12
 1050       CONTINUE
              HB(NS+IU,J+ID)=-HILFC(NS-1)*ONE12-HILFC(NS)*TEN12
              HB(NS+ID,J+IU)=-HILFD(NS-1)*ONE12-HILFD(NS)*TEN12
 1030     CONTINUE
             HB(NS+IU,MS+ID) = -VNLO(MS,IA)*ONE12*GN1(IK+1)
             HB(NS+ID,MS+IU) = -VNLO(MS,IA)*ONE12*GN1(IC)
         IA = IA+1
 1020   CONTINUE
 1010 CONTINUE
C--------------------------------------------------------C
C divide HB by HBAR^2/(2*Nucleon mass) (=H2M=20.735)     C
C define vector BB                                       C
C--------------------------------------------------------C
              DO 1066 M=1,NDIM
              DO 1066 N=1,NDIM
                HB(M,N) = HB(M,N)/H2M
 1066         CONTINUE
C
          CALL F03AFF(NDIM,EPS,HB,NDIM,DETR,ID,RINT,0)
             DET = DETR*2.D0**ID
                IF(E.LT.0.D0 ) THEN
                                 WRITE(10,533) E,DET
                                       RETURN
                               ELSE
                                 WRITE(10,534) DET
                ENDIF
C
          DO 124 IC = 1,K
               XX(IC) = 0.D0
  124     CONTINUE
C
          DO  2010  M = 1,K
             J = MS*(M-1)+NS
             L = MS*(M-1)+MS
C
               DO 2011 I=1,NDIM
                 BB(I,1) = 0.D0
 2011          CONTINUE
               BB(J,1) = RSN(M)
               BB(L,1) = RSM(M)
C
            CALL F04AJF(NDIM,1,HB,NDIM,RINT,BB,NDIM)
            DO 1388 IC = 1,K
               XKMAT(M,IC) = BB(IC*MS,1)
 1388       CONTINUE
C
            DO 1311 N=1,K
             DO 1377 I=1,NS
               VECT(I,N,M) = BB((N-1)*MS+I,1)
 1377        CONTINUE
 1311       CONTINUE
C
 2010    CONTINUE
C
                       WRITE(10,536)
C
  533    FORMAT(/,1X,'ENERGY : ',F12.6,'*MEV',3X,' DETERMINANT OF'
     &          ,1X,'HB : ',E11.5)
  534    FORMAT(/,1X,'DETERMINANT OF HB (A * 2**R) : ',2X,E11.5,//)
  536    FORMAT(1X,//)
C-------------------------------------------C
C calculate S-matrix from K-matrix using    C
C                                           C
C                1 + i * K                  C
C           S =  ---------                  C
C                1 - i * K                  C
C-------------------------------------------C
         DO 2111 M=1,K
          DO 2111 N=1,K
             QMAT(M,N) = DKRDEL(M,N)-III*XKMAT(M,N)
             EMAT(M,N) = DKRDEL(M,N)
 2111    CONTINUE
C
           CALL F03AHF(K,QMAT,PK,DETR,DETI,ID,SINT,0)
           CALL F04AKF(K,K,QMAT,PK,SINT,EMAT,PK)

         DO 2112 M=1,K
          DO 2112 N=1,K
                 SUMME = 0.D0
           DO 2113 J=1,K
                 SUMME = SUMME+EMAT(M,J)*(DKRDEL(J,N)+III*XKMAT(J,N))
 2113     CONTINUE
                        SMAT(M,N) = SUMME
 2112    CONTINUE

C---------------------------------------------------------------------C
C  calculate absolute value (R) and phase (PHI) of S-matrix elements  C
C               A + B*I = R * EXP(I*PHI)                              C
C  scattering phase is defined in degrees ( 0 - 180*DEG)              C
C---------------------------------------------------------------------C
                     WRITE(10,140)
         DO 1201 M = 1,K
         DO 1201 N = 1,K
                    SBTRAG = CDABS(SMAT(M,N))
                 IF(M.EQ.N) THEN
                     PHASE = 0.5D0*DIMAG(CDLOG((SMAT(M,N)+1.D-60)/
     &                                               (SBTRAG+1.D-60)))
                            ELSE
                     PHASE = DIMAG(CDLOG((SMAT(M,N)+1.D-60)/
     &                                               (SBTRAG+1.D-60)))
                 ENDIF
C
                           DEL = (PHASE*180.D0)/3.141592654D0
C
               IF(DEL.LT.0.D0) DEL = DEL+180.D0
                        IF(CHNSWI(M)*CHNSWI(N).NE.0.D0) THEN
                                     WRITE(10,141) M,N,SBTRAG,DEL
                        ENDIF
 1201    CONTINUE
                     WRITE(10,140)
C-------------------------------------C
C calculate astrophysical S-factors   C
C-------------------------------------C
        DO 1500 I=1,K
         DO 1501 J=1,K
          IF(ER(I).LT.0.D0) THEN
                        SIG = 0.D0
                       SFAK = 0.D0
                           ELSE
        QMAT(I,J) = CDEXP(III*(SIGMA(I)+SIGMA(J)))*SMAT(I,J)
              SIG = PI/Q(I)**2*(2.D0*TOTJ+1.D0)/CHSPIN(I)
     &                        *CDABS(DKRDEL(I,J)-QMAT(I,J))**2/100.D0
             SFAK = ER(I)*SIG*DEXP(0.989D0*DFLOAT(NZ(I))
     &                                         *DSQRT(AMUE(I)/ER(I)))
          ENDIF
                  IF(CHNSWI(I)*CHNSWI(J).NE.0.D0) THEN
                                  WRITE(10,143) I,J,SIG,I,J,SFAK
                  ENDIF
 1501    CONTINUE
                          IF(CHNSWI(I).NE.0.D0) WRITE(10,*)
 1500   CONTINUE
C-------------------------------C
C calculate wave functions      C
C-------------------------------C
      IF(MAU.EQ.1) THEN
           JJ = 1
        DO 1405 J=1,K
         IF(ER(J).GT.0.D0) THEN
                        MBX=MBU(JJ)
                         JJ=JJ+1
                           ELSE
                        GOTO 1405
         ENDIF
C
         DO 1505 I=1,K
             DO 1605 M=1,NS
              SUMME = 0.D0
                DO 1705 IX=1,K
                 SUMME = SUMME -III*(DKRDEL(J,IX)+SMAT(J,IX))*
     1                   VECT(M,I,IX)
 1705           CONTINUE
C
C ---------   WAVE(M,OUT,IN)
C
              FLUSS       = 2.D0*PI*III**(LMIN+1)*CDEXP(III*SIGMA(J))*
     1                                               DSQRT(Q(J)/AMUE(J))
              WAVE(M,I,J) = FLUSS*SUMME
 1605        CONTINUE
          IF(MBX.NE.0)  WRITE(MBX) IDNINT(DKRDEL(I,J)),SMAT(I,J),
     1                                 (WAVE(M,I,J),M=1,NS)
 1505    CONTINUE
 1405   CONTINUE
      ENDIF
C--------------------------------------------------------------C
C calculate determinant of S-matrix; test of unitarity         C
C--------------------------------------------------------------C
            HDIA = 0.D0
            RDIA = 0.D0
         DO 1401 M=1,K
         DO 1401 N=1,K
              SUMME = 0.D0
                   DO 1403 J= 1,K
                    SUMME = SUMME+SMAT(M,J)*DCONJG(SMAT(N,J))
 1403              CONTINUE
           IF(M.EQ.N) THEN
                         ZXZY =DABS(CDABS(SUMME)-1.D0)
                         IF(ZXZY.GT.HDIA) HDIA = ZXZY
                      ELSE
                         ZXZY =CDABS(SUMME)
                         IF(ZXZY.GT.RDIA) RDIA = ZXZY
           ENDIF
 1401    CONTINUE
C
         DO 1300 I=1,K
           IF(SMAT(I,I).EQ.(0.D0,0.D0)) SMAT(I,I) = 1.D0
 1300    CONTINUE
                       CALL F03AHF(K,SMAT,PK,DETR,DETI,ID,SINT,0)
                       SMADET = (DETR+III*DETI)*2**ID
                       ABSDET = CDABS(SMADET)
C
C----  DET(S) = EXP(2*I*PHASE)
C
                        ARGUM = DATAN(DABS(DETI/DETR))
                        ARGUM = (180.D0*ARGUM)/3.141592654D0
                                 IF(DETI.GT.0.D0 ) THEN
                                   IF(DETR.LT.0.D0 ) THEN
                                            ARGUM = 180.D0-ARGUM
                                   ENDIF
                                                   ELSE
                                   IF(DETR.LT.0.D0) THEN
                                            ARGUM = 180.D0+ARGUM
                                                    ELSE
                                            ARGUM = 360.D0-ARGUM
                                   ENDIF
                                  ENDIF
      WRITE(10,142) ABSDET,ARGUM,HDIA,RDIA
C
  140 FORMAT(7X,60('*'),/)
  141 FORMAT(7X,'*',1X,'SMAT(',I1,',',I1,')',T22,'VALUE:',F13.10,3X,
     &       'PHASE :',F9.4,'*DEG *',/)
  142 FORMAT(7X,'DETERMINANT of SMATRIX : ',E13.8,3X,F9.4,'*DEG',//
     &       ,20X,'///  UNITARITY of SMATRIX  ///',/,19X,'HDIA : '
     &       ,E10.5,4X,'NDIA : ',E10.5,/,1X,80('='),/)
  143 FORMAT(6X,' SIG(',I1,'->',I1,') : ',E11.5,' *B ',6X,'S-FAK',
     1       I1,I1,' : ',E11.5,' MEV*B')
      RETURN
      END
C---------------------------------------------C
C  calculate Coulomb-functions at XNS, XNS+1  C
C---------------------------------------------C
      SUBROUTINE   C O U F C T  (L,LS,Q,AMUE,NZ,XN,HSTEP,FN,FN1,GN,GN1,
     1                           PHASE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(31),G(31),FP(31),GP(31),IEXP(31),SIGMA(31)
C
       H2M       = 20.735D0
       GAMMA     = DFLOAT(NZ)*AMUE
       XN1       = XN+HSTEP
       ETA       = GAMMA*1.4409D0/(41.47D0*Q)
       IEXP(L+1) = 0
C
       CALL   D F C O U L ( L,ETA,Q*XN,F,FP,G,GP,IEXP,SIGMA )
C
C ----------- consider factor H2M
C
       FAK = 1.D0
       IF(L.NE.LS) FAK =-1.D0
C
       FN  = DSQRT(AMUE/Q)*FAK*F(L+1)/H2M
       GN  = DSQRT(AMUE/Q)*FAK*G(L+1)
C
       CALL   D F C O U L ( L,ETA,Q*XN1,F,FP,G,GP,IEXP,SIGMA )
C
       FN1 = DSQRT(AMUE/Q)*FAK*F(L+1)/H2M
       GN1 = DSQRT(AMUE/Q)*FAK*G(L+1)
C
       PHASE = SIGMA(L+1)
      RETURN
      END
