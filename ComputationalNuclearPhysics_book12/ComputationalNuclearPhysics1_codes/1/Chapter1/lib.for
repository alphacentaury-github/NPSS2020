C
C========================================== fitting  version  ========C
C                                                                     C
C    FDU0 Package                                                     C
C    ---- -------                                                     C
C                                                                     C
C    Utility subroutines                                              C
C    ------- -----------                                              C
C                                                                     C
C=====================================================================C
C
      FUNCTION ISECH(INDEXTAB,IEND,INDEX,IRC)
C=====================================================================C
C     Index Search                                                    C
C=====================================================================C
      INTEGER INDEXTAB(1) !table of index in increasing order
      INTEGER IEND        !number of index in the table
      INTEGER INDEX       !index to be searched
      INTEGER IRC         !return control when index does not exist.
C                          irc.lt.0 --> return -1000
C                          irc.ge.0 --> inform and stop
      IP1=1
      IP2=IEND
  1   IF((INDEXTAB(IP2).LT.INDEX).OR.(INDEXTAB(IP1).GT.INDEX)) THEN
       IF(IRC.LT.0) THEN
        ISECH=-1000
        RETURN
       ELSE
        PRINT*,'    Data search failed for index=',index
        PAUSE '***Pause at subroutine ISECH***'
       END IF
      END IF
      ISECH=(IP1+IP2)/2
      IF(INDEXTAB(ISECH).LT.INDEX) THEN
       IP1=ISECH+1
       GOTO 1
      ELSE IF(INDEXTAB(ISECH).GT.INDEX) THEN
       IP2=ISECH-1
       GOTO 1
      END IF
      RETURN
      END

      FUNCTION SIXJ ( JA, JB, JC,
     >                LA, LB, LC  )
C======================================================================C
C  SIXJ IS A 6J-SYMBOL (RACAH-COEFF.); ALL ARGUMENTS ARE GIVEN AS 2J+1 C
C======================================================================C
      REAL*8           GMA(150), SQG(150), TERM, TR, SUM,sixj
      LOGICAL          SET
      SAVE             GMA, SET
      DATA SET /.TRUE./
      TR(N1,N2,N3,N4) = SQG(MA-N2) + SQG(MB-N3) +
     >         SQG(MC-N4) - SQG(N1)
      IF (SET) THEN
         GMA(1) = 0.
         SQG(1) = 0.
         DO 10 I = 3, 149, 2
            GMA(I) = DLOG(DBLE((I-1)/2.0))+GMA(I-2)
            SQG(I) = GMA(I)/2.
   10    CONTINUE
         SET = .FALSE.
      END IF
      N1 = JA + JB + JC
      N2 = JA + LB + LC
      N3 = LA + JB + LC
      N4 = LA + LB + JC
      MA = JB + LB + JC + LC
      MB = JA + LA + JC + LC
      MC = JA + LA + JB + LB
      K1 = MAX (N1, N2, N3, N4)
      KL = MIN (MA, MB, MC)
      SIXJ = 0.
C                              triangular and integer-half integer rules
      IF ( K1.GE.KL ) RETURN
      IF ( MOD(JA+JB+JC-3,2).NE.0 .OR. MOD(JA+LB+LC-3,2).NE.0 .OR.
     >     MOD(LA+JB+LC-3,2).NE.0 .OR. MOD(LA+LB+JC-3,2).NE.0 )
     >            RETURN
      IF(MOD(K1,2).EQ.0) THEN
            SUM=0.
            TERM=0.
      ELSE
            TERM =GMA(K1)-(GMA(MA-K1)+GMA(MB-K1)+GMA(MC-K1)+
     >           GMA(K1+1-N1)+GMA(K1+1-N2)+GMA(K1+1-N3)+GMA(K1+1-N4))
            SUM =DEXP(TERM)*(MOD(K1,4)-2)
            TERM=SUM
      ENDIF
      DO 20 KJ = K1+2, KL, 2
      TERM = -(KJ-1)*(MA-(KJ-1))*(MB-(KJ-1))*(MC-(KJ-1))*TERM/
     >                             ((KJ-N1)*(KJ-N2)*(KJ-N3)*(KJ-N4))
      SUM = SUM + TERM
   20 CONTINUE
      SIXJ = DEXP ( TR(N1, N2, N3, N4) +
     >             TR(N2, N1, N4, N3) +
     >             TR(N3, N4, N1, N2) +
     >             TR(N4, N3, N2, N1) ) * SUM

      RETURN
      END

      SUBROUTINE COPY(FILE1,FILE2)
C=====================================================================C
C     To copy a file which is in logical name file1                   C
C         to logical name file2                                       C
C=====================================================================C
      INTEGER FILE1,FILE2
      CHARACTER*80 LINE
 1    READ(FILE1,11,END=100) LINE
      WRITE(FILE2,11) LINE
 11   FORMAT(80A)
      GOTO 1
 100  REWIND FILE1
      RETURN
      END


      SUBROUTINE EIGRL(A,INDEX,NT,N,LAN,NS,JOBN,D,Z,IZ,WK,IER)
C=====================================================================C
C     LANCZOS diagonalization routine                                 C
C=====================================================================C
      PARAMETER (MAXOTH=15,EPS=1E-3)
      REAL A(1)        !input non-zero matrix elements
      INTEGER INDEX(1) !index(k)=i*(n+1)+j for a(k)
      INTEGER NT       !number of non-zero elements
      INTEGER N        !dimension of a
      INTEGER LAN      !dimension of lanzos matrix
      INTEGER NS       !number of eigenvalues desired(ns.le.lan)
      INTEGER JOBN     !controler 0-->eigenvalues only,1-->both
      REAL    D(1)     !output eigenvalues, dim>=lan
      REAL    Z(IZ,1)  !output eigenvectors,n by lan,z(1:n,k)-->k th eig
      INTEGER IZ       !row dim of z in calling program
      REAL    WK(1)    !working space. >=n+lan*(lan+1)
      INTEGER IER      !error cord
      INTEGER BETA,ZLAN
      BETA(I)=N+I
      ZLAN(I,J)=N+J*LAN+I
C                                           Lanczos matrix
      IF(NS.GT.LAN) STOP 'NS.GT.LAN IN EIGRL'
      IRAN=200001                           !initialize random
      NL=1                                  !number generalator
      ANORM=1.0/SQRT(REAL(N))
      DO 101 IC=1,N
       Z(IC,1)=ANORM
 101  CONTINUE
 89   DO 102 IC=1,N
       WK(IC)=0.0
 102  CONTINUE
      DO 103 IC=1,NT
       I=INDEX(IC)/(N+1)
       J=MOD(INDEX(IC),N+1)
       H=A(IC)
       IF(I.EQ.J) THEN
        WK(I)=WK(I)+Z(I,NL)*H
       ELSE
        WK(J)=WK(J)+Z(I,NL)*H
        WK(I)=WK(I)+Z(J,NL)*H
       END IF
  103 CONTINUE
      OVERLAP=0.
      DO I=1,N
       OVERLAP=OVERLAP+Z(I,NL)*WK(I)
      END DO
      D(NL)=OVERLAP
      IF(NL.EQ.1) THEN
       DO I=1,N
        WK(I)=WK(I)-OVERLAP*Z(I,1)
       END DO
      ELSE
       OVERLAP=0.
       DO I=1,N
        OVERLAP=OVERLAP+Z(I,NL-1)*WK(I)
       END DO
       DO I=1,N
        WK(I)=WK(I)-D(NL)*Z(I,NL)-OVERLAP*Z(I,NL-1)
       END DO
      END IF
      ANORM=0.0
      DO 107 I=1,N
       ANORM=ANORM+WK(I)**2
  107 CONTINUE
      IF(NL.EQ.LAN) GOTO 87
      IF(NL.EQ.1) THEN
       COMPAIR=ANORM/( D(NL)**2 )
      ELSE
       COMPAIR=ANORM/( D(NL)**2+WK(BETA(NL-1))**2 )
      END IF
      ANORM=SQRT(ANORM)
      WK( BETA(NL) )=ANORM
      IF(COMPAIR.LT.1E-10) THEN
C                                           Chose a random state
  85   DO 108 I=1,N
        WK(I)=RAN(IRAM)
  108  CONTINUE
       ANORM=0.0
       DO 113 I=1,N
        ANORM=ANORM+WK(I)**2
 113   CONTINUE
       ANORM=SQRT(ANORM)
      END IF
C                                          generate nl+1 state
      DO I=1,N
       WK(I)=WK(I)/ANORM
      END DO
      IOTH=0
 401  IREP=0
      IOTH=IOTH+1
      IF(IOTH.GT.MAXOTH) THEN
       PRINT*,' STOP AT EIGRL,WHEN NL=',NL
       STOP
      END IF
      DO II=NL,1,-1
       OVERLAP=0.0
       DO I=1,N
        OVERLAP=OVERLAP+Z(I,II)*WK(I)
       END DO
       DO I=1,N
        WK(I)=WK(I)-OVERLAP*Z(I,II)
       END DO
       IF(ABS(OVERLAP).GT.EPS) IREP=1
      END DO
      ANORM=0.
      DO I=1,N
       ANORM=ANORM+WK(I)**2
      END DO
      ANORM=SQRT(ANORM)
      IF(IREP.EQ.1) THEN
       DO I=1,N
        WK(I)=WK(I)/ANORM
       END DO
       GOTO 401
      ELSE
       DO I=1,N
        Z(I,NL+1)=WK(I)/ANORM
       END DO
      END IF
      NL=NL+1
      GOTO 89
C                                            Diagonalize Lanczos matrix
 87   CONTINUE
      DO I=2,LAN
       WK(2*(I-1)+1)=WK(BETA(I-1))
      END DO
      DO I=1,LAN
       WK(2*(I -1)+2)=D(I)
      END DO
      IF(JOBN.EQ.1) THEN
       CALL EVCSB(LAN,WK,2  ,1,D,WK(ZLAN(1,1)),LAN)
       DO J=1,LAN
       SUM=0.
        DO I=1,LAN
         SUM=SUM+WK(ZLAN(I,J))**2
        END DO
        SUM=SQRT(SUM)
        DO I=1,LAN
         WK(ZLAN(I,J))=WK(ZLAN(I,J))/SUM
        END DO
       END DO
C      CALL EQRT2S(D,WK(BETA(1)),LAN,WK(ZLAN(1,1)),LAN,IER)
C
C
C        <--LAN-->
C   ^    Õ       å              Õ      å
C   ]    Õ       å              Õ      å
C   ]    Õ       å Õ      å     Õ      å
C   N    Õ  Z    å*Õ ZLAN å ==> Õ  Z   å
C   ]    Õ       å Õ      å     Õ      å
C   ]    Õ       å              Õ      å
C   V    Õ       å              Õ      å
C
C
       DO 118 I=1,N
        DO 119 II=1,LAN
         WK( BETA(II) )=Z(I,II)
 119    CONTINUE
        DO 120 II=1,NS
         Z(I,II)=0.
         DO 121 III=1,LAN
          Z(I,II)=Z(I,II)+WK( BETA(III) )*WK( ZLAN(III,II) )
 121     CONTINUE
 120    CONTINUE
 118   CONTINUE
      ELSE                       !JOBN=0
C      CALL EQRT1S(D,WK(BETA(1)),LAN,NS,0,IER)
      CALL EVLSB(LAN,WK,2 ,1,D)
      END IF                     !JOBN
      RETURN
      END

      FUNCTION RAN(I)
C=====================================================================C
C     Random Number                                                   C
C=====================================================================C
      A=1234567.0*I+765432
      A=A-1E7*INT(A/1E7)
      RAN=A/1E7
      I=RAN
      RETURN
      END

      SUBROUTINE ERR ( MESG )
C=====================================================================C
C     Error handler                                                   C
C=====================================================================C
      CHARACTER MESG*(*)
C
      WRITE ( *, '(//1X,A//1X,A//)' )
     :  '*** Fatal error encountered ***',MESG
      STOP
      END

      SUBROUTINE SLNN
C=====================================================================C
C     Computes array SLNFISLNJ, SLNF2 -- CG                           C
C=====================================================================C
      COMMON/LNFJ/ SLNF(300),SLNJ(300)
      SAVE /LNFJ/
C
      DO 1 N=2,300
    1 SLNJ(N)=ALOG(FLOAT(N-1))
      SLNJ(1)=0.
      FLN=0.
      DO 2 N=2,300
      FLN=FLN+SLNJ(N)
    2 SLNF(N)=FLN
      SLNF(1)=0.
      RETURN
      END

      FUNCTION CG(J1,M1,J2,M2,J,M)
C=====================================================================C
C     Clebsh-Gordan Coefficient                                       C
C=====================================================================C
      REAL    J1,M1,J2,M2,J,M
      INTEGER X(3),Y(3)
      COMMON/LNFJ/ ALNF(300),ALNJ(300)
      SAVE /LNFJ/,KCALL
      DATA KCALL/1/
C
      GOTO(3333,3334) KCALL
 3333 CALL SLNN
      KCALL=2
 3334 IF((NINT(J1+J2-J).GE.0).AND.(NINT(J-ABS(J1-J2)).GE.0)
     &.AND.(NINT(M1+M2-M).EQ.0).AND.
     *(NINT(ABS(M1)-J1).LE.0).AND.(NINT(ABS(M2)-J2).LE.0))GO TO 10
   30 CG=0.0
      GO TO 50
   10 IF(J1.NE.0.0.AND.J2.NE.0.0) GO TO 20
      CG=1.0
      GO TO 50
   20 X(1)=NINT(J1+J2-J)
      X(2)=NINT(J1-M1)
      X(3)=NINT(J2+M2)
      Y(1)=NINT(J2-M1-J)
      Y(2)=NINT(J1+M2-J)
      Y(3)=0
      IZ2=MIN0(X(1),X(2),X(3))
      IZ1=MAX0(Y(1),Y(2),Y(3))
      IF(IZ1.GT.IZ2) GO TO 30
      K1=NINT(2*J+2)
      K2=NINT(J1+J2-J+1)
      K3=NINT(J1-J2+J+1)
      K4=NINT(-J1+J2+J+1)
      K5=NINT(J1+J2+J+2)
      K6=NINT(J1+M1+1)
      K7=NINT(J1-M1+1)
      K8=NINT(J2+M2+1)
      K9=NINT(J2-M2+1)
      K10=NINT(J+M+1)
      K11=NINT(J-M+1)
      A=(ALNJ(K1)+ALNF(K2)+ALNF(K3)+ALNF(K4)-ALNF(K5)+ALNF(K6)
     *+ALNF(K7)+ALNF(K8)+ALNF(K9)+ALNF(K10)+ALNF(K11))*0.5
      CG=0.0
      IZ2=IZ2+1
      IZ1=IZ1+1
      DO 100 I=IZ1,IZ2
      E=0.0
      DO 200 N=1,3
      N1=(X(N)-I+2)
      N2=(I-Y(N))
  200 E=E+ALNF(N1)+ALNF(N2)
      C=A-E
  100 CG=CG+(-1)**(I-1)*EXP(C)
   50 RETURN
      END

      FUNCTION SLND(A,B,C)
C=====================================================================C
C     Used in UXSH                                                    C
C=====================================================================C
      COMMON/LNFJ/SLNF(300),SLNJ(300)
      SAVE /LNFJ/
C
      I=NINT(A+B-C+1)
      J=NINT(A-B+C+1)
      K=NINT(B-A+C+1)
      L=NINT(A+B+C+2)
      SLND=SLNF(I)+SLNF(J)+SLNF(K)-SLNF(L)
      RETURN
      END

      FUNCTION UXSH(J1,J2,J,J3,J12,J23)
C=====================================================================C
C     Computes Racah coeffients                                       C
C=====================================================================C
      COMMON/LNFJ/ SLNF(300),SLNJ(300)
      SAVE /LNFJ/,KCALL
      REAL    J1,J2,J,J3,J12,J23
      INTEGER X(3),Y(4)
      DATA KCALL/1/
C
      GOTO(3333,3334) KCALL
 3333 CALL SLNN
      KCALL=2
 3334 IF(NINT(J1+J2-J12).LT.0) GO TO 2
      IF(NINT(ABS(J1-J2)-J12).GT.0) GO TO 2
      IF(NINT(J2+J3-J23).LT.0) GO TO 2
      IF(NINT(ABS(J2-J3)-J23).GT.0) GO TO 2
      IF(NINT(J1+J23-J).LT.0) GO TO 2
      IF(NINT(ABS(J1-J23)-J).GT.0) GO TO 2
      IF(NINT(J3+J12-J).LT.0) GO TO 2
      IF(NINT(ABS(J3-J12)-J).GT.0) GO TO 2
      IF(J1.EQ.0.) GO TO 11
      IF(J2.EQ.0.) GO TO 11
      IF(J3.EQ.0.) GO TO 11
      IF(J.EQ.0.) GO TO 11
      X(1)=NINT(J1+J2+J3+J)
      X(2)=NINT(J2+J12+J23+J)
      X(3)=NINT(J12+J1+J23+J3)
      Y(1)=NINT(J1+J2+J12)
      Y(2)=NINT(J1+J+J23)
      Y(3)=NINT(J3+J2+J23)
      Y(4)=NINT(J3+J+J12)
      IZ2=MIN0(X(1),X(2),X(3))
      IZ1=MAX0(Y(1),Y(2),Y(3),Y(4))
      IF(IZ1.GT.IZ2) GO TO 2
      D=SLND(J1,J2,J12)
      D=D+SLND(J2,J3,J23)
      D=SLND(J1,J2,J12)+SLND(J2,J3,J23)+SLND(J1,J23,J)+SLND(J12,J3,J)
      D=(D+SLNJ(NINT(2*J12+2))+SLNJ(NINT(2*J23+2)))*.5
      S=0.
      DO 1000 I=IZ1,IZ2
      A=0.
      DO 200 N=1,3
  200 A=A+SLNF(X(N)-I+1)+SLNF(I-Y(N)+1)
      A=D+SLNF(I+2)-A-SLNF(I-Y(4)+1)
 1000 S=S+(-1)**I*EXP(A)
      UXSH=S*(-1)**X(1)
      RETURN
   11 UXSH=1.
      RETURN
   2  UXSH=0.
      END

      SUBROUTINE GETFILE(SYMMETRY,TYPE,OMI,N,LENGTH,LIMIT,INDEX,
     &    TAB,TAB8)
C=====================================================================C
C     To read in the Matrix Elements files                            C
C=====================================================================C
      CHARACTER SYMMETRY*3,TYPE*(*),FNAME*30
      INTEGER   OMI,INDEX(1),OMI_READ
      REAL      TAB(1)
      REAL*8    TAB8(1)
C
C     FNAME='/ '//SYMMETRY//TYPE//' DAT'     ! IBM type name convention
      FNAME=SYMMETRY//TYPE//'.DAT'           ! VAX type name convention
      OPEN(1,FILE=FNAME,STATUS='OLD',FORM='UNFORMATTED')
      OMI_READ=-1
      N_READ=-1
      DO WHILE((OMI_READ.NE.OMI).OR.(N.NE.N_READ))
       READ(1,END=333) OMI_READ,N_READ
       READ(1)
      END DO
      BACKSPACE(1)
      IF(TYPE.EQ.'JSIZE') THEN
       READ(1) (INDEX(I),I=1,2*N+1)
      ELSE
       BACKSPACE(1)
       READ(1) OMI,N,LENGTH
       IF(LENGTH.GT.LIMIT) THEN
        PRINT*,'LENGTH=',LENGTH,'  FOR FILE '//FNAME
        STOP
       END IF
       IF(TYPE.EQ.'P') THEN
        READ(1) (INDEX(I),TAB8(I),I=1,LENGTH)
       ELSE
        READ(1) (INDEX(I),TAB(I),I=1,LENGTH)
       END IF
      END IF
      CLOSE(1)
      RETURN
 333  PRINT*, 'CAN NOT FIND FILE FOR OMI=',OMI,' N1=',N
      STOP
      END
      SUBROUTINE DIMERR (NAME,I)
C=====================================================================C
C     DIMENSION Error handler                                         C
C=====================================================================C
      CHARACTER NAME*(*)
      WRITE(*,'(1X,A,A,I6)') NAME,' SHOULD>',I
      STOP 'Dimension error detected'
      END
