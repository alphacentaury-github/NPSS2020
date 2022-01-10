CCCCC&C   V U U - AUXILIARIES TO IMPLEMENT IN VUU.FOR
C   CONTAINS FUNCTION IAND
C            FUNCTION IOR
C            SUBROUTINE PAULI
C            SUBROUTINE DENS
C            SUBROUTINE TRSPRT
C
C     IAND  SIMULATES BINARY AND IN 31 BITS OF TWO INTEGER*4 VARIABLES
C
      FUNCTION IAND(I1IN,I2IN)
      INTEGER*4 I1IN,I2IN,I1,I2,I1D,I2D,I1R,I2R,IS,IR,IM,IAND
      I1=IABS(I1IN)
      I2=IABS(I2IN)
      IS=0
      IM=1
      DO 1 I=1,30
        I1D=I1/IM
        I1R=I1-I1D*IM
        I2D=I2/IM
        I2R=I2-I2D*IM
        IR=MIN(I1R,I2R)
        IS=IS+IR
        I1=I1-I1R
        I2=I2-I2R
        IM=IM*2
 1    CONTINUE
        I1D=I1/IM
        I1R=I1-I1D*IM
        I2D=I2/IM
        I2R=I2-I2D*IM
        IR=MIN(I1R,I2R)
        IS=IS+IR
      IAND=IS
      RETURN
      END
C
C     IOR  SIMULATES BINARY OR IN 31 BITS OF TWO INTEGER*4 VARIABLES
C
      FUNCTION IOR(I1IN,I2IN)
      INTEGER*4 I1IN,I2IN,I1,I2,I1D,I2D,I1R,I2R,IS,IR,IM,IOR
      I1=IABS(I1IN)
      I2=IABS(I2IN)
      IS=0
      IM=1
      DO 1 I=1,30
        I1D=I1/IM
        I1R=I1-I1D*IM
        I2D=I2/IM
        I2R=I2-I2D*IM
        IR=MIN(I1R,I2R)
        IS=IS+IR
        I1=I1-I1R
        I2=I2-I2R
        IM=IM*2
 1    CONTINUE
        I1D=I1/IM
        I1R=I1-I1D*IM
        I2D=I2/IM
        I2R=I2-I2D*IM
        IR=MAX(I1R,I2R)
        IS=IS+IR
      IOR=IS
      RETURN
      END
C
C   SUBROUTINES PAULI,DENS,TRSPRT WITHOUT LARGE VECTORS  DENSPP, XYZSQ
C  WHEN IMPLEMENTING THESE ROUTINES, COMMON BLOCK NOINFB CAN BE DELETED
C
C
C*****************  P A U L I  ************************************
C
      SUBROUTINE PAULI(I,NTAG)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C--
CD     DIMENSION XMESSAGE(10)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      COMMON /PLIBLK/VEFFR(100),VEFFP(100),CONVR,CONVP
      DIMENSION RCM(3),PCM(3)
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 25-JAN-85,12:02:30'/
      DATA PRONAM/'P A U L I  -- PAULI BLOCKING'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
CHK      LOGICAL PROXIM
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (TTVERS) THEN
        WRITE(6,*)'---------------- VERSION NUMBER REPORT -------------'
        WRITE(6,*)PRONAM
        WRITE(6,*)VERNUM
        WRITE(6,*)'----------------------------------------------------'
        TTVERS=.FALSE.
      END IF
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      NTAG=0
      IF (ITYP(I) .GT. 2) RETURN
C--
CD     NCHECK=0
C--
      NCUM=-1
      IO=0
      DO 100 L=1,3
        PCM(L)=0.0
        RCM(L)=0.0
 100  CONTINUE
      DO 101 KKK=1,NSIMUL
CHK        IO=(KKK-1)*(NT+NPIONM)
        DO 102 J=1+IO,NT+IO
          XYZ= RTEST2-(X(J)-X(I))**2
     &            -(Y(J)-Y(I))**2-(Z(J)-Z(I))**2
          IF (XYZ.GT. 0) THEN
           IF (ITYP(J) .LE. 2) THEN
             PXYZSQ= PTEST2-(PX(J)-PX(I))**2
     &            -(PY(J)-PY(I))**2-(PZ(J)-PZ(I))**2
             IF (PXYZSQ.GT. 0) THEN
CHK            IF ((IAND(LOCX(I),LOCX(J)) .NE. 0) .AND.
CHK     *          (IAND(LOCY(I),LOCY(J)) .NE. 0) .AND.
CHK     *          (IAND(LOCZ(I),LOCZ(J)) .NE. 0))       THEN
C--
CD             NCHECK=NCHECK+1
C--
                  NCUM=NCUM+1
                  PCM(1)=PCM(1)+PX(J)
                  PCM(2)=PCM(2)+PY(J)
                  PCM(3)=PCM(3)+PZ(J)
                  RCM(1)=RCM(1)+X(J)
                  RCM(2)=RCM(2)+Y(J)
                  RCM(3)=RCM(3)+Z(J)
              END IF
            END IF
CHK            END IF
          END IF
 102    CONTINUE
        IO=IO+NT+NPIONM
 101  CONTINUE
C
      DO 104 L=1,3
        PCM(L)=PCM(L)/FLOAT(NCUM+1)
        RCM(L)=RCM(L)/FLOAT(NCUM+1)
 104  CONTINUE
      FIR=CONVR*((RCM(1)-X(I))**2+(RCM(2)-Y(I))**2+(RCM(3)-Z(I))**2)
      FIP=CONVP*((PCM(1)-PX(I))**2+(PCM(2)-PY(I))**2+(PCM(3)-PZ(I))**2)
      IRK=INT(FIR+0.5)
      IPK=INT(FIP+0.5)
      F=FLOAT(NCUM)/(ATEST*VEFFR(MIN0(100,IRK+1))*
     * VEFFP(MIN0(100,IPK+1)))
      IF (F.GT.RANF(0)) NTAG=-1
C--
CD     XMESSAGE(1)=NCUM
CD     XMESSAGE(2)=F
CD     XMESSAGE(3)=NCHECK
CD     XMESSAGE(4)=NTAG
CD     WRITE(13)1,XMESSAGE
C--
      RETURN
      END
C
C
      SUBROUTINE DENS(NOUTSI)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
      INTEGER SPHX(2),SPHY(2),SPHZ(2)
C     DIMENSION XMESSAGE(10)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 25-JAN-85,11:52:40'/
      DATA PRONAM/'DENS  -- FORCES AND ACCELERATIONS'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (TTVERS) THEN
        WRITE(6,*)'---------------- VERSION NUMBER REPORT -------------'
        WRITE(6,*)PRONAM
        WRITE(6,*)VERNUM
        WRITE(6,*)'----------------------------------------------------'
        TTVERS=.FALSE.
      END IF
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C  ...  LOOP OVER ALL ENSEMBLES AND PARTICLES
C
      DO 100 KKK=1,NSIMUL
        IO=(KKK-1)*(NT+NPIONM)
        DO 101 IE=1,NT
          I=IE+IO
C
C  ...  ZERO THE COUNTING SPHERE AROUND THE PARTICLE
C
          SPHX(1)=0
          SPHY(1)=0
          SPHZ(1)=0
          SPHX(2)=0
          SPHY(2)=0
          SPHZ(2)=0
C
C  ...  LOOP OVER ALL PARTICLES <> I
C
          DO 102 LLL=1,NSIMUL
            JO=(LLL-1)*(NT+NPIONM)
            DO 103 JE=1,NT
              J=JE+JO
C
C  ...  PRESELECT USING THE INTEGER LOCATOR FUNCTION
C
              IF ( (J .NE. I) .AND. (IAND(LOCX(I),LOCX(J)) .NE. 0) .AND.
     &             (IAND(LOCY(I),LOCY(J)) .NE. 0) .AND.
     &             (IAND(LOCZ(I),LOCZ(J)) .NE. 0) ) THEN
                XIJ=X(I)-X(J)
                YIJ=Y(I)-Y(J)
                ZIJ=Z(I)-Z(J)
                R2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
C
C  ...  IF THE PARTICLE IS IN THE TEST-SPHERE:
C
                IF (R2 .LE. FRAD2) THEN
C
C  ...  FIND WHETHER IT IS ON THE LEFT (L=1) OR RIGHT (L=2) OF THE  I
C       PARTICLE, ON EACH COORDINTAE AXIS AND
C
C
C  ...  COUNT THE PARTICLE IN THE APPROPRIATE HALF-SPHERE
C
                  IF (XIJ .GE. 0.0) THEN
                    SPHX(1)=SPHX(1)+1
                  ELSE
                    SPHX(2)=SPHX(2)+1
                  END IF
                  IF (YIJ .GE. 0.0) THEN
                    SPHY(1)=SPHY(1)+1
                  ELSE
                    SPHY(2)=SPHY(2)+1
                  END IF
                  IF (ZIJ .GE. 0.0) THEN
                    SPHZ(1)=SPHZ(1)+1
                  ELSE
                    SPHZ(2)=SPHZ(2)+1
                  END IF
                END IF
              END IF
 103        CONTINUE
 102      CONTINUE
C
C  ...  COMPUTE THE DENSIT AND THE DERIVATIVE OF THE POTENTIAL
C       (NOTE THAT ANY PAIR OF HALF-SPHERES CONTAINS ALL PARTICLES)
C
          DENSIT=FLOAT(SPHX(1)+SPHX(2))/VOLUME
          III=NINT(DENSIT/DELRHO)
          TEMP=DVRTBL(MIN0(NTABLE,III))
C         XMESSAGE(8)=DENSIT
C         XMESSAGE(9)=TEMP
C
C  ...  COMPUTE THE FORCES AND ACCELERATIONS
C
          TEMP=TEMP*GRAD
          FX(I)=TEMP*FLOAT(SPHX(1)-SPHX(2))
          FY(I)=TEMP*FLOAT(SPHY(1)-SPHY(2))
          FZ(I)=TEMP*FLOAT(SPHZ(1)-SPHZ(2))
          CALL ACCEL(I)
 101    CONTINUE
 100  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE TRSPRT(TADV,BOX)
C
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      LOGICAL CHGED
C  SAVE STORAGE: REPLACE E1 BY A SCALAR
C<>VERNUM<>
      DATA VERNUM/'VER   1.006 -- 25-JAN-85,12:05:30'/
      DATA PRONAM/'TRSPRT -- POSSIBLE CURVED TRAJECTORIES'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (TTVERS) THEN
        WRITE(6,*)'---------------- VERSION NUMBER REPORT -------------'
        WRITE(6,*)PRONAM
        WRITE(6,*)VERNUM
        WRITE(6,*)'----------------------------------------------------'
        TTVERS=.FALSE.
      END IF
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C  ...  TRSPRT ALL PARTICLES AHEAD IN TIME BY TADV
C
      GOTO LABEL
C
C ...         SWITCH(6)=.TRUE.    BOX .GT. 0.0
C          SWITCH(6) CHANGED TO MFIELD
 1000 CONTINUE
C     IF (SWITCH(6)) THEN
        DO 100 I1E=1,NT+NPION(KSIMUL)
          I1=I1E+IOFSET
C
C  ...  TAKE CARE OF  HARD WALL  RUNS
C
C         IF (BOX .GT. 0.0) THEN
            CHGED=.FALSE.
            IF (X(I1) .LT. -BOX) THEN
              PX(I1)=ABS(PX(I1))
              CHGED=.TRUE.
            END IF
            IF (X(I1) .GT. BOX) THEN
              PX(I1)=-ABS(PX(I1))
              CHGED=.TRUE.
            END IF
            IF (Y(I1) .LT. -BOX) THEN
              PY(I1)=ABS(PY(I1))
              CHGED=.TRUE.
            END IF
            IF (Y(I1) .GT. BOX) THEN
              PY(I1)=-ABS(PY(I1))
              CHGED=.TRUE.
            END IF
            IF (Z(I1) .LT. -BOX) THEN
              PZ(I1)=ABS(PZ(I1))
              CHGED=.TRUE.
            END IF
            IF (Z(I1) .GT. BOX) THEN
              PZ(I1)=-ABS(PZ(I1))
              CHGED=.TRUE.
            END IF
            IF (CHGED) CALL CUPDAT(I1,-1)
C         END IF
 100    CONTINUE
C
C ...     SWITCH(6)=.TRUE.  BOX.GT.0.0 AND BOX.LE.0.0
 2000   CONTINUE
        NTNPIO=NT+NPION(KSIMUL)+IOFSET
        DO 200 I1=1+IOFSET,NTNPIO
          E1=SQRT(EM(I1)**2+PZ(I1)**2+PX(I1)**2+PY(I1)**2)
C???          IF (E1 .LE. 0.0) WRITE(6,*)I1E,EM(I1),PX(I1),PY(I1),PZ(I1)
          X(I1)=X(I1)+TADV*PX(I1)/E1
          Y(I1)=Y(I1)+TADV*PY(I1)/E1
          Z(I1)=Z(I1)+TADV*PZ(I1)/E1
 200    CONTINUE
C
C ...  UPDATE THE LOCATOR ARRAYS
C
C CH    DO 201 I1=1+IOFSET,NTNPIO
C DISABEL LOCX(I1)=LOCATE(X(I1),FRAD)
C         LOCY(I1)=LOCATE(Y(I1),FRAD)
C         LOCZ(I1)=LOCATE(Z(I1),FRAD)
C
 201    CONTINUE
      GOTO 4000
C     ELSE
C          SWITCH(6)=.FALSE.
 3000   CONTINUE
C       IF(NPION(KSIMUL).LE.NPIONM)THEN
C        NOFSET=(NT+NPIONM)*(KSIMUL-1)
         NTNPIO=IOFSET+NT+NPION(KSIMUL)
C       ELSE
C        NTNPIO=IOFSET+NT+NPIONM
C        WRITE(6,*) 'PIONNUMBER TOO LARGE ENSEMBLE', KSIMUL
C       END IF
C       WRITE(6,*) KSIMUL,NOFSET,IOFSET,NTNPIO,NPION(KSIMUL),NPIONM
        DO 300 I1=1+IOFSET,NTNPIO
C
          E1=SQRT(EM(I1)**2+PZ(I1)**2+PX(I1)**2+PY(I1)**2)
C
          X(I1)=X(I1)+TADV*(PX(I1)/E1+.5*AX(I1)*TADV)
          Y(I1)=Y(I1)+TADV*(PY(I1)/E1+.5*AY(I1)*TADV)
          Z(I1)=Z(I1)+TADV*(PZ(I1)/E1+.5*AZ(I1)*TADV)
C         PX(I1)=PX(I1)+TADV*AX(I1)
C         PY(I1)=PY(I1)+TADV*AY(I1)
C         PZ(I1)=PZ(I1)+TADV*AZ(I1)
C         CHANGE AX TO FX
          PX(I1)=PX(I1)+TADV*FX(I1)
          PY(I1)=PY(I1)+TADV*FY(I1)
          PZ(I1)=PZ(I1)+TADV*FZ(I1)
C
C ...  UPDATE THE ACCELERATION (FORCES REMAIN UNCHGED)
C
C         CALL ACCEL(I1)
          E1= EM(I1)*EM(I1)+PX(I1)*PX(I1)+PY(I1)*PY(I1)+PZ(I1)*PZ(I1)
          AZ(I1)=PX(I1)*FX(I1)+PY(I1)*FY(I1)+PZ(I1)*FZ(I1)
          AZ(I1)=AZ(I1)/E1
          E1=SQRT(E1)
          AX(I1)=(FX(I1)-PX(I1)*AZ(I1))/E1
          AY(I1)=(FY(I1)-PY(I1)*AZ(I1))/E1
          AZ(I1)=(FZ(I1)-PZ(I1)*AZ(I1))/E1
 300    CONTINUE
C
C ...  UPDATE THE LOCATOR ARRAYS
C
C CH    DO 302 I1=1+IOFSET,NTNPIO
C DISABEL LOCX(I1)=LOCATE(X(I1),FRAD)
C         LOCY(I1)=LOCATE(Y(I1),FRAD)
C         LOCZ(I1)=LOCATE(Z(I1),FRAD)
C
C302    CONTINUE
C     END IF
4000  TIME=TIME+TADV
      RETURN
      ENTRY TSPLBL(BOX)
      IF(.NOT.MFIELD)THEN
       IF(BOX.GT.0.0)THEN
        ASSIGN 1000 TO LABEL
       ELSE
        ASSIGN 2000 TO LABEL
       END IF
      ELSE
        ASSIGN 3000 TO LABEL
      END IF
      RETURN
      END
C     ********************************************************
C     *****************************************************************
