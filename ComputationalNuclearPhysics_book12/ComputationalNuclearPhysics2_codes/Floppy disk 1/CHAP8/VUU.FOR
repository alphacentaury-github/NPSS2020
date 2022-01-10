CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PROGRAM VUU
C
C MODIFIED 11-JAN-83   HANS KRUSE
C MODIFIED 28-JAN-83   TO SHIFT NUCLEUS CENTER-OF-MASS INTO
C  THE CORRECT POSITION. HANS KRUSE
C
C----------------------------------------------------------------
C VERSION 2: 11-AUG-1983
C
C  BOX  ADDED 16-AUG-1983
C
C VERSION 3: 17-AUG-1983
C
C       VERSION 4: 28-AUG-1983
C
C       VERSION 6: 17-SEP-1983 ( INCLUDES ALL VERSION 5 MODIFICATIONS)
C
C       VERSION 7: 25-JAN-85 REVISED VERSION TO RUN ON FPS164
C
C CH   CHANGES 16/01/87 IN SUBROUTINES PAULI, DENS, TRSPRT(AX->FX)
C      FINAL CHANGES IN PREPARING FILE VUU.FOR    C.H.
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
        INCLUDE (VUUPRM)
      CHARACTER*4 HEADER
      CHARACTER*13 FILE
      CHARACTER*3 REF
      CHARACTER*10 DANAME(22)
      LOGICAL FEX
      LOGICAL HIT
      CHARACTER*1 WHAT
      COMMON/AAASK/N1,N2,IZ2,IZ1
      DATA DANAME/
     & 'A-PROJ=N1 ','Z-PROJ=IZ1',' R0-RADIUS','X3-FCT.PR.','PZ-PROJ.  ',
     & 'A-TARG=N2 ','Z-TARG=IZ2','X-SECT-FCT','X3-FCT.TA.','PZ-TARG.  ',
     & 'IMPACTPAR.','TIMEST.SIZ','DELWIDTH-F','SWITCHES  ','PFERMI-FCT',
     & 'EOS-PARM. ','X-SCT-XCUT','EKIMIN-CUT','BOXSIZE   ','SEED-RAND.',
     & 'CHEM.POT. ','TEMP-F.GAS'/
C
        INCLUDE (VUUDIM)
C
      DIMENSION NCOLLT(NSIMUL),NPAUBL(NSIMUL),NDELDC(NSIMUL)
C
        INCLUDE (VUUSIG)
C
        WRITE(6,*)'----------------------------------------------------'
        WRITE(6,*)'  VUU  -  VLASOV UEHLING UHLENBECK  APPROACH        '
        WRITE(6,*)'----------------------------------------------------'
        WRITE(6,*)
C
C<>VERNUM<>
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DATA VERNUM/'VER   7.010 -- 25-JAN-85,12:16:20'/
      DATA PRONAM /'VUU/SCATTER -MAIN -- '/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DATA TTVERS/.TRUE./
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
C  THIS IS PIECE OF CODE TO BE INCLUDED IN THE CASCADE
C
      CONST1=ATAN(CONST)
      CONST2=PI/2.-CONST1
C
C  ...  NUMBER OF NUCLEONS DEFINING THE TEST VOLUME FOR THE FERMI BLOCK
C
      ATEST=4
C
C  ...  INITIALIZE THE  GAP  BETWEEN SIMULTANEOUS RUNS FOR PION STORAGE
C
      NPIONM=NPIONI
C
C  ...  SET DEFAULTS IN THE DA ARRAY
C
      DA(3) = 1.12
C RADIUS PARAMETER, CORRESPONDS TO 0.17 PART/FM**3
      DA(4) = 1.0
C Z-COORDINATE PROJECTILE
      DA(8) = 1.0
C CROSS SECTION PARAMETER - ONLY USED IN VERSION >6
      DA(9) = -1.0
C Z-COORDINATE TARGET
      DA(11)= 0.0
C IMPACT PARAMETER
      DA(12)= 0.25
C TIME STEP
      DA(13)= 1.0
C DELTA DECAY WIDTH MODIFIER
      DA(14)= 0.0
C ALL SWITCHES ARE OFF
      DA(15)= 1.0
C NORMAL FERMI MOMENTUM
      DA(16)=-2.0
C COMPRESSIBILITY PARAMETER
      DA(17)= 1.0
C COMP. FRAME DISTANCE CUTOFF -- MULTIPLIER
      DA(18)= EKININ
C NON-PAULI KINETIC ENERGY CUTOFF
      DA(19)= 0.0
C CONFINED SYSTEM -- BOX SIZE
      DA(20)= 0.0
C RANDOM NUMBER GEN. SEED -- AUTOMATIC
      DA(21)= 0.0
C CHEMICAL POTENTIAL
      DA(22)= 0.0
C TEMPERATURE
C
C ... COLLISION SEARCH CUTOFF PARAMETERS
C
      DO 10 I=1,MAXTYP
        SEACUT(I)=0.0
 10   CONTINUE
C
      DO 20 I=1,3
        SEACUT(I)=1.7
C           NN-COLLISIONS
 20   CONTINUE
C
      SEACUT(6)=1.7
C DELTA-DELTA COLLISIONS
      SEACUT(10)=1.7
      SEACUT(15)=1.7
      SEACUT(21)=1.7
      SEACUT(9)=1.7
      SEACUT(13)=1.7
      SEACUT(14)=1.7
      SEACUT(18)=1.7
      SEACUT(19)=1.7
      SEACUT(20)=1.7
C
      SEACUT(4)=1.7
C DELTA-NUCLEON ELASTIC
      SEACUT(17)=1.7
      SEACUT(5)=1.7
C DELTA-NUCLEON WITH RECOMBINATION
      SEACUT(7)=1.7
      SEACUT(8)=1.7
      SEACUT(11)=1.7
      SEACUT(12)=1.7
      SEACUT(16)=1.7
C
      SEACUT(22)=3.5
C  PION-NUCLEON
      SEACUT(23)=3.5
      SEACUT(29)=3.5
      SEACUT(30)=3.5
      SEACUT(37)=3.5
      SEACUT(38)=3.5
C
C ...  LOG PARAMETERS AND CONSTANTS
C
      WRITE(6,8900)EMNUC,EMPION,EMDELT,
     + WIDDLT,
     + CONST,CONST1,CONST2,
     + DLTLIM
C
 8900 FORMAT('0',
     + ' MASSES (GEV) : NUCLEON = ',F7.3,', PION = ',F7.3,', DELTA = ',
     + F7.3,/,
     + ' DELTA-WIDTH (GEV) = ',F7.3,/,
     + ' DELTA SELECTION CONSTANTS : ',3F9.5,/,
     + ' UPPER LIMIT FOR DELTA RESONANCE SCATTERING = ',F7.3,' GEV',
     + //)
C
      JSIZE=INT(TBLSIZ)
      DO 30 INDEX=1,JSIZE
        SNNTOT(INDEX)=SPPTOT(INDEX)
        SNNEL(INDEX)=SPPEL(INDEX)
 30   CONTINUE
 8910 FORMAT(//,1X,1A,/,1X,'     S (GEV)','       TOTAL','     ELASTIC',
     * '   INELASTIC',/,(1X,1P,4E12.4))
C
C     END DEFAULT VALUES
C
C  ...  SET UP TIMERS
C
      FILE=' '
C
C................................................BEGIN INPUT SECTION
C
      DO 9993 L=1,NSIMUL
      DO 9993 K=1,MAXTYP
      DO 9993 J=1,MXCHAN
      ICOUNT(J,K,L)=0
 9993 CONTINUE
C
C  ...  BEGIN  LINE  INPUT
      WRITE(6,*) 'INPUT OF HEADER (4 CHARACTERS) FORMAT 1A4 E.G. TEST'
      READ(5,8000)HEADER
 8000 FORMAT(1A4)
      WRITE(6,*)'FILE HEADER - ',HEADER
      DO 100 L=1,4
        IF (HEADER(L:L) .EQ. ' ') HEADER(L:L)='X'
 100  CONTINUE
C
C       DEFAULTS FOR TEST RUN
C
        IF (HEADER(1:4) .EQ. 'TEST') THEN
           N1=4
           IZ1=2
           N2=4
           IZ2=2
           NDT=100
           NPRINT=100
           KQT=2
           ELMEV=800.
           REF='NN '
            DA(1)=N1
            DA(2)=IZ1
            DA(6)=N2
            DA(7)=IZ2
            ELAB=ELMEV/1000.
            PPROJ=SQRT(ELAB*(2.*EMNUC+ELAB))
            S2=(ELAB+2*EMNUC)**2-PPROJ**2
            PST=SQRT(S2/4.-EMNUC**2)
            DA(5)=-PST
            DA(10)=PST
           GOTO 3
        END IF
C
      WRITE(6,*) 'INPUT OF A,Z PROJ., A,Z TARG.  FMT=*'
      READ(5,*)N1,IZ1,N2,IZ2
      DA(1)=N1
      DA(2)=IZ1
      DA(6)=N2
      DA(7)=IZ2
      WRITE(6,*) 'INPUT OF #TIMEST.,PRINT-OUT PERIOD,#RUNS,ELAB(MEV)'
      READ(5,*)NDT,NPRINT,KQT,ELMEV
C  !!!  NOTE THAT THE VALUE FOR REF MUST BE ENCLOSED IN QUOTES
C
C  ...  USE THE VALUE OF ELMEV TO PRESET THE MOMENTA DA(5) (PROJ.<0)
C       AND DA(10) (TARG.>0).
      ELAB=ELMEV/1000.
      PPROJ=SQRT(ELAB*(2.*EMNUC+ELAB))
        S2=(ELAB+2*EMNUC)**2-PPROJ**2
        PST=SQRT(S2/4.-EMNUC**2)
        DA(5)=-PST
        DA(10)=PST
C
C  ...  OPTIONAL INPUT CYCLE
C
      WRITE(6,8101)
 8101 FORMAT(1X,80('-'))
      WRITE(6,*) 'INPUT CYCLE FOR DA-BLOCK: INDEX  VALUE  - IND=-1:LIST'
    4 READ(5,*)NA,E
      IF(NA.EQ.-1) GOTO 31
      IF(NA.LT.1) GOTO 3
      IF(NA.GT.22) GOTO 99
      WRITE(6,*)NA,E, DANAME(NA)
      DA(NA)=E
      GOTO 4
 31   CONTINUE
      WRITE(6,*) 'DA-BLOCK : SHORTLIST'
      WRITE(6,32) (I,DANAME(I), I=1,22)
 32   FORMAT(4(5(I3,'=',A10),/),2(I3,'=',A10),16X,'0=END DALOOP',
     & ' -1=THIS LIST')
      GOTO 4
C
    3 CONTINUE
C
C................................................INIT, LOG, PREP. FILES
C
C
C  ...  PARTICLE NUMBERS
C
      WRITE(6,8101)
      DA(6)=AMAX1(DA(6),0.0)
      NT=DA(1)+DA(6)
      N1=DA(1)
      IZ1=DA(2)
C  !!!  NOTE THAT THE VALUE FOR REF MUST BE ENCLOSED IN QUOTES
      N2=DA(6)
      IZ2=DA(7)
      WRITE(6,8100)N1,IZ1,N2,IZ2,NT
 8100 FORMAT(/,' COLLISION OF A,Z=(',2I4,') ON A,Z=(',2I4,'); TOTAL=',
     *       I4)
C
C     CHECK STORAGE OF VECTORS
C
      IF(NT+NPIONM.GT.NPEVMA)THEN
        WRITE(6,*) ' ******* TOO MUCH PARTICLES PER RUN ***********'
        WRITE(6,*) ' NUMBER OF BARYONS',NT
        WRITE(6,*) ' MAX. NUMBER OF PIONS',NPIONM
        WRITE(6,*) ' MAX NUM. TOTAL',NPEVMA
        IF(NT.LT.NPEVMA)THEN
           IF(ELMEV.LE.200.) THEN
C       MOREORLESS NO PIONS
             WRITE(6,*) ' FOR LOW ENERGIES REDUCE NPIONI IN VUUDIM'
     &       ,' (E<200MEV)'
           ELSEIF((NPEVMA*5/6.GT.NT).AND.(ELMEV.LT.1000.)) THEN
C       NPI/A LE 0.2 FOR ELAB LT 1000
             WRITE(6,*) ' FOR SMALL SYSTEMS REDUCE NPIONI IN VUUDIM'
     &       ,' (E<1GEV:NPI/A < 0.2-0.25)'
           END IF
           CALL ERREX(' INSUFFICIENT VECTOR STORAGE')
        END IF
      END IF
C
C  ...  MOMENTA AND GAMMA FACTORS
C
      P01=DA(5)
      P02=DA(10)
      GAM1=SQRT(EMNUC**2+P01*P01)/EMNUC
      GAM2=SQRT(EMNUC**2+P02*P02)/EMNUC
      WRITE(6,8110)P01,GAM1,P02,GAM2
 8110 FORMAT(' PROJ. P,GAMMA: ',2F8.4,'; TARG. P,GAMMA: ',2F8.4)
C
C  ...  RADII
C
      RR1=DA(3)*DA(1)**0.333333333
      RR2=DA(3)*DA(6)**0.333333333
      WRITE(6,9902)RR1,RR2
 9902 FORMAT('RADIUS PROJ, TARG:',2F8.4)
C
C  ...  TIME STEP
C
      DT=DA(12)
      WRITE(6,9901)DT
 9901 FORMAT('TIMESTEP ',1F8.4)
C
C  ...  POSITIONING OF THE NUCLEI
C
      X0=DA(11)/2.
      IF (N2.GT.0) THEN
       Z01=DA(4)*RR1
      ELSE
       Z01=0.0
      END IF
      Z02=DA(9)*RR2
      WRITE(6,9912)Z01,Z02
 9912 FORMAT('Z-COORDINATE PROJ., TARG.: ',2F8.4)
      IF (X0 .GE. 0.0) WRITE(6,9911)X0
 9911 FORMAT('FIXED IMPACT PARAMETER = ',1F6.4)
C
C  ...  CROSS SECTION MODIFIER
C
        IF(DA(8).GT.0.0001) THEN
          DOCOLL=.TRUE.
          XM=DA(8)
        ELSE
C         NO COLLISIONS
          DOCOLL=.FALSE.
          XM=0.0
        END IF
C
C  ...  LAB SYSTEM CUTOFF, ( TIME-SAVER )
C
      XM=XM*DA(17)
      DO 101 I=1,MAXTYP
        SEACUT(I)=SEACUT(I)*XM
        LOCCUT(I)=INT(SEACUT(I)+1.)
 101  CONTINUE
      WRITE(6,8120)(SEACUT(I),LOCCUT(I),I=1,MAXTYP)
 8120 FORMAT(' COMP. SYSTEM SEARCH CUTOFFS',
     * /,(5(1X,F7.1,'(',I5,'),')))
C
C  ...  BOX SIZE FOR THERMAL EQUILIBRIUM RUNS
C
      BOX = DA(19)
      IF (BOX .GT. 0.) WRITE(6,9921)BOX
 9921 FORMAT('SYSTEM CONFINED IN BOX ',1F8.4)
      CALL TSPLBL(BOX)
C
C  ...  SEED FOR RANDOM NUMBERS
C
      ISEED = INT(DA(20)+0.5)
      CALL SEED(ISEED)
C
C  ...  FERMI MOMENTUM AND FERMI TEST RADII
C
C *** DEGENERACY * (4*PI/3)**2 ***
      TEMP=4.*(4.*PI/3.)**2
      PFM=DA(15)*HC/(DA(3)*TEMP**.333333333)
      IF (PFM .GT. 0.) THEN
        RTEST=SQRT(HC*DA(3)*FLOAT(NT)**.333333333/PFM)*
     &         (ATEST/TEMP)**.166666666
C *** USING RADIUS=R0 * A-TOTAL**1/3 ***
        PTEST=PFM*RTEST/(DA(3)*FLOAT(NT)**.333333333)
        RTEST2=RTEST*RTEST
        LOCRAD=INT(RTEST+1.)
        PTEST2=PTEST*PTEST
        ATEST=ATEST*NSIMUL
      ELSE
C *** SELECT RTEST FOR FIELD CALCULATION ***
        RTEST=1.
        RTEST2=1.
      END IF
C
      CALL INITP(RTEST,PTEST)
C
      WRITE(6,9931)PFM
 9931 FORMAT(' FERMI MOMENTUM: ',1F8.4)
      WRITE(6,9904)ATEST,RTEST,PTEST,LOCRAD
 9904 FORMAT(' FERMI TEST A,R,P,LOC:',3F8.4,1I5)
C
C  ...  CONSTANTS FOR THE FIELD CALCULATION;
C
      WRITE(6,*)'RADIUS OF DENSIT-SPHERE : FRAD = ',FRAD
      WRITE(6,9922)VOLUME,GRAD
 9922 FORMAT(' VOLUME, GRAD : ',2F12.6)
C
C  ...  DECODE SWITCHES FROM DA(14)
C
      ISW=INT(ABS(DA(14))+0.5)
      MASK=1
      DO 102 I=1,NUMSWI
        SWITCH(I)=.FALSE.
        IF (IAND(MASK,ISW) .NE. 0) SWITCH(I)=.TRUE.
        MASK=MASK*2
 102  CONTINUE
C
C  ...  SET EKINMI ACCORDING TO SWITCH#4
C
      IF(SWITCH(1)) WRITE(6,*) ' FIXED DELTA DECAY WIDTH '
      IF (SWITCH(4)) THEN
        EKINMI=DA(18)
      ELSE
        EKINMI=.0001
        IF (PFM .LE. 0.0) WRITE(6,*)'??? WARNING, PFM = 0, SW4 OFF'
      END IF
      WRITE(6,*)'MIN. KINETIC ENERGY = ',EKINMI
C
C  ...  SET LOW-KIN.-ENERGY CUTOFF
C
      DO 103 I=1,MAXTYP
        SM2CUT(I)=0.0
 103  CONTINUE
      DO 104 I=1,3
        SM2CUT(I)=(2.*EMNUC+EKINMI)**2
 104  CONTINUE
      XXX=0.0
      SM2CUT(4)=XXX
      SM2CUT(5)=XXX
      SM2CUT(7)=XXX
      SM2CUT(8)=XXX
      SM2CUT(11)=XXX
      SM2CUT(12)=XXX
      SM2CUT(16)=XXX
      SM2CUT(17)=XXX
C
      WRITE(6,8920)SM2CUT
 8920 FORMAT(' LOW ENERGY CUTOFF (SQUARED):',/,9(5F6.1/))
C
C  ...  INIT EQUATION OF STATE
C
      IF(DA(16).EQ.0.0) THEN
       WRITE(6,*) ' CASCADE MODE '
      ELSE
        IF (SWITCH(6)) THEN
C CH ***  CAUTION WHEN CHANGING:  MOMENTA STORED IN FX,FY,FZ
          WRITE(6,*) ' FROZEN FERMI WORKS ONLY IN CASCADE MODE '
          CALL ERREX('DA(16).NE.0  AND SWITCH(6)')
        END IF
        CALL SETEQS(DA(16))
      END IF
      CALL TSPLBL(BOX)
C       ANNOUNCE SET SWITCHES
      IF(SWITCH(1)) WRITE(6,*) ' FIXED DELTA DECAY WIDTH '
      IF(SWITCH(2)) WRITE(6,*) ' NO DELTA DECAY UP TO END OF CALC.'
      IF(SWITCH(3)) WRITE(6,*) ' OUTPUT OF AFTER FIRST TIMESTEP'
      IF(SWITCH(4)) WRITE(6,*) ' NO PAULI BLOCKING '
      IF(SWITCH(5)) WRITE(6,*) ' FIXED DELTA MASS '
      IF(SWITCH(6)) WRITE(6,*) ' FROZEN FERMI DISTRIBUTION'
      IF(SWITCH(7)) WRITE(6,*) ' ALL DELTA DECAY BEFORE FINAL OUTPUT'
      IF(SWITCH(8)) WRITE(6,*) ' REDUCED COLLISIONS '
C
      WRITE(6,*)'INPUT COMPLETE'
      WRITE(6,903)DA
      WRITE(6,*)'SWITCHES :',SWITCH
      WRITE(6,*)'+++++++++++++++++++++++++++++++++++++++++++++++++'
C
C ...  PREPARE THE OUTPUT FILES
C
      IUNIT=98
      TIME=0.0
      DO 4000 K=1,NDT
        TIME=TIME+DT
        KZ=NPRINT*(K/NPRINT)
        IF ( (K.EQ.NDT) .OR.
     &       ( ((K .EQ. 1).AND.SWITCH(3)) .OR. (K .EQ. KZ) ) ) THEN
C
C ........  OUTPUT OF TIMESTEP INFORMATION. FILE NAME IS
C           <HEADER> FOLLOWED BY THE STEP NUMBER, EXTENSION  OUT .
C           THE FIRST RECORD IS CHECKED FOR COMPATIBILITY, THE NEW
C           INFORMATION IS THEN APPENDED.
C
C
C ...  CREATE FILE NAME
C
          WRITE(FILE,8010)HEADER,K
 8010     FORMAT(1A4,I3,'.VUU')
C   HERE A VERSION NUMBER COULD BE IMPLEMENTED
          DO 4010 I=1,7
            IF(FILE(I:I) .EQ. ' ') FILE(I:I)='0'
 4010     CONTINUE
C     WRITE HEADERS WITH APPROPRIATE UNIT AND FILE ID TO MAIN OUTPUT FILE.
          WRITE(7,913) IUNIT,FILE,TIME,KQT
          WRITE(7,903) DA
 903      FORMAT((1X,5F10.4))
          WRITE(7,913) IUNIT,FILE,TIME,KQT
 913      FORMAT(1I3,A13,1E14.7,1I3)
          GOTO 4500
 4500     IUNIT=IUNIT-1
        END IF
 4000 CONTINUE
 4001 CONTINUE
      WRITE(7,913) IUNIT-1,FILE,-999.0,KQT
C
C
      W=1.0
C
C ...  STOP IF THERE ARE NO PARTICLES
      IF (NT .LE. 0)
     *  CALL ERREX( 'NO PARTICLES')
C ...  IF N2 < 0 A THERMAL SOURCE IS TO BE SET UP
      IF (N2 .LT. 0) THEN
        CALL PTINIT(P02)
        WRITE(6,8040)P02
 8040   FORMAT(' THERMAL SOURCE WITH T = ',1P,E12.4,' GEV')
      END IF
C
C................................................LOOP OVER RUNS
C
C
      XBLOCK=0
      XBLCUT=0
      XCOLL=0
      TMFCE=0
      TMCOLL=0
      TMTPRT=0
      TMSCTR=0
      TMPALI=0
C
      DO 1010 KQ=1,KQT
C
C  FIND IMPACT PARAMETER
C
      IF (DA(11) .LT. 0.) X0 = -DA(11)*SQRT(RANF(W))/2.
C
C................................................BEGIN INIT  SIMUL  LOOP
C
      DO 105 KSIMUL=1,NSIMUL
        IOFSET=(KSIMUL-1)*(NT+NPIONM)
C
C  SET INITIAL CONDITIONS
C  PROJECTILE AND TARGET SEPARATELY, REFERENCE FRAME SELECTED BY INPUT
C
C
C................................................INIT PROJECTILE
C
        DO 106 IE=1,N1
          I=IOFSET+IE
          LASTCL(I)=-1
          NCLCNT(I)=0
          DECAY(I)=1.E35
C CH      LOCK(I)=.FALSE.
          R=RR1*RANF(W)**0.3333
          CTH=1.-2.*RANF(W)
          STH=SQRT(1.-CTH**2)
          PHI=PI2*RANF(W)
          Z(I)=(R*CTH+Z01)/GAM1
          X(I)=R*STH*COS(PHI)+X0
          Y(I)=R*STH*SIN(PHI)
C
C
          IF ((N2 .LT. 0).AND.(DA(22).EQ.0.0)) THEN
            CALL PTHERM(PX(I),PY(I),PZ(I))
           ELSE IF ((N2.EQ.0).AND.(DA(22).GT.0.0)) THEN
            CALL FTHERM(PX(I),PY(I),PZ(I),DA(21),DA(22))
          ELSE
            PR=PFM*RANF(0)**(1./3.)
            CTH=1.-2.*RANF(W)
            STH=SQRT(1.-CTH**2)
            PHI=PI2*RANF(W)
            PZ(I)=P01+PR*CTH
            PX(I)=PR*STH*COS(PHI)
            PY(I)=PR*STH*SIN(PHI)
          END IF
          IF (IE .LE. IZ1) THEN
            ITYP(I)=2
            EM(I)=EMNUC
          ELSE
            ITYP(I)=1
            EM(I)=EMNUC
          END IF
 106    CONTINUE
C
        CALL SHIFT(IOFSET,X0,0.,Z01/GAM1,P01,N1)
C
        IF (N2 .GT. 0) GOTO 3800
      IF ((N2.LT.0).AND.(DA(22).EQ.0.0)) THEN
        PXT=0.0
        PYT=0.0
        PZT=0.0
        DO 3000 IE=1,N1
          I=IOFSET+IE
          PXT=PXT+PX(I)
          PYT=PYT+PY(I)
          PZT=PZT+PZ(I)
 3000   CONTINUE
        PXT=PXT/N1
        PYT=PYT/N1
        PZT=PZT/N1
        DO 3010 IE=1,N1
          I=IE+IOFSET
          PX(I)=PX(I)-PXT
          PY(I)=PY(I)-PYT
          PZ(I)=PZ(I)-PZT
 3010   CONTINUE
        WRITE(6,8300)PXT,PYT,PZT
 8300   FORMAT(' C.M.-CORRECTION : ',1P,3E12.4)
      END IF
        GOTO 3900
C
C
C................................................INIT TARGET
C
 3800   DO 107 IE=N1+1,NT
          I=IOFSET+IE
          LASTCL(I)=-1
          NCLCNT(I)=0
          DECAY(I)=1.E35
C  CH     LOCK(I)=.FALSE.
          R=RR2*RANF(W)**0.3333
          CTH=1.-2.*RANF(W)
          STH=SQRT(1.-CTH**2)
          PHI=PI2*RANF(W)
          Z(I)=(R*CTH+Z02)/GAM2
          X(I)=R*STH*COS(PHI)-X0
          Y(I)=R*STH*SIN(PHI)
C
C
          PR=PFM*RANF(0)**(1./3.)
          CTH=1.-2.*RANF(W)
          STH=SQRT(1.-CTH**2)
          PHI=PI2*RANF(W)
          PZ(I)=P02+PR*CTH
          PX(I)=PR*STH*COS(PHI)
          PY(I)=PR*STH*SIN(PHI)
          IF ((IE-N1) .LE. IZ2) THEN
            ITYP(I)=2
            EM(I)=EMNUC
          ELSE
            ITYP(I)=1
            EM(I)=EMNUC
          END IF
 107    CONTINUE
C
        CALL SHIFT(IOFSET+N1,-X0,0.,Z02/GAM2,P02,N2)
C
C CH         FROZEN FERMI - ONLY IN CASCADE MODE
            IF(SWITCH(6))THEN
C            CASCADE: FORCES NOT USED FX,FY,FZ= FROZEN FERMI MOMENTA
              DO 116 IE=1,N1
                I=IOFSET+IE
                 FX(I)=PX(I)
                 FY(I)=PY(I)
                 FZ(I)=PZ(I)-P01
                 PX(I)=0.0
                 PY(I)=0.0
                 PZ(I)=P01
  116         CONTINUE
              DO 117 IE=N1+1,NT
                I=IOFSET+IE
                 FX(I)=PX(I)
                 FY(I)=PY(I)
                 FZ(I)=PZ(I)-P02
                 PX(I)=0.0
                 PY(I)=0.0
                 PZ(I)=P02
  117         CONTINUE
            END IF
C
 3900   NCOLLT(KSIMUL)=0
        NPAUBL(KSIMUL)=0
        NDELDC(KSIMUL)=0
        NPION(KSIMUL)=0
 105  CONTINUE
C
C................................................CLOSE INIT  SIMUL  LOOP
C
      IUNIT=98
      TIME=0.0
C
C................................................LOOP OVER TIME STEPS
C
C
      DO 26 K=1,NDT
C
        TIMEHO=TIME
        TIMMAX=TIME+DT
C
C  ...  UPDATE DENSIT ARRAY
C
        IF (MFIELD) THEN
          CALL DENS(NOUTSI)
        END IF
C
C..............................BEGIN LOOP OVER SIMULTANEOUS ENSEMBLES
C
        DO 108 KSIMUL=1,NSIMUL
          IOFSET=(KSIMUL-1)*(NT+NPIONM)
          NCOLL=0
          TIME=TIMEHO
C
C          COLLISION PART
C
         IF(DOCOLL)THEN
C
C  ...  LOAD COLLISION ARRAYS
C
          CALL COLOAD
C
C  ...  SEARCH FOR NEXT EVENT
C
 1000     CALL FNDNXT(IND,WHAT)
          IF (WHAT .EQ. 'E') GOTO 1999
C *** FINISH TIME-STEP ***
          IF (WHAT .EQ. 'C') THEN
            TNEXEV=CTIME(IND)
          ELSE
            TNEXEV=DECAY(IND)
          END IF
C
          IF (TNEXEV .GT. TIMMAX) THEN
            WRITE(6,*)'??? TIME>TIMMAX',TNEXEV,TIMMAX,WHAT
            GOTO 1999
          END IF
C
C  ...  TRSPRT TO TIME OF COLLISION OR DECAY
C
          TADV=TNEXEV-TIME
          CALL TRSPRT(TADV,BOX)
C
C  ...  COLLISION
C
          IF (WHAT .EQ. 'C') THEN
            NCOLL=NCOLL+1
            I1=I1C(IND)
            I2=I2C(IND)
            SRT=CS(IND)
            SIG=CSIG(IND)
            SIGEL=CSIGEL(IND)
C
            III=MAX0(ITYP(I1),ITYP(I2))
            JJJ=MIN0(ITYP(I1),ITYP(I2))
            ICLTYP=(III*(III-1))/2+JJJ
            IF ( (ICLTYP .LE. 0) .OR. (ICLTYP .GT. MAXTYP) ) THEN
              WRITE(6,*)'ICLTYP,1,2',ICLTYP,I1,ITYP(I1),I2,ITYP(I2)
              GOTO 1000
            END IF
C CH
C   REDUCED COLLISIONS-SWITCH(8): FIRST COLLISION OF A PARTICLE
C     HAS TO BE A COLLISION WITH A PART. FROM OTHER NUCLEUS OR
C     WITH A PART. THAT ALREADY COLLIDED.
            IF(SWITCH(8).AND.(ICLTYP.LE.3))THEN
              IF(LASTCL(I1)+LASTCL(I2).EQ.-2)THEN
                IC12MI=MIN(I1,I2)-IOFSET
                IC12MA=MAX(I1,I2)-IOFSET
                IF((IC12MI.GT.N1).OR.(IC12MA.LE.N1))THEN
C                   BOTH TARGET   OR  PROJECTILE PARTICLES
                   GOTO 1000
                END IF
              END IF
            END IF
C
C  ...  SAVE INFORMATION IN CASE COLLISION WILL BE BLOCKED
C   PUT HERE  BECAUSE OF FROZEN FERMI   CH
              CALL SVEIFO(I1,1)
              CALL SVEIFO(I2,2)
C CH         FROZEN FERMI - ONLY IN CASCADE MODE
            IF(SWITCH(6))THEN
C            CASCADE: FORCES NOT USED FX,FY,FZ= FROZEN FERMI MOMENTA
               IF(LASTCL(I1).EQ.-1)THEN
                  PX(I1)=PX(I1)+FX(I1)
                  PY(I1)=PY(I1)+FY(I1)
                  PZ(I1)=PZ(I1)+FZ(I1)
                  FX(I1)=0.0
                  FY(I1)=0.0
                  FZ(I1)=0.0
               END IF
               IF(LASTCL(I2).EQ.-1)THEN
                  PX(I2)=PX(I2)+FX(I2)
                  PY(I2)=PY(I2)+FY(I2)
                  PZ(I2)=PZ(I2)+FZ(I2)
                  FX(I2)=0.0
                  FY(I2)=0.0
                  FZ(I2)=0.0
               END IF
            END IF
C CH  UPDATE COLLISION NUMBERS
              NCLCNT(I1)=NCLCNT(I1)+1
              NCLCNT(I2)=NCLCNT(I2)+1
C  ...  STORE PARTICLE PROPERTIES
            X1=X(I1)
            Y1=Y(I1)
            Z1=Z(I1)
            X2=X(I2)
            Y2=Y(I2)
            Z2=Z(I2)
            PZ1=PZ(I1)
            PZ2=PZ(I2)
            PX1=PX(I1)
            PX2=PX(I2)
            PY1=PY(I1)
            PY2=PY(I2)
C  ...  FIND ENERGIES AND BETA OF 2-BODY SYSTEM
            E1=SQRT(EM(I1)**2+PZ1**2+PX1**2+PY1**2)
            E2=SQRT(EM(I2)**2+PZ2**2+PX2**2+PY2**2)
            BX=(PX1+PX2)/(E1+E2)
            BY=(PY1+PY2)/(E1+E2)
            BZ=(PZ1+PZ2)/(E1+E2)
            GAM=1./SQRT(1.-BX**2-BY**2-BZ**2)
C  ...  FIND RELATIVE MOMENTUM IN THE 2-BODY CM FRAME
            B2=.0001+BX**2+BY**2+BZ**2
            PB=BX*PX1+BY*PY1+BZ*PZ1
            GG=(GAM-1.)*PB/B2-GAM*E1
            PZR=GG*BZ+PZ1
            PXR=GG*BX+PX1
            PYR=GG*BY+PY1
C
C  ...  COLLISIONS WITH TWO-BODY EXIT CHANNELS
C
            IF (ICLTYP .LE. 21) THEN
C
              CALL SCATTR(PXR,PYR,PZR,SRT,I1,I2,
     &                    EM(I1),EM(I2),ITYP(I1),ITYP(I2),
     &                    ICLTYP,SIG,SIGEL,ICHANL)
C
              PB=BX*PXR+BY*PYR+BZ*PZR
              EA=SQRT(EM(I1)**2+PXR**2+PYR**2+PZR**2)
              GG=(GAM-1.)*PB/B2+GAM*EA
              PX(I1)=GG*BX+PXR
              PY(I1)=GG*BY+PYR
              PZ(I1)=GG*BZ+PZR
              EB=SQRT(EM(I2)**2+PXR**2+PYR**2+PZR**2)
              GG=(GAM-1.)*(-PB)/B2+GAM*EB
              PX(I2)=GG*BX-PXR
              PY(I2)=GG*BY-PYR
              PZ(I2)=GG*BZ-PZR
C
C  ...  CHECK ON PAULI PRINCIPLE
C
              IF (SWITCH(4)) THEN
                NTAG1=0
                NTAG2=0
              ELSE
                CALL PAULI(I1,NTAG1)
                CALL PAULI(I2,NTAG2)
              END IF
             IF(ICHANL.GT.4) WRITE(6,*)' *** WARNING CANNEL NO',ICHANL
              IF ((NTAG1 .LT. 0) .OR. (NTAG2 .LT. 0)) THEN
                NPAUBL(KSIMUL)=NPAUBL(KSIMUL)+1
                XBLOCK=XBLOCK+1
                CALL RSTIFO(I1,1)
                CALL RSTIFO(I2,2)
                I1C(IND)=-1
C CH COUNT PAULIBLOCKS IN CHANNEL 5 -FREE BECAUSE OF CHANGE OF SWITCH(2)
                ICOUNT(5,ICLTYP,KSIMUL)=
     &            ICOUNT(5,ICLTYP,KSIMUL)+1
              ELSE IF (ICHANL .GT. 0) THEN
                ICOUNT(ICHANL,ICLTYP,KSIMUL)=
     &               ICOUNT(ICHANL,ICLTYP,KSIMUL)+1
C
C  ...  VARIOUS BOOKEEPING FUNCTIONS OCCUR HERE
C
C  ...  RECOMPUTE THE ACCELERATIONS FOR NUCLEONS AND DELTAS
C
                IF (MFIELD) THEN
                  IF (ITYP(I1) .LE. 6) CALL ACCEL(I1)
                  IF (ITYP(I2) .LE. 6) CALL ACCEL(I2)
                END IF
C
C  ...    SEE WHETHER WEE MUST UPDATE DECAY
C
                CALL DUPDAT(I1)
                CALL DUPDAT(I2)
C
C  ...  UPDATE THE COLLISION ARRAY
C
                CALL CUPDAT(I1,I2)
              END IF
            ELSE
C
C  ...  ONE BODY EXIT CHANNELS, ACCEL. WILL BE UPDATED IN  ABSORB
C
              CALL ABSORB(I1,I2,SRT,ICLTYP)
C
C  ...  UPDATE COLLISION ARRAY
C
              CALL CUPDAT(I1,I2)
C

            END IF
C
          END IF
C
C  ...  DECAY EVENTS
C
          IF (WHAT .EQ. 'D') THEN
            CALL DELDEC(IND)
            NDELDC(KSIMUL)=NDELDC(KSIMUL)+1
C
C  ...  UPDATE THE COLLISION ARRAY
C
            IF (MFIELD) CALL ACCEL(IND)
            CALL CUPDAT(IND,IOFSET+NT+NPION(KSIMUL))
C
          END IF
C
C................................................COLLISION LOOPS CLOSED
C
          GOTO 1000
C
C ...  RECOVER ABSORBED PIONS
C
 1999     IF (NPION(KSIMUL) .GT. 0) THEN
            NPOINT=NPION(KSIMUL)
  109       IF  ((NPION(KSIMUL) .GT. 0) .AND. (NPOINT .GT. 0)) THEN
              IF (ITYP(NT+NPOINT+IOFSET) .LT. 0) THEN
                IF (NPOINT .LT. NPION(KSIMUL)) THEN
                  DO 110 I=NT+NPOINT+1+IOFSET,NT+NPION(KSIMUL)+IOFSET
                    J=I-1
                    X(J)=X(I)
                    Y(J)=Y(I)
                    Z(J)=Z(I)
                    PX(J)=PX(I)
                    PY(J)=PY(I)
                    PZ(J)=PZ(I)
                    EM(J)=EM(I)
                    ITYP(J)=ITYP(I)
                    DECAY(J)=DECAY(I)
                    LASTCL(J)=LASTCL(I)
                    NCLCNT(J)=NCLCNT(I)
 110              CONTINUE
                END IF
                NPION(KSIMUL)=NPION(KSIMUL)-1
              END IF
              NPOINT=NPOINT-1
              GOTO 109
          END IF
          END IF
C *** CLOSE PION STORAGE COMPRESSION
C
         END IF
C
C        END OF COLLISION PART
C
C
C  ...  TRSPRT TO END OF TIME-SLICE
C
          IF (TIME .LT. TIMMAX) THEN
            CALL TRSPRT(TIMMAX-TIME,BOX)
          END IF
C ................................................STATS AND OUTPUT
C
          IF ( SWITCH(7).AND.(K.EQ.NDT) )THEN
             DO 118 I=IOFSET+1,IOFSET+NT
               IF(ITYP(I).GT.2) THEN
                 CALL DELDEC(I)
                 NDELDC(KSIMUL)=NDELDC(KSIMUL)+1
               END IF
  118        CONTINUE
          END IF
C
          NCOLLT(KSIMUL)=NCOLLT(KSIMUL)+NCOLL
          XCOLL=XCOLL+NCOLL
          KZ=NPRINT*(K/NPRINT)
          IF ( (K.EQ.NDT) .OR.
     &      ( ((K .EQ. 1).AND.SWITCH(3)) .OR. (K .EQ. KZ) ) ) THEN
C
C  ...   WRITE NEW INFORMATION
C
            WRITE(7,912)IUNIT,MXCHAN,MAXTYP,KSIMUL
            WRITE(7,905)IUNIT,NT,NPION(KSIMUL),2.0*X0
          WRITE(7,906)IUNIT,NCOLLT(KSIMUL),NPAUBL(KSIMUL),NDELDC(KSIMUL)
            WRITE(7,908) ((ICOUNT(I,J,KSIMUL),J=1,21),I=1,MXCHAN)
            WRITE(7,909) ((ICOUNT(I,J,KSIMUL),J=22,MAXTYP),I=1,MXCHAN)
C
            IF(SWITCH(6))THEN
C            CASCADE: FORCES NOT USED FX,FY,FZ= FROZEN FERMI MOMENTA
            DO 4058 I=IOFSET+1,IOFSET+NT+NPION(KSIMUL)
            IF(LASTCL(I).GT.0)THEN
              LSTCL=LASTCL(I)-IOFSET
            WRITE(7,907) X(I),Y(I),Z(I),PX(I),PY(I),PZ(I),EM(I),ITYP(I),
     &           LSTCL,NCLCNT(I)
            ELSE
              LSTCL=LASTCL(I)
C WRITE OUT MOMENTA WITH FROZEN FERMI PART TO CONSERVE TOTAL MOMENTUM
            WRITE(7,907) X(I),Y(I),Z(I),PX(I)+FX(I),PY(I)+FY(I),
     &           PZ(I)+FZ(I),EM(I),ITYP(I),LSTCL,NCLCNT(I)
            END IF
 4058       CONTINUE
            ELSE
            DO 4059 I=IOFSET+1,IOFSET+NT+NPION(KSIMUL)
            IF(LASTCL(I).GT.0)THEN
              LSTCL=LASTCL(I)-IOFSET
            ELSE
              LSTCL=LASTCL(I)
            END IF
            WRITE(7,907) X(I),Y(I),Z(I),PX(I),PY(I),PZ(I),EM(I),ITYP(I),
     &           LSTCL,NCLCNT(I)
 4059       CONTINUE
            END IF
 912  FORMAT(4I4)
 908  FORMAT((3I5,18I3))
 909  FORMAT((24I3))
 906  FORMAT(1I3,3I8)
 905  FORMAT(1I3,1I5,1I6,1E14.7)
 907  FORMAT(3E13.5,3F8.4,F7.4,I2,I4,I3)
            IF (KSIMUL .EQ. NSIMUL) IUNIT=IUNIT-1
          END IF
C ------  END OF OUTPUT SECTION
 108    CONTINUE
C ***  ...............  SIMUL  BLOCK CLOSED
  26  CONTINUE
C
C................................................TIME LOOP CLOSED
C
      CALL SCATTR(0,0,0,-1.,0,0,0,0,0,0)
C
C  ...  TIMING UPDATE
C
 1010 CONTINUE
C
C................................................ENSEMBLE LOOP CLOSED
C
C ...  TIMING INFORMATION
C
 99   WRITE(6,8093) XCOLL/FLOAT(NSIMUL*KQT),
     * XBLOCK/FLOAT(NSIMUL*KQT),XBLCUT/FLOAT(NSIMUL*KQT)
 8093 FORMAT(1X,/,'                      PER ENSEMBLE',/,
     * '-----------------------------------',/,
     * ' COLLISIONS      = ',F10.4,/,
     * ' FERMI BLOCKS    = ',F10.4,/,
     * ' LOW ENERGY CUTS = ',F10.4,/)
C
      STOP
      END
C
C*********************    S C A T T R   *******************************
C
      SUBROUTINE SCATTR(PXR,PYR,PZR,SRT,IND1,IND2,
     &  EM1,EM2,ITYP1,ITYP2,ICLTYP,
     &  SIG,SIGEL,ICHANL)
C<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        INCLUDE (VUUPRM)
        INCLUDE (VUUSIG)
        INCLUDE (VUUDIM)
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      DIMENSION SIGIN(MXCHAN)
C<>VERNUM<>
      DATA VERNUM/'VER   6.020 -- 25-JAN-85,13:22:00'/
      DATA PRONAM/'MODULE: S C A T T R'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
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
      IF (SRT .LE. 0) GOTO 9700
      IF (ICLTYP .GT. MAXTYP) THEN
        WRITE(6,*)' INVALID TYPE IN  SCATTR  : ',ICLTYP
        CALL ERREX( 'FATAL ERROR')
      END IF
C
      S=SRT
      PR=SQRT(PXR*PXR+PYR*PYR+PZR*PZR)
      C2=PZR/PR
C
C ...  DETERMINE ELASTIC/INELASTIC CROSSECTIONS AND SELECT BRANCH
C
        RAND=RANF(0)
        IF (RAND .GT. SIGEL/SIG) GOTO 3000
C
C ...  ELASTIC SCATTERING BRANCH STARTS HERE
C
 2000 LASTCL(IND1)=IND2
      LASTCL(IND2)=IND1
C
      ICHANL=1
C
C ...  GET ANGULAR DISTRIBUTION
C
      CALL ANGEL(S,PR,ICLTYP,T1,C1)
      GOTO 5000
C
C ...  INELASTIC SCATTERING BRANCH
C
 3000 CALL CROSS2(S,PR,EM1,EM2,ICLTYP,SIG,SIGEL,SIGIN,ICNUM)
C
C ...  SELECT AN INELASTIC EXIT CHANNEL
C
C
      IF (ICNUM .LE. 1) GOTO 2000
C
      SUM=SIGEL
      DO 10 I=2,ICNUM
      IC=I-1
      SUM=SUM+SIGIN(IC)
      IF (RAND .LE. SUM/SIG) GOTO 3100
   10 CONTINUE
      ICHANL=0
      RETURN
C
C ...  CHANNEL SELECTED, FIND NEW MASSES (AND TYPES)
C
 3100 IF (SWITCH(5) .AND. (S .LT. EMNUC+EMDELT))
     *   GOTO 2000
      CALL CMASS(S,PR,EM1,EM2,ITYP1,ITYP2,ICLTYP,IC)
      IF (IC .EQ. 0) GOTO 2000
C
C ...  COMPUTE MAGNITUDE OF REL. MOMENTUM USING THE NEW MASSES
C
C
C  ...  TAG BOTH PARTICLES USING NEGATIVE TYPES; THE MAIN PROGRAM
C       MUST STILL UPDATE THE DECAY ARRAY
C
      ITYP1=-ITYP1
      ITYP2=-ITYP2
C
      LASTCL(IND1)=IND2
      LASTCL(IND2)=IND1
C
      ICHANL=IC+1
      Y1=EM1*EM1
      Y2=EM2*EM2
      Y3=S*S-Y1-Y2
      PR=SQRT(Y3*Y3-4.*Y1*Y2)/(2.*S)
C
C ...  GET THE ANGULAR DISTRIBUTION
C
C ...  SUBROUTINE CALL DELETED UNTIL SUCH TIME WHEN WE WILL HAVE
C      ANG. DISTR. OTHER THAN ISOTROPIC ONES.
C
C      CALL ANGIN(S,PR,ICLTYP,T1,C1)
      C1=1.-2.*RANF(0)
      T1=PI2*RANF(0)
C
C ...  ROTATE MOMENTUM VECTOR TO NEW POSITION AND RETREAT
C
 5000 S2=SQRT(1.-C2**2)
      IF(PXR .EQ. 0.0 .AND. PYR .EQ. 0.0 ) THEN
        T2=0.0
      ELSE
        T2=ATAN2(PYR,PXR)
      END IF
      S1=SQRT(1.-C1**2)
      CT1=COS(T1)
      ST1=SIN(T1)
      CT2=COS(T2)
      ST2=SIN(T2)
      PZR=PR*(C1*C2-S1*S2*CT1)
      SS=C2*S1*CT1+S2*C1
      PXR=PR*(SS*CT2-S1*ST1*ST2)
      PYR=PR*(SS*ST2+S1*ST1*CT2)
      RETURN
C
C ...  RESET ALL COUNTERS
C
 9700 DO 9710 K=1,NSIMUL
      DO 9710 I=1,MAXTYP
      DO 9710 J=1,MXCHAN
      ICOUNT(J,I,K)=0
 9710 CONTINUE
      RETURN
      END
C
C************************   S A V E  I N F O    *****************
C
      SUBROUTINE SVEIFO(IP,IS)
C
      DIMENSION PKEEP(3,2),EKEEP(2),ITKEEP(2),DECEEP(2),
     &          LCKEEP(2),AKEEP(3,2),FKEEP(3,2),NCKEEP(2),
     &          LKKEEP(2)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
      IF ((IS .LE. 0) .OR. (IS .GT. 2)) CALL ERREX('IS -- SAVE')
C
      PKEEP(1,IS)=PX(IP)
      PKEEP(2,IS)=PY(IP)
      PKEEP(3,IS)=PZ(IP)
      EKEEP(IS)=EM(IP)
      ITKEEP(IS)=ITYP(IP)
      DECEEP(IS)=DECAY(IP)
      LCKEEP(IS)=LASTCL(IP)
      NCKEEP(IS)=NCLCNT(IP)
      AKEEP(1,IS)=AX(IP)
      AKEEP(2,IS)=AY(IP)
      AKEEP(3,IS)=AZ(IP)
      FKEEP(1,IS)=FX(IP)
      FKEEP(2,IS)=FY(IP)
      FKEEP(3,IS)=FZ(IP)
C
      RETURN
C
C-------------------  R E S T O R E  I N F O   ---------------------
C
      ENTRY RSTIFO(IP,IS)
C
      IF ((IS .LE. 0) .OR. (IS .GT. 2)) CALL ERREX('IS -- RESTORE')
C
      PX(IP)=PKEEP(1,IS)
      PY(IP)=PKEEP(2,IS)
      PZ(IP)=PKEEP(3,IS)
      EM(IP)=EKEEP(IS)
      ITYP(IP)=ITKEEP(IS)
      DECAY(IP)=DECEEP(IS)
      LASTCL(IP)=LCKEEP(IS)
      NCLCNT(IP)=NCKEEP(IS)
      AX(IP)=AKEEP(1,IS)
      AY(IP)=AKEEP(2,IS)
      AZ(IP)=AKEEP(3,IS)
      FX(IP)=FKEEP(1,IS)
      FY(IP)=FKEEP(2,IS)
      FZ(IP)=FKEEP(3,IS)
C
      RETURN
      END
C
C**********************   P R O X I M   ****************************
C
C%%%%%%%%%%%%%%%%  CURRENTLY NOT USED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION PROXIM(M1,M2,IW)
C
      INTEGER SHIFT
      M=0
      DO 100 I=-IW,IW,1
        M=IOR(M,SHIFT(M1,I))
 100  CONTINUE
      PROXIM=(IAND(M,M2) .NE. 0)
C
      RETURN
      END
C*********************   E R R E X  ******************************
C
      SUBROUTINE ERREX(TEXT)
      CHARACTER*(*) TEXT
C
C  ...  TRIP ERROR CONDITION  NON-FORTRAN SPECIFIC ERROR
C
      WRITE(6,*)TEXT
      STOP
      RETURN
      END
C
C *************  F T H E R M   ************************************
C
      SUBROUTINE FTHERM(PX,PY,PZ,U,T)
C
C THIS SUBROUTINE INITIALIZES A FINITE TEMPERATURE FERMI GAS
C GIVEN THE CHEMICAL POTENTIAL U AND TEMPERATURE T AS INPUTS
C
      PARAMETER (PI=3.1415927)
      REAL M
      M=938.0
      DP=0.1
      NUP=12000
      P=0.0
      Y=0.0
      DO 100 I=1,NUP
       P=P+DP
       RE=SQRT(P**2+M**2)
       IF (RE/T-U/T.LE.40.0) Y=Y+DP*(P**2)/(EXP(RE/T-U/T)+1.0)
 100  CONTINUE
      PABS=RANF(1.0)
      P=0.0
      X=0.0
      DO 200 J=1,NUP
       P=P+DP
       RE=SQRT(P**2+M**2)
       IF (RE/T-U/T.LE.40.0) X=X+DP*(P**2)/(EXP(RE/T-U/T)+1.0)
       IF (X/Y .GE. PABS) GOTO 1000
 200  CONTINUE
 1000 P=P/1000.
      CTH=1.-2.*RANF(1.)
      STH=SQRT(1.-CTH*CTH)
      PHI=2.*PI*RANF(1.)
      PX=P*STH*COS(PHI)
      PY=P*STH*SIN(PHI)
      PZ=P*CTH
      RETURN
      END
C
C *************  P T H E R M   AND    P T H E R M  I N I T   ***********
C
C%%%%%%%%%%%%%%% TEMPORARILY DISABLED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE PTHERM(PX,PY,PZ)
        INCLUDE (VUUPRM)
      PARAMETER (NBIN=1000,TMULT=8.,TMIN=0.01)
      COMMON /STORE/ BIN(NBIN),DE
      IF (III .EQ. 0) CALL ERREX( 'PTHERM NOT IMPLEMENTED')
      ITRY=0
 1111 X=RANF(1.)
      DO 10 I=1,NBIN
      K=I
      IF (X .LT. BIN(I)) GOTO 1000
   10 CONTINUE
      WRITE(6,8900)X
 8900 FORMAT(' BIN SELECTION ERROR :',1P,E12.4)
      ITRY=ITRY+1
      IF (ITRY .LT. 3) GOTO 1111
      CALL ERREX( 'RECOVERY FAILED')
 1000 EN=(FLOAT(K)-0.5)*DE
      P=SQRT(2.*0.938*EN)
      CTH=1.-2.*RANF(1.)
      STH=1.-CTH*CTH
      PHI=PI2*RANF(1.)
      PX=P*STH*COS(PHI)
      PY=P*STH*SIN(PHI)
      PZ=P*CTH
      RETURN
C ------------------------------
      ENTRY PTINIT(T)
      IF (III .EQ. 0) CALL ERREX( 'PTHERM NOT IMPLEMENTED')
C ------------------------------
      EMAX=TMULT*MAX(T,TMIN)
      DE=EMAX/FLOAT(NBIN)
      CON=1.128379167*DE/T**1.5
      SUM=0.0
      DO 20 I=1,NBIN
      EN=(FLOAT(I)-0.5)*DE
      S=CON*SQRT(EN)*EXP(-EN/T)
      SUM=SUM+S
      BIN(I)=SUM
   20 CONTINUE
      IF (SUM .GT. 1.0) THEN
        S=1./SUM
        WRITE(6,8000)S
 8000   FORMAT(' PTHERM -- ADJUSTMENT MULT. : ',1P,E12.4)
        DO 100 I=1,NBIN
          BIN(I)=BIN(I)*S
 100    CONTINUE
      ELSE IF (SUM .LT. 1.0) THEN
        S=1.0-BIN(NBIN)
        WRITE(6,8010)S
 8010   FORMAT(' PTHERM -- ADJUSTMENT LAST CELL ADD : ',1P,E12.4)
        BIN(NBIN)=1.0
      END IF
      RETURN
      END
      SUBROUTINE SHIFT(IND,X0,Y0,Z0,PZ0,N)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C SHIFT THE N NUCLEONS (X,Y,Z,EM) SUCH THAT THEIR CENTER OF
C MASS IS AT (X0,Y0,Z0).
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C CREATED 28-JAN-83    HANS KRUSE
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C 9/19/83 ADDED SHIFTS OF MOMENTA
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   3.006 -- 25-JAN-85,14:03:00'/
      DATA PRONAM/'MODULE: S H I F T '/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (TTVERS) THEN
        WRITE(6,*)'---------------- VERSION NUMBER REPORT -------------'
        WRITE(6,*)PRONAM
        WRITE(6,*)VERNUM
        WRITE(6,*)'----------------------------------------------------'
        TTVERS=.FALSE.
      END IF
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CX=0.
      CY=0.
      CZ=0.
      S=0.
      CPX=0.
      CPY=0.
      CPZ=0.
      DO 100 IE=1,N
        I=IND+IE
        CX=CX+X(I)*EM(I)
        CY=CY+Y(I)*EM(I)
        CZ=CZ+Z(I)*EM(I)
        S=S+EM(I)
        CPX=CPX+PX(I)
        CPY=CPY+PY(I)
        CPZ=CPZ+PZ(I)
 100  CONTINUE
      CX=X0-CX/S
      CY=Y0-CY/S
      CZ=Z0-CZ/S
      CPX=-CPX/FLOAT(N)
      CPY=-CPY/FLOAT(N)
      CPZ=PZ0-CPZ/FLOAT(N)
      DO 101 IE=1,N
        I=IND+IE
        X(I)=X(I)+CX
        Y(I)=Y(I)+CY
        Z(I)=Z(I)+CZ
        PX(I)=PX(I)+CPX
        PY(I)=PY(I)+CPY
        PZ(I)=PZ(I)+CPZ
 101  CONTINUE
      RETURN
      END
C************************ V U U C O L L ***************************
C********************   A B S O R B   *********************************
C
      SUBROUTINE ABSORB(I1,I2,S,ICLTYP)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
C  ... PERFORM ABSORPTION OF A PION ON A NUCLEON TO MAKE A DELTA
C
C      I1,I2   : PARTICLE INDICES
C      S       : CM ENERGY
C      ICLTYP  : COLLISION TYPE
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 16-JAN-85,19:13:00'/
      DATA PRONAM/'A B S O R B  -- PION ABSORPTION'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
      IF (ITYP(I1) .LE. 6) THEN
        ID=I1
        IPI=I2
      ELSE
        ID=I2
        IPI=I1
      END IF
C
      EM(ID)=S
      PX(ID)=PX(ID)+PX(IPI)
      PY(ID)=PY(ID)+PY(IPI)
      PZ(ID)=PZ(ID)+PZ(IPI)
C
C  ...  RECOMPUTE THE ACCELERATION
C
      IF (MFIELD) CALL ACCEL(ID)
C
      IF (ICLTYP .EQ. 22) ITYP(ID)=3
      IF ( (ICLTYP .EQ. 23) .OR. (ICLTYP .EQ. 29) ) ITYP(ID)=4
      IF ( (ICLTYP .EQ. 30) .OR. (ICLTYP .EQ. 37) ) ITYP(ID)=5
      IF (ICLTYP .EQ. 38) ITYP(ID)=6
      DECAY(ID)=TIME+DECTIM(EM(ID),PX(ID),PY(ID),PZ(ID))
      LASTCL(ID)=0
C
      ITYP(IPI)=-9999
C
      ICOUNT(2,ICLTYP,KSIMUL)=ICOUNT(2,ICLTYP,KSIMUL)+1
C
      RETURN
      END
C*******************   C O L L  ***************************************
C*******************   C O L L  L O A D    ****************************
C
      SUBROUTINE COLOAD
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
C  ...  LOAD INITIAL VALUES INTO COLLISION PREDICTOR ARRAYS
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      LOGICAL HIT
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 16-JAN-85,19:18:25'/
      DATA PRONAM /'COLOAD -- COLLISION ARRAY REBUILD'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
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
C  ...  LOOP OVER ALL PARTICLE PAIRS AND STORE ALL PREDICTED COLLISIONS,
C       ALL TIMES ARE IN THE COMP. FRAME.
C
      NUMENT=0
      NTOTAL=NT+NPION(KSIMUL)
C ..............................FIND COLLISIONS INVOLVING PARTICLE I1
      DO 21 I1E=2,NTOTAL
        I1=I1E+IOFSET
        IP=I1E-1
C ..............................FIND COLLISIONS BETWEEN I1 AND I2
        DO 22 I2E=1,IP
          I2=I2E+IOFSET
          IF (HIT(I1,I2,SRT,SIG,SIGEL,THIT)) THEN
            NUMENT=NUMENT+1
            IF (NUMENT .GT. NUMEMA)
     *        CALL ERREX( 'COLOAD, COLLISION STORAGE EXCEEDED')
            I1C(NUMENT)=I1
            I2C(NUMENT)=I2
            CTIME(NUMENT)=THIT
            CS(NUMENT)=SRT
            CSIG(NUMENT)=SIG
            CSIGEL(NUMENT)=SIGEL
          END IF
   22   CONTINUE
   21 CONTINUE
      RETURN
      END
C
C*******************   I M P A C T  *******************************
C
C
      SUBROUTINE IMPACT(I1,I2,E1,E2,BSQ,TMIN,TMAX,GOTEST)
C
C  ...  RETURNS THE IMPACT PARAMETER COMPUTED BY LBL CODE
C  ...  TMIN IS RELATIVE TO  TIME
C
C%%%%%%%%%%%%%%%%%%%%%%  MODIFICATIONS  %%%%%%%%%%%%%%%%%%%%%
C
C HK 12/15/83 ADDED ALTERNATE RETURN FOR TMIN OUT-OF-RANGE
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
      PARAMETER (A2SMAL = .04)
      LOGICAL QUBIC,GOTEST
C
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 16-JAN-85,19:28:15'/
      DATA PRONAM /'IMPACT -- IMP. PARAM. AND COLL. TIME'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (TTVERS) THEN
        WRITE(6,*)'---------------- VERSION NUMBER REPORT -------------'
        WRITE(6,*)PRONAM
        WRITE(6,*)VERNUM
        WRITE(6,*)'----------------------------------------------------'
        TTVERS=.FALSE.
      END IF
      GOTEST=.TRUE.
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DBX=PX(I1)/E1-PX(I2)/E2
      DBY=PY(I1)/E1-PY(I2)/E2
      DBZ=PZ(I1)/E1-PZ(I2)/E2
      DX=X(I1)-X(I2)
      DY=Y(I1)-Y(I2)
      DZ=Z(I1)-Z(I2)
      DDB=DBX*DX+DBY*DY+DBZ*DZ
C
C  ...  DECIDE WHETHER THE ACCELERATION IS LARGE ENOUGH TO REQUIRE
C       THE FULL TRAJECTORY CALCULATION.
C
      IF (MFIELD) THEN
        DAX=AX(I1)-AX(I2)
        DAY=AY(I1)-AY(I2)
        DAZ=AZ(I1)-AZ(I2)
        DA2=DAX*DAX+DAY*DAY+DAZ*DAZ
        QUBIC = (DA2 .GT. A2SMAL)
      ELSE
        DA2=0.0
        QUBIC =.FALSE.
      END IF
C
C     QUBIC = (DA2 .GT. A2SMAL)
C
C
      IF (QUBIC) THEN
C   MUST DO FULL CALCULATION
        DB2=DBX*DBX+DBY*DBY+DBZ*DBZ
        DBDA=DBX*DAX+DBY*DAY+DBZ*DAZ
        DDA=DX*DAX+DY*DAY+DZ*DAZ
        DA2I=1./DA2
        AA3=DBDA*DA2I
        AA=3.*AA3
        BB=2.*(DB2+DDA)*DA2I
        CC=2.*DDB*DA2I
        PP=-AA*AA3+BB
        QQ=AA3*(2.*AA3*AA3-BB)+CC
        Q=(.333333*PP)**3+(.5*QQ)**2
        IF (Q) 300,200,100
C               SELECT CASE OF 3,2,1 REAL ROOTS
  300     IF (PP .GE. 0.0)
     &               CALL ERREX('IMPACT, LABEL 300, PP>0 @ Q<0')
          PP3=PP/3.
          CAL=-QQ*.5/SQRT(-PP3*PP3*PP3)
          AL3=ACOS(CAL)/3.
C
          AL3P=AL3+1.047198
          AL3M=AL3-1.047198
          FAC=2.*SQRT(-PP3)
          T1=COS(AL3)*FAC-AA3
          T2=-COS(AL3P)*FAC-AA3
          T3=-COS(AL3M)*FAC-AA3
        GOTO 1000
  200     YY=-SIGN((.5*ABS(QQ))**.333333,QQ)
          T1=2.*YY-AA3
          T2=-YY-AA3
          T3=1.E35
        GOTO 1000
  100     SQ=SQRT(Q)
          R1=-.5*QQ+SQ
          R2=-.5*QQ-SQ
          YY=SIGN(ABS(R1)**.333333,R1)
          ZZ=SIGN(ABS(R2)**.333333,R2)
          T1=YY+ZZ-AA3
          T2=1.E35
          T3=1.E35
C       ...  SOLUTIONS TO THE QUBIC EQUATIONS HAVE NOW BEEN FOUND
C
 1000   CONTINUE
C
C
        IF (T1 .LT. 0.0) T1=1.E35
        IF (T2 .LT. 0.0) T2=1.E35
        IF (T3 .LT. 0.0) T3=1.E35
        TMIN=AMIN1(T1,T2,T3)
        IF (TMIN .GT. TMAX) RETURN
C
C  ...  FIND THE VELOCITY AT THE FLY-BY AND SET UP FOR THE DISTANCE CALC.
C
        DXN=DX+TMIN*(DBX+.5*TMIN*DAX)
        DYN=DY+TMIN*(DBY+.5*TMIN*DAY)
        DZN=DZ+TMIN*(DBZ+.5*TMIN*DAZ)
C       CHANGE AX TO FX
        PX1N=PX(I1)+TMIN*FX(I1)
        PY1N=PY(I1)+TMIN*FY(I1)
        PZ1N=PZ(I1)+TMIN*FZ(I1)
        PX2N=PX(I2)+TMIN*FX(I2)
        PY2N=PY(I2)+TMIN*FY(I2)
        PZ2N=PZ(I2)+TMIN*FZ(I2)
        E1N=SQRT(EM(I1)*EM(I1)+PX1N*PX1N+PY1N*PY1N+PZ1N*PZ1N)
        E2N=SQRT(EM(I2)*EM(I2)+PX2N*PX2N+PY2N*PY2N+PZ2N*PZ2N)
        DENOM=(E1N+E2N)
        BX=(PX1N+PX2N)/DENOM
        BY=(PY1N+PY2N)/DENOM
        BZ=(PZ1N+PZ2N)/DENOM
      ELSE
C
        IF (DDB .GT. 0) RETURN
C
        DB2=DBX*DBX+DBY*DBY+DBZ*DBZ
        TMIN=-DDB/DB2
C
        IF (TMIN .GT. TMAX) RETURN
C
        BX=(PX(I1)+PX(I2))/(E1+E2)
        BY=(PY(I1)+PY(I2))/(E1+E2)
        BZ=(PZ(I1)+PZ(I2))/(E1+E2)
        DXN=DX+TMIN*DBX
        DYN=DY+TMIN*DBY
        DZN=DZ+TMIN*DBZ
      END IF
C
C  ...  FIND THE SQUARE OF THE DISTANCE OF CLOSEST APPROACH, IN THE
C       CM FRAME OF THE TWO PARTICLES AT THE TIME OF THE FLY-BY.
C
      GAM2=1./(1.-BX*BX-BY*BY-BZ*BZ)
      DN2=DXN*DXN+DYN*DYN+DZN*DZN
      BDN=BX*DXN+BY*DYN+BZ*DZN
C
      BSQ=DN2+GAM2*BDN*BDN
C
      GOTEST=.FALSE.
      RETURN
      END
C*********************  H I T  *************************************
C
      LOGICAL FUNCTION HIT(I1,I2,SRT,SIG,SIGEL,THIT)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      LOGICAL   GOTEST
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 16-JAN-85,19:35:50'/
      DATA PRONAM /'HIT -- DETERMINE WHETHER COLLISION OCCURS'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C      DIMENSION XMESSAGE(10)
C      LOGICAL PROXIM
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (TTVERS) THEN
        WRITE(6,*)'---------------- VERSION NUMBER REPORT -------------'
        WRITE(6,*)PRONAM
        WRITE(6,*)VERNUM
        WRITE(6,*)'----------------------------------------------------'
        TTVERS=.FALSE.
C        SET SOME DEFAULT VALUES
         IF(DA(8).GT.0.0001)THEN
C          CROSS SECTION FACTOR  FOR  10.*PI*BSQ .GT. DA(8)*SIG
           CRSFCT=10.*PI/DA(8)
           WRITE(6,*) ' CROSS SECTION FACTOR OF ',DA(8)
         ELSE
           CALL ERREX(' ZERO X SECTION IN FUNCTION HIT')
         END IF
      END IF
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      IF (I1 .EQ. I2) GOTO 22
C      AVOID COLLISIONS WITH PIONS WHICH ARE ALREADY ABSORBED
      IF (ITYP(I1) .LT. 0) GOTO 22
      IF (ITYP(I2) .LT. 0) GOTO 22
C ...  AVOID SECOND COLLISION BETWEEN TWO PARTICLES THAT JUST COLLIDED
C
      IF ( (LASTCL(I1) .EQ. I2) .AND. (LASTCL(I2) .EQ. I1) ) GOTO 22
C
C
C ..   ALSO, DONT LET PION COLLIDE WITH ITS DECAY NUCLEON
C
        IF ((ITYP(I1) .GE. 7) .AND. (LASTCL(I1) .EQ. I2)) GOTO 22
        IF ((ITYP(I2) .GE. 7) .AND. (LASTCL(I2) .EQ. I1)) GOTO 22
C
C ...  IDENTIFY TYPE OF COLLISION
C
      III=MAX0(ITYP(I1),ITYP(I2))
      JJJ=MIN0(ITYP(I1),ITYP(I2))
      ICLTYP=(III*(III-1))/2+JJJ
      IF ( (ICLTYP .LT. 0) .OR. (ICLTYP .GT. MAXTYP) ) THEN
        WRITE(6,*)'ICLTYP,1,2',ICLTYP,I1,ITYP(I1),I2,ITYP(I2)
        GOTO 22
      END IF
C
C  ...  THIS IS A TEMPORARY SAVINGS, SINCE CURRENTLY ONLY INELASTIC
C       PI-N COLLISIONS ARE  INCLUDED (VERSION 3)
C
      IF ( (ICLTYP .GT. 21) .AND.
     &     (ICLTYP .NE. 22) .AND.
     &     (ICLTYP .NE. 23) .AND.
     &     (ICLTYP .NE. 29) .AND.
     &     (ICLTYP .NE. 30) .AND.
     &     (ICLTYP .NE. 37) .AND.
     &     (ICLTYP .NE. 38)       )
     &                               GOTO 22
C
      XM=SEACUT(ICLTYP)
C
C ...  PRESELCT AND STORE PROPERTIES IN SCALARS
C
      X1=X(I1)
      Y1=Y(I1)
      Z1=Z(I1)
      X2=X(I2)
      IF(ABS(X1-X2).GT.XM) GO TO 22
      Y2=Y(I2)
      IF(ABS(Y1-Y2).GT.XM) GO TO 22
      Z2=Z(I2)
      IF(ABS(Z1-Z2).GT.XM) GO TO 22
      R2=(X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2
      XM2=XM*XM
      IF (R2 .GT. XM2) GO TO 22
      PZ1=PZ(I1)
      PZ2=PZ(I2)
      PX1=PX(I1)
      PX2=PX(I2)
      PY1=PY(I1)
      PY2=PY(I2)
C
C  ...  4-COMPONENT OF THE MOMENTUM VECTOR
C
      E1=SQRT(EM(I1)**2+PZ1**2+PX1**2+PY1**2)
      E2=SQRT(EM(I2)**2+PZ2**2+PX2**2+PY2**2)
C
C  ...   FERMI  BLOCK, I.E. LOW-ENERGY CUTOFF
C
C
C  ...  TOTAL ENERGY SQUARED
C
      S=(E1+E2)**2-(PZ1+PZ2)**2-(PX1+PX2)**2-(PY1+PY2)**2
      IF (S .LT. SM2CUT(ICLTYP) ) THEN
        XBLCUT=XBLCUT+1
        GOTO 22
      END IF
C
C  ...  GET TIME OF CLOSEST APPROACH AND IMPACT PARAMETER-SQUARED
C
      CALL IMPACT(I1,I2,E1,E2,BSQ,TMIN,TIMMAX-TIME,GOTEST)
      IF (GOTEST) GOTO 22
      IF (BSQ .GT. XM2) GOTO 22
C
C  ...  FOR ALL COLLISIONS WHICH DO NOT NEED DETAILED BALANCE,
C       CHECK CROSSECTION NOW
C
      IF ( (ICLTYP .NE. 5) .AND. (ICLTYP .NE. 7) .AND.
     *     (ICLTYP .NE. 8) .AND. (ICLTYP .NE. 11) .AND.
     *     (ICLTYP .NE. 12) .AND. (ICLTYP .NE. 16)) THEN
        PR=-999.
C
C  ...  GET CROSSECTION
C
        SRT=SQRT(S)
        CALL CROSS1(SRT,PR,ICLTYP,SIG,SIGEL)
        IF (CRSFCT*BSQ .GT. SIG) GOTO 22
C  HIT OUTSIDE CROSSECTION
      ELSE
C
C  ...  VELOCITY OF TWO-BODY SYSTEM
C
        BX=(PX1+PX2)/(E1+E2)
        BY=(PY1+PY2)/(E1+E2)
        BZ=(PZ1+PZ2)/(E1+E2)
        GAM=1./SQRT(1.-BX**2-BY**2-BZ**2)
        B2=.0001+BX**2+BY**2+BZ**2
        PB=BX*PX1+BY*PY1+BZ*PZ1
        GG=(GAM-1.)*PB/B2-GAM*E1
C
C  ...  RELATIVE MOMENTUM
C
        PZR=GG*BZ+PZ1
        PXR=GG*BX+PX1
        PYR=GG*BY+PY1
        PR=SQRT(PXR**2+PYR**2+PZR**2)
C
C  ...  GET CROSSECTION
C
        SRT=SQRT(S)
        CALL CROSS1(SRT,PR,ICLTYP,SIG,SIGEL)
        IF (CRSFCT*BSQ .GT. SIG) GOTO 22
C  HIT OUTSIDE CROSSECTION
      END IF
C
C  ...  COLLISION WILL OCCUR:
C
      HIT=.TRUE.
      THIT=TIME+TMIN
      RETURN
C
C  ...  COLLISION WILL NOT OCCUR
C
   22 HIT=.FALSE.
      RETURN
      END
C***************   F I N D  N E X T   ************************
C
      SUBROUTINE FNDNXT(INDEX,WHAT)
C
C  ...  SEARCH DECAY AND COLLISION ARRAYS TO FIND NEXT EVENT
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
      CHARACTER*1 WHAT
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 19-JAN-85,19:42:20'/
      DATA PRONAM/'FNDNXT'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
      TMIN=TIMMAX
      IMIN=0
      IF (NUMENT .LE. 0) GOTO 2000
C
C  ...  SEARCH COLLISION ARRAY
C
      DO 100 I=1,NUMENT
        IF (I1C(I) .GT. 0) THEN
          IF (CTIME(I) .LE. TMIN) THEN
            TMIN=CTIME(I)
            IMIN=I
          END IF
        END IF
  100 CONTINUE
C
C  ...  SEARCH DECAY ARRAY FOR NEXT DECAY EVENT
C     SWITCH(2) KEEPS THE DELTAS STABLE
C
 2000 TDEC=TIMMAX
      JMIN=0
      IF(.NOT.SWITCH(2))THEN
      DO 101 JE=1,NT
        J=JE+IOFSET
        IF (DECAY(J) .LE. TDEC) THEN
          TDEC=DECAY(J)
          JMIN=J
        END IF
  101 CONTINUE
      END IF
C
      IF (TMIN .LE. TDEC) THEN
        INDEX=IMIN
        WHAT='C'
      ELSE
        INDEX=JMIN
        WHAT='D'
      END IF
      IF (INDEX .EQ. 0) WHAT='E'
      RETURN
      END
C ----------------------  CROSES.FOR  ----------------------------
C
      BLOCK DATA
C
C CROSSECTIONS FOR THE CASCADE CODE  VUU/SCATTER ; NOTE THAT THE NN
C
C
C CROSSECTIONS WILL BE EQUATED TO THE PP ONES USING  CINIT
C
        INCLUDE  (VUUSIG)
C
C PN - TOTAL -- SMOOTHED AND EXTRAPOLATED
C
      DATA SPNTOT/
     @ 44.50, 43.35, 42.16, 40.96, 39.73, 38.48, 37.20, 35.98,
     @ 35.02, 34.47, 34.37, 34.67, 35.23, 35.97, 36.75, 37.37,
     @ 37.77, 38.03, 38.40, 38.83, 39.26, 39.67, 40.06, 40.45,
     @ 40.79, 41.06, 41.31, 41.52, 41.70, 41.81, 41.87, 41.98,
     @ 42.12, 42.29, 42.55, 42.82, 43.01, 43.12, 43.16, 43.14,
     @ 43.06, 42.95, 42.81, 42.67, 42.54, 42.45, 42.38, 42.33,
     @ 42.30, 42.29, 42.28, 42.26, 42.24, 42.21, 42.17, 42.14,
     @ 42.10, 42.07, 42.06, 42.05, 42.04, 42.03, 42.02, 42.00,
     @ 41.97, 41.94, 41.89, 41.84, 41.79, 41.73, 41.67, 41.61,
     @ 41.55, 41.49, 41.44, 41.38, 41.34, 41.31, 41.29, 41.28,
     @ 41.27, 41.28, 41.30, 41.33, 41.36, 41.40, 41.44, 41.49
     @/
C
C  PN - ELASTIC FROM SUBTRACTION, TOTAL-SCALED PP->PPPI0 DATA
C
      DATA SPNEL /
     @ 44.50, 43.35, 42.16, 40.96, 39.73, 38.48, 37.20, 35.98,
     @ 35.02, 34.30, 32.77, 31.58, 30.61, 29.77, 29.04, 28.38,
     @ 27.76, 27.25, 27.05, 27.10, 27.26, 27.48, 27.75, 28.08,
     @ 28.41, 28.72, 29.03, 29.32, 29.59, 29.77, 29.89, 30.06,
     @ 30.27, 30.51, 30.86, 31.24, 31.59, 31.91, 32.18, 32.43,
     @ 32.65, 32.85, 33.05, 33.25, 33.46, 33.70, 33.94, 34.17,
     @ 34.39, 34.57, 34.70, 34.78, 34.80, 34.77, 34.71, 34.60,
     @ 34.47, 34.32, 34.17, 34.03, 33.88, 33.76, 33.65, 33.56,
     @ 33.50, 33.48, 33.50, 33.60, 33.83, 34.27, 34.91, 35.56,
     @ 36.03, 36.31, 36.41, 36.40, 36.29, 36.16, 36.04, 35.97,
     @ 35.97, 36.02, 36.12, 36.27, 36.47, 36.69, 36.95, 37.22
     @/
C
C PP-TOTAL SMOOTHED
C
      DATA SPPTOT/
     @ 39.48, 31.76, 26.26, 24.05, 23.94, 23.77, 23.74, 24.13,
     @ 25.09, 26.23, 28.04, 32.30, 34.93, 37.85, 41.47, 44.73,
     @ 46.73, 47.33, 47.29, 47.49, 47.50, 47.37, 47.54, 47.39,
     @ 47.38, 47.38, 47.22, 47.14, 46.91, 46.60, 46.39, 45.80,
     @ 45.47, 45.54, 45.28, 45.00, 44.95, 44.74, 44.32, 44.21,
     @ 43.86, 43.63, 43.50, 43.04, 43.01, 42.79, 42.65, 42.44,
     @ 42.20, 42.02, 41.87, 41.75, 41.64, 41.53, 41.43, 41.33,
     @ 41.24, 41.16, 41.09, 41.03, 40.97, 40.90, 40.83, 40.75,
     @ 40.67, 40.59, 40.51, 40.45, 40.39, 40.34, 40.29, 40.25,
     @ 40.22, 40.19, 40.17, 40.15, 40.13, 40.11, 40.10, 40.08,
     @ 40.06, 40.04, 40.02, 40.00, 39.98, 39.95, 39.93, 39.91
     @/
C
C PP - ELASTIC -- SMOOTHED
C
      DATA SPPEL /
     @ 41.59, 35.72, 30.85, 27.18, 24.72, 23.29, 22.65, 22.57,
     @ 22.83, 23.26, 23.74, 24.22, 24.58, 24.75, 24.74, 24.61,
     @ 24.43, 24.29, 24.21, 24.10, 23.86, 23.48, 23.01, 22.49,
     @ 21.94, 21.39, 20.87, 20.42, 20.09, 19.92, 19.88, 19.92,
     @ 19.95, 19.91, 19.70, 19.25, 18.73, 18.27, 17.85, 17.46,
     @ 17.09, 16.73, 16.37, 15.99, 15.58, 15.13, 14.63, 14.11,
     @ 13.63, 13.24, 12.94, 12.74, 12.62, 12.56, 12.56, 12.59,
     @ 12.64, 12.69, 12.71, 12.69, 12.64, 12.56, 12.46, 12.33,
     @ 12.19, 12.04, 11.88, 11.73, 11.58, 11.46, 11.36, 11.28,
     @ 11.22, 11.16, 11.11, 11.06, 10.99, 10.90, 10.77, 10.61,
     @ 10.42, 10.25, 10.10, 10.02, 10.00, 10.03, 10.09, 10.17
     @/
C
C PP - INELASTIC, FROM PI0 DATA, SMOOTHED AND EXTRAPOLATED
C
      DATA SNNND /
     @  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @  0.00,  0.35,  3.21,  6.17,  9.25, 12.41, 15.41, 17.98,
     @ 20.03, 21.58, 22.69, 23.46, 24.00, 24.38, 24.62, 24.74,
     @ 24.76, 24.69, 24.56, 24.40, 24.23, 24.08, 23.96, 23.84,
     @ 23.71, 23.57, 23.39, 23.15, 22.83, 22.44, 21.97, 21.43,
     @ 20.84, 20.20, 19.53, 18.85, 18.16, 17.50, 16.88, 16.32,
     @ 15.83, 15.44, 15.16, 14.97, 14.88, 14.87, 14.94, 15.07,
     @ 15.26, 15.50, 15.77, 16.04, 16.32, 16.56, 16.75, 16.88,
     @ 16.94, 16.91, 16.79, 16.49, 15.92, 14.93, 13.53, 12.10,
     @ 11.04, 10.37, 10.04,  9.98, 10.09, 10.30, 10.50, 10.61,
     @ 10.62, 10.53, 10.36, 10.11,  9.79,  9.41,  8.98,  8.53
     @/
C
C   P+PI+ -> DELTA, SIGMAX =   200.0000
C
      DATA SNPID /
     @  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, 32.58,
     @ 42.46, 57.10, 79.25,112.34,156.30,194.35,192.95,152.69,
     @107.71, 74.46, 52.55, 38.28, 28.77, 22.24, 17.60, 14.22,
     @ 11.68,  9.74,  8.23,  7.03,  6.06,  5.28,  4.62,  4.08,
     @  3.62,  3.23,  2.90,  2.62,  2.37,  2.15,  1.96,  1.79,
     @  1.64,  1.51,  1.39,  1.29,  1.19,  1.11,  1.03,  0.96,
     @  0.90,  0.84,  0.78,  0.73,  0.69,  0.65,  0.61,  0.58,
     @  0.54,  0.51,  0.48,  0.46,  0.43,  0.41,  0.39,  0.37,
     @  0.35,  0.34,  0.32,  0.30,  0.29,  0.28,  0.26,  0.25,
     @  0.24,  0.23,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
     @  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00
     @/
      END
C*******************    C R O S S    *******************************
C*******************   C R O S S 1   *******************************
C
C   UPDATED 7/26/83 BVJ.
C   CHANGES:      -- CALL  INCLUDE FILE TO DEFINE CROSS SECTION
C   COMMON BLOCK.
C       -- MAKE CUTOFFS, OFFSETS INTO TABLES, AND TABLE
C          LOOKUPS DEPEND ON TABLE PARAMETERS.
C
C   UPDATED VERSION 9/6/82.    HK
C
C   CHANGES  INCLUDE: -- END-OF-TABLE ENERGY IS 19.11 GEV
C                    -- A SLIGHT CHANGE IN CORRECTING OUT-OF-BOUNDS
C                       INDEX VALUES.
C                    -- A CHANGE IN COMPUTING THE ELASTIC CROSSECTION
C                       FOR ENERGIES BEYOND THE END-OF-TABLE VALUE.
C                    -- THE USE OF OFFSETS AND ASSIGNED GOTO S FOR
C                       THE TABLE LOOK-UP.
C                    -- DELTA-NUCLEON CROSSECTIONS USING DETAILED
C                       BALANCE FOR THE INELASTIC PORTION. THIS
C                       REQUIRES THE USE OF THE RELATIVE MOMENTUM  PR
C
C   CHANGES 8/11/83  HK:
C                    -- ADDED COINEL TO CORRECTLY DETERMINE
C                       THE INELASTIC CROSSECTION NEEDED IN THE
C                       DETAILED BALANCE.
C
C   CHGED 8/18/83 HK: UPDATE TO VERSION 3
C                       NEW CROSSECTION COMMON BLOCKS
C                       PION ABSORPTION CROSSECTIONS
C                       THESE REQUIRE A SHIFT PARAMETER TO USE THE
C                       LOOKUP SEQUENCE
C
      SUBROUTINE CROSS1(S,PR,ICLTYP,SIG,SIGEL)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUSIG)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      LOGICAL COINEL
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 16-JAN-85,20:01:30'/
      DATA PRONAM /'CROSS1 -- TOTAL AND ELASTIC CROSSECTIONS'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C      COMPUTE TOTAL CROSSECTION AND ELASTIC PARTIAL CROSSECTION
C      FOR COLLISION TYPE  ICLTYP  AT C.M. ENERGY S.
C
C       S      : C.M. ENERGY                    (INPUT)
C       PR     : RELATIVE MOMENTUM              (INPUT)
C       ICLTYP : COLLISION TYPE                 (INPUT)
C       SIG    : TOTAL CROSSECTION              (OUTPUT)
C       SIGEL  : ELASTIC CROSSECTION            (OUTPUT)
C
C  ...  COINEL CONTROLS WHETHER THE TABLE LOOKUP COMPUTES THE
C       INELASTIC CROSSECTION (TRUE), OR TOTAL AND ELASTIC (FALSE)
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
      COINEL=.FALSE.
      TBLSFT=0.0
      GOTO (
     >  10, 20, 30, 90, 50,
     >  60, 50, 50, 90, 60,
     >  50, 50, 90, 90, 60,
     >  50, 90, 90, 90, 90,
     >  60, 80, 80, 70, 70,
     >  70, 70, 70, 80, 80,
     >  70, 70, 70, 70, 70,
     >  70, 80, 80, 70, 70,
     >  70, 70, 70, 70, 70
     & ) , ICLTYP
C
      WRITE(6,*)' ICLTYP INVALID IN  CROSS1   ',ICLTYP
      CALL ERREX( 'FATAL ERROR')
C  N-N COLLISION
 10   CONTINUE
      IF(S.GE.TBLHIG)THEN
        SIG=0.0645*(S-19.11)+40.
        SIGEL=0.0215*(S-19.11)+6.7
        RETURN
      END IF
      IOFSET=0
      ASSIGN 100 TO LABEL
      GOTO 1000
  100 RETURN
C
C  P-N COLLISION
C
 20   CONTINUE
      IF(S.GE.TBLHIG)THEN
        SIG=0.0645*(S-19.11)+40.
        SIGEL=0.0
        RETURN
      END IF
      IOFSET=TBLSIZ * 2
      ASSIGN 200 TO LABEL
      GOTO 1000
  200 RETURN
C
C  P-P COLLISION
C
 30   CONTINUE
      IF(S.GE.TBLHIG)THEN
        SIG=0.0645*(S-19.11)+40.
        SIGEL=0.0215*(S-19.11)+6.7
        RETURN
      END IF
      IOFSET=TBLSIZ * 4
      ASSIGN 300 TO LABEL
      GOTO 1000
  300 RETURN
C
C DELTA-NUCLEON SCATTERING WITH POSSIBLE DELTA ABSORPTION
C
   50 CONTINUE
C
C ...  THE CUT-OFF IN THE FOLLOWING STATEMENT SHOULD REFLECT THE
C      LARGEST CM ENERGY AT WHICH DELTA-DOMINANCE IS STILL VALID,
C      SINCE WE WILL DETERMINE THE DELTA ABSORPTION FROM THE
C      INELASTIC NUCLEON-NUCLEON CROSSECTION.
      IF (S .GT. DLTLIM) GOTO 90
C
C ...  ELASTIC CROSSECTION FROM P-N
C
      IOFSET=TBLSIZ * 2
      ASSIGN 500 TO LABEL
      GOTO 1000
  500 SIGSAV=SIGEL
C
C ...  EXIT CHANNEL MOMENTUM FOR DETAILED BALANCE
C
      PR2NEW=.25*S*S-.880
      IOFSET=TBLSIZ * 6
      IF (ICLTYP.EQ.5) THEN
        WEIGHT=.75
      ELSE IF (ICLTYP.EQ.7) THEN
        WEIGHT=.25
      ELSE IF ((ICLTYP.EQ.11).OR.(ICLTYP.EQ.8)) THEN
        WEIGHT=.5
      ELSE IF (ICLTYP.EQ.16) THEN
        WEIGHT=.75
      ELSE IF (ICLTYP.EQ.12) THEN
        WEIGHT=.25
      ELSE
        CALL ERREX( 'FATAL ERROR -- SUBSELECT, CROSS1, AT 50')
      END IF
      COINEL=.TRUE.
      ASSIGN 510 TO LABEL
      GOTO 1000

  510 SIGIN=WEIGHT*SIGIN
C   SIGIN WAS COMPUTED DUE TO COINEL
      COINEL=.FALSE.
      SIGEL=SIGSAV
      SIG=SIGEL
C
C ... CHECK THAT  PR  IS VALID
C
      IF (PR .LT. 0) CALL ERREX(' CROSS1, AT 510, PR INVALID')
      SIG=SIGEL+.5*PR2NEW*SIGIN/(PR*PR)
      RETURN
C
C IDENTICAL DELTA-DELTA SCATTERING, ELASTIC ONLY.
C
   60 CONTINUE
      IF(S.GE.TBLHIG)THEN
        SIGEL=0.0215*(S-19.11)+6.7
        SIG=SIGEL
        RETURN
      END IF
      IOFSET=0
      ASSIGN 600 TO LABEL
      GOTO 1000
  600 SIG=SIGEL
      RETURN
C
C PION SCATTERING WITHOUT ABSORPTION POSSIBILITY, CURRENTLY 0
C
   70 CONTINUE
      SIG=0.0
      SIGEL=0.0
      RETURN
C
C PION ABSORPTION SCATTERING
C
   80 CONTINUE
      COINEL=.TRUE.
      TBLSFT=EMNUC
      IF ( (S+TBLSFT .LT. TBLLOW) .OR.
     *     (S+TBLSFT .GT. TBLHIG) ) THEN
        SIG=0.0
        SIGEL=0.0
        RETURN
      END IF
      ASSIGN 800 TO LABEL
      IOFSET = TBLSIZ * 7
      GOTO 1000
  800 COINEL=.FALSE.
      IF ( (ICLTYP .EQ. 22) .OR. (ICLTYP .EQ. 38) ) WEIGHT=1.0
      IF ( (ICLTYP .EQ. 37) .OR. (ICLTYP .EQ. 23) ) WEIGHT=0.33333333
      IF ( (ICLTYP .EQ. 29) .OR. (ICLTYP .EQ. 30) ) WEIGHT=0.66666666
      SIG=SIGIN*WEIGHT
      SIGEL=0.
      RETURN
C
C NON-IDENTICAL DELTA-DELTA SCATTERING AS WELL AS (D- + N), (D++ + P),
C               ELASTIC ONLY
   90 CONTINUE
      IF (S .GE. TBLHIG) THEN
        SIG=0.0
        SIGEL=0.0
        RETURN
      END IF
      IOFSET=TBLSIZ * 2
      ASSIGN 900 TO LABEL
      GOTO 1000
  900 SIG=SIGEL
      RETURN
C
C  ...  TABLE LOOK-UP SEQUENCE.  IOFSET   DETERMINES THE CROSSECTIONS
C       PICKED UP FROM THE COMMON BLOCK, THE SEGMENT RETURNS TO THE
C       LOCATION ASSIGNED TO  LABEL .
C  ...  IF  COINEL     IS TRUE, THE INELASTIC CROSSECTION IS RETURNED
C       IN  SIGIN , OTHERWISE THE TOTAL AND ELASTIC CROSSECTIONS ARE
C       FOUND AND RETURNED IN  SIG  AND  SIGEL .  HK
C
C  ...  THE TABLE INDEX IS FOUND AS LOG(CMENERGY + TBLSFT)
C
C  FIND INDEX OF ELEMENT BELOW ENERGY OF INTEREST
 1000 XPT=ALOG(S+TBLSFT)
      TLOW=ALOG(TBLLOW)
      MAX=INT(TBLSIZ)-1
      INDEX=INT((XPT-TLOW)/TBLSTP)+1
      IF(INDEX.LE.0)INDEX=1
      IF(INDEX.GT.MAX)INDEX=MAX
C  FIND SLOPES AND CROSSECTIONS
      INDTOT=IOFSET+INDEX
      INDEL=INDTOT+TBLSIZ
      INDIN=IOFSET+INDEX
C  THE INELASTIC CROSSECTION HAS ITS OWN OFFSET
      X1=(INDEX-1)*TBLSTP      + TLOW
C  CORR-HK 10/5/83
      X2=X1 + TBLSTP
      IF (COINEL) THEN
        Y1IN=SIGCOM(INDIN)
        Y2IN=SIGCOM(INDIN+1)
        SLIN=(Y2IN-Y1IN)/TBLSTP
C  UPDATED-HK 10/5/83
        SIGIN=SLIN*(XPT-X1) + Y1IN
      ELSE
        Y1TOT=SIGCOM(INDTOT)
        Y2TOT=SIGCOM(INDTOT+1)
        Y1EL=SIGCOM(INDEL)
        Y2EL=SIGCOM(INDEL+1)
        SLTOT=(Y2TOT-Y1TOT)/TBLSTP
C  UPDATED-HK 10/5/83
        SLEL=(Y2EL-Y1EL)/TBLSTP
C  UPDATED-HK 10/5/83
        SIG=SLTOT*(XPT-X1) + Y1TOT
        SIGEL=SLEL*(XPT-X1) + Y1EL
      END IF
C
      GOTO LABEL
C
      END
C
C*****************   A N G E L   **********************************
C
      SUBROUTINE ANGEL(S,PR,ICLTYP,T1,C1)
C
        INCLUDE (VUUPRM)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   4.006 -- 16-JAN-85,20:09:20'/
      DATA PRONAM/'MODULE: A N G E L'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C      DETERMINE FOR ELASTIC SCATTERING THE NEW MOMENTUM
C      DIRECTION (RELATIVE TO THE ORIGINAL DIRECTION). RETURN
C      THE AZIMUTHAL ANGLE, T1, AND THE COSINE OF THE POLAR
C      ANGLE, C1.
C
C        S      : C.M. ENERGY                      (INPUT)
C        PR     : MAGNITUDE OF REL. MOMENTUM       (INPUT)
C        ICLTYP : COLLSION TYPE                    (INPUT)
C        T1     : AZIMUTHAL ANGLE                  (OUTPUT)
C        C1     : COSINE(POLAR ANGLE)              (OUTPUT)
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C ERROR IN C1= . . .      CORRECTED  2/10/83  HK
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
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
      GOTO (
     >   1,  2,  1,  2,  2,
     >   1,  2,  2,  2,  1,
     >   2,  2,  2,  2,  1,
     >   2,  2,  2,  2,  2,
     >   1
     & ) , ICLTYP
C
      WRITE(6,*)' ICLTYP INVALID IN  ANGEL   ',ICLTYP
      CALL ERREX( 'FATAL ERROR')
C
C ...  ELASTIC SCATTERING OF IDENTICAL PARTICLES
C
    1 TA=-2.*PR*PR
  101 AS=(3.65*(S-1.8766))**6
C CORRECTED 9/13/83; HK
      A=6.*AS/(1.+AS)
      Y3=RANF(0)
      PT=ALOG((1.-Y3)*EXP(A*TA)+Y3)/A
      C1=1.+PT/(2*PR*PR)
      IF (ABS(C1) .GT. 1.0) C1=SIGN(1.0,C1)
      T1=PI2*RANF(0)
      RETURN
C
C ...  ELASTIC SCATTERING OF NON-IDENTICAL PARTICLES
C
    2 TA=-4.*PR*PR
      GOTO 101
      END
C
C****************    C R O S S 2   **************************************
C
C  CHANGES BY B.JACAK:
C                  -- TABLE LOOKUP GOVERNED BY THE PARAMETERS IN
C                      VUUSIG.FOR .
C
C  CHANGES 8/11/83 HK:
C                  -- NUCL.-NUCL. CROSSECTIONS EXEEDING THE
C                     LOOKUP TABLE ARE SET TO THE HIGHEST TABLE
C                     VALUE - THIS IS A TEMPORARY FIX!
C
C          8/18/83 -- LOOKUP CHGED TO NEW COMMON BLOCKS
C
      SUBROUTINE CROSS2(S,PR,EM1,EM2,ICLTYP,SIG,SIGEL,SIGIN,ICNUM)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUSIG)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      DIMENSION SIGIN(MXCHAN)
C<>VERNUM<>
      DATA VERNUM/'VER   4.006 -- 16-JAN-85,20:12:35'/
      DATA PRONAM/'MODULE: C R O S S 2'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C      DETERMINE THE NUMBER OF TOTAL CHANNELS (# OF INELASTIC + 1)
C      FOR THIS COLLISION AND RETURN THE (ICNUM-1) INELASTIC
C      CROSSECTIONS IN SIGIN.
C
C        S      : C.M. ENERGY                           (INPUT)
C        PR     : REL. MOMENTUM                         (INPUT)
C        EM1    : MASS PART. 1                          (INPUT)
C        EM2    : MASS PART. 2                          (INPUT)
C        ICLTYP : COLLISION TYPE                        (INPUT)
C        SIG    : TOTAL CROSSECTION FROM CROSS1         (INPUT)
C        SIGEL  : ELASTIC CROSSETION FROM CROSS1        (INPUT)
C        SIGIN  : (MXCHAN) INEL. CROSSECTIONS           (OUTPUT)
C        ICNUM  : # OF TOTAL EXIT CHANNELS              (OUTPUT)
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
      GOTO (
     >  10, 20, 30, 60, 50,
     >  60, 50, 80, 60, 60,
     >  80, 50, 60, 60, 60,
     >  50, 60, 60, 60, 60,
     >  60
     & ) , ICLTYP
C
      WRITE(6,*)' ICLTYP INVALID IN  CROSS2   ',ICLTYP
      CALL ERREX( 'FATAL ERROR')
C
C ---   N+N -> N+DELTA
C
C  N-N COLLISION
 10   CONTINUE
C
C ...  MAKE SURE THERE IS ENOUGH ENERGY TO MAKE A DELTA, S>2*M(N)+M(PI)
C      ALSO DO NOT ALLOW INELASTIC SCATTERING ABOVE THE DELTA DOMINANC E
C      REGION, THIS WILL BE LIFTED AS SOON AS PIONS ARE AVAILABLE.
C
      IF ( (S .LT. (2.*EMNUC+EMPION+EKINMI)) .OR.
     *     (S .GT. DLTLIM) ) THEN
        ICNUM=1
        RETURN
      END IF
        IF(S.GT.TBLHIG)THEN
        SINTOT=SNNND(INT(TBLSIZ+0.5))
        GOTO 100
      END IF
      IOFSET=TBLSIZ * 6
      ASSIGN 100 TO LABEL
      GOTO 1000
  100 ICNUM=3
      SIGIN(1)=.25*SINTOT
      SIGIN(2)=.75*SINTOT
      RETURN
C
C  P-N COLLISION
 20   CONTINUE
      IF ( (S .LT. (2.*EMNUC+EMPION+EKINMI)) .OR.
     *     (S .GT. DLTLIM) ) THEN
        ICNUM=1
        RETURN
      END IF
        IF(S.GT.TBLHIG)THEN
          SINTOT=SNNND(INT(TBLSIZ+0.5))
          GOTO 200
      END IF
      IOFSET=TBLSIZ * 6
      ASSIGN 200 TO LABEL
      GOTO 1000
  200 ICNUM=3
      SIGIN(1)=.25*SINTOT
C  CHGED FROM VERSION 2, THE FACTOR OF 2 WAS
C  HIDDEN IN THE DIFFERENT INEL.  CROSSECITONS
      SIGIN(2)=SIGIN(1)
      RETURN
C
C  P-P COLLISION
 30   CONTINUE
      IF ( (S .LT. (2.*EMNUC+EMPION+EKINMI)) .OR.
     *     (S .GT. DLTLIM) ) THEN
        ICNUM=1
        RETURN
      END IF
      IF(S.GT.TBLHIG)THEN
        SINTOT=SNNND(INT(TBLSIZ+0.5))
        GOTO 300
      END IF
      IOFSET=TBLSIZ * 6
      ASSIGN 300 TO LABEL
      GOTO 1000
  300 ICNUM=3
      SIGIN(1)=.25*SINTOT
      SIGIN(2)=.75*SINTOT
      RETURN
C
C ...  COMBINATIONS WITHOUT INELASTIC CHANNEL
C
   60 ICNUM=1
      RETURN
C
C ---  DELTA+N -> N+N   ONE EXIT CHANNEL
C
   50 IF (S .GT. DLTLIM) THEN
        ICNUM=1
        RETURN
      END IF
      ICNUM=2
      SIGIN(1)=SIG-SIGEL
C THIS IS CORRECT, SINCE CROSS1 COMPUTED SIG
C AS SIGEL+SIGIN.
      RETURN
C
C ---  DELTA+N -> N+N   TWO EXIT CHANNELS, DUMMY SINCE EXIT CHANNELS
C                       CURRENTLY CARRY THE SAME WEIGHT.
C
   80 GOTO 50
C
C ...  TABLE LOOK-UP (SEE CROSS1)
C
C  FIND INDEX OF ELEMENT BELOW ENERGY OF INTEREST
 1000 XPT=ALOG(S)
      TLOW=ALOG(TBLLOW)
      MAX=INT(TBLSIZ)-1
      INDEX=INT((XPT-TLOW)/TBLSTP)+1
      IF(INDEX.LE.0)INDEX=1
      IF(INDEX.GT.MAX)INDEX=MAX
C  FIND SLOPES AND CROSSECTIONS
      X1=(INDEX-1)*TBLSTP + TLOW
C  CORR-HK 10/5/83
      X2=X1+TBLSTP
      INDEX=INDEX+IOFSET
C INELASTIC CROSSECTION HAS ITS OWN OFFSET
      Y1IN=SIGCOM(INDEX)
      Y2IN=SIGCOM(INDEX+1)
      SLIN=(Y2IN-Y1IN)/TBLSTP
C  UPDATED-HK 10/5/83
      SINTOT=SLIN*(XPT-X1) + Y1IN
C
      GOTO LABEL
C
      END
C
C*****************   A N G I N   *********************************
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CURRENTLY NOT USED  %%%%%%%%%%%%%%
C
      SUBROUTINE ANGIN(S,PR,ICLTYP,T1,C1)
C
C      ANGULAR DISTRIBUTION FOR INELASTIC SCATTERING
C      (SEE ANGEL)
C
C        S      : C.M. ENERGY                          (INPUT)
C        PR     : RELATIVE MOMENTUM (INITIAL)          (INPUT)
C        ICLTYP : COLLISION TYPE                       (INPUT)
C        T1     : AZIMUTHAL ANGLE                      (OUTPUT)
C        C1     : COSINE(POLAR ANGLE)                  (OUTPUT)
C
C
        INCLUDE (VUUPRM)
C
      GOTO (
     >   1,  1,  1,  4,  4,
     >   4,  1,  1,  4,  4,
     >   1,  1,  4,  4,  4,
     >   4,  4,  4,  4,  4,
     >   4
     & ) , ICLTYP
C
      WRITE(6,*)' ICLTYP INVALID IN  ANGIN   ',ICLTYP
      CALL ERREX( 'FATAL ERROR')
C
    1 C1=1.-2.*RANF(0)
      T1=PI2*RANF(0)
      RETURN
    4 CALL ERREX( ' ANGIN -- 4')
      END
C
C***********************    C M A S S    ******************************
C
      SUBROUTINE CMASS(S,PR,EM1,EM2,ITYP1,ITYP2,ICLTYP,IC)
C
        INCLUDE (VUUPRM)
C
C---------------  UPDATES:  ------------------------------------------
C
C  8/24/83  SELECT MASSES USING 1ST AND 2ND CHOICE METHOD, TAG
C           2ND CHOICE DELTA WITH TYPE=100+DELTA-TYPE IF SWITCH
C           NUMBER 2 IS SET. LOGICAL INDICATOR SECOND ADDED IN SELECTION
C           SEQUENCE.
C  MEANING OF SW2 CHANGED -> NO DELTA DECAY
C
C  9/2/83   SW5 ADDED, FORCES MASS(DELTA)=1232
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      LOGICAL SECOND
C<>VERNUM<>
      DATA VERNUM/'VER   4.006 -- 16-JAN-85,20:19:00'/
      DATA PRONAM/'MODULE: C M A S S '/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C      FIND TYPES AND MASSES OF THE PARTICLES IN A BINARY
C      INELASTIC EXIT CHANNEL.
C
C        S      : C.M. ENERGY                           (INPUT)
C        PR     : REL. MOMENTUM (INITIAL)               (INPUT)
C        EM1    : MASS PART. 1                          (INPUT/OUTPUT)
C        EM2    : MASS PART. 2                          (INPUT/OUTPUT)
C        ITYP1  : TYPE PART. 1                          (INPUT/OUTPUT)
C        ITYP2  : TYPE PART. 2                          (INPUT/OUTPUT)
C        ICLTYP : COLLISION TYPE                        (INPUT)
C        IC     : EXIT CHANNEL SELECTOR                 (INPUT)
C
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
C ...  MAKE A RANDOM CHOICE FOR THE NUCLEON THAT WILL BECOME A DELTA
C
      IF (ICLTYP .LE. 3) THEN
        IF (RANF(0.) .GT. .5) THEN
          TMP=EM1
          EM1=EM2
          EM2=TMP
          I=ITYP1
          ITYP1=ITYP2
          ITYP2=I
        END IF
C
C ...  OTHERWISE FLIP A DELTA INTO THE  1  POSITION
C
      ELSE IF (ITYP1 .LE. 2) THEN
        TMP=EM1
        EM1=EM2
        EM2=TMP
        I=ITYP1
        ITYP1=ITYP2
        ITYP2=I
      END IF
C
      GOTO (
     >   1,  2,  3,  4,  5,
     >   4,  7,  8,  4,  4,
     >  11, 12,  4,  4,  4,
     >  16,  4,  4,  4,  4,
     >   4
     & ) , ICLTYP
C
      WRITE(6,*)' ICLTYP INVALID IN  CMASS   ',ICLTYP
      CALL ERREX( 'FATAL ERROR')
C
C ---  N+N -> N+DELTA
C
    1 ASSIGN 100 TO LABEL
      GOTO 1000
  100 IF (IC .EQ. 1) THEN
        EM1=EM
        ITYP1=4
      ELSE IF (IC .EQ. 2) THEN
        EM1=EM
        ITYP1=3
        ITYP2=2
      ELSE
        WRITE(6,*)'IC-ERROR CMASS-100'
      END IF
      RETURN
C
C ---  P+N -> P+DELTA // N+DELTA
C
    2 ASSIGN 200 TO LABEL
      GOTO 1000
  200 IF (IC .EQ. 1) THEN
        EM1=EM
        ITYP1=5
        ITYP2=1
      ELSE IF (IC .EQ. 2) THEN
        EM1=EM
        ITYP1=4
        ITYP2=2
      ELSE
        WRITE(6,*)'IC-ERROR CMASS-200'
      END IF
      RETURN
C
C ---  P+P -> P+DELTA
C
    3 ASSIGN 300 TO LABEL
      GOTO 1000
  300 IF (IC .EQ. 1) THEN
        EM1=EM
        ITYP1=5
      ELSE IF (IC .EQ. 2) THEN
        EM1=EM
        ITYP1=6
        ITYP2=1
      ELSE
        WRITE(6,*)'IC-ERROR CMASS-300'
      END IF
      RETURN
C
C --- INVALID COMBINATION
C
    4 CALL ERREX( ' CMASS -- 4')
C
C ---  DELTA ABSORPTION
C ...  D-,P
C
    5 EM1=EMNUC
      ITYP1=1
      ITYP2=1
      RETURN
C
C ...  D0,N
C
    7 EM1=EMNUC
      ITYP1=1
      RETURN
C
C ...  D0,P
C
    8 EM1=EMNUC
      IF (RANF(0.) .GT. .5) THEN
        ITYP1=1
      ELSE
        ITYP1=2
        ITYP2=1
      END IF
      RETURN
C
C ...  D+,N
   11 EM1=EMNUC
      IF (RANF(0.) .GT. .5) THEN
        ITYP1=2
      ELSE
        ITYP1=1
        ITYP2=2
      END IF
      RETURN
C
C ...  D+,P
C
   12 EM1=EMNUC
      ITYP1=2
      RETURN
C
C ...  D++,N
C
   16 EM1=EMNUC
      ITYP1=2
      ITYP2=2
      RETURN
C
C ...  MASS SELECTION SEQUENCE
C
 1000 EMMAX = S - EMNUC - EKINMI
C
      IF (SWITCH(5)) THEN
        EM=AMIN1(EMMAX,EMDELT)
        SECOND=.FALSE.
        GOTO LABEL
      END IF
C
      XMAX=(ATAN(2.*(EMMAX-EMDELT)/WIDDLT)-CONST1)/
     # CONST2
C
      XMAX=AMIN1(1.,XMAX)
      IF (XMAX .LE. 0.) THEN
        IC=0
        RETURN
      END IF
      SECOND=.FALSE.
      X=RANF(0.)
C        FIRST CHANCE DELTA
      IF (X .GT. XMAX) THEN
C  MISSED,  MAKE A SECOND-CHANCE DELTA
        SECOND=.TRUE.
        X=XMAX*RANF(0.)
      END IF
C
      T=TAN(X*CONST2)
      EM=EMDELT+.5*WIDDLT*(CONST+T)/(1.-CONST*T)
      GOTO LABEL
C
      END
C*********************   C  U P D A T E   **************************
C
      SUBROUTINE CUPDAT(I1,I2)
C
C
C
C NOTE: CUPDAT CAN UPDATE THE COLLISIONS FOR JUST ONE PARTICLE
C       IF ONE USES I2<-1. IN THAT CASE I1 MUST (!) REPRESENT A
C       LEGAL ENTRY (ITYP(I1)>0).
C
C
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
      LOGICAL HIT
C  ...  UPDATE THE COLLISION ARRAY
C
C       (A) REMOVE COLLSIONS INCLUDING I1,I2
C
      DO 100 JJJ=1,NUMENT
        IF ((I1C(JJJ) .EQ. I1) .OR. (I2C(JJJ) .EQ. I1) .OR.
     *      (I1C(JJJ) .EQ. I2) .OR. (I2C(JJJ) .EQ. I2))
     *     I1C(JJJ)=-1
  100 CONTINUE
C
C       (B) REVISE COLLISIONS WITH NEW I1,I2 PARTICLES
C
      IF (ITYP(I1) .GT. 0) THEN
        I=I1
        J=I2
      ELSE
        I=I2
        J=I1
      END IF
C
      DO 101 JJJE=1,NT+NPION(KSIMUL)
        JJJ=JJJE+IOFSET
        IF (HIT(I,JJJ,SRT,SIG,SIGEL,THIT)) THEN
C
C  ...  FIND A SLOT
C
          DO 102 L=1,NUMENT
            IF (I1C(L) .LE. 0) THEN
              JENTRY=L
              GOTO 1000
            END IF
 102      CONTINUE
          NUMENT=NUMENT+1
          IF (NUMENT .GT. NUMEMA)
     *      CALL ERREX( 'CUPDAT COLL. STORAGE')
          JENTRY=NUMENT
 1000     I1C(JENTRY)=I
          I2C(JENTRY)=JJJ
          CTIME(JENTRY)=THIT
          CS(JENTRY)=SRT
          CSIG(JENTRY)=SIG
          CSIGEL(JENTRY)=SIGEL
        END IF
        IF (J .GT. 0) THEN
          IF ((JJJ .NE. I) .AND. (ITYP(J) .GT. 0) .AND.
     *        HIT(J,JJJ,SRT,SIG,SIGEL,THIT)) THEN
C
C  ...  FIND A SLOT
C
            DO 103 L=1,NUMENT
              IF (I1C(L) .LE. 0) THEN
                JENTRY=L
                GOTO 1100
              END IF
 103        CONTINUE
            NUMENT=NUMENT+1
            IF (NUMENT .GT. NUMEMA)
     *         CALL ERREX( 'CUPDAT COLL. STORAGE')
            JENTRY=NUMENT
1100        I1C(JENTRY)=J
            I2C(JENTRY)=JJJ
            CTIME(JENTRY)=THIT
            CS(JENTRY)=SRT
            CSIG(JENTRY)=SIG
            CSIGEL(JENTRY)=SIGEL
          END IF
        END IF
 101  CONTINUE
      RETURN
      END
C**********************   D  U P D A T E   ****************************
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C 10/21/83 HK : DUPDAT IS CALLED FOR EVERY BINARY COLLSION,
C               IT NOW ACCEPTS ANY PARTICLE, NUCLEONS AND PIONS
C               ARE GIVEN A LIFETIME OF 1.E35,  TYPE-2  DELTAS
C               WILL HAVE A LIFETIME OF 1.E34.
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE DUPDAT(I1)
C
C  ...  UPDATE THE DECAY TIMES FOR A NEWLY CREATED PARTICLE
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
C  ...  DON'T TOUCH TYPE-2 DELTAS IF THIS IS NOT A NEW PARTICLE.
C
      IF ((ITYP(I1) .GT. 0) .AND. (DECAY(I1) .EQ. 1.E34)) RETURN
C
      ITYP(I1)=IABS(ITYP(I1))
C REMOVE THE  NEW PARTICLE  TAG
C
      IF (ITYP(I1) .GT. 100) THEN
C   TYPE 2  NON-DECAYING DELTA
        ITYP(I1)=ITYP(I1)-100
        DECAY(I1)=1.E34
C
      ELSE IF ((ITYP(I1) .LE. 2) .OR. (ITYP(I1) .GT. 6)) THEN
C NUCLEONS
        DECAY(I1)=1.E35
C + PIONS
      ELSE
        DECAY(I1)=TIME+DECTIM(EM(I1),PX(I1),PY(I1),PZ(I1))
      END IF
C
      RETURN
      END
C******************  D E L T A D E C  ******************************
C****************   D M O M , D W I D T H   ************************
C
C  ...  THIS FUNCTION SETS UP THE CALL TO DWIDTH, IT DETERMINES
C       THE MOMENTUM OF THE DECAY PION IN THE DELTA REST FRAME
C       FROM THE DELTA MASS.
C
C       EM     : MASS OF THE DECAYING DELTA (INPUT)
C       DMOM   : MOMENTUM OF THE DECAY PION IN THE DELTA RET FRAME
C                (OUTPUT)
C
      FUNCTION DMOM(EM)
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
        INCLUDE (VUUPRM)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 25-JAN-85,11:41:00'/
      DATA PRONAM/'D M O M -- DELTA DECAY MOMENTUM'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (TTVERS) THEN
        WRITE(6,*)'---------------- VERSION NUMBER REPORT -------------'
        WRITE(6,*)PRONAM
        WRITE(6,*)VERNUM
        WRITE(6,*)'----------------------------------------------------'
        TTVERS=.FALSE.
      END IF
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      EM2=EM*EM
      PSQ=((.5*(EM2-S12M22))**2-P12M22)/EM2
      IF (PSQ .LT. 0) THEN
        WRITE(6,*)'PSQ<0, DMOM, ',EM,PSQ
        DMOM=0.
      ELSE
        DMOM=SQRT(PSQ)
      END IF
      RETURN
      END
C---------------------------------------------------------------------
C
C  ...  THIS FUNCTION COMPUTES THE MASS-DEPENDENT DELTA DECAY WIDTH
C       REF: RANDRUP, NP A314 (1979) 429.
C            RITTENBERG, REV.MOD.PHYS. 43 (1971) S1.
C
C  ...  PPI   : PION MOMENTUM IN THE DELTA REST FRAME, IN GEV/C (INPUT)
C       DWIDTH: DELTA DECAY WIDTH
C
      REAL FUNCTION DWIDTH(PPI)
C
C-------------------   UPDATES:  -------------------------------------
C
C  08/24/83  IF SWITCH NUMBER 1 IS SET, USE A CONSTANT DECAY WIDTH
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
        INCLUDE (VUUPRM)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 25-JAN-85,11:42:30'/
      DATA PRONAM/'D W I D T H  -- DELTA DECAY WIDTH'/
      DATA TTVERS/.TRUE./
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
      IF (SWITCH(1)) THEN
       WIDTH=WIDDLT
      ELSE
       WIDTH=WIDDLT*(PPI**3/(1.+(PPI/PSCAL1)**2+(PPI/PSCAL2)**4))
     @ /(P0REF**3/(1.+(P0REF/PSCAL1)**2+(P0REF/PSCAL2)**4))
      END IF
C
C  ...  CONVERT FROM GEV TO (FM/C)**-1
C
      DWIDTH=5.0674*WIDTH
C
      RETURN
      END
C
C***********************   D E C T I M E   ************************
C
C  ...  COMPUTE A RANDOM CHOICE FOR THE TIME AT WHICH A DELTA
C       WILL DECAY; TRANSFORM TO THE COMPUTATIONAL FRAME.
C
C  ...  EMD      :  DELTA MASS (INPUT)
C       PXD      :
C       PYD      : ' DELTA MOMENTUM COMPONENTS (INPUT)
C       PXD      :
C       DECTIM   :  TIME UNTIL DECAY IN THE COMP. FRAME (OUTPUT)
C
C
      REAL FUNCTION DECTIM(EMD,PXD,PYD,PZD)
C
        INCLUDE (VUUPRM)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 25-JAN-85,11:45:00'/
      DATA PRONAM/'D E C T I M E  -- DELTA DECAY TIME'/
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
C  ...  GAMMA FOR THE DELTA RESTFRAME <-> COMP. FRAME TRANSFORMATION
      ED=SQRT(EMD*EMD+PXD*PXD+PYD*PYD+PZD*PZD)
      GG=ED/EMD
C
C  ...  GET THE WIDTH
      PPI=DMOM(EMD)
      WID=DWIDTH(PPI)*DA(13)
C
C  ...  REST FRAME DECAY TIME
      IF (WID .GT. 0.) THEN
        TAU=-LOG(AMAX1(1.-RANF(0),1.E-35))/WID
      ELSE
        DECTIM=1.E35
        RETURN
      END IF
C
C  ...  APPLY TIME DILATION
      DECTIM=TAU*GG
C
      RETURN
      END
C
C****************   D E L D E C  ******************************
C
C  ...  ALLOW A DELTA TO DECAY, CHANGE THE DELTA TO AN APPROPRIATE
C       NUCLEON AND ADD THE PION TO THE LIST.
C
C  ...  IND   : ARRAY INDEX OF THE DECAYING DELTA (INPUT)
C               NO OUTPUT
C  ...  NOTE - THE COMMON BLOCK /PARTICLEBLK/ WILL BE ALTERED.
C
      SUBROUTINE DELDEC(IND)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
      LOGICAL THIRD
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 25-JAN-85,11:46:15'/
      DATA PRONAM/'D E L D E C  -- DELTA DECAY SEQUENCE'/
      DATA TTVERS/.TRUE./
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C  ...  THIS FUNCTION IS TRUE IF THE RANDOM SELECTION PICKS THE
C       1/3 BRANCH
      THIRD=(RANF(0) .GT. .666666666)
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
C  ...  BRANCH ON DELTA TYPE
      GOTO (103,104,105,106),ITYP(IND)-2
      WRITE(6,*)'ERROR, SELECT,  DELDEC , IND,ITYP:',IND,ITYP(IND)
      CALL ERREX('FATAL ERROR')
C
C DL-
  103 IPI=7
      IN =1
      IC=3
      GOTO 1000
C
C DL0
  104 IF (THIRD) THEN
        IPI=7
        IN= 2
        IC=3
      ELSE
        IPI=8
        IN =1
        IC=4
      END IF
      GOTO 1000
C
C DL+
  105 IF (THIRD) THEN
        IPI=9
        IN =1
        IC=3
      ELSE
        IPI=8
        IN =2
        IC=4
      END IF
      GOTO 1000
C
C DL++
  106 IPI=9
      IN =2
      IC=3
      GOTO 1000
C
 1000 ICL=(ITYP(IND)*(ITYP(IND)+1))/2
      ICOUNT(IC,ICL,KSIMUL)=ICOUNT(IC,ICL,KSIMUL)+1
C
      PREL=DMOM(EM(IND))
C
C  ...  PICK RANDOM DIRECTION
C
      CT=1.-2.*RANF(0)
      ST=SQRT(1.-CT*CT)
      PHI=RANF(0)*PI2
      PXREL=PREL*COS(PHI)*ST
      PYREL=PREL*SIN(PHI)*ST
      PZREL=PREL*CT
C
C  ...  PION:
      IF (NPION(KSIMUL) .GE. NPIONM)
     &                CALL ERREX( ' DELDEC -- NO MORE STORAGE')
      NPION(KSIMUL)=NPION(KSIMUL)+1
      L=IOFSET+NT+NPION(KSIMUL)
      X(L)=X(IND)
      Y(L)=Y(IND)
      Z(L)=Z(IND)
      ITYP(L)=IPI
      EM(L)=EMPION
      DECAY(L)=1.E35
      LASTCL(L)=IND
      NCLCNT(L)=1
      CALL TRANS1(PXREL,PYREL,PZREL,EM(L),
     &            PX(L),PY(L),PZ(L),0,
     &            PX(IND),PY(IND),PZ(IND),EM(IND),
     &            +1.0,.FALSE.)
      AX(L)=0.0
      AY(L)=0.0
      AZ(L)=0.0
      FX(L)=0.0
      FY(L)=0.0
      FZ(L)=0.0
C
C  ...  NUCLEON
      ITYP(IND)=IN
      CALL TRANS1(-PXREL,-PYREL,-PZREL,EMNUC,
     &            PX(IND),PY(IND),PZ(IND),0,
     &            PX(IND),PY(IND),PZ(IND),EM(IND),
     &            +1.0,.FALSE.)
      EM(IND)=EMNUC
      DECAY(IND)=1.E35
      LASTCL(IND)=0
C
      RETURN
      END
C
C************************  T R A N S 1  ***********************
C
C  ...  LORENTZ TRANSFORMATION:
C       ARGS 1-4 DESCRIBE A 4-VECTOR IN THE REST FRAME OF THE PARTICLE
C       DEFINED BY ARGS 9-12. THE SUBROUTINE COMPUTES THE 4-VECTOR IN
C       THE COMP. FRAME AND STORES IT IN ARGS 5-8. THIS TRANSFORMATION
C       CAN BE REVERSED BY SETTING SIGN=-1. THE ENERGY COMPONENT IS
C       COMPUTED ONLY IF ECOMP=.TRUE.
C       NOTE - ALL 3 4-VECTORS CAN BE IDENTICAL
C
      SUBROUTINE TRANS1(PXP,PYP,PZP,EMP,
     &                  PX,PY,PZ,E,
     &                  P1,P2,P3,XM,
     &                  SIGN,ECOMP)
      LOGICAL ECOMP
C
C  ...  PXP,PYP,PZP  :  3-MOMENTUM IN THE PRIMED FRAME (INPUT)
C       EMP          :  MASS (!) IN THE PRIMED FRAME (INPUT)
C       PX,PY,PZ     :  3-MOMENTUM IN THE UNPRIMED FRAME (OUTPUT)
C       E            :  ENERGY (!) IN THE UNPRIMED FRAME
C                       (CONDITIONAL OUTPUT)
C       P1,P2,P3     :  3-MOMENTUM OF THE PARTICLE DEFINING THE REST
C                       FRAME (INPUT)
C       XM           :  MASS (!) OF THE REST-FRAME PARTICLE (INPUT)
C       SIGN         :  +1 : PRIMED FRAME IS THE REST FRAME, TRANSFORM
C                            TO THE COMP. FRAME
C                       -1 : PRIMED FRAME IS THE COMP. FRAME, TRANSFORM
C                            TO THE REST FRAME
C                       (INPUT)
C       ECOMP        :  LOGICAL SWITCH FOR COMP. OF E (INPUT)
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      CHARACTER*33 VERNUM
      CHARACTER*46 PRONAM
      LOGICAL*1 TTVERS
C<>VERNUM<>
      DATA VERNUM/'VER   7.006 -- 25-JAN-85,11:49:15'/
      DATA PRONAM/'T R A N S 1  -- LORENTZ TRANSFORMATION'/
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
      EP=SQRT(EMP*EMP+PXP*PXP+PYP*PYP+PZP*PZP)
      ET=SQRT(XM*XM+P1*P1+P2*P2+P3*P3)
      B1=P1/ET
      B2=P2/ET
      B3=P3/ET
      G=ET/XM
      BP=B1*PXP+B2*PYP+B3*PZP
      BPG=G*BP/(G+1.)
      BB=G*(BPG+EP*SIGN)
C
      PX=PXP+B1*BB
      PY=PYP+B2*BB
      PZ=PZP+B3*BB
C
      IF (ECOMP) E=G*(EP+SIGN*BP)
C
      RETURN
      END
C********************* V U U V L A S *********************************
C********************* C O M P ***************************************
C********************* COMP ******************************************
      FUNCTION COMP(R)
C
      INTEGER*4 IEQ
      REAL*4 P(12)
      COMMON /EDATA/P,IEQ
      IF (R .LE. 0.) THEN
        COMP=0.0
        RETURN
      END IF
      COMP=P(1)*R+(P(2)*R**(P(3)+1.))/(P(3)+1.)
      RETURN
      END
C
C****************** DVDRHO ******************************************
C
      FUNCTION DVDRHO(R)
C
C  ...  FIRST DERIVATIVE OF U W.R.T. RHO.
C
      INTEGER*4 IEQ
      REAL*4 P(12)
      COMMON /EDATA/P,IEQ
      IF (R .LE. 0.) THEN
        DVDRHO=0.0
        RETURN
      END IF
      IF (P(3) .EQ. 1.) THEN
        DVDRHO=P(1)+R*P(2)
      ELSE
        DVDRHO=P(1)+P(2)*R**P(3)
      END IF
      RETURN
      END
C
C****************************** INPUT ********************************
      SUBROUTINE SETEQS(COMINI)
C
        INCLUDE (VUUPRM)
C
      REAL*4 P(12)
      REAL*4 REQ,BEQ,COM,PMASS,RSH,COND,DEL
      INTEGER*4 IEQ,ISEL,IU
      COMMON /EDATA/P,IEQ
      COMMON /IOPRMS/REQ,BEQ,COM,PMASS,RSH,COND,DEL,ISEL
C     ... READ FIRST RECORD : IEQ,REQ,BEQ
      IEQ=IEQINI
      REQ=REQINI
      BEQ=BEQINI
      IF ((IEQ .LT. 1) .OR. (IEQ .GT. 3)) IEQ=1
C     ... SECOND RECORD : COM,PMASS,ISEL
C         ( ISEL MUST BE SPECIFIED BUT IS USED ONLY FOR IEQ=3 )
      COM=COMINI

      PMASS=EMNUC*1000.
C  EVERYTHING HERE IS IN MEV
      ISEL=ISEINI
      IF ((ISEL .LT. 1) .OR. (ISEL .GT. 2)) ISEL=1
      CALL CMPUTE
      RETURN
      END
C********************** CMPUTE ****************************************
      SUBROUTINE CMPUTE
C
        INCLUDE (VUUPRM)
C
C     ...  CMPUTE  PROCESSES THE GIVEN PARAMETERS
      REAL*4 P(12),R43,R53,R23,REQ1,EEQ1,COM1
      REAL*4 REQ,BEQ,COM,PMASS,RSH,COND,DEL,EF,EDGE
      INTEGER*4 IEQ,ISEL
      LOGICAL SUSOFT
      COMMON /EDATA/P,IEQ
      COMMON /IOPRMS/REQ,BEQ,COM,PMASS,RSH,COND,DEL,ISEL
      COMMON /VALUES/EF,EDGE
      DATA SUSOFT/.FALSE./
C
C
C     ... STORE DOUBLE PRECISION VALUES
      P(6)=REQ
      P(7)=-BEQ
      P(8)=COM
      P(9)=PMASS
C     ... FERMI ENERGY AND EDGE PARAMETERS
      P(5)=1.173742583E+05/P(9)
      EDGE=P(5)*REQ**(2./3.)
      EF=0.6*EDGE
C
      IF(COM.GT.0.0)THEN
C        CALCULATE PARAMETERS FROM THE COMPRESSIBILITY
C        ALL VALUES SHOULD BE IN GEV
C CHTEST
        IF(EF.LT.1)WRITE(6,*) '*** ATTENTION EF IN GEV'
        EFGEV=EF/1000.
        EB=-1.E-3 *BEQ
        C9=COM/9000.
C CHTEST
        WRITE(6,*) 'EB,EFGEV,C9',EB,EFGEV,C9
        GAM=(C9+EFGEV*2./9.)/(EFGEV/3. - EB)
        BET=(EFGEV/3.-EB)*(GAM+1.)/(GAM-1.)
        ALP=2.*(EB-EFGEV-BET/(GAM+1))
        WRITE(6,*)' POTENTIAL U=A*RHO/RHO0+B*(RHO/RHO0)**GAMMA'
        WRITE(6,*)' A,B,GAMMA (IN GEV) =',ALP,BET,GAM
        P(1)=ALP/REQ
        P(2)=BET*GAM/REQ**GAM
        P(3)=GAM-1.
        P(8)=COM
      ELSE
C        USE PARAMETER SETS
      ICOM=INT(ABS(COM)+0.5)+1
C       SELECT EQUATION OF STATE
      GOTO (1000,2000,3000,4000,1400,1500,1600),ICOM
 1000 CONTINUE
      WRITE(6,1100)ICOM
 1100 FORMAT(1X,'COM= ',I2)
      CALL ERREX('K=290 UNIMPLEMENTED EQUATION OF STATE')
      GOTO 5000
 2000 CONTINUE
      P(1)=-.356/REQ
      P(2)=.3535/REQ**1.1666666
      P(3)=.1666666
      P(8)=200.
      COM=200.
      WRITE(6,*) ' SOFT EOS, COM=',COM
      GOTO 5000
 3000 CONTINUE
      P(1)=-.124/REQ
      P(2)=.141/REQ**2
      P(3)=1.
      P(8)=380.
      COM=380.
      WRITE(6,*) ' HARD EOS, COM=',COM
      GOTO 5000
 4000 CONTINUE
      P(1)=-.124/REQ
      P(2)=.141/REQ**2
      P(3)=1.
      P(8)=380.
      COM=0.
      SUSOFT=.TRUE.
      WRITE(6,*) ' SUPERSOFT EOS, COM=',COM
      GOTO 5000
 1400 CONTINUE
      P(1)=.1552/REQ
      P(2)=-.1633/REQ**0.80
      P(3)=-0.20
      P(8)=121.
      COM=121.
      WRITE(6,*) ' VERY SOFT EOS, COM=',COM
      GOTO 5000
 1500 CONTINUE
      P(1)=-.1573/REQ
      P(2)=.1633/REQ**1.58
      P(3)=0.58
      P(8)=284.
      COM=284.
      WRITE(6,*) ' INTERMED. EOS, COM=',COM
      GOTO 5000
 1600 CONTINUE
      P(1)=-.1104/REQ
      P(2)=.1355/REQ**2.40
      P(3)=1.40
      P(8)=456.
      COM=456.
      WRITE(6,*) ' VERY HARD EOS, COM=',COM
      END IF
C
 5000 CONTINUE
      CALL EQSOUT
C
C... PRESET DERIVATIVE TABLE
C
      DELRHO=REQ*RHOMAX/FLOAT(NTABLE-1)
      VAL=0.0
      DVRTBL(0)=0.0
      DO 100 I=1,NTABLE
        DVRTBL(I)=DVDRHO(VAL)
        IF(SUSOFT.AND.(VAL.GT.REQ))DVRTBL(I)=0.
        VAL=VAL+DELRHO
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE EQSOUT
C     ... PRINT ON UNIT 6 A SUMMARY OF THE EQUATION OF STATE USED.
      REAL*4 P(12)
      REAL*4 REQ,BEQ,COM,PMASS,RSH,COND,DEL,EF,EDGE
      INTEGER*4 IEQ,ISEL
      COMMON /EDATA/P,IEQ
      COMMON /IOPRMS/REQ,BEQ,COM,PMASS,RSH,COND,DEL,ISEL
      COMMON /VALUES/EF,EDGE
      WRITE(6,8000)IEQ,ISEL,REQ,BEQ,COM,PMASS,RSH,COND,DEL
      WRITE(6,8010)EF,EDGE,P
      WRITE(6,8900)
      RETURN
 8000 FORMAT(1X,30(1H*),' EQUATION OF STATE ',30(1H*)//
     * ' VALUES SPECIFIED :'/
     * '   EQUATION TYPE ',I2,' (SELECT SWITCH FOR TYPE 3 IS ',I1,')'/
     * '   EQUILIBRIUM DENSIT : ',F8.4,' 1/FM**3'/
     * '   BINDING ENERGY AT EQUILIBRIUM : ',F8.2,' MEV'/
     * '   INCOMPRESSIBILITY : ',F9.2,' MEV'/
     * '   PARTICLE MASS : ',F9.2,' MEV'/
     * '   CONDENSATION ENERGY :'/
     * '      SHIFT DENSIT : ',F7.2/
     * '      STRENGTH      : ',F9.2,' MEV'/
     * '      WIDTH         : ',F8.4/1X,79(1H-))
 8010 FORMAT(' COMPUTED VALUES:'/
     * '   FERMI ENERGY AT T=0 IS ',F8.2,' MEV'/
     * '   FERMI EDGE AT T=0 IS AT ',F8.2,' MEV'/
     * ' THE PARAMETER FIELD CONTAINS :'/
     * 2(1X,1P,6E12.4/),1X,79(1H-))
 8900 FORMAT(1X,79(1H*))
      END
C******************* F I E L D  *************************************
C********************  D E N S  *************************************
C      MOVE TO END OF FILE
C***********************   ACCEL  *******************************
C
      SUBROUTINE ACCEL(I)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
C  ...  COMPUTE ACCELERATION FROM FORCES AND MOMENTA
C
      E2=EM(I)*EM(I)+PX(I)*PX(I)+PY(I)*PY(I)+PZ(I)*PZ(I)
      E=SQRT(E2)
      BX=PX(I)/E
      BY=PY(I)/E
      BZ=PZ(I)/E
      FXE=FX(I)/E
      FYE=FY(I)/E
      FZE=FZ(I)/E
      BFE=BX*FXE+BY*FYE+BZ*FZE
      AX(I)=FXE-BX*BFE
      AY(I)=FYE-BY*BFE
      AZ(I)=FZE-BZ*BFE
      RETURN
      END
C*********************   INITP        **********************************
C
      SUBROUTINE INITP(RR,RP)
C
      COMMON /PLIBLK/VEFFR(100),VEFFP(100),CONVR,CONVP
C
      WRITE(6,*)'EFFECTIVE VOLUME, R-SPACE'
      CALL EFFSET(RR,CONVR,VEFFR)
      WRITE(6,*)'EFFECTIVE VOLUME, P-SPACE'
      CALL EFFSET(RP,CONVP,VEFFP)
      WRITE(6,*)'------------------------------------------------------'
C
      RETURN
      END
C....................   EFFSET   ........................
C
      SUBROUTINE EFFSET(R,C,VEFF)
      DIMENSION VEFF(100)
C
C... COMPUTE THE EFFECTIVE VOLUME AS A FUNCTION OF THE CM SHIFT
C    BY THE WAY, THIS IS A TERRIBLE PIECE OF CODE
C
      CMM2=(3.*R/8.)**2
      D=CMM2/99.
C
      DZM=R/9999.
      VEFF(1)=1.0
      VEFF(100)=0.5
C
      I=2
      ZM=R-DZM
 1000 CM=SQRT(D*(I-1))
 1100 CMTEST=0.25*((R*R-ZM*ZM)**2)/(ZM*(R*R-ZM*ZM/3.)+2.*R*R*R/3.)
      IF (CMTEST .GE. CM) GOTO 1200
      ZM=ZM-DZM
      GOTO 1100
 1200 VEFF(I)=.5+.75*ZM/R-.25*ZM*ZM*ZM/R/R/R
      I=I+1
      IF (I .LT. 100) GOTO 1000
C
      C=1./D
C
      WRITE(6,1300)R,C
 1300 FORMAT(2F15.6)
      WRITE(6,1400)VEFF
 1400 FORMAT(1X,(5F10.5))
C
      RETURN
      END
C
C SAVE STORAGE: REPLACE THESE SUBROUTINES BY SUBROUTINES FROM VUUAUX
C               AND DELETE COMMON BLOCK NOINFB IN VUUDIM
C*****************  P A U L I  ************************************
C
      SUBROUTINE PAULI(I,NTAG)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C--
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
      NCUM=-1
      IO=0
      DO 100 L=1,3
        PCM(L)=0.0
        RCM(L)=0.0
 100  CONTINUE
      DO 101 KKK=1,NSIMUL
        DO 102 J=1+IO,NT+IO
          XYZSQ(J)= RTEST2-(X(J)-X(I))**2
     &            -(Y(J)-Y(I))**2-(Z(J)-Z(I))**2
 102    CONTINUE
        DO 103 JE=1,NT
          J=IO+JE
          IF (XYZSQ(J).GT. 0) THEN
           IF (ITYP(J) .LE. 2) THEN
             PXYZSQ= PTEST2-(PX(J)-PX(I))**2
     &            -(PY(J)-PY(I))**2-(PZ(J)-PZ(I))**2
             IF (PXYZSQ.GT. 0) THEN
                  NCUM=NCUM+1
                  PCM(1)=PCM(1)+PX(J)
                  PCM(2)=PCM(2)+PY(J)
                  PCM(3)=PCM(3)+PZ(J)
                  RCM(1)=RCM(1)+X(J)
                  RCM(2)=RCM(2)+Y(J)
                  RCM(3)=RCM(3)+Z(J)
              END IF
            END IF
          END IF
 103    CONTINUE
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
      RETURN
      END
C********************  D E N S  *************************************
      SUBROUTINE DENS(NOUTSI)
C
        INCLUDE (VUUPRM)
        INCLUDE (VUUDIM)
C
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
 10   NTPIMX=NT+NPIONM
      LOOPMX=NSIMUL*NTPIMX
      DO 80 I=1,LOOPMX
       DENSPP(I)=0.
       FX(I)=0.
       FY(I)=0.
       FZ(I)=0.
 80   CONTINUE
      DO 100 KKK=1,NSIMUL
        IO=(KKK-1)*NTPIMX
        DO 101 IE=1,NT
          I=IE+IO
C
C  ...  ZERO THE COUNTING SPHERE AROUND THE PARTICLE
C
          XII=X(I)
          YII=Y(I)
          ZII=Z(I)
C
C  ...  LOOP OVER ALL PARTICLES <> I
C
          JEE=I
          JO=IO
          DO 140 LLL=KKK,NSIMUL
            DO 141 J=JEE+1,JO+NT
C  ...  IF THE PARTICLE IS IN THE TEST-SPHERE:
C
              AX(J)=XII-X(J)
              AY(J)=YII-Y(J)
              AZ(J)=ZII-Z(J)
              XYZSQ(J)=FRAD2-AX(J)*AX(J)-AY(J)*AY(J)-AZ(J)*AZ(J)
 141        CONTINUE
            DO 142 J=JEE+1,JO+NT
                IF (XYZSQ(J).GT.0) THEN
C
C  ...  COUNT THE PARTICLE IN THE APPROPRIATE HALF-SPHERE
C
                 AX(J)=SIGN(1.,AX(J))
                 AY(J)=SIGN(1.,AY(J))
                 AZ(J)=SIGN(1.,AZ(J))
                 FX(I)=FX(I)+AX(J)
                 FY(I)=FY(I)+AY(J)
                 FZ(I)=FZ(I)+AZ(J)
                 DENSPP(I)=DENSPP(I)+1.
                 FX(J)=FX(J)-AX(J)
                 FY(J)=FY(J)-AY(J)
                 FZ(J)=FZ(J)-AZ(J)
                 DENSPP(J)=DENSPP(J)+1.
C               ELSE
C                AX(J)=0.
C                AY(J)=0.
C                AZ(J)=0.
                END IF
 142        CONTINUE
 143        CONTINUE
            JEE=LLL*NTPIMX
            JO=JEE
 140      CONTINUE
C
C  ...  COMPUTE THE DENSIT AND THE DERIVATIVE OF THE POTENTIAL
C       (NOTE THAT ANY PAIR OF HALF-SPHERES CONTAINS ALL PARTICLES)
C
          DENSPP(I)=DENSPP(I)/VOLUME
          III=INT(DENSPP(I)/DELRHO +0.5)
C        USE AX ONLY FOR INTERMEDIATE STORAGE
          AX(I)=DVRTBL(MIN0(NTABLE,III))
C
C  ...  COMPUTE THE FORCES AND ACCELERATIONS
C
          AX(I)=GRAD*AX(I)
 101    CONTINUE
 100  CONTINUE
      DO 170 KKK=1,NSIMUL
         IO=(KKK-1)*NTPIMX
         DO 171 I=IO+1,IO+NT
          FX(I)=AX(I)*FX(I)
          FY(I)=AX(I)*FY(I)
          FZ(I)=AX(I)*FZ(I)
 171     CONTINUE
         DO 172 I=IO+1,IO+NT
          XYZSQ(I)= EM(I)*EM(I)+PX(I)*PX(I)+PY(I)*PY(I)+PZ(I)*PZ(I)
          AZ(I)=PX(I)*FX(I)+PY(I)*FY(I)+PZ(I)*FZ(I)
          AZ(I)=AZ(I)/XYZSQ(I)
          XYZSQ(I)=SQRT(XYZSQ(I))
          AX(I)=(FX(I)-PX(I)*AZ(I))/XYZSQ(I)
          AY(I)=(FY(I)-PY(I)*AZ(I))/XYZSQ(I)
          AZ(I)=(FZ(I)-PZ(I)*AZ(I))/XYZSQ(I)
 172     CONTINUE
 170  CONTINUE
      RETURN
      END
C
C SAVE STORAGE: END OF REPLACE
C
C**************   T R A N S P O R T   ********************************
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
      DIMENSION E1(NTOTMA)
      EQUIVALENCE (XYZSQ,E1)
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
          E1(I1)=SQRT(EM(I1)**2+PZ(I1)**2+PX(I1)**2+PY(I1)**2)
          X(I1)=X(I1)+TADV*PX(I1)/E1(I1)
          Y(I1)=Y(I1)+TADV*PY(I1)/E1(I1)
          Z(I1)=Z(I1)+TADV*PZ(I1)/E1(I1)
 200    CONTINUE
 201    CONTINUE
      GOTO 4000
 3000   CONTINUE
         NTNPIO=IOFSET+NT+NPION(KSIMUL)
        DO 300 I1=1+IOFSET,NTNPIO
C
          E1(I1)=SQRT(EM(I1)**2+PZ(I1)**2+PX(I1)**2+PY(I1)**2)
C
          X(I1)=X(I1)+TADV*(PX(I1)/E1(I1)+.5*AX(I1)*TADV)
          Y(I1)=Y(I1)+TADV*(PY(I1)/E1(I1)+.5*AY(I1)*TADV)
          Z(I1)=Z(I1)+TADV*(PZ(I1)/E1(I1)+.5*AZ(I1)*TADV)
C         CHANGE AX TO FX
          PX(I1)=PX(I1)+TADV*FX(I1)
          PY(I1)=PY(I1)+TADV*FY(I1)
          PZ(I1)=PZ(I1)+TADV*FZ(I1)
 300    CONTINUE
C
C ...  UPDATE THE ACCELERATION (FORCES REMAIN UNCHGED)
C
C         CALL ACCEL(I1)
        DO 301 I=1+IOFSET,NTNPIO
          E1(I)= EM(I)*EM(I)+PX(I)*PX(I)+PY(I)*PY(I)+PZ(I)*PZ(I)
          AZ(I)=PX(I)*FX(I)+PY(I)*FY(I)+PZ(I)*FZ(I)
          AZ(I)=AZ(I)/E1(I)
          E1(I)=SQRT(E1(I))
          AX(I)=(FX(I)-PX(I)*AZ(I))/E1(I)
          AY(I)=(FY(I)-PY(I)*AZ(I))/E1(I)
          AZ(I)=(FZ(I)-PZ(I)*AZ(I))/E1(I)
 301    CONTINUE
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
C
C SAVE STORAGE: END OF REPLACE
C
C
C************************ RANF ***********************************
C
      FUNCTION RANF(X)
      RANF=DRAND48()
      RETURN
C           INITIALIZATION OF THE SEED WITH AN ODD INTEGER
      ENTRY SEED(I)
      ISEED=I+1-MOD(I,2)
      CALL SRAND48(ISEED)
      WRITE(6,*)'FIXED SEED = ',ISEED
      RETURN
C
      END
