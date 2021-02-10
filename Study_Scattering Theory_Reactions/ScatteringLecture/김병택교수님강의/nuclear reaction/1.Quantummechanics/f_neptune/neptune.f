      PROGRAM NEPTUNE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC   NEW VERSION OF NEPTUNE                          CCCCCCCCCCCCCCC
CCCCCC         (SPRING OF 2012 BY B. T. KIM)             CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CONST /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/BSX/   URSAVE(500),ENEPT,VNEPT,VSOR,DFNR,DFNSO,RZR,RZSO,
     1              RZC,XMES2,Q,FTR,NODER,LBTR,JBTRTW,NXRAWF,KTRL4,
     2              KOPOT,LBTRY,JBTWY
      COMMON/UNCPSA/KTRL3,ITBEMX,ACURCY,AMUPMU,
     1              KTRL2,KTRL8,KEX2,KEX4,KEX40,KEX41,KEX42,KEX43,
     2              TMAS,PMAS,ZZT,ZZP,RMAS,ZZ,XMES1,PERCNT,VSX,
     3              ISTW,NXRA,NXRM,NXRMP1,NXRMP2,NXRMP3,NXRMP4,NXRMP5,
     4              NODE,KGES,EGESRD,EGES,EGEST,DELGES,FKAPPA,FKAPIN,
     5              URRMIN,URRMEX,RGDLIN(3),RGDLEX(3),
     6              XMEM(514),VCENTR(514),VSPIN(514),VCOULM(514),
     7              PFORM1(2),PFORM2(2),PFORM3(2),GESMEM(20),NODMEM(10)
      COMMON/CNTRL/ KTRL(24),KEXCOM(24),KTLOUT(24),EXTCOM(50)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC      IF KTRL2 = 0; B.E. IS SEARCHED, = 1; VSX IS SEARCHED
CCCCCC      IF KTRL3 = 1; R=RZ*(A1**(1/3)+A2**(1/3))
CCCCCC      IF KTRL4 = 0; WAVE FUNCTION ITSELF IS STORED IN URSAVE.
CCCCCC               = 1; WAVE FT TIMES POTENTIAL
CCCCCC               = 2; WAVE FT TIMES R**(-L) 
CCCCCC               = 3; WAVE FT TIMES POTENTIAL TIMES R**(-L)
CCCCCC
CCCCCC      KEX2=NXCPL2 IF NON-ZERO
CCCCCC      KEX4+NXCPL2=NXRA IF KEX4 IS NON-ZERO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      open(unit=7, file='nep.DAT',status='old')
      open(unit=8, file='nep.OUT',status='unknown')
CCCCCC
CCCCCC     SECTION 1. DEFINE CONSTANTS 
CCCCCC
      FACLOG(1)=0.0D0
      FACLOG(2)=0.0D0
      FN=1.D0
      DO 5 N=3,500
      FN=FN+1.D0
    5 FACLOG(N)=FACLOG(N-1)+DLOG(FN)
      PI  =4.0*DATAN(1.D0) 
      HBAR=197.327053D0
      AMAS=931.49432D0
      WNUNIT=DSQRT(2.D0*AMAS)/HBAR
      FINE=1.0D0/137.0359896D0
CCCCCC
CCCCCC     SECTION 2. READ INITIAL DATA      CCCCCCCCC
CCCCCC
   10 FORMAT(14I5)
   11 FORMAT(10F7.3)
      READ(7,10) (KTRL(N),N=1,14)
	READ(7,10) (KEXCOM(N),N=1,14)
      READ(7,10) (KTLOUT(N),N=1,14)
	READ(7,11) (EXTCOM(N),N=1,10)
	READ(7,10) ITBEMX,JBTRTW,LBTR,NODER
      READ(7,11) TMAS,PMAS,ZZT,ZZP,ACURCY,PERCNT,EGES
      READ(7,11) VSX,VSOR,DFNR,DFNSO,RZR,RZSO,RZC
CCCCCC
CCCCCCCCCCCCC    END OF INPUT DATA             CCCCCCCCCCCCCCCCCCCC           
CCCCCC
CCCCCC
CCCCCC     SECTION 3. PRINT OUT OF INITIAL INPUT 
CCCCCC
CCCCCC
CCCCCC     SECTION 4. CALCULATES PH STATES IN THE TARGET SYSTEM     
CCCCCC 
CCCCCC
CCCCCC   EXCHANGE OF PARTICEL AND HOLE FOR NUCLEON-NUCLEUS CASE
CCCCCC
      ISTW=1
      KOPOT=KTLOUT(2)
      KTRL2=KTRL(2)
      KTRL4=KTRL(4)
      KTRL8=KTRL(8)
      KEX2=KEXCOM(2)
      KEX4=KEXCOM(4)
      AMUPMU=0.0
      XMES2=EXTCOM(1)
      CALL BSAXON
	STOP
	END
CCCCCC
CCCCCC
C======================================================================CCCN52470
      SUBROUTINE BSAXON                                                 CCN52480
C======================================================================CCCN52490
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CONST /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/BSX/   URSAVE(500),ENEPT,VNEPT,VSOR,DFNR,DFNSO,RZR,RZSO,   CCN52500
     1              RZC,XMES2,Q,FTR,NODER,LBTR,JBTRTW,NXRAWF,KTRL4,     CCN52510
     2              KOPOT,LBTRY,JBTWY                                   CCN52520
      COMMON/UNCPSA/KTRL3,ITBEMX,ACURCY,AMUPMU,                         CCN52530
     1              KTRL2,KTRL8,KEX2,KEX4,KEX40,KEX4L,KEX42,KEX43,      CCN52540
     2              TMAS,PMAS,ZZT,ZZP,RMAS,ZZ,XMES1,PERCNT,VSX,         CCN52550
     3              ISTW,NXRA,NXRM,NXRMP1,NXRMP2,NXRMP3,NXRMP4,NXRMP5,  CCN52560
     4              NODE,KGES,EGESRD,EGES,EGEST,DELGES,FKAPPA,FKAPIN,   CCN52570
     5              URRMIN,URRMEX,RGDLIN(3),RGDLEX(3),                  CCN52580
     6              XMEM(514),VCENTR(514),VSPIN(514),VCOULM(514),       CCN52590
     7              PFORM1(2),PFORM2(2),PFORM3(2),GESMEM(20),NODMEM(10) CCN52600
      COMMON/CNTRL/ KTRL(24),KEXCOM(24),KTLOUT(24),EXTCOM(50)
      DATA          IQI,IQJ/4HB.E.,4HVSX  /                             CCN52610
      DATA          KK,KL,KM,KN / 1HX,4HVCEN,4HVSPI,4HCOUL /            CCN52620
      PAI=PI                                                            CCN52630
      ITENOD=0                                                          CCN52640
      PERCNT=0.2                                                        CCN52650
      UNITOK=0.2187066  *(1.  -AMUPMU)+.21953760  *AMUPMU               CCN52660
      EGEST=EGES                                                        CCN52670
CCCCCC  ******                INITIAL OUTPUT                      ******CCN52680
      KOUT=KTLOUT(1)
      IF (KOPOT.NE.0) KOUT=1
      IF (KOUT.NE.0) THEN
      WRITE(8,61)                                                       CCN52690
   61 FORMAT(1H ,55X,17H BSAXON IS CALLED)                              CCN52700
      WRITE(8,63) KTRL2,KTRL3,KTRL4,KTRL8,KEX2,KEX4                     CCN52710
   63 FORMAT(//10X,12HKTRL2,3,4,8=,4I3,10H   KEX2,4=,2I5)               CCN52720
      WRITE(8,65) PMAS,ZZP,TMAS,ZZT                                     CCN52730
   65 FORMAT(/10X,34HMOTION BETWEEN NUCLEI WITH (A,Z)=(,2F7.2,7H) AND (,CCN52740
     1        2F7.2,1H))                                                CCN52750
      WRITE(8,67) EGES,LBTR,ISTW,JBTRTW,NODER                           CCN52760
   67 FORMAT(10X,10HWITH EGES=,F7.3,6H  LTR=,I1,5H  2S=,I1,5H  2J=,I2,  CCN52770
     1      8H  NODER=,I1)                                              CCN52780
      WRITE(8,69)                                                       CCN52790
   69 FORMAT(/17X,3HVSX,6X,4HVSOR,6X,4HDFNR,5X,5HDFNSO,7X,3HRZR,6X,     CCN52800
     1     4HRZSO,7X,3HRZC)                                             CCN52810
      WRITE(8,71) VSX,VSOR,DFNR,DFNSO,RZR,RZSO,RZC                      CCN52820
   71 FORMAT(10X,7F10.3)                                                CCN52830
      WRITE(8,73) AMUPMU,ACURCY                                         CCN52840
   73 FORMAT(/10X,7HAMUPMU=,F4.1,11H  ACCURACY=,F8.6)                   CCN52850
      ENDIF

      RMAS=(TMAS*PMAS)/(TMAS+PMAS)                                      CCN52860
      ZZ=ZZT*ZZP                                                        CCN52870
      NODE=NODER                                                        CCN52880
      FKTRL2=KTRL2                                                      CCN52890
      XMES1=XMES2*.125                                                  CCN52900
      XBFAC=TMAS**.3333333333                                           CCN52910
      IF(KTRL3.EQ.1) XBFAC=XBFAC+PMAS**.3333333333                      CCN52920
      XBAR=RZR*XBFAC                                                    CCN52930
      XRA=XBAR+10.0  *DFNR                                              CCN52940
      NXRAWF=XRA/XMES2                                                  CCN52950
      IF(KEX2.EQ.0) GO TO 113                                           CCN52960
      NXRAWF=KEX2                                                       CCN52970
  113 IF(KEX4.EQ.0) GO TO 117                                           CCN52980
      NXRA=NXRAWF+KEX4                                                  CCN52990
      GO TO 119                                                         CCN53000
  117 NXRA=NXRAWF+30                                                    CCN53010
119   NXRM=(XBAR+.5*DFNR*DFLOAT(NODE))/XMES2                            CCN53020
      FNXRM=NXRM                                                        CCN53030
      NXRMP1=NXRM+1                                                     CCN53040
      NXRMP2=NXRM+2                                                     CCN53050
      NXRMP3=NXRM+3                                                     CCN53060
      NXRMP4=NXRM+4                                                     CCN53070
      NXRMP5=NXRM+5                                                     CCN53080
      NXRA12=NXRA+12                                                    CCN53090
      IF (KOUT.NE.0) THEN
      WRITE(8,135) NXRAWF,NXRA,NXRM,XMES2,XBAR                          CCN53100
135   FORMAT( 10X,7HNXRAWF=,I5,3X,5HNXRA=,I5,3X,5HNXRM=,I5,             CCN53110
     1      8H  XMES2=,F7.4,7H  XBAR=,F7.4)                             CCN53120
      ENDIF
CCCCCC  ******    CONSTRUCTION OF THE BINDING POTENTIAL(S)        ******CCN53130
      NTTLMS=NXRA+12                                                    CCN53140
      XBARSO=XBFAC*RZSO                                                 CCN53150
      XBARC =XBFAC*RZC                                                  CCN53160
      VSPFC=2.0  *VSOR/DFNSO                                            CCN53170
      VCLFC2=1.4398650  *ZZ                                             CCN53180
      VCLFC1=VCLFC2*0.5  /XBARC                                         CCN53190
      DO 180 ND=1,4                                                     CCN53200
      FND=2.0  **(ND-1)                                                 CCN53210
      DX=XMES1*FND                                                      CCN53220
      IF(ND.NE.1) GO TO 163                                             CCN53230
      X=0.0                                                             CCN53240
      NXIMIN=1                                                          CCN53250
      NXIMAX=8                                                          CCN53260
      GO TO 165                                                         CCN53270
  163 NXIMIN=NXIMAX+1                                                   CCN53280
      NXIMAX=NXIMAX+4                                                   CCN53290
      IF(ND.EQ.4) NXIMAX=NTTLMS+2                                       CCN53300
      IF(NXIMAX.LE.1500) GO TO 164                                      CCN53310
      WRITE(8,725) NXIMAX                                               CCN53320
  725 FORMAT(2X,'NXIMAX.GT.514,=',I5)                                   CCN53330
      STOP                                                              CCN53340
  164 CONTINUE                                                          CCN53350
  165 DO 175 NX=NXIMIN,NXIMAX                                           CCN53360
      X=X+DX                                                            CCN53370
      XMEM(NX)=X                                                        CCN53380
      PFORM1(1)= DEXP((X-XBAR  )/DFNR )                                 CCN53390
      PFORM1(2)= DEXP((X-XBARSO)/DFNSO)                                 CCN53400
      DO 167 N=1,2                                                      CCN53410
      PFORM2(N)=1.0  /(1.0  +PFORM1(N))                                 CCN53420
      PFORM3(N)=PFORM1(N)*PFORM2(N)*PFORM2(N)                           CCN53430
  167 CONTINUE                                                          CCN53440
      VCENTR(NX)=-VSX*PFORM2(1)                                         CCN53450
      VSPIN(NX)=(VSPFC*PFORM3(2))/X                                     CCN53460
      IF(X.GT.XBARC) GO TO 172                                          CCN53470
      VCOULM(NX)=VCLFC1*(3.0  -((X/XBARC)**2))                          CCN53480
      GO TO 175                                                         CCN53490
  172 VCOULM(NX)=VCLFC2/X                                               CCN53500
  175 CONTINUE                                                          CCN53510
  180 CONTINUE                                                          CCN53520
      IF(KOPOT.EQ.0) GO TO 190                                          CCN53530
      NXWRIT=NTTLMS+2                                                   CCN53540
      WRITE(8,182) KK                                                   CCN53550
      WRITE(8,184)  (XMEM  (NX),NX=1,20)                                CCN53560
      WRITE(8,184)  (XMEM  (NX),NX=21,NXWRIT,KOPOT)                     CCN53570
      WRITE(8,182) KL                                                   CCN53580
      WRITE(8,184)  (VCENTR(NX),NX=1,NXWRIT,KOPOT)                      CCN53590
      WRITE(8,182) KM                                                   CCN53600
      WRITE(8,184)  (VSPIN (NX),NX=1,NXWRIT,KOPOT)                      CCN53610
      WRITE(8,182) KN                                                   CCN53620
      WRITE(8,184)  (VCOULM(NX),NX=1,NXWRIT,KOPOT)                      CCN53630
  182 FORMAT(1X,A4)                                                     CCN53640
  184 FORMAT(10E12.4)                                                   CCN53650
CCCCCC  ******    BEGINNING OF ITERATIVE DETERMINATION OF BINDING ******CCN53660
CCCCCC  ******    ENERGY OR DEPTH OF THE SAXON POTENTIAL          ******CCN53670
  190 KEX41=0                                                           CCN53680
      KEX42=0                                                           CCN53690
  205 ITEBE=0                                                           CCN53700
  210 ITEBE=ITEBE+1                                                     CCN53710
      GESMEM(ITEBE)=EGEST*(1.0  -FKTRL2)+VSX*FKTRL2                     CCN53720
CCCCCC  ******    K40=1 AND 2 ARE INTEGRATION OF THE INTERNAL     ******CCN53730
CCCCCC  ******    SOLUTION OUTWARDS AND THE EXTERNAL SOLUTION     ******CCN53740
CCCCCC  ******    INWARDS RESPECTIVELY TO BE MATCHED SMOOTHLY     ******CCN53750
CCCCCC  ******    AT THE NUCLEAR RADIUS.                          ******CCN53760
      DO 400 K40=1,2                                                    CCN53770
      KEX40=K40-1                                                       CCN53780
      IF(KTRL2.EQ.1) GO TO 213                                          CCN53790
      EGESDL=0.005  *EGEST                                              CCN53800
      EGES=EGEST-2.0  *EGESDL                                           CCN53810
      GO TO 214                                                         CCN53820
  213 EGES=EGEST                                                        CCN53830
  214 DO 390 KGES=1,3                                                   CCN53840
      FKGES=KGES                                                        CCN53850
      IF(KTRL2.EQ.1) GO TO 232                                          CCN53860
      EGES=EGES+EGESDL                                                  CCN53870
      GO TO 235                                                         CCN53880
  232 VCOREC=1.0  +0.005  *(FKGES-2.0  )                                CCN53890
      VCORIN=1.0  /VCOREC                                               CCN53900
      IF(KGES.EQ.2) GO TO 235                                           CCN53910
      DO 234 NX=1,NXRA12                                                CCN53920
      VCENTR(NX)=VCENTR(NX)*VCOREC                                      CCN53930
  234 CONTINUE                                                          CCN53940
  235 EGESIN=VSX-EGES                                                   CCN53950
      FKAPPA=UNITOK   * DSQRT(RMAS*EGES)                                CCN53960
      FKAPIN=UNITOK   * DSQRT(RMAS*EGESIN)                              CCN53970
      CALL UNCPST                                                       CCN53980
      IF(KGES.EQ.2.OR.KTRL2.EQ.0) GO TO 315                             CCN53990
      DO 300 NX=1,NXRA12                                                CCN54000
      VCENTR(NX)=VCENTR(NX)*VCORIN                                      CCN54010
300   CONTINUE                                                          CCN54020
315   IF(KEX42) 317,390,316                                             CCN54030
  316 CORNOD=1.1  *(1.0  -FKTRL2)+0.9  *FKTRL2                          CCN54040
      GO TO 318                                                         CCN54050
  317 CORNOD=0.9  *(1.0  -FKTRL2)+1.1  *FKTRL2                          CCN54060
  318 IF(KTRL2.EQ.1) GO TO 323                                          CCN54070
      EGEST=EGEST*CORNOD                                                CCN54080
      GO TO 341                                                         CCN54090
  323 VSX=VSX*CORNOD                                                    CCN54100
      DO 335 NX=1,NXRA12                                                CCN54110
      VCENTR(NX)=VCENTR(NX)*CORNOD                                      CCN54120
  335 CONTINUE                                                          CCN54130
  341 ITENOD=ITENOD+1                                                   CCN54140
      NODMEM(ITENOD)=KEX42+NODE                                         CCN54150
      KEX42=0                                                           CCN54160
      IF(ITENOD.LE.10) GO TO 205                                        CCN54170
      WRITE(8,354) (NODMEM(I1),I1=1,ITENOD)                             CCN54180
  354 FORMAT(10H ITENOD=10,10I3)                                        CCN54190
CCCCCC  ******        CONVERGENCE HAS NOT BEEN ACHIEVED AND       ******CCN54200
CCCCCC  ******        THE CALCULATION IS GIVEN UP                 ******CCN54210
      STOP                                                              CCN54220
  390 CONTINUE                                                          CCN54230
  400 CONTINUE                                                          CCN54240
      GEST=EGEST*(1.0  -FKTRL2)+VSX*FKTRL2                              CCN54250
      IF(DABS(DELGES/GEST).GT. .01) DELGES=.7*DELGES                    CCN54260
      IF(DABS(DELGES/GEST).GT. .05) DELGES=.7*DELGES                    CCN54270
      IF( DABS(DELGES/GEST).LE.PERCNT) GO TO 417                        CCN54280
      DELGES=   (GEST*DELGES/ DABS(DELGES))*PERCNT                      CCN54290
  417 GEST=GEST-DELGES                                                  CCN54300
      IF(KTRL2.EQ.1) GO TO 427                                          CCN54310
      EGEST=GEST                                                        CCN54320
      GO TO 440                                                         CCN54330
  427 VCOREC=GEST/(GEST+DELGES)                                         CCN54340
      VSX=VSX*VCOREC                                                    CCN54350
      DO 430 NX=1,NXRA12                                                CCN54360
      VCENTR(NX)=VCENTR(NX)*VCOREC                                      CCN54370
  430 CONTINUE                                                          CCN54380
  440 IF( DABS(DELGES/GEST).LE.ACURCY) GO TO 505                        CCN54390
      IF(KOPOT.NE.0)                                                    CCN54400
     1WRITE (8,447)  ITEBE,VSX,EGEST,EGES,EGESIN,DELGES,GEST            CCN54410
  447 FORMAT (I10, 8F10.5)                                              CCN54420
      IF(ITEBE.LE.ITBEMX) GO TO 210                                     CCN54430
      IF (KOUT.NE.0) THEN
      WRITE(8,455)                                                      CCN54440
  455 FORMAT(//15H0ITEBE=ITBEMX+1)                                      CCN54450
      ENDIF
      IF(KTRL2.EQ.1) GO TO 472                                          CCN54460
      ISERCH=IQI                                                        CCN54470
      GO TO 475                                                         CCN54480
  472 ISERCH=IQJ                                                        CCN54490
  475 WRITE(8,520) ISERCH,(GESMEM(I),I=1,ITEBE)                         CCN54500
      STOP                                                              CCN54510
  505 EGES=EGEST                                                        CCN54520
CCCCCC  ******             CONVERGENCE HAS BEEN ACHIEVED          ******CCN54530
      KGES=2                                                            CCN54540
      EGESIN=VSX-EGES                                                   CCN54550
      EGAB=DABS(EGES)                                                   CCN54560
      FKAPPA=UNITOK*DSQRT(RMAS*EGAB)                                    CCN54570
      FKAPIN=UNITOK*DSQRT(RMAS*EGESIN)                                  CCN54580
      KEX41=1                                                           CCN54590
      KEX40=0                                                           CCN54600
      CALL UNCPST                                                       CCN54610
      ENEPT=EGEST                                                       CCN54620
      VNEPT=VSX                                                         CCN54630
      KEX40=1                                                           CCN54640
      CALL UNCPST                                                       CCN54650
      IF (KOUT.NE.0) THEN
      WRITE(8,510) EGEST,VSX                                            CCN54660
  510 FORMAT( /  10X,36HFINAL VALUE OF THE PARAMETERS. B.E.=, F10.5,    CCN54670
     1  5H VSX=F10.5)    
      ENDIF                                                             CCN54680
      IF(KTRL2.EQ.1) GO TO 512                                          CCN54690
      ISERCH=IQI                                                        CCN54700
      GO TO 515                                                         CCN54710
  512 ISERCH=IQJ                                                        CCN54720
  515 IF (KOUT.NE.0) THEN
      WRITE(8,520) ISERCH,(GESMEM(I),I=1,ITEBE)                         CCN54730
  520 FORMAT(10X,A4,15H HAS VARIES AS ,10F9.4/29X,10F9.4)               CCN54740
      ENDIF
      IF(ITENOD.EQ.0) GO TO 527
      IF (KOUT.NE.0) THEN
      WRITE(8,525) ITENOD,(NODMEM(I1),I1=1,ITENOD)                      CCN54760
  525 FORMAT(/  15X,7HITENOD=,I2,11H,WITH NODE=,10I3)                   CCN54770
      ENDIF
  527 NXRMWF=NXRMP1                                                     CCN54780
      ABURX=DABS(URRMEX)                                                CCN54790
      ABURM=DABS(URRMIN)                                                CCN54800
      IF(ABURX.GT.1.D-20) GO TO 533                                     CCN54810
      IF(ABURM.GT.1.D-20) GO TO 545                                     CCN54820
      T1=0.0                                                            CCN54830
      GO TO 534                                                         CCN54840
  533 T1=URRMIN/URRMEX                                                  CCN54850
  534 DO 540 NX=NXRMWF,NXRAWF                                           CCN54860
  540 URSAVE(NX   )=URSAVE(NX   )*T1                                    CCN54870
  545 IIPA=0                                                            CCN54880
      SS=0.0                                                            CCN54890
      DO 575 NX=1,NXRAWF                                                CCN54900
      IF(IIPA.EQ.1) GO TO 552                                           CCN54910
      F1=4.0                                                            CCN54920
      IIPA=1                                                            CCN54930
      GO TO 553                                                         CCN54940
  552 F1=2.0                                                            CCN54950
      IIPA=0                                                            CCN54960
  553 SS=SS+(URSAVE(NX)**2)*F1                                          CCN54970
  575 CONTINUE                                                          CCN54980
      SS=SS*XMES2*.3333333333                                           CCN54990
      SQ1= DSQRT(1.  /SS)                                               CCN55000
      DO 580 NX=1,NXRAWF                                                CCN55010
  580 URSAVE(NX)=URSAVE(NX)*SQ1                                         CCN55020
      X=0.0                                                             CCN55030
      DO 647 NX=1,NXRAWF                                                CCN55040
      X=X+XMES2                                                         CCN55050
  647 URSAVE(NX   )=URSAVE(NX   )/X                                     CCN55060
      NXSTEP=1
      IF (KOUT.NE.0) THEN
      IF(KOPOT.EQ.0) NXSTEP=10                                          CCN55080
      XMOUT=DFLOAT(NXSTEP)*XMES2                                        CCN55090
      WRITE(8,649) XMOUT                                                CCN55100
  649 FORMAT( //,18H WAVE FUNCTIONS AT,F7.2,15H FERMI INTERVAL)         CCN55110
      WRITE(8,655) (URSAVE(NX),NX=NXSTEP,NXRAWF,NXSTEP)                 CCN55120
  655 FORMAT(10E13.4)
c      x=9.9
c      do 999 nx=100,130
c	x=x+xmes2
c	ursave(nx) = dlog(ursave(nx))
c	write(8,998) x, ursave(nx)
c999   continue
c998   format(2f10.5)    
      ENDIF
      IF(MOD(KTRL4,2).EQ.0) GO TO 667                                   CCN55140
CCCCCC  ******         MULTIPLICATION BY THE BINDING POTENTIAL    ******CCN55150
CCCCCC  ******          (FOR KTRL4 = 1 OR 3 )                     ******CCN55160
      LBTRTW=2*LBTR                                                     CCN55170
      LLSQ=LBTRTW*(LBTRTW+2)                                            CCN55180
      ISSQ=ISTW*(ISTW+2)                                                CCN55190
      JJSQ=JBTRTW*(JBTRTW+2)                                            CCN55200
      SPFAC=JJSQ-LLSQ-ISSQ                                              CCN55210
      SPFF=0.250*SPFAC                                                  CCN55220
      SPFAC=SPFF                                                        CCN55230
      IF(ISTW.EQ.0) GO TO 657                                           CCN55240
      FISTW=ISTW                                                        CCN55250
      SPFAC=SPFAC/FISTW                                                 CCN55260
  657 URSAVE(1)=URSAVE(1 )*(VCENTR(8)-VSPIN(8)*SPFAC)                   CCN55270
      URSAVE(2)=URSAVE(2 )*(VCENTR(12)-VSPIN(12)*SPFAC)                 CCN55280
      URSAVE(3)=URSAVE(3 )*(VCENTR(14)-VSPIN(14)*SPFAC)                 CCN55290
      DO 661 NX=4,NXRAWF                                                CCN55300
      MX=NX+12                                                          CCN55310
  661 URSAVE(NX)=URSAVE(NX)*(VCENTR(MX)-VSPIN(MX)*SPFAC)                CCN55320
      IF (KOUT.NE.0) THEN
      WRITE(8,662)                                                      CCN55330
  662 FORMAT( //,33H WAVE FUNCTION TIMES VCENTR+VSPIN)                  CCN55340
      WRITE(8,655) (URSAVE(NX),NX=1,NXRAWF,NXSTEP)
      ENDIF
      IF(TMAS.GE.4.0) GO TO 667                                         CCN55360
      S=0.0                                                             CCN55370
      IIPA=0                                                            CCN55380
      VFAI=0.0                                                          CCN55390
      DO 665 NX=1,NXRAWF                                                CCN55400
      S=S+XMES2                                                         CCN55410
      IF(IIPA.EQ.1) GO TO 663                                           CCN55420
      F1=4.0                                                            CCN55430
      IIPA=1                                                            CCN55440
      GO TO 664                                                         CCN55450
  663 F1=2.0                                                            CCN55460
      IIPA=0                                                            CCN55470
  664 VFAI=VFAI+F1*S**2*URSAVE(NX)                                      CCN55480
  665 CONTINUE                                                          CCN55490
      VFAI=VFAI*XMES2*0.333333333333                                    CCN55500
      VFF=VFAI**2*4.0*PAI                                               CCN55510
      VFAI=VFF
      IF (KOUT.NE.0) THEN
      WRITE(8,666)VFAI                                                  CCN55530
  666 FORMAT(//9X,49HCORRESPONDING ZERO RANGE NORMALIZATION D SQUARE =, CCN55540
     1    E20.5)         
      ENDIF
667   IF(KTRL4.LE.1) GO TO 700                                          CCN55560
CCCCCC  ******        DIVISION OF WAVE FUNCTION BY R**LBTR        ******CCN55570
CCCCCC  ******        IS MADE ( FOR KTRL4 = 2 OR 3 )              ******CCN55580
      IF(LBTR.EQ.0) GO TO 669                                           CCN55590
      R=0.                                                              CCN55600
      DO 668 NX=1,NXRAWF                                                CCN55610
      R=R+XMES2                                                         CCN55620
      URSAVE(NX)=URSAVE(NX)*(R**(-LBTR))                                CCN55630
  668 CONTINUE                                                          CCN55640
  669 IF (KOUT.NE.0) THEN
      IF(KTRL4.EQ.2) WRITE(8,672)                                       CCN55650
      IF(KTRL4.EQ.3) WRITE(8,673)                                       CCN55660
  672 FORMAT( //,31H WAVE FUNCTION TIMES R**(-LBTR))                    CCN55670
  673 FORMAT( //,50H WAVE FUNCTION TIMES VCENTR+VSPIN TIMES R**(-LBTR)) CCN55680
      WRITE(8,655) (URSAVE(NX),NX=1,NXRAWF,NXSTEP)
      ENDIF
  700 RETURN                                                            CCN55700
      END                                                               CCN55710
CCCCCC
CCCCCC
C====================================================================== CCN55730
      SUBROUTINE UNCPST                                                 CCN55740
C====================================================================== CCN55750
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CONST /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/BSX/   URSAVE(500),ENEPT,VNEPT,VSOR,DFNR,DFNSO,RZR,RZSO,   CCN55770
     1              RZC,XMES2,Q,FTR,NODER,LBTR,JBTRTW,NXRAWF,KTRL4,     CCN55780
     2              KOPOT,LBTRY,JBTWY                                   CCN55790
      COMMON/UNCPSA/KTRL3,ITBEMX,ACURCY,AMUPMU,                         CCN55800
     1              KTRL2,KTRL8,KEX2,KEX4,KEX40,KEX4L,KEX42,KEX43,      CCN55810
     2              TMAS,PMAS,ZZT,ZZP,RMAS,ZZ,XMES1,PERCNT,VSX,         CCN55820
     3              ISTW,NXRA,NXRM,NXRMP1,NXRMP2,NXRMP3,NXRMP4,NXRMP5,  CCN55830
     4              NODE,KGES,EGESRD,EGES,EGEST,DELGES,FKAPPA,FKAPIN,   CCN55840
     5              URRMIN,URRMEX,RGDLIN(3),RGDLEX(3),                  CCN55850
     6              XMEM(514),VCENTR(514),VSPIN(514),VCOULM(514),       CCN55860
     7              PFORM1(2),PFORM2(2),PFORM3(2),GESMEM(20),NODMEM(10) CCN55870
      COMMON/CNTRL/ KTRL(24),KEXCOM(24),KTLOUT(24),EXTCOM(50)
      DIMENSION     UR(4),FPRER(5),FPRERM(5)                            CCN55880
     
      PRE41=19.  /240.                                                  CCN55890
      PRE42=-2.  /5.                                                    CCN55900
      PRE43=97.  /120.                                                  CCN55910
      PRE44=-11. /15.                                                   CCN55920
      PRE45=299.  /240.                                                 CCN55930
      ECM=-EGES                                                         CCN55940
      FKEX40=KEX40                                                      CCN55950
      VENSFC=EGES/FKAPPA**2                                             CCN55960
      WNI=FKAPIN                                                        CCN55970
      MULSPN=ISTW+1                                                     CCN55980
      LBTRTW=2*LBTR                                                     CCN55990
      LLSQ=LBTRTW*(LBTRTW+2)                                            CCN56000
      ISSQ=ISTW*(ISTW+2)                                                CCN56010
      JJSQ=JBTRTW*(JBTRTW+2)                                            CCN56020
      LLTR=LBTRTW/2                                                     CCN56030
      LL=LLTR+1                                                         CCN56040
      ENSFAC=DFLOAT(LLTR*(LLTR+1))*VENSFC                               CCN56050
      SPFAC=JJSQ-LLSQ-ISSQ                                              CCN56060
      SPFF=0.250*SPFAC                                                  CCN56070
      SPFAC=SPFF                                                        CCN56080
      IF(ISTW.EQ.0) GO TO 425                                           CCN56090
      FISTW=ISTW                                                        CCN56100
      SPFAC=SPFAC/FISTW                                                 CCN56110
  425 NDFREP=(1-KEX40)*4+KEX40                                          CCN56120
      UR1M=1.0                                                          CCN56130
CCCCCC  ******       BEGINS NUMERICAL INTEGRATION WITH            ******CCN56140
CCCCCC  ******       FOUR STEP OR SINGLE STEP STORMER METHOD      ******CCN56150
      DO 730 ND=1,NDFREP                                                CCN56160
      FND=2.0  **(ND-1)                                                 CCN56170
      DDX=XMES1*FND*(1.  -FKEX40)-XMES2*FKEX40                          CCN56180
      DR=DDX*FKAPPA                                                     CCN56190
      DRSQ=(DR*DR)/EGES                                                 CCN56200
      DRRSQ=DABS(DRSQ)                                                  CCN56210
      DRSQ=DRRSQ                                                        CCN56220
      IF(KEX40.EQ.1) GO TO 510                                          CCN56230
CCCCCC  ******    PREPARATION FOR INTEGRATING INTERNAL SOLUTION   ******CCN56240
CCCCCC  ******    OUTWARDS WITH FOUR STEP STORMER METHOD          ******CCN56250
      K4COR2=4*ND-5                                                     CCN56260
      IF(ND.NE.1) GO TO 485                                             CCN56270
      UR(1)=.0                                                          CCN56280
      FPRERM(1)=.0                                                      CCN56290
      FPRER(1)=.0                                                       CCN56300
      UR(2)=.5  *((2.  *DDX*WNI)**LL)* DEXP(FACLOG(LL)-FACLOG(2*LL))    CCN56310
      IF(LL.NE.2) GO TO 473                                             CCN56320
      FPRER(1)=VENSFC*0.666666666666  *WNI*WNI*DRSQ                     CCN56330
  473 KEX42=0                                                           CCN56340
      NODCAL=0                                                          CCN56350
      FNODFC=1.0                                                        CCN56360
      IF(MOD(NODE,2).EQ.0) GO TO 478                                    CCN56370
      FNODFC=-1.0                                                       CCN56380
  478 UR(2)=UR(2)*FNODFC                                                CCN56390
      URSAVE(1)=UR(2)                                                   CCN56400
      FPRERM(1)=FPRER(1)*FNODFC                                         CCN56410
      FPRER(1)=FPRERM(1)                                                CCN56420
      X=0.0                                                             CCN56430
      NXIMIN=2                                                          CCN56440
      NXIMAX=8                                                          CCN56450
      NXPRCH=3                                                          CCN56460
      NNX=2                                                             CCN56470
      N3COR=0                                                           CCN56480
      GO TO 540                                                         CCN56490
  485 NXIMIN=5                                                          CCN56500
      X=3.  *DDX                                                        CCN56510
      NXIMAX=8                                                          CCN56520
      IF(ND.NE.4) GO TO 491                                             CCN56530
      NXIMAX=NXRMP5                                                     CCN56540
  491 DO 495 NQ=1,4                                                     CCN56550
      URSAVE(NQ)=URSAVE(2*NQ)                                           CCN56560
  495 FPRER(NQ)=FPRERM(NQ)*4.0                                          CCN56570
      FPRERM(1)=FPRER(1)                                                CCN56580
      FPRERM(2)=FPRER(3)                                                CCN56590
      UR(1)=UR1M                                                        CCN56600
      NXPRCH=5                                                          CCN56610
      NNX=3                                                             CCN56620
      GO TO 540                                                         CCN56630
  510 XRA=NXRA                                                          CCN56640
CCCCCC  ******    PREPARATION FOR INTEGRATING EXTERNAL SOLUTION   ******CCN56650
CCCCCC  ******     INWARDS WITH SINGLE STEP STORMER METHOD        ******CCN56660
      XRRA=XMES2*XRA                                                    CCN56670
      XRA=XRRA                                                          CCN56680
      NXRAM1=NXRA-1                                                     CCN56690
      FPRER(1)=0.0                                                      CCN56700
      UR(1)=FPRER(1)                                                    CCN56710
      URSAVE(NXRA)=UR(1)                                                CCN56720
      UR(2)=0.0001                                                      CCN56730
      URSAVE(NXRAM1)=UR(2)                                              CCN56740
      X=XRA                                                             CCN56750
      NXIMIN=2                                                          CCN56760
      NXIMAX=NXRA-NXRM-1                                                CCN56770
      N3COR =NXRA+13                                                    CCN56780
  540 DO 700 NX=NXIMIN,NXIMAX                                           CCN56790
      X=X+DDX                                                           CCN56800
      K2=(NX+K4COR2)*(1-KEX40)+(N3COR-NX)*KEX40                         CCN56810
      K3=NX        *(1-KEX40)+(NXRA -NX)*KEX40                          CCN56820
      AR2=VCENTR(K2)+VCOULM(K2)+ENSFAC/(X*X)-SPFAC*VSPIN(K2)-ECM        CCN56830
      BR2=AR2*DRSQ*UR(2)                                                CCN56840
      IF(ND-1) 551,551,554                                              CCN56850
  551 IF(NX-4) 552,552,554                                              CCN56860
  552 TERMR=BR2                                                         CCN56870
      FPRER(NX)=TERMR                                                   CCN56880
      GO TO 555                                                         CCN56890
  554 FPRER(5)=BR2                                                      CCN56900
      TERMR=PRE41*FPRER(1)+PRE42*FPRER(2)+PRE43*FPRER(3)                CCN56910
     1     +PRE44*FPRER(4)+PRE45*FPRER(5)                               CCN56920
  555 UR(3)=2.0  *UR(2)-UR(1)+TERMR                                     CCN56930
      URSAVE(K3)=UR(3)                                                  CCN56940
      UR(1)=UR(2)                                                       CCN56950
      UR(2)=UR(3)                                                       CCN56960
      IF(ND.EQ.NDFREP) GO TO 570                                        CCN56970
      IF(NX.NE.NXPRCH) GO TO 570                                        CCN56980
      FPRERM(NNX)=BR2                                                   CCN56990
      NXPRCH=NXPRCH+2                                                   CCN57000
      NNX=NNX+1                                                         CCN57010
      IF(NX.NE.7) GO TO 570                                             CCN57020
      UR1M=UR(1)                                                        CCN57030
  570 IF(ND.NE.1) GO TO 580                                             CCN57040
      IF(NX.LE.4) GO TO 585                                             CCN57050
  580 DO 583 NQ=1,4                                                     CCN57060
  583 FPRER(NQ)=FPRER(NQ+1)                                             CCN57070
  585 IF(KEX40.EQ.1) GO TO 700                                          CCN57080
      IF(UR(1)*UR(2).LT.0.0  ) NODCAL=NODCAL+1                          CCN57090
  700 CONTINUE                                                          CCN57100
  730 CONTINUE                                                          CCN57110
CCCCCC  ******      CALCULATE LOGARITHMIC DERIVATIVE OF WAVE      ******CCN57120
CCCCCC  ******      FUNCTIONS AT THE MATCHING RADIUS AND          ******CCN57130
CCCCCC  ******      ESTIMATE NEW SEARCH PARAMETER.                ******CCN57140
      IF(KEX40.EQ.1) GO TO 760                                          CCN57150
      KEX42=NODCAL-NODE                                                 CCN57160
      IF(KEX42.NE.0) GO TO 1000                                         CCN57170
  760 URRM=URSAVE(NXRMP3)                                               CCN57180
      RGDL=(8.0  *(URSAVE(NXRMP4)-URSAVE(NXRMP2))                       CCN57190
     1           -(URSAVE(NXRMP5)-URSAVE(NXRMP1)))/(12.  *XMES2*URRM)   CCN57200
      IF(KEX40.NE.0) GO TO 770                                          CCN57210
      RGDLIN(KGES)=RGDL                                                 CCN57220
      GO TO 773                                                         CCN57230
  770 RGDLEX(KGES)=RGDL                                                 CCN57240
  773 IF(KGES.NE.2) GO TO 910                                           CCN57250
      IF(KEX40.NE.0) GO TO 777                                          CCN57260
      URRMIN=URRM                                                       CCN57270
      GO TO 910                                                         CCN57280
  777 URRMEX=URRM                                                       CCN57290
  910 IF(KEX40.EQ.0) GO TO 1000                                         CCN57300
      IF(KGES.NE.3) GO TO 1000                                          CCN57310
      IF(KTRL2.EQ.1) GO TO 923                                          CCN57320
      DEN=.005  *EGEST                                                  CCN57330
      GO TO 925                                                         CCN57340
  923 DEN=.005  *VSX                                                    CCN57350
  925 DLEVIN=(RGDLIN(3)-RGDLIN(1))/DEN                                  CCN57360
      DLEVEX=(RGDLEX(3)-RGDLEX(1))/DEN                                  CCN57370
      DELGES=(RGDLIN(2)-RGDLEX(2))/(DLEVIN-DLEVEX)                      CCN57380
 1000 RETURN                                                            CCN57390
      END                                                               CCN57400
