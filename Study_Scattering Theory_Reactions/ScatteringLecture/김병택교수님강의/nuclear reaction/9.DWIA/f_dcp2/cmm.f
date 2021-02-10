      SUBROUTINE BSAXON                                                 CCN52480
CCCCCC
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
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
      COMMON/CNTRL /KTRLD(9),KTLOUT(24),NERUN,KRTYPE
      DATA          IQI,IQJ/4HB.E.,4HVSX  /                             CCN52610
      DATA          KK,KL,KM,KN / 1HX,4HVCEN,4HVSPI,4HCOUL /            CCN52620
CCCCCC
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
119   NXRM=(XBAR+.5*DFNR*FLOAT(NODE))/XMES2                             CCN53020
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
      PRINT 725, NXIMAX                                                 CCN53320
  725 FORMAT(2X,'NXIMAX.GT.514,=',I5)                                   CCN53330
      STOP                                                              CCN53340
  164 CONTINUE                                                          CCN53350
  165 DO 175 NX=NXIMIN,NXIMAX                                           CCN53360
      X=X+DX                                                            CCN53370
      XMEM(NX)=X                                                        CCN53380
      PFORM1(1)= EXP((X-XBAR  )/DFNR )                                  CCN53390
      PFORM1(2)= EXP((X-XBARSO)/DFNSO)                                  CCN53400
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
      FKAPPA=UNITOK   * SQRT(RMAS*EGES)                                 CCN53960
      FKAPIN=UNITOK   * SQRT(RMAS*EGESIN)                               CCN53970
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
      IF(ABS(DELGES/GEST).GT. .01) DELGES=.7*DELGES                     CCN54260
      IF(ABS(DELGES/GEST).GT. .05) DELGES=.7*DELGES                     CCN54270
      IF( ABS(DELGES/GEST).LE.PERCNT) GO TO 417                         CCN54280
      DELGES=   (GEST*DELGES/ ABS(DELGES))*PERCNT                       CCN54290
  417 GEST=GEST-DELGES                                                  CCN54300
      IF(KTRL2.EQ.1) GO TO 427                                          CCN54310
      EGEST=GEST                                                        CCN54320
      GO TO 440                                                         CCN54330
  427 VCOREC=GEST/(GEST+DELGES)                                         CCN54340
      VSX=VSX*VCOREC                                                    CCN54350
      DO 430 NX=1,NXRA12                                                CCN54360
      VCENTR(NX)=VCENTR(NX)*VCOREC                                      CCN54370
  430 CONTINUE                                                          CCN54380
  440 IF( ABS(DELGES/GEST).LE.ACURCY) GO TO 505                         CCN54390
      IF(KOPOT.NE.0)                                                    CCN54400
     1WRITE (6,447)  ITEBE,VSX,EGEST,EGES,EGESIN,DELGES,GEST            CCN54410
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
      EGAB=ABS(EGES)                                                    CCN54560
      FKAPPA=UNITOK*SQRT(RMAS*EGAB)                                     CCN54570
      FKAPIN=UNITOK*SQRT(RMAS*EGESIN)                                   CCN54580
      KEX41=1                                                           CCN54590
      KEX40=0                                                           CCN54600
      CALL UNCPST                                                       CCN54610
      ENEPT=EGEST                                                       CCN54620
      VNEPT=VSX                                                         CCN54630
      KEX40=1                                                           CCN54640
      CALL UNCPST                                                       CCN54650
      IF (KOUT.NE.0) THEN
      WRITE(8,510) EGEST,VSX                                            CCN54660
  510 FORMAT( /  10X,36HFINAL VALUE OF THE PARAMETERS. B.E.=,F10.5,      CCN54670
     1  5H VSX=,F10.5)    
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
  525 FORMAT(/  15X,7HITENOD=,I2,11H,WITH NODE=,10I3)                     CCN54770
      ENDIF
  527 NXRMWF=NXRMP1                                                     CCN54780
      ABURX=ABS(URRMEX)                                                 CCN54790
      ABURM=ABS(URRMIN)                                                 CCN54800
      IF(ABURX.GT.1.E-20) GO TO 533                                     CCN54810
      IF(ABURM.GT.1.E-20) GO TO 545                                     CCN54820
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
      SQ1= SQRT(1.  /SS)                                                CCN55000
      DO 580 NX=1,NXRAWF                                                CCN55010
  580 URSAVE(NX)=URSAVE(NX)*SQ1                                         CCN55020
      X=0.0                                                             CCN55030
      DO 647 NX=1,NXRAWF                                                CCN55040
      X=X+XMES2                                                         CCN55050
  647 URSAVE(NX   )=URSAVE(NX   )/X                                     CCN55060
      NXSTEP=1
      IF (KOUT.NE.0) THEN
      IF(KOPOT.EQ.0) NXSTEP=10                                          CCN55080
      XMOUT=FLOAT(NXSTEP)*XMES2                                         CCN55090
      WRITE(8,649) XMOUT                                                CCN55100
  649 FORMAT( //,18H WAVE FUNCTIONS AT,F7.2,15H FERMI INTERVAL)         CCN55110
      WRITE(8,655) (URSAVE(NX),NX=NXSTEP,NXRAWF,NXSTEP)                 CCN55120
  655 FORMAT(10E13.4)    
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


C====================================================================== CCN55730
      SUBROUTINE UNCPST                                                 CCN55740
C====================================================================== CCN55750
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
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
      COMMON/CNTRL /KTRLD(9),KTLOUT(24),NERUN,KRTYPE
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
      ENSFAC=FLOAT(LLTR*(LLTR+1))*VENSFC                                CCN56050
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
      DRRSQ=ABS(DRSQ)                                                   CCN56210
      DRSQ=DRRSQ                                                        CCN56220
      IF(KEX40.EQ.1) GO TO 510                                          CCN56230
CCCCCC  ******    PREPARATION FOR INTEGRATING INTERNAL SOLUTION   ******CCN56240
CCCCCC  ******    OUTWARDS WITH FOUR STEP STORMER METHOD          ******CCN56250
      K4COR2=4*ND-5                                                     CCN56260
      IF(ND.NE.1) GO TO 485                                             CCN56270
      UR(1)=.0                                                          CCN56280
      FPRERM(1)=.0                                                      CCN56290
      FPRER(1)=.0                                                       CCN56300
      UR(2)=.5  *((2.  *DDX*WNI)**LL)* EXP(FACLOG(LL)-FACLOG(2*LL))     CCN56310
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


C======================================================================CCCN40160
      SUBROUTINE OPT (IDCHNL,NEXPT)                                     CCN40170
C======================================================================CCCN40180
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4)
      PARAMETER(NXA=420,NXB=140,NIN=300)
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/COUCC/ FC(130),GC(130),FCP(130),GCP(130),EXSGRI(130),L5,L6,CCN40250
     1              ETA,SIGMAZ,RD,Z,KTOUT7,K9,STP,ACCR                  CCN40260
      COMPLEX       EXSGRI 
      COMMON/DWCC/  LDWMIR(4),LDWMXR(4),WN(4),LDWSTR(4),CE(4),CFUNIR(4),CCN40270
     1              ECM(4),XMESR(4),NXMXR(4),NXMIR(4),MXMAX,MXSTEP,     CCN40280
     2              TMASD(4),PMASD(4),TZD(4),PZD(4),XBARD(4),XBARID(4)  CCN40290
     3              ,WNXD(20),EXPD(20)                                  CCN40300
      COMMON/OPTL/  THMIN,THMAX,THINC                                   CCN40310
      COMMON/POTCC/ VD(4),WD(4),WSD(4),ARD(4),AID(4),AISD(4),RZRD(4),   CCN40330
     1              RZID(4),RZISD(4),RZCD(4),IDCH(4),VRIT(900),NXMN     CCN40340
     2              ,NXMX      
      COMMON/CNTRL/ KTRLD(9),KTLOUT(24),NERUN,KRTYPE                    CCN40360
      COMMON/CBST/  DISWY(NPS*NXA),USAVH(NHS*NXA),CY(LXA)
      COMPLEX       CY
      COMMON/DISW / DISWB(LXB*NXB),DISWA(LXA*NXA)
C      COMMON        DISWB(LXB*NXB),DISWA(LXA*NXA)
      COMPLEX       DISWA,DISWB                                         CCN40220
      COMMON        U(5),P(130),F1(180),CL(130),SIG(130),TETADG(180),   CCN40410
     1              SGRUTH(180),XRATIO(130),SMATX(130),NID(4),          CCN40420
     2              DISWT(900),DUMMY1(6578)                             CCN40430
      COMPLEX       SMATX,COR1,COR,B1,AA,CAB,D,WR,FCC,CC,CON1,DURM,     CCN40190
     1              CL,U,ZERO,TTI,SMAT,DISWT,
     2              ARI,URI1,URI2,URI3,BRI,VRIT                         CCN40210
      NEXPTA=NEXPT
      HBSQ=HBAR**2                                                      CCN40490
      D=1.0E-10                                                         CCN40510
      ZERO=(0.,0.)                                                      CCN40520
      TTI=(0.,1.)                                                       CCN40530
      STEP=0.0                                                          CCN40540
      ACCUR=0.0                                                         CCN40550
      STP=STEP                                                          CCN40560
      ACCR=ACCUR                                                        CCN40570
      ID=IDCHNL                                                         CCN40580
      NID(1)=84000                                                      CCN40590
      NID(2)=28000                                                      CCN40600
      LP1MI=LDWMIR(ID)+1                                                CCN40610
      LP1MX=LDWMXR(ID)+1                                                CCN40620
      LPSTP=LDWSTR(ID)                                                  CCN40630
      NMX=NXMXR(ID)                                                     CCN40640
      NMI=NXMIR(ID)                                                     CCN40650
      H=XMESR(ID)                                                       CCN40660
      RM=H*FLOAT(NMX)                                                   CCN40670
      NMX1=NMX+1                                                        CCN40680
      NMX2=NMX+2                                                        CCN40690
      NMXM2=NMX-2                                                       CCN40700
      TMST=TMASD(ID)                                                    CCN40710
      Z1=TZD(ID)                                                        CCN40720
      PMST=PMASD(ID)                                                    CCN40730
      Z2=PZD(ID)                                                        CCN40740
      IF(IDCH(ID).EQ.0) TMPM=TMST**.3333333333                          CCN40750
      IF(IDCH(ID).EQ.1) TMPM=TMST**.3333333333 + PMST**.3333333333      CCN40760
      L4=LP1MX                                                          CCN40770
      RO=TMPM*RZRD(ID)                                                  CCN40780
      RI=TMPM*RZID(ID)                                                  CCN40790
      RMP2=RM+2.*XMESR(ID)                                              CCN40800
      NONX=NMX-NMI+1                                                    CCN40810
      TM=TMST*AMAS                                                      CCN40820
      PM=PMST*AMAS                                                      CCN40830
      E=ECM(ID)                                                         CCN40850
      C=PM*TM/(PM+TM)                                                   CCN40860
      Q=WN(ID)                                                          CCN40870
      ETA=C*Z2*Z1*HBAR*FINE/(Q*HBSQ)                                    CCN40880
      ETA=CE(ID) 
      IF (KTRLD(2).EQ.1) ETA=0.0
CCCCCC     FOR THE PLANE WAVE APPROXIMATION
      Z=RM*Q                                                            CCN40900
      L5=L4+2                                                           CCN40910
      RD=H*Q                                                            CCN40920
      KTOUT7=KTLOUT(7)                                                  CCN40930
      CALL FLGLCH                                                       CCN40940
      SIG(1)=SIGMAZ                                                     CCN40950
      CE(ID)=ETA                                                        CCN40960
      XBARD(ID)=RO                                                      CCN40970
      XBARID(ID)=RI                                                     CCN40980
      DO 310 L=1,L4                                                     CCN40990
      SS=L                                                              CCN41000
310   SIG(L+1)=SIG(L)+ATAN(ETA/SS)                                      CCN41010
      DO 340 L=1,L4                                                     CCN41020
      FCP(L)=Q*FCP(L)                                                   CCN41030
340   GCP(L)=Q*GCP(L)                                                   CCN41040
      NXMN=1                                                            CCN41050
      NXMX=NMX2                                                         CCN41060
      IF(KTRLD(3).EQ.1) GO TO 341
      CALL POTEN(ID,1)                                                  CCN41070
  341	XSTEP=XMESR(ID)*10.0
	WRITE(8,342) XSTEP
	WRITE(8,343) (VRIT(NX),NX=1,NXMX,10)            
  342 FORMAT(/1H ,'  OM POTENTIAL IN STEPS OF', F8.3/)
  343 FORMAT(1H ,3X,10E10.3)
      TMU=-(Q**2)/E                                                     CCN41100
      LM=0                                                              CCN41150
      NXMF=0                                                            CCN41160
      LMTL=(LP1MX-LP1MI)/LPSTP+1                                        CCN41170
      NXMTL=LMTL*NONX                                                   CCN41180
      IF(ID.EQ.1.AND.NXMTL.LE.NID(1)) GO TO 350                         CCN41190
      IF(ID.EQ.2.AND.NXMTL.LE.NID(2)) GO TO 350                         CCN41200
      IF(ID.EQ.3.AND.NXMTL.LE.NID(1)) GO TO 350                         CCN41210
      IF(ID.EQ.4.AND.NXMTL.LE.NID(1)) GO TO 350                         CCN41220
      WRITE(8,650) ID,NXMTL,NID(ID)                                     CCN41230
  650 FORMAT(2X,'FOR CHANNEL',I3,2X,'NXMTL=',I6,2X,'.GT. THAN',I6)      CCN41240
      STOP                                                              CCN41250
  350 CONTINUE                                                          CCN41260
      DO 220 LP1=LP1MI,LP1MX,LPSTP                                      CCN41270
      LM=LM+1                                                           CCN41280
      LP=LP1-1                                                          CCN41290
      CALL DISWAVE(ID,LM,LP,NMI,NMX,TMU,E,H)
      NXMBAS=(LM-1)*NONX
      NXM=NXMBAS                                                        CCN41610
      R=0.0                                                             CCN41620
      DO 210 NX=NMI,NMX                                                 CCN41630
      R=H*FLOAT(NX)                                                     CCN41640
      NXM=NXM+1                                                         CCN41650
      IF(ID.EQ.2) DISWB(NXM)=DISWT(NX)*R*R                              CCN41670
      IF(ID.EQ.1) DISWA(NXM)=DISWT(NX)*R*R                              CCN41680
210   CONTINUE                                                          CCN41710
      IF(ID.NE.4) GO TO 214                                             CCN41720
      NXMT=NXMBAS+NONX-2                                                CCN41730
      URI1=DISWT(NXMT)                                                  CCN41740
      COR1=URI1/CABS(URI1)                                              CCN41750
      CY(LM)=COR1                                                       CCN41760
      NXMT=NXMBAS                                                       CCN41770
      DO 212 NX=1,NONX                                                  CCN41780
      NXMT=NXMT+1                                                       CCN41790
      DISWY(NXMT)=DISWT(NXMT)/COR1                                      CCN41800
  212 CONTINUE                                                          CCN41810
  214 IF(KTLOUT(2).EQ.0) GO TO 220                                      CCN41820
      KTL2=KTLOUT(2)                                                    CCN41880
      NXM=NXMBAS                                                        CCN41890
      DO 215 NX=KTL2,NONX,KTL2                                          CCN41900
      NXM=NXM+KTL2                                                      CCN41910
      ARI=DISWT(NX)                                                     CCN41940
      DISAB=CABS(ARI)                                                   CCN41950
      AR=REAL(ARI)                                                      CCN41960
      AI=AIMAG(ARI)                                                     CCN41970
  215 CONTINUE                                                          CCN42010
      NX1=NXMBAS+KTL2                                                   CCN42020
      NX2=NXMBAS+NONX                                                   CCN42030
      WRITE(8,232) LP,KTL2,(DISWT(NX),NX=KTL2,NONX,KTL2)                CCN42040
      NXM=NXMBAS                                                        CCN42050
  220 CONTINUE 
      IF(KTLOUT(5).EQ.0) GO TO 420                                      CCN42190
      WRITE(8,610)                                                      CCN42200
      WRITE(8,620) (CL(L),L=1,L4)                                       CCN42210
      WRITE(8,480)                                                      CCN42240
      DO 370 L=1,L4                                                     CCN42250
      LL=L-1                                                            CCN42260
      SMATX(L)=1.+2.*TTI*CL(L)                                          CCN42270
      SMAT=SMATX(L)                                                     CCN42280
      X1=REAL(SMAT)                                                     CCN42290
      X2=AIMAG(SMAT)                                                    CCN42300
      PHASE=ATAN2(X2,X1)                                                CCN42310
      ABSM=(X1**2)+(X2**2)                                              CCN42320
      ABSMAT=SQRT(ABSM)                                                 CCN42330
      WRITE(8,490) LL,SMAT,ABSMAT,PHASE                                 CCN42340
370   CONTINUE                                                          CCN42350
      THETA=THMIN-THINC                                                 CCN42380
      NTH=(THMAX-THMIN)/THINC + 1.00001                                 CCN42390
      NHAF=NTH/2                                                        CCN42400
      DO 410 M=1,NTH                                                    CCN42410
      THETA=THETA+THINC                                                 CCN42420
      TETADG(M)=THETA                                                   CCN42430
      TH=THETA*(PI/180.)                                                CCN42440
      P(1)=1.                                                           CCN42450
      P(2)=COS(TH)                                                      CCN42460
      DO 390 L=3,L4                                                     CCN42470
      S=L                                                               CCN42480
390   P(L)=(1./(S-1.))*((2.*S-3.)*COS(TH)*P(L-1)-(S-2.)*P(L-2))         CCN42490
      WR=0.                                                             CCN42500
      DO 400 N=1,L4                                                     CCN42510
      S=N                                                               CCN42520
      CC=CMPLX(0.,2.*SIG(N))                                            CCN42530
400   WR=(1./Q)*(2.*S-1.)*CL(N)*P(N)*EXP(CC)+WR                         CCN42540
      AB=(SIN(TH/2.))**2                                                CCN42550
      FCC=(ETA/(2.*Q*AB))                                               CCN42560
      AA=CMPLX(0.,-ETA*LOG(AB)+2.*SIGMAZ)                               CCN42570
      FCC=FCC*EXP(AA)                                                   CCN42580
      SGRUTH(M)=10.*FCC*CONJG(FCC)                                      CCN42590
      F1(M)=10.*(WR-FCC)*CONJG(WR-FCC)                                  CCN42600
      XRATIO(M)=F1(M)/SGRUTH(M)                                         CCN42630
410   CONTINUE                                                          CCN42640
      WRITE(8,630) E                                                    CCN42650
      WRITE(8,640) (TETADG(M),F1(M),XRATIO(M),TETADG(M+NHAF),F1(M+NHAF),CCN42660
     1 XRATIO(M+NHAF),M=1,NHAF) 
  420 CONTINUE                                                          CCN42670
                                                                        CCN42680
      RETURN                                                            CCN42690
  232 FORMAT(/' DW (AMPL) FOR LA=',I2,' IN EVERY NX=',I2,/(5X,6E10.3))  CCN42700
  233 FORMAT(/' DW (AMPL) FOR LB=',I2,' IN EVERY NX=',I2,/(5X,6E10.2))  CCN42710
  480 FORMAT(2(/),4X,'L',21X,'S MATRIX',23X,'ABSMAT',14X,'PHASE',/)     CCN42720
  490 FORMAT(I5,4E20.5)                                                 CCN42730
  610 FORMAT(2(/),3X,29HTHE C MATRIX FOR L=0 TO LMAX,,1(/))              CCN42760
  620 FORMAT(2E15.5)                                                    CCN42770
  630 FORMAT(2(/),3X,30HTHE ELASTIC CROSS SECTION     ,4HECM=,
     1 E15.5,//, 10X,'ANGLE',10X,'SIGMA',10X,'RATIO',
     2 10X,'ANGLE',10X,'SIGMA',10X,'RATIO'/)                            
  640 FORMAT(6E15.5)                                                    CCN42800
      END                                                               CCN42810


                                                                        CCN48660
C======================================================================CCCN48670
      SUBROUTINE POTEN (IDCHNL,IDC)                                     CCN48680
C======================================================================CCCN48690
      PARAMETER(NXA=420,NXB=140,NIN=300)
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/POTCC/ VD(4),WD(4),WSD(4),ARD(4),AID(4),AISD(4),RZRD(4),   CCN48730
     1              RZID(4),RZISD(4),RZCD(4),IDCH(4),VRIT(900),NXMN,    CCN48740
     2              NXMX                                                CCN48750
      COMPLEX       VRIT
      COMMON/DWCC/  LDWMIR(4),LDWMXR(4),WN(4),LDWSTR(4),CE(4),DFUNIR(4),CCN48760
     1              ECM(4),XMESR(4),NXMXR(4),NXMIR(4),MXMAX,MXSTEP,     CCN48770
     2              TMASD(4),PMASD(4),TZD(4),PZD(4),XBARD(4),XBARID(4)
     3       	    ,EXTRA(40)
      COMMON/OPTO/  FOLD(NXA)                                           CCN48820
      COMMON/CNTRL/ KTRLD(9),KTLOUT(24),NERUN,KRTYPE                    CCN48720
      DIMENSION     RTS(96),WGT(96)
      DIMENSION     VRITR(900),VRITI(900),VRITWC(900),                  CCN48830
     1              OPTEXR(900),OPTEXI(900)                             CCN48840
      COMPLEX       TTI,U1,OPTEX(900),VRITWC,TTR                        CCN48850
      TTI=(0.,1.)                                                       CCN48860
      TTR=(1.,0.)                                                       CCN48870
      RENO=1.00                                                         CCN48880
      IF (ID.EQ.1) GO TO 770
      NNN=NXMXR(1)+2                                                    CCN48900
      DO 878 J=1,900                                                    CCN48910
      VRITWC(J)=(0.,0.)                                                 CCN48920
      VRITR(J)=0.                                                       CCN48930
  878 VRIT(J)=(0.,0.)                                                   CCN48940
  770 CONTINUE
      ID=IDCHNL                                                         CCN48950
      TMST=TMASD(ID)                                                    CCN48960
      PMST=PMASD(ID)                                                    CCN48970
      H=XMESR(ID)                                                       CCN48980
      IF(IDCH(ID).EQ.0) TMPM=TMST**.3333333333                          CCN48990
      IF(IDCH(ID).EQ.1) TMPM=TMST**.3333333333 + PMST**.3333333333      CCN49000
      VO=VD(ID)                                                         CCN49010
      WO=WD(ID)                                                         CCN49020
      RO=TMPM*RZRD(ID)                                                  CCN49030
      RI=TMPM*RZID(ID)                                                  CCN49040
      RC=TMPM*RZCD(ID)                                                  CCN49050
      A=ARD(ID)                                                         CCN49060
      B=AID(ID)                                                         CCN49070
      TMST=TMASD(ID)                                                    CCN49100
      Z1=TZD(ID)                                                        CCN49110
      PMST=PMASD(ID)                                                    CCN49120
      Z2=PZD(ID)                                                        CCN49130
      RN=TMPM*RZCD(ID)                                                  CCN49140
      VOS=VD(ID+2)           ! TEMP
      ROS=TMPM*RZRD(ID+2)    ! TEMP
      AOS=ARD(ID+2)          ! TEMP
      RZIS=TMPM*RZISD(ID)                                               CCN49150
      AIS=AISD(ID)                                                      CCN49160
      WS=WSD(ID)                                                        CCN49170
      CH=Z1*Z2                                                          CCN49180
      IF(IDC.EQ.0) CH=0.                                                CCN49190
      EE=HBAR*FINE                                                      CCN49200
                                                                        CCN49210
      VCC=.5*CH*EE/RN                                                   CCN49220
      IK=0                                                              CCN49230
      RMIN=H*FLOAT(NXMN)                                                CCN49240
      R=RMIN-H                                                          CCN49250
      DO 270 NX=NXMN,NXMX                                               CCN49260
      IK=IK+1                                                           CCN49270
      R=R+H                                                             CCN49280
      VR=-VO/(1.+EXP((R-RO)/A))                                         CCN49290
      VR1=VOS/(1.+EXP((R-ROS)/AOS))**2   ! TEMP
      VR=VR+VR1                          ! TEMP
      VI1=-WO/(1.+EXP((R-RI)/B))                                        CCN49300
      EXPD=EXP((R-RZIS)/AIS)                                            CCN49310
      EX1=1./(1.+EXPD)                                                  CCN49320
      VI2=-4.*WS*(EX1**2)*EXPD                                          CCN49330
      VI=VI1+VI2                                                        CCN49340
      IF (R-RN) 230,230,240                                             CCN49350
  230 VVC=VCC*(3.-((R/RN)**2))                                          CCN49360
      GO TO 250                                                         CCN49370
  240 VVC=CH*EE/R                                                       CCN49380
  250 CONTINUE                                                          CCN49390
      VRITR(IK)=VVC                                                     CCN49400
      VRITI(IK)=VI                                                      CCN49410
      VRIT(IK)=VR+VVC+VI*TTI                                            CCN49420
      IF(KTRLD(2).EQ.1) VRIT(IK)=(0.0,0.0)              
      VRITWC(IK)=(VR+VI*TTI)*RENO                                       CCN49450
  270 CONTINUE                                                          CCN49460
      KMAS=PMST+0.1                                                     CCN49610
      IF(KMAS.LE.1) GO TO 266                                           CCN49620
      IF(KTRLD(2).GE.1) GO TO 266                                       CCN49630
      IF(KTRLD(3).EQ.0) GO TO 266
      N2MAX=NXMXR(3) 
      NGAUS=24                                                          CCN50010
      CALL GAUSF(NGAUS,RTS,WGT)                                         CCN50020
      N1MAX=NXMXR(4)
      KEX1M1=NXMX-1                                                     CCN50030
      NXMAX=NXMXR(1)
      DO 4 J=NXMN,NXMX                                                  CCN50040
      OPTEXR(J)=0.                                                      CCN50050
      OPTEXI(J)=0.                                                      CCN50060
    4 OPTEX(J)=(0.,0.)                                                  CCN50070
      F=2.0
      DO 1 N=NXMN,NXMX                                                  CCN50080
      R=XMESR(ID)*FLOAT(N)                                              CCN50090
      RSQ=R*R                                                           CCN50100
      F=2.                                                              CCN50110
      DO 2 N2=1,N2MAX                                                   CCN50120
      F=6.-F                                                            CCN50130
      R2=XMESR(3)*FLOAT(N2)                                             CCN50140
      R2SQ=R2*R2                                                        CCN50150
      R2R=R2*R*2.                                                       CCN50160
      DO 3 NG=1,NGAUS                                                   CCN50170
      ROOT=RTS(NG)                                                      CCN50180
      R1P=SQRT(RSQ+R2SQ-R2R*ROOT)                                       CCN50190
      X1P=R1P/XMESR(ID)                                                 CCN50200
      N1P=X1P                                                           CCN50210
      IF(N1P.LT.1) N1P=1                                                CCN50220
      IF(N1P.GT.KEX1M1) GOTO 3                                          CCN50230
      FN1P=N1P                                                          CCN50240
      U1=VRITWC(N1P)+(VRITWC(N1P+1)-VRITWC(N1P))*(X1P-FN1P)             CCN50250
      U1=U1*WGT(NG)*FOLD(N2)*R2SQ*XMESR(3)*F/(3.0*2.0)                  CCN50260
      OPTEX(N)=OPTEX(N)+U1                                              CCN50270
    3 CONTINUE                                                          CCN50280
    2 CONTINUE                                                          CCN50290
    1 CONTINUE                                                          CCN50300
      DO 18 J=NXMN,NXMX                                                 CCN50310
      VRIT(J)=(OPTEX(J)+VRITR(J))                                       CCN50320
   18 CONTINUE                                                          CCN50330
  266 CONTINUE                                                          CCN50380
      RETURN                                                            CCN50430
      END


C======================================================================CCCN50460
      SUBROUTINE FLGLCH                                                 CCN50470
C======================================================================CCCN50480
      COMPLEX       EXSGRI                                              CCN50490
      REAL          K,K1,K2,K3,K4,M1,M2,M3,M4                           CCN50500
      COMMON/CNTRL/ KTRLD(9),KTLOUT(24),NERUN,KRTYPE                    CCN50510
      COMMON/COUCC/ FC(130),GC(130),FCP(130),GCP(130),EXSGRI(130),      CCN50520
     1              MAXL,N1,EETA,SIGMAX,RD,RHO,KTOUT7,K9,STEP,ACCUR     CCN50530
      ETA=EETA                                                          CCN50540
      IF(KTRLD(2).EQ.1) ETA=0.0                                         CCN50560
      IF(ETA-10.) 280,280,290                                           CCN50580
280   ETA2=ETA**2                                                       CCN50590
      ETA2A=2.*ETA                                                      CCN50600
      ETA6=ETA2+16.                                                     CCN50610
      SIGMAO=-(ETA/(12.*ETA6))*(1.+(ETA2-48.)/(30.*(ETA6**2))+((ETA2-   CCN50620
     1 160.)*ETA2+1280.)/(105.*(ETA6**4)))-ETA+(ETA/2.)*LOG(ETA6)+3.5*  CCN50630
     2 ATAN(.25*ETA)-(ATAN(ETA)+ATAN(.5*ETA)+ATAN(ETA/3.))              CCN50640
      GO TO 300                                                         CCN50650
290   EINV1=1./ETA                                                      CCN50660
      EINV2=EINV1*EINV1                                                 CCN50670
      EINV3=EINV1*EINV2                                                 CCN50680
      EINV5=EINV3*EINV2                                                 CCN50690
      EINV7=EINV5*EINV2                                                 CCN50700
      EINV9=EINV7*EINV2                                                 CCN50710
      SIGMAO=.7853981634+ETA*LOG(ETA)-ETA-(.08333333333*EINV1+          CCN50720
     1  .00277777777*EINV3+.00079365079*EINV5+.00059523810*EINV7+       CCN50730
     2  .00084175084*EINV9)                                             CCN50740
300   MODTPI=SIGMAO/6.2831853072                                        CCN50750
      FMODTP=MODTPI                                                     CCN50760
      SIGMAZ=SIGMAO-6.2831853072*FMODTP                                 CCN50770
      SIGMAX=SIGMAZ                                                     CCN50780
      MINL = 0                                                          CCN50790
      F=0.0                                                             CCN50800
      FP=0.0                                                            CCN50810
      PACE=STEP                                                         CCN50820
      ACC=ACCUR                                                         CCN50830
      IF(PACE.LT.100.0)PACE=100.0                                       CCN50840
      IF(ACC.LT.1.0E-15.OR.ACC.GT.1.0E-6) ACC = 1.0E-6                  CCN50850
      R=RHO                                                             CCN50860
      KTR=1                                                             CCN50870
      LMAX=MAXL                                                         CCN50880
      LMIN1=MINL+1                                                      CCN50890
      XLL1=FLOAT(MINL*LMIN1)                                            CCN50900
      ETA2=ETA*ETA                                                      CCN50910
      TURN=ETA+SQRT(ETA2+XLL1)                                          CCN50920
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.0E-6)KTR=-1                        CCN50930
      KTRP=KTR                                                          CCN50940
      GOTO2                                                             CCN50950
1     R=TURN                                                            CCN50960
      TF=F                                                              CCN50970
      TFP=FP                                                            CCN50980
      LMAX=MINL                                                         CCN50990
      KTRP=1                                                            CCN51000
2     ETAR=ETA*R                                                        CCN51010
      RHO2=R*R                                                          CCN51020
      PL=FLOAT(LMAX+1)                                                  CCN51030
      PMX=PL+0.5                                                        CCN51040
      FP=ETA/PL+PL/R                                                    CCN51050
      DK=ETAR*2.00                                                      CCN51060
      DEL=0.0                                                           CCN51070
      D=0.0                                                             CCN51080
      F=1.0                                                             CCN51090
      K=(PL*PL-PL+ETAR)*(2.0*PL-1.0)                                    CCN51100
      ABPL=ABS(PL**2+PL+ETAR)                                           CCN51110
      IF(ABPL.GT.1.0E-30) GO TO 3                                       CCN51120
      R=R+1.0E-6                                                        CCN51130
      GOTO2                                                             CCN51140
3     H=(PL*PL+ETA2)*(1.0-PL*PL)*RHO2                                   CCN51150
      K=K+DK+PL*PL*6.0                                                  CCN51160
      D=1.0/(D*H+K)                                                     CCN51170
      DEL=DEL*(D*K-1.0)                                                 CCN51180
      IF(PL.LT.PMX)DEL=-R*(PL*PL+ETA2)*(PL+1.0)*D/PL                    CCN51190
      PL=PL+1.0                                                         CCN51200
      FP=FP+DEL                                                         CCN51210
      IF(D.LT.0.0)F=-F                                                  CCN51220
      IF(PL.GT.20000.)GOTO11                                            CCN51230
      IF(ABS(DEL/FP).GE.ACC)GOTO3                                       CCN51240
      FP=F*FP                                                           CCN51250
      IF(LMAX.EQ.MINL)GOTO5                                             CCN51260
      FC(LMAX+1)=F                                                      CCN51270
      FCP(LMAX+1)=FP                                                    CCN51280
      L=LMAX                                                            CCN51290
      DO 4 LP=LMIN1,LMAX                                                CCN51300
      PL = FLOAT(L)                                                     CCN51310
      GC(L+1)=ETA/PL+PL/R                                               CCN51320
      GCP(L+1)=SQRT(ETA2+PL*PL)/PL                                      CCN51330
      FC(L)=(GC(L+1)*FC(L+1)+FCP(L+1))/GCP(L+1)                         CCN51340
      FCP(L)=GC(L+1)*FC(L)-GCP(L+1)*FC(L+1)                             CCN51350
4     L=L-1                                                             CCN51360
      F=FC(LMIN1)                                                       CCN51370
      FP=FCP(LMIN1)                                                     CCN51380
5     IF(KTRP.EQ.-1)GOTO1                                               CCN51390
      P=0.0                                                             CCN51400
      Q=R-ETA                                                           CCN51410
      PL=0.0                                                            CCN51420
      AR=-(ETA2+XLL1)                                                   CCN51430
      AI=ETA                                                            CCN51440
      BR=2.0*Q                                                          CCN51450
      BI=2.0                                                            CCN51460
      WI=2.0*ETA                                                        CCN51470
      DR=BR/(BR*BR+BI*BI)                                               CCN51480
      DI=-BI/(BR*BR+BI*BI)                                              CCN51490
      DP=-(AR*DI+AI*DR)                                                 CCN51500
      DQ=(AR*DR-AI*DI)                                                  CCN51510
6     P=P+DP                                                            CCN51520
      Q=Q+DQ                                                            CCN51530
      PL=PL+2.0                                                         CCN51540
      AR=AR+PL                                                          CCN51550
      AI=AI+WI                                                          CCN51560
      BI=BI+2.0                                                         CCN51570
      D=AR*DR-AI*DI+BR                                                  CCN51580
      DI=AI*DR+AR*DI+BI                                                 CCN51590
      T=1.0/(D*D+DI*DI)                                                 CCN51600
      DR=T*D                                                            CCN51610
      DI=-T*DI                                                          CCN51620
      H=BR*DR-BI*DI-1.00                                                CCN51630
      K=BI*DR+BR*DI                                                     CCN51640
      T=DP*H-DQ*K                                                       CCN51650
      DQ=DP*K+DQ*H                                                      CCN51660
      DP=T                                                              CCN51670
      IF(PL.GT.46000.)GOTO11                                            CCN51680
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC)GOTO6                   CCN51690
      P=P/R                                                             CCN51700
      Q=Q/R                                                             CCN51710
      G=(FP-P*F)/Q                                                      CCN51720
      GP=P*G-Q*F                                                        CCN51730
      W=1.0/SQRT(FP*G-F*GP)                                             CCN51740
      G=W*G                                                             CCN51750
      GP=W*GP                                                           CCN51760
      IF(KTR.EQ.1)GOTO8                                                 CCN51770
      F=TF                                                              CCN51780
      FP=TFP                                                            CCN51790
      LMAX=MAXL                                                         CCN51800
      IF(RHO.LT.0.2*TURN)PACE=999.00                                    CCN51810
      R3=1.0/3.0                                                        CCN51820
      H=(RHO-TURN)/(PACE+1.0)                                           CCN51830
      H2=0.5*H                                                          CCN51840
      I2=PACE+0.001                                                     CCN51850
      ETAH=ETA*H                                                        CCN51860
      H2LL=H2*XLL1                                                      CCN51870
      S=(ETAH+H2LL/R)/R-H2                                              CCN51880
7     RH2=R+H2                                                          CCN51890
      T=(ETAH+H2LL/RH2)/RH2-H2                                          CCN51900
      K1=H2*GP                                                          CCN51910
      M1=S*G                                                            CCN51920
      K2=H2*(GP+M1)                                                     CCN51930
      M2=T*(G+K1)                                                       CCN51940
      K3=H*(GP+M2)                                                      CCN51950
      M3=T*(G+K2)                                                       CCN51960
      RM3=2.0*M3                                                        CCN51970
      M3=RM3                                                            CCN51980
      K4=H2*(GP+M3)                                                     CCN51990
      RH=R+H                                                            CCN52000
      S=(ETAH+H2LL/RH)/RH-H2                                            CCN52010
      M4=S*(G+K3)                                                       CCN52020
      G=G+(K1+K2+K3+K2+K4)*R3                                           CCN52030
      GP=GP+(M1+M2+M2+M3+M4)*R3                                         CCN52040
      IF(G.GT.1.0E30.OR.GP.GT.1.0E30)GOTO11                             CCN52050
      R=RH                                                              CCN52060
      I2=I2-1                                                           CCN52070
      IF(I2.GE.0)GOTO7                                                  CCN52080
      W=1.0/(FP*G-F*GP)                                                 CCN52090
8     GC(LMIN1)=G                                                       CCN52100
      GCP(LMIN1)=GP                                                     CCN52110
      IF(LMAX.EQ.MINL) GO TO 10                                         CCN52120
      IF(LMAX.EQ.MINL) GO TO 10                                         CCN52130
      DO 9 L=LMIN1,LMAX                                                 CCN52140
      T=GC(L+1)                                                         CCN52150
      GC(L+1)=(GC(L)*GC(L+1)-GCP(L))/GCP(L+1)                           CCN52160
      GCP(L+1)=GC(L)*GCP(L+1)-GC(L+1)*T                                 CCN52170
      FC(L+1)=W*FC(L+1)                                                 CCN52180
9     FCP(L+1)=W*FCP(L+1)                                               CCN52190
      FC(LMIN1)=FC(LMIN1)*W                                             CCN52200
      FCP(LMIN1)=FCP(LMIN1)*W                                           CCN52210
      GO TO 107                                                         CCN52220
10    FC(LMIN1)=W*F                                                     CCN52230
      FCP(LMIN1)=W*FP                                                   CCN52240
      GO TO 107                                                         CCN52250
11    W=0.0                                                             CCN52260
      G=0.0                                                             CCN52270
      GP=0.0                                                            CCN52280
      GOTO8                                                             CCN52290
107   CONTINUE                                                          CCN52300
      IF(KTOUT7.EQ.0) GO TO 330                                         CCN52310
      L4=LMAX                                                           CCN52320
      WRITE(8,560)                                                      CCN52330
      WRITE(8,570) (FC(L),L=1,L4)                                       CCN52340
      WRITE(8,580) (FCP(L),L=1,L4)                                      CCN52350
      WRITE(8,590) (GC(L),L=1,L4)                                       CCN52360
      WRITE(8,600) (GCP(L),L=1,L4)                                      CCN52370
560   FORMAT(//,5X,43H THE COULOMB WAVE FUNCTIONS FOR L=0 TO LMAX,2(/)) CCN52380
570   FORMAT(6HFC(L)=,5E15.5)                                           CCN52390
580   FORMAT(7HFCP(L)=,5E15.5)                                          CCN52400
590   FORMAT(6HGC(L)=,5E15.5)                                           CCN52410
600   FORMAT(7HGCP(L)=,5E15.5)                                          CCN52420
330   CONTINUE                                                          CCN52430
      RETURN                                                            CCN52440
      END                                                               CCN52450


      SUBROUTINE DISWAVE(ID,LM,LP,NMI,NMX,TMU,E,H)
CCCCCC
      PARAMETER(NHS=9,NPS=9,LXA=150,LXB=150,LPH=10,LMM=4)
      PARAMETER(NXA=420,NXB=140,NIN=300)
CCCCCC
      COMMON/CNST  /FACLOG(500),PI,HBAR,AMAS,WNUNIT,FINE
      COMMON/COUCC/ FC(130),GC(130),FCP(130),GCP(130),EXSGRI(130),L5,L6,
     1              ETA,SIGMAZ,RD,Z,KTOUT7,K9,STP,ACCR                  
      COMPLEX       EXSGRI 
      COMMON/POTCC/ VD(4),WD(4),WSD(4),ARD(4),AID(4),AISD(4),RZRD(4),   CCN40330
     1              RZID(4),RZISD(4),RZCD(4),IDCH(4),VRIT(900),NXMN     CCN40340
     2              ,NXMX      
      COMMON/DISW / DISWB(LXB*NXB),DISWA(LXA*NXA)
      COMPLEX       DISWA,DISWB                                         CCN40220
      COMMON        U(5),P(130),F1(180),CL(130),SIG(130),TETADG(180),   CCN40410
     1              SGRUTH(180),XRATIO(130),SMATX(130),NID(4),          CCN40420
     2              DISWT(900),DUMMY1(6578) 
      COMPLEX       SMATX,COR1,COR,B1,AA,CAB,D,WR,FCC,CC,CON1,DURM,     CCN40190
     1              CL,U,ZERO,TTI,SMAT,DISWT,
     2              ARI,URI1,URI2,URI3,BRI,VRIT
      DOUBLE PRECISION  TTMU,EE,R,RL,YMG,YMP2,HH,HSQ,HSQ12,SG,
     1                  FFC,FFCP,GGC,GGCP,FF1
      COMPLEX*16    AARI,BBRI,UURI1,UURI2,UURI3,DDURM,CCON1,BB1,CCAB,
     1              CCOR,CCL,VVRIT,TTTI,DD,
     2              UU(5),DDISWT(900)
CCCCCC
	ZERO=CMPLX(0.0,0.0)
      TTTI=DCMPLX(0.0,1.0)
      NMX2=NMX+2
      NMXM2=NMX-2
      NONX=NMX-NMI+1
      TTMU=TMU
      EE=E
      HH=H
      YMG=1.E-18 
      YMP2=-.2*YMG 
      HSQ=HH**2 
      HSQ12=HSQ/12.0 
      RL=LP         
      UURI1=DCMPLX(0.0,0.0)
      IF(LP.EQ.1) UURI1=YMP2
      UURI2=YMG            
      R=0.0              
      NXMBAS=(LM-1)*NONX
      NXM=NXMBAS       
      DO 205 NX=1,NMX2
      R=R+HH
      VVRIT=VRIT(NX)         
      AARI=(EE-VVRIT)*TTMU+RL*(RL+1.)/(R**2)
      BBRI=1.-HSQ12*AARI                        
      UURI3=(2.+(HSQ*AARI/BBRI))*UURI2-UURI1      
      IF(NX.LT.NMI) GO TO 204               
      DDISWT(NX)=UURI2/(BBRI*R)               
  204 UURI1=UURI2                           
      UURI2=UURI3                          
      IF(NX.LT.NMXM2) GO TO 205         
      NX1=NX-NMXM2+1                   
      UU(NX1)=UURI1/BBRI                 
  205 CONTINUE                       
      DDURM=(1./(12.*HH))*(8.*(UU(4)-UU(2))-UU(5)+UU(1))
      CCON1=-DDURM/UU(3)                
      L=LP+1
      FFC=FC(L)
      FFCP=FCP(L)
      GGC=GC(L)
      GGCP=GCP(L)                         
      BB1=CCON1*FFC+FFCP         
      CCAB=-CCON1*DCMPLX(GGC,FFC)
      DD=-DCMPLX(GGCP,FFCP)    
      CCL=BB1/(CCAB+DD)
      CL(L)=CCL          
      SG=SIG(L)                
      FF1=CDABS(1.+2.*TTTI*CCL)
      F1(L)=FF1
      CCON1=FFC+CCL*(GGC+TTTI*FFC) 
      CCOR=CDEXP(SG*TTTI)*CCON1/UU(3)
      DO 210 NX=NMI,NMX
      DDISWT(NX)=DDISWT(NX)*CCOR
      IF (CDABS(DDISWT(NX)).LE.1.0E-30) DDISWT(NX)=ZERO
      DISWT(NX)=DDISWT(NX)
210   CONTINUE 
CCCCCC
      RETURN
      END


C======================================================================CCCN58490
      SUBROUTINE RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)                    CCN58500
C======================================================================CCCN58510
      DIMENSION FACLOG(500)                                             CCN58520
      DIMENSION LT(6)                                                   CCN58530
      K1=IA+IB-IE                                                       CCN58540
      K3=IC+ID-IE                                                       CCN58550
      K5=IA+IC-IFF                                                      CCN58560
      K7=IB+ID-IFF                                                      CCN58570
      K2=IE-IABS(IA-IB)                                                 CCN58580
      K4=IE-IABS(IC-ID)                                                 CCN58590
      K6=IFF-IABS(IA-IC)                                                CCN58600
      K8=IFF-IABS(IB-ID)                                                CCN58610
      K9=MIN0(K1,K2,K3,K4,K5,K6,K7,K8)                                  CCN58620
      RAC=0.0                                                           CCN58630
      IF (K9) 4000,20,20                                                CCN58640
   20 K2=K1-2*(K1/2)                                                    CCN58650
      K4=K3-2*(K3/2)                                                    CCN58660
      K6=K5-2*(K5/2)                                                    CCN58670
      K8=K7-2*(K7/2)                                                    CCN58680
      IF(MAX0(K2,K4,K6,K8)) 4000,25,4000                                CCN58690
   25 LTMIN=MIN0(IA,IB,IC,ID,IE,IFF)                                    CCN58700
      IF(LTMIN) 4000,30,150                                             CCN58710
   30 LT(1)=IA                                                          CCN58720
      LT(2)=IB                                                          CCN58730
      LT(3)=IC                                                          CCN58740
      LT(4)=ID                                                          CCN58750
      LT(5)=IE                                                          CCN58760
      LT(6)=IFF                                                         CCN58770
      LTMIN=LT(1)                                                       CCN58780
      KMIN=1                                                            CCN58790
      DO 40 N=2,6                                                       CCN58800
      IF(LT(N)-LTMIN)35,40,40                                           CCN58810
   35 LTMIN=LT(N)                                                       CCN58820
      KMIN=N                                                            CCN58830
   40 CONTINUE                                                          CCN58840
      S1=1.0                                                            CCN58850
      F1=IE                                                             CCN58860
      F2=IFF                                                            CCN58870
      GO TO (55,55,55,55,45,50),KMIN                                    CCN58880
   45 F1=IA                                                             CCN58890
      F2=IC                                                             CCN58900
      S1=1-2*MOD(K5/2,2)                                                CCN58910
      GO TO 55                                                          CCN58920
   50 F1=IA                                                             CCN58930
      F2=IB                                                             CCN58940
      S1=1-2*MOD(K1/2,2)                                                CCN58950
   55 RAC=S1/ SQRT((F1+1.  )*(F2+1.  ))                                 CCN58960
      GO TO 4000                                                        CCN58970
  150 IABEP=(IA+IB+IE)/2+1                                              CCN58980
      ICDEP=(IC+ID+IE)/2+1                                              CCN58990
      IACFP=(IA+IC+IFF)/2+1                                             CCN59000
      IBDFP=(IB+ID+IFF)/2+1                                             CCN59010
      IABE=IABEP-IE                                                     CCN59020
      IEAB=IABEP-IB                                                     CCN59030
      IBEA=IABEP-IA                                                     CCN59040
      ICDE=ICDEP-IE                                                     CCN59050
      IECD=ICDEP-ID                                                     CCN59060
      IDEC=ICDEP-IC                                                     CCN59070
      IACF=IACFP-IFF                                                    CCN59080
      IFAC=IACFP-IC                                                     CCN59090
      ICFA=IACFP-IA                                                     CCN59100
      IBDF=IBDFP-IFF                                                    CCN59110
      IFBD=IBDFP-ID                                                     CCN59120
      IDFB=IBDFP-IB                                                     CCN59130
      NZMAX=MIN0(IABE,ICDE,IACF,IBDF)                                   CCN59140
      IABCD1=(IA+IB+IC+ID+4)/2                                          CCN59150
      IEFMAD=(IE+IFF-IA-ID)/2                                           CCN59160
      IEFMBC=(IE+IFF-IB-IC)/2                                           CCN59170
      NZMI1=-IEFMAD                                                     CCN59180
      NZMI2=-IEFMBC                                                     CCN59190
      NZMIN=MAX0(0,NZMI1,NZMI2)+1                                       CCN59200
      SQLOG=0.5  *(FACLOG(IABE)+FACLOG(IEAB)+FACLOG(IBEA)+FACLOG(ICDE)  CCN59210
     1          +FACLOG(IECD)+FACLOG(IDEC)+FACLOG(IACF)+FACLOG(IFAC)    CCN59220
     2          +FACLOG(ICFA)+FACLOG(IBDF)+FACLOG(IFBD)+FACLOG(IDFB)    CCN59230
     3 -FACLOG(IABEP+1)-FACLOG(ICDEP+1)-FACLOG(IACFP+1)-FACLOG(IBDFP+1))CCN59240
      DO 200 NZ=NZMIN,NZMAX                                             CCN59250
      NZM1=NZ-1                                                         CCN59260
      S1=1-2*MOD(NZM1,2)                                                CCN59270
      K1=IABCD1-NZM1                                                    CCN59280
      K2=IABE-NZM1                                                      CCN59290
      K3=ICDE-NZM1                                                      CCN59300
      K4=IACF-NZM1                                                      CCN59310
      K5=IBDF-NZM1                                                      CCN59320
      K6=NZ                                                             CCN59330
      K7=IEFMAD+NZ                                                      CCN59340
      K8=IEFMBC+NZ                                                      CCN59350
      SSLOG=SQLOG+FACLOG(K1)-FACLOG(K2)-FACLOG(K3)-FACLOG(K4)           CCN59360
     1           -FACLOG(K5)-FACLOG(K6)-FACLOG(K7)-FACLOG(K8)           CCN59370
      SSTERM=S1* EXP(SSLOG)                                             CCN59380
      RAC=RAC+SSTERM                                                    CCN59390
  200 CONTINUE                                                          CCN59400
      IF( ABS(RAC).LT.1.0E-10 ) RAC=0.0                                 CCN59410
 4000 RETURN                                                            CCN59420
      END                                                               CCN59430

                                                                        CCN59440
C======================================================================CCCN59450
      SUBROUTINE CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)                    CCN59460
C======================================================================CCCN59470
      DIMENSION FACLOG(500)                                             CCN59480
      RAC=0.0                                                           CCN59490
      IF(ID+IE-IFF)   1000,105,1000                                     CCN59500
  105 K1=IA+IB+IC                                                       CCN59510
      IF(K1-2*(K1/2)) 1000,110,1000                                     CCN59520
  110 K1=IA+IB-IC                                                       CCN59530
      K2=IC-IABS(IB-IA)                                                 CCN59540
      K3=MIN0(K1,K2)                                                    CCN59550
      IF(K3) 1000,130,130                                               CCN59560
  130 IF((-1)**(IB+IE))  1000,1000,140                                  CCN59570
  140 IF((-1)**(IC+IFF)) 1000,1000,150                                  CCN59580
  150 IF(IA-IABS(ID))  1000,152,152                                     CCN59590
  152 IF(IB-IABS(IE))  1000,154,154                                     CCN59600
  154 IF(IC-IABS(IFF)) 1000,160,160                                     CCN59610
  160 IF(IA) 1000,175,165                                               CCN59620
  165 IF(IB) 1000,175,170                                               CCN59630
  170 IF(IC) 1000,180,250                                               CCN59640
  175 RAC=1.0                                                           CCN59650
      GO TO 1000                                                        CCN59660
  180 FB=IB+1                                                           CCN59670
      RAC=((-1.0  )**((IA-ID)/2))/ SQRT(FB)                             CCN59680
      GO TO 1000                                                        CCN59690
  250 FC2=IC+1                                                          CCN59700
      IABCP=(IA+IB+IC)/2+1                                              CCN59710
      IABC=IABCP-IC                                                     CCN59720
      ICAB=IABCP-IB                                                     CCN59730
      IBCA=IABCP-IA                                                     CCN59740
      IAPD=(IA+ID)/2+1                                                  CCN59750
      IAMD=IAPD-ID                                                      CCN59760
      IBPE=(IB+IE)/2+1                                                  CCN59770
      IBME=IBPE-IE                                                      CCN59780
      ICPF=(IC+IFF)/2+1                                                 CCN59790
      ICMF=ICPF-IFF                                                     CCN59800
      SQFCLG=0.5  *(LOG(FC2)-FACLOG(IABCP+1)                            CCN59810
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)                     CCN59820
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)                     CCN59830
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))                    CCN59840
      NZMIC2=(IB-IC-ID)/2                                               CCN59850
      NZMIC3=(IA-IC+IE)/2                                               CCN59860
      NZMI=MAX0(0,NZMIC2,NZMIC3)+1                                      CCN59870
      NZMX=MIN0(IABC,IAMD,IBPE)                                         CCN59880
      S1=1-2*MOD(NZMI-1,2)                                              CCN59890
      DO 400 NZ=NZMI,NZMX                                               CCN59900
      NZM1=NZ-1                                                         CCN59910
      NZT1=IABC-NZM1                                                    CCN59920
      NZT2=IAMD-NZM1                                                    CCN59930
      NZT3=IBPE-NZM1                                                    CCN59940
      NZT4=NZ-NZMIC2                                                    CCN59950
      NZT5=NZ-NZMIC3                                                    CCN59960
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)                CCN59970
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)                CCN59980
      SSTERM=S1*EXP(TERMLG)                                             CCN59990
      RAC=RAC+SSTERM                                                    CCN60000
  400 S1=-S1                                                            CCN60010
      IF(ABS(RAC).LT.1.0E-10) RAC=0.0                                   CCN60020
 1000 RETURN                                                            CCN60030
      END                                                               CCN60040

                                                                        CCN60050
C======================================================================CCCN60060
      SUBROUTINE CLEBZ(IA,IB,IC,FACLOG,RAC)                             CCN60070
C======================================================================CCCN60080
      DIMENSION FACLOG(500)                                             CCN60090
      RAC=0.0                                                           CCN60100
      IGTW=(IA+IB+IC)/2                                                 CCN60110
      IF(MOD(IGTW,2).NE.0) GO TO 1000                                   CCN60120
      K1=IA+IB-IC                                                       CCN60130
      K2=IC-IABS(IB-IA)                                                 CCN60140
      K3=MIN0(K1,K2)                                                    CCN60150
      IF(K3.LT.0)  GO TO 1000                                           CCN60160
      IG=IGTW/2                                                         CCN60170
      IAHF=IA/2                                                         CCN60180
      IBHF=IB/2                                                         CCN60190
      ICHF=IC/2                                                         CCN60200
      S1=1-2*MOD(IG+ICHF,2)                                             CCN60210
      IABC=IGTW+2                                                       CCN60220
      IABMC=IAHF+IBHF-ICHF+1                                            CCN60230
      ICAMB=ICHF+IAHF-IBHF+1                                            CCN60240
      IBCMA=IBHF+ICHF-IAHF+1                                            CCN60250
      IGMA=IG-IAHF+1                                                    CCN60260
      IGMB=IG-IBHF+1                                                    CCN60270
      IGMC=IG-ICHF+1                                                    CCN60280
      R1=0.5  *(FACLOG(IABMC)+FACLOG(ICAMB)+FACLOG(IBCMA)-FACLOG(IABC)) CCN60290
     1         +FACLOG(IG+1)-FACLOG(IGMA)-FACLOG(IGMB)-FACLOG(IGMC)     CCN60300
      H1=FLOAT(IC+1)                                                    CCN60310
      RAC=S1* SQRT( H1 )* EXP(R1)                                       CCN60320
 1000 RETURN                                                            CCN60330
      END                                                               CCN60340


C======================================================================CCCN57420
      SUBROUTINE NINEJ(L9,FACLOG,U9)                                    CCN57430
C======================================================================CCCN57440
      DIMENSION FACLOG(500),L9(10)                                      CCN57450
      DIMENSION LT(9)                                                   CCN57460
      U9=0.0                                                            CCN57470
      KEX=0                                                             CCN57480
      LT(1)=L9(1)                                                       CCN57490
      LMIN=LT(1)                                                        CCN57500
      KP=LMIN                                                           CCN57510
      IMIN=1                                                            CCN57520
      DO 20 I=2,9                                                       CCN57530
      LT(I)=L9(I)                                                       CCN57540
      KP=KP+L9(I)                                                       CCN57550
      IF(LT(I)-LMIN) 15,15,20                                           CCN57560
   15 LMIN=LT(I)                                                        CCN57570
      IMIN=I                                                            CCN57580
   20 CONTINUE                                                          CCN57590
      IF(LMIN) 1000,70,50                                               CCN57600
   50 GO TO (110,300,300,300,300,150,300,170,300),IMIN                  CCN57610
   70 GO TO (110,110,110,110,150,150,170,170,190),IMIN                  CCN57620
  110 MM=(IMIN-1)/2+1                                                   CCN57630
      M1=MM+MM-1                                                        CCN57640
      M2=M1+1                                                           CCN57650
      M3=MM+4                                                           CCN57660
      L1=LT(7)                                                          CCN57670
      LT(7)=LT(M1)                                                      CCN57680
      LT(M1)=L1                                                         CCN57690
      L1=LT(8)                                                          CCN57700
      LT(8)=LT(M2)                                                      CCN57710
      LT(M2)=L1                                                         CCN57720
      L1=LT(9)                                                          CCN57730
      LT(9)=LT(M3)                                                      CCN57740
      LT(M3)=L1                                                         CCN57750
      IMIN=IMIN+(7-M1)                                                  CCN57760
      GO TO 175                                                         CCN57770
  150 KEX=1                                                             CCN57780
      M1=7                                                              CCN57790
      M2=8                                                              CCN57800
      M3=IMIN+IMIN-9                                                    CCN57810
      M4=M3+1                                                           CCN57820
      GO TO 180                                                         CCN57830
  170 KEX=1                                                             CCN57840
  175 M1=5                                                              CCN57850
      M2=6                                                              CCN57860
      M3=IMIN-6                                                         CCN57870
      M4=M3+2                                                           CCN57880
  180 L1=LT(M1)                                                         CCN57890
      LT(M1)=LT(M3)                                                     CCN57900
      LT(M3)=L1                                                         CCN57910
      L1=LT(M2)                                                         CCN57920
      LT(M2)=LT(M4)                                                     CCN57930
      LT(M4)=L1                                                         CCN57940
      L1=LT(9)                                                          CCN57950
      LT(9)=LT(IMIN)                                                    CCN57960
      LT(IMIN)=L1                                                       CCN57970
  190 IF(LT(9))1000,200,300                                             CCN57980
  200 IF(LT(5)-LT(6)) 1000,210,1000                                     CCN57990
  210 IF(LT(7)-LT(8)) 1000,220,1000                                     CCN58000
  220 IA=LT(1)                                                          CCN58010
      IB=LT(2)                                                          CCN58020
      IC=LT(3)                                                          CCN58030
      ID=LT(4)                                                          CCN58040
      IE=LT(5)                                                          CCN58050
      IFF=LT(7)                                                         CCN58060
      R1=(IE+1)*(IFF+1)                                                 CCN58070
      R1= SQRT(R1)                                                      CCN58080
      S1=(-1.0  )**((IE+IFF-IA-ID)/2)                                   CCN58090
      CALL RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)                          CCN58100
      U9=RAC*S1/R1                                                      CCN58110
      GO TO 370                                                         CCN58120
  300 K1=IABS (LT(2)-LT(7))                                             CCN58130
      K2=IABS (LT(3)-LT(5))                                             CCN58140
      K3=IABS (LT(4)-LT(9))                                             CCN58150
      MINRDA=MAX0  (K1,K2,K3)                                           CCN58160
      K1=LT(2)+LT(7)                                                    CCN58170
      K2=LT(3)+LT(5)                                                    CCN58180
      K3=LT(4)+LT(9)                                                    CCN58190
      MAXRDA=MIN0  (K1,K2,K3)                                           CCN58200
      IF(MINRDA-MAXRDA) 320,320,1000                                    CCN58210
  320 DO 350 N1=MINRDA,MAXRDA,2                                         CCN58220
      RAMDA2=N1                                                         CCN58230
      IA=LT(2)                                                          CCN58240
      IB=LT(5)                                                          CCN58250
      IC=LT(7)                                                          CCN58260
      ID=LT(3)                                                          CCN58270
      IE=LT(1)                                                          CCN58280
      IFF=N1                                                            CCN58290
      CALL RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)                          CCN58300
      W1=(RAMDA2+1.0  )*RAC                                             CCN58310
      IB=LT(4)                                                          CCN58320
      ID=LT(9)                                                          CCN58330
      IE=LT(8)                                                          CCN58340
      CALL RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)                          CCN58350
      W1=W1*RAC                                                         CCN58360
      IA=LT(3)                                                          CCN58370
      IC=LT(5)                                                          CCN58380
      IE=LT(6)                                                          CCN58390
      CALL RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)                          CCN58400
      W1=W1*RAC                                                         CCN58410
  350 U9=U9+W1                                                          CCN58420
      IF( ABS(U9).LT.1.0E-10) U9=0.0                                    CCN58430
  370 IF(KEX) 400,1000,400                                              CCN58440
  400 U9=U9*((-1.0  )**(KP/2))                                          CCN58450
 1000 RETURN                                                            CCN58460
      END                                                               CCN58470


C======================================================================CCCN47310
      SUBROUTINE GAUSF(NGAUS,RTS,WGT)                                   CCN47320
C======================================================================CCCN47330
      DIMENSION     RTS(96),WGT(96)
      DIMENSION     FLDR(500)                                           CCN47360
      EPS=1.0E-12                                                       CCN47370
      PAI=4.0*ATAN(1.)                                                  CCN47380
      FN=NGAUS                                                          CCN47390
      DS=PAI/(FN+0.5)                                                   CCN47400
      S1=-0.25*DS                                                       CCN47410
      L1MAX=NGAUS+1                                                     CCN47420
      NGAUS2=NGAUS-NGAUS/2                                              CCN47430
      FLDR(1)=1.0                                                       CCN47440
      DO 50 K=1,NGAUS2                                                  CCN47450
      S1=S1+DS                                                          CCN47460
      X0=COS(S1)                                                        CCN47470
      FLDR(2)=X0                                                        CCN47480
      DO 45 I=1,999                                                     CCN47490
      IF(L1MAX.LE.2) GO TO 42                                           CCN47500
      ALDR0=1.0                                                         CCN47510
      ALDR1=X0                                                          CCN47520
      DO 41 L1=3,L1MAX                                                  CCN47530
      L=L1-1                                                            CCN47540
      FL=L                                                              CCN47550
      ALDR2=((2.0*FL-1.0)*X0*ALDR1-(FL-1.0)*ALDR0)/FL                   CCN47560
      FLDR(L1)=ALDR2                                                    CCN47570
      ALDR0=ALDR1                                                       CCN47580
      ALDR1=ALDR2                                                       CCN47590
   41 CONTINUE                                                          CCN47600
   42 PNM1=FLDR(NGAUS)                                                  CCN47610
      PN=FLDR(L1MAX)                                                    CCN47620
      PND=FN*(PNM1-X0*PN)/(1.0-X0*X0)                                   CCN47630
      X1=X0-PN/PND                                                      CCN47640
      IF(ABS(X0-X1).LT.EPS) GO TO 48                                    CCN47650
      X0=X1                                                             CCN47660
   45 CONTINUE                                                          CCN47670
   48 WGHT=2.0*(1.0-X1*X1)/(FN*FN*PNM1*PNM1)                            CCN47680
      KMG=NGAUS-K+1                                                     CCN47690
      RTS(KMG)=X1                                                       CCN47700
      RTS(K)=-X1                                                        CCN47710
      WGT(KMG)=WGHT                                                     CCN47720
      WGT(K)=WGHT                                                       CCN47730
   50 CONTINUE                                                          CCN47740
      RETURN                                                            CCN47750
      END                                                               CCN47760


C==========================================================================
      SUBROUTINE BESSEL(X,LP1MX,AK,AI)                                
C==========================================================================
      DIMENSION AK(500),AI(500)
      DIMENSION F(51),G(51)                                             CCN45060
                                                                        CCN45070
C     THE FOLLOWING IS FOR SPHERICAL BESSEL AND NEUMANN                 CCN45080
      AK(1)=SIN(X)/X                                                    CCN45090
      AK(2)=(AK(1)-COS(X))/X                                            CCN45100
      AI(1)=-COS(X)/X                                                   CCN45110
      AI(2)=(AI(1)-SIN(X))/X                                            CCN45120
      DO 100 I=3,LP1MX                                                  CCN45130
      L=I-1                                                             CCN45140
      T1=(L*2-1)                                                        CCN45150
      AK(I)=AK(I-1)*T1/X-AK(I-2)                                        CCN45160
      AI(I)=AI(I-1)*T1/X-AI(I-2)                                        CCN45170
  100 CONTINUE                                                          CCN45180
 7810 FORMAT(1X,F9.4)                                                   CCN45210
 7800 FORMAT(1X,6E12.4)                                                 CCN45220
C     THE FOLLOWING IS ANOTHER VERSION CALCULATING SPHERICAL BESSEL     CCN45240
      R=X                                                               CCN45250
      F(2)=0.00                                                         CCN45260
      DO 2000 IL=1,LP1MX                                                CCN45270
      IF(IL.GE.18) GO TO 1900
      L=IL-1                                                            CCN45280
      LL=L+1                                                            CCN45290
      IF(LL-50)20,45,45                                                 CCN45300
   20 ELP=50.                                                           CCN45310
      J=50                                                              CCN45320
   45 DEL1=R                                                            CCN45330
      DEL=R-10.00                                                       CCN45340
      X=R                                                               CCN45350
      I=2                                                               CCN45360
      T2=2.*X                                                           CCN45370
      T3=X                                                              CCN45380
      SS=1.                                                             CCN45390
      TS=0.                                                             CCN45400
      SL=0.                                                             CCN45410
      TL=1.                                                             CCN45420
      SSS=1.                                                            CCN45430
      STS=0.                                                            CCN45440
      SSL=0.                                                            CCN45450
      STL=TL                                                            CCN45460
      EN=0.                                                             CCN45470
      DO 620 K=1,25                                                     CCN45480
      T5=EN+1.                                                          CCN45490
      T6=T5+EN                                                          CCN45500
      T7=EN*T5                                                          CCN45510
      T9=-T7/(T2*T5)                                                    CCN45520
      T5=-T9*TS                                                         CCN45530
      TS=T9*SS                                                          CCN45540
      SS=T5                                                             CCN45550
      IF(ABS (SS/SSS)-1.0E-10)630,630,540                               CCN45560
  540 T5=-T9*TL-SS/X                                                    CCN45570
      TL=T9*SL-TS/X                                                     CCN45580
      SL=T5                                                             CCN45590
      SSS=SSS+SS                                                        CCN45600
      STS=STS+TS                                                        CCN45610
      SSL=SSL+SL                                                        CCN45620
      STL=STL+TL                                                        CCN45630
      EN=EN+1.                                                          CCN45640
  620 CONTINUE                                                          CCN45650
  630 T8=SIN (T3)                                                       CCN45660
      T9=COS (T3)                                                       CCN45670
      Y=SSS*T9-STS*T8                                                   CCN45680
      Z=SSL*T9-STL*T8                                                   CCN45690
      G(1)=Y                                                            CCN45700
      G(2)=(1./R*G(1)-Z)                                                CCN45710
      I=L+11                                                            CCN45720
      N=2.*R+11.                                                        CCN45730
      IF(I-N)960,960,950                                                CCN45740
  950 N=I                                                               CCN45750
  960 Y=1.0E-20                                                         CCN45760
      X=Y                                                               CCN45770
      Z=0.                                                              CCN45780
      T2=N                                                              CCN45790
 1000 T3=T2+1.                                                          CCN45800
      T4=T2+T3                                                          CCN45810
      Y =(T4/R*Y -Z )                                                   CCN45820
      IF(N-LL)1060,1060,1080                                            CCN45830
 1060 F(N)=Y                                                            CCN45840
      GO TO 1120                                                        CCN45850
 1080 IF(1.-ABS (Y))1090,1090,1120                                      CCN45860
 1090 CONTINUE                                                          CCN45870
      Y =Y *1.0E-20                                                     CCN45880
      X =X *1.0E-20                                                     CCN45890
 1120 N=N-1                                                             CCN45900
      Z=X                                                               CCN45910
      X=Y                                                               CCN45920
      T2=T2-1.                                                          CCN45930
      IF(N)1150,1150,1000                                               CCN45940
 1150 Y=F(1)*G(2)-F(2)*G(1)                                             CCN45950
      Z=1./Y                                                            CCN45960
      AK(IL)=F(LL)*Z/R                                                  CCN45970
      IF(L.EQ.0) AK(1)=SIN(R)/R                                         CCN45980
      GO TO 2000
 1900 IF(R.LT.5) THEN
        T1=IL-1
        IF(IL.GT.1) THEN
          T5=1.0
          DO 1910 N1=1,IL
          T4=FLOAT(IL+N1-1)
          T5=T5/T4 
 1910     CONTINUE         
          T3=ALOG(T5)
        ELSE
          T3=0.0
        ENDIF    
        T2=T1*LOG(R*2.0)+T3
        AK(IL)=EXP(T2)
      ELSE
        T1=(IL-1)*2-1
        AK(IL)=AK(IL-1)*T1/R-AK(IL-2)
      ENDIF
 2000 CONTINUE                                                          CCN45990
      RETURN                                                            CCN46010
      END                                                               CCN46020


      SUBROUTINE YLCAL(A,LLL,MMM,PKK)
      DIMENSION  PKK(130,30),PK(130,10)
      DOUBLE PRECISION A,EPS,PAI4,SQR4PI,ROOT,PK,Y,FK,TFK,FK1,
     1				 Z,TM1,TM2,TM3,FKP,FAC,T1,T2,H1
      EPS=1.0D-30
      PAI4=12.5663706
      SQR4PI=1.0/DSQRT(PAI4)
      ROOT=A
      LP1MAX=LLL
      MP1MAX=MMM
      PK(1,1)=SQR4PI
      PK(1,2)=0.0
      DO 620 M=3,5
      PK(1,M)=0.0
  620 PK(2,M)=0.0
      PK(2,1)=ROOT*SQR4PI
      Y=1.0-ROOT*ROOT
      IF(Y.LT.0.0) Y=0.0
      Y= DSQRT(Y)
      PK(2,2)=Y*SQR4PI
      IF(LP1MAX.LE.2) GO TO 660
      FK=1.0
      TFK=1.0
      DO 650 K=3,LP1MAX
      FK=FK+1.0
      FK1=FK-1.0
      TFK=TFK+2.0
      PK(K,1)=(TFK*ROOT*PK(K-1,1)-FK1*PK(K-2,1))/FK
      PK(K,2)=(TFK*ROOT*PK(K-1,2)-FK *PK(K-2,2))/FK1
CCCCC
      IF(DABS(FK).LT.EPS.OR.DABS(FK1).LT.EPS) 
     1        WRITE(6,998) ROOT,LP1MAX,MP1MAX
  998 FORMAT(1H   ,'DIVID CHECK FOR ROOT, LP1MAX, MP1MAX=',E10.2,2I5)
CCCCC
	if (y.lt.1.0d-7) go to 640
      Z=ROOT/Y
      TM1=0.0
      TM2=FK+1.0
      TM3=FK
      DO 635 MP=3,MP1MAX
      IF(MP.GT.K) GO TO 630
      TM1=TM1+2.0
      TM2=TM2-1.0
      TM3=TM3+1.0
      PK(K,MP)=TM1*Z*PK(K,MP-1)-TM2*TM3*PK(K,MP-2)
      GO TO 635
  630 PK(K,MP)=0.0
  635 CONTINUE
      GO TO 650
  640 DO 645 MP=3,MP1MAX
  645 PK(K,MP)=0.
  650 CONTINUE
  660 IF(LP1MAX.EQ.1) GO TO 1000 !RETURN
      DO 420 LP1=2,LP1MAX
      FK=LP1-1
      FKP=FK+1.0
      FAC=1.0
      H1=DSQRT(FK*2.0+1.0) 
      PK(LP1,1)=PK(LP1,1)*H1
      MP1X=MIN0(LP1,MP1MAX)
      DO 410 MP1=2,MP1X
      T1=-FAC/DSQRT(FK*FKP) 
      FAC=T1
      IF(DABS(FK).LT.EPS.OR.DABS(FKP).LT.EPS)
     1        WRITE(6,999) ROOT,LP1,MP1
  999 FORMAT(1H  ,'DIVID CHECK FOR ROOT, LP1,MP1=',E10.2,2I5)
      T2=PK(LP1,MP1)*FAC*H1
      PK(LP1,MP1)=T2
      FK=FK-1.0
      FKP=FKP+1.0
  410 CONTINUE
  420 CONTINUE
 1000 DO 800 LP1=1,LP1MAX
      DO 800 MP1=1,MP1MAX
	PKK(LP1,MP1)=PK(LP1,MP1)
  800 CONTINUE 
      RETURN
      END
                                                                        CCN55720
                                                                        CCN47770
C====================================================================== CCN47780
      SUBROUTINE DSPLS3(X,Y,N,XI,FI,M,Q,AU,IGO)                         CCN47790
C====================================================================== CCN47800
C                                                                       CCN47810
C     ******************************************************************CCN47820
C                                                                       CCN47830
C     CUBIC SPLINE, STARTING WITH ZERO SECOND DERIVATIVES AT THE        CCN47840
C     BOUNDARIES OF THE APPROXIMATION INTERVAL.                         CCN47850
C                                                                       CCN47860
C     IGO = 0      BUILD UP SPLINE ONLY.                                CCN47870
C     IGO = 1      BUILD UP SPLINE AND INTERPOLATE FOR M POINTS.        CCN47880
C                                                                       CCN47890
C     DOUBLE PRECISION VERSION.        J.GALONSKA, 3.11.1970            CCN47900
C                                                                       CCN47910
C     ******************************************************************CCN47920
C                                                                       CCN47930
      DIMENSION X(500), Y(500), XI(500), Q(500), AU(500), FI(500)       CCN47940
C                                                                       CCN47950
      SIX = 6.                                                          CCN47960
      FACT = 0.16666666666666667                                        CCN47970
      ZERO = 0.                                                         CCN47980
C                                                                       CCN47990
      Q(1) = ZERO                                                       CCN48000
      AU(1) = ZERO                                                      CCN48010
      AU(N) = ZERO                                                      CCN48020
      HK = X(2) - X(1)                                                  CCN48030
      YSAVE = (Y(2)-Y(1)) / HK                                          CCN48040
      AUX = ZERO                                                        CCN48050
      NN = N - 1                                                        CCN48060
      DO 10  K = 2,NN                                                   CCN48070
        HX = X(K+1) - X(K-1)                                            CCN48080
        DIVQ = (HK*Q(K-1)+HX+HX)                                        CCN48090
        HK = X(K+1) - X(K)                                              CCN48100
        YK = (Y(K+1)-Y(K)) / HK                                         CCN48110
        Q(K) = - HK / DIVQ                                              CCN48120
        AU(K) = (SIX*(YK-YSAVE)-AUX) / DIVQ                             CCN48130
        YSAVE = YK                                                      CCN48140
        AUX = AU(K) * HK                                                CCN48150
   10 CONTINUE                                                          CCN48160
C                                                                       CCN48170
      NN2 = NN + 2                                                      CCN48180
      DO 20  KK = 2,NN                                                  CCN48190
        K = NN2 - KK                                                    CCN48200
        AU(K) = Q(K) * AU(K+1) + AU(K)                                  CCN48210
   20 CONTINUE                                                          CCN48220
C                                                                       CCN48230
      IF (IGO.EQ.0) RETURN                                              CCN48240
C                                                                       CCN48250
C     ******************************************************************CCN48260
C                                                                       CCN48270
      ENTRY APPROX(XI,FI,M)                                             CCN48280
C                                                                       CCN48290
C     ******************************************************************CCN48300
C                                                                       CCN48310
C     FOR INTERPOLATION ONLY.                                           CCN48320
C                                                                       CCN48330
C     WRITE(*,*)X(1),X(10),XI(10),Y(10),M,N
      DO 100  J = 1,M                                                   CCN48340
C     WRITE(*,*)X(1),X(10),XI(10),Y(10),J,XI(J),N,X(N)
        IF(X(1).GT.XI(J)) THEN                                          CCN48350
          M1 = 1                                                        CCN48360
          M2 = 2                                                        CCN48370
        ELSE IF (XI(J).GT.X(N)) THEN                                    CCN48380
          M1 = N - 1                                                    CCN48390
          M2 = N                                                        CCN48400
        ELSE                                                            CCN48410
          M1 = 1                                                        CCN48420
          M2 = N                                                        CCN48430
          DO 50 IKI=1,20                                                CCN48440
            M3 = (M2+M1)/2                                              CCN48450
            IF (XI(J).LT.X(M3)) THEN                                    CCN48460
              M2 = M3                                                   CCN48470
            ELSE                                                        CCN48480
              M1 = M3                                                   CCN48490
            END IF                                                      CCN48500
            IF (M1+1.EQ.M2) GOTO 90                                     CCN48510
   50     CONTINUE                                                      CCN48520
   90     CONTINUE                                                      CCN48530
        END IF                                                          CCN48540
        DIJ = X(M2) - XI(J)                                             CCN48550
        DIJ3 = DIJ * DIJ * DIJ                                          CCN48560
        HI = X(M2) - X(M1)                                              CCN48570
        HI2 = HI * HI                                                   CCN48580
        DIM1J = X(M1) - XI(J)                                           CCN48590
        DIM1J3 = DIM1J * DIM1J * DIM1J                                  CCN48600
        FI(J) = FACT * (AU(M1)*DIJ3-AU(M2)*DIM1J3+(SIX*Y(M1)-HI2*AU(M1))CCN48610
     1               *DIJ-(SIX*Y(M2)-HI2*AU(M2))*DIM1J) / HI            CCN48620
  100 CONTINUE                                                          CCN48630
      RETURN                                                            CCN48640
      END                                                               CCN48650