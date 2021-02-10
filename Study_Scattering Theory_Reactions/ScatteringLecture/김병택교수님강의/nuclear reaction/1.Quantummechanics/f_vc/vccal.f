      PROGRAM MAIN 
	IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FACLOG(500), L9(10)

      open(unit=5, file='vccal.DAT',status='old')
      open(unit=6, file='vccal.out',status='unknown')

      FACLOG(1)=0.0D0
      FACLOG(2)=0.0D0
      FN=1.D0
      DO 5 N=3,500
      FN=FN+1.D0
    5 FACLOG(N)=FACLOG(N-1)+DLOG(FN)

      READ(5,10) KTRL, NCAL
	IF (KTRL.LE.2) WRITE(6,100) NCAL
	IF (KTRL.EQ.3) WRITE(6,101) NCAL
	IF (KTRL.EQ.4) WRITE(6,102) NCAL  

	IF (KTRL.LE.2) THEN
	   IF (KTRL.EQ.2) THEN
	   DO 50 N=1,NCAL
	   READ(5,11) AIA,AIB,AIC,AID,AIE,AIF 
	   IA = AIA*2.0D0+0.001D0
	   IB = AIB*2.0D0+0.001D0
	   IC = AIC*2.0D0+0.001D0
	   ID = AID*2.0D0+0.001D0
	   IF (AID.LT.0.0D0) ID= AID*2.0D0-0.001D0
	   IE = AIE*2.0D0+0.001D0
	   IF (AIE.LT.0.0) IE= AIE*2D0-0.001D0
	   IFF = AIF*2.0D0+0.001D0
	   IF (AIF.LT.0.0) IFF= AIF*2.0D0-0.001D0
	   CALL CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
	   WRITE(6,110) AIA,AID,AIB,AIE,AIC,AIF, RAC
50       CONTINUE
	   ELSE
	   DO 60 N=1,NCAL
	   READ(5,11) AIA,AIB,AIC
	   AID=0.0D0
	   AIE=0.0D0
	   AIF=0.0D0
	   IA = AIA*2.0D0+0.001D0
	   IB = AIB*2.0D0+0.001D0
	   IC = AIC*2.0D0+0.001D0
	   CALL CLEBZ(IA,IB,IC,FACLOG,RAC) 
	   WRITE(6,110) AIA,AID,AIB,AIE,AIC,AIF, RAC
60       CONTINUE
         ENDIF
      ELSE
	   IF (KTRL.EQ.3) THEN
	   DO 70 N=1,NCAL
	   READ(5,11) AIA,AIB,AIC,AID,AIE,AIF 
	   IA = AIA*2.0D0+0.001D0
	   IB = AIB*2.0D0+0.001D0
	   IC = AIC*2.0D0+0.001D0
	   ID = AID*2.0D0+0.001D0
	   IE = AIE*2.0D0+0.001D0
	   IFF = AIF*2.0D0+0.001D0
	   CALL RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC) 
	   WRITE(6,120) AIA,AIB,AIC,AID,AIE,AIF, RAC
70       CONTINUE
	   ELSE
	   DO 80 N=1,NCAL
	   READ(5,10) (L9(I),I=1,9)
	   CALL NINEJ(L9,FACLOG,U9) 
	   WRITE(6,130) (L9(I),I=1,9), U9
80       CONTINUE
         ENDIF
	ENDIF
10    FORMAT(10I5)
11    FORMAT(10F7.3)
100   FORMAT(//,1H ,'CLEBSCH-GORDAN COEFFICIENTS,   NCAL=',I3,/)
101   FORMAT(//,1H ,'RACAH COEFFICIENTS,   NCAL=',I3,/)
102   FORMAT(//,1H ,'9-J COEFFICIENTS,   NCAL=',I3,/)
110   FORMAT(1H ,'(',4F5.1,' | ',2F5.1,' ) =',F12.8)
120   FORMAT(1H ,'W(',4F5.1,' ; ',2F5.1,' ) =',F12.8)
130   FORMAT(1H ,'9-J FOR (',9I5,' ) =',F12.8)
      STOP
      END
CCCCCC
CCCCCC
C======================================================================CCCN58490
      SUBROUTINE RAC7(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)                    CCN58500
C======================================================================CCCN58510
	IMPLICIT REAL*8(A-H,O-Z)
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
      RAC=0.0D0                                                         CCN58630
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
      S1=1.0D0                                                          CCN58850
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
   55 RAC=S1/ DSQRT((F1+1.D0  )*(F2+1.D0  ))                            CCN58960
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
      SQLOG=0.5D0*(FACLOG(IABE)+FACLOG(IEAB)+FACLOG(IBEA)+FACLOG(ICDE)  CCN59210
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
      SSTERM=S1* DEXP(SSLOG)                                            CCN59380
      RAC=RAC+SSTERM                                                    CCN59390
  200 CONTINUE                                                          CCN59400
      IF( DABS(RAC).LT.1.0D-10 ) RAC=0.0D0                              CCN59410
 4000 RETURN                                                            CCN59420
      END                                                               CCN59430
CCCCCC
CCCCCC                                                                  CCN59440
C======================================================================CCCN59450
      SUBROUTINE CLEB(IA,IB,IC,ID,IE,IFF,FACLOG,RAC)                    CCN59460
C======================================================================CCCN59470
	IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FACLOG(500)                                             CCN59480
      RAC=0.0D0                                                         CCN59490
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
  175 RAC=1.0D0                                                         CCN59650
      GO TO 1000                                                        CCN59660
  180 FB=IB+1                                                           CCN59670
      RAC=((-1.0D0  )**((IA-ID)/2))/ DSQRT(FB)                          CCN59680
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
      SQFCLG=0.5D0  *(LOG(FC2)-FACLOG(IABCP+1)                          CCN59810
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
      SSTERM=S1*DEXP(TERMLG)                                            CCN59990
      RAC=RAC+SSTERM                                                    CCN60000
  400 S1=-S1                                                            CCN60010
      IF(DABS(RAC).LT.1.0D-10) RAC=0.0D0                                CCN60020
 1000 RETURN                                                            CCN60030
      END                                                               CCN60040
CCCCCC
CCCCCC                                                                  CCN60050
C======================================================================CCCN60060
      SUBROUTINE CLEBZ(IA,IB,IC,FACLOG,RAC)                             CCN60070
C======================================================================CCCN60080
	IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FACLOG(500)                                             CCN60090
      RAC=0.0D0                                                         CCN60100
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
      R1=0.5D0*(FACLOG(IABMC)+FACLOG(ICAMB)+FACLOG(IBCMA)-FACLOG(IABC)) CCN60290
     1         +FACLOG(IG+1)-FACLOG(IGMA)-FACLOG(IGMB)-FACLOG(IGMC)     CCN60300
      H1=FLOAT(IC+1)                                                    CCN60310
      RAC=S1* DSQRT( H1 )* DEXP(R1)                                     CCN60320
 1000 RETURN                                                            CCN60330
      END                                                               CCN60340
CCCCCC
CCCCCC
C======================================================================CCCN57420
      SUBROUTINE NINEJ(L9,FACLOG,U9)                                    CCN57430
C======================================================================CCCN57440
	IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FACLOG(500),L9(10)                                      CCN57450
      DIMENSION LT(9)                                                   CCN57460
      U9=0.0D0                                                          CCN57470
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
      R1= DSQRT(R1)                                                     CCN58080
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
      W1=(RAMDA2+1.0D0  )*RAC                                           CCN58310
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
      IF( DABS(U9).LT.1.0D-10) U9=0.0D0                                 CCN58430
  370 IF(KEX) 400,1000,400                                              CCN58440
  400 U9=U9*((-1.0D0  )**(KP/2))                                        CCN58450
 1000 RETURN                                                            CCN58460
      END                                                               CCN58470
