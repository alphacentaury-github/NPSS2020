      PROGRAM MAIN 
	IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16       EXSGRI                                              
      REAL          K,K1,K2,K3,K4,M1,M2,M3,M4                           
      COMMON/COUCC/ FC(130),GC(130),FCP(130),GCP(130),EXSGRI(130),      
     1              MAXL,N1,EETA,SIGMAX,RD,RHO,KTOUT7,K9,STEP,ACCUR     

      open(unit=5, file='flglcal.dat',status='old')
      open(unit=6, file='flglcal.out',status='unknown')

      PI  =4.0*DATAN(1.D0) 
      HBAR=197.327053
      AMAS=931.49432
      WNUNIT=DSQRT(2.*AMAS)/HBAR
      FINE=1.0/137.0359896

      READ(5,10) TMAS,ZZT,PMAS,ZZP,ELAB
      READ(5,10) XMIN,XSTEP
      READ(5,11) NXMAX,MAXL
10    FORMAT(10F7.3)
11    FORMAT(10I5)  

      TM=TMAS*AMAS                                                      
      PM=PMAS*AMAS                                                                                                                 
      C=PM*TM/((PM+TM)*AMAS)                                                   
      ECM=ELAB*TM/(PM+TM)
      WN=WNUNIT*DSQRT(C*ECM)
	EETA=C*ZZT*ZZP*HBAR*FINE/(WN*HBAR*HBAR)                 
      STEP=0.0                                                          
      ACCUR=0.0 
	RHOT=EETA+DSQRT(EETA**2+DFLOAT(MAXL*(MAXL-1))) 
	      
      WRITE(6,15) TMAS,ZZT,PMAS,ZZP,C
	WRITE(6,16) ELAB,ECM,WN,EETA,RHOT    
15    FORMAT(1H ,'COULOMB WAVE FUNCTIONS FOR'/
     15X,'TMAS=',F7.3,' ZZT=',F7.3,' PMAS=',F7.3,' ZZP=',F7.3,' RMAS='
     2F7.3)
16    FORMAT(5X,'ELAB=',F7.3,' ECM=',F7.3,' WN  =',F7.3,' ETA=',F7.3,
     2 ' TURN=',F7.3)                                            
      LL=MAXL-1
      WRITE(6,20) LL
20    FORMAT(//1H ,'REGULAR COULOMB WAVE FUNCTIONS UP TO L=', I3, /
     1 1H ,'   X         L=0            L=1           L=2',/) 

      DO 100 N=1,NXMAX
	XX =XMIN+ DFLOAT(N-1)*XSTEP
	X=XX
	RHO=X*WN
	IF (RHO.LT.RHOT*0.2) GO TO 100
	CALL FLGLCH
      WRITE(6,21) XX,(FC(L),L=1,MAXL)                                     
100   CONTINUE
21    FORMAT(F7.3,10E14.6)

      WRITE(6,22) LL
22    FORMAT(//,1H ,'IRREGULAR COULOMB WAVE FUNCTIONS UP TO L=', I3, /
     1 1H ,'   X         L=0            L=1           L=2',/) 
      DO 200 N=1,NXMAX
	XX =XMIN+ DFLOAT(N-1)*XSTEP
	X=XX
	RHO=X*WN
	IF (RHO.LT.RHOT*0.2) GO TO 200
	CALL FLGLCH
      WRITE(6,21) XX,(GC(L),L=1,MAXL)                                     
200   CONTINUE
      STOP
      END
CCCCCC
CCCCCC
C======================================================================CCCN50460
      SUBROUTINE FLGLCH                                                 CCN50470
C======================================================================CCCN50480
	IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16       EXSGRI                                           CCN50490
      REAL          K,K1,K2,K3,K4,M1,M2,M3,M4                           CCN50500
      COMMON/COUCC/ FC(130),GC(130),FCP(130),GCP(130),EXSGRI(130),      CCN50520
     1              MAXL,N1,EETA,SIGMAX,RD,RHO,KTOUT7,K9,STEP,ACCUR     CCN50530
      ETA=EETA                                                          CCN50540
      IF(ETA-10.) 280,280,290                                           CCN50580
280   ETA2=ETA**2                                                       CCN50590
      ETA2A=2.*ETA                                                      CCN50600
      ETA6=ETA2+16.                                                     CCN50610
      SIGMAO=-(ETA/(12.*ETA6))*(1.+(ETA2-48.)/(30.*(ETA6**2))+((ETA2-   CCN50620
     1 160.)*ETA2+1280.)/(105.*(ETA6**4)))-ETA+(ETA/2.)*DLOG(ETA6)+3.5* CCN50630
     2 DATAN(.25*ETA)-(DATAN(ETA)+DATAN(.5*ETA)+DATAN(ETA/3.))          CCN50640
      GO TO 300                                                         CCN50650
290   EINV1=1./ETA                                                      CCN50660
      EINV2=EINV1*EINV1                                                 CCN50670
      EINV3=EINV1*EINV2                                                 CCN50680
      EINV5=EINV3*EINV2                                                 CCN50690
      EINV7=EINV5*EINV2                                                 CCN50700
      EINV9=EINV7*EINV2                                                 CCN50710
      SIGMAO=.7853981634+ETA*DLOG(ETA)-ETA-(.08333333333*EINV1+         CCN50720
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
      IF(ACC.LT.1.0D-15.OR.ACC.GT.1.0D-6) ACC = 1.0E-6                  CCN50850
      R=RHO                                                             CCN50860
      KTR=1                                                             CCN50870
      LMAX=MAXL                                                         CCN50880
      LMIN1=MINL+1                                                      CCN50890
      XLL1=DFLOAT(MINL*LMIN1)                                           CCN50900
      ETA2=ETA*ETA                                                      CCN50910
      TURN=ETA+DSQRT(ETA2+XLL1)                                         CCN50920
      IF(R.LT.TURN.AND.DABS(ETA).GE.1.0D-6)KTR=-1                       CCN50930
      KTRP=KTR                                                          CCN50940
      GOTO2                                                             CCN50950
1     R=TURN                                                            CCN50960
      TF=F                                                              CCN50970
      TFP=FP                                                            CCN50980
      LMAX=MINL                                                         CCN50990
      KTRP=1                                                            CCN51000
2     ETAR=ETA*R                                                        CCN51010
      RHO2=R*R                                                          CCN51020
      PL=DFLOAT(LMAX+1)                                                 CCN51030
      PMX=PL+0.5                                                        CCN51040
      FP=ETA/PL+PL/R                                                    CCN51050
      DK=ETAR*2.00                                                      CCN51060
      DEL=0.0                                                           CCN51070
      D=0.0                                                             CCN51080
      F=1.0                                                             CCN51090
      K=(PL*PL-PL+ETAR)*(2.0*PL-1.0)                                    CCN51100
      ABPL=DABS(PL**2+PL+ETAR)                                          CCN51110
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
      IF(PL.GT.20000.D0)GOTO11                                          CCN51230
      IF(DABS(DEL/FP).GE.ACC)GOTO3                                      CCN51240
      FP=F*FP                                                           CCN51250
      IF(LMAX.EQ.MINL)GOTO5                                             CCN51260
      FC(LMAX+1)=F                                                      CCN51270
      FCP(LMAX+1)=FP                                                    CCN51280
      L=LMAX                                                            CCN51290
      DO 4 LP=LMIN1,LMAX                                                CCN51300
      PL = DFLOAT(L)                                                    CCN51310
      GC(L+1)=ETA/PL+PL/R                                               CCN51320
      GCP(L+1)=DSQRT(ETA2+PL*PL)/PL                                     CCN51330
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
      IF(DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q))*ACC)GOTO6               CCN51690
      P=P/R                                                             CCN51700
      Q=Q/R                                                             CCN51710
      G=(FP-P*F)/Q                                                      CCN51720
      GP=P*G-Q*F                                                        CCN51730
      W=1.0/DSQRT(FP*G-F*GP)                                            CCN51740
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
      WRITE(6,560)                                                      CCN52330
      WRITE(6,570) (FC(L),L=1,L4)                                       CCN52340
      WRITE(6,580) (FCP(L),L=1,L4)                                      CCN52350
      WRITE(6,590) (GC(L),L=1,L4)                                       CCN52360
      WRITE(6,600) (GCP(L),L=1,L4)                                      CCN52370
560   FORMAT(//,5X,43H THE COULOMB WAVE FUNCTIONS FOR L=0 TO LMAX,2(/)) CCN52380
570   FORMAT(6HFC(L)=,5E15.5)                                           CCN52390
580   FORMAT(7HFCP(L)=,5E15.5)                                          CCN52400
590   FORMAT(6HGC(L)=,5E15.5)                                           CCN52410
600   FORMAT(7HGCP(L)=,5E15.5)                                          CCN52420
330   CONTINUE                                                          CCN52430
      RETURN                                                            CCN52440
      END                                                               CCN52450

                                                                        
