      OVERLAY(JPWKB,0,0)                                                        
      PROGRAM JP1CTRL(INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT,TAPE7,TAPE8)        
      COMMON /POTPAR/VSXD,WSXD,WSFD,DFND,DFNWD,DFNSD,RZED,RZEWD,RZESD,          
     1               VSXR,WSXC,WSFC,DFNC,DFNWC,DFNSC,RZEC,RZEWC,RZESC,          
     2               RZECD,DFNLD,DFNWLD,ALD,AWLD,                               
     3               RZECC,DFNLC,DFNWLC,ALC,AWLC                                
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC      ******      COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES   *****
      COMMON      EXTRA4(17090)                                         12/1/75 
      COMMON      AMS(3,3,3),BMTR(10,10,3),BMTI(10,10,3),NENSBP(3),             
     1            NENSBT(3),SGMEXP(100,6),POLEXP(100,2),FAI(4),NFAI,            
     2            NPOLST                                                        
C     COMMON     DUMMY(10000)                                                   
C     COMMON     DUMMY(7500)                                            12/1/75 
      COMMON     DUMY(20000)                                            10/31/77
C     CALL REINDEF(FACLOG(1),DUMMY(10000))                                      
      CALL OVERLAY(5HJPWKB,1,0,6HRECALL)                                        
      CALL SECOND(T1)                                                           
      T0=T1                                                                     
      CALL OVERLAY(5HJPWKB,2,0,6HRECALL)                                        
      CALL SECOND(T1)                                                           
      TIME=T1-T0                                                                
      WRITE (6,50) TIME                                                         
   50 FORMAT(///1H ,10H**********,22H TIME SPENT IN COUPLE=F10.5,               
     1              10H**********)                                              
      CALL OVERLAY(5HJPWKB,3,0,6HRECALL)                                        
      END                                                                       
      SUBROUTINE NLJJJK                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
      KTRL1T=KTRL(1)*(1-KTRL(7))  $  JJ=KEXCOM(45)  $  K=KEXCOM(46)             
      KT1PIS=KTRL1T+ISTRTW  $  KTISCK=KT1PIS-2*(KT1PIS/2)                       
      JJTRTW=2*JJ-2+KTISCK  $  NROW=0                                           
      KEXCOM(47)=JJTRTW                                                         
      DO 600 I1=1,IICPLE  $  IF(KTRL(7)) 210,230,210                            
  210 K1=K+1  $  J1TWMI=JJTRTW+2  $  J1TWMX=2*JJJMAX+2-KTISCK $GO TO 250        
  230 K1=K+KPRITR(I1)  $  IITRTW=2*IIREAD(I1)-KTRL(1)                           
      J1TWMI=IABS(IITRTW-JJTRTW)+2  $  J1TWMX=IITRTW+JJTRTW+2                   
  250 KPCHEK=K1-2*(K1/2)                                                        
      DO 600 J1TW=J1TWMI,J1TWMX,2 $ J1TRTW=J1TWMI+J1TWMX-J1TW-2                 
      L1TWMI=IABS(J1TRTW-ISTRTW)+2  $  L1TWMX=J1TRTW+ISTRTW+2                   
      DO 500 L1TW=L1TWMI,L1TWMX,2                                               
      L1TRTW=L1TWMI+L1TWMX-L1TW-2  $  L1TR=L1TRTW/2                             
      L1PRCK=L1TR-2*(L1TR/2)  $  IF(L1PRCK-KPCHEK) 500,410,500                  
  410 IF(L1TR-(KEXCOM(14)-1)) 420,420,500                                       
  420 NROW=NROW+1 $ NNROW(NROW)=I1 $ LLROW(NROW)=L1TRTW                         
      JJROW(NROW)=J1TRTW                                                        
  500 CONTINUE                                                                  
  600 CONTINUE                                                                  
      MXROW=NROW                                                                
      IF(NNROW(1)-1) 610,620,610                                                
  610 MXROW=0                                                                   
  620 RETURN                                                                    
      END                                                                       
      SUBROUTINE CLEB                                                           
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      RAC=0.0                                                                   
      IF(ID+IE-IF) 1000,105,1000                                                
  105 K1=IA+IB+IC  $  IF(K1-2*(K1/2)) 1000,110,1000                             
  110 K1=IA+IB-IC  $  K2=IC-IABS(IB-IA)  $  K3=MIN0(K1,K2)                      
      IF(K3) 1000,130,130                                                       
  130 IF((-1)**(IB+IE)) 1000,1000,140                                           
  140 IF((-1)**(IC+IF)) 1000,1000,150                                           
  150 IF(IA-IABS(ID)) 1000,152,152                                              
  152 IF(IB-IABS(IE)) 1000,154,154                                              
  154 IF(IC-IABS(IF)) 1000,160,160                                              
  160 IF(IA) 1000,175,165                                                       
  165 IF(IB) 1000,175,170                                                       
  170 IF(IC) 1000,180,250                                                       
  175 RAC=1.0  $  GO TO 1000                                                    
  180 FB=IB+1  $  RAC=((-1.0)**((IA-ID)/2))/SQRT(FB)  $  GO TO 1000             
  250 FC2=IC+1  $  IABCP=(IA+IB+IC)/2+1                                         
      IABC=IABCP-IC  $  ICAB=IABCP-IB  $  IBCA=IABCP-IA                         
      IAPD=(IA+ID)/2+1 $ IAMD=IAPD-ID $ IBPE=(IB+IE)/2+1 $ IBME=IBPE-IE         
      ICPF=(IC+IF)/2+1 $ ICMF=ICPF-IF                                           
      SQFCLG=0.5*(ALOG(FC2)-FACLOG(IABCP+1)                                     
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)                             
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)                             
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))                            
      NZMIC2 = (IB-IC-ID)/2  $  NZMIC3 = (IA-IC+IE)/2                           
      NZMI=MAX0(0,NZMIC2,NZMIC3)+1  $  NZMX=MIN0(IABC,IAMD,IBPE)                
      S1=(-1.0)**(NZMI-1)                                                       
      DO 400 NZ=NZMI,NZMX                                                       
      NZM1=NZ-1  $  NZT1=IABC-NZM1  $  NZT2=IAMD-NZM1                           
      NZT3=IBPE-NZM1  $  NZT4=NZ-NZMIC2  $  NZT5=NZ-NZMIC3                      
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)                        
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)                        
      SSTERM=S1*EXP(TERMLG)  $  RAC=RAC+SSTERM                                  
  400 S1=-S1                                                                    
 1000 RETURN                                                                    
      END                                                                       
      SUBROUTINE RAC7                                                           
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      DIMENSION      LT(6)                                              10/21/75
      K1=IA+IB-IE  $  K3=IC+ID-IE  $  K5=IA+IC-IF  $  K7=IB+ID-IF               
      K2=IE-IABS(IA-IB)  $  K4=IE-IABS(IC-ID)  $  K6=IF-IABS(IA-IC)             
      K8=IF-IABS(IB-ID)  $  K9=MIN0(K1,K2,K3,K4,K5,K6,K7,K8)                    
      RAC=0.0  $  IF (K9) 4000,20,20                                            
   20 K2=K1-2*(K1/2)  $  K4=K3-2*(K3/2)  $  K6=K5-2*(K5/2)                      
      K8=K7-2*(K7/2)  $  IF(MAX0(K2,K4,K6,K8)) 4000,25,4000                     
   25 LTMIN=MIN0(IA,IB,IC,ID,IE,IF)  $  IF(LTMIN) 4000,30,150                   
   30 LT(1)=IA $ LT(2)=IB $ LT(3)=IC $ LT(4)=ID $ LT(5)=IE $ LT(6)=IF           
      LTMIN=LT(1)  $  KMIN=1  $  DO 40 N=2,6  $  IF(LT(N)-LTMIN)35,40,40        
   35 LTMIN=LT(N)  $  KMIN=N                                                    
   40 CONTINUE                                                                  
      S1=1.0  $  F1=IE  $  F2=IF  $  GO TO (55,55,55,55,45,50),KMIN             
   45 F1=IA  $  F2=IC  $  S1=(-1.0)**(K5/2)  $  GO TO 55                        
   50 F1=IA  $  F2=IB  $  S1=(-1.0)**(K1/2)                                     
   55 RAC=S1/SQRT((F1+1.)*(F2+1.))  $  GO TO 4000                               
  150 IABEP=(IA+IB+IE)/2+1  $  ICDEP=(IC+ID+IE)/2+1                             
      IACFP=(IA+IC+IF)/2+1  $  IBDFP=(IB+ID+IF)/2+1                             
      IABE=IABEP-IE  $  IEAB=IABEP-IB  $  IBEA=IABEP-IA                         
      ICDE=ICDEP-IE  $  IECD=ICDEP-ID  $  IDEC=ICDEP-IC                         
      IACF=IACFP-IF  $  IFAC=IACFP-IC  $  ICFA=IACFP-IA                         
      IBDF=IBDFP-IF  $  IFBD=IBDFP-ID  $  IDFB=IBDFP-IB                         
      NZMAX=MIN0(IABE,ICDE,IACF,IBDF)  $  IABCD1=(IA+IB+IC+ID+4)/2              
      IEFMAD=(IE+IF-IA-ID)/2  $  IEFMBC=(IE+IF-IB-IC)/2                         
      NZMI1=-IEFMAD  $  NZMI2=-IEFMBC  $  NZMIN=MAX0(0,NZMI1,NZMI2)+1           
      SQLOG=0.5*(FACLOG(IABE)+FACLOG(IEAB)+FACLOG(IBEA)+FACLOG(ICDE)            
     1          +FACLOG(IECD)+FACLOG(IDEC)+FACLOG(IACF)+FACLOG(IFAC)            
     2          +FACLOG(ICFA)+FACLOG(IBDF)+FACLOG(IFBD)+FACLOG(IDFB)            
     3 -FACLOG(IABEP+1)-FACLOG(ICDEP+1)-FACLOG(IACFP+1)-FACLOG(IBDFP+1))        
      DO 200 NZ=NZMIN,NZMAX  $  NZM1=NZ-1                                       
      K1=IABCD1-NZM1  $  K2=IABE-NZM1  $  K3=ICDE-NZM1  $  K4=IACF-NZM1         
      K5=IBDF-NZM1  $  K6=NZ  $  K7=IEFMAD+NZ  $  K8=IEFMBC+NZ                  
      SSLOG=SQLOG+FACLOG(K1)-FACLOG(K2)-FACLOG(K3)-FACLOG(K4)                   
     1           -FACLOG(K5)-FACLOG(K6)-FACLOG(K7)-FACLOG(K8)                   
      SSTERM=((-1.0)**NZM1)*EXP(SSLOG)  $  RAC=RAC+SSTERM                       
  200 CONTINUE                                                                  
 4000 RETURN                                                                    
      END                                                                       
      SUBROUTINE FLGLCH(F,G,FP,GP,L)                                            
C                                                                               
CCCC   COULOMB FUNCTION BY BARNETT,FENG AND STEED                               
C                                                                               
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON       KXXX(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                    
      COMMON FC(142),GC(142),FCP(140),GCP(140),ETA,SIGMAZ,RHO,RD,MAXL           
      REAL K,K1,K2,K3,K4,M1,M2,M3,M4                                            
      MINL=L  $  MAXL=L                                                         
      ACCUR=1.0E-6                                                              
      ACCUR=1.0E-10                                                             
      STEP=100.0                                                                
      STEP=1000.0                                                               
      PACE=STEP                                                                 
      ACC=ACCUR                                                                 
      IF(PACE.LT.100.0)PACE=100.0                                               
      IF(ACC.LT.1.0E-15.OR.ACC.GT.1.0E-6) ACC = 1.0E-6                          
      R=RHO                                                                     
      KTR=1                                                                     
      LMAX=MAXL                                                                 
      LMIN1=MINL+1                                                              
      XLL1=FLOAT(MINL*LMIN1)                                                    
      ETA2=ETA*ETA                                                              
      TURN=ETA+SQRT(ETA2+XLL1)                                                  
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.0E-6)KTR=-1                                
      KTRP=KTR                                                                  
      GOTO2                                                                     
1     R=TURN                                                                    
      TF=F                                                                      
      TFP=FP                                                                    
      LMAX=MINL                                                                 
      KTRP=1                                                                    
2     ETAR=ETA*R                                                                
      RHO2=R*R                                                                  
      PL=FLOAT(LMAX+1)                                                          
      PMX=PL+0.5                                                                
      FP=ETA/PL+PL/R                                                            
      DK=ETAR*2.00                                                              
      DEL=0.0                                                                   
      D=0.0                                                                     
      F=1.0                                                                     
      K=(PL*PL-PL+ETAR)*(2.0*PL-1.0)                                            
      IF(PL*PL+PL+ETAR.NE.0.0)GOTO3                                             
      R=R+1.0E-6                                                                
      GOTO2                                                                     
3     H=(PL*PL+ETA2)*(1.0-PL*PL)*RHO2                                           
      K=K+DK+PL*PL*6.0                                                          
      D=1.0/(D*H+K)                                                             
      DEL=DEL*(D*K-1.0)                                                         
      IF(PL.LT.PMX)DEL=-R*(PL*PL+ETA2)*(PL+1.0)*D/PL                            
      PL=PL+1.0                                                                 
      FP=FP+DEL                                                                 
      IF(D.LT.0.0)F=-F                                                          
      IF(PL.GT.20000.)GOTO11                                                    
      IF(ABS(DEL/FP).GE.ACC)GOTO3                                               
      FP=F*FP                                                                   
      IF(LMAX.EQ.MINL)GOTO5                                                     
5     IF(KTRP.EQ.-1)GOTO1                                                       
      P=0.0                                                                     
      Q=R-ETA                                                                   
      PL=0.0                                                                    
      AR=-(ETA2+XLL1)                                                           
      AI=ETA                                                                    
      BR=2.0*Q                                                                  
      BI=2.0                                                                    
      WI=2.0*ETA                                                                
      DR=BR/(BR*BR+BI*BI)                                                       
      DI=-BI/(BR*BR+BI*BI)                                                      
      DP=-(AR*DI+AI*DR)                                                         
      DQ=(AR*DR-AI*DI)                                                          
6     P=P+DP                                                                    
      Q=Q+DQ                                                                    
      PL=PL+2.0                                                                 
      AR=AR+PL                                                                  
      AI=AI+WI                                                                  
      BI=BI+2.0                                                                 
      D=AR*DR-AI*DI+BR                                                          
      DI=AI*DR+AR*DI+BI                                                         
      T=1.0/(D*D+DI*DI)                                                         
      DR=T*D                                                                    
      DI=-T*DI                                                                  
      H=BR*DR-BI*DI-1.00                                                        
      K=BI*DR+BR*DI                                                             
      T=DP*H-DQ*K                                                               
      DQ=DP*K+DQ*H                                                              
      DP=T                                                                      
      IF(PL.GT.46000.)GOTO11                                                    
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC)GOTO6                           
      P=P/R                                                                     
      Q=Q/R                                                                     
      G=(FP-P*F)/Q                                                              
      GP=P*G-Q*F                                                                
      W=1.0/SQRT(FP*G-F*GP)                                                     
      G=W*G                                                                     
      GP=W*GP                                                                   
      IF(KTR.EQ.1)GOTO8                                                         
      F=TF                                                                      
      FP=TFP                                                                    
      LMAX=MAXL                                                                 
      IF(RHO.LT.0.2*TURN)PACE=999.00                                            
      R3=1.0/3.0D0                                                              
      H=(RHO-TURN)/(PACE+1.0)                                                   
      H2=0.5*H                                                                  
      I2=IFIX(PACE+0.001)                                                       
      ETAH=ETA*H                                                                
      H2LL=H2*XLL1                                                              
      S=(ETAH+H2LL/R)/R-H2                                                      
7     RH2=R+H2                                                                  
      T=(ETAH+H2LL/RH2)/RH2-H2                                                  
      K1=H2*GP                                                                  
      M1=S*G                                                                    
      K2=H2*(GP+M1)                                                             
      M2=T*(G+K1)                                                               
      K3=H*(GP+M2)                                                              
      M3=T*(G+K2)                                                               
      M3=M3+M3                                                                  
      K4=H2*(GP+M3)                                                             
      RH=R+H                                                                    
      S=(ETAH+H2LL/RH)/RH-H2                                                    
      M4=S*(G+K3)                                                               
      G=G+(K1+K2+K3+K2+K4)*R3                                                   
      GP=GP+(M1+M2+M2+M3+M4)*R3                                                 
      IF(G.GT.1.0E65.OR.GP.GT.1.0E65)GOTO11                                     
      R=RH                                                                      
      I2=I2-1                                                                   
      IF(I2.GE.0)GOTO7                                                          
      W=1.0/(FP*G-F*GP)                                                         
    8 F=W*F                                                             9/18/75 
      FP=W*FP                                                           9/18/75 
      IF (KTLOUT(5).EQ.1) WRITE(6,4000) L,F,FP,G,GP                     9/21/75 
 4000 FORMAT(1H,*L,F,FP,G,GP=*,I3,4E13.4)                               9/21/75 
      RETURN                                                                    
11     W=0.0                                                                    
      G=0.0                                                                     
      GP=0.0                                                                    
      GOTO8                                                                     
      END                                                                       
      SUBROUTINE FLGLNG                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      LMX=LMAX+2                                                                
      R=RHOMX  $  RI=1.0/R  $  EX=EXP(-R)                                       
      G(1)=EX*RI  $  G(2)=EX*(RI+RI*RI)                                         
      FL=1.0                                                                    
      DO 50 L=1,LMX                                                             
      TWELP1=2.0*FL+1.0                                                         
      G(L+2)=TWELP1*RI*G(L+1)+G(L)                                              
      GD(L)=-G(L+1)+RI*G(L)                                                     
      F(L)=FD(L)=0.0                                                            
   50 FL=FL+1.0                                                                 
      GD(1)=-G(2)                                                               
      IF(KTLOUT(5)) 4100,4190,4100                                              
 4100 WRITE(6,4110) R,RI,EX                                                     
 4110 FORMAT(19H0IN FLGLNG R,RI,EX=3E15.7)                                      
      WRITE(6,4120) (G (L),L=1,LMX)                                             
      WRITE(6,4130) (GD(L),L=1,LMX)                                             
 4120 FORMAT(4H G =,7E15.7)                                                     
 4130 FORMAT(4H GD=,7E15.7)                                                     
 4190 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      OVERLAY(JPWKB,1,0)                                                        
      PROGRAM JP1WKB                                                            
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON /POTPAR/VSXD,WSXD,WSFD,DFND,DFNWD,DFNSD,RZED,RZEWD,RZESD,          
     1               VSXR,WSXC,WSFC,DFNC,DFNWC,DFNSC,RZEC,RZEWC,RZESC,          
     2               RZECD,DFNLD,DFNWLD,ALD,AWLD,                               
     3               RZECC,DFNLC,DFNWLC,ALC,AWLC                                
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC      ******      COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES   *****
      COMMON      EXTRA1(10931),LAMD(5),EXTRA2(6154)                    12/1/75 
      COMMON      AMS(3,3,3),BMTR(10,10,3),BMTI(10,10,3),NENSBP(3),             
     1            NENSBT(3),SGMEXP(100,6),POLEXP(100,2),FAI(4),NFAI,            
     2            NPOLST                                                        
CCCCCC      ******      DIMENSION FIELDS USED ONLY IN THIS ROUTINE      ********
      DIMENSION   IIOT(6),IIBB(6),KPBB(6)                                       
      TYPE INTEGER ISBB,IIBB,KPBB,AMPMWR                                        
      EQUIVALENCE (AR,ISBB),(AR(2),IIBB),(AR(8),KPBB),(AR(14),IIOT)             
CCCCCC      IF KTRL(1)=0 TARGET IS OF EVEN A, IF=1 TARGET IS OF ODD A.          
CCCCCC      IF KTRL(2)=1 TMAS**(1/3)+PMAS**(1/3)                                
CCCCCC      IF KTRL(3)=1 LEE POTENTIAL IS USED                                  
CCCCCC      IF KTRL(4)=1 FOLDING POTENTIAL IS USED (USED IN POTENT)             
CCCCCC      IF KTRL(5)=1 DIAGONAL ELEMENT IN BMATRX (INTYPE=3) IS ADDED         
CCCCCC      IF KTRL(9)=1 RANGES OF JJJ AND K ARE RESTRICTED BY KEXCOM(5,6,7).   
CCCCCC      IF KTRL(9)=2 RANGES OF JJJ AND K ARE RESTRICTED BY KEXCOM(5,6).     
CCCCCC      IF KTRL(10)=1 DERIVATIVE WITH RESPECT TO SMALL R IS USED            
CCCCCC                    FOR RADIAL FORM FACTOR                                
CCCCCC      IF KTRL(10)=2 FOR GDR                                               
CCCCCC      IF KTRL(11)=1 EXPANSION OF V IN LEGENDRE POLINOMIAL IS MADE.        
CCCCCC      IF KTRL(12)=1 COMPLEX COUPLING IS TAKEN.                            
CCCCCC      IF KTRL(13)=1 COULOMB EXCITATION IS ADDED.                          
CCCCCC      IF KTRL(14)=1 WKB INTEGRATION IS MADE                               
CCCCCC      IF KTRL(15).NE.0 PADE APPROXIMATION IS USED                         
CCCCCC      IF KTRL(16).GT.0 INTERPOLATIONIS MADE                       10/16/75
CCCCCC      IF KTRL(27)=1 IT IS ASSUMED THAT C-MARICES WERE ALREADY COMPUTED.   
CCCCCC      KEXCOM(1)=NXIMIN                                            10/17/75
CCCCCC      KEXCOM(2)=NXMAX=NXCPLE IF NON-ZERO (CCCTRL,SN-181)                  
CCCCCC      KEXCOM(3)=NLCUT (LOWER CUT-OFF POINT)                               
CCCCCC      EXTCOM(1)=RWKBMX                                                    
CCCCCC      EXTCOM(5)=RFTPMI                                                    
CCCCCC      EXTCOM(6)=RFTPMX                                                    
CCCCCC      KEXCOM(4)=LLCUT (LOWER L-CUT OFF)                                   
CCCCCC      KEXCOM(5,6,7)=JMIN,JMAX,KPARITY USED WHEN KTRL(9) IS NON-ZERO.      
CCCCCC      KEXCOM(14)=LLMAX IN NLJJJK.                                         
CCCCCC      KEXCOM(15)=NUMBER OF MIN ANGLE OF TESTING PADE COVERGENCE           
CCCCCC      KEXCOM(16)=LLWKB                                                    
CCCCCC      KEXCOM(28)=KIREPT (BMATRX,AMATRX)                                   
CCCCCC      KEXCOM(41)=KBOUND  (C.F.COUPLE)                                     
CCCCCC      KEXCOM(42)=NBOUND  (C.F.COUPLE)                                     
CCCCCC      KEXCOM(43)=MXROWI  (C.F.COUPLE)                                     
CCCCCC      KEXCOM(45,46)=JJ,K IN NLJJJK.                                       
CCCCCC      KEXCOM(47)=JJTRTW IN NLJJJK                                         
CCCCCC      KEXCOM(48)=ISMTXZ IN (JJKLOOP,WKBCTR)                               
CCCCCC           NAMES OF SUBROUTINE WHERE KTLOUT(N) IS USED                    
CCCCCC      KTLOUT(3)    JKLOOP                                                 
CCCCCC      KTLOUT(5)    FLGLCH AND FLGLNG                                      
CCCCCC      KTLOUT(6)    CCCTRL                                                 
CCCCCC      KTLOUT(7)    POTENT                                                 
CCCCCC      KTLOUT(9)    AMATRIX                                                
CCCCCC      KTLOUT(10)   COUPLE                                                 
CCCCCC      KTLOUT(13)   SCATSG                                                 
CCCCCC      KTLOUT(14)   WKBCTR                                                 
CCCCCC      KTLOUT(15)   WKBINT                                                 
CCCCCC      KTLOUT(16)   CMATRX                                                 
CCCCCC      KTLOUT(17)   AMPDER                                                 
CCCCCC      KTLOUT(18)   RMATRX                                                 
CCCCCC      KTLOUT(20)   JP1CRS (PLOT)                                          
      TTR=(1.0,0.0)  $  TTI=(0.0,1.0)  $  ZERO=(0.0,0.0)                        
  701 FORMAT(14I5)                                                              
  702 FORMAT(10F7.2)                                                            
  703 FORMAT(3(2I5,F10.3))                                                      
  704 FORMAT(I5,4F7.2)                                                          
  705 FORMAT(7I10)                                                      11/21/75
  706 FORMAT(10F7.3,I2)                                                 5/1/76  
 4000 FORMAT(1H1)                                                               
 4100 FORMAT(1H0)                                                               
      FACLOG(1)=FACLOG(2)=0.0  $  FN=1.0                                        
      DO 10 N=3,1300                                                    12/1/75 
      FN=FN+1.0                                                                 
   10 FACLOG(N)=FACLOG(N-1)+ALOG(FN)                                            
  100 READ(5,701) (KTRL(N),N=1,28),(KEXCOM(N),N=1,14),(KTLOUT(N),N=1,28)CARD 1-5
      READ(5,702) (EXTCOM(N),N=1,10)                                            
      READ(5,701) IICPLE,INTYPE,INTMAX,ISTRTW,NANGLR,IIXCAL,IIXPLT,     CARD   7
     1            IIPCAL,IIPPLT,KANGRD                                          
      IIRMAX=IICPLE                                                     4/28/76 
      READ(5,703) (IIREAD(N),KPRITR(N),QVALUE(N),N=1,IIRMAX)            CARD   8
  105 READ(5,702) ELAB,PMAS,TMAS,CHARGE,XMES1,XMES2,AMUPMU              CARD   9
      IF(KANGRD) 112,111,112                                                    
  111 READ(5,702) (ANGLER(N),N=1,NANGLR)  $  GO TO 117                          
  112 READ(5,702) ANGMIN,ANGDEL  $  ANG=ANGMIN                          CARD  10
      DO 113 N=1,NANGLR  $  ANGLER(N)=ANG                                       
  113 ANG=ANG+ANGDEL                                                            
      IF(IIXPLT)115,117,115                                                     
  115 DO 116 I1=1,IIXPLT                                                        
  116 READ(5,702) (SGMEXP(N,I1),N=1,NANGLR)                                     
  117 DO 138 N=1,NANGLR                                                 4/28/76 
  138 THETA(N)=ANGLER(N)*0.0174532925                                           
CCCCCC                                                                          
CCCCCC      OUTPUT OF BASIC INPUT DATA                                          
CCCCCC                                                                          
  400 IF(KTRL(3).EQ.1) GO TO 410                                                
      READ(5,702) VSXD,WSXD,WSFD,DFND,DFNWD,DFNSD,RZED,RZEWD,RZESD,RZECD        
      READ(5,702) VSXC,WSXC,WSFC,DFNC,DFNWC,DFNSC,RZEC,RZEWC,RZESC,RZECC        
      READ(5,702) (WC(I),I=1,6)                                                 
      GO TO 450                                                                 
  410 READ(5,702) VSXD,WSXD,DFND,DFNLD,DFNWD,DFNWLD,                            
     1            RZED,RZEWD,RZECD                                              
      READ(5,702) VSXC,WSXC,DFNC,DFNLC,DFNWC,DFNWLC,                            
     1            RZEC,RZEWC,RZECC                                              
      READ(5,702) ALD,AWLD,ALC,AWLC,(WC(I),I=1,6)                               
      DFNSD=DFNSC=0.5                                                           
      RZESD=RZESC=RZED                                                          
      DFNW=1.0  $  DFNS=1.0                                                     
  450 RZERO=RZED  $  DFN=DFND                                                   
      VSX=VSXD                                                                  
      READ(5,706)(VCOUPL(I),I=1,10),KREAD1                                      
      IF (INTYPE.EQ.3) READ (5,702) (VCOUPL(I),I=11,20)                         
      READ(5,701) KPPADE,KTCONV,LINTPO                                  6/8/76  
      READ(5,702) RFTPMI,RFTPMX,RWKBMX                                  6/8/76  
  500 WRITE(6,501)                                                              
  501 FORMAT(1H1,//////////)                                                    
      WRITE(6,502)                                                              
  502 FORMAT(31X,    10H**********,6X,27HCOUPLED CHANNEL CALCULATION,5X,        
     110H**********)                                                            
      WRITE(6,503)                                                              
  503 FORMAT(1H0,48X,22H(ON PROGRAM JUPITOR-1))                                 
      WRITE(6,504) ELAB,CHARGE,TMAS,PMAS                                        
  504 FORMAT(////21X,5HELAB=F6.2,10H,  CHARGE=F5.1,8H,  TMAS=F7.3,8H,  P        
     1MAS=F6.3)                                                                 
      IF(KTRL(3).EQ.1) GO TO 510                                                
      WRITE(6,506)                                                              
     1       VSXD,WSXD,WSFD,DFND,DFNWD,DFNSD,RZED,RZEWD,RZESD,RZECD,            
     2       VSXC,WSXC,WSFC,DFNC,DFNWC,DFNSC,RZEC,RZEWC,RZESC,RZECC             
  506 FORMAT(///38X,3HVSX,5X,3HWSX,5X,3HWSF,5X,3HDFN,4X,4HDFNW,4X,4HDFNS        
     1          ,5X,3HRZE,4X4HRZEW,4X,4HRZES,4X,4HRZEC/                         
     2          21X,13HDIADONAL POT=,10F8.3/                                    
     3          21X,13HCOUPLING POT=,10F8.3)                                    
      WRITE (6,508) (WC(I),I=1,6)                                               
  508 FORMAT(21X,13HWC          =,6F8.4)                                        
      GO TO 514                                                                 
  510 WRITE (6,511)                                                             
     1      VSXD,WSXD,DFND,DFNLD,DFNWD,DFNWLD,ALD,AWLD,RZED,RZEWD,RZESD,        
     2      VSXC,WSXC,DFNC,DFNLC,DFNWC,DFNWLC,ALC,AWLC,RZEC,RZEWC,RZESC         
  511 FORMAT(///38X,3HVSX,5X,3HWSX,5X,3HDFN,4X,4HDFNL,4X,4HDFNW                 
     1       ,3X,5HDFNWL,6X,2HAL,5X,3HAWL,5X,3HRZE,4X,4HRZEW,4X,4HRZEC/         
     2          21X,13HDIADONAL POT=,11F8.3/                                    
     3          21X,13HCOUPLING POT=,11F8.3)                                    
      WRITE (6,508) (WC(I),I=1,6)                                               
  514 ISOE=ISTRTW-2*(ISTRTW/2)  $  IF(ISOE) 516,515,516                         
  515 ISTROT=ISTRTW/2  $  ISBB=2H    $  GO TO 517                               
  516 ISTROT=ISTRTW    $  ISBB=2H/2                                             
  517 WRITE(6,518) ISTROT,ISBB                                                  
  518 FORMAT(1H0,20X,16HPROJECTILE SPIN=I1,A2)                                  
      WRITE(6,520) IIRMAX,INTYPE,INTMAX,(VCOUPL(N),N=1,10)                      
  520 FORMAT(1H0,20X,7HIIRMAX=I1,8H,INTYPE=I1,8H,INTMAX=I1,6H,BETA=10F6.        
     13)                                                                        
      IF(INTYPE-3) 530,524,530                                          4/28/76 
  524 WRITE(6,527) (VCOUPL(N),N=11,20)                                          
  527 FORMAT(48X,5HBETA=10F6.3)                                                 
  530 DO 538 I1=1,IIRMAX  $  IF(KTRL(1)) 533,532,533                            
  532 IIOT(I1)=IIREAD(I1)      $  IIBB(I1)=2H    $  GO TO 534                   
  533 IIOT(I1)=IIREAD(I1)*2-1  $  IIBB(I1)=2H/2                                 
  534 IF(KPRITR(I1)-1) 536,535,536                                              
  535 KPBB(I1)=3H(+)  $  GO TO 538                                              
  536 KPBB(I1)=3H(-)                                                            
  538 CONTINUE                                                                  
      WRITE(6,540) (IIOT(I),IIBB(I),KPBB(I),QVALUE(I),I=1,IIRMAX)               
  540 FORMAT(1H0,20X,14HTARGET STATES,3(2X,I1,A2,A3,4H AT F6.3,5H MEV,))        
      WRITE(6,4100)  $  WRITE(6,562) (KTRL(N),N=1,28)                   4/28/76 
      WRITE(6,563) (KEXCOM(N),N=1,28) $ WRITE(6,564) (KTLOUT(N),N=1,28)         
      WRITE(6,565) (EXTCOM(N),N=1,10)                                           
      AMPMWR=3HAMU  $  IF(AMUPMU) 557,558,557                                   
  557 AMPMWR=3HPMU                                                              
  558 WRITE(6,566) XMES1,XMES2,AMPMWR                                           
  561 FORMAT(/21X,31HNO. OF STATES IN VARIOUS BANDS=7I5)                        
  562 FORMAT(///7H KTRL  ,28I4)                                                 
  563 FORMAT(7H KEXCOM,28I4)                                                    
  564 FORMAT(7H KTLOUT,28I4)                                                    
  565 FORMAT(7H EXTCOM,10F10.2)                                                 
  566 FORMAT(7H XMES1=F7.5,3X,7H XMES2=F7.5,3X,6H UNIT=A3)                      
      KTRL(15)=KPPADE                                                   4/28/76 
      KTRL(16)=LINTPO                                                   4/28/76 
      KEXCOM(15)=KTCONV                                                 4/28/76 
      EXTCOM(1)=RWKBMX                                                  4/28/76 
      EXTCOM(5)=RFTPMI                                                  4/28/76 
      EXTCOM(6)=RFTPMX                                                  4/28/76 
      WRITE(6,570) KTRL(15),KEXCOM(15),KTRL(16),                        6/8/76  
     1             EXTCOM(5),EXTCOM(6),EXTCOM(1)                        6/8/76  
  570 FORMAT (///,21X,25HKPPADE,KTCONV,LINTPO    =,3I5                  6/8/76  
     1          /,21X,25HRFTPMI,RFTPMX,RWKBMX    =,3F8.2)               6/8/76  
  690 WRITE(6,4000)                                                             
      CALL CCCTRL                                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CCCTRL                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(420,3)  $  TYPE COMPLEX EXSGRI                  9/18/75 
      LMAX=KEXCOM(14)  $  LMAXM1=LMAX-1                                         
      DX=XMES2                                                                  
      DO 130 I1=1,IIXCAL                                                        
  130 IIMULT(I1)=2*IIREAD(I1)+1-KTRL(1)                                         
      RMAS=PMAS*TMAS/(PMAS+TMAS)  $  ECM(1)=ELAB*RMAS/PMAS                      
      ETUNIT=0.15745400*(1.0-AMUPMU)+0.15805086*AMUPMU                          
      WNUNIT=0.21870660*(1.0-AMUPMU)+0.21953760*AMUPMU                          
      CFUNIT=1.0/(RMAS*WNUNIT*WNUNIT)                                           
      IIMAX=0                                                           10/18/75
      DO 150 I1=1,IICPLE                                                        
      IF(IIREAD(I1).GT.IIMAX)   IIMAX=IIREAD(I1)                        10/18/75
      ECM(I1)=ECM(1)-QVALUE(I1)  $  E1=ABS(ECM(I1))                             
      CE   (I1)=ETUNIT*CHARGE*SQRT(RMAS/E1)                                     
      WN   (I1)=WNUNIT*SQRT(RMAS*E1)                                            
      WNINI(I1)=WNUNIT*SQRT(RMAS*(ECM(I1)+VSX))                                 
      ETA=CE(I1)                                                                
      IF(ETA-10.0) 131,133,133                                                  
  131 ETA2=ETA*ETA  $  ETA2A=2.0*ETA  $  ETA6=ETA2+16.0                         
      SIGMA0=-(ETA/(12.*ETA6))*(1.+(ETA2-48.)/(30.*(ETA6**2))+((ETA2            
     1-160.)*ETA2+1280.)/(105.*(ETA6**4)))-ETA+(ETA/2.)*ALOG(ETA6)              
     2+3.5*ATAN(0.25*ETA)-(ATAN(ETA)+ATAN(0.5*ETA)+ATAN(ETA/3.))                
      GO TO 135                                                                 
  133 EINV1=1.0/ETA  $  EINV2=EINV1*EINV1  $  EINV3=EINV1*EINV2                 
      EINV5=EINV3*EINV2  $  EINV7=EINV5*EINV2  $  EINV9=EINV7*EINV2             
      SIGMA0=0.7853981634+ETA*ALOG(ETA)-ETA                                     
     1      -(0.08333333333*EINV1+0.00277777777*EINV3                           
     2       +0.00079365079*EINV5+0.00059523810*EINV7                           
     3       +0.00084175084*EINV9)                                              
  135 MODTPI=SIGMA0/6.2831853072                                                
      FMODTP=MODTPI                                                             
      SGMAZZ(I1)=SIGMA0-6.2831853072*FMODTP                                     
  150 CONTINUE                                                                  
      NDFMES=1  $  TX=2.0*XMES2                                                 
      DO 165 N=1,10                                                             
      TX=0.5*TX                                                                 
      IF(TX-0.00001-XMES1) 167,167,165                                          
  165 NDFMES=NDFMES+1                                                           
  167 XMES1=TX                                                                  
      A1=TMAS**0.333333333333                                                   
      IF(KTRL(2).EQ.1) A1=A1+PMAS**0.3333333333                         5/28/76 
      XBAR=RZERO*A1                                                             
      NXMAX=KEXCOM(2)                                                   10/8/75 
      NXCPLE=NXMAX                                                      10/8/75 
      FNXMAX=NXMAX  $  XMAX=DX*FNXMAX  $  KEX14=KEXCOM(14)              9/18/75 
      IIMAX=IIMAX*2                                                     10/18/75
      JJJMAX=KEX14-(IIMAX+ISTRTW-KTRL(1))/2                             10/18/75
      WRITE(6,263) (ECM(I1),I1=1,IICPLE)                                        
      WRITE(6,265) (WN(I1),I1=1,IICPLE)                                         
      WRITE(6,267) (WNINI(I1),I1=1,IICPLE)                                      
      WRITE(6,269) (CE(I1),I1=1,IICPLE)                                         
      WRITE(6,271) (SGMAZZ(I1),I1=1,IICPLE)                                     
      WRITE(6,273) JJJMAX,KEXCOM(14),NXCPLE,NXMAX,XMAX,XBAR                     
  263 FORMAT(7H ECM   6E15.5)                                                   
  265 FORMAT(7H WN    6E15.5)                                                   
  267 FORMAT(7H WNINI 6E15.5)                                                   
  269 FORMAT(7H ETA   6E15.5)                                                   
  271 FORMAT(7H SIGM0 6E15.5)                                                   
  273 FORMAT(27H0JJJMAX,KEX14,NXCPLE,NXMAX=4I5,11H XMAX,XBAR=2E15.7)            
      DO 510 I1=1,IICPLE  $  ETA=CE(I1)  $  ETASQ=ETA*ETA                       
      SG=SGMAZZ(I1)  $  EXSGRI(1,I1)=COS(SG)+TTI*SIN(SG)  $  FL=1.0             
      DO 510 L=1,LMAXM1                                                         
      DENOM=SQRT(1.0/(ETASQ+FL*FL))                                             
      EXSGRI(L+1,I1)=(FL+TTI*ETA)*EXSGRI(L,I1)*DENOM                            
  510 FL=FL+1.0                                                                 
      IF(KTLOUT(6)) 4110,4190,4110                                              
 4110 DO 4120 I1=1,IICPLE  $  WRITE(6,4115) I1                                  
 4115 FORMAT(15H0EXSGRI FOR I1=I1)                                              
      WRITE(6,4060) (EXSGRI(L,I1),L=1,LMAX)                                     
 4120 CONTINUE                                                                  
 4190 CONTINUE                                                                  
 4060 FORMAT(4(1X,1H(E13.6,E14.6,1H)))                                          
      CALL POTENT                                                               
 1000 RETURN                                                                    
      END                                                                       
      SUBROUTINE POTENT                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON /POTPAR/VSXD,WSXD,WSFD,DFND,DFNWD,DFNSD,RZED,RZEWD,RZESD,          
     1               VSXR,WSXC,WSFC,DFNC,DFNWC,DFNSC,RZEC,RZEWC,RZESC,          
     2               RZECD,DFNLD,DFNWLD,ALD,AWLD,                               
     3               RZECC,DFNLC,DFNWLC,ALC,AWLC                                
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(420,3)  $  TYPE COMPLEX EXSGRI                  10/21/75
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(5),             
     2           VFACI(5),VFACO(5)                                              
CCCCCC      ******     DIMENSION FIELDS USED ONLY IN THIS ROUTINE     **********
      COMMON    PFORM1(4),PFORM2(4),PFORM3(4),PFORM4(4),PFORM5(4),              
     1          VLAMD1(3,3),VLAMD2(3,3),VTERM1(3,3),VTERM2(3,3)                 
     2         ,XMEM(1200)                                              12/1/75 
      TYPE INTEGER VWRITE                                                       
      VSX=VSXD  $  WSX=WSXD  $  WSF=WSFD  $  DFN=DFND  $  DFNW=DFNWD            
      DFNS=DFNSD  $  RZERO=RZED  $  RZEROW=RZEWD  $  RZEROS=RZESD               
      RZEROC=RZECD                                                              
      IF (INTYPE.EQ.4) GO TO 40                                         4/28/76 
      BETA0=VCOUPL(1)                                                   4/28/76 
      IF (BETA0.EQ.0.0) BETA0=VCOUPL(3)                                 4/28/76 
      IF (BETA0.EQ.0.0) GO TO 40                                        4/28/76 
      DO 30 I=10,20                                                     4/28/76 
   30 VCOUPL(I)=VCOUPL(I)/BETA0                                         4/28/76 
      DO 35 I=1,4                                                       4/28/76 
   35 VCOUPL(I+4)=VCOUPL(I)/BETA0                                       4/28/76 
   40 NTTOTL=KEXCOM(2)+4*(NDFMES-1)                                             
      DO 88 ND=1,NDFMES                                                 4/28/76 
      FND=2.0**(ND-1)                                                           
      DX=XMES1*FND                                                              
      IF(ND-1) 84,82,84                                                 4/28/76 
   82 X=0.0                                                                     
      NXIMIN=1                                                                  
      NXIMAX=8                                                                  
      IF (NDFMES-1)86,83,86                                             4/28/76 
   83 NXIMIN=NXMAX+2                                                            
      GO TO 86                                                                  
   84 NXIMIN=NXIMAX+1                                                           
      NXIMAX=NXIMAX+4                                                           
      IF (ND-NDFMES)86,85,86                                                    
   85 NXIMAX=NTTOTL                                                     5/14/76 
   86 DO 88 NX=NXIMIN,NXIMAX                                            4/28/76 
      X=X+DX                                                                    
   88 XMEM(NX)=X                                                                
      NXMAX2=NTTOTL                                                     5/14/76 
      XBFAC=TMAS**0.333333333                                                   
      XBFAC1=XBFAC                                                      9/29/75 
      IF(KTRL(2).EQ.1) XBFAC=XBFAC+PMAS**0.333333333333                 5/28/76 
      RAT1=XBFAC1/XBFAC                                                 10/1/75 
      RAT2=RAT1*RAT1                                                    10/21/75
      GO TO (150,120,120,140,150), INTYPE                               5/5/76  
  120 VCOUPL(1)=VCOUPL(1)*RAT1                                          5/5/76  
      VCOUPL(2)=VCOUPL(2)*RAT1                                          5/5/76  
      VCOUPL(3)=VCOUPL(3)*RAT2                                          10/21/75
      VCOUPL(4)=VCOUPL(4)*RAT2                                          10/21/75
      GO TO 150                                                         5/5/76  
  140 DO 145 I1=1,INTMAX                                                5/5/76  
      I2=I1+INTMAX                                                      5/5/76  
      VCOUPL(I1)=VCOUPL(I1)*RAT1                                        5/5/76  
  145 VCOUPL(I2)=VCOUPL(I2)*RAT2                                        5/5/76  
  150 IF(KTRL(4).EQ.0) GO TO 200                                                
      READ(5,10) (VNUCLE(NX),NX=1,NTTOTL)                                       
      READ(5,10) (VCPL2R(NX),NX=1,NTTOTL)                                       
      READ(5,10) (VCPL1R(NX),NX=1,NTTOTL)                                       
      READ(5,10) (VCPL2I(NX),NX=1,NTTOTL)                                       
      VFAC1=VSX*VCOUPL(1)                                                       
      WFAC1=WSX*VCOUPL(1)                                                       
   10 FORMAT(5E14.6)                                                            
      DO 180 NX=1,NTTOTL                                                        
      VCENTR(NX)=-VNUCLE(NX)*VSX+VCPL2R(NX)                                     
      VCENTI(NX)=-VNUCLE(NX)*WSX                                                
      VCPL1R(NX)=-VCPL1R(NX)                                                    
      VCPL1I(NX)= VCPL1R(NX)*WFAC1                                              
      VCPL1R(NX)= VCPL1R(NX)*VFAC1+VCPL2I(NX)*VCOUPL(3)                         
  180 CONTINUE                                                                  
  200 XBARW=RZEROW*XBFAC                                                4/19/75 
      XBARS=RZEROS*XBFAC  $  XBARC=RZEROC*XBFAC                                 
      VCLFC2=1.4398650*CHARGE  $  VCLFC1=VCLFC2*0.5/XBARC                       
      VFAC1=-VSX*XBAR/DFN                                                       
      WFAC1=-WSX*XBARW/DFNW                                                     
      WDFC1=-4.0*WSF*XBARS/DFNS                                                 
      VFAC2=VFAC1*XBAR*0.5/DFN                                                  
      WFAC2=WFAC1*XBARW*0.5/DFNW                                                
      WDFC2=WDFC1*XBARS*0.5/DFNS                                                
      VFAC1=VFAC1*VCOUPL(1)                                             10/1/75 
      WFAC1=WFAC1*VCOUPL(1)                                             10/1/75 
      WDFC1=WDFC1*VCOUPL(1)                                             10/1/75 
      IF(KTRL(10).EQ.2) FACGDR=1.02333/XBAR                             10/9/80 
      IF(KTRL(10).EQ.2)  VFAC2=VFAC2*VCOUPL(1)                                  
      IF(KTRL(10).EQ.2) GO TO 215                                               
      VFAC2=VFAC2*(VCOUPL(1)**2)                                        5/19/76 
      WFAC2=WFAC2*(VCOUPL(1)**2)                                        5/19/76 
      WDFC2=WDFC2*(VCOUPL(1)**2)                                        5/19/76 
      IF (KTRL(4).NE.0) GO TO 410                                               
  215 IF(KTRL(11)) 257,258,257                                                  
  257 MESIZE=31  $  LAMAX=3  $  FN=30.0                                         
      DCOSI=1.0/FN   $  DCOSI3=DCOSI/3.0                                        
      BETA5=0.6307831*VCOUPL(1)  $  BETA9=0.846284367*VCOUPL(2)                 
      ALFC0= 3.5449078*DCOSI3  $  ALFC2= 7.9266549*DCOSI3                       
      ALFC4=10.6347234*DCOSI3  $  ALFC6=12.7813476*DCOSI3                       
  258 DO 400 NX=1,NXIMAX                                                        
      X=XMEM(NX)                                                                
      IF (KTRL(11).EQ.1) GO TO 285                                      4/28/76 
      PFORM1(1)=EXP((X-XBAR )/DFN )                                             
      PFORM1(2)=EXP((X-XBARW)/DFNW)                                             
      PFORM1(3)=EXP((X-XBARS)/DFNS)                                             
C     PFORM1(4)=EXP((X-XBARSP)/DFNSP)                                           
      DO 267 N=1,3                                                              
      PFORM2(N)=1.0/(1.0+PFORM1(N))                                             
      PFORM3(N)=PFORM1(N)*PFORM2(N)*PFORM2(N)                                   
      PFORM4(N)=(PFORM1(N)-1.0)*PFORM2(N)*PFORM3(N)                             
      PFORM5(N)=(1.0-PFORM1(N)*(4.0-PFORM1(N)))*PFORM2(N)                       
     1          *PFORM2(N)*PFORM3(N)                                            
  267 CONTINUE                                                                  
      VNUCLE(NX)=-VSX*PFORM2(1)                                                 
      VCENTR(NX)=VNUCLE(NX)                                             5/16/76 
      VCENTI(NX)=-WSX*PFORM2(2)-4.0*WSF*PFORM3(3)                               
      IF(KTRL(10).EQ.1)  GO TO 275                                              
      VCPL1R(NX)=VFAC1*PFORM3(1)                                                
      VCPL2R(NX)=VFAC2*PFORM4(1)                                                
      VCPL1I(NX)=WFAC1*PFORM3(2)+WDFC1*PFORM4(3)                                
      VCPL2I(NX)=WFAC2*PFORM4(2)+WDFC2*PFORM5(3)                                
      IF(KTRL(10).EQ.2) VCPL1R(NX)=(VCPL1R(NX)-VCPL2R(NX)*2./3.)*FACGDR         
      IF(KTRL(11).EQ.0) GO TO 400                                               
      GO TO 285                                                                 
  275 XR1=X/XBAR                                                                
      XW1=X/XBARW                                                               
      XS1=X/XBARS                                                               
      XR2=XR1*XR1                                                               
      XW2=XW1*XW1                                                               
      XS2=XS1*XS1                                                               
      VFC1X=VFAC1*XR1                                                           
      VFC2X=VFAC2*XR2                                                           
      WFC1X=WFAC1*XW1                                                           
      WDF1X=WDFC1*XS1                                                           
      WFC2X=WFAC2*XW2                                                           
      WDF2X=WDFC2*XS2                                                           
      VCPL1R(NX)=VFC1X*PFORM3(1)                                                
      VCPL2R(NX)=VFC2X*PFORM4(1)                                                
      VCPL1I(NX)=WFC1X*PFORM3(2)+WDF1X*PFORM4(3)                                
      VCPL2I(NX)=WFC2X*PFORM4(2)+WDF2X*PFORM5(3)                                
      IF( KTRL(3).EQ.0) GO TO 400                                               
      IF(NX.NE.1) GO TO 281                                             10/31/77
      EX1FC=  VCOUPL(1)*(ALD/DFNLD)                                     10/31/77
      EX2FC=  VCOUPL(1)/DFN                                             10/31/77
  281 EX1=EXP((X-XBAR)/DFNLD)                                           10/31/77
      EX2=EXP((X-XBAR)/DFN)                                             10/31/77
      EXFC1=1./(1.+ALD*EX1+EX2)                                         10/31/77
      EXFC2=(EX1FC*EX1+EX2FC*EX2)*(EXFC1**2)*X                          10/31/77
      VCENTR(NX)=VNUCLE(NX)=-VSX*EXFC1                                  10/31/77
      VCENTI(NX)           =-WSX*EXFC1                                  10/31/77
      VCPL1R(NX)=VSX*EXFC2                                              10/31/77
      VCPL1I(NX)=WSX*EXFC2                                              10/31/77
      GO TO 400                                                                 
  285 DO 287 LA=1,LAMAX                                                 4/28/76 
      DO 287 N=1,3                                                      4/28/76 
      VLAMD1(N,LA)=0.0  $  VLAMD2(N,LA)=0.0                                     
  287 CONTINUE                                                                  
      IIPA=1  $  COSI=-DCOSI                                                    
      DO 300 ME=1,MESIZE  $  COSI=COSI+DCOSI                                    
      COSI2=COSI*COSI $ COSI4=COSI2*COSI2 $ COSI6=COSI2*COSI4                   
      BR (1)=1.0  $  BR (2)=1.5*COSI2-0.5                                       
       BR (3)=4.375*COSI4-3.75*COSI2+0.375                                      
      BR (4)=14.4375*COSI6-19.6875*COSI4+6.5625*COSI2-0.3125                    
      BETAFC=1.0+BETA5*BR(2)+BETA9*BR(3)                                        
      PFORM1(1)=EXP((X-XBAR *BETAFC)/DFN )                                      
      PFORM1(2)=EXP((X-XBARW*BETAFC)/DFNW)                                      
      PFORM1(3)=EXP((X-XBARS*BETAFC)/DFNS)                                      
      DO 288 N=1,3                                                              
      PFORM2(N)=1.0/(1.0+PFORM1(N))                                             
      PFORM3(N)=PFORM1(N)*PFORM2(N)*PFORM2(N)                                   
      PFORM4(N)=(PFORM1(N)-1.0)*PFORM2(N)*PFORM3(N)                             
  288 CONTINUE                                                                  
      DO 289 LA=1,LAMAX                                                         
      VTERM1(1,LA)=PFORM2(1)*BR(LA)  $  VTERM2(1,LA)=PFORM3(1)*BR(LA)           
      VTERM1(2,LA)=PFORM2(2)*BR(LA)  $  VTERM2(2,LA)=PFORM3(2)*BR(LA)           
      VTERM1(3,LA)=PFORM3(3)*BR(LA)  $  VTERM2(3,LA)=PFORM4(3)*BR(LA)           
  289 CONTINUE                                                                  
      IF(ME-1) 291,295,291                                                      
  291 IF(ME-MESIZE) 292,295,292                                                 
  292 IF(IIPA) 293,294,293                                                      
  293 TERMFC=4.0  $  IIPA=0  $  GO TO 296                                       
  294 TERMFC=2.0  $  IIPA=1  $  GO TO 296                                       
  295 TERMFC=1.0                                                                
  296 DO 297 LA=1,LAMAX  $  DO 297 N=1,3                                        
      VLAMD1(N,LA)=VLAMD1(N,LA)+VTERM1(N,LA)*TERMFC                             
      VLAMD2(N,LA)=VLAMD2(N,LA)+VTERM2(N,LA)*TERMFC                             
  297 CONTINUE                                                                  
  300 CONTINUE                                                                  
      IF(KTRL(11).EQ.2)  GO TO 302                                              
      VNUCLE(NX)= -VSX*VLAMD1(1,1)*DCOSI3                               10/22/75
      VCENTR(NX)=VNUCLE(NX)                                             5/16/76 
      VCENTI(NX)=-(WSX*VLAMD1(2,1)+4.0*WSF*VLAMD1(3,1))*DCOSI3                  
  302 IF(KTRL(11).EQ.3)  GO TO 400                                              
      VCPL1R(NX)=-VSX*VLAMD1(1,2)*ALFC2                                         
      VCPL2R(NX)=-VSX*VLAMD1(1,3)*ALFC4                                         
      VCPL1I(NX)=(-WSX*VLAMD1(2,2)-4.0*WSF*VLAMD1(3,2))*ALFC2                   
      VCPL2I(NX)=(-WSX*VLAMD1(2,3)-4.0*WSF*VLAMD1(3,3))*ALFC4                   
  400 CONTINUE                                                                  
  410 VFACI(1)=0.0  $  VFACI(2)=0.0  $  VFACO(1)=0.0  $  VFACO(2)=0.0   4/5/75  
      XBARC2=XBARC*XBARC  $  XBARC3=XBARC*XBARC2                        4/4/75  
      XBARC4=XBARC*XBARC3                                               4/4/75  
      VCOUVX=4.319595*CHARGE                                            4/4/75  
      JLSCMX=1                                                          5/6/76  
      JLS=1                                                             5/1/76  
      LAMD(1)=IIREAD(2)                                                 5/1/76  
      NVJLS(1)=1                                                        5/1/76  
      GO TO (420,420,412,416,420,420),INTYPE                            5/1/76  
  412 NTPS=IICPLE-2                                                     5/1/76  
      DO 414 I1=3,IICPLE                                                5/1/76  
      JLS=JLS+1                                                         5/1/76  
      LAMD(JLS)=IIREAD(I1)                                              5/1/76  
      NVJLS(JLS)=1                                                      5/1/76  
      JLS1=JLS+NTPS                                                     5/1/76  
      LAMD(JLS1)=IIREAD(I1)                                             5/1/76  
  414 NVJLS(JLS1)=2                                                     5/1/76  
      JLSCMX=JLS                                                        5/6/76  
      JLS=NTPS*2+1                                                      5/1/76  
      GO TO 420                                                         5/1/76  
  416 IF (INTMAX.EQ.1) GO TO 420                                        5/1/76  
      DO 418 I=2,INTMAX                                                 5/1/76  
      JLS=JLS+1                                                         5/1/76  
      LAMD(JLS)=IIREAD(I)                                               5/1/76  
  418 NVJLS(JLS)=I                                                      5/1/76  
      JLSCMX=JLS                                                        5/19/76 
  420 JLSMAX=JLS                                                        5/1/76  
      IF (KTRL(13).EQ.0) GO TO 600                                      5/19/76 
      DO 430 I=1,JLSMAX                                                 5/6/76  
      II=LAMD(I)                                                        11/14/75
      XBARCL=XBARC**II                                                  11/17/75
      FACL=1.0/FLOAT(II*2+1)                                            11/14/75
      VCXL=VCOUVX*XBARCL*FACL                                           11/14/75
      VCXIL=VCOUVX*FACL/(XBARCL*XBARC)                                  11/14/75
      GO TO (430,421,422,423,430), INTYPE                               5/19/76 
  421 V1=VCOUPL(3)                                                      5/19/76 
      GO TO 426                                                         5/19/76 
  422 V1=VCOUPL(3)                                                      5/19/76 
      IF (NVJLS(I).EQ.1) GO TO 426                                      5/19/76 
      V1=VCOUPL(II+12)**2                                               5/19/76 
      VCXL=VCXL*0.5*(II+2)                                              5/19/76 
      VCXIL=VCXIL*0.5*(1-II)                                            5/19/76 
      GO TO 426                                                         5/19/76 
  423 I2=INTMAX+I                                                       5/19/76 
      V1=VCOUPL(I2)                                                     5/19/76 
  426 VFACI(I)=VCXIL*V1                                                 5/19/76 
      VFACO(I)=VCXL*V1                                                  5/19/76 
  430 CONTINUE                                                          5/19/76 
  600 IF(KTLOUT(7)) 4510,4590,4510                                      4/7/75  
 4510 WRITE(6,4515)                                                     6/22/76 
 4515 FORMAT(//,1X,21HWOODS-SAXON POTENTIAL,/)                                  
      KT7=KTLOUT(7)                                                     5/5/76  
      VWRITE=6HX       $  WRITE(6,4520) VWRITE                          5/5/76  
      WRITE(6,4525) (XMEM  (NX),NX=2,NXMAX2,KT7)                        5/5/76  
      VWRITE=6HVCENTR  $  WRITE(6,4520) VWRITE                                  
      WRITE(6,4525) (VCENTR(NX),NX=2,NXMAX2,KT7)                        5/5/76  
      VWRITE=6HVCENTI  $  WRITE(6,4520) VWRITE                                  
      WRITE(6,4525) (VCENTI(NX),NX=2,NXMAX2,KT7)                        5/5/76  
      VWRITE=6HVCPL1R  $  WRITE(6,4520) VWRITE                                  
      WRITE(6,4525) (VCPL1R(NX),NX=2,NXMAX2,KT7)                        5/5/76  
      VWRITE=6HVCPL1I  $  WRITE(6,4520) VWRITE                                  
      WRITE(6,4525) (VCPL1I(NX),NX=2,NXMAX2,KT7)                        5/5/76  
      VWRITE=6HVCPL2R  $  WRITE(6,4520) VWRITE                                  
      WRITE(6,4525) (VCPL2R(NX),NX=2,NXMAX2,KT7)                        5/5/76  
      VWRITE=6HVCPL2I  $  WRITE(6,4520) VWRITE                                  
      WRITE(6,4525) (VCPL2I(NX),NX=2,NXMAX2,KT7)                        5/5/76  
 4520 FORMAT(1X,A6)                                                             
 4525 FORMAT(10E12.4)                                                           
 4590 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      OVERLAY(JPWKB,2,0)                                                        
      PROGRAM JKLOOP                                                            
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ** 
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     FC(140,3),FDC(140,3),GC(140,3),GDC(140,3)                      
      COMMON     EXSGRI(140,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON     H(8400)                                                        
      COMMON CMAT(20,20)                                                        
      COMPLEX CMAT                                                              
      KEXCOM(48)=0                                                      5/2 /75 
      IF(KTRL(27)) 560,580,560                                                  
  560 REWIND 7  $  REWIND 8  $  GO TO 900                                       
  580 REWIND 7                                                                  
      WRITE(6,500)                                                      6/22/76 
  500 FORMAT(1H1)                                                               
      JJREPX=1                                                          4/7/75  
      IF (                  KTRL(16).GE.1)JJREPX=2                      10/18/75
      DO 850 JJREP=1,JJREPX                                             4/7/75  
      IF(JJREPX.EQ.1)  GO TO 590                                        10/17/75
      JJ1=2-JJREP                                                       10/18/75
      JJ2=JJREP-1                                                       10/18/75
      JMIN=JJ1+JJ2*(KTRL(16)+1)                                         10/18/75
      JMAX=JJ1* KTRL(16)   +JJ2*JJJMAX                                  10/18/75
      JINC=JJ1+JJ2*10                                                           
      GO TO 605                                                         10/17/75
  590 JMIN=1                                                            10/17/75
      JMAX=JJJMAX                                                       10/17/75
      JINC=1                                                            10/17/75
      IF(KTRL(9)) 601,605,601                                           4/7/75  
  601 JMIN=KEXCOM(5)  $  JMAX=KEXCOM(6)                                         
  605 DO 800 JJ=JMIN,JMAX,JINC                                                  
      KEXCOM(45)=JJ                                                             
      DO 790 K=1,2  $  IF(KTRL(9)-1) 608,606,608                                
  606 IF(K-KEXCOM(7)) 790,608,790                                               
  608 KEXCOM(46)=K  $  CALL NLJJJK  $  IF(MXROW) 790,790,720                    
  720 REWIND 8                                                                  
      IF (KTLOUT(3)-1) 780,725,723                                      6/22/76 
  723 IF(JJ.EQ.1.OR.JJ.EQ.JJJMAX) GO TO 725                             6/22/76 
      JJA=JJ/25                                                                 
      JJA=JJA*25+1                                                              
      IF(JJ.EQ.JJA) GO TO 725                                                   
      GO TO 780                                                                 
  725 J1CHEK=JJROW(1)-2*(JJROW(1)/2)  $  J1CKDV=2-J1CHEK                        
      DO 750 N=1,MXROW  $  LLROW(N)=LLROW(N)/2                                  
  750 JJROW(N)=JJROW(N)/J1CKDV                                                  
      WRITE(6,761) (NNROW(N),N=1,MXROW)                                         
      WRITE(6,762) (LLROW(N),N=1,MXROW)                                         
      IF(J1CHEK) 756,754,756                                                    
  754 WRITE(6,763) (JJROW(N),N=1,MXROW)   $  GO TO 770                          
  756 WRITE(6,765) (JJROW(N),N=1,MXROW)                                         
  761 FORMAT(3H0N=,23I5)                                                        
  762 FORMAT(3H L=,23I5)                                                        
  763 FORMAT(3H J=,23I5)                                                        
  765 FORMAT(4H J= ,23(I3,2H/2))                                                
  770 DO 775 N=1,MXROW  $  LLROW(N)=LLROW(N)*2                                  
  775 JJROW(N)=JJROW(N)*J1CKDV                                                  
  780 IF(KEXCOM(13).GE.JJ)  GO TO 785                                           
      CALL COUPLE                                                               
      GO TO 790                                                                 
  785 MXROWI=0                                                                  
      DO 786 M1=1,MXROW                                                         
      IF(NNROW(M1).GE.2)  GO TO 787                                             
      DO 786 M2=1,MXROW                                                         
      CMAT(M2,M1)=ZERO                                                          
  786 CONTINUE                                                                  
  787 MXROWI=M1-1                                                               
      DO 788 N1=1,MXROWI                                                        
      WRITE (7) (CMAT(M1,N1),M1=1,MXROW)                                        
  788 CONTINUE                                                                  
  790 CONTINUE                                                                  
  800 CONTINUE                                                                  
  850 CONTINUE                                                          4/7/75  
  900 RETURN                                                                    
      END                                                                       
      SUBROUTINE COUPLE                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),PR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXMACH,NANGLR,NDFMES,       9/18/75 
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(420,3)  $  TYPE COMPLEX EXSGRI                  10/21/75
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPLR(1200,2),VCPLI(1200,2)                           12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(5),             
     2           VFACI(5),VFACO(5)                                              
      COMMON     UCSAVE(20,20),UDSAVE(20,20)                            4/6/75  
      COMPLEX    UCSAVE,UDSAVE,EX                                       4/6/75  
      COMMON     UCR1(30),UCI1(30),UCR2(30),UCI2(30),                           
     1           URFORD(30,5),UIFORD(30,5),SPFACT(30),                          
     2           FPRER(5,30),FPREI(5,30),FPRERM(4,30),FPREIM(4,30),             
     3           UCR1M(30),UCI1M(30)                                            
      COMMON     DUMMY(4245)                                                    
      DIMENSION  FCRA(1,1,1 ),FDCRA(1,1,1), GCRA(1,1,1 ),GDCRA(1,1,1 )  5/1/76  
      DIMENSION VFACT(5),VC(5),WLDEP(500)                                       
      DATA KT/0/                                                                
      KWCTRL=0                                                                  
      IF(EXTCOM(10).NE.0.0) KWCRTL=1                                            
      IF(KT.EQ.0)  GO TO 10                                                     
      IF(KWCTRL.EQ.0) GO TO 10                                                  
      KT=1                                                                      
      DO 5 L1=1,500                                                             
      FL=L1-1                                                                   
      T1=1.0+EXP((FL-EXTCOM(9))/EXTCOM(10))                                     
      WLDEP(L1)=1.0/T1                                                          
    5 CONTINUE                                                                  
   10 MXROWI=0                                                          9/18/75 
      DO 20 M1=1,MXROW                                                          
      N1=NNROW(M1)                                                      9/18/75 
   14 IF(N1-1) 20,17,20                                                         
   17 MXROWI=MXROWI+1                                                           
   20 CONTINUE                                                                  
      KEXCOM(43)=MXROWI                                                 9/18/75 
      NXMACH=NXMAX                                                      9/18/75 
      NXRMM=NXMAX                                                       10/8/75 
      PRE41=19./240.  $  PRE42=-2./5.  $  PRE43=97./120.                        
      PRE44=-11./15.  $  PRE45=299./240.                                        
      NQPR=NQPRED=4  $  FNQPR=4.0                                               
      REWIND 8 $  CALL BMATRX                                           4/4/75  
       XBARC=XBAR*RZEROC/RZERO                                          5/14/76 
      VCLFC2=1.4398650*CHARGE                                           5/14/76 
      VCLFC1=VCLFC2*0.5/XBARC                                           5/14/76 
      M3MAX=MXROW                                                       9/18/75 
      INEX=1                                                            9/19/75 
      DO 1900 M3=1,M3MAX                                                        
      NDMAX=NDFMES                                                      9/18/75 
      DO 1890 ND=1,NDMAX                                                        
      FND=2.0**(ND-1)  $  DX=XMES1*FND                                  9/18/75 
      DR=DX*WN(1)   $  DRSQ=(DR*DR)/ECM(1)  $  K4COR2=4*ND-5                    
      IF(KEXCOM(1).EQ.0)  GO TO 810                                     10/8/75 
      IF(ND.LT.NDFMES)  GO TO 1890                                      5/5/75  
      DO 750 M1=1,MXROW                                                 5/5/75  
      UCR1(M1)=0.0                                                      5/5/75  
      UCR2(M1)=0.0                                                      5/5/75  
      UCI1(M1)=0.0                                                      5/5/75  
      UCI2(M1)=0.0                                                      5/5/75  
      DO 750 N=1,5                                                      5/5/75  
      FPRER (N,M1)=0.0                                                  5/5/75  
      FPREI (N,M1)=0.0                                                  5/5/75  
      FPRERM(N,M1)=0.0                                                  5/5/75  
      FPREIM(N,M1)=0.0                                                  5/5/75  
  750 CONTINUE                                                          5/5/75  
      FL1=LLROW(1)/2+1                                                  9/19/75 
      RT=(CE(1)+SQRT(CE(1)**2+FL1**2))/WN(1)                            9/19/75 
      R1=RT+EXTCOM(5)                                                   9/19/75 
      NXIMAX=R1/XMES2                                                   9/18/75 
      NXIMAX=MAX0(NXIMAX,KEXCOM(2))                                     10/18/75
      NXMACH=NXIMAX-2                                                   9/25/75 
      NXIMIN=(RT-EXTCOM(6))/XMES2                                       10/8/75 
      NXIMIN=MAX0(NXIMIN,KEXCOM(1),5)                                   10/8/75 
      NXMAX3=NXIMAX-5                                                   9/18/75 
      X=XMES2*(NXIMIN-2)                                                5/5/75  
      N1=NNROW(M3)  $  WNI=WN(N1)                                       5/5/75  
      LL=LLROW(M3)/2+1                                                  5/5/75  
      ETA=CE(N1)  $  RHOMX =X*WNI  $  L=LL-1                            9/29/75 
      CALL FLGLCH(F1,G1,FD1,GD1,L)                                      9/29/75 
      UCR2(M3)=F1                                                       9/29/7K 
      GO TO 1225                                                        5/5/75  
  810 IF(ND-1) 860,822,860                                                      
  822 DO 825 M1=1,MXROW                                                         
      UCR1(M1)=UCI1(M1)=UCR2(M1)=UCI2(M1)=0.0                                   
      DO 825 N=1,5                                                              
      FPRER(N,M1)=FPREI(N,M1)=FPRERM(N,M1)=FPREIM(N,M1)=0.0                     
  825 CONTINUE                                                                  
      IF(INEX-1) 925,830,925                                                    
  830 N1=NNROW(M3)  $  J1=JJROW(M3)  $  L1=LLROW(M3)                            
      LL=(L1/2)+1  $  WNI=WNINI(N1)                                             
      UCR2(M3)=0.5*((2.0*DX*WNI)**LL)*EXP(FACLOG(LL)-FACLOG(2*LL))              
      IF(LL-2) 840,834,840                                                      
  834 FPRER(1,M3)=FPRERM(1,M3)=CFUNIT*0.6666666666*WNI*WNI*DRSQ                 
  840 X=0.0  $  NXIMIN=2  $  NXIMAX=8  $  NXPRCH=3  $  NNX=2                    
      IF(NDFMES-1) 1225,857,1225                                                
  857 NXIMAX=NXRMM+2  $  NXMAX3=NXIMAX-5  $  GO TO 1225                         
  860 NXIMIN=5  $  X=3.0*DX                                                     
      IF(ND-NDFMES) 862,863,862                                                 
  862 NXIMAX=8  $  GO TO 865                                                    
  863 NXIMAX=NXRMM+2   $  NXMAX3=NXIMAX-5                                       
  865 DO 875 M1=1,MXROW                                                         
      DO 870 NQ=1,NQPR                                                          
      FPRER(NQ,M1)=FPRERM(NQ,M1)*4.0                                            
  870 FPREI(NQ,M1)=FPREIM(NQ,M1)*4.0                                            
      FPRERM(1,M1)=FPRER(1,M1)  $  FPRERM(2,M1)=FPRER(3,M1)                     
      FPREIM(1,M1)=FPREI(1,M1)  $  FPREIM(2,M1)=FPREI(3,M1)                     
      UCR1(M1)=UCR1M(M1)  $  UCI1(M1)=UCI1M(M1)                                 
  875 CONTINUE                                                                  
      NXPRCH=5  $  NNX=3                                                        
      GO TO 1225                                                                
  925 GR1=FI1=1.0  $  GI1=FR1=0.0  $  M4=M3                                     
      IF(M3-MXROW) 935,935,930                                                  
  930 GR1=GI1=FI1=0.0  $  FR1=1.0  $  M4=M3-MXROW                               
  935 N1=NNROW(M4)  $  L1=LLROW(M4)  $  J1=JJROW(M4)  $  LL=(L1/2)+1            
      FL=LL  $  FNX=NXMAX+3                                                     
      X=      FNX*XMES2  $  CFFCLL=FL*(FL-1)*CFUNIT                             
      FAC=L1+1  $   T1=SQRT(FAC*((WN(1)/WN(N1))**3))                            
      FAC=T1                                                                    
      DO 950 NQ=1,5                                                             
      X=X-XMES2                                                                 
      V1=(1.4398650*CHARGE/X+CFFCLL/(X*X)-ECM(N1))*DRSQ                         
      FRR   =GR1*GCRA(LL,NQ,N1)+FR1*FCRA(LL,NQ,N1)                              
      FII   =GI1*GCRA(LL,NQ,N1)+FI1*FCRA(LL,NQ,N1)                              
      EX=(FRR+TTI*FII)*FAC*EXSGRI(LL,N1)                                        
      FPRERR=EX  $  FPREII=AIMAG(EX)                                            
      FPRER(NQ,M4)=FPRERR*V1  $  FPREI(NQ,M4)=FPREII*V1                         
      IF(NQ-4) 950,945,947                                                      
  945 UCR1(M4)=FPRERR  $  UCI1(M4)=FPREII  $  GO TO 950                         
  947 UCR2(M4)=FPRERR  $  UCI2(M4)=FPREII                                       
  950 CONTINUE                                                                  
      NXIMIN=1  $  NXIMAX=NXMAX-NXRM  $  NXMAX3=NXIMAX-5                        
      X=(NXMAX-1)*XMES2  $  K4COR2=NXMAX+4*NDFMES-5                             
 1225 DO 1800 NX=NXIMIN,NXIMAX                                                  
      N3=(NX+K4COR2)*(2-INEX)+(K4COR2-NX)*(INEX-1)                              
      X=X+DX  $  XX2=X*X  $  XX2INV=1.0/XX2                                     
      VI1=0.0                                                           5/14/76 
      IF (X-XBARC) 1230,1230,1232                                       5/14/76 
 1230 VR1=VCLFC1*     (3.0-((X/XBARC)**2))                              5/14/76 
      GO TO 1235                                                        5/14/76 
 1232 VR1=VCLFC2/X                                                      5/14/76 
 1235 IF (NX.GE.KEXCOM(2)) GO TO 1240                                   5/14/76 
      VR1=VR1+VCENTR(N3)                                                5/14/76 
      IF(KTRL(4).NE.0)VR1=VCENTR(N3)                                            
      VI1=VCENTI(N3)                                                    5/14/76 
 1240 DO 1246 NL=1,JLSMAX                                               5/14/76 
      XX1=X**LAMD(NL)                                                   5/14/76 
      XX1INV=1.0/XX1/X                                                  5/14/76 
      IF (X-XBARC) 1242,1242,1244                                       5/14/76 
 1242 VC(NL)=VFACI(NL)*XX1                                              5/14/76 
      GO TO 1246                                                        5/14/76 
 1244 VC(NL)=VFACO(NL)*XX1INV                                           5/14/76 
 1246 CONTINUE                                                          5/14/76 
      DO 1700 M1=1,MXROW                                                        
      L1=LLROW(M1) $ LL=(L1/2)+1 $ CFFCLL=(LL-1)*LL*CFUNIT $N1=NNROW(M1)        
      WCOREC=WC(N1)                                                     10/17/75
      IF(KWCTRL.EQ.1) WCOREC=WCOREC*WLDEP(LL)                                   
      AR2=VR1       +CFFCLL*XX2INV-ECM(N1)                              5/14/76 
      AI2=VI1       *WCOREC                                             5/14/76 
      BR2=AR2*UCR2(M1)-AI2*UCI2(M1)                                     10/17/75
      BI2=AR2*UCI2(M1)+AI2*UCR2(M1)                                             
      BR(M1)=BR2*DRSQ  $  BI(M1)=BI2*DRSQ                                       
CCCCCC                                                                          
CCCCCC      ******      COUPLING          ***********************               
CCCCCC                                                                          
 1330 IF(NX.LE.KEXCOM(3)) GO TO 1700                                            
      DO 1500 NL=1,JLSMAX                                                       
      NV=NVJLS(NL)                                                      4/5/75  
      VI2=0.0                                                           5/19/76 
      VR2=VC(NL)                                                        5/19/76 
      IF (NX.GE.KEXCOM(2)) GO TO 1335                                   5/14/76 
      VR2=VR2+VCPLR(N3,NV)                                              5/19/76 
      IF(KTRL(4).NE.0)VR2=VCPLR(N3,NV)                                          
      VI2=VCPLI(N3,NV)                                                  5/19/76 
 1335 AFAC1R=DRSQ*VR2                                                   5/19/76 
      AFAC1I=DRSQ*VI2                                                   5/19/76 
      AR1=AFAC1R*AMAT(M1,NL)                                            4/4/75  
      AI1=AFAC1I*AMAT(M1,NL)                                            4/4/75  
      BR(M1)=BR(M1)+AR1*UCR2(M1)                                        4/5/75  
      BI(M1)=BI(M1)+AR1*UCI2(M1)                                        4/5/75  
      IF(KTRL(12).EQ.0)  GO TO 1340                                     4/5/75  
      BR(M1)=BR(M1)-AI1*UCI2(M1)                                        4/5/75  
      BI(M1)=BI(M1)+AI1*UCR2(M1)                                        4/5/75  
 1340 IF(M1.EQ.1)  GO TO 1500                                           4/4/75  
      M12MIN=M12MX1(M1-1,NL)+1                                          4/4/75  
      M12MAX=M12MX1(M1,NL)                                              4/4/75  
      DO 1460 M12=M12MIN,M12MAX                                         4/4/75  
      M2=MCOL(M12,NL)                                                   4/4/75  
      AR1=AFAC1R*AMAT(M12,NL)                                           4/4/75  
      AI1=AFAC1I*AMAT(M12,NL)                                           4/4/75  
      BR(M1)=BR(M1)+AR1*UCR2(M2)                                        4/5/75  
      BI(M1)=BI(M1)+AR1*UCI2(M2)                                        4/5/75  
      BR(M2)=BR(M2)+AR1*UCR2(M1)                                        4/5/75  
      BI(M2)=BI(M2)+AR1*UCI2(M1)                                        4/5/75  
      IF(KTRL(12).EQ.0)  GO TO 1410                                     4/5/75  
      BR(M1)=BR(M1)-AI1*UCI2(M2)                                        4/5/75  
      BI(M1)=BI(M1)+AI1*UCR2(M2)                                        4/5/75  
      BR(M2)=BR(M2)-AI1*UCI2(M1)                                        4/5/75  
      BI(M2)=BI(M2)+AI1*UCR2(M1)                                        4/5/75  
 1410 CONTINUE                                                          4/4/75  
 1460 CONTINUE                                                          4/4/75  
 1500 CONTINUE                                                          4/4/75  
 1700 CONTINUE                                                                  
CCCCCC                                                                          
CCCCCC                                                                          
CCCCCC                                                                          
      IF(KTLOUT(10)-2)4190,4110,4190                                            
 4110 IF(KTRL(9).NE.0)  GO TO 4120                                      5/5/75  
      IF(LLROW(1)-2)  4120,4120,4190                                    5/5/75  
 4120 WRITE(6,4125) NX,X,(UCR2(M),UCI2(M),UCR1(M),UCI1(M),M=1,2)                
 4125 FORMAT(4H NX=I3,3H X=F8.4,8E12.4)                                         
 4190 CONTINUE                                                                  
      DO 1740 M1=1,MXROW                                                        
      IF(INEX-1) 1714,1711,1714                                                 
 1711 IF(ND-1) 1714,1712,1714                                                   
 1712 IF(NX-4) 1713,1713,1714                                                   
 1713 FPRER(NX,M1)=TERMR=BR(M1)  $  FPREI(NX,M1)=TERMI=BI(M1)                   
      GO TO 1715                                                                
 1714 FPRER(5,M1)=BR(M1)  $  FPREI(5,M1)=BI(M1)                                 
      TERMR =PRE41*FPRER (1,M1)+PRE42*FPRER (2,M1)+PRE43*FPRER (3,M1)           
     1      +PRE44*FPRER (4,M1)+PRE45*FPRER (5,M1)                              
      TERMI =PRE41*FPREI (1,M1)+PRE42*FPREI (2,M1)+PRE43*FPREI (3,M1)           
     1      +PRE44*FPREI (4,M1)+PRE45*FPREI (5,M1)                              
 1715 ARR=2.0*UCR2(M1)-UCR1(M1)+TERMR                                           
      AII=2.0*UCI2(M1)-UCI1(M1)+TERMI                                           
      UCR1(M1)=UCR2(M1)  $  UCR2(M1)=ARR                                        
      UCI1(M1)=UCI2(M1)  $  UCI2(M1)=AII                                        
      IF(INEX-1) 1735,1717,1735                                                 
 1717 IF(ND-NDFMES) 1718,1729,1718                                              
 1718 IF(NX-NXPRCH) 1729,1719,1729                                              
 1719 FPRERM(NNX,M1)=BR(M1)  $  FPREIM(NNX,M1)=BI(M1)                           
      IF(M1-MXROW) 1725,1723,1725                                               
 1723 NXPRCH=NXPRCH+2  $  NNX=NNX+1                                             
 1725 IF(NX-(NXIMAX-1)) 1729,1728,1729                                          
 1728 UCR1M(M1)=UCR1(M1)  $  UCI1M(M1)=UCI1(M1)                                 
 1729 IF(ND-1) 1735,1730,1735                                                   
 1730 IF(NX-4) 1740,1740,1735                                                   
 1735 DO 1738 NQ=1,NQPR                                                         
      FPRER(NQ,M1)=FPRER(NQ+1,M1)                                               
 1738 FPREI(NQ,M1)=FPREI(NQ+1,M1)                                               
 1740 CONTINUE                                                                  
      IF(ND-NDMAX ) 1800,1760,1800                                              
 1760 NXCH2=NX-NXMAX3  $  IF(NXCH2)1800,1800,1765                               
 1765 DO 1770 M1=1,MXROW                                                        
      URFORD(M1,NXCH2)=UCR2(M1)  $  UIFORD(M1,NXCH2)=UCI2(M1)                   
 1770 CONTINUE                                                                  
 1800 CONTINUE                                                                  
CCCCCC      ******      CALCULATION OF THE DERIVATIVES      ******              
      IF(ND-NDMAX ) 1890,1805,1890                                              
 1805 DO 1850 M1=1,MXROW                                                        
      UCR1(M1)=URFORD(M1,3)  $  UCI1(M1)=UIFORD(M1,3)                           
      N1=NNROW(M1)  $  DENOM=1.0/(12.0*DX*WN(N1))                               
      UCR2(M1)   =(8.0*(URFORD(M1,4)-URFORD(M1,2))                              
     1                -(URFORD(M1,5)-URFORD(M1,1)))*DENOM                       
      UCI2(M1)   =(8.0*(UIFORD(M1,4)-UIFORD(M1,2))                              
     1                -(UIFORD(M1,5)-UIFORD(M1,1)))*DENOM                       
      UCSAVE(M1,M3)=UCR1(M1)+UCI1(M1)*TTI                               5/1/76  
      UDSAVE(M1,M3)=UCR2(M1)+UCI2(M1)*TTI                               5/1/76  
 1850 CONTINUE                                                                  
 1890 CONTINUE                                                                  
 1900 CONTINUE                                                                  
 1950 CONTINUE                                                                  
      CALL WKBCTR                                                               
 2000 RETURN                                                                    
      END                                                                       
      SUBROUTINE BMATRX                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(140,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON     FC(140,3),FDC(140,3),GC(140,3),GDC(140,3)                      
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(5),             
     2           VFACI(5),VFACO(5)                                              
      COMMON     B(10,10,10),KII(24),JLSD(24)                           4/4/75  
      COMMON     IITYPE(10),KQNOTW(7),ALP23P(7,3),ALSQ3P(7,3)           4/4/75  
      JJ=KEXCOM(45)  $ KTRL1T=KTRL(1)                                   4/4/75  
      KT1PIS=KTRL1T+ISTRTW  $  KTISCK=KT1PIS-2*(KT1PIS/2)                       
      JJTRTW=2*JJ-2+KTISCK                                                      
      DO 50 N1=1,6  $  DO 50 N2=1,6  $  DO 50 K1=1,10                           
   50 B(N1,N2,K1)=0.0                                                           
      COMFAC=0.2820948                                                          
      I=5                                                               5/1/76  
      IF (VCOUPL(1).EQ.0.0) I=7                                         5/1/76  
      J=I+1                                                             5/1/76  
      B(1,2,1)=COMFAC*VCOUPL(I)                                         4/28/76 
      B(2,2,1)=COMFAC*VCOUPL(J)                                         4/28/76 
      GO TO (1000,1000,310,410,1000),INTYPE                             4/28/76 
CCCCCC                                                                          
CCCCCC     **********     INTYPE=3     *****************************            
CCCCCC                                                                          
  310 COMFC2=COMFAC/2.5066282                                                   
      NTPHST=IICPLE-2                                                   5/6/76  
      DO 340 N1=3,IICPLE                                                5/1/76  
      I1=IIREAD(N1)                                                     4/28/76 
      N2=I1+11                                                          4/28/76 
      N3=N2+1                                                           4/28/76 
      N4=(I1/2)+2                                                       4/28/76 
      HAT=I1+I1+1                                                       4/28/76 
      HAT=SQRT(HAT)                                                             
      B(2,N1,1)=COMFAC*VCOUPL(N2)*HAT*0.6324555                                 
      N22=(I1/2)+8                                                              
      N44=N22-3                                                                 
      B(1,N1,N44)=COMFAC*VCOUPL(N22+ 9)                                 4/28/76 
      I1M2=I1-2  $  IF(I1M2)331,333,335                                         
  331 CLEB1=0.4472136  $  GO TO 337                                             
  333 CLEB1=-0.53452243  $  GO TO 337                                           
  335 CLEB1=0.71713714                                                          
  337 VC12=VCOUPL(N3)                                                           
      B(1,N1,N4)=(VC12**2)*COMFC2*CLEB1                                         
  340 CONTINUE                                                                  
      INTMAX=7 $  IF(KTRL(5)) 351,375,351                               5/1/76  
  351 ALFADG=VCOUPL(1)**2  $  ALFADG=ALFADG*COMFC2                              
      B(1,1,2)=ALFADG*0.7070107  $  B(2,2,2)=ALFADG*2.2135946                   
      B(2,2,3)=ALFADG*0.7559290  $  B(2,2,4)=ALFADG*1.01418507                  
      DO 370 K1=2,4                                                             
      K1M1=K1-1  $  GO TO (353,355,357),K1M1                                    
  353 FKFAC=1.118034  $  FKHAT=1.0  $  IF=0  $  GO TO 360                       
  355 FKFAC=0.59761437  $  FKHAT=2.236068  $  IF=4  $  GO TO 360                
  357 FKFAC=0.59761437  $  FKHAT=3.0  $  IF=8                                   
  360 IA=4  $  IC=4  $  IE=4                                                    
      DO  370 N1=1,NTPHST  $  I1=IIREAD(N1+2)  $  IB=I1+I1                      
      HAT1SQ=IB+1  $  HAT1=SQRT(0.2*HAT1SQ)                                     
      DO  370 N2=N1,NTPHST  $  IF(N1-N2) 369,365,369                            
  365 IF(K1-2) 369,367,369                                                      
  367 B(N1+2,N2+2,2)=HAT1                                                       
  369 I2=IIREAD(N2+2) $ ID=I2+I2 $ HAT2SQ=ID+1 $ HAT2=SQRT(3.2*HAT2SQ)          
      CALL RAC7  $  SNDTRM=HAT1*HAT2*FKHAT*RAC                                  
      B(N1+2,N2+2,K1)=FKFAC*VCOUPL(2*N1+10)*VCOUPL(2*N2+10)             4/28/76 
     1        *(B(N1+2,N2+2,K1)+SNDTRM)                                         
  370 CONTINUE                                                                  
  375 GO TO 1000                                                                
CCCCCC                                                                          
CCCCCC     **********     INTYPE=4     *****************************            
CCCCCC                                                                          
  410 ROTFAC=0.2820948                                                          
      ID=IF=2*IIREAD(1)-KTRL(1)  $  IE=0                                        
      DO 450 N3=1,INTMAX                                                        
      IB=2*LAMD(N3)                                                     5/1/76  
      DO 430 N1= 1,IICPLE  $  I1=IIREAD(N1)  $  IC=2*I1-KTRL(1)                 
      DO 430 N2=N1,IICPLE  $  I2=IIREAD(N2)  $  IA=2*I2-KTRL(1)                 
      HATA=IA+1  $  HATA=SQRT(HATA)                                             
      CALL CLEB                                                                 
  430 B(N1,N2,N3)=ROTFAC*HATA*RAC                                               
  450 CONTINUE                                                          5/1/76  
 1000 CALL AMATRX                                                               
 2000 RETURN                                                                    
      END                                                                       
      SUBROUTINE AMATRX                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(140,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON     FC(140,3),FDC(140,3),GC(140,3),GDC(140,3)                      
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(5),             
     2           VFACI(5),VFACO(5)                                              
      COMMON     B(10,10,10),KII(24),JLSD(24)                           4/4/75  
      KIREPT=KEXCOM(28)  $  JJ=KEXCOM(45)                                       
      KTRL1T=KTRL(1)              $  KT1PIS=KTRL1T+ISTRTW                       
      KTISCK=KT1PIS-2*(KT1PIS/2)  $  JJTRTW=2*JJ-2+KTISCK                       
      MJBRTW=-JJTRTW                                                            
      DO 950 JLS=1,JLSMAX                                               4/4/75  
      NL=JLS                                                            4/4/75  
      K1=JLS                                                                    
      KI=LAMD(JLS)                                                              
      IF (INTYPE.NE. 3) GO TO 100                                       5/6/76  
      IF (JLS.EQ.1) GO TO 100                                           5/6/76  
      NV1=NVJLS(JLS)                                                    5/6/76  
      K1=KI/2+2+(2-NV1)*3                                               5/6/76  
  100 DO 900 M1=1,MXROW                                                         
      IU=0                                                                      
      N1=NNROW(M1) $ I1=IIREAD(N1) $ J1=JJROW(M1) $ L1=LLROW(M1)                
      DO 890 M2=1,M1                                                            
      N2=NNROW(M2) $ I2=IIREAD(N2) $ J2=JJROW(M2) $ L2=LLROW(M2)                
      T1=0.0                                                                    
      B1=B(N2,N1,K1)*((-1.0)**(I1-I2+KI))                               4/5/75  
      IF(B1.EQ.0.0)  GO TO 750                                          4/5/75  
      KITRTW=KI+KI $ IA=I1+I1-KTRL(1) $ IB=J1 $ IC=I2+I2-KTRL(1)                
      ID=J2 $ IE=JJTRTW $ IF=KITRTW $ CALL RAC7 $ RAC1=RAC                      
      IA=L1 $ IB=L2 $ IC=KITRTW $ ID=0 $ IE=0 $ IF=0 $ CALL CLEB                
      RAC2=RAC                                                                  
      IA=L1 $ IB=J1 $ IC=L2 $ ID=J2 $ IE=ISTRTW $ IF=KITRTW                     
      CALL RAC7 $  RAC3=RAC                                                     
      HAT=(L1+1)*(L2+1)*(J1+1)*(J2+1)  $  HAT=SQRT(HAT)                         
      S1=(-1.0)**((JJTRTW-ISTRTW-I2*2+KTRL(1)+L1+L2)/2                          
     1          +((L1-L2)/2+KPRITR(N1)-KPRITR(N2))/2)                           
      T2=S1*HAT*RAC1*RAC2*RAC3                                          4/4/74  
      T1=T1+T2*B1                                                       4/4/75  
  750 IF(M1.NE.M2)  GO TO 810                                                   
      AMAT(M1,NL)=T1                                                            
      MCOL(M1,NL)=M1                                                    4/5/75  
      IF(M1.EQ.1)  M12=MXROW                                                    
      GO TO 890                                                                 
  810 IF(T1.EQ.0.0)  GO TO 890                                                  
      IU=1                                                                      
      M12=M12+1                                                                 
      MCOL(M12,NL)=M2                                                           
      AMAT(M12,NL)=T1                                                           
  890 CONTINUE                                                                  
      IF(IU.EQ.1)  GO TO 895                                                    
      IF(M1.EQ.1)  GO TO 900                                            4/5/75  
      M12=M12+1                                                                 
      MCOL(M12,NL)=1                                                            
      AMAT(M12,NL)=0.0                                                          
  895 M12MX1(M1,NL)=M12                                                         
  900 CONTINUE                                                                  
      M12MX1(1,NL)=MXROW                                                4/5/75  
      M12MX2(NL)=M12                                                            
  950 CONTINUE                                                                  
 4050 FORMAT(8E15.5)                                                            
 4060 FORMAT(1H0)                                                               
      IF(KTLOUT(9).EQ.0)  GO TO 4270                                    5/4/76  
      WRITE(6,4230) INTYPE,INTMAX,(LAMD(N),N=1,JLSMAX)                  5/3/76  
 4230 FORMAT(21H0  BMATRX FOR INTYPE=I1,8H INTMAX=I1,7H LAMDA=10I3)             
      DO 4231 K1=1,JLSMAX                                               5/3/76  
      WRITE(6,4060)                                                             
      DO 4231 N1=1,IICPLE                                                       
      WRITE(6,4050) (B(N1,N2,K1),N2=1,IICPLE)                                   
 4231 CONTINUE                                                                  
      DO 4260 NL=1,JLSMAX                                               5/4/76  
      M12MX=M12MX2(NL)                                                          
      M12MI3=1                                                          5/3/76  
      M12MX3=8                                                          5/3/76  
      KREP=(M12MX-1)/8+1                                                5/13/76 
      DO 4245 K=1,KREP                                                  5/3/76  
      IF (K.EQ.KREP) M12MX3=M12MX                                       5/3/76  
      IF (K.GT.1) GO TO 4240                                            5/3/76  
      WRITE(6,4235) NL, (MCOL(M12,NL),AMAT(M12,NL),M12=M12MI3,M12MX3)   5/3/76  
 4235 FORMAT(1H0,21HMCOL AND AMAT FOR NL=, I2,4H ARE,8(I3,E10.3))       5/3/76  
 4236 FORMAT(1H,27X, 8(I3,E10.3))                                       5/3/76  
      GO TO 4242                                                        5/3/76  
 4240 WRITE(6,4236) (MCOL(M12,NL),AMAT(M12,NL),M12=M12MI3,M12MX3)       5/3/76  
 4242 M12MI3=M12MI3+8                                                   5/3/76  
      M12MX3=M12MX3+8                                                   5/3/76  
 4245 CONTINUE                                                          5/3/76  
 4260 CONTINUE                                                                  
 4270 I1=I1                                                             5/4/76  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WKBCTR                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXMACH,NANGLR,NDFMES,       9/18/75 
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(420,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(5),             
     2           VFACI(5),VFACO(5)                                              
      COMMON      URI(20,20),UDRI(20,20),Q(20,40),R(20,40),QRI(100,5),  4/29/76 
     1            FRI(20,20,4),SMAT(20,20)                              4/29/76 
      COMPLEX     URI,UDRI,FRI,QRI                                      4/29/76 
      COMPLEX    Q,R,HRI1,HRI2,HRI3,HRIINV,ERI1,ERI2,ERI3                       
      DIMENSION   V2(100,5),WNRATO(70),CMAT(20,20)                      4/28/76 
      COMPLEX SMAT,CMAT                                                         
      EQUIVALENCE (CMAT(1,1),FRI(1,1,1))                                5/4/76  
      EQUIVALENCE (V2(1,1),AMAT(1,1))                                   4/28/76 
      EQUIVALENCE (WNRATO(1),AR(1))                                     4/28/76 
CCCCCC                                                                          
CCCCCC                                                                          
CCCCCC                                                                          
      KBOUND=0  $  NBOUND=KEXCOM(42)  $  MXROWI=KEXCOM(43)                      
      JJ=KEXCOM(45)  $  K=KEXCOM(46)                                            
      JJTRTW=KEXCOM(47)                                                         
      JJRT=JJTRTW/2                                                             
      X=XMES2*FLOAT(NXMACH)                                             9/18/75 
      W=X*WN(1)                                                                 
      JLSCMX=JLSMAX                                                     4/29/76 
      NTPS=0                                                            4/29/76 
      IF (INTYPE.EQ.3) NTPS=IICPLE-2                                    4/29/76 
      JLSCMX=JLSMAX-NTPS                                                4/29/76 
      MXROW1=MXROW+1                                                            
      DO   5 I1=1,IICPLE                                                        
      WNR=WN(I1)/WN(1)  $  WNRATO(I1)=WNR                                       
    5 CONTINUE                                                                  
      KWRIT2=0                                                                  
      XBARC=XBAR*RZEROC/RZERO                                                   
      VCOUVX=4.319595*CHARGE*WN(1)/ECM(1)                               11/17/75
      DO  30 NL=1,JLSCMX                                                5/1/76  
      II=LAMD(NL)                                                       11/17/75
      XBARCL=XBARC**II                                                  11/17/75
      FACL=VCOUVX/FLOAT(II*2+1)                                         11/17/75
      BETA=VCOUPL(3)                                                    5/1/76  
      IF (INTYPE.EQ.4) BETA=VCOUPL(INTMAX+NL)                           5/1/76  
      AR1=FACL*XBARCL*BETA                                              5/1/76  
      M12MAX=M12MX2(NL)                                                 4/5/75  
      DO  20 M12=1,M12MAX                                               4/5/75  
      V2(M12,NL)=AR1   *AMAT(M12,NL)                                    5/1/76  
      IF (INTYPE.NE.3) GO TO 20                                         5/1/76  
      IF (NL.EQ.1) GO TO 20                                             5/12/76 
      NL1=NL+NTPS                                                       5/1/76  
      V2(M12,NL)=V2(M12,NL)+AR1*AMAT(M12,NL1)*(LAMD(NL1)+2)/2.0         5/1/76  
   20 CONTINUE                                                                  
      IF(KTLOUT(14).EQ.0)   GO TO  30                                           
      M12MI3=1                                                          5/3/76  
      M12MX3=8                                                          5/3/76  
      KREP=(M12MAX-1)/8+1                                               5/13/76 
      KREP=M12MAX/8+1                                                   5/3/76  
      DO 28 K=1,KREP                                                    5/3/76  
      IF (K.EQ.KREP) M12MX3=M12MAX                                      5/3/76  
      IF(K.GT.1)  GO TO 24                                              5/3/76  
      WRITE(6,4000) NL, (MCOL(M12,NL),AMAT(M12,NL),M12=M12MI3,M12MX3)   5/3/76  
 4000 FORMAT(1H0,21HMCOL AND AMAT FOR NL=, I2,4H ARE,8(I3,E10.3))       5/3/76  
      GO TO 26                                                          5/3/76  
   24 WRITE(6,4010) (MCOL(M12,NL),AMAT(M12,NL),M12=M12MI3,M12MX3)       5/3/76  
 4010 FORMAT(1H,27X, 8(I3,E10.3))                                       5/3/76  
   26 M12MI3=M12MI3+8                                                   5/3/76  
      M12MX3=M12MX3+8                                                   5/3/76  
   28 CONTINUE                                                          5/3/76  
   30 CONTINUE                                                                  
CCCCCC                                                                          
CCCCCC     *****     CALCULATION OF SMAT AND CMAT     ****************          
CCCCCC                                                                          
      RHO  =XMES2*NXMACH                                                9/22/75 
      DO 40 M1=1,MXROW                                                  4/28/76 
      N=NNROW(M1)                                                       4/28/76 
      ETA=CE(N)                                                         9/18/75 
      W2=WN(N)  $  RHOMX =W2*RHO                                        9/22/75 
      L=LLROW(M1)/2                                                     4/28/76 
      CALL FLGLCH(F1,G1,FD1,GD1,L)                                      9/18/75 
      F(M1)=F1                                                          4/28/76 
      G(M1)=G1                                                          4/28/76 
      FD(M1)=FD1                                                        4/28/76 
   40 GD(M1)=GD1                                                        4/28/76 
      IF(KTLOUT(14)) 43,49,43                                                   
   43 IF(KTRL(9)) 48,44,48                                                      
   44 IF(JJ-3) 48,48,49                                                         
   48 KWRIT2=1                                                                  
   49 NMATCH=(1-KBOUND)+KBOUND*MXROWI                                           
      IM=MXROW1                                                                 
 4050 FORMAT(1H ,6X,8E15.6)                                                     
      IF(KWRIT2) 4120,4190,4120                                                 
 4120 WRITE(6,4150)                                                             
 4150 FORMAT(8H0URC ETC)                                                        
      DO 4160 M3=1,MXROW                                                        
      WRITE(6,4050) (URI(M1,M3),UDRI(M1,M3),M1=1,MXROW)                 5/1/76  
 4160 CONTINUE                                                                  
 4190 CONTINUE                                                                  
      DO 90 I=1,MXROW                                                           
      HRI1=G(I)+F(I)*TTI                                                        
      HRI2=GD(I)+FD(I)*TTI                                                      
      DO 80 J=1,MXROW                                                           
      Q(I,J)=UDRI(I,J)*HRI1-URI(I,J)*HRI2                               5/1/76  
      JM=J+MXROW                                                                
      Q(I,JM)=0.0                                                               
      IF(I.EQ.J)  Q(I,JM)=1.0                                                   
   80 CONTINUE                                                                  
   90 CONTINUE                                                                  
      CALL SIMLEQ                                                       4/9/75  
      DO 96 I=1,MXROW                                                           
      DO 94 J=1,MXROW                                                           
      ERI1=0.0                                                                  
      DO 92 L=1,MXROW                                                           
      ERI1=ERI1+(URI(I,L)*FD(I)-UDRI(I,L)*F(I))*R(L,J)                  5/1/76  
   92 CONTINUE                                                                  
      Q(I,J)=ERI1                                                               
   94 CONTINUE                                                                  
   96 CONTINUE                                                                  
      IF(KWRIT2.EQ.0) GO TO 4250                                                
      WRITE(6,4200)                                                             
 4200 FORMAT(17H0INVERSION MATRIX//1H )                                         
      DO 4220 M3=1,MXROW                                                        
      WRITE(6,4050) (Q(M1,M3),M1=1,MXROW)                                       
 4220 CONTINUE                                                                  
 4250 CONTINUE                                                                  
      ERI3=TTI*2.0                                                              
      DO 120 M1=1,MXROW                                                         
      N1=NNROW(M1)                                                              
      L1=LLROW(M1)/2+1                                                          
      HRI1=EXSGRI(L1,N1)                                                        
      FL1=LLROW(M1)+1                                                           
      T1=SQRT(WN(N1)**3/FL1)                                                    
      ERI1=T1/HRI1                                                              
      DO 110 M2=1,MXROW                                                         
      N2=NNROW(M2)                                                              
      L2=LLROW(M2)/2+1                                                          
      HRI2=EXSGRI(L2,N2)                                                        
      FL2=LLROW(M2)+1                                                           
      T2=SQRT(WN(N2)**3/FL2)                                                    
      ERI2=HRI2/T2                                                              
      CMAT(M1,M2)=Q(M1,M2)*ERI1*ERI2                                            
      SMAT(M1,M2)=Q(M1,M2)*ERI3*SQRT(WN(N1)/WN(N2))                     4/21/75 
  110 CONTINUE                                                                  
      SMAT(M1,M1)=SMAT(M1,M1)+TTR                                               
  120 CONTINUE                                                                  
CCCCCC                                                                          
CCCCCC       ********     INITIAL AMPLITUDE     **********                      
CCCCCC                                                                          
      DO 150 M1=1,MXROW                                                 5/3/76  
      DO 150 M3=1,MXROW                                                 5/3/76  
      URI(M1,M3)=0.0                                                    4/28/76 
      IF(M1.EQ.M3) URI(M1,M3)=1.0                                       4/28/76 
  150 CONTINUE                                                          4/28/76 
      IF(KTRL(14).EQ.1) GO TO 180                                               
      DO 170 I=1,JLSMAX                                                 5/3/76  
      NF=MXROW                                                                  
      V2 (1,I)=V2(1,I)*0.5                                              4/28/76 
      DO 160 M1=2,MXROW                                                         
      N1=NNROW(M1)                                                              
      V2(M1,I)=V2(M1,I)/(WNRATO (N1)*2.0)                                       
      NA=NF+1                                                                   
      NF=M12MX1(M1,I)                                                   4/5/75  
      DO 160  J=NA,NF                                                           
      M2=MCOL(J,I)                                                      4/5/75  
      N2=NNROW(M2)                                                              
      V2(J,I)=V2(J,I)*SQRT(CE(N1)*CE(N2))/(CE(1)*2.0)                           
  160 CONTINUE                                                                  
  170 CONTINUE                                                                  
CCCCCC                                                                          
CCCCCC     **********     OUT PUT       ***************************             
CCCCCC                                                                          
  180 IF(KTLOUT(14).EQ.0)  GO TO 225                                            
      WRITE (6,4310)                                                    5/3/76  
 4310 FORMAT(1H0,34HC-MATRIX FOR THE INTERIOR SOLUTION)                 5/5/76  
      DO 4320 M2=1,MXROW                                                5/3/76  
      WRITE (6,4300) (CMAT(M1,M2),M1=1,MXROW)                           5/3/76  
 4320 CONTINUE                                                          5/3/76  
      WRITE (6,4330)                                                    5/3/76  
 4330 FORMAT(1H0,34HS-MATRIX FOR THE INTERIOR SOLUTION)                 5/5/76  
      DO 4340 M2=1,MXROW                                                5/3/76  
      WRITE (6,4300) (SMAT(M1,M2),M1=1,MXROW)                           5/3/76  
 4340 CONTINUE                                                          5/3/76  
  225 CONTINUE                                                                  
      IF(KTRL(14).EQ.0) GO TO 260                                       5/28/76 
      NLJMAX=KEXCOM(43)                                                 4/7/75  
      DO 250 NLJ=1,NLJMAX                                               4/7/75  
      WRITE(7) (CMAT(M1,NLJ),M1=1,MXROW)                                4/7/75  
  250 CONTINUE                                                          4/7/75  
      IF (KTLOUT(3).LE.1) WRITE(6,255) (CMAT(M1,1),M1=1,MXROW)                  
  255 FORMAT(1H0,9HC-MATRIX,10E11.3/10X,10E11.3/10X,10E11.3)            5/4/76  
      GO TO 3000                                                        4/7/75  
  260 MXRW2=MXROW*2                                                     4/7/75  
      IF(KTLOUT(14).EQ.0)  GO TO 270                                    4/12/75 
      WRITE (6,4350)                                                    5/3/76  
 4350 FORMAT(1H0,18HINITIAL AMPLITUDES)                                 5/3/76  
      WRITE (6,4300)  ((URI(I,J  ),I=1,MXROW),J=1,MXROW)                5/1/76  
 4300 FORMAT(1H ,6X,10E11.3)                                            5/3/76  
  270 CALL WKBINT                                                       4/12/75 
 3000 RETURN                                                            4/7/75  
      END                                                               ABOY0751
      SUBROUTINE WKBINT                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON     FH(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXMACH,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(140,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON     FC(140,3),FDC(140,3),GC(140,3),GDC(140,3)                      
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(15)             
      COMMON      URI(20,20),UDRI(20,20),Q(20,40),R(20,40),QRI(100,5),  4/29/76 
     1            FRI(20,20,4)                                          4/29/76 
      COMPLEX     URI,UDRI,FRI,QRI                                      4/29/76 
      DIMENSION   Q1(20,20)                                             4/29/76 
      EQUIVALENCE (Q1(1,1),FRI(1,1,4))                                  4/29/76 
      COMPLEX    Q,R,ARI ,HRI2,HRI3,HRIINV,ERI1,ERI2,ERI3               4/29/76 
      COMPLEX    Q1,RK1,RK2,RK3                                         4/29/76 
      DATA RB1,RB2,C1,C2/ 0.5857864,3.4142136,0.1213204,-4.1213204/     10/22/75
      DATA G1,G2,ACCUR,ACCO50/3.0,2.0,196.E-12,196.E-16/                10/22/75
      DATA ITEREP/20/                                                   10/22/75
      RB23=RB2*0.66666666667                                            10/22/75
      C23=C2*0.66666666667                                              10/22/75
      JJTRTW=KEXCOM(47)                                                 10/22/75
      FNXMX=NXMACH                                                              
      FACTO=1.0                                                                 
      ICHECK=0                                                                  
      IRUN=0                                                                    
      ITE=0                                                                     
      RWKBMX=EXTCOM(1)                                                  4/29/76 
      ROMAWB=200.0                                                      6/8/76  
      DX1=50.0                                                          6/8/76  
      X1=ROMAWB                                                         4/28/75 
      DRMAWB=DX1                                                        4/28/75 
      DX=XMES2                                                          4/28/75 
      X=XMES2*FNXMX                                                     4/28/75 
      MXRW2=MXROW+MXROW                                                         
      NAM=0                                                                     
      IF(KTLOUT(15).EQ.0)  GO TO 5                                              
      WRITE (6,4010) RWKBMX  ,ROMAWB,X1,DX1                             4/29/76 
 4010 FORMAT(1H0,* RWKBMX  ,ROMAWB,X1,DX1=*,5F10.5)                     4/29/76 
CCCCCC                                                                          
C     THE RUNGE-KUTTA-GILL INTEGRATION PROCEDURE                        ABOY0776
CCCCCC                                                                          
    5 CALL CMATRX(X)                                                    4/28/75 
      CALL AMPDER                                                               
   10 DX=DX*FACTO                                                       4/28/75 
      D2X=DX*2.0                                                        4/28/75 
      DXD3=DX*0.33333333                                                4/28/75 
      DX3D8=DX*0.75                                                     4/28/75 
      DO 20 L=1,MXROW                                                   4/29/76 
      DO 20 M1=1,MXROW                                                  4/29/76 
      FRI(M1,L,1)=UDRI(M1,L)                                            4/29/76 
   20 CONTINUE                                                                  
      DO 120 NAM=2,4                                                            
      DO 40 L=1,MXROW                                                   4/29/76 
      DO 30 M1=1,MXROW                                                  4/29/76 
      Q1(M1,L)=DX*UDRI(M1,L)                                            4/29/76 
      URI(M1,L)=URI(M1,L)+Q1(M1,L)                                      4/29/76 
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
      IF(KTLOUT(15).EQ.0.OR .NAM.NE.2)  GO TO 42                                
      WRITE (6,4000)  NAM,X ,DX                                         4/28/75 
      WRITE (6,4100)  ((URI(I,J  ),I=1,MXROW),J=1,MXROW)                5/1/76  
 4000 FORMAT(1H0,*   NAM,X,DX=*,I5,2F10.5)                              4/28/75 
 4100 FORMAT(1H ,6X,8E15.6)                                                     
   42 X=X+DX                                                            4/28/75 
      CALL CMATRX(X)                                                    4/28/75 
      CALL AMPDER                                                               
      DO 60 L=1,MXROW                                                   4/29/76 
      DO 50 M1=1,MXROW                                                  4/29/76 
      RK1=DX*UDRI(M1,L)                                                 4/29/76 
      URI(M1,L)=URI(M1,L)+RB1*(RK1-Q1(M1,L))                            4/29/76 
      Q1(M1,L)=Q1(M1,L)*C1 +RB1*RK1                                     4/29/76 
   45 CONTINUE                                                                  
   50 CONTINUE                                                                  
   60 CONTINUE                                                                  
      CALL AMPDER                                                               
      DO 80 L=1,MXROW                                                   4/29/76 
      DO 70 M1=1,MXROW                                                  4/29/76 
      RK2=DX*UDRI(M1,L)                                                 4/29/76 
      URI(M1,L)=URI(M1,L)+RB2*(RK2-Q1(M1,L))                            4/29/76 
      Q1(M1,L)=RB23*RK2+C23*Q1(M1,L)                                    4/29/76 
   65 CONTINUE                                                                  
   70 CONTINUE                                                                  
   80 CONTINUE                                                                  
   82 X=X+DX                                                            4/28/75 
      CALL CMATRX(X)                                                    4/28/75 
      CALL AMPDER                                                               
      DO 100 L=1,MXROW                                                  4/29/76 
      DO 90 M1=1,MXROW                                                  4/29/76 
      RK3=DXD3*UDRI(M1,L)                                               4/29/76 
      URI(M1,L)=URI(M1,L)    +(RK3-Q1(M1,L))                            4/29/76 
   85 CONTINUE                                                                  
   90 CONTINUE                                                                  
  100 CONTINUE                                                                  
      CALL AMPDER                                                               
      IF (X.GT.RWKBMX) GO TO 210                                        4/29/76 
      DO 110 L=1,MXROW                                                  4/29/76 
      DO 110 M1=1,MXROW                                                 4/29/76 
      FRI(M1,L,NAM)=UDRI(M1,L)                                          4/29/76 
  105 CONTINUE                                                                  
  110 CONTINUE                                                                  
  112 IF(X.LT.X1)  GO TO 120                                            4/28/75 
      X1=X1+DX1                                                         4/28/75 
  120 CONTINUE                                                                  
      GO TO 150                                                         ABOY0832
CCCCCC                                                                          
CCCCCC                                                                          
CCCCCC                                                                          
  130 ITE=ITE+1                                                                 
      DO 140 L=1,MXROW                                                  4/29/76 
      DO 140 M1=1,MXROW                                                 4/29/76 
      DO 138 N=1,3                                                      4/29/76 
  138 FRI(M1,L,N)=FRI(M1,L,N+1)                                         4/29/76 
      FRI(M1,L,4)=UDRI(M1,L)                                            4/29/76 
  140 CONTINUE                                                                  
  150 DO 160 L=1,MXROW                                                  4/29/76 
      DO 160 M1=1,MXROW                                                 4/29/76 
      URI(M1,L)=URI(M1,L)+DX3D8*(6.111111111*FRI(M1,L,4)-6.55555555556  4/29/76 
     1          *FRI(M1,L,3)+4.111111111*FRI(M1,L,2)-FRI(M1,L,1))       4/29/76 
  160 CONTINUE                                                                  
      X=X+D2X                                                           4/28/75 
      CALL CMATRX(X)                                                    4/28/75 
      AAPMAX=0.                                                         ABOY0849
      CALL AMPDER                                                               
      DO 170 L=1,MXROW                                                  4/29/76 
      DO 170 M1=1,MXROW                                                 4/29/76 
      ARI=DX3D8*(UDRI(M1,L)-4.0*(FRI(M1,L,4)+FRI(M1,L,2))               4/29/76 
     1    +6.0*FRI(M1,L,3)+FRI(M1,L,1))                                 4/29/76 
      AA1=CABS(ARI) **2                                                 4/29/76 
      IF (AA1.GT.AAPMAX)  AAPMAX=AA1                                            
      URI(M1,L)=URI(M1,L)+ARI                                           4/29/76 
  170 CONTINUE                                                                  
      CALL AMPDER                                                               
      IF (X.GT.RWKBMX) GO TO 210                                        4/29/76 
      FACTO=1.                                                          ABOY0865
      IF (AAPMAX.LT.ACCO50)  FACTO=2.                                   ABOY0866
      IF (AAPMAX.GT.ACCUR)  FACTO=0.5                                   ABOY0867
      IF(FACTO.EQ.1.0)   GO TO 178                                              
      DW=DX*WN(1)                                                       4/28/75 
      IF(DW*FACTO.GE.2.*CE(1)) GO TO 178                                        
      IF(KTLOUT(15).GE.1)  WRITE(6,4200)  X,DX,AAPMAX                   4/28/75 
 4200 FORMAT(1H0,*  X,DX,AAPMAX=*,2F10.5,E10.3)                         4/28/75 
      GO TO 10                                                                  
  178 IF(X.LT.X1)  GO TO 180                                            4/28/75 
      CALL RMATRX(0,1,X)                                                4/28/75 
      CALL CHEKRM(ICHECK,IRUN,IOK)                                              
      IF(IOK.EQ.1)  GO TO 210                                                   
      IF(IRUN.NE.0)  GO TO 130                                                  
      X1=X1+DX1                                                         4/28/75 
  180 IF(X.LT.RWKBMX) GO TO 130                                         4/29/76 
  210 IOK=1                                                                     
      EXTCOM(50)=DX                                                     4/28/75 
      CALL RMATRX(1,1,X)                                                4/28/75 
      A1=X                                                              4/28/75 
      IF (JJTRTW.NE.0) RETURN                                                   
      ROMAWB=X-3.*DRMAWB                                                4/28/75 
      RETURN                                                            ABOY0884
      END                                                               ABOY0885
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      SUBROUTINE CMATRX(X)                                                      
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON     FH(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    ****** 
      COMMON     EXSGRI(140,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON     FC(140,3),FDC(140,3),GC(140,3),GDC(140,3)                      
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(15)             
      COMMON      URI(20,20),UDRI(20,20),Q(20,40),R(20,40),QRI(100,5),  4/29/76 
     1            FRI(20,20,4)                                          4/29/76 
      COMPLEX     URI,UDRI,FRI,QRI                                      4/29/76 
      DIMENSION  FAK(20),PHAS(20)                                       5/1/76  
      COMPLEX    Q,R,HRI1,HRI2,HRI3,HRIINV,ERI1,ERI2,ERI3                       
      DATA PIH/1.5707963/                                               ABOY0894
      W=X*WN(1)                                                         4/28/75 
      W1=CE(1)*W                                                                
      DO 10 M1=1,MXROW                                                          
      N1=NNROW(M1)                                                              
      WX=W1/CE(N1)                                                              
      LLTR=LLROW(M1)/2                                                          
      ALLL=(LLTR+1)*LLTR                                                        
      S1=    (WX**2-(W1+W1)-ALLL)                                       5/3/75  
      IF(S1.GT.0.0)  GO TO 8                                            5/3/75  
      FAK(M1)=0.0                                                       5/3/75  
      PHAS(M1)=0.0                                                      5/3/75  
      GO TO 10                                                          5/3/75  
    8 QS=SQRT(S1)                                                       5/3/75  
      ETAL=SQRT(CE(N1)**2+ALLL)                                                 
      FAK(M1)=SQRT(WX/QS)                                                       
      ARG=(W1+ALLL)/(WX*ETAL)                                           ABOY0905
      IF (ARG.GT.1.0)  ARG=1.                                           ABOY0906
      ALFSQ=SQRT(ALLL)                                                          
      PHASIS=QS-ALFSQ    *ACOS(ARG)                                             
   10 PHAS(M1)=PHASIS+CE(N1)*ALOG(ETAL/(WX+QS-CE(N1)))                          
      IF(KTLOUT(16).EQ.1) WRITE (6,4000)                                        
 4000 FORMAT(1H0,2X,9HQRI(I,NL))                                                
      DO 30 NL=1,JLSCMX                                                 5/20/76 
      NE=MXROW                                                                  
      I1=LAMD(NL)                                                       5/20/76 
      W1=X**(I1+1)                                                      5/20/76 
      W3=1.0/W1                                                         5/29/76 
      QRI(1,NL)=FAK(1)**2*W3* AMAT(1,NL)*TTI                            4/29/76 
      DO 20 I=2,MXROW                                                           
      NA=NE+1                                                           ABOY0913
      QRI(I,NL)=FAK(I)**2*W3*AMAT(I,NL) *TTI                            4/29/76 
      NE=M12MX1(I,NL)                                                   4/5/75  
      DO 20 J=NA,NE                                                             
      IND=MCOL(J,NL)                                                    4/5/75  
      VV=W3*AMAT(J,NL)                                                  5/1/76  
      PHASE=PHAS(IND)-PHAS(I)                                           ABOY0919
      TREA=COS(PHASE)                                                   ABOY0920
      TIMA=SIN(PHASE)                                                   ABOY0921
      FAKT=FAK(IND)*FAK(I)*VV                                           ABOY0922
      QRI(J,NL)=(TIMA+TREA*TTI)*FAKT                                    4/29/76 
   20 CONTINUE                                                                  
      IF(KTLOUT(16).EQ.0)  GO TO 30                                             
      WRITE (6,4100) NL,(QRI(I,NL),I=1,NE)                                      
 4100 FORMAT(1H ,3HNL=,I5,12E10.2/6X,12E10.2)                           5/4/76  
   30 CONTINUE                                                                  
      RETURN                                                            ABOY0927
      END                                                               ABOY1179
      SUBROUTINE RMATRX(ICC,IWK,X)                                              
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON     FH(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(420,3)  $  TYPE COMPLEX EXSGRI                  9/23/75 
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(15)             
      COMMON      URI(20,20),UDRI(20,20),Q(20,40),R(20,40),QRI(100,5),  4/29/76 
     1            FRI(20,20,4),SMAT(20,20)                              4/29/76 
      COMPLEX     URI,UDRI,FRI,QRI                                      4/29/76 
      DIMENSION   CMAT(20,20),REF(20),ABREF(20),PHREF(20)                       
      COMPLEX      SMAT ,CMAT,REF                                               
      EQUIVALENCE  (CMAT(1,1),URI(1,1))                                 5/4/76  
      COMPLEX    Q,R,HRI1,HRI2,HRI3,HRIINV,ERI1,ERI2,ERI3                       
      DATA ONE/1./                                                      ABOY0983
CCCCCC                                                                          
CCCCCC   *******     CALCULATION OF S-2 MATRIX     ***************              
CCCCCC                                                                          
      JJ=KEXCOM(45)                                                             
      J=0                                                                       
      W=X*WN(1)                                                         4/28/75 
      ISW=0                                                                     
      JJTRTW=KEXCOM(47)                                                         
      IF(JJTRTW.EQ.0)  IUNIT=0                                                  
      DO 10 I=1,MXROW                                                           
      DO 10 L=1,MXROW                                                           
      LS=L+MXROW                                                                
      Q(I,L)=URI(I,L)                                                   5/1/76  
      Q(I,LS)=0.0                                                               
      IF(I.EQ.L)  Q(I,LS)=1.0                                                   
   10 CONTINUE                                                                  
      CALL SIMLEQ                                                               
      DO 20 I=1,MXROW                                                           
      DO 20 J=1,MXROW                                                           
      Q(J,I)=CONJG(URI(I,J))                                            5/1/76  
   20 CONTINUE                                                                  
      IF(KTLOUT(18).EQ.0.OR. ICC.EQ .0)  GO TO 25                               
      WRITE (6,4010)                                                            
 4010 FORMAT(1H0,*  RMAT-2*)                                                    
      DO 4020 M2=1,MXROW                                                        
      WRITE (6,4100) (R   (M1,M2),M1=1,MXROW)                                   
 4020 CONTINUE                                                                  
 4100 FORMAT(1H ,6X,10E11.3)                                            9/25/75 
      WRITE (6,4025)                                                            
 4025 FORMAT(1H0,*  R-ADJ *)                                            9/25/75 
      DO 4030 M2=1,MXROW                                                        
      WRITE (6,4100) (Q(M1,M2),M1=1,MXROW)                                      
 4030 CONTINUE                                                                  
   25 CONTINUE                                                                  
      IF(ICC.EQ.0)  GO TO 2000                                                  
CCCCCC                                                                          
CCCCCC     **********     CALCULATION OF FINAL C-MATRX     **********           
CCCCCC                                                                          
      ERI2= TTI*0.5                                                             
      DO 60 I=1,MXROW                                                           
      DO 55 J=1,MXROW                                                           
      ERI1=ZERO                                                                 
      DO 52 K=1,MXROW                                                   9/25/75 
      DO 52 L=1,MXROW                                                           
   52 ERI1=ERI1+SMAT(K,L)*R(L,J)*R(I,K)                                 9/27/75 
      IF(JJ.LE.KEXCOM(4))  ERI1=0.0                                             
      Q(I,J)=ERI1                                                       9/25/75 
      CMAT(I,J)=-ERI1*ERI2                                                      
   55 CONTINUE                                                                  
      CMAT(I,I)=CMAT(I,I)+ERI2                                                  
   60 CONTINUE                                                                  
      DO 220 M1=1,MXROW                                                         
      N1=NNROW(M1)                                                              
      L1=LLROW(M1)/2+1                                                          
      LTRAN=2*IIREAD(N1)                                                        
      HRI1=EXSGRI(L1,N1)                                                        
      FL1=LLROW(M1)+1                                                           
      T1=WN(N1)/SQRT(FL1)                                                       
      ERI1=T1/HRI1                                                              
      DO 210 M2=1,MXROW                                                         
      N2=NNROW(M2)                                                              
      L2=LLROW(M2)/2+1                                                          
      HRI2=EXSGRI(L2,N2)                                                        
      FL2=LLROW(M2)+1                                                           
      T2=WN(N2)/SQRT(FL2)                                                       
      ERI2=HRI2/T2                                                              
      CMAT(M1,M2)=CMAT(M1,M2)*ERI1*ERI2                                 4/21/75 
  210 CONTINUE                                                                  
  220 CONTINUE                                                                  
      ERQ=REAL(Q(1,1))                                                          
      EIQ=AIMAG(Q(1,1))                                                         
      PHASEQ=ATAN2(EIQ,ERQ)                                                     
      ABSQ=CABS(Q(1,1))                                                         
      NLJMAX=KEXCOM(43)                                                         
      DO 100 NLJ=1,NLJMAX                                                       
      WRITE (7) (CMAT(M1,NLJ),M1=1,MXROW)                                       
  100 CONTINUE                                                                  
      IF (KTLOUT(3)-2) 230,228,240                                              
  228 IF(JJ.EQ.1.OR.JJ.EQ.JJJMAX) GO TO 230                             6/22/76 
      JJA=JJ/25                                                                 
      JJA=JJA*25+1                                                              
      IF(JJ.EQ.JJA) GO TO 230                                                   
      GO TO 2000                                                                
  230 WRITE (6,4200) X, (CMAT(M1,1),M1=1,MXROW)                         5/4/76  
 4200 FORMAT(1H0,11HC-MATRIX(X=,F5.1,1H),10E11.3/18X,10E11.3/18X,       5/4/76  
     1       10E11.3)                                                   5/4/76  
      DO 235 M1=1,MXROW                                                         
      A=0.0                                                                     
      IF (M1.EQ.1) A=1.0                                                        
      REF(M1)=2.0*TTI*CMAT(M1,1)+A                                              
      ER=REAL(REF(M1))                                                          
      EI=AIMAG(REF(M1))                                                         
      PHREF(M1)=ATAN2(EI,ER)                                                    
  235 ABREF(M1)=CABS(REF(M1))                                                   
      WRITE (6,4205) (ABREF(M1),PHREF(M1),M1=1,MXROW)                           
 4205 FORMAT(1H0,11HS-MATRIX   ,6X,10E11.3/18X,10E11.3/18X,10E11.3)             
      GO TO 2000                                                        10/20/75
  240 DO 300 M2=1,MXROW                                                 9/25/75 
      WRITE(6,4100) (CMAT(M1,M2),M1=1,MXROW)                            9/25/75 
  300 CONTINUE                                                          9/25/75 
      WRITE(6,4210)                                                     9/25/75 
 4210 FORMAT(1H0,*   FINAL S-MATRX *)                                   9/25/75 
      DO 301 M2=1,MXROW                                                 9/25/75 
      WRITE(6,4100) (Q   (M1,M2),M1=1,MXROW)                            9/25/75 
  301 CONTINUE                                                          9/25/75 
 2000 RETURN                                                            ABOY1066
      END                                                                       
      SUBROUTINE AMPDER                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON     FH(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(140,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON     FC(140,3),FDC(140,3),GC(140,3),GDC(140,3)                      
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(15)             
      COMMON      URI(20,20),UDRI(20,20),Q(20,40),R(20,40),QRI(100,5),  4/29/76 
     1            FRI(20,20,4)                                          4/29/76 
      COMPLEX     URI,UDRI,FRI,QRI                                      4/29/76 
      COMPLEX    Q,R,HRI1,HRI2,HRI3,HRIINV,ERI1,ERI2,ERI3                       
      DO 20 M1=1,MXROW                                                  4/29/76 
      DO 20 L=1,MXROW                                                           
      HRI1=ZERO                                                         4/29/76 
      DO 10 NL=1,JLSCMX                                                 5/20/76 
      HRI1=HRI1+QRI(M1,NL)*URI(M1,L)                                    4/29/76 
   10 CONTINUE                                                                  
      UDRI(M1,L)=HRI1                                                   4/29/76 
   15 CONTINUE                                                                  
   20 CONTINUE                                                                  
      DO 50 NL=1,JLSCMX                                                 5/20/76 
      NE=MXROW                                                                  
      DO 40 M1=2,MXROW                                                  4/29/76 
      NA=NE+1                                                           ABOY0951
      NE=M12MX1(M1,NL)                                                  4/29/76 
      DO 35 IS=NA,NE                                                            
      M2 =MCOL(IS,NL)                                                   4/29/76 
      DO 30 L=1,MXROW                                                           
      UDRI(M1,L)=UDRI(M1,L)+QRI(IS,NL)*URI(M2,L)                        4/29/76 
      UDRI(M2,L)=UDRI(M2,L)-CONJG(QRI(IS,NL))*URI(M1,L)                 5/4/76  
   30 CONTINUE                                                                  
   35 CONTINUE                                                                  
   40 CONTINUE                                                                  
   50 CONTINUE                                                                  
      IF(KTLOUT(17).EQ.0)  GO TO 100                                            
      WRITE (6,4000)                                                            
 4000 FORMAT(1H0,*  AMPDOT  *)                                                  
      WRITE (6,4100)  ((URI(I,J  ),I=1,MXROW),J=1,MXROW)                4/29/76 
 4100 FORMAT(1H ,10E10.3)                                                       
  100 RETURN                                                                    
      END                                                               ABOY0968
      SUBROUTINE SIMLEQ                                                         
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON     FH(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(140,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON     FC(140,3),FDC(140,3),GC(140,3),GDC(140,3)                      
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(15)             
      COMMON      URI(20,20),UDRI(20,20),Q(20,40),R(20,40),QRI(100,5),  4/29/76 
     1            FRI(20,20,4)                                          4/29/76 
      COMPLEX     URI,UDRI,FRI,QRI                                      4/29/76 
      DIMENSION ARI(70),BRI(70) $ EQUIVALENCE (AR,ARI),(BR,BRI)                 
      COMPLEX    Q,R,ARI,BRI,EX,DIVQ,SUMRI                                      
      MXROWM=MXROW                                                      9/22/75 
      IM=MXROW+MXROW                                                            
      DO 870 I=1,MXROW                                                          
      BIGQ=CABS(Q(I,I))  $  MXI=I                                               
      DO 816 I1=I,MXROW                                                         
      B1=CABS(Q(I1,I))  $  IF(BIGQ-B1) 815,816,816                              
  815 BIGQ=B1  $  MXI=I1                                                        
  816 CONTINUE                                                                  
      IF(MXI-I) 818,822,818                                                     
  818 DO 820 J=I,IM                                                             
      EX=Q(I,J)  $  Q(I,J)=Q(MXI,J)                                             
  820 Q(MXI,J)=EX                                                               
  822 DIVQ=Q(I,I)  $  DIVQAB=CABS(DIVQ)  $  IF(DIVQAB) 840,830,840              
  830 WRITE(6,835)                                                              
  835 FORMAT(17HODIVQ=0 IN SIMLEQ)                                              
      GO TO 1000                                                                
  840 DIVQ=1.0/DIVQ                                                             
      DO 850 J=I,IM                                                             
  850 Q(I,J)=Q(I,J)*DIVQ                                                        
      I1=I+1                                                                    
      IF(I.EQ.MXROW)  GO TO 880                                                 
      DO 870 K1=I1,MXROW                                                        
      EX=Q(K1,I)                                                                
      DO 870 K2=I,IM                                                            
      Q(K1,K2)=Q(K1,K2)-EX*Q(I,K2)                                              
  870 CONTINUE                                                                  
  880 DO 882 J=1,MXROW                                                          
      JM=J+MXROW                                                                
  882 R(MXROW,J)=Q(MXROW,JM)                                                    
  890 I=MXROW                                                                   
      IF(MXROW-1) 892,910,892                                                   
  892 SUMRI=ZERO                                                                
      DO 895 J=1,MXROW                                                          
      JM=MXROW+J                                                        4/6/75  
      SUMRI=ZERO                                                        4/6/75  
      DO 894 K2=I,MXROW                                                         
      SUMRI=SUMRI+Q(I-1,K2)*R(K2,J)                                             
  894 CONTINUE                                                                  
      L=I-1                                                                     
      R(L,J)=Q(L,JM)-SUMRI                                                      
  895 CONTINUE                                                                  
      I=I-1                                                                     
      IF(I-1) 910,910,892                                                       
  910 I1=I1                                                                     
 1000 MXROW=MXROWM                                                              
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CHEKRM(ICHECK,IRUN,IOK)                                        
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES   ***
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON     FH(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      TYPE COMPLEX TTR,TTI,ZERO                                                 
CCCCCC     *****     COMMON FIELDS USED IN COMMON BY SEVERAL ROUTINES    *******
      COMMON     EXSGRI(140,3)  $  TYPE COMPLEX EXSGRI                          
      COMMON     FC(140,3),FDC(140,3),GC(140,3),GDC(140,3)                      
      COMMON      VCENTR(1200),VCENTI(1200),VNUCLE(1200),               12/1/75 
     1            VCPL1R(1200),VCPL2R(1200),VCPL1I(1200),VCPL2I(1200)   12/1/75 
      COMMON     JLSMAX, VFACN(5),VFACC(5),LAMD(5),M12MX1(20,5),                
     1           M12MX2(5),AMAT(100,5),MCOL(100,5),JLSCMX,NVJLS(15)             
      COMMON      URI(20,20),UDRI(20,20),Q(20,40),R(20,40),QRI(100,5),  4/29/76 
     1            FRI(20,20,4)                                          4/29/76 
      COMPLEX     URI,UDRI,FRI,QRI                                      4/29/76 
      COMPLEX    Q,R,HRI1,HRI2,HRI3,HRIINV,ERI1,ERI2,ERI3                       
      DIMENSION RVOR(25), AMAXI(5) ,TEST(20,1)                                  
      DATA AMAXI/0.,0.,0.,0.,0./                                        ABOY1158
      IOK=0                                                                     
      IF(ICHECK.EQ.1)  GO TO 1                                                  
      ICHECK=1                                                                  
      DO 20 M1=1,MXROW                                                          
   20 RVOR(M1)=0.0                                                              
      RETURN                                                            ABOY1161
    1 DO 2 M=1,MXROW                                                            
      N=NNROW(M)                                                                
      TEST(M,1)=CABS(R(M,1))                                                    
      IF (TEST(M,1).GT.AMAXI(N))  AMAXI(N)=TEST(M,1)                    ABOY1165
    2 CONTINUE                                                          ABOY1166
      DO 3 M=1,MXROW                                                            
      N=NNROW(M)                                                                
      AMAX=ABS(RVOR(M)-TEST(M,1))/AMAXI(N)                              ABOY1169
      IF (AMAX.GT.1.E-4)  GO TO 4                                       ABOY1170
    3 CONTINUE                                                          ABOY1171
      IRUN=IRUN+1                                                               
      IF(IRUN.EQ.3)  GO TO 10                                                   
      GO TO 5                                                           ABOY1174
    4 IRUN=0                                                                    
    5 DO 6 M=1,MXROW                                                            
    6 RVOR(M)=TEST(M,1)                                                 ABOY1177
      RETURN                                                            ABOY1178
   10 IOK=1                                                                     
      RETURN                                                                    
       END                                                               ABOY092
      OVERLAY(JPWKB,3,0)                                                        
      PROGRAM JP1AMP                                                            
      COMMON/LEGENC/LCALTR,MMXTR,RADIAN,PLM20,PLM10,PL(15)                      
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTIN        
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      COMPLEX    TTR,TTI,ZERO                                                   
      COMMON     EXSGRI(420,3)                                          9/23/75 
      COMMON HBAR(2000),OVE(50),ABAMP(200),AMPL(4),FASE(4),CJLS(100),           
     1        CJLSR(100),HBART(500)                                             
      COMPLEX HBAR,OVE,EX,CJLS,EXSGRI,CJLSR,HBART                               
      COMMON CJLS1(100),CJLS2(100),CJLS3(100),CJLS4(100)                10/17/75
      COMPLEX CJLS1,CJLS2,CJLS3,CJLS4                                   10/17/75
      COMMON LLCUMD(10)                                                         
      EQUIVALENCE (JJ,KEXCOM(45)),(K,KEXCOM(46))                                
      DIMENSION FLA(20),FGR(20),FAN(20),FGE(20)                         4/7/75  
      IF(KTLOUT(12).NE.0)  PRINT 101                                    2/7/75  
  101 FORMAT(1H1,20X,10(1H*),10X,21HDWAMP HAS BEEN CALLED,10X,10(1H*))          
      REWIND 7                                                          2/7/75  
      SQRT10=SQRT(10.0)                                                 2/7/75  
      PAI=3.141592653589793                                                     
      TWOPAI=2.*PAI                                                             
      LLCMMX=0                                                                  
      JJ1=JJJMAX                                                        10/18/75
      IF (KTRL(16).EQ.0) GO TO 150                                              
      DO 140 N=1,10                                                             
      FN=N                                                                      
      DEL=FN*0.1                                                                
      DELP1=DEL+1.0                                                             
      DELP2=DEL+2.0                                                             
      DELM1=DEL-1.0                                                             
      DELM2=DEL-2.0                                                             
      FLA(N)=-DEL*DELM1*DELM2/6.0                                               
      FGR(N)=0.5*DELP1*DELM1*DELM2                                              
      FAN(N)=-.5*DELP1*DEL*DELM2                                                
  140 FGE(N)=DELP1*DEL*DELM1/6.0                                                
      JJINT=KTRL(16)+1                                                          
      JJ1=JJINT-1                                                               
      JJM10=JJINT-10                                                            
      KEXCOM(14)=KEXCOM(14)-10                                                  
      JJJMAX=JJJMAX-10                                                          
  150 DO 155 II=1,IICPLE                                                2/7/75  
      LLTW=IIREAD(II)*2                                                 4/9/75  
      LLCUMD(II)=LLCMMX                                                 2/7/75  
      LLCMMX=LLCMMX+(LLTW/2)+1                                          2/7/75  
  155 CONTINUE                                                                  
      LEP1MX=KEXCOM(14)                                                 2/7/75  
      NHBMAX=LEP1MX*LLCMMX                                                      
      DO 710 NHB=1,NHBMAX                                                       
  710 HBAR(NHB)=ZERO                                                            
      JMIN=1  $  JMAX=JJJMAX  $  IF(KTRL(9)) 601,605,601                2/7/75  
  601 JMIN=KEXCOM(5)                                                    2/7/75  
      JMAX=KEXCOM(6)                                                    2/7/75  
  605 DO 800   JJ=JMIN,JMAX                                             2/7/75  
      DO 790 K=1,2                                                      2/7/75  
      IF(KTRL(9)-1)  608,606,608                                        2/7/75  
  606 IF(K-KEXCOM(7)) 790,608,790                                       2/7/75  
  608 CALL NLJJJK  $  IF(MXROW.EQ.0)  GO TO 790                         2/7/75  
      LITW=LLROW(1)  $  LITR=LITW/2                                     2/7/75  
      MXROWE=MXROW                                                              
      MXROWI=1                                                          2/7/75  
      IF (JJ.GT.JJ1) GO TO 210                                                  
      READ(7) (CJLS(M1),M1=1,MXROWE)                                            
      IF(JJ.NE.JJM10)  GO TO 290                                        10/18/75
      DO 205 M1=1,MXROW                                                         
  205 CJLS1(M1)=CJLS(M1)                                                        
      GO TO 290                                                                 
  210 IF (JJ.NE.JJINT) GO TO 220                                                
      READ(7) (CJLS2(M1),M1=1,MXROW)                                            
      READ(7) (CJLS3(M1),M1=1,MXROW)                                            
      READ(7) (CJLS4(M1),M1=1,MXROW)                                            
      DO 215 M1=1,MXROW                                                         
  215 CJLS(M1)=CJLS2(M1)                                                        
      GO TO 290                                                                 
  220 J1=JJ-JJINT                                                               
       J1M=MOD(J1,10)                                                           
      IF (J1M.NE.0) GO TO 230                                                   
      DO 225 M1=1,MXROW                                                         
      CJLS1(M1)=CJLS2(M1)                                                       
      CJLS2(M1)=CJLS3(M1)                                                       
      CJLS3(M1)=CJLS4(M1)                                                       
  225 CJLS(M1)=CJLS2(M1)                                                        
      READ(7) (CJLS4(M1),M1=1,MXROW)                                            
      GO TO 290                                                                 
  230 I=J1M                                                                     
      DO 235 M1=1,MXROW                                                         
  235 CJLS(M1)=FLA(I)*CJLS1(M1)+FGR(I)*CJLS2(M1)                                
     1        +FAN(I)*CJLS3(M1)+FGE(I)*CJLS4(M1)                                
  290 IF (KTLOUT(12).EQ.0) GO TO 620                                            
      WRITE (6,610)  JJ,K                                               2/5/75  
  610 FORMAT(/* JJ,K=*,2I4)                                             2/5/75  
      WRITE (6,611)  (CJLS(M1),M1=1,MXROW)                              2/5/75  
  611 FORMAT(*  CJLS=*,4(2E11.3,3X))                                    2/5/75  
  620 CONTINUE                                                          2/5/75  
      DO 770 NLJE=1,MXROWE                                                      
      LETW=LLROW(NLJE)                                                          
      LETR=LETW/2                                                               
      HAT2=(LETW+1)*(LITW+1)                                            4/9/75  
      II=NNROW(NLJE)   $   IETR=IIREAD(II)                              4/10/75 
      LFOMTW=IIREAD(II)*2                                               4/9/75  
      HAT1=SQRT(HAT2/FLOAT(LFOMTW+1))                                   4/9/75  
      K1=LFOMTW-IABS(LITW-LETW)                                                 
      K2=LITW+LETW-LFOMTW                                                       
      K3=MIN0(K1,K2)                                                            
      IF(K3.GE.0) GO TO 719                                                     
      GO TO 760                                                                 
  719 EX=CJLS(NLJE)/WN(II)                                              2/7/75  
      EX=EX*SQRT10                                                      2/7/75  
      LLCUM=LLCUMD(II)                                                  2/7/75  
      LETR1=LETR+1  $  EX=EX*(EXSGRI(LETR1,II)**2)                      2/7/75  
      MLE1MX=LFOMTW/2+1                                                         
      DO 750 MLE1=1,MLE1MX                                                      
      LLCUM=LLCUM+1                                                             
      MLETW=2*(MLE1-1)                                                          
      IA=LITW                                                                   
      IB=LETW                                                                   
      IC=LFOMTW                                                                 
      ID=0                                                                      
      IE=MLETW                                                                  
      IF=MLETW                                                                  
      CALL CLEB                                                                 
      C4=RAC                                                                    
      IF(C4.EQ.0.  ) GO TO 750                                                  
      C5=(-1.0)**(LITR-IETR+MLETW/2)                                    4/13/75 
      F1=C4*C5*HAT1                                                     4/9/75  
      N1=LETR*LLCMMX+LLCUM                                              2/7/75  
      HBAR(N1)=HBAR(N1)+F1*EX                                                   
  750 CONTINUE                                                                  
  760 CONTINUE                                                                  
  770 CONTINUE                                                                  
  790 CONTINUE                                                          2/7/75  
  800 CONTINUE                                                                  
      IK=IICPLE                                                         2/8/75  
      MLE1MX=IIREAD(IK)/2+1                                             2/8/75  
      N1MAX=MLE1MX*LEP1MX                                               2/8/75  
      IF(KEXCOM(12).EQ.0) GO TO 816                                             
      LLCMX2=2*LLCMMX                                                           
      LE1MX=LEP1MX-1                                                            
      DO 815 LEP1=1,LEP1MX                                                      
      IF(MOD(LEP1,2).NE.0) GO TO 815                                            
      IF(LEP1.EQ.2) GO TO 805                                                   
      IF(LEP1.EQ.LE1MX) GO TO 805                                               
      N0BAS=(LEP1-4)*LLCMMX                                                     
      NMAX=4                                                                    
      GO TO 808                                                                 
  805 IF(LEP1.EQ.2) N0BAS=0                                                     
      IF(LEP1.EQ.LE1MX) N0BAS=(LEP1MX-3)*LLCMMX                                 
      NMAX=2                                                                    
  808 DO 814 M1=1,LLCMMX                                                        
      N0=N0BAS+M1                                                               
      N1=N0+LLCMMX*(NMAX-1)                                                     
      DO 811 N=1,NMAX                                                           
      FASE(N)=0.                                                                
      NT=N0+(N-1)*LLCMX2                                                        
      HREAL=HBAR(NT)                                                            
      HIMAG=-TTI*HBAR(NT)                                                       
      AMPL(N)= SQRT(HREAL**2+HIMAG**2)                                          
      IF(AMPL(N).EQ.0.0  ) GO TO 811                                            
      PHASE= ATAN2(HIMAG,HREAL)                                                 
      IF(PHASE.LT.0.0  ) PHASE=PHASE+TWOPAI                                     
      FASE(N)=PHASE                                                             
      IF(N.EQ.1) GO TO 811                                                      
  809 IF(FASE(N).GE.FASE(N-1)) GO TO 811                                        
      FASE(N)=FASE(N)+TWOPAI                                                    
      GO TO 809                                                                 
  811 CONTINUE                                                                  
      IF(NMAX.EQ.4) GO TO 812                                                   
      AMP1=(AMPL(1)+AMPL(2))/2.0                                                
      FAS1=(FASE(1)+FASE(2))/2.0                                                
      GO TO 813                                                                 
  812 AMP1=.0625  *(-(AMPL(1)+AMPL(4))+9.  *(AMPL(2)+AMPL(3)))                  
      FAS1=.0625  *(-(FASE(1)+FASE(4))+9.  *(FASE(2)+FASE(3)))                  
  813 HBAR(N1)=AMP1*( COS(FAS1)+TTI* SIN(FAS1))                                 
  814 CONTINUE                                                                  
  815 CONTINUE                                                                  
  816 IF(KTLOUT(12)-1)  900,823,817                                     2/7/75  
  817 WRITE(6,818)                                                              
  818 FORMAT(16H OVLAP4-818.HBAR)                                               
      DO 820  LE1=1,LEP1MX                                                      
      N1MIN=(LE1-1)*LLCMMX+1                                                    
      N1MAX=N1MIN+LLCMMX-1                                                      
      WRITE(6,819) (HBAR(N1),N1=N1MIN,N1MAX)                                    
  819 FORMAT(5(E11.3,E11.3,2X))                                                 
  820 CONTINUE                                                                  
  823 LLC=0                                                             2/7/75  
      DO 850 I1=1,IICPLE                                                        
      MLE1MX=IIREAD(I1)+1                                               4/10/75 
      DO 850 MLE1=1,MLE1MX  $  LLC=LLC+1                                2/7/75  
      DO 830 LE1=1,LEP1MX                                                       
      N1=(LE1-1)*LLCMMX+LLC                                                     
      EX=HBAR(N1)                                                       4/9/75  
      EXR=EX                                                                    
      EXI=-TTI*EX                                                               
      ABAMP(LE1)=EXR*EXR+EXI*EXI                                                
  830 CONTINUE                                                                  
      WRITE(6,840) (ABAMP(LE1),LE1=1,LEP1MX)                            2/7/75  
  840 FORMAT(1X,10E13.5)                                                        
  850 CONTINUE                                                                  
  900 CONTINUE                                                                  
      CALL JP1CRS                                                               
      END                                                                       
      SUBROUTINE JP1CRS                                                         
      COMMON/LEGENC/LCALTR,MMXTR,RADIAN,PLM20,PLM10,PL(15)                      
CCCCCC      ******   COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTIN        
      COMMON      FACLOG(1300),RAC,IA,IB,IC,ID,IE,IF,L9(9),U9           12/1/75 
      COMMON      KTRL(30),KEXCOM(50),EXTCOM(50),KTLOUT(30)                     
      COMMON      F(142),G(142),FD(140),GD(140),ETA,SIGMAZ,RHOMX,RD,LMAX        
      COMMON      ANGLER(100),AR(70),AI(70),BR(70),BI(70),CE(10),DR(10),        
     1            ECM(10),                                                      
     2            IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),         
     3            NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),        
     4            WN(10),WNINI(10),WC(10)                                       
      COMMON      ISTRTW,IICPLE,INTYPE,INTMAX,IIXCAL,IIXPLT,IIPCAL,             
     1            IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,               
     2            AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT,         
     3            PMAS,RMAS,RZERO,RZEROC,RZEROS,RZEROW,RZROSP,                  
     4            RMAX,RBAR,SGMAR,TMAS,VSX,VSF,VSO,WSX,WSF,XMAX,XBAR,           
     5            XMES1,XMES2,WNUNIT,TTR,TTI,ZERO                               
      COMPLEX    TTR,TTI,ZERO,EXSGRI                                            
      COMMON     EXSGRI(420,3)                                          9/23/75 
      COMMON HBAR(2000),XAMP(90,40),CROSEX(90)   ,VAPEX(90),                    
     1        PM(90,20),PMST(2,90),LLOUT(7),FCRI(100),SGRUTH(100)               
      COMPLEX    HBAR,XAMP,EX,SUM,HBFAC,FCRI,ERI                                
      COMMON     GCROS(90),PCROS(90,5),LLCUMD(20)                       4/9/75  
      COMMON     CROSOT(90,5)                                           4/14/75 
      DIMENSION  HBART(2000),ACOEFF(90,5)                                       
      DIMENSION XX(100,2),YY(100,2),NDATA(2),XLINE(5),YLINE(5),         12/23/75
     1AREA( 13,301),ISYMBL(2),YSCALE(301),NDAIN(2)                      12/23/75
      COMPLEX HRI1,HRI2,HRI3,HRI4,HBART                                         
      EQUIVALENCE (XAMP(1,1),HBART(1)),(ACOEFF(1,1),HBAR(1))                    
      DATA NEFR,NNR,NZR,IINC,IEXIT / 3HEFR,3HNR ,3HZR ,3HIN ,3HOUT /            
      NTHEMI=KEXCOM(15)                                                 9/10/75K
      KTRL15=KTLOUT(15)     $     KTRL16=1                              2/7/75  
      SQRT10=SQRT(10.0)                                                 2/7/75  
      WN1=WN(1)                                                         2/7/75  
      ETA=CE(1)                                                         2/7/75  
      NTHETA=NANGLR                                                             
      DO 120 NA = 1,NTHETA                                              2/7/75  
      SN = SIN(0.5*THETA(NA))                                           2/7/75  
      SN=SN*SN                                                          2/7/75  
      FLN=ETA*ALOG(SN)-2.0*SGMAZZ(1)                                    2/7/75  
      FNO=ETA/(2.00*WN1*SN)                                             2/7/75  
      ERI=-FNO*SQRT10*(COS(FLN)-TTI*SIN(FLN))                           2/7/75  
      ER=ERI                                                            2/7/75  
      EI=-TTI*ERI                                                       2/7/75  
      SGRUTH(NA)=ER*ER+EI*EI                                            2/7/75  
      FCRI(NA)=ERI                                                      2/7/75  
  120 CONTINUE                                                          2/7/75  
      LLTWMX=0                                                                  
      LLCMMX=0                                                                  
      DO 155 II=1,IICPLE                                                2/7/75  
      LLCUMD(II)=LLCMMX                                                 4/14/75 
      LLTW=IIREAD(II)*2                                                 4/9/75  
      IF(LLTW.GT.LLTWMX) LLTWMX=LLTW                                            
      LLCMMX=LLCMMX+(LLTW/2)+1                                                  
      LLOUT(II)=LLTW/2                                                  2/7/75  
  155 CONTINUE                                                                  
      MMXTR1=LLTWMX/2+1                                                         
      MMXTR=MMXTR1-1                                                            
      LEP1MX=KEXCOM(14)+1                                               4/22/75 
      NBARSU=LLCMMX*LLTWMX                                              4/30/75 
      NBARMX=LEP1MX*LLCMMX                                              4/9/75  
      NITEMX=1                                                          4/9/75  
      IF(KTRL(15).NE.0)  NITEMX=KTRL(15)                                4/14/75 
      DO 2000 NITE1=1,NITEMX                                            4/9/75  
      IF (KTRL(15).EQ.0) GO TO 156                                      10/17/75
      LEP1MX=LEP1MX-1                                                   4/22/75 
      NBARMX=NBARMX-LLCMMX                                              4/22/75 
  156 IF(NITE1.EQ.1)  GO TO 770                                         4/15/75 
CCCCCC                                                                          
CCCCCC     *********     PADE PROCEDURE        *************                    
      DO 750 LEP1=1,LEP1MX                                              4/9/75  
      LETR=LEP1-1                                                       4/15/75 
      FLE2=(LEP1-1)*2                                                   4/9/75  
      FLE3=(FLE2+1.0)*(FLE2-1.0)                                        4/9/75  
      FLE4=(FLE2+1.0)*(FLE2+3.0)                                        4/9/75  
      T3=0.0                                                            4/9/75  
      IF(LEP1.GE.2)  T3=1.0/FLE3                                        4/9/75  
      T4=1.0/FLE4                                                       4/9/75  
      DO 745 II=1,IICPLE                                                4/9/75  
      LLCUM=LLCUMD(II)                                                  4/9/75  
      MLE1MX=IIREAD(II)+1                                               4/9/75  
      FLMP=LEP1-2                                                       4/9/75  
      FLMM=LEP1                                                         4/9/75  
      DO 740 MLE1=1,MLE1MX                                              4/9/75  
      HRI4=0.0                                                          4/15/75 
      LLCUM=LLCUM+1                                                     4/14/75 
      FLMP=FLMP+1.0                                                     4/9/75  
      FLMM=FLMM-1.0                                                     4/9/75  
      N1=LETR*LLCMMX+LLCUM                                              4/9/75  
      IF(MLE1.GT.LEP1)  GO TO 730                                       4/15/75 
      N1P=N1+LLCMMX                                                     4/9/75  
      N1M=N1-LLCMMX                                                     4/9/75  
      HRI2=HBAR(N1)                                                     4/9/75  
      HRI1=0.0                                                          4/9/75  
      HRI3=0.0                                                          4/9/75  
      FLP1=FLMP*FLMM                                                    4/14/75 
      FLP2=(FLMP+1.0)*(FLMM+1.0)                                        4/14/75 
      IF(LEP1.EQ.1)  GO TO 710                                          4/9/75  
      IF(FLP1.LT.0.0)  FLP1=0.0                                         4/14/75 
      HRI1=HBAR(N1M)                                                    4/9/75  
  710 HRI3=HBAR(N1P)                                                    4/9/75  
  720 HRI4=HRI2-SQRT(FLP1*T3  )*HRI1   -SQRT(FLP2*T4)*HRI3              4/14/75 
  730 HBART(N1)=HRI4                                                    4/9/75  
  740 CONTINUE                                                          4/9/75  
  745 CONTINUE                                                          4/9/75  
  750 CONTINUE                                                          4/9/75  
      DO 760 N1=1,NBARMX                                                4/9/75  
  760 HBAR(N1)=HBART(N1)                                                4/9/75  
  770 IF(KTLOUT(13).EQ.0) GO TO 790                                     4/14/75 
      NBART=NBARMX-NBARSU                                               4/30/75 
      WRITE(6,4100) NITE1                                               4/14/75 
 4100 FORMAT(1H0,*HBAR FOR ITE1=*,I2)                                   4/14/75 
      WRITE(6,4200) (HBAR(N1),N1=1,NBART)                               4/30/75 
 4200 FORMAT(1H ,3X,10E11.3)                                            4/14/75 
  790 CONTINUE                                                          4/14/75 
CCCCCC                                                                  4/14/75 
CCCCCC   ***********CALCULATION OF XAMP*****************                4/14/75 
CCCCCC                                                                  4/14/75 
      NTHE=NTHETA                                                       2/7/75  
      DO 805 NA=1,NTHE                                                          
      DO 805 LLC=1,LLCMMX                                                       
  805 XAMP(NA,LLC)=ZERO                                                         
      LEP1TX=LEP1MX                                                     10/17/75
      IF(KTRL(15).EQ.0)  GO TO 806                                      10/17/75
      LEP1TX=LEP1MX-LLTWMX                                              4/30/75 
  806 DO 850 LEP1=1,LEP1TX                                              4/30/75 
      IF(LEP1.GT.KEXCOM(14))  GO TO 850                                 2/7/75  
      LETR=LEP1-1                                                               
      LETW=LETR*2                                                       4/9/75  
      HAT1=SQRT(FLOAT(LETW+1))                                                  
      LCALTR=LETR                                                               
      IF(LETR-1) 808,811,813                                                    
  808 DO 809 IAN=1,NTHE                                                         
      PMST(1,IAN)=1.                                                            
  809 PM(IAN,1)=1.                                                              
      GO TO 825                                                                 
  811 DO 812 IAN=1,NTHE                                                         
      TH=THETA(IAN)                                                     2/7/75  
      PM(IAN,1)= COS(TH)                                                        
      PMST(2,IAN)=PM(IAN,1)                                                     
  812 PM(IAN,2)= SIN(TH)                                                        
      GO TO 825                                                                 
  813 DO 820 IAN=1,NTHE                                                         
      RADIAN=THETA(IAN)                                                 2/7/75  
      PLM20=PMST(1,IAN)                                                         
      PLM10=PMST(2,IAN)                                                         
      CALL LEGNDR                                                               
      DO 815 M=1,MMXTR1                                                         
  815 PM(IAN,M)=PL(M)                                                           
      PMST(1,IAN)=PMST(2,IAN)                                                   
      PMST(2,IAN)=PL(1)                                                         
  820 CONTINUE                                                                  
  825 LLCUM=0                                                                   
      NHBBAS=(LEP1-1)*LLCMMX                                                    
      DO 840 II=1,IICPLE                                                2/7/75  
      MLE1MX=IIREAD(II)+1                                               4/9/75  
      DO 835 MLE1=1,MLE1MX                                                      
      LLCUM=LLCUM+1                                                             
      NHB=NHBBAS+LLCUM                                                          
      MLETAB=IABS((MLE1-1)*2)                                           4/12/75 
      K1=(LETW-MLETAB)/2+1                                                      
      K2=(LETW+MLETAB)/2+1                                                      
      G1= EXP(0.5  *(FACLOG(K1)-FACLOG(K2)))                                    
      HBFAC=HBAR(NHB)*G1*   HAT1                                        4/9/75  
      IF(LEP1.LT.MLE1) GO TO 835                                                
      DO 830 NA=1,NTHE                                                          
      XAMP(NA,LLCUM)=XAMP(NA,LLCUM)+HBFAC*PM(NA,MLE1)                           
  830 CONTINUE                                                                  
  835 CONTINUE                                                                  
  840 CONTINUE                                                                  
  850 CONTINUE                                                                  
      IF(NITE1.EQ.1) GO TO 1001                                         4/14/75 
      NITE2=NITE1-1                                                     4/14/75 
      DO 852 NA=1,NTHE                                                  4/14/75 
      TH=THETA(NA)                                                      4/14/75 
      T1=(1.0-COS(TH))**NITE2                                           4/14/75 
      T2=1.0/T1                                                         4/14/75 
      DO 851 LLC=1,LLCMMX                                               4/14/75 
  851 XAMP(NA,LLC)=XAMP(NA,LLC)*T2                                      4/14/75 
  852 CONTINUE                                                          4/14/75 
CCCCCC                                                                          
CCCCCC     **********     TEST OF CONVERGENCE       **************              
CCCCCC                                                                          
 1001 ITEST=0                                                           4/14/75 
      LLCUM=0                                                           4/9/75  
      IF(NITE1.EQ.NITEMX)  GO TO 1200                                   4/9/75  
      DO 1100 II=1,IICPLE                                               4/9/75  
      LLTW=IIREAD(II)*2                                                 4/9/75  
      MLE1MX=IIREAD(II)+1                                               4/9/75  
      DO 1002 NA=NTHEMI,NTHE                                                    
 1002 CROSEX(NA)=0.0                                                    4/9/75  
      DO 1010 MLE1=1,MLE1MX                                             4/9/75  
      LLCUM=LLCUM+1                                                             
      DO 1005 NA=NTHEMI,NTHE                                                    
      EX=XAMP(NA,LLCUM)                                                 4/9/75  
      IF(II.EQ.1)  EX=EX+FCRI(NA)                                       4/9/75  
      F1=EX*CONJG(EX)                                                   4/9/75  
      IF(MLE1.NE.1)  F1=F1*2.0                                          4/9/75  
      CROSEX(NA)=CROSEX(NA)+F1                                          4/9/75  
 1005 CONTINUE                                                          4/9/75  
 1010 CONTINUE                                                          4/9/75  
      IF(II.NE.1)  GO TO 1020                                           4/9/75  
      DO 1015 NA=NTHEMI,NTHE                                                    
 1015 GCROS(NA)=CROSEX(NA)                                              4/9/75  
      GO TO 1100                                                        4/9/75  
 1020 AMAX=0.0                                                          4/9/75  
      DO 1030 NA=NTHEMI,NTHE                                                    
      P1=CROSEX(NA)/GCROS (NA)                                                  
      IF(NITE1.EQ.1)  GO TO 1025                                        4/9/75  
      P2=PCROS(NA,II)                                                   4/9/75  
      R1=ABS(P1/P2-1.0)                                                 4/9/75  
      IF(R1.GT.AMAX)  AMAX=R1                                           4/9/75  
 1025 PCROS(NA,II)=P1                                                   4/9/75  
 1030 CONTINUE                                                          4/9/75  
      IF(AMAX.LE.1.0E-03)  GO TO 1100                                   4/9/75  
      ITEST=1                                                           4/9/75  
      WRITE (6,4000) AMAX,II,NITE1                                      4/9/75  
 4000 FORMAT(1H0,* AMAX=*,E10.3,*  FOR II,NITE1=*,2I4)                  4/9/75  
 1100 CONTINUE                                                          4/9/75  
      IF(KTLOUT(13).EQ.1)  GO TO 1200                                   5/3/75  
      IF(NITE1.EQ.1)  GO TO 2000                                        4/14/75 
      IF(ITEST.EQ.0)  GO TO 1200                                        4/9/75  
      GO TO 2000                                                        4/9/75  
 1200 CONTINUE                                                          4/9/75  
CCCCCC                                                                          
CCCCCC     **********     CALCULATION OF CROSS SECTION     **********           
CCCCCC                                                                          
      LLCUM=0                                                                   
      IICPL1=IICPLE+1                                                   4/14/75 
      DO 880 II=1,IICPLE                                                4/14/75 
      MLE1MX=IIREAD(II)+1                                               4/12/75 
      II1=II+1                                                          4/15/75 
      DO 855 NA=1,NTHE                                                          
  855 CROSEX(NA)=0.0                                                    4/9/75  
      DO 860 MLE1=1,MLE1MX                                                      
      LLCUM=LLCUM+1                                                             
      DO 857 NA=1,NTHE                                                          
      EX=XAMP(NA,LLCUM)                                                         
      IF(II.EQ.1)  EX=EX+FCRI(NA)                                       2/7/75  
      F1=EX* CONJG(EX)                                                          
      IF(MLE1.NE.1) F1=2.  *F1                                                  
      CROSEX(NA)=CROSEX(NA)+F1                                          4/9/75  
  857 CONTINUE                                                                  
  860 CONTINUE                                                                  
      DO 870 NA=1,NTHE                                                  4/14/75 
      CROSOT(NA,II1)=CROSEX(NA)                                         4/14/75 
      IF(II.EQ.1)  CROSOT(NA,1)=CROSEX(NA)/SGRUTH(NA)                   4/14/75 
  870 CONTINUE                                                          4/14/75 
  880 CONTINUE                                                          4/14/75 
CCCCCC                                                                          
CCCCCC     **********   OUTPUT OF CROSS SECTIONS     **************             
CCCCCC                                                                          
      WRITE (6,962) ELAB                                                4/14/75 
  962 FORMAT(1H1,/1X,*    COUPLED CHANNEL CROSS SECTION FOR  ELAB=*,    2/7/75  
     1        F7.3,4H-MEV)                                              2/7/75  
      WRITE (6,982) (IIREAD(N),N=1,IICPLE)                              4/14/75 
  982 FORMAT(/2X,5HANGLE,3X,8HSIG/RUTH,5X,2HI=,7(I1,8X,2HI=))           4/14/75 
      IICPL1 =IICPLE+1                                                  2/7/75  
      DO 985 NA=1,NTHETA                                                        
      IF(MOD(NA-1,5).EQ.0.AND.NA.NE.1) WRITE(6,983)                             
  983 FORMAT( )                                                                 
      PRINT 984,ANGLER(NA),  (CROSOT(NA,I),I=1,IICPL1)                          
  984 FORMAT(1X,F6.2,9E11.3)                                            2/7/75  
  985 CONTINUE                                                                  
 2000 CONTINUE                                                          4/9/75  
CCCCCC                                                                          
CCCCCC     *****   BETWEEN 2000 - 2500 IS TEMPORARAL ***********                
      IF(IIREAD(2).NE.2) GO TO 2500                                             
      WRITE (6,2210)                                                            
 2210 FORMAT(1H1,*   COEFFICIENTS OF ANGULAR CORR  */1H )                       
      SQRT1=SQRT(1.5)                                                           
      DO 2250 NA=1,NTHE                                                         
      HRI4=-0.5*(XAMP(NA,4)-SQRT1*XAMP(NA,2))                                   
      HRI2=-(SQRT1*XAMP(NA,4)+0.5*XAMP(NA,2))                                   
      HRI1=HRI4+TTI*XAMP(NA,3)                                                  
      HRI3=HRI4-TTI*XAMP(NA,3)                                                  
      HRI1=SQRT1*HRI1  $  HR1=HRI1  $  HI1=AIMAG(HRI1)                          
      HRI2=-HRI2       $  HR2=HRI2  $  HI2=AIMAG(HRI2)                          
      HRI3=SQRT1*HRI3  $  HR3=HRI3  $  HI3=AIMAG(HRI3)                          
      TH1=ATAN2(HI1,HR1)  $  AB1=CABS(HRI1)                                     
      TH2=ATAN2(HI2,HR2)  $  AB2=CABS(HRI2)                                     
      TH3=ATAN2(HI3,HR3)  $  AB3=CABS(HRI3)                                     
      HRI4=AB1*CEXP(-TTI*(TH1-TH2))+AB3*CEXP(TTI*(TH3-TH2))                     
      HR4=HRI4  $  HI4=AIMAG(HRI4)                                              
      TH4=ATAN2(HI4,HR4)  $  AB4=CABS(HRI4)                                     
      T1=AB1**2+AB3**2                                                          
      T2=AB2**2                                                                 
      A=T1+T2                                                                   
      FN=T1*0.6666667+T2                                                        
      B=2.0*AB2*AB4/A                                                           
      C=2.0*AB1*AB3/A                                                           
      THB=0.5*TH4                                                               
      THC=(TH3-TH1)*0.25                                                        
      A=(A/FN)*0.099471841                                                      
      ACOEFF(NA,1)=CROSOT(NA,3)*A                                               
      ACOEFF(NA,2)=B                                                            
      ACOEFF(NA,3)=C                                                            
      ACOEFF(NA,4)=THB                                                          
      ACOEFF(NA,5)=THC                                                          
      WRITE (6,2242) (ACOEFF(NA,N),N=1,5),HRI1,HRI2,HRI3                        
 2242 FORMAT(1H ,5X,5E10.2,2X,6E10.2)                                           
 2250 CONTINUE                                                                  
 2500 CONTINUE                                                                  
      REWIND 8                                                                  
      WRITE (8) NTHE,((ACOEFF(NA,N),N=1,5),NA=1,NTHE)                           
      IF (KTLOUT(20).EQ.0) GO TO 1000                                   5/4/76  
CCCCCC                                                                          
CCCCCC       *****   PLOT    ******                                             
CCCCCC                                                                          
      PRINT 6300                                                        12/30/75
 6300 FORMAT(1H1)                                                       12/30/75
      READ(5,5000) NDATA(1),NDATA(2)                                    12/23/75
      ND1=NDATA(1)                                                      12/23/75
      ND2=NDATA(2)                                                      12/23/75
      DO 6000 J=1,ND1                                                   12/23/75
      READ(5,5100) XX(J,1),YY(J,1)                                      12/23/75
 6000 CONTINUE                                                          12/23/75
 5000 FORMAT(2I5)                                                       12/23/75
 5100 FORMAT(2F7.3)                                                     12/23/75
      DEV=0.0                                                                   
      DO 6100 K=1,ND2                                                   12/23/75
      XX(K,2)=ANGLER(K)                                                 12/23/75
      YY(K,2)=CROSOT(K,1)                                               12/23/75
      DEV=DEV+((YY(K ,2)-YY(K ,1))/YY(K ,2))**2                                 
 6100 CONTINUE                                                          12/23/75
      NDMAX=100                                                         12/23/75
      ISYMBL(1)=1H*                                                     12/23/75
      ISYMBL(2)=1HO                                                     12/23/75
      NF=2                                                              12/23/75
      MX=0                                                              12/23/75
      MY=1                                                              12/23/75
      YLINE(1)=1.0                                                      12/23/75
      NLS=50                                                            12/23/75
      NCL=100                                                           12/23/75
      MM=1                                                              12/23/75
      LL=2                                                              12/23/75
      PRINT 6200,DEV                                                    12/23/75
 6200 FORMAT(*ANGULAR DISTRIBUTION OF ELASTIC SCATTERING*,E12.3)                
      CALL PLOTT(XX,YY,NDATA,NDMAX,ISYMBL,NF,XLINE,MX,YLINE,MY,NLS,     12/23/75
     1NCL,MM,LL,AREA,YSCALE)                                            12/23/75
      PRINT 6300                                                        12/30/7P
      READ(5,5000) NDAIN(1),NDAIN(2)                                    12/30/75
      NDA1=NDAIN(1)                                                     12/30/75
      NDA2=NDAIN(2)                                                     12/30/75
      DO 7000 JI=1,NDA1                                                 12/30/75
      READ(5,5100) XX(JI,1),YY(JI,1)                                    12/30/75
 7000 CONTINUE                                                          12/30/75
      DEV=0.0                                                                   
      DO 7100 KI=1,NDA2                                                 12/30/75
      XX(KI,2)=ANGLER(KI)                                               12/30/75
      YY(KI,2)=CROSOT(KI,3)                                             12/30/75
      DEV=DEV+((YY(KI,2)-YY(KI,1))/YY(KI,2))**2                                 
 7100 CONTINUE                                                          12/30/75
      MY=2                                                              12/30/75
      YLINE(1)=10.0                                                     12/30/75
      YLINE(2)=1.0                                                      12/30/75
      PRINT 7200,DEV                                                    12/30/75
 7200 FORMAT(*ANGULAR DISTRIBUTION OF INELASTIC SCATTERING*,E12.3)              
      CALL PLOTT(XX,YY,NDAIN,NDMAX,ISYMBL,NF,XLINE,MX,YLINE,MY,NLS,     12/23/75
     1NCL,MM,LL,AREA,YSCALE)                                            12/23/75
 1000 RETURN                                                                    
      END                                                                       
      SUBROUTINE LEGNDR                                                         
      COMMON/LEGENC/LCALTR,MMXTR,RADIAN,PLM20,PLM10,PL(15)                      
      FL=LCALTR                                                                 
      FLM1=FL-1                                                                 
      FLP1=FL+1                                                                 
      TWLM1=FL+FLM1                                                             
      TWLP1=FL+FLP1                                                             
      CO= COS(RADIAN)                                                           
      COSAB= ABS(CO)                                                            
      IF( ABS(COSAB-1.  ).LT.1.E-7) GO TO 221                                   
      SI=1.  / SIN(RADIAN)                                                      
      CT=2.  *CO*SI                                                             
      PL0  =(TWLM1*CO*PLM10-FLM1*PLM20)/FL                                      
      PLP10=(TWLP1*CO*PL0  -FL  *PLM10)/FLP1                                    
      PL(1)=PL0                                                                 
      PL(2)=FLP1*SI*(CO*PL0-PLP10)                                              
      MMAX=MIN0(LCALTR,MMXTR)+1                                                 
      IF(MMAX.LT.3) GO TO 200                                                   
      FM=0.                                                                     
      DO 150 M=3,MMAX                                                           
      FMP1=FM+1.                                                                
      PL(M)=CT*FMP1*PL(M-1)-(FL-FM)*(FL+FMP1)*PL(M-2)                           
      FM=FMP1                                                                   
  150 CONTINUE                                                                  
      GO TO 200                                                                 
  221 K1=MOD(LCALTR,2)                                                          
      PL(1)=CO**K1                                                              
      DO 225 M=2,MMAX                                                           
  225 PL(M)=0.                                                                  
  200 RETURN                                                                    
      END                                                                       
