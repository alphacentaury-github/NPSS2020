      PROGRAM OPTXSQ
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                       
      COMMON/COUCC/ FC(200),GC(200),FCP(200),GCP(200),SIG(200),ETA,
     1              SIGMAZ,RD,Z,KTOUT7,L5                            
      COMMON/POTPA/PARAM(39),VREAL(850),VREPC(850),WVOL(850),WSUR1(850),       
     1             WSUR2(850),WTOT(850),WF,AF,RF                               
      COMMON/DWAVE/ TMST,PMST,TZ,PZ,XMES,LMAX,MMAX,NXMX                        
      COMMON/PLMCC/ FACLOG(500),PK(200,200,1)                                  
      COMMON/SGM/   SGMT5(400),SGMEXP(400),DSGMEX(400),THETAD(400)             
      COMMON/CONSTANT/HBAR,HBRSQ,CONST,TTI,FAC,FAC1,FAC2,NONX,NRCUT,
     1                WNUNIT        
      COMMON/INPUT/ E,ELAB,Q,H,WN,CFUNIT,RMST,RCUT,THMIN,THMAX,THINC,
     1              TZPZ,ALI,DLI,ALR,DLR,RLI,RLR,NANGLR,NRCMX,NMX2,NMXM2
      COMMON/SRH/ AKOFF,B(10,300),CHIT,DY(39),HA,PARAMI(39),FCTR,
     1            SGMST5(400),FLAMDA,CHIN(6),KTRL5,JDATA,KPLT(10),
     2            IPA(9),JMAX,KW,KM,NDIM,NOVAR,KTR(24),KOFF,KPRNT(10),
     3            MANGLR,JJXCAL,NCUT                                           
      COMMON/FCRS/ FCROS,FCROSL,FCROSL2                                        
      COMMON/BASIC/ ECM
      COMMON/ARRAY/ SGMEXP1(400),THETAD1(400)
      COMPLEX*16  TTI,SMATX,DISTORWAVE              
      COMMON/WAVE/  DISTORWAVE(400,3000),ABSDISTOR(200),SQDISTOR(200)
	COMMON/TEST/  CROSSRATIO,EXPDATA
      COMMON/SMAT/ SMATX(400),ABSMATL(400),PHASEL(400),TCOER(400)
	COMMON/DIFF/  DIFR,DIF,FUSTL,DRSTL,RECROSL
      DIMENSION  DIRL(200),RTH(200),FLL(200),PARDIR(200),EXPDATA(200),
     1           ANGFUS(200),PARFUS(200),REL(200),CROSSRATIO(200),
     2           DIF(200),RFUS(200),RDIR(200),REDD(200),RDELSIG(200),
     3           DIFR(200),DIFCR(200),FUSTL(200),DRSTL(200),RECROSL(200)                                  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C           KTRLD(1)=1 ECM IS USED                                             
C           KTRLD(5)=1 DIFFERENTIAL ELASTIC CROSS SECTION IS CALCULATED        
C           KTRLD(10)=1 CHI-SQUARE FITTING IS DONE.                             
C           KTRLD(13)=0 DIFF/RUTHERFORD ARE USED IN SGMT5(NTH)          
C                    =1 DIFFERENTIAL X-SECTIONS USED IN SGMT5(NTH) 
C           KTR(1)   =1 RSO=RO, ASO=AO                                          
C           KTR(2)   =1 RI=RIS=RO, AI=AIS=AO                                    
C           KTR(9)   =0 NO OUTPUT OF      DETAILS OF CHI-SQUARE FITTING         
C           KTR(10)  =0 NO OUTPUT OF MORE DETAILS OF CHI-SQUARE FITTING         
C           KTLOUT(1)   PARTIAL WAVE FUSION CROSS SECTION                       
C           KTLOUT(2)   DISTORTED WAVE, SQUARE OF DISTORTED WAVE TIMES W        
C           KTLOUT(3)   ELASTIC DIFFERENTIAL CROSS SECTIONS                     
C           KTLOUT(4)   C-MATRIX                                                
C           KTLOUT(5)   S-MATRIX                                                
C           KTLOUT(6)   TRANSMISSION COEFFICIENTS                               
C           KTLOUT(7)   COULOMB WAVE FN AT MATCHING RADIUS                      
C           KTLOUT(8)   LEGENDRE FUNCTION ( L=0,LMAX ) AT THETA                 
C           KTLOUT(9)   RADIAL WAVE AROUND MATCHING RADIUS                      
C           KTLOUT(10)  GRAPHS OF POTENTIALS                                    
C           KTLOUT(11)  FIVE GRAPHS OF S-MATRIX                                 
C           KTLOUT(12)  THREE GRAPHS OF DISTORTED WAVES                         
C           KTLOUT(13)  PLOT DIFF ELASTIC CROSS SECTION                         
C           KTLOUT(20)  USED FOR POTENTIAL OUTPUT                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C           TARGET MASS AND CHARGE, PROJECTILE MASS AND CHARGE.                 
C           ELABI=STARTING ENERGY. (ECM OR ELAB DEPENDING ON KTRLD(1))          
C           XMES =MESH SIZE.                                                    
C           THETA=ANGLE WHERE DIFF. ELASTIC CROSS SECTION IS CALCULATED         
C           LMAX =NUMBER OF PARTIAL WAVES. LESS THAN 100                        
C           NXMX =NUMBER OF MESH POINTS.   LESS THAN 200                        
C           MMAX =M VALUE FOR ASSOCIATED LEGENDRE. SET MMAX=0 FOR LEGEND        
C           USUAL OMP VALUES                                                    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCC
           OPEN(5,FILE='optxsq.DAT',STATUS='OLD')                                                                    
           OPEN(6,FILE='optxsq.OUT',STATUS='unknown')                                                                     
CCCCCC
C     ***************     READ INITIAL DATA     ***************             
10    FORMAT(24I3)                                           
11    FORMAT(14I5)                                       
15    FORMAT(10F7.4)                                         
16    FORMAT(10F9.5)                                                 
17    FORMAT(18F4.1)
18    FORMAT(10F8.4)                                           
CCCCCC
      READ(5,10) (KTRLD(N),N=1,24)                                              
      READ(5,10) (KTLOUT(N),N=1,24)                                             
      READ(5,15) TMST,TZ,PMST,PZ                           
      READ(5,15) ELABI,RCUT,XMES,THMIN,THMAX,THINC                
      READ(5,11) LMAX,MMAX,NXMX,NANGLR                    
      READ(5,15) VO,WO,WS1,AO,AI,AIS1,RO,RI,RIS1,RC                             
      READ(5,15) WS2,AIS2,RIS2,VS1,VV1,VV2  
      IF(KTRLD(10).EQ.0) GO TO 60                                               
      NOTH=NANGLR                                                               
      MANGLR=NANGLR                                                             
      READ(5,10) (KTR(N),N=1,24)                                                
      READ(5,11) NOVAR,(IPA(I),I=1,NOVAR)                                       
      READ(5,16) FKM,AKOFF,FCTR,FLAMDA                                          
      READ(5,17) (DY(N), N=1,18)
C     TAKE CARE IN CHOOSING FCTR. TOO BIG = NO ACCUR. TOO SMALL = SLOW          
      IF(FCTR.LE.0.0D00) FCTR=0.01D00
      IF(FLAMDA.LE.0.0D00) FLAMDA=0.01D00                               
      KM=FKM+0.1D00                                                      
      READ(5,18) (THETAD(N),N=1,NANGLR)                                         
      READ(5,18) (SGMEXP(N),N=1,NANGLR)                                         
      READ(5,18) (DSGMEX(N),N=1,NANGLR) 
   60 CONTINUE
CCCCCC                                        
C     *******************      CONSTANTS      ****************                  
      PI=3.141592653D0 
      HBAR=197.3285851D0                                             
      CONST=931.501626D0                                              
      WNUNIT=DSQRT(2.0D0*CONST)/HBAR                                   
      HBRSQ =HBAR*HBAR
CCCCCC
      IF (KTRLD(5).EQ.0) GO TO 80                      
      FACLOG(1)=0.0D00                                             
      FACLOG(2)=0.0D00                                                  
      FON=1.D00                                                         
      DO 70 N=3,400                                                         
      FON=FON+1.D00                                                  	
   70 FACLOG(N)=FACLOG(N-1)+DLOG(FON)                                           
   80 CONTINUE                                                  
CCCCCC
      NDIM=18                                                                   
      PARAM(1)=VO                                                               
      PARAM(2)=WO                                                               
      PARAM(3)=WS1                                                              
      PARAM(4)=WS2                                                              
      PARAM(5)=RO                                                               
      PARAM(6)=RI                                                               
      PARAM(7)=RIS1                                                             
      PARAM(8)=RIS2                                                             
      PARAM(9)=RCUT                                                             
      PARAM(10)=AO                                                              
      PARAM(11)=AI                                                              
      PARAM(12)=AIS1                                                            
      PARAM(13)=AIS2                                                            
      PARAM(14)=RC                                                              
      PARAM(15)=VS1                                                            
      PARAM(16)=VV1                                         
      PARAM(17)=VV2                                                     
      JDATA=NOVAR                                                               
      JMAX=NANGLR+KTRLD(11)                 
CCCCCCCCCCCCCCCCCCCC   BASIC CONSTANT     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RMST=TMST*PMST/(TMST+PMST)                                                
      TZPZ=TZ*PZ                                                                
      ELAB=ELABI                                                                
      ECM=ELAB*TMST/(TMST+PMST)                                                 
CCCCCC
CCCCCC
C     ***************     WRITE INPUT DATA     ****************                 
      WRITE(6,99)                                                               
      WRITE(6,100)                                                              
      WRITE(6,101) (KTRLD(N),N=1,14)                                            
      WRITE(6,102) (KTLOUT(N),N=1,14)                                           
      WRITE(6,103) TMST,TZ,PMST,PZ                           
      WRITE(6,104)                                                              
      WRITE(6,105) VO,AO,RO                                                     
      WRITE(6,106) VV1,WO,AI,RI                           
      WRITE(6,107) VS1,WS1,AIS1,RIS1                                 
      WRITE(6,1071)VV2,WS2,AIS2,RIS2                                   
      WRITE(6,108) RC                                                           
      WRITE(6,109) LMAX,MMAX,NXMX,XMES                                          
      IF (KTRLD(1).EQ.1) GO TO 90                                               
      WRITE(6,110) ELABI                                                        
      GO TO 91                                                                  
   90 WRITE(6,111) ELABI                                                        
   91 CONTINUE
                                                            
CCCCCC 
      TMPM=TMST**.33333333333D00+PMST**.333333333333D00               
      IF(PMST.LT.4.0D00) TMPM=TMST**.3333333333333D00                  
      WRITE(6,113) TMPM                                                         
      IF (KTRLD(5).EQ.0) GO TO 95                                             
      WRITE(6,115) THMIN,THMAX,THINC,NANGLR                          
      WRITE(6,117) (KTR(N),N=1,24)                                              
      WRITE(6,118) NOVAR,(IPA(I),I=1,NOVAR)                                     
      WRITE(6,119) FKM,AKOFF,FCTR,FLAMDA                                        
      WRITE(6,125) (DY(N), N=1,18)
      WRITE(6,120) (THETAD(N),N=1,NANGLR)                                       
      WRITE(6,121) (SGMEXP(N),N=1,NANGLR)                                       
      WRITE(6,122) (DSGMEX(N),N=1,NANGLR)
   95 CONTINUE                                       
CCCCCC
   99 FORMAT(///,1X,'OM PARAMETERS SEARCH PROGRAM FOR ELASTIC CROSS S
     1ECTIONS BY CHI SQUARE FITTING',//)                              
  100 FORMAT(5X,'       N =  1  2  3  4  5  6  7  8  9 10 11 12 13 14')         
  101 FORMAT(5X,'KTRLD (N)=',14I3)                                              
  102 FORMAT(5X,'KTLOUT(N)=',14I3,//)                                           
  103 FORMAT(1X, 'TARGET     MASS=',F8.2,' CHARGE=',F8.2/                       
     1       1X, 'PROJECTILE MASS=',F8.2,' CHARGE=',F8.2//)                     
  104 FORMAT(1X, 'INITIAL OPTICAL MODEL PARAMETERS',/,                                  
     1       11X,' REAL(DEPTH)   DEPTH     DIFF   RADIUS')             
  105 FORMAT(1X, 'REAL     =          ',3F10.5)                          
  106 FORMAT(1X, 'IMAGINARY=',4F10.5)                                           
  107 FORMAT(1X, 'SUR.IMAG1=',4F10.5)                                           
 1071 FORMAT(1X, 'SUR.IMAG2=',4F10.5)                                           
  108 FORMAT(1X, 'COULOMB  =',30X,F10.5,//)                                        
  109 FORMAT(1X, 'LMAX, MMAX, NXMX, XMES =',3I6,F8.3)                        
  110 FORMAT(1X, 'ELABI       =',F8.2)                                       
  111 FORMAT(1X, 'ECM         =',F8.2)                                       
  113 FORMAT(1X, 'TM1/3+PM1/3 =',F8.3,//)                                   
  115 FORMAT(1X, 'THMIN,THMAX,THINC,NANGLR  =',3F8.2,I5,//)          
  117 FORMAT(1X, 'INPUT DATA FOR CHI-SQUARE FITTING',//,                   
     1' KTR(I)=              ',24I3,/)                                          
  118 FORMAT(' NOVAR,IPA(I=1,NOVAR) = ',14I5)                                 
  119 FORMAT(' FKM,AKOFF,FCTR,FLAMDA= ',4F14.5)                          
  125 FORMAT(' DY(N=1,18)           = ',10F7.2)
  120 FORMAT(' THETAD   (I=1,NANGLR)= ',10F7.2)             
  121 FORMAT(' SGMEXP   (I=1,NANGLR)= ',10F7.2)                   
  122 FORMAT(' DSGMEX   (I=1,NANGLR)= ',10F7.2)                  
CCCCCC
C  ******************************************************************        
      TTI=(0.D00,1.D00)
      NRCUT=(RCUT*TMPM)/(XMES*10000.D0)*10000.D0+0.0001D0              
      NRCMX=NRCUT                                       
      IF (KTRLD(1).EQ.1) THEN
      ECM=ELAB                                               
      ELAB=ECM*(TMST+PMST)/TMST                              
      ENDIF
      WN=WNUNIT*DSQRT(ECM*RMST)                          
      CFUNIT=1.0D00/(RMST*WNUNIT*WNUNIT)                               
      FAC=40.0D00*3.14159265359D00*RMST*WNUNIT**2/(WN**3)        
      FAC1=4.0D00*RMST*WNUNIT**2/WN                      
      FAC2=40.0D00*3.14159265359D00/(WN*WN)                          
      E=ECM                                                                     
      C=RMST                                                                    
      Q=WN                                                                      
      Z1=TZ                                                                     
      Z2=PZ                                                                     
      L4=LMAX+1                                                                 
      NONX=NXMX-KTRLD(3)                                                        
      H=XMES                                                                    
      RM=H*DFLOAT(NXMX)                       
      NMI=1                                                                     
      NMX2=NXMX+2                                                               
      NMXM2=NXMX-2                                                              
      ETA=C*Z2*Z1*HBAR*1.0D00/137.0359896D00*CONST/(Q*HBRSQ)       
      Z=RM*Q                                                                    
      L5=L4+2                                                                   
      RD=H*Q                                                                    
      KTOUT7=KTLOUT(7)                                 
CCCCCC
  305 CALL FLGLCH                                                               
CCCCCC
      WRITE(6,827)ECM,ELAB,WN,ETA                                
  827 FORMAT(//,1H ,'ECM,ELAB,WN,ETA=',4F8.3)                                   
      SIG(1)=SIGMAZ                                                             
      DO 310 L=1,L4                                                             
      SS=L                                                                      
 310  SIG(L+1)=SIG(L)+DATAN(ETA/SS)                            
      DO 340 L=1,L4                                                             
      FCP(L)=Q*FCP(L)                                                           
 340  GCP(L)=Q*GCP(L)                                                           
CCCCCC
      CALL AUTO                                                                 
CCCCCC
      STOP                                                                      
      END 
CCCCCC
CCCCCC
C     ***************     SUBROUTINE AUTO     ***************                   
                                                                                
      SUBROUTINE AUTO                                                           
	IMPLICIT REAL*8(A-H,O-Z)
                                                                                
C   KW=0,1  PARABOLIC APPROX                                                    
C   KW=2    CENTRAL RUN OF SEARCH                                               
C   KW=4    DERIV RUN                                                           
C   KW=3    OUTPUT RUN                                                          
                                                                                
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                       
      COMMON/CALC/ NCALCUL                                                      
      COMMON/SGM/ SGMT5(400),SGMEXP(400),DSGMEX(400),THETAD(400)              
      COMMON/SRH/ AKOFF,B(10,300),CHIT,DY(39),HA,PARAMI(39),FCTR,
     1            SGMST5(400),FLAMDA,CHIN(6),KTRL5,JDATA,KPLT(10),
     2            IPA(9),JMAX,KW,KM,NDIM,NOVAR,KTR(24),KOFF,KPRNT(10),
     3            MANGLR,JJXCAL,NCUT                                           
      COMMON/POTPA/PARAM(39),VREAL(850),VREPC(850),WVOL(850),WSUR1(850),        
     1             WSUR2(850),WTOT(850),WF,AF,RF                                
      DIMENSION   CP(39),C(9,10),DXPR(39),DELX(39),                             
     1            PRPAR(9),PARAMZ (9),SGMZ(400),IPAP(9),Y(39)                
                                                                                
	IF(KTR(10).EQ.0) GO TO 11                                                 
      WRITE(6,10)                                                                  
   11 KOFF=0                                                                    
      DO 15  I=1,NDIM                                                           
   15 CP(I)=PARAM(I)                                                            
      NOVAR1=NOVAR+1                                                            
      NOVAR2=NOVAR-1                                                            
      LC=0                                                                      
      NCUT=1                                                                    
      FLAM=FLAMDA                                                               
      DO 380 I=1,NDIM                                                           
  380 DXPR(I)=DY(I)                                                             
                                                                                
C-------------    CENTRAL GUESS  (KOFF NUMBER OF GUESS)    ------------------   
                                                                                
  400 DO 420 I=1,NDIM                                                           
      Y(I)=0.D00                                                                   
      DELX(I)=0.D00                                                                
  420 PARAM(I)=CP(I)                                                            
      DO 430 J=1,300                                                            
      DO 430 I=1,10                                                             
  430 B(I,J)=0.0D00                                                                
      KW=2                                                                      
      NCALCUL=1                                                                 
                                                                                
      CALL CALCUL                                                        1ST CAL
                                                                                
      KW=0                                                                      
  449 DEL0=CHIT                                                                 
  450 IF(KOFF) 520,520,460                                                      
  460 LC=0                                                                      
      IF(CHIT-CHIZ) 468,468,461                                                 
  461 DEL0=CHIZ                                                                 
      CHIT=CHIZ                                                                 
      DO 462 J=1,JMAX                                                           
  462 SGMST5(J)=SGMZ(J)                                                         
  465 DO 466 I=1,NOVAR                                                          
      IB=IPA(I)                                                                 
  466 CP(IB)=PARAMZ(I)                                                          
C     CALL CLOCK(KLK)                                                           
C     KLK=1
      WRITE(6,467) CHIT,ZM,ZA,ZB                                            
      DO 464 I=1,NOVAR                                                          
      IPAP(I)=IPA(I)                                                            
      IB=IPA(I)                                                                 
  464 PARAMI(I)=CP(IB)                                                          
      WRITE(6,2000) (IPAP(I),PARAMI(I),I=1,NOVAR)                               
  468  CONTINUE                                                                 

  470 IF(CHIT.EQ.0.D00)  GO TO 480                                                 
  500 IF(DABS((PRCHI-CHIT)/CHIT).LE.AKOFF)  GO TO 540                            
      IF(CHIT.GT.PRCHI)  GO TO 550                                              

  510 IF(KOFF-KM) 520,520,620                                                   
  520 DO 530 I=1,NOVAR                                                          
      IB=IPA(I)                                                                 
  530 PRPAR(I )=CP(IB)                                                          
      PRCHI=CHIT                                                                
      GO TO 630                                                                 
  540 LC=3                                                                      
      GO TO 1320                                                                
  550 IF(NCUT-1) 560,560,610                                                    
  560 DO 570 I=1,NOVAR                                                          
      IB=IPA(I)                                                                 
  570 CP(IB)=PRPAR(I)                                                           
  580 KOFF=0                                                                    
      NCUT=2                                                                    
      WRITE(6,590)                                                              
      DO 600 I=1,NDIM                                                           
  600 DY(I)=DY(I)*.2D00                                                            
      GO TO 400                                                                 
C                    " REDO THE CENTAL GUESS ]"   ------------------            
                                                                                
  610 LC=4                                                                      
  531 CHIT=PRCHI                                                                
      DO 532 I=1,NOVAR                                                          
      IB=IPA(I)                                                                 
      PARAM(IB)=PRPAR(I)                                                        
  532 CP(IB)=PRPAR(I)                                                           
      GO TO 1320                                                                
  620 LC=1                                                                      
  621 GO TO 1320                                                                
  630 KW=4                                                                      
      IF(KTR(9).EQ.0) GO TO 671                                                 
      WRITE(6,7621) KOFF                                                           
  671 NCALCUL=2                                                                 
                                                                                
      CALL CALCUL                                                        2ND CAL
                                                                                
      DO 640 I=1,JMAX                                                           
      IF(SGMEXP(I).LE.0.D00)  GO TO 640   
C      WRITE(6,1278) DSGMEX(I)
 1278 FORMAT(' DSGMEX ',F14.4)                                         
      B(NOVAR1,I)=(SGMEXP(I)-SGMST5(I))/DSGMEX(I)                               
  640 CONTINUE                                                                  
      DO 440 I=1,9                                                              
      DO 440 J=1,10                                                             
  440 C(I,J)=0.0D00                                                                
      DO 740 I=1,NOVAR                                                          
      DO 740 J=1,NOVAR1                                                         
      DO 731 K=1,JMAX                                                           
      IF(SGMEXP(K).LE.0.D00)  GO TO 731                                            
      C(I,J)=C(I,J)+B(I,K)*B(J,K)                                               
  731 CONTINUE                                                                  
      IF(KTR(9).EQ.0) GO TO 740                                                 
      WRITE(6,4457) C(I,J),I,J                                                     
  740 CONTINUE                                                                  
      DO 742 ILM=1,NOVAR                                                        
  742 C(ILM,ILM)=C(ILM,ILM)*(1.D00+FLAM)                                           
      FLAM=FLAM*0.5D00                                                             
      IF(NOVAR.GT.1)  GO TO 750                                                 
      IF(C(1,1).EQ.0.D00)  GO TO 743                                               
      Y(1)=C(1,2)/C(1,1)                                                        
      GO TO 880                                                                 
  750 DO 820 NL=1,NOVAR2                                                        
      AC=DABS(C(NL,NL))                                                          
      IC=NL                                                                     
      I1=NL+1                                                                   
      DO 770 I=I1,NOVAR                                                         
      BC=DABS(C(I,NL))                                                           
      IF(AC-BC) 760,770,770                                                     
  760 AC=BC                                                                     
      IC=I                                                                      
  770 CONTINUE                                                                  
      DO 780 K=NL,NOVAR1                                                        
      DK=C(NL,K)                                                                
      C(NL,K)=C(IC,K)                                                           
  780 C(IC,K)=DK                                                                
      IF(C(NL,NL).EQ.0.D00)  GO TO 800                                             
      DO 790 J=I1,NOVAR1                                                        
      C(NL,J)=C(NL,J)/C(NL,NL)                                                  
      DO 785 I=I1,NOVAR                                                         
  785 C(I,J)=C(I,J)-C(I,NL)*C(NL,J)                                             
  790 CONTINUE                                                                  
  820 CONTINUE                                                                  
      IF(C(NOVAR,NOVAR).EQ.0.D00)  GO TO 830                                       
      Y(NOVAR)=C(NOVAR,NOVAR1)/C(NOVAR,NOVAR)                                   
      DO 870 I=2,NOVAR                                                          
      IS=NOVAR-I+2                                                              
      SUML=0.0D00                                                                  
      DO 860 J=IS,NOVAR                                                         
  860 SUML=SUML+C(IS-1,J)*Y(J)                                                  
CX  870 Y(IS-1)=(C(IS-1,NOVAR1)-SUML)/C(IS-1,IS-1)
  870 Y(IS-1)=(C(IS-1,NOVAR1)-SUML)
  880 PHYM=0.D00                                                                   
      DO 910 I=1,NOVAR                                                          
      IB=IPA(I)                                                                 
      DELX(IB)=Y(I)                                                             
      PHYN=Y(I)/DY(IB)                                                          
      PHYM=DMAX1(PHYM,DABS(PHYN))                                                
  910 CONTINUE                                                                  
                                                                                
      IF(PHYM-0.1D00) 980,940,940                                                 
  940 IF(PHYM-0.3D00) 990,950,950                                                 
  950 IF(PHYM-0.5D00) 1000,960,960                                                
  960 IF(PHYM-1.0D00) 1010,970,970                                                
  970 ZA=.05D00/PHYM                                                               
      GO TO 1020                                                                
  980 ZA=1.0D00                                                                    
      GO TO 1020                                                                
  990 ZA=0.5D00                                                                    
      GO TO 1020                                                                
 1000 ZA=0.2D00                                                                    
      GO TO 1020                                                                
 1010 ZA=0.1D00                                                                    
 1020 DO 1030 I=1,NDIM                                                          
 1030 PARAM(I)=CP(I)+ZA*DELX(I)                                                 
      IF(KTR(9).EQ.0) GO TO 739                                                 
      WRITE(6,7654) (Y(I),I=1,NOVAR)                                               
      WRITE(6,7655) (PARAM(I),I=1,NDIM)
      	                                            
  739 KW=1                                                                      
      NCALCUL=3                                                                 
                                                                                
      CALL CALCUL                                                        3RD CAL
                                                                                
      DEL1=CHIT                                                                 
      KRAZY=0                                                                   
C---------------------------------------------------                            
 1031 CHIZ=CHIT                                                                 
      DO 1032 I=1,NOVAR                                                         
      IB=IPA(I)                                                                 
 1032 PARAMZ(I )=PARAM(IB)                                                      
      DO 1033 J=1,JMAX                                                          
 1033 SGMZ(J)=SGMT5(J)                                                          
 1036 IF(KRAZY-1) 1037,1080,1080                                                
 1037 IF(DEL0-DEL1) 1060,1040,1040                                              
 1040 ZB=2.0D00*ZA                                                                 
      DO 1050 I=1,NDIM                                                          
 1050 PARAM(I)=CP(I)+ZB*DELX(I)                                                 

      KW=1                                                                      
      NCALCUL=4                                                                 
                                                                                
      CALL CALCUL                                                       4TH CALC
                                                                                
      DEL2=CHIT                                                                 
      GO TO 1071                                                                
 1060 DEL2=DEL1                                                                 
      ZB=ZA                                                                     
      ZA=0.5D00*ZA                                                                 
      DO 1070 I=1,NDIM                                                          
 1070 PARAM(I)=CP(I)+ZA*DELX(I)                                                 
      KW=1                                                                      
      NCALCUL=5                                                                 
                                                                                
      CALL CALCUL                                                       5TH CALC
                                                                                
      DEL1=CHIT                                                                 
 1071 IF(CHIZ-CHIT) 1080,1080,1072                                              
 1072 KRAZY=1                                                                   
      GO TO 1031                                                                
C------------------------------------------------------                         
                                                                                
 1080 AC=DEL0-2.0D00*DEL1+DEL2                                                     
      IF(AC) 1150,1090,1110                                                     
 1090 WRITE(6,1100)                                                             
      LC=2                                                                      
      GO TO 1320                                                                
 1110 ZM=ZB*((DEL0-DEL2)/(4.0D00*AC)+0.5D00)                                          
      IF(ZM+ZB) 1120,1130,1130                                                  
 1120 ZM=-ZB                                                                    
      GO TO 1180                                                                
 1130 IF(ZM-ZB) 1180,1140,1140                                                  
 1140 ZM=2.0D00*ZB                                                                 
      GO TO 1180                                                                
 1150 IF(DEL2-DEL0) 1160,1170,1170                                              
 1160 ZM=3.0D00*ZA                                                                 
      GO TO 1180                                                                
 1170 ZM=-ZA                                                                    
 1180 DO 1190 I=1,NDIM                                                          
 1190 CP(I)=CP(I)+ZM*DELX(I)                                                    
      KOFF=KOFF+1                                                               
      GO TO 400                                                                 
C                  " REDO THE CENTRAL GUESS ]"       --------------             
                                                                                
 1320 NCUT=1                                                                    
      DO 1370 I=1,NDIM                                                          
      DY(I)=DXPR(I)                                                             
 1370 PARAM(I)=CP(I)                                                            
      KW=3                                                                      
C     KTLOUT(20)=0                                                             
      NCALCUL=6                                                                 
                                                                                
      CALL CALCUL                                                          6TH C
                                                                                
 1440 GO TO (1450,1760,1510,1490),LC                                            
 1450 WRITE(6,1460) KM                                                          
      GO TO 1760                                                                
 1490 WRITE(6,1500) NCUT                                                        
      GO TO 1760                                                                
 1510 WRITE(6,1520) AKOFF                                                       
 1760 GO TO 200                                                                 
  480 KCHECK=2                                                                  
      GO TO 100                                                                 
  743 KCHECK=4                                                                  
      GO TO 100                                                                 
  890 KCHECK=7                                                                  
      GO TO 100                                                                 
  830 KCHECK=6                                                                  
      GO TO 100                                                                 
  800 KCHECK=5                                                                  
  100 WRITE(6,110) KCHECK                                                       
  200 RETURN                                                                    
  467 FORMAT('---BEST IN PARABOLIC APPROX,---','CHIT=',E20.8,2X,'ZM=',        
     1 F8.3,'  ZA=',E12.5,' ZB=',E12.5)                                   
  590 FORMAT(13H0---CUT DX---)                                                  
 1100 FORMAT(42H0DEL0 DEL1 AND DEL2 ARE ON A STRAIGHT LINE)                     
 1460 FORMAT(' SEARCH LIMITED BY KM=',I3,3X,'NEED LARGER FKM?')               
 1500 FORMAT(' SEARCH LIMITED BY NCUT=',I5)                                    
 1520 FORMAT(' SEARCH LIMITED BY AKOFF=',F11.7)                                
  110 FORMAT(' MAIN KCHECK=',I3)                                               
 2000 FORMAT(4(' PARAM(',I2,')=',F11.5))                                        
C     FOLLOWINGS ARE ZUNC FORMATS.                                             
   10 FORMAT('1',' AUTO IS CALLED',///)                                        
 7621 FORMAT(/,'    KOFF= ',I5)                                                
 4457 FORMAT('    C(I,J) ',E12.5,',  WHEN I AND J = ',2I3)                     
 7654 FORMAT(/,'    DELX(I) = ',9E14.5)                                        
 7655 FORMAT(  '    PARAM(I)= ',9F9.4)                                         
      END                                                                      
CCCCCC
CCCCCC
C     ***************      SUBROUTINE CALCUL     ****************              
                                                                               
      SUBROUTINE CALCUL                                                        
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                       
      COMMON/SGM/   SGMT5(400),SGMEXP(400),DSGMEX(400),THETAD(400)             
      COMMON/SRH/ AKOFF,B(10,300),CHIT,DY(39),HA,PARAMI(39),FCTR,
     1            SGMST5(400),FLAMDA,CHIN(6),KTRL5,JDATA,KPLT(10),
     2            IPA(9),JMAX,KW,KM,NDIM,NOVAR,KTR(24),KOFF,KPRNT(10),
     3            MANGLR,JJXCAL,NCUT                                           
      COMMON/POTPA/PARAM(39),VREAL(850),VREPC(850),WVOL(850),
     1             WSUR1(850),WSUR2(850),WTOT(850),WF,AF,RF                               
      COMMON/CALC/ NCALCUL                                                     
      DIMENSION IBP(9)                                                         
                                                                               
      IF(KTR(10).EQ.0) GO TO 23                                                
      WRITE(6,25) NCALCUL                                                         
   23 IF(KW.EQ.4)  GO TO 13                                                     
                                                                                
      CALL CCCTRL                                                               
      
      CALL CHISQ                                                                
                                                                                
      IF(KTR(10).EQ.0) GO TO 24                                                 
      WRITE(6,3377) KW                                                             
   24 IF(KW.NE.2)  GO TO 90                                                     
      DO 12 J=1,JMAX                                                            
   12 SGMST5(J)=SGMT5(J)                                                        
      GO TO 130                                                                 
                                                                                
C     **********          BEGINNING OF WHEN KW=4          *****************     
   13 DO 30 I=1,NOVAR                                                     *     
      KB=IPA(I)                                                           *     
      PRMIB=PARAM(KB)                                                     *     
      DPRM=PRMIB*FCTR                                                     *     
      PARAM(KB)=PRMIB+DPRM                                                *     
                                                                          *     
      CALL CCCTRL                                                         *     
                                                                          *     
      PARAM(KB)=PRMIB                                                     *     
      DO 20 J=1,JMAX                                                      *     
      IF(SGMEXP(J).LE.0.D00)  GO TO 20                                    *     
CX      write(6,*) ipa(i),param(kb)
      IF(KTR(9).EQ.0) GO TO 19                                            *     
 
C      IF(DPRM.EQ.0.0) GOTO 20
   19 B(I,J)=(SGMT5(J)-SGMST5(J))/(DPRM*DSGMEX(J))                        *     
   20 CONTINUE                                                            *     
   30 CONTINUE                                                            *     
      GO TO 190                                                           *     
C     **********          END OF WHEN KW=4 AND END OF CALCUL         ******     
                                                                                
   90 IF(KW-2) 100,130,190                                                      
  100 IF(KW-1) 190,120,130                                                      
C 120 CALL CLOCK(KLK)                                                           
C     KLK=1
  120 CONTINUE
      IF(KTR(9).EQ.0) GO TO 160                                                 
C     WRITE(6,600) KLK                                                          
      GO TO 160                                                                 
C 130 CALL CLOCK(KLK)                                                           
C     KLK=1
  130 CONTINUE
C      WRITE(6,602) KOFF,KLK                                                  
  160 DO 161 I=1,NOVAR                                                          
      IBP(I)=IPA(I)                                                             
      KB=IPA(I)                                                                 
  161 PARAMI(I)=PARAM(KB)                                                       
      WRITE(6,604) CHIT                                                         
      WRITE(6,1000) (IBP(I),PARAMI(I),I=1,NOVAR)                                

  190 RETURN                                                                    
  600 FORMAT(25H0   FOR PARABOLIC APPROX.,4X,A10)                                
  602 FORMAT(/,17H0   CENTRAL GUESS,5X,I5,5X,A10)                                 
  604 FORMAT( '    CHIT=',1E14.5)                                                
 1000 FORMAT(4('   PARAM(',I2,')=',F11.5))                                      
C     FOLLOWINGS ARE ZUNC FORMATS.                                              
   25 FORMAT(////,I5,'TH CALCUL IS CALLED',/)                                   
 3377 FORMAT('    KW IN CALCUL = ',I2)                                          
      END                                                                       
CCCCCC
CCCCCC
C     ***************     SUBROUTINE CHISQ     ***************                  
                                                                                
      SUBROUTINE CHISQ                                                          
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                        
      COMMON/SGM/   SGMT5(400),SGMEXP(400),DSGMEX(400),THETAD(400)              
      COMMON/SRH/ AKOFF,B(10,300),CHIT,DY(39),HA,PARAMI(39),FCTR,
     1            SGMST5(400),FLAMDA,CHIN(6),KTRL5,JDATA,KPLT(10),
     2            IPA(9),JMAX,KW,KM,NDIM,NOVAR,KTR(24),KOFF,KPRNT(10),
     3            MANGLR,JJXCAL,NCUT                                           
      COMMON/POTPA/PARAM(39),VREAL(850),VREPC(850),WVOL(850),WSUR1(850),        
     1             WSUR2(850),WTOT(850),WF,AF,RF                                
      
      IF(JDATA.EQ.0)  GO TO 900                                                 
      LDATA=0                                                                   
      CHIS=0.D00                                                                   
      DO 20 J=1,JMAX                                                            
      IF(SGMEXP(J).LE.0.D00)  GO TO 20                                             
      LDATA=LDATA+1                                                             
      CHIS=CHIS+((SGMT5(J)-SGMEXP(J))/DSGMEX(J))**2                          
   20 CONTINUE                                                                  
c     CHIN=CHIS/DFLOAT(LDATA)                                                 
      CHIT=CHIS/DFLOAT(LDATA)                                                    
  900 RETURN                                                                    
  
      END                                                                       
CCCCCC
CCCCCC
C     ***************     SUBROUTINE CCCTRL     ***************                 
                                                                                
      SUBROUTINE CCCTRL                                                         
	IMPLICIT REAL*8(A-H,O-Z)
CCCCCC ****  COMMON STATEMENTS THAT ARE COMMON TO ALL THE ROUTINES ****         
                                                 
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                       
      COMMON/AA/     FACLOG(500),RAC,IA,IB,IC,ID,IE,IG,L9(10),U9                    
      COMMON/BB/ ANGLER(100),AR(60),AI(60),BC(60),CE(10),DR(10),ECM(10),        
     2           IIMULT(10),IIREAD(10),KPRITR(10),JJROW(30),LLROW(30),          
     3           NNROW(30),QVALUE(10),SGMAZZ(10),THETA(100),VCOUPL(20),         
     4           WNA(10),WNINI(10),WC(10),CSTH(100),SNTH(100)          
      COMMON/INPUT/ E,ELAB,Q,H,WN,CFUNIT,RMST,RCUT,THMIN,THMAX,THINC,
     1              TZPZ,ALI,DLI,ALR,DLR,RLI,RLR,NANGLR,NRCMX,NMX2,NMXM2
      COMMON/INT/ISTRTW,IICPLE,INTYPE,INTMAX,IIXPLT,IIPCAL,                     
CX     1           IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NANGLR,NDFMES,NYMAX,          
     1           IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NDFMES,          
CX     2           AMUPMU,CHARGE,CFUNIT,DFN,DFNS,DFNW,DFNSP,ELAB,ETUNIT          
     2           AMUPMU,CHARGE,DFN,DFNS,DFNW,DFNSP,ETUNIT          

      COMMON/CONSTANT/HBAR,HBRSQ,CONST,TTI,FAC,FAC1,FAC2,NONX,NRCUT,
     1                WNUNIT        
      COMMON/REA/PMAS,RMAS,RMAX,RBAR,SGMAR,VSF,XMAX,XMES1,FB1                                   

      COMMON/DWAVE/ TMST,PMST,TZ,PZ,XMES,LMAX,MMAX,NXMX                         
      COMMON/SRH/ AKOFF,B(10,300),CHIT,DY(39),HA,PARAMI(39),FCTR,
     1            SGMST5(400),FLAMDA,CHIN(6),KTRL5,JDATA,KPLT(10),
     2            IPA(9),JMAX,KW,KM,NDIM,NOVAR,KTR(24),KOFF,KPRNT(10),
     3            MANGLR,JJXCAL,NCUT                                           
      COMMON/POTPA/PARAM(39),VREAL(850),VREPC(850),WVOL(850),WSUR1(850),        
     1             WSUR2(850),WTOT(850),WF,AF,RF                                
      COMPLEX*16 TTI,BC                                                                            
      DIMENSION   WPOT(850)                                

      IF(KW.EQ.3)  GO TO 141                                                    
CX      DO 142 I=1,30                                                             
CX  142 KTLOUT(I)=0                                                               

  141 IF(KTR(1).EQ.0)  GO TO 143                                                
      PARAM(8)=PARAM(5)                                                         
      PARAM(13)=PARAM(10)                                                       
  143 IF(KTR(2).EQ.0)  GO TO 144                                                
      PARAM(6)=PARAM(5)                                                         
      PARAM(7)=PARAM(5)                                                         
      PARAM(11)=PARAM(10)                                                       
      PARAM(12)=PARAM(10)                                                       
  144 DO 3 I=1,6                                                                
    3 WC(I)=PARAM(I+13)                                                         
C      IF(KTR(6)) 7,8,7                                                         
C    7 DO 9 I=22,26,2                                                           
C    9 PARAM(I)=SQRT(PARAM(I-1)*PARAM(20))                                      
C    8 DO 4 I=1,20                                                              
C    4 VCOUPL(I)=PARAM(I+19)                                                    
      IF(KTR(4).EQ.0)  GO TO 6                                                  
      FNCUT=1.D00                                                                  
      IF(NCUT.NE.1)  FNCUT=.2D00                                                   
      DO 5 I=1,39                                                               
    5 DY(I)=PARAM(I)*FNCUT                                                      
    6 DX=XMES                                                                   
      A1=TMST**0.333333333333D00                                                   
      XCPLE=0.D00                                                                  
      DO 151 I=5,8                                                              
      T1=PARAM(I)*A1                                                            
  151 XCPLE=DMAX1(XCPLE,(T1+10.D00*PARAM(I+5)))                                    
      NYMAX=XCPLE/DX                                                            

      CALL POTEN(WPOT)                                                                

      IMTA=1                                                                    
      JMIN=1                                                                    
                                                                                
      CALL OPT(WPOT)                                                                  
                                                                          
      RETURN                                                                    
      END                                                                       
CCCCCC
CCCCCC
C     ***************     SUBROUTINE OPT     ***************                    
                                                                                
      SUBROUTINE OPT(WPOT)                                                            
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                        
      COMMON/COUCC/ FC(200),GC(200),FCP(200),GCP(200),SIG(200),ETA,
     1              SIGMAZ,RD,Z,KTOUT7,L5                            
      COMMON/POTPA/PARAM(39),VREAL(850),VREPC(850),WVOL(850),WSUR1(850),        
     1             WSUR2(850),WTOT(850),WF,AF,RF                                
      COMMON/DWAVE/ TMST,PMST,TZ,PZ,XMES,LMAX,MMAX,NXMX                         
      COMMON/PLMCC/ FACLOG(500),PK(200,200,1)                                   
      COMMON/SGM/   SGMT5(400),SGMEXP(400),DSGMEX(400),THETAD(400)              
      COMMON/CONSTANT/HBAR,HBRSQ,CONST,TTI,FAC,FAC1,FAC2,NONX,NRCUT,
     1                WNUNIT        
      COMMON/INPUT/ E,ELAB,Q,H,WN,CFUNIT,RMST,RCUT,THMIN,THMAX,THINC,
     1              TZPZ,ALI,DLI,ALR,DLR,RLI,RLR,NANGLR,NRCMX,NMX2,NMXM2
      COMMON/FCRS/  FCROS,FCROSL,FCROSL2                                         
	COMMON/DIFF/  DIFR,DIF,FUSTL,DRSTL,RECROSL
      DIMENSION     U(5),CL(200)                          
      DIMENSION     REALSM(200),RIMASM(200)                        
      DIMENSION     DISWT(850),FUST(200),FUSTL(200),DRST(200),DRSTL(200)                                         
      DIMENSION     TCOE(200),RECROS(200),RECROSL(200)                                                  
      DIMENSION     DIFR(200),DIF(200),DIFCR(200)                               
      COMPLEX*16    COR,B1,CAB,D,CON1,DURM,
     1              CL,U,TTI,SMAT,DISAB,DISPH,
     2              ARI,URI1,URI2,URI3,BRI,SMATX,DISWT,CLL,	  
     3              EXPSIG,CPH,FCC,FN,FNPC,VOFR,SMM1
      DIMENSION     WPOT(850),RDIR(200),ANGDIR(200),RATIO(200),DRCM(200)
     1              ,D1(200),SINCM(200),DELTHE(200),THECM(200),RFUS(200)                                
      COMMON/ARRAY/ SGMEXP1(400),THETAD1(400)
      DIMENSION    RTH(200),FLL(200),PARDIR(200),ANGFUS(200),PARFUS(200)
      COMMON/FNN/   F_L(200),FN_L,FNL_L,PARF_L(200)
      COMPLEX*16    F_L,FN_L,FNL_L

	COMMON/TEST/  CROSSRATIO,EXPDATA
      DIMENSION     CROSSRATIO(200),EXPDATA(200)
c      COMPLEX*16    DISTORWAVE              
c      COMMON/WAVE/  DISTORWAVE(400,3000),ABSDISTOR(400),SQDISTOR(400)
      COMMON/SMAT/  SMATX(400),ABSMATL(400),PHASEL(400),TCOER(400)
	PI=3.141592653D00
	FINE=137.0359896D00
      ECM=E                                                                     
	GAMMA=PMST/TMST
      A_1=TZPZ/FINE*HBAR/ECM
	B=A_1/(16.0D00*PI)
      TMPM=TMST**.33333333333D00+PMST**.333333333333D00                               
      TMU=-(Q**2)/E
	YMG=1.D00*DEXP(-90.0D00)                                                     
      YMP2=-.2D00*YMG                                                              
      HSQ=H**2                                                                  
      HSQ12=HSQ/12.0D00                                                            
      RCUT=PARAM(9)                                                             
                                                                                
      DO 1769 L=1,LMAX                                                          
      TCOE(L)=0.0D00                                                               
      FUST(L)=0.0D00                                                               
 1769 DRST(L)=0.0D00                        !13/06/2005                                       
      LM=0                                                                      

      DO 220 LP1=1,LMAX                                                         
      LM=LM+1                                                                   
      LP=LP1-1                                                                  
      FL1=DFLOAT(LP*2+1)                                                         
      RL=LP                                                                   
                                                                                
                          
  190 URI1=0.D00                                                                   
      IF(LP.EQ.1) URI1=YMP2                                                     
      URI2=YMG                                                                  
      IF (KTRLD(6).EQ.0) GO TO 7675                                             
      EXXPR=DEXP((RL-RLR)/DLR)                                                   
      EXXPI=DEXP((RL-RLI)/DLI)                                                   
      FACV =1.0D00+ALR/(1.0D00+EXXPR)                                                 
      FACW =1.0D00+ALI/(1.0D00+EXXPI)                                                 
 7675 R=0.0D00                                                                     
      NMI=1                                                                     
  200 NXM=0  
  
      DO 205 NX=NMI,NMX2                                                        
	R=R+H 
      VOFR=VREPC(NX)+WTOT(NX)*TTI                                               
      IF (KTRLD(6).EQ.0) GO TO 7676                                             
      VOFR =VREAL(NX)*(FACV-1.0D00)+VREPC(NX)+WTOT(NX)*FACW*TTI                    
 7676 ARI=(E-VOFR)*TMU+RL*(RL+1.D00)/(R**2)                                        
      BRI=1.D00-HSQ12*ARI                                                          
      URI3=(2.D00+(HSQ*ARI/BRI))*URI2-URI1                                         
      NXM=NXM+1                                                                 
      DISWT(NXM)=URI2/BRI                                                       
      URI1=URI2                                                                 
      URI2=URI3                                                                 
      IF(NX.LT.NMXM2) GO TO 205                                                 
      NX1=NX-NMXM2+1
      U(NX1)=URI1/BRI                                                           

  205 CONTINUE                                                                  

      DURM=(1.D00/(12.D00*H))*(8.D00*(U(4)-U(2))-U(5)+U(1))                              
      CON1=-DURM/U(3)                                                           
      L=LP1                                                                     
      B1=CON1*FC(L)+FCP(L)                                                      
      CAB=-CON1*DCMPLX(GC(L),FC(L))                                              
      D=-DCMPLX(GCP(L),FCP(L))                                                   
      CL(L)=B1/(CAB+D)                                                          
      SG=SIG(L)                                                                 
      CON1=FC(L)+CL(L)*(GC(L)+TTI*FC(L))                                        
      COR=CDEXP(SG*TTI)*CON1/U(3)                                                
      NXM=0                                                                     
      DO 210 NX=1,NONX
      NXM=NXM+1                                                                 
      DISWT(NXM)=DISWT(NXM)*COR                                                 
      RADIUS=H*DFLOAT(NXM-1)
      DIST=CDABS(DISWT(NXM))
 210   CONTINUE                                                                  
                                                                                
  450 NRC=1                                                                     
  451 DO 1790 NX=1,NXMX                                                        
      IF(KTRLD(7).EQ.0) VI1=WVOL(NX)                                           
      IF(KTRLD(7).EQ.1) VI1=WVOL(NX)+WSUR2(NX)                                  
      IF(KTRLD(7).EQ.2) VI1=WTOT(NX)                                            
      IF(KTRLD(2).EQ.1) VI1=WPOT(NX)                                            
      IF(KTRLD(8).EQ.0) VI2=WSUR1(NX)                                           
      OVE=-DISWT(NX)*DCONJG(DISWT(NX))*VI1*XMES                                  
      OVE2=-DISWT(NX)*DCONJG(DISWT(NX))*VI2*XMES                                  
      TCOE(LP1)=TCOE(LP1)+OVE*FAC1                                              
C      FUST(LP1)=FUST(LP1)+OVE*FAC*FL1             !07/20/01                              
C      FUSTL(LP1)=FUST(LP1)                        !07/20/01
C      DRST(LP1)=DRST(LP1)+OVE2*FAC*FL1            !12/06/2005                             
C      DRSTL(LP1)=DRST(LP1)                        !12/06/2005
C      DISTORWAVE(LP1,NX)=DISWT(NX)                !2003/08/04                                                  
 1790 CONTINUE                                                                  


  470 KTL2=KTLOUT(2)                                                            
      IF(KTL2.EQ.0) KTL2=1                                                      
      NXM=0                                                                     
      NX1=KTL2                                                                  
      NX2=NONX                                                                  
      IF(KTLOUT(2).EQ.0) GO TO 216
      IF(MOD(LP,50).NE.0) GO TO 216                                             
  216 NXM=0                                                                     
      DO 217 NX=KTL2,NONX,KTL2                                                  
      NXM=NXM+KTL2                                                              
      ARI=DISWT(NXM)                                                            
      DISAB=ARI                                                                 
      DISPH=-TTI*ARI                                                            
      DISWT(NXM)=DISAB*CDEXP(TTI*DISPH)                                          
  217 CONTINUE     
  220 CONTINUE                                                                  
      
	FCROSL=0.0D00                                                                
      FCROSL2=0.0D00                                                               
      FCROS=0.0D00                                                                 
	DRCROSL=0.0D00                                                                
      DRCROSL2=0.0D00                                                               
      DRCROS=0.0D00                                                                 
      DO 1795 L=1,LMAX                                                          
      FCROSL=FCROSL+FUST(L)*DFLOAT(L-1)                                          
      FCROSL2=DRCROSL2+FUST(L)*(DFLOAT(L-1)**2)                                   
C     DRCROSL=DRCROSL+DRST(L)*DFLOAT(L-1)         !12/06/2005                                 
C      DRCROSL2=FCROSL2+DRST(L)*(DFLOAT(L-1)**2)   !12/06/2005                                
      FCROS=FCROS+FUST(L)                                                       
 1795 DRCROS=DRCROS+DRST(L)                    !12/06/2005                                     
C      FCROSL=FCROSL/FCROS                      !18/07/01
C      FCROSL2=FCROSL2/FCROS                    !18/07/01
C      DRCROSL=DRCROSL/DRCROS                   !12/06/2005 
C      DRCROSL2=DRCROSL2/DRCROS                 !12/06/2005 
      
C      IF (KTLOUT(1).EQ.0) GO TO 211                                             
C      IF (KTLOUT(20).EQ.0) GO TO 211                                            
C      WRITE(6,123) RCUT                                                         
C      WRITE(6,122) (FUST(L),L=1,LMAX)                                           
  122 FORMAT(1H ,10E12.5)                                                       
C  123 FORMAT(/,1H ,'PARTIAL WAVE FUSION CROSS SECTION AT RCUT=',F8.3)           
C  211 IF(KTLOUT(1).EQ.0) GO TO 230                                              
C      IF(KTLOUT(20).EQ.0) GO TO 230                                             
C      IF(ETA.EQ.0.0D00) GO TO 229                                                  
C      DO 224 L=1,LMAX                                                           
C      IF (TCOE(L).LT.0.5D00) GO TO 224                                             
C      LH=L                                                                      
C  224 CONTINUE                                                                  
C      LHH=LH-1                                                                  
C      HL=DFLOAT(LH)-0.5D00                                                          
C      A=1.0D00+(HL/ETA)**2                                                         
C      SAR=(ETA/WN)*(1.0D00+DSQRT(A))                                                
C      WRITE(6,124) RCUT                                                         
C      WRITE(6,122) (TCOE(L),L=1,LMAX)                                           
C  229 CONTINUE                                                                  
C  124 FORMAT(/,1H ,'FUSION TRANSMISSION COEFFICIENTS AT RCUT=                   
C     1   ',F8.3)                                                                
C                                                                               
  230 IF(KTLOUT(9).LT.2) GO TO 270                                              
      WRITE(6,520)                                                              
  270 CONTINUE                                                                  
      IF(KTLOUT(4).EQ.0) GO TO 362                                              
      WRITE(6,610)                                                              
      DO 360 L=1,LMAX                                                           
      CLL=CL(L)         !CL(L)=DSIN(PHASEC)*(DCOS(PHASEC)+TTI*DSIN(PHASEC))                                                        
      IF (CLL.EQ.(0.0D00,0.0D00)) GO TO 360                                           
      CRR=DREAL(CLL)                                                             
      CII=DIMAG(CLL)                                                            
      PHASEC=DATAN2(CII,CRR)                                                     
      ABCL2=CRR**2+CII**2                                                       
      ABCLL=DSQRT(ABCL2)                                                         
  360 CONTINUE                                                                  

 362  IF(KTLOUT(5).EQ.0) GO TO 371                                              
      WRITE(6,480)                                   
 371  CROSS=0.0D00                                                                 
      ELCROSS=0.0D00
      AREA=PI/(WN*WN) 
      DO 370 L=1,LMAX                                                           
      LL=L-1                                                                    
      FL1=DFLOAT(LL*2+1)                                                         
      SMATX(L)=1.D00+2.D00*TTI*CL(L)                                                  
      SMAT=SMATX(L) 
      TCOER(L)=1.D00-SMATX(L)*DCONJG(SMATX(L))                                       
      RECROS(L)=10.D00*TCOER(L)*FL1*AREA 
      RECROSL(L)=RECROS(L)                                                     
      CROSS=CROSS+RECROS(L)                                                     
clhb  IF(TZPZ.NE.0.0) GO TO 298                                                 
      SMM1=1.0D00-SMAT                                                             
      SMM12=SMM1*DCONJG(SMM1)                                                    
      ELCROSS=ELCROSS+10.0D00*SMM12*FL1*AREA 
  298 X1=DREAL(SMAT)                                                             
      X2=DIMAG(SMAT)                                                            
      REALSM(L)=X1                                                              
      RIMASM(L)=X2                                                              
      PHASE=DATAN2(X2,X1)                                                        
      PHASEL(L)=PHASE                                                           
      ABSM=(X1**2)+(X2**2)                                                      
      ABSMAT=DSQRT(ABSM) 
      ABSMATL(L)=ABSMAT 
      IF (KTLOUT(5).EQ.0) GO TO 370                                             
      WRITE(6,490) L,SMAT,ABSMAT,PHASE,TCOER(L)                                  
  480 FORMAT(4(/),9X,'L',21X,'S MATRIX',23X,'ABSMAT',15X,'PHASE',15X,           
     1       'TCOE',/)                                                          
  490 FORMAT(5X,I5,10X,2E13.6,12X,1E13.6,7X,1E13.6,7X,1E13.6)                                                        
  370 CONTINUE                                                                  
      TOTCROS=CROSS+ELCROSS
                                                                                
C      NAN1=NANGLR+NUMB+1                                                             
C      NAN2=NANGLR+NUMB+2                                                             
C      NAN3=NANGLR+NUMB+3
C      IF(KTRLD(11).EQ.1 .AND. KTRLD(12).EQ.1) SGMT5(NAN1)=FCROS                 
C      IF(KTRLD(11).EQ.1 .AND. KTRLD(12).EQ.2) SGMT5(NAN1)=CROSS                 
C      IF(KTRLD(11).GE.2) SGMT5(NAN1)=FCROS                                      
C      IF(KTRLD(11).GE.2) SGMT5(NAN2)=DRCROS                                      
C      IF(KTRLD(11).GE.2) SGMT5(NAN3)=ELCROSS                                                                          
C  616 IF(KTLOUT(6).EQ.0) GO TO 282                                              
C      IF(KTLOUT(20).EQ.0) GO TO 282                                             
C      IF(ETA.EQ.0.0D00) GO TO 282                                                  
C      DO 284 L=1,LMAX                                                           
C      IF (TCOER(L).LT.0.5D00) GO TO 284                                            
C      LH=L                                                                      
C  284 CONTINUE                                                                  
C      LHH=LH-1                                                                  
C      HL=DFLOAT(LH)-0.5D00                                                          
C      A=1.0D00+(HL/ETA)**2                                                         
C      SAR=(ETA/WN)*(1.0D00+DSQRT(A))                                                
C      WRITE(6,725) SAR,LHH                                                      
C      WRITE(6,724)                                                              
C      WRITE(6,722) (TCOER(L),L=1,LMAX)                                          
C  282 CONTINUE                                                                  
C  722 FORMAT(1H ,10E12.5)                                                       
C  724 FORMAT(/,1H ,'TRANSMISSION COEFFICIENTS        ')                         
C  725 FORMAT(/,1H ,'STRONG ABSORPTION RADIUS IS    ',F8.3,                      
C     1             '     AND L IS', I5)                                         
                                                                                
clhb  IF(TZPZ.NE.0.0) GO TO 299                                                 
      WRITE(6,499)ELCROSS,TOTCROS                                                         
  499 FORMAT(1X,///,' TOTAL ELASTIC & TOTAL CROSS SECTION IS =',2E12.4)                  
C  299 IF(KTLOUT(20).EQ.0) GO TO 4000                                            
C      WRITE(6,497)                                                              
C      WRITE(6,122) (RECROS(L),L=1,LMAX)                                         
C     WRITE(6,491) CROSS                                                           
C  491 FORMAT(1X,//,' TOTAL REACTION CROSS SECTION IS =',E14.6,//)               
C  497 FORMAT(1X,//,' PARTIAL REACTION CROSS SECTIONS',/)                        
C 4000 CONTINUE                                                                  
                                                                                

CCCCCC      DIFFERENTIAL ELASTIC CROSS SECTION      CCCCCC                      
      IF (KTRLD(5).EQ.0) GO TO 9010                                             

      THETA=THMIN-THINC                                                         
      NOTH=(THMAX-THMIN)/THINC + 1.00001D00  
	IF(KTRLD(10).EQ.0) NANGLR=NOTH                                      
      IF(KTRLD(10).EQ.1) NOTH=NANGLR !+NUMB    !2002/07/22
      DO 489 KTH=1,NOTH                                                         
      IF(KTRLD(10).NE.0) GO TO 476                                              
      THETA=THETA+THINC                                                         
      THETAD(KTH)=THETA                                                         
      TH=THETA*3.14159265359D00/180.0D00                                              
      GO TO 477                                                                 
  476 TH=THETAD(KTH)*3.14159265359D00/180.0D00                                        

  477 CALL LEGNDR(TH,LMAX,MMAX)                                                 
      LMAX1=LMAX+1                                                              
      FN=(0.0D00,0.0D00)                                                              
      DO 486 L1=1,LMAX                                                          
      L=L1-1                                                                    
      L21=2*L+1                                                                 
      EXPSIG=CDEXP(2.D00*TTI*SIG(L1))
	F_L(L1)=1.00D00/(2.00D00*TTI*Q)*L21*EXPSIG*
     1        (SMATX(L1)-1.00D00)          !2003/02/13
      FN=FN + L21*EXPSIG*CL(L1)*PK(L1,1,1)                                      
  486 CONTINUE                                                                  
      FN=FN/Q                                                                   
      SIN2=(DSIN(0.5D00*TH))**2                                                     
      FC0=ETA/(2.D00*Q*SIN2)                                                       
      DIFR(KTH)=10.D00*FC0**2                                                      
      CPH=DCMPLX(0.D00,-ETA*DLOG(SIN2)+2.D00*SIGMAZ)                                   
      FCC=FC0*CDEXP(CPH)                                                         
      FNPC=FN-FCC                                                               
      DIF(KTH)=10.D00*FNPC*DCONJG(FNPC)                                             
      IF(TZPZ.EQ.0.0D00) DIFCR(KTH)=DIF(KTH)                                       
      IF(TZPZ.EQ.0.0D00) GO TO 489                                                 
      DIFCR(KTH)=DIF(KTH)/DIFR(KTH)                                             
      
	IF(KTH.LE.NANGLR) THEN                                        
         IF(KTRLD(13).EQ.0) THEN
            SGMT5(KTH)=DIFCR(KTH)                                                     
         ELSE
            SGMT5(KTH)=DIF(KTH)
         ENDIF
      ELSE
            SGMT5(KTH)=0.0D00
	END IF

  489 CONTINUE 

C      IF(NUMB.EQ.0) GO TO 2024

CCCCCCCCCCCCCCCCCC   CALCULATION OF DIRECT REACTION & FUSION CROSS SECTION    CCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

C	DO L=1,LMAX,1
C	DRSTL(L)=RECROSL(L)-FUSTL(L)
C	END DO

C         WRITE(6,103)
C  103    FORMAT(//,1H ,13X,'      THETA',7X,'STR-RADIUS',7X,'ANGFUS',9X,       
C     1          'ANGDIR',8X,'FUS-RATIO',6X,'DIR-RATIO')                  

C	   DO N=1,NANGLR,1
C           SINTH=DSIN(THETAD(N)*PI/180.D00/2.D00)
C	     COSTH=DCOS(THETAD(N)*PI/180.D00/2.D00)
C
C           RTH(N)=A_1/2.D00*(1.D00+(1.D00/SINTH))
C
C	     FLL(N)=WN*DSQRT(RTH(N)*(RTH(N)-A_1))
C
C	     MFLL=FLL(N)
C	     FLLX1=DFLOAT(MFLL)
C	     FLLX2=DFLOAT(MFLL+1)
C
C	     DFLL=(FUSTL(MFLL+1)-FUSTL(MFLL))*(FLL(N)-FLLX1)/(FLLX2-FLLX1)
C	     DDRL=(DRSTL(MFLL+1)-DRSTL(MFLL))*(FLL(N)-FLLX1)/(FLLX2-FLLX1)
C           PARFUS(N)=FUSTL(MFLL)+DFLL
C	     PARDIR(N)=DRSTL(MFLL)+DDRL
C	     ANGFUS(N)=WN*B*(1.D00/(COSTH*SINTH**3))*PARFUS(N)
C	     ANGDIR(N)=WN*B*(1.D00/(COSTH*SINTH**3))*PARDIR(N)
C           RFUS(N)=ANGFUS(N)/DIFR(N)
C           RDIR(N)=ANGDIR(N)/DIFR(N)
C          END DO

CCCCCCCCCCCCCCCCCC   CALCULATION OF DIRECT REACTION & FUSION CROSS SECTION    CCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

 2024 CONTINUE

      WRITE(6,621)                                                              
      WRITE(6,622)(THETAD(K),DIFR(K),DIF(K),DIFCR(K),SGMEXP(K),THETAD(K)        
     1             ,K=1,NANGLR) 
        

  234 FORMAT(/'DW FOR LX=',I3/(4X,10E11.3))                                     
  235 FORMAT(/8H(DW**)*W/(4X,10E11.3))                                          
  520 FORMAT(/,5X,80H THE RADIAL WAVE AT RM-2*H, RM-H, RM, RM+H, AND RM+        
     12*H- AND ITS DERIVATIVE AT RM,2(/))                                       
  610 FORMAT(5(/),5X,29H THE C MATRIX FOR L=0 TO LMAX,2(/))                     
  621 FORMAT(1H ,'  THETA',20X,'RUTHERFORD',13X,'ELASTIC',8X,                   
     1       'DIFFERENTIAL',8X,'SGM EXPERIMENT',9X,'THETA',/)                   
  622 FORMAT(1H ,F7.2,10X,4E20.5,F14.2)                                         
  623 FORMAT('SWY',1H ,F7.2,10X,3E20.5,F14.2)                                         
  626 FORMAT('SWY1 ','FUSION ',1H ,6E15.6,' REACTION')                                         
  624 FORMAT(/,1H ,'     THETA',19X,'CALCULATION',9X,'EXPERIMENT',13X,                   
     1       'ERROR',11X,'THETA')                   

 9010 RETURN                                                                    
      END                                                                       
CCCCCC
CCCCCC
C     ***************     SUBROUTINE POTEN     ***************                  
                                                                                
      SUBROUTINE POTEN(WPOT)                                                          
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                        
      COMMON/POTPA/PARAM(39),VREAL(850),VREPC(850),WVOL(850),WSUR1(850),        
     1             WSUR2(850),WTOT(850),WF,AF,RFR                               
      COMMON/DWAVE/ TMST,PMST,TZ,PZ,XMES,LMAX,MMAX,NXMX                         
      DIMENSION   WPOT(850),VF(850),VF2(850),VD(850)                                
      COMPLEX*16  TTI                                                      
      COMMON/INPUT/ E,ELAB,Q,H,WN,CFUNIT,RMST,RCUT,THMIN,THMAX,THINC,
     1              TZPZ,ALI,DLI,ALR,DLR,RLI,RLR,NANGLR,NRCMX,NMX2,NMXM2
      COMMON/INT/ISTRTW,IICPLE,INTYPE,INTMAX,IIXPLT,IIPCAL,                     
     1           IIPPLT,JJJMAX,MXROW,NXMAX,NXCPLE,NDFMES,          
     2           AMUPMU,CHARGE,DFN,DFNS,DFNW,DFNSP,ETUNIT          
      COMMON/CONSTANT/HBAR,HBRSQ,CONST,TTI,FAC,FAC1,FAC2,NONX,NRCUT,
     1                WNUNIT        
      COMMON/REA/PMAS,RMAS,RMAX,RBAR,SGMAR,VSF,XMAX,XMES1,FB1                                  
                                                                               
      TTI=(0.0D00,1.0D00)                                                             
      VO=PARAM(1)                                                               
      WO=PARAM(2)                                                               
      WS1=PARAM(3)                                                              
      WS2=PARAM(4)                                                              
      ROR=PARAM(5)                                                              
      RIR=PARAM(6)                                                              
      RIS1R=PARAM(7)                                                            
      RIS2R=PARAM(8)                                                            
      RCUT=PARAM(9)                                                             
      AO=PARAM(10)                                                              
      AI=PARAM(11)                                                              
      AIS1=PARAM(12)                                                            
      AIS2=PARAM(13)                                                            
      RC=PARAM(14)                                                              
      VS1=PARAM(15)                                                            
      VR1=PARAM(16)                                                              
      VV2=PARAM(17)                                                              
      NXMN=1                                                                    
      NXMX2=NXMX+2                                                              
      TMPM=TMST**0.3333333333333D00+PMST**0.33333333333333D00                         
      IF(PMST.LT.4.0D00) TMPM=TMST**.3333333333333D00                                 
	RO=ROR* TMPM                                                              
      RI=RIR* TMPM                                                              
	RIS1=RIS1R*TMPM                                                           
      RIS2=RIS2R*TMPM                                                           
      RN=RC*TMPM                                                                
      IF (KTRLD(2).EQ.1) RF=RFR*TMPM                                            
      CH=TZ*PZ                                                                  
CCCC  FOR THE NEUTRON SCATTERING      
      IF(CH.LT.1.0D00) CH=0.0D00
CCCC      
      H=XMES                                                                    
      EE=1.4398650D00                                                              
      VCC=.5D00*CH*EE/RN                                                           
      IK=0                                                                      
      RMIN=H*DFLOAT(NXMN)                                                        
      R=RMIN-H         
	
C      WRITE(7,112)
C 112 FORMAT('        R','         VREAL','       VREPC','      WTOT')
      DO 270 NX=NXMN,NXMX2                                                      
      IK=IK+1                                                                   
      R=R+H                                                                     
      VR=-VO/(1.D00+DEXP((R-RO)/AO))                                                
      T1=1.0D00/(1.D00+DEXP((R-RI)/AI))                                                
      VRV=-VR1*T1                                               
      VIV=-WO*T1                                               
      IF (KTRLD(2).EQ.1) VIF=VIV/(1.D00+DEXP((R-RF)/AF))                            
       EXPA1=DEXP((R-RIS1) /AIS1)                                                
       EXPA2=DEXP((R-RIS2) /AIS2)                                                
       EX1=1.D00/(1.D00+EXPA1)                                                         
       EX2=1.D00/(1.D00+EXPA2)                                                         
       VRS1=-4.D00*VS1*(EX1**2)*EXPA1                                               
       VRV2=-VV2*EX2                                       ! VOLUME PART  2003/08/14                               
       VIS1=-4.D00*WS1*(EX1**2)*EXPA1                                               
       VIS2=-WS2*EX2                                       ! VOLUME PART  2003/08/14                                        
       VIS12=VIS1+VIS2                                                           
       VIVS=VIV+VIS12
	 VEFFECT=0*(1+1)*CFUNIT/R**2                                                            
       IF (R-RN) 230,230,240                                                     
  230 VVC=VCC*(3.D00-((R/RN)**2))                                                  
      GO TO 250                                                                 
  240 VVC=CH*EE/R                                                               
  250 CONTINUE                                                                  
      WPOT(IK)=VIVS                                                             
      IF (KTRLD(2).EQ.1) WPOT(IK)=VIF                                           
      VREAL(IK)=VR+VRV+VRV2+VRS1                                                              
      VF(IK)=VRV                                                         
      VD(IK)=VRS1                                                          
      VREPC(IK)=VR+VRV+VRV2+VRS1+VVC                                                          
      WVOL(IK)=VIV                                                              
      WSUR1(IK)=VIS1                                                            
      WSUR2(IK)=VIS2                                                            
      WTOT(IK)=VIVS                                                             
  270 CONTINUE                                                                  

      IF(KTLOUT(20).EQ.0) RETURN                                                
      WRITE(6,300) H*10.0D00                                                            
      WRITE(6,301) (VREAL(NX),NX=1,NXMX,10)                                        
      WRITE(6,302) H*10.0D00                                                            
      WRITE(6,301) (VREPC(NX),NX=1,NXMX,10)                                        
      WRITE(6,303) H*10.0D00                                                            
      WRITE(6,301) (WTOT(NX),NX=1,NXMX,10)                                         
                                                                                
CCCCCC      TWO-PARAMETER THEORY IS APPLIED       CCCCCC                        
                                                                                
      IF (KTRLD(2).EQ.0) GO TO 280                                              
      WRITE(6,304) H*10D00                                                            
      WRITE(6,301) (WPOT (NX),NX=1,NXMX)                                        
  280 CONTINUE                                                                  
                                                                                
CCCCCC      TWO-PARAMETER THEORY IS APPLIED.     CCCCCC                         
                                                                                
  300 FORMAT(1H ,'REAL WOODS-SAXON POT IN STEPS OF',F7.3)                       
  302 FORMAT(1H ,'REAL VNUC + VCOU POT IN STEPS OF',F7.3)                       
  303 FORMAT(1H ,'IMAGINARY POTENTIAL  IN STEPS OF',F7.3)                       
  304 FORMAT(1H ,'FUSION IMG POTENTIAL IN STEPS OF',F7.3)                       
  301 FORMAT(1H ,10E12.3)                                                       
      RETURN                                                                    
      END                                                                       
CCCCCC
CCCCCC
C     ***************     SUBROUTINE FLGLCH     ***************                 
                                                                                
      SUBROUTINE FLGLCH                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      REAL*8          K,K1,K2,K3,K4,M1,M2,M3,M4                                   
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                       
      COMMON/COUCC/ FC(200),GC(200),FCP(200),GCP(200),SIG(200),EETA,
     1              SIGMAX,RD,RHO,KTOUT7,MAXL                            
CX      COMMON/COUCC/ FC(200),GC(200),FCP(200),GCP(200),MAXL,                          
CX     1              EETA,SIGMAX,RD,RHO,KTOUT7,SIG(200)                     
                                                                                
      STEP=0.0D00                                                                  
      ACCUR=0.0D00                                                                 
      STP=STEP                                                                  
      ACCR=ACCUR                                                                
      ETA=EETA                                                                  
      IF(ETA-10.D00) 280,280,290                                                   
280   ETA2=ETA**2                                                               
CCCCC ETA2A=2.*ETA                                                              
      ETA6=ETA2+16.D00                                                             
      SIGMAO=-(ETA/(12.D00*ETA6))*(1.D00+(ETA2-48.D00)/           
     1 (30.D00*(ETA6**2))+((ETA2-160.D00)*ETA2+1280.D00)/         
     2 (105.D00*(ETA6**4)))-ETA+(ETA/2.D00)*DLOG(ETA6)+
     3 3.5D00*DTAN(.25D00*ETA)-(DTAN(ETA)+DTAN(.5D00*ETA)+
     4 DTAN(ETA/3.D00))                      
      GO TO 300                                                                 
290   EINV1=1.D00/ETA                                                              
      EINV2=EINV1*EINV1                                                         
      EINV3=EINV1*EINV2                                                         
      EINV5=EINV3*EINV2                                                         
      EINV7=EINV5*EINV2                                                         
      EINV9=EINV7*EINV2                                                         
      SIGMAO=.7853981634D00+ETA*DLOG(ETA)-ETA-(.08333333333D00*EINV1+                 
     1 .00277777777D00*EINV3+.00079365079D00*EINV5+.00059523810D00*EINV7               
     2 +.00084175084D00*EINV9)                                                     
300   MODTPI=SIGMAO/6.2831853072D00                                                
      FMODTP=MODTPI                                                             
      SIGMAZ=SIGMAO-6.2831853072D00*FMODTP                                         
      SIGMAX=SIGMAZ                                                             
      MINL = 0                                                                  
      F=0.0D00                                                                     
      FP=0.0D00                                                                    
      PACE=STEP                                                                 
      ACC=ACCUR                                                                 
      IF(PACE.LT.100.0D00)PACE=100.0D00                                               
      IF(ACC.LT.1.0E-15.OR.ACC.GT.1.0E-6) ACC = 1.0E-6                          
      R=RHO                                                                     
      KTR=1                                                                     
      LMAX=MAXL                                                                 
      LMIN1=MINL+1                                                              
      XLL1=DFLOAT(MINL*LMIN1)                                                    
      ETA2=ETA*ETA                                                              
      TURN=ETA+DSQRT(ETA2+XLL1)                                                  
      IF(R.LT.TURN.AND.DABS(ETA).GE.1.0E-6)KTR=-1                                
      KTRP=KTR                                                                  
      GOTO2                                                                     
1     R=TURN                                                                    
      TF=F                                                                      
      TFP=FP                                                                    
      LMAX=MINL                                                                 
      KTRP=1                                                                    
2     ETAR=ETA*R                                                                
      RHO2=R*R                                                                  
      PL=DFLOAT(LMAX+1)                                                          
      PMX=PL+0.5D00                                                                
      FP=ETA/PL+PL/R                                                            
      DK=ETAR*2.00D00                                                              
      DEL=0.0D00                                                                   
      D=0.0D00                                                                     
      F=1.0D00                                                                     
      K=(PL*PL-PL+ETAR)*(2.0D00*PL-1.0D00)                                            
      ABPL=DABS(PL**2+PL+ETAR)                                                   
      IF(ABPL.GT.1.0E-30) GO TO 3                                               
      R=R+1.0E-6                                                                
      GOTO2                                                                     
3     H=(PL*PL+ETA2)*(1.0D00-PL*PL)*RHO2                                           
      K=K+DK+PL*PL*6.0D00                                                          
      D=1.0D00/(D*H+K)                                                             
      DEL=DEL*(D*K-1.0D00)                                                         
      IF(PL.LT.PMX)DEL=-R*(PL*PL+ETA2)*(PL+1.0D00)*D/PL                            
      PL=PL+1.0D00                                                                 
      FP=FP+DEL                                                                 
      IF(D.LT.0.0D00)F=-F                                                          
      IF(PL.GT.20000.D00)GOTO11                                                    
      IF(DABS(DEL/FP).GE.ACC)GOTO3                                               
      FP=F*FP                                                                   
      IF(LMAX.EQ.MINL)GOTO5                                                     
      FC(LMAX+1)=F                                                              
      FCP(LMAX+1)=FP                                                            
      L=LMAX                                                                    
      DO 4 LP=LMIN1,LMAX                                                        
      PL =DFLOAT(L)                                                             
      GC(L+1)=ETA/PL+PL/R                                                       
      GCP(L+1)=DSQRT(ETA2+PL*PL)/PL                                              
      FC(L)=(GC(L+1)*FC(L+1)+FCP(L+1))/GCP(L+1)                                 
      FCP(L)=GC(L+1)*FC(L)-GCP(L+1)*FC(L+1)                                     
4     L=L-1                                                                     
      F=FC(LMIN1)                                                               
      FP=FCP(LMIN1)                                                             
5     IF(KTRP.EQ.-1)GOTO1                                                       
      P=0.0D00                                                                     
      Q=R-ETA                                                                   
      PL=0.0D00                                                                    
      AR=-(ETA2+XLL1)                                                           
      AI=ETA                                                                    
      BR=2.0D00*Q                                                                  
      BI=2.0D00                                                                    
      WI=2.0D00*ETA                                                                
      DR=BR/(BR*BR+BI*BI)                                                       
      DI=-BI/(BR*BR+BI*BI)                                                      
      DP=-(AR*DI+AI*DR)                                                         
      DQ=(AR*DR-AI*DI)                                                          
6     P=P+DP                                                                    
      Q=Q+DQ                                                                    
      PL=PL+2.0D00                                                                 
      AR=AR+PL                                                                  
      AI=AI+WI                                                                  
      BI=BI+2.0D00                                                                 
      D=AR*DR-AI*DI+BR                                                          
      DI=AI*DR+AR*DI+BI                                                         
      T=1.0D00/(D*D+DI*DI)                                                         
      DR=T*D                                                                    
      DI=-T*DI                                                                  
      H=BR*DR-BI*DI-1.00D00                                                        
      K=BI*DR+BR*DI                                                             
      T=DP*H-DQ*K                                                               
      DQ=DP*K+DQ*H                                                              
      DP=T                                                                      
      IF(PL.GT.46000.D00)GOTO11                                                    
      IF(DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q))*ACC)GOTO6                           
      P=P/R                                                                     
      Q=Q/R                                                                     
      G=(FP-P*F)/Q                                                              
      GP=P*G-Q*F                                                                
      W=1.0D00/DSQRT(FP*G-F*GP)                                                     
      G=W*G                                                                     
      GP=W*GP                                                                   
      IF(KTR.EQ.1)GOTO8                                                         
      F=TF                                                                      
      FP=TFP                                                                    
      LMAX=MAXL                                                                 
      IF(RHO.LT.0.2D00*TURN)PACE=999.00D00                                            
      R3=1.0D00/3.0D00                                                                
      H=(RHO-TURN)/(PACE+1.0D00)                                                   
      H2=0.5D00*H                                                                  
      I2=IDINT(PACE+0.001D00)                                                       
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
      RM3=2.0D00*M3                                                                
      M3=RM3                                                                    
      K4=H2*(GP+M3)                                                             
      RH=R+H                                                                    
      S=(ETAH+H2LL/RH)/RH-H2                                                    
      M4=S*(G+K3)                                                               
      G=G+(K1+K2+K3+K2+K4)*R3                                                   
      GP=GP+(M1+M2+M2+M3+M4)*R3
	BBB=65.0D00                                                 
      IF(G.GT.1.0D00*DEXP(BBB).OR.GP.GT.1.0D00*DEXP(BBB))GOTO11                         
      R=RH                                                                      
      I2=I2-1                                                                   
      IF(I2.GE.0)GOTO7                                                          
      W=1.0D00/(FP*G-F*GP)                                                         
8     GC(LMIN1)=G                                                               
      GCP(LMIN1)=GP                                                             
      IF(LMAX.EQ.MINL) GO TO 10                                                 
      IF(LMAX.EQ.MINL) GO TO 10                                                 
      DO 9 L=LMIN1,LMAX                                                         
      T=GC(L+1)                                                                 
      GC(L+1)=(GC(L)*GC(L+1)-GCP(L))/GCP(L+1)                                   
      GCP(L+1)=GC(L)*GCP(L+1)-GC(L+1)*T                                         
      FC(L+1)=W*FC(L+1)                                                         
9     FCP(L+1)=W*FCP(L+1)                                                       
      FC(LMIN1)=FC(LMIN1)*W                                                     
      FCP(LMIN1)=FCP(LMIN1)*W                                                   
      GO TO 107                                                                 
10    FC(LMIN1)=W*F                                                             
      FCP(LMIN1)=W*FP                                                           
      GO TO 107                                                                 
11    W=0.0D00                                                                    
      G=0.0D00                                                                     
      GP=0.0D00                                                                    
      GOTO8                                                                     
107   CONTINUE                                                                  
      IF(KTOUT7.EQ.0) GO TO 330                                                 
      L4=LMAX                                                                   
      WRITE(6,560)                                                              
      WRITE(6,570) (FC(L),L=1,L4)                                               
      WRITE(6,580) (FCP(L),L=1,L4)                                              
      WRITE(6,590) (GC(L),L=1,L4)                                               
      WRITE(6,600) (GCP(L),L=1,L4)                                              
560   FORMAT(//,5X,43H THE COULOMB WAVE FUNCTIONS FOR L=0 TO LMAX,2(/))         
570   FORMAT(6HFC(L)=,5E15.5)                                                   
580   FORMAT(7HFCP(L)=,5E15.5)                                                  
590   FORMAT(6HGC(L)=,5E15.5)                                                   
600   FORMAT(7HGCP(L)=,5E15.5)                                                  
330   CONTINUE                                                                  
      RETURN                                                                    
      END                   
CCCCCC
CCCCCC
C     ***************     SUBROUTINE LEGNDR     ***************                 
                                                                                
      SUBROUTINE LEGNDR(TH,LMAX,MMAX)                                           
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CNTRL/ KTRLD(24),KTLOUT(24)                                        
CCCCC SUBROUTINE FOR D-FUCTION                                                  
CCCCC LEGENDRE FUCTION=PK(K,MP,IN)                                              
      COMMON/PLMCC/ FACLOG(500),PK(200,200,1)                                   
                                                                                
      IND=1                                                                     
      THETA=TH                                                                  
      LMAX1=LMAX+1                                                              
      LMAX2=LMAX+2                                                              
      MMAX1=MMAX+1                                                              
      DO 920 L=1,LMAX2                                                          
      DO 920 M=1,MMAX1                                                          
      PK(L,M,1)=0.0D00                                                             
  920 CONTINUE                                                                  
                                                                                
                                                                                
      IF (THETA.NE.0.0D00) GO TO 950                                               
      DO 940 L=1,LMAX1                                                          
      PK(L,1,1)=1.0D00                                                             
  940 CONTINUE                                                                  
      RETURN                                                                    
  950 CONTINUE                                                                  
                                                                                
      DO 8500 L=1,LMAX1                                                         
      IL=LMAX1-L+1                                                              
      ILT=IL-1                                                                  
      APH=-1.0D00 
      IF (MOD(ILT,2).EQ.0) APH=1.0D00                                              
      IL2=ILT*2+1                                                               
      ASI=DSIN(THETA)/2                                                          
      ASI=DLOG(ASI)                                                             
      DFT=FACLOG(IL2)/2-FACLOG(IL)+ASI*DFLOAT(ILT)
      PK(IL,IL,1)=DEXP(DFT)                                                      
      PK(IL,IL,1)=APH*PK(IL,IL,1)                                               
      PK(IL+1,IL,1)=PK(IL,IL,1)*DCOS(THETA)*DSQRT(DFLOAT(IL2))                     
 8500 CONTINUE                                                                  
                                                                                
      LMAXM1=LMAX-1                                                             
CCCCCC      SET LMAXM1=1, IF YOU WANT ONLY LEGENDRE. YOU CAN SAVE TIME.         
      LMAXM1=1                                                                  
      DO 8550 ML=1,LMAXM1                                                       
      AML=DFLOAT(ML-1)                                                           
      A2=PK(ML+1,ML,1)                                                          
      A1=PK(ML,ML,1)                                                            
      LMN=ML+2                                                                  
      DO 8600 L=LMN,LMAX1                                                       
      LT=L-2                                                                    
      ALT=DFLOAT(LT)                                                             
      B1=DSQRT(ALT**2-AML**2)*A1                                                 
      B2=(2.0D00*ALT+1.0D00)*DCOS(THETA)*A2                                            
      B3=DSQRT((ALT+AML+1.D00)*(ALT+1.D00-AML))                                        
      PK(L,ML,1)=(B2-B1)/B3                                                     
      A1=A2                                                                     
      A2=PK(L,ML,1)                                                             
CCCCCCPRINT 9998,A2,ML,L                                                        
 9998 FORMAT(1H ,'PK=',E12.5,'  ML=',I5,'  L=',I5)                              
 8600 CONTINUE                                                                  
 8550 CONTINUE                                                                  
CCCCCCPRINT 888,KTLOUT(8)                                                       
  888 FORMAT(1H ,'KTLOUT(8)=',I7)                                               
      IF (KTLOUT(8).EQ.0) GO TO 8553                                            
      WRITE(6,662) (PK(I,1,1),I=1,LMAX1)                                        
  662 FORMAT(1H ,10E12.4)                                                       
 8553 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
