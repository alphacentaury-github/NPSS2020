C ABPASATURN-1-FOR-EFR-DWBA. EXACT FINITE RANGE DWBA CALCULATIONS FOR           
C 1   HEAVY-ION INDUCED NUCLEAR REACTION.  TAMURA,T., LOW,K.S.                  
C REF. IN COMP. PHYS. COMMUN. 8 (1974) 349                                      
C PLEASE NOTE: UNDEFINED VARIABLES IN THE PRG. NEED TO BE SET ZERO.             
C      PROGRAM SATURN (INPUT,OUTPUT,PUNCH,TAPE5=INPUT,TAPE6=OUTPUT,              
C     1              TAPE11,TAPE12,TAPE13)                                       
      PROGRAM SATURN                                     
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC,L9(10),U9                   
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/RARB/NATOTL,NBMIN,NBMAX,NBREAD,NBTAIL,MESFCA,MESFCB,               
     1            N1MAX,N2MAX,MGAUS,NGAUS                   
      COMMON/RARB1/XMESA,XMESB,XMES(2)                   
      COMMON/GAUS/ABSCIS(222),WEIGHT(222),WWL1(500),WWL2(500),                  
     1            BSPAR(46),IBSPR(10)                                           
      COMMON USINGL (500),ISTORE(4),STORE(11)                                   
CCCCCC    /CRAC/   APPEARS IN SATURN,CFCCAL,FFRARB,NRFMFC.                      
CCCCCC    /SATN/   APPEARS IN SATURN,MONITR,KERNEL,GKRARB,FRFMFC,CFCCAL,        
CCCCCC                        FFRARB,NRFMFC,INTERP.                             
CCCCCC    /BIND/   APPEARS IN BSAXON,UNCPST.                                    
CCCCCC    /RARB/   APPEARS IN SATURN,MONITR,KERNEL,GKRARB,FRFMFC,FFRARB,        
CCCCCC                        NRFMFC,INTERP.                                    
CCCCCC    /GAUS/   APPEARS IN SATURN,MONITR,KERNEL,GKRARB,CFCCAL,FFRARB,        
CCCCCC                        NRFMFC.                                           
CCCCCC    /ST12/   APPEARS IN KERNEL,GKRARB,CFCCAL.                             
CCCCCC    /GKAB/   APPEARS IN KERNEL,GKRARB.                                    
CCCCCC      IF KTRLST(1)=0 EXACT FINITE RANGE CALCULATION IS MADE               
CCCCCC                  =1 NO-RECOIL APPROXIMATION IS MADE.                     
CCCCCC      IF KEXCST(1)=NON-ZERO LAMIN=KEXCST(1).                              
CCCCCC      IF KEXCST(2)=NON-ZERO LASTEP=KEXCST(2). (IF KEXCST(2)=0,LAST        
CCCCCC      KEXCST(14)=KPUNCH  (SEE INPUT OF KPUNCH)                            
CCCCCC      KEXCST(15)=NGRNMX (SEE MONITR,CFCCAL,FFRARB).                       
C      DATA PLUS,MINS,EVEN,ODD/ 3H(+),3H(-),2H  ,2H/2 /,                         
C     1     NNR,NEFR / 4HNR) ,4HEFR) /                                           
C      INTEGER PLUS,MINS,EVEN,ODD                                                
C      DIMENSION NAT(100)                                                        
C      DATA NEUTRN,NAT / 2H N,                                                   
C     *        2H H,2HHE,2HLI,2HBE,2H B,2H C,2H N,2H O,2H F,2HNE,2HNA,           
C     *        2HMG,2HAL,2HSI,2H P,2H S,2HCL,2HAR,2H K,2HCA,2HSC,2HTI,           
C     *        2H V,2HCR,2HMN,2HFE,2HCO,2HNI,2HCU,2HZN,2HGA,2HGE,2HAS,           
C     *        2HSE,2HBR,2HKR,2HRB,2HSR,2H Y,2HZR,2HNB,2HMO,2HTC,2HRU,           
C     *        2HRH,2HPD,2HAG,2HCD,2HIN,2HSN,2HSB,2HTE,2H I,2HXE,2HCS,           
C     *        2HBA,2HLA,2HCE,2HPR,2HND,2HPR,2HSM,2HEU,2HGD,2HTB,2HDY,           
C     *        2HHO,2HER,2HTM,2HYB,2HLU,2HHF,2HTA,2H W,2HRE,2HOS,2HIR,           
C     *        2HPT,2HAU,2HHG,2HTL,2HPB,2HBI,2HPO,2HAT,2HRN,2HFR,2HRA,           
C     *        2HAC,2HTH,2HPA,2H U,2HNP,2HPU,2HAM,2HCM,2HBK,2HCF,2HES,           
C     *        2HFM/                                                             
CCCCCC
      OPEN(UNIT=5, FILE='saturn.DAT',STATUS='OLD')
      OPEN(UNIT=6, FILE='saturn.OUT',STATUS='UNKNOWN')
CCCCCC                                      
      FACLOG(1)=0.D0                                                              
      FACLOG(2)=0.D0                                                              
      FN=1.D0                                                                     
      DO 5 N=3,500                                                              
      FN=FN+1.D0                                                                  
    5 FACLOG(N)=FACLOG(N-1)+DLOG(FN)                                            
   10 FORMAT(14I5)                                                              
   15 FORMAT(24I3)                                                              
   20 FORMAT(10F7.2)                                                            
CCCCCC  ******               INPUT DATA FOR SATURN                ******        
      READ(5,15) (KTRLST(N),N=1,10),(KEXCST(N),N=1, 6),(KOUTST(N),N=1,8)        
      READ(5,15) MASCA,NZCA,KPCA,ICATW,MASCB,NZCB,KPCB,ICBTW,                   
     1           MASSA,NZSA,KPSA,ISATW,MASSB,NZSB,KPSB,ISBTW                    
      READ(5,15) JLSMAX,(JJTWR(N),LLTWR(N),ISTWR(N),N=1,JLSMAX)                 
      IF(KTRLST(1).EQ.0) READ(5,10) NATOTL,NBMIN,NBMAX,                         
     1                              MESFCA,MESFCB,MGAUS,LDWMAX                  
      IF(KTRLST(1).EQ.1) READ(5,10) NBMIN,NBMAX,MESFCB,MGAUS,KPUNCH             
      READ(5,20) XMESB,XMESA                                                    
      IF(KTRLST(1).EQ.1) KEXCST(14)=KPUNCH                                      
CCCCCC  ******           PREPARATION OF INITIAL OUTPUT            ******        
      K1=KPCA+KPCB+KPSA+KPSB                                                    
      KPCHAG=MOD(K1,2)                                                          
C      NATCA=NAT(NZCA)                                                           
C      NATCB=NAT(NZCB)                                                           
C      NATSA=NAT(NZSA)                                                           
C      NATSB=NAT(NZSB)                                                           
C      ICAWR=ICATW                                                               
C      KCAWR=EVEN                                                                
C      KCAPA=PLUS                                                                
C      K1=MOD(ICATW,2)                                                           
C      IF(K1.EQ.0) ICAWR=ICATW/2                                                 
C      IF(K1.EQ.1) KCAWR=ODD                                                     
C      IF(KPCA.EQ.1) KCAPA=MINS                                                  
C      ICBWR=ICBTW                                                               
C      KCBWR=EVEN                                                                
C      KCBPA=PLUS                                                                
C      K1=MOD(ICBTW,2)                                                           
C      IF(K1.EQ.0) ICBWR=ICBTW/2                                                 
C      IF(K1.EQ.1) KCBWR=ODD                                                     
C      IF(KPCB.EQ.1) KCBPA=MINS                                                  
C      ISAWR=ISATW                                                               
C      KSAWR=EVEN                                                                
C      KSAPA=PLUS                                                                
C      K1=MOD(ISATW,2)                                                           
C      IF(K1.EQ.0) ISAWR=ISATW/2                                                 
C      IF(K1.EQ.1) KSAWR=ODD                                                     
C      IF(KPSA.EQ.1) KSAPA=MINS                                                  
C      ISBWR=ISBTW                                                               
C      KSBWR=EVEN                                                                
C      KSBPA=PLUS                                                                
C      K1=MOD(ISBTW,2)                                                           
C      IF(K1.EQ.0) ISBWR=ISBTW/2                                                 
C      IF(K1.EQ.1) KSBWR=ODD                                                     
C      IF(KPSB.EQ.1) KSBPA=MINS                                                  
CCCCCC  ******                    INITIAL OUTPUT                  ******        
C      NREFR=NEFR                                                                
C      IF(KTRLST(1).EQ.1) NREFR=NNR                                              
C      WRITE(6,115) NREFR                                                        
      IF(KTRLST(1).EQ.0) WRITE(6,115)                                                         
      IF(KTRLST(1).EQ.1) WRITE(6,116)                                                         
  115 FORMAT(1H1,///,55X,20(1H*),///,45X,43HFORM FACTOR CALCULATION ON P        
     1ROGRAM SATURN (,4HEFR),//)                                                    
  116 FORMAT(1H1,///,55X,20(1H*),///,45X,43HFORM FACTOR CALCULATION ON P        
     1ROGRAM SATURN (,3HNR),//)                                                    
      WRITE(6,124) NZCA,MASCA,NZSA,MASSA,NZSB,MASSB,NZCB,MASCB              
  124 FORMAT(42X,5H FOR ,1H(,I3,2H -,I3,4H) ((,I3,2H -,I3,3H),(,
     1      I3,2H -,I3,4H)) (,I3,2H -,I3,17H)   (CHARGE-MASS))                       
      WRITE(6,126) ICATW,KPCA,ISATW,KPSA,                      
     1             ISBTW,KPSB,ICBTW,KPCB                          
  126 FORMAT( 49X,1H(,I1,2H -,I2,1H),5X,1H(,I1,2H -,I2,1H),
     1  4X,1H(,I1,2H -,I2,1H), 5X,1H(,I1,2H -,I2,1H),
     2  4X,27H(2*SPIN - PARITY, 0=+, 1=-)   )       
      WRITE(6,132)(KTRLST(N),N=1,10),(KEXCST(N),N=1,6),(KOUTST(N),N=1,8)        
  132 FORMAT(55X,20(1H*),///,50X,7HKTRLST=,10I3,/,50X,7HKEXCST=,6I3,/,          
     1   50X,7HKOUTST=,8I3)                                                     
      WRITE(6,170)(JJTWR(N),LLTWR(N),ISTWR(N),N=1,JLSMAX)                       
  170 FORMAT(//,30X,12H(2J,2L,2S)= ,5(1H(,I2,1H,,I2,1H,,I2,4H) , ))             
      IF(KTRLST(1).EQ.0) WRITE(6,230) NATOTL,NBMIN,NBMAX,                       
     1                         MESFCA,MESFCB,MGAUS,LDWMAX,XMESA,XMESB           
  230 FORMAT(//,30X,28HNATOTL,NBMIN,NBMAX         =,3I5,                        
     1      /30X,28HMESFCA,MESFCB,MGAUS,LDWMAX =4I5,                            
     2       /,30X,6HXMESA=,F8.3,5X,6HXMESB=,F8.3)                              
      IF(KTRLST(1).EQ.1) WRITE(6,240) NBMIN,NBMAX,MESFCB,MGAUS,KPUNCH,          
     1                   XMESB                                                  
  240 FORMAT(//,30X,34H NBMIN,NBMAX,MESFCB,MGAUS,KPUNCH =,5I5,                  
     1      /,30X,6H XMESB,27X    ,1H=,F8.3)                                    
CCCCCC  ******    CALCULATION OF BOUND STATE WAVE FUNCTIONS       ******        
      DO 635 N=1,2                                                              
      II=JJTWR(1)                                                               
      IF(N.EQ.1) GO TO 604                                                      
      II=ISTWR(1)                                                               
  604 CALL BSAXON                                                               
      XMES(N)=STORE(9)                                                          
      I1=(N-1)*3                                                                
      DO 608 I=1,3                                                              
      I1=I1+1                                                                   
  608 IBSPR(I1)=ISTORE(I)                                                       
      I1=(N-1)*9                                                                
      DO 609 I=1,9                                                              
      I1=I1+1                                                                   
  609 BSPAR(I1)=STORE(I)                                                        
      IF(ISTORE(3).EQ.II) GO TO 611                                             
      WRITE(6,610)                                                              
  610 FORMAT(1X,10(1H*),61H ERROR, CHECK INPUT OF(J,L,S) WITH INPUT OF          
     1JBTRTW IN BSAXON)                                                         
      STOP                                                                      
  611 IF(N.EQ.2) GO TO 621                                                      
      N1MAX=ISTORE(4)                                                           
      DO 612 NX=1,N1MAX                                                         
  612 WWL1(NX) =USINGL(NX)                                                      
      DO 614 JLS=1,JLSMAX                                                       
  614 L1TRTW(JLS)=2*ISTORE(2)                                                   
      GO TO 635                                                                 
  621 N2MAX=ISTORE(4)                                                           
      DO 622 NX=1,N2MAX                                                         
  622 WWL2(NX) =USINGL(NX)                                                      
      DO 624 JLS=1,JLSMAX                                                       
  624 L2TRTW(JLS)=2*ISTORE(2)                                                   
  635 CONTINUE                                                                  
CCCCCC  ******    MONITORING ROUTINE CHECKS FOR PROPER DIMENSION  ******        
CCCCCC  ******    ALLOCATION IN SATURN AND IN MARS.               ******        
      CALL MONITR                                                               
      IF(KTRLST(1).EQ.1) GO TO 720                                              
CCCCCC  ******         CALCULATION OF EFR FORM FACTOR BEGINS      ******        
      CALL KERNEL                                                               
      CALL FRFMFC                                                               
      GO TO 740                                                                 
CCCCCC  ******         CALCULATION OF  NR FORM FACTOR BEGINS      ******        
  720 CALL NRFMFC                                                               
  740 CONTINUE                                                                  
      STOP                                                                   
      END
                                                                      
      SUBROUTINE BSAXON                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/BIND/KTRL2,KTRL8,KTRL10,KEX2,KEX4,KEX40,KEX41,KEX42,               
     1            TMAS,PMAS,ZZT,ZZP,RMAS,ZZ,XMES1,ACURCY,PERCNT,VSX,            
     2            ISTW,NXRA,NXRM,NXRMP1,NXRMP2,NXRMP3,NXRMP4,NXRMP5,            
     3            NODE,KGES,EGESRD,EGES,EGEST,DELGES,FKAPPA,FKAPIN,             
     4            URRMIN,URRMEX,RGDLIN(3),RGDLEX(3)                             
      COMMON URSAVE(500),NODER,LBTR,JBTRTW,NXRAWF,ENEPT,VNEPT,                  
     1       VSOR,DFNR,DFNSO,RZR,RZSO,RZC,XMES2,                                
     2       XMEM(514),VCENTR(514),VSPIN(514),VCOULM(514)                       
C      DATA IQI,IQJ/4HB.E.,4HVSX  /                                              
C      DATA KK,KL,KM,KN / 1HX,4HVCEN,4HVSPI,4HCOUL /                             
CCCCCC      IF KTRL2=0 B.E. IS SEARCHED, =1 VSX IS SEARCHED.                    
CCCCCC      IF KTRL3=1 R=RZ*(A1**(1/3)+A2**(1/3))                               
CCCCCC      IF KTRL4=0 WAVE FUNCTION ITSELF IS STORED IN URSAVE.                
CCCCCC              =1 WAVE FUNCTION TIMES POTENTIAL IS STORED IN URSAVE        
CCCCCC              =2 WAVE FUNCTION TIMES R**(-L) IS STORED IN URSAVE.         
CCCCCC              =3 WAVE FUNCTION TIMES POTENTIAL TIMES R**(-L) IS           
CCCCCC                 STORED.                                                  
CCCCCC      IF KTRL8=1 URSAVE IS PUNCHED OUT.                                   
CCCCCC      KEX2=NXCPL2 IF NON-ZERO.                                            
CCCCCC      KEX4+NXCPL2=NXRA IF KEX4 IS NON-ZERO.                               
      DIMENSION PFORM1(2),PFORM2(2),PFORM3(2),GESMEM(20),NODMEM(10)             
      PAI=3.1415926535D0                                                          
      ITENOD=0                                                                  
      PERCNT=0.2D0                                                                
CCCCCC  ******              INPUT DATA FOR BSAXON                 ******        
      READ(5,51) KTRL2,KTRL3,KTRL4,KTRL8,KEX2,KEX4,NODER,LBTR,ISTW,             
     1           JBTRTW,ITBEMX,KOPOT                                            
      READ(5,53) EGES,TMAS,PMAS,ZZT,ZZP,XMES2,ACURCY,AMUPMU                     
      READ(5,53) VSX,VSOR,DFNR,DFNSO,RZR,RZSO,RZC                               
   51 FORMAT(14I5)                                                              
   53 FORMAT(10F7.4)                                                            
      UNITOK=0.2187066D0  *(1.  -AMUPMU)+.21953760D0  *AMUPMU                       
      EGEST=EGES                                                                
CCCCCC  ******                INITIAL OUTPUT                      ******        
      WRITE(6,61)                                                               
   61 FORMAT(1H1,55X,17H BSAXON IS CALLED)                                      
      WRITE(6,63) KTRL2,KTRL3,KTRL4,KTRL8,KEX2,KEX4                             
   63 FORMAT(//10X,12HKTRL2,3,4,8=,4I3,10H   KEX2,4=,2I5)                       
      WRITE(6,65) PMAS,ZZP,TMAS,ZZT                                             
   65 FORMAT(/10X34HMOTION BETWEEN NUCLEI WITH (A,Z)=(,2F7.2,7H) AND (,         
     1        2F7.2,1H))                                                        
      WRITE(6,67) EGES,LBTR,ISTW,JBTRTW,NODER                                   
   67 FORMAT(10X10HWITH EGES=,F7.3,6H  LTR=,I1,5H  2S=,I1,5H  2J=,I2,           
     1      8H  NODER=,I1)                                                      
      WRITE(6,69)                                                               
   69 FORMAT(/17X,3HVSX,6X,4HVSOR,6X,4HDFNR,5X,5HDFNSO,7X,3HRZR,6X,             
     1     4HRZSO,7X,3HRZC)                                                     
      WRITE(6,71) VSX,VSOR,DFNR,DFNSO,RZR,RZSO,RZC                              
   71 FORMAT(10X,7F10.3)                                                        
      WRITE(6,73) AMUPMU,ACURCY                                                 
   73 FORMAT(/10X,7HAMUPMU=,F4.1,11H  ACCURACY=,F8.6)                           
      RMAS=(TMAS*PMAS)/(TMAS+PMAS)                                              
      ZZ=ZZT*ZZP                                                                
      NODE=NODER                                                                
      FKTRL2=KTRL2                                                              
      XMES1=XMES2*.125D0                                                          
      XBFAC=TMAS**.3333333333                                                   
      IF(KTRL3.EQ.1) XBFAC=XBFAC+PMAS**.3333333333                              
      XBAR=RZR*XBFAC                                                            
      XRA=XBAR+10.0D0  *DFNR                                                      
      NXRAWF=XRA/XMES2                                                          
      IF(KEX2.EQ.0) GO TO 113                                                   
      NXRAWF=KEX2                                                               
  113 IF(KEX4.EQ.0) GO TO 117                                                   
      NXRA=NXRAWF+KEX4                                                          
      GO TO 119                                                                 
  117 NXRA=NXRAWF+30                                                            
  119 NXRM=XBAR/XMES2                                                           
      FNXRM=NXRM                                                                
      NXRMP1=NXRM+1                                                             
      NXRMP2=NXRM+2                                                             
      NXRMP3=NXRM+3                                                             
      NXRMP4=NXRM+4                                                             
      NXRMP5=NXRM+5                                                             
      NXRA12=NXRA+12                                                            
      WRITE(6,135) NXRAWF,NXRA,NXRM,XMES2,XBAR                                  
  135 FORMAT( 10X,7HNXRAWF=,I3,3X,5HNXRA=,I3,3X,5HNXRM=,I3,                     
     1      8H  XMES2=,F7.4,7H  XBAR=,F7.4)                                     
CCCCCC  ******    CONSTRUCTION OF THE BINDING POTENTIAL(S)        ******        
      NTTLMS=NXRA+12                                                            
      XBARSO=XBFAC*RZSO                                                         
      XBARC =XBFAC*RZC                                                          
      VSPFC=2.0D0  *VSOR/DFNSO                                                    
      VCLFC2=1.4398650D0  *ZZ                                                     
      VCLFC1=VCLFC2*0.5D0  /XBARC                                                 
      DO 180 ND=1,4                                                             
      FND=2.0D0  **(ND-1)                                                         
      DX=XMES1*FND                                                              
      IF(ND.NE.1) GO TO 163                                                     
      X=0.0D0                                                                    
      NXIMIN=1                                                                  
      NXIMAX=8                                                                  
      GO TO 165                                                                 
  163 NXIMIN=NXIMAX+1                                                           
      NXIMAX=NXIMAX+4                                                           
      IF(ND.EQ.4) NXIMAX=NTTLMS+2                                               
  165 DO 175 NX=NXIMIN,NXIMAX                                                   
      X=X+DX                                                                    
      XMEM(NX)=X                                                                
      PFORM1(1)= DEXP((X-XBAR  )/DFNR )                                          
      PFORM1(2)= DEXP((X-XBARSO)/DFNSO)                                          
      DO 167 N=1,2                                                              
      PFORM2(N)=1.0D0 /(1.0D0  +PFORM1(N))                                         
      PFORM3(N)=PFORM1(N)*PFORM2(N)*PFORM2(N)                                   
  167 CONTINUE                                                                  
      VCENTR(NX)=-VSX*PFORM2(1)                                                 
      VSPIN(NX)=(VSPFC*PFORM3(2))/X                                             
      IF(X.GT.XBARC) GO TO 172                                                  
      VCOULM(NX)=VCLFC1*(3.0D0  -((X/XBARC)**2))                                  
      GO TO 175                                                                 
  172 VCOULM(NX)=VCLFC2/X                                                       
  175 CONTINUE                                                                  
  180 CONTINUE                                                                  
      IF(KOPOT.EQ.0) GO TO 190                                                  
      NXWRIT=NTTLMS+2                                                           
      WRITE(6,181)                                                            
      WRITE(6,185)  (XMEM  (NX),NX=1,20)                                        
      WRITE(6,185)  (XMEM  (NX),NX=21,NXWRIT,KOPOT)                             
      WRITE(6,182)                                                            
      WRITE(6,185)  (VCENTR(NX),NX=1,NXWRIT,KOPOT)                              
      WRITE(6,183)                                                            
      WRITE(6,185)  (VSPIN (NX),NX=1,NXWRIT,KOPOT)                              
      WRITE(6,184)                                                            
      WRITE(6,185)  (VCOULM(NX),NX=1,NXWRIT,KOPOT)                              
  181 FORMAT(1X,'X')                                                             
  182 FORMAT(1X,'VCENTRAL')                                                             
  183 FORMAT(1X,'VSPIN')                                                             
  184 FORMAT(1X,'VCOULOMB')                                                             
  185 FORMAT(10E12.4)                                                           
CCCCCC  ******    BEGINNING OF ITERATIVE DETERMINATION OF BINDING ******        
CCCCCC  ******    ENERGY OR DEPTH OF THE SAXON POTENTIAL          ******        
  190 KEX41=0                                                                   
      KEX42=0                                                                   
  205 ITEBE=0                                                                   
  210 ITEBE=ITEBE+1                                                             
      GESMEM(ITEBE)=EGEST*(1.0D0  -FKTRL2)+VSX*FKTRL2                             
CCCCCC  ******    K40=1 AND 2 ARE INTEGRATION OF THE INTERNAL     ******        
CCCCCC  ******    SOLUTION OUTWARDS AND THE EXTERNAL SOLUTION     ******        
CCCCCC  ******    INWARDS RESPECTIVELY TO BE MATCHED SMOOTHLY     ******        
CCCCCC  ******    AT THE NUCLEAR RADIUS.                          ******        
      DO 400 K40=1,2                                                            
      KEX40=K40-1                                                               
      IF(KTRL2.EQ.1) GO TO 213                                                  
      EGESDL=0.005D0  *EGEST                                                      
      EGES=EGEST-2.0D0  *EGESDL                                                   
      GO TO 214                                                                 
  213 EGES=EGEST                                                                
  214 DO 390 KGES=1,3                                                           
      FKGES=KGES                                                                
      IF(KTRL2.EQ.1) GO TO 232                                                  
      EGES=EGES+EGESDL                                                          
      GO TO 235                                                                 
  232 VCOREC=1.0D0  +0.005D0  *(FKGES-2.0D0  )                                        
      VCORIN=1.0D0  /VCOREC                                                       
      IF(KGES.EQ.2) GO TO 235                                                   
      DO 234 NX=1,NXRA12                                                        
      VCENTR(NX)=VCENTR(NX)*VCOREC                                              
  234 CONTINUE                                                                  
  235 EGESIN=VSX-EGES                                                           
      FKAPPA=UNITOK   * DSQRT(RMAS*EGES)                                         
      FKAPIN=UNITOK   * DSQRT(RMAS*EGESIN)                                       
      CALL UNCPST                                                               
      IF(KEX42) 317,380,316                                                     
  316 CORNOD=1.1D0 *(1.0D0  -FKTRL2)+0.9D0  *FKTRL2                                  
      GO TO 318                                                                 
  317 CORNOD=0.9D0  *(1.0D0  -FKTRL2)+1.1D0  *FKTRL2                                  
  318 IF(KTRL2.EQ.1) GO TO 323                                                  
      EGEST=EGEST*CORNOD                                                        
      GO TO 341                                                                 
  323 VSX=VSX*CORNOD                                                            
      DO 335 NX=1,NXRA12                                                        
      VCENTR(NX)=VCENTR(NX)*CORNOD                                              
  335 CONTINUE                                                                  
  341 ITENOD=ITENOD+1                                                           
      NODMEM(ITENOD)=KEX42+NODE                                                 
      KEX42=0                                                                   
      IF(ITENOD.LE.10) GO TO 205                                                
      WRITE(6,354) (NODMEM(I1),I1=1,ITENOD)                                     
  354 FORMAT(10H ITENOD=10,10I3)                                                
CCCCCC  ******        CONVERGENCE HAS NOT BEEN ACHIEVED AND       ******        
CCCCCC  ******        THE CALCULATION IS GIVEN UP                 ******        
      STOP                                                                      
  380 IF(KTRL2.EQ.0) GO TO 390                                                  
      IF(KGES.EQ.2) GO TO 390                                                   
      DO 385 NX=1,NXRA12                                                        
      VCENTR(NX)=VCENTR(NX)*VCORIN                                              
  385 CONTINUE                                                                  
  390 CONTINUE                                                                  
  400 CONTINUE                                                                  
      GEST=EGEST*(1.0D0  -FKTRL2)+VSX*FKTRL2                                      
      IF( DABS(DELGES/GEST).LE.PERCNT) GO TO 417                                 
      DELGES=   (GEST*DELGES/ DABS(DELGES))*PERCNT                               
  417 GEST=GEST-DELGES                                                          
      IF(KTRL2.EQ.1) GO TO 427                                                  
      EGEST=GEST                                                                
      GO TO 440                                                                 
  427 VCOREC=GEST/(GEST+DELGES)                                                 
      VSX=VSX*VCOREC                                                            
      DO 430 NX=1,NXRA12                                                        
      VCENTR(NX)=VCENTR(NX)*VCOREC                                              
  430 CONTINUE                                                                  
  440 IF( ABS(DELGES/GEST).LE.ACURCY) GO TO 505                                 
      IF(KOPOT.NE.0)                                                            
     1WRITE (6,447)  ITEBE,VSX,EGEST,EGES,EGESIN,DELGES,GEST                    
  447 FORMAT (I10, 8F10.5)                                                      
      IF(ITEBE.LE.ITBEMX) GO TO 210                                             
      WRITE(6,455)                                                              
  455 FORMAT(//15H0ITEBE=ITBEMX+1)                                              
C      IF(KTRL2.EQ.1) GO TO 472                                                  
C      ISERCH=IQI                                                                
C      GO TO 475                                                                 
C  472 ISERCH=IQJ                                                                
C  475 WRITE(6,520) ISERCH,(GESMEM(I),I=1,ITEBE)                                 
      IF(KTRL2.EQ.0) WRITE(6,520) (GESMEM(I),I=1,ITEBE)                                                   
      IF(KTRL2.EQ.1) WRITE(6,521) (GESMEM(I),I=1,ITEBE)                                                   
      STOP                                                                      
  505 EGES=EGEST                                                                
CCCCCC  ******             CONVERGENCE HAS BEEN ACHIEVED          ******        
      KGES=2                                                                    
      EGESIN=VSX-EGES                                                           
      FKAPPA=UNITOK   * DSQRT(RMAS*EGES)                                         
      FKAPIN=UNITOK   * DSQRT(RMAS*EGESIN)                                       
      KEX41=1                                                                   
      KEX40=0                                                                   
      CALL UNCPST                                                               
      KEX40=1                                                                   
      CALL UNCPST                                                               
      WRITE(6,510) EGEST,VSX                                                    
  510 FORMAT( /  10X,36HFINAL VALUE OF THE PARAMETERS. B.E.=F10.5,              
     1  5H VSX=F10.5)                                                           
      ENEPT=EGEST                                                               
      VNEPT=VSX                                                                 
C      IF(KTRL2.EQ.1) GO TO 512                                                  
C      ISERCH=IQI                                                                
C      GO TO 515                                                                 
C  512 ISERCH=IQJ                                                                
  515 IF(KTRL2.EQ.0) WRITE(6,520) (GESMEM(I),I=1,ITEBE)                                 
      IF(KTRL2.EQ.1) WRITE(6,521) (GESMEM(I),I=1,ITEBE)                                 
  520 FORMAT( 10X,19HB.E. HAS VARIED AS 10F8.4/29X,10F8.4)                       
  521 FORMAT( 10X,19HVSX  HAS VARIED AS 10F8.4/29X,10F8.4)                       
      IF(ITENOD.EQ.0) GO TO 527                                                 
      WRITE(6,525) ITENOD,(NODMEM(I1),I1=1,ITENOD)                              
  525 FORMAT(/  15X,7HITENOD=I2,11H,WITH NODE=10I3)                             
  527 NXRMWF=NXRMP1                                                             
      IF(URRMEX.NE.0.D0  ) GO TO 533                                              
      IF(URRMIN.NE.0.D0  ) GO TO 545                                              
      T1=0.0D0                                                                    
      GO TO 534                                                                 
  533 T1=URRMIN/URRMEX                                                          
  534 DO 540 NX=NXRMWF,NXRAWF                                                   
  540 URSAVE(NX   )=URSAVE(NX   )*T1                                            
  545 IIPA=0                                                                    
      SS=0.0D0                                                                    
      DO 575 NX=1,NXRAWF                                                        
      IF(IIPA.EQ.1) GO TO 552                                                   
      F1=4.0D0                                                                    
      IIPA=1                                                                    
      GO TO 553                                                                 
  552 F1=2.0D0                                                                    
      IIPA=0                                                                    
  553 SS=SS+(URSAVE(NX)**2)*F1                                                  
  575 CONTINUE                                                                  
      SS=SS*XMES2*.3333333333D0                                                   
      SQ1= DSQRT(1.D0  /SS)                                                        
      DO 580 NX=1,NXRAWF                                                        
  580 URSAVE(NX)=URSAVE(NX)*SQ1                                                 
      X=0.D0                                                                      
      DO 647 NX=1,NXRAWF                                                        
      X=X+XMES2                                                                 
  647 URSAVE(NX   )=URSAVE(NX   )/X                                             
      NXSTEP=1                                                                  
      IF(KOPOT.EQ.0) NXSTEP=10                                                  
      XMOUT=NXSTEP*XMES2                                                        
      WRITE(6,649) XMOUT                                                        
  649 FORMAT( //,18H WAVE FUNCTIONS AT,F7.2,15H FERMI INTERVAL)                 
      WRITE(6,655) (URSAVE(NX),NX=1,NXRAWF,NXSTEP)                              
  655 FORMAT(10E13.4)                                                           
      IF(MOD(KTRL4,2).EQ.0) GO TO 667                                           
CCCCCC  ******         MULTIPLICATION BY THE BINDING POTENTIAL    ******        
CCCCCC  ******          (FOR KTRL4 = 1 OR 3 )                     ******        
      LBTRTW=2*LBTR                                                             
      LLSQ=LBTRTW*(LBTRTW+2)                                                    
      ISSQ=ISTW*(ISTW+2)                                                        
      JJSQ=JBTRTW*(JBTRTW+2)                                                    
      SPFAC=JJSQ-LLSQ-ISSQ                                                      
      SPFAC=0.250D0  *SPFAC                                                       
      IF(ISTW.EQ.0) GO TO 657                                                   
      FISTW=ISTW                                                                
      SPFAC=SPFAC/FISTW                                                         
  657 URSAVE(1)=URSAVE(1 )*(VCENTR(8)-VSPIN(8)*SPFAC)                           
      URSAVE(2)=URSAVE(2 )*(VCENTR(12)-VSPIN(12)*SPFAC)                         
      URSAVE(3)=URSAVE(3 )*(VCENTR(14)-VSPIN(14)*SPFAC)                         
      DO 661 NX=4,NXRAWF                                                        
      MX=NX+12                                                                  
  661 URSAVE(NX)=URSAVE(NX)*(VCENTR(MX)-VSPIN(MX)*SPFAC)                        
      WRITE(6,662)                                                              
  662 FORMAT( //,33H WAVE FUNCTION TIMES VCENTR+VSPIN)                          
      WRITE(6,655) (URSAVE(NX),NX=1,NXRAWF,NXSTEP)                              
      IF(TMAS.GE.4.0) GO TO 667                                                 
      S=0.0D0                                                                     
      IIPA=0                                                                    
      VFAI=0.0D0                                                                  
      DO 665 NX=1,NXRAWF                                                        
      S=S+XMES2                                                                 
      IF(IIPA.EQ.1) GO TO 663                                                   
      F1=4.0D0                                                                    
      IIPA=1                                                                    
      GO TO 664                                                                 
  663 F1=2.0D0                                                                    
      IIPA=0                                                                    
  664 VFAI=VFAI+F1*S**2*URSAVE(NX)                                              
  665 CONTINUE                                                                  
      VFAI=VFAI*XMES2*0.333333333333D0                                            
      VFAI=VFAI**2*4.0D0*PAI                                                      
      WRITE(6,666)VFAI                                                          
  666 FORMAT(//9X,49HCORRESPONDING ZERO RANGE NORMALIZATION D SQUARE =,         
     1    E20.5)                                                                
  667 IF(KTRL4.LE.1) GO TO 674                                                  
CCCCCC  ******        DIVISION OF WAVE FUNCTION BY R**LBTR        ******        
CCCCCC  ******        IS MADE ( FOR KTRL4 = 2 OR 3 )              ******        
      IF(LBTR.EQ.0) GO TO 669                                                   
      R=0.                                                                      
      DO 668 NX=1,NXRAWF                                                        
      R=R+XMES2                                                                 
      URSAVE(NX)=URSAVE(NX)*(R**(-LBTR))                                        
  668 CONTINUE                                                                  
  669 IF(KTRL4.EQ.2) WRITE(6,672)                                               
      IF(KTRL4.EQ.3) WRITE(6,673)                                               
  672 FORMAT( //,31H WAVE FUNCTION TIMES R**(-LBTR))                            
  673 FORMAT( //,50H WAVE FUNCTION TIMES VCENTR+VSPIN TIMES R**(-LBTR))         
      WRITE(6,655) (URSAVE(NX),NX=1,NXRAWF,NXSTEP)                              
  674 IF(KTRL8.EQ.0) GO TO 700                                                  
      NXREPT=(NXRAWF-1)/5+1                                                     
      DO 678 NXR=1,NXREPT                                                       
      NXI=5*(NXR-1)+1                                                           
      NXM=NXI+4                                                                 
C      PUNCH 676,(URSAVE(NX   ),NX=NXI,NXM),NXR                                  
C  676 FORMAT(5E14.6,7X,I3)                                                      
  678 CONTINUE                                                                  
  700 RETURN                                                                    
      END                                                                       
      SUBROUTINE UNCPST                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC,L9(10),U9                   
      COMMON/BIND/KTRL2,KTRL8,KTRL10,KEX2,KEX4,KEX40,KEX41,KEX42,               
     1            TMAS,PMAS,ZZT,ZZP,RMAS,ZZ,XMES1,ACURCY,PERCNT,VSX,            
     2            ISTW,NXRA,NXRM,NXRMP1,NXRMP2,NXRMP3,NXRMP4,NXRMP5,            
     3            NODE,KGES,EGESRD,EGES,EGEST,DELGES,FKAPPA,FKAPIN,             
     4            URRMIN,URRMEX,RGDLIN(3),RGDLEX(3)                             
      COMMON URSAVE(500),NODER,LBTR,JBTRTW,NSPAC2,ENEPT,VNEPT,                  
     1       VSOR,DFNR,DFNSO,RZR,RZSO,RZC,XMES2,                                
     2       XMEM(514),VCENTR(514),VSPIN(514),VCOULM(514),                      
     3          PFORM1(3),PFORM2(3),PFORM3(3),PFORM4(3),PFORM5(3),              
     4          VLAMD1(3,3),VLAMD2(3,3),VTERM1(3,3),VTERM2(3,3),                
     5          UR(4),UI(4),BR(4),FPRER(5),FPRERM(5),URSVIN(16)                 
CCCCCC  ******    STORMER COEFFICIENTS AND SPIN ORBIT FACTORS     ******        
      PRE41=19.D0  /240.D0                                                          
      PRE42=-2.D0  /5.D0                                                            
      PRE43=97.D0  /120.D0                                                          
      PRE44=-11.D0  /15.D0                                                          
      PRE45=299.D0  /240.D0                                                         
      FKEX40=KEX40                                                              
      VENSFC=1.D0  /(RMAS*0.04783258D0  )                                           
      WNI=FKAPIN                                                                
      MULSPN=ISTW+1                                                             
      LBTRTW=2*LBTR                                                             
      LLSQ=LBTRTW*(LBTRTW+2)                                                    
      ISSQ=ISTW*(ISTW+2)                                                        
      JJSQ=JBTRTW*(JBTRTW+2)                                                    
      LLTR=LBTRTW/2                                                             
      LL=LLTR+1                                                                 
      ENSFAC=LLTR*(LLTR+1)*VENSFC                                               
      SPFAC=JJSQ-LLSQ-ISSQ                                                      
      SPFAC=0.250D0  *SPFAC                                                       
      IF(ISTW.EQ.0) GO TO 425                                                   
      FISTW=ISTW                                                                
      SPFAC=SPFAC/FISTW                                                         
  425 NDFREP=(1-KEX40)*4+KEX40                                                  
      UR1M=1.0D0                                                                  
CCCCCC  ******       BEGINS NUMERICAL INTEGRATION WITH            ******        
CCCCCC  ******       FOUR STEP OR SINGLE STEP STORMER METHOD      ******        
      DO 730 ND=1,NDFREP                                                        
      FND=2.0D0  **(ND-1)                                                         
      DDX=XMES1*FND*(1.D0  -FKEX40)-XMES2*FKEX40                                  
      DR=DDX*FKAPPA                                                             
      DRSQ=(DR*DR)/EGES                                                         
      IF(KEX40.EQ.1) GO TO 510                                                  
CCCCCC  ******    PREPARATION FOR INTEGRATING INTERNAL SOLUTION   ******        
CCCCCC  ******    OUTWARDS WITH FOUR STEP STORMER METHOD          ******        
      K4COR2=4*ND-5                                                             
      IF(ND.NE.1) GO TO 485                                                     
      UR(1)=.0D0                                                                  
      FPRERM(1)=.0D0                                                              
      FPRER(1)=.0D0                                                               
      UR(2)=.5D0 *((2.D0*DDX*WNI)**LL)* DEXP(FACLOG(LL)-FACLOG(2*LL))             
      IF(LL.NE.2) GO TO 473                                                     
      FPRER(1)=VENSFC*0.666666666666D0  *WNI*WNI*DRSQ                             
  473 KEX42=0                                                                   
      NODCAL=0                                                                  
      FNODFC=1.0D0                                                                
      IF(MOD(NODE,2).EQ.0) GO TO 478                                            
      FNODFC=-1.0D0                                                               
  478 UR(2)=UR(2)*FNODFC                                                        
      URSAVE(1)=UR(2)                                                           
      FPRERM(1)=FPRER(1)*FNODFC                                                 
      FPRER(1)=FPRERM(1)                                                        
      X=0.0D0                                                                     
      NXIMIN=2                                                                  
      NXIMAX=8                                                                  
      NXPRCH=3                                                                  
      NNX=2                                                                     
      N3COR=0                                                                   
      GO TO 540                                                                 
  485 NXIMIN=5                                                                  
      X=3.D0  *DDX                                                                
      NXIMAX=8                                                                  
      IF(ND.NE.4) GO TO 491                                                     
      NXIMAX=NXRMP5                                                             
  491 DO 495 NQ=1,4                                                             
      URSAVE(NQ)=URSAVE(2*NQ)                                                   
  495 FPRER(NQ)=FPRERM(NQ)*4.0D0                                                  
      FPRERM(1)=FPRER(1)                                                        
      FPRERM(2)=FPRER(3)                                                        
      UR(1)=UR1M                                                                
      NXPRCH=5                                                                  
      NNX=3                                                                     
      N3COR=0                                                                   
      GO TO 540                                                                 
  510 XRA=NXRA                                                                  
CCCCCC  ******    PREPARATION FOR INTEGRATING EXTERNAL SOLUTION   ******        
CCCCCC  ******     INWARDS WITH SINGLE STEP STORMER METHOD        ******        
      XRA=XMES2*XRA                                                             
      NXRAM1=NXRA-1                                                             
      FPRER(1)=0.0D0                                                              
      UR(1)=FPRER(1)                                                            
      URSAVE(NXRA)=UR(1)                                                        
      UR(2)=0.0001D0                                                              
      URSAVE(NXRAM1)=UR(2)                                                      
      X=XRA                                                                     
      NXIMIN=2                                                                  
      NXIMAX=NXRA-NXRM-1                                                        
      N3COR =NXRA+13                                                            
  540 DO 700 NX=NXIMIN,NXIMAX                                                   
      X=X+DDX                                                                   
      K2=(NX+K4COR2)*(1-KEX40)+(N3COR-NX)*KEX40                                 
      K3=NX        *(1-KEX40)+(NXRA -NX)*KEX40                                  
      AR2=VCENTR(K2)+VCOULM(K2)+ENSFAC/(X*X)-SPFAC*VSPIN(K2)+EGES               
      BR2=AR2*DRSQ*UR(2)                                                        
      IF(ND-1) 551,551,554                                                      
  551 IF(NX-4) 552,552,554                                                      
  552 TERMR=BR2                                                                 
      FPRER(NX)=TERMR                                                           
      GO TO 555                                                                 
  554 FPRER(5)=BR2                                                              
      TERMR=PRE41*FPRER(1)+PRE42*FPRER(2)+PRE43*FPRER(3)                        
     1     +PRE44*FPRER(4)+PRE45*FPRER(5)                                       
  555 UR(3)=2.0D0  *UR(2)-UR(1)+TERMR                                             
      URSAVE(K3)=UR(3)                                                          
      UR(1)=UR(2)                                                               
      UR(2)=UR(3)                                                               
      IF(ND.EQ.NDFREP) GO TO 570                                                
      IF(NX.NE.NXPRCH) GO TO 570                                                
      FPRERM(NNX)=BR2                                                           
      NXPRCH=NXPRCH+2                                                           
      NNX=NNX+1                                                                 
      IF(NX.NE.7) GO TO 570                                                     
      UR1M=UR(1)                                                                
  570 IF(ND.NE.1) GO TO 580                                                     
      IF(NX.LE.4) GO TO 585                                                     
  580 DO 583 NQ=1,4                                                             
  583 FPRER(NQ)=FPRER(NQ+1)                                                     
  585 IF(KEX40.EQ.1) GO TO 700                                                  
      IF(UR(1)*UR(2).LT.0.0D0  ) NODCAL=NODCAL+1                                  
  700 CONTINUE                                                                  
  730 CONTINUE                                                                  
CCCCCC  ******      CALCULATE LOGARITHMIC DERIVATIVE OF WAVE      ******        
CCCCCC  ******      FUNCTIONS AT THE MATCHING RADIUS AND          ******        
CCCCCC  ******      ESTIMATE NEW SEARCH PARAMETER.                ******        
      IF(KEX40.EQ.1) GO TO 760                                                  
      KEX42=NODCAL-NODE                                                         
      IF(KEX42.NE.0) GO TO 1000                                                 
  760 URRM=URSAVE(NXRMP3)                                                       
      RGDL=(8.0D0*(URSAVE(NXRMP4)-URSAVE(NXRMP2))                               
     1           -(URSAVE(NXRMP5)-URSAVE(NXRMP1)))/(12.D0*XMES2*URRM)           
      IF(KEX40.NE.0) GO TO 770                                                  
      RGDLIN(KGES)=RGDL                                                         
      GO TO 773                                                                 
  770 RGDLEX(KGES)=RGDL                                                         
  773 IF(KGES.NE.2) GO TO 910                                                   
      IF(KEX40.NE.0) GO TO 777                                                  
      URRMIN=URRM                                                               
      GO TO 910                                                                 
  777 URRMEX=URRM                                                               
  910 IF(KEX40.EQ.0) GO TO 1000                                                 
      IF(KGES.NE.3) GO TO 1000                                                  
      IF(KTRL2.EQ.1) GO TO 923                                                  
      DEN=.005D0  *EGEST                                                          
      GO TO 925                                                                 
  923 DEN=.005D0  *VSX                                                            
  925 DLEVIN=(RGDLIN(3)-RGDLIN(1))/DEN                                          
      DLEVEX=(RGDLEX(3)-RGDLEX(1))/DEN                                          
      DELGES=(RGDLIN(2)-RGDLEX(2))/(DLEVIN-DLEVEX)                              
 1000 RETURN                                                                    
      END                                                                       
      SUBROUTINE MONITR                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/RARB/NATOTL,NBMIN,NBMAX,NBREAD,NBTAIL,MESFCA,MESFCB,               
     1            N1MAX,N2MAX,MGAUS,NGAUS                   
      COMMON/RARB1/XMESA,XMESB,XMES(2)                   
      COMMON/GAUS/ABSCIS(222),WEIGHT(222),WWL1(500),WWL2(500),                  
     1            BSPAR(46),IBSPR(10)                                           
      WRITE(6,201)                                                              
  201 FORMAT(1H1,20X,47HTHE FOLLOWING STATEMENTS ARE ISSUED FROM MONITR)        
      IF(JLSMAX.GT.5) WRITE(6,202)                                              
  202 FORMAT(12H JLSMAX GT 5)                                                   
      IF(JLSMAX.GT.5) STOP                                                      
      KSTOP=0                                                                   
      QVALGR=BSPAR(10)-BSPAR(1)                                                 
      IF(MASSA.GT.MASSB) QVALGR=-QVALGR                                         
      IF(KTRLST(1).EQ.1) GO TO 501                                              
      LDARMX=(L1TRTW(1)+L2TRTW(1))/2+1                                          
      NGRNMX=0                                                                  
      L12MOD= MOD(LDARMX-1,2)                                                   
      LBSTTL=0                                                                  
      LLMAX=0                                                                   
      DO 140 JLS=1,JLSMAX                                                       
      NGRUN=0                                                                   
      LL=LLTWR(JLS)/2                                                           
      LBSTTL=LBSTTL+(LL+1-MOD(LL+KPCHAG,2))                                     
      LLMAX=MAX0(LL,LLMAX)                                                      
      DO 130 LDAR=1,LDARMX                                                      
      LDA=LDAR-1                                                                
      LDBRMI=IABS(LDA-LL)+1                                                     
      LDBRMX=MIN0(LDA+LL,LDARMX-1)+1                                            
      IF(LDBRMI.GT.LDBRMX) GO TO 130                                            
      DO 125 LDBR=LDBRMI,LDBRMX                                                 
      IF(MOD(LDAR+LDBR,2).NE.L12MOD) GO TO 125                                  
      NGRUN=NGRUN+1                                                             
  125 CONTINUE                                                                  
  130 CONTINUE                                                                  
      NGRNMX=MAX0(NGRNMX,NGRUN)                                                 
  140 CONTINUE                                                                  
      KEXCST(15)=NGRNMX                                                         
      NN=NGRNMX*LDARMX*JLSMAX                                                   
      IF(NN.LE.1500) GO TO 161                                                  
      KSTOP=1                                                                   
      WRITE(6,145)NN                                                            
  145 FORMAT(///1X,20(1H*),5X, 5H  NN=,I5,14H EXCEEDS 1500.,/26X,               
     1 51H INCREASE DIMENSION OF CFAC(1500) IN CFCCAL, FFRARB)                  
  161 CONTINUE                                                                  
      NBTOTL=NBMAX-NBMIN+1                                                      
      KMAX=LDWMAX+(L1TRTW(1)+L2TRTW(1))/2                                       
      IF(KMAX.LT.250) GO TO 220                                                 
      KSTOP=1                                                                   
      WRITE(6,215)KMAX                                                          
  215 FORMAT(///1X,20(1H*),6H KMAX=,I5,50H GT. 250. INCREASE DIMENSION          
     1OF PK(250) IN KERNEL )                                                    
  220 MSFCRT=(XMESB/XMESA)+0.1D0                                                  
      IF(MOD(MESFCB,MESFCA).NE.0) GO TO 310                                     
      IF(MESFCB/MESFCA .EQ.MSFCRT) GO TO 320                                    
  310 KSTOP=1                                                                   
      WRITE(6,315)MESFCA,MESFCB,XMESA,XMESB                                     
  315 FORMAT(///1X,20(1H*),5X,27H MESFCA,MESFCB,XMESA,XMESB=,2I3,2F7.3,         
     1 17H ARE INCOMPATIBLE)                                                    
  320 NARTOT=(NBTOTL-1)*MSFCRT+NATOTL                                           
      L1L2TR=(L1TRTW(1)+L2TRTW(1))/2                                            
      NRAPWC=NARTOT*(L1L2TR+1)                                                  
      IF(NRAPWC.LE.3600) GO TO 330                                              
      KSTOP=1                                                                   
      WRITE(6,325)NRAPWC                                                        
  325 FORMAT(///1X,20(1H*),5X,8H NRAPWC=,I5,14H EXCEEDS 3600.,/26X,             
     1 53H INCREASE DIMENSION OF RAPOWR(3600) IN FRFMFC, FFRARB)                
  330 NRBPWC=NBTOTL*(L1L2TR+1)                                                  
      IF(NRBPWC.LE.500) GO TO 340                                               
      KSTOP=1                                                                   
      WRITE(6,335) NRBPWC                                                       
  335 FORMAT(///1X,20(1H*),5X,8H NRBPWC=,I5,13H EXCEEDS 500.,/26X,              
     1 53H INCREASE DIMENSION OF RBPOWR(500)  IN FRFMFC, FFRARB)                
  340 KRANGE=L1TRTW(1)+L2TRTW(1)+1                                              
      NNMAX=LDARMX*LBSTTL*KRANGE                                                
      IF(NNMAX.LE.2000) GO TO 360                                               
      KSTOP=1                                                                   
      WRITE(6,351)NNMAX                                                         
  351 FORMAT(//1X,10(1H*),7H NNMAX=,I5,13H EXCEEDS 2000,                        
     1 50H INCREASE DIMENSION OF GEOK(2000) IN FRFMFC,FFRARB)                   
  360 NGUSE=NATOTL*NBTOTL*KRANGE                                                
      IF(NGUSE.LE.4000) GO TO 420                                               
      KSTOP=2                                                                   
      WRITE(6,415) NGUSE                                                        
  415 FORMAT(///1X,20(1H*),5X,7H NGUSE=,I5,14H EXCEEDS 4000.,/26X,              
     162HINCREASE DIMENSION OF GTAPE(4000) IN FRFMFC, FFRARB, OR SIMPLY,        
     2 36H INCREASE THE REQUESTED FIELD LENGTH)                                 
  420 NAFMES=(NATOTL-1)*MESFCA+1                                                
      NAFBLK=LBSTTL*NAFMES                                                      
      NABLOK=LBSTTL*NATOTL                                                      
      MSTOT=(NBMAX-NBMIN-2)*MESFCB                                              
      IF(NAFBLK.LE.650) GO TO 440                                               
      KSTOP=2                                                                   
      WRITE(6,435)NAFBLK                                                        
  435 FORMAT(///1X,20(1H*),5X,14HNAFMES*LBSTTL=,I5,13H EXCEEDS 650.,/26X        
     1,49H INCREASE DIMENSION OF FFORA(650) IN MARS, INTPOL)                    
  440 IF(NABLOK.LE.650) GO TO 450                                               
      KSTOP=2                                                                   
      WRITE(6,445)NABLOK                                                        
  445 FORMAT(///1X,20(1H*),5X,14HNATOTL*LBSTTL=,I5,13H EXCEEDS 650.,/26X        
     1,61H INCREASE DIMENSION OF FB1,FB2,FB3,FB4 AND FFA(650) IN INTPOL)        
  450 NBRD1=500/NATOTL                                                          
      NBRD2=2500/(NATOTL*LBSTTL)                                                
      NBRD3=(1500-MSTOT)/(LBSTTL*MESFCB)+3                                      
      NBREAD=MIN0(NBRD1,NBRD2,NBRD3,NBTOTL)                                     
      IF(NBREAD.GE.4) GO TO 460                                                 
      KSTOP=1                                                                   
      WRITE(6,455)NBRD1,NBRD2,NBRD3                                             
  455 FORMAT(///1X,20(1H*),19H NBRD1,NBRD2,NBRD3=,3I5,/,20X,                    
     1      50H ALL THE ABOVE MUST BE GREATER THAN OR EQUAL TO 4 ,/20X,         
     2      50H INCREASE THE NECESSARY DIMENSION AS DEFINED BELOW,/20X,         
     3      55H NBRD1=500/NATOTL. INCREASE FTEM1(500) IN FRFMFC,FFRARB,         
     4 /20X,50H NBRD2=2500/(NATOTL*LBSTTL). INCREASE FTAPE(2500) ,              
     5 /20X,50H        IN FRFMFC,FFRARB AND ALSO FORM(2500),     ,              
     6 /20X,50H        FTAPE(2500) IN MARS,INTGMS,INTPOL.        ,              
     7 /20X,50H NBRD3=(1500-MSTOT)/(LBSTTL*MESFCB)+3.            ,              
     8 /20X,50H       INCREASE DIMENSION OF URISV(1500) IN INTGMS,              
     9 /20X,50H       AND INTPOL                                 )              
  460 IF(LBSTTL.LE.40) GO TO 470                                                
      KSTOP=1                                                                   
      WRITE(6,465) LBSTTL                                                       
  465 FORMAT(///1X,20(1H*),7HLBSTTL=,I3,12H EXCEEDS 40.,/30X,                   
     1 60HINCREASE ALL THOSE DIMENSIONED QUANTITIES IN FFRARB FROM 40.          
     2   ,/30X,51HAND THOSE IN OVLAP,INTGMS,INTPOL,URENOM AND DWAMP. )          
  470 LROWMX=LLMAX+1-MOD(LLMAX+KPCHAG,2)                                        
      IF(LROWMX.LE.15) GO TO 480                                                
      KSTOP=2                                                                   
      WRITE(6,475) LROWMX                                                       
  475 FORMAT(///1X,20(1H*),7HLROWMX=,I3,12H EXCEEDS 15.,/30X,                   
     1 62HINCREASE LLROW(15) IN THE BLOCK COMMON /JP1/ AND ALL THE BLANK        
     2 ,/30X,49H  COMMON RELATED TO THIS DIMENSION OF 15 IN MARS.)              
  480 IF(KSTOP.EQ.1) STOP                                                       
CCCCCC                                                                          
      IF(KSTOP.NE.0) WRITE(6,492)                                               
  492 FORMAT(///1X,20(1H*),41HCHECK THE ABOVE STATEMENT, IF ANY, BEFORE,        
     1 13H RUNNING MARS)                                                        
  501 WRITE(6,502)                                                              
  502 FORMAT(///1X,20(1H*),40H THE FOLLOWING QUANTITIES ARE WRITTEN ON,         
     1 8H TAPE-12)                                                              
CCCCCC  ******             STORE BASIC QUANTITIES ON TAPE12       ******        
CCCCCC  ******             TO BE READ BACK IN MARS                ******        
      WRITE(12)  MASCA,NZCA,KPCA,ICATW,MASCB,NZCB,KPCB,ICBTW,                   
     1           MASSA,NZSA,KPSA,ISATW,MASSB,NZSB,KPSB,ISBTW,                   
     2           JLSMAX,LDWMAX,NBREAD,LBSTTL,KEXCST(1),KEXCST(2)                
      WRITE(12) (JJTWR(N),LLTWR(N),ISTWR(N),N=1,JLSMAX)                         
      IF(KTRLST(1).EQ.0)                                                        
     1WRITE(12) NATOTL,NBMIN,NBMAX,MESFCA,MESFCB,XMESA,XMESB,QVALGR             
      KEXCB1=(NBMIN-1)*MESFCB                                                   
      KEXCB2=NBMAX*MESFCB                                                       
      MST=KEXCB2-KEXCB1                                                         
      RMES=XMESB/MESFCB                                                         
      MESFAC=1                                                                  
      IF(KTRLST(1).EQ.1) WRITE(12) KEXCB1,KEXCB2,MST,MESFAC,RMES ,QVALGR        
      WRITE(12) (IBSPR(N),N=1,6),(BSPAR(N),N=1,18)                              
      IF(KTRLST(1)*KEXCST(14).NE.1) GO TO 513                                   
  504 FORMAT(24I3)                                                              
  505 FORMAT(4I5,2F10.4)                                                        
  506 FORMAT(6I5,/,9F7.2,/,9F7.2)                                               
C      PUNCH 504, MASCA,NZCA,KPCA,ICATW,MASCB,NZCB,KPCB,ICBTW,                   
C     1           MASSA,NZSA,KPSA,ISATW,MASSB,NZSB,KPSB,ISBTW,                   
C     2           JLSMAX,LDWMAX,NBREAD,LBSTTL,KEXCST(1),KEXCST(2)                
C      PUNCH 504,(JJTWR(N),LLTWR(N),ISTWR(N),N=1,JLSMAX)                         
C      PUNCH 505, KEXCB1,KEXCB2,MST,MESFAC,RMES ,QVALGR                          
C      PUNCH 506,(IBSPR(N),N=1,6),(BSPAR(N),N=1,18)                              
  513 WRITE(6,515) MASCA,NZCA,KPCA,ICATW,MASCB,NZCB,KPCB,ICBTW,                 
     1           MASSA,NZSA,KPSA,ISATW,MASSB,NZSB,KPSB,ISBTW,                   
     2           JLSMAX,LDWMAX,NBREAD,LBSTTL,KEXCST(1),KEXCST(2)                
  515 FORMAT(//45H MASCA,NZCA,KPCA,ICATW,MASCB,NZCB,KPCB,ICBTW=,8I5,/           
     1         45H MASCA,NZSA,KPSA,ISATW,MASSB,NZSB,KPSB,ISBTW=,8I5,/           
     2         45H JLSMAX,LDWMAX,NBREAD,LBSTTL,KEX1,KEX2      =,6I5)            
      WRITE(6,517)(JJTWR(N),LLTWR(N),ISTWR(N),N=1,JLSMAX)                       
  517 FORMAT(//28H JJTWR(N),LLTWR(N),ISTWR(N)=,15I5)                            
      WRITE(6,519)NATOTL,NBMIN,NBMAX,MESFCA,MESFCB,XMESA,XMESB,QVALGR           
  519 FORMAT(/53H NATOTL,NBMIN,NBMAX,MESFCA,MESFCB,XMESA,XMESB,QVALGR=,         
     1 5I5,3F13.5)                                                              
      WRITE(6,521)(BSPAR(N),N=1,18)                                             
  521 FORMAT(47H ENEPT,VNEPT,VSOR,DFNR,DFNSO,RZR,RZSO,RZC,XMES=,9F7.2)          
      RETURN                                                                    
      END                                                                       
      SUBROUTINE NRFMFC                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC,L9(10),U9                   
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/RARB/NATOTL,NBMIN,NBMAX,NBREAD,NBTAIL,MESFCA,MESFCB,               
     1            N1MAX,N2MAX,MGAUS,NGAUS                   
      COMMON/RARB1/XMESA,XMESB,XMES(2)                   
      COMMON/GAUS/ABSCIS(222),WEIGHT(222),WWL1(500),WWL2(500),                  
     1            BSPAR(46),IBSPR(10)                                           
      COMMON FORMFC(400,5),GEOFC(200),JLSMEM(200),LD1MEM(200),                  
     1       LD1PME(200),KP1MEM(200),RAMM(200),WFACMM(200),WL2MM(200),          
     2       NRTMEM(10),SUMALL(10),RTS(96),WGT(96),PK(96),                      
     3       RAD(2),W(2),XMAX(2),DE(4)                                          
      DIMENSION NUMG(12),IPOS(12)                                               
      DATA NUMG / 4,8,12,16,20,24,32,40,48, 64, 80, 96 /                        
      DATA IPOS / 1,3, 7,13,21,31,43,59,79,103,135,175 /                        
      KOUT=KOUTST(5)                                                            
      KPUNCH=KEXCST(14)                                                         
      SQRTPI= DSQRT(3.14159265359D0  )                                             
      WRITE(6,5)                                                                
    5 FORMAT(1H1,40X,17HNRFMFC IS CALLED,//)                                    
      XMAX(1)=XMES(1)*(N1MAX-2)                                                 
      XMAX(2)=XMES(2)*(N2MAX-2)                                                 
      NGAUS=NUMG(MGAUS)                                                         
      IP=IPOS(MGAUS)                                                            
      NG2=NGAUS/2                                                               
      DO 20 NG=1,NG2                                                            
      RTS(NG2-NG+1)=-ABSCIS(IP)                                                 
      WGT(NG2-NG+1)= WEIGHT(IP)                                                 
      RTS(NG+NG2)=ABSCIS(IP)                                                    
      WGT(NG+NG2)=WEIGHT(IP)                                                    
   20 IP=IP+1                                                                   
      IF(KOUT.LT.2) GO TO 25                                                    
      WRITE(6,22)NGAUS                                                          
   22 FORMAT(10X,23HINTEGRATING WITH NGAUS=,I3)                                 
      WRITE(6,24)(RTS(NG),NG=1,NGAUS)                                           
      WRITE(6,23)                                                               
   23 FORMAT(/)                                                                 
      WRITE(6,24)(WGT(NG),NG=1,NGAUS)                                           
   24 FORMAT(5(1X,E21.14))                                                      
   25 KMAX=0                                                                    
      TM1=MAX0(MASCA,MASCB)                                                     
      R1=BSPAR(6)*TM1**0.3333333333333333                                       
      TM1=MAX0(MASSA,MASSB)                                                     
      R12=R1+BSPAR(15)*TM1**0.333333333333333333                                
      R12=R12+3.0D0*(BSPAR(4)+BSPAR(13))                                          
      NBMIN=1                                                                   
      NBINT=R12/XMESB                                                           
      NBTAIL=NBMAX-NBINT+3                                                      
      WRITE(6,27)NBMIN,NBMAX,NBINT,NBTAIL,MGAUS,NGAUS,KPUNCH,N1MAX,N2MAX        
   27 FORMAT(20X,24HNBMIN,NBMAX,NBINT,NBTAIL,9X,1H=,4I5,                        
     1      /20X,34HMGAUS,NGAUS,KPUNCH,N1MAX,N2MAX   =,5I5)                     
      WRITE(6,29)XMESB,XMES(1),XMES(2)                                          
   29 FORMAT(20X,17HXMESB,XMES1,XMES2,16X,1H=,3F9.4)                            
      NRUNT=1                                                                   
      IXTW=MOD(IABS(MASCA-MASCB),2)                                             
      DO 50 JLS=1,JLSMAX                                                        
      LLTW=LLTWR(JLS)                                                           
      LL=LLTW/2                                                                 
      JJTW=JJTWR(JLS)                                                           
      ISTW=ISTWR(JLS)                                                           
      L1TR=L1TRTW(JLS)/2                                                        
      L2TR=L2TRTW(JLS)/2                                                        
      L1TW=2*L1TR                                                               
      L2TW=2*L2TR                                                               
      K1=IABS(ISTW+L2TW-IXTW)/2                                                 
      SSD=1-2*MOD(K1,2)                                                         
      HATD=1.0D0                                                                  
      IA=L1TW                                                                   
      IB=L2TW                                                                   
      IC=JJTW                                                                   
      ID=ISTW                                                                   
      IE=LLTW                                                                   
      IFF=IXTW                                                                   
      CALL RAC7                                                                 
      DCOE(JLS)=SSD*HATD*RAC                                                    
      L1P1=L1TR+1                                                               
      NRTMEM(JLS   )=NRUNT                                                      
      DO 40 LD1R=1,L1P1                                                         
      LD1=LD1R-1                                                                
      LD1P=L1TR-LD1                                                             
      LD1TW=2*LD1                                                               
      LD1PTW=2*LD1P                                                             
      CB1= DEXP(0.5D0*(FACLOG(L1TW+2)-FACLOG(LD1TW+2)-FACLOG(LD1PTW+2)))         
      HAT=(L1TW+1)*(L2TW+1)*(LD1TW+1)*(LD1PTW+1)                                
      HAT= DSQRT(HAT)                                                            
      CC=   0.5D0*CB1*HAT                                                         
      KMIN=MAX0(IABS(LL-LD1),IABS(L2TR-LD1P))+1                                 
      KMMX=MIN0(    (LL+LD1),    (L2TR+LD1P))+1                                 
      KMAX=MAX0(KMAX,KMMX)                                                      
      IF (KMIN.GT.KMMX) GO TO 40                                                
      DO 35 KP1=KMIN,KMMX,2                                                     
      K=KP1-1                                                                   
      K1=(L1TR+L2TR-LL)/2+K                                                     
      SS=1-2*MOD(K1,2)                                                          
      KTW=2*(K)                                                                 
      IA=LD1TW                                                                  
      IB=LLTW                                                                   
      IC=KTW                                                                    
      ID=0                                                                      
      IE=0                                                                      
      IFF=0                                                                      
      CALL CLEB                                                                 
      C1=RAC                                                                    
      IA=LD1PTW                                                                 
      IB=L2TW                                                                   
      CALL CLEB                                                                 
      C2=RAC                                                                    
      IA=L1TW                                                                   
      IB=L2TW                                                                   
      IC=LD1TW                                                                  
      ID=KTW                                                                    
      IE=LLTW                                                                   
      IFF=LD1PTW                                                                 
      CALL RAC7                                                                 
      R1=RAC                                                                    
      GEOFC(NRUNT)=SS*CC*C1*C2*R1                                               
      JLSMEM(NRUNT)=JLS                                                         
      LD1MEM(NRUNT)=LD1                                                         
      LD1PME(NRUNT)=LD1P                                                        
      KP1MEM(NRUNT)=KP1                                                         
      IF(KOUT.LT.2) GO TO 34                                                    
      WRITE(6,33)JLS,     L1TR,L2TR,LD1,K,NRUNT,CB1,HAT,CC,SS,C1,C2,R1,         
     1      GEOFC(NRUNT)                                                        
   33 FORMAT(42H NRFMFC-33,JLS,     L1TR,L2TR,LD1,K,NRUNT=,6I5,/,               
     1    30H CB1,HAT,CC,SS,C1,C2,R1,GEOFC=,8E11.3)                             
   34 NRUNT=NRUNT+1                                                             
   35 CONTINUE                                                                  
   40 CONTINUE                                                                  
   50 CONTINUE                                                                  
      NRUNTL=NRUNT-1                                                            
      NGAUS2=2*NGAUS                                                            
      AMASS=MAX0(MASSA,MASSB)                                                   
      R21=BSPAR(15)*AMASS**0.333333+1.0D0                                         
      R22=R21+6.0D0                                                               
      DO 90 I=1,2                                                               
      IF(I.EQ.2) GO TO 57                                                       
      RUP=R21                                                                   
      RANGE=R21                                                                 
      GO TO 61                                                                  
   57 RUP=R22                                                                   
      RANGE=RUP-R21                                                             
   61 DO 90 N=1,NGAUS                                                           
      NG2=(I-1)*NGAUS+N                                                         
      RA=RUP-(1.D0-RTS(N ))*RANGE/2.0D0                                             
      RAMM(NG2)=RA                                                              
      WFACMM(NG2)=WGT(N)*RANGE/2.D0                                               
      IR2=RA/XMES(2)                                                            
      IF(IR2.GE.2) GO TO 63                                                     
      DELX=RA/XMES(2)                                                           
      INMODE=1                                                                  
      GO TO 67                                                                  
   63 INMODE=2                                                                  
      IR2M2=IR2-2                                                               
      DELX=RA    /XMES(2)-IR2                                                   
      DEP1=DELX+1.D0                                                              
      DEM1=DELX-1.D0                                                            
      DEM2=DELX-2.D0                                                              
      DEM12=DEM1*DEM2/6.D0                                                        
      DEP1X=DEP1*DELX/6.D0                                                        
      DE(1)=-DELX*DEM12                                                         
      DE(2)= DEP1*DEM12*3.D0                                                      
      DE(3)=-DEM2*DEP1X*3.D0                                                      
      DE(4)= DEM1*DEP1X                                                         
   67 DO 80 JLS=1,JLSMAX                                                        
      IF(RA    .GT.XMAX(2)) GO TO 77                                            
      W(2)=-(DELX-2.D0  )*WWL2(1)   +(DELX-1.D0  )*WWL2(2)                          
      IF(INMODE.EQ.1) GO TO 76                                                  
      W(2)=0.D0                                                                   
      DO  75 IRR=1,4                                                            
      IR=IR2M2+IRR                                                              
   75 W(2)=W(2)+DE(IRR)*WWL2(IR)                                                
   76 WL2MM(NG2)=W(2)                                                           
      GO TO 80                                                                  
   77 WL2MM(NG2)=0.0D0                                                            
   80 CONTINUE                                                                  
   90 CONTINUE                                                                  
      NTOTAL=NBMAX+1                                                            
      DO 500 NB=NBMIN,NTOTAL                                                    
      RB=NB*XMESB                                                               
      DO 112 JLS=1,JLSMAX                                                       
  112 FORMFC(NB,JLS)=0.D0                                                         
      IF(KOUT.LT.2) GO TO 118                                                   
      WRITE(6,116)RB                                                            
  116 FORMAT(//,10X,3HRB=,F7.4)                                                 
  118 DO 400 NG2=1,NGAUS2                                                       
      RA=RAMM(NG2)                                                              
      FFFF=WFACMM(NG2)                                                          
      DO 122 JLS=1,JLSMAX                                                       
  122 SUMALL(JLS)=0.                                                            
      DO 200 NG=1,NGAUS                                                         
      ROOT=RTS(NG)                                                              
      PK(1)=1.D0                                                                  
      PK(2)=ROOT                                                                
      IF(KMAX.LT.3) GO TO 128                                                   
      DO 126 K=3,KMAX                                                           
      PK(K)=((2*K-3)*ROOT*PK(K-1)-(K-2)*PK(K-2))/(K-1)                          
  126 CONTINUE                                                                  
  128 WFAC=WGT(NG)                                                              
      DO 200 NRUNT=1,NRUNTL                                                     
      JLS=JLSMEM(NRUNT)                                                         
      LD1=LD1MEM(NRUNT)                                                         
      LD1P=LD1PME(NRUNT)                                                        
      KP1=KP1MEM(NRUNT)                                                         
      IF(NRUNT.NE.NRTMEM(JLS   )) GO TO 172                                     
      DC=DCOE(JLS)                                                              
      RAD(1)= SQRT(RA**2+RB**2+2.D0  *RA*RB*ROOT)                                 
      IF(RAD(1).GT.XMAX(1)) GO TO 200                                           
      IR2=RAD(1)/XMES(1)                                                        
      IF(IR2.GE.2) GO TO 135                                                    
      DELX=RAD(1)/XMES(1)                                                       
      W(1)=-(DELX-2.D0  )*WWL1(1)   +(DELX-1.D0  )*WWL1(2)                          
      GO TO 141                                                                 
  135 IR2M2=IR2-2                                                               
      DELX=RAD(1)/XMES(1)-IR2                                                   
      DEP1=DELX+1.D0                                                              
      DEM1=DELX-1.D0                                                              
      DEM2=DELX-2.D0                                                              
      DEM12=DEM1*DEM2/6.D0                                                        
      DEP1X=DEP1*DELX/6.D0                                                        
      DE(1)=-DELX*DEM12                                                         
      DE(2)= DEP1*DEM12*3.D0                                                      
      DE(3)=-DEM2*DEP1X*3.D0                                                      
      DE(4)= DEM1*DEP1X                                                         
      W(1)=0.D0                                                                   
      DO 137 IRR=1,4                                                            
      IR=IR2M2+IRR                                                              
  137 W(1)=W(1)+DE(IRR)*WWL1(IR)                                                
  141 W(2)=WL2MM(NG2)                                                           
      W12=WFAC*W(1)*W(2)*DC                                                     
      IF(KOUT.LT.3) GO TO 172                                                   
      WRITE(6,152)RAD(1),W(1),RAD(2),W(2),DC,W12                                
  152 FORMAT(27H RAD1,W1,RAD2,W2,DCOE,W12 =,7E12.4)                             
  172 RPOWER=RA**(LD1P+2)*RB**LD1                                               
      TERM=GEOFC(NRUNT)*PK(KP1)*RPOWER*W12                                      
      SUMALL(JLS)=SUMALL(JLS)+TERM                                              
      IF(KOUT.LT.2) GO TO 200                                                   
      N=KP1                                                                     
      WRITE(6,175)N,NRUNT,PK(N),GEOFC(NRUNT),RPOWER,TERM,SUMALL(JLS)            
  175 FORMAT(9H K,NRUNT=,2I5,30H PKK,GEOFC,RPOWER,TERM,SUMALL=,5E12.4)          
  200 CONTINUE                                                                  
      DO 395 JLS=1,JLSMAX                                                       
  395 FORMFC(NB,JLS)=FORMFC(NB,JLS)+FFFF*SUMALL(JLS)                            
  400 CONTINUE                                                                  
      IF(KOUT.LT.1) GO TO 500                                                   
      WRITE(6,418)RB,   (FORMFC(NB,JLS),JLS=1,JLSMAX)                           
  418 FORMAT(  10X,3HRB=,F7.4,5X,(7HFORMFC=,3E20.10))                           
  500 CONTINUE                                                                  
      DO 530 JLS=1,JLSMAX                                                       
      WRITE(6,510)JLS,NBMIN,NBMAX,XMESB                                         
  510 FORMAT(//20X,4HJLS=,I5,5X,12HNBMIN,NBMAX=,2I5,5X,6HXMESB=,F8.4,/)         
      WRITE(6,520) (FORMFC(N,JLS),N=NBMIN,NTOTAL)                               
  520 FORMAT(10E13.5)                                                           
  530 CONTINUE                                                                  
      CALL INTERP                                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE INTERP                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/RARB/NATOTL,NBMIN,NBMAX,NBREAD,NBTAIL,MESFCA,MESFCB,               
     1            N1MAX,N2MAX,MGAUS,NGAUS                   
      COMMON/RARB1/XMESA,XMESB,XMES(2)                   
      COMMON FORMFC(400,5),FTEM(400),FINAL(800),                                
     1       DEL(50),FC1(50),FC2(50),FC3(50),FC4(50)                            
      KPUNCH=KEXCST(14)                                                         
      NTIME1=MESFCB                                                             
      NBINT=NBMAX-NBTAIL+3                                                      
      NBP1=NBINT+1                                                              
      NM1=NTIME1-1                                                              
      DO 60 NN=1,NM1                                                            
      FNN=NN                                                                    
      FNTIME=NTIME1                                                             
      DEL(NN)=FNN/FNTIME                                                        
      DELN=DEL(NN)                                                              
      FC1(NN)=-DELN*(DELN-1.D0  )*(DELN-2.D0  )/6.D0                                  
      FC2(NN)=0.5D0  *(DELN+1.D0  )*(DELN-1.D0  )*(DELN-2.D0  )                         
      FC3(NN)=-0.5D0  *(DELN+1.D0  )*DELN*(DELN-2.D0  )                               
      FC4(NN)=(DELN+1.D0  )*DELN*(DELN-1.D0  )/6.D0                                   
   60 CONTINUE                                                                  
      DO 500 JLS=1,JLSMAX                                                       
      DO 70 NN=NBMIN,NBP1                                                       
      FTEM(NN)=FORMFC(NN,JLS)                                                   
      N=NTIME1*NN                                                               
      FINAL(N)=FTEM(NN)                                                         
   70 CONTINUE                                                                  
      NINT=0                                                                    
      F1=FTEM(1)                                                                
      F2=FTEM(2)                                                                
      F3=FTEM(3)                                                                
      DO 80 NN=1,NM1                                                            
      DELN=DEL(NN)                                                              
      FAC1=0.5D0*(DELN-2.D0  )*(DELN-3.D0  )                                          
      FAC2=-(DELN-1.D0  )*(DELN-3.D0  )                                             
      FAC3=0.5D0*(DELN-1.D0  )*(DELN-2.D0  )                                          
      FINAL(NN)=FAC1*F1+FAC2*F2+FAC3*F3                                         
   80 CONTINUE                                                                  
      NINT=NTIME1                                                               
  110 DO 150 NN=1,NM1                                                           
      DELN=DEL(NN)                                                              
      FAC1=0.5D0*(DELN-1.D0  )*(DELN-2.D0  )                                          
      FAC2=-DELN*(DELN-2.D0  )                                                    
      FAC3=0.5D0*DELN*(DELN-1.D0  )                                                 
      FINAL(NN+NINT)=FAC1*F1+FAC2*F2+FAC3*F3                                    
  150 CONTINUE                                                                  
      NINT=NINT+NTIME1                                                          
      NP2=NBMIN+2                                                               
      DO 300 NB=NP2,NBINT                                                       
      F4=FTEM(NB+1)                                                             
      DO 250 NN=1,NM1                                                           
      FINAL(NN+NINT)=FC1(NN)*F1+FC2(NN)*F2+FC3(NN)*F3+FC4(NN)*F4                
  250 CONTINUE                                                                  
      F1=F2                                                                     
      F2=F3                                                                     
      F3=F4                                                                     
      NINT=NINT+NTIME1                                                          
  300 CONTINUE                                                                  
CCCCC        ASYMPTOTIC FORM FACTOR INTERPOLATION                               
      NBM2=NBINT-2                                                              
      DSIG= DSIGN(1.D0  ,FORMFC(NBP1,JLS))                                         
      DO 310 N=1,NBTAIL                                                         
  310 FTEM(N)=DLOG( DABS(FORMFC(NBM2+N,JLS)))                                    
      F1=FTEM(1)                                                                
      F2=FTEM(2)                                                                
      F3=FTEM(3)                                                                
      NTM1=NBTAIL-1                                                             
      DO 350 NB=3,NTM1                                                          
      F4=FTEM(NB+1)                                                             
      DO 340 NN=1,NM1                                                           
      FINAL(NN+NINT)=FC1(NN)*F1+FC2(NN)*F2+FC3(NN)*F3+FC4(NN)*F4                
  340 CONTINUE                                                                  
      F1=F2                                                                     
      F2=F3                                                                     
      F3=F4                                                                     
      NINT=NINT+NTIME1                                                          
      FINAL(NINT)=FTEM(NB)                                                      
  350 CONTINUE                                                                  
      NP1=NBINT*NTIME1+1                                                        
      NTOTAL=NBMAX*NTIME1                                                       
      DO 360 N=NP1,NTOTAL                                                       
  360 FINAL(N)= DEXP(FINAL(N))*DSIG                                              
      XMESA=XMESB/MESFCB                                                        
      WRITE(6,420)JLS,XMESA                                                     
  420 FORMAT(///,32H NO RECOIL FORM FACTOR FOR JLS =,I5,5X,3H AT,F7.3,          
     1   11H FERMI MESH,/)                                                      
      WRITE(6,430)(FINAL(N),N=1,NTOTAL)                                         
  430 FORMAT(10E13.5)                                                           
      IF(KPUNCH-1) 450,431,445                                                  
  431 NXREPT=(NTOTAL-1)/5+1                                                     
      DO 436 NXR=1,NXREPT                                                       
      NXI=5*(NXR-1)+1                                                           
      NXM=NXI+4                                                                 
      IF(NXR.EQ.NXREPT) NXM=NTOTAL                                              
C      PUNCH 434,(FINAL(N),N=NXI,NXM),NXR                                        
C  434 FORMAT(5E14.6,7X,I3)                                                      
  436 CONTINUE                                                                  
      GO TO 450                                                                 
  445 WRITE (12) NTOTAL,(FINAL(N),N=1,NTOTAL)                                   
      IF(KPUNCH.EQ.3) GO TO 431                                                 
  450 CONTINUE                                                                  
  500 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE KERNEL                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/RARB/NATOTL,NBMIN,NBMAX,NBREAD,NBTAIL,MESFCA,MESFCB,               
     1            N1MAX,N2MAX,MGAUS,NGAUS                   
      COMMON/RARB1/XMESA,XMESB,XMES(2)                   
      COMMON/GAUS/ABSCIS(222),WEIGHT(222),WWL1(500),WWL2(500),                  
     1            BSPAR(46),IBSPR(10)                                           
      COMMON/ST12/S1,S2,T1,T2                                                   
      COMMON/GKAB/RA,RB,WGT(96),RTS(96),XMAX(2),COSAB2                          
      COMMON GTAPE( 8192),GTEMP(4096),PK(250)                                   
      DIMENSION NUMG(12),IPOS(12),DE(4) ,RAS(2),RBT(2),RAD(2),W(2)              
      DATA NUMG / 4,8,12,16,20,24,32,40,48, 64, 80, 96 /                        
      DATA IPOS / 1,3, 7,13,21,31,43,59,79,103,135,175 /                        
      WRITE(6,3)                                                                
    3 FORMAT(1H1,30X,15(1H*),5X,29HFINITE RANGE KERNEL IS CALLED,5X,            
     1 15(1H*))                                                                 
      KOUT=KOUTST(5)                                                            
CCCCCC    ******   THE NUMBER 1024 CORRESPONDS TO 2K OCTAL BUFFER SIZE.         
      NBUFER=4*1024                                                             
      NGTPMX=8*1024                                                             
CCCCCC  ******    GAUSSIAN QUADRATURE ROOT AND WEIGHT FACTORS     ******        
      NGAUS=NUMG(MGAUS)                                                         
      IP=IPOS(MGAUS)                                                            
      NG2=NGAUS/2                                                               
      DO 40 NG=1,NG2                                                            
      RTS(NG2-NG+1)=-ABSCIS(IP)                                                 
      WGT(NG2-NG+1)= WEIGHT(IP)                                                 
      RTS(NG+NG2)=ABSCIS(IP)                                                    
      WGT(NG+NG2)=WEIGHT(IP)                                                    
   40 IP=IP+1                                                                   
CCCCCC  ******    PREPARE COORDINATE TRANSFORMATION COEFFICIENTS  ******        
      AC=MASCA                                                                  
      BC=MASCB                                                                  
      AS=MASSA                                                                  
      BS=MASSB                                                                  
      TM=AC+AS                                                                  
      XM= DABS(AC-BC)                                                            
      IF(MASCA.GT.MASCB) GO TO  10                                              
      S1=AS*BC/(XM*TM)                                                          
      S2=(AC/BC)*S1                                                             
      T1=-(BS/AS)*S1                                                            
      T2=-S1                                                                    
      GO TO 20                                                                  
   10 T1=BS*AC/(XM*TM)                                                          
      T2=(BC/AC)*T1                                                             
      S1=-(AS/BS)*T1                                                            
      S2=-T1                                                                    
   20 RZR=1.2D0                                                                   
      XBARSB=RZR*(BS**.333333333333  )                                          
      XMAX(1)=XMES(1)*(N1MAX-2)                                                 
      XMAX(2)=XMES(2)*(N2MAX-2)                                                 
      NBTOTL=NBMAX-NBMIN+1                                                      
      WRITE(6,31)NBMIN,NBMAX,NBREAD,NBTOTL,NATOTL,MESFCA,MESFCB,                
     1         KMAX,N1MAX,N2MAX,NGAUS,MGAUS,NGTPMX,NBUFER,                      
     2         XMESA,XMESB,XMES(1),XMES(2),                                     
     3         S1,S2,T1,T2,XBARSB                                               
   31 FORMAT(//,20X,47HNBMIN,NBMAX,NBREAD,NBTOTL,NATOTL,MESFCA,MESFCB=,         
     1      7I5/20X,47HKMAX,N1MAX,N2MAX,NGAUS,MGAUS,NGTPMX,NBUFER    =,         
     2   6I5,I7/20X,28HXMESA,XMESB,XMES(1),XMES(2)=,                            
     3    4F9.4/20X,28HS1,S2,T1,T2,XBARSB         =,5F9.4/)                     
  127 NATLHF=NATOTL/2                                                           
      NATRMI=NBMIN*XMESB/XMESA+0.1D0                                              
      NATRMX=NBMAX*XMESB/XMESA+0.1D0                                              
      NAGBLK=NATOTL*NGAUS                                                       
      NABBLK=NATOTL*NBTOTL                                                      
      NAKBLK=NATOTL*KMAX                                                        
      NTOTAL=NATOTL*NBTOTL*NGAUS                                                
      ITAP13=0                                                                  
      IF(NTOTAL.GT.NGTPMX) ITAP13=1                                             
      NBSTEP=NBUFER/NAGBLK                                                      
      IF(ITAP13.EQ.0) NBSTEP=NBTOTL                                             
      NBSTEP=MIN0(NBSTEP,NBTOTL)                                                
      NBREPT=(NBTOTL+NBSTEP-1)/NBSTEP                                           
      NBSTPF=NBTOTL-NBSTEP*(NBREPT-1)                                           
      WRITE(6,131)  ITAP13,     NBREPT,NBSTEP,NBSTPF                            
  131 FORMAT(20X,33HITAP13,     NBREPT,NBSTEP,NBSTPF=,5I5)                      
      KSTEP=NGTPMX/NABBLK                                                       
      IF(ITAP13.EQ.0) KSTEP=NBUFER/NABBLK                                       
      KSTEP=MIN0(KSTEP,KMAX)                                                    
      KREPT=(KMAX+KSTEP-1)/KSTEP                                                
      KSTPF=KMAX-KSTEP*(KREPT-1)                                                
      WRITE(6,133)KREPT,KSTEP,KSTPF                                             
  133 FORMAT(20X,18HKREPT,KSTEP,KSTPF=,3I5)                                     
      IF(NBSTEP.EQ.0.OR.KSTEP.EQ.0) WRITE(6,135)NBSTEP,KSTEP                    
      IF(NBSTEP.EQ.0.OR.KSTEP.EQ.0) STOP                                        
  135 FORMAT(//,1X,10(1H*),13HNBSTEP,KSTEP=,2I5,78HINCREASE DIMENSION OF        
     1 EITHER NBUFER(=8K OCTAL) OR NGTPMX(=16K OCTAL) IN KERNEL)                
      DO 500 NBREP=1,NBREPT                                                     
      NGTAPE=0                                                                  
      NBMINT=NBMIN+(NBREP-1)*NBSTEP                                             
      NBMAXT=NBMINT+NBSTEP-1                                                    
      NBMAXT=MIN0(NBMAX,NBMAXT)                                                 
      DO 450 NB=NBMINT,NBMAXT                                                   
      RB=NB*XMESB                                                               
      RA=RB                                                                     
      RAS(2)=RA*S2                                                              
      RBT(2)=RB*T2                                                              
      SUM2=RAS(2)**2+RBT(2)**2                                                  
      DENOM=2.D0  *RAS(2)*RBT(2)                                                  
      R2MIN= DABS(RAS(2)+RBT(2))                                                 
      R2MAX= DABS(RAS(2)-RBT(2))                                                 
      R2MP6=R2MIN+6.D0                                                            
      COSAB2  =(R2MP6   **2-SUM2)/DENOM                                         
      IF(R2MAX.LT.R2MP6   ) COSAB2  =-1.D0                                        
      IF(KOUT.LT.3) GO TO 141                                                   
      WRITE(6,138)RB                                                            
  138 FORMAT(10X,3HRB=,F7.4)                                                    
  141 RA=RB-NATLHF*XMESA-XMESA                                                  
      NATR=NB*XMESB/XMESA+0.1D0                                                   
      NATR=NATR-NATLHF                                                          
CCCCCC  ******            NA LOOP BEGINS FOR A GIVEN NB           ******        
      DO 400 NA=1,NATOTL                                                        
      RA=RA+XMESA                                                               
      IF(NATR.LT.NATRMI) GO TO 151                                              
      IF(NATR.GT.NATRMX) GO TO 151                                              
      CALL GKRARB                                                               
      GO TO 361                                                                 
  151 DO 157 NG=1,NGAUS                                                         
  157 GTEMP(NG)=0.0D0                                                             
  361 DO 365 NG=1,NGAUS                                                         
      NGTAPE=NGTAPE+1                                                           
      GTAPE(NGTAPE)=GTEMP(NG)                                                   
  365 CONTINUE                                                                  
      NATR=NATR+1                                                               
      IF(KOUT.LT.3) GO TO 400                                                   
      WRITE(6,372)RA,(GTEMP(NG),NG=1,NGAUS)                                     
  372 FORMAT(4H RA=,F6.3,(1X,10E11.3))                                          
  400 CONTINUE                                                                  
  450 CONTINUE                                                                  
      IF(ITAP13.EQ.1) WRITE(13) (GTAPE(N),N=1,NGTAPE)                           
  500 CONTINUE                                                                  
      KK=0                                                                      
      DO 900 KR=1,KREPT                                                         
      KSTEPT=KSTEP                                                              
      IF (KR.EQ.KREPT)  KSTEPT=KSTPF                                            
      KKI=KK+1                                                                  
      KKF=KK+KSTEPT                                                             
      KK=KKF                                                                    
      IF( ITAP13-1) 610,620,620                                                 
  610 DO 615 N=1,NBUFER                                                         
  615 GTEMP(N)=0.D0                                                               
      GO TO 630                                                                 
  620 DO 625 N=1,NGTPMX                                                         
  625 GTAPE(N)=0.D0                                                               
  630 NTAPE=NAGBLK*NBSTEP                                                       
      IF(ITAP13.EQ.1) REWIND 13                                                 
      NBRUN=0                                                                   
      RB=(NBMIN-1)*XMESB                                                        
      DO 800 NBREP=1,NBREPT                                                     
      NBST=NBSTEP                                                               
      IF(NBREP.EQ.NBREPT) NBST=NBSTPF                                           
      IF(NBREP.EQ.NBREPT) NTAPE=NAGBLK*NBSTPF                                   
      IF(ITAP13.EQ.1) READ(13) (GTEMP(N),N=1,NTAPE)                             
      DO 750 NB=1,NBST                                                          
      NBRUN=NBRUN+1                                                             
      RB=RB+XMESB                                                               
      RA=RB                                                                     
      RAS(2)=RA*S2                                                              
      RBT(2)=RB*T2                                                              
      SUM2=RAS(2)**2+RBT(2)**2                                                  
      DENOM=2.  *RAS(2)*RBT(2)                                                  
      R2MIN= DABS(RAS(2)+RBT(2))                                                 
      R2MAX= DABS(RAS(2)-RBT(2))                                                 
      R2MP6=R2MIN+6.D0                                                            
      COSAB2  =(R2MP6   **2-SUM2)/DENOM                                         
      IF(R2MAX.LT.R2MP6   ) COSAB2  =-1.D0                                        
      RANGE=1.0D0-COSAB2                                                          
      DO 700 NG=1,NGAUS                                                         
      ROOT=1.0D0-(1.D0  -RTS(NG))*RANGE/2.D0                                          
      PK(1)=1.D0                                                                  
      PK(2)=ROOT                                                                
      FK=1.0D0                                                                    
      TFK=1.0D0                                                                   
      DO 635 K=3,KKF                                                            
      FK=FK+1.D0                                                                  
      TFK=TFK+2.D0                                                                
      PK(K)=(TFK*ROOT*PK(K-1)-(FK-1.D0  )*PK(K-2))/FK                             
  635 CONTINUE                                                                  
      RA=RB-NATLHF*XMESA-XMESA                                                  
      DO 650 NA=1,NATOTL                                                        
      RA=RA+XMESA                                                               
      NT13=(NB-1)*NAGBLK+(NA-1)*NGAUS+NG                                        
      IF(ITAP13.EQ.0) W12=GTAPE(NT13)                                           
      IF(ITAP13.EQ.1) W12=GTEMP(NT13)                                           
      DO 640 KS=1,KSTEPT                                                        
      KP1=KKI+KS-1                                                              
      NT11=(KS-1)*NABBLK+(NBRUN-1)*NATOTL+NA                                    
      IF(ITAP13.EQ.0) GTEMP(NT11)=GTEMP(NT11)+W12*PK(KP1)                       
      IF(ITAP13.EQ.1) GTAPE(NT11)=GTAPE(NT11)+W12*PK(KP1)                       
  640 CONTINUE                                                                  
  650 CONTINUE                                                                  
  700 CONTINUE                                                                  
  750 CONTINUE                                                                  
  800 CONTINUE                                                                  
      NT1=1-NABBLK                                                              
      DO 850 KS=1,KSTEPT                                                        
      NT1=NT1+NABBLK                                                            
      NT2=NT1+NABBLK-1                                                          
      IF(              ITAP13.EQ.0) WRITE(11) (GTEMP(N),N=NT1,NT2)              
      IF(              ITAP13.EQ.1) WRITE(11) (GTAPE(N),N=NT1,NT2)              
  850 CONTINUE                                                                  
  851 IF(KOUT.EQ.0.AND.MOD(KR,5).NE.1) GO TO 900                                
      KPRINK=1                                                                  
      IF(KOUT.EQ.0) KPRINK=10                                                   
      KPRINA=1                                                                  
      IF(KOUT.LE.1) KPRINA=(NATOTL+8)/9                                         
      NBPR1=(NBMIN+9)/10                                                        
      NBPR2=NBMAX/10                                                            
      DO 870 NBPR=NBPR1,NBPR2                                                   
      RB=XMESB*NBPR*10                                                          
      NB1=NBPR*10-NBMIN                                                         
      NT1=NB1*NATOTL+1-NABBLK                                                   
      NT2=NT1+NATOTL                                                            
      DO 870 KS=1,KSTEPT                                                        
      K=KKI+KS-2                                                                
      NT1=NT1+NABBLK                                                            
      NT2=NT2+NABBLK                                                            
      IF(MOD(K,KPRINK).NE.0) GO TO 870                                          
      IF(ITAP13.EQ.0) WRITE(6,860)K,RB,(GTEMP(N),N=NT1,NT2,KPRINA)              
      IF(ITAP13.EQ.1) WRITE(6,860)K,RB,(GTAPE(N),N=NT1,NT2,KPRINA)              
  860 FORMAT(3H K=,I3,2X,3HRB=,F5.1,9E12.4)                                     
  870 CONTINUE                                                                  
  900 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE GKRARB                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/RARB/NATOTL,NBMIN,NBMAX,NBREAD,NBTAIL,MESFCA,MESFCB,               
     1            N1MAX,N2MAX,MGAUS,NGAUS                   
      COMMON/RARB1/XMESA,XMESB,XMES(2)                   
      COMMON/GAUS/ABSCIS(222),WEIGHT(222),WWL1(500),WWL2(500),                  
     1            BSPAR(46),IBSPR(10)                                           
      COMMON/ST12/S1,S2,T1,T2                                                   
      COMMON/GKAB/RA,RB,WGT(96),RTS(96),XMAX(2),COSAB2                          
      COMMON GTAPE( 8192),GTEMP(4096),PK(250)                                   
      DIMENSION RAS(2),RBT(2),RAD(2),W(2),DE(4)                                 
CCCCCC   ******   EVALUATE THE NGAUS INTEGRAND OF THE KERNEL FUNCTION **        
CCCCCC            GK FOR A GIVEN (RA,RB) IN THE GAUSSIAN QUADRATURE.  **        
      RAS(1)=RA*S1                                                              
      RAS(2)=RA*S2                                                              
      RBT(1)=RB*T1                                                              
      RBT(2)=RB*T2                                                              
      DO 85 NG=1,NGAUS                                                          
   85 GTEMP(NG)=0.0D0                                                             
      RANGE=1.0D0-COSAB2                                                          
      DO 290 NG=1,NGAUS                                                         
      ROOT=1.0D0-(1.0D0 -RTS(NG))*RANGE/2.D0                                          
      WFAC=WGT(NG)*RANGE/2.D0                                                     
      DO 150 I=1,2                                                              
      RAD(I)= DSQRT(RAS(I)**2+RBT(I)**2+2.D0  *RAS(I)*RBT(I)*ROOT)                 
      IF(RAD(I).GT.XMAX(I)) GO TO 290                                           
CCCCCC  ******    INTERPOLATION OF BOUND STATE WAVE FUNCTION      ******        
CCCCCC  ******    USING 2 OR 4 POINT LAGRANGIAN METHOD.           ******        
      IR2=RAD(I)/XMES(I)                                                        
      IF(IR2.GE.2) GO TO 115                                                    
      DELX=RAD(I)/XMES(I)                                                       
      DEL1=DELX-1.D0                                                              
      DEL2=DELX-2.D0                                                              
      DE(1)= DEL1*DEL2*0.5D0                                                      
      DE(2)=-DELX*DEL2                                                          
      DE(3)= DELX*DEL1*0.5D0                                                      
      W(I)=0.D0                                                                   
      IF(I.EQ.2) GO TO 112                                                      
      DO 111 IRR=1,3                                                            
  111 W(I)=W(I)+DE(IRR)*WWL1(IRR)                                               
      GO TO 150                                                                 
  112 DO 113 IRR=1,3                                                            
  113 W(I)=W(I)+DE(IRR)*WWL2(IRR)                                               
      GO TO 150                                                                 
  115 IR2M2=IR2-2                                                               
      DELX=RAD(I)/XMES(I)-IR2                                                   
      DEP1=DELX+1.D0                                                              
      DEM1=DELX-1.D0                                                              
      DEM2=DELX-2.D0                                                              
      DEM12=DEM1*DEM2/6.D0                                                        
      DEP1X=DEP1*DELX/6.D0                                                        
      DE(1)=-DELX*DEM12                                                         
      DE(2)= DEP1*DEM12*3.D0                                                      
      DE(3)=-DEM2*DEP1X*3.D0                                                      
      DE(4)= DEM1*DEP1X                                                         
      W(I)=0.D0                                                                   
      IF(I.EQ.2) GO TO 122                                                      
      DO 121 IRR=1,4                                                            
      IR=IR2M2+IRR                                                              
  121 W(I)=W(I)+DE(IRR)*WWL1(IR)                                                
      GO TO 150                                                                 
  122 DO 123 IRR=1,4                                                            
      IR=IR2M2+IRR                                                              
  123 W(I)=W(I)+DE(IRR)*WWL2(IR)                                                
  150 CONTINUE                                                                  
      GTEMP(NG)=WFAC*W(1)*W(2)                                                  
  290 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE FRFMFC                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/RARB/NATOTL,NBMIN,NBMAX,NBREAD,NBTAIL,MESFCA,MESFCB,               
     1            N1MAX,N2MAX,MGAUS,NGAUS                   
      COMMON/RARB1/XMESA,XMESB,XMES(2)                   
      COMMON/FRFM/LA,LAMIN,KBASE,NGAB,KRANGE
      COMMON RAPOWR(3600),RBPOWR(500),               
     1       GEOK(2000),FTEM1(500),FTAPE(2500),GTAPE(4000)                      
CCCCCC  ***  FOLLOWING QUANTITIES DECIDE REQUIRED DIMENSIONS OF FTAPE,GT        
CCCCCC      LBSTTL=TOTAL NUMBER OF POSSIBLE LB FOR A GIVEN LA FOR ALL JL        
CCCCCC      KRANGE=2*(L1TR+L2TR)+1.                                             
CCCCCC      NGAB=NATOTL*NBTOTL                                                  
CCCCCC      NGUSE=NGAB*KRANGE           ( THIS MUST BE LT. NGTPMX ).            
CCCCCC      NFBLOK=LBSTTL*NBREAD*NATOTL ( THIS MUST BE LT. NFTPMX ).            
      KOUT=KOUTST(6)                                                            
      KPRINT=1                                                                  
      IF(KOUT.EQ.0) KPRINT=100                                                  
      WRITE(6,15)                                                               
   15 FORMAT(1H1,36X,20(1H*),20H  FRFMFC IS CALLED  ,20(1H*))                   
      REWIND 11                                                                 
      NBRED3=NBREAD-3                                                           
      NATLHF=(NATOTL-1)/2                                                       
      NBTOTL=(NBMAX-NBMIN)+1                                                    
      KRANGE=L1TRTW(1)+L2TRTW(1)+1                                              
      KRANHF=(KRANGE-1)/2                                                       
      NGAB=NATOTL*NBTOTL                                                        
      NGUSE=NGAB*KRANGE                                                         
      NABLOK=LBSTTL*NATOTL                                                      
      NFBLOK=NABLOK*NBREAD                                                      
      NGBLOK=KMAX                                                               
      WRITE(6,31)NBMIN,NBMAX,NBREAD,NBTOTL,NATOTL,MESFCA,MESFCB,                
     1         NFBLOK,NGBLOK,KRANGE,NGAB,NGUSE,KMAX ,                           
     2         LBSTTL,              XMESA,XMESB                                 
   31 FORMAT(//,20X,47HNBMIN,NBMAX,NBREAD,NBTOTL,NATOTL,MESFCA,MESFCB=,         
     1    7I5,/,20X,47HNFBLOK,NGBLOK,KRANGE,NGAB,NGUSE,KMAX          =,         
     2    6I5,/,20X,32HLBSTTL               XMESA,XMESB ,14X,1H=, I5,           
     3      2F7.3,//)                                                           
CCCCCC  ******        CFCCAL IS CALLED TO CALCULATE GEOMETRICAL   ******        
CCCCCC  ******        FACTORS DFAC AND CFAC.                      ******        
      CALL CFCCAL                                                               
      LAMIN=1                                                                   
      LAMAX=LDWMAX                                                              
      LASTEP=1                                                                  
      IF(KEXCST(2).NE.0) LASTEP=KEXCST(2)                                       
      IF(KEXCST(2).NE.0) LAMIN=1+MOD(KPCHAG,2)                                  
      IF(KEXCST(1).NE.0) LAMIN=KEXCST(1)+MOD(KEXCST(1)+KPCHAG+1,2)              
      KRMAX=LAMIN-(KRANHF+1)                                                    
      IF(KRMAX.LT.1) GO TO 276                                                  
      DO 275 K=1,KRMAX                                                          
      READ(11) (GTAPE(N),N=1,NGAB)                                              
  275 CONTINUE                                                                  
  276 N1=1-NGAB                                                                 
      DO 278 K=1,KRANGE                                                         
      N1=N1+NGAB                                                                
      N2=N1+NGAB-1                                                              
      READ(11) (GTAPE(N),N=N1,N2)                                               
  278 CONTINUE                                                                  
      KBASE=0                                                                   
      IF(KRMAX.GE.1) KBASE=KRMAX                                                
      IF(KOUT.LT.4) GO TO 305                                                   
      WRITE(6,281)KRMAX,KBASE,(GTAPE(N),N=1,NGAB)                               
  281 FORMAT(19H KRMAX,KBASE,GTAPE=,2I5,/,(1X,10E12.5))                         
  305 DO 800 LAP1=LAMIN,LAMAX,LASTEP                                            
      LA=LAP1-1                                                                 
      LATW=2*LA                                                                 
      IF(LAP1.EQ.LAMIN) GO TO 331                                               
      IF(LAP1.LE.KRANHF+1) GO TO 331                                            
      DO 320 LREAD=1,LASTEP                                                     
      NGABM1=(KRANGE-1)*NGAB                                                    
CCCCC     NOTE  ---   CONSIDER THE CASE WITH KRANGE=1                           
      DO 316 N=1,NGABM1                                                         
  316 GTAPE(N)=GTAPE(N+NGAB)                                                    
      KBASE=KBASE+1                                                             
      NGABM2=NGABM1+NGAB                                                        
      NGABM1=NGABM1+1                                                           
      READ(11) (GTAPE(N),N=NGABM1,NGABM2)                                       
  320 CONTINUE                                                                  
CCCCCC  ******            FOR EACH LA, FFRARB IS CALLED.          ******        
  331 CALL FFRARB                                                               
  800 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CFCCAL                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC,L9(10),U9                   
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/CFCC/CFAC(1500)                                                    
      COMMON/ST12/S1,S2,T1,T2                                                   
      DIMENSION SUM(20)                                                         
      KOUT=KOUTST(6)                                                            
      IXTW=MOD(JJTWR(1),2)                                                      
      L1TR=L1TRTW(1)/2                                                          
      L2TR=L2TRTW(1)/2                                                          
      LD12PM=L1TR+L2TR+1                                                        
      NGRNMX=KEXCST(15)                                                         
      NN=LD12PM*JLSMAX*NGRNMX                                                   
      DO 125 N=1,NN                                                             
  125 CFAC(N)=0.0D0                                                               
      L12MOD=MOD(L1TR+L2TR,2)                                                   
      L1TW=2*L1TR                                                               
      L2TW=2*L2TR                                                               
      L1P1=L1TR+1                                                               
      L2P1=L2TR+1                                                               
      L12TR=L1TR+L2TR                                                           
      LDARMX=L12TR+1                                                            
      DO 300 JLS=1,JLSMAX                                                       
      NRUNG=0                                                                   
      LLTW=LLTWR(JLS)                                                           
      LL=LLTW/2                                                                 
      JJTW=JJTWR(JLS)                                                           
      ISTW=ISTWR(JLS)                                                           
CCCCCC  ******    CALCULATION OF DFAC WITH UNIT SPEC. FACTOR.     ******        
      K1=IABS(ISTW+L2TW-IXTW)/2                                                 
      SSD=1-2*MOD(K1,2)                                                         
      IA=L1TW                                                                   
      IB=L2TW                                                                   
      IC=JJTW                                                                   
      ID=ISTW                                                                   
      IE=LLTW                                                                   
      IFF=IXTW                                                                   
      CALL RAC7                                                                 
      DCOE(JLS)=SSD*RAC                                                         
CCCCCC  ******         CALCULATION OF GEOMETRICAL FACTOR CFAC     ******        
CCCCCC  ******                     LAMBDA A LOOP                  ******        
      DO 280 LDAR=1,LDARMX                                                      
      LDA=LDAR-1                                                                
      LDATW=2*LDA                                                               
      LDBRMI=IABS(LDA-LL)+1                                                     
      LDBRMX=MIN0(LDA+LL,LDARMX-1)+1                                            
      IF(LDBRMI.GT.LDBRMX) GO TO 280                                            
CCCCCC  ******                     LAMBDA B LOOP                  ******        
      DO 270 LDBR=LDBRMI,LDBRMX                                                 
      IF(MOD(LDAR+LDBR,2).NE.L12MOD) GO TO 270                                  
      NRUNG=NRUNG+1                                                             
      LDSTBS=LD12PM*((JLS-1)*NGRNMX+(NRUNG-1))                                  
      LDB=LDBR-1                                                                
      LDBTW=2*LDB                                                               
      DO 220 N=1,LD12PM                                                         
  220 SUM(N)=0.D0                                                                 
CCCCCC  ******                     LAMBDA 1 LOOP                  ******        
      DO 260 LD1R=1,L1P1                                                        
      LD1=LD1R-1                                                                
      LD1P=L1TR-LD1                                                             
      LD1TW=2*LD1                                                               
      LD1PTW=2*LD1P                                                             
      CB1= DEXP(0.5D0*(FACLOG(L1TW+2)-FACLOG(LD1TW+2)-FACLOG(LD1PTW+2)))         
     1   *(S1**LD1)*(T1**LD1P)                                                  
      LD2RMI=IABS(LDA-LD1)+1                                                    
      LD2RMX=MIN0(L2TR,LDA+LD1)+1                                               
      IF(LD2RMI.GT.LD2RMX) GO TO 260                                            
CCCCCC  ******                     LAMBDA 2 LOOP                  ******        
      DO 250 LD2R=LD2RMI,LD2RMX                                                 
      LD2=LD2R-1                                                                
      LD2P=L2TR-LD2                                                             
      IF(MOD(LDA+LD1+LD2,2).NE.0) GO TO 250                                     
      IF(LDB-IABS(LD1P-LD2P).LT.0) GO TO 250                                    
      IF(LD1P+LD2P-LDB.LT.0) GO TO 250                                          
      LD2TW=2*LD2                                                               
      LD2PTW=2*LD2P                                                             
      CB2= DEXP(0.5D0*(FACLOG(L2TW+2)-FACLOG(LD2TW+2)-FACLOG(LD2PTW+2)))         
     1   *(S2**LD2)*(T2**LD2P)                                                  
      IA=LD1TW                                                                  
      IB=LD2TW                                                                  
      IC=LDATW                                                                  
      ID=0                                                                      
      IE=0                                                                      
      IFF=0                                                                      
      CALL CLEB                                                                 
      C1=RAC                                                                    
      IF(C1.EQ.0.0D0  ) GO TO 250                                                 
      IA=LD1PTW                                                                 
      IB=LD2PTW                                                                 
      IC=LDBTW                                                                  
      ID=0                                                                      
      IE=0                                                                      
      IFF=0                                                                      
      CALL CLEB                                                                 
      C2=RAC                                                                    
      IF(C2.EQ.0.0D0  ) GO TO 250                                                 
      L9(1)=LD1TW                                                               
      L9(2)=LD2TW                                                               
      L9(3)=LD1PTW                                                              
      L9(4)=LD2PTW                                                              
      L9(5)=LDATW                                                               
      L9(6)=LDBTW                                                               
      L9(7)=L1TW                                                                
      L9(8)=L2TW                                                                
      L9(9)=LLTW                                                                
      CALL NINEJ                                                                
      U91=U9                                                                    
      IF(U91.EQ.0.0D0  ) GO TO 250                                                
      HAT=(LD1TW+1)*(LD2TW+1)*(LD1PTW+1)*(LD2PTW+1)*(L1TW+1)*(L2TW+1)           
     1   *(LDATW+1)*(LDBTW+1)                                                   
      HAT= DSQRT(HAT)                                                            
      TS=CB1*CB2*HAT*C1*C2*U91*DCOE(JLS)                                        
      LD12P1=LD1+LD2+1                                                          
      SUM(LD12P1)=SUM(LD12P1)+TS                                                
      IF(KOUTST(6).LT.1) GO TO 250                                              
      IF(NRUNG.GT.30) GO TO 250                                                 
      WRITE(6,242) JLS,L1TR,L2TR,LDA,LDB,LD1,LD2,LD1P,LD2P                      
  242 FORMAT(/,52H CFCCAL-242.JLS,L1TR,L2TR,LDA,LDB,LD1,LD2,LD1P,LD2P=,         
     1       10I4)                                                              
      WRITE(6,244) CB1,CB2,HAT,C1,C2,U91,TS                                     
  244 FORMAT(26H CB1,CB2,HAT,C1,C2,U91,TS=,7E14.5)                              
  250 CONTINUE                                                                  
  260 CONTINUE                                                                  
      DO 265 LD12=1,LD12PM                                                      
      LDST=LDSTBS+LD12                                                          
  265 CFAC(LDST)=SUM(LD12)                                                      
CCCCCC  ******               OUTPUT OF CFAC AND DFAC              ******        
      IF(KOUT.LT.1) GO TO 270                                                   
      NNMI=LDSTBS+1                                                             
      NNMX=LDSTBS+LD12PM                                                        
      WRITE(6,267)(CFAC(NN),NN=NNMI,NNMX)                                       
  267 FORMAT(6H CFAC=,10E12.4)                                                  
  270 CONTINUE                                                                  
  280 CONTINUE                                                                  
  300 CONTINUE                                                                  
      WRITE(6,310) (DCOE(JLS),JLS=1,JLSMAX)                                     
  310 FORMAT(//,6H DCOE=,8E14.6)                                                
      RETURN                                                                    
      END                                                                       
      SUBROUTINE FFRARB                                                         
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC,L9(10),U9                   
      COMMON/SATN/KTRLST(15),KEXCST(15),KOUTST(15),                             
     1            MASCA,MASCB,MASSA,MASSB,ICATW,ICBTW,ISATW,ISBTW,              
     2            NZCA,NZCB,NZSA,NZSB,KPCA,KPCB,KPSA,KPSB,KPCHAG,               
     3            JLSMAX,JJTWR(5),LLTWR(5),ISTWR(5),L1TRTW(5),L2TRTW(5),        
     4            LDWMAX,KMAX,LLMAX,LBSTTL,DCOE(5)                              
      COMMON/RARB/NATOTL,NBMIN,NBMAX,NBREAD,NBTAIL,MESFCA,MESFCB,               
     1            N1MAX,N2MAX,MGAUS,NGAUS                   
      COMMON/RARB1/XMESA,XMESB,XMES(2)                   
      COMMON/CFCC/CFAC(1500)                                                    
      COMMON/FRFM/LA,LAMIN,KBASE,NGAB,KRANGE
      COMMON RAPOWR(3600),RBPOWR(500),               
     1       GEOK(2000),FTEM1(500),FTAPE(2500),GTAPE(4000)                      
      DIMENSION LLMEM(40),LBMEM(40),LBJLME(40),LLOUT(40),LBOUT(40),             
     1          LBJLMM(5),LBMINM(5),LBMAXM(5)                                   
      KOUT=KOUTST(6)                                                            
      LAP1=LA+1                                                                 
      IF(LAP1.NE.LAMIN) GO TO 301                                               
      NBTOTL=(NBMAX-NBMIN)+1                                                    
      NATLHF=(NATOTL-1)/2                                                       
      LDNB=0                                                                    
      NANB=0                                                                    
      L1TR=L1TRTW(1)/2                                                          
      L2TR=L2TRTW(1)/2                                                          
      L12TR=L1TR+L2TR                                                           
      L12MOD=MOD(L12TR,2)                                                       
      LD12PM=L12TR+1                                                            
      NGRNMX=KEXCST(15)                                                         
CCCCCC  ******            CALCULATE POWERS OF RB AND RA           ******        
      DO 250 NB=NBMIN,NBMAX                                                     
      RB=NB*XMESB                                                               
      RBP=RB**LD12PM                                                            
      DO 240 LD12P1=1,LD12PM                                                    
      RBP=RBP/RB                                                                
      LDNB=LDNB+1                                                               
      RBPOWR(LDNB)=RBP                                                          
  240 CONTINUE                                                                  
  250 CONTINUE                                                                  
      NARMIN=NBMIN*XMESB/XMESA+0.1D0                                              
      NARMIN=NARMIN-NATLHF                                                      
      NARMAX=NBMAX*XMESB/XMESA+0.1D0+NATLHF                                       
      NARTOT=NARMAX-NARMIN+1                                                    
      RAMM=XMESA *(NARMIN-1)                                                    
      RA=RAMM                                                                   
      NARZ=0                                                                    
      DO 255 NAR=1,NARTOT                                                       
      RA=RA+XMESA                                                               
      IF(DABS(RA).LT.0.00001D0) RA=0.0D0                                             
      IF(RA.EQ.0.0D0) NARZ=NAR                                                    
      RAPOWR(NAR)=1.0D0                                                           
      IF(RA.LT.0.0D0) RAPOWR(NAR)=0.0D0                                             
  255 CONTINUE                                                                  
      IF(LD12PM.EQ.1) GO TO 301                                                 
      LDNA=NARTOT                                                               
      DO 260 LD12P1=2,LD12PM                                                    
      RA=RAMM                                                                   
      DO 257 NAR=1,NARTOT                                                       
      RA=RA+XMESA                                                               
      LDNA=LDNA+1                                                               
      RAPOWR(LDNA)=RAPOWR(LDNA-NARTOT)*RA                                       
  257 CONTINUE                                                                  
      IF(LD12P1.NE.2) GO TO 260                                                 
      IF(NARZ.NE.0) RAPOWR(NARTOT+NARZ)=0.0D0                                     
  260 CONTINUE                                                                  
      IF(KOUT.LT.4) GO TO 301                                                   
      WRITE(6,271) (RAPOWR(NA),NA=1,LDNA)                                       
      WRITE(6,273) (RBPOWR(NB),NB=1,LDNB)                                       
  271 FORMAT(8H RAPOWR= ,10E12.4)                                               
  273 FORMAT(8H RBPOWR= ,10E12.4)                                               
CCCCCC  ******         REARRANGING OF LB AND LL FOR A GIVEN LA    ******        
  301 LATW=2*LA                                                                 
      LBJLMM(1)=0                                                               
      LLMIN=100                                                                 
      LLMAX=0                                                                   
      LBP1MI=1000                                                               
      LBP1MX=0                                                                  
      LBJL=0                                                                    
      DO 310 JLS=1,JLSMAX                                                       
      LL=LLTWR(JLS)/2                                                           
      LBMIN=IABS(LA-LL)+1+MOD(KPCHAG+LL,2)                                      
      LBMAX=LA+LL+1                                                             
      LBMINM(JLS)=LBMIN                                                         
      LBMAXM(JLS)=LBMAX                                                         
      IF(LBMIN.GT.LBMAX) GO TO 306                                              
      IF(LL.LT.LLMIN) LLMIN=LL                                                  
      IF(LL.GT.LLMAX) LLMAX=LL                                                  
      IF(LBMIN.LT.LBP1MI) LBP1MI=LBMIN                                          
      IF(LBMAX.GT.LBP1MX) LBP1MX=LBMAX                                          
      DO 305 LBP1=LBMIN,LBMAX,2                                                 
      LBJL=LBJL+1                                                               
      LLMEM(LBJL)=LL                                                            
      LBMEM(LBJL)=LBP1-1                                                        
  305 CONTINUE                                                                  
  306 IF(JLS.NE.JLSMAX) LBJLMM(JLS+1)=LBJLMM(JLS)+(LBMAX-LBMIN)/2+1             
      IF(JLS.NE.JLSMAX.AND.LBMIN.GT.LBMAX) LBJLMM(JLS+1)=LBJLMM(JLS)            
      IF(KOUT.LT.4) GO TO 310                                                   
      WRITE(6,307)JLS,LL,KPCHAG,LBMIN,LBMAX,LBJL                                
  307 FORMAT(32H JLS,LL,KPCHAG,LBMIN,LBMAX,LBJL=,6I5)                           
  310 CONTINUE                                                                  
      LBJLMX=LBJL                                                               
      LBP1RG=(LBP1MX-LBP1MI)+1                                                  
      LLP1RG=(LLMAX-LLMIN)+1                                                    
      LBJLS=0                                                                   
      LB=LBP1MX                                                                 
      DO 340 LBP1=1,LBP1RG                                                      
      LB=LB-1                                                                   
      LL=LLMAX+1                                                                
      DO 330 LLP1=1,LLP1RG                                                      
      LL=LL-1                                                                   
      DO 325 LBJL=1,LBJLMX                                                      
      IF(LBMEM(LBJL).NE.LB) GO TO 325                                           
      IF(LLMEM(LBJL).NE.LL) GO TO 325                                           
      LBJLS=LBJLS+1                                                             
      LBJLME(LBJLS)=LBJL                                                        
      LBOUT(LBJLS)=LB                                                           
      LLOUT(LBJLS)=LL                                                           
      GO TO 330                                                                 
  325 CONTINUE                                                                  
  330 CONTINUE                                                                  
  340 CONTINUE                                                                  
      NMAX=LD12PM*LBJLMX                                                        
      NNMAX=NMAX*KRANGE                                                         
CCCCCC  ******       CALCULATION OF GEOMETRICAL FACTOR GEOK       ******        
  353 DO 355 NN=1,NNMAX                                                         
  355 GEOK(NN)=0.0D0                                                              
      DO 500 JLS=1,JLSMAX                                                       
      NGRUN=0                                                                   
      LLTW=LLTWR(JLS)                                                           
      LL=LLTW/2                                                                 
      LBMIN=LBMINM(JLS)                                                         
      LBMAX=LBMAXM(JLS)                                                         
      IF(LBMIN.GT.LBMAX) GO TO 500                                              
      LDARMX=LD12PM                                                             
CCCCCC  ******                      LAMBDA A LOOP                 ******        
      DO 480 LDAR=1,LDARMX                                                      
      LDA=LDAR-1                                                                
      LDATW=2*LDA                                                               
      LDBRMI=IABS(LDA-LL)+1                                                     
      LDBRMX=MIN0(LDA+LL,LDARMX-1)+1                                            
      IF(LDBRMI.GT.LDBRMX) GO TO 480                                            
CCCCCC  ******                      LAMBDA B LOOP                 ******        
      DO 470 LDBR=LDBRMI,LDBRMX                                                 
      IF(MOD(LDAR+LDBR,2).NE.L12MOD) GO TO 470                                  
      NGRUN=NGRUN+1                                                             
      LDSTBS=LD12PM*((JLS-1)*NGRNMX+(NGRUN-1))                                  
      LDB=LDBR-1                                                                
      LDBTW=2*LDB                                                               
      LBJL=LBJLMM(JLS)                                                          
CCCCCC  ******                         LB LOOP                    ******        
      DO 460 LBP1=LBMIN,LBMAX,2                                                 
      LBJL=LBJL+1                                                               
      LDLBBS=(LBJL-1)*LD12PM                                                    
      LB=LBP1-1                                                                 
      LBTW=2*LB                                                                 
      LCHANG=IABS(L1TR+L2TR-LA-LB)/2                                            
      KMIN=MAX0(IABS(LDA-LA),IABS(LDB-LB))+1                                    
      KMMX=MIN0(LDA+LA,LDB+LB)+1                                                
      IF(KMIN.GT.KMMX) GO TO 460                                                
      IF(KMIN-KBASE.LT.0) GO TO 423                                             
      IF(KMMX-KBASE-1.GT.KRANGE) GO TO 423                                      
      GO TO 426                                                                 
  423 WRITE(6,424)LA,LB,LDA,LDB,JLS,KBASE,KMIN,KMMX                             
  424 FORMAT(//,1X,10(1H*),38HFOR LA,LB,LDA,LDB,JLS,KBASE,KMIN,KMMX=,           
     1  8I4,20H RANGE OF K IS WRONG)                                            
      STOP                                                                      
  426 IF(LDA+LDB.GT.L12TR) GO TO 460                                            
CCCCCC   ******                       K LOOP                      ******        
      DO 450 KP1=KMIN,KMMX,2                                                    
      KRUN=KP1-KBASE                                                            
      NGOKBS=(KRUN-1)*NMAX+LDLBBS                                               
      K=KP1-1                                                                   
      KTW=2*K                                                                   
      IA=LDATW                                                                  
      IB=KTW                                                                    
      IC=LATW                                                                   
      ID=0                                                                      
      IE=0                                                                      
      IFF=0                                                                      
      CALL CLEBZ                                                                
      C3=RAC                                                                    
      IF(C3.EQ.0.0D0  ) GO TO 450                                                 
      IA=LDBTW                                                                  
      IC=LBTW                                                                   
      CALL CLEBZ                                                                
      C4=RAC                                                                    
      IF(C4.EQ.0.0D0  ) GO TO 450                                                 
      IA=LATW                                                                   
      IB=LDATW                                                                  
      IC=LBTW                                                                   
      ID=LDBTW                                                                  
      IE=KTW                                                                    
      IFF=LLTW                                                                   
      CALL RAC7                                                                 
      R1=RAC                                                                    
      IF(R1.EQ.0.0D0  ) GO TO 450                                                 
      K1=K+LL+LCHANG                                                            
      SS=1-2*MOD(K1,2)                                                          
      TS=SS*C3*C4*R1*(KTW+1)*0.5D0                                                
      DO 430 LD12P1=1,LD12PM                                                    
      LDST=LDSTBS+LD12P1                                                        
      CCC=CFAC(LDST)*TS                                                         
      NGOK=NGOKBS+LD12P1                                                        
      GEOK(NGOK)=GEOK(NGOK)+CCC                                                 
  430 CONTINUE                                                                  
  450 CONTINUE                                                                  
  460 CONTINUE                                                                  
  470 CONTINUE                                                                  
  480 CONTINUE                                                                  
  490 CONTINUE                                                                  
  500 CONTINUE                                                                  
      IF(KOUT.LT.4) GO TO 515                                                   
      WRITE(6,505)(GEOK(NN),NN=1,NNMAX)                                         
  505 FORMAT((6H GEOK=,10E12.4))                                                
CCCCCC  ******    END OF PREPARATION OF GEOMETRICAL FACTORS       ******        
  515 NBRED3=NBREAD-3                                                           
      NABLOK=LBSTTL*NATOTL                                                      
CCCCCC  ******              NB LOOP BROKEN INTO NRR TIMES         ******        
      NRR=(NBTOTL-4)/NBRED3+1                                                   
      DO 800 NR=1,NRR                                                           
      IF(NR.NE.1) GO TO 551                                                     
      NBLMCM=0                                                                  
      NANBMX=NATOTL*MIN0 (NBREAD,NBTOTL)                                        
      NBMINT=NBMIN                                                              
      NBMAXT=MIN0(NBMAX,NBREAD+NBMIN-1)                                         
      NFBLOK=NABLOK*MIN0(NBREAD,NBTOTL)                                         
      GO TO 553                                                                 
  551 NANBMX=NATOTL*MIN0(NBRED3,(NBTOTL-NBLMCM))                                
      NBMINT=NBMAXT+1                                                           
      NBMAXT=MIN0(NBMAXT+NBRED3,NBMAX)                                          
      NFBLOK=NABLOK*MIN0(NBRED3,(NBTOTL-NBLMCM))                                
  553 NBLMCM=NBLMCM+(NFBLOK/NABLOK)                                             
      DO 555 N=1,NFBLOK                                                         
  555 FTAPE(N)=0.0D0                                                              
CCCCCC  ******   LB,K,LAMBDA 1 AND 2, NB PARTIAL AND NA LOOPS BEGIN  ***        
      DO 700 LBJLS=1,LBJLMX                                                     
      LBJL=LBJLME(LBJLS)                                                        
      LDLBBS=(LBJL-1)*LD12PM                                                    
      DO 650 KRUN=1,KRANGE                                                      
      NGOKBS=(KRUN-1)*NMAX+LDLBBS                                               
      DO 607 NANB=1,NANBMX                                                      
  607 FTEM1(NANB)=0.D0                                                            
      DO 640 LD12P1=1,LD12PM                                                    
      NGOK=NGOKBS+LD12P1                                                        
      G1=GEOK(NGOK)                                                             
      IF(G1.EQ.0.D0  ) GO TO 640                                                  
      NANB=0                                                                    
      LDNB=LD12P1+(NBMINT-NBMIN-1)*LD12PM                                       
      LDNABS=(LD12P1-1)*NARTOT+(NBMINT-NBMIN-1)*MESFCB/MESFCA                   
      DO 630 NB=NBMINT,NBMAXT                                                   
      LDNB=LDNB+LD12PM                                                          
      G2=G1*RBPOWR(LDNB)                                                        
      LDNABS=LDNABS+MESFCB/MESFCA                                               
      DO 630 NA=1,NATOTL                                                        
      NANB=NANB+1                                                               
      LDNA=LDNABS+NA                                                            
      FTEM1(NANB)=FTEM1(NANB)+G2*RAPOWR(LDNA)                                   
  630 CONTINUE                                                                  
  640 CONTINUE                                                                  
      NFTP=LBJLS-LBSTTL                                                         
      NGLOCT=(KRUN-1)*NGAB+NATOTL*(NBMINT-NBMIN)                                
      DO 645 NANB=1,NANBMX                                                      
      NFTP=NFTP+LBSTTL                                                          
      NGLOCT=NGLOCT+1                                                           
      FTAPE(NFTP)=FTAPE(NFTP)+FTEM1(NANB)*GTAPE(NGLOCT)                         
  645 CONTINUE                                                                  
  650 CONTINUE                                                                  
  700 CONTINUE                                                                  
CCCCCC  ******    WRITE ON TAPE12 AND OUTPUT OF FORM FACTOR.      ******        
      WRITE(12) (FTAPE(N),N=1,NFBLOK)                                           
      LASKIP=LDWMAX/5+5-MOD(LDWMAX/5,5)                                         
      IF(KOUTST(1).NE.0) LASKIP=KOUTST(1)                                       
      IF(KOUT.NE.0) LASKIP=1                                                    
      IF(MOD(LAP1-LAMIN,LASKIP).NE.0) GO TO 800                                 
      NBSKIP=5                                                                  
      IF(KOUT.NE.0) NBSKIP=1                                                    
      NW3=(NATOTL+10)/11                                                        
      IF(KOUT.NE.0) NW3=1                                                       
      NAW3=LBSTTL*NW3                                                           
      IF(NR.NE.1) GO TO 729                                                     
      IF(LAP1.EQ.LAMIN) WRITE(6,725)NW3                                         
      WRITE(6,721)LA,(LBOUT(L),LLOUT(L),L=1,LBJLMX)                             
  721 FORMAT(/10X,3HLA=,I5,8X,7H(LB,L)=,(1X,10(1H(,I3,1H,,I2,1H))))             
  725 FORMAT(20X,52HOUTPUT OF FORM FACTOR FOR A GIVEN RB AND (LB,L) FROM        
     1,27H NA=1 TO NATOTL AT STEPS OF,I3)                                       
  729 DO 740 NB=NBMINT,NBMAXT                                                   
      MODNB=MOD((NB-NBMIN),NBSKIP)                                              
      IF(MODNB.NE.0) GO TO 740                                                  
      RB=NB*XMESB                                                               
      DO 735 LBST=1,LBSTTL                                                      
      NAW1=(NB-NBMINT)*NABLOK+LBST                                              
      NAW2=NAW1+(NABLOK-LBSTTL)                                                 
      WRITE(6,731)RB,(FTAPE(NA),NA=NAW1,NAW2,NAW3)                              
  731 FORMAT(5H  RB=,F5.1,(1X,11E11.3))                                         
  735 CONTINUE                                                                  
  740 CONTINUE                                                                  
  800 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      BLOCK DATA                                                                
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/GAUS/ABSCI1,ABSCI2,ABSCI3,ABSCI4,ABSCI5,ABSCI6,                    
     1            WEIGH1,WEIGH2,WEIGH3,WEIGH4,WEIGH5,WEIGH6                     
      DIMENSION ABSCI1(39),ABSCI2(39),ABSCI3(39),ABSCI4(39),ABSCI5(39),         
     1          ABSCI6(27)                                                      
      DIMENSION WEIGH1(39),WEIGH2(39),WEIGH3(39),WEIGH4(39),WEIGH5(39),         
     1          WEIGH6(27)                                                      
      DATA ABSCI1 /                                                             
     1   0.339981043584856  , 0.861136311594053  , 0.183434642495650  ,         
     2   0.525532409916329  , 0.796666477413627  , 0.960289856497536  ,         
     3   0.125233408511469  , 0.367831498998180  , 0.587317954286617  ,         
     4   0.769902674194305  , 0.904117256370475  , 0.981560634246719  ,         
     5   0.095012509837637  , 0.281603550779259  , 0.458016777657227  ,         
     6   0.617876244402644  , 0.755404408355003  , 0.865631202387832  ,         
     7   0.944575023073233  , 0.989400934991650  , 0.993128599185094  ,         
     8   0.963971927277913  , 0.912234428251325  , 0.839116971822218  ,         
     9   0.746331906460150  , 0.636053680726515  , 0.510867001950827  ,         
     A   0.373706088715419  , 0.227785851141645  , 0.076526521133497  ,         
     B   0.995187219997021  , 0.974728555971309  , 0.938274552002732  ,         
     C   0.886415527004401  , 0.820001985973902  , 0.740124191578554  ,         
     D   0.648093651936975  , 0.545421471388839  , 0.433793507626045   /        
      DATA ABSCI2 /                                                             
     1   0.315042679696163  , 0.191118867473616  , 0.064056892862605  ,         
     2   0.997263861849481  , 0.985611511545268  , 0.964762255587506  ,         
     3   0.934906075937739  , 0.896321155766052  , 0.849367613732569  ,         
     4   0.794483795967942  , 0.732182118740289  , 0.663044266930215  ,         
     5   0.587715757240762  , 0.506899908932229  , 0.421351276130635  ,         
     6   0.331868602282127  , 0.239287362252137  , 0.144471961582796  ,         
     7   0.048307665687738  , 0.998237709710559  , 0.990726238699457  ,         
     8   0.977259949983774  , 0.957916819213791  , 0.932812808278676  ,         
     9   0.902098806968874  , 0.865959503212259  , 0.824612230833311  ,         
     A   0.778305651426519  , 0.727318255189927  , 0.671956684614179  ,         
     B   0.612553889667980  , 0.549467125095128  , 0.483075801686178  ,         
     C   0.413779204371605  , 0.341994090825758  , 0.268152185007253  ,         
     D   0.192697580701371  , 0.116084070675255  , 0.038772417506050   /        
      DATA ABSCI3 /                                                             
     1   0.998771007252426  , 0.993530172266350  , 0.984124583722826  ,         
     2   0.970591592546247  , 0.952987703160430  , 0.931386690706554  ,         
     3   0.905879136715569  , 0.876572020274247  , 0.843588261624393  ,         
     4   0.807066204029442  , 0.767159032515740  , 0.724034130923814  ,         
     5   0.677872379632663  , 0.628867396776513  , 0.577224726083972  ,         
     6   0.523160974722233  , 0.466902904750958  , 0.408686481990716  ,         
     7   0.348755886292160  , 0.287362487355455  , 0.224763790394689  ,         
     8   0.161222356068891  , 0.097004699209462  , 0.032380170962869  ,         
     9   0.999305041735772  , 0.996340116771955  , 0.991013371476744  ,         
     A   0.983336253884625  , 0.973326827789910  , 0.961008799652053  ,         
     B   0.946411374858402  , 0.929569172131939  , 0.910522137078502  ,         
     C   0.889315445995114  , 0.865999398154092  , 0.840629296252580  ,         
     D   0.813265315122797  , 0.783972358943341  , 0.752819907260531   /        
      DATA ABSCI4 /                                                             
     1   0.719881850171610  , 0.685236313054233  , 0.648965471254657  ,         
     2   0.611155355172393  , 0.571895646202634  , 0.531279464019894  ,         
     3   0.489403145707052  , 0.446366017253464  , 0.402270157963991  ,         
     4   0.357220158337668  , 0.311322871990210  , 0.264687162208767  ,         
     5   0.217423643740007  , 0.169644420423992  , 0.121462819296120  ,         
     6   0.072993121787799  , 0.024350292663424  , 0.999553822651630  ,         
     7   0.997649864398237  , 0.994227540965688  , 0.989291302499755  ,         
     8   0.982848572738629  , 0.974909140585727  , 0.965485089043799  ,         
     9   0.954590766343634  , 0.942242761309872  , 0.928459877172445  ,         
     A   0.913263102571757  , 0.896675579438770  , 0.878722567678213  ,         
     B   0.859431406663111  , 0.838831473580255  , 0.816954138681463  ,         
     C   0.793832717504605  , 0.769502420135041  , 0.744000297583597  ,         
     D   0.717365185362099  , 0.689637644342027  , 0.660859898986119   /        
      DATA ABSCI5 /                                                             
     1   0.631075773046871  , 0.600330622829751  , 0.568671268122709  ,         
     2   0.536145920897131  , 0.502804111888784  , 0.468696615170544  ,         
     3   0.433875370831756  , 0.398393405881969  , 0.362304753499487  ,         
     4   0.325664370747701  , 0.288528054884511  , 0.250952358392272  ,         
     5   0.212994502857666  , 0.174712291832646  , 0.136164022809143  ,         
     6   0.097408398441584  , 0.058504437152420  , 0.019511383256793  ,         
     7   0.999689503883230  , 0.998364375863181  , 0.995981842987209  ,         
     8   0.992543900323762  , 0.988054126329623  , 0.982517263563014  ,         
     9   0.975939174585136  , 0.968326828463264  , 0.959688291448742  ,         
     A   0.950032717784437  , 0.939370339752755  , 0.927712456722308  ,         
     B   0.915071423120898  , 0.901460635315852  , 0.886894517402420  ,         
     C   0.871388505909296  , 0.854959033434601  , 0.837623511228187  ,         
     D   0.819400310737931  , 0.800308744139140  , 0.780369043867433   /        
      DATA ABSCI6 /                                                             
     1   0.759602341176647  , 0.738030643744400  , 0.715676812348967  ,         
     2   0.692564536642171  , 0.668718310043916  , 0.644163403784967  ,         
     3   0.618925840125468  , 0.593032364777572  , 0.566510418561397  ,         
     4   0.539388108324357  , 0.511694177154667  , 0.483457973920596  ,         
     5   0.454709422167743  , 0.425478988407300  , 0.395797649828908  ,         
     6   0.365696861472313  , 0.335208522892625  , 0.304364944354496  ,         
     7   0.273198812591049  , 0.241743156163840  , 0.210031310460567  ,         
     8   0.178096882367618  , 0.145973714654896  , 0.113695850110665  ,         
     9   0.081297495464425  , 0.048812985136049  , 0.016276744849602   /        
      DATA WEIGH1 /                                                             
     1   0.652145154862546  , 0.347854845137454  , 0.362683783378362  ,         
     2   0.313706645877887  , 0.222381034453374  , 0.101228536290376  ,         
     3   0.249147045813403  , 0.233492536538355  , 0.203167426723066  ,         
     4   0.160078328543346  , 0.106939325995318  , 0.047175336386512  ,         
     5   0.189450610455069  , 0.182603415044924  , 0.169156519395003  ,         
     6   0.149595988816577  , 0.124628971255534  , 0.095158511682493  ,         
     7   0.062253523938648  , 0.027152459411754  , 0.017614007139152  ,         
     8   0.040601429800386  , 0.062672048334109  , 0.083276741576704  ,         
     9   0.101930119817240  , 0.118194531961518  , 0.131688638449176  ,         
     A   0.142096109318382  , 0.149172986472603  , 0.152753387130725  ,         
     B   0.012341229799987  , 0.028531388628933  , 0.044277438817419  ,         
     C   0.059298584915436  , 0.073346481411080  , 0.086190161531953  ,         
     D   0.097618652104113  , 0.107444270115965  , 0.115505668053725   /        
      DATA WEIGH2 /                                                             
     1   0.121670472927803  , 0.125837456346828  , 0.127938195346752  ,         
     2   0.007018610009470  , 0.016274394730905  , 0.025392065309262  ,         
     3   0.034273862913021  , 0.042835898022226  , 0.050998059262376  ,         
     4   0.058684093478535  , 0.065822222776361  , 0.072345794108848  ,         
     5   0.078193895787070  , 0.083311924226946  , 0.087652093004403  ,         
     6   0.091173878695763  , 0.093844399080804  , 0.095638720079274  ,         
     7   0.096540088514727  , 0.004521277098533  , 0.010498284531152  ,         
     8   0.016421058381907  , 0.022245849194166  , 0.027937006980023  ,         
     9   0.033460195282547  , 0.038782167974472  , 0.043870908185673  ,         
     A   0.048695807635072  , 0.053227846983936  , 0.057439769099391  ,         
     B   0.061306242492928  , 0.064804013456601  , 0.067912045815233  ,         
     C   0.070611647391286  , 0.072886582395804  , 0.074723169057968  ,         
     D   0.076110361900626  , 0.077039818164247  , 0.077505947978424   /        
      DATA WEIGH3 /                                                             
     1   0.003153346052305  , 0.007327553901276  , 0.011477234579234  ,         
     2   0.015579315722943  , 0.019616160457355  , 0.023570760839324  ,         
     3   0.027426509708356  , 0.031167227832798  , 0.034777222564770  ,         
     4   0.038241351065830  , 0.041545082943464  , 0.044674560856694  ,         
     5   0.047616658492490  , 0.050359035553854  , 0.052890189485193  ,         
     6   0.055199503699984  , 0.057277292100403  , 0.059114839698395  ,         
     7   0.060704439165893  , 0.062039423159892  , 0.063114192286254  ,         
     8   0.063924238584648  , 0.064466164435950  , 0.064737696812683  ,         
     9   0.001783280721696  , 0.004147033260562  , 0.006504457968978  ,         
     A   0.008846759826363  , 0.011168139460131  , 0.013463047896718  ,         
     B   0.015726030476024  , 0.017951715775697  , 0.020134823153530  ,         
     C   0.022270173808383  , 0.024352702568710  , 0.026377469715054  ,         
     D   0.028339672614259  , 0.030234657072402  , 0.032057928354851   /        
      DATA WEIGH4 /                                                             
     1   0.033805161837141  , 0.035472213256882  , 0.037055128540240  ,         
     2   0.038550153178615  , 0.039953741132720  , 0.041262563242623  ,         
     3   0.042473515123653  , 0.043583724529323  , 0.044590558163756  ,         
     4   0.045491627927418  , 0.046284796581314  , 0.046968182816210  ,         
     5   0.047540165714830  , 0.047999388596458  , 0.048344762234802  ,         
     6   0.048575467441503  , 0.048690957009139  , 0.001144950003186  ,         
     7   0.002663533589512  , 0.004180313124694  , 0.005690922451403  ,         
     8   0.007192904768117  , 0.008683945269260  , 0.010161766041103  ,         
     9   0.011624114120797  , 0.013068761592401  , 0.014493508040509  ,         
     A   0.015896183583725  , 0.017274652056269  , 0.018626814208299  ,         
     B   0.019950610878141  , 0.021244026115782  , 0.022505090246332  ,         
     C   0.023731882865930  , 0.024922535764115  , 0.026075235767565  ,         
     D   0.027188227500486  , 0.028259816057276  , 0.029288369583267   /        
      DATA WEIGH5 /                                                             
     1   0.030272321759557  , 0.031210174188114  , 0.032100498673487  ,         
     2   0.032941939397645  , 0.033733214984611  , 0.034473120451753  ,         
     3   0.035160529044747  , 0.035794393953416  , 0.036373749905835  ,         
     4   0.036897714638276  , 0.037365490238730  , 0.037776364362001  ,         
     5   0.038129711314477  , 0.038424993006959  , 0.038661759774076  ,         
     6   0.038839651059051  , 0.038958395962769  , 0.039017813656306  ,         
     7   0.000796792065552  , 0.001853960788946  , 0.002910731817934  ,         
     8   0.003964554338444  , 0.005014202742927  , 0.006058545504235  ,         
     9   0.007096470791153  , 0.008126876925698  , 0.009148671230783  ,         
     A   0.010160770535008  , 0.011162102099838  , 0.012151604671088  ,         
     B   0.013128229566961  , 0.014090941772314  , 0.015038721026994  ,         
     C   0.015970562902562  , 0.016885479864245  , 0.017782502316045  ,         
     D   0.018660679627411  , 0.019519081140145  , 0.020356797154333   /        
      DATA WEIGH6 /                                                             
     1   0.021172939892191  , 0.021966644438744  , 0.022737069658329  ,         
     2   0.023483399085926  , 0.024204841792364  , 0.024900633222483  ,         
     3   0.025570036005349  , 0.026212340735672  , 0.026826866725591  ,         
     4   0.027412962726029  , 0.027970007616848  , 0.028497411065085  ,         
     5   0.028994614150555  , 0.029461089958167  , 0.029896344136328  ,         
     6   0.030299915420827  , 0.030671376123669  , 0.031010332586313  ,         
     7   0.031316425596861  , 0.031589330770727  , 0.031828758894411  ,         
     8   0.032034456231992  , 0.032206204794030  , 0.032343822568575  ,         
     9   0.032447163714064  , 0.032516118713868  , 0.032550614492363   /        
      END                                                                       
      SUBROUTINE NINEJ                                                          
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CRAC/FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC,L9(10),U9                   
      DIMENSION LT(9)                                                           
      U9=0.0D0                                                                    
      KEX=0                                                                     
      LT(1)=L9(1)                                                               
      LMIN=LT(1)                                                                
      KP=LMIN                                                                   
      IMIN=1                                                                    
      DO 20 I=2,9                                                               
      LT(I)=L9(I)                                                               
      KP=KP+L9(I)                                                               
      IF(LT(I)-LMIN) 15,15,20                                                   
   15 LMIN=LT(I)                                                                
      IMIN=I                                                                    
   20 CONTINUE                                                                  
      IF(LMIN) 1000,70,50                                                       
   50 GO TO (110,300,300,300,300,150,300,170,300),IMIN                          
   70 GO TO (110,110,110,110,150,150,170,170,190),IMIN                          
  110 MM=(IMIN-1)/2+1                                                           
      M1=MM+MM-1                                                                
      M2=M1+1                                                                   
      M3=MM+4                                                                   
      L1=LT(7)                                                                  
      LT(7)=LT(M1)                                                              
      LT(M1)=L1                                                                 
      L1=LT(8)                                                                  
      LT(8)=LT(M2)                                                              
      LT(M2)=L1                                                                 
      L1=LT(9)                                                                  
      LT(9)=LT(M3)                                                              
      LT(M3)=L1                                                                 
      IMIN=IMIN+(7-M1)                                                          
      GO TO 175                                                                 
  150 KEX=1                                                                     
      M1=7                                                                      
      M2=8                                                                      
      M3=IMIN+IMIN-9                                                            
      M4=M3+1                                                                   
      GO TO 180                                                                 
  170 KEX=1                                                                     
  175 M1=5                                                                      
      M2=6                                                                      
      M3=IMIN-6                                                                 
      M4=M3+2                                                                   
  180 L1=LT(M1)                                                                 
      LT(M1)=LT(M3)                                                             
      LT(M3)=L1                                                                 
      L1=LT(M2)                                                                 
      LT(M2)=LT(M4)                                                             
      LT(M4)=L1                                                                 
      L1=LT(9)                                                                  
      LT(9)=LT(IMIN)                                                            
      LT(IMIN)=L1                                                               
  190 IF(LT(9))1000,200,300                                                     
  200 IF(LT(5)-LT(6)) 1000,210,1000                                             
  210 IF(LT(7)-LT(8)) 1000,220,1000                                             
  220 IA=LT(1)                                                                  
      IB=LT(2)                                                                  
      IC=LT(3)                                                                  
      ID=LT(4)                                                                  
      IE=LT(5)                                                                  
      IFF=LT(7)                                                                  
      R1=(IE+1)*(IFF+1)                                                          
      R1= DSQRT(R1)                                                              
      S1=(-1.0D0  )**((IE+IFF-IA-ID)/2)                                            
      CALL RAC7                                                                 
      U9=RAC*S1/R1                                                              
      GO TO 370                                                                 
  300 K1=IABS (LT(2)-LT(7))                                                     
      K2=IABS (LT(3)-LT(5))                                                     
      K3=IABS (LT(4)-LT(9))                                                     
      MINRDA=MAX0  (K1,K2,K3)                                                   
      K1=LT(2)+LT(7)                                                            
      K2=LT(3)+LT(5)                                                            
      K3=LT(4)+LT(9)                                                            
      MAXRDA=MIN0  (K1,K2,K3)                                                   
      IF(MINRDA-MAXRDA) 320,320,1000                                            
  320 DO 350 N1=MINRDA,MAXRDA,2                                                 
      RAMDA2=N1                                                                 
      IA=LT(2)                                                                  
      IB=LT(5)                                                                  
      IC=LT(7)                                                                  
      ID=LT(3)                                                                  
      IE=LT(1)                                                                  
      IFF=N1                                                                     
      CALL RAC7                                                                 
      W1=(RAMDA2+1.0D0  )*RAC                                                     
      IB=LT(4)                                                                  
      ID=LT(9)                                                                  
      IE=LT(8)                                                                  
      CALL RAC7                                                                 
      W1=W1*RAC                                                                 
      IA=LT(3)                                                                  
      IC=LT(5)                                                                  
      IE=LT(6)                                                                  
      CALL RAC7                                                                 
      W1=W1*RAC                                                                 
  350 U9=U9+W1                                                                  
      IF( DABS(U9).LT.1.0D-10) U9=0.0                                            
  370 IF(KEX) 400,1000,400                                                      
  400 U9=U9*((-1.0D0  )**(KP/2))                                                  
 1000 RETURN                                                                    
      END                                                                       
      SUBROUTINE RAC7                                                           
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON / CRAC / FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC                         
      DIMENSION LT(6)                                                           
      K1=IA+IB-IE                                                               
      K3=IC+ID-IE                                                               
      K5=IA+IC-IFF                                                               
      K7=IB+ID-IFF                                                               
      K2=IE-IABS(IA-IB)                                                         
      K4=IE-IABS(IC-ID)                                                         
      K6=IFF-IABS(IA-IC)                                                         
      K8=IFF-IABS(IB-ID)                                                         
      K9=MIN0(K1,K2,K3,K4,K5,K6,K7,K8)                                          
      RAC=0.0D0                                                                   
      IF (K9) 4000,20,20                                                        
   20 K2=K1-2*(K1/2)                                                            
      K4=K3-2*(K3/2)                                                            
      K6=K5-2*(K5/2)                                                            
      K8=K7-2*(K7/2)                                                            
      IF(MAX0(K2,K4,K6,K8)) 4000,25,4000                                        
   25 LTMIN=MIN0(IA,IB,IC,ID,IE,IFF)                                             
      IF(LTMIN) 4000,30,150                                                     
   30 LT(1)=IA                                                                  
      LT(2)=IB                                                                  
      LT(3)=IC                                                                  
      LT(4)=ID                                                                  
      LT(5)=IE                                                                  
      LT(6)=IFF                                                                  
      LTMIN=LT(1)                                                               
      KMIN=1                                                                    
      DO 40 N=2,6                                                               
      IF(LT(N)-LTMIN)35,40,40                                                   
   35 LTMIN=LT(N)                                                               
      KMIN=N                                                                    
   40 CONTINUE                                                                  
      S1=1.0D0                                                                    
      F1=IE                                                                     
      F2=IFF                                                                     
      GO TO (55,55,55,55,45,50),KMIN                                            
   45 F1=IA                                                                     
      F2=IC                                                                     
      S1=1-2*MOD(K5/2,2)                                                        
      GO TO 55                                                                  
   50 F1=IA                                                                     
      F2=IB                                                                     
      S1=1-2*MOD(K1/2,2)                                                        
   55 RAC=S1/ DSQRT((F1+1.  )*(F2+1.  ))                                         
      GO TO 4000                                                                
  150 IABEP=(IA+IB+IE)/2+1                                                      
      ICDEP=(IC+ID+IE)/2+1                                                      
      IACFP=(IA+IC+IFF)/2+1                                                      
      IBDFP=(IB+ID+IFF)/2+1                                                      
      IABE=IABEP-IE                                                             
      IEAB=IABEP-IB                                                             
      IBEA=IABEP-IA                                                             
      ICDE=ICDEP-IE                                                             
      IECD=ICDEP-ID                                                             
      IDEC=ICDEP-IC                                                             
      IACF=IACFP-IFF                                                             
      IFAC=IACFP-IC                                                             
      ICFA=IACFP-IA                                                             
      IBDF=IBDFP-IFF                                                             
      IFBD=IBDFP-ID                                                             
      IDFB=IBDFP-IB                                                             
      NZMAX=MIN0(IABE,ICDE,IACF,IBDF)                                           
      IABCD1=(IA+IB+IC+ID+4)/2                                                  
      IEFMAD=(IE+IFF-IA-ID)/2                                                    
      IEFMBC=(IE+IFF-IB-IC)/2                                                    
      NZMI1=-IEFMAD                                                             
      NZMI2=-IEFMBC                                                             
      NZMIN=MAX0(0,NZMI1,NZMI2)+1                                               
      SQLOG=0.5D0*(FACLOG(IABE)+FACLOG(IEAB)+FACLOG(IBEA)+FACLOG(ICDE)          
     1          +FACLOG(IECD)+FACLOG(IDEC)+FACLOG(IACF)+FACLOG(IFAC)            
     2          +FACLOG(ICFA)+FACLOG(IBDF)+FACLOG(IFBD)+FACLOG(IDFB)            
     3 -FACLOG(IABEP+1)-FACLOG(ICDEP+1)-FACLOG(IACFP+1)-FACLOG(IBDFP+1))        
      DO 200 NZ=NZMIN,NZMAX                                                     
      NZM1=NZ-1                                                                 
      S1=1-2*MOD(NZM1,2)                                                        
      K1=IABCD1-NZM1                                                            
      K2=IABE-NZM1                                                              
      K3=ICDE-NZM1                                                              
      K4=IACF-NZM1                                                              
      K5=IBDF-NZM1                                                              
      K6=NZ                                                                     
      K7=IEFMAD+NZ                                                              
      K8=IEFMBC+NZ                                                              
      SSLOG=SQLOG+FACLOG(K1)-FACLOG(K2)-FACLOG(K3)-FACLOG(K4)                   
     1           -FACLOG(K5)-FACLOG(K6)-FACLOG(K7)-FACLOG(K8)                   
      SSTERM=S1* DEXP(SSLOG)                                                     
      RAC=RAC+SSTERM                                                            
  200 CONTINUE                                                                  
      IF( DABS(RAC).LT.1.0D-10 ) RAC=0.0D0                                         
 4000 RETURN                                                                    
      END                                                                       
      SUBROUTINE CLEB                                                           
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON / CRAC / FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC                         
      RAC=0.0D0                                                                   
      IF(ID+IE-IFF) 1000,105,1000                                                
  105 K1=IA+IB+IC                                                               
      IF(K1-2*(K1/2)) 1000,110,1000                                             
  110 K1=IA+IB-IC                                                               
      K2=IC-IABS(IB-IA)                                                         
      K3=MIN0(K1,K2)                                                            
      IF(K3) 1000,130,130                                                       
  130 IF((-1)**(IB+IE)) 1000,1000,140                                           
  140 IF((-1)**(IC+IFF)) 1000,1000,150                                           
  150 IF(IA-IABS(ID)) 1000,152,152                                              
  152 IF(IB-IABS(IE)) 1000,154,154                                              
  154 IF(IC-IABS(IFF)) 1000,160,160                                              
  160 IF(IA) 1000,175,165                                                       
  165 IF(IB) 1000,175,170                                                       
  170 IF(IC) 1000,180,250                                                       
  175 RAC=1.0D0                                                                   
      GO TO 1000                                                                
  180 FB=IB+1                                                                   
      RAC=((-1.0D0  )**((IA-ID)/2))/ DSQRT(FB)                                     
      GO TO 1000                                                                
  250 FC2=IC+1                                                                  
      IABCP=(IA+IB+IC)/2+1                                                      
      IABC=IABCP-IC                                                             
      ICAB=IABCP-IB                                                             
      IBCA=IABCP-IA                                                             
      IAPD=(IA+ID)/2+1                                                          
      IAMD=IAPD-ID                                                              
      IBPE=(IB+IE)/2+1                                                          
      IBME=IBPE-IE                                                              
      ICPF=(IC+IFF)/2+1                                                          
      ICMF=ICPF-IFF                                                              
      SQFCLG=0.5D0  *(DLOG(FC2)-FACLOG(IABCP+1)                                   
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)                             
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)                             
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))                            
      NZMIC2=(IB-IC-ID)/2                                                       
      NZMIC3=(IA-IC+IE)/2                                                       
      NZMI=MAX0(0,NZMIC2,NZMIC3)+1                                              
      NZMX=MIN0(IABC,IAMD,IBPE)                                                 
      S1=1-2*MOD(NZMI-1,2)                                                      
      DO 400 NZ=NZMI,NZMX                                                       
      NZM1=NZ-1                                                                 
      NZT1=IABC-NZM1                                                            
      NZT2=IAMD-NZM1                                                            
      NZT3=IBPE-NZM1                                                            
      NZT4=NZ-NZMIC2                                                            
      NZT5=NZ-NZMIC3                                                            
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)                        
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)                        
      SSTERM=S1* DEXP(TERMLG)                                                    
      RAC=RAC+SSTERM                                                            
  400 S1=-S1                                                                    
      IF( DABS(RAC).LT.1.0D-10) RAC=0.0D0                                          
 1000 RETURN                                                                    
      END                                                                       
      SUBROUTINE CLEBZ                                                          
	IMPLICIT REAL*8(A-H,O-Z)
      COMMON / CRAC / FACLOG(500),IA,IB,IC,ID,IE,IFF,RAC                         
      RAC=0.0D0                                                                   
      IGTW=(IA+IB+IC)/2                                                         
      IF(MOD(IGTW,2).NE.0) GO TO 1000                                           
      IG=IGTW/2                                                                 
      IAHF=IA/2                                                                 
      IBHF=IB/2                                                                 
      ICHF=IC/2                                                                 
      S1=1-2*MOD(IG+ICHF,2)                                                     
      IABC=IGTW+2                                                               
      IABMC=IAHF+IBHF-ICHF+1                                                    
      ICAMB=ICHF+IAHF-IBHF+1                                                    
      IBCMA=IBHF+ICHF-IAHF+1                                                    
      IGMA=IG-IAHF+1                                                            
      IGMB=IG-IBHF+1                                                            
      IGMC=IG-ICHF+1                                                            
      R1=0.5D0*(FACLOG(IABMC)+FACLOG(ICAMB)+FACLOG(IBCMA)-FACLOG(IABC))         
     1         +FACLOG(IG+1)-FACLOG(IGMA)-FACLOG(IGMB)-FACLOG(IGMC)             
      RAC=S1* DSQRT( DFLOAT(IC+1))* DEXP(R1)                                       
 1000 RETURN                                                                    
      END                                                                       
