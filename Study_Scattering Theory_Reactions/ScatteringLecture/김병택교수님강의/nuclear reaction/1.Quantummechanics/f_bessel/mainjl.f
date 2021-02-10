      PROGRAM MAIN 
	IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AK(500),AI(500)

      open(unit=5, file='jlcal.DAT',status='old')
      open(unit=6, file='jlcal.OUT',status='unknown')

      PI  =4.0*ATAN(1.D0) 
      HBAR=197.327053
      AMAS=931.49432
      WNUNIT=SQRT(2.*AMAS)/HBAR
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
15    FORMAT(1H ,'SPHERICAL BESSEL FUNCTIONS FOR'/
     15X,'TMAS=',F7.3,' ZZT=',F7.3,' PMAS=',F7.3,' ZZP=',F7.3,' RMAS='
     2F7.3)
16    FORMAT(5X,'ELAB=',F7.3,' ECM=',F7.3,' WN  =',F7.3,' ETA=',F7.3,
     2 ' TURN=',F7.3)                                            

      LL=MAXL-1
      WRITE(6,20) LL
20    FORMAT(1H ,'REGULAR SPHERICAL BESSEL FUNCTIONS UP TO L=', I3, /
     1 1H ,'   X         L=0            L=1           L=2',/) 

      DO 100 N=1,NXMAX
	XX =XMIN+ DFLOAT(N-1)*XSTEP
	X=XX*WN
	CALL BESSEL(X,MAXL,AK,AI)
	WRITE(6,21) XX, (AK(L),L=1,MAXL)
100   CONTINUE
21    FORMAT(F7.3,10E14.6)
      STOP
      END                                                                        
CCCCCC
CCCCCC
C==========================================================================
      SUBROUTINE BESSEL(X,LP1MX,AK,AI)                                
C==========================================================================
	IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AK(500),AI(500)
      DIMENSION F(51),G(51)                                             CCN45060
                                                                        CCN45070
C     THE FOLLOWING IS FOR SPHERICAL BESSEL AND NEUMANN                 CCN45080
      AK(1)=DSIN(X)/X                                                   CCN45090
      AK(2)=(AK(1)-DCOS(X))/X                                           CCN45100
      AI(1)=-DCOS(X)/X                                                  CCN45110
      AI(2)=(AI(1)-DSIN(X))/X                                           CCN45120
      DO 100 I=3,LP1MX                                                  CCN45130
      L=I-1                                                             CCN45140
      T1=(L*2-1)                                                        CCN45150
      AK(I)=AK(I-1)*T1/X-AK(I-2)                                        CCN45160
      AI(I)=AI(I-1)*T1/X-AI(I-2)                                        CCN45170
  100 CONTINUE                                                          CCN45180
 7810 FORMAT(1X,F9.4)                                                   CCN45210
 7800 FORMAT(1X,6E12.4)                                                 CCN45220
C     THE FOLLOWING IS ANOTHER VERSION CALCULATING SPHERICAL BESSEL     CCN45240
      IF (X.EQ.0.0) THEN
	AK(1)=1.0
	DO 10 IL=2,LP1MX
	AK(IL)=0.0
10    CONTINUE
      RETURN
	ENDIF
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
      IF(DABS (SS/SSS)-1.0E-10)630,630,540                              CCN45560
  540 T5=-T9*TL-SS/X                                                    CCN45570
      TL=T9*SL-TS/X                                                     CCN45580
      SL=T5                                                             CCN45590
      SSS=SSS+SS                                                        CCN45600
      STS=STS+TS                                                        CCN45610
      SSL=SSL+SL                                                        CCN45620
      STL=STL+TL                                                        CCN45630
      EN=EN+1.                                                          CCN45640
  620 CONTINUE                                                          CCN45650
  630 T8=DSIN (T3)                                                      CCN45660
      T9=DCOS (T3)                                                      CCN45670
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
 1080 IF(1.-DABS (Y))1090,1090,1120                                     CCN45860
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
      IF(L.EQ.0) AK(1)=DSIN(R)/R                                        CCN45980
      GO TO 2000
 1900 IF(R.LT.5) THEN
        T1=IL-1
        IF(IL.GT.1) THEN
          T5=1.0
          DO 1910 N1=1,IL
          T4=DFLOAT(IL+N1-1)
          T5=T5/T4 
 1910     CONTINUE         
          T3=DLOG(T5)
        ELSE
          T3=0.0
        ENDIF    
        T2=T1*DLOG(R*2.0)+T3
        AK(IL)=DEXP(T2)
      ELSE
        T1=(IL-1)*2-1
        AK(IL)=AK(IL-1)*T1/R-AK(IL-2)
      ENDIF
 2000 CONTINUE                                                          CCN45990
      RETURN                                                            CCN46010
      END                                                               CCN46020
