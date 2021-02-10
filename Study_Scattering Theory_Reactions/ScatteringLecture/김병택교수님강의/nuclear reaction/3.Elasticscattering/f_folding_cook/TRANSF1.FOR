      SUBROUTINE TRANSF1 (L,ITYPE,AIN,AOUT,DR,IMAX,DQ,KMAX)             ABQP0708
C                                                                       ABQP0709
C     L IS THE MULTIPOLARITY                                            ABQP0710
C                                                                       ABQP0711
C     ITYPE = 0    CALCULATES FOURIER TRANSFORM                         ABQP0712
C     ITYPE = 1    CALCULATES INVERSE TRANSFORM                         ABQP0713
C                                                                       ABQP0714
C     BASED ON FILON'S FORMULA IN M ABRAMOWITZ AND I A STEGUN (1965)    ABQP0715
C     HANDBOOK OF MATHEMATICAL FUNCTIONS, DOVER, NEW YORK.              ABQP0716
C                                                                       ABQP0717
      IMPLICIT REAL*8 (A-H,O-Z)                                         ABQP0718
CX      COMMON/ALL/PI,DR,DQ,IMAX,KMAX,KPRINT                              ABQP0719
      DIMENSION    AIN(1000),AOUT(1000)                                 ABQP0720
C                                                                       ABQP0721
C     AIN   CONTAINS THE FUNCTION                                       ABQP0722
C     AOUT  CONTAINS THE FOURIER TRANSFORM                              ABQP0723
C                                                                       ABQP0724
      REAL*8 K,KH,KH2,KH3                                               ABQP0725
      PI=3.14159265358979D0                                             
      PIPI=2.0D0*PI                                                     ABQP0726
      IF(ITYPE.EQ.1) GOTO 5                                             ABQP0727
C                                                                       ABQP0728
C     SETS INTERVALS AND LIMITS FOR FOURIER TRANSFORM                   ABQP0729
C                                                                       ABQP0730
      HA=DR                                                             ABQP0731
      HB=DQ                                                             ABQP0732
      NA=IMAX                                                           ABQP0733
      NB=KMAX                                                           ABQP0734
      CONST=4.0D0*PI                                                    ABQP0735
      GOTO 15                                                           ABQP0736
C                                                                       ABQP0737
C     SETS INTERVALS AND LIMITS FOR INVERSE TRANSFORM                   ABQP0738
C                                                                       ABQP0739
    5 HA=DQ                                                             ABQP0740
      HB=DR                                                             ABQP0741
      NA=KMAX                                                           ABQP0742
      NB=IMAX                                                           ABQP0743
      CONST=1.0D0/(2.0D0*PI**2)                                         ABQP0744
   15 IF(L.GT.0) GOTO 110                                               ABQP0745
C                                                                       ABQP0746
C     INTEGRATION FOR L.EQ.0                                            ABQP0747
C                                                                       ABQP0748
      DO 100 J=1,NB                                                     ABQP0749
      K=HB*DFLOAT(J)                                                    ABQP0750
      KH=K*HA                                                           ABQP0751
      KH2=KH*KH                                                         ABQP0752
      KH3=KH2*KH                                                        ABQP0753
      SINKH=DSIN(DMOD(KH,PIPI))                                         ABQP0754
      COSKH=DCOS(DMOD(KH,PIPI))                                         ABQP0755
      A=1.0D0/KH+SINKH*COSKH/KH2-2.0D0*SINKH*SINKH/KH3                  ABQP0756
      B=(1.0D0+COSKH*COSKH)/KH2-2.0D0*SINKH*COSKH/KH3                   ABQP0757
      C=4.0D0*(SINKH/KH3-COSKH/KH2)                                     ABQP0758
      EVEN=0.0D0                                                        ABQP0759
      ODD=0.0                                                           ABQP0760
      DO 20 I=2,NA,2                                                    ABQP0761
      F=DFLOAT(I)                                                       ABQP0762
      X=F*KH                                                            ABQP0763
      R=F*HA                                                            ABQP0764
   20 EVEN=EVEN+AIN(I)*R*DSIN(DMOD(X,PIPI))                             ABQP0765
      DO 25 I=1,NA,2                                                    ABQP0766
      F=DFLOAT(I)                                                       ABQP0767
      X=F*KH                                                            ABQP0768
      R=F*HA                                                            ABQP0769
   25 ODD=ODD+AIN(I)*R*DSIN(DMOD(X,PIPI))                               ABQP0770
      R=HA*DFLOAT(NA)                                                   ABQP0771
      X=R*K                                                             ABQP0772
      AA=AIN(NA)*R                                                      ABQP0773
      S=-AA*DSIN(DMOD(X,PIPI))+2.0D0*EVEN                               ABQP0774
  100 AOUT(J)=CONST*HA*(-A*AA*DCOS(DMOD(X,PIPI))+B*S+C*ODD)/K           ABQP0775
      RETURN                                                            ABQP0776
C                                                                       ABQP0777
C     INTEGRATION FOR L.GT.0                                            ABQP0778
C                                                                       ABQP0779
  110 DO 200 J=1,NB                                                     ABQP0780
      K=HB*DFLOAT(J)                                                    ABQP0781
      KH=K*HA                                                           ABQP0782
      KH2=KH*KH                                                         ABQP0783
      KH3=KH2*KH                                                        ABQP0784
      SINKH=DSIN(DMOD(KH,PIPI))                                         ABQP0785
      COSKH=DCOS(DMOD(KH,PIPI))                                         ABQP0786
      A=1.0D0/KH+SINKH*COSKH/KH2-2.0D0*SINKH*SINKH/KH3                  ABQP0787
      B=(1.0D0+COSKH*COSKH)/KH2-2.0D0*SINKH*COSKH/KH3                   ABQP0788
      C=4.0D0*(SINKH/KH3-COSKH/KH2)                                     ABQP0789
      EVEN=0.0D0                                                        ABQP0790
      ODD=0.0D0                                                         ABQP0791
      DO 240 I=2,NA,2                                                   ABQP0792
      F=DFLOAT(I)                                                       ABQP0793
      X=F*KH                                                            ABQP0794
      R=F*HA                                                            ABQP0795
      AA=AIN(I)*R*R                                                     ABQP0796
      U=AA*SBES(L,X)                                                    ABQP0797
      V=AA*CBES(L,X)                                                    ABQP0798
  240 EVEN=EVEN+U*DSIN(DMOD(X,PIPI))+V*DCOS(DMOD(X,PIPI))               ABQP0799
      DO 250 I=1,NA,2                                                   ABQP0800
      F=DFLOAT(I)                                                       ABQP0801
      X=F*KH                                                            ABQP0802
      R=F*HA                                                            ABQP0803
      AA=AIN(I)*R*R                                                     ABQP0804
      U=AA*SBES(L,X)                                                    ABQP0805
      V=AA*CBES(L,X)                                                    ABQP0806
  250 ODD=ODD+U*DSIN(DMOD(X,PIPI))+V*DCOS(DMOD(X,PIPI))                 ABQP0807
      R=HA*DFLOAT(NA)                                                   ABQP0808
      X=R*K                                                             ABQP0809
      AA=AIN(NA)*R*R                                                    ABQP0810
      UA=AA*SBES(L,X)                                                   ABQP0811
      VA=AA*CBES(L,X)                                                   ABQP0812
      S=-UA*DSIN(DMOD(X,PIPI))-VA*DCOS(DMOD(X,PIPI))+2.0D0*EVEN         ABQP0813
  200 AOUT(J)=CONST*HA*(-A*UA*DCOS(DMOD(X,PIPI))+A*VA*DSIN(DMOD(X,PIPI))ABQP0814
     *       +B*S+C*ODD)                                                ABQP0815
      RETURN                                                            ABQP0816
      END                                                               ABQP0817
