$LARGE
C==============================================================C
C   CALCULATION OF THE RGM-KERNEL FOR THE D - ALPHA - SYSTEM   C
C==============================================================C
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 COU,DN,J,J2,L,PO,PZN,PSN,STANZ,DIS,SDUMMY
C
      COMMON /BREITE/ A,B(3,3),C(3),DIS,DN
      COMMON /POTPAR/ H2M,POT(8,3),PZN,PSN,COU
      COMMON /DISK/ STEP,STANZ
C
      H2M   = 20.735D0
      PI    = 3.141592653589793D0
      ELEML = DSQRT(1.4409D0)
C
      OPEN(5,FILE='ALPHAD.INP')
C -----------------------------------C
C READ IN THE VARIOUS PARAMETERS     C
C -----------------------------------C
      READ(5,100)
      READ(5,*)    PO
      READ(5,100)
      READ(5,*)    DIS

      READ(5,*)    J2,L,SDUMMY
                    J = J2/2

      READ(5,*)
      READ(5,*)    U
      READ(5,*)    POT(8,1),POT(7,1)
              pot(8,1) = - pot(8,1)
             PSN     =  1
      READ(5,*)    STEP,STANZ
         COU = 1
C
        CLOSE(5)
C
100    FORMAT(1X,//)

C
      IF (PO.EQ.1) THEN
C         CENTRAL POTENTIAL
             PZN= 2
             POT(1,1)=     U/4.D0
             POT(2,1)=     U/4.D0
             POT(3,1)= (2.D0-U)/4.D0
             POT(4,1)= (2.D0-U)/4.D0
             POT(5,1)=    0.415D0
             POT(6,1)=   -66.92D0
             POT(1,2)=     U/4.D0
             POT(2,2)=    -U/4.D0
             POT(3,2)=-(2.D0-U)/4.D0
             POT(4,2)= (2.D0-U)/4.D0
             POT(5,2)=    0.292D0
             POT(6,2)=   -29.05D0
             PSN     =  1

C  WIDTH OF ALPHA-PARTICLE WAVE FUNCTION
             A = 0.514d0
C
C  CONSTANTS IN THE DEUTERON WAVE FUNCTION
               DN  = 3
               C(1)= 0.07284D0
               C(2)= 0.36570D0
               C(3)= 1.46960D0
               B(1,1)=  0.438511973721659740E-01
               B(1,2)=  0.159209549369291978
               B(1,3)=  0.251955354241984342
               B(2,1)=  0.121414399053524857
               B(2,2)= -0.319359254404235446
               B(2,3)= -0.119867858528174068
               B(3,1)=  0.504353095551319033E-01
               B(3,2)= -0.529019766071656655
               B(3,3)=  1.41741354474524539
      ELSE
         IF(PO.EQ.2) THEN
              PZN= 2
                YY  =  0.94D0
                MM  =  U
              POT(6,1) = -72.98D0
              POT(6,2) = -72.98D0
              POT(5,1) =  0.46D0
              POT(5,2) =  0.46D0
C
              POT(1,1) =             YY*(1.D0+MM)/4.D0
              POT(4,1) =             YY*(1.D0+MM)/4.D0
              POT(2,1) =             YY*(1.D0-MM)/4.D0
              POT(3,1) =             YY*(1.D0-MM)/4.D0
C
              POT(1,2) = (1.D0-YY)*(1.D0-3.D0*MM)/6.D0
              POT(4,2) = (1.D0-YY)*(1.D0+3.D0*MM)/3.D0
              POT(2,2) = (1.D0-YY)*(1.D0+3.D0*MM)/6.D0
              POT(3,2) = (1.D0-YY)*(1.D0-3.D0*MM)/3.D0
C
C  WIDTH OF ALPHA-PARTICLE WAVE FUNCTION
              A       = 0.540D0
C
C  CONSTANTS IN THE DEUTERON WAVE FUNCTION
               DN  = 3

                 C(1) = 0.320D0
                 C(2) = 0.440D0
                 C(3) = 0.600D0
C
              B(1,1) = 1.1130824469D0
              B(1,2) =-2.2102402241D0
              B(1,3) = 1.5368531890D0
C
              B(2,1) =-2.4016415961D0
              B(2,2) = 4.2766276372D0
              B(2,3) =-1.5562673117D0
C
              B(3,1) = 3.8234300316D0
              B(3,2) =-10.5292546988D0
              B(3,3) = 7.3266860773D0
C
                  ELSE
              PZN= 3
              MM = U
              POT(6,1) =   -6.000D0
C             POT(6,2) = -546.000D0
              POT(6,2) = -474.200d0
CC see Kanada,Kaneko,Saito,Tang, Nucl. Phys. A444 (1985) 209.
              POT(6,3) = 1655.000D0
              POT(5,1) =    0.160D0
              POT(5,2) =    1.127D0
              POT(5,3) =    3.400D0
C
              POT(1,1) =-0.2361D0
              POT(4,1) = 1.1528D0
              POT(2,1) = 0.5972D0
              POT(3,1) =-0.5139D0
C
              POT(1,2) = 0.8297D0-MM
              POT(4,2) =      MM
              POT(2,2) = 0.1401D0
              POT(3,2) = 0.0302D0
C
              POT(1,3) = 0.4474D0
              POT(4,3) = 0.3985D0
              POT(2,3) = 0.1015D0
              POT(3,3) = 0.0526D0
C
C  WIDTH OF ALPHA-PARTICLE WAVE FUNCTION
              A       =     0.580D0

C  CONSTANTS IN THE DEUTERON WAVE FUNCTION
               DN  = 3
C
                 C(1) = 0.120D0
                 C(2) = 1.200D0
                 C(3) = 1.400D0
C
              B(1,1) = 0.0924494391D0
              B(1,2) = 1.7517299019D0
              B(1,3) =-1.5374092840D0
C
              B(2,1) = 0.1468048247D0
              B(2,2) =-3.9229456793D0
              B(2,3) = 3.4827875051D0
C
              B(3,1) = 0.0502491686D0
              B(3,2)= -8.9267477273D0
              B(3,3) =10.2133481498D0
         ENDIF
      ENDIF

C
C     **************
      CALL KERN(J,L)
C     **************
C
C
      STOP
      END
C
C     ********************
      SUBROUTINE KERN(J,L)
C     ********************
c
      IMPLICIT REAL*8 (A-Z)
c For PC under DOS: Maximum number of meshpoints is 50.
      PARAMETER (MAXPNT=50)
      INTEGER*4  COU,I,J,K,L,M,PZN,PSN,Q,R,STANZ,BRA,KET,DIS,FIN,
     1           START,DN,KK
      REAL*8 KMAT(MAXPNT,MAXPNT,3,3),KNORM(MAXPNT,MAXPNT,3,3),
     1       KVEC(MAXPNT+1,3,3),KENE(3,3),KENED(3,3),NA0(3,3)
      CHARACTER NUMBER(9) /'1','2','3','4','5','6','7','8','9'/
      CHARACTER*13 OUTNAME(6)
C
      COMMON /BREITE/ A,B(3,3),C(3),DIS,DN
      COMMON /POTPAR/ H2M,POT(8,3),PZN,PSN,COU
      COMMON /DISK/ STEP,STANZ
C
      common /tkcom/tkb1,tkb2,tkb3,tkb4,tkb5,tkb6,tkb7,
     &              tkc1,tkc2,tkc3,tkc4,tkc5,tkc6,tkc7
      common /vzcom/vzb0,vzb1,vzb2,vzb3,
     &              vzc0,vzc1,vzc2,vzc3,
     &              vzd0,vzd1,vzd2,vzd3,
     &              vze0,vze1,vze2,vze3,
     &              vzg0,vzg1,vzg2,vzg3,
     &              vzh0,vzh1,vzh2,vzh3,
     &              vzi0,vzi1,vzi2,vzi3,
     &              vzj0,vzj1,vzj2,vzj3,
     &              vzk0,vzk1,vzk2,vzk3,
     &              vzl0,vzl1,vzl2,vzl3,
     &              vzm0,vzm1,vzm2,vzm3,
     &              vzo0,vzo1,vzo2,vzo3,
     &              vzp0,vzp1,vzp2,vzp3
      common /vscom/vsb0,vsb1,vsb2,vsb3,
     &              vsc0,vsc1,vsc2,vsc3,
     &              vsd0,vsd1,vsd2,vsd3
C
      PI    = 3.141592653589793D0
      E2    = 1.4409D0
C
      IF (STANZ.GT.MAXPNT) THEN
        WRITE(6,*)'NUMBER OF MESHPOINTS TOO BIG! PROGRAM TERMINATES!'
      ENDIF
C
C     LS-COUPLING
C
      IF(J.EQ.L+1) THEN
            LKOFF = DFLOAT(L)
      ELSE
         IF(J.EQ.L) THEN
            LKOFF =-1.D0
         ELSE
            LKOFF =-(DFLOAT(L)+1.D0)
         ENDIF
      ENDIF
      LKOFF=2.D0*LKOFF


C
C     COULOMB-PREFAKTORS
C
      FVC1=1.D0
      FVC2=1.D0
      FVC3=0.D0
      FVC4=0.D0
C
C  *****  SKALAR PART  *****
C
C     BINDING ENERGY OF THE ALPHA-PARTICLE
C
      TKA2A=H2M*9.D0*A/2.D0
      VCA0A=E2*DSQRT(2.D0*A/PI)
      VZA0 =0.D0
      DO 100 M=1,PZN,1
           PW=POT(1,M)
           PM=POT(4,M)
           LAM=POT(5,M)
           V0=POT(6,M)
           VZA0=VZA0+V0*6.D0*(PW+PM)*DSQRT((A/(A+2.D0*LAM))**3)
100   CONTINUE
      ENEA=TKA2A+VCA0A+VZA0


C
C     BINDING ENERGY OF THE DEUTERON
C
      DO 201 BRA=1,DIS,1
      DO 201 KET=1,DIS,1
         NA0D=0.D0
         TKA2D=0.D0
         VZN0=0.D0
         VCMD=0.D0
         DO 1000 I=1,DN,1
            AI=C(I)
            BI=B(BRA,I)
         DO 1000 K=1,DN,1
            AJ=C(K)
            BJ=B(KET,K)
C
            NA0D=NA0D+BI*BJ/DSQRT((AI+AJ)**3)
            TKA2D=TKA2D+3.D0*BI*BJ*AI*AJ/DSQRT((AI+AJ)**5)
            VCMD=VCMD+FVC3*BI*BJ/(DSQRT(PI)*(AI+AJ))
C
            DO 110 M=1,PZN,1
               PW=POT(1,M)
               PB=POT(2,M)
               PH=POT(3,M)
               PM=POT(4,M)
               LAM=POT(5,M)
               V0=POT(6,M)
               VZN0=VZN0+V0*BI*BJ*(PW+PB+PH+PM)/
     1              DSQRT((AI+AJ+4.D0*LAM)**3)
110         CONTINUE
1000     CONTINUE
         TKA2D=H2M*TKA2D/NA0D
         VZN0=VZN0/NA0D
         VCMD=E2*VCMD/NA0D
         KENED(BRA,KET)=TKA2D+VZN0+VCMD
201   CONTINUE
C


C     NORMALIZATION CONSTANTS
C
      DO 202 BRA=1,DIS,1
      DO 202 KET=1,DIS,1
         NA0(BRA,KET)=0.D0
202   CONTINUE
C
      DO 1010 I=1,DN,1
         AI=C(I)
      DO 1010 K=1,DN,1
         AJ=C(K)
         H=(64.*SQRT((AI+AJ)*A*PI)*PI**6)/(SQRT(PI)*A**5*(AI**2+2.*
     1     AI*AJ+AJ**2))
         DO 203 BRA=1,DIS,1
            BI=B(BRA,I)
         DO 203 KET=1,DIS,1
            BJ=B(KET,K)
            NA0(BRA,KET)=NA0(BRA,KET)+BI*BJ*H
203      CONTINUE
1010  CONTINUE


C
C  *****  VECTOR PART  *****
C
      DO 204 BRA=1,DIS,1
      DO 204 KET=1,DIS,1
         DO 10 Q=1,STANZ+1,1
            KVEC(Q,BRA,KET)=0.D0
10       CONTINUE
204   CONTINUE
C
      DO 1020 I=1,DN,1
         AI=C(I)
      DO 1020 K=1,DN,1
         AJ=C(K)
C
C     COULOMBPOTENTIAL
C
             H=4.D0*A+3.D0*(AI+AJ)
             HD=AI+AJ
             HE=DSQRT(DSQRT(AI*AJ))**3.D0
       VCE0=SQRT((16384.*PI**12)/(A**9*(AI**3+3.*AI**2*AJ+3.*
     1      AI*AJ**2+AJ**3)))
       VCE0=FVC1*VCE0
       VCE1=SQRT((4.*A*(AI+AJ))/(4.*A+3.*AI+3.*AJ))
         DO 11 Q=1,STANZ+1,1
            X=STEP*DFLOAT(Q)
            HC=E2*VCE0*DERF(VCE1*X)/X
            DO 205 BRA=1,DIS,1
               BI=B(BRA,I)
            DO 205 KET=1,DIS,1
               BJ=B(KET,K)
               KVEC(Q,BRA,KET)=KVEC(Q,BRA,KET)+BI*BJ*HC
205         CONTINUE
11       CONTINUE


C
C     CENTRAL POTENTIAL
C
         DO 120 M=1,PZN,1
            PW=POT(1,M)
            PM=POT(4,M)
            LAM=POT(5,M)
            V0=POT(6,M)
            PB=POT(2,M)
            PH=POT(3,M)
C
          H=(512.*SQRT((AI+AJ)*(A*AI+2*A*LAM+A*AJ+AI*LAM+LAM*AJ)
     1      *A*PI)*SQRT(4*A+AI+AJ)*PI**6)/(SQRT(PI)*SQRT(A)*SQRT(AI+
     2      AJ)*SQRT(16*A**2*AI+16*A**2*LAM+16*A**2*AJ+4*A*AI**2+16*A
     3      *AI*LAM+8*A*AI*AJ+16*A*LAM*AJ+4*A*AJ**2+3*AI**2*LAM+6*AI*
     4      LAM*AJ+3*LAM*AJ**2)*SQRT(A*AI+2*A*LAM+A*AJ+AI*LAM+LAM*AJ)
     5      *A**3*(4.*A*AI+4.*A*LAM+4.*A*AJ+3.*AI*LAM+3.*LAM*AJ))
       VZF0=(8.D0*PW+4.D0*PB-4.D0*PH-2.D0*PM)*H
       VZF1=(-(4.*A*LAM)*(AI+AJ))/(4.*A*AI+4.*A*LAM+4.*A*AJ+3.*
     1      AI*LAM+3.*LAM*AJ)
            DO 12 Q=1,STANZ+1,1
               X=STEP*DFLOAT(Q)
               HF=V0*VZF0*DEXP(VZF1*X*X)
               DO 206 BRA=1,DIS,1
                  BI=B(BRA,I)
               DO 206 KET=1,DIS,1
                  BJ=B(KET,K)
                  KVEC(Q,BRA,KET)=KVEC(Q,BRA,KET)+BI*BJ*HF
206            CONTINUE
12          CONTINUE
120      CONTINUE
C
C     SPIN-ORBIT-POTENTIAL
C
         IF (L.NE.0) THEN
            DO 130 M=1,PSN,1
               LAM=POT(7,M)
               V0=POT(8,M)
      VSA0=(768.*SQRT(2.)*SQRT((AI+AJ)*(A*AI+2.*A*LAM+A*AJ+AI*
     1     LAM+LAM*AJ)*A)*SQRT(4.*A+AI+AJ)*PI**6*(AI+AJ))/(SQRT(A)*
     2     SQRT((A*AI+2.*A*LAM+A*AJ+AI*LAM+LAM*AJ)*2.)*SQRT(AI+AJ)*
     3     SQRT(16.*A**2*AI+16.*A**2*LAM+16.*A**2*AJ+
     3     4.*A*AI**2+16.*A*AI*LAM+8.*A*AI*AJ+
     4     16.*A*LAM*AJ+4.*A*AJ**2+3.*AI**2*LAM+6.*AI*LAM*
     5     AJ+3.*LAM*AJ**2)*A**2*(16.*A**2*AI**2+32.*A**2*AI*LAM+32.*
     6     A**2*AI*AJ+16.*A**2*LAM**2+32.*A**2*LAM*AJ+16.*A**2*AJ**2
     7     +24.*A*AI**2*LAM+24.*A*AI*LAM**2+48.*A*AI*LAM*AJ+24.*A*
     8     LAM**2*AJ+24.*A*LAM*AJ**2+9.*AI**2*LAM**2+18.*AI*LAM**2*
     9     AJ+9.*LAM**2*AJ**2))
      VSA1=(-(4.*A*LAM)*(AI+AJ))/(4.*A*AI+4.*A*LAM+4.*A*AJ+
     1     3.*AI*LAM+3.*LAM*AJ)
      VSA0= 8.D0*VSA0
C
               DO 31 Q=1,STANZ+1,1
                  X=STEP*DFLOAT(Q)
                  HA=V0*VSA0*LKOFF*DEXP(VSA1*X*X)
                  DO 207 BRA=1,DIS,1
                     BI=B(BRA,I)
                  DO 207 KET=1,DIS,1
                     BJ=B(KET,K)
                     KVEC(Q,BRA,KET)=KVEC(Q,BRA,KET)+BI*BJ*HA
207               CONTINUE
31             CONTINUE
130         CONTINUE
         ELSE
         ENDIF
1020  CONTINUE
C
C  *****  MATRIX PART  *****
C
      DO 13 Q=1,STANZ,1
      DO 13 R=1,STANZ,1
         DO 208 BRA=1,DIS,1
         DO 208 KET=1,DIS,1
            KMAT(Q,R,BRA,KET)=0.D0
            KNORM(Q,R,BRA,KET)=0.D0
208      CONTINUE
13    CONTINUE
C
      DO 1030 I=1,DN,1
         AI=C(I)
      DO 1030 K=1,DN,1
         AJ=C(K)
C
C     NORM
C
      NB0=(-65536.*SQRT(PI)*PI**4)/(27.*SQRT(4.*A+3.*AI+3.*AJ)*A
     1    **3*(4.*A+3.*AI+3.*AJ))
      NB1=(2.*(-6.*A**2-5.*A*AI-17.*A*AJ-6.*AI*AJ))/(9.*(4.*A
     1    +3.*AI+3.*AJ))
      NB2=(2.*(-6.*A**2-17.*A*AI-5.*A*AJ-6.*AI*AJ))/(9.*(4.*A
     1    +3.*AI+3.*AJ))
      NB3=(8.*(3.*A**2+A*AI+A*AJ+3.*AI*AJ))/(9.*(4.*A+3.*AI+
     1    3.*AJ))
C
      NC0=(8192.D0*SQRT(PI*2.D0)*PI**4)/(27.D0*SQRT((A**2+A*AI+A*AJ+
     1    AI*AJ)*A)*A*(A**2+A*AI+A*AJ+AI*AJ))
      NC1=(-10.D0*A)/9.D0
      NC2=(-10.D0*A)/9.D0
      NC3=(-16.D0*A)/9.D0
C
         DO 14 Q=1,STANZ,1
            X=STEP*DFLOAT(Q)
         DO 14 R=1,STANZ,1
            Y=STEP*DFLOAT(R)
C
            HB=NB0*S(-NB1,-NB2,-NB3,L,X,Y)
            HC=NC0*S(-NC1,-NC2,-NC3,L,X,Y)
            HILF=HB+HC
C
            DO 209 BRA=1,DIS,1
               BI=B(BRA,I)
            DO 209 KET=1,DIS,1
               BJ=B(KET,K)
      KNORM(Q,R,BRA,KET)=KNORM(Q,R,BRA,KET)+BI*BJ*HILF
209         CONTINUE
14       CONTINUE
C
C     KINETIC ENERGY
C
      call tk(a,ai,aj)
C
         DO 15 Q=1,STANZ,1
            X=STEP*DFLOAT(Q)
         DO 15 R=1,STANZ,1
            Y=STEP*DFLOAT(R)
C
            HB=((TKB1*X*X+TKB3*Y*Y+TKB4)*
     1         S(-TKB5,-TKB6,-TKB7,L,X,Y)+
     2         TKB2*X*Y*T(-TKB5,-TKB6,-TKB7,L,X,Y))
            HC=((TKC1*X*X+TKC3*Y*Y+TKC4)*
     1         S(-TKC5,-TKC6,-TKC7,L,X,Y)+
     2         TKC2*X*Y*T(-TKC5,-TKC6,-TKC7,L,X,Y))
            HILF=HB+HC
C
            DO 210 BRA=1,DIS,1
               BI=B(BRA,I)
            DO 210 KET=1,DIS,1
               BJ=B(KET,K)
               KMAT(Q,R,BRA,KET)=KMAT(Q,R,BRA,KET)-H2M*BI*BJ*HILF
210         CONTINUE
15       CONTINUE
C
C     COULOMBPOTENTIAL
      IF (COU.EQ.1) THEN
C
          H=4.D0*A+3.D0*(AI+AJ)
          HA=4.D0*A+AI+AJ
          HB=A+AI
          HC=A+AJ
          HD=AI+AJ
          HE=DSQRT(PI**9.D0)/(A**3.D0)
          HF=DSQRT((PI**9.D0)/(A**3.D0))
      VCB0=-8192.D0*HE/(9.D0*DSQRT(H))
      VCB1=2.D0*DSQRT(2.D0*A/(HA*H))/3.D0
      VCB2=2.D0*A+AI-AJ
      VCB3=-2.D0*(HB+HD)
      VCB4=-2.D0*(A*(6.D0*A+5.D0*AI+17.D0*AJ)+6.D0*AI*AJ)/(9.D0*H)
      VCB5=-2.D0*(A*(6.D0*A+17.D0*AI+5.D0*AJ)+6.D0*AI*AJ)/(9.D0*H)
      VCB6=8.D0*(3.D0*(A*A+AI*AJ)+A*HD)/(9.D0*H)
          HG=(PI**4.D0)/DSQRT(A**5.D0)
      VCC0=-32768.D0*DSQRT(2.D0)*HG/(27.D0*DSQRT(H)**3.D0)
      VCC4=VCB4
      VCC5=VCB5
      VCC6=VCB6
      VCD0=4096.D0*HF*DSQRT(2.D0/(HB*HC)**3.D0)/9.D0
      VCD1=2.D0*DSQRT(2.D0*A*HC/(3.D0*A+AJ))/3.D0
      VCD2=1.D0
      VCD3=2.D0
      VCD4=-10.D0*A/9.D0
      VCD5=VCD4
      VCD6=-16.D0*A/9.D0
      VCF0=-8192.D0*HE/(DSQRT(H**3.D0)*9.D0)
      VCF4=VCB4
      VCF5=VCB5
      VCF6=VCB6
      VCG0=-16384.D0*HE/(9.D0*DSQRT(H))
      VCG1=2.D0*DSQRT(2.D0*A/((HB+HC)*H))/3.D0
      VCG2=3.D0*A+2.D0*AI+AJ
      VCG3=3.D0*A+AI+2.D0*AJ
      VCG4=VCB4
      VCG5=VCB5
      VCG6=VCB6
      VCH0=-8192.D0*HE/(9.D0*DSQRT(H))
      VCH1=2.D0*DSQRT(2.D0*A/(HA*H))/3.D0
      VCH2=-2.D0*(HC+HD)
      VCH3=2.D0*A-AI+AJ
      VCH4=VCB4
      VCH5=VCB5
      VCH6=VCB6
      VCI0=4096.D0*HF*DSQRT(2.D0)/(9.D0*DSQRT((HB*HC)**3.D0))
      VCI1=2.D0*DSQRT(2.D0*A*HB/(3.D0*A+AI))/3.D0
      VCI2=2.D0
      VCI3=1.D0
      VCI4=VCD4
      VCI5=VCD5
      VCI6=VCD6
      VCJ0=4096.D0*HF*DSQRT(2.D0)/(9.D0*DSQRT((HB*HC)**3.D0))
      VCJ1=2.D0*DSQRT(HB*HC/(HB+HC))/3.D0
      VCJ2=1.D0
      VCJ3=-1.D0
      VCJ4=VCD4
      VCJ5=VCD5
      VCJ6=VCD6
        HG=PI**4.D0/(A*HB*HC)
      VCK0=8192.D0*HG*DSQRT(2.D0/(A*HB))/27.D0
      VCK4=VCD4
      VCK5=VCD5
      VCK6=VCD6
      VCL0=8192*HE/(9.D0*DSQRT(H))
      VCL0=-2.D0*VCL0
      VCL1=2.D0/(3.D0*DSQRT(3.D0*H))
      VCL2=5.D0*A+3.D0*AI
      VCL3=A-3.D0*AI
      VCL4=VCB4
      VCL5=VCB5
      VCL6=VCB6
      VCN0=8192*HE/(9.D0*DSQRT(H))
      VCN0=-2.D0*VCN0
      VCN1=2.D0/(3.D0*DSQRT(3.D0*H))
      VCN2=A-3.D0*AJ
      VCN3=5.D0*A+3.D0*AJ
      VCN4=VCB4
      VCN5=VCB5
      VCN6=VCB6
      VCO0=8192.D0*HG*DSQRT(2.D0/(A*HC))/27.D0
      VCO4=VCD4
      VCO5=VCD5
      VCO6=VCD6
      VCP0=16384.D0*PI**4.D0/(27.D0*A*DSQRT(HB*HC)**3.D0)
      VCP4=VCD4
      VCP5=VCD5
      VCP6=VCD6
C
      VCB0=VCB0*FVC1
      VCC0=VCC0*FVC2
      VCD0=VCD0*FVC1*FVC2
      VCE0=VCE0*FVC1
      VCF0=VCF0*FVC1
      VCG0=VCG0*(FVC1+FVC2)/2.D0
      VCH0=VCH0*FVC1
      VCI0=VCI0*FVC1*FVC2
      VCJ0=VCJ0*FVC1*FVC1
      VCK0=VCK0*FVC3
      VCL0=VCL0*FVC3
      VCN0=VCN0*FVC3
      VCO0=VCO0*FVC3
      VCP0=VCP0*FVC4
C
         DO 16 Q=1,STANZ,1
            X=STEP*DFLOAT(Q)
         DO 16 R=1,STANZ,1
            Y=STEP*DFLOAT(R)
C
            HB=E2*2.D0*PI*VCB0*VCB1*
     1         COULI(-VCB4,-VCB5,VCB1*VCB2,VCB1*VCB3,-VCB6,L,X,Y,1)
            HC=E2*VCC0*S(-VCC4,-VCC5,-VCC6,L,X,Y)
            HD=E2*2.D0*PI*VCD0*VCD1*
     1         COULI(-VCD4,-VCD5,VCD1*VCD2,VCD1*VCD3,-VCD6,L,X,Y,1)
            HF=E2*2.D0*PI*VCF0*
     1         COULI(-VCF4,-VCF5,1.D0,-1.D0,-VCF6,L,X,Y,2)
            HG=E2*2.D0*PI*VCG0*VCG1*
     1         COULI(-VCG4,-VCG5,VCG1*VCG2,VCG1*VCG3,-VCG6,L,X,Y,1)
            HH=E2*2.D0*PI*VCH0*VCH1*
     1         COULI(-VCH4,-VCH5,VCH1*VCH2,VCH1*VCH3,-VCH6,L,X,Y,1)
            HI=E2*2.D0*PI*VCI0*VCI1*
     1         COULI(-VCI4,-VCI5,VCI1*VCI2,VCI1*VCI3,-VCI6,L,X,Y,1)
            HJ=E2*2.D0*PI*VCJ0*VCJ1*
     1         COULI(-VCJ4,-VCJ5,VCJ1*VCJ2,VCJ1*VCJ3,-VCJ6,L,X,Y,1)
            HK=E2*VCK0*S(-VCK4,-VCK5,-VCK6,L,X,Y)
            HL=E2*2.D0*PI*VCL0*VCL1*
     1         COULI(-VCL4,-VCL5,VCL1*VCL2,VCL1*VCL3,-VCL6,L,X,Y,1)
            HN=E2*2.D0*PI*VCN0*VCN1*
     1         COULI(-VCN4,-VCN5,VCN1*VCN2,VCN1*VCN3,-VCN6,L,X,Y,1)
            HO=E2*VCO0*S(-VCO4,-VCO5,-VCO6,L,X,Y)
            HP=E2*VCP0*S(-VCP4,-VCP5,-VCP6,L,X,Y)
            HILF=HB+HC+HD+HF+HG+HH+HI+HJ+HK+HL+HN+HO+HP
C
            DO 211 BRA=1,DIS,1
               BI=B(BRA,I)
            DO 211 KET=1,DIS,1
               BJ=B(KET,K)
               KMAT(Q,R,BRA,KET)=KMAT(Q,R,BRA,KET)+BI*BJ*HILF
211         CONTINUE
16       CONTINUE
      ELSE
      ENDIF
C
C     CENTRAL POTENTIAL
C
         DO 121 M=1,PZN,1
            PW=POT(1,M)
            PM=POT(4,M)
            LAM=POT(5,M)
            V0=POT(6,M)
            PB=POT(2,M)
            PH=POT(3,M)
C
      call vz1(pb,pm,ph,pw,a,ai,aj,lam)
      call vz2(pb,pm,ph,pw,a,ai,aj,lam)
C
            DO 17 Q=1,STANZ,1
               X=STEP*DFLOAT(Q)
            DO 17 R=1,STANZ,1
               Y=STEP*DFLOAT(R)
C
               HB=V0*VZB0*S(-VZB1,-VZB2,-VZB3,L,X,Y)
               HC=V0*VZC0*S(-VZC1,-VZC2,-VZC3,L,X,Y)
               HD=V0*VZD0*S(-VZD1,-VZD2,-VZD3,L,X,Y)
               HE=V0*VZE0*S(-VZE1,-VZE2,-VZE3,L,X,Y)
               HG=V0*VZG0*S(-VZG1,-VZG2,-VZG3,L,X,Y)
               HH=V0*VZH0*S(-VZH1,-VZH2,-VZH3,L,X,Y)
               HI=V0*VZI0*S(-VZI1,-VZI2,-VZI3,L,X,Y)
               HJ=V0*VZJ0*S(-VZJ1,-VZJ2,-VZJ3,L,X,Y)
               HK=V0*VZK0*S(-VZK1,-VZK2,-VZK3,L,X,Y)
               HL=V0*VZL0*S(-VZL1,-VZL2,-VZL3,L,X,Y)
               HM=V0*VZM0*S(-VZM1,-VZM2,-VZM3,L,X,Y)
               HO=V0*VZO0*S(-VZO1,-VZO2,-VZO3,L,X,Y)
               HP=V0*VZP0*S(-VZP1,-VZP2,-VZP3,L,X,Y)
               HILF=HB+HC+HD+HE+HG+HH+HI+HJ+HK+HL+HM+HO+HP
C
               DO 212 BRA=1,DIS,1
                  BI=B(BRA,I)
               DO 212 KET=1,DIS,1
                  BJ=B(KET,K)
      KMAT(Q,R,BRA,KET)=KMAT(Q,R,BRA,KET)+BI*BJ*HILF
C
212            CONTINUE
17          CONTINUE
121      CONTINUE
C
C     SPIN-ORBIT-POTENTIAL
C
         IF (L.NE.0) THEN
            DO 131 M=1,PSN,1
               LAM=POT(7,M)
               V0=POT(8,M)
C
      call vs(a,ai,aj,lam)
C
               DO 32 Q=1,STANZ,1
                  X=STEP*DFLOAT(Q)
               DO 32 R=1,STANZ,1
                  Y=STEP*DFLOAT(R)
C
               HB=V0*LKOFF*VSB0/VSB3*S(-VSB1,-VSB2,-VSB3,L,X,Y)
               HC=V0*LKOFF*VSC0/VSC3*S(-VSC1,-VSC2,-VSC3,L,X,Y)
               HD=V0*LKOFF*VSD0/VSD3*S(-VSD1,-VSD2,-VSD3,L,X,Y)
               HILF=HB+HC+HD
C
                  DO 213 BRA=1,DIS,1
                     BI=B(BRA,I)
                  DO 213 KET=1,DIS,1
                     BJ=B(KET,K)
                     KMAT(Q,R,BRA,KET)=KMAT(Q,R,BRA,KET)+BI*BJ*HILF
213               CONTINUE
32             CONTINUE
131         CONTINUE
         ELSE
         ENDIF
1030  CONTINUE
C
C *** DIVISION BY NORMALIZATION FACTOR ***
C
      DO 214 BRA=1,DIS,1
      DO 214 KET=1,DIS,1
         KNA0=DSQRT(NA0(BRA,BRA)*NA0(KET,KET))
         DO 18 Q=1,STANZ+1,1
            KVEC(Q,BRA,KET)=KVEC(Q,BRA,KET)/KNA0
18       CONTINUE
         DO 19 Q=1,STANZ,1
            DO 20 R=1,STANZ,1
               KMAT(Q,R,BRA,KET)=KMAT(Q,R,BRA,KET)/KNA0
               KNORM(Q,R,BRA,KET)=KNORM(Q,R,BRA,KET)/KNA0
20          CONTINUE
19       CONTINUE
214   CONTINUE
C
C *** STORAGE OF THE MATRIXELEMENTS ***
C
C
       OUTNAME(1)='alphad.ou1'
       OUTNAME(2)='alphad.ou2'
       OUTNAME(3)='alphad.ou3'
       OUTNAME(4)='alphad.ou4'
       OUTNAME(5)='alphad.ou5'
       OUTNAME(6)='alphad.ou6'
C
      KK=1
      FIN=20
      DO 200 BRA=1,DIS
      DO 200 KET=BRA,DIS
         KENE(BRA,KET)=ENEA+KENED(BRA,KET)
         FIN=FIN+KK
C
         OPEN(FIN,FILE=OUTNAME(KK),
     &        FORM='UNFORMATTED')
C
         IF (BRA.EQ.KET) THEN
C
C *** PRINT OUT THE BINDING ENERGIES ***
C
            WRITE(FIN) KENE(BRA,BRA)
            WRITE(6,*) 'BINDING ENERGIES:'
            WRITE(6,*) 'CHANNEL NR. ',BRA
            WRITE(6,*) ' TOTAL          ',KENE(BRA,BRA)
            WRITE(6,*) ' ALPHA-PARTICLE ',ENEA
            WRITE(6,*) ' DEUTERON       ',KENED(BRA,BRA)
         ENDIF
C
         DO 40 Q=1,STANZ+1
            WRITE(FIN) KVEC(Q,BRA,KET)
40       CONTINUE
C
         START=1
         DO 41 Q=1,STANZ
         IF (BRA.EQ.KET) START=Q
         DO 41 R=START,STANZ
            WRITE(FIN) KMAT(Q,R,BRA,KET),KNORM(Q,R,BRA,KET)
41       CONTINUE
C
      KK=KK+1
C
200   CONTINUE
      RETURN
      END
C
      subroutine tk(a,ai,aj)
C  Subroutine for the kinetic energy
      IMPLICIT REAL*8 (A-Z)
      common /tkcom/tkb1,tkb2,tkb3,tkb4,tkb5,tkb6,tkb7,
     &              tkc1,tkc2,tkc3,tkc4,tkc5,tkc6,tkc7
      PI    = 3.141592653589793D0
C
      TKB1=(524288.*SQRT(4*A+AI+AJ)*SQRT(A**2+A*AI+A*AJ+AI*
     1      AJ)*PI**5*(-36.*A**4-48.*A**3*AI-84.*A**3*AJ-20.*A**2*AI
     2      **2-112.*A**2*AI*AJ-149.*A**2*AJ**2-12.*A*AI**2*AJ-102.*A
     3      *AI*AJ**2-45.*AI**2*AJ**2))/(729.*SQRT(PI)*SQRT(A+AJ)*
     4      SQRT(A+AI)*SQRT(16*A**2+16*A*AI+16*A*AJ+3*AI**2+6*AI*AJ+3
     5      *AJ**2)*A**3*(64.*A**3+144.*A**2*AI+144.*A**2*AJ+108.*A*
     6      AI**2+216.*A*AI*AJ+108.*A*AJ**2+27.*AI**3+81.*AI**2*AJ+
     7      81.*AI*AJ**2+27.*AJ**3))
      TKB2=(1048576.*SQRT(4*A+AI+AJ)*SQRT(A**2+A*AI+A*AJ+AI
     1      *AJ)*PI**5*(36.*A**4+66.*A**3*AI+66.*A**3*AJ+8.*A**2*AI**
     2      2+49.*A**2*AI*AJ+8.*A**2*AJ**2+57.*A*AI**2*AJ+57.*A*AI*AJ
     3      **2+45.*AI**2*AJ**2))/(729.*SQRT(PI)*SQRT(A+AJ)*SQRT(A+AI
     4      )*SQRT(16*A**2+16*A*AI+16*A*AJ+3*AI**2+6*AI*AJ+3*AJ**2)*A
     5      **3*(64.*A**3+144.*A**2*AI+144.*A**2*AJ+108.*A*AI**2+216.
     6      *A*AI*AJ+108.*A*AJ**2+27.*AI**3+81.*AI**2*AJ+81.*AI*AJ**2
     7      +27.*AJ**3))
      TKB3=(524288.*SQRT(4*A+AI+AJ)*SQRT(A**2+A*AI+A*AJ+AI*
     1      AJ)*PI**5*(-36.*A**4-84.*A**3*AI-48.*A**3*AJ-149.*A**2*AI
     2      **2-112.*A**2*AI*AJ-20.*A**2*AJ**2-102.*A*AI**2*AJ-12.*A*
     3      AI*AJ**2-45.*AI**2*AJ**2))/(729.*SQRT(PI)*SQRT(A+AJ)*SQRT
     4      (A+AI)*SQRT(16*A**2+16*A*AI+16*A*AJ+3*AI**2+6*AI*AJ+3*AJ
     5      **2)*A**3*(64.*A**3+144.*A**2*AI+144.*A**2*AJ+108.*A*AI**
     6      2+216.*A*AI*AJ+108.*A*AJ**2+27.*AI**3+81.*AI**2*AJ+81.*AI
     7      *AJ**2+27.*AJ**3))
      TKB4=(65536.*SQRT(4*A+AI+AJ)*SQRT(A**2+A*AI+A*AJ+AI*
     1      AJ)*PI**5*(24.*A**2+35.*A*AI+35.*A*AJ+15.*AI*AJ))/(27.*
     2      SQRT(PI)*SQRT(A+AJ)*SQRT(A+AI)*SQRT(16*A**2+16*A*AI+16*A*
     3      AJ+3*AI**2+6*AI*AJ+3*AJ**2)*A**3*(16.*A**2+24.*A*AI+24.*A
     4      *AJ+9.*AI**2+18.*AI*AJ+9.*AJ**2))
      TKB5=(2.*(-6.*A**2-5.*A*AI-17.*A*AJ-6.*AI*AJ))/(9.*(
     1      4.*A+3.*AI+3.*AJ))
      TKB6=(2.*(-6.*A**2-17.*A*AI-5.*A*AJ-6.*AI*AJ))/(9.*(
     1      4.*A+3.*AI+3.*AJ))
      TKB7=(8.*(3.*A**2+A*AI+A*AJ+3.*AI*AJ))/(9.*(4.*A+3.*
     1      AI+3.*AJ))
C
      TKC1=(1835008.*SQRT(A)*SQRT(A**2+A*AI+A*AJ+AI*AJ)*PI
     1      **5)/(729.*SQRT(2.)*SQRT(PI)*(A**4+2.*A**3*AI+2.*A**3*AJ+A
     2      **2*AI**2+4.*A**2*AI*AJ+A**2*AJ**2+2.*A*AI**2*AJ+2.*A*AI*
     3      AJ**2+AI**2*AJ**2))
      TKC2=(3407872.*SQRT(A)*SQRT(A**2+A*AI+A*AJ+AI*AJ)*PI
     1      **5)/(729.*SQRT(2.)*SQRT(PI)*(A**4+2.*A**3*AI+2.*A**3*AJ+A
     2      **2*AI**2+4.*A**2*AI*AJ+A**2*AJ**2+2.*A*AI**2*AJ+2.*A*AI*
     3      AJ**2+AI**2*AJ**2))
      TKC3=(1835008.*SQRT(A)*SQRT(A**2+A*AI+A*AJ+AI*AJ)*PI
     1      **5)/(729.*SQRT(2.)*SQRT(PI)*(A**4+2.*A**3*AI+2.*A**3*AJ+A
     2      **2*AI**2+4.*A**2*AI*AJ+A**2*AJ**2+2.*A*AI**2*AJ+2.*A*AI*
     3      AJ**2+AI**2*AJ**2))
      TKC4=(8192.*PI**5*(-19.*SQRT(A)*SQRT((4*A+AI+AJ)*2)*
     1      SQRT(A**2+A*AI+A*AJ+AI*AJ)*A**2-23.*SQRT(A)*SQRT((4*A+AI+
     2      AJ)*2)*SQRT(A**2+A*AI+A*AJ+AI*AJ)*A*AI-25.*SQRT(A)*SQRT((
     3      4*A+AI+AJ)*2)*SQRT(A**2+A*AI+A*AJ+AI*AJ)*A*AJ-29.*SQRT(A)
     4      *SQRT((4*A+AI+AJ)*2)*SQRT(A**2+A*AI+A*AJ+AI*AJ)*AI*AJ-2.*
     5      SQRT(2.)*SQRT((A**2+A*AI+A*AJ+AI*AJ)*A)*SQRT(4*A+AI+AJ)*A*
     6      AI-2.*SQRT(2.)*SQRT((A**2+A*AI+A*AJ+AI*AJ)*A)*SQRT(4*A+AI+
     7      AJ)*AI*AJ))/(27.*SQRT(2.)*SQRT(PI)*SQRT((4*A+AI+AJ)*2)*A*(
     8      A**6+3.*A**5*AI+3.*A**5*AJ+3.*A**4*AI**2+9.*A**4*AI*AJ+3.
     9      *A**4*AJ**2+A**3*AI**3+9.*A**3*AI**2*AJ+9.*A**3*AI*AJ**2+
     1      A**3*AJ**3+3.*A**2*AI**3*AJ+9.*A**2*AI**2*AJ**2+3.*A**2*
     2      AI*AJ**3+3.*A*AI**3*AJ**2+3.*A*AI**2*AJ**3+AI**3*AJ**3))
      TKC5=(-10.*A)/9.
      TKC6=(-10.*A)/9.
      TKC7=(-16.*A)/9.
      return
      end
C
      subroutine vz1(pb,pm,ph,pw,a,ai,aj,lam)
C  subroutine for the central potential
      IMPLICIT REAL*8 (A-Z)
      common /vzcom/vzb0,vzb1,vzb2,vzb3,
     &              vzc0,vzc1,vzc2,vzc3,
     &              vzd0,vzd1,vzd2,vzd3,
     &              vze0,vze1,vze2,vze3,
     &              vzg0,vzg1,vzg2,vzg3,
     &              vzh0,vzh1,vzh2,vzh3,
     &              vzi0,vzi1,vzi2,vzi3,
     &              vzj0,vzj1,vzj2,vzj3,
     &              vzk0,vzk1,vzk2,vzk3,
     &              vzl0,vzl1,vzl2,vzl3,
     &              vzm0,vzm1,vzm2,vzm3,
     &              vzo0,vzo1,vzo2,vzo3,
     &              vzp0,vzp1,vzp2,vzp3
      PI    = 3.141592653589793D0
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT((4*A**2+3*A*AI+8*
     1      A*LAM+3*A*AJ+2*AI*LAM+2*LAM*AJ)*A)*A*(4.*A**2+3.*A*AI+8.*
     2      A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZB0=-6.D0*(PW+PM)*H
      VZB1=(2.*(-6.*A**3-5.*A**2*AI-20.*A**2*LAM-17.*A**2*AJ-
     1      6.*A*AI*LAM-6.*A*AI*AJ-14.*A*LAM*AJ-4.*AI*LAM*AJ))/(9.*(4.
     2      *A**2+3.*A*AI+8.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZB2=(2.*(-6.*A**3-17.*A**2*AI-20.*A**2*LAM-5.*A**2*AJ-
     1      54.*A*AI*LAM-6.*A*AI*AJ-14.*A*LAM*AJ-4.*AI*LAM*AJ))/(9.*(
     2      4.*A**2+3.*A*AI+8.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZB3=(8.*(3.*A**3+A**2*AI+10.*A**2*LAM+A**2*AJ+6.*A*AI*
     1      LAM+3.*A*AI*AJ-2.*A*LAM*AJ+2.*AI*LAM*AJ))/(9.*(4.*A**2+3.
     2      *A*AI+8.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
C
          H=(8192.*SQRT(PI*2)*PI**4)/(27.*SQRT((A**2+A*AI+4*A*
     1      LAM+A*AJ+4*AI*LAM+AI*AJ)*A)*A*(A**2+A*AI+4.*A*LAM+A*AJ+4.
     2      *AI*LAM+AI*AJ))
      VZC0=(PW+PB+PH+PM)*H
      VZC1=(-10.*A)/9.
      VZC2=(-10.*A)/9.
      VZC3=(-16.*A)/9.
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT((4*A**2+3*A*AI+8*
     1      A*LAM+3*A*AJ+6*AI*LAM+6*LAM*AJ)*A)*A*(4.*A**2+3.*A*AI+8.*
     2      A*LAM+3.*A*AJ+6.*AI*LAM+6.*LAM*AJ))
      VZD0=-6.D0*(PW+PM)*H
      VZD1=(2.*(-6.*A**2-5.*A*AI-17.*A*AJ-6.*AI*AJ))/(9.*(4.*A
     1      +3.*AI+3.*AJ))
      VZD2=(2.*(-6.*A**2-17.*A*AI-5.*A*AJ-6.*AI*AJ))/(9.*(4.*A
     1      +3.*AI+3.*AJ))
      VZD3=(8.*(3.*A**2+A*AI+A*AJ+3.*AI*AJ))/(9.*(4.*A+3.*AI+
     1      3.*AJ))
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT(2*A**3+2*A**2*AI+
     1      3*A**2*LAM+2*A**2*AJ+3*A*AI*LAM+2*A*AI*AJ+A*LAM*AJ+AI*LAM
     2      *AJ)*(2.*A**3+2.*A**2*AI+3.*A**2*LAM+2.*A**2*AJ+3.*A*AI*
     3      LAM+2.*A*AI*AJ+A*LAM*AJ+AI*LAM*AJ))
      VZE0=(4.D0*PW-2.D0*PB-2.D0*PH+4.D0*PM)*H
      VZE1=(2.*A*(-10.*A**2-19.*A*LAM-10.*A*AJ-9.*LAM*AJ))/(9.
     1      *(2.*A**2+3.*A*LAM+2.*A*AJ+LAM*AJ))
      VZE2=(2.*A*(-10.*A**2-31.*A*LAM-10.*A*AJ-21.*LAM*AJ))/(
     1      9.*(2.*A**2+3.*A*LAM+2.*A*AJ+LAM*AJ))
      VZE3=(16.*A*(-2.*A**2-5.*A*LAM-2.*A*AJ-3.*LAM*AJ))/(9.*(
     1      2.*A**2+3.*A*LAM+2.*A*AJ+LAM*AJ))
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT((4*A**2+3*A*AI+4*
     1      A*LAM+3*A*AJ+2*AI*LAM+2*LAM*AJ)*A)*A*(4.*A**2+3.*A*AI+4.*
     2      A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZG0=(-6.D0*PW-2.D0*PB+6.D0*PH+4.D0*PM)*H
      VZG1=(2.*(-6.*A**3-5.*A**2*AI-24.*A**2*LAM-17.*A**2*AJ-
     1      14.*A*AI*LAM-6.*A*AI*AJ-14.*A*LAM*AJ-4.*AI*LAM*AJ))/(9.*(
     2      4.*A**2+3.*A*AI+4.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZG2=(2.*(-6.*A**3-17.*A**2*AI-24.*A**2*LAM-5.*A**2*AJ-
     1      14.*A*AI*LAM-6.*A*AI*AJ-14.*A*LAM*AJ-4.*AI*LAM*AJ))/(9.*(
     2      4.*A**2+3.*A*AI+4.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZG3=(8.*(3.*A**3+A**2*AI-6.*A**2*LAM+A**2*AJ-2.*A*AI*
     1      LAM+3.*A*AI*AJ-2.*A*LAM*AJ+2.*AI*LAM*AJ))/(9.*(4.*A**2+3.
     2      *A*AI+4.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT(4*A+3*AI+3*AJ)*A
     1      **3*(4.*A+3.*AI+3.*AJ))
      VZH0=(-2.D0*PW-4.D0*PB+4.D0*PH+8.D0*PM)*H
      VZH1=(2.*(-6.*A**2-5.*A*AI-32.*A*LAM-17.*A*AJ-24.*AI*LAM
     1      -6.*AI*AJ-24.*LAM*AJ))/(9.*(4.*A+3.*AI+3.*AJ))
      VZH2=(2.*(-6.*A**2-17.*A*AI-32.*A*LAM-5.*A*AJ-24.*AI*LAM
     1      -6.*AI*AJ-24.*LAM*AJ))/(9.*(4.*A+3.*AI+3.*AJ))
      VZH3=(8.*(3.*A**2+A*AI+16.*A*LAM+A*AJ+12.*AI*LAM+3.*AI*
     1      AJ+12.*LAM*AJ))/(9.*(4.*A+3.*AI+3.*AJ))
      return
      end
C
      subroutine vz2(pb,pm,ph,pw,a,ai,aj,lam)
C  subroutine for the central potential
      IMPLICIT REAL*8 (A-Z)
      common /vzcom/vzb0,vzb1,vzb2,vzb3,
     &              vzc0,vzc1,vzc2,vzc3,
     &              vzd0,vzd1,vzd2,vzd3,
     &              vze0,vze1,vze2,vze3,
     &              vzg0,vzg1,vzg2,vzg3,
     &              vzh0,vzh1,vzh2,vzh3,
     &              vzi0,vzi1,vzi2,vzi3,
     &              vzj0,vzj1,vzj2,vzj3,
     &              vzk0,vzk1,vzk2,vzk3,
     &              vzl0,vzl1,vzl2,vzl3,
     &              vzm0,vzm1,vzm2,vzm3,
     &              vzo0,vzo1,vzo2,vzo3,
     &              vzp0,vzp1,vzp2,vzp3
      PI    = 3.141592653589793D0
C
          H=(8192.*SQRT(PI*2)*PI**4)/(27.*SQRT((A**2+A*AI+2*A*
     1      LAM+A*AJ+AI*LAM+AI*AJ+LAM*AJ)*A)*A*(A**2+A*AI+2.*A*LAM+A*
     2      AJ+AI*LAM+AI*AJ+LAM*AJ))
      VZI0=(4.D0*PW+6.D0*PB-2.D0*PH-6.D0*PM)*H
      VZI1=(2.*(-5.*A**3-5.*A**2*AI-12.*A**2*LAM-5.*A**2*AJ-7.
     1      *A*AI*LAM-5.*A*AI*AJ-7.*A*LAM*AJ-2.*AI*LAM*AJ))/(9.*(A**2
     2      +A*AI+2.*A*LAM+A*AJ+AI*LAM+AI*AJ+LAM*AJ))
      VZI2=(2.*(-5.*A**3-5.*A**2*AI-12.*A**2*LAM-5.*A**2*AJ-7.
     1      *A*AI*LAM-5.*A*AI*AJ-7.*A*LAM*AJ-2.*AI*LAM*AJ))/(9.*(A**2
     2      +A*AI+2.*A*LAM+A*AJ+AI*LAM+AI*AJ+LAM*AJ))
      VZI3=(8.*(-2.*A**3-2.*A**2*AI-3.*A**2*LAM-2.*A**2*AJ-A*
     1      AI*LAM-2.*A*AI*AJ-A*LAM*AJ+AI*LAM*AJ))/(9.*(A**2+A*AI+2.*
     2      A*LAM+A*AJ+AI*LAM+AI*AJ+LAM*AJ))
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT((4*A**2+3*A*AI+8*
     1      A*LAM+3*A*AJ+2*AI*LAM+2*LAM*AJ)*A)*A*(4.*A**2+3.*A*AI+8.*
     2      A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZJ0=-6.D0*(PW+PM)*H
      VZJ1=(2.*(-6.*A**3-5.*A**2*AI-20.*A**2*LAM-17.*A**2*AJ-
     1      14.*A*AI*LAM-6.*A*AI*AJ-54.*A*LAM*AJ-4.*AI*LAM*AJ))/(9.*(
     2      4.*A**2+3.*A*AI+8.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZJ2=(2.*(-6.*A**3-17.*A**2*AI-20.*A**2*LAM-5.*A**2*AJ-
     1      14.*A*AI*LAM-6.*A*AI*AJ-6.*A*LAM*AJ-4.*AI*LAM*AJ))/(9.*(4.
     2      *A**2+3.*A*AI+8.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VZJ3=(8.*(3.*A**3+A**2*AI+10.*A**2*LAM+A**2*AJ-2.*A*AI*
     1      LAM+3.*A*AI*AJ+6.*A*LAM*AJ+2.*AI*LAM*AJ))/(9.*(4.*A**2+3.
     2      *A*AI+8.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT(4*A+3*AI+12*LAM+3
     1      *AJ)*A**3*(4.*A+3.*AI+12.*LAM+3.*AJ))
      VZK0=-2.D0*(PW+PB+PH+PM)*H
      VZK1=(2.*(-6.*A**2-5.*A*AI-68.*A*LAM-17.*A*AJ-24.*AI*LAM
     1      -6.*AI*AJ))/(9.*(4.*A+3.*AI+12.*LAM+3.*AJ))
      VZK2=(2.*(-6.*A**2-17.*A*AI-20.*A*LAM-5.*A*AJ-24.*AI*LAM
     1      -6.*AI*AJ))/(9.*(4.*A+3.*AI+12.*LAM+3.*AJ))
      VZK3=(8.*(3.*A**2+A*AI+4.*A*LAM+A*AJ+12.*AI*LAM+3.*AI*AJ
     1      ))/(9.*(4.*A+3.*AI+12.*LAM+3.*AJ))
C
          H=(8192.*SQRT(PI*2)*PI**4)/(27.*SQRT(A**3+A**2*AI+2*A
     1      **2*LAM+A**2*AJ+2*A*AI*LAM+A*AI*AJ+2*A*LAM*AJ+2*AI*LAM*AJ
     2      )*(A**3+A**2*AI+2.*A**2*LAM+A**2*AJ+2.*A*AI*LAM+A*AI*AJ+
     3      2.*A*LAM*AJ+2.*AI*LAM*AJ))
      VZL0=(PW+PB+PH+PM)*H
      VZL1=(-10.*A)/9.
      VZL2=(-10.*A)/9.
      VZL3=(-16.*A)/9.
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT(2*A**3+2*A**2*AI+
     1      3*A**2*LAM+2*A**2*AJ+A*AI*LAM+2*A*AI*AJ+3*A*LAM*AJ+AI*LAM
     2      *AJ)*(2.*A**3+2.*A**2*AI+3.*A**2*LAM+2.*A**2*AJ+A*AI*LAM+
     3      2.*A*AI*AJ+3.*A*LAM*AJ+AI*LAM*AJ))
      VZM0=(4.D0*PW-2.D0*PB-2.D0*PH+4.D0*PM)*H
      VZM1=(2.*A*(-10.*A**2-10.*A*AI-31.*A*LAM-21.*AI*LAM))/(
     1      9.*(2.*A**2+2.*A*AI+3.*A*LAM+AI*LAM))
      VZM2=(2.*A*(-10.*A**2-10.*A*AI-19.*A*LAM-9.*AI*LAM))/(9.
     1      *(2.*A**2+2.*A*AI+3.*A*LAM+AI*LAM))
      VZM3=(16.*A*(-2.*A**2-2.*A*AI-5.*A*LAM-3.*AI*LAM))/(9.*(
     1      2.*A**2+2.*A*AI+3.*A*LAM+AI*LAM))
C
          H=(32768.*SQRT(PI)*PI**4)/(27.*SQRT(4*A+3*AI+12*LAM+3
     1      *AJ)*A**3*(4.*A+3.*AI+12.*LAM+3.*AJ))
      VZO0=-2.D0*(PW+PB+PH+PM)*H
      VZO1=(2.*(-6.*A**2-5.*A*AI-20.*A*LAM-17.*A*AJ-6.*AI*AJ-
     1      24.*LAM*AJ))/(9.*(4.*A+3.*AI+12.*LAM+3.*AJ))
      VZO2=(2.*(-6.*A**2-17.*A*AI-68.*A*LAM-5.*A*AJ-6.*AI*AJ-
     1      24.*LAM*AJ))/(9.*(4.*A+3.*AI+12.*LAM+3.*AJ))
      VZO3=(8.*(3.*A**2+A*AI+4.*A*LAM+A*AJ+3.*AI*AJ+12.*LAM*AJ
     1      ))/(9.*(4.*A+3.*AI+12.*LAM+3.*AJ))
C
          H=(8192.*SQRT(PI*2)*PI**4)/(27.*SQRT((A**2+A*AI+4*A*
     1      LAM+A*AJ+AI*AJ+4*LAM*AJ)*A)*A*(A**2+A*AI+4.*A*LAM+A*AJ+AI
     2      *AJ+4.*LAM*AJ))
      VZP0=(PW+PB+PH+PM)*H
      VZP1=(-10.*A)/9.
      VZP2=(-10.*A)/9.
      VZP3=(-16.*A)/9.
      return
      end
C
      subroutine vs(a,ai,aj,lam)
C  subroutine for the spin-orbit-potential
      IMPLICIT REAL*8 (A-Z)
      common /vscom/vsb0,vsb1,vsb2,vsb3,
     &              vsc0,vsc1,vsc2,vsc3,
     &              vsd0,vsd1,vsd2,vsd3
      PI    = 3.141592653589793D0
C
      VSB0=(262144.*SQRT(PI)*PI**4*(A**2+A*AI+A*AJ+AI*AJ))/(
     1     81.*SQRT((4.*A**2+3.*A*AI+4.*A*LAM+3.*A*AJ+
     1     2.*AI*LAM+2.*LAM*AJ)*A)*
     2     (16.*A**4+24.*A**3*AI+32.*A**3*LAM+24.*A**3*AJ+9.*A**2
     3     *AI**2+40.*A**2*AI*LAM+18.*A**2*AI*AJ+16.*A**2*LAM**2+40.
     4     *A**2*LAM*AJ+9.*A**2*AJ**2+12.*A*AI**2*LAM+16.*A*AI*LAM**
     5     2+24.*A*AI*LAM*AJ+16.*A*LAM**2*AJ+12.*A*LAM*AJ**2+4.*AI**
     6     2*LAM**2+8.*AI*LAM**2*AJ+4.*LAM**2*AJ**2))
      VSB1=(2.*(-6.*A**3-5.*A**2*AI-24.*A**2*LAM-17.*A**2*AJ
     1     -14.*A*AI*LAM-6.*A*AI*AJ-14.*A*LAM*AJ-4.*AI*LAM*AJ))/(9.*
     2     (4.*A**2+3.*A*AI+4.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VSB2=(2.*(-6.*A**3-17.*A**2*AI-24.*A**2*LAM-5.*A**2*AJ
     1     -14.*A*AI*LAM-6.*A*AI*AJ-14.*A*LAM*AJ-4.*AI*LAM*AJ))/(9.*
     2     (4.*A**2+3.*A*AI+4.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
      VSB3=(8.*(3.*A**3+A**2*AI-6.*A**2*LAM+A**2*AJ-2.*A*AI*
     1     LAM+3.*A*AI*AJ-2.*A*LAM*AJ+2.*AI*LAM*AJ))/(9.*(4.*A**2+3.
     2     *A*AI+4.*A*LAM+3.*A*AJ+2.*AI*LAM+2.*LAM*AJ))
C
      VSC0=-(524288.*SQRT(PI)*PI**4*(AI+AJ))/(81.*SQRT(4.*A+3.*
     1     AI+3.*AJ)*A**2*(16.*A**2+24.*A*AI+24.*A*AJ+9.*AI**2+18.*AI
     2     *AJ+9.*AJ**2))
      VSC1=(2.*(-6.*A**2-5.*A*AI-32.*A*LAM-17.*A*AJ-24.*AI*
     1     LAM-6.*AI*AJ-24.*LAM*AJ))/(9.*(4.*A+3.*AI+3.*AJ))
      VSC2=(2.*(-6.*A**2-17.*A*AI-32.*A*LAM-5.*A*AJ-24.*AI*
     1     LAM-6.*AI*AJ-24.*LAM*AJ))/(9.*(4.*A+3.*AI+3.*AJ))
      VSC3=(8.*(3.*A**2+A*AI+16.*A*LAM+A*AJ+12.*AI*LAM+3.*AI
     1     *AJ+12.*LAM*AJ))/(9.*(4.*A+3.*AI+3.*AJ))
C
      VSD0=-(32768.*SQRT(PI*2.)*PI**4*(A**2+A*AI+A*AJ+AI*AJ))/
     1     (81.*SQRT((A**2+A*AI+2.*A*LAM+A*AJ+AI*LAM+AI*AJ+LAM*AJ)*A)
     2     *(A**4+2.*A**3*AI+4.*A**3*LAM+2.*A**3*AJ+A**2*AI**2+6.*A
     3     **2*AI*LAM+4.*A**2*AI*AJ+4.*A**2*LAM**2+6.*A**2*LAM*AJ+A
     4     **2*AJ**2+2.*A*AI**2*LAM+2.*A*AI**2*AJ+4.*A*AI*LAM**2+8.*
     5     A*AI*LAM*AJ+2.*A*AI*AJ**2+4.*A*LAM**2*AJ+2.*A*LAM*AJ**2+
     6     AI**2*LAM**2+2.*AI**2*LAM*AJ+AI**2*AJ**2+2.*AI*LAM**2*AJ+
     7     2.*AI*LAM*AJ**2+LAM**2*AJ**2))
      VSD1=(2.*(-5.*A**3-5.*A**2*AI-12.*A**2*LAM-5.*A**2*AJ-
     1     7.*A*AI*LAM-5.*A*AI*AJ-7.*A*LAM*AJ-2.*AI*LAM*AJ))/(9.*(A**
     2     2+A*AI+2.*A*LAM+A*AJ+AI*LAM+AI*AJ+LAM*AJ))
      VSD2=(2.*(-5.*A**3-5.*A**2*AI-12.*A**2*LAM-5.*A**2*AJ-
     1     7.*A*AI*LAM-5.*A*AI*AJ-7.*A*LAM*AJ-2.*AI*LAM*AJ))/(9.*(A**
     2     2+A*AI+2.*A*LAM+A*AJ+AI*LAM+AI*AJ+LAM*AJ))
      VSD3=(8.*(-2.*A**3-2.*A**2*AI-3.*A**2*LAM-2.*A**2*AJ-A
     1     *AI*LAM-2.*A*AI*AJ-A*LAM*AJ+AI*LAM*AJ))/(9.*(A**2+A*AI+2.
     2     *A*LAM+A*AJ+AI*LAM+AI*AJ+LAM*AJ))
C
      VSB0=-4.D0*VSB0
      VSC0=-4.D0*VSC0
      VSD0= 8.D0*VSD0
      return
      end
C
      DOUBLE PRECISION FUNCTION DERF(X)
      real*8 s15aef,x
      integer ifail
      IFAIL=0
      DERF=s15aef(x,ifail)
      RETURN
      END
