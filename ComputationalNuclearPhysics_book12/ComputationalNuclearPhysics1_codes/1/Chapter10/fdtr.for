      PROGRAM FDTR
C
C========================================== fitting  version  ========C
C                                                                     C
C    FDU0 Package                                                     C
C    ---- -------                                                     C
C                                                                     C
C    Transitions and Moments                                          C
C    -----------     -------
C                                                                     C
C=====================================================================C
C
      PARAMETER (MAXNS=4200)              !maximun dimension for any j
      PARAMETER (MAXP=11221)
      NAMELIST /INP/NEUTRON_SYMMETRY,
     &              PROTON_SYMMETRY,
     &              OMP,OMN,              !basis parameter,shell sizes
     &
     &              N1P,N1N,              !pair numbers
     &
     &              EN,B2N,B3N,G0N,G2N,   !hamiltonian parameters
     &              EP,B2P,B3P,G0P,G2P,
     &              ENP,B2NP,B3NP,
     &
     &              KEY,                  !input control
     &                                    ! key.eq.0, the above two sets of
     &                                    !  parameters will not be used
     &                                    !  their values will be read
     &                                    !  from fdsm.wav.
     &                                    ! key.ne.0, input match
     &              EFFN,EFFP,GRN,GRP     !effective charges and g-factor
      INTEGER   JNFROM(MAXNS),ANFROM(MAXNS)
      INTEGER   JPFROM(MAXNS),APFROM(MAXNS)
      INTEGER   JNTO(MAXNS),ANTO(MAXNS)
      INTEGER   JPTO(MAXNS),APTO(MAXNS)
      INTEGER   OMN,OMP,OMN2,OMP2,AN,AN2,AP,AP2,A,R
      REAL      ZFROM(MAXNS),ZTO(MAXNS),RESULT
      CHARACTER FNAME*20,LINE*80,NEUTRON_SYMMETRY*3,PROTON_SYMMETRY*3
      CHARACTER NEUTRON_SYMMETRY2*3,PROTON_SYMMETRY2*3
      LOGICAL   EXIST,READWAVE,LTRIANGLE
      REAL*8    PTABN(MAXP),PTABP(MAXP),SIXJ
      INTEGER   PENDN,PENDP,PINDEXN(MAXP),PINDEXP(MAXP)
C
      INDP(AP,JP,R,A,J)=20000000*J+200000*JP+10000*R+100*A+AP
      LTRIANGLE(J1,J2,J3)= (ABS(J1-J2).LE.J3).AND.(J1+J2.GE.J3)
      PHASE(I)=1-2*MOD(ABS(I),2)
C
      OPEN(11,FILE='FDTR.INP',STATUS='OLD')
C                                                 set defaults
      EFFN=1.0
      EFFP=1.0
      GRN=1.0
      GRP=1.0
      KEY=1
      NEUTRON_SYMMETRY='SP6'
      PROTON_SYMMETRY='SP8'
      OMP=10
      OMN=15
      N1P=0
      N1N=0
      EP=0.
      B2P=0.
      B3P=0.
      G0P=0.
      G2P=0.
      EN=0.
      B2N=0.
      B3N=0.
      G0N=0.
      G2N=0.
      B2NP=0.
      B3NP=0.
      ENP=0.
      READ(11,INP)
      READ(11,*)
      OPEN(3,FILE='FDSM.WAV',STATUS='OLD',FORM='UNFORMATTED')
      IF(KEY.EQ.0) THEN
       READ(3) NEUTRON_SYMMETRY,PROTON_SYMMETRY,OMN,OMP,N1N,N1P,
     &           EN,B2N,B3N,G0N,G2N,
     &           EP,B2P,B3P,G0P,G2P,
     &           ENP,B2NP,B3NP
      ELSE
       READ(3) NEUTRON_SYMMETRY2,PROTON_SYMMETRY2,OMN2,OMP2,N1N2,N1P2,
     &           EN2,B2N2,B3N2,G0N2,G2N2,
     &           EP2,B2P2,B3P2,G0P2,G2P2,
     &           ENP2,B2NP2,B3NP2
       EXIST=(NEUTRON_SYMMETRY.EQ.NEUTRON_SYMMETRY2)
     &  .AND.(PROTON_SYMMETRY.EQ.PROTON_SYMMETRY2)
     &  .AND.(OMN.EQ.OMN2).AND.(OMP.EQ.OMP2)
     &  .AND.(N1N.EQ.N1N2).AND.(N1P.EQ.N1P2)
     &  .AND.(EN.EQ.EN2).AND.(B2N.EQ.B2N2).AND.(B2N.EQ.B3N2)
     &  .AND.(G0N.EQ.G0N2).AND.(G2N.EQ.G2N2)
     &  .AND.(EP.EQ.EP2).AND.(B2P.EQ.B2P2).AND.(B2P.EQ.B3P2)
     &  .AND.(G0P.EQ.G0P2).AND.(G2P.EQ.G2P2)
     &  .AND.(ENP.EQ.ENP2).AND.(B2NP.EQ.B2NP2).AND.(B3NP.EQ.B3NP2)
       IF(.NOT.EXIST) THEN
        PRINT*,'****Input wavefunction error************'
        PRINT*,'Solution 1: If you have wavefunctions in file with'
        PRINT*,'     name other than FDSM.WAV, rename it to FDSM.WAV'
        PRINT*,'Solution 2: Modify FDSM.INP so that it agree with'
        PRINT*,'      FDSMTR.INP, and run FDSMU0'
        PRINT*,'Solution 3: Set key=0 in FDSMTR.INP'
        STOP '**********************************************'
       END IF
      END IF
      IF(N1N.NE.0) THEN
       CALL GETFILE(NEUTRON_SYMMETRY,'P',OMN,N1N,PENDN,MAXP,
     &                             PINDEXN,ANY,PTABN)
      ELSE
       PENDN=1
       PINDEXN(1)=INDP(1,0,0,1,0)
       PTABN(1)=0.D0
      END IF
      IF(N1P.NE.0) THEN
       CALL GETFILE(PROTON_SYMMETRY,'P',OMP,N1P,PENDP,MAXP,
     &                                   PINDEXP,ANY,PTABP)
      ELSE
       PENDP=1
       PINDEXP(1)=INDP(1,0,0,1,0)
       PTABP(1)=0.D0
      END IF
      OPEN(2,FILE='FDTR.OUT',STATUS='NEW')
      DO I=1,80
       LINE(I:I)='*'
      END DO
      WRITE(2,'(A,/,A,/,1X,2(A,A,4X),/,1X,4(A,I2,4X),//
     &,1X,5(A,F8.4,3X),/
     &1X,5(A,F8.4,3X),/,1X,3(A,F8.4,3X),//,2(A,F8.4,A,F8.4,A))')
     &' Input parameters:',
     &LINE,
     &'NEUTRON_SYMMETRY=',NEUTRON_SYMMETRY,
     &'PROTON_SYMMETRY=',PROTON_SYMMETRY,
     &'OMN=',OMN,'OMP=',OMP,'N1N=',N1N,'N1P=',N1P,
     &'EN=',EN,'B2N=',B2N,'B3N=',B3N,'G0N=',G0N,'G2N=',G2N,
     &'EP=',EP,'B2P=',B2P,'B3P=',B3P,'G0P=',G0P,'G2P=',G2P,
     &'ENP=',ENP,'B2NP=',B2NP,'B3NP=',B3NP,
     &' Q2=',EFFN,'*P2(N)+',EFFP,'*P2(P)    ',
     &' M1=',GRN,'*J(N)+',GRP,'*J(P)'
      WRITE(2,'(A)') LINE
      WRITE(2,'(A,A)') ' LM   FROM             TO',
     &'         SN        SP        RESULT'
      PRINT '(A,A)', ' LM   FROM             TO',
     &'         SN        SP        RESULT'
C                                                        loop begins
  10   READ(11,*,END=100) LM,FROM,TO
       J1=FROM
       J2=TO
       IF((LM.NE.1).AND.(LM.NE.2)) GOTO 10
       EXIST=READWAVE(FROM,NS1,JNFROM,ANFROM,JPFROM,APFROM,E1,ZFROM)
       IF(.NOT.EXIST) GOTO 10
       EXIST=READWAVE(TO,  NS2,JNTO,ANTO,JPTO,APTO,E2,ZTO)
       IF(.NOT.EXIST) GOTO 10
       IF(.NOT.LTRIANGLE(J1,J2,LM)) THEN
        SN=0.
        SP=0.
        RESULT=0.0
        CALL OUTPUT(LM,FROM,TO,E1,E2,SN,SP,RESULT)
        GOTO 10
       END IF
       IF(LM.EQ.1) THEN
C                                                 summation for dipole
        SN=0.
        SP=0.
        DO 3 I1=1,NS1
         IF(ABS(ZFROM(I1)).LT.1E-20) GOTO 3
         JN=JNFROM(I1)
         AN=ANFROM(I1)
         JP=JPFROM(I1)
         AP=APFROM(I1)
         DO 4 I2=1,NS2
          IF(ABS(ZTO(I2)).LT.1E-20) GOTO 4
          IF(JNTO(I2).NE.JN) GOTO 4
          IF(ANTO(I2).NE.AN) GOTO 4
          IF(JPTO(I2).NE.JP) GOTO 4
          IF(APTO(I2).NE.AP) GOTO 4
          SN=SN+ZFROM(I1)*ZTO(I2)*
     &      REAL( SIXJ(2*JP+1,2*JN+1,2*J1+1,
     &                 3     ,2*J2+1,2*JN+1) )*
     &      PHASE(JP+JN+1-J1)*
     &      SQRT( (2*J1+1)*(2*JN+1)*JN*(JN+1.0) )
C     &      UXSH(REAL(JP),REAL(JN),REAL(J2),1.0,REAL(J1),REAL(JN))*
C     &      PHASE(J2-J1)*SQRT(JN*(JN+1.0))
          SP=SP+ZFROM(I1)*ZTO(I2)*
     &      REAL( SIXJ(2*JN+1,2*JP+1,2*J1+1,
     &                 3     ,2*J2+1,2*JP+1) )*
     &      PHASE(JN+JP+J2+1)*
     &      SQRT( (2*J1+1)*(2*JP+1)*JP*(JP+1.0) )
C     &      UXSH(REAL(JN),REAL(JP),REAL(J2),1.0,REAL(J1),REAL(JP))*
C     &      SQRT(JP*(JP+1.0))
  4      END DO  !I2
  3     END DO !I1
        S=GRN*SN+GRP*SP
        IF(FROM.EQ.TO) THEN
         RESULT=SQRT(4*3.1415926/3.0)*S*
     &   SQRT(J1/(J1+1.))
        ELSE
         RESULT=(2*J2+1)/(2.0*J1+1)*S*S
        END IF
        CALL OUTPUT(LM,FROM,TO,E1,E2,SN,SP,RESULT)
       ELSE                   !LM=2
C                                                 summation for quadrupole
        SN=0.
        SP=0.
        DO 13 I1=1,NS1
         IF(ABS(ZFROM(I1)).LT.1E-20) GOTO 13
         JN=JNFROM(I1)
         AN=ANFROM(I1)
         JP=JPFROM(I1)
         AP=APFROM(I1)
         DO 14 I2=1,NS2
          IF(ABS(ZTO(I2)).LT.1E-20) GOTO 14
          JN2=JNTO(I2)
          AN2=ANTO(I2)
          JP2=JPTO(I2)
          AP2=APTO(I2)
          IF( (JP.EQ.JP2).AND.(AP.EQ.AP2) )THEN
           IF(.NOT.LTRIANGLE(JN,JN2,2)) THEN
            P2N=0.
           ELSE
            INDEX=INDP(AN2,JN2,2,AN,JN)
            P2N=PTABN(ISECH(PINDEXN,PENDN,INDEX,1))
           END IF
           SN=SN+ZFROM(I1)*ZTO(I2)*
     &      REAL(SIXJ(2*JP+1,2*JN+1,2*J1+1,
     &                5     ,2*J2+1,2*JN2+1) )*
     &      PHASE(-JN2-J1+JP)*
     &      SQRT( (2*J1+1)*(2*JN2+1.0) )*
     &      P2N
C     &      UXSH(REAL(JP),REAL(JN),REAL(J2),2.0,REAL(J1),REAL(JN2))*
C     &      PHASE(JN-JN2+J2-J1)*P2N
          END IF
          IF( (JN.EQ.JN2).AND.(AN.EQ.AN2) )THEN
           IF(.NOT.LTRIANGLE(JP,JP2,2)) THEN
            P2P=0.
           ELSE
            INDEX=INDP(AP2,JP2,2,AP,JP)
            P2P=PTABP(ISECH(PINDEXP,PENDP,INDEX,1))
           END IF
           SP=SP+ZFROM(I1)*ZTO(I2)*
     &      REAL( SIXJ(2*JN+1,2*JP+1,2*J1+1,
     &                 5,     2*J2+1,2*JP2+1) )*
     &      PHASE(JN+JP+J2+2)*
     &      SQRT( (2*J1+1)*(2*JP2+1.0) )*
     &      P2P
C     &      UXSH(REAL(JN),REAL(JP),REAL(J2),2.0,REAL(J1),REAL(JP2))*
C     &      P2P
          END IF
 14      END DO  !I2
 13     END DO !I1
        S=EFFN*SN+EFFP*SP
        IF(FROM.EQ.TO) THEN
         RESULT=SQRT(16*3.1415926/5.0)*S*
     &   SQRT( J1*(2*J1-1.)/(J1+1.)/(2*J1+3.) )
        ELSE
         RESULT=(2*J2+1)/(2.0*J1+1)*S*S
        END IF
        CALL OUTPUT(LM,FROM,TO,E1,E2,SN,SP,RESULT)
       END IF
      GOTO 10
 100  END

      SUBROUTINE OUTPUT(LM,FROM,TO,E1,E2,SN,SP,RESULT)
C=====================================================================C
C     Output the results                                              C
C=====================================================================C
      CHARACTER C*15
C
      IF(LM.EQ.1) THEN
       IF(FROM.EQ.TO) THEN
        C='(dipole)'
       ELSE
        C='(BM1)'
       END IF
      ELSE IF(LM.EQ.2) THEN
       IF(FROM.EQ.TO) THEN
        C='(qardrupole)'
       ELSE
        C='(BE2)'
       END IF
      ELSE
       C='(not valid)'
      END IF
      PRINT 10,LM,FROM,E1,TO,E2,SN,SP,RESULT,C
      WRITE(2,10) LM,FROM,E1,TO,E2,SN,SP,RESULT,C
 10   FORMAT(1X,I1,2(1X,F4.1,'(',F7.3,')'),3(1X,F8.4),A)
      RETURN
      END

      LOGICAL FUNCTION
     &      READWAVE(STATE,NS,JN,AN,JP,AP,E,Z)
C=====================================================================C
C     Reads the wavefunctions                                         C
C=====================================================================C
      REAL     STATE,E,Z(1)
      INTEGER  NS,JN(1),AN(1),JP(1),AP(1)
C
      READWAVE=.TRUE.
      J=STATE
      N=NINT(10*(STATE-J))
      J_KEEP=-100
      READ(3,END=10) J_READ,NS,NJ
      GOTO 20
 10   REWIND(3)
      READ(3)
  1   READ(3,END=30) J_READ,NS,NJ
 20   IF(J_READ.LT.J) THEN
       READ(3,END=30)
       READ(3,END=30)
       GOTO 1
      ELSE IF(J_READ.GT.J) THEN
       IF(J_KEEP.EQ.J_READ) THEN
        BACKSPACE(3,ERR=30)
        GOTO 30
       END IF
       J_KEEP=J_READ
       DO I=1,4
        BACKSPACE(3,ERR=30)
       END DO
       GOTO 1
      ELSE
       IF(N.GT.NJ) THEN
        BACKSPACE(3,ERR=30)
        GOTO 30
       END IF
       READ(3,END=30) (JN(I),AN(I),JP(I),AP(I),I=1,NS)
       READ(3,END=30) ( E,(Z(I),I=1,NS) ,II=1,N)
       RETURN
      END IF
 30   READWAVE=.FALSE.
      RETURN
      END
