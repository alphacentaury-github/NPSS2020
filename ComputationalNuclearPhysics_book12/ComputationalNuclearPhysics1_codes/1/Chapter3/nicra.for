C***********************************************************************
C****  main        N I L S S O N     C R A N K E R             *********
C***********************************************************************
C*********************Mathematical Physics, Lund 1989*******************
C***********************************************************************
C ********** parameters controlling the dimension of matrices:
C
C MAXN= max oscillator shell in calcs. MAXDIM= the resulting dimension
C MAXN    4   5   6   7   8   9   10  11  12
C MAXDIM  22 34  50  70  95  125  161 203 252  for MSYM=2
C MSYM=# of symmetry groups for one kind of particle, i.e. parity
C MOMEGA = max # of omegas (omega = rot. frequency)
C MISO = max no. of particle types= 2! (prot & neu)
C MDEF = max. no. of deformations for one deformation parameter
C MERR = max. no. of errors to be stored in ERR_LOG
C
C  *** contents of communication common blocks:
C     (other common blocks are for memory reduction only,
C      values can be regarded as unimportant outside the calling
C      sequence of the subroutines in which they appear)
C
C  /HAMPAR/ parameters for the present diagonalization
C     i.e. deformation, frequency, delta, lambda, dimension & isospin
C     note:GAM is in radians and gamma+120
C         isospin ISO =1 for protons, 2 for neutrons
C     NUMBER(ISO) is particle number for protons and neutrons
C     ISO = present isospin  IANT = present storing point
C     MOVE = present isospin storing offset, NEUMOV = offset storing
C         point for neutrons
C
C  /NILSPAR/ parameters for Nilsson Hamiltonian
C     kappa's & my's for 0-MAXN shells and proton, neutrons
C     RHOMEGA(ISO) = oscillator strength, without volume conservation
C     OMOM as usual the volume conservation term
C     NUU max no. of shells coupled, if <0 => JCO=0 i.e. no coupling
C         between J-shells
C     NMAX(ISO) max N-shell used
C
C  /HAMILTON/ matrices for the construction of Hamiltonian matrix
C     EIJ single particle energy matrix
C     RJXIJ jx matrix    Q2IJ Q2 matrix
C
C  /QUANTA/ quantum number vectors for the matrix elements in /HAMILTON/
C     NVAL main shell quanta, LVAL orbital angular momentum
C     JVAL =2*angular momentum ,MVAL = 2* projection of angular momentum
C         on cranking axis
C
C  /ISOS/ calculated total quantities for protons and neutrons
C      at different omega
C     ICONFI = signature and parity for optimal state = 100*is+ip
C         where is=2*alpha(=0, 1, 2, 3), ip= 1(-) or 0(+)
C     VAR = vector keeping the main variable as selected by NVAR
C         usually omega
C
C  /SETTINGS/ various settings for the run:
C     IF_SPOUT: sp-levels on file for plotting if 1
C     IF_STRUT: 0 = no strutinsky corr OR liquid drop calc
C               1 = stutinsky corr for omega=0
C               2 = 1 + renormalization of moment of inertia
C     LEV_PRINT: contolls size of printout, the larger the more
C         = 20 (21) total (proton&neutron) quantities outside omega-loop
C         > 50 quantities inside omega loop
C         output from specific subroutines is obtained if
C         LEV_PRINT = routine #
C     IF_LEVSAV:
C         if +1/-1 saves/reads(if possible) sp-levels from savefiles
C         if 0 no savefiles used - all calculations from scratch!
C     IF_ISOENEROUT: 1 gives printing on file of proton and neutron
C                      total quantities (LEV_PRINT .LE. 21)
C     IF_TOTEOUT: 1 gives printing on file of nuclear total quant
C                    for calc frequencies
C     MODELTYP: =1 (modified oscillator)
C     ISYM1,ISYM2: =1,2 parities '+' and '-'
C     ISO1,ISO2: =1,2 if both protons(=1) and neutrons(=2) studied
C         1,1 => only protons,  2,2 => only neutrons
C     NVAR: # of main variable for calculation (used for output)
C         NVAR=1 => EPS, 2=>GAM, 3=>OMEG, 4=>EPS4,  >4 unspec
C         should be 3 if rotation is considered!
C     IF_MEV: =1 if /HAMILTON/ matrices are calculated in MeV
C     SPEMIN,SPEMAX: max & min sp energy for printout
C     IF_Q2 : =1 if electric quadrupole moment wanted
C             =2 if total nuclear density quadrupole moment
C     IF_NIVPRI= 0 no orbital values printed, only summed quantities
C              = 1 also orbital values printed
C              = 2 no summed quantities calculated (or printed!)
C     IF_DROPCALC = 1 if liquid drop parameters calculated
C           ONLY if IF_STRUT>0
C     IF_BCS = 1 if BCS calculation performed for omega=0 and BCS energy
C         added to total energies at omega=I=0.
C     IN_LEV input complexity level selector, see description of input
C     ISTEP > 0 gives printing for interpolated spin values in steps
C                     of ISTEP
C           < 0 total routhians, R(w)=E(w)-wI, as output
C
C  /FILES/ output file numbers used
C  /SAVFIL/ various parameters initially read from savefiles, for check
C     IPSAV=max no of (eps,gamma,eps4,omega) points expected to be
C     stored in SAVFILE system
C
C  /ORBITALS/ single-particle energies E, <jx> RJX, quadrupole moment Q2
C     and wave functions VEC
C     signature +1/2 saved in 1 to NDIM, -1/2 in  NDIM+1 to 2*NDIM
C     for both isospins if NEUMOV > 0
C
C  /STRUTI/ strutinsky, liquid drop and BCS results, parameters for
C         BCS calculation G0, G1 and KOEFF_RANGE (orbitals within
C         sqrt(KOEFF_RANGE*NUMBER(ISO)) from fermi energy included in
C         BCS calc.
C
C  /LOOPS/ LOOP1 parameters:
C     IDEF_TYP = -1 x(=eps) & y(=gamma) is input as a rectangular
C                  cartesian mesh. loop1 is over eps4,y,x in that order
C                                (i.e. inner loop over eps4)
C              = 0 standard input, loop1 is over eps4,gam,eps
C              = 1 eps and gam vary simultaneously, loop1: eps4,(eps,gam)
C              = 2 eps,gam and eps4 vary simultaneously
C              = 3 eps,gam,eps4 and omega vary simultaneously
C                  for IDEF_TYP < 3, omega varies faster than def.
C     FI = tilt of cartesian mesh (see FUN2 and XYCOORD)
C     NEPS,  EPSI(MDEF) epsilon (or x if IDEF_TYP=-1) deformations
C     NGAM,  GAMI(MDEF) gamma (or y)
C     NEPS4, EPS4I(MDEF) epsilon 4, hexadecapole
C     NOMEGA, OMEGA(MOMEGA) cranking frequencies
C     IE,IG,I4 = the deformation mesh points
C
C  /TOTALS/ accumulated total quantities
C
C  /ERR_LOG/ just that!
C
C  /NUC_TEXT/ character variables
C
C Note: variable names may be truncated to 6 characters and underscores
C     omitted
*-------------------------------------
C     input structure:
C     note: each input unit has a unique number first.
C     this number is for error detection and convenience
C     as a specific READ statement goes into a loop until it has read
C     a unit number higher or equal to that it is supposed to read!
C  1,NUMBER(1),NUMBER(2),LEV_PRINT,IN_LEV       !Z,N,write read level
C  2,NEPS,(EPSI(I),I=1,NEPS)      !* alternative input form
C  3,NGAM,(GAMI(I),I=1,NGAM)      !*    for defs & omega:
C  4,NEPS4,(EPS4I(I),I=1,NEPS4)   !* ###,-3,start,stop,step
C  5,NOMEGA,(OMEGA(I),I=1,NOMEGA)
*     if(IN_LEV.le.0) ready with reading, read 1000
C  6,IDEF_TYP,NVAR,FI
C  7,'filename',IF_SPOUT,SPERANGE
C  8,'filename',IF_LEVSAV
C  9,'filename',IF_TOTEOUT,ISTEP
*     if(IN_LEV.le.1) ready with reading, read 1000
C  10,IF_MEV,IF_STRUT,IF_Q2,IF_BCS
C  11,IF_DROPCALC,ISO1,ISO2,IF_ISOENEROUT
*     if(IN_LEV.le.2) ready with reading, read 1000
C  12,NUU,NMAX(pro),NMAX(neu),G0,G1,KOEFF_RANGE
*  do shell=0,11
C     13+shell,(RKAPPA(shell,iso),RMY(shell,iso),iso=pro,neu)
*  enddo
C 1000,IF_GO,(0=>stop)
C
C  end of input
C NOTE: three files used if IF_LEVSAV.ne.0
C       1. 'filename' for information data
C       2. 'EN-filename(1:13)' for E & JX
C       3. 'LQ-filename(1:13)' for liquid drop parrameters
*
      CHARACTER SPFILE*24,SAVFILE*24,ENERFILE*24
      COMMON/TEST/ CCCC20,CCCC22
* VAX FORTRAN include statements!
      INCLUDE 'NICRAINC.FOR/LIST'
      CALL CHECK_READ(1,0)
C  read particle numbers Z & N
1     READ(*,*)INI,NUMBER(1),NUMBER(2),LEV_PRINT,IN_LEV
      CALL CHECK_READ(INI,1,*1)
      CALL LOOP_READER(NDEFS)
      IF(IN_LEV.GT.0)THEN
7         READ(*,*)INI,SPFILE,IF_SPOUT,SPERANGE
          CALL CHECK_READ(INI,7,*7)
          IF(IF_SPOUT.EQ.1)THEN
             OPEN(NSPFILE,FILE=SPFILE,ACCESS='SEQUENTIAL',
     &       STATUS='UNKNOWN')
          ENDIF
8         READ(*,*)INI,SAVFILE,IF_LEVSAV
          CALL CHECK_READ(INI,8,*8)
C  file opened in NILS_READER!
9         READ(*,*)INI,ENERFILE,IF_TOTEOUT,ISTEP
          CALL CHECK_READ(INI,9,*9)
      ENDIF
      IF(IN_LEV.GT.1)THEN
10        READ(*,*)INI,IF_MEV,IF_STRUT,IF_Q2,IF_BCS
          CALL CHECK_READ(INI,10,*10)
11        READ(*,*)INI,IF_DROPCALC,ISO1,ISO2,IF_ISOENEROUT
          CALL CHECK_READ(INI,11,*11)
      ENDIF
      IF(IF_TOTEOUT+IF_ISOENEROUT.GE.1)THEN
            OPEN(NTOTFILE,FILE=ENERFILE,ACCESS='SEQUENTIAL',
     &      STATUS='UNKNOWN')
      ENDIF
* find out if offset for s.p. values storing is needed & possible:
      NEUMOV=0
      IF(NVAR.NE.3.AND.ISO1.NE.ISO2)THEN
          IF(2*NDEFS.LE.MOMEGA)NEUMOV=NDEFS
      ENDIF
      IF(MODELTYP.EQ.1)CALL NILS_READER(SAVFILE)
      CALL EXPLAIN_RUN(SPFILE,SAVFILE,ENERFILE)
1000  READ(*,*)INI,IF_GO
      CALL CHECK_READ(INI,1000,*1000)
      IF(IF_GO.EQ.0)THEN
          WRITE(*,*)' OH, you only wanted to test-read input.'
          WRITE(*,*)'  Then I''ll stop now'
          GO  TO 999
      ENDIF
      CALL CHECK_READ(1,-1,*1001)
*
* no more reading after this point!
*
1001  IANT=1
      TID1=SECOND(U)
      DO 100 ILOOP1=1,NDEFS
          ILOOP=ILOOP1
          TID2=SECOND(U)
          CALL LOOP1_SET
          CALL DEFSET
          IF((LEV_PRINT.GE.20.AND.NVAR.EQ.3).OR.LEV_PRINT.GE.50)THEN
               IF(ILOOP1.GT.1)WRITE(*,'(1H1)')
               WRITE(*,1958)ILOOP1
1958           FORMAT(35('-'),'ILOOP1=',I3,35('-'))
               WRITE(*,1957)NUMBER(1)+NUMBER(2),ANUC(NUMBER(1))
1957           FORMAT(24X,I4,A,' calculations with:')
               WRITE(*,1955)EPS,57.29*GAM-120,EPS4,IE,IG,I4
1955           FORMAT(' deformation:  E P S =',F7.3,'  G A M M A =',
     &         F7.1,'   E P S 4 =',
     &         F7.3/20X,' mesh point (Ieps,Igam,Ieps4)=(',3(i2,','),')')
               WRITE(*,1956)OMOM,DROPENERGY,DROPINERTIA
1956           FORMAT(12X,' w0/w00=',f8.4,' L.D. energy=',f7.3,
     &     ' inertia=',f8.3)
          ENDIF
          DO 103 ISO=ISO1,ISO2
              IF(ISO.EQ.ISO1)THEN
                  MOVE=0
              ELSE
                  MOVE=NEUMOV
              ENDIF
              IF(MODELTYP.EQ.1)CALL SETMAT_NILS
              CALL ORBITS
103       CONTINUE
          IF_OUT=1
          IF(NVAR.EQ.4.AND.I4.LT.NEPS4)IF_OUT=0
          IF(NVAR.EQ.2.AND.IG.LT.NGAM)IF_OUT=0
          IF(NVAR.EQ.1.AND.IE.LT.NEPS)IF_OUT=0
          IF(NVAR.GT.4)IF_OUT=0
          IF(IF_OUT.EQ.1.AND.LEV_PRINT.GE.21.AND.NVAR.NE.3)
     &     CALL PRINIV(IANT)
          IF(ISO2-ISO1.GT.0.AND.IF_NIVPRI.LE.1)THEN
              CALL SUM_CONFS
              IF(IF_STRUT.GE.2.AND.LEV_PRINT.GE.20)
     &        WRITE(*,1954)STRUTINERTIA,DROPINERTIA
1954          FORMAT(' strutinsky inertia=',F8.2,' L.D. inertia=',
     &        F8.2)
              IF(IF_OUT.GT.0)THEN
                  WRITE(*,'(1H1)')
                  WRITE(*,'(80(''=''))')
                  WRITE(*,*)'the final results:'
                  CALL OUT_WITHIT
                  IANT=0
              ENDIF
          ENDIF
          TID3=SECOND(U)
          IF((LEV_PRINT.GE.20.AND.NVAR.EQ.3).OR.LEV_PRINT.GE.50)THEN
              TID=TID3-TID2
              WRITE(*,'(A,F10.1,A)')' time for deformation:',TID,'s'
          ENDIF
          NANT=IANT
          IANT=IANT+1
          IF(IANT.GT.MOMEGA)THEN
              WRITE(ERR_TEX(NERR+1),1800)ILOOP1
1800          FORMAT(' LOOP1 error: IANT>MOMEGA at ILOOP1=',I3)
              CALL ERR_OUT(1)
          ENDIF
100   CONTINUE
      IF(NVAR.GT.4)THEN
          IF(LEV_PRINT.GE.21)CALL PRINIV(NANT)
          IF(IF_NIVPRI.LE.1)THEN
              WRITE(*,'(1H1)')
              WRITE(*,'(80(''=''))')
              WRITE(*,*)'the final results:'
              CALL OUT_WITHIT
          ENDIF
      ENDIF
999   TID3=SECOND(U)
      TID=(TID3-TID1)/60.
      CALL ERR_OUT(0)
      WRITE(*,'(A,F9.2,A)')' total time:',TID,' minutes'
      WRITE(*,*)' thanks for now'
      STOP
      END
      SUBROUTINE EXPLAIN_RUN(SPFILE,SAVFILE,ENERFILE)
*********************************************************************
* utility routine in NILSSON-CRANKER                                *
* purpose: printing and checking the parameters for the present run *
* called by: MAIN                                                   *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      CHARACTER SPFILE*24,SAVFILE*24,ENERFILE*24
*  first: checking
      MOD_INPUT=0
      IF(IF_MEV.EQ.0)THEN
          IF(IF_NIVPRI.LT.2)THEN
              WRITE(*,*)' with orbitals calculated in hw0, I cannot'
              WRITE(*,*)' do a decent energy calculation, will do none'
C              WRITE(*,*)' IF_NIVPRI reset to 2'
              IF_NIVPRI=2
              MOD_INPUT=MOD_INPUT+1
          ENDIF
      ENDIF
      IF(IF_NIVPRI.EQ.2)THEN
          IF_ISOENEROUT=0
          IF_TOTEOUT=0
      ENDIF
* check IDEF_TYP/variable consistency
      IF(IDEF_TYP.EQ.-1.AND.(NVAR.LE.2.OR.NVAR.GT.4))THEN
          WRITE(*,*)' for a cartesian grid the variable could not be'
          WRITE(*,*)' eps or gamma or undefined, variable is put omega'
          NVAR=3
          MOD_INPUT=MOD_INPUT+1
      ENDIF
      IF(NVAR.LE.2)THEN
          IF(NEPS4.GT.1.AND.IDEF_TYP.LT.2)THEN
              IF(NOMEGA.GT.1)THEN
                  WRITE(*,*)' when main variable is eps or gam, you',
     &            ' cannot have both eps4 and omega loops',
     &            ', no eps4 loop'
                  NEPS4=1
                  MOD_INPUT=MOD_INPUT+1
              ELSE
                  WRITE(*,*)' when main variable is eps or gam, you',
     &            ' cannot have eps4 loop, no eps4 loop'
                  NEPS4=1
                  MOD_INPUT=MOD_INPUT+1
              ENDIF
         ELSEIF(NOMEGA.GT.1.AND.IDEF_TYP.LT.3)THEN
             WRITE(*,*)' if Nomega >1 variable should be omega'
             NVAR=3
             MOD_INPUT=MOD_INPUT+1
         ENDIF
      ENDIF
* check consistency of ISTEP variable
      IF (ISTEP .GE. 1) THEN
         IF (NVAR .NE. 3) THEN
            WRITE(*,*) ' no interpolation to fixed spins if omega is'
            WRITE(*,*) ' not main variable,  ISTEP reset to 0'
            ISTEP = 0
            MOD_INPUT=MOD_INPUT+1
         ENDIF
         IF (IDEF_TYP .GT. 2) THEN
            WRITE(*,*) ' no interpolation to fixed spins if omega and'
            WRITE(*,*) ' def vary simultaneously; ISTEP reset to 0'
            ISTEP = 0
            MOD_INPUT=MOD_INPUT+1
         ENDIF
      ENDIF
      IF(IF_STRUT.GE.2.AND.(NVAR.NE.3.OR.NOMEGA.EQ.1))THEN
          IF_STRUT=1
          WRITE(*,*)' IF_STRUT reset to 1, no Strut. renorm.of inertia'
          WRITE(*,*)' when variable not omega or only one omega'
      ENDIF
      IF(IF_STRUT.GE.2.AND.IF_DROPCALC.LT.1)THEN
          IF_DROPCALC=1
          WRITE(*,*)' A liquid drop calculation is neccesary to do'
          WRITE(*,*)' a Struinsky renormalisation of inertia'
          MOD_INPUT=MOD_INPUT+1
      ENDIF
      IF(IF_LEVSAV.LT.0.AND.IF_Q2.GT.0)THEN
          IF_Q2=0
          WRITE(*,*)' As wavefunctions are not available from'
          WRITE(*,*)' saver files, Q2''s cannot be calculated'
          MOD_INPUT=MOD_INPUT+1
      ENDIF
      IF(MOD_INPUT.NE.0)THEN
          WRITE(ERR_TEX(NERR+1),1800)MOD_INPUT
1800      FORMAT(I4,' input inconsistencies modified')
          CALL ERR_OUT(1)
      ENDIF
*  then : printing
      WRITE(*,'(1H1/1X,80(''*''))')
      WRITE(*,1)
1     FORMAT(' ***',7X,'     about this run of         ',36X,'***'/
     &       ' ***',7X,'N I L S S O N     C R A N K E R',36X,'***')
      WRITE(*,'(1X,50(''*''),A,''******'')')'Math Phys Lund 1989*****'
      WRITE(*,*)
      WRITE(*,2)NUMBER(1)+NUMBER(2),ANUC(NUMBER(1)),NUMBER(1),NUMBER(2)
2     FORMAT(' The nucleus calculated for is',I4,A,' with Z=',I3,
     &' N=',I4)
      IF(ISO1.EQ.ISO2)WRITE(*,*)' only ',ISOTEXT(ISO1),' calculated'
*
      WRITE(*,'(/1X,35(''-''),''potential'',36(''-''))')
      ENT=' MeV'
      IF(MODELTYP.EQ.1)THEN
          WRITE(*,*)' Nilsson (modified oscillator) potential used'
          WRITE(*,*)'    with the parameters:'
          WRITE(*,3)(ITTE-1,RKAPPA(MIN(ITTE,15),1),RMY(MIN(ITTE,15),1),
     &    RKAPPA(MIN(ITTE,15),2),RMY(MIN(ITTE,15),2),
     &    ITTE=1,MAX(NMAX(1),NMAX(2))+1)
 3        FORMAT(/'  N    kappaP  myP     kappaN  myN'/13(I3,4F8.4/))
          WRITE(*,4)(ISOTEXT(IS),NMAX(IS),(NDIM(IY,IS),IY=1,2),
     &     IS=ISO1,ISO2)
4         FORMAT(' max shell for ',A,' is',I2,' =>',
     &    ' dimension for parity - is:',I4,', parity + is:',I4)
          WRITE(*,5)NUU,RHOMEGA(1),RHOMEGA(2)
5         FORMAT(I5,' shells are coupled, hw00(P)=',F8.4,
     &    'MeV, hw00(N)=',F8.4,'MeV')
          IF(IF_MEV.EQ.0)WRITE(*,*)' orbital energies in hw0'
          IF(IF_MEV.EQ.0)ENT=' hw0'
      ENDIF
      IF(IF_MEV.EQ.1)WRITE(*,*)' energies given in MeV'
      IF(IF_STRUT.EQ.0)WRITE(*,*)' no strutinsky or L.D. calculations'
      IF(IF_BCS.GT.0)THEN
          WRITE(*,*)' BCS energy at omega=0 added'
          WRITE(*,*)' the pairing force constant G=(g0+-g1(N-Z)/A)/A '
          WRITE(*,'(A,F6.2,A,F6.2,A)')
     &     ' with g0=',G0,'MeV, g1=',G1,'MeV.'
          WRITE(*,'(A,I2,A)')' Orbitals within sqrt(',KOEFF_RANGE,
     &     '*N or Z) from fermi level, included in BCS calc.'
      ENDIF
      IF(IF_STRUT.GE.1)THEN
          WRITE(*,*)' strutinsky correction at omega=0'
          IF(IF_DROPCALC.EQ.1)WRITE(*,*)' and liquid drop energy added'
          IF(IF_STRUT.GE.2) WRITE(*,*)' and Strut. renormalization',
     &    ' of moment of inertia'
      ENDIF
*
      WRITE(*,'(/1X,34(''+''),''deformations'',34(''+''))')
      IF(IDEF_TYP.EQ.-1)THEN
          WRITE(*,600)FI
600       FORMAT(' eps,gamma deformations vary in a cartesian grid'/
     & ' gamma=arctan(y/x)+',F5.0,', eps=Sqrt(x**2+y**2)')
          WRITE(*,601)' x= ',(EPSI(I),I=1,NEPS)
          WRITE(*,601)' y= ',(GAMI(I),I=1,NGAM)
601       FORMAT(A,20F8.3)
604       FORMAT(A,6(F8.3,A)/(7X,6(F8.3,A)))
      ELSE
          WRITE(*,*)' the deformations used are:'
          WRITE(*,601)' eps= ',(EPSI(I),I=1,NEPS)
          WRITE(*,601)' gam= ',(GAMI(I),I=1,NGAM)
      ENDIF
      WRITE(*,601)' eps4=',(EPS4I(I),I=1,NEPS4)
      WRITE(*,*)' cranking frequencies used are:'
      IF(IDEF_TYP.NE.3)THEN
          I1=1
          NOM=NOMEGA
      ELSE
          I1=2
          NOM=2
          DO 700 I=MOMEGA,2,-1
              IF(OMEGA(I).LT.-10)NOM=I-1
700       CONTINUE
      ENDIF
      WRITE(*,604)' omega= ',(OMEGA(I),ENT,I=I1,NOM)
      IF(IDEF_TYP.EQ.1)WRITE(*,*)' eps and gamma vary',
     &' simultaneously'
      IF(IDEF_TYP.EQ.2)WRITE(*,*)' eps, gamma and eps4 vary',
     &' simultaneously'
      IF(IDEF_TYP.EQ.3)WRITE(*,*)' eps, gamma,eps4 and omega vary',
     &' simultaneously'
      IF(NVAR.EQ.1)WRITE(*,*)' main output variable is eps'
      IF(NVAR.EQ.2)WRITE(*,*)' main output variable is gamma'
      IF(NVAR.EQ.3)WRITE(*,*)' main output variable is omega'
      IF(NVAR.EQ.4)WRITE(*,*)' main output variable is eps4'
      IF(NVAR.GT.4)WRITE(*,*)' all parameter variations traced'
      WRITE(*,*)
      IF(IF_Q2.GT.0)THEN
          IF(IF_MEV.GT.0)THEN
              WRITE(*,*)' <Q20> calculated in fm**2'
          ELSE
              WRITE(*,*)' <Q20> calculated in h**2/Mw0'
          ENDIF
          IF(IF_Q2.EQ.1)THEN
              WRITE(*,*)' electric <Q20> from protons calculated'
          ELSEIF(IF_Q2.EQ.2)THEN
              WRITE(*,*)' total nuclear density <Q20> calculated'
          ENDIF
      ELSE
          WRITE(*,*)' <Q20> not calculated'
      ENDIF
*
      WRITE(*,'(/1X,36(''+''),''printout'',36(''+''))')
      WRITE(*,6)LEV_PRINT
6     FORMAT('  printout level is',I3,' i.e.')
      IF(LEV_PRINT.GT.54.OR.LEV_PRINT.GT.21.AND.LEV_PRINT.LT.50)
     & WRITE(*,*)' debugging printout'
      IF(LEV_PRINT.GE.50)WRITE(*,*)' printout from innermost parameter',
     &' loop'
      IF(LEV_PRINT.GE.21)WRITE(*,*)
     &' proton and neutron printout from deformation loop'
      IF(LEV_PRINT.GE.20)WRITE(*,*)' printout from deformation loop'
      WRITE(*,*)
      WRITE(*,*)' file output:'
      IF(IF_SPOUT.EQ.1.AND.LEV_PRINT.GE.21)THEN
          IF(IF_NIVPRI.EQ.0)IF_NIVPRI=1
          WRITE(*,7)SPERANGE,ENT,SPFILE
7         FORMAT('  orbital energies and <jx> in the range',F5.2,
     &    A,' from fermi energy,'/'  written to file ',A)
          IF(NVAR.NE.3)THEN
               IF(NEUMOV.EQ.0)
     &         WRITE(*,*)' only orbitals for ',ISOTEXT(ISO2),
     &        ' written, due to memory arrangement'
          ENDIF
      ELSE
          WRITE(*,*)' orbitals not written to file'
          IF(IF_NIVPRI.EQ.0)WRITE(*,*)' and not printed out either'
      ENDIF
      IF(IF_TOTEOUT+IF_ISOENEROUT.GE.1)THEN
          WRITE(*,*)
          IF(IF_TOTEOUT.GE.1)THEN
            WRITE(*,*)' total nuclear quantities',
     &      ' for calculated frequencies'
            IF(ISTEP.GT.0)WRITE(*,*)' and interpolated spins'
            IF(ISTEP.LT.0)WRITE(*,*)
     &       ' total routhian: R = E - wI printed instead for E'
          ENDIF
          IF(IF_ISOENEROUT.GE.1)WRITE(*,*)' proton and neutron total'
          WRITE(*,*)' quantities written to file :',ENERFILE
      ELSE
          WRITE(*,*)' no total quantities written to file'
          IF(IF_NIVPRI.GE.2)WRITE(*,*)' and not calculated either'
      ENDIF
      IF(IF_LEVSAV.LT.0)THEN
          WRITE(*,*)' will try to read sp-values from file system:',
     &     SAVFILE
          WRITE(*,*)' if deemed possible'
      ENDIF
      IF(IF_LEVSAV.GT.0)THEN
          WRITE(*,*)' will save sp-values in file system:',SAVFILE
      ENDIF
      WRITE(*,'(/1X,80(''-''))')
      WRITE(*,*)' OK, let''s try it'
      WRITE(*,'(1H1)')
      RETURN
      END
      SUBROUTINE LOOP_READER(NDEF)
*********************************************************************
* input routine in NILSSON-CRANKER                                  *
* purpose: reading deformation definitions                          *
* called by: MAIN                                                   *
* common blocks changed: LOOPS                                      *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
* reading deformation mesh i.e. LOOP1
2     READ(*,*)INI,NEPS,(EPSI(I),I=1,MIN(ABS(NEPS),MDEF))
      CALL CHECK_READ(INI,2,*2)
      IF(NEPS.EQ.-3)THEN
          START=EPSI(1)
          STEP=EPSI(3)
          NEPS=(EPSI(2)-START)/STEP+1.5
          NEPS=MIN(NEPS,2*MDEF)
          DO 101 I=1,NEPS
101           EPSI(I)=START+(I-1)*STEP
      ENDIF
      NEPS=MIN(NEPS,MDEF)
3     READ(*,*)INI,NGAM,(GAMI(I),I=1,MIN(ABS(NGAM),MDEF))
      CALL CHECK_READ(INI,3,*3)
      IF(NGAM.EQ.-3)THEN
          START=GAMI(1)
          STEP=GAMI(3)
          NGAM=(GAMI(2)-START)/STEP+1.5
          NGAM=MIN(NGAM,2*MDEF)
          DO 102 I=1,NGAM
102           GAMI(I)=START+(I-1)*STEP
      ENDIF
      NGAM=MIN(NGAM,MDEF)
4     READ(*,*)INI,NEPS4,(EPS4I(I),I=1,MIN(ABS(NEPS4),MDEF))
      CALL CHECK_READ(INI,4,*4)
      IF(NEPS4.EQ.-3)THEN
          START=EPS4I(1)
          STEP=EPS4I(3)
          NEPS4=(EPS4I(2)-START)/STEP+1.5
          NEPS4=MIN(NEPS4,2*MDEF)
          DO 103 I=1,NEPS4
103           EPS4I(I)=START+(I-1)*STEP
      ENDIF
      NEPS4=MIN(NEPS4,MDEF)
* cranking frequency
5     READ(*,*)INI,NOMEGA,(OMEGA(I),I=1,MIN(ABS(NOMEGA),MOMEGA))
      CALL CHECK_READ(INI,5,*5)
      IF(NOMEGA.EQ.-3)THEN
          START=OMEGA(1)
          STEP=OMEGA(3)
          NOMEGA=(OMEGA(2)-START)/STEP+1.5
          NOMEGA=MIN(NOMEGA,MOMEGA)
          DO 104 I=1,NOMEGA
104           OMEGA(I)=START+(I-1)*STEP
      ENDIF
      NOMEGA=MIN(NOMEGA,MOMEGA)
      NDEF=NEPS*NGAM*NEPS4
      IF(IN_LEV.GT.0)THEN
* type of deformation variation
6         READ(*,*)INI,IDEF_TYP,NVAR,FI
          CALL CHECK_READ(INI,6,*6)
      ENDIF
      IF(IDEF_TYP.EQ.1)NDEF=NEPS*NEPS4
      IF(IDEF_TYP.GE.2)NDEF=NEPS
      IF(IDEF_TYP.EQ.3)THEN
          DO 400 I=NOMEGA,1,-1
              OMEGA(I+1)=OMEGA(I)
400       CONTINUE
          OMEGA(NOMEGA+2)=-20
          NOMEGA=1
      ENDIF
      RETURN
      END
      SUBROUTINE NILS_READER(SAVFILE)
*********************************************************************
* input routine in NILSSON-CRANKER                                  *
* purpose:reading or setting default values of                      *
*       Nilsson (Modified Oscillator) potential parameters          *
* called by: MAIN     calling: SETQN_NILS                           *
* common blocks changed: NILSPAR                                    *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      CHARACTER SAVFILE*(*),ENSAV*16,LQSAV*16
      NMAX(1)=MAXN
      NMAX(2)=MAXN
      IF(IN_LEV.LE.2)GO TO 501
12    READ(*,*)INI,NUU,NMAX(1),NMAX(2),G0,G1,KOEFF_RANGE
      CALL CHECK_READ(INI,12,*12)
      IF(NUU.LT.0)JCO=0
      NANT=MAX(NMAX(1),NMAX(2))+1
      DO 500 I=1,NANT
13        READ(*,*)INI,(RKAPPA(I,J),RMY(I,J),J=1,2)
          CALL CHECK_READ(INI,13+I-1,*13)
 500  CONTINUE
 501  CONTINUE
      IF(NMAX(1).LT.INT((3.*FLOAT(NUMBER(1)))**0.33333+1.2))THEN
          WRITE(ERR_TEX(NERR+1),1800)ISOTEXT(1)
1800      FORMAT('warning: too few ',A,' shells used, really')
          CALL ERR_OUT(1)
      ENDIF
      IF(NMAX(2).LT.((3.*FLOAT(NUMBER(2)))**0.33333+1.2))THEN
          WRITE(ERR_TEX(NERR+1),1800)ISOTEXT(2)
          CALL ERR_OUT(1)
      ENDIF
      BB=FLOAT(NUMBER(2)-NUMBER(1))/3./FLOAT(NUMBER(2)+NUMBER(1))
      AA=41./(NUMBER(2)+NUMBER(1))**(1./3.)
      RHOMEGA(1)=AA*(1-BB)
      RHOMEGA(2)=AA*(1+BB)
      CALL SETQN_NILS
      IF(IF_LEVSAV.NE.0)THEN
          IF(IF_LEVSAV.LT.0)THEN
              IERR=9
              OPEN(NSAV1FIL,FILE=SAVFILE,ACCESS='DIRECT',
     &        STATUS='OLD',RECL=1*16,ERR=900)
* VAX open statements! F77 standard should have 1* replaced by 4*
*    record length in bytes/block! VAX wants words/block!
              IERR=8
              READ(NSAV1FIL,REC=1,ERR=900)IPSAV,ENSAV,ITEST,IBLK
              IERR=7
              IF(IBLK.LT.
     &        2*MAX(NDIM(1,1)+NDIM(2,1),NDIM(1,2)+NDIM(2,2)))GO TO 900
              IERR=6
              OPEN(NSAV2FIL,FILE=ENSAV,ACCESS='DIRECT',STATUS='OLD',
     &        RECL=1*IBLK,ERR=900)
              ISL=LEN(SAVFILE)
              LQSAV='LQ-'//SAVFILE(1:MIN(13,ISL))
              IF_SAVDROP=+1
              OPEN(NSAV3FIL,FILE=LQSAV,STATUS='OLD',
     &        ACCESS='DIRECT',RECL=1*6,ERR=899)
              IF_SAVDROP=-1
              GO TO 899
900           CONTINUE
              WRITE(ERR_TEX(NERR+1),1801)IERR,SAVFILE(1:16)
1801          FORMAT(i4,' error opening ',A,' saving instead')
              CALL ERR_OUT(1)
              IF_LEVSAV=+1
              IF_SAVDROP=+1
899           CONTINUE
          ENDIF
          IF(IF_LEVSAV.GT.0)THEN
* VAX open statements! F77 standard should have 1* replaced by 4*
              OPEN(NSAV1FIL,FILE=SAVFILE,ACCESS='DIRECT',
     &        STATUS='UNKNOWN',RECL=1*16)
              ISL=LEN(SAVFILE)
              ENSAV='EN-'//SAVFILE(1:MIN(13,ISL))
              LQSAV='LQ-'//SAVFILE(1:MIN(13,ISL))
              ITEST=3
              IF(ISO1.EQ.ISO2)ITEST=2
              IBLK=2*MAX(NDIM(1,1)+NDIM(2,1),NDIM(1,2)+NDIM(2,2))
              NPROT=NUMBER(1)
              NNEU=NUMBER(2)
              WRITE(NSAV1FIL,REC=1)IPSAV,ENSAV,ITEST,IBLK
              OPEN(NSAV2FIL,FILE=ENSAV,ACCESS='DIRECT',STATUS='UNKNOWN',
     &        RECL=1*IBLK)
              OPEN(NSAV3FIL,FILE=LQSAV,STATUS='UNKNOWN',
     &        ACCESS='DIRECT',RECL=1*6)
              IF_SAVDROP=+1
          ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE LOOP1_SET
*********************************************************************
* routine in NILSSON-CRANKER                                        *
* purpose: setting parameters for loop1 (main deformation loop)     *
* called by: MAIN                                                   *
* common blocks changed: HAMPAR                                     *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      I4=MOD(ILOOP-1,NEPS4)+1
      KLOOP=(ILOOP-1)/NEPS4+1
      EPS4=EPS4I(I4)
      IF(IDEF_TYP.EQ.0)THEN
          IG=MOD(KLOOP-1,NGAM)+1
          KLOOP=(KLOOP-1)/NGAM+1
          IE=MOD(KLOOP-1,NEPS)+1
          EPS=EPSI(IE)
          GAM=GAMI(IG)
      ELSEIF(IDEF_TYP.EQ.-1)THEN
          IG=MOD(KLOOP-1,NGAM)+1
          KLOOP=(KLOOP-1)/NGAM+1
          IE=MOD(KLOOP-1,NEPS)+1
          X=EPSI(IE)
          Y=GAMI(IG)
          CALL FUN2(X,Y,EPS,GAM)
      ELSEIF(IDEF_TYP.EQ.1)THEN
          IE=MIN(NEPS,KLOOP)
          EPS=EPSI(IE)
          IG=MIN(NGAM,KLOOP)
          GAM=GAMI(IG)
      ELSEIF(IDEF_TYP.GE.2)THEN
          I4=MIN(NEPS4,ILOOP)
          EPS4=EPS4I(I4)
          IE=MIN(NEPS,ILOOP)
          EPS=EPSI(IE)
          IG=MIN(NGAM,ILOOP)
          GAM=GAMI(IG)
          IF(IDEF_TYP.EQ.3.AND.OMEGA(ILOOP+1).GT.-20)
     &     OMEGA(1)=OMEGA(ILOOP+1)
      ENDIF
      GAM=(GAM+120.)*.017453292
      RETURN
      END
      SUBROUTINE ORBITS
*********************************************************************
* routine number 21  at level 50  in NILSSON-CRANKER                *
* purpose: driver for a omega loop for one kind of particles        *
* called by: MAIN         calling: PREDIA, SPDIA,                   *
*    ORBITVAL, PRIORB, PRINIV, TOTISO, LAST_STRUT, SPROUT           *
* common blocks changed: TOTAL                                      *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      TID1=SECOND(U)
      IF(LEV_PRINT.GE.50)THEN
        WRITE(*,*)'*********entered ORBITS**************************'
      ENDIF
      IF(NVAR.EQ.3)IANT=1
      IF(NVAR.EQ.1)VAR(IANT)=EPS
      IF(NVAR.EQ.2)VAR(IANT)=57.2957795*GAM-120.
      IF(NVAR.EQ.4)VAR(IANT)=EPS4
      CALL PREDIA
      IF((LEV_PRINT.GE.21.AND.NVAR.EQ.3).OR.LEV_PRINT.GE.50)THEN
          IF(ISO.GT.ISO1)THEN
              WRITE(*,'(1H0)')
C             WRITE(*,'(1H1)')
          ELSE
              WRITE(*,'(80(''-''))')
          ENDIF
          WRITE(*,1954)ISOTEXT(ISO),NUMBER(1)+NUMBER(2),ANUC(NUMBER(1)),
     &    NUMBER(ISO)
1954      FORMAT(12X,A,I4,A,', particle number is',I4)
          WRITE(*,1955)EPS,57.2957*GAM-120,EPS4
1955      FORMAT(' deformation: E P S =',F7.3,'   G A M M A =',F8.2,
     &'   E P S 4 =',F7.3/)
          WRITE(*,1956)ESTRUT(ISO),ENT,OMOM*RHOMEGA(ISO)
1956      FORMAT(5X,' smoothed energy=',F8.3,A,5X,' hw0=',F8.5,'MeV')
          IF(IF_DROPCALC.GE.1.AND.ISO2-ISO1.NE.0)THEN
          WRITE(*,1957)DROPENERGY,DROPINERTIA
1957      FORMAT(8X,' L.D. energy=',F8.3,'MeV , L.D. inertia=',F8.2,
     &    ' 1/MeV')
          ENDIF
          IF(IF_BCS.GT.0)THEN
          WRITE(*,1958)BCS_EN(ISO),BCS_DEL(ISO),GPAIR(ISO)
1958      FORMAT(8X,' BCS energy:',F8.3,'MeV, Delta:',F5.3,'MeV, G=',
     &    F5.3,'MeV')
          ENDIF
      ENDIF
*
      DO 1 IOM=1,NOMEGA
*
          TID2=SECOND(U)
          IF(LEV_PRINT.GE.50.AND.IOM.GT.1)WRITE(*,'(1H1)')
          OMEG=OMEGA(IOM)
          IF(NVAR.EQ.3)VAR(IANT)=OMEG
          IF(NVAR.GE.4)IANT=1
          IF(OMEG.GT.0)CALL SPDIA
          IF(IF_LEVSAV.GE.0)CALL ORBITVAL
          IF(IF_LEVSAV.GT.0)CALL SPROUT
          IF(IF_STRUT.GE.2.AND.IOM.EQ.NOMEGA) CALL LAST_STRUT
          IF(LEV_PRINT.GE.53)WRITE(*,1962)EPS,57.29*GAM-120,EPS4,
     &    OMEG,IANT
1962      FORMAT(' eps=',F5.3,' gam=',F5.0,' eps4=',F5.3,' omega=',
     &    F5.3,' #',I3/)
          IF(IF_NIVPRI.LT.2)CALL TOTISO
          IF(LEV_PRINT.GE.50)CALL PRIORB
          IF(NVAR.EQ.3.AND.IOM.NE.NOMEGA)IANT=IANT+1
          IF(IANT.GT.MOMEGA)THEN
              WRITE(ERR_TEX(NERR+1),1800)IOM
1800          FORMAT(' omega loop error: IANT>MOMEGA at IOM=',I3)
              CALL ERR_OUT(1)
          ENDIF
          IF(LEV_PRINT.GE.50)THEN
              TID3=SECOND(U)
              WRITE(*,1961)' time for omega :',TID3-TID2
          ENDIF
1     CONTINUE
      INUM=IANT
      IF(LEV_PRINT.GE.21.AND.NVAR.EQ.3)CALL PRINIV(INUM)
      IF((LEV_PRINT.GE.21.AND.NVAR.EQ.3).OR.LEV_PRINT.GE.50)THEN
        TID2=SECOND(U)
        WRITE(*,1961)' time in ORBITS :',TID2-TID1
1961    FORMAT(A,F10.1,'s')
        WRITE(*,*)'I-----------leaving ORBITS----------------------I'
      ENDIF
      RETURN
      END
      SUBROUTINE OUT_WITHIT
*********************************************************************
* main output routine in NILSSON-CRANKER                            *
* called by: MAIN                                                   *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      CHARACTER HEAD1*100,HEAD2*100,HEAD3*100,CONFTEX*16
      NANT=IANT
      RENO=0
      IF(IF_STRUT.GE.2)THEN
          RENO=1./DROPINERTIA - 1./STRUTINERTIA
          WRITE(*,'(A,F6.4,A)')
     &     ' renormalized hw printed out. corr=',RENO,'*I MeV'
      ENDIF
      IF(MAX_SPIN.GT.0)NANT=MAX_SPIN+1
      DO 1 I=1,10
          HEAD1((I-1)*10+1:10*I)='          '
          HEAD2((I-1)*10+1:10*I)='          '
          HEAD3((I-1)*10+1:10*I)='          '
1     CONTINUE
      IA=NUMBER(1)+NUMBER(2)
      IF(ISTEP.GE.0)THEN
         WRITE(HEAD2,112)
      ELSE
         WRITE(HEAD2,114)
      ENDIF
      WRITE(HEAD3,113)
112   FORMAT(2X,'I    eps    gam   eps4     Etot    hw     dE      Q2',
     &'   Pcon.Ncon')
114   FORMAT(2X,'I    eps    gam   eps4     Rout    hw     dE      Q2',
     &'   Pcon.Ncon')
113   FORMAT(2X,'I    eps    gam   eps4     Etot    hw      Q2')
      IF(ANUC(NUMBER(1)).NE.ANUC(100))THEN
          WRITE(HEAD1,100)IA,ANUC(NUMBER(1)),NUMBER(2)
100       FORMAT(I5,A,' N=',I3,' optimal states ')
      ELSE
          WRITE(HEAD1,101)NUMBER(1),NUMBER(2)
101       FORMAT(' Z=  ',I2,' N=',I3,' optimal states ')
      ENDIF
      ISTART=1
      ISTOP=NANT
      WRITE(*,*)HEAD1
      WRITE(*,*)HEAD2
      WRITE(NTOTFILE,'(A)')HEAD1
      WRITE(NTOTFILE,108)NANT,HEAD2
108   FORMAT(I4,',',A)
      DO 6 IO=ISTART,ISTOP
          DE=0
          IF(IO.GT.ISTART)DE=UT_TOT(IO)-UT_TOT(IO-1)
          EPSU=UT_EPS(IO)
          GAMU=UT_GAM(IO)*57.2957-120
          IF(IDEF_TYP.LT.0)THEN
              CALL FUN2(UT_EPS(IO),UT_GAM(IO),EPSU,GAMU)
          ENDIF
          OMUT=UT_OM(IO)+RENO*UT_SPIN(IO)
          TOTUT=UT_TOT(IO)
          IF(ISTEP.LT.0)TOTUT=TOTUT-UT_OM(IO)*UT_SPIN(IO)
          IF(UT_OM(IO) .GT. 1.E-5)TOTUT=TOTUT-BCS_EN(1)-BCS_EN(2)
          WRITE(*,103)UT_SPIN(IO),EPSU,GAMU,
     &    UT_EP4(IO),TOTUT,OMUT,DE,UT_Q2(IO),
     &    UT_CONF(IO)
          IF(IF_TOTEOUT.EQ.1)THEN
             WRITE(NTOTFILE,109)UT_SPIN(IO),EPSU,GAMU,
     &       UT_EP4(IO),TOTUT,OMUT,DE,UT_Q2(IO),
     &       UT_CONF(IO)
          ENDIF
6     CONTINUE
103   FORMAT(1X,F7.2,F6.3,F7.1,F6.3,F9.4,F6.3,F7.3,F8.2,F9.3)
109   FORMAT(F9.2,',',F6.3,',',F5.0,',',F6.3,',',F9.4,',',
     &F5.3,',',F7.3,',',F6.0,',',F8.3,',')

      IF (ISTEP .GT. 0) THEN
         IMAX=UT_SPIN(ISTOP)
         NANTI = IMAX/ISTEP+1
         WRITE(NTOTFILE,'(A)')HEAD1
         WRITE(NTOTFILE,108)NANTI,HEAD3
         WRITE(*,'(//,2A,/)') ' CRUDE LINEAR INTERPOLATION TO ',
     &   'PRESCRIBED SPINS DISREGARDING BAND CROSSINGS'
         WRITE(*,*) HEAD3
         IO2 = ISTART+1
         DO 4 ISPIN = 0,IMAX,ISTEP
            DO 7 IO = IO2,ISTOP
               IF (UT_SPIN(IO) .GT. FLOAT(ISPIN)) THEN
                  IO2 = IO
                  GO TO 8
               ENDIF
    7       CONTINUE
    8       CONTINUE
            IO1=IO2-1
            OMISPIN = UT_OM(IO1) + (UT_OM(IO2)-UT_OM(IO1))*
     &      (ISPIN-UT_SPIN(IO1))/(UT_SPIN(IO2)-UT_SPIN(IO1))
            OMISPIN=OMISPIN+RENO*ISPIN
            EISPIN = UT_TOT(IO1) + (UT_TOT(IO2)-UT_TOT(IO1))*
     &      (ISPIN**2-UT_SPIN(IO1)**2)/(UT_SPIN(IO2)**2-UT_SPIN(IO1)**2)
            IF (ISPIN .NE. 0) EISPIN = EISPIN-BCS_EN(1)-BCS_EN(2)
            EPSU1=UT_EPS(IO1)
            GAMU1=UT_GAM(IO1)*57.2957-120
            IF (IDEF_TYP .LT. 0)
     &      CALL FUN2(UT_EPS(IO1),UT_GAM(IO1),EPSU1,GAMU1)
            Q2ISPIN = UT_Q2(IO1) + (UT_Q2(IO2)-UT_Q2(IO1))*
     &      (ISPIN-UT_SPIN(IO1))/(UT_SPIN(IO2)-UT_SPIN(IO1))
            WRITE(NTOTFILE,107) FLOAT(ISPIN),EPSU1,GAMU1,
     &      UT_EP4(IO1),EISPIN,OMISPIN,Q2ISPIN
107         FORMAT(F9.2,',',F6.3,',',F5.0,',',F6.3,',',F9.4,',',
     &      F5.3,',',F6.0,',')
            WRITE(*,104)FLOAT(ISPIN),EPSU1,GAMU1,
     &      UT_EP4(IO1),EISPIN,OMISPIN,Q2ISPIN
104         FORMAT(1X,F7.2,F6.3,F7.1,F6.3,F9.4,F6.3,F8.2)
    4    CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE TOTISO
*********************************************************************
* routine number 53 in NILSSON-CRANKER                              *
* purpose: calculating total quantities for one isospin             *
*  NOTE: this is the simplest variant of total energy possible      *
* called by:ORBITS                                                  *
* common blocks changed: ISOS                                       *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /HELPER/ EHELP(2*MAXD2),KNDX(2*MAXD2),RJHELP(2*MAXD2)
      DIMENSION ESUM(MSYM),RJSUM(MSYM),NSUM(MSYM)
      TID0=SECOND(U)
      IF(LEV_PRINT.GE.53)THEN
          WRITE(*,*)'******************TOTISO**************************'
      ENDIF
*
      ETOTI(IANT,ISO)=-ESTRUT(ISO)+BCS_EN(ISO)
      SPINI(IANT,ISO)=0
      Q2TOTI(IANT,ISO)=0
      ICONFI(IANT,ISO)=0
* unpaired  orbitals
      IHELP=0
      DO 100 ISYM=ISYM1,ISYM2
          RJSUM(ISYM)=0
          ESUM(ISYM)=0
          NSUM(ISYM)=0
          NDIM2=NDIM(ISYM,ISO)
          DO 101 II=1,NDIM2
               IHELP=IHELP+1
               EHELP(IHELP)=E(II,ISYM,IANT+MOVE)
               RJHELP(IHELP)=RJX(II,ISYM,IANT+MOVE)
               KNDX(IHELP)=(2*ISYM-3)*(II+2000)
               IHELP=IHELP+1
               EHELP(IHELP)=E(II+NDIM2,ISYM,IANT+MOVE)
               RJHELP(IHELP)=RJX(II+NDIM2,ISYM,IANT+MOVE)
               KNDX(IHELP)=(2*ISYM-3)*(1000+NDIM2+II)
101       CONTINUE
100   CONTINUE
      CALL SORT3R(EHELP,RJHELP,KNDX,IHELP)
      ISIG_SUM=0
      IPAR_SUM=0
      DO 103 I=1,NUMBER(ISO)
          INDEX=MOD(ABS(KNDX(I)),1000)
          ISYM=(3+SIGN(1,KNDX(I)))/2
          ISIG=2*ABS(KNDX(I))/1000-3
          ISIG_SUM=ISIG_SUM+ISIG
          IPAR_SUM=IPAR_SUM+ISYM
          ETOTI(IANT,ISO)=ETOTI(IANT,ISO)+EHELP(I)
          SPINI(IANT,ISO)=SPINI(IANT,ISO)+RJHELP(I)
          Q2TOTI(IANT,ISO)=Q2TOTI(IANT,ISO)+
     &    Q2(INDEX,ISYM,IANT+MOVE)
          RJSUM(ISYM)=RJSUM(ISYM)+RJHELP(I)
          ESUM(ISYM)= ESUM(ISYM)+EHELP(I)
          NSUM(ISYM)= NSUM(ISYM)+1
103   CONTINUE
      IPAR=MOD(IPAR_SUM,2)
      ISIG=MOD(ISIG_SUM+100,4)
      ICONFI(IANT,ISO)=ISIG*100+IPAR
      ETOTI(IANT,ISO)=ETOTI(IANT,ISO)+OMEG*SPINI(IANT,ISO)
      IF(IDEF_TYP.GE.0)THEN
          UT_EPS(IANT)=EPS
          UT_GAM(IANT)=GAM
      ELSE
          CALL XYCOORD(EPS,GAM,X,Y,IX,IY)
          UT_EPS(IANT)=X
          UT_GAM(IANT)=Y
      ENDIF
      UT_EP4(IANT)=EPS4
      IOTILL=0
      IF(IDEF_TYP.EQ.3)IOTILL=1
      UT_OM(IANT)=OMEGA(IANT+IOTILL)
      TID2=SECOND(U)
      TIME=TID2-TID0
      IF(LEV_PRINT.GE.53)THEN
          DO 12 ISYM=ISYM1,ISYM2
          WRITE(*,11)' parity',2*ISYM-3,' :',ESUM(ISYM),
     &     RJSUM(ISYM),NSUM(ISYM)
12        CONTINUE
11        FORMAT(A,I2,A,'E=',F10.3,' <jx>=',F7.2,' N=',I4)
10        FORMAT(A,'E=',F10.3,' <jx>=',F7.2,' <Q2>=',F8.2,' sigpar=',I5)
          WRITE(*,10)' total  quantities:',ETOTI(IANT,ISO),
     &    SPINI(IANT,ISO),Q2TOTI(IANT,ISO),ICONFI(IANT,ISO)
          WRITE(*,*)'**************************************************'
      ENDIF
      IF(LEV_PRINT.GE.53)WRITE(*,'(A,F8.1,A)')
     &' time in TOTISO: ',TIME,'s'
      RETURN
      END
      SUBROUTINE DEFSET
*********************************************************************
* routine number 22 at level 70  in NILSSON-CRANKER                 *
* purpose: stetting various parameters common to both particle      *
*       types prior to omega loop                                   *
* called by: MAIN         calling: OMO, DROP                        *
* common blocks changed: NILSPAR, FYRCOM                            *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /FYRCOM/ D40,D42,D44,EPS2,EPS22,OMX,OMY,OMZ,FAC,HAC
      DATA PI/3.14159265/
      TID1=SECOND(U)
      EPS2=EPS*COS(GAM)
      EPS22=EPS*SIN(GAM)/SQRT(0.75)
      D40=(5.*COS(GAM)**2+1)/6.
      D42=-SQRT(30.)/12.*SIN(2.*GAM)
      D44=SQRT(70.)/12.*SIN(GAM)**2
C this recipie is applicable in entire gamma-plane
      OMX=1.-2./3.*EPS*COS(GAM+2.*PI/3.)
      OMY=1.-2./3.*EPS*COS(GAM-2.*PI/3.)
      OMZ=1.-2./3.*EPS*COS(GAM)
      FAC=0.5*(SQRT(OMY/OMX)+SQRT(OMX/OMY))-1.
      HAC=0.5*(SQRT(OMY/OMX)-SQRT(OMX/OMY))
* volume conservation
      CALL OMO(OMOM)
      IF(IF_DROPCALC.NE.0.AND.IF_STRUT.GT.0)CALL DROP
      IF(LEV_PRINT.EQ.22.OR.LEV_PRINT.GE.70)THEN
        TID2=SECOND(U)
        WRITE(*,1961)' time in DEFSET :',TID2-TID1
1961    FORMAT(A,F8.1,'s')
      ENDIF
      RETURN
      END
      SUBROUTINE PRINIV(NANT)
*********************************************************************
* output routine in NILSSON-CRANKER                                 *
* purpose: printing q.p. energies, <jx> and total quantities        *
*          for a omega sequence, eventually on file                 *
* called by: ORBITS      calling: SORT                              *
* common blocks changed: none                                       *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /HELPER/ EHELP2(2*MAXD2),KNDX(2*MAXD2),RJHELP2(2*MAXD2)
      CHARACTER IDENT*30,SH_OUT*8,BCS_OUT*8
      GAMMA=57.29557*GAM-120
      IF(NVAR.GT.2)THEN
          WRITE(IDENT,500)EPS,GAMMA
500       FORMAT(' eps=',F5.3,' gam=',F5.0)
      ELSE
          WRITE(IDENT,501)
501       FORMAT('  deformation path            ')
      ENDIF
      IF(LEV_PRINT.GE.50)WRITE(*,'(1H1)')
      ISOS=ISO
      IF(NVAR.NE.3)THEN
          ISOS=ISO2
          IF(NEUMOV.GT.0)ISOS=ISO1
      ENDIF
607   CONTINUE
      MOVE=0
      SPEMIN=ELAM(ISOS)-SPERANGE
      SPEMAX=ELAM(ISOS)+SPERANGE
      IF(ISOS.EQ.ISO2)MOVE=NEUMOV
      IF(IF_NIVPRI.EQ.0)GO TO 600
      WRITE(*,300)EPS,GAMMA,EPS4,NUMBER(ISOS),ISOTEXT(ISOS)
 300  FORMAT(' eps=',F5.3,' gam=',F6.1,' eps4=',F6.3,
     &' N=',I4,1X,A)
 311  FORMAT(I10)
      KOMEGA=NANT
      MO=1
      IF(NANT.GT.10)THEN
        KOMEGA=NANT/2
        MO=2
      ENDIF
      IF(IF_SPOUT.EQ.1)THEN
 307  FORMAT(I4,1X,A,A)
          WRITE(NSPFILE,307)NUMBER(ISOS),ISOTEXT(ISOS),
     &    IDENT//' unrenorm. hw!'
          WRITE(NSPFILE,311)NANT
          WRITE(NSPFILE,310)(VAR(IOMEGA),IOMEGA=1,NANT)
      ENDIF
 310  FORMAT(15F8.3)
 391  FORMAT(4X,' omega=',10F12.3)
 301  FORMAT(4X,' vari =',10F12.3)
      KANT=0
      LOW_SKIP=0
      DO 100 ISYM=ISYM1,ISYM2
          NDIM2=NDIM(ISYM,ISOS)
          DO 957 IN=1,NDIM2
              LSKIP=0
              DO 955 IO=1,NANT
                  IF(E(IN,ISYM,IO+MOVE).LT.SPEMIN)THEN
                     LSKIP=1
                     GO TO 955
                  ENDIF
                  IF(E(IN,ISYM,IO+MOVE).GT.SPEMAX)GO TO 955
                  GO TO 954
955           CONTINUE
              IF(LSKIP.GT.0)LOW_SKIP=LOW_SKIP+1
              GO TO 953
954           KANT=KANT+1
              EHELP2(KANT)=E(IN,ISYM,1+MOVE)
              KNDX(KANT)=IN+1000*ISYM
953           CONTINUE
              DO 952 IO=1,NANT
                  IF(E(IN+NDIM2,ISYM,IO+MOVE).LT.SPEMIN)THEN
                     LSKIP=1
                     GO TO 952
                  ENDIF
                  IF(E(IN+NDIM2,ISYM,IO+MOVE).GT.SPEMAX)GO TO 952
                  GO TO 951
952           CONTINUE
              IF(LSKIP.GT.0)LOW_SKIP=LOW_SKIP+1
              GO TO 957
951           KANT=KANT+1
              EHELP2(KANT)=E(IN+NDIM2,ISYM,1+MOVE)
              KNDX(KANT)=-(IN+NDIM2+1000*ISYM)
 957      CONTINUE
      WRITE(*,*)
 100  CONTINUE
      IF(KANT.EQ.0)GO TO 101
      CALL SORT(EHELP2,KNDX,KANT)
      IF(NVAR.NE.3)THEN
      WRITE(*,301)(VAR(IOMEGA),IOMEGA=1,NANT,MO)
      ELSE
      WRITE(*,391)(VAR(IOMEGA),IOMEGA=1,NANT,MO)
      ENDIF
      WRITE(*,303)
  303 FORMAT(3X,'#   ',4X,8('    E   <jx>'))
      DO 956 INN=1,KANT
         IN=ABS(KNDX(INN))
         ISYM=IN/1000
         IN=MOD(IN,1000)
         IKOD=(2*ISYM-3)*((SIGN(1,KNDX(INN))+3)/2*10000+IN)
         WRITE(*,302)INN+LOW_SKIP,IKOD,(E(IN,ISYM,IO+MOVE),
     &   RJX(IN,ISYM,IO+MOVE),IO=1,NANT,MO)
         IF(INN+LOW_SKIP.EQ.NUMBER(ISOS))WRITE(*,*)' fermi level '
         IF(IF_SPOUT.EQ.1)THEN
         WRITE(NSPFILE,312)IN,IKOD,(E(IN,ISYM,IO+MOVE),
     &   RJX(IN,ISYM,IO+MOVE),IO=1,NANT)
         ENDIF
 312     FORMAT(1X,I3,I8,(8(F7.3,F7.2)))
 302     FORMAT(1X,I3,I8,10(F6.2,F6.2))
 956  CONTINUE
 101  CONTINUE
      NOLL=0
      IF(IF_SPOUT.EQ.1)WRITE(NSPFILE,312)NOLL,IDUM,
     & (EPS,EPS,IOP=1,NANT)
600   CONTINUE
      IF(IF_NIVPRI.GE.2)GO TO 608
      IF(NVAR.NE.3)THEN
          IF(NEUMOV.EQ.0)ISOS=ISO1
          WRITE(*,503)ISOTEXT(ISOS)
503       FORMAT(2X,A,' total quantities'/)
      ELSE
          WRITE(*,*)
          WRITE(*,*)'  ',ISOTEXT(ISOS),' total quantities:'
      ENDIF
606   CONTINUE
304   FORMAT(1X,F7.2,F6.3,F7.1,F6.3,F9.4,F6.3,F7.3,F8.2,I6,2X,2A)
309   FORMAT(4X,'   I    eps  gam   eps4   Etot    hw      dE     Q2',
     &'   sigpar  R-shell  E(BCS)')
504   FORMAT(I3,',I    eps      gam   eps4     Etot   ',
     &' hw    dE    Q2  sigpar  R-shell E(BCS)')
      WRITE(*,309)
      IF(IF_ISOENEROUT.EQ.1)THEN
          WRITE(NTOTFILE,307)NUMBER(ISOS),ISOTEXT(ISOS),IDENT
          WRITE(NTOTFILE,504)NANT
      ENDIF
305       FORMAT(F9.2,' ',F6.3,' ',F5.0,' ',F6.3,' ',F8.4,' ',
     &    F5.3,' ',F6.3,' ',F5.0,I5,2A)
      DO 308 IO=1,NANT
          DIF_E=0
          IF(IO.GT.1)DIF_E=ETOTI(IO,ISOS)-ETOTI(IO-1,ISOS)
          EPSU=EPS
          GAMU=GAM*57.2957-120
          EPS4U=EPS4
          SH_OUT = '        '
          BCS_OUT ='        '
          IF(NVAR.NE.3)THEN
              EPSU=UT_EPS(IO)
              GAMU=UT_GAM(IO)*57.2957-120
              EPS4U=UT_EP4(IO)
              IF(IO.EQ.NANT)THEN
                  EPSU=EPS
                  GAMU=GAM*57.2957-120
                  EPS4U=EPS4
              ENDIF
              IF(IDEF_TYP.LT.0)THEN
                  CALL FUN2(UT_EPS(IO),UT_GAM(IO),EPSU,GAMU)
              ENDIF
              IF (OMEGA(IO) .LT. 1.E-5)THEN
                 WRITE(SH_OUT,'(F8.3)') ESHELL(IO,ISOS)
                 IF (IF_BCS.EQ.1) WRITE(BCS_OUT,'(F8.3)') EBCS(IO,ISOS)
              ENDIF
          ELSE
            IF (IF_STRUT .GE. 2)THEN
              SUMD=ETOTI(IO,ISOS)-OMEGA(IO)*SPINI(IO,ISOS)-BCS_EN(ISOS)
              SUMS=(ESTRU_LAST(ISOS)-OMEGA(NOMEGA)*SPIN_STRU(ISOS)-
     &           ESTRUT(ISOS))*OMEGA(IO)**2/OMEGA(NOMEGA)**2
              WRITE(SH_OUT,'(F8.3)') SUMD-SUMS
              IF (OMEGA(IO) .LT. 1.E-5  .AND. IF_BCS .EQ. 1)
     &           WRITE(BCS_OUT,'(F8.3)') EBCS(1,ISOS)
            ENDIF
          ENDIF
          EOUT=ETOTI(IO,ISOS)
          IF (OMEGA(IO).GT.1.E-5) EOUT=ETOTI(IO,ISOS)-EBCS(1,ISOS)
          WRITE(*,304)SPINI(IO,ISOS),EPSU,GAMU,EPS4U,
     &    EOUT,OMEGA(IO),DIF_E,Q2TOTI(IO,ISOS),
     &    ICONFI(IO,ISOS),SH_OUT,BCS_OUT
          IF(IF_ISOENEROUT.EQ.1)THEN
              WRITE(NTOTFILE,305)SPINI(IO,ISOS),EPSU,GAMU,EPS4U,
     &        EOUT,OMEGA(IO),DIF_E,Q2TOTI(IO,ISOS),
     &        ICONFI(IO,ISOS),SH_OUT,BCS_OUT
          ENDIF
308   CONTINUE
608   IF((NVAR.NE.3).AND.ISOS.NE.ISO2)THEN
          ISOS=ISOS+1
          IF(IF_NIVPRI.LT.2) WRITE(*,503)ISOTEXT(ISOS)
          IF(NEUMOV.GT.0)GO TO 607
          GO TO 606
      ENDIF
      RETURN
      END
      SUBROUTINE PRIORB
*********************************************************************
* output routine in NILSSON-CRANKER                                 *
* purpose: printing q.p. energies, <jx> and total quantities        *
*          for a fixed omega                                        *
* called by: ORBITS      calling: SORT                              *
* common blocks changed: none                                       *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /HELPER/ EHELP1(MAXD2),EHELP2(MAXD2),KNDX1(MAXD2),
     &KNDX2(MAXD2),RJHELP2(2*MAXD2)
      WRITE(*,1962)EPS,57.29*GAM-120,EPS4,OMEG,IANT
1962  FORMAT(' eps=',F5.3,' gam=',F5.0,' eps4=',F6.3,' omega=',
     &F5.3,' #',I3/)
      DO 2 I=1,NDIM(1,ISO)*2
        EHELP1(I)=E(I,1,IANT+MOVE)
2       KNDX1(I)=I
      CALL SORT(EHELP1,KNDX1,NDIM(1,ISO)*2)
      DO 22 I=1,NDIM(2,ISO)*2
        EHELP2(I)=E(I,2,IANT+MOVE)
22      KNDX2(I)=I
      CALL SORT(EHELP2,KNDX2,NDIM(2,ISO)*2)
      WRITE(*,209)
 209  FORMAT(12X,'parity=-',42x,'parity=+'/
     &2('   #   sign      energy    <jx>   <Q2>  ',13X))
 208  FORMAT(I5,I7,3x,2F7.3,F7.2,17X,I4,I7,6x,2F7.3,F7.2)
      DO 3 I=1,30
          II=MAX(NUMBER(ISO)/2+15-I,0)
          IF(II.EQ.0)GO TO 3
          IN1=KNDX1(II)
          IN2=KNDX2(II)
          INA1=IN1
          INA2=IN2
          IF(IN1.GT.NDIM(1,ISO))INA1=-(IN1-NDIM(1,ISO))
          IF(IN2.GT.NDIM(2,ISO))INA2=-(IN2-NDIM(2,ISO))
          WRITE(*,208)II,INA1,E(IN1,1,IANT+MOVE),
     &    RJX(IN1,1,IANT+MOVE),Q2(IN1,1,IANT+MOVE),II,INA2,
     &    E(IN2,2,IANT+MOVE),RJX(IN2,2,IANT+MOVE),Q2(IN2,2,IANT+MOVE)
          IF(I.EQ.14)WRITE(*,*)
 3    CONTINUE
      WRITE(*,206)ETOTI(IANT,ISO),SPINI(IANT,ISO),
     &Q2TOTI(IANT,ISO)
206   FORMAT(' E=',F8.3,'MeV, <I>=',F8.2,' <Q2>=',F8.2)
      RETURN
      END
      SUBROUTINE ORBITVAL
*********************************************************************
* routine number 72  at level 80  in NILSSON-CRANKER                *
* purpose: calculating various expectation values for individual    *
*         orbitals                                                  *
* called by: ORBITS                                                 *
* common blocks changed: ORBITALS                                   *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      TID1=SECOND(U)
      IF(LEV_PRINT.EQ.72.OR.LEV_PRINT.GE.80)THEN
        WRITE(*,*)'************ ORBITVAL output*************'
      ENDIF
      DO 1 ISYM=ISYM1,ISYM2
      NDIM2=NDIM(ISYM,ISO)
      DO 2 I=1,2*NDIM2
        RJX(I,ISYM,IANT+MOVE)=0.
        Q2(I,ISYM,IANT+MOVE)=0
        DO 3 J=1,NDIM2
          RJSUM=0
          Q2SUM=0
          DO 4 K=1,NDIM2
              RJSUM=RJSUM+VEC(K,I,ISYM)*RJXIJ(J,K,ISYM)
              IF(IF_Q2.LT.ISO)GO TO 4
              Q2SUM=Q2SUM+VEC(K,I,ISYM)*Q2IJ(J,K,ISYM)
4         CONTINUE
          RJX(I,ISYM,IANT+MOVE)=RJX(I,ISYM,IANT+MOVE)-
     &    VEC(J,I,ISYM)*RJSUM
          Q2(I,ISYM,IANT+MOVE)=Q2(I,ISYM,IANT+MOVE)+
     &    VEC(J,I,ISYM)*Q2SUM
3       CONTINUE
        IF(I.GT.NDIM2)RJX(I,ISYM,IANT+MOVE)=-RJX(I,ISYM,IANT+MOVE)
2     CONTINUE
1     CONTINUE
      IF(LEV_PRINT.EQ.72.OR.LEV_PRINT.GE.80)THEN
          DO 32 ISYM=ISYM1,ISYM2
          NDIM2=NDIM(ISYM,ISO)
          WRITE(*,*)' parity=',2*ISYM-3
          WRITE(*,30)'      ',(I,I=NDIM2-4,NDIM2+5)
30        FORMAT(A,10I10)
          WRITE(*,31)' e    ',(E(I,ISYM,IANT+MOVE),I=NDIM2-4,NDIM2+5)
31        FORMAT(A,10F10.3)
          WRITE(*,31)' <jx> ',(RJX(I,ISYM,IANT+MOVE),I=NDIM2-4,NDIM2+5)
32        CONTINUE
      ENDIF
      IF(LEV_PRINT.GE.52)THEN
        TID2=SECOND(U)
        WRITE(*,1961)' time in ORBITVAL :',TID2-TID1
1961    FORMAT(A,F8.1,'s')
      ENDIF
      IF(LEV_PRINT.EQ.72.OR.LEV_PRINT.GE.80)THEN
        WRITE(*,*)'I-----------------leaving ORBITVAL--------------I'
      ENDIF
      RETURN
      END
      SUBROUTINE SPDIA
*********************************************************************
* routine number 73  at level 80  in NILSSON-CRANKER                *
* purpose: performing an unpaired diagonalization at some omega     *
* called by: ORBITS   calling: MATDIA, OVER, SPREAD                 *
* common blocks changed: ORBITALS                                   *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /DIAGVEC/DMAT(MAXDIM,MAXDIM),EIGENE(MAXDIM),
     &EIGENV(MAXDIM,MAXDIM)
      IF(LEV_PRINT.EQ.73.OR.LEV_PRINT.GE.80)THEN
          WRITE(*,*)'*********entered SP dia*********************'
      ENDIF
      TID1=SECOND(U)
      TID2=TID1
*
*  trying to read saved values
      IF(IF_LEVSAV.LT.0)CALL SPREAD
      IF(IF_LEVSAV.LT.0)GO TO 101
*
      DO 100 ISYM=ISYM1,ISYM2
      NDIM2=NDIM(ISYM,ISO)
      IF(LEV_PRINT.EQ.73.OR.LEV_PRINT.GE.80)write(*,*)' isym,ndim2',
     &ISYM,NDIM2
      DO 400 I=1,NDIM2
      DO 399 J=1,NDIM2
        DMAT(I,J)=0
399   CONTINUE
400   CONTINUE
*
      TID3=SECOND(U)
      TIMEZ=TID3-TID2
*
50    CONTINUE
      DO 407 ISIG=1,2
        TECK=3-2*ISIG
      DO 401 I=1,NDIM2
      DO 402 J=1,NDIM2
          DMAT(J,I)=EIJ(J,I,ISYM)+TECK*OMEG*RJXIJ(J,I,ISYM)
*   first signature calculated by h0-w jx
*   other signature calculated by h0+w jx
402   CONTINUE
401   CONTINUE
*
      TID2=SECOND(U)
      TIMES=TID2-TID3
*
      CALL MATDIA(DMAT,EIGENE,NDIM2,EIGENV)
*
      TID3=SECOND(U)
      TIMED=TID2-TID3
      ITILL=(ISIG-1)*NDIM2
*
*     place crosser call here
*
      DO 403 I=1,NDIM2
          E(I+ITILL,ISYM,IANT+MOVE)=EIGENE(I)
          DO 404 J=1,NDIM2
              VEC(J,I+ITILL,ISYM)=EIGENV(J,I)
404       CONTINUE
403   CONTINUE
*
      TID2=SECOND(U)
      TIMESS=TID2-TID3
      IF(LEV_PRINT.EQ.73.OR.LEV_PRINT.GE.80)
     &WRITE(*,1970)TIMEZ,TIMES,TIMED,TIMESS
1970  FORMAT(' times: zero:',f6.1,' set:',f6.1,' diag:',f8.1,
     &' store :',f6.1)
407   CONTINUE
*
100   CONTINUE
*
101   CONTINUE
      TID3=SECOND(U)
      IF(LEV_PRINT.GE.52)WRITE(*,1969)' SPDIAG time:',TID3-TID1
1969  FORMAT(A,F8.1,'s')
      IF(LEV_PRINT.EQ.73.OR.LEV_PRINT.GE.80)THEN
       WRITE(*,*)'I-----------leaving SPDIA------------------I'
      ENDIF
      RETURN
      END
      SUBROUTINE SORTIND(EHELP,KNDX,IHELP)
*********************************************************************
* routine number 82 at level 90 in NILSSON-CRANKER                  *
* purpose: sorting both symmetry groups of one-particle energies    *
*       and returning them in vectors *HELP. original index         *
*       in vector KNDX                                              *
* called by: SPDIA     calling: SORT                                *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      DIMENSION EHELP(*),KNDX(*)
      TID1=SECOND(U)
      IHELP=0
      DO 1 ISYM=ISYM1,ISYM2
      NDIM2=NDIM(ISYM,ISO)
      DO 1 II=1,NDIM2
          IHELP=IHELP+1
          EHELP(IHELP)=E(II,ISYM,IANT+MOVE)
          KNDX(IHELP)=(2*ISYM-3)*(II+2000)
          IHELP=IHELP+1
          EHELP(IHELP)=E(II+NDIM2,ISYM,IANT+MOVE)
          KNDX(IHELP)=(2*ISYM-3)*(1000+NDIM2+II)
1     CONTINUE
      CALL SORT(EHELP,KNDX,IHELP)
      IF(LEV_PRINT.EQ.82.OR.LEV_PRINT.GE.90)THEN
        TID2=SECOND(U)
        WRITE(*,1961)' time in SORTIND:',TID2-TID1
1961    FORMAT(A,F6.1,'s')
      ENDIF
      RETURN
      END
      SUBROUTINE SETMAT_NILS
*********************************************************************
* routine number 81  at level 90  in NILSSON-CRANKER                *
* purpose: setting matrices for a Nilsson model diagonalization     *
* called by: ORBITS  calling: CLEBI, MRHOT                          *
* common blocks changed: HAMILTON                                   *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /CMN1/L2F(18)
      COMMON /FYRCOM/ D40,D42,D44,EPS2,EPS22,OMX,OMY,OMZ,FAC,HAC
      INTEGER RNP,RN,RLP,RL,RJP,RJ,OMP,OM
      TID1=SECOND(U)
      DATA PI/3.14159265/
      GAC=(1.+EPS2/3.+EPS22/2.)*(1.+EPS2/3.-EPS22/2.)*(1.-2.*EPS2/3.)
      RMEV=1
      RFM2=1
      IF(IF_MEV.GT.0)THEN
          RMEV=RHOMEGA(ISO)*OMOM
          RFM2=41.5/RMEV
      ENDIF
*     Transition quadrupole moment, Q22(x) (x=rotation axis):
*      (static moment, Q20(x), :SQ20=0.0,   CQ20=1.0
*                      Q20(z)  :SQ20=SQRT(1.5),   CQ20=-0.5)
          SQ20=1.0
          CQ20=0
      DO 211 ISYM=ISYM1,ISYM2
      NDIM2=NDIM(ISYM,ISO)
      DO 210 IN=1,NDIM2
        RNP=NVAL(IN,ISYM)
        RLP=LVAL(IN,ISYM)
        RJP=JVAL(IN,ISYM)
        OMP=MVAL(IN,ISYM)
        DO 210 IP=1,NDIM2
          RN=NVAL(IP,ISYM)
          RL=LVAL(IP,ISYM)
          RJ=JVAL(IP,ISYM)
          OM=MVAL(IP,ISYM)
          IF(ABS(RNP-RN).GT.NUU.OR.ABS(OMP-OM).GT.8) GO TO 303
          IF(ABS(RJ-RJP).GT.JCO) GO TO 303
          E20=0.0
          E22=0.0
          ELZ=0.0
          END=0.0
          E40=0.0
          E42=0.0
          E44=0.0
          Q20S=0
          Q22S=0
          RHO2S=0
          IF(ABS(OMP-OM).EQ.8) GO TO 305
          IF(ABS(OMP-OM).EQ.4) GO TO 306
*    OM=OMP
          IF(RNP.EQ.RN)THEN
            E20=-2.0*EPS2/3.0*CLEBI(RJ,4,RJP,-1,0,-1)*
     $      CLEBI(RJ,4,RJP,OM,0,OMP)
          ENDIF
          IF(ABS(RNP-RN).LE.2 .AND. OM.EQ.OMP)THEN
            Q20S=CLEBI(RJ,4,RJP,-1,0,-1)*CLEBI(RJ,4,RJP,OM,0,OMP)
          ENDIF
          E40=EPS4*D40*
     $    CLEBI(RJ,8,RJP,-1,0,-1)*CLEBI(RJ,8,RJP,OM,0,OMP)
          IF(RNP.EQ.RN.AND.RLP.EQ.RL) ELZ=-FAC*
     $    ((OM-1.)/2.*CLEBI(2*RL,1,RJP,OM-1,1,OM)*
     $    CLEBI(2*RL,1,RJ,OM-1,1,OM)+
     $    (OM+1.)/2.*CLEBI(2*RL,1,RJP,OM+1,-1,OM)*
     $    CLEBI(2*RL,1,RJ,OM+1,-1,OM))
          GO TO 304
  306     CONTINUE
* OM=OMP+-4  ie K'=K+-2
          IF(RNP.EQ.RN)THEN
            E22=EPS22/SQRT(6.0)*CLEBI(RJ,4,RJP,-1,0,-1)*
     $      (CLEBI(RJ,4,RJP,OM,4,OMP)+CLEBI(RJ,4,RJP,OM,-4,OMP))
          ENDIF
          IF(ABS(RNP-RN).LE.2)THEN
            Q22S=CLEBI(RJ,4,RJP,-1,0,-1)*
     $      (CLEBI(RJ,4,RJP,OM,4,OMP)+CLEBI(RJ,4,RJP,OM,-4,OMP))
          ENDIF
          IF(ABS(RN-RNP).EQ.2) END= 1./SQRT(1.5)*(RNP-RN)/2.*
     &    HAC*CLEBI(RJ,4,RJP,-1,0,-1)*
     $    (CLEBI(RJ,4,RJP,OM,4,OMP)-CLEBI(RJ,4,RJP,OM,-4,OMP))
          E42=EPS4*D42*
     $    CLEBI(RJ,8,RJP,-1,0,-1)*
     $    (CLEBI(RJ,8,RJP,OM,4,OMP)+CLEBI(RJ,8,RJP,OM,-4,OMP))
          GO TO 304
  305     CONTINUE
* omp=om+-8
          E44=EPS4*D44*
     $    CLEBI(RJ,8,RJP,-1,0,-1)*
     $    (CLEBI(RJ,8,RJP,OM,8,OMP)+CLEBI(RJ,8,RJP,OM,-8,OMP))
  304     CONTINUE
* common finale:
          CALL MRHOT(RNP,RLP,0,RN,RL,0,2,RME)
          EIJ(IN,IP,ISYM)=RME*SQRT((RJ+1.)/(RJP+1.))*
     &     (E20+E22+E40+E42+E44)*RMEV
          IF(RJ.EQ.RJP .AND. OM.EQ.OMP) RHO2S=RME*RFM2
          Q20S=2.*RME*SQRT((RJ+1.)/(RJP+1.))*Q20S*RFM2
          Q22S=RME*SQRT((RJ+1.)/(RJP+1.))*Q22S*RFM2
          Q20=((1.+EPS2/3.-EPS22*EPS22/6.)*Q20S+
     $    EPS22/SQRT(6.)*(1.-2.*EPS2/3.)*Q22S+
     $    (2.*EPS2/3.*(1.+EPS2/3.)-EPS22*EPS22/6.)*RHO2S)/GAC
          Q22=((1.+EPS2/3.)*Q22S+
     $    EPS22/SQRT(24.)*Q20S-
     $    EPS22/SQRT(6.)*RHO2S)/GAC*(1.-2.*EPS2/3.)
          Q2IJ(IN,IP,ISYM)=CQ20*Q20+SQ20*Q22
          RJXIJ(IN,IP,ISYM)=ELZ+ RME*SQRT((RJ+1.)/(RJP+1.))*END
          N1=MIN(RN+1,15)
          IF(IN.EQ.IP)THEN
          RJXIJ(IN,IP,ISYM)=RJXIJ(IN,IP,ISYM)-FLOAT(OM)/2
          EIJ(IN,IP,ISYM)=EIJ(IN,IP,ISYM)+(FLOAT(RN)+1.5-
     $RKAPPA(N1,ISO)/OMOM*(RJ*(RJ+2.)/4.-RL*(RL+1.)-0.75)-
     $RKAPPA(N1,ISO)/OMOM*RMY(N1,ISO)*(RL*(RL+1.)-
     &FLOAT(L2F(RN+1))))*RMEV
          ENDIF
  303 CONTINUE
  210 CONTINUE
      IF(LEV_PRINT.EQ.81.OR.LEV_PRINT.GE.90)THEN
        WRITE(*,*)'************SETMAT_NILS output*******************'
      WRITE(*,'(A,I2)')' parity=',2*ISYM-3
      WRITE(*,*)' single particle energy matrix:'
         WRITE(*,203)'   #   ',(IN,IN=1,NDIM2)
         WRITE(*,203)'   N   ',(NVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,203)'   L   ',(LVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,201)'   J   ',(JVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,201)'   K   ',(MVAL(IN,ISYM),IN=1,NDIM2)
         DO 975 IN=1,NDIM2
             WRITE(*,202)MVAL(IN,ISYM),(EIJ(IN,INP,ISYM),INP=1,NDIM2)
 975     CONTINUE
      WRITE(*,*)
      WRITE(*,*)' angular momentum matrix:'
         WRITE(*,203)'   #   ',(IN,IN=1,NDIM2)
         WRITE(*,203)'   N   ',(NVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,203)'   L   ',(LVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,201)'   J   ',(JVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,201)'   K   ',(MVAL(IN,ISYM),IN=1,NDIM2)
         DO 9765 IN=1,NDIM2
             WRITE(*,202)MVAL(IN,ISYM),(RJXIJ(IN,INP,ISYM),INP=1,NDIM2)
 9765    CONTINUE
      WRITE(*,*)
      WRITE(*,*)' quadrupole moment matrix:'
         WRITE(*,203)'   #   ',(IN,IN=1,NDIM2)
         WRITE(*,203)'   N   ',(NVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,203)'   L   ',(LVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,201)'   J   ',(JVAL(IN,ISYM),IN=1,NDIM2)
         WRITE(*,201)'   K   ',(MVAL(IN,ISYM),IN=1,NDIM2)
         DO 8765 IN=1,NDIM2
             WRITE(*,202)MVAL(IN,ISYM),(Q2IJ(IN,INP,ISYM),INP=1,NDIM2)
 8765    CONTINUE
203     FORMAT(A,(20I6))
201     FORMAT(A,(20(I4,'/2')))
202     FORMAT(I5,'/2',(20F6.2))
        TID2=SECOND(U)
        WRITE(*,1961)' time in SETMAT_NILS:',TID2-TID1
1961    FORMAT(A,F8.1,'s')
        WRITE(*,*)'I-----------------------------------------------I'
      ENDIF
211   CONTINUE
      RETURN
      END
      SUBROUTINE SETQN_NILS
*********************************************************************
* routine in NILSSON-CRANKER                                        *
* purpose: setting up quantum number vectors                        *
* called by: NILS_READER                                            *
* common blocks affected: QUANTA, HAMPAR                            *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      NMAXA=MAX(NMAX(1),NMAX(2))
      IF(NMAXA.GT.MAXN)THEN
          WRITE(ERR_TEX(NERR+1),1800)
1800      FORMAT(' MORE SHELLS WANTED THAN DECLARED FOR ')
          CALL ERR_OUT(1)
          WRITE(*,*) 'change the MAXN &MAXDIM parameters!'
          WRITE(*,*) ' I WILL NOW REDUCE MAXIMAL SHELL USED!'
          WRITE(*,*) ' B E W A R E!!!!'
          NMAX(1)=MIN(NMAX(1),MAXN)
          NMAX(2)=MIN(NMAX(2),MAXN)
          NMAXA=MAX(NMAX(1),NMAX(2))
      ENDIF
      CALL HELP
      DO 1004 ISYM=ISYM1,ISYM2
      NSHM=NMAXA
      IF(MOD(NSHM-ISYM,2).EQ.1)NSHM=NSHM-1
      II=0
      NDIM(ISYM,1)=0
      NDIM(ISYM,2)=0
      IST=2-ISYM
      DO 1005 NSH=IST,NSHM,2
      DO 1001 LSH=NSH,0,-2
        ITEC=+1
1002    JSH=2*LSH+ITEC
        IOM=JSH
        IF(MOD((JSH-1)/2,2).EQ.1)IOM=JSH-2
        DO 1003 LOM=IOM,-JSH,-4
          II=II+1
          IF(NSH.LE.NMAX(1))NDIM(ISYM,1)=NDIM(ISYM,1)+1
          IF(NSH.LE.NMAX(2))NDIM(ISYM,2)=NDIM(ISYM,2)+1
          NVAL(II,ISYM)=NSH
          JVAL(II,ISYM)=JSH
          LVAL(II,ISYM)=LSH
          MVAL(II,ISYM)=LOM
1003    CONTINUE
        IF(ITEC.GT.0)THEN
          ITEC=-1
          GO TO 1002
        ENDIF
1001  CONTINUE
1005  CONTINUE
1004  CONTINUE
      RETURN
      END
      SUBROUTINE PREDIA
*********************************************************************
* routine number 41  at level 70  in NILSSON-CRANKER                *
* purpose: setting inititial (omega=0) values for strutinsky        *
*       energy                                                      *
* called by : ORBITS   calling: SPDIA, SORTIND, STRUTS, BCS         *
* common blocks changed: STRUTI                                     *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /HELPER/ EHELP2(2*MAXD2),KNDX(2*MAXD2),RJHELP2(2*MAXD2)
      DATA GAMST/1.0/
* GAMST is Strutinsky smearing width, also present in LAST_STRUT
      IF(LEV_PRINT.EQ.41.OR.LEV_PRINT.GE.70)THEN
        WRITE(*,*)'***************PREDIA output*********************'
      ENDIF
      TID1=SECOND(U)
      OMEG=0
      CALL SPDIA
      CALL SORTIND(EHELP2,KNDX,J)
      IF(IF_STRUT.GT.0)THEN
          GAMSTR=GAMST*RHOMEGA(ISO)*OMOM
          IF(IF_MEV.EQ.0)GAMSTR=GAMST
          CALL STRUTS(EHELP2,RJHELP2,OMEG,J,NUMBER(ISO),GAMSTR,
     &    ESMOTH,AISMOT,SUM,EF,AG1,AG2,SPIN)
          ESTRUT(ISO)=ESMOTH
          ELAM(ISO)=EF
          ESHELL(IANT,ISO) = SUM-ESTRUT(ISO)
      ELSE
          ESUM=0
          ELAM(ISO)=EHELP2(NUMBER(ISO))
          DO 101 I=1,NUMBER(ISO)
              ESUM=ESUM+EHELP2(I)
101       CONTINUE
          ESTRUT(ISO)=ESUM
      ENDIF
      IF(IF_BCS.GE.1)THEN
          A=NUMBER(1)+NUMBER(2)
          IA=A
          GPAIR(ISO)=(G0-(2*ISO-3)*G1*(NUMBER(2)-NUMBER(1))/A)/A
          G=GPAIR(ISO)
          NN=NUMBER(ISO)
          RMEV=OMOM*RHOMEGA(ISO)
          IF(IF_MEV.NE.0)THEN
              RMEV=1
          ENDIF
          CALL BCS(EHELP2,J,NN,IA,DEL,BCSLAM,DEN,RMEV,G,BCSDNDL)
          BCS_EN(ISO)=DEN
          BCS_DEL(ISO)=DEL
      ELSE
          BCS_EN(ISO)=0.
          BCS_DEL(ISO)=0.
      ENDIF
      EBCS(IANT,ISO) = BCS_EN(ISO)
      IF(LEV_PRINT.EQ.41.OR.LEV_PRINT.GE.70)THEN
        WRITE(*,1966)ESTRUT(ISO),OMOM,SUM,SUM-ESTRUT(ISO)
1966    FORMAT(' smoothed energy:',F10.3,' omom=',F8.4,
     &  ' sum energy:',F10.3,' shell energy:',F10.3)
        TID2=SECOND(U)
        WRITE(*,1961)' time in PREDIA :',TID2-TID1
1961    FORMAT(A,F8.1,'s')
        WRITE(*,*)'I----------leaving PREDIA-------------------------I'
      ENDIF
      RETURN
      END
      SUBROUTINE LAST_STRUT
*********************************************************************
* routine number 41  at level 80  in NILSSON-CRANKER                *
* purpose: performing the last-omega strutinsky calculation needed  *
*     for renomalization om the moment of inertia                   *
* called by: ORBITS  calling: SORT3R, STRUTS                        *
* common blocks changed:  STRUTI                                    *
*                                                                   *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /HELPER/ EHELP(2*MAXD2),KNDX(2*MAXD2),RJHELP(2*MAXD2)
      DATA GAMST/1.0/
* GAMST is Strutinsky smearing width, also present in PREDIA
      TID1=SECOND(U)
      IHELP=0
      DO 1 ISYM=ISYM1,ISYM2
      NDIM2=NDIM(ISYM,ISO)
      DO 2 II=1,NDIM2
          IHELP=IHELP+1
          EHELP(IHELP)=E(II,ISYM,IANT+MOVE)
          RJHELP(IHELP)=RJX(II,ISYM,IANT+MOVE)
          KNDX(IHELP)=(2*ISYM-3)*(II+2000)
          IHELP=IHELP+1
          EHELP(IHELP)=E(II+NDIM2,ISYM,IANT+MOVE)
          RJHELP(IHELP)=RJX(II+NDIM2,ISYM,IANT+MOVE)
          KNDX(IHELP)=(2*ISYM-3)*(1000+NDIM2+II)
2     CONTINUE
1     CONTINUE
      CALL SORT3R(EHELP,RJHELP,KNDX,IHELP)
      GAMSTR=GAMST*RHOMEGA(ISO)*OMOM
      IF(IF_MEV.EQ.0)GAMSTR=GAMST
      CALL STRUTS(EHELP,RJHELP,OMEG,IHELP,NUMBER(ISO),GAMSTR,
     &ESMOTH,AISMOT,SUM,EF,AG1,AG2,SPIN)
      ESTRU_LAST(ISO)=ESMOTH
      SPIN_STRU(ISO)=AISMOT
      IF(ISO.EQ.2)THEN
          SPI=SPIN_STRU(1)+SPIN_STRU(2)
          EDIF=ESTRU_LAST(1)+ESTRU_LAST(2)-ESTRUT(1)-ESTRUT(2)
          STRUTINERTIA=SPI**2/EDIF*0.5
      ENDIF
      TID2=SECOND(U)
      IF(LEV_PRINT.EQ.41.OR.LEV_PRINT.GE.70)THEN
          WRITE(*,*)'****************LAST-STRUT************************'
          SUM=SUM-OMEG*SPIN
          ESMOTH=ESMOTH-OMEG*AISMOT
          WRITE(*,10)ESMOTH,AISMOT,SUM,SPIN,SUM-ESMOTH,TID2-TID1
10        FORMAT(' smoothed energy and spin in rotating system:',
     &       F9.3,F6.2/' discrete energy and spin in rotating system:',
     &       F9.3,F6.2/
     &    ' quasi-shell energy:',F9.4/
     &    ' time in LAST_STRUT:',F8.1,'s')
          IF(ISO.EQ.2)THEN
            WRITE(*,11)STRUTINERTIA
11          FORMAT(' Strutinsky inertia from spins:',F7.2)
          ENDIF
          WRITE(*,*)'**************************************************'
      ENDIF
      RETURN
      END
      SUBROUTINE DROP
*********************************************************************
* routine number 27  at level 70  in NILSSON-CRANKER                *
* purpose: calculating liquid drop energy and inertia               *
* called by: DEFSET  calling: BSC24G, ROT2G, LDROP                  *
* common blocks changed:  STRUTI                                    *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      COMMON /FYRFAC/DA,DB,DC
      REAL LDROP
      DOUBLE PRECISION DA,DB,DC,GAMMAR
      TID1=SECOND(U)
      IF(IF_LEVSAV.LT.0.AND.IF_SAVDROP.LT.0)THEN
          IPOST=ILOOP
          READ(NSAV3FIL,REC=IPOST,ERR=99)BS,BC,XJ,YJ,ZJ,OMOMA
          GO TO 98
99        IF_SAVDROP=+1
98        CONTINUE
      ENDIF
      IF(IF_SAVDROP.LT.0)GO TO 97
      GAMMAR=DBLE(GAM-2.094395)
      DA=(5.D0*DCOS(GAMMAR)**2+1.D0)/6.D0
      DB=-DSQRT(30.D0)/12.D0*DSIN(2.D0*GAMMAR)
      DC=DSQRT(70.D0)/12.D0*DSIN(GAMMAR)**2
      D40=DA
      D42=DB
      D44=DC
      GUMMA=GAM*57.2957-120.
      CALL BSC24G(EPS,GUMMA,EPS4,OMOMA,BS,BC,12,12)
      CALL ROT2G(16,XJ,YJ,ZJ)
* replace single precision OMOM with OMOMA as we have now have it:
      OMOM=OMOMA
97    IA=NUMBER(1)+NUMBER(2)
      IZ=NUMBER(1)
      DROPENERGY=LDROP(IA,IZ,BS,BC)
      DROPINERTIA=FLOAT(IA)**(5./3.)*XJ/72.
      IF(IF_LEVSAV.LT.0.AND.IF_SAVDROP.GT.0)THEN
          IPOST=ILOOP
          WRITE(NSAV3FIL,REC=IPOST)BS,BC,XJ,YJ,ZJ,OMOMA
      ENDIF
      TID2=SECOND(U)
      IF(LEV_PRINT.EQ.27.OR.LEV_PRINT.GE.70)THEN
          WRITE(*,*)'*******************DROP***************************'
          IF(IF_SAVDROP.NE.0)WRITE(*,*)' liquid-drop saving status:',
     &    IF_SAVDROP
          WRITE(*,7)
 7        FORMAT('  OMOM',8X,'BS',9X,'BC',9X,' ENERGY',
     &    4X,'JX',9X,'JY',9X,'JZ')
          WRITE(*,8)OMOMA,BS,BC,DROPENERGY,XJ,YJ,ZJ
 8        FORMAT(F11.6,2F11.6,F7.3,3F11.6)
          WRITE(*,9)DROPINERTIA,TID2-TID1
 9        FORMAT(' L.D. inertia:',F7.2,' time in DROP:',F8.1,'s')
          WRITE(*,*)'**********leaving  DROP***************************'
      ENDIF
      RETURN
      END
      SUBROUTINE SUM_CONFS
*********************************************************************
* routine number 23 at level 70  in NILSSON-CRANKER                 *
* purpose: summing proton and neutron total quantities              *
* NOTE: this is nothing more than the simplest possible variant     *
* called by: MAIN                                                   *
* common blocks changed: TOTALS                                     *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      IF(IF_STRUT.GE.2)THEN
          CORRINERTIA=.5/DROPINERTIA-.5/STRUTINERTIA
      ELSE
          CORRINERTIA=0
      ENDIF
      MAXER=IANT
      MINER=1
      IF(NVAR.NE.3)THEN
          MINER=MAXER
      ENDIF
      IOTILL=0
      IF(IDEF_TYP.EQ.3)IOTILL=1
      DO 300 I=MINER,MAXER
         UT_TOT(I)=ETOTI(I,1)+ETOTI(I,2)+DROPENERGY
         UT_SPIN(I)=SPINI(I,1)+SPINI(I,2)
         IF(IF_STRUT.GE.2)UT_TOT(I)=UT_TOT(I)+
     &    CORRINERTIA*UT_SPIN(I)**2
         IF(IF_Q2.GE.1) UT_Q2(I)=Q2TOTI(I,1)
         IF(IF_Q2.GE.2) UT_Q2(I)=UT_Q2(I)+Q2TOTI(I,2)
         UT_OM(I)=OMEGA(I+IOTILL)
         IF(IDEF_TYP.GE.0)THEN
            UT_EPS(I)=EPS
            UT_GAM(I)=GAM
         ELSE
            CALL XYCOORD(EPS,GAM,X,Y,IX,IY)
            UT_EPS(I)=X
            UT_GAM(I)=Y
         ENDIF
         UT_EP4(I)=EPS4
         UT_CONF(I)=ICONFI(I,1)+ICONFI(I,2)/1000.
300   CONTINUE
      RETURN
      END
      SUBROUTINE SPROUT
*********************************************************************
* routine in NILSSON-CRANKER                                        *
* purpose: saving sp-values on file in GENITZ package format        *
* called by: SPDIA                                                  *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      DIMENSION KVANT(4)
      DATA KVANT/201,101,202,102/,RZERO/0.0/
      IPOST=(MOD(ITEST,2)+1)*IANT
      IF(NVAR.EQ.2)THEN
         IPOST=(MOD(ITEST,2)+1)*(IANT+(IE-1)*NGAM)
      ELSEIF(NVAR.EQ.3)THEN
         IPOST=(MOD(ITEST,2)+1)*(IANT+(ILOOP-1)*NOMEGA)
      ELSEIF(NVAR.EQ.4)THEN
         IPOST=(MOD(ITEST,2)+1)*(IANT+(IG-1)*NEPS4+(IE-1)*NGAM*NEPS4)
      ENDIF
      IF(MOD(ITEST,2).EQ.0)IPOST=IPOST+1
      IF(ISO1.NE.ISO2)IPOST=IPOST+ISO-1
      IPSST=2*(IPOST-2)+1
      IHEL=2*(NDIM(1,ISO)+NDIM(2,ISO))
      RMEV=1
      IF(IF_MEV.GT.0)THEN
          RMEV=RHOMEGA(ISO)*OMOM
      ENDIF
      WRITE(NSAV1FIL,REC=IPOST)EPS,GAM,EPS4,RZERO,OMEG/RMEV,OMOM,IHEL,
     &(NDIM(I,ISO),NDIM(I,ISO),I=1,2),(KVANT(I),I=1,4),ISO
      WRITE(NSAV2FIL,REC=IPSST)((E(I,ISYM,IANT+MOVE)/RMEV
     &,I=1,2*NDIM(ISYM,ISO)),ISYM=1,2)
      IPSST=IPSST+1
      WRITE(NSAV2FIL,REC=IPSST)((RJX(I,ISYM,IANT+MOVE),
     & I=1,2*NDIM(ISYM,ISO)),ISYM=1,2)
      RETURN
      END
      SUBROUTINE SPREAD
*********************************************************************
* routine in NILSSON-CRANKER                                        *
* purpose: reading sp-values from file in GENITZ package format     *
*     if failure - resetting to saving mode                         *
* called by: SPDIA                                                  *
* common blocks changed: ORBITALS                                   *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      DIMENSION ND(2)
      IERR=10
      RMEV=1
      IF(IF_MEV.GT.0)THEN
          RMEV=RHOMEGA(ISO)*OMOM
      ENDIF
      IPOST=(MOD(ITEST,2)+1)*IANT
      IF(NVAR.EQ.2)THEN
         IPOST=(MOD(ITEST,2)+1)*(IANT+(IE-1)*NGAM)
      ELSEIF(NVAR.EQ.3)THEN
         IPOST=(MOD(ITEST,2)+1)*(IANT+(ILOOP-1)*NOMEGA)
      ELSEIF(NVAR.EQ.4)THEN
         IPOST=(MOD(ITEST,2)+1)*(IANT+(IG-1)*NEPS4+(IE-1)*NGAM*NEPS4)
      ENDIF
      IF(MOD(ITEST,2).EQ.0)IPOST=IPOST+1
      IF(ISO1.NE.ISO2)IPOST=IPOST+ISO-1
      IPSST=2*(IPOST-2)+1
      IHEL=2*(NDIM(1,ISO)+NDIM(2,ISO))
      READ(NSAV1FIL,REC=IPOST,ERR=900)EPSA,GAMA,EPS4A,RZERO,
     & OMEGB,OMOMA,IHELA,(ND(I),ND(I),I=1,2),KV1,KV2,KV3,KV4,ISA
      IERR=9
      IF(ABS(EPS-EPSA).GT.0.001)GO TO 900
      IF(ABS(GAM-GAMA).GT.0.001)GO TO 900
      IF(ABS(EPS4-EPS4A).GT.0.001)GO TO 900
      IF(ABS(OMEGB-OMEG/RMEV).GT.0.001)GO TO 900
      IERR=8
      IF(IHELA.LT.IHEL)GO TO 900
      IERR=7
      IF(ISO.NE.ISA)GO TO 900
      IERR=6
      READ(NSAV2FIL,REC=IPSST,ERR=900)
     & ((E(I,ISYM,IANT+MOVE),I=1,2*ND(ISYM)),ISYM=1,2)
      IERR=5
      IPSST=IPSST+1
      READ(NSAV2FIL,REC=IPSST,ERR=900)
     & ((RJX(I,ISYM,IANT+MOVE),I=1,2*ND(ISYM)),ISYM=1,2)
      IF(IF_MEV.GT.0)THEN
          DO 1 ISYM=1,2
              DO 2 I=1,2*ND(ISYM)
                  E(I,ISYM,IANT+MOVE)=RMEV*E(I,ISYM,IANT+MOVE)
2             CONTINUE
1         CONTINUE
      ENDIF
      RETURN
900   CONTINUE
      WRITE(ERR_TEX(NERR+1),'(A,I2)')
     & ' saver file read error, type:',IERR
      CALL ERR_OUT(1)
      IF_LEVSAV=+1
CCCCC      IF_SAVDROP=+1
      WRITE(*,*)' shifting into saving mode!'
      RETURN
      END
      BLOCK DATA SETTER
      INCLUDE 'NICRAINC.FOR'
      DATA ISOTEXT(1)/'protons '/,ISOTEXT(2)/'neutrons'/
      DATA ANUC(14),ANUC(15),ANUC(16),ANUC(17),ANUC(18),ANUC(19)
     &/'Si','P ','S ','Cl','Ar','K '/
      DATA ANUC(20),ANUC(21),ANUC(22),ANUC(23)
     &/'Ca','Sc','Ti','V '/
      DATA ANUC(24),ANUC(25),ANUC(26),ANUC(27),ANUC(28),ANUC(29)
     &/'Cr','Mn','Fe','Co','Ni','Cu'/
      DATA ANUC(30),ANUC(31),ANUC(32),ANUC(33)
     &/'Zn','Ga','Ge','As'/
      DATA ANUC(34),ANUC(35),ANUC(36),ANUC(37),ANUC(38),ANUC(39)
     &/'Se','Br','Kr','Rb','Sr','Y '/
      DATA ANUC(40),ANUC(41),ANUC(42),ANUC(43)
     &/'Zr','Nb','Mo','Tc'/
      DATA ANUC(44),ANUC(45),ANUC(46),ANUC(47),ANUC(48),ANUC(49)
     &/'Ru','Rh','Pd','Ag','Cd','In'/
      DATA ANUC(50),ANUC(51),ANUC(52),ANUC(53)
     &/'Sn','Sb','Te','I '/
      DATA ANUC(54),ANUC(55),ANUC(56),ANUC(57),ANUC(58),ANUC(59)
     &/'Xe','Cs','Ba','La','Ce','Pr'/
      DATA ANUC(60),ANUC(61),ANUC(62),ANUC(63)
     &/'Nd','Pm','Sm','Eu'/
      DATA ANUC(64),ANUC(65),ANUC(66),ANUC(67),ANUC(68),ANUC(69)
     &/'Gd','Tb','Dy','Ho','Er','Tm'/
      DATA ANUC(70),ANUC(71),ANUC(72),ANUC(73)
     &/'Yb','Lu','Hf','Ta'/
      DATA ANUC(74),ANUC(75),ANUC(76),ANUC(77),ANUC(78),ANUC(79)
     &/'W ','Re','Os','Ir','Pt','Au'/
      DATA ANUC(80),ANUC(81),ANUC(82),ANUC(83)
     &/'Hg','Tl','Pb','Bi'/
      DATA ANUC(84),ANUC(85),ANUC(86),ANUC(87),ANUC(88),ANUC(89)
     &/'Po','At','Rn','Fr','Ra','Ac'/
      DATA ANUC(90),ANUC(91),ANUC(92),ANUC(93)
     &/'Th','Pa','U ','Np'/
      DATA ANUC(94),ANUC(95),ANUC(96),ANUC(97),ANUC(98),ANUC(99)
     &/'Pu','Am','Cm','Bk','Cf','Es'/
      DATA NSPFILE,NSAV1FIL,NSAV2FIL,NSAV3FIL,NTOTFILE/13,14,15,16,17/
      DATA IF_SPOUT/0/,IF_STRUT/2/,LEV_PRINT/21/,IF_LEVSAV/0/,
     &IF_ISOENEROUT/0/,IF_TOTEOUT/0/,MODELTYP/1/,ISYM1/1/,ISYM2/2/,
     &ISO1/1/,ISO2/2/,NVAR/3/,IF_MEV/1/,IF_SAVDROP/0/,IN_LEV/0/
     &SPERANGE/3/,IF_Q2/1/,NERR/0/,IF_NIVPRI/1/,IF_DROPCALC/1/
      DATA IPSAV/100/,IF_BCS/1/,G0,G1/19.2,7.4/,KOEFF_RANGE/15/
      DATA NUU,JCO/4,100/,ISTEP/0/
      DATA ((RKAPPA(I,J),RMY(I,J),J=1,2),I=1,15)
     &/ 0.12,0.0,  0.12,0.0,
     &  0.12,0.0,  0.12,0.0,
     &  0.105,0.0, 0.105,0.0,
     &  0.09,0.30, 0.09,0.25,
     &  0.065,0.57,0.070,0.39,
     &  0.060,0.65,0.062,0.43,
     &  0.054,0.69,0.062,0.34,
     &  0.054,0.69,0.062,0.26,
     &  0.054,0.69,0.062,0.26,
     &  0.054,0.69,0.062,0.26,
     &  0.054,0.69,0.062,0.26,
     &  0.054,0.69,0.062,0.26,
     &  0.054,0.69,0.062,0.26,
     &  0.054,0.69,0.062,0.26,
     &  0.000,0.00,0.000,0.00/
*
      END
      SUBROUTINE BCS(E1,Q,Z,A,D,RL,DEP,HOM,G,DNDL)
*********************************************************************
* routine number 42  at level 70  in NILSSON-CRANKER                *
* purpose: perfoming a BCS calculation at omega=0                   *
* called by: PREDIA  calling: BARDEE                                *
* no external common blocks affected                                *
*********************************************************************
      INCLUDE 'NICRAINC.FOR'
      DIMENSION E1(2*MAXD2),V2B(MAXD2),EB(MAXD2)
      INTEGER Q,Z,A
      COMMON /AREA/K1,K2,BDNDL
      TID1=SECOND(U)
      D=0.7
      NQ=Q/2
      I=0
      DO 1 II=1,NQ
          I=I+1
          EB(I)=HOM*E1(2*II-1)
1     CONTINUE
      NZ=Z/2
      KZ=2*NZ
      RL=(EB(NZ)+EB(NZ+1))/2.
      IF(MOD(Z,2).EQ.1)RL=EB(NZ+1)
      CALL BARDEE(Z,EB,G,0.001,D,RL,V2B,IFLAG,KOEFF_RANGE)
      D=SQRT(D)
      IF(IFLAG.NE.0)THEN
          WRITE(ERR_TEX(NERR+1),1800)IE,IG,I4
1800      FORMAT('BCS calculation error at Ieps,Igam,Ieps4=',3I3)
          CALL ERR_OUT(1)
      ENDIF
      DNDL=BDNDL
      SIP=0.
      IF(IFLAG.NE.0) GO TO 999
      DO 2 I=K1,K2
 2        SIP=SIP+(2*EB(I)-G*V2B(I))*V2B(I)
      SIP=SIP-D**2/G
      KM1=K1-1
      DO 3 I=1,KM1
 3        SIP=SIP+2*EB(I)
      SIP=SIP-G*KM1
 999  SUPP=0.
      DO 4 I=1,NZ
 4        SUPP=SUPP+2*EB(I)
      SUP=SUPP-G*NZ
      IF(MOD(Z,2).EQ.1)SUP=SUP+E1(IODD)*HOM
      SIPP=SIP+G*NZ
      DEP=SIP-SUP
      IF(IFLAG.NE.0)DEP=0.
      IF(LEV_PRINT.EQ.42.OR.LEV_PRINT.GE.70)THEN
        WRITE(*,*)'**************BCS output*********************'
         WRITE(*,*)' orbitals:'
          WRITE(*,1962)'  #  ',(2*I,I=NZ-5,NZ+5)
          WRITE(*,1963)' e sp',(EB(I),
     &     I=NZ-5,NZ+5)
          WRITE(*,1963)' V**2 ',(V2B(I),
     &    I=NZ-5,NZ+5)
1962      FORMAT(A,11I10)
1963      FORMAT(A,11F10.3)
          IF(MOD(Z,2).EQ.1)WRITE(*,'(A,I3,F10.3)')' odd # & energy',
     &IODD,E1(IODD)
         WRITE(*,100)Z,A,G,HOM,D,RL,SUPP,SIPP,DEP,IFLAG
  100    FORMAT(' Z=',I3,' A=',I3,' G=',F7.4,' hw=',F7.4/
     &   ' delta,lambda=',2F8.4,' sum,Ebcs=',2F10.4,' dE=',F7.4/
     &   ' flag=',I3)
        TID2= SECOND(U)
        WRITE(*,1961)' time in BCS:',TID2-TID1
1961    FORMAT(A,F6.3,'s')
        WRITE(*,*)'I-----------------------------------------------I'
      ENDIF
      D=D/HOM
      RL=RL/HOM
      DEP=DEP/HOM
      RETURN
      END
      SUBROUTINE HELP
* service routine for CLEBI and MRHOT
      COMMON /PTR/ F(57)
      COMMON /GM/ GM(150)
      COMMON /CMN1/L2F(18)
C*****Computation of f(n+1)=ln(n!)
      F(1) = 0.
      DO 101 I=1,56
      AII = FLOAT(I)
      F(I+1) = F(I) + ALOG(AII)
  101 CONTINUE
C*****Computation of (l2)av = l2f(n+1)
      L2F(1) = 0
      DO 106 I = 1,17
      L2F(I+1) = L2F(I) + I + 1
  106 CONTINUE
C*****Computation of gamma-function in integer and half-integer points.
      GM(52)=1.0
      GM(51)=1.77245385
      DO 34 I=1,28
      GM(2*I+52)=I*GM(2*I+50)
   34 GM(2*I+51)=(I-0.5)*GM(2*I+49)
      GM(49)=-2.0*GM(51)
      DO 35 I=1,20
      IR=-2*I+49
      IIR=IR+2
   35 GM(IR)=GM(IIR)/(-I-0.5)
      RETURN
      END
      SUBROUTINE OMO(OMOM)
* calculates volume conservation factor OMOM
      INCLUDE 'NICRAPAR.FOR'
      COMMON /HAMPAR/EPS,GAM,EPS4,OMEG,DELT,RLAM,NDIM(MSYM,MISO),ISO
      DIMENSION Z(9),WZ(9)
      EXTERNAL FVCON
      IF(ISSE.NE.1267)CALL LGAUSS(Z,WZ,9)
      ISSE=1267
      EPS2=EPS*COS(GAM)
      EPS22=EPS*SIN(GAM)/SQRT(0.75)
      IF(ABS(EPS4).LT.0.0001) GO TO 4
      PI=3.1415926536
      NX=2
      NY=3
      A=BVCON(FVCON,-1.0,1.0,0.0,2.*PI,9,9,Z,WZ,Z,WZ,NX,NY)
      OMOM=A**(1.0/3.0)
      GO TO 5
    4 OMOM=((1.+EPS2/3.+EPS22/2.)*(1.+EPS2/3.-EPS22/2.)*(1.-2.*EPS2/3.))
     $**(-1.0/3.0)
    5 CONTINUE
      GO TO 6
    6 CONTINUE
      RETURN
      END
      FUNCTION FVCON(XX,YY)
* subroutine for OMO
      INCLUDE 'NICRAPAR.FOR'
      COMMON /HAMPAR/EPSP,GAM,EPS4,OMEG,DELT,RLAM,NDIM(MSYM,MISO),ISO
      COMMON /FYRCOM/ D40,D42,D44,EPS2,EPS22,OMX,OMY,OMZ,FAC,HAC
      X=XX
      Y=YY
      PI=3.1415926536
      ANM=1.0-EPS2/3.0*(3.*X*X-1.)+0.5*EPS22*(1.-X*X)*COS(2.*Y)+
     $EPS4/4.0*D40*
     $(35.*X*X*X*X-30.*X*X+3.)+5.*EPS4*D42/SQRT(10.)*(7.*X*X-1.)
     $*(1.-X*X)*COS(2.*Y)+35.*EPS4*D44/2./SQRT(70.)*(1.-X*X)*
     $(1.-X*X)*COS(4.*Y)
      FVCON=1./(4.*PI)/SQRT((1.+EPS2/3.+EPS22/2.)*(1.+EPS2/3.-EPS22/2.)*
     $(1.-2.*EPS2/3.))/(ANM**(1.5))
      RETURN
      END
*ELEMENT LGAUSS
      SUBROUTINE LGAUSS(X,W,N)
C     subroutine lgauss computes the zeroes of the legendre polynomial
C     and their associated weights for a gaussian quadrature
      DIMENSION X(1),W(1)
      IF(N-1) 1,2,3
C     request for a zero point formula is meaningless so  stop
   1  WRITE(6,20)
   20 FORMAT(1X,34HA ZERO POINT FORMULA WAS REQUESTED)
      STOP
  2   X(1)=0.0
      W(1)=2.0
      RETURN
C     for a one point formula send back results without computing
  3   R=N
      G=-1.
C     the initial guess for the smallest root of p(n) is taken as -1.
      DO 147 I=1,N
      TEST=-2.
      IC=N+1-I
C     whenever we find a root of the polynomial, its negative is also a
C     the index ic tells where to store the other root
      IF(IC.LT.I) GO TO 150
  4   S=G
      T=1.
      U=1.
      V=0.
C     evaluation of the n-th legendre polynomial and its first derivativ
C     where   U=DS/DX
C             V=DT/DX
C             DP=DP/DX
      DO 50 K=2,N
      A=K
      P=((2.0*A-1.0)*S*G-(A-1.0)*T)/A
      DP=((2.0*A-1.0)*(S+G*U)-(A-1.0)*V)/A
      V=U
      U=DP
      T=S
  50  S=P
      IF(ABS(TEST+G).LT.00.00001)GO TO 100
      IF (ABS((TEST-G)/(TEST+G)).LT.0.0000005) GO TO 100
      SUM=0.
      IF(I.EQ.1) GO TO 52
C     the following computes the reduced legendre polynomial and its
C     derivative:
      DO 51 K=2,I
  51  SUM=SUM+1./(G-X(K-1))
  52  TEST=G
      G=G-P/(DP-P*SUM)
      GO TO 4
  100 X(IC)=-G
      X(I)=G
      W(I)=2./(R*T*DP)
      W(IC)=W(I)
  147 G=G-R*T/((R+2.)*G*DP+R*V-2.*R*T*SUM)
  150 RETURN
      END
*ELEMENT BVCON
      FUNCTION BVCON(FVCON,XA,XB,YA,YB,NX,NY,X,WX,Y,WY,IX,IY)
C     double integration subroutine for volume conservation condition.
C     fvcon   is the external function to be integrated.
C     X,WX    are the arrays of normalized points and weights for x-vari
C     Y,WY    are the arrays of normalized points and weights for y-vari
      DIMENSION X(1),WX(1),Y(1),WY(1)
C     calculation of the denormalizing x and y slopes
      HX=(XB-XA)/(2.*IX)
      HY=(YB-YA)/(2.*IY)
      BVCON=0.0
      DO 10 L=1,IY
C     calculation of the denormalizing y intercept
      BY=YA+(2.*L-1.)*HY
      DO 10 K=1,IX
C     calculation of the denormalizing x intercept
      BX=XA+(2.*K-1.)*HX
      DO 10 J=1,NY
C     calculation of the denormalized point y
      A=HY*Y(J)+BY
      DO 10 I=1,NX
C     calculation of the denormalized point x and summation of the
C     approximating sum
   10 BVCON=BVCON+WY(J)*WX(I)*FVCON(HX*X(I)+BX,A)
C     adjustment of the calculated sum.
      BVCON=HX*HY*BVCON
      RETURN
      END
      FUNCTION SECOND(XX)
* gives cpu time used in seconds, can return 0 if cpu-time unobtainable
* this one is for VAX VMS systems
      SECOND=0.0
C     LOGICAL STATUS
C     DATA I0,IPID,JPI_CPU/0,0,1031/
C     IF(IPID.EQ.0)THEN
C         STATUS = LIB$GETJPI(I0,IPID,,ISEX,,)
C         IF(.NOT.STATUS)THEN
C             WRITE(*,*)' something wrong with lib$getjpi!'
C         ENDIF
C     ENDIF
C     STATUS = LIB$GETJPI(JPI_CPU,IPID,,ISEX,,)
C     IF(.NOT.STATUS)THEN
C         WRITE(*,*)' something wrong with lib$getjpi!'
C     ENDIF
C     SECOND=FLOAT(ISEX)/100.
      RETURN
      END
      SUBROUTINE SORT3R(E,KV1,KV2,NI)
* quicksort for 3 entities
      DIMENSION E(1),KV1(1),KV2(1)
      REAL KV1
      INTEGER Y2
      N=NI
      M=2**INT(1.4427*ALOG(FLOAT(N-1)))-1
    1 IF(M.LE.0) RETURN
      I=1
    2 J=I
      K=I+M
      Y=E(K)
      Y1=KV1(K)
      Y2=KV2(K)
    3 IF(Y.LT.E(J)) GO TO 5
    4 L=J+M
      E(L)=Y
      KV1(L)=Y1
      KV2(L)=Y2
      I=I+1
      IF(I+M.LE.N) GO TO 2
      M=(M-1)/2
      GO TO 1
    5 L=J+M
      E(L)=E(J)
      KV1(L)=KV1(J)
      KV2(L)=KV2(J)
      J=J-M
      IF(J.GE.1) GO TO 3
      GO TO 4
      END
      FUNCTION CLEBI(I1,I2,I3,N1,N2,N3)
* Clepsh-Gordan coefficients, all input parameters have twice their
*  physical value
      COMMON /PTR/ FCT(57)
      INTEGER Z,ZMIN,ZMAX
      J1=I1
      J2=I2
      J =I3
      N=57
      M1=N1
      M2=N2
      M=-N3
      CC=0.0
      JSUM=J1+J2+J
      JM1 =J1-IABS(M1)
      JM2 =J2-IABS(M2)
      JM3 =J -IABS(M )
      IF((MOD(JSUM,2).NE.0).OR.(MOD(JM1,2).NE.0).OR.(MOD(JM2,2).NE.0)
     1.OR.(MOD(JM3,2).NE.0).OR.(JM1.LT.0).OR.(JM2.LT.0).OR.(JM3.LT.0))
     2GO TO 1
      IF((M1+M2+M.NE.0).OR.(J.GT.J1+J2).OR.(J.LT.IABS(J1-J2))) GO TO 1
      ZMIN=0
      IF(J-J2+M1.LT.0) ZMIN=-J+J2-M1
      IF(J-J1-M2+ZMIN.LT.0) ZMIN=-J+J1+M2
      ZMAX=J1+J2-J
      IF(J2+M2-ZMAX.LT.0) ZMAX=J2+M2
      IF(J1-M1-ZMAX.LT.0) ZMAX=J1-M1
      JA=(J1+M1)/2+1
      JB=JA-M1
      JC=(J2+M2)/2+1
      JD=JC-M2
      JE=(J +M )/2+1
      JF=JE-M
      JG=(J1+J2-J)/2+1
      JH=JA+JB-JG
      JI=JC+JD-JG
      JJ=JE+JF+JG-1
      IF(JJ.GT.N) GO TO 5
      IA=ZMIN/2
      IB=JG-IA+1
      IC=JB-IA+1
      ID=JC-IA+1
      IE=JA-JG+IA
      IF=JD-JG+IA
      FASE=1.0
      IF(MOD(IA,2).EQ.0) FASE=-FASE
      Z =ZMIN
      ARII=FLOAT(J+1)
      ARI=(FCT(JA)+FCT(JB)+FCT(JC)+FCT(JD)+FCT(JE)+FCT(JF)+FCT(JG)
     1 +FCT(JH)+FCT(JI)-FCT(JJ)+ALOG(ARII))/2.0
    2 IA=IA+1
      IB=IB-1
      IC=IC-1
      ID=ID-1
      IE=IE+1
      IF=IF+1
      FASE=-FASE
      ARP=-FCT(IA)-FCT(IB)-FCT(IC)-FCT(ID)-FCT(IE)-FCT(IF)
      ARS=ARI+ARP
      IF(ARS.GT.80.0) GO TO 5
      IF(ARS.LT.-80.0) GO TO 10
      CC=CC+FASE*EXP(ARS)
   10 CONTINUE
      Z=Z+2
      IF(Z.LE.ZMAX) GO TO 2
    1 CLEBI=CC
      RETURN
C   3 WRITE(*,4)
C   4 FORMAT(26H0 ERROR IN ARGUMENT OF VCC)
C     GO TO 7
    5 WRITE(*,6)
      STOP 3
    6 FORMAT(29H0 ERROR - FACTORIALS EXCEEDED)
C   7 WRITE(61,8)I1,I2,I3,N1,N2,N3
C   8 FORMAT(10I10)
C     RETURN
      END
      SUBROUTINE SORT(E,KV1,NI)
* quicksort for 2 entities
      DIMENSION E(1),KV1(1)
      INTEGER Y1
      N=NI
      M=2**INT(1.4427*ALOG(FLOAT(N-1)))-1
    1 IF(M.LE.0) RETURN
      I=1
    2 J=I
      K=I+M
      Y=E(K)
      Y1=KV1(K)
    3 IF(Y.LT.E(J)) GO TO 5
    4 L=J+M
      E(L)=Y
      KV1(L)=Y1
      I=I+1
      IF(I+M.LE.N) GO TO 2
      M=(M-1)/2
      GO TO 1
    5 L=J+M
      E(L)=E(J)
      KV1(L)=KV1(J)
      J=J-M
      IF(J.GE.1) GO TO 3
      GO TO 4
      END
      SUBROUTINE MRHOT(NNP,LLP,LAMP,NN,LL,LAM,T,ME)
* matrix elements ME=<NNP LLP LAMP/ r**T /NN LL LAM>
      COMMON /LU/ LU,MAXBL
      COMMON /GM/ GM(150)
      REAL ME
      INTEGER T,TA,E,C,D,SLUT
    1 FORMAT(12I6)
  300 FORMAT(7H DANGER)
  301 FORMAT(31H TERMINATION SOMETHING IS WRONG)
  302 FORMAT(33H IS NOT THIS A WONDERFUL COMPILER)
      NANP=NNP
      LALP=LLP
      LAAMP=LAMP
      NAN=NN
      LAL=LL
      LAAM=LAM
      TA=T
      IF(LAAM-LAAMP) 231,201,231
  201 CONTINUE
      NP=(NANP-LALP)/2+1
      N=(NAN-LAL)/2+1
      IF(NP-N)202,203,203
  202 NPNP=NP
      NP=N
      N=NPNP
      L=LALP
      LALP=LAL
      LAL=L
  203 CONTINUE
      C=(LAL-LALP-TA)/2
      D=(LALP-LAL-TA)/2
      E=D+NP-N-1
      IF(2*C.EQ.(LAL-LALP-TA)) GO TO 401
      ME=0.0
      DO 3 K=1,N
      K1=2*(N-K+1)+50
      K2=LAL-LALP-TA+2*(K-1)+50
      K3=LALP-LAL-TA+2*(K+NP-N-1)+50
      K4=2*(NP-N+K)+50
      K5=2*(N-K+1)+LAL+LALP+1+TA+50
      K6=LAL-LALP-TA+50
      K7=LALP-LAL+50-TA
    3 ME=ME+GM(2*N+50)*GM(2*NP+50)/GM(K1)*GM(K2)/GM(2*K+50)*GM(K3)
     */GM(K4)*GM(K5)/GM(K6)/GM(K7)
      K1=2*N+2*LAL+51
      K2=2*NP+2*LALP+51
      ME=ME/SQRT(GM(2*N+50)*GM(K1)*GM(2*NP+50)*GM(K2))
      IF(MOD(ABS(NP-N),2).NE.0)ME=-ME
      GO TO 230
  401 CONTINUE
      IF(C.LE.0.AND.D.LE.0)  GO TO 204
      GO TO 208
  204 CONTINUE
      IF(E.GT.-1) GO TO 206
      GO TO 209
  206 SLUT=0
      GO TO  220
  209 CONTINUE
      IF(C.LE.E+1) GO TO 207
      GO TO 221
  207 SLUT=-E
      GO TO 220
  221 SLUT=1-C
      GO TO 210
  208 CONTINUE
      IF(C.LE.0.AND.D.GT.0) GO TO  211
      GO TO 212
  211 SLUT=1-C
      GO TO 213
  212 CONTINUE
      IF(C.GT.0.AND.D.LE.0) GO TO 214
      GO TO 215
  214 CONTINUE
      IF(E.GT.-1) GO TO 216
      GO TO 217
  216 SLUT=0
      GO TO 218
  217 SLUT=-E
  218 CONTINUE
      GO TO 219
  215 CONTINUE
      WRITE(LU,300)
      WRITE(LU,301)
      WRITE(LU,302)
      STOP
  219 CONTINUE
  210 CONTINUE
  220 CONTINUE
  213 CONTINUE
      IF(SLUT.EQ.0) GO TO 231
      IF(SLUT.LT.0) GO TO 222
      GO TO 223
  222 WRITE(LU,300)
      WRITE(LU,301)
      WRITE(LU,302)
      STOP
  223 ME=0.0
      IF(SLUT.GT.N)  SLUT=N
      DO 234 K=1,SLUT
      KK1=2*(N-K+1)+50
      KK2=2*(N-K+1)+LAL+LALP+1+TA+50
      KK3=2*(NP-N+K)+50
      KK4=2-2*C+50
      KK5=2*(2-K-C)+50
      KK6=2*(C+K-1)+50
      KK7=2*C+50
      IF(C.LE.0) GO TO 224
      GO TO 225
  224 GLI=(-1)**(K+1)*GM(KK4)/GM(KK5)
      GO TO 226
  225 GLI=GM(KK6)/GM(KK7)
  226 CONTINUE
      KK8=2-2*D+50
      KK9=2*(2-D-K-NP+N)+50
      KK10=2*(D+K+NP-N-1)+50
      KK11=2*D+50
      IF(D.LE.0.) GO TO 227
      GO TO 228
  227 GLA=GM(KK8)/GM(KK9)*(-1)**(1+K+NP-N)
      GO TO 229
  228 GLA=GM(KK10)/GM(KK11)
  229 CONTINUE
      ME=ME+GM(2*N+50)/GM(2*K+50)*GM(2*NP+50)/GM(KK1)*GM(KK2)/GM(KK3)*
     *GLI*GLA
  234 CONTINUE
      KK12=2*(N+LAL)+51
      KK13=2*(NP+LALP)+51
      ME=ME/SQRT(GM(2*N+50)*GM(2*NP+50)*GM(KK12)*GM(KK13))
      IF(MOD(ABS(NP-N),2).NE.0)ME=-ME
      GO TO 230
  231 ME=0.0
  230 CONTINUE
      RETURN
      END
      SUBROUTINE MATDIA(DMAT,EIG,IDIM,VEC)
      INCLUDE 'NICRAPAR.FOR'
* driver for diagonaization routine
* DMAT matrix to diagonalize, IDIM actual size of problem
* EIG & VEC eigenvalues & vectors
* NOTE!!! eigenvalues &-vectors assumed to be in decending
*  energy order!!!!!!
      DIMENSION DMAT(MAXDIM,MAXDIM),EIG(MAXDIM),
     &VEC(MAXDIM,MAXDIM),VW(MAXDIM,6)
      MUDIM=MAXDIM
      IF(IDIM.LT.1)RETURN
*
      CALL DNICES(IDIM,IDIM,IDIM,DMAT,MUDIM,1.E-6,
     & IDIM,EIG,VEC,MUDIM,ICON,
     & VW(1,1),VW(1,2),VW(1,3),VW(1,4),VW(1,5),VW(1,6))
*
      IF(ICON.NE.0)THEN
              CALL ERR_SERV(' diagonalization error')
      ENDIF
      RETURN
      END
      SUBROUTINE CHECK_READ(N1,N2,*)
* input unit checker
      IF(N1.LE.0)RETURN 1
      IF(N2.EQ.0)THEN
          IERR=0
          NERR=0
          RETURN
      ELSEIF(N2.LT.0)THEN
          IF(IERR.GE.-N2)STOP 'Too many input errors'
          RETURN
      ENDIF
      IF(NERR.NE.0)THEN
          IF(N2.EQ.NERR)THEN
              IF(N1.NE.N2)THEN
                  IF(N1.LT.N2)THEN
                      RETURN 1
                  ELSE
                      IERR=IERR+1
                      NERR=0
                      RETURN
                  ENDIF
              ELSE
                  NERR=0
                  WRITE(*,2)N2
                  RETURN
              ENDIF
          ELSE
              NERR=0
          ENDIF
      ENDIF
      IF(NERR.EQ.0)THEN
          IF(N1.NE.N2)THEN
              WRITE(*,1)N1,N2
              IF(N1.LT.N2)THEN
                  NERR=N2
                  RETURN 1
              ELSE
                  IERR=IERR+1
                  RETURN
              ENDIF
          ELSE
              RETURN
          ENDIF
      ENDIF
1     FORMAT(' input unit:',I4,' found when unit:',I4,' expected')
2     FORMAT(' but no worry, now unit:',I4,' is found')
      END
      SUBROUTINE ERR_OUT(NNN)
      INCLUDE 'NICRAINC.FOR'
* error handler
      IF(NNN.EQ.1)THEN
          WRITE(*,*)' **** ERROR ****'
          WRITE(*,10)ERR_TEX(NERR+1)
          IF(NERR+1.LT.MERR)NERR=NERR+1
          IF(NERR.EQ.1)KERR=0
          KERR=KERR+1
          RETURN
      ENDIF
      IF(NNN.EQ.0)THEN
          IF(NERR.EQ.0)THEN
              WRITE(*,*)'No logged errors has occured'
              RETURN
          ENDIF
          WRITE(*,*)
          WRITE(*,*)'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
          WRITE(*,*)'X  THERE HAS BEEN SOME ERRORS HERE     YYYYYYYYYYY'
          WRITE(*,*)'ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ'
          IF(NERR.EQ.MERR-1)WRITE(*,*)' error buffer became filled'
          WRITE(*,*)' THE FOLLOWING ERROR MESSAGES HAS BEEN LOGGED:'
          DO 1 I=1,MIN(KERR,MERR)
              WRITE(*,10)ERR_TEX(I)
1         CONTINUE
10        FORMAT(1X,A)
      ENDIF
      RETURN
      END
      SUBROUTINE ERR_SERV(ERRLINE)
      INCLUDE 'NICRAINC.FOR'
* error handler service
      CHARACTER ERRLINE*(*)
      ILEN=LEN(ERRLINE)
      ERR_TEX(NERR+1)=ERRLINE(1:ILEN)
      CALL ERR_OUT(1)
      RETURN
      END
      SUBROUTINE STRUTS(E,RM,OMEGA,K,NN,GAMMA,ES,AIS,
     &SUM,EF,G1,G2,SPIN)
* strutinsky calculation
C  E=array of single particle energies
C  RM=array of spin components
C  OMEGA=rotation frequency
C  K=number of calculated energies in e
C  NN=nucleon number
C  GAMMA=shell smearing parameter
C  ES=smoothed energy
C  SUM=sum of single particle energies
C  AIS=smoothed spin
C  SPIN=angular momentum for tilted fermi surface
C  not more than 10 iterations are allowed
C  EF=fermi level
C  G1=level density at fermi level
C  G2=spin  density at fermi level
C
      DIMENSION E(K),RM(K)
      DATA  SIPI/.5641896/
      DATA XMIN,XMAX/3.4,3.5/
      T(X)=(((((4.30638E-5*X+2.765672E-4)*X+1.520143E-4)*X+9.2705272E-3)
     1*X+4.22820123E-2)*X+7.05230784E-2)*X
      TT(Z)=1./(T(Z)+1.)
      ERF(Y)=1.-(TT(Y)*TT(Y)*TT(Y)*TT(Y)*TT(Y)*TT(Y)*TT(Y)*TT(Y)*TT(Y)*
     *TT(Y)*TT(Y)*TT(Y)*TT(Y)*TT(Y)*TT(Y)*TT(Y))
      N=NN
      EF=0.5*(E(NN)+E(NN+1))
      IT=0
 100  CONTINUE
      IF(IT-20)111,111,110
 110  PRINT 9,K,N
 9    FORMAT(' ****** ERROR IN STRUTS** K=',I4,' N=',I4)
      STOP
 111  CONTINUE
      AN1=0.
      EN1=0.
      AI1=0.
      DO 50 I=1,K
      X=(EF-E(I))/GAMMA
      IF(X.LE.XMAX) GO TO 50
      AN1=AN1+1.0
      EN1=EN1+E(I)+OMEGA*RM(I)
      AI1=AI1+RM(I)
 50   CONTINUE
      DN1=0.
      AN=0.
      EN=0.
      AI=0.
      G1=0.
      G2=0.
      DO 200 I=1,K
      X=(EF-E(I))/GAMMA
      IF(X)1,2,3
 1    IF(-X.GT.XMIN)GO TO 200
      PHI=-ERF(-X)
      GO TO 4
 2    PHI=0.
      GO TO 4
 3    IF(X.GT.XMAX) GO TO 200
      PHI=ERF(X)
 4    CONTINUE
      X2=X*X
      EX=SIPI*EXP(-X2)
      EXG=EX/GAMMA
      ANY=(1.+PHI)/2.+(1./12.*X2*X2-2./3.*X2+19./16.)*X*EX
C4    ANY=(1.+PHI)/2.-(1./4.*X2-7./8.)*X*EX
      ENY=(1./12.*X2*X2*X2-5./8.*X2*X2+15./16.*X2-5./32.)*EX*GAMMA+
     &E(I)*ANY
C4    ENY=(-1./4.*X2*X2+3./4.*X2-3./16.)*EX*GAMMA+E(I)*ANY
      AIY=ANY*RM(I)
      DNY=(35./16.-35./8.*X2+7./4.*X2*X2-1./6.*X2*X2*X2)*EXG
C4    DNY=(15./8.-5./2.*X2+1./2.*X2*X2)*EXG
      DN1=DN1+DNY
      AN=AN+ANY
      EN=EN+ENY
      AI=AI+AIY
      G1=G1+DNY
      G2=G2+DNY*RM(I)
 200  CONTINUE
      AIS=AI+AI1
      ES=EN+EN1+OMEGA*AI
      AS=AN+AN1
      IF(IODD.EQ.0) GO TO 202
C
C         odd particles treated as even
C
 202  CONTINUE
      DNN=FLOAT(N)-AS
      DF=DNN/DN1
      IF(ABS(DNN).LT.1.E-3) GO TO 300
      EF=EF+DF
      IT=IT+1
      GO TO 100
 300  CONTINUE
 1111 CONTINUE
      EF=E(NN)+0.001
      N=NN
      SUM=0.
      SPIN=0.
      IJ=0.
      DO 301 I=1,N
      SPIN=SPIN+RM(I)
      SUM=SUM+E(I)+OMEGA*RM(I)
      X=(EF-E(I))/GAMMA
      IF(X.LE.0.0)IJ=IJ+1
  301 CONTINUE
      RETURN
      END
      SUBROUTINE BSC24G(EPS,GAMMA,EPS4,OTMO,BS,BC,NTH,NFI)
* calculation of deformed liquid drop surface and
* Columb corrections BS and BC
C
C***** indata: EPS,GAMMA,EPS4 - deformation
C              NTH, NFI - number of lattice points per interval pi/2
C                         TH-wise and FI-wise respectively. limitations
C                         are in array dimensions below and in abscwt.
C                         recommended values 8 and 8.
C     EXTERNAL: NRDS,ABSCWT
C
      DOUBLE PRECISION EPS2,EPS22,EPS4D,DS,R,SINTH,AR(3),AN(3),PI,PIQ,
     $   ATH(16),WTH(16),AFI(16),WFI(16),OMOMD,BSD,TH,FI,WT,SP1,SP2
     $,RVEKT(16,16,3),SNVEKT(16,16,3),BCD,SUMBCD,AVST2,R2MR1(3),AVST,W
     $,RNORM,SNORM
      COMMON /ECOM/ EPS2,EPS22,EPS4D,OMOM,AR,AN,DS,R,SINTH
      DATA PI /3.14159265358979324D0/
      EPS2=DBLE(EPS)*DCOS(DBLE(GAMMA)*PI/180.D0)
      EPS22=DBLE(EPS)*DSIN(DBLE(GAMMA)*PI/180.D0)/DSQRT(.75D0)
      EPS4D=DBLE(EPS4)
      PIQ=PI/4.D0
      CALL ABSCWT(ATH,WTH,NTH)
      CALL ABSCWT(AFI,WFI,NFI)
      OMOMD=0.D0
      BSD=0.D0
      DO 1 L5=1,NTH
      TH=PIQ*(ATH(L5)+1.D0)
      DO 1 L6=1,NFI
      FI=PIQ*(AFI(L6)+1.D0)
      CALL NRDS(TH,FI)
      SP1=0.D0
      SP2=0.D0
      DO 2 L9=1,3
      SP1=SP1+AR(L9)*AN(L9)
      SP2=SP2+AN(L9)**2
      RVEKT(L5,L6,L9)=AR(L9)
      SNVEKT(L5,L6,L9)=AN(L9)
    2 CONTINUE
      WT=WTH(L5)*WFI(L6)
      OMOMD=OMOMD+SP1*WT
      BSD=BSD+DSQRT(SP2)*WT
    1 CONTINUE
      OMOMD=OMOMD*8.D0*PIQ**2
      OMOMD=DCBRT(OMOMD/4.D0/PI)
      BSD=BSD*8.D0*PIQ**2
      BSD=BSD/4.D0/PI/OMOMD**2
      OMOM=OMOMD
      BS=BSD
      OTMO=OMOM
      RNORM=OMOMD
      SNORM=OMOMD**2
      DO 3 L5=1,NTH
      DO 3 L6=1,NFI
      DO 3 L9=1,3
      RVEKT(L5,L6,L9)=RVEKT(L5,L6,L9)/RNORM
      SNVEKT(L5,L6,L9)=SNVEKT(L5,L6,L9)/SNORM
 3    CONTINUE
      BCD=0.0
      DO 5 L5=1,NTH
      DO 5 L6=1,NFI
      W=WTH(L5)*WFI(L6)
      DO 6 L9=1,3
      AR(L9)=RVEKT(L5,L6,L9)
      AN(L9)=SNVEKT(L5,L6,L9)
 6    CONTINUE
      DO 5 L7=1,NTH
      DO 5 L8=1,NFI
      SUMBCD=0.D0
      ASSIGN 11 TO LABEL
      GO TO 7
 11   CONTINUE
      AR(3)=-AR(3)
      AN(3)=-AN(3)
      ASSIGN 12 TO LABEL
      GO TO 7
 12   CONTINUE
      AR(1)=-AR(1)
      AN(1)=-AN(1)
      ASSIGN 13 TO LABEL
      GO TO 7
 13   CONTINUE
      AR(3)=-AR(3)
      AN(3)=-AN(3)
      ASSIGN 14 TO LABEL
      GO TO 7
 14   CONTINUE
      AR(2)=-AR(2)
      AN(2)=-AN(2)
      ASSIGN 15 TO LABEL
      GO TO 7
 15   CONTINUE
      AR(3)=-AR(3)
      AN(3)=-AN(3)
      ASSIGN 16 TO LABEL
      GO TO 7
 16   CONTINUE
      AR(1)=-AR(1)
      AN(1)=-AN(1)
      ASSIGN 17 TO LABEL
      GO TO 7
 17   CONTINUE
      AR(3)=-AR(3)
      AN(3)=-AN(3)
      ASSIGN 18 TO LABEL
 7    AVST2=0.D0
      DO 8 L9=1,3
      R2MR1(L9)=AR(L9)-RVEKT(L7,L8,L9)
      AVST2=AVST2+R2MR1(L9)**2
 8    CONTINUE
      AVST=DSQRT(AVST2)
      IF(AVST.LT.1.D-15) GO TO 71
      SP1=0.D0
      SP2=0.D0
      DO 9 L9=1,3
      SP1=SP1+R2MR1(L9)*AN(L9)
      SP2=SP2+R2MR1(L9)*SNVEKT(L7,L8,L9)
 9    CONTINUE
      TH=SP1*SP2
      SUMBCD=SUMBCD+TH/AVST
 71   CONTINUE
      GO TO LABEL
 18   CONTINUE
      WT=W*WTH(L7)*WFI(L8)
      BCD=BCD-SUMBCD*WT
 5    CONTINUE
      BCD=BCD*PIQ**4*8.D0
      BCD=BCD*5.D0/64.D0/PI**2
      BC=BCD
      RETURN
      END
      SUBROUTINE NRDS(TH,FI)
C
C        the radius vector AR(3), surface normal vector AN(3)  and their
C        lengths R and DS at point TH, FI of a surface with deformation
C        EPS, EPS22, EPS4. for volume normalization AR, R and AN, ds
C        must be divided by OMOM and OMOM2 respectively.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL OMOM
      COMMON /ECOM/ EPS2,EPS22,EPS4,OMOM,AR(3),AN(3),DS,R,SINTH
      COMMON /FYRFAC/DA,DB,DC
C EPS4-POTENTIAL  EPS4*(DA*Y40 + DB*(Y42+Y4-2) +DC*(Y44+Y4-4))
      EA=EPS2/3.D0
      EB=EPS22/2.D0
      EC=EA+1.D0
      ED=1.D0-EA-EA
      AAAA=DSQRT(10.D0)/2.D0
      BBBB=DSQRT(70.D0)/4.D0
      COSTH=DCOS(TH)
      COS2TH=COSTH**2
      SIN2TH=1.D0-COS2TH
      SINTH=DSQRT(SIN2TH)
      COSFI=DCOS(FI)
      SINFI=DSIN(FI)
      COS2FI=COSFI**2-SINFI**2
      SIN2FI=2.D0*SINFI*COSFI
      GA=EC+EB*COS2FI
      GB=EC-EPS2*COS2TH+EB*SIN2TH*COS2FI
      T=ED/GB*COS2TH
      TA=1.D0-T
      GC=T*(7.D0*T-8.D0)+1.D0
      V=(EC*COS2FI+EB)/GA
      VA=2.D0*V**2-1.D0
      BA=1.D0-EA-2.D0*EA**2+(EPS2*EC-EB**2)*T-EB*ED*TA*V
      BB=EC-EPS2*T+EB*TA*V+EPS4*(DA*(.75D0-T*(7.5D0-T*8.75D0))
     $ -DB*GC*V*AAAA+DC*TA**2*VA*BBBB)
      DBADT=EPS2*EC-EB*(EB-ED*V)
      DBBDT=-EPS2-EB*V+EPS4*(DA*35.D0*(T-3.D0/7.D0)/2.D0
     $  -DB*AAAA*(14.D0*T-8.D0)*V-2.D0*TA*VA*DC*BBBB)
      DBADV=-EB*ED*TA
      DBBDV=EB*TA+EPS4*(4.D0*TA**2*V*DC*BBBB-DB*AAAA*GC)
      R=DSQRT(BA/(BB*(EC+EB)*(EC-EB)*ED))
      DRDTP=R*(DBADT/BA-DBBDT/BB)
      DRDTH=-DRDTP*ED*COSTH*SINTH*GA/GB**2
      DRDFI=DRDTP*ED*COS2TH*EB*SIN2TH*SIN2FI/GB**2-R*(DBADV/BA-DBBDV/BB)
     $   *SIN2FI*(EC**2-EB**2)/GA**2
      AR(1)=R*SINTH*COSFI
      AR(2)=R*SINTH*SINFI
      AR(3)=R*COSTH
      AN(1)=R*DRDFI*SINFI+AR(1)*(R*SINTH-DRDTH*COSTH)
      AN(2)=AR(2)*(R*SINTH-DRDTH*COSTH)-R*DRDFI*COSFI
      AN(3)=R*SINTH*(DRDTH*SINTH+R*COSTH)
      DS=DSQRT(AN(1)**2+AN(2)**2+AN(3)**2)
      RETURN
      END
      REAL FUNCTION LDROP(AA,ZZ,BS,BC)
* deformed liquid drop energy
      INTEGER AA,ZZ
      A=AA
      Z=ZZ
      AS=BS-1.0
      AC=BC-1.0
      LDROP=AS*(17.9439*A**(2.0/3.0)-31.9868*(A-2.0*Z)**2.0/
     F      A**(4.0/3.0))+
     F      AC*(0.70531*Z*Z/A**(1.0/3.0))
      RETURN
      END
      SUBROUTINE ABSCWT(ABSC,WT,NG)
C     sets up the abscissas and weights for gaussian quadrature of a
C     given order ng.
C
      DIMENSION ABSC(NG),WT(NG)
      DOUBLE PRECISION ABSC,WT
      IF (NG.EQ.4) GO TO 4
      IF (NG.EQ.8) GO TO 8
      IF (NG.EQ.12) GO TO 12
      IF (NG.EQ.16) GO TO 16
      IF (NG.EQ.32) GO TO 32
      IF (NG.EQ.64) GO TO 64
        WRITE(6,999)NG
  999 FORMAT(3HNG=,I3,' IS NOT AVAILABLE AT PRESENT. NG IS SET EQUAL',
     $'TO 16.'/)
      NG=16
      GO TO 16
    4 CONTINUE
      ABSC( 1)=  .339981043584856D0
      WT( 1)=  .652145154862546D0
      ABSC( 2)=  .861136311594053D0
      WT( 2)=  .347854845137454D0
      GO TO 1
    8 CONTINUE
      ABSC( 1)=  .18343464249564981D0
      WT( 1)=  .36268378337836199D0
      ABSC( 2)=  .52553240991632899D0
      WT( 2)=  .31370664587788729D0
      ABSC( 3)=  .79666647741362674D0
      WT( 3)=  .22238103445337447D0
      ABSC( 4)=  .96028985649753623D0
      WT( 4)=  .10122853629037627D0
      GO TO 1
   12 CONTINUE
      ABSC( 1)=  .12523340851146892D0
      WT( 1)=  .24914704581340279D0
      ABSC( 2)=  .36783149899818019D0
      WT( 2)=  .23349253653835481D0
      ABSC( 3)=  .58731795428661745D0
      WT( 3)=  .20316742672306592D0
      ABSC( 4)=  .76990267419430469D0
      WT( 4)=  .16007832854334623D0
      ABSC( 5)=  .90411725637047486D0
      WT( 5)=  .10693932599531843D0
      ABSC( 6)=  .98156063424671925D0
      WT( 6)=  .04717533638651182D0
      GO TO 1
   16 CONTINUE
      ABSC( 1)=  .09501250983763744D0
      WT( 1)=  .18945061045506850D0
      ABSC( 2)=  .28160355077925891D0
      WT( 2)=  .18260341504492359D0
      ABSC( 3)=  .45801677765722739D0
      WT( 3)=  .16915651939500254D0
      ABSC( 4)=  .61787624440264375D0
      WT( 4)=  .14959598881657673D0
      ABSC( 5)=  .75540440835500304D0
      WT( 5)=  .12462897125553387D0
      ABSC( 6)=  .86563120238783174D0
      WT( 6)=  .09515851168249278D0
      ABSC( 7)=  .94457502307323258D0
      WT( 7)=  .06225352393864789D0
      ABSC( 8)=  .98940093499164993D0
      WT( 8)=  .02715245941175408D0
      GO TO 1
   32 CONTINUE
      ABSC( 1)=  .04830766568773832D0
      WT( 1)=  .09654008851472780D0
      ABSC( 2)=  .14447196158279649D0
      WT( 2)=  .09563872007927486D0
      ABSC( 3)=  .23928736225213708D0
      WT( 3)=  .09384439908080457D0
      ABSC( 4)=  .33186860228212765D0
      WT( 4)=  .09117387869576389D0
      ABSC( 5)=  .42135127613063535D0
      WT( 5)=  .08765209300440381D0
      ABSC( 6)=  .50689990893222939D0
      WT( 6)=  .08331192422694676D0
      ABSC( 7)=  .58771575724076233D0
      WT( 7)=  .07819389578707031D0
      ABSC( 8)=  .66304426693021520D0
      WT( 8)=  .07234579410884851D0
      ABSC( 9)=  .73218211874028968D0
      WT( 9)=  .06582222277636185D0
      ABSC(10)=  .79448379596794241D0
      WT(10)=  .05868409347853555D0
      ABSC(11)=  .84936761373256997D0
      WT(11)=  .05099805926237618D0
      ABSC(12)=  .89632115576605213D0
      WT(12)=  .04283589802222668D0
      ABSC(13)=  .93490607593773969D0
      WT(13)=  .03427386291302143D0
      ABSC(14)=  .96476225558750643D0
      WT(14)=  .02539206530926206D0
      ABSC(15)=  .98561151154526834D0
      WT(15)=  .01627439473090567D0
      ABSC(16)=  .99726386184948156D0
      WT(16)=  .00701861000947009D0
      GO TO 1
   64 CONTINUE
      ABSC( 1)=  .02435029266342443D0
      WT( 1)=  .04869095700913972D0
      ABSC( 2)=  .07299312178779904D0
      WT( 2)=  .04857546744150343D0
      ABSC( 3)=  .12146281929612055D0
      WT( 3)=  .04834476223480296D0
      ABSC( 4)=  .16964442042399282D0
      WT( 4)=  .04799938859645831D0
      ABSC( 5)=  .21742364374000708D0
      WT( 5)=  .04754016571483031D0
      ABSC( 6)=  .26468716220876742D0
      WT( 6)=  .04696818281621002D0
      ABSC( 7)=  .31132287199021096D0
      WT( 7)=  .04628479658131442D0
      ABSC( 8)=  .35722015833766812D0
      WT( 8)=  .04549162792741815D0
      ABSC( 9)=  .40227015796399160D0
      WT( 9)=  .04459055816375657D0
      ABSC(10)=  .44636601725346409D0
      WT(10)=  .04358372452932346D0
      ABSC(11)=  .48940314570705296D0
      WT(11)=  .04247351512365359D0
      ABSC(12)=  .53127946401989455D0
      WT(12)=  .04126256324262353D0
      ABSC(13)=  .57189564620263404D0
      WT(13)=  .03995374113272034D0
      ABSC(14)=  .61115535517239325D0
      WT(14)=  .03855015317861563D0
      ABSC(15)=  .64896547125465734D0
      WT(15)=  .03705512854024005D0
      ABSC(16)=  .68523631305423324D0
      WT(16)=  .03547221325688239D0
      ABSC(17)=  .71988185017161083D0
      WT(17)=  .03380516183714161D0
      ABSC(18)=  .75281990726053190D0
      WT(18)=  .03205792835485155D0
      ABSC(19)=  .78397235894334141D0
      WT(19)=  .03023465707240248D0
      ABSC(20)=  .81326531512279756D0
      WT(20)=  .02833967261425948D0
      ABSC(21)=  .84062929625258036D0
      WT(21)=  .02637746971505466D0
      ABSC(22)=  .86599939815409282D0
      WT(22)=  .02435270256871087D0
      ABSC(23)=  .88931544599511411D0
      WT(23)=  .02227017380838325D0
      ABSC(24)=  .91052213707850281D0
      WT(24)=  .02013482315353021D0
      ABSC(25)=  .92956917213193958D0
      WT(25)=  .01795171577569734D0
      ABSC(26)=  .94641137485840282D0
      WT(26)=  .01572603047602472D0
      ABSC(27)=  .96100879965205372D0
      WT(27)=  .01346304789671864D0
      ABSC(28)=  .97332682778991097D0
      WT(28)=  .01116813946013112D0
      ABSC(29)=  .98333625388462596D0
      WT(29)=  .00884675982636395D0
      ABSC(30)=  .99101337147674432D0
      WT(30)=  .00650445796897837D0
      ABSC(31)=  .99634011677195528D0
      WT(31)=  .00414703326056246D0
      ABSC(32)=  .99930504173577214D0
      WT(32)=  .00178328072169639D0
      GO TO 1
    1 NGD2=NG/2
      DO 2 I=1,NGD2
      IP=NG+1-I
      ABSC(IP)=-ABSC(I)
    2 WT(IP)=WT(I)
      RETURN
      END
*ELEMENT ROT2G
       SUBROUTINE ROT2G(N,XJ,YJ,ZJ)
C     the three moments of inertia, relative to that of a sphere, of
C     a triaxial homogeneous mass distribution. note, however, that when
C     used. the quantity (hbar**2/2J) is easily obtained as
C     (36.*IA**(-5./3.)/XJ). recommended value: N=16.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL XJ,YJ,ZJ
      REAL OMOM
      DIMENSION AR(3),AN(3),X(64),W(64)
      COMMON /ECOM/EPS2,EPS22,EPS4D,OMOM,AR,AN,DS,R,SINTHA
      DATA PI /3.14159265358979324D0/
      CALL ABSCWT(X,W,N)
      H=PI/4.D0
      XJD=0.D0
      YJD=0.D0
      ZJD=0.D0
      DO 4 I=1,N
      TH=H*(X(I)+1.D0)
      TH2=TH+2*H
      SINTH=DSIN(TH)
      SINTH2=DSIN(TH2)
      DO 4 J=1,N
      FI =H*(X(J)+1.D0)
      COSFI=DCOS(FI)
      CALL NRDS(TH,FI)
      R5=R**5*W(I)*W(J)
      CALL NRDS(TH2,FI)
      R52=R**5*W(I)*W(J)
      XJD=XJD+R5*(1.D0-(SINTH*COSFI)**2)*SINTH
     ++R52*(1.D0-(SINTH2*COSFI)**2)*SINTH2
      YJD=YJD+R5*(1.D0-SINTH**2*(1.D0-COSFI**2))*SINTH
     ++R52*(1.D0-SINTH2**2*(1.D0-COSFI**2))*SINTH2
      ZJD=ZJD+R5*SINTH**3+R52*SINTH2**3
    4 CONTINUE
      XJD=8.D0*XJD/10.D0*(PI/4.D0)**2/(8.D0*PI/15.D0)
      YJD=8.D0*YJD/10.D0*(PI/4.D0)**2/(8.D0*PI/15.D0)
      ZJD=8.D0*ZJD/10.D0*(PI/4.D0)**2/(8.D0*PI/15.D0)
      XJ=XJD/OMOM**5
      YJ=YJD/OMOM**5
      ZJ=ZJD/OMOM**5
      RETURN
      END
      SUBROUTINE DNICES(NN,NNE,NNV,A,NMA,EPS,IORD,E,V,NMV,ILL,
     &                  W1,W2,W3,W4,W5,W6)
C **** DNICES(NN,NNE,NNV,A,NMA,EPS,IORD,E,V,NMV,ILL,W1,W2,W3,W4,W5,W6) *
C * DIAGNALIZATION OF SYMMETRIC MATRIX                                 *
C *                 ORIGINAL BY Y. BEPPU AND I. NINOMIYA, 1980         *
C *I------------------------------------------------------------------I*
C *I        ::: "NSHOUD" IN LIBRARY 'NICER' :::                       I*
C *I SUBPROGRAM FOR STANDARD EIGEN-PROBLEM , A*V=V*E ,   BY           I*
C *I(HOUSEHOLDER)-(BISECTION & NO-ROOT-QR)-(INVERSE-ITERATION) METHOD I*
C *I---------------------------------------------(VERSION-2,LEVEL-1)--I*
C * ::: INPUT :::                                                      *
C * NN       --- DIM. OF DIAGONALIZATION PROBLEM                       *
C * NNE,NNV  --- NUMBER OF EIGEN-VALUES , VECTORS TO BE CALCULATED     *
C * A,NMA    --- INPUT SYMMETRIC MATRIX AND ITS ADJUSTABLE ARREY DIM.  *
C *                  NOTE : ONLY THE UPPER-HALF PART IS USED           *
C * EPS      --- ACCURACY PARAMETER FOR CONVERGENCE   ((STD.:1.D-10))  *
C * IORD     --- < 0 , > 0 FOR ASCENDING , DESCENDING ORDER            *
C *                          OF EIGEN-VALUES                           *
C * ::: OUTPUT :::                                                     *
C * E        --- EIGEN-VALUES                                          *
C * V,NMV    --- EIGEN-VETORS MATRIX AND ITS ADJUSTABLE ARREY DIM.     *
C * ILL      --- CONDITION CODE = 0 FOR NORMAL END                     *
C * ::: DUMMY :::                                                      *
C * W1-W6    --- WORKING SPACE-ARREYS                                  *
C **********************************************************************
C     IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NMA,NN),E(NN),V(NMV,NN)        ,LS(20),KS(20)
      DIMENSION W1(NN),W2(NN),W3(NN),W4(NN),W5(NN),W6(NN)
C     DATA ZERO,HALF,FIBONA,ONE /0.D 0,0.5D 0,0.6180339D 0,1.D 0/
C     DATA EXPM30,EXPM20,EXPM6  /1.D-30,1.D-20,1.D-6/
C     DATA EXPP6,EXPP12,EXPP18  /1.D 6,1.D12,1.D18/
      DATA ZERO,HALF,FIBONA,ONE /0.,0.5,0.6180339,1./
      DATA EXPM30,EXPM20,EXPM6  /1.E-30,1.E-20,1.E-6/
      DATA EXPP6,EXPP12,EXPP18  /1.E+6,1.E12,1.E18/
      DATA LOOPMX/100/
C-----------------------------------------------------------------------
C        INITIALIZATION
C-----------------------------------------------------------------------
      NMAX=NMA
      N=NN
      NE=NNE
      NV=NNV
      ILL=300
      IF(NMAX.LT.2.OR.N.GT.NMAX.OR.N.LT.2.OR.NE.LT.1.OR.NV.LT.0) RETURN
      IF(NE.GT.N.OR.NV.GT.NE.OR.EPS.LE.ZERO.OR.EPS.GT.ONE)       RETURN
      NM1=N-1
      NM2=N-2
      NE8=NE*8
      IF(N.EQ.2)                  GO TO 130
C-----------------------------------------------------------------------
C        TRI-DIAGONALIZATION                      ( BY A.HOUSEHOLDER )
C-----------------------------------------------------------------------
      DO 120 K=1,NM2
      KP1=K+1
      E(K)=A(K,K)
      SUM=ZERO
      DO 10 J=KP1,N
      E(J)=A(K,J)
   10 SUM=E(J)*E(J)+SUM
C     S=DSIGN(DSQRT(SUM),E(KP1))
      S=SIGN(SQRT(SUM),E(KP1))
      W1(K)=-S
      E(KP1)=E(KP1)+S
      A(K,KP1)=E(KP1)
      H=E(KP1)*S
      IF(H.EQ.ZERO)     GO TO 110
      SUMM=ZERO
      DO 40 I=KP1,NM1
      SUM=ZERO
      DO 20 J=KP1,I
   20 SUM=A(J,I)*E(J)+SUM
      IP1=I+1
      DO 30 J=IP1,N
   30 SUM=A(I,J)*E(J)+SUM
      W1(I)=SUM/H
      SUMM=W1(I)*E(I)+SUMM
   40 CONTINUE
      SUM=ZERO
      DO 50 J=KP1,N
   50 SUM=A(J,N)*E(J)+SUM
      W1(N)=SUM/H
      SUMM=W1(N)*E(N)+SUMM
      U=SUMM*HALF/H
C------- FOR A SCALAR MACHINE ------------------------------------------
61    DO 60 J=KP1,N
      W1(J)=E(J)*U-W1(J)
      DO 60 I=KP1,J
   60 A(I,J)=W1(I)*E(J)+W1(J)*E(I)+A(I,J)
C------- FOR A VECTOR MACHINE ------------------------------------------
C65    CONTINUE
C      DO 100 J=KP1,N
C      W1(J)=E(J)*U-W1(J)
C      EJ=E(J)
C      DO 80 I=KP1,J
C   80 A(I,J)=W1(I)*EJ+A(I,J)
C      W1J=W1(J)
C      DO 90 I=KP1,J
C   90 A(I,J)=E(I)*W1J+A(I,J)
C  100 CONTINUE
C-----------------------------------------------------------------------
  110 A(K,K)=H
  120 CONTINUE
C
  130 CONTINUE
      E(NM1)=A(NM1,NM1)
      E(N)=A(N,N)
      W1(NM1)=A(NM1,N)
      W1(N)=ZERO
C     GERSCH=DABS(E(1))+DABS(W1(1))
      GERSCH=ABS(E(1))+ABS(W1(1))
      DO 140 I=1,NM1
C     SUM=DABS(E(I+1))+DABS(W1(I))+DABS(W1(I+1))
      SUM=ABS(E(I+1))+ABS(W1(I))+ABS(W1(I+1))
      IF(SUM.GT.GERSCH) GERSCH=SUM
      IF(NV.NE.0)  V(I,NV)=E(I)
  140 CONTINUE
      IF(NV.NE.0)  V(N,NV)=E(N)
      DEL=EPS*GERSCH
      IF(DEL.EQ.ZERO)             RETURN
      ILL=0
      IF(NE8.LT.N)                GO TO 360
C-----------------------------------------------------------------------
C        NO-ROOT-QR METHOD FOR EIGENVALUES
C                    (BY PAL-WALKER-KAHAN & M.SHIMASAKI)
C-----------------------------------------------------------------------
      DD=DEL*DEL
      DO 150 I=1,NM1
  150 W3(I+1)=W1(I)*W1(I)
      K=N
      LOO=0
  160 KM1=K-1
      LOO=LOO+1
      IF(W3(K).LT.DD) GO TO 190
      IF(LOO.GT.LOOPMX)THEN
          ILL=ILL+1
          GO TO 190
      ENDIF
      EPE=(E(KM1)+E(K))*HALF
      EME=E(K)-EPE
C      QRSHIF=EPE-DSIGN(DSQRT(W3(K)+EME*EME),EPE)
      QRSHIF=EPE-SIGN(SQRT(W3(K)+EME*EME),EPE)
      CC=ONE
      SS=ZERO
      G=E(1)-QRSHIF
      PP=G*G
      DO 180 I=1,KM1
      BB=W3(I+1)
      TT=PP+BB
      W3(I)=SS*TT
      OLDCC=CC
      SS=BB/TT
      CC=PP/TT
      OLDG=G
      IF(CC.EQ.ZERO) GO TO 170
      G=(E(I+1)-QRSHIF)*CC-OLDG*SS
      E(I)=E(I+1)+OLDG-G
      PP=(G*G)/CC
      GO TO 180
  170 G=-OLDG
      E(I)=E(I+1)+OLDG+OLDG
      PP=BB*OLDCC
  180 CONTINUE
      E(K)=G+QRSHIF
      W3(K)=SS*PP
      GO TO 160
  190 K=K-1
      LOO=0
      IF(K.GT.1) GO TO 160
C-----------------------------------------------------------------------
C        QUICK SORT OF EIGENVALUES                ( BY C.HOARE )
C-----------------------------------------------------------------------
      ISP=0
      L=1
      K=N
  200 IF(K-L.LT.16)          GO TO 290
      M=(K+L)/2
      MAX=K
      IF(E(M).GT.E(K))   MAX=M
      IF(E(L).GT.E(MAX)) MAX=L
      IF(MAX.EQ.K)           GO TO 210
      BK=E(MAX)
      E(MAX)=E(K)
      E(K)=BK
  210 IF(E(L).GE.E(M))       GO TO 220
      BK=E(L)
      E(L)=E(M)
      E(M)=BK
  220 BK=E(L)
      I=L
      J=K
      GO TO 250
  230 E(J)=E(I)
      E(I)=BK
  240 J=J-1
  250 IF(BK.LT.E(J))         GO TO 240
      IF(J.LE.I)             GO TO 270
      E(I)=E(J)
      E(J)=BK
  260 I=I+1
      IF(E(I).LT.BK)         GO TO 260
      IF(J.GT.I)             GO TO 230
  270 ISP=ISP+1
      IF(K-I.GE.I-L)         GO TO 280
      LS(ISP)=L
      KS(ISP)=I-1
      L=I+1
      GO TO 200
  280 LS(ISP)=I+1
      KS(ISP)=K
      K=I-1
      GO TO 200
  290 IF(K-L.LT.1)           GO TO 330
      J=K
  300 BK=E(J-1)
      I=J
  310 IF(E(I).GE.BK)         GO TO 320
      E(I-1)=E(I)
      I=I+1
      IF(I.LE.K)             GO TO 310
  320 E(I-1)=BK
      J=J-1
      IF(J.GT.L)             GO TO 300
  330 IF(ISP.EQ.0)           GO TO 340
      L=LS(ISP)
      K=KS(ISP)
      ISP=ISP-1
      GO TO 200
  340 IF(IORD.LT.0)               GO TO 500
      NP1=N+1
      NE2=N/2
      DO 350 K=1,NE2
      TEMP=E(K)
      E(K)=E(NP1-K)
  350 E(NP1-K)=TEMP
      GO TO 500
C-----------------------------------------------------------------------
C        BISECTION METHOD FOR EIGENVALUES         ( BY J.GIVENS)
C-----------------------------------------------------------------------
  360 CONTINUE
      DO 370 I=1,NE
      W2(I)=E(I)
      E(I)=-GERSCH
      W4(I)=GERSCH
  370 W3(I)=-W1(I)*W1(I)
      NEP1=NE+1
      DO 380 I=NEP1,N
      W2(I)=E(I)
  380 W3(I)=-W1(I)*W1(I)
      IF(IORD.GT.0)               GO TO 400
      DO 390 I=1,N
  390 W2(I)=-W2(I)
  400 CONTINUE
      DO 490 K=1,NE
      LOO=0
  410 X=(E(K)+W4(K))*HALF
      LOO=LOO+1
C     IF(DABS(W4(K)-X).LE.DEL)    GO TO 490
      IF(ABS(W4(K)-X).LE.DEL)    GO TO 490
      IF(LOO.GT.LOOPMX)THEN
          ILL=ILL+1
          GO TO 490
      ENDIF
      NAG=0
      I=1
  420 S=W2(I)-X
  430 IF(S.GE.ZERO) NAG=NAG+1
C     IF(DABS(S).LT.EXPM30)       GO TO 440
      IF(ABS(S).LT.EXPM30)       GO TO 440
      I=I+1
      IF(I.GT.N)                  GO TO 450
      S=W3(I-1)/S+W2(I)-X
      GO TO 430
  440 I=I+2
      IF(I.LE.N)                  GO TO 420
  450 IF(NAG.GE.K)                GO TO 470
      DO 460 J=K,NE
      IF(X.LT.W4(J)) W4(J)=X
  460 CONTINUE
      GO TO 410
  470 MG=NAG
      IF(NE.LT.MG) MG=NE
      DO 480 J=K,MG
  480 E(J)=X
      GO TO 410
  490 E(K)=X
C-----------------------------------------------------------------------
C        INVERSE ITERATION FOR EIGENVECTORS      ( BY H.WIELANDT )
C-----------------------------------------------------------------------
  500 CONTINUE
      IF(NV.EQ.0)                 GO TO 810
      IF(NE8.GE.N .OR. IORD.GT.0) GO TO 520
      DO 510 I=1,N
      W1(I)=-W1(I)
  510 V(I,NV)=-V(I,NV)
C 520 FN=DFLOAT(N)
  520 FN=FLOAT(N)
C     SN=DSQRT(FN)
      SN=SQRT(FN)
C     SEPS=DSQRT(EPS)
      SEPS=SQRT(EPS)
      EPS1=(GERSCH*EXPM6)/(FN*SEPS)
      EPS2=EXPP12*FN/SEPS
      TM6N=EXPM6*SN
      TP6N=EXPP6*SN
      TP18N=EXPP18*SN
C
      RN=ZERO
      RA=EPS*FIBONA
      IG=1
      DO 750 I=1,NV
      IM1=I-1
      DO 530 J=1,N
      W3(J)=ZERO
      W4(J)=W1(J)
      W5(J)=V(J,NV)-E(I)
      RN=RN+RA
      IF(RN.GE.EPS) RN=RN-EPS
  530 W6(J)=RN
      DO 580 J=1,NM1
C     IF(DABS(W5(J)).GE.DABS(W1(J)))   GO TO 560
      IF(ABS(W5(J)).GE.ABS(W1(J)))   GO TO 560
      W2(J)=-W5(J)/W1(J)
      W5(J)=W1(J)
      T=W5(J+1)
      W5(J+1)=W4(J)
      W4(J)=T
      W3(J)=W4(J+1)
      IF(W3(J).EQ.ZERO) W3(J)=DEL
      W4(J+1)=ZERO
      GO TO 570
  560 IF(W5(J).EQ.ZERO) W5(J)=DEL
      W2(J)=-W1(J)/W5(J)
  570 W4(J+1)=W3(J)*W2(J)+W4(J+1)
  580 W5(J+1)=W4(J)*W2(J)+W5(J+1)
      IF(W5(N).EQ.ZERO) W5(N)=DEL
      ITELIM=2
      IF(EPS.LT.EXPM20)                          ITELIM=3
C     IF(I.GT.1.AND.DABS(E(I)-E(IM1)).LT.EPS1)   ITELIM=3
      IF(I.GT.1.AND.ABS(E(I)-E(IM1)).LT.EPS1)   ITELIM=3
      DO 640 IT=1,ITELIM
      IF(IT.EQ.1)       GO TO 600
      DO 590 J=1,NM1
      IF(W3(J).EQ.ZERO) GO TO 590
      T=W6(J)
      W6(J)=W6(J+1)
      W6(J+1)=T
  590 W6(J+1)=W6(J)*W2(J)+W6(J+1)
  600 W6(N)=W6(N)/W5(N)
      W6(NM1)=(W6(NM1)-W6(N)*W4(NM1))/W5(NM1)
C     VN=DABS(W6(N))+DABS(W6(NM1))
      VN=ABS(W6(N))+ABS(W6(NM1))
      IF(N.EQ.2)                       GO TO 620
      DO 610 KK=2,NM1
      K=N-KK
      W6(K)=(W6(K)-W6(K+1)*W4(K)-W6(K+2)*W3(K))/W5(K)
C 610 VN=DABS(W6(K))+VN
  610 VN=ABS(W6(K))+VN
  620 IF(VN.GT.TM6N.AND.VN.LT.TP6N)    GO TO 640
      IF(IT.EQ.ITELIM)                 GO TO 650
      IF(VN.GT.EPS2)                   GO TO 650
      DUMP=FN/VN
      DO 630 J=1,N
  630 W6(J)=W6(J)*DUMP
  640 CONTINUE
C
C------- RE-ORTHOGONALIZATION ------------------ ( BY GRAM & SCHMIDT )
C
  650 CONTINUE
      DO 660 J=IG,I
C     IF(DABS(E(J)-E(I)).LT.EPS1)      GO TO 670
      IF(ABS(E(J)-E(I)).LT.EPS1)      GO TO 670
  660 CONTINUE
      J=I
  670 IG=J
      IF(IG.EQ.I.AND.VN.LT.TP18N)      GO TO 720
      DUMP=FN/VN
      DO 680 J=1,N
  680 W6(J)=W6(J)*DUMP
      IF(IG.EQ.I)                      GO TO 720
      DO 710 K=IG,IM1
      SUM=ZERO
      DO 690 J=1,N
  690 SUM=V(J,K)*W6(J)+SUM
      S=-SUM
      DO 700 J=1,N
  700 W6(J)=V(J,K)*S+W6(J)
  710 CONTINUE
C
C------- NORMALIZATION ------------------------------------------------
C
  720 SUM=ZERO
      DO 730 J=1,N
  730 SUM=W6(J)*W6(J)+SUM
C     SINV=ONE/DSQRT(SUM)
      SINV=ONE/SQRT(SUM)
      DO 740 J=1,N
  740 V(J,I)=W6(J)*SINV
  750 CONTINUE
C-----------------------------------------------------------------------
C        BACK-TRANSFORMATION OF EIGEN-VECTORS
C-----------------------------------------------------------------------
      IF(N.EQ.2)                       GO TO 810
      DO 800 J=1,NM2
      K=N-J-1
      IF(A(K,K).EQ.ZERO)               GO TO 800
      DO 760 KK=K,N
  760 W1(KK)=A(K,KK)
      KP1=K+1
      DO 790 I=1,NV
      SUM=ZERO
      DO 770 KK=KP1,N
  770 SUM=V(KK,I)*W1(KK)+SUM
      S=-SUM/W1(K)
      DO 780 KK=KP1,N
  780 V(KK,I)=W1(KK)*S+V(KK,I)
  790 CONTINUE
  800 CONTINUE
C-----------------------------------------------------------------------
C        ENDING
C-----------------------------------------------------------------------
  810 CONTINUE
C     ILL=0
      IF(NE8.GE.N .OR. IORD.GT.0) RETURN
      DO 820 I=1,NE
  820 E(I)=-E(I)
      RETURN
      END
      SUBROUTINE XYCOORD(EPSAI,GAMMA,AX,AY,NX,NY)
* polar to cartesian transformation
      INCLUDE 'NICRAINC.FOR'
      DATA PIFAC/57.29577/
      AY=EPSAI*SIN(GAMMA-(120+FI)/PIFAC)
      AX=EPSAI*COS(GAMMA-(120+FI)/PIFAC)
      NY=AY+0.1
      NX=AX+0.1
      RETURN
      END
      SUBROUTINE FUN2(X,Y,R,THETA)
* cartesian to polar transformation
      INCLUDE 'NICRAINC.FOR'
      DATA PIFAC/57.29577/
      AX=X
      AY=Y
      R=SQRT(AX*AX+AY*AY)
      IF(ABS(AX).GE.0.001)ALFA=ATAN(AY/AX)*PIFAC
      IF(ABS(AX).LT.0.001.AND.AY.GT.0)ALFA=90.
      IF(ABS(AX).LT.0.001.AND.AY.LT.0)ALFA=-90.
      IF(AX.LT.-.001)ALFA=ALFA-180.
      THETA=ALFA+FI
      RETURN
      END
      DOUBLE PRECISION FUNCTION DCBRT(ARG)
      DOUBLE PRECISION ARG
      DCBRT=ARG**(1./3.)
      RETURN
      END
      SUBROUTINE BARDEE (NOI,E,GOI,TAI,D,L,V2,IFLAG,KOEFF)
* performs iteration in BCS-calc
      REAL L
      DIMENSION V2(1),E(1)
      REAL EGI,EI,A,B,G,NR,S,DG,DN,Z,EIL,DL,DD
      INTEGER I
      COMMON /AREA/K1,K2,BDNDL
      NO=NOI
      GO=GOI
      TA=TAI
      IF (NO .LT.1141) GO TO 350
  100 FORMAT(1X,43HFIRST PARAMETER IN BARDEE GREATER THAN 1141)
      WRITE(6,100)
      STOP
  350 CONTINUE
      IFLAG=0
  351 CONTINUE
      NR=0.0
      S=NR
      B=S
      A=B
      K1=NO/2-IFIX(SQRT(FLOAT(KOEFF*NO)))
      IF(K1.GE.1) GO TO 249
  101 FORMAT(1X,29HK1 IN BARDEE IS LESS THAN ONE)
C     WRITE(6,101)
C     STOP
       K1=1
  249 CONTINUE
      K2=NO/2+IFIX(SQRT(FLOAT(KOEFF*NO)))
       IF(K1.EQ.1)K2=NO
CB    WRITE(6,*)' K1,K2=',K1,K2
      DO 352 I=K1,K2
      EI=E(I)
      EIL=EI-L
      EGI=SQRT(EIL*EIL+ABS(D))
      EGI3=1.0/EGI**3
      A=A+EIL*EGI3
      B=B+EGI3
      S=S+1.0/EGI
      NR=NR+1.0-EIL/EGI
  352 CONTINUE
      G=2.0/S
      DG=G-GO
      DN=NR-NO+2.0*(K1-1.0)
CB    WRITE(6,10)G,DG,DN
CB10  FORMAT(' G,DG,DN=',3F8.4)
      IF (ABS(DG)+ABS(DN) .LT. TA) GO TO 353
      Z=A*A+D*B*B
      DL=(-2.0*A/G/G*DG+B*DN)/Z
      DD=2.0*(A*DN+2.0*D*B/G/G*DG)/Z
      IF(DN/DL.GT.0)BDNDL=DN/DL
      L=L-DL
      D=D-DD
CB    WRITE(6,11)Z,L,DL,D,DD
CB11  FORMAT(' Z=',F11.3,' L,DL=',2F7.3,' D,DD=',2F7.3)
      IF ((D .GE. -3.0) .AND. (D .LE. 100.0)) GO TO 351
      D=0.04
      IFLAG=1
      RETURN
  353 IF (D.GE.0.04) GO TO 355
      D=0.04
      IFLAG=1
      RETURN
  355 DO 354 I=K1,K2
      EIL=E(I)-L
      EGI=SQRT(EIL*EIL+D)
      V2(I)=0.5*(1.0-(E(I)-L)/EGI)
  354 CONTINUE
      RETURN
      END
