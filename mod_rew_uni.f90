!---------------------------------------------------------------------------
module reweight
    use expect 
    use auxiliary 
    use variables

    implicit none
    !-----waves and determinants 
    integer*4,parameter :: nf_rew = n_f-1 ,np_rew = n_p 
    !---Slater determinants on the left/right hand side
    complex*16 zwaveL_rew(L3, 0:1, 0:1, nf_rew)
    complex*16 zwaveR_rew(L3, 0:1, 0:1, nf_rew)

    !---inverse Z^-1 and determinant det(Z) of the correlation matrix
    !---probabilty weight is absorbed into zDet 
    complex*16 zcorrinv_rew(nf_rew, nf_rew) ! For S_n  
    complex*16 zcorrmatrix_rew(nf_rew, nf_rew)
    complex*16 zDet_rew

    complex*16 zcorrinv2_rew(n_f-2, n_f-2)  ! For S_2n 
    complex*16 zcorrmatrix2_rew(n_f-2, n_f-2)
    complex*16 zDet2_rew

    complex*16 zphase_rew(3)
    
    !---Bins for collecting data
    complex*16 zIabsBin_rew(6)         !<det(Z)> and <|det(Z)|> for phase
    complex*16 zIyeBin_rew(dimB)       !<det(Z)>, common denominator
    complex*16 zIyeplusBin_rew(dimB)       !<det(Z)>, common denominator
    complex*16 zIyemnusBin_rew(dimB)       !<det(Z)>, common denominator
    
    complex*16 zIabsBin2_rew(6)         !<det(Z)> and <|det(Z)|> for phase
    complex*16 zIyeBin2_rew(dimB)       !<det(Z)>, common denominator
    !---leading order energy
    complex*16 zKinBin_rew(dimB)       !kinetic energy <K>
    complex*16 zQ0SBin_rew(dimB, 3)    !Q0 SU(4) conserving term
    complex*16 zQSIBin_rew(dimB, 3)    !Q0 SU(4) conserving term
    complex*16 zVLOLBin_rew(dimB, 100) !Q0 SU(4) conserving term
    complex*16 zVLONLBin_rew(dimB, 100)!Q0 SU(4) conserving term
    complex*16 zOPEBin_rew(dimB, 2)    !LO OPEP term <V_pion>
    complex*16 zOPWBin_rew(dimB)       !LO OPEP weakened <V_pion_weak>
    complex*16 zQ0VBin_rew(dimB)       !Q0 SU(4) breaking term
    complex*16 zVdiffN3LOBin_rew(dimB,nopr)
    complex*16 zVdiffN3LOGIRBin_rew(dimB,3*nopr)
    
    complex*16 zTMXBin_rew(dimB)       !transfer matrix
    !---higher order energies
    complex*16 zGIRBin_rew(dimB, 3)         !Galilean Invariance Restoration
    complex*16 zCouBin_rew(dimB, 2)         !COUlomb energy
    
    complex*16 zTBFBin_rew(dimB, 20)         !ThreeBody Forces
    complex*16 zVGIR0QNBin_rew(dimB, nopr)  !higher order contact terms
    complex*16 zVGIR1QNBin_rew(dimB, nopr)  !higher order contact terms
    complex*16 zVGIR2QNBin_rew(dimB, nopr)  !higher order contact terms
    complex*16 zA3NBin_rew(dimB, 20)        !additional 3N forces

    !---density related observables
    complex*16 zrhoBin_rew(dimB, 4)    !local density correlation \rho^i, i=1:4
    
    complex*16 zRSQBin_rew(dimB, 2)    !r. m. s. radii
    complex*16 zCNMBin_rew(dimB, 3)    !charge radii
    
    !---for thermodynamics
    
    !---for pinhole algorithm
    
    !---for quick test
    complex*16 zXXXBin_rew(dimB, 9)    !test measurement
    
    !------printing  
contains
    subroutine initialize_rew()
      ! At initialization stage,
      ! prepare initial wave functions 
      ! 
      ! initialization of waves 
      implicit none  
      if (n_c < 1) stop 'Error: Re-weight for multi-channel not available yet.'

      if ((nf_rew .ne. n_f).or.(np_rew .ne. n_p)) then 
        write(*,*) 'Initialize waves for ReWeighting calculation.'
        if ((nf_rew < n_f).and.(np_rew==n_p)) then 
          write(*,*) 'No need of new waves. Proceed.'
        else
          write(*,*) 'Need to setup new waves. Under Construction' 
          stop  
        end if 
      end if  
      !---initialize counters
      zIabsBin_rew = Z0
      zIyeBin_rew = Z0 
      zIyeplusBin_rew = Z0 
      zIyemnusBin_rew = Z0 
      zKinBin_rew = Z0 
      zQ0SBin_rew = Z0 
      zQSIBin_rew = Z0 
      zVLOLBin_rew = Z0 
      zVLONLBin_rew = Z0 
      zOPEBin_rew = Z0 
      zOPWBin_rew = Z0 
      zQ0VBin_rew = Z0 
      zVdiffN3LOBin_rew = Z0 
      zVdiffN3LOGIRBin_rew = Z0 
      zTMXBin_rew = Z0 
      zGIRBin_rew = Z0 
      zCouBin_rew = Z0 
      zTBFBin_rew = Z0 
      zVGIR0QNBin_rew = Z0 
      zVGIR1QNBin_rew = Z0 
      zVGIR2QNBin_rew = Z0 
      zA3NBin_rew = Z0 
      zrhoBin_rew = Z0 
      zRSQBin_rew = Z0 
      zCNMBin_rew = Z0 
      zXXXBin_rew = Z0 

      zIabsBin2_rew = Z0
      zIyeBin2_rew = Z0 

    end subroutine initialize_rew

    subroutine prepare_expect()
      ! 
      ! INPUT : existing zwaveL, zwaveR, zDet, zcorrmatrix, zcorrinv 
      ! OUTOUT: new      zwaveL_rew ... 
      !
      ! This should be called everytime before measurement
      ! when zves are updated.                
      implicit none   
      complex*16 cph(L3)
      integer*4 n, ni, ns, np
      !------from here---------------------------------
      ! define zDet, zcorrmatrix, zwaveL, zwaveR
      !       for the REW calculation 
      !       in each problems 
      zcorrmatrix_rew = zcorrmatrix_A(1:nf_rew,1:nf_rew) 
      zcorrinv_rew = inv(zcorrmatrix_rew)
      !---zDet = det(Z)/|det(W)|,  zDet is not det(Z)
      zDet_rew = exp(getLogDet(zcorrmatrix_rew)-dble(getLogDet(zcorrmatrix_A)))
      zDet_rew = zDet_rew * cdexp( -sum( cu_A * cu_i ) )       

      !----For S_2n 
      zcorrmatrix2_rew = zcorrmatrix_A(1:(n_f-2),1:(n_f-2)) 
      zcorrinv2_rew = inv(zcorrmatrix2_rew)
      !---zDet = det(Z)/|det(W)|,  zDet is not det(Z)
      zDet2_rew = exp(getLogDet(zcorrmatrix2_rew)-dble(getLogDet(zcorrmatrix_A)))
      zDet2_rew = zDet2_rew * cdexp( -sum( cu_A * cu_i ) )  

      !---<bra| and |ket> for calculating obeservables
      zwaveL_rew = zvecs(:, :, :, 1:nf_rew, ntsweep) 
      call MulbyTM(zvecs(:, :, :, 1:nf_rew, ntsweep - 1), zwaveR_rew, ntsweep, +1)

      !-------end here--------------------------------- 
      !---multiply CM phase factor 
      cph = cdexp( Zi * ( nPcm(1) * qx + nPcm(2) * qy + nPcm(3) * qz ) )
      do np = 1, nf_rew 
          do ni = 0, 1; do ns = 0, 1
              zwaveL_rew(:, ns, ni, np) = zwaveL_rew(:, ns, ni, np) * conjg(cph)
              zwaveR_rew(:, ns, ni, np) = zwaveR_rew(:, ns, ni, np) * cph
          enddo; enddo
      enddo      
      return 
    end subroutine prepare_expect 
    

    !--------------------subroutines for expect-------------------------------------
    function getTMX_test(zwaveL_,zwaveR_,zDet_,zcorrmatrix_,nf_) result(res)
        !---expectation value of transfer matrix <M> with importance sampling.
        !---
        !---Be careful that all inputs are zcorrmatrix_ is for nf_ 
        !---     zDet_ = zDet_(nf_)/|det Z(n_f)| 
        !---However, on the other hand, we now want to sample auxiliary fields
        !---according to < zwaveL_|M(s)| zwaveR_> for nf_ particles.
        !---  Note zDet_ is not used here actually. 
        implicit none 
        complex*16 res
        integer*4,intent(in) :: nf_ 
        complex*16,intent(in) :: zwaveL_(L3,0:1,0:1,nf_)
        complex*16,intent(in) :: zwaveR_(L3,0:1,0:1,nf_)
        complex*16,intent(in) :: zDet_ 
        complex*16,intent(in) :: zcorrmatrix_(nf_,nf_) !!! Be careful. This time input is zcorrmatrix not it's inverse. 

        real*8 s0(L3), sI(L3, 3), sv(L3, 0:3), pn(L3, 3), cu(L3), &
            Act, logProb, u(3), f(L3)
        real*8 dVds0(L3), dVdsI(L3, 3), dVdsv(L3, 0:3), dVdpn(L3, 3), dVdcu(L3)
        complex*16 cdVdcu(L3), cdVdsI(L3, 3)
        complex*16 zwave(L3, 0:1, 0:1, nf_), zcorrmatrix(nf_, nf_), &
            zcorrinv(nf_, nf_), zLogDet, cs
        integer*4 n, nRepeat
        integer*4,parameter  :: MeasureRepeat = 16

        res = Z0

        !---Act0 = -log(|det(W)|)
        Act = -dble(getLogDet(zcorrmatrix_A)) !!! not zcorrmatrix_ !!!  
        
        !---calculate derivative d/ds(logP(s)) |s = 0
        zcorrinv(:, :) = inv(zcorrmatrix_)   

        !call getdiff_test('MID', zwaveL_, zwaveR_, zcorrinv, &
        !    dVds0, cdVdsI, dVdsv, dVdpn, cdVdcu,nf_) 
        call getdiff_test('MID', zwaveL_, zwaveR_, zcorrinv, &
            dVds0, dVdsI, dVdsv, dVdpn, cdVdcu,nf_)     
        dVdcu = dble(cdVdcu)
        ! dVdsI = dble(cdVdsI)

        !---normalization factor
        do n = 1, L3
            u(:) = -dVds0(n) * svalue(:)
            f(n) = dlog( sum( sweigh(:) * &
                dexp(u - maxval(u)) ) ) + maxval(u)
        enddo

        cs = Z0
        do nRepeat = 1, MeasureRepeat 

            !---normalize probability to avoid overflow
            do n = 1, L3
                u(:) = -dVds0(n) * svalue(:)
                s0(n) = svalue( randint( 3, &
                    sweigh(:) * dexp(u - maxval(u)) ) )
            enddo

            call randfield(sI, 3 * L3, zero, one)
            sI = sI - dVdsI

            call randfield(sv, 4 * L3, zero, one)
            sv = sv - dVdsv

            call randfield(pn, 3 * L3, zero, one)
            pn = pn - dVdpn

            call randfield(cu, L3, zero, one)
            cu = cu - dVdcu

            !---normalized probability correction factor
            logProb = -sum(dVds0 * s0 + f) &
                -sum(dVdsI * sI + dVdsI**2 / 2) &
                -sum(dVdsv * sv + dVdsv**2 / 2) &
                -sum(dVdpn * pn + dVdpn**2 / 2) &
                -sum(dVdcu * cu + dVdcu**2 / 2)

            call trans('MID', s0, sI, sv, pn, cu + Z0, zwaveR_, zwave, nf_, +1)
            call getcorr(zwaveL_, zwave, zcorrmatrix) !new zcorrmatrix

            cs = cs + cdexp( getLogDet(zcorrmatrix) - logProb + Act ) &
                      * cdexp( -sum( cu_A * cu_i ) )

        enddo

        res = res + cs / MeasureRepeat

        return
    end function getTMX_test

    function getKin_test(zwaveL_,zwaveR_,zDet_,zcorrinv_,nf_) result(res)
        !---measure the kinetic energy <K>
        implicit none 
        complex*16 res
        integer*4,intent(in) :: nf_ 
        complex*16,intent(in) :: zwaveL_(L3,0:1,0:1,nf_)
        complex*16,intent(in) :: zwaveR_(L3,0:1,0:1,nf_)
        complex*16,intent(in) :: zDet_ 
        complex*16,intent(in) :: zcorrinv_(nf_,nf_) 
    
        complex*16 zKin(nf_, nf_), zwave(L3, 0:1, 0:1, nf_)
        integer*4 np, ni, ns 
        complex*16 :: c
    
        !---multiply by the kinetic energy operator
        zwave = KINETICW2(zwaveR_, +1)
    
        call getcorr(zwaveL_, zwave, zKin)
    
        zKin = zcorrinv_ .x. zKin
        res = zDet_ * getTrace(nf_,zKin)
    
        return
    end function getKin_test    
    !---------------------------------------------------------------------------
    function getQ0S_test(zwaveL_,zwaveR_,zDet_,zcorrinv_,nf_) result(res)
        !---measure the SU(4) contact terms, 
        !---including two-body, three-body and four-body terms.

        implicit none  
        complex*16 res(3)
        integer*4,intent(in) :: nf_ 
        complex*16,intent(in) :: zwaveL_(L3,0:1,0:1,nf_)
        complex*16,intent(in) :: zwaveR_(L3,0:1,0:1,nf_)
        complex*16,intent(in) :: zDet_ 
        complex*16,intent(in) :: zcorrinv_(nf_,nf_) 
        complex*16 :: zsmearL(L3, 0:1, 0:1, nf_), zsmearR(L3, 0:1, 0:1, nf_)
        complex*16 :: zrho(nf_, nf_), zsum(L3, nf_, nf_), c2, c3, c4
        integer*4 np, np1, np2, ns, ni, n
        real*8 :: fp(L3)

        !---wave function smearing
        zsmearL(:, :, :, :) = &
            SMEARW(zwaveL_(:, :, :, :), smearNL, -1)
        zsmearR(:, :, :, :) = &
            SMEARW(zwaveR_(:, :, :, :), smearNL, +1)

        !---density operator smearing
        !fp  = dexp(-qr**2 / (2 * Lam**2))
        zsum(:, :, :) = Z0
        do np2 = 1, nf_; do np1 = 1, nf_
            do ni = 0, 1; do ns = 0, 1
                zsum(:, np1, np2) = zsum(:, np1, np2) + &
                    zsmearL(:, ns, ni, np1) * zsmearR(:, ns, ni, np2)
            enddo; enddo
        enddo; enddo
        do np2 = 1, nf_; do np1 = 1, nf_
            zsum(:, np1, np2) = SMEAR(zsum(:, np1, np2), smearL1, smearL2, smearL3)
            !zsum(:, np1, np2) = zFourier(zFourier(zsum(:, np1, np2), +1) * fp, -1)        
        enddo; enddo
        
        c2 = Z0; c3 = Z0; c4 = Z0
        do n = 1, L3
            zrho(:, :) = zsum(n, :, :)
            zrho = zcorrinv_ .x. zrho

            c2 = c2 + zDet_ * getTrace(nf_,zrho, 2)
            if (cc3 /= zero) &
                c3 = c3 + zDet_ * getTrace(nf_,zrho, 3)
            if (cc4 /= zero) &
                c4 = c4 + zDet_ * getTrace(nf_,zrho, 4)
        enddo
    
        res(1) = c2 * (cc2/2)
        res(2) = c3 * (cc3/6)
        res(3) = c4 * (cc4/24)

        return
    end function getQ0S_test

    function getQSI_test(zwaveL_,zwaveR_,zDet_,zcorrinv_,nf_) result(res)
        !---measure the SU(4) contact terms, 
        !---including two-body, three-body and four-body terms.

        use variables
        implicit none 
        complex*16 res(3)
        integer*4,intent(in) :: nf_ 
        complex*16,intent(in) :: zwaveL_(L3,0:1,0:1,nf_)
        complex*16,intent(in) :: zwaveR_(L3,0:1,0:1,nf_)
        complex*16,intent(in) :: zDet_ 
        complex*16,intent(in) :: zcorrinv_(nf_,nf_) 
        complex*16 :: zsmearL(L3, 0:1, 0:1, nf_), zsmearR(L3, 0:1, 0:1, nf_)
        complex*16 :: zrho(nf_, nf_), zsum(L3, nf_, nf_), c2, c3, c4
        integer*4 np, np1, np2, ns, ni, nii, n, iso
        real*8 :: fp(L3)

        res = Z0

        !---wave function smearing
        zsmearL(:, :, :, :) = &
            SMEARW(zwaveL_(:, :, :, :), smearNL, -1)
        zsmearR(:, :, :, :) = &
            SMEARW(zwaveR_(:, :, :, :), smearNL, +1)

        do iso = 1, 3
           !---density operator smearing
           !fp  = dexp(-qr**2 / (2 * Lam**2))
           zsum(:, :, :) = Z0
           do np2 = 1, nf_; do np1 = 1, nf_
              do ni = 0, 1; do nii = 0, 1; do ns = 0, 1
                 zsum(:, np1, np2)                 &
                    = zsum(:, np1, np2)            &
                    + ztau2x2(ni, nii, iso)        &
                      * zsmearL(:,  ns, ni,  np1)  &
                       * zsmearR(:, ns, nii, np2)
              enddo; enddo; enddo
           enddo; enddo
           do np2 = 1, nf_; do np1 = 1, nf_
              zsum(:, np1, np2) = SMEAR(zsum(:, np1, np2), smearL1, smearL2, smearL3)
              !zsum(:, np1, np2) = zFourier(zFourier(zsum(:, np1, np2), +1) * fp, -1)        
           enddo; enddo
        
           c2 = Z0; c3 = Z0; c4 = Z0
           do n = 1, L3
              zrho(:, :) = zsum(n, :, :)
              zrho = zcorrinv_ .x. zrho

               c2 = c2 + zDet_ * getTrace(nf_,zrho, 2)
           enddo

           res(1) = res(1) + c2 * (cc2_I/2)
        
        enddo
        
        return

    end function getQSI_test

    !=======================subroutines for measurements===========================

    subroutine Measure_REW()
        ! energy(and other) measurements for nf_rew particle system  
        ! This have to be called after normal measurement where mconfig is updated. 
        implicit none 
        real*8 Pcm_phys, Pcm
        complex*16 cph(L3)
        integer*4 n, ni, ns, np
        integer*4 nfq_Vdiff  !--yhs 
        complex*16 temp, temp2(2),temp3(3),temp100(100)
        complex*16 temp20(20),temp4(4) 
        complex*16 tempnopr(nopr),temp3nopr(3*nopr) 
        

        call prepare_expect() 
        call accum(zDet_rew, zIyeBin_rew)
        call accum(zDet2_rew, zIyeBin2_rew)

        if(dble(zDet_rew) .gt. 0.D0)then
           call accum(zDet_rew, zIyeplusBin_rew)
        else
           call accum(zDet_rew, zIyemnusBin_rew)
        endif
        !---averaged phase <e^i\theta>
        zIabsBin_rew(1) = zIabsBin_rew(1) + zDet_rew
        zIabsBin_rew(2) = zIabsBin_rew(2) + abs(zDet_rew)
        !---tranfer matrix
        temp = getTMX_test(zwaveL_rew,zwaveR_rew,zDet_rew,zcorrmatrix_rew, nf_rew) ! A-1 case         
        call accum(temp, zTMXBin_rew)
        !---kinetic energy
        temp = getKin_test(zwaveL_rew,zwaveR_rew,zDet_rew,zcorrinv_rew,nf_rew)  ! A case   
        call accum(temp, zKinBin_rew, vrot_A, n_cubic, npar_A, n_reflc ) ! with projection
        !---LO contact term
        temp3 = getQ0S_test(zwaveL_rew,zwaveR_rew,zDet_rew,zcorrinv_rew,nf_rew) ! for A-1
        call accum(temp3, zQ0SBin_rew, vrot_A, n_cubic, npar_A, n_reflc )  
        temp3 = getQSI_test(zwaveL_rew,zwaveR_rew,zDet_rew,zcorrinv_rew,nf_rew) ! for A-1
        call accum(temp3, zQSIBin_rew, vrot_A, n_cubic, npar_A, n_reflc )
        if (Same_Random_Flag==1) call sgrnd(myseed + 10 * myid)
       
        return   
    end subroutine Measure_REW 
    !--------------------subroutines for statistics---------------------------------
    subroutine statistics_rew()
         
        use mpi
        use cubic
        use variables
        use auxiliary
        implicit none
        integer*4 i, k, n, nmin, ndim, ip, rep, npar, nrot, nf, &
        nw, nwmin, ierr, ni, nrow, ncol, nQ
        complex*16 zIabs(6), &
            zIye(dimB),zIyeplus(dimB),zIyemnus(dimB), &
            zVGIR0QN(dimB, nopr), zVGIR1QN(dimB, nopr), zVGIR2QN(dimB, nopr), &
            zQ0S(dimB, 3), zQSI(dimB,3), zVLOL(dimB, 100), zVLONL(dimB, 100), zQ0V(dimB), zISB(dimB, 2),  &
            zKin(dimB), zCou(dimB, 2), &
            zOPE(dimB, 2), zOPW(dimB), zVdiffN3LO(dimB, nopr), zVdiffN3LOGIR(dimB, 3*nopr), &
            zTBF(dimB, 20), &
            zPGF(dimB, dimR), zMDF(dimB, dimR),               &
            zGIR(dimB, 3), zDPF(dimB, dimR*2),                &
            zDPFplus(dimB, dimR*2), zDPFmnus(dimB, dimR*2),   &
            zIBP(dimB), zXXX(dimB, 9),                        &
            zA3N(dimB, 20), zrho(dimB, 4),                    &
            zdIy(dimB, 0:(nQord+2)), zdDP(dimB, dimR*2, 0:(nQord+2)), &
            zdIyplus(dimB, 0:(nQord+2)), zdIymnus(dimB, 0:(nQord+2)), &
            zdIyplus3N(dimB, 1:30),                           &
            zdIymnus3N(dimB, 1:30),                           &
            zdDPplus(dimB, dimR*2, 0:(nQord+2)),              &
            zdDPplus3N(dimB, dimR*2, 1:30),                   &
            zdDPmnus3N(dimB, dimR*2, 1:30),                   &
            zdDPmnus(dimB, dimR*2, 0:(nQord+2)), &
            zTMX(dimB), zdTM(dimB), &
            zRSQ(dimB, 2), zCNM(dimB, 3), zCLU(dimB, n_f), &
            zMXC(dimB, n_f), zMUL(dimB), &
            zDCF(dimB, dimR), &
            zImu(dimB, dimN), &
            zOcc(dimB, 50*2), &
            zCOR(dimB, n_f, n_f)
        complex*16 zIye0(dimB),zIye2(dimB)
        real*8 r(dimR), f(4), t, f0(dimR*2), f1(dimR*2), d1,      &
            f0plus(dimR*2), f1plus(dimR*2), f1plus3N(dimR*2),  &
            f0mnus(dimR*2), f1mnus(dimR*2), f1mnus3N(dimR*2),  &
            d0plus, d0mnus, d1plus, d1mnus, d1plus3N, d1mnus3N
        complex*16 hh0(n_c, n_c), hh(n_c, n_c), kk(n_c, n_c), v(n_c, n_c)
        real*8 x(n_c)
        real*8 cc0, cc1, epsln

        if (mconfig .eq. 0) return
    
        !---collect data from all threads, self-thread
        !---ommitted for jackknife
        call analyze(zIabsBin_rew / mconfig, zIabs)
        call analyze(zKinBin_rew / mconfig, zKin)
        call analyze(zQ0SBin_rew / mconfig, zQ0S)
        call analyze(zQSIBin_rew / mconfig, zQSI)
        call analyze(zQ0VBin_rew / mconfig, zQ0V)

        call analyze(zXXXBin_rew / mconfig, zXXX)

        call analyze(zIyeBin / mconfig, zIye0)
        call analyze(zIyeBin_rew / mconfig, zIye)
        call analyze(zIyeBin2_rew / mconfig, zIye2)

        call analyze(zIyeplusBin_rew / mconfig, zIyeplus)
        call analyze(zIyemnusBin_rew / mconfig, zIyemnus)

        call analyze(zTMXBin_rew / mconfig, zTMX)
        
        zphase_rew(1) = zIabs(1) / zIabs(2)
        zphase_rew(2) = zIabs(3) / zIabs(2)
        zphase_rew(3) = zIabs(5) / zIabs(2)
        
        !------diagonalize LO Hamiltonian in multi-channel. 
        !      Hamiltonian here only have K+ V_{0S}+OPW +V_{0V} 
        !----------
        ! Caution!  
        ! We do not assume multi-channel in REW calculation.   
        !----------
        if ((n_c > 1).or.(n_wf > 1)) stop 'Error: Re-Weighting is not available for multi-channel yet!'
        
        kk = reshape(zIye(:), (/n_c, n_c/))
        hh = reshape(zKin(:) + sum(zQ0S(:, :), 2) + zOPW(:) + zQ0V , (/n_c, n_c/))
        
        call gdiag(n_c, hh, kk, x, v, n_wf, 1.d-4) 
        
        do nf = 1, n_wf
            Qn_rew(nf)%v(:) = v(:, nf)
        enddo
 
        !______________________operator insertions____________________________________
        do nf = 1, n_wf

            Qn_rew(nf)%P = n_reflc
            Qn_rew(nf)%O = n_cubic

            !______________________energies in transfer matrix formalism___________________
            !
            call getexpect(zTMX, Qn_rew(nf)%v, Qn_rew(nf)%EM0) 
            Qn_rew(nf)%EM0 = -dlog( Qn_rew(nf)%EM0 ) / atovera

            call getexpect(zKin, Qn_rew(nf)%v, Qn_rew(nf)%Kin) 
            call getexpect(zOPE, Qn_rew(nf)%v, Qn_rew(nf)%OPE)
            call getexpect(zQ0V, Qn_rew(nf)%v, Qn_rew(nf)%Q0V)
            call getexpect(zQ0S, Qn_rew(nf)%v, Qn_rew(nf)%Q0S)
            call getexpect(zQSI, Qn_rew(nf)%v, Qn_rew(nf)%QSI)
    
            call getexpect(zXXX, Qn_rew(nf)%v, Qn_rew(nf)%XXX)

            !___________________ energies of the wave function Hamiltonian  ________________
            Qn_rew(nf)%EHWF = Qn_rew(nf)%Kin          &
                        + sum(Qn_rew(nf)%Q0S(:))  &
                        + sum(Qn_rew(nf)%QSI(:))  &
                        + Qn_rew(nf)%OPW

            !______________________energies in Hamiltonian formalism_______________________
            Qn_rew(nf)%EH0 = Qn_rew(nf)%Kin + sum(Qn_rew(nf)%OPE(:)) 

        enddo
        !----------------------------         
        ! Let XXX(3) = EHWF(A-1)-EHWF(A) 
        !     XXX(4) = Z_{A-1}/Z_A
        !     XXX(5) = Z_{A-2}/Z_A 
        !----------------------------  
        Qn_rew(1)%XXX(3) = Qn_rew(1)%EHWF - Qn(1)%EHWF  
        Qn_rew(1)%XXX(4) = zIye(1)/zIye0(1) 
        Qn_rew(1)%XXX(5) = zIye2(1)/zIye0(1) 

        !---re-collect data and do jackknife
        do nf = 1, n_wf
            call geterror(Qn_rew(nf)%EHWF, Qn_rew(nf)%EHWF_err)

            call geterror(Qn_rew(nf)%Kin, Qn_rew(nf)%Kin_err)
            call geterror(Qn_rew(nf)%Q0S, Qn_rew(nf)%Q0S_err)
            call geterror(Qn_rew(nf)%QSI, Qn_rew(nf)%QSI_err)
            call geterror(Qn_rew(nf)%Q0V, Qn_rew(nf)%Q0V_err)

            call geterror(Qn_rew(nf)%EM0, Qn_rew(nf)%EM0_err)

            call geterror(Qn_rew(nf)%XXX, Qn_rew(nf)%XXX_err)

        enddo
        return 
    end subroutine statistics_rew
    !--------------------subroutines for printer------------------------------------  
    subroutine print_rew() 
        !---print observables
        !   in Qn_rew 
        use variables
        use auxiliary
        use cubic
        
        implicit none
        integer*4 n, m, i, j, ni, nf, nQ, ord(4), madd(4), ih
        real*8 r(dimR), dtau, ptau, t, u, acpt
        logical mask(Lt), Rmask(dimR)
        character*3 TMtype
        character*40,parameter :: pad = repeat(' ', 40)
        integer*4,parameter :: pdst = 24
        integer*4 :: mm, nn, Lamb
        character*40 text
        
        write(pout,*) '=================RE-Weight Results================='
        write(pout,*) '!!! Results for A=',nf_rew,' Z=', np_rew, 'myid=', myid            
        !---print averaged sign
        write(pout,*) ' '
        write(pout,'(1x,a20,5x,f10.7,SP,f10.7,"j")') &
                'PHASE <e^{i*theta}>:', zphase_rew(1)
        write(pout,'(1x,a20,5x,f10.7,SP,f10.7,"j")') &
                'PHASE <e^{i*theta}>:', zphase_rew(2)
        write(pout,'(1x,a20,5x,f10.7,SP,f10.7,"j")') &
                'PHASE <e^{i*theta}>:', zphase_rew(3)
        
        !---print energies
        write(pout, *) '[[[[[[ ENERGIES ]]]]]]'
        
        write(pout, *) ' '
        write(pout, *) 'NON-PERTURBATIVE ENERGY <H0> (MeV):'
        call TableHead()
        call list('EHWF', Qn_rew(:)%EHWF, Qn_rew(:)%EHWF_err, 'MeV')
        call list('  > KIN', Qn_rew(:)%KIN, Qn_rew(:)%KIN_err, 'MeV')
        call list('  > Q0S_2N', Qn_rew(:)%Q0S(1), Qn_rew(:)%Q0S_err(1), 'MeV')
        call list('  > QSI', Qn_rew(:)%QSI(1), Qn_rew(:)%QSI_err(1), 'MeV')
        call list('  > Q0V', Qn_rew(:)%Q0V, Qn_rew(:)%Q0V_err, 'MeV')
        
        IF(1 .EQ. 1)THEN
        
            if (Mode .eq. 'NML' .or. Mode .eq. 'HOT') then
                
                write(pout, *) ' '
                write(pout, *) 'TRANSFER MATRIX ENERGY -1/a_t * <:exp(-a_t * H0):> (MeV):'
                call TableHead()
                call list('EM0', Qn_rew(:)%EM0, Qn_rew(:)%EM0_err, 'MeV')
        
           endif
        
        ENDIF        
        
        !---print testing measurement
        write(pout, *) ' '
        write(pout, *) '[[[[[[ TESTING MEASUREMENT ]]]]]]'
        write(pout, *) ' '
        call TableHead()
        do n = 1, 9
                call list('XXX('//i2c(n)//')', Qn_rew(:)%XXX(n), Qn_rew(:)%XXX_err(n),'l.u.')
        enddo
        
        write(pout,*) '=================End of RE-Weight Results================='
            
        return

        contains
            !-----------------------------------------------------------------------
            subroutine TableHead(xname)
                !---print a table head
        
                character(len=*),optional :: xname
                character*14 name
                integer*4 i
        
                name = 'name'
                if (present(xname)) name = xname
                write(pout, '(a1, 4x, a14, 31x, 100(7x, a19))') &
                    '*', name, ( 'w.f. ' // i2c(i) // '   ' // &
                    'w.f. ' // i2c(i) // ' err' // pad, i = 1, n_wf)
                write(pout, *) repeat('-', 51 + 2 * n_wf * 12)
        
            end subroutine TableHead
            !-----------------------------------------------------------------------
            subroutine makeTable(name, x, values, errors, unit, mask)
                !---print a table of function values y = f(x) + error
                !---values are in l.u. convert to n.u. if unit is present
        
                real*8,intent(in) :: x(:), &
                    values(size(x), n_wf), errors(size(x), n_wf)
                character(len=3),intent(in) :: name
                character(len=*),intent(in),optional :: unit
                logical,intent(in),optional :: mask(:)
                real*8 u
                character(len=100) ustring
                integer*4 i, n
                logical msk(size(x))
        
                ustring = ''
                if (present(unit)) ustring = unit
        
                msk = .true.
                if (present(mask)) msk = mask
        
                u = ParseUnit(ustring)
        
                m = 1
                do n = 1, size(x)
                    if (.not. msk(n)) cycle
                    write(pout, '(5x, a5, f12.4, 28x, 100e13.5)') name//i2c(m), x(n), & 
                        ((/values(n, i) * u, errors(n, i) * u/), i=1, n_wf)
                    m = m + 1
                enddo
        
                return
            end subroutine makeTable
            !-----------------------------------------------------------------------
            subroutine list(names, values, errors, units)
                !---print the entity with values and errors
                !---values are in l.u. convert to n.u. if unit is present
        
                character(len=*),intent(in) :: names
                real*8,intent(in) :: values(:), errors(:)
                character(len=*),intent(in) :: units
                real*8 u
                character(len=100) ustring
        
                ustring = units
        
                u = ParseUnit(ustring)
                ustring = ' [' // trim(ustring) // ']'
        
                write(pout, '(5x, a31, 9x, 100e13.5)') &
                    names // trim(ustring) // pad, & 
                    ((/values(i) * u, errors(i) * u/), i=1, n_wf)
        
            end subroutine list
            !----------------------------------------------------------------------        
    end subroutine print_rew 

    !===============================================================================
    subroutine getdiff_test(TMtype, zwaveL_, zwaveR_, zcorrinv_, &
        dVds0, dVdsI, dVdsv, dVdpn, cdVdcu, nf_)
    !---calculate the derivative -d/ds(log(|detZ|)),  Z is the nf 
    !---by nf correlation matrix Z_{i, j} = <i|j> and s is auxiliary fields.
    
    use auxiliary
    
    implicit none
    character*3 :: TMtype
    integer*4,intent(in) :: nf_
    complex*16,intent(in) :: zwaveL_(L3, 0:1, 0:1, nf_),  &
        zwaveR_(L3, 0:1, 0:1, nf_)
    complex*16,intent(in) :: zcorrinv_(nf_, nf_)
    real*8,intent(out) :: dVds0(L3), dVdsI(L3, 3),    &
                          dVdpn(L3, 3), dVdsv(L3, 0:3)
    complex*16,intent(out) :: cdVdcu(L3)
    real*8 cc, ccs, ccv, cpion, ccou
    complex*16 zdV(L3, 0:3, 0:3), zsum(L3), c, ccI
    complex*16 zsmearL(L3, 0:1, 0:1, nf_),  &
        zsmearR(L3, 0:1, 0:1, nf_)
    integer*4 np, np1, np2, ns, ns1, ns2, nss, ni, nii, ni1, ni2, iso, vec
    real*8 :: fGsmear(L3)
    
    !---different parameters for outer time slices
    select case (TMtype)
    case ('FUL')
        ccs = ccs_INN; ccv = ccv_INN
        cpion = cpion_INN; ccou = ccou_INN
    case ('OUT')
        ccs = ccs_OUT; ccv = ccv_OUT
        cpion = cpion_OUT; ccou = ccou_OUT
    case ('MID')
        ccs = ccs_H; ccv = ccv_H
        cpion = cpion_H; ccou = ccou_H
    end select
    
    dVds0 = zero; dVdsI = zero; dVdsv = zero; dVdpn = zero; cdVdcu = Z0
    
    if (atovera < 1.d-14) return
    
    !-------doubly smeared contact interaction------
    if (cc2 /= zero) then
        zsmearL(:, :, :, :) = &
            SMEARW(zwaveL_(:, :, :, :), smearNL, -1)
        zsmearR(:, :, :, :) = &
            SMEARW(zwaveR_(:, :, :, :), smearNL, +1)
        
        zsum(:) = Z0
        do np2 = 1, nf_; do np1 = 1, nf_
            do ni = 0, 1; do ns = 0, 1
                zsum(:) = zsum(:) + &
                    zcorrinv_(np2, np1) &
                    * zsmearL(:, ns, ni, np1) &
                    * zsmearR(:, ns, ni, np2)
            enddo; enddo
        enddo; enddo
    
        dVds0(:) = dVds0(:) -  dsqrt(-atovera * cc2) &
            * SMEAR( dble(zsum(:)), smearL1, smearL2, smearL3)
    
    endif
    
    if (cc2_I /= zero ) then
       ccI = Z1 * dsqrt(atovera * abs(cc2_I))
       if(cc2_I > zero) ccI = Zi * ccI
        
       !---wave function smearing
       zsmearL(:, :, :, :) = &
           SMEARW(zwaveL_(:, :, :, :), smearNL, -1)
       zsmearR(:, :, :, :) = &
           SMEARW(zwaveR_(:, :, :, :), smearNL, +1)
    
        do iso = 1, 3
           zsum(:) = Z0
           do np2 = 1, nf_; do np1 = 1, nf_
              do ni = 0, 1; do nii = 0, 1
                 do ns = 0, 1
                    zsum(:) = zsum(:) + zcorrinv_(np2, np1) * &
                              ztau2x2(ni, nii, iso) *        &
                              zsmearL(:, ns, ni , np1) *     &
                              zsmearR(:, ns, nii, np2)
                 enddo
              enddo; enddo
           enddo; enddo
    
           dVdsI(:, iso) = dVdsI(:, iso)                         &
             - dble(ccI * SMEAR( zsum(:), smearL1, smearL2, smearL3))
    
    
        enddo
    endif
    
    
    if (ccs /= zero .or. ccv /= zero) then
        !---wave function smearing
        zsmearL(:, :, :, :) = &
            SMEARW(zwaveL_(:, :, :, :), smearHO, -1)
        zsmearR(:, :, :, :) = &
            SMEARW(zwaveR_(:, :, :, :), smearHO, +1)
    
        do vec = 0, 3
    
            zsum(:) = Z0
            do np2 = 1, nf_; do np1 = 1, nf_
                do ni = 0, 1
                    do ns = 0, 1; do nss = 0, 1
                        zsum(:) = zsum(:) + zcorrinv_(np2, np1) * &
                            ztau2x2(ns, nss, vec) * &
                            zsmearL(:,  ns, ni, np1) * &
                            zsmearR(:, nss, ni, np2)
                    enddo; enddo
                enddo
            enddo; enddo
    
            if (vec .eq. 0) cc = ccs
            if (vec /=   0) cc = ccv
            dVdsv(:, vec) = dVdsv(:, vec) - &
                dsqrt(-atovera * cc) * dble(zsum(:))
    
        enddo
    endif
    
    !------1 pion exchange potential------
    if (cpion /= zero) then
        zdV = Z0
    
        do ni1 = 0, 1; do ni2 = 0, 1
            do ns1 = 0, 1; do ns2 = 0, 1
    
                zsum(:) = Z0
                do np1 = 1, nf_; do np2 = 1, nf_
                    zsum(:) = zsum(:) + &
                        zcorrinv_(np2, np1) &
                        * zwaveL_(:, ns1, ni1, np1) &
                        * zwaveR_(:, ns2, ni2, np2)
                enddo; enddo
    
                do iso = 1, 3; do vec = 1, 3
                    if ( ztau2x2(ns1, ns2, vec) .eq. Z0 .or. &
                        ztau2x2(ni1, ni2, iso) .eq. Z0 ) cycle
                    zdV(:, vec, iso) = zdV(:, vec, iso) + &
                        ztau2x2(ns1, ns2, vec) * &
                        ztau2x2(ni1, ni2, iso) * zsum(:)
                enddo; enddo
    
            enddo; enddo
        enddo; enddo;
    
        do iso = 1, 3; do vec = 1, 3
            dVdpn(:, iso) = dVdpn(:, iso) &
                +dsqrt(-atovera * cpion) &
                *FourierP2X(zFourierX2P(dble(zdV(:, vec, iso))) &
                *zp2dtp(:, :, :, vec)) 
        enddo; enddo
    endif
    
    if (ccou /= zero) then
        zsum(:) = Z0
        do np2 = 1, nf_; do np1 = 1, nf_
            do ns = 0, 1
                zsum(:) = zsum(:) + &
                    zcorrinv_(np2, np1) &
                    * zwaveL_(:, ns, 0, np1) &
                    * zwaveR_(:, ns, 0, np2)
            enddo
        enddo; enddo
        c = Z1 * dsqrt(atovera * abs(ccou))
        if (ccou > zero) c = Zi * c
        cdVdcu(:) = cdVdcu(:) - c * &
            zFourier(zFourier(zsum(:), +1) * czc2dtc, -1)
    endif
    
    return
    end subroutine getdiff_test 
    
end module reweight    
