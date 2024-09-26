module Mod_Household
    use PARAMS

    implicit none
    
    real(8), dimension(na, n_ofsh, ntheta, nkappa, nz, Tret) :: Vfun_ret, cfun_ret, afun_ret, dVfun_ret
    real(8), dimension(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork) :: Vfun, cfun, lfun, dVfun, afun, offshoring
    real(8) :: xsols(2, na, n_ofsh, ntheta, nkappa, nz, nxi, Twork)
    real(8), dimension(na, ntheta, nkappa, nz, nxi) :: dV_next, V_next
    
    !real(8), allocatable, dimension(:,:,:,:,:,:) :: Vfun_ret, cfun_ret, afun_ret, dVfun_ret 
    !real(8), allocatable, dimension(:,:,:,:,:,:,:) :: Vfun, cfun, lfun, dVfun, afun, offshoring
    !real(8), allocatable :: xsols(:,:,:,:,:,:,:,:)
    
    real(8) :: offshoring_ret(na, n_ofsh, ntheta, nkappa, nz)
    !real(8), allocatable, dimension(:,:,:,:,:) :: offshoring_ret, dV_next, V_next
    real(8) :: pretax_inc(na)

    real(8) :: c_mod, lambda_mod, nonlabinc_mod, w_mod, theta_mod, kappa_mod, ap_mod, a_mod, rcur_mod, psi_mod
    

contains

    !subroutine allocate_policy_fns()
    !    implicit none
    !    
    !    allocate(Vfun(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork))
    !    allocate(cfun(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork))
    !    allocate(afun(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork))
    !    allocate(lfun(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork))
    !    allocate(dVfun(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork))
    !    allocate(offshoring(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork))
    !    
    !    allocate(Vfun_ret(na, n_ofsh, ntheta, nkappa, nz, Tret))
    !    allocate(cfun_ret(na, n_ofsh, ntheta, nkappa, nz, Tret))
    !    allocate(afun_ret(na, n_ofsh, ntheta, nkappa, nz, Tret))
    !    allocate(dVfun_ret(na, n_ofsh, ntheta, nkappa, nz, Tret))
    !    
    !    allocate(xsols(2, na, n_ofsh, ntheta, nkappa, nz, nxi, Twork))
    !    
    !    allocate(dV_next(na, ntheta, nkappa, nz, nxi))
    !    allocate(V_next(na, ntheta, nkappa, nz, nxi))
    !    
    !    allocate(offshoring_ret(na, n_ofsh, ntheta, nkappa, nz))
    !    
    !end subroutine allocate_policy_fns
    
    
    subroutine SolveHH(save_res)
    use params
    use MyLinInterp 
    !use zbren_int
    use NEQNF_INT
    !use int_tictoc
    use root_module
    use focs_mod
    use glob_vars_mod
        logical :: save_res
        integer :: ia, nq, scp
        real(8) :: cons, lab
        real(8) :: vp, Vnext
        real(8) :: acutoff(n_ofsh), mcutoff(n_ofsh)
        real(8) :: m0
        real(8) :: gridm(na)
        real(8) :: gridm_ret(na, n_ofsh, ntheta, nkappa, nz)
        !real(8), allocatable :: gridm_ret(:,:,:,:,:)
        real(8) :: mtmp(na), ctmp(na), htmp(na), vtmp(na), atmp(na), aptmp(na)
        real(8) :: ctest(na), vtest(na)
        real(8) :: mdiff(na-1)
        logical :: mdiff_neg(na-1)
        integer :: nkinks
        real(8) :: aprime
        integer :: fne
        real(8) :: hlo, hhi, flo, fhi
        real(8) :: hrs, labinc
        integer :: IERSVR, IPACT, ISACT
        integer :: ikink, ij, ii
        integer :: fall(na), increase(na)
        real(8) :: VV
        real(8), allocatable :: CommonVt(:,:), CommonC(:,:), CommonH(:,:), CommonA(:,:)
        integer :: i1, i2
        integer :: idmax(na)
        !real(8) :: test_bc(nty,ns,na,J)
        !real(8) :: test_hrs(nty,ns,na,J)
        real(8) :: max_test_bc
        real(8) :: vals(2), testV, testH, testC, testA
        integer :: inds(2)
        integer :: jj
        real(8) :: fff_static, fff_helper_func1 !, fff_static_BL_h
        
        real(8) :: acur, mcur
        integer :: itheta, ikappa, iz, ithetap, ixi
        real(8) :: rcur, drcur
        real(8) :: gross_inc, net_inc, net_inc_ofsh, net_inc_nofsh
        real(8) :: vp_theta(ntheta), Vnext_theta(ntheta), Vnext_theta0(ntheta), Vnext0
        real(8) :: dDa(na), ra(na), dD(na, n_ofsh, ntheta, nkappa, nz), dD_
        real(8) :: beta_surv, Ucur
        real(8) :: xsol(2), xguess(2), fnorm, fvals(2), fnorm_test
        integer :: izp, ixip, ikappap
        real(8) :: hsol, fval
        integer :: ic
        real(8) :: ap_test(na)
        integer, parameter :: na_test = 50, nh_test = 50
        real(8) :: GridATest(na_test), GridHTest(nh_test) 
        real(8) :: da_test, dh_test
        integer :: ih
        logical :: solver_failed    
        real(8) :: asol, rtmp, ytmp !, net_inc_ofsh, net_inc_nofsh 

        external static_focs
        external static_focs_ofsh
        external static_focs_nofsh

        
        IERSVR = 0
        IPACT = 0 !1
        ISACT = 0
        call erset(IERSVR, IPACT, ISACT)  
        
        !allocate(gridm_ret(na, n_ofsh, ntheta, nkappa, nz))
        
        !da_test = 1.0d0/dble(na_test-1)
        !dh_test = 1.0d0/dble(nh_test-1)
        !GridATest = [(ia*da_test, ia=0,na_test-1)]
        !GridHTest = [(ih*dh_test, ih=0,nh_test-1)]
        !GridATest = 1.2d0 - 1.5d0/2d0 + GridATest*(1.5d0)
        !dh_test = 0.2d0*0.8d0
        !GridHTest = 0.2d0 - dh_test/2d0 + GridHTest*(dh_test)
        
        yb_cutoff = ( theta0*(1.0d0-theta1)/(1.0d0-tau_max) )**(1d0/theta1)
        
        ! Last period
        offshoring_ret = 0.0d0
        do jj = 1, n_ofsh
            do itheta = 1, ntheta
                do ikappa = 1, nkappa
                    do iz = 1, nz
                        !open(21, file='test.txt')
                        do ia = 1, na
                            acur = grida(ia)
                            rcur = rfunc(acur, thetas(itheta), kappas(ikappa))
                            drcur = Da_rfunc(acur, thetas(itheta), kappas(ikappa))
                            !ra(ia) = rcur
                            gross_inc = b_ret(iz) + rcur*acur
                            net_inc_ofsh = frac_ofsh*gross_inc + after_tax_income((1d0-frac_ofsh)*gross_inc) - psi_vals(jj)
                            net_inc_nofsh = after_tax_income(gross_inc)  
                            net_inc = max(net_inc_ofsh, net_inc_nofsh) 
                            mcur = acur + TrB + net_inc 
                            gridm_ret(ia, jj, itheta, ikappa, iz) = mcur
                            cons = mcur/(1.0d0+tauc) 
                            Vfun_ret(ia, jj, itheta, ikappa, iz, Tret) = ( cons**(1.0d0-sig1) - 1.0d0 )/(1.0d0-sig1) !U(cons,lab) 
                            cfun_ret(ia, jj, itheta, ikappa, iz, Tret) = cons
                            afun_ret(ia, jj, itheta, ikappa, iz, Tret) = 0d0
                            if (net_inc_ofsh >= net_inc_nofsh) then
                                offshoring_ret(ia, jj, itheta, ikappa, iz) = 1d0
                                dD_ = 1.0d0 + (frac_ofsh + D_after_tax_income((1d0-frac_ofsh)*gross_inc)*(1d0-frac_ofsh))*(rcur + acur*drcur)
                                dD(ia, jj, itheta, ikappa, iz) = dD_
                                dVfun_ret(ia, jj, itheta, ikappa, iz, Tret) = cons**(-sig1)*dD_
                            else
                                !dD_ = rcur*D_after_tax_income(gross_inc)
                                dD_ = (1.0d0 + D_after_tax_income(gross_inc)*(rcur + acur*drcur))
                                dD(ia, jj, itheta, ikappa, iz) = dD_
                                dVfun_ret(ia, jj, itheta, ikappa, iz, Tret) = cons**(-sig1)*dD_
                            end if    
                            !write(21, '(2f16.6)') grida(ia), DDa(ia)
                        end do
                        !close(21)
                    end do
                end do            
            end do
        end do

        
        do jj = 1, n_ofsh
            psi_mod = psi_vals(jj)
            psi_glob = psi_mod
            ! Retirement part
            do jc = Tret-1, 1, -1 
                do iz = 1, nz
                    beta_surv = bbeta*surv(jr+jc-1)
                    do itheta = 1, ntheta
                        do ikappa = 1, nkappa   
                            gridm = gridm_ret(:, jj, itheta, ikappa, iz)
                            dDa = dD(:, jj, itheta, ikappa, iz)
                            ! Invert Euler equation to find current consumption
                            do ia = 1, na    
                                do ithetap = 1, ntheta
                                    vp_theta(ithetap) = sum(pi_kappa*dVfun_ret(ia, jj, ithetap, :, iz, jc+1))
                                    Vnext_theta(ithetap) = sum(pi_kappa*Vfun_ret(ia, jj, ithetap, :, iz, jc+1))
                                end do
                                vp = beta_surv*sum(pi_theta(itheta,:)*vp_theta)                               
                                Vnext = sum(pi_theta(itheta,:)*Vnext_theta)
                                cons = vp**(-1.0d0/sig1)
                                mcur = grida(ia) + cons*(1.0d0+tauc)
                                mtmp(ia) = mcur
                                ctmp(ia) = cons
                                Ucur = ( cons**(1.0d0-sig1) - 1.0d0 )/(1.0d0-sig1) 
                                !vtmp(ia) = -1d0 / (Ucur + beta_surv*Vnext)
                                vtmp(ia) = ( Ucur + beta_surv*Vnext )
                                aptmp(ia) = grida(ia)
                            end do


                            if (save_res) then
                                if (jj == 1 .and. itheta == 1 .and. ikappa == 1 .and. iz == 1) then
                                    open(20, file='mtmp.txt')
                                    do ia=1,na
                                        write(20,'(5f16.6)') grida(ia), mtmp(ia), ctmp(ia), vtmp(ia), offshoring_ret(ia, jj, itheta, ikappa, iz) 
                                    end do
                                    close(20)
                                end if
                            end if

                            mdiff = mtmp(2:na) - mtmp(1:na-1)
                            nkinks = count( mdiff < 0d0 )

                            !call tic()
                            do ithetap = 1, ntheta
                                Vnext_theta0(ithetap) = sum(pi_kappa*Vfun_ret(1, jj, ithetap, :, iz, jc+1))
                            end do
                            Vnext0 = sum(pi_theta(itheta,:)*Vnext_theta0)  
                            
                            if (nkinks > 0) then
                                !print *, 'NKINKS = ', nkinks
                                increase(1) = 1
                                ikink = 0
                                ii = 2
                                do while (ii < na)
                                    if (mtmp(ii+1) < mtmp(ii) .or. (mtmp(ii) > mtmp(ii-1) .and. vtmp(ii) < vtmp(ii-1))) then
                                        ikink = ikink + 1
                                        fall(ikink) = ii
                                        ij = ii
                                        do while (mtmp(ij+1) < mtmp(ij))
                                            ij = ij + 1
                                        end do
                                        increase(ikink+1) = ij
                                        ii = ij
                                    end if
                                    ii = ii + 1
                                end do
                                nkinks = ikink + 1
                                fall(nkinks) = na
                                allocate(CommonVt(na,nkinks))
                                allocate(CommonC(na,nkinks))
                                allocate(CommonA(na,nkinks))
                                CommonVt = -10d0**8d0
                                do ia = 1, na
                                    mcur = gridm(ia)
                                    if (mcur <= mtmp(1)) then
                                        aprime = grida(1)
                                        cons = (mcur - aprime)/(1.0d0+tauc)
                                        Ucur = ( cons**(1.0d0-sig1) - 1.0d0 )/(1.0d0-sig1)
                                        !VV = -1d0/( Ucur + beta_surv*Vnext0 )
                                        VV = Ucur + beta_surv*Vnext0 
                                        CommonVt(ia,:) = VV
                                        CommonC(ia,:) = cons
                                        CommonA(ia,:) = aprime
                                    else 
                                        do ii = 1, nkinks
                                            if (mcur > mtmp(na)) then
                                                ! If the point on the exo grid "gridm" falls outside the last point of the endo grid "mtmp", extrapolate using the previous points
                                                CommonC(ia,ii) = LinInterp_1d(mcur,gridm(1:ia-1),CommonC(1:ia-1,ii),ia-1)   
                                                CommonA(ia,ii) = LinInterp_1d(mcur,gridm(1:ia-1),CommonA(1:ia-1,ii),ia-1) 
                                                CommonVt(ia,ii) = LinInterp_1d(mcur,gridm(1:ia-1),CommonVt(1:ia-1,ii),ia-1) 
                                            else
                                                i1 = increase(ii)
                                                i2 = fall(ii)
                                                if (mcur >= mtmp(i1) .and. mcur <= mtmp(i2)) then
                                                    call basefun(mtmp(i1:i2),i2-i1+1,mcur,vals,inds)
                                                    CommonC(ia,ii) = vals(1)*ctmp(i1-1+inds(1))+vals(2)*ctmp(i1-1+inds(2))
                                                    CommonVt(ia,ii) = vals(1)*vtmp(i1-1+inds(1))+vals(2)*vtmp(i1-1+inds(2))
                                                    CommonA(ia,ii) = vals(1)*aptmp(i1-1+inds(1))+vals(2)*aptmp(i1-1+inds(2))    
                                                end if
                                            end if
                                        end do
                                    end if
                                end do

                                idmax = maxloc(CommonVt, dim=2) 
                                do ia = 1, na
                                    cons = CommonC(ia,idmax(ia))
                                    cfun_ret(ia, jj, itheta, ikappa, iz, jc) = cons
                                    afun_ret(ia, jj, itheta, ikappa, iz, jc) = CommonA(ia,idmax(ia))
                                    Vfun_ret(ia, jj, itheta, ikappa, iz, jc) = CommonVt(ia,idmax(ia))
                                    dVfun_ret(ia, jj, itheta, ikappa, iz, jc) = cons**(-sig1)*dDa(ia)
                                end do
                                deallocate(CommonVt, CommonC, CommonA)
                            else
                                ! Interpolation step
                                do ia = 1, na
                                    mcur = gridm(ia)
                                    if (mcur <= mtmp(1)) then
                                        aprime = grida(1)
                                        cons = (mcur - aprime)/(1.0d0+tauc)
                                        Ucur = ( cons**(1.0d0-sig1) - 1.0d0 )/(1.0d0-sig1) 
                                        !VV = -1d0/( Ucur + beta_surv*Vnext0 )
                                        VV = Ucur + beta_surv*Vnext0 
                                    else
                                        call basefun(mtmp,na,mcur,vals,inds)

                                        !if (inds(2) == na) then
                                        !    print *, 'WARNING: grida too small?'
                                        !end if

                                        cons = vals(1)*ctmp(inds(1))+vals(2)*ctmp(inds(2))
                                        aprime = vals(1)*aptmp(inds(1))+vals(2)*aptmp(inds(2))
                                        !VV = -1d0/(vals(1)*vtmp(inds(1))+vals(2)*vtmp(inds(2)))
                                        VV = vals(1)*vtmp(inds(1))+vals(2)*vtmp(inds(2))

                                    end if
                                    afun_ret(ia, jj, itheta, ikappa, iz, jc) = aprime
                                    cfun_ret(ia, jj, itheta, ikappa, iz, jc) = cons
                                    Vfun_ret(ia, jj, itheta, ikappa, iz, jc) = VV
                                    dVfun_ret(ia, jj, itheta, ikappa, iz, jc) = cons**(-sig1)*dDa(ia)   
                                end do
                            end if

                            !call toc()
                            if (save_res == .true.) then
                                if (jj == 1 .and. itheta == 1 .and. ikappa == 1 .and. iz == 1) then
                                    open(21, file='gridm.txt')
                                    open(22, file='afun.txt')
                                    open(23, file='cfun.txt')
                                    open(24, file='Vfun.txt')
                                    open(25, file='dVfun.txt')
                                    open(26, file='test_plot_data.txt')
                                    do ia=1,na
                                        write(21,'(f26.6)') gridm(ia)
                                        write(22,'(f26.6)') afun_ret(ia, jj, itheta, ikappa, iz, jc)
                                        write(23,'(f26.6)') cfun_ret(ia, jj, itheta, ikappa, iz, jc)
                                        write(24,'(f26.6)') Vfun_ret(ia, jj, itheta, ikappa, iz, jc)
                                        write(25,'(f26.12)') dVfun_ret(ia, jj, itheta, ikappa, iz, jc)
                                    end do
                                    !write(26, '(2f16.6, i4)') acutoff(jj), mcutoff(jj), jc
                                    write(26, '(i4)') jc
                                    close(21); close(22); close(23); close(24); close(25) ; close(26)
                                end if
                            end if


                        end do
                    end do
                end do
            end do
        
        
            ! Working age part
            !dV_next(na, ntheta, nkappa, nz, nxi)
            !v_ret(na, n_ofsh, ntheta, nkappa, nz, Tret)
            do ixi = 1, nxi
                dV_next(:, :, :, :, ixi) = dVfun_ret(:, jj, :, :, :, 1)
                V_next(:, :, :, :, ixi) = Vfun_ret(:, jj, :, :, :, 1)
            end do
            
            !gridm = grida + TrB*(1d0+r*(1d0-tk))
            !print *, jr
            !do jc = jr-1, 1, -1
            do jc = Twork, 1, -1
                beta_surv = bbeta*surv(jc)
                do itheta = 1, ntheta
                    theta_mod = thetas(itheta)
                    do ikappa = 1, nkappa  
                        kappa_mod = kappas(ikappa)
                        do iz = 1, nz
                            do ixi = 1, nxi
                                w_mod = eta(iz)*xi(ixi)*ep(1,jc)*w
                                
                                solver_failed = .false.
                                ! Backward solution step
                                do ia = 1, na
                                    ap_mod = grida(ia)
                                    !nonlabinc_mod = grida(ia) + TrB
                                    vp = 0.0d0
                                    Vnext = 0.0d0
                                    do ithetap = 1, ntheta
                                        do ikappap = 1, nkappa
                                            do izp = 1, nz
                                                do ixip = 1, nxi
                                                    vp = vp + pi_xi(ixip)*pi(iz,izp)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*dV_next(ia, ithetap, ikappap, izp, ixip)
                                                    Vnext = Vnext + pi_xi(ixip)*pi(iz,izp)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*V_next(ia, ithetap, ikappap, izp, ixip)
                                                end do
                                            end do
                                        end do
                                    end do
                                    vp = beta_surv*vp
                                    cons = vp**(-1.0d0/sig1)
                                    c_mod = cons
                                    lambda_mod = cons**(-sig1)/(1d0+tauc)
                                    
                                    if (ia == 1) then
                                        xguess = [grida(ia), 0.5d0]
                                    else
                                        xguess = xsols(:, ia-1, jj, itheta, ikappa, iz, ixi, jc)
                                    end if
                    
                                    call static_focs(xguess, fvals, 2)
                                    call d_NEQNF (static_focs, xsol, xguess=xguess, fnorm=fnorm)
                                    !fnorm_test = sum(fvals**2d0)/2d0

                                    if (fnorm > 1d-6) then
                                    !if (ia == 103 .and. iz == 8) then
                                        !print *, 'WARNING: solver failed'
                                        solver_failed = .true.
                                        call d_NEQNF (static_focs_nofsh, xsol, xguess=xguess, fnorm=fnorm)
                                        if (fnorm > 1d-6) then
                                            print *, 'WARNING: solver failed in no-offshoring'
                                        end if
                                        asol = xsol(1)
                                        hsol = xsol(2)
                                        rtmp = rfunc(asol, theta_mod, kappa_mod)
                                        ytmp = rcur*asol + w_mod*hsol
                                        net_inc_ofsh = frac_ofsh*ytmp + after_tax_income((1d0-frac_ofsh)*ytmp) - psi_mod
                                        net_inc_nofsh = after_tax_income(ytmp)  
                                        if (net_inc_nofsh < net_inc_ofsh) then
                                            call d_NEQNF (static_focs_ofsh, xsol, xguess=xguess, fnorm=fnorm)    
                                            if (fnorm > 1d-6) then
                                                print *, 'WARNING: solver failed in offshoring'
                                            end if     
                                            asol = xsol(1)
                                            hsol = xsol(2)
                                            rtmp = rfunc(asol, theta_mod, kappa_mod)
                                            ytmp = rcur*asol + w_mod*hsol
                                            net_inc_ofsh = frac_ofsh*ytmp + after_tax_income((1d0-frac_ofsh)*ytmp) - psi_mod
                                            net_inc_nofsh = after_tax_income(ytmp)  
                                            if (net_inc_nofsh > net_inc_ofsh) then
                                                print *, 'WARNING: something is wrong with offshoring / no-offshoring cases'
                                            end if
                                        end if
                                    end if
                                    xsols(:, ia, jj, itheta, ikappa, iz, ixi, jc) = xsol
                                                                
                                    atmp(ia) = xsol(1)
                                    hrs = xsol(2)
                                    ctmp(ia) = cons
                                    htmp(ia) = hrs
                                    aptmp(ia) = grida(ia)
                                    vtmp(ia) = U(cons,hrs) + beta_surv*Vnext
                                    !vtmp(ia) = -1d0 / (U(cons,hrs) + beta_surv*Vnext)
                                end do
                                !open(20, file='mtmp.txt')
                                !do ia=1,na
                                    !write(20,'(5f20.9)') grida(ia,jj), mtmp(ia), ctmp(ia), htmp(ia), vtmp(ia)
                                !end do
                                !close(20)
                        
                                ! Interpolation step

                                mdiff = atmp(2:na) - atmp(1:na-1)
                                nkinks = count( mdiff < 0d0 )
                                !call tic()
                                
                                Vnext0 = 0.0d0
                                do ithetap = 1, ntheta
                                    do ikappap = 1, nkappa
                                        do izp = 1, nz
                                            do ixip = 1, nxi
                                                Vnext0 = Vnext0 + pi_xi(ixi)*pi(iz,izp)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*V_next(1, ithetap, ikappap, izp, ixip)
                                            end do
                                        end do
                                    end do
                                end do                                
                        
                                if (nkinks > 0) then
                                    ! Upper envelope
                                    !print *, 'NKINKS = ', nkinks
                                    increase(1) = 1
                                    ikink = 0
                                    ii = 2
                                    do while (ii < na)
                                        if (atmp(ii+1) < atmp(ii) .or. (atmp(ii) > atmp(ii-1) .and. vtmp(ii) < vtmp(ii-1))) then
                                            ikink = ikink + 1
                                            fall(ikink) = ii
                                            ij = ii
                                            do while (atmp(ij+1) < atmp(ij))
                                                ij = ij + 1
                                            end do
                                            increase(ikink+1) = ij
                                            ii = ij
                                        end if
                                        ii = ii + 1
                                    end do
                                    nkinks = ikink + 1
                                    fall(nkinks) = na
                                    allocate(CommonVt(na,nkinks))
                                    allocate(CommonH(na,nkinks))
                                    allocate(CommonC(na,nkinks))
                                    allocate(CommonA(na,nkinks))
                                    CommonVt = -10d0**8d0
                                    do ia = 1, na
                                        if (grida(ia) <= atmp(1)) then
                                            
                                            r_glob = rfunc(grida(ia), theta_mod, kappa_mod)
                                            a_glob = grida(ia)
                                            aprime = grida(1)
                                            ap_glob = aprime
                                            w_glob = w_mod                                            
                                            
                                            hlo = 0.000001d0
                                            hhi = 1.0d0
                                            
                                            flo = static_focs_BL_h(hlo)
                                            fhi = static_focs_BL_h(hhi)
                            
                                            if (flo*fhi <= 0d0) then
                                                call zeroin(static_focs_BL_h, hlo, hhi, tol=1d-8, xzero=hsol, fzero=fval, iflag=ic)
                                            else
                                                print *, 'WARNING: flo and fhi must have opposite signs!'
                                                do while (flo*fhi > 0d0)
                                                    hhi = hhi*2d0
                                                    fhi = static_focs_BL_h(hhi)
                                                end do
                                                call zeroin(static_focs_BL_h, hlo, hhi, tol=1d-8, xzero=hsol, fzero=fval, iflag=ic)
                                            end if                                            
                                                                                   
                                            cons = c_glob
                                            hrs = hsol
                                            
                                            VV = U(cons,hrs) + beta_surv*Vnext0
                                            CommonVt(ia,:) = VV
                                            CommonH(ia,:) = hrs
                                            CommonC(ia,:) = cons
                                            CommonA(ia,:) = aprime
                                        else 
                                            do ii = 1, nkinks
                                                if (grida(ia) > atmp(na)) then
                                                    ! If the point on the exo grid "gridm" falls outside the last point of the endo grid "mtmp", extrapolate using the previous points
                                                    CommonH(ia,ii) = LinInterp_1d(grida(ia),grida(1:ia-1),CommonH(1:ia-1,ii),ia-1)   
                                                    CommonC(ia,ii) = LinInterp_1d(grida(ia),grida(1:ia-1),CommonC(1:ia-1,ii),ia-1)   
                                                    CommonA(ia,ii) = LinInterp_1d(grida(ia),grida(1:ia-1),CommonA(1:ia-1,ii),ia-1) 
                                                    CommonVt(ia,ii) = LinInterp_1d(grida(ia),grida(1:ia-1),CommonVt(1:ia-1,ii),ia-1) 
                                                else
                                                    i1 = increase(ii)
                                                    i2 = fall(ii)
                                                    if ((grida(ia) >= atmp(i1) .and. grida(ia) <= atmp(i2)) .or. (ii == nkinks .and. grida(ia) > atmp(i2))) then                                                
                                                        call basefun(atmp(i1:i2),i2-i1+1,grida(ia),vals,inds)
                                                        CommonH(ia,ii) = vals(1)*htmp(i1-1+inds(1))+vals(2)*htmp(i1-1+inds(2))
                                                        CommonC(ia,ii) = vals(1)*ctmp(i1-1+inds(1))+vals(2)*ctmp(i1-1+inds(2))
                                                        !CommonVt(ia,ii) = -1d0/(vals(1)*vtmp(i1-1+inds(1))+vals(2)*vtmp(i1-1+inds(2)))
                                                        CommonVt(ia,ii) = (vals(1)*vtmp(i1-1+inds(1))+vals(2)*vtmp(i1-1+inds(2)))
                                                        CommonA(ia,ii) = vals(1)*aptmp(i1-1+inds(1))+vals(2)*aptmp(i1-1+inds(2))
                                                    end if
                                                end if
                                            end do
                                        end if
                                    end do
                            
                                    idmax = maxloc(CommonVt, dim=2) 
                                    do ia = 1, na
                                        acur = grida(ia)
                                        hrs = CommonH(ia,idmax(ia))
                                        lfun(ia,jj,itheta,ikappa,iz,ixi,jc) = hrs
                                        cons = CommonC(ia,idmax(ia))
                                        cfun(ia,jj,itheta,ikappa,iz,ixi,jc) = cons
                                        afun(ia,jj,itheta,ikappa,iz,ixi,jc) = CommonA(ia,idmax(ia))
                                        Vfun(ia,jj,itheta,ikappa,iz,ixi,jc) = CommonVt(ia,idmax(ia))

                                        rcur = rfunc(acur, thetas(itheta), kappas(ikappa))
                                        drcur = Da_rfunc(acur, thetas(itheta), kappas(ikappa))
                                        gross_inc = rcur*acur + w_mod*hrs
                                        net_inc_ofsh = frac_ofsh*gross_inc + after_tax_income((1d0-frac_ofsh)*gross_inc) - psi_vals(jj)
                                        net_inc_nofsh = after_tax_income(gross_inc)  
                                        net_inc = max(net_inc_ofsh, net_inc_nofsh) 
                                        if (net_inc_ofsh >= net_inc_nofsh) then
                                            offshoring(ia,jj,itheta,ikappa,iz,ixi,jc) = 1d0
                                            dD_ = 1.0d0 + (frac_ofsh + D_after_tax_income((1d0-frac_ofsh)*gross_inc)*(1d0-frac_ofsh))*(rcur + acur*drcur)
                                            dVfun(ia,jj,itheta,ikappa,iz,ixi,jc) = cons**(-sig1)*dD_
                                        else
                                            offshoring(ia,jj,itheta,ikappa,iz,ixi,jc) = 0d0
                                            dD_ = (1.0d0 + D_after_tax_income(gross_inc)*(rcur + acur*drcur))
                                            dVfun(ia,jj,itheta,ikappa,iz,ixi,jc) = cons**(-sig1)*dD_
                                        end if    
                                        ap_test(ia) = ap_bc(acur, cons, hrs, rcur, w_mod)
                                        !dVfun(ia,jj,itheta,ikappa,iz,ixi,jc) = marginal_utility(cons,hrs)
                                    end do
                                    deallocate(CommonVt, CommonH, CommonC, CommonA)

                                else                

                                    do ia = 1, na
                                        acur = grida(ia)
                                        rcur = rfunc(acur, thetas(itheta), kappas(ikappa))
                                        if (acur <= atmp(1)) then
                                            
                                            r_glob = rfunc(grida(ia), theta_mod, kappa_mod)
                                            a_glob = grida(ia)
                                            aprime = grida(1)
                                            ap_glob = aprime
                                            w_glob = w_mod                                            
                                            
                                            hlo = 0.000001d0
                                            hhi = 1.0d0
                                            
                                            flo = static_focs_BL_h(hlo)
                                            fhi = static_focs_BL_h(hhi)
                            
                                            if (flo*fhi <= 0d0) then
                                                call zeroin(static_focs_BL_h, hlo, hhi, tol=1d-8, xzero=hsol, fzero=fval, iflag=ic)
                                            else
                                                print *, 'WARNING: flo and fhi must have opposite signs!'
                                                do while (flo*fhi > 0d0)
                                                    hhi = hhi*2d0
                                                    fhi = static_focs_BL_h(hhi)
                                                end do
                                                call zeroin(static_focs_BL_h, hlo, hhi, tol=1d-8, xzero=hsol, fzero=fval, iflag=ic)
                                            end if                                            
                                                                                   
                                            cons = c_glob
                                            hrs = hsol                                            
                                            
                                            VV = U(cons,hrs) + beta_surv*Vnext0
                                        else
                                            !aprime = LinInterp_1d(gridm(ia),mtmp,atmp,na)
                                            !hrs = LinInterp_1d(gridm(ia),mtmp,htmp,na)
                                            !cons = LinInterp_1d(gridm(ia),mtmp,ctmp,na)
                                            !VV = -1d0/LinInterp_1d(gridm(ia),mtmp,vtmp,na)
                                            !!VV = LinInterp_1d(gridm(ia),mtmp,vtmp,na)
                                    
                                            call basefun(atmp,na,grida(ia),vals,inds)
                                            hrs = vals(1)*htmp(inds(1))+vals(2)*htmp(inds(2))
                                            cons = vals(1)*ctmp(inds(1))+vals(2)*ctmp(inds(2))
                                            !VV = -1d0/(vals(1)*vtmp(inds(1))+vals(2)*vtmp(inds(2)))
                                            VV = (vals(1)*vtmp(inds(1))+vals(2)*vtmp(inds(2)))
                                            aprime = vals(1)*aptmp(inds(1))+vals(2)*aptmp(inds(2))                                
                                    
                                        end if
                                        afun(ia,jj,itheta,ikappa,iz,ixi,jc) = aprime
                                        lfun(ia,jj,itheta,ikappa,iz,ixi,jc) = hrs
                                        cfun(ia,jj,itheta,ikappa,iz,ixi,jc) = cons
                                        Vfun(ia,jj,itheta,ikappa,iz,ixi,jc) = VV                             
                                        
                                        ap_test(ia) = ap_bc(acur, cons, hrs, rcur, w_mod)
                                        
                                        gross_inc = rcur*acur + w_mod*hrs
                                        pretax_inc(ia) = gross_inc
                                        net_inc_ofsh = frac_ofsh*gross_inc + after_tax_income((1d0-frac_ofsh)*gross_inc) - psi_vals(jj)
                                        net_inc_nofsh = after_tax_income(gross_inc)  
                                        net_inc = max(net_inc_ofsh, net_inc_nofsh) 
                                        if (net_inc_ofsh >= net_inc_nofsh) then
                                            offshoring(ia,jj,itheta,ikappa,iz,ixi,jc) = 1d0
                                            dD_ = 1.0d0 + (frac_ofsh + D_after_tax_income((1d0-frac_ofsh)*gross_inc)*(1d0-frac_ofsh))*(rcur + acur*drcur)
                                            dVfun(ia,jj,itheta,ikappa,iz,ixi,jc) = cons**(-sig1)*dD_
                                        else
                                            offshoring(ia,jj,itheta,ikappa,iz,ixi,jc) = 0d0
                                            dD_ = (1.0d0 + D_after_tax_income(gross_inc)*(rcur + acur*drcur))
                                            dVfun(ia,jj,itheta,ikappa,iz,ixi,jc) = cons**(-sig1)*dD_
                                        end if    
                                        
                                        
                                        !dVpun(ia,jj,itheta,ikappa,iz,ixi,jc) = marginal_utility(cons,hrs)
                                    end do
                                end if

                                !call toc()
                                if (save_res == .true.) then
                                    !if (solver_failed) then
                                    if (itheta == 2 .and. iz == 1) then
                                        open(21, file='gridm.txt')
                                        open(22, file='afun.txt')
                                        open(23, file='cfun.txt')
                                        open(24, file='Vfun.txt')
                                        open(25, file='dVfun.txt')
                                        open(26, file='hfun.txt')
                                        open(27, file='test_plot_data.txt')
                                        open(28, file='ap_test.txt')
                                        open(29, file='offshoring.txt')
                                        open(30, file='pretax_inc.txt')
                                        open(31, file='capinc.txt')
                                        open(32, file='labinc.txt')
                                        open(33, file='constax.txt')
                                        !open(28, file='test_bc.txt')
                                        do ia=1,na
                                            write(21,'(f16.6)') grida(ia)
                                            write(22,'(f16.6)') afun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                            write(23,'(f16.6)') cfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                            write(24,'(f16.6)') Vfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                            write(25,'(f16.6)') dVfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                            write(26,'(f16.6)') lfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                            !write(28,'(f16.6)') test_bc(tyc,sc,ia,jc)
                                            write(28, '(f16.6)') ap_test(ia)
                                            if (jc == Twork) then
                                                write(29, '(2f16.6)') offshoring(ia,jj,itheta,ikappa,iz,ixi,jc), offshoring_ret(ia, jj, itheta, ikappa, iz)
                                            else
                                                write(29, '(2f16.6)') offshoring(ia,jj,itheta,ikappa,iz,ixi,jc), offshoring(ia,jj,itheta,ikappa,iz,ixi,jc+1) 
                                            end if
                                            write(30, '(f16.6)') pretax_inc(ia)
                                            write(31, '(f16.6)') grida(ia)*rfunc(grida(ia), thetas(itheta), kappas(ikappa))
                                            write(32, '(f16.6)') lfun(ia,jj,itheta,ikappa,iz,ixi,jc)*w_mod
                                            write(33, '(f16.6)') cfun(ia,jj,itheta,ikappa,iz,ixi,jc)*tauc
                                        end do
                                        ! write(27, '(2f16.6, i4)') acutoff(jj), mcutoff(jj), jc
                                        write(27, '(i4)') jc

                                        close(21); close(22); close(23); close(24); close(25); close(26); close(27); close(28); close(29); close(30); close(31); close(32); close(33) 
                                    end if
                                end if
                        
                            end do
                        end do
                    end do
                end do
                
                !max_test_bc = maxval(abs(test_bc))
                
                !print *, max_test_bc
                
                do ixi = 1, nxi
                    dV_next(:, :, :, :, ixi) = dVfun(:, jj, :, :, :, ixi, jc)
                    V_next(:, :, :, :, ixi) = Vfun(:, jj, :, :, :, ixi, jc)
                end do                
                
            end do
            
            print *, 'jj = ', jj
            
        end do
        
        !deallocate(gridm_ret)
        
        print *, 'HH optimisation done'
    end subroutine SolveHH
        
    subroutine check_focs_on_grid_2d(focs, GridX, GridY, nx, ny)
        integer :: nx, ny
        real(8) :: GridX(nx)
        real(8) :: GridY(ny)
        integer :: ix, iy
        real(8) :: Errs(2, nx, ny)
        real(8) :: err(2), x(2)
        
        interface 
            subroutine focs(x, f, n)
                integer :: n
                real(8) :: x(n), f(n)
            end subroutine focs
        end interface
        
        open(71, file='GridX.txt')
        open(73, file='ErrsOnGrid.txt')
        do ix = 1, nx
            write(71, '(f16.8)') GridX(ix)
            if (ix == 1) open(72, file='GridY.txt')
            do iy = 1, ny
                x(1) = GridX(ix)
                x(2) = GridY(iy)
                call focs(x, err, 2)
                Errs(:, ix, iy) = err
                write(72, '(f16.8)') GridY(iy)
                write(73, '(2f16.8)') Errs(:, ix, iy)
            end do
            if (ix == 1) close(72)
        end do  
        close(71); close(73)
    
    end subroutine check_focs_on_grid_2d 
    
    function ap_bc(a, c, h, r, w) result(res)
        real(8), intent(in) :: a
        real(8), intent(in) :: c
        real(8), intent(in) :: h
        real(8), intent(in) :: r
        real(8), intent(in) :: w
        real(8) :: rcur, y, net_inc_ofsh, net_inc_nofsh, net_inc, dyh
        real(8) :: res
        
        y = r*a + w*h
        net_inc_ofsh = frac_ofsh*y + after_tax_income((1d0-frac_ofsh)*y) - psi_mod
        net_inc_nofsh = after_tax_income(y)  
        net_inc = max(net_inc_ofsh, net_inc_nofsh)     
    
        if (net_inc_ofsh >= net_inc_nofsh) then
            dyh = D_after_tax_income((1d0-frac_ofsh)*y)*(1d0-frac_ofsh)
        else
            dyh = D_after_tax_income(y)
        end if        
        
        res = net_inc + a + TrB - c*(1.0d0+tauc)
        
    end function ap_bc

end module Mod_Household
