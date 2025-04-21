module Mod_Distribution
    use PARAMS
    use Het_returns, only: nkappa, ntheta
    implicit none
    
    real(8), allocatable :: Phi(:, :, :, :, :, :, :)
    !real(8), allocatable :: Tot_income(:, :, :, :, :, :, :)
    real(8), allocatable :: Phi_ret(:, :, :, :, :, :)
    !dimension(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork) :: Phi
    integer, parameter :: na_small = 50
    real(8) :: grida_small(na_small)
    integer, parameter :: ntot_work = na*n_ofsh*ntheta*nkappa*nz*nxi*Twork
    integer, parameter :: ntot = ntot_work + na*n_ofsh*ntheta*nkappa*nz*Tret
    real(8), allocatable :: PhiAll_1d(:), TotInc_1d(:), Wealth_1d(:), Lorenz_x(:), Lorenz_y(:), &
        TotInc_1d_ord(:), Wealth_1d_ord(:), PE_share_in_totinc_1d(:), Offsh_1d(:), Wages_1d(:), Hours_1d(:), &
        Wages_1d_ord(:), Hours_1d_ord(:), LabInc_1d(:), LabInc_1d_ord(:), Lorenz_x_work(:), Lorenz_y_work(:), &
        PhiAll_1d_work(:), CapInc_1d(:), CapInc_1d_ord(:), PE_share_in_capinc_1d(:), PE_inc_1d(:), tmp_arr_sorted(:), &
        CapInc_share_in_totinc_1d(:), Omega_1d(:), Theta_1d(:), ZProd_1d(:), tmp_arr_work_sorted(:)
    integer, allocatable :: sort_key(:), sort_key_work(:)
    
contains
    
    subroutine init_distr()
        integer :: ierr
        allocate(Phi(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork), stat=ierr)
        if (ierr /= 0) then
            print *, 'Error allocating Phi'
            stop
        end if
        allocate(Phi_ret(na, n_ofsh, ntheta, nkappa, nz, Tret), stat=ierr)
        if (ierr /= 0) then
            print *, 'Error allocating Phi_ret'
            stop
        end if 
        
        !allocate(Tot_income(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork), stat=ierr)
        allocate(PhiAll_1d(ntot), stat=ierr)
        allocate(Wealth_1d(ntot), stat=ierr)
        allocate(TotInc_1d(ntot), stat=ierr)
        allocate(Lorenz_x(ntot), stat=ierr)
        allocate(Lorenz_y(ntot), stat=ierr)
        allocate(TotInc_1d_ord(ntot), stat=ierr)
        allocate(Wealth_1d_ord(ntot), stat=ierr)
        allocate(PE_share_in_totinc_1d(ntot), stat=ierr)
        allocate(Offsh_1d(ntot), stat=ierr)
        allocate(tmp_arr_sorted(ntot), stat=ierr)
        allocate(sort_key(ntot), stat=ierr)
        allocate(Wages_1d(ntot_work), stat=ierr)
        allocate(tmp_arr_work_sorted(ntot_work), stat=ierr)
        allocate(ZProd_1d(ntot_work), stat=ierr)
        allocate(Hours_1d(ntot_work), stat=ierr)
        allocate(Wages_1d_ord(ntot_work), stat=ierr)
        allocate(Hours_1d_ord(ntot_work), stat=ierr)
        allocate(LabInc_1d(ntot_work), stat=ierr)  
        allocate(LabInc_1d_ord(ntot_work), stat=ierr)
        allocate(Lorenz_x_work(ntot_work), stat=ierr)
        allocate(Lorenz_y_work(ntot_work), stat=ierr) 
        allocate(PhiAll_1d_work(ntot_work), stat=ierr)
        allocate(sort_key_work(ntot_work), stat=ierr)
        allocate(CapInc_1d(ntot), stat=ierr)
        allocate(CapInc_share_in_totinc_1d(ntot), stat=ierr)
        allocate(CapInc_1d_ord(ntot), stat=ierr)
        allocate(PE_share_in_capinc_1d(ntot), stat=ierr)
        allocate(PE_inc_1d(ntot), stat=ierr)
        allocate(Omega_1d(ntot), stat=ierr)
        allocate(Theta_1d(ntot), stat=ierr)
        
    end subroutine init_distr
    
    subroutine GetDistribution(save_res)
        ! THIS SUBROUTINE COMPUTES STEADY STATE DISTRIBUTION OF ASSETS
        use params
        use Het_returns
        use int_tictoc
        !use svrgp_int
        !use toolbox
        use MyLinInterp
        use Mod_Household, only: afun, afun_ret 
        use io, only: save_array, read_array
        use ogpf

        IMPLICIT NONE
        
        integer :: ia, itheta, ikappa, iz,  ixi, jj
        integer :: ithetap, ikappap, izp, ixip
        real(8) :: test 
        integer :: inds(2)
        real(8) :: vals(2)
        real(8) :: TT1, TT2
        !integer :: iunit_lc
        CHARACTER (LEN=*), PARAMETER :: outDir = "tmp/"
        INTEGER :: iunit_phi        
        integer :: get_phi, get_phi_ret
        real(8) :: help(Twork+Tret)
        type(gpf):: gp
        integer :: id_tmp
        logical :: save_res
        real(8) :: Phi_prev
        real(8) :: Phi_next(na, n_ofsh, ntheta, nkappa, nz, nxi)
        real(8) :: Phi_ret_next(na, n_ofsh, ntheta, nkappa, nz)
        
        ! Initialize Distribution by Computing Distribution for first Generation
        get_phi = 1
        get_phi_ret = 1
        
        print *, 'ntot = ', ntot
        
        if (get_phi == 1) then
            Phi=0.0d0
            !$OMP PARALLEL PRIVATE(jj, itheta, ikappa, iz, ixi)
            !$OMP DO collapse(5)
            !!$OMP DO collapse(5) SCHEDULE(DYNAMIC)
            do jj = 1, n_ofsh
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        do iz = 1, nz
                            do ixi = 1, nxi
                                Phi(1,jj,itheta,ikappa,iz,ixi,1) = pi_kappa(ikappa)*pi_xi(ixi)*pini(iz)*psi_prob(jj)*pi_theta_stat(itheta)   
                            end do
                        end do
                    end do
                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL
                
            test = sum(Phi(:,:,:,:,:,:,1))
            if (abs(test-1d0) > 1d-6) then
                print *, 'Error in initial distribution'
                print *, 1, test 
            end if
            
            
            ! Loop to find distributions for ages 2 to J
            !!call tic
            
            do jc=2,Twork
                Phi_next = 0d0
                !$OMP PARALLEL PRIVATE(jj, itheta, ikappa, iz, ixi, ia, vals, inds, TT1, TT2, izp, ithetap, ikappap, ixip, Phi_prev) SHARED(jc)
                !$OMP DO collapse(6) reduction(+:Phi_next)
                !!$OMP DO SCHEDULE(DYNAMIC)
                do jj=1,n_ofsh
                    do itheta=1,ntheta
                        do ikappa=1,nkappa
                            do iz=1,nz
                                do ixi=1,nxi
                                    do ia=1,na
                                        Phi_prev = Phi(ia,jj,itheta,ikappa,iz,ixi,jc-1)
                                        call basefun(grida(1:na),na,afun(ia,jj,itheta,ikappa,iz,ixi,jc-1),vals,inds)
                                        do izp=1,nz
                                            do ithetap=1,ntheta
                                                do ikappap=1,nkappa
                                                    do ixip=1,nxi
                                                        TT1 = vals(1)*pi(iz,izp)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)*pi_xi(ixip)
                                                        TT2 = vals(2)*pi(iz,izp)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)*pi_xi(ixip)
                                                        Phi_next(inds(1),jj,ithetap,ikappap,izp,ixip)=Phi_next(inds(1),jj,ithetap,ikappap,izp,ixip)+Phi_prev*TT1
                                                        Phi_next(inds(2),jj,ithetap,ikappap,izp,ixip)=Phi_next(inds(2),jj,ithetap,ikappap,izp,ixip)+Phi_prev*TT2
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL
                Phi(:,:,:,:,:,:,jc)=Phi_next
                test = sum(Phi(:,:,:,:,:,:,jc))
                if (abs(test-1d0) > 1d-6) then
                    print *, 'Error in distribution'
                    print *, jc, test
                end if
            end do
            
            ! call save_array(Phi, outDir // "Phi.bin")
        
        else 
            call read_array(Phi, outDir // "Phi.bin")
            test = sum(Phi(:,:,:,:,:,:,Twork))
            print *, jc, test            
        end if
        
        
        if (get_phi_ret == 1) then
            Phi_ret = 0d0
            Phi_ret_next = 0d0
            !print *, 'Retirement: '
            ! First period of retirement
            !$OMP PARALLEL PRIVATE(jj,itheta,ikappa,iz,ixi,ia,ithetap,ikappap,vals,inds,TT1,TT2,phi_prev)
            !$OMP DO collapse(6) reduction(+:Phi_ret_next)
            !!$OMP DO SCHEDULE(DYNAMIC)
            do jj = 1, n_ofsh
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        do iz = 1, nz
                            do ixi = 1, nxi
                                do ia = 1,na
                                    Phi_prev = Phi(ia,jj,itheta,ikappa,iz,ixi,Twork)
                                    call basefun(grida(1:na),na,afun(ia,jj,itheta,ikappa,iz,ixi,Twork),vals,inds)
                                    do ithetap=1,ntheta
                                        do ikappap=1,nkappa
                                            TT1 = vals(1)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                            TT2 = vals(2)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                            Phi_ret_next(inds(1),jj,ithetap,ikappap,iz) = Phi_ret_next(inds(1),jj,ithetap,ikappap,iz) + Phi_prev*TT1
                                            Phi_ret_next(inds(2),jj,ithetap,ikappap,iz) = Phi_ret_next(inds(2),jj,ithetap,ikappap,iz) + Phi_prev*TT2
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL
            Phi_ret(:,:,:,:,:,1) = Phi_ret_next
            test = sum(Phi_ret(:,:,:,:,:,1))
            if (abs(test-1d0) > 1d-6) then
                 print *, Twork+1, test
             end if
            
            ! All other retirement periods
            do jc = 2, Tret
                Phi_ret_next = 0d0
                !$OMP PARALLEL PRIVATE(jj,itheta,ikappa,iz,ixi,ia,vals,inds,TT1,TT2,ithetap,ikappap,phi_prev)
                !$OMP DO collapse(5) reduction(+:Phi_ret_next)
                !!$OMP DO SCHEDULE(DYNAMIC)
                do jj = 1, n_ofsh
                    do itheta = 1, ntheta
                        do ikappa = 1, nkappa
                            do iz = 1, nz
                                do ia = 1,na
                                    Phi_prev = Phi_ret(ia,jj,itheta,ikappa,iz,jc-1)
                                    call basefun(grida(1:na),na,afun_ret(ia,jj,itheta,ikappa,iz,jc-1),vals,inds)
                                    do ithetap=1,ntheta
                                        do ikappap=1,nkappa
                                            TT1 = vals(1)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                            TT2 = vals(2)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                            Phi_ret_next(inds(1),jj,ithetap,ikappap,iz) = Phi_ret_next(inds(1),jj,ithetap,ikappap,iz) + Phi_prev*TT1
                                            Phi_ret_next(inds(2),jj,ithetap,ikappap,iz) = Phi_ret_next(inds(2),jj,ithetap,ikappap,iz) + Phi_prev*TT2                                            
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do 
                !$OMP END DO
                !$OMP END PARALLEL        
                Phi_ret(:,:,:,:,:,jc) = Phi_ret_next
                test = sum(Phi_ret(:,:,:,:,:,jc))
                if (abs(test-1d0) > 1d-6) then
                     print *, Twork+jc, test
                end if
                !print *, Twork+jc, test            
            end do
            ! call save_array(Phi_ret, outDir // "Phi_ret.dat")
        else
            call read_array(Phi_ret, outDir // "Phi_ret.dat")    
        end if
        
    end subroutine GetDistribution
    
    subroutine SummarizeDistribution()
        use params
        use Taxes
        use Het_returns
        use Mod_Household, only: afun, cfun, lfun, afun_ret, cfun_ret, offshoring, offshoring_ret
        use MyLinInterp
        use toolbox
        use ogpf
        implicit none
        !!call tic
        integer :: iunit_lc
        real(8) :: rtmp
        real(8), dimension(Twork+Tret) :: abar, cbar, inctaxbar 
        real(8), dimension(Twork) :: lbar, labar        
        integer :: ii, ia, iz, ixi, itheta, ikappa, ithetap, ikappap
        integer :: jj
        integer :: get_phi, get_phi_ret
        real(8) :: help(Twork+Tret)
        type(gpf):: gp
        integer :: id_tmp
        logical :: save_res
        real(8) :: labinc_tmp, totinc_tmp, tax_tmp 
        real(8) :: aprime 
        real(8) :: TaxIncTest
        real(8) :: RaggTmp, RauxaggTmp, TrBnTmp, TaxETmp
        real(8) :: atmp, ctmp, aprimetmp, testbc_tmp, aftertaxtmp, aftertaxtest
        real(8) :: rptmp
        real(8) :: AsOffshore
        real(8) :: gini_inc, gini_wealth, gini_test
        real(8) :: LorenzFx(ntot), LorenzFx_work(ntot_work)
        real(8) :: frac_below
        real(8) :: PE_share, r_PE, pe_share_tmp
        real(8) :: cut_1pct, cut_0_1pct, cut_0_01pct, cut_10pct
        real(8) :: pe_share_1pct, pe_share_bot90pct, pe_share_0_1pct, pe_share_0_01pct
        real(8) :: offsh_share_1pct, offsh_share_0_1pct, offsh_share_0_01pct
        real(8) :: offsh_share
        real(8) :: d_test
        real(8) :: corr_hw
        real(8) :: test
        integer :: ntest
        integer :: i, iu
        real(8), allocatable :: LorenzTest(:,:), tmp_arr(:)
        real(8) :: tmp
        integer, parameter :: nplot = 1000
        real(8) :: plot_x(nplot), plot_y(nplot)
        real(8) :: dx_plot
        real(8) :: num_, den_
        logical:: make_plots
        logical:: report_ineq_stats
        integer:: iunit_ineq
        
        make_plots = .false.
        report_ineq_stats = .true.
        
        ! Find Stationary Distribution over Asset Holdings
        
        open(newunit=iunit_lc,file='LifeCycle.txt')
        open(newunit=iunit_ineq,file='InequalityStats.txt')
        RetAgg = 0d0
        TrBn = 0d0
        TaxE = 0d0
        As = 0d0
        AsOffshore = 0d0
        RAgg = 0d0
        RauxAgg = 0d0
        LAgg = 0d0
        TaxIncAboveYb = 0d0
        YfBelowYb = 0d0
        DBelowYb = 0d0
        TotOffshCost = 0d0
        PE_share = 0d0
        ii=1
        do jc = 1, Twork
            abar(jc) = sum(Phi(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,1:nxi,jc)*afun(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,1:nxi,jc))
            cbar(jc) = sum(Phi(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,1:nxi,jc)*cfun(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,1:nxi,jc))
            lbar(jc) = 0d0
            do iz = 1, nz
                do ixi = 1, nxi
                    lbar(jc) = lbar(jc) + sum(Phi(1:na,1:n_ofsh,1:ntheta,1:nkappa,iz,ixi,jc)*eta(iz)*xi(ixi)*ep(1,jc)*lfun(1:na,1:n_ofsh,1:ntheta,1:nkappa,iz,ixi,jc))
                end do
            end do
            labar(jc) = sum(Phi(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,1:nxi,jc)*lfun(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,1:nxi,jc))
            inctaxbar(jc) = 0d0
            do ia = 1, na
                atmp = grida(ia)
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        rtmp = rfunc(grida(ia), thetas(itheta), Kappas(ikappa))
                        r_pe = rfunc_PE(grida(ia), thetas(itheta), Kappas(ikappa))
                        do iz = 1, nz
                            do ixi = 1, nxi
                                do jj = 1, n_ofsh
                                    labinc_tmp = w*eta(iz)*xi(ixi)*ep(1,jc)*lfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    totinc_tmp = labinc_tmp + rtmp*grida(ia)
                                    pe_share_tmp = r_pe*grida(ia)/totinc_tmp
                                    PhiAll_1d(ii) = Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*mu(jc)
                                    PhiAll_1d_work(ii) = Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*mu(jc)
                                    Omega_1d(ii) = Omega(grida(ia), thetas(itheta))
                                    Theta_1d(ii) = thetas(itheta)
                                    TotInc_1d(ii) = totinc_tmp
                                    LabInc_1d(ii) = labinc_tmp
                                    CapInc_1d(ii) = rtmp*grida(ia)
                                    Wages_1d(ii) = w*eta(iz)*xi(ixi)*ep(1,jc)
                                    ZProd_1d(ii) = eta(iz)
                                    Hours_1d(ii) = lfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    ctmp = cfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    aprimetmp = afun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    Wealth_1d(ii) = grida(ia) 
                                    ii = ii + 1
                                    
                                    
                                    if (offshoring(ia,jj,itheta,ikappa,iz,ixi,jc) < 0.5d0) then
                                        tax_tmp = tax_income(totinc_tmp)
                                        YfBelowYb = YfBelowYb + min(yb_cutoff, totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        DBelowYb = DBelowYb + after_tax_income_aux(totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, totinc_tmp-yb_cutoff)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                    
                                        aftertaxtmp = totinc_tmp - tax_tmp  
                                        Offsh_1d(ii) = 0d0
                                    else
                                        tax_tmp = tax_income((1d0-frac_ofsh)*totinc_tmp)
                                        YfBelowYb = YfBelowYb + min(yb_cutoff, (1d0-frac_ofsh)*totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        DBelowYb = DBelowYb + after_tax_income_aux((1d0-frac_ofsh)*totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, ((1d0-frac_ofsh)*totinc_tmp-yb_cutoff))*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)                                       
                                        TotOffshCost = TotOffshCost + psi_vals(jj)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc) 
                                        
                                        aftertaxtmp = (1d0-frac_ofsh)*totinc_tmp - tax_tmp + frac_ofsh*totinc_tmp - psi_vals(jj)
                                        AsOffshore = AsOffshore + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*frac_ofsh*grida(ia)*Nu(jc)
                                        !OffshWealth_1d(ii) = frac_ofsh*grida(ia)
                                        Offsh_1d(ii) = frac_ofsh
                                    end if
                                    
                                    
                                    inctaxbar(jc) = inctaxbar(jc) + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*tax_tmp
                                    
                                    aprime = afun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    
                                    As = As + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*aprime*Nu(jc)
                                    PE_share = PE_share + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*pe_share_tmp*Nu(jc)
                                    PE_share_in_totinc_1d(ii) = pe_share_tmp
                                    CapInc_share_in_totinc_1d(ii) = rtmp*grida(ia)/totinc_tmp
                                    PE_share_in_capinc_1d(ii) = r_pe/rtmp !*grida(ia)
                                    LAgg = LAgg + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*eta(iz)*xi(ixi)*ep(1,jc)*lfun(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                    ! PE_inc_1d(ii) = r_pe*grida(ia)
                                    PE_inc_1d(ii) = r_pe*grida(ia)
                                    
                                    RaggTmp = 0d0
                                    RauxaggTmp = 0d0 
                                    TrBnTmp = 0d0
                                    TaxETmp = 0d0
                                    do ithetap = 1, ntheta
                                        do ikappap = 1, nkappa
                                            rptmp = rfunc(aprime, thetas(ithetap), Kappas(ikappap))
                                            RaggTmp = RaggTmp + pi_kappa(ikappap)*pi_theta(itheta,ithetap)*rptmp*aprime !*Nu(jc)    
                                            RauxaggTmp = RauxaggTmp + pi_kappa(ikappap)*pi_theta(itheta,ithetap)*(rptmp-rbar)*aprime !*Nu(jc)
                                                                                      
                                            TrBnTmp = TrBnTmp + pi_kappa(ikappap)*pi_theta(itheta,ithetap)*(aprime*(1d0 + rptmp)) !*Nu(jc)*(1d0-surv(jc))
                                            TaxETmp = TaxETmp + pi_kappa(ikappap)*pi_theta(itheta,ithetap)*tau_estate*max(aprime*(1d0 + rptmp) - a_estate, 0d0)
                                        end do
                                    end do
                                    RAgg = RAgg + RaggTmp*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                    RauxAgg = RauxAgg + RauxaggTmp*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                    TrBn = TrBn + TrBnTmp*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)*(1d0-surv(jc))
                                    TaxE = TaxE + TaxETmp*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)*(1d0-surv(jc))
                                end do    
                            end do
                        end do
                    end do
                end do
            end do
            write(iunit_lc, '(i0, 4f9.6)') jc, abar(jc), labar(jc), lbar(jc), cbar(jc) !, inctaxbar(jc)
        end do
        do jc = 1, Tret
            abar(Twork+jc) = sum(Phi_ret(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc)*afun_ret(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc))
            cbar(Twork+jc) = sum(Phi_ret(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc)*cfun_ret(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc))
            inctaxbar(Twork+jc) = 0d0
            
            do ia = 1, na
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        rtmp = rfunc(grida(ia), thetas(itheta), Kappas(ikappa))
                        r_pe = rfunc_PE(grida(ia), thetas(itheta), Kappas(ikappa))
                        do iz = 1, nz
                            totinc_tmp = b_ret(iz) + rtmp*grida(ia)
                            pe_share_tmp = r_pe*grida(ia)/totinc_tmp
                            do jj = 1, n_ofsh
                                
                                PhiAll_1d(ii) = Phi_ret(ia,jj,itheta,ikappa,iz,jc)*mu(Twork+jc)
                                TotInc_1d(ii) = totinc_tmp     
                                CapInc_1d(ii) = rtmp*grida(ia)
                                CapInc_share_in_totinc_1d(ii) = rtmp*grida(ia)/totinc_tmp
                                
                                if (offshoring_ret(ia, jj, itheta, ikappa, iz) < 0.5d0) then
                                    tax_tmp = tax_income(totinc_tmp)
                                    YfBelowYb = YfBelowYb + min(yb_cutoff, totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    DBelowYb = DBelowYb + after_tax_income_aux(totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, totinc_tmp-yb_cutoff)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    Offsh_1d(ii) = 0d0
                                else
                                    tax_tmp = tax_income((1d0-frac_ofsh)*totinc_tmp)
                                    YfBelowYb = YfBelowYb + min(yb_cutoff, (1d0-frac_ofsh)*totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    DBelowYb = DBelowYb + after_tax_income_aux((1d0-frac_ofsh)*totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, ((1d0-frac_ofsh)*totinc_tmp-yb_cutoff))*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)                                    
                                    TotOffshCost = TotOffshCost + psi_vals(jj)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    
                                    AsOffshore = AsOffshore + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*frac_ofsh*grida(ia)*Nu(Twork+jc)
                                    !OffshWealth_1d(ii) = frac_ofsh*grida(ia)
                                    Offsh_1d(ii) = frac_ofsh
                                end if
                                                                     
                                RetAgg = RetAgg + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*b_ret(iz)*Nu(Twork+jc)
                                
                                inctaxbar(Twork+jc) = inctaxbar(Twork+jc) + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*tax_tmp
                                
                                aprime = afun_ret(ia,jj,itheta,ikappa,iz,jc)
                                Wealth_1d(ii) = grida(ia) 
                                !Wealth_1d(ii) = aprimetmp
                                PE_share_in_totinc_1d(ii) = pe_share_tmp
                                PE_share_in_capinc_1d(ii) = r_pe/rtmp
                                PE_inc_1d(ii) = r_pe*grida(ia)
                                Omega_1d(ii) = Omega(grida(ia), thetas(itheta))
                                Theta_1d(ii) = thetas(itheta)
                                
                                ii = ii + 1 
                                
                                As = As + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*aprime*Nu(Twork+jc)
                                PE_share = PE_share + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*pe_share_tmp*Nu(Twork+jc)
                                
                                RaggTmp = 0d0
                                RauxaggTmp = 0d0     
                                TrBnTmp = 0d0
                                TaxETmp = 0d0                                
                                do ithetap = 1, ntheta
                                    do ikappap = 1, nkappa
                                        rptmp = rfunc(aprime, thetas(ithetap), Kappas(ikappap))
                                        RaggTmp = RaggTmp + pi_kappa(ikappap)*pi_theta(itheta,ithetap)*rptmp*aprime
                                        RauxaggTmp = RauxaggTmp + pi_kappa(ikappap)*pi_theta(itheta,ithetap)*(rptmp-rbar)*aprime
                                        
                                        TrBnTmp = TrBnTmp + pi_kappa(ikappap)*pi_theta(itheta,ithetap)*(aprime*(1d0 + rptmp))
                                        TaxETmp = TaxETmp + pi_kappa(ikappap)*pi_theta(itheta,ithetap)*tau_estate*max(aprime*(1d0 + rptmp) - a_estate, 0d0)
                                    end do
                                end do
                                RAgg = RAgg + RaggTmp*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                RauxAgg = RauxAgg + RauxaggTmp*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                TrBn = TrBn + TrBnTmp*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)*(1d0-surv(Twork+jc))
                                TaxE = TaxE + TaxETmp*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)*(1d0-surv(Twork+jc))
                            end do
                        end do
                    end do
                end do
            end do
            write(iunit_lc, '(i0, 4f9.6)') Twork+jc, abar(Twork+jc), 0d0, 0d0, cbar(Twork+jc) !, inctaxbar(Twork+jc)
        end do
        close(iunit_lc)
        
        TrBn = TrBn - TaxE
        TrBn = TrBn/sum(Nu)
        
        HrsAgg = sum( Nu(1:Twork)*lbar(1:Twork) ) / sum( Nu(1:Twork) )
        print *, 'Av Hours worked: ', HrsAgg
        
        AAgg = sum( Nu(1:J)*abar(1:J) )
        print *, 'Test A: ', AAgg - As
        LabS = sum( Nu(1:Twork)*lbar(1:Twork) )
        print *, 'Test L: ', LAgg - LabS
        CAgg = sum( Nu(1:J)*cbar(1:J) )
        print *, 'Share of wealth in offshore assets: ', AsOffshore/As
        TaxIncTest = sum( Nu(1:J)*inctaxbar(1:J) )
        TaxInc = TaxIncAboveYb + YfBelowYb - theta0*DBelowYb
        print *, 'Test TaxInc: ', TaxInc - TaxIncTest
        TaxC = tauc*CAgg
        TaxTot = TaxInc + TaxC + TaxE
        
        do ia = 1, na
            do jc = 1, Twork
                help(jc) = sum(Phi(ia,1:n_ofsh,1:ntheta,1:nkappa,1:nz,1:nxi,jc))
            end do
            do jc = 1, Tret
                help(Twork+jc) = sum(Phi_ret(ia,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc))
            end do
            ADis(ia) = sum(help*mu)
        end do
        
        if ( ADis(na) > 0.0 ) then
            print*,'Enlarge Grid', ADis(na)
        end if
        
        test = sum(ADis(1:na))
        print *, test
        if ( (test > 1.01) .or. (test<0.99 )) then
            print*,'Should equal 1', sum(ADis(1:na))
            print *, 'ADis in DISTRIB.f90.'
            !pause
        else
            if (make_plots == .true.) then
                ! Annotation: set title, xlabel, ylabel, line specification
                call gp%title('Distribution over assets')
                call gp%xlabel('assets')
                call gp%ylabel('f')

                !Call Plot to draw a vector against a vector of data
                !The last argument defines the line specification
                !call gp%plot(grida,ADis,'with linespoints lt 2 pt 4')    
                id_tmp = maxloc(grida, dim=1, mask=(ADis > 1d-6))
                !call gp%plot(grida(1:id_tmp),ADis(1:id_tmp),'with lines') 
            end if
        end if
        
        
        ! Earnings and wealth distribution
        test = sum(PhiAll_1d) 
        print *, 'sum(PhiAll_1)= ', test
        PhiAll_1d = PhiAll_1d/test
        
        test = sum(PhiAll_1d_work)
        print *, 'sum(PhiAll_1d_work)= ', test
        PhiAll_1d_work = PhiAll_1d_work/test
        
        ! Total Income
        write(iunit_ineq, '(a)') 'Total Income Distribution: *************************************************'
        call lorenz_s(PhiAll_1d, TotInc_1d, LorenzFx, TotInc_1d_ord, Lorenz_x, Lorenz_y, gini_inc, sort_key)
        write(iunit_ineq, '(a60,f10.4)') 'Income Gini = ', gini_inc
        !call gp%title('Lorenz curve (total income)')
        !call gp%xlabel('Population share')
        !call gp%ylabel('Income share')
        !call gp%plot(Lorenz_x, Lorenz_y, 'with lines') 
        
        ! cut_10pct = LinInterp_1d(0.9d0,Lorenz_x,TotInc_1d_ord,ntot)
        ! cut_1pct = LinInterp_1d(0.99d0,Lorenz_x,TotInc_1d_ord,ntot)
        ! cut_0_1pct = LinInterp_1d(0.999d0,Lorenz_x,TotInc_1d_ord,ntot)
        ! cut_0_01pct = LinInterp_1d(0.9999d0,Lorenz_x,TotInc_1d_ord,ntot)
        
        tmp_arr_sorted = Offsh_1d(sort_key)
        !offsh_share_1pct = sum( Offsh_1d*LorenzFx, mask=(TotInc_1d_ord >= cut_1pct) ) / 0.01d0
        offsh_share = sum(tmp_arr_sorted*LorenzFx, mask=(Lorenz_x >= 0.99d0) ) / 0.01d0
        write(iunit_ineq,'(a60,f10.4,a3)') 'Share of offshored income among top 1% in tot. income = ', offsh_share*100.0d0, '%'
        ! offsh_share_0_1pct = sum( Offsh_1d*LorenzFx, mask=(TotInc_1d_ord >= cut_0_1pct) ) / 0.001d0   
        offsh_share = sum(tmp_arr_sorted*LorenzFx, mask=(Lorenz_x >= 0.999d0) ) / 0.001d0
        write(iunit_ineq,'(a60,f10.4,a3)') 'Share of offshored income among top 0.1% in tot. income = ', offsh_share*100.0d0, '%'
        ! offsh_share_0_01pct = sum( Offsh_1d*LorenzFx, mask=(TotInc_1d_ord >= cut_0_01pct) ) / 0.0001d0   
        offsh_share = sum(tmp_arr_sorted*LorenzFx, mask=(Lorenz_x >= 0.9999d0) ) / 0.0001d0
        write(iunit_ineq,'(a60,f10.4,a3)') 'Share of offshored income among top 0.01% in tot. income = ', offsh_share*100.0d0, '%'        
        
        
        test = sum( LorenzFx, mask=(Lorenz_x > 0.99d0) ) 
        d_test = abs(test - 0.01d0)/0.01d0
        if (d_test > 1d-2) then
            write(iunit_ineq,'(a,f12.5,a,f12.5)') 'Error in the mass of income share of top 1%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if        
        test = sum( LorenzFx, mask=(Lorenz_x > 0.999d0) ) 
        d_test = abs(test - 0.001d0)/0.001d0
        if (d_test > 1d-2) then
            write(iunit_ineq,'(a,f12.5,a,f12.5)') 'Error in the mass of income share of top 0.1%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if        
        test = sum( LorenzFx, mask=(Lorenz_x > 0.9999d0) ) 
        d_test = abs(test - 0.0001d0)/0.0001d0
        if (d_test > 1d-2) then
            write(iunit_ineq,'(a,f12.5,a,f12.5)') 'Error in the mass of income share of top 0.01%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if        

        tmp_arr_sorted = PE_share_in_totinc_1d(sort_key)
        pe_share = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x <= 0.90d0) )/0.90d0
        write(iunit_ineq,'(a60,f10.4)') 'PE/TotInc among bottom 90% in tot. income = ', pe_share
        pe_share = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.99d0) )/0.01d0
        write(iunit_ineq,'(a60,f10.4)') 'PE/TotInc among top 1% in tot. income = ', pe_share
        pe_share = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.999d0) )/0.001d0
        write(iunit_ineq,'(a60,f10.4)') 'PE/TotInc among top 0.1% in tot. income = ', pe_share
        pe_share = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9999d0) )/0.0001d0
        write(iunit_ineq,'(a60,f10.4)') 'PE/TotInc among top 0.01% in tot. income = ', pe_share
        
        tmp_arr_sorted = CapInc_share_in_totinc_1d(sort_key)
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x <= 0.90d0) )/0.90d0
        write(iunit_ineq,'(a60,f10.4)') 'CapInc/TotInc among bottom 90% in tot. income = ', tmp
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.99d0) )/0.01d0
        write(iunit_ineq,'(a60,f10.4)') 'CapInc/TotInc among top 1% of total income = ', tmp
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.999d0) )/0.001d0
        write(iunit_ineq,'(a60,f10.4)') 'CapInc/TotInc among top 0.1% of total income = ', tmp
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9999d0) )/0.0001d0
        write(iunit_ineq,'(a60,f10.4)') 'CapInc/TotInc among top 0.01% of total income = ', tmp
        
        frac_below = LinInterp_1d(0.9d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Income share of top 10% = ', (1d0-frac_below)*100.0d0, '%'        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Income share of top 1% = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Income share of top 0.1% = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Income share of top 0.01% = ', (1d0-frac_below)*100.0d0, '%'          
                
        ! Wealth
        write(iunit_ineq, '(a)') 'Wealth Distribution: *****************************************************'
        call lorenz_s(PhiAll_1d, Wealth_1d, LorenzFx, Wealth_1d_ord, Lorenz_x, Lorenz_y, gini_wealth, sort_key)
        write(iunit_ineq,'(a60,f10.4)') 'Wealth Gini = ', gini_wealth
        
        cut_10pct= LinInterp_1d(0.9d0, Lorenz_x, Wealth_1d_ord, ntot)
        write(iunit_ineq,'(a60,f10.4)') 'Cutoff for 90 percentile of wealth distribution = ', cut_10pct
        
        tmp_arr_sorted = PE_share_in_capinc_1d(sort_key)
        
        if (make_plots == .true.) then
            dx_plot = 1d0/(nplot-1)
            plot_x = [ (dx_plot*(i-1), i=1,nplot) ]
            plot_y = [ (LinInterp_1d(plot_x(i),Lorenz_x,tmp_arr_sorted,ntot), i=1,nplot) ]
            !plot_y = [ (LinInterp_1d(plot_x(i),Lorenz_x,Wealth_1d_ord,ntot), i=1,nplot) ]
            call gp%title('PE share along the wealth distribution')
            call gp%xlabel('Wealth share')
            !call gp%ylabel('Wealth share')
            call gp%ylabel('PE share')
            !call gp%plot(Lorenz_x, Lorenz_y, 'with lines')  
            call gp%plot(plot_x, plot_y, 'with lines') 
        
            plot_y = [ (LinInterp_1d(plot_x(i),Lorenz_x,Wealth_1d_ord,ntot), i=1,nplot) ]
            call gp%title('Wealth level along the wealth distribution')
            call gp%xlabel('Wealth share')
            !call gp%ylabel('Wealth share')
            call gp%ylabel('Wealth')
            !call gp%plot(Lorenz_x, Lorenz_y, 'with lines')  
            call gp%plot(plot_x, plot_y, 'with lines') 
        
            tmp_arr_sorted = Omega_1d(sort_key)
            !plot_y = [ (LinInterp_1d(plot_x(i),Lorenz_x,tmp_arr_sorted,ntot), i=1,nplot) ]
            do i = 1, nplot
                plot_y(i) = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > plot_x(i) ) )/(1d0-plot_x(i))
            end do        
            call gp%title('Omega along the wealth distribution')
            call gp%xlabel('Wealth share')
            !call gp%ylabel('Wealth share')
            call gp%ylabel('E(Omega), w > wbar')
            !call gp%plot(Lorenz_x, Lorenz_y, 'with lines')  
            call gp%plot(plot_x, plot_y, 'with lines')     
        end if
        
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x <= 0.90d0) )/0.90d0
        write(iunit_ineq,'(a60,f10.4)') 'Omega among bottom 90% in wealth = ', tmp
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.90d0) )/0.10d0
        write(iunit_ineq,'(a60,f10.4)') 'Omega among top 10% in wealth = ', tmp 
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.95d0) )/0.05d0
        write(iunit_ineq,'(a60,f10.4)') 'Omega among top 5% in wealth = ', tmp         
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.99d0) )/0.01d0
        write(iunit_ineq,'(a60,f10.4)') 'Omega among top 1% of wealth = ', tmp
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.999d0) )/0.001d0
        write(iunit_ineq,'(a60,f10.4)') 'Omega among top 0.1% of wealth = ', tmp
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9995d0) )/0.0005d0
        write(iunit_ineq,'(a60,f10.4)') 'Omega among top 0.05% of wealth = ', tmp          
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9999d0) )/0.0001d0
        write(iunit_ineq,'(a60,f10.4)') 'Omega among top 0.01% of wealth = ', tmp        
        
        
        tmp_arr_sorted = Theta_1d(sort_key)
        if (make_plots == .true.) then
            !plot_y = [ (LinInterp_1d(plot_x(i),Lorenz_x,tmp_arr_sorted,ntot), i=1,nplot) ]
            do i = 1, nplot
                plot_y(i) = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > plot_x(i) ) )/(1d0-plot_x(i))
            end do
            call gp%title('Theta along the wealth distribution')
            call gp%xlabel('Wealth share')
            !call gp%ylabel('Wealth share')
            call gp%ylabel('E(Theta), w > wbar')
            !call gp%plot(Lorenz_x, Lorenz_y, 'with lines')  
            call gp%plot(plot_x, plot_y, 'with lines') 
        end if
        
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x <= 0.90d0) )/0.90d0
        write(iunit_ineq,'(a60,f10.4)') 'Theta among bottom 90% in wealth = ', tmp
        !tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.90d0) )/0.10d0
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.90d0) )/(1d0-0.90d0)
        write(iunit_ineq,'(a60,f10.4)') 'Theta among top 10% in wealth = ', tmp 
        !tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.95d0) )/0.05d0
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.95d0) )/(1d0-0.95d0)
        write(iunit_ineq,'(a60,f10.4)') 'Theta among top 5% in wealth = ', tmp         
        !tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.99d0) )/0.01d0
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.99d0) )/(1d0-0.99d0)
        write(iunit_ineq,'(a60,f10.4)') 'Theta among top 1% of wealth = ', tmp
        !tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.999d0) )/0.001d0
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.999d0) )/(1d0-0.999d0)
        write(iunit_ineq,'(a60,f10.4)') 'Theta among top 0.1% of wealth = ', tmp
        !tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9995d0) )/0.0005d0
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9995d0) )/(1d0-0.9995d0)
        write(iunit_ineq,'(a60,f10.4)') 'Theta among top 0.01% of wealth = ', tmp        
        !tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9999d0) )/0.0001d0
        tmp = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9999d0) )/(1d0-0.9999d0)
        write(iunit_ineq,'(a60,f10.4)') 'Theta among top 0.01% of wealth = ', tmp
        
        
        
        !call gp%title('Lorenz curve (wealth)')
        !call gp%xlabel('Population share')
        !!call gp%ylabel('Wealth share')
        !call gp%ylabel('PE share')
        !!call gp%plot(Lorenz_x, Lorenz_y, 'with lines')  
        !call gp%plot(Lorenz_x, tmp_arr_sorted, 'with lines') 
        
        test = sum( LorenzFx, mask=(Lorenz_x > 0.99d0) ) 
        d_test = abs(test - 0.01d0)/0.01d0
        if (d_test > 1d-2) then
            write(iunit_ineq,'(a,f12.5,a,f12.5)') 'Error in the mass of wealth share of top 1%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if        
        test = sum( LorenzFx, mask=(Lorenz_x > 0.999d0) ) 
        d_test = abs(test - 0.001d0)/0.001d0
        if (d_test > 1d-2) then
            write(iunit_ineq,'(a,f12.5,a,f12.5)') 'Error in the mass of wealth share of top 0.1%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if        
        test = sum( LorenzFx, mask=(Lorenz_x > 0.9999d0) ) 
        d_test = abs(test - 0.0001d0)/0.0001d0
        if (d_test > 1d-2) then
            write(iunit_ineq,'(a,f12.5,a,f12.5)') 'Error in the mass of wealth share of top 0.01%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if        

        frac_below = LinInterp_1d(0.9d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Wealth share of top 10% = ', (1d0-frac_below)*100.0d0, '%'        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Wealth share of top 1% individuals = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Wealth share of top 0.1% individuals = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Wealth share of top 0.01% individuals = ', (1d0-frac_below)*100.0d0, '%'    
        
        
        write(iunit_ineq, '(a)') 'Labor income distribution: *********************************************'
        call lorenz_s(PhiAll_1d_work, LabInc_1d, LorenzFx_work, LabInc_1d_ord, Lorenz_x_work, Lorenz_y_work, gini_wealth, sort_key_work)
        frac_below = LinInterp_1d(0.9d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Labor inc share of top 10% of workers = ', (1d0-frac_below)*100.0d0, '%'        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Labor inc share of top 1% of workers = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Labor inc share of top 0.1% of workers = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Labor inc share of top 0.01% of workers = ', (1d0-frac_below)*100.0d0, '%'    
        
        write(iunit_ineq, '(a)') 'Wages distribution: *********************************************'
        call lorenz_s(PhiAll_1d_work, Wages_1d, LorenzFx_work, Wages_1d_ord, Lorenz_x_work, Lorenz_y_work, gini_wealth, sort_key_work)
        frac_below = LinInterp_1d(0.9d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Wages share of top 10% of workers = ', (1d0-frac_below)*100.0d0, '%'        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Wages share of top 1% of workers = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Wages share of top 0.1% of workers = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Wages share of top 0.01% of workers = ', (1d0-frac_below)*100.0d0, '%'    

        write(iunit_ineq, '(a)') 'Z-productivity distribution: *************************************'
        call lorenz_s(PhiAll_1d_work, ZProd_1d, LorenzFx_work, tmp_arr_work_sorted, Lorenz_x_work, Lorenz_y_work, gini_wealth, sort_key_work)
        frac_below = LinInterp_1d(0.9d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'ZProd share of top 10% of workers = ', (1d0-frac_below)*100.0d0, '%'        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'ZProd of top 1% of workers = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'ZProd of top 0.1% of workers = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'ZProd of top 0.01% of workers = ', (1d0-frac_below)*100.0d0, '%'    
        
        write(iunit_ineq, '(a)') 'Hours distribution: *********************************************'
        call lorenz_s(PhiAll_1d_work, Hours_1d, LorenzFx_work, Hours_1d_ord, Lorenz_x_work, Lorenz_y_work, gini_wealth, sort_key_work)
        frac_below = LinInterp_1d(0.9d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Hours share of top 10% of workers = ', (1d0-frac_below)*100.0d0, '%'        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Hours share of top 1% of workers = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Hours share of top 0.1% of workers = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x_work,Lorenz_y_work,ntot_work)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Hours share of top 0.01% of workers = ', (1d0-frac_below)*100.0d0, '%'   

        corr_hw = covar(Hours_1d, Wages_1d, PhiAll_1d_work)/(stdev(Hours_1d, PhiAll_1d_work)*stdev(Wages_1d, PhiAll_1d_work))
        write(iunit_ineq,'(a60,f10.4)') 'Corr between hours and wages = ', corr_hw
        
        write(iunit_ineq, '(a)') 'Capital income distribution: *********************************************'
        call lorenz_s(PhiAll_1d, CapInc_1d, LorenzFx, CapInc_1d_ord, Lorenz_x, Lorenz_y, gini_inc, sort_key)
        frac_below = LinInterp_1d(0.9d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Capital income share of top 10% of individuals = ', (1d0-frac_below)*100.0d0, '%'        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Capital income share of top 1% of individuals = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Capital income share of top 0.1% of individuals = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x,Lorenz_y,ntot)
        write(iunit_ineq,'(a60,f10.4,a3)') 'Capital income share of top 0.01% of individuals = ', (1d0-frac_below)*100.0d0, '%'   
        
        ! pe_share = sum( PE_inc_1d*PhiAll_1d, mask=(CapInc_1d >= cut_1pct) ) / sum( CapInc_1d*PhiAll_1d, mask=(CapInc_1d >= cut_1pct) )
        tmp_arr_sorted = PE_share_in_capinc_1d(sort_key)
        pe_share = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.99d0) )/0.01d0
        write(iunit_ineq,'(a60,f10.4)') 'PE/CapInc among top 1% of capital income = ', pe_share
        ! pe_share = sum( PE_inc_1d*PhiAll_1d, mask=(CapInc_1d >= cut_0_1pct) ) / sum( CapInc_1d*PhiAll_1d, mask=(CapInc_1d >= cut_0_1pct) )
        pe_share = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.999d0) )/0.001d0
        write(iunit_ineq,'(a60,f10.4)') 'PE/CapInc among top 0.1% of capital income = ', pe_share        
        ! pe_share = sum( PE_inc_1d*PhiAll_1d, mask=(CapInc_1d >= cut_0_01pct) ) / sum( CapInc_1d*PhiAll_1d, mask=(CapInc_1d >= cut_0_01pct) )
        pe_share = sum( tmp_arr_sorted*LorenzFx, mask=(Lorenz_x > 0.9999d0) )/0.0001d0
        write(iunit_ineq,'(a60,f10.4)') 'PE/CapInc among top 0.01% of capital income = ', pe_share         

        close(iunit_ineq)
        close(iunit_lc)

    end subroutine SummarizeDistribution
    
    function stdev(x, f)
    use mod_types, only: dp
        real(dp), intent(in) :: x(:), f(:)
        real(dp) :: stdev
        real(dp) :: mean
        integer :: i
        !mean = sum(x)/size(x)
        mean = sum(x*f)
        stdev = 0d0
        do i = 1, size(x)
            !stdev = stdev + (x(i)-mean)**2d0
            stdev = stdev + (x(i)-mean)**2d0*f(i)
        end do
        stdev = sqrt(stdev/size(x))
    end function stdev
    
    function covar(x, y, f)
    use mod_types, only: dp
        real(dp), intent(in) :: x(:), y(:), f(:)
        real(dp) :: covar
        real(dp) :: mean_x, mean_y
        integer :: i
        !mean_x = sum(x)/size(x)
        mean_x = sum(x*f)
        !mean_y = sum(y)/size(y)
        mean_y = sum(y*f)
        covar = 0d0
        do i = 1, size(x)
            !covar = covar + (x(i)-mean_x)*(y(i)-mean_y)
            covar = covar + (x(i)-mean_x)*(y(i)-mean_y)*f(i)
        end do
        covar = covar/size(x)
    end function covar
    
end module Mod_Distribution
