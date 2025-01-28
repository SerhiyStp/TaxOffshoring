module Mod_Distribution
    use PARAMS
    implicit none
    
    real(8), allocatable :: Phi(:, :, :, :, :, :, :)
    !real(8), allocatable :: Tot_income(:, :, :, :, :, :, :)
    real(8), allocatable :: Phi_ret(:, :, :, :, :, :)
    !dimension(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork) :: Phi
    integer, parameter :: na_small = 50
    real(8) :: grida_small(na_small)
    integer, parameter :: ntot_work = na_small*n_ofsh*ntheta*nkappa*nz*nxi*Twork
    integer, parameter :: ntot = ntot_work + na*n_ofsh*ntheta*nkappa*nz*Tret
    real(8), allocatable :: PhiAll_1d(:), TotInc_1d(:), Wealth_1d(:), Lorenz_x(:), Lorenz_y(:), &
        TotInc_1d_ord(:), Wealth_1d_ord(:), PE_share_1d(:), OffshWealth_1d(:), Wages_1d(:), Hours_1d(:), &
        Wages_1d_ord(:), Hours_1d_ord(:), LabInc_1d(:), LabInc_1d_ord(:)
    
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
        allocate(PE_share_1d(ntot), stat=ierr)
        allocate(OffshWealth_1d(ntot), stat=ierr)
        allocate(Wages_1d(ntot), stat=ierr)
        allocate(Hours_1d(ntot), stat=ierr)
        allocate(Wages_1d_ord(ntot), stat=ierr)
        allocate(Hours_1d_ord(ntot), stat=ierr)
        allocate(LabInc_1d(ntot), stat=ierr)  
        allocate(LabInc_1d_ord(ntot), stat=ierr)
        
    end subroutine init_distr
    
    subroutine Distribution(save_res)
        ! THIS SUBROUTINE COMPUTES STEADY STATE DISTRIBUTION OF ASSETS

        use params
        use int_tictoc
        !use svrgp_int
        use toolbox
        use MyLinInterp
        use Mod_Household, only: afun, lfun, cfun, afun_ret, cfun_ret, offshoring, offshoring_ret
        use io, only: save_array, read_array
        use ogpf
        !use moments, only: sim_moms, sim_moms_aux, sim_moms_klp

        IMPLICIT NONE
        
        integer :: ia, itheta, ikappa, iz,  ixi, jj
        integer :: ithetap, ikappap, izp, ixip
        real(8) :: test, test2, test3
        integer :: inds(2)
        real(8) :: vals(2)
        real(8) :: TT1, TT2
        real(8) :: rtmp
        real(8), dimension(Twork+Tret) :: abar, cbar, inctaxbar !, yauxbar, afttaxauxbar !, inctaxbar_aboveyb
        real(8), dimension(Twork) :: lbar, labar
        integer :: iunit_lc
        CHARACTER (LEN=*), PARAMETER :: outDir = "tmp/"
        INTEGER :: iunit_phi        
        integer :: get_phi, get_phi_ret
        real(8) :: help(Twork+Tret)
        type(gpf):: gp
        integer :: id_tmp
        logical :: save_res
        real(8) :: labinc_tmp, totinc_tmp, tax_tmp !, yaux_tmp, aftertaxaux_tmp, taxaboveyb_aux
        real(8) :: aprime !, TaxE
        real(8) :: TaxIncTest
        real(8) :: RaggTmp, RauxaggTmp, TrBnTmp, TaxETmp
        real(8) :: atmp, ctmp, aprimetmp, testbc_tmp, aftertaxtmp, aftertaxtest
        real(8) :: rptmp
        real(8) :: AsOffshore
        integer :: ii
        real(8) :: gini_inc, gini_wealth, gini_test
        real(8) :: LorenzFx(ntot)
        real(8) :: frac_below
        real(8) :: PE_share, r_PE, pe_share_tmp
        real(8) :: inc_cut_1pct, inc_cut_0_1pct, inc_cut_0_01pct, inc_cut_10pct
        real(8) :: wealth_cut_1pct, wealth_cut_0_1pct, wealth_cut_0_01pct
        real(8) :: pe_share_1pct, pe_share_bot90pct, pe_share_0_1pct, pe_share_0_01pct
        real(8) :: offsh_share_1pct, offsh_share_0_1pct, offsh_share_0_01pct
        real(8) :: d_test
        
        ! Initialize Distribution by Computing Distribution for first Generation
        get_phi = 1
        get_phi_ret = 1
        
        
        if (get_phi == 1) then
            Phi=0.0d0
            !$OMP PARALLEL PRIVATE(jj, itheta, ikappa, iz, ixi)
            !$OMP DO SCHEDULE(DYNAMIC)
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
                
            ! test = sum(Phi(:,:,:,:,:,:,1))
            
            ! Loop to find distributions for ages 2 to J
            !!call tic
            
            do jc=2,Twork
                !$OMP PARALLEL PRIVATE(jj, itheta, ikappa, iz, ixi, ia, vals, inds, TT1, TT2) SHARED(jc)
                !$OMP DO SCHEDULE(DYNAMIC)
                do jj=1,n_ofsh
                    do itheta=1,ntheta
                        do ikappa=1,nkappa
                            do iz=1,nz
                                do ixi=1,nxi
                                    do ia=1,na
                                        call basefun(grida(1:na),na,afun(ia,jj,itheta,ikappa,iz,ixi,jc-1),vals,inds)
                                        do izp=1,nz
                                            do ithetap=1,ntheta
                                                do ikappap=1,nkappa
                                                    do ixip=1,nxi
                                                        TT1 = vals(1)*pi(iz,izp)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)*pi_xi(ixip)
                                                        TT2 = vals(2)*pi(iz,izp)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)*pi_xi(ixip)
                                                        Phi(inds(1),jj,ithetap,ikappap,izp,ixip,jc)=Phi(inds(1),jj,ithetap,ikappap,izp,ixip,jc)+Phi(ia,jj,itheta,ikappa,iz,ixi,jc-1)*TT1
                                                        Phi(inds(2),jj,ithetap,ikappap,izp,ixip,jc)=Phi(inds(2),jj,ithetap,ikappap,izp,ixip,jc)+Phi(ia,jj,itheta,ikappa,iz,ixi,jc-1)*TT2
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
                ! test = sum(Phi(:,:,:,:,:,:,jc))
                !print *, jc, test 
            end do
            
            ! call save_array(Phi, outDir // "Phi.bin")
        
        else 
            call read_array(Phi, outDir // "Phi.bin")
            test = sum(Phi(:,:,:,:,:,:,Twork))
            print *, jc, test            
        end if
        
        
        if (get_phi_ret == 1) then
            Phi_ret = 0d0
            !print *, 'Retirement: '
            ! First period of retirement
            !$OMP PARALLEL PRIVATE(jj,itheta,ikappa,iz,ixi,ia,ithetap,ikappap,vals,inds,TT1,TT2)
            !$OMP DO SCHEDULE(DYNAMIC)
            do jj = 1, n_ofsh
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        do iz = 1, nz
                            do ixi = 1, nxi
                                do ia = 1,na
                                    call basefun(grida(1:na),na,afun(ia,jj,itheta,ikappa,iz,ixi,Twork),vals,inds)
                                    do ithetap=1,ntheta
                                        do ikappap=1,nkappa
                                            TT1 = vals(1)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                            TT2 = vals(2)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                            Phi_ret(inds(1),jj,ithetap,ikappap,iz,1) = Phi_ret(inds(1),jj,ithetap,ikappap,iz,1) +  Phi(ia,jj,itheta,ikappa,iz,ixi,Twork)*TT1 
                                            Phi_ret(inds(2),jj,ithetap,ikappap,iz,1) = Phi_ret(inds(2),jj,ithetap,ikappap,iz,1) +  Phi(ia,jj,itheta,ikappa,iz,ixi,Twork)*TT2
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
            ! test = sum(Phi_ret(:,:,:,:,:,1))
            ! if (abs(test-1d0) > 1d-9) then
            !     print *, Twork+1, test
            ! end if
            !print *, Twork+1, test 
            
            ! All other retirement periods
            do jc = 2, Tret
                !$OMP PARALLEL PRIVATE(jj,itheta,ikappa,iz,ixi,ia,vals,inds,TT1,TT2)
                !$OMP DO SCHEDULE(DYNAMIC)
                do jj = 1, n_ofsh
                    do itheta = 1, ntheta
                        do ikappa = 1, nkappa
                            do iz = 1, nz
                                do ia = 1,na
                                    call basefun(grida(1:na),na,afun_ret(ia,jj,itheta,ikappa,iz,jc-1),vals,inds)
                                    do ithetap=1,ntheta
                                        do ikappap=1,nkappa
                                            TT1 = vals(1)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                            TT2 = vals(2)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                            Phi_ret(inds(1),jj,ithetap,ikappap,iz,jc) = Phi_ret(inds(1),jj,ithetap,ikappap,iz,jc) +  Phi_ret(ia,jj,itheta,ikappa,iz,jc-1)*TT1 
                                            Phi_ret(inds(2),jj,ithetap,ikappap,iz,jc) = Phi_ret(inds(2),jj,ithetap,ikappap,iz,jc) +  Phi_ret(ia,jj,itheta,ikappa,iz,jc-1)*TT2
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do 
                !$OMP END DO
                !$OMP END PARALLEL                
                ! test = sum(Phi_ret(:,:,:,:,:,jc))
                ! if (abs(test-1d0) > 1d-9) then
                !     print *, Twork+jc, test
                ! end if
                !print *, Twork+jc, test            
            end do
            ! call save_array(Phi_ret, outDir // "Phi_ret.dat")
        else
            call read_array(Phi_ret, outDir // "Phi_ret.dat")    
        end if
        
        !!call tic
        
        ! Find Stationary Distribution over Asset Holdings
        
        open(newunit=iunit_lc,file='LifeCycle.txt')
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
                                    TotInc_1d(ii) = totinc_tmp
                                    LabInc_1d(ii) = labinc_tmp
                                    Wages_1d(ii) = w*eta(iz)*xi(ixi)*ep(1,jc)
                                    Hours_1d(ii) = lfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    ctmp = cfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    aprimetmp = afun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    Wealth_1d(ii) = grida(ia) 
                                    !Wealth_1d(ii) = aprimetmp
                                    ii = ii + 1
                                    
                                    
                                    if (offshoring(ia,jj,itheta,ikappa,iz,ixi,jc) < 0.5d0) then
                                        tax_tmp = tax_income(totinc_tmp)
                                        YfBelowYb = YfBelowYb + min(yb_cutoff, totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        DBelowYb = DBelowYb + after_tax_income_aux(totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, totinc_tmp-yb_cutoff)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                    
                                        aftertaxtmp = totinc_tmp - tax_tmp  
                                        OffshWealth_1d(ii) = 0d0
                                    else
                                        tax_tmp = tax_income((1d0-frac_ofsh)*totinc_tmp)
                                        YfBelowYb = YfBelowYb + min(yb_cutoff, (1d0-frac_ofsh)*totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        DBelowYb = DBelowYb + after_tax_income_aux((1d0-frac_ofsh)*totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, ((1d0-frac_ofsh)*totinc_tmp-yb_cutoff))*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)                                       
                                        TotOffshCost = TotOffshCost + psi_vals(jj)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc) 
                                        
                                        aftertaxtmp = (1d0-frac_ofsh)*totinc_tmp - tax_tmp + frac_ofsh*totinc_tmp - psi_vals(jj)
                                        AsOffshore = AsOffshore + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*frac_ofsh*grida(ia)*Nu(jc)
                                        !OffshWealth_1d(ii) = frac_ofsh*grida(ia)
                                        OffshWealth_1d(ii) = frac_ofsh
                                    end if
                                    
                                    
                                    inctaxbar(jc) = inctaxbar(jc) + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*tax_tmp
                                    
                                    aprime = afun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    
                                    As = As + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*aprime*Nu(jc)
                                    PE_share = PE_share + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*pe_share_tmp*Nu(jc)
                                    PE_share_1d(ii) = pe_share_tmp
                                    LAgg = LAgg + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*eta(iz)*xi(ixi)*ep(1,jc)*lfun(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                    
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
                                LabInc_1d(ii) = 0d0
                                Wages_1d(ii) = 0d0
                                Hours_1d(ii) = 0d0
                                                               
                                
                                if (offshoring_ret(ia, jj, itheta, ikappa, iz) < 0.5d0) then
                                    tax_tmp = tax_income(totinc_tmp)
                                    YfBelowYb = YfBelowYb + min(yb_cutoff, totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    DBelowYb = DBelowYb + after_tax_income_aux(totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, totinc_tmp-yb_cutoff)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    OffshWealth_1d(ii) = 0d0
                                else
                                    tax_tmp = tax_income((1d0-frac_ofsh)*totinc_tmp)
                                    YfBelowYb = YfBelowYb + min(yb_cutoff, (1d0-frac_ofsh)*totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    DBelowYb = DBelowYb + after_tax_income_aux((1d0-frac_ofsh)*totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, ((1d0-frac_ofsh)*totinc_tmp-yb_cutoff))*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)                                    
                                    TotOffshCost = TotOffshCost + psi_vals(jj)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    
                                    AsOffshore = AsOffshore + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*frac_ofsh*grida(ia)*Nu(Twork+jc)
                                    !OffshWealth_1d(ii) = frac_ofsh*grida(ia)
                                    OffshWealth_1d(ii) = frac_ofsh
                                end if
                                                                     
                                RetAgg = RetAgg + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*b_ret(iz)*Nu(Twork+jc)
                                
                                inctaxbar(Twork+jc) = inctaxbar(Twork+jc) + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*tax_tmp
                                
                                aprime = afun_ret(ia,jj,itheta,ikappa,iz,jc)
                                Wealth_1d(ii) = grida(ia) 
                                !Wealth_1d(ii) = aprimetmp
                                PE_share_1d(ii) = pe_share_tmp
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
        
        
        ! Earnings and wealth distribution
        test = sum(PhiAll_1d) 
        print *, test
        PhiAll_1d = PhiAll_1d/test
        
        call lorenz_s(PhiAll_1d, TotInc_1d, LorenzFx, gini_inc, TotInc_1d_ord, Lorenz_x, Lorenz_y)
        print '(a30,f10.4)', 'Income Gini = ', gini_inc
        !call gp%title('Lorenz curve (total income)')
        !call gp%xlabel('Population share')
        !call gp%ylabel('Income share')
        !call gp%plot(Lorenz_x, Lorenz_y, 'with lines') 
        
        inc_cut_10pct = LinInterp_1d(0.9d0,Lorenz_x,TotInc_1d_ord,ntot)
        inc_cut_1pct = LinInterp_1d(0.99d0,Lorenz_x,TotInc_1d_ord,ntot)
        inc_cut_0_1pct = LinInterp_1d(0.999d0,Lorenz_x,TotInc_1d_ord,ntot)
        inc_cut_0_01pct = LinInterp_1d(0.9999d0,Lorenz_x,TotInc_1d_ord,ntot)
        
        test = sum( LorenzFx, mask=(TotInc_1d_ord > inc_cut_1pct) )
        d_test = abs(test - 0.01d0)/0.01d0
        if (d_test > 1d-2) then
            print '(a,f8.5,a,f8.5)', 'Error in the mass of income share of top 1%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if        
        pe_share_1pct = sum( PE_share_1d*LorenzFx, mask=(TotInc_1d_ord > inc_cut_1pct) ) / test
        print '(a30,f10.4,a3)', 'PE share among top 1% = ', pe_share_1pct*100.0d0, '%'
        
        test = sum( LorenzFx, mask=(TotInc_1d_ord > inc_cut_0_1pct) )
        d_test = abs(test - 0.001d0)/0.001d0
        if (d_test > 1d-2) then
            print '(a,f8.5,a,f8.5)', 'Error in the mass of income share of top 0.1%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if           
        pe_share_0_1pct = sum( PE_share_1d*LorenzFx, mask=(TotInc_1d_ord > inc_cut_0_1pct) ) / test        
        print '(a30,f10.4,a3)', 'PE share among top 0.1% = ', pe_share_0_1pct*100.0d0, '%'
        
        test = sum( LorenzFx, mask=(TotInc_1d_ord > inc_cut_0_01pct) )
        d_test = abs(test - 0.0001d0)/0.0001d0
        if (d_test > 1d-2) then
            print '(a,f8.5,a,f8.5)', 'Error in the mass of income share of top 0.01%, test = ', test*100.0d0, ', d_test = ', d_test*100.0d0
        end if           
        pe_share_0_01pct = sum( PE_share_1d*LorenzFx, mask=(TotInc_1d_ord > inc_cut_0_01pct) ) / test     
        print '(a30,f10.4,a3)', 'PE share among top 0.01% = ', pe_share_0_01pct*100.0d0, '%'
        
        test = sum( LorenzFx, mask=(TotInc_1d_ord <= inc_cut_10pct) )
        pe_share_bot90pct = sum( PE_share_1d*LorenzFx, mask=(TotInc_1d_ord <= inc_cut_10pct) ) / test
        print '(a30,f10.4,a3)', 'PE share among bottom 90% = ', pe_share_bot90pct*100.0d0, '%'
        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x,Lorenz_y,ntot)
        print '(a30,f10.4,a3)', 'Income share of top 1% = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x,Lorenz_y,ntot)
        print '(a30,f10.4,a3)', 'Income share of top 0.1% = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x,Lorenz_y,ntot)
        print '(a30,f10.4,a3)', 'Income share of top 0.01% = ', (1d0-frac_below)*100.0d0, '%'          
                
        call lorenz_s(PhiAll_1d, Wealth_1d, LorenzFx, gini_wealth, Wealth_1d_ord, Lorenz_x, Lorenz_y)
        print '(a30,f10.4)', 'Wealth Gini = ', gini_wealth
        !call gp%title('Lorenz curve (total income)')
        !call gp%xlabel('Population share')
        !call gp%ylabel('Income share')
        !call gp%plot(Lorenz_x, Lorenz_y, 'with lines')    
        
        wealth_cut_1pct = LinInterp_1d(0.99d0,Lorenz_x,Wealth_1d_ord,ntot)
        wealth_cut_0_1pct = LinInterp_1d(0.999d0,Lorenz_x,Wealth_1d_ord,ntot)
        wealth_cut_0_01pct = LinInterp_1d(0.9999d0,Lorenz_x,Wealth_1d_ord,ntot)
        
        test = sum( LorenzFx, mask=(Wealth_1d_ord > wealth_cut_1pct) )
        d_test = abs(test - 0.01d0)/0.01d0
        if (d_test > 1d-2) then
            print '(a,f8.5,a,f8.5)', 'Error in the mass of wealth share of top 1%, test = ', test*100.0d0, ', d_test = ', d_test
        end if
        offsh_share_1pct = sum( OffshWealth_1d*LorenzFx, mask=(Wealth_1d_ord > wealth_cut_1pct) ) / test
        print '(a45,f10.4,a3)', 'Share of wealth offshored among top 1% = ', offsh_share_1pct*100.0d0, '%'
        test = sum( LorenzFx, mask=(Wealth_1d_ord > wealth_cut_0_1pct) )
        d_test = abs(test - 0.001d0)/0.001d0
        if (d_test > 1d-2) then
            print '(a,f8.5,a,f8.5)', 'Error in the mass of wealth share of top 0.1%, test = ', test*100.0d0, ', d_test = ', d_test
        end if
        offsh_share_0_1pct = sum( OffshWealth_1d*LorenzFx, mask=(Wealth_1d_ord > wealth_cut_0_1pct) ) / test   
        print '(a45,f10.4,a3)', 'Share of wealth offshored among top 0.1% = ', offsh_share_0_1pct*100.0d0, '%'
        test = sum( LorenzFx, mask=(Wealth_1d_ord > wealth_cut_0_01pct) )
        d_test = abs(test - 0.0001d0)/0.0001d0
        if (d_test > 1d-2) then
            print '(a,f8.5,a,f8.5)', 'Error in the mass of wealth share of top 0.01%, test = ', test*100.0d0, ', d_test = ', d_test
        end if
        offsh_share_0_01pct = sum( OffshWealth_1d*LorenzFx, mask=(Wealth_1d_ord > wealth_cut_0_01pct) ) / test   
        print '(a45,f10.4,a3)', 'Share of wealth offshored among top 0.01% = ', offsh_share_0_01pct*100.0d0, '%'
        
        frac_below = LinInterp_1d(0.99d0,Lorenz_x,Lorenz_y,ntot)
        print '(a30,f10.4,a3)', 'Wealth share of top 1% = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.999d0,Lorenz_x,Lorenz_y,ntot)
        print '(a30,f10.4,a3)', 'Wealth share of top 0.1% = ', (1d0-frac_below)*100.0d0, '%'
        frac_below = LinInterp_1d(0.9999d0,Lorenz_x,Lorenz_y,ntot)
        print '(a30,f10.4,a3)', 'Wealth share of top 0.01% = ', (1d0-frac_below)*100.0d0, '%'    
        
        
        
        
    end subroutine Distribution
    
end module Mod_Distribution
