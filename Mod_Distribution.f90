module Mod_Distribution
    use PARAMS
    implicit none
    
    real(8), allocatable :: Phi(:, :, :, :, :, :, :)
    real(8), allocatable :: Phi_ret(:, :, :, :, :, :)
    !dimension(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork) :: Phi
    integer, parameter :: na_small = 50
    real(8) :: grida_small(na_small)
    
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
    end subroutine init_distr
    
    subroutine Distribution(save_res)
        ! THIS SUBROUTINE COMPUTES STEADY STATE DISTRIBUTION OF ASSETS

        use params
        use int_tictoc
        use svrgp_int
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
        real(8), dimension(Twork+Tret) :: abar, cbar !, inctaxbar !, yauxbar, afttaxauxbar !, inctaxbar_aboveyb
        real(8), dimension(Twork) :: lbar, labar
        integer :: iunit_lc
        CHARACTER (LEN=*), PARAMETER :: outDir = "tmp/"
        INTEGER :: iunit_phi        
        integer :: get_phi, get_phi_ret
        real(8) :: help(Twork+Tret)
        type(gpf):: gp
        integer :: id_tmp
        logical :: save_res
        real(8) :: labinc_tmp, totinc_tmp !, tax_tmp !, yaux_tmp, aftertaxaux_tmp, taxaboveyb_aux
        real(8) :: aprime !, TaxE
        real(8) :: TaxIncTest
        
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
        RAgg = 0d0
        RauxAgg = 0d0
        LAgg = 0d0
        TaxIncAboveYb = 0d0
        YfBelowYb = 0d0
        DBelowYb = 0d0
        TotOffshCost = 0d0
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
            do ia = 1, na
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        rtmp = rfunc(grida(ia), thetas(itheta), Kappas(ikappa))
                        do iz = 1, nz
                            do ixi = 1, nxi
                                do jj = 1, n_ofsh
                                    labinc_tmp = w*eta(iz)*xi(ixi)*ep(1,jc)*lfun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    totinc_tmp = labinc_tmp + rtmp*grida(ia)
                                    if (offshoring(ia,jj,itheta,ikappa,iz,ixi,jc) < 0.5d0) then
                                        YfBelowYb = YfBelowYb + min(yb_cutoff, totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        DBelowYb = DBelowYb + after_tax_income_aux(totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, totinc_tmp-yb_cutoff)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                    else
                                        YfBelowYb = YfBelowYb + min(yb_cutoff, frac_ofsh*totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        DBelowYb = DBelowYb + after_tax_income_aux(frac_ofsh*totinc_tmp)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                        TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, (frac_ofsh*totinc_tmp-yb_cutoff))*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)                                       
                                        TotOffshCost = TotOffshCost + psi_vals(jj)*Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc) 
                                    end if
                                    
                                    aprime = afun(ia,jj,itheta,ikappa,iz,ixi,jc)
                                    
                                    As = As + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*aprime*Nu(jc)
                                    LAgg = LAgg + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*eta(iz)*xi(ixi)*ep(1,jc)*lfun(ia,jj,itheta,ikappa,iz,ixi,jc)*Nu(jc)
                                    do ithetap = 1, ntheta
                                        do ikappap = 1, nkappa
                                            rtmp = rfunc(aprime, thetas(ithetap), Kappas(ikappap))
                                            RAgg = RAgg + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*rtmp*aprime*Nu(jc)    
                                            RauxAgg = RauxAgg + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*(rtmp-rbar)*aprime*Nu(jc)
                                            
                                            TrBn = TrBn + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*(aprime*(1d0 + rtmp))*Nu(jc)*(1d0-surv(jc))
                                            TaxE = TaxE + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*tau_estate*max(aprime*(1d0 + rtmp) - a_estate, 0d0)*Nu(jc)*(1d0-surv(jc))                                            
                                            
                                        end do
                                    end do

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
            do ia = 1, na
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        rtmp = rfunc(grida(ia), thetas(itheta), Kappas(ikappa))
                        do iz = 1, nz
                            totinc_tmp = b_ret(iz) + rtmp*grida(ia)
                            do jj = 1, n_ofsh
                                if (offshoring_ret(ia, jj, itheta, ikappa, iz) < 0.5d0) then
                                    YfBelowYb = YfBelowYb + min(yb_cutoff, totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    DBelowYb = DBelowYb + after_tax_income_aux(totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, totinc_tmp-yb_cutoff)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                else
                                    YfBelowYb = YfBelowYb + min(yb_cutoff, frac_ofsh*totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    DBelowYb = DBelowYb + after_tax_income_aux(frac_ofsh*totinc_tmp)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)
                                    TaxIncAboveYb = TaxIncAboveYb + tau_max*max(0d0, (frac_ofsh*totinc_tmp-yb_cutoff))*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc)                                    
                                    TotOffshCost = TotOffshCost + psi_vals(jj)*Phi_ret(ia,jj,itheta,ikappa,iz,jc)*Nu(Twork+jc) 
                                end if
                                                                     
                                RetAgg = RetAgg + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*b_ret(iz)*Nu(Twork+jc)
                                
                                aprime = afun_ret(ia,jj,itheta,ikappa,iz,jc)
                                As = As + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*aprime*Nu(Twork+jc)
                                do ithetap = 1, ntheta
                                    do ikappap = 1, nkappa
                                        rtmp = rfunc(aprime, thetas(ithetap), Kappas(ikappap))
                                        RAgg = RAgg + Phi_ret(ia,jj,ithetap,ikappap,iz,jc)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*rtmp*aprime*Nu(Twork+jc)
                                        RauxAgg = RauxAgg + Phi_ret(ia,jj,ithetap,ikappap,iz,jc)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*(rtmp-rbar)*aprime*Nu(Twork+jc)
                                        TrBn = TrBn + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*(aprime*(1d0 + rtmp))*Nu(Twork+jc)*(1d0-surv(Twork+jc))
                                        TaxE = TaxE + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*pi_kappa(ikappap)*pi_theta(itheta,ithetap)*tau_estate*max(aprime*(1d0 + rtmp) - a_estate, 0d0)*Nu(Twork+jc)*(1d0-surv(Twork+jc))
                                    end do
                                end do
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
        
        AAgg = sum( Nu(1:J)*abar(1:J) )
        print *, 'Test A: ', AAgg - As
        LabS = sum( Nu(1:Twork)*lbar(1:Twork) )
        print *, 'Test L: ', LAgg - LabS
        CAgg = sum( Nu(1:J)*cbar(1:J) )
        !Rs = sum( Nu(1:J)*rabar(1:J) )
        !print *, 'Test R: ', RAgg - Rs
        !Rs_aux = sum( Nu(1:J)*rabar_aux(1:J) )
        !TaxInc = sum( Nu(1:J)*inctaxbar(1:J) )
        !TaxIncTest = TaxIncAboveYb + YfBelowYb - theta0*DBelowYb
        TaxInc = TaxIncAboveYb + YfBelowYb - theta0*DBelowYb
        !print *, 'Test TaxInc: ', TaxInc - TaxIncTest
        !TaxS = sum( Nu(1:J)*(inctaxbar(1:J) + tauc*cbar(1:J)) ) + TaxE
        TaxC = tauc*CAgg
        TaxTot = TaxInc + TaxC + TaxE
        !YauxS = sum( Nu(1:J)*yauxbar(1:J) )
        !AftTaxauxS = sum( Nu(1:J)*afttaxauxbar(1:J) )
        !TaxaboveybS = sum( Nu(1:J)*inctaxbar_aboveyb(1:J) )
        !print *, 'Test Tax above Yb: ', TaxaboveybS - TaxIncAboveYb
        !print *, 'Test Yf below yb: ', YauxS - YfBelowYb
        !print *, 'Test NetY below yb: ', AftTaxauxS - DBelowYb
        
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

    end subroutine Distribution
    
end module Mod_Distribution
