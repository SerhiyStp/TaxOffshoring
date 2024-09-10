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
        real(8), dimension(Twork+Tret) :: abar, rabar, cbar, totincbar, inctaxbar, rabar_aux, yauxbar, afttaxauxbar, inctaxbar_aboveyb
        real(8), dimension(Twork) :: lbar, labar
        integer :: iunit_lc
        CHARACTER (LEN=*), PARAMETER :: outDir = "tmp/"
        INTEGER :: iunit_phi        
        integer :: get_phi, get_phi_ret
        real(8) :: help(Twork+Tret)
        type(gpf):: gp
        integer :: id_tmp
        logical :: save_res
        real(8) :: labinc_tmp, totinc_tmp, tax_tmp, yaux_tmp, aftertaxaux_tmp, taxaboveyb_aux
        
        ! Initialize Distribution by Computing Distribution for first Generation
        get_phi = 1
        get_phi_ret = 1
        
        
        if (get_phi == 1) then
            
        Phi=0.0d0
        !(na, n_ofsh, ntheta, nkappa, nz, nxi, Twork)
        !do ia = 1, na
            do jj = 1, n_ofsh
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        do iz = 1, nz
                            do ixi = 1, nxi
                                !Phi(ia,jj,itheta,ikappa,iz,ixi,1)
                                Phi(1,jj,itheta,ikappa,iz,ixi,1) = pi_kappa(ikappa)*pi_xi(ixi)*pini(iz)*psi_prob(jj)*pi_theta_stat(itheta)   
                            end do
                        end do
                    end do
                end do
            end do
        !end do
            
        test = sum(Phi(:,:,:,:,:,:,1))
        
        ! Loop to find distributions for ages 2 to J
        !!call tic
        do jc=2,Twork
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
            test = sum(Phi(:,:,:,:,:,:,jc))
            print *, jc, test
        end do
        
        call save_array(Phi, outDir // "Phi.bin")
        
        else 
            call read_array(Phi, outDir // "Phi.bin")
            test = sum(Phi(:,:,:,:,:,:,Twork))
            print *, jc, test            
        end if
        
        
        if (get_phi_ret == 1) then
        Phi_ret = 0d0
        print *, 'Retirement: '
        ! First period of retirement
        do jj = 1, n_ofsh
            do itheta = 1, ntheta
                do ikappa = 1, nkappa
                    do iz = 1, nz
                        do ixi = 1, nxi
                            do ia = 1,na
                                call basefun(grida(1:na),na,afun(ia,jj,itheta,ikappa,iz,ixi,Twork),vals,inds)
                                test = 0d0
                                do ithetap=1,ntheta
                                    do ikappap=1,nkappa
                                        TT1 = vals(1)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                        TT2 = vals(2)*pi_theta(itheta,ithetap)*pi_kappa(ikappap)
                                        !test = test + TT1 + TT2
                                        Phi_ret(inds(1),jj,ithetap,ikappap,iz,1) = Phi_ret(inds(1),jj,ithetap,ikappap,iz,1) +  Phi(ia,jj,itheta,ikappa,iz,ixi,Twork)*TT1 
                                        Phi_ret(inds(2),jj,ithetap,ikappap,iz,1) = Phi_ret(inds(2),jj,ithetap,ikappap,iz,1) +  Phi(ia,jj,itheta,ikappa,iz,ixi,Twork)*TT2
                                    end do
                                end do
                                !if (abs(test-1d0) > 1d-9) then
                                !    print *, test
                                !end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
        test = sum(Phi_ret(:,:,:,:,:,1))
        print *, Twork+1, test 
        
        ! All other retirement periods
        do jc = 2, Tret
            do jj = 1, n_ofsh
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        do iz = 1, nz
                            do ia = 1,na
                                !(na, n_ofsh, ntheta, nkappa, nz, Tret)
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
            test = sum(Phi_ret(:,:,:,:,:,jc))
            print *, Twork+jc, test            
        end do
        
        call save_array(Phi_ret, outDir // "Phi_ret.dat")
        
        else
            call read_array(Phi_ret, outDir // "Phi_ret.dat")    
        end if
        
        !!call tic
        
        ! Find Stationary Distribution over Asset Holdings
        
        open(newunit=iunit_lc,file='LifeCycle.txt')
        RetS = 0d0
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
            rabar(jc) = 0d0
            rabar_aux(jc) = 0d0
            do ia = 1, na
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        rtmp = rfunc(grida(ia), thetas(itheta), Kappas(ikappa))
                        rabar(jc) = rabar(jc) + sum(Phi(ia,1:n_ofsh,itheta,ikappa,1:nz,1:nxi,jc))*rtmp*grida(ia)
                        rabar_aux(jc) = rabar_aux(jc) + sum(Phi(ia,1:n_ofsh,itheta,ikappa,1:nz,1:nxi,jc))*(rtmp-rbar)*grida(ia)
                    end do
                end do
            end do
            totincbar(jc) = 0d0
            inctaxbar(jc) = 0d0
            yauxbar(jc) = 0d0
            afttaxauxbar(jc) = 0d0
            inctaxbar_aboveyb(jc) = 0d0
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
                                        tax_tmp = tax_income(totinc_tmp)
                                        yaux_tmp = min(yb_cutoff, totinc_tmp)
                                        aftertaxaux_tmp = after_tax_income_aux(totinc_tmp)
                                        taxaboveyb_aux = tau_max*max(0d0, totinc_tmp-yb_cutoff)
                                    else
                                        tax_tmp = tax_income(frac_ofsh*totinc_tmp)
                                        yaux_tmp = min(yb_cutoff, frac_ofsh*totinc_tmp)
                                        aftertaxaux_tmp = after_tax_income_aux(frac_ofsh*totinc_tmp)
                                        test2 = (1d0-tau_max)*max(0d0, (frac_ofsh*totinc_tmp-yb_cutoff))
                                        taxaboveyb_aux = tau_max*max(0d0, (frac_ofsh*totinc_tmp-yb_cutoff))
                                        !test = taxaboveyb_aux + yaux_tmp - aftertaxaux_tmp
                                        !if (abs(test - tax_tmp) > 1d-8) then
                                        !    print *, 'WARNING'
                                        !    tax_tmp = tax_income(frac_ofsh*totinc_tmp)
                                        !    yaux_tmp = min(yb_cutoff, frac_ofsh*totinc_tmp)
                                        !    aftertaxaux_tmp = after_tax_income_aux(frac_ofsh*totinc_tmp)
                                        !    test2 = (1d0-tau_max)*max(0d0, (frac_ofsh*totinc_tmp-yb_cutoff))
                                        !    taxaboveyb_aux = tau_max*max(0d0, (frac_ofsh*totinc_tmp-yb_cutoff))
                                        !    test3 = taxaboveyb_aux + yaux_tmp - aftertaxaux_tmp
                                        !end if                                         
                                    end if
                                    
                                    !test = taxaboveyb_aux + yaux_tmp - aftertaxaux_tmp
                                    !if (abs(test - tax_tmp) > 1d-8) then
                                    !    print *, 'WARNING'
                                    !end if
                                    
                                    totincbar(jc) = totincbar(jc) + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*totinc_tmp
                                    inctaxbar(jc) = inctaxbar(jc) + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*tax_tmp
                                    yauxbar(jc) = yauxbar(jc) + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*yaux_tmp
                                    afttaxauxbar(jc) = afttaxauxbar(jc) + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*aftertaxaux_tmp
                                    inctaxbar_aboveyb(jc) = inctaxbar_aboveyb(jc) + Phi(ia,jj,itheta,ikappa,iz,ixi,jc)*taxaboveyb_aux
                                end do    
                            end do
                        end do
                    end do
                end do
            end do
            write(iunit_lc, '(i0, 7f9.6)') jc, abar(jc), labar(jc), lbar(jc), rabar(jc), cbar(jc), totincbar(jc), inctaxbar(jc)
        end do
        do jc = 1, Tret
            abar(Twork+jc) = sum(Phi_ret(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc)*afun_ret(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc))
            cbar(Twork+jc) = sum(Phi_ret(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc)*cfun_ret(1:na,1:n_ofsh,1:ntheta,1:nkappa,1:nz,jc))
            rabar(Twork+jc) = 0d0
            totincbar(Twork+jc) = 0d0
            inctaxbar(Twork+jc) = 0d0    
            rabar_aux(Twork+jc) = 0d0
            yauxbar(Twork+jc) = 0d0
            afttaxauxbar(Twork+jc) = 0d0
            inctaxbar_aboveyb(Twork+jc) = 0d0
            do ia = 1, na
                do itheta = 1, ntheta
                    do ikappa = 1, nkappa
                        rtmp = rfunc(grida(ia), thetas(itheta), Kappas(ikappa))
                        rabar(Twork+jc) = rabar(Twork+jc) + sum(Phi_ret(ia,1:n_ofsh,itheta,ikappa,1:nz,jc))*rtmp*grida(ia)
                        rabar_aux(Twork+jc) = rabar_aux(Twork+jc) + sum(Phi_ret(ia,1:n_ofsh,itheta,ikappa,1:nz,jc))*(rtmp-rbar)*grida(ia)
                        do iz = 1, nz
                            totinc_tmp = b_ret(iz) + rtmp*grida(ia)
                            do jj = 1, n_ofsh
                                if (offshoring_ret(ia, jj, itheta, ikappa, iz) < 0.5d0) then
                                    tax_tmp = tax_income(totinc_tmp)
                                    yaux_tmp = min(yb_cutoff, totinc_tmp)
                                    aftertaxaux_tmp = after_tax_income_aux(totinc_tmp)
                                    taxaboveyb_aux = tau_max*max(0d0, totinc_tmp-yb_cutoff)
                                else
                                    tax_tmp = tax_income(frac_ofsh*totinc_tmp)
                                    yaux_tmp = min(yb_cutoff, frac_ofsh*totinc_tmp)
                                    aftertaxaux_tmp = after_tax_income_aux(frac_ofsh*totinc_tmp)
                                    taxaboveyb_aux = tau_max*max(0d0, (frac_ofsh*totinc_tmp-yb_cutoff))
                                    !test = taxaboveyb_aux + yaux_tmp - aftertaxaux_tmp
                                    !if (abs(test - tax_tmp) > 1d-8) then
                                    !    print *, 'WARNING'
                                    !end if                                      
                                end if
                                
                                !test = taxaboveyb_aux + yaux_tmp - aftertaxaux_tmp
                                !if (abs(test - tax_tmp) > 1d-8) then
                                !    print *, 'WARNING'
                                !end if                                
                                
                                totincbar(Twork+jc) = totincbar(Twork+jc) + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*totinc_tmp
                                inctaxbar(Twork+jc) = inctaxbar(Twork+jc) + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*tax_tmp
                                yauxbar(Twork+jc) = yauxbar(Twork+jc) + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*yaux_tmp
                                afttaxauxbar(Twork+jc) = afttaxauxbar(Twork+jc) + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*aftertaxaux_tmp
                                inctaxbar_aboveyb(Twork+jc) = inctaxbar_aboveyb(Twork+jc) + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*taxaboveyb_aux
                                RetS = RetS + Phi_ret(ia,jj,itheta,ikappa,iz,jc)*b_ret(iz)*Nu(Twork+jc)
                            end do
                        end do
                    end do
                end do
            end do
            write(iunit_lc, '(i0, 7f9.6)') Twork+jc, abar(Twork+jc), 0d0, 0d0, rabar(Twork+jc), cbar(Twork+jc), totincbar(Twork+jc), inctaxbar(Twork+jc)
        end do
        close(iunit_lc)
        
        As = sum( Nu(1:J)*abar(1:J) )
        LabS = sum( Nu(1:Twork)*lbar(1:Twork) )
        Rs = sum( Nu(1:J)*rabar(1:J) )
        Rs_aux = sum( Nu(1:J)*rabar_aux(1:J) )
        TaxS = sum( Nu(1:J)*(inctaxbar(1:J) + tauc*cbar(1:J)) )
        TaxCS = sum( Nu(1:J)*(tauc*cbar(1:J)) )
        YauxS = sum( Nu(1:J)*yauxbar(1:J) )
        AftTaxauxS = sum( Nu(1:J)*afttaxauxbar(1:J) )
        TaxaboveybS = sum( Nu(1:J)*inctaxbar_aboveyb(1:J) )
        
        !test = (YauxS - AftTaxauxS)
        !test2 = TaxaboveybS + test + TaxCS
        !test3 = TaxS - test2
        
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
            pause
        else
            ! Annotation: set title, xlabel, ylabel, line specification
            call gp%title('Distribution over assets')
            call gp%xlabel('assets')
            call gp%ylabel('f')

            !Call Plot to draw a vector against a vector of data
            !The last argument defines the line specification
            !call gp%plot(grida,ADis,'with linespoints lt 2 pt 4')    
            id_tmp = maxloc(grida, dim=1, mask=(ADis > 1d-6))
            call gp%plot(grida(1:id_tmp),ADis(1:id_tmp),'with lines') 
        end if
        !
        !do jc=1,J
        !    test = sum(Phi(1:nty,1:ns,1:na,jc,1:n_ofsh))
        !    if ( (sum(Phi(1:nty,1:ns,1:na,jc,1:n_ofsh)) > 1.001) .or. (sum(Phi(1:nty,1:ns,1:na,jc,1:n_ofsh))<0.999 )) then
        !        print*,'Should equal 1', sum(Phi(1:nty,1:ns,1:na,jc,1:n_ofsh))
        !        print *, 'Phi in DISTRIB.f90.'
        !        !pause
        !    end if
        !end do
        !
        !! Find Aggregate Asset Holdings (end of Period)
        !
        !do jc=1,J
        !    abar(jc)=0.0
        !    astartbar(jc)=0.0
        !    do tyc=1,nty
        !        do jj=1,n_ofsh
        !            abartype(jc,tyc,jj)=0.0
        !            astartbartype(jc,tyc,jj)=0.0
        !            do sc=1,ns
        !                !abar(jc)=abar(jc)+ sum( Phi(tyc,sc,1:na,jc,jj)*afun(tyc,sc,1:na,jc,jj), mask=afun(tyc,sc,1:na,jc,jj) < a_cutoff(jj) ) 
        !                abar(jc)=abar(jc)+ sum( Phi(tyc,sc,1:na,jc,jj)*afun(tyc,sc,1:na,jc,jj) )
        !                astartbar(jc)=astartbar(jc)+ sum( Phi(tyc,sc,1:na,jc,jj)*grida(1:na) )
        !                abartype(jc,tyc,jj)=abartype(jc,tyc,jj)+ sum( Phi(tyc,sc,1:na,jc,jj)*afun(tyc,sc,1:na,jc,jj) )/measty(tyc)/psi_prob(jj) 			
        !                astartbartype(jc,tyc,jj)=astartbartype(jc,tyc,jj)+ sum( Phi(tyc,sc,1:na,jc,jj)*grida(1:na) )/measty(tyc)/psi_prob(jj)
        !                
        !                if (Phi(tyc,sc,na,jc,jj) > 1d-9) then 
        !                    print *, Phi(tyc,sc,na,jc,jj)
        !                    print *, 'WARNING: positive probability of hitting Abar at age ', jc, ', productivity state ', sc
        !                end if  
        !                
        !            end do
        !        end do
        !    end do
        !end do
        !
        !
        !
        !
        !! Computing aggregate labor supply
        !
        !do jc=1,J
        !    lbar(jc)=0.0
        !    do tyc=1,nty
        !        do jj=1,n_ofsh
        !            lbartype(jc,tyc,jj)=0.0
        !            do sc=1,ns
        !                lbar(jc)=lbar(jc)+ sum(eta(sc)*ep(tyc,jc)*Phi(tyc,sc,1:na,jc,jj)*lfun(tyc,sc,1:na,jc,jj))
        !                lbartype(jc,tyc,jj)=lbartype(jc,tyc,jj) + sum(eta(sc)*ep(tyc,jc)*Phi(tyc,sc,1:na,jc,jj)*lfun(tyc,sc,1:na,jc,jj))/measty(tyc)/psi_prob(jj)
        !            end do
        !        end do
        !    end do
        !end do
        !
        !!open(15, file='long_hours.txt')
        !do jc=1,J
        !    labar(jc)=0.0
        !    lbar_test(jc,1) = sum(Phi(1:nty,1:ns,1:na,jc,1:n_ofsh), mask=(lfun(1:nty,1:ns,1:na,jc,1:n_ofsh)>1.0d0) )
        !    lbar_test(jc,2) = sum(Phi(1:nty,1:ns,1:na,jc,1:n_ofsh), mask=(lfun(1:nty,1:ns,1:na,jc,1:n_ofsh)>2.0d0) )
        !    lbar_test(jc,3) = sum(Phi(1:nty,1:ns,1:na,jc,1:n_ofsh), mask=(lfun(1:nty,1:ns,1:na,jc,1:n_ofsh)>5.0d0) )
        !    lbar_test(jc,4) = sum(Phi(1:nty,1:ns,1:na,jc,1:n_ofsh), mask=(lfun(1:nty,1:ns,1:na,jc,1:n_ofsh)>10.0d0) )
        !    !write(15, '(i3, 4f9.6)') jc, lbar_test(jc, :)
        !    do tyc=1,nty
        !        do jj=1,n_ofsh
        !            labartype(jc,tyc,jj)=0.0
        !            do sc=1,ns
        !                labar(jc)=labar(jc)+ sum(Phi(tyc,sc,1:na,jc,jj)*lfun(tyc,sc,1:na,jc,jj))
        !                labartype(jc,tyc,jj)=labartype(jc,tyc,jj)+ sum(Phi(tyc,sc,1:na,jc,jj)*lfun(tyc,sc,1:na,jc,jj))/measty(tyc)
        !            end do
        !        end do
        !    end do
        !end do
        !!close(15)
        !
        !
        !! Compute average hours worked
        !
        !hours=sum(Nu(1:jr-1)*labar(1:jr-1))/sum(Nu(1:jr-1))
        !
        !! Compute Income Statistics over the Life Cycle
        !
        !!open(20, file='inc_distibution.txt')
        !
        !do jc=1,jr-1
        !
        !    probb(jc)=0.0
        !    meanearn(jc)=0.0
        !    logmean(jc)=0.0
        !    varlogearn(jc)=0.0
        !
        !    do tyc=1,nty
        !        do jj=1,n_ofsh
        !            meanearntype(jc,tyc,jj)=0.0
        !            logmeantype(jc,tyc,jj)=0.0
        !
        !            do sc=1,ns
        !                do ac=1,na
        !                    ! Labor Income before taxes (this is what STY estimate their process on)
        !                    earnings(tyc,sc,ac,jc,jj)=eta(sc)*ep(tyc,jc)*lfun(tyc,sc,ac,jc,jj)
        !                    meanearn(jc)=meanearn(jc)+earnings(tyc,sc,ac,jc,jj)*Phi(tyc,sc,ac,jc,jj)
        !                    meanearntype(jc,tyc,jj)=meanearntype(jc,tyc,jj)+earnings(tyc,sc,ac,jc,jj)*Phi(tyc,sc,ac,jc,jj)/measty(tyc)/psi_prob(jj)
        !
        !                    if ( earnings(tyc,sc,ac,jc,jj) > 0.00001 ) then 
        !
        !                        logearnings(tyc,sc,ac,jc,jj)=log(earnings(tyc,sc,ac,jc,jj))
        !                        logmean(jc)=logmean(jc)+logearnings(tyc,sc,ac,jc,jj)*Phi(tyc,sc,ac,jc,jj)
        !                        logmeantype(jc,tyc,jj)=logmeantype(jc,tyc,jj)+logearnings(tyc,sc,ac,jc,jj)*Phi(tyc,sc,ac,jc,jj)/measty(tyc)/psi_prob(jj)
        !                        probb(jc)=probb(jc)+Phi(tyc,sc,ac,jc,jj)
        !
        !                    endif
        !
        !                end do
        !            end do
        !        end do
        !    end do
        !
        !    logmean(jc)=logmean(jc)/probb(jc)
        !
        !    ! Variance of Log-Earnings
        !
        !    do tyc=1,nty
        !        do jj=1,n_ofsh
        !            varlogearntype(jc,tyc,jj)=0.0
        !            do sc=1,ns
        !                do ac=1,na
        !                    if ( ( earnings(tyc,sc,ac,jc,jj) > 0.00001 ) .and. (Phi(tyc,sc,ac,jc,jj) > 0.0) )then 
        !                        varlogearn(jc)=varlogearn(jc)+ Phi(tyc,sc,ac,jc,jj)*(logearnings(tyc,sc,ac,jc,jj)-logmean(jc))**2.0
        !                        varlogearntype(jc,tyc,jj)=varlogearntype(jc,tyc,jj)+ (Phi(tyc,sc,ac,jc,jj)/measty(tyc)/psi_prob(jj))*(logearnings(tyc,sc,ac,jc,jj)-logmeantype(jc,tyc,jj))**2.0
        !                    endif
        !                end do
        !            end do
        !        end do
        !    end do
        !
        !    varlogearn(jc)=varlogearn(jc)/probb(jc)
        !    !write(20, '(i3, 3f15.9)') jc, probb(jc), logmean(jc), varlogearn(jc)
        !
        !end do
        !
        !!close(20)
        !
        !! Computing aggregate Consumption
        !
        !do jc=1,J
        !    cbar(jc)=0.0
        !    do tyc=1,nty
        !        do jj=1,n_ofsh
        !            cbartype(jc,tyc,jj)=0.0
        !            do sc=1,ns
        !                cbar(jc)=cbar(jc)+ sum(Phi(tyc,sc,1:na,jc,jj)*cfun(tyc,sc,1:na,jc,jj))
        !                cbartype(jc,tyc,jj)=cbartype(jc,tyc,jj)+ sum(Phi(tyc,sc,1:na,jc,jj)*cfun(tyc,sc,1:na,jc,jj))/measty(tyc)
        !            end do
        !        end do
        !    end do
        !end do
        !
        !! Computing total taxes paid
        !
        !Totinctax=0.0
        !TotSStax=0.0
        !
        !do jc=1,J
        !    do tyc=1,nty
        !        do jj=1,n_ofsh
        !            do sc=1,ns
        !                do ac=1,na
        !                    if (jc<jr) then
        !                        earnl0=w*eta(sc)*ep(tyc,jc)*lfun(tyc,sc,ac,jc,jj)    
        !                        !earnl=earnl0-0.5*taup*MIN(maxSS,earnl0)  
        !                        earnl=earnl0  !-0.5*taup*MIN(maxSS,earnl0)  
        !
        !                    elseif (jc>=jr) then
        !                        earnl0=0.0
        !                        earnl=0.0 
        !                    end if
        !
        !                    acur = grida(ac)
        !                    if (acur < a_cutoff(jj)) then
        !                        atilde = acur*qd 
        !                    else
        !                        atilde = acur*qd*(1d0-frac_ofsh)
        !                    end if
        !                    earncap = (atilde + TrB)*r
        !                    ap = afun(tyc,sc,ac,jc,jj)
        !                    !wealth_tot(tyc,sc,ac,jc,jj) = ainit + TrB
        !                    wealth_tot(tyc,sc,ac,jc,jj) = ap + TrB
        !                    if (ap < a_cutoff(jj)) then
        !                        wealth_obs(tyc,sc,ac,jc,jj) = ap + TrB
        !                    else
        !                        wealth_obs(tyc,sc,ac,jc,jj) = ap*(1d0-frac_ofsh) + TrB
        !                    end if
        !                    tot_income(tyc,sc,ac,jc,jj) = earnl + earncap
        !                    Totinctax=Totinctax+Nu(jc)*Phi(tyc,sc,ac,jc,jj)*Tax4(earncap,earnl,jj)
        !
        !                    !TotSStax=TotSStax+Nu(jc)*Phi(tyc,sc,ac,jc)*taup*MIN(maxSS,earnl0)  
        !                    TotSStax=TotSStax+Nu(jc)*Phi(tyc,sc,ac,jc,jj)*taup*earnl0  
        !
        !                end do
        !            end do
        !        end do
        !    end do
        !end do
        !
        !
        !LabS = sum( Nu(1:J)*lbar(1:J) )
        !
        !As = sum( Nu(1:J)*abar(1:J) )
        !
        !Astart = sum(Nu(1:J)*astartbar(1:J) )
        !
        !Trstart = sum(Nu(1:J))*TrB
        !
        !C = sum( Nu(1:J)*cbar(1:J) )
        !
        !
        !do jc=1,J
        !    do jj=1,n_ofsh
        !        do sc=1,ns
        !            !Trn(sc,jc,jj)=Nu(jc)*(1.0-surv(jc))*sum(Phi(1:nty,sc,1:na,jc,jj)*afun(1:nty,sc,1:na,jc,jj), mask=afun(1:nty,sc,1:na,jc,jj) < a_cutoff(jj))
        !            Trn(sc,jc,jj)=Nu(jc)*(1.0-surv(jc))*sum(Phi(1:nty,sc,1:na,jc,jj)*afun(1:nty,sc,1:na,jc,jj))
        !        end do
        !    end do
        !end do
        !
        !Tr=0.0
        !Tr=(1.0+r)*sum(Trn)/(1.0+nn) 
        !Transagg=(1.0+r)*sum(Trn)
        !
        !TrBn=(sum(Trn)/(1.0+nn))/sum(Nu) ! accidental bequests lump-sum transferred to each agent
        !
        !
        !SSn=TotSStax/sum(NU(jr:J))
        !
        !!call toc
        !
        !
        !!! Save Distribution
        !!open(unit=81,file='ADis.txt')
        !!do ac=1,na
        !    !write(81,fmt='(2f12.6)') grida(ac,1), ADis(ac) !, a_cutoff
        !!end do
        !!close(81)
        !
        !
        !! Earnings distribution details
        !ntot = nty*ns*na*J*n_ofsh
        !
        !do jc = 1, Jr-1
        !    Phi_work(1:nty,1:ns,1:na,jc,1:n_ofsh) = Phi(1:nty,1:ns,1:na,jc,1:n_ofsh)*Nu(jc)/sum(Nu(1:Jr-1))
        !end do
        !
        !testPhi = sum(Phi_work)
        !print *, testPhi
        !
        !do sc = 1, ns
        !    DistXW(sc) = sum(Phi_work(1:nty,sc,1:na,1:Jr-1,1:n_ofsh))
        !end do
        !
        !print *, sum(DistXW)         
        !
        !do jc = 1, J
        !    Phi(1:nty,1:ns,1:na,jc,1:n_ofsh) = Phi(1:nty,1:ns,1:na,jc,1:n_ofsh)*Nu(jc)/sum(Nu)
        !end do
        !
        !testPhi = sum(Phi)
        !print *, testPhi
        !
        !do sc = 1, ns
        !    DistX(sc) = sum(Phi(1:nty,sc,1:na,1:J,1:n_ofsh))
        !end do
        !
        !print *, sum(DistX)
        !
        !
        !
        !!test_phi = sum( Phi )
        !Share_Offshoring = 0d0
        !do jj=1,n_ofsh
        !    Share_Offshoring  = Share_Offshoring + sum( Phi(1:nty,1:ns,1:na,1:J, jj), mask=afun(1:nty,1:ns,1:na,1:J,jj) >= a_cutoff(jj))
        !end do
        !
        !print *, 'Share offshoring: ', Share_Offshoring
        !
        !! Total income distribution 
        !tot_income_1d = reshape(tot_income, [ntot])
        !phi_1d = reshape(phi, [ntot])
        !lorenz_inc = tot_income_1d
        !call lorenz(phi_1d, lorenz_inc, fx, gini_inc)  
        !frac_below = LinInterp_1d(99.99d0/100d0,fx,lorenz_inc,ntot)
        !inc_top_001 = 1d0-frac_below
        !frac_below = LinInterp_1d(99.9d0/100d0,fx,lorenz_inc,ntot)
        !frac_above = inc_top_001
        !inc_top_01_001 = 1d0-frac_above-frac_below
        !frac_below = LinInterp_1d(99d0/100d0,fx,lorenz_inc,ntot)
        !frac_above = inc_top_001 + inc_top_01_001
        !inc_top_1_01 = 1d0-frac_above-frac_below
        !frac_below = LinInterp_1d(90d0/100d0,fx,lorenz_inc,ntot)
        !frac_above = inc_top_001 + inc_top_01_001 + inc_top_1_01
        !inc_top_10_1 = 1d0-frac_above-frac_below
        !inc_0_90 = frac_below
        !test1 = inc_0_90 + inc_top_10_1 + inc_top_1_01 + inc_top_01_001 + inc_top_001 
        !!open(21, file='inc_distr.txt')
        !!do i = 1, ntot
        !!    write(21, '(2f12.6)') fx(i), lorenz_inc(i)
        !!end do
        !!close(21)
        !
        !! Labor income distribution 
        !lab_income_1d = reshape(earnings, [ntot])
        !!phi_1d = reshape(phi, [ntot])
        !lorenz_inc = lab_income_1d
        !call lorenz(phi_1d, lorenz_inc, fx, gini_lab_inc)  
        !frac_below = LinInterp_1d(99.99d0/100d0,fx,lorenz_inc,ntot)
        !lab_inc_top_001 = 1d0-frac_below
        !frac_below = LinInterp_1d(99.9d0/100d0,fx,lorenz_inc,ntot)
        !frac_above = inc_top_001
        !lab_inc_top_01_001 = 1d0-frac_above-frac_below
        !frac_below = LinInterp_1d(99d0/100d0,fx,lorenz_inc,ntot)
        !frac_above = inc_top_001 + inc_top_01_001
        !lab_inc_top_1_01 = 1d0-frac_above-frac_below
        !frac_below = LinInterp_1d(90d0/100d0,fx,lorenz_inc,ntot)
        !frac_above = inc_top_001 + inc_top_01_001 + inc_top_1_01
        !lab_inc_top_10_1 = 1d0-frac_above-frac_below
        !lab_inc_0_90 = frac_below
        !test1 = lab_inc_0_90 + lab_inc_top_10_1 + lab_inc_top_1_01 + lab_inc_top_01_001 + lab_inc_top_001        
        !
        !! Observed wealth distribution
        !wealth_obs_1d = reshape(wealth_obs, [ntot])
        !!phi_1d = reshape(phi, [ntot])
        !lorenz_wealth = wealth_obs_1d
        !call lorenz(phi_1d, lorenz_wealth, fx, gini_wealth_obs)  
        !frac_below = LinInterp_1d(0.9999d0,fx,lorenz_wealth,ntot)
        !wealth_obs_top_001 = 1d0-frac_below
        !frac_below = LinInterp_1d(0.999d0,fx,lorenz_wealth,ntot)
        !frac_above = wealth_obs_top_001
        !wealth_obs_top_01_001 = 1d0-frac_above-frac_below
        !frac_below = LinInterp_1d(0.99d0,fx,lorenz_wealth,ntot)
        !frac_above = wealth_obs_top_001 + wealth_obs_top_01_001
        !wealth_obs_top_1_01 = 1d0-frac_above-frac_below
        !frac_below = LinInterp_1d(0.90d0,fx,lorenz_wealth,ntot)
        !frac_above = wealth_obs_top_001 + wealth_obs_top_01_001 + wealth_obs_top_1_01
        !wealth_obs_top_10_1 = 1d0-frac_above-frac_below
        !wealth_obs_0_90 = frac_below
        !test2 = wealth_obs_0_90 + wealth_obs_top_10_1 + wealth_obs_top_1_01 +  wealth_obs_top_001 + wealth_obs_top_01_001
        !!open(22, file='wealth_obs_distr.txt')
        !!do i = 1, ntot
        !!    write(22, '(2f12.6)') fx(i), lorenz_wealth(i)
        !!end do
        !!close(22)        
        !
        !! Total wealth distribution
        !wealth_tot_1d = reshape(wealth_tot, [ntot])
        !!phi_1d = reshape(phi, [ntot])
        !lorenz_wealth = wealth_tot_1d
        !call lorenz_mod(phi_1d, lorenz_wealth, fx, gini_wealth_tot, wealth_tot_sorted_1d, phi_1d_sorted)  
        !frac_below = LinInterp_1d(0.9999d0,fx,lorenz_wealth,ntot)
        !if (frac_below > 1d0) then
        !    if (save_res) then
        !        open(61, file='Lorenz_wealth_test.txt')
        !        do i = 1, ntot
        !            write(61, '(4f35.16)') fx(i), lorenz_wealth(i), wealth_tot_sorted_1d(i), phi_1d_sorted(i)
        !        end do
        !        close(61)
        !    end if
        !    frac_below = LinInterp_1d(0.9999d0,fx,lorenz_wealth,ntot)
        !end if
        !wealth_tot_top_001 = 1d0-frac_below
        !frac_below = LinInterp_1d(0.999d0,fx,lorenz_wealth,ntot)
        !frac_above = wealth_tot_top_001
        !wealth_tot_top_01_001 = 1d0-frac_above-frac_below
        !frac_below = LinInterp_1d(0.99d0,fx,lorenz_wealth,ntot)
        !frac_above = wealth_tot_top_001 + wealth_tot_top_01_001
        !wealth_tot_top_1_01 = 1d0-frac_above-frac_below
        !frac_below = LinInterp_1d(0.90d0,fx,lorenz_wealth,ntot)
        !frac_above = wealth_tot_top_001 + wealth_tot_top_01_001 + wealth_tot_top_1_01
        !wealth_tot_top_10_1 = 1d0-frac_above-frac_below
        !wealth_tot_0_90 = frac_below
        !test3 = wealth_tot_0_90 + wealth_tot_top_10_1 + wealth_tot_top_1_01 +  wealth_tot_top_001 + wealth_tot_top_01_001
        !!open(22, file='wealth_tot_distr.txt')
        !!do i = 1, ntot
        !!    write(22, '(2f12.6)') fx(i), lorenz_wealth(i)
        !!end do
        !!close(22)         
        !
        !wealth_tot_cutoff_top_001 = LinInterp_1d(0.9999d0,fx,wealth_tot_sorted_1d,ntot)
        !wealth_tot_cutoff_top_01 = LinInterp_1d(0.999d0,fx,wealth_tot_sorted_1d,ntot)
        !wealth_tot_cutoff_top_1 = LinInterp_1d(0.99d0,fx,wealth_tot_sorted_1d,ntot)
        !wealth_tot_cutoff_top_10 = LinInterp_1d(0.90d0,fx,wealth_tot_sorted_1d,ntot)
        !
        !sort_key=[ (i,i=1,ntot) ]
        !CALL dlasrt2('I',ntot,wealth_tot_1d,sort_key,info)
        !phi_1d_sorted = phi_1d(sort_key)
        !!open(22, file='wealth_tot_sorted.txt')
        !!do i = 1, ntot
        !!    write(22, '(f10.6, f25.15)') phi_1d_sorted(i), wealth_tot_1d(i)
        !!end do
        !!close(22)          
        !
        !i = 1
        !phi_tot = phi_1d_sorted(1)
        !do while (phi_tot <= 0.9d0)
        !    i = i + 1
        !    phi_tot = phi_tot + phi_1d_sorted(i)
        !end do
        !
        !test1 = (wealth_tot_1d(i) + wealth_tot_1d(i-1))/2d0
        !
        !
        !
        !test2 = sum(phi_1d_sorted, mask=(wealth_tot_1d < test1))
        !test3 = sum(phi, mask=(wealth_tot < test1))
        !test4 = sum(phi_1d)
        !
        !do while (phi_tot <= 0.99d0)
        !    i = i + 1
        !    phi_tot = phi_tot + phi_1d_sorted(i)
        !end do
        !test2 = (wealth_tot_1d(i) + wealth_tot_1d(i-1))/2d0
        !
        !test3 = sum(phi, mask=(wealth_tot > test1 .and. wealth_tot <= test2))
        !
        !do while (phi_tot <= 0.999d0)
        !    i= i + 1
        !    phi_tot = phi_tot + phi_1d_sorted(i)
        !end do
        !test3 = (wealth_tot_1d(i) + wealth_tot_1d(i-1))/2d0
        !
        !test4 = sum(phi, mask=(wealth_tot > test2 .and. wealth_tot <= test3))
        !
        !do while (phi_tot <= 0.9999d0)
        !    i= i + 1
        !    phi_tot = phi_tot + phi_1d_sorted(i)
        !end do
        !test4 = (wealth_tot_1d(i) + wealth_tot_1d(i-1))/2d0   
        !
        !test5 = sum(phi, mask=(wealth_tot > test4))
        !
        !
        !! Share offshoring
        !share_off_top_001 = 0d0
        !share_off_top_01_001 = 0d0
        !share_off_top_1_01 = 0d0
        !share_off_top_10_1 = 0d0
        !share_off_0_90 = 0d0
        !prob_top_001 = 0d0
        !prob_top_01_001 = 0d0
        !prob_top_1_01 = 0d0
        !prob_top_10_1 = 0d0
        !prob_0_90 = 0d0  
        !
        !test1 = 0d0
        !do i = 1, ntot
        !    if (wealth_tot_1d(i) <= wealth_tot_cutoff_top_10) then
        !        test1 = test1 + phi_1d_sorted(i)
        !    end if
        !end do
        !print *, test1
        !
        !
        !do jc=1,J
        !    do tyc=1,nty
        !        do jj=1,n_ofsh
        !            do sc=1,ns
        !                do ac=1,na
        !                    ap = afun(tyc,sc,ac,jc,jj)
        !                    !if (ap > a_cutoff(jj)) then
        !                        if (wealth_tot(tyc,sc,ac,jc,jj) <= wealth_tot_cutoff_top_10) then
        !                            prob_0_90 = prob_0_90 + phi(tyc,sc,ac,jc,jj)
        !                            if (ap > a_cutoff(jj)) then
        !                                share_off_0_90 = share_off_0_90 + phi(tyc,sc,ac,jc,jj)
        !                            end if
        !                        else if (wealth_tot(tyc,sc,ac,jc,jj) > wealth_tot_cutoff_top_10 .and. wealth_tot(tyc,sc,ac,jc,jj) <= wealth_tot_cutoff_top_1) then
        !                            prob_top_10_1 = prob_top_10_1 + phi(tyc,sc,ac,jc,jj)
        !                            if (ap > a_cutoff(jj)) then
        !                                share_off_top_10_1 = share_off_top_10_1 + phi(tyc,sc,ac,jc,jj) 
        !                            end if
        !                        else if (wealth_tot(tyc,sc,ac,jc,jj) > wealth_tot_cutoff_top_1 .and. wealth_tot(tyc,sc,ac,jc,jj) <= wealth_tot_cutoff_top_01) then
        !                            prob_top_1_01 = prob_top_1_01 + phi(tyc,sc,ac,jc,jj)
        !                            if (ap > a_cutoff(jj)) then
        !                                share_off_top_1_01 = share_off_top_1_01 + phi(tyc,sc,ac,jc,jj)  
        !                            end if
        !                        else if (wealth_tot(tyc,sc,ac,jc,jj) > wealth_tot_cutoff_top_01 .and. wealth_tot(tyc,sc,ac,jc,jj) <= wealth_tot_cutoff_top_001) then
        !                            prob_top_01_001 = prob_top_01_001 + phi(tyc,sc,ac,jc,jj)
        !                            if (ap > a_cutoff(jj)) then
        !                                share_off_top_01_001 = share_off_top_01_001 + phi(tyc,sc,ac,jc,jj)
        !                            end if
        !                        else 
        !                            prob_top_001 = prob_top_001 + phi(tyc,sc,ac,jc,jj)
        !                            if (ap > a_cutoff(jj)) then
        !                                share_off_top_001 = share_off_top_001 + phi(tyc,sc,ac,jc,jj) 
        !                            end if
        !                        end if
        !                    !end if
        !                end do
        !            end do
        !        end do
        !    end do
        !end do
        !test1 = prob_top_001 - 0.0001d0
        !test2 = prob_top_01_001 - (0.9999d0 - 0.999d0)
        !test3 = prob_top_1_01 - (0.999d0 - 0.99d0)
        !test4 = prob_0_90 - 0.90d0
        !
        !share_off_top_001 = share_off_top_001/prob_top_001
        !share_off_top_01_001 = share_off_top_01_001/prob_top_01_001
        !share_off_top_1_01 = share_off_top_1_01/prob_top_1_01
        !share_off_top_10_1 = share_off_top_10_1/prob_top_10_1
        !share_off_0_90 = share_off_0_90/prob_0_90              
        !
        !share_off_dist=[share_off_0_90, share_off_top_10_1, share_off_top_1_01, share_off_top_01_001, share_off_top_001]
        !wealth_obs_dist=[wealth_obs_0_90, wealth_obs_top_10_1, wealth_obs_top_1_01, wealth_obs_top_01_001, wealth_obs_top_001, gini_wealth_obs]
        !wealth_tot_dist=[wealth_tot_0_90, wealth_tot_top_10_1, wealth_tot_top_1_01, wealth_tot_top_01_001, wealth_tot_top_001, gini_wealth_tot]
        !income_dist=[inc_0_90, inc_top_10_1, inc_top_1_01, inc_top_01_001, inc_top_001, gini_inc]
        !lab_income_dist=[lab_inc_0_90, lab_inc_top_10_1, lab_inc_top_1_01, lab_inc_top_01_001, lab_inc_top_001, gini_lab_inc]
        !
        !open(71, file='income_wealth_ExpKLP.tex')
        !write(71, *) '\documentclass{article}'
        !write(71, *) '\usepackage{amsmath}'
        !write(71, *) '\begin{document}'
        !write(71, *) '\begin{table}[h!]'
        !write(71, *) '\centering'
        !write(71, *) '\begin{tabular}{l|c|c|c|c|c}'
        !write(71, *) '\hline'
        !write(71, *) ' Moment & Income & Labor income & Wealth (obs) & Wealth (tot) & Offshoring \\'
        !write(71, *) '\hline'
        !write(71, '(a15, 4(f8.3, a3), f8.3, a5)') 'P0-90 & ', inc_0_90, '&', lab_inc_0_90, '&', wealth_obs_0_90, '&', wealth_tot_0_90, '&', share_off_0_90, ' \\'
        !write(71, '(a15, 4(f8.3, a3), f8.3, a5)') 'P90-99 & ', inc_top_10_1, '&', lab_inc_top_10_1, '&', wealth_obs_top_10_1, '&', wealth_tot_top_10_1, '&', share_off_top_10_1, ' \\'
        !write(71, '(a15, 4(f8.3, a3), f8.3, a5)') 'P99-99.9 & ', inc_top_1_01, '&', lab_inc_top_1_01, '&', wealth_obs_top_1_01, '&', wealth_tot_top_1_01, '&', share_off_top_1_01,  ' \\'
        !write(71, '(a15, 4(f8.3, a3), f8.3, a5)') 'P99.9-99.99 & ', inc_top_01_001, '&', lab_inc_top_01_001, '&', wealth_obs_top_01_001, '&', wealth_tot_top_01_001, '&', share_off_top_01_001, ' \\'
        !write(71, '(a15, 4(f8.3, a3), f8.3, a5)') 'P99.99-100 & ', inc_top_001, '&', lab_inc_top_001, '&', wealth_obs_top_001, '&', wealth_tot_top_001, '&', share_off_top_001, ' \\'
        !write(71, '(a15, 3(f8.3, a3), f8.3, a9)') 'Gini & ', gini_inc, '&',  gini_lab_inc, '&', gini_wealth_obs, '&', gini_wealth_tot, '& NA \\'
        !write(71, *) '\hline'
        !write(71, *) '\end{tabular}'
        !write(71, *) '\caption{Earnings/income/wealth distributions}'
        !write(71, *) '\end{table}'
        !!
        !!write(71, *) '\begin{table}[h!]'
        !!write(71, *) '\centering'
        !!write(71, *) '\begin{tabular}{lcccccccc}'
        !!write(71, *) '\hline'
        !!write(71, *) ' & $z_1$ & $z_2$ & $z_3$ & $z_4$ & $z_5$ & $z_6$ & $z_7$ & $z_8$ \\'
        !!write(71, *) '\hline'
        !!do i = 1, ns
        !!    write(71, '(a, f8.2, <ns>(a, f8.3), a)') '$z_1$ = ', eta(i), ' & ', pi(i,1) , ' & ', pi(i,2), ' & ', pi(i,3), ' & ', pi(i,4), ' & ', pi(i,5), ' & ', pi(i,6), ' & ', pi(i,7), ' & ', pi(i,8), ' \\' 
        !!end do
        !!write(71, *) '\hline'
        !!write(71, '(<ns>(a, f8.3), a)') 'stationary distr. &', pstat(1,1), ' & ', pstat(1,2), ' & ', pstat(1,3), ' & ', pstat(1,4), ' & ', pstat(1,5), ' & ', pstat(1,6), ' & ', pstat(1,7), ' & ', pstat(1,8), ' \\'
        !!write(71, *) '\hline'
        !!write(71, *) '\end{tabular}'
        !!write(71, *) '\caption{Shock to productivity process}'
        !!write(71, *) '\end{table}'
        !!write(71, *) '\end{document}'
        !
        !
        !close(71)
        !
        !!sim_moms = [wealth_obs_0_90, wealth_obs_top_10_1, wealth_obs_top_1_01, wealth_obs_top_01_001, wealth_obs_top_001, gini_wealth_obs, &
        !!            inc_0_90, inc_top_10_1, inc_top_1_01, inc_top_01_001, inc_top_001, gini_inc, &
        !!            share_off_0_90, share_off_top_10_1, share_off_top_1_01, share_off_top_01_001, share_off_top_001, &
        !!            As/((1.0+nn)*Y), hours]
        !sim_moms = [wealth_obs_top_01_001, wealth_obs_top_001, gini_wealth_obs, inc_top_01_001, inc_top_001, gini_inc]     
        !sim_moms_aux = [share_off_0_90, share_off_top_10_1, share_off_top_1_01, share_off_top_01_001, share_off_top_001, As/((1.0+nn)*Y), hours]
        !sim_moms_klp = [wealth_obs_top_1_01 + wealth_obs_top_01_001 + wealth_obs_top_001,  wealth_obs_top_01_001 + wealth_obs_top_001, gini_wealth_obs, &
        !                inc_top_1_01 + inc_top_01_001 + inc_top_001, inc_top_01_001 + inc_top_001, gini_inc]

    end subroutine Distribution
    
end module Mod_Distribution
