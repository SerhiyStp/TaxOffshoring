module Mod_Household

    implicit none

    real(8) :: cmod, lambda_mod, nonlabinc_mod

contains

    subroutine SolveHH(save_res)
    use params
    use MyLinInterp 
    use zbren_int
    use int_tictoc
        logical :: save_res
        integer :: ia, nq, scp
        real(8) :: cons, lab
        real(8) :: vp, Vnext
        real(8) :: acutoff(n_ofsh), mcutoff(n_ofsh)
        real(8) :: m0
        real(8) :: gridm(na,n_ofsh)
        real(8) :: mtmp(na), ctmp(na), htmp(na), vtmp(na), atmp(na)
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
        real(8) :: test_bc(nty,ns,na,J)
        real(8) :: test_hrs(nty,ns,na,J)
        real(8) :: max_test_bc
        real(8) :: vals(2), testV, testH, testC, testA
        integer :: inds(2)
        integer :: jj
        real(8) :: fff_static, fff_helper_func1
        
        !external foc_static
        !external helper_func1
        external fff_static
        external fff_helper_func1
        
        IERSVR = 0
        IPACT = 1
        ISACT = 0
        call erset(IERSVR, IPACT, ISACT)  
        
        !acutoff = psi_vals/frac_ofsh/(qd-qw)
        acutoff = psi_vals/(qd - qwd)
        mcutoff = SS + acutoff + TrB*(1d0+r*(1d0-tk))
        
        a_cutoff = acutoff
        
        do jj = 1, n_ofsh
            if (acutoff(jj) < grida_tmp1(na-2)) then
                ia = 1
                do while (grida_tmp1(ia) <= acutoff(jj))
                    ia = ia + 1
                end do
                ia = ia - 1
                grida(1:ia,jj) = grida_tmp1(1:ia)
                grida(ia+1,jj) = acutoff(jj)-1d-3
                grida(ia+2,jj) = acutoff(jj)+1d-3
                grida(ia+3:na,jj) = grida_tmp1(ia+1:na-2)   
            else
                print *, 'WARNING: a_cutoff is not within grida'
                grida(:,jj) = grida_tmp2
            end if
        end do
        

        gridm = SS + grida + TrB*(1d0+r*(1d0-tk))
        ! Last generation
        lab = 0.0d0
        do jj = 1, n_ofsh
            do ia = 1, na
                cons = gridm(ia,jj)/(1.0d0+tauc)
                Vfun(1:nty,1:ns,ia,J,jj)=U(cons,lab)
                lfun(1:nty,1:ns,ia,J,jj)=lab
                cfun(1:nty,1:ns,ia,J,jj)=cons
                afun(1:nty,1:ns,ia,J,jj)=0.0d0
                vpfun(1:nty,1:ns,ia,J,jj)=marginal_utility(cons,lab)
            end do
        end do

    
        do jj = 1, n_ofsh

            ! Retirement part
            do jc=J-1,jr,-1 
                ! Invert Euler equation to find current consumption
                do ia = 1, na
                    vp = bbeta*surv(jc)*vpfun(1,1,ia,jc+1,jj)
                    Vnext = Vfun(1,1,ia,jc+1,jj)
                    if (grida(ia,jj) <= acutoff(jj)) then
                        ! Not offshoring
                        cons = (vp/qd)**(-1.0d0/sig1)
                        m0 = grida(ia,jj)*qd + cons*(1.0d0+tauc)
                    else
                        ! Offshoring
                        cons = (vp/qwd)**(-1.0d0/sig1)
                        m0 = grida(ia,jj)*qwd + psi_offshore + cons*(1.0d0+tauc)
                    end if
                    mtmp(ia) = m0
                    ctmp(ia) = cons
                    vtmp(ia) = -1d0 / (U(cons,0d0) + bbeta*surv(jc)*Vnext)
                    !vtmp(ia) = (U(cons,0d0) + bbeta*surv(jc)*Vnext)
                    atmp(ia) = grida(ia,jj)
                end do
                
                if (save_res) then
                    open(20, file='mtmp.txt')
                    do ia=1,na
                        write(20,'(4f16.6)') grida(ia,jj), mtmp(ia), ctmp(ia), vtmp(ia)
                    end do
                    close(20)
                end if
                
                mdiff = mtmp(2:na) - mtmp(1:na-1)
                nkinks = count( mdiff < 0d0 )
                
                !call tic()
                
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
                        if (gridm(ia,jj) <= mtmp(1)) then
                            aprime = grida(1,jj)
                            if (aprime <= acutoff(jj)) then
                                cons = (gridm(ia,jj) - aprime*qd)/(1.0d0+tauc)
                            else
                                cons = (gridm(ia,jj) - aprime*qwd - psi_offshore)/(1.0d0+tauc)
                            end if
                            VV = U(cons,0d0) + bbeta*surv(jc)*Vfun(1,1,1,jc+1,jj)
                            CommonVt(ia,:) = VV
                            CommonC(ia,:) = cons
                            CommonA(ia,:) = aprime
                        else 
                            do ii = 1, nkinks
                                if (gridm(ia,jj) > mtmp(na)) then
                                    ! If the point on the exo grid "gridm" falls outside the last point of the endo grid "mtmp", extrapolate using the previous points
                                    CommonC(ia,ii) = LinInterp_1d(gridm(ia,jj),gridm(1:ia-1,jj),CommonC(1:ia-1,ii),ia-1)   
                                    CommonA(ia,ii) = LinInterp_1d(gridm(ia,jj),gridm(1:ia-1,jj),CommonA(1:ia-1,ii),ia-1) 
                                    CommonVt(ia,ii) = LinInterp_1d(gridm(ia,jj),gridm(1:ia-1,jj),CommonVt(1:ia-1,ii),ia-1) 
                                else
                                    i1 = increase(ii)
                                    i2 = fall(ii)
                                    if (gridm(ia,jj) >= mtmp(i1) .and. gridm(ia,jj) <= mtmp(i2)) then
                                        !!CommonVt(ia,ii) = LinInterp_1d(gridm(ia),mtmp(i1:i2),vtmp(i1:i2),i2-i1+1)
                                        !CommonVt(ia,ii) = -1d0/LinInterp_1d(gridm(ia),mtmp(i1:i2),vtmp(i1:i2),i2-i1+1)
                                        !CommonC(ia,ii) = LinInterp_1d(gridm(ia),mtmp(i1:i2),ctmp(i1:i2),i2-i1+1)
                                        !CommonA(ia,ii) = LinInterp_1d(gridm(ia),mtmp(i1:i2),atmp(i1:i2),i2-i1+1)
                                    
                                        call basefun(mtmp(i1:i2),i2-i1+1,gridm(ia,jj),vals,inds)
                                        CommonC(ia,ii) = vals(1)*ctmp(i1-1+inds(1))+vals(2)*ctmp(i1-1+inds(2))
                                        CommonVt(ia,ii) = -1d0/(vals(1)*vtmp(i1-1+inds(1))+vals(2)*vtmp(i1-1+inds(2)))
                                        CommonA(ia,ii) = vals(1)*atmp(i1-1+inds(1))+vals(2)*atmp(i1-1+inds(2))    
                                    end if
                                end if
                            end do
                        end if
                    end do

                    idmax = maxloc(CommonVt, dim=2) 
                    do ia = 1, na
                        hrs = 0d0
                        lfun(1:nty,1:ns,ia,jc,jj) = hrs
                        cons = CommonC(ia,idmax(ia))
                        cfun(1:nty,1:ns,ia,jc,jj) = cons
                        afun(1:nty,1:ns,ia,jc,jj) = CommonA(ia,idmax(ia))
                        vpfun(1:nty,1:ns,ia,jc,jj) = marginal_utility(cons,hrs)
                        Vfun(1:nty,1:ns,ia,jc,jj) = CommonVt(ia,idmax(ia))
                    end do
                    deallocate(CommonVt, CommonC, CommonA)
                else
                    ! Interpolation step
                    do ia = 1, na
                        if (gridm(ia,jj) <= mtmp(1)) then
                            aprime = grida(1,jj)
                            if (aprime <= acutoff(jj)) then
                                cons = (gridm(ia,jj) - aprime*qd)/(1.0d0+tauc)
                            else
                                cons = (gridm(ia,jj) - aprime*qwd - psi_vals(jj))/(1.0d0+tauc)
                            end if
                            VV = U(cons,0d0) + bbeta*surv(jc)*Vfun(1,1,1,jc+1,jj)
                        else
                            !aprime = LinInterp_1d(gridm(ia),mtmp,grida,na)
                            !cons = LinInterp_1d(gridm(ia),mtmp,ctmp,na)
                            !VV = -1d0/LinInterp_1d(gridm(ia),mtmp,vtmp,na)
                            !!VV = LinInterp_1d(gridm(ia),mtmp,vtmp,na)
                            
                            call basefun(mtmp,na,gridm(ia,jj),vals,inds)
                            
                            !if (inds(2) == na) then
                            !    print *, 'WARNING: grida too small?'
                            !end if
                            
                            cons = vals(1)*ctmp(inds(1))+vals(2)*ctmp(inds(2))
                            VV = -1d0/(vals(1)*vtmp(inds(1))+vals(2)*vtmp(inds(2)))
                            aprime = vals(1)*atmp(inds(1))+vals(2)*atmp(inds(2))                           
                            
                        end if
                        afun(1:nty,1:ns,ia,jc,jj) = aprime
                        cfun(1:nty,1:ns,ia,jc,jj)=cons
                        lfun(1:nty,1:ns,ia,jc,jj)=0d0
                        vpfun(1:nty,1:ns,ia,jc,jj)=marginal_utility(cons,0d0)
                        Vfun(1:nty,1:ns,ia,jc,jj) = VV
                    end do
                end if
                
                !call toc()
                if (save_res == .true.) then
                    open(21, file='gridm.txt')
                    open(22, file='afun.txt')
                    open(23, file='cfun.txt')
                    open(24, file='Vfun.txt')
                    open(25, file='dVfun.txt')
                    open(26, file='am_cutoff.txt')
                    do ia=1,na
                        write(21,'(f26.6)') gridm(ia,jj)
                        write(22,'(f26.6)') afun(1,1,ia,jc,jj)
                        write(23,'(f26.6)') cfun(1,1,ia,jc,jj)
                        write(24,'(f26.6)') Vfun(1,1,ia,jc,jj)
                        write(25,'(f26.6)') vpfun(1,1,ia,jc,jj)
                    end do
                    write(26, '(2f16.6, i4)') acutoff(jj), mcutoff(jj), jc
                    close(21); close(22); close(23); close(24); close(25); close(26)
                end if
            end do
        
        
            ! Working age part
            gridm = grida + TrB*(1d0+r*(1d0-tk))
            mcutoff = acutoff + TrB*(1d0+r*(1d0-tk))
            !open(26, file='am_cutoff.txt')
            !write(26, '(2f16.6)') acutoff, mcutoff
            !close(26)
            do jc = jr-1, 1, -1
                !print *, 'jc = ', jc
                do tyc = 1, nty
                    do sc = 1, ns

                        ! Backward solution step
                        do ia = 1, na
                            vp = bbeta*surv(jc)*sum(pi(sc,:)*vpfun(tyc,:,ia,jc+1,jj))
                            Vnext = sum(pi(sc,:)*Vfun(tyc,:,ia,jc+1,jj))

                            if (grida(ia,jj) <= acutoff(jj)) then
                                ! Not offshoring
                                cons = (vp/qd)**(-1d0/sig1)
                                m0 = grida(ia,jj)*qd + cons*(1.0d0+tauc)
                            else
                                ! Offshoring
                                cons = (vp/qwd)**(-1d0/sig1)
                                m0 = grida(ia,jj)*qwd + psi_offshore + cons*(1.0d0+tauc)
                            end if
                            !cmod = cons
                            lambda_mod = cons**(-sig1)/(1d0+tauc)
                            fne = 100
                            hlo = 1d-9
                            if (sc < 6) then
                                hhi = 5d0
                            else if (sc == 6) then
                                hhi = 10d0
                            else 
                                hhi = 20d0
                            end if
                            flo = foc_static(hlo)
                            fhi = foc_static(hhi)
                            
                            if (flo*fhi <= 0d0) then
                                !call zbren(foc_static, hlo, hhi, errabs=1d-9, errrel=1d-9, maxfn=fne) 
                                call zbren(fff_static, hlo, hhi, errabs=1d-9, errrel=1d-9, maxfn=fne)
                            else
                                print *, 'WARNING: flo and fhi must have opposite signs!'
                                do while (flo*fhi > 0d0)
                                    hhi = hhi*2d0
                                    fhi = foc_static(hhi)
                                end do
                                !call zbren(foc_static, hlo, hhi, errabs=1d-9, errrel=1d-9, maxfn=fne)
                                call zbren(fff_static, hlo, hhi, errabs=1d-9, errrel=1d-9, maxfn=fne)
                            end if
                            hrs = hhi
                            labinc = w*eta(sc)*ep(tyc,jc)*hrs
                            m0 = m0 + taup*labinc - theta0*min(yb_cutoff, labinc)**(1d0-theta1) - (1.0d0-tau_max)*max(0d0, labinc-yb_cutoff)
                            mtmp(ia) = m0
                            ctmp(ia) = cons
                            htmp(ia) = hrs
                            atmp(ia) = grida(ia,jj)
                            !vtmp(ia) = U(cons,hrs) + bbeta*surv(jc)*Vnext
                            vtmp(ia) = -1d0 / (U(cons,hrs) + bbeta*surv(jc)*Vnext)
                        end do
                        !open(20, file='mtmp.txt')
                        !do ia=1,na
                            !write(20,'(5f20.9)') grida(ia,jj), mtmp(ia), ctmp(ia), htmp(ia), vtmp(ia)
                        !end do
                        !close(20)
                        
                        ! Interpolation step

                        mdiff = mtmp(2:na) - mtmp(1:na-1)
                        nkinks = count( mdiff < 0d0 )
                        !call tic()
                        
                        if (nkinks > 0) then
                            ! Upper envelope
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
                            allocate(CommonH(na,nkinks))
                            allocate(CommonC(na,nkinks))
                            allocate(CommonA(na,nkinks))
                            CommonVt = -10d0**8d0
                            do ia = 1, na
                                if (gridm(ia,jj) <= mtmp(1)) then
                                    aprime = grida(1,jj)
                                    if (aprime <= acutoff(jj)) then
                                        nonlabinc_mod = gridm(ia,jj) - aprime*qd
                                    else
                                        nonlabinc_mod = gridm(ia,jj) - aprime*qwd - psi_offshore
                                    end if
                                    hlo = 1d-4
                                    hhi = 5d0
                                    fne = 100
                                    flo = helper_func1(hlo)
                                    fhi = helper_func1(hhi)
                                    if (flo*fhi < 0d0) then
                                        !call zbren(helper_func1, hlo, hhi, errabs=1d-9, errrel=1d-9, maxfn=fne) 
                                        call zbren(fff_helper_func1, hlo, hhi, errabs=1d-9, errrel=1d-9, maxfn=fne) 
                                        hrs = hhi
                                    else
                                        print *, 'WARNING: helper_func1 must have different signs!'
                                    end if                        
                                    labinc = w*eta(sc)*ep(tyc,jc)*hrs
                                    if (labinc < yb_cutoff) then
                                        cons = ( theta0*(1d0-theta1)*(w*eta(sc)*ep(tyc,jc))**(1d0-theta1)*hrs**(-theta1) - taup*w*eta(sc)*ep(tyc,jc) )/(1d0+tauc)/chi/(hrs**sig2)
                                    else
                                        cons = ( (1.0d0-tau_max)*(w*eta(sc)*ep(tyc,jc)) - taup*w*eta(sc)*ep(tyc,jc) )/(1d0+tauc)/chi/(hrs**sig2)
                                    end if
                                    cons = cons**(1d0/sig1)
                                    VV = U(cons,hrs) + bbeta*surv(jc)*sum(pi(sc,:)*Vfun(tyc,:,1,jc+1,jj))
                                    CommonVt(ia,:) = VV
                                    CommonH(ia,:) = hrs
                                    CommonC(ia,:) = cons
                                    CommonA(ia,:) = aprime
                                else 
                                    do ii = 1, nkinks
                                        if (gridm(ia,jj) > mtmp(na)) then
                                            ! If the point on the exo grid "gridm" falls outside the last point of the endo grid "mtmp", extrapolate using the previous points
                                            CommonH(ia,ii) = LinInterp_1d(gridm(ia,jj),gridm(1:ia-1,jj),CommonH(1:ia-1,ii),ia-1)   
                                            CommonC(ia,ii) = LinInterp_1d(gridm(ia,jj),gridm(1:ia-1,jj),CommonC(1:ia-1,ii),ia-1)   
                                            CommonA(ia,ii) = LinInterp_1d(gridm(ia,jj),gridm(1:ia-1,jj),CommonA(1:ia-1,ii),ia-1) 
                                            CommonVt(ia,ii) = LinInterp_1d(gridm(ia,jj),gridm(1:ia-1,jj),CommonVt(1:ia-1,ii),ia-1) 
                                        else
                                            i1 = increase(ii)
                                            i2 = fall(ii)
                                            if ((gridm(ia,jj) >= mtmp(i1) .and. gridm(ia,jj) <= mtmp(i2)) .or. (ii == nkinks .and. gridm(ia,jj) > mtmp(i2))) then
                                                !!CommonVt(ia,ii) = LinInterp_1d(gridm(ia),mtmp(i1:i2),vtmp(i1:i2),i2-i1+1)
                                                !CommonVt(ia,ii) = -1d0/LinInterp_1d(gridm(ia),mtmp(i1:i2),vtmp(i1:i2),i2-i1+1)
                                                !CommonH(ia,ii) = LinInterp_1d(gridm(ia),mtmp(i1:i2),htmp(i1:i2),i2-i1+1)
                                                !CommonC(ia,ii) = LinInterp_1d(gridm(ia),mtmp(i1:i2),ctmp(i1:i2),i2-i1+1)
                                                !CommonA(ia,ii) = LinInterp_1d(gridm(ia),mtmp(i1:i2),atmp(i1:i2),i2-i1+1)
                                                
                                                call basefun(mtmp(i1:i2),i2-i1+1,gridm(ia,jj),vals,inds)
                                                CommonH(ia,ii) = vals(1)*htmp(i1-1+inds(1))+vals(2)*htmp(i1-1+inds(2))
                                                CommonC(ia,ii) = vals(1)*ctmp(i1-1+inds(1))+vals(2)*ctmp(i1-1+inds(2))
                                                CommonVt(ia,ii) = -1d0/(vals(1)*vtmp(i1-1+inds(1))+vals(2)*vtmp(i1-1+inds(2)))
                                                CommonA(ia,ii) = vals(1)*atmp(i1-1+inds(1))+vals(2)*atmp(i1-1+inds(2))
                                            end if
                                        end if
                                    end do
                                end if
                            end do
                            
                            idmax = maxloc(CommonVt, dim=2) 
                            do ia = 1, na
                                hrs = CommonH(ia,idmax(ia))
                                lfun(tyc,sc,ia,jc,jj) = hrs
                                cons = CommonC(ia,idmax(ia))
                                cfun(tyc,sc,ia,jc,jj) = cons
                                afun(tyc,sc,ia,jc,jj) = CommonA(ia,idmax(ia))
                                vpfun(tyc,sc,ia,jc,jj) = marginal_utility(cons,hrs)
                                Vfun(tyc,sc,ia,jc,jj) = CommonVt(ia,idmax(ia))
                                
                                !if (afun(tyc,sc,ia,jc,jj) >= grida(na,jj)) then
                                !    print *, 'WARNING: grida too small?'
                                !end if

                                labinc = w*eta(sc)*ep(tyc,jc)*hrs
                                test_bc(tyc,sc,ia,jc) = gridm(ia,jj) + theta0*min(yb_cutoff,labinc)**(1d0-theta1) + (1.0d0-tau_max)*max(0d0,labinc-yb_cutoff)  - cons*(1d0+tauc) - taup*labinc
                                if (afun(tyc,sc,ia,jc,jj) <= acutoff(jj)) then 
                                    test_bc(tyc,sc,ia,jc) = test_bc(tyc,sc,ia,jc) - afun(tyc,sc,ia,jc,jj)*qd 
                                else
                                    test_bc(tyc,sc,ia,jc) = test_bc(tyc,sc,ia,jc) - afun(tyc,sc,ia,jc,jj)*qwd - psi_offshore
                                end if

                            end do
                            deallocate(CommonVt, CommonH, CommonC, CommonA)

                        else                

                            do ia = 1, na
                                if (gridm(ia,jj) <= mtmp(1)) then
                                    aprime = grida(1,jj)
                                    if (aprime <= acutoff(jj)) then
                                        nonlabinc_mod = gridm(ia,jj) - aprime*qd
                                    else
                                        nonlabinc_mod = gridm(ia,jj) - aprime*qwd - psi_offshore
                                    end if
                                    hlo = 1d-4
                                    hhi = 5d0
                                    fne = 100
                                    flo = helper_func1(hlo)
                                    fhi = helper_func1(hhi)
                                    if (flo*fhi < 0d0) then
                                        !call zbren(helper_func1, hlo, hhi, errabs=1d-9, errrel=1d-9, maxfn=fne) 
                                        call zbren(fff_helper_func1, hlo, hhi, errabs=1d-9, errrel=1d-9, maxfn=fne) 
                                        hrs = hhi
                                    else
                                        print *, 'WARNING: helper_func1 must have different signs!'
                                    end if                        
                                    labinc = w*eta(sc)*ep(tyc,jc)*hrs
                                    if (labinc < yb_cutoff) then
                                        cons = ( theta0*(1d0-theta1)*(w*eta(sc)*ep(tyc,jc))**(1d0-theta1)*hrs**(-theta1) - taup*w*eta(sc)*ep(tyc,jc) )/(1d0+tauc)/chi/(hrs**sig2)
                                    else
                                        cons = ( (1.0d0-tau_max)*(w*eta(sc)*ep(tyc,jc)) - taup*w*eta(sc)*ep(tyc,jc) )/(1d0+tauc)/chi/(hrs**sig2)
                                    end if
                                    cons = cons**(1d0/sig1)
                                    VV = U(cons,hrs) + bbeta*surv(jc)*sum(pi(sc,:)*Vfun(tyc,:,1,jc+1,jj))
                                else
                                    !aprime = LinInterp_1d(gridm(ia),mtmp,atmp,na)
                                    !hrs = LinInterp_1d(gridm(ia),mtmp,htmp,na)
                                    !cons = LinInterp_1d(gridm(ia),mtmp,ctmp,na)
                                    !VV = -1d0/LinInterp_1d(gridm(ia),mtmp,vtmp,na)
                                    !!VV = LinInterp_1d(gridm(ia),mtmp,vtmp,na)
                                    
                                    call basefun(mtmp,na,gridm(ia,jj),vals,inds)
                                    hrs = vals(1)*htmp(inds(1))+vals(2)*htmp(inds(2))
                                    cons = vals(1)*ctmp(inds(1))+vals(2)*ctmp(inds(2))
                                    VV = -1d0/(vals(1)*vtmp(inds(1))+vals(2)*vtmp(inds(2)))
                                    aprime = vals(1)*atmp(inds(1))+vals(2)*atmp(inds(2))                                
                                    
                                end if
                                afun(tyc,sc,ia,jc,jj) = aprime
                                lfun(tyc,sc,ia,jc,jj) = hrs
                                cfun(tyc,sc,ia,jc,jj) = cons
                                vpfun(tyc,sc,ia,jc,jj) = marginal_utility(cons,hrs)
                                Vfun(tyc,sc,ia,jc,jj) = VV
                                
                                labinc = w*eta(sc)*ep(tyc,jc)*hrs
                                test_bc(tyc,sc,ia,jc) = gridm(ia,jj) + theta0*min(labinc,yb_cutoff)**(1d0-theta1) + (1d0-tau_max)*max(0d0,labinc-yb_cutoff) - cons*(1d0+tauc) - taup*labinc
                                if (afun(tyc,sc,ia,jc,jj) <= acutoff(jj)) then 
                                    test_bc(tyc,sc,ia,jc) = test_bc(tyc,sc,ia,jc) - aprime*qd 
                                else
                                    test_bc(tyc,sc,ia,jc) = test_bc(tyc,sc,ia,jc) - aprime*qwd - psi_offshore 
                                end if  
                                
                                !if (afun(tyc,sc,ia,jc,jj) >= grida(na,jj)) then
                                !    print *, 'WARNING: grida too small?'
                                !end if                                
                                
                            end do
                        end if

                        !call toc()
                        if (save_res == .true.) then
                            open(21, file='gridm.txt')
                            open(22, file='afun.txt')
                            open(23, file='cfun.txt')
                            open(24, file='Vfun.txt')
                            open(25, file='dVfun.txt')
                            open(26, file='hfun.txt')
                            open(27, file='am_cutoff.txt')
                            open(28, file='test_bc.txt')
                            do ia=1,na
                                write(21,'(f16.6)') gridm(ia,jj)
                                write(22,'(f16.6)') afun(tyc,sc,ia,jc,jj)
                                write(23,'(f16.6)') cfun(tyc,sc,ia,jc,jj)
                                write(24,'(f16.6)') Vfun(tyc,sc,ia,jc,jj)
                                write(25,'(f16.6)') vpfun(tyc,sc,ia,jc,jj)
                                write(26,'(f16.6)') lfun(tyc,sc,ia,jc,jj)
                                write(28,'(f16.6)') test_bc(tyc,sc,ia,jc)
                            end do
                            write(27, '(2f16.6, i4)') acutoff(jj), mcutoff(jj), jc
                            close(21); close(22); close(23); close(24); close(25); close(26); close(27); close(28)    
                        end if
                    end do
                end do
                
                !max_test_bc = maxval(abs(test_bc))
                
                !print *, max_test_bc
            end do
        end do

    end subroutine SolveHH
    
    function foc_static(hrs)
        real(8) :: hrs
        real(8) :: foc_static
        real(8) :: fff_static
        external fff_static
        
        foc_static = fff_static(hrs)
    
    end function foc_static
    
    !function foc_static(hrs) 
    !use params, only: chi, sig2, ep, eta, jc, sc, tyc, w, theta0, theta1, taup, tau_max, yb_cutoff
    !    real(8) :: hrs
    !    real(8) :: wtot
    !    real(8) :: labinc
    !    !real(8) :: err
    !    real(8) :: foc_static
    !    
    !    wtot = w*eta(sc)*ep(tyc,jc)
    !    labinc = wtot*hrs
    !    if (labinc < yb_cutoff) then
    !        foc_static = chi*hrs**sig2 - lambda_mod*( theta0*(1d0-theta1)*(wtot)**(1d0-theta1)*hrs**(-theta1) - taup*wtot)
    !    else
    !        foc_static = chi*hrs**sig2 - lambda_mod*( (1-tau_max-taup)*wtot )
    !    end if
    !    !err = chi*hrs**sig2 - lambda_mod*(w*eta(sc)*ep(tyc,jc))
    !
    !end function foc_static    
    
    function helper_func1(hrs)
        real(8) :: hrs
        real(8) :: helper_func1
        real(8) :: fff_helper_func1
        external fff_helper_func1
        
        helper_func1 = fff_helper_func1(hrs)
    
    end function helper_func1    

end module Mod_Household
