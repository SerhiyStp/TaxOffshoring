MODULE PARAMS

    ! THIS IS THE MODULE FOR PARAMETERS AND VARIABLES DEFINITION

    implicit none

    integer,parameter:: prec=selected_real_kind(15,307)

    integer :: nm_iter

    integer, parameter :: file_res_id = 52

    ! Offshoring
    real(prec),parameter:: psi_offshore=2.5d0 !1.5d0 !0.5d0
    integer,parameter:: n_ofsh=3 !1 !3
    !real(prec),parameter:: psi_vals(n_ofsh)=[psi_offshore-0.5d0, psi_offshore, psi_offshore+3.5d0]
    real(prec) :: psi_vals(n_ofsh)
    real(8) :: psi_val_top, psi_val_mid
    real(8) :: p_out_h
    real(prec) :: psi_prob(n_ofsh) !=[1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0]
    real(prec),parameter:: frac_ofsh=0.4d0

    real(8) :: a_cutoff(n_ofsh)


    ! Grids
    !integer,parameter:: ns=7, na=501, nl=66, nty=2, maxit=10000  ! Size of grids
    !integer,parameter:: ns=8, na=501, nl=66, nty=1, maxit=10000  ! Size of grids
    !integer,parameter:: ns=8, na=3001, nl=1, nty=1, maxit=10000  ! Size of grids
    !integer, parameter :: ns=9 ! Persistent wage shocks - old
    integer, parameter :: nz=9 ! Persistent wage shocks
    integer, parameter :: na=1001 !201 !401 ! Assets
    integer, parameter :: nl=1 !
    integer, parameter :: nty=1
    integer, parameter :: maxit=10000
    integer, parameter :: nxi=3  ! Temporary wage shocks

    integer:: na0,na1,ntauk

    !real(8) :: DistX(ns), DistXW(ns)


    ! Population
    integer,parameter:: Jr = 46
    integer,parameter:: J = 81
    integer,parameter:: Tret = J-Jr+1
    integer,parameter:: Twork = Jr-1
    real(prec),parameter:: nn=0.0d0 !0.011


    ! Indicators
    real(prec),parameter:: delta_S	= 0.0833d0  !0.045d0 !0.0780 !0.0833 

    ! (2) parameters for separable preference
    !real(prec),parameter:: chi	= 17.5d0 !3.70 !3.55 !3.20 !2.15 !2.05 !1.95 !1.92 ! to be calibrated  
    real(prec) :: chi
    real(prec),parameter:: sig1	= 2.0d0 !1.509 !2.0  
    real(prec),parameter:: frisch=0.6d0
    real(prec),parameter:: sig2	= 1d0/frisch !2d0 !10d0/6d0 !3.0


    ! Production function
    real(prec),parameter:: alpha= 0.36
    real(prec),parameter:: TFP	= 1.0

    ! Value of the borrowing constraint
    real(prec),parameter:: blimit=0.0

    ! Wage shocks and transition probabilities
    real(8) :: sig_z = 0.50d0 !0.58d0 !0.05d0 !0.02d0
    real(8) :: rho_z = 0.9d0 !0.973d0 !0.95d0 !0.8d0
    real(8) :: sig_xi = 0.25d0
    real(prec),dimension(nz,nz)::pi
    real(8) :: eta(nz)
    real(8) :: pi_xi(nxi)
    real(8) :: xi(nxi)
    real(prec)::avgeta
    real(prec),dimension(nz)::pistat, pini
    
    
    ! Retirement
    real(8) :: b_ret(nz)
    
    ! Points in grid of assets 
    real(prec),dimension(na)::grida


    ! Demographic Parameters
    real(prec),dimension(J):: ephansen,surv,mu,Nu
    real(prec),dimension(nty,J):: ep
    real(prec),dimension(nty):: measty
    real(prec):: bbeta,delta,topop,pop(81)


    ! Counters
    integer :: sc,ac,lc,jc,tc,date,tyc,agec

    ! Asset distribution
    real(prec),dimension(na):: Adis

    ! Average variables by age and type
    real(prec),dimension(J,nty,n_ofsh):: abartype,astartbartype,lbartype,labartype
    real(prec),dimension(J,nty,n_ofsh)::meanearntype,logmeantype,varlogearntype,cbartype


    ! Prices of capital and labor; labor supply and capital stock, and other 
    real(prec):: r, w, N, LabS, K, As, Astart, Y, C, Tr, exdem, Totinctax, hours, Transagg, stdle, stdleini
    real(8) :: YauxS, AftTaxauxS, TaxCS, TaxE
    real(8) :: AAgg, RetAgg, HrsAgg
    real(8) :: RAgg, RauxAgg, LAgg, CAgg
    real(8) :: TaxInc, TaxC, TaxTot
    real(8) :: TaxIncAboveYb, YfBelowYb, DBelowYb
    real(8) :: TotOffshCost

    ! Bequest 
    real(prec):: TrB,TrBn 

    !===========================================================================
    ! FUNCTIONS
    !===========================================================================

CONTAINS
    
    
    subroutine SetParams()
        !===========================================================================
        ! This subroutine sets the parameters for the model
        !===========================================================================    
        use Het_returns, only: set_returns
        
        implicit none

        call GRID
        call PREFERENCE
        call DEMOGRAPHICS
        call LABOR
        call SET_RETURNS()
        
        ! Offshoring
        psi_vals = [0.1d0, 0.5d0, 1.5d0] 
        psi_prob = [1d0/3d0, 1d0/3d0, 1d0/3d0]
        
    end subroutine SetParams
        
    
    subroutine LABOR
        ! THIS SUBROUTINE DEFINES THE STOCHASTIC PROCESS FOR LABOR PRODUCTIVITY
        !use params
        use TAUCHEN_mod, only: discretize_w_pareto

        implicit none
        
        ! endowment process parameters
        integer::i
        real(8)::p0(1,nz), p1(1,nz) 
        real(8)::dist
        
        real(8) :: m_tauch, pareto_cutoff
        real(8) :: alpha_pareto        


        !yb_cutoff = ( theta0*(1.0d0-theta1)/(1.0d0-tau_max) )**(1d0/theta1)

        if (nz==1) then
            eta=1.0d0
            pi=1.0d0
            pini=1.0d0
        end if
        
        pareto_cutoff = 0.9d0
        m_tauch = 2.7d0  
        alpha_pareto = 1.9d0
        call discretize_w_pareto(sig_z, rho_z, nz, m_tauch, pareto_cutoff, alpha_pareto, eta, pi)
        !call read_mc(eta, pi)
        p0(1,:) = 1d0/dble(nz) ![1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns)]

        dist = 1d0
        do i = 1, 500
            p1 = matmul(p0, pi)
            dist = sum( (p1 - p0)**2d0 )
            p0 = p1
            if (dist < 1d-12) exit
        end do             
        pini = p0(1,:)
        
        xi = exp([-sig_xi, 0.0d0, sig_xi]) ! exp(0.0d0) !
        pi_xi = [0.25d0, 0.5d0, 0.25d0] !1d0 !

        ! Determine the initial Distribution of Labor Productivities

        ! A: Fixed Effects: sigma**2_alpha=0.247
        !ep(1,1:J)=exp(-sqrt(sigma2alpha))*ephansen(1:J)
        !ep(2,1:J)=exp(sqrt(sigma2alpha))*ephansen(1:J)
        
        ep(1,1:J)=ephansen(1:J)
        
        b_ret = 0.4d0*eta

    end subroutine LABOR     
    
    subroutine GRID
        ! THIS SUBROUTINE DEFINES GRID FOR INDIVIDUAL ASSET HOLDING
        !use params
        implicit none

        !real(prec),parameter::scale=75.0,curv=1.1 !2.5
        !real(prec),parameter::scale=105.0,curv=1.1 !2.5
        !real(prec),parameter::scale=1800.0,curv=1.1 !2.5
        !real(prec),parameter::scale=2000.0,curv=1.1 !2.5
        !real(prec),parameter::scale=7500.0,curv=1.1 !2.5
        !real(prec),parameter::scale=505000.0,curv=1.1 !2.5
        !real(prec),parameter::scale=400000.0,curv=1.1 !2.5
        !real(prec),parameter::scale=6000.0,curv=1.1 !2.5
        !real(prec),parameter::scale=40000.0,curv=1.1 !2.5
        !real(prec),parameter::scale=2000.0,curv=1.2
        real(prec),parameter::scale=15000.0d0 !800.0d0 !500.0d0 !150.0
        real(prec),parameter::curv=1.0 !2
        real(prec)::step

        grida(1)=blimit
        do ac=2,na
            grida(ac)=grida(1)+scale*((ac-1.0)/(na-1.0))**curv
        end do

    end subroutine GRID 
    
    subroutine PREFERENCE
        ! THIS SUBROUTINE DEFINES THE PREFERENCE SPECIFICATION
        !use params

        implicit none

        !bbeta	= beta_S
        delta	= delta_S
        
        bbeta = 0.945d0 !0.980d0 !0.989d0!0.959d0
        chi = 500.0d0 !250.0d0 !45.0d0 !25.0d0 !20.0d0 !17.4d0        

    end subroutine PREFERENCE     
    
    subroutine DEMOGRAPHICS
        ! THIS SUBROUTINE BUILDS GRIDS FOR LABOR EFFICIENCY UNITS, SURVIVAL RATES
        ! AND AGE-DISTRIBUTION
        !use params
        implicit none

        integer::i

        ! Age-Efficiency Units from Hansen (1993)
        ephansen(1)=1.0000
        ephansen(2)=1.0719
        ephansen(3)=1.1438
        ephansen(4)=1.2158
        ephansen(5)=1.2842
        ephansen(6)=1.3527
        ephansen(7)=1.4212
        ephansen(8)=1.4897
        ephansen(9)=1.5582
        ephansen(10)=1.6267
        ephansen(11)=1.6952
        ephansen(12)=1.7217
        ephansen(13)=1.7438
        ephansen(14)=1.7748
        ephansen(15)=1.8014
        ephansen(16)=1.8279
        ephansen(17)=1.8545
        ephansen(18)=1.8810
        ephansen(19)=1.9075
        ephansen(20)=1.9341
        ephansen(21)=1.9606
        ephansen(22)=1.9623
        ephansen(23)=1.9640
        ephansen(24)=1.9658
        ephansen(25)=1.9675
        ephansen(26)=1.9692
        ephansen(27)=1.9709
        ephansen(28)=1.9726
        ephansen(29)=1.9743
        ephansen(30)=1.9760
        ephansen(31)=1.9777
        ephansen(32)=1.9700
        ephansen(33)=1.9623
        ephansen(34)=1.9546
        ephansen(35)=1.9469
        ephansen(36)=1.9392
        ephansen(37)=1.9315
        ephansen(38)=1.9238
        ephansen(39)=1.9161
        ephansen(40)=1.9084
        ephansen(41)=1.9007
        ephansen(42)=1.8354
        ephansen(43)=1.7701
        ephansen(44)=1.7048 
        ephansen(45)=1.6396

        do i=jr,J
            ephansen(i)=0.0
        end do         


        do tyc=1,nty
            ep(tyc,1:J)=ephansen(1:J)
        end do

        measty(1:nty)=1.0/nty

        ! Population Numbers from Bell and Miller (2002)

        pop(1)=	197316
        pop(2)=	197141
        pop(3)=	196959
        pop(4)=	196770
        pop(5)=	196580
        pop(6)=	196392
        pop(7)=	196205
        pop(8)=	196019
        pop(9)=	195830
        pop(10)=195634
        pop(11)=195429
        pop(12)=195211
        pop(13)=194982
        pop(14)=194739
        pop(15)=194482
        pop(16)=194211
        pop(17)=193924
        pop(18)=193619
        pop(19)=193294
        pop(20)=192945
        pop(21)=192571
        pop(22)=192169
        pop(23)=191736
        pop(24)=191271
        pop(25)=190774
        pop(26)=190243
        pop(27)=189673
        pop(28)=189060
        pop(29)=188402
        pop(30)=187699
        pop(31)=186944
        pop(32)=186133
        pop(33)=185258
        pop(34)=184313
        pop(35)=183290
        pop(36)=182181
        pop(37)=180976
        pop(38)=179665
        pop(39)=178238
        pop(40)=176689
        pop(41)=175009
        pop(42)=173187
        pop(43)=171214
        pop(44)=169064
        pop(45)=166714
        pop(46)=164147
        pop(47)=161343
        pop(48)=158304
        pop(49)=155048
        pop(50)=151604
        pop(51)=147990
        pop(52)=144189
        pop(53)=140180
        pop(54)=135960
        pop(55)=131532
        pop(56)=126888
        pop(57)=122012
        pop(58)=116888
        pop(59)=111506
        pop(60)=105861
        pop(61)=99957
        pop(62)=93806
        pop(63)=87434
        pop(64)=80882
        pop(65)=74204
        pop(66)=67462
        pop(67)=60721
        pop(68)=54053
        pop(69)=47533
        pop(70)=41241
        pop(71)=35259
        pop(72)=29663
        pop(73)=24522
        pop(74)=19890
        pop(75)=15805
        pop(76)=12284
        pop(77)=9331
        pop(78)=6924
        pop(79)=5016
        pop(80)=3550
        pop(81)=2454


        ! Survival probabilities: surv(i)=prob(alive in i+1|alive in i)

        do i = 1,J-1
            surv(i) = pop(i+1)/pop(i)
        end do

        surv(J) = 0.0

        ! Number of Agents in population

        Nu(1) = 1.0

        do i = 2,J
            Nu(i) = surv(i-1)*Nu(i-1)/(1.0+nn)	  
        end do

        ! Fraction of agents in population

        do i = 1,J
            mu(i) = Nu(i)/sum(Nu)
        end do
        
        print *, 'Sum(mu) = ', sum(mu)

        ! open(unit=32,file='measpop.txt')
        ! rewind(32)
        ! write(32,fmt=*) mu
        ! rewind(32)
        ! close(32)

        topop=sum(Nu)

    end subroutine DEMOGRAPHICS !    
    
    !==========================================================================
    function U(c,l)
        !==========================================================================

        implicit none
        real(prec):: U,c,l

        U = ( c**(1.0d0-sig1) - 1.0d0 )/(1.0d0-sig1) - chi * ( l**(1.0d0+sig2) )/(1.d0+sig2) !- 10d0

        !IF (ind_pref.eq.1) THEN 
        !
        !	if (c<=0.0) then
        !		U=umin-penscale*c**2.0
        !	else if (l<0.0) then
        !		U=umin-penscale*l**2.0
        !	else if (l>=1.0) then
        !		U=umin-penscale*(l-1.0)**2.0
        !	else 
        !		U= (1.0/(1.0-sigma))*( ( (c**gamma)*((1.0-l)**(1.0-gamma)) )**(1.0-sigma) )
        !	end if
        !
        !ELSE
        !
        !	if (c<=0.0) then
        !		U=umin-penscale*c**2.0
        !	else if (l<0.0) then
        !		U=umin-penscale*l**2.0
        !	else if (l>=1.0) then
        !		U=umin-penscale*(l-1.0)**2.0
        !	else 
        !		U=( c**(1.0-sig1) )/(1.0-sig1) + chi * ( (1.0-l)**(1.0-sig2) )/(1.0-sig2) 
        !	end if
        !
        !ENDIF

    end function U
    !==========================================================================


    !==========================================================================
    function marginal_utility(c, l)
        !==========================================================================

        real (prec), intent (in) :: c, l
        real (prec) :: marginal_utility

        !IF (ind_pref.eq.1) THEN 
        !
        !	if (c>0.0) then
        !		marginal_utility = ( gamma*( c**gamma*(1.0-l)**(1.0-gamma) )**(1.0-sigma) )/c
        !	else
        !		marginal_utility = 1000000.0+ abs(c)**2.0
        !	endif
        !
        !ELSE
        !
        if (c>0.0d0) then
            marginal_utility = c**(-sig1)
        else
            marginal_utility = 1000000.0+ abs(c)**2.0d0
        endif
        !
        !ENDIF

    end function marginal_utility
    !==========================================================================

    
    subroutine read_mc(zvals, prob)
        !use PARAMS, only: nz
        integer :: i, iu
        real(8) :: prob(nz,nz)
        real(8) :: zvals(nz), logzvals(nz)
        real(8) :: tmp(nz)
        real(8) :: s_z, m_z
        
        open(newunit=iu, file='MC_b20.txt')
        do i = 1, nz
            read(iu, '(<nz+2>f26.16)'), tmp(i), logzvals(i), prob(i,:)   
            print *, sum(prob(i, :))
        end do
        close(iu)
        s_z = tmp(1)
        m_z = tmp(3)   
        zvals = exp(logzvals)
    end subroutine read_mc    
    
    
END MODULE PARAMS
