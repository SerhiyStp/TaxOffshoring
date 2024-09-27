MODULE PARAMS

    ! THIS IS THE MODULE FOR PARAMETERS AND VARIABLES DEFINITION

    implicit none

    integer,parameter:: prec=selected_real_kind(15,307)

    integer, parameter :: n_moms = 6 !19
    integer, parameter :: klp_n_moms = 6
    integer :: nm_iter

    integer, parameter :: file_res_id = 52

    ! Offshoring
    real(8) :: Share_Offshoring, share_off_dist(5), wealth_obs_dist(6), wealth_tot_dist(6),income_dist(6), lab_income_dist(6)
    real(prec):: qd, qw, qwd
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
    integer, parameter :: ns=9 ! Persistent wage shocks - old
    integer, parameter :: nz=9 ! Persistent wage shocks
    integer, parameter :: na=1001 !201 !401 ! Assets
    integer, parameter :: nl=1 !
    integer, parameter :: nty=1
    integer, parameter :: maxit=10000
    integer, parameter :: nxi=3  ! Temporary wage shocks

    integer:: na0,na1,ntauk

    real(8) :: DistX(ns), DistXW(ns)


    ! Population
    integer,parameter:: jr = 46
    integer,parameter:: J = 81
    integer,parameter:: Tret = J-Jr+1
    integer,parameter:: Twork = Jr-1
    real(prec),parameter:: nn=0.0d0 !0.011


    ! Indicators
    integer:: ind_pref 
    ! = 1 non-separable 
    ! = 2 separable preference
    integer:: indext 
    ! = 1 for benchmark with Gouveia-Strauss tax function

    ! Parameters to be calibrated for each preference specs (to be calibrated)
    !real(prec),parameter:: beta_S	= 0.959d0 !0.98d0 !0.974 !0.9750 !0.97172 !1.2d0 
    real(prec) :: beta_S
    real(prec),parameter:: delta_S	= 0.0833d0  !0.045d0 !0.0780 !0.0833 

    real(prec),parameter:: beta_NS	= 1.00093 
    real(prec),parameter:: delta_NS	= 0.0833  

    ! (1) parameters for NON-separable preference
    real(prec),parameter:: gamma	= 0.377 
    real(prec),parameter:: sigma	= 4.0		

    ! (2) parameters for separable preference
    !real(prec),parameter:: chi	= 17.5d0 !3.70 !3.55 !3.20 !2.15 !2.05 !1.95 !1.92 ! to be calibrated  
    real(prec) :: chi
    real(prec),parameter:: sig1	= 1.509 !2.0  
    real(prec),parameter:: frisch=0.6d0
    real(prec),parameter:: sig2	= 1d0/frisch !2d0 !10d0/6d0 !3.0


    ! Production function
    real(prec),parameter:: alpha= 0.36
    real(prec),parameter:: TFP	= 1.0

    ! Government Policies
    real(prec),parameter:: govconsNS	= 6.670 ! Gov't consumption Non Separable utility
    real(prec),parameter:: govconsS		= 5.870 ! Gov't consumption Separable utility

    real(prec),parameter:: tauc=0.05 	!Consumption tax
    real(prec),parameter:: b=0.5		!Social Security Replacement Rate


    real(prec),parameter::umin=-1.0E+2
    real(prec),parameter::penscale=10000000
    real(prec),parameter::maxl=0.99

    ! Value of the borrowing constraint
    real(prec),parameter:: blimit=0.0

    ! Wage shocks and transition probabilities
    real(8) :: sig_z = 0.02d0
    real(8) :: rho_z = 0.8d0
    real(8) :: sig_xi = 0.01d0
    real(prec),dimension(nz,nz)::pi
    real(8) :: eta(nz)
    real(8) :: pi_xi(nxi)
    real(8) :: xi(nxi)
    real(prec)::avgeta
    real(prec),dimension(nz)::pistat, pini
    
    
    ! Retirement
    real(8) :: b_ret(nz)

    
    real(8)::pstat(1,ns)
    
    ! Heterogeneous returns
    real(8), parameter :: omega1 = 0.072d0
    real(8), parameter :: omega2 = 0.20d0
    real(8), parameter :: gamma_omega = 0.30d0
    real(8), parameter :: abar_omega = 0.0d0
    real(8), parameter :: omegabar = 0.4d0
    real(8), parameter :: rF = 0.03d0 !0.01d0
    real(8), parameter :: rR = 0.09d0 !0.06d0
    real(8), parameter :: sig_kappa = 0.05d0
    integer, parameter :: nkappa = 3 !1 !3 ! Rate of return temporary shocks
    real(8) :: Kappas(nkappa)
    real(8) :: pi_kappa(nkappa)
    integer, parameter :: ntheta=2 !1 !2 ! Rate of return persistent types    
    real(8) :: thetas(ntheta)  
    real(8) :: pi_theta(ntheta,ntheta)    
    real(8) :: pi_theta_stat(ntheta)
    

    ! Variables needed for the Tauchen routine
    integer :: NVAR,nlag=1
    integer :: nval(10)
    integer :: nsm
    real(prec) :: theta(100)		   ! 1:nvar = constant
    ! nvar*nvar = autoregressive matrix
    ! nvar*nvar = epsilons var cov matrix

    ! Parameters derived from Tauchen procedure (Markov)
    real(prec),allocatable :: mstae(:)			 ! Stationary distribution of e
    real(prec),allocatable :: mstates(:)		 ! Matrix of states e1,e2
    real(prec),allocatable :: mprobs(:,:)      ! Transition probabilities


    ! Points in grid of assets 
    real(prec),dimension(na)::grida


    ! Demographic Parameters
    real(prec),dimension(J):: ephansen,surv,mu,Nu
    real(prec),dimension(nty,J):: ep
    real(prec),dimension(nty):: measty
    real(prec):: bbeta,delta,topop,pop(81)


    ! Counters
    integer :: sc,ac,lc,jc,tc,date,tyc,agec


    ! Output
    real(prec),allocatable:: matresul(:,:,:,:)


    ! Tax code
    real(prec):: a2
    real(prec),allocatable:: a0(:)
    real(prec),allocatable:: a1(:)
    real(prec),allocatable:: tauk(:)

    real(8),parameter:: tk=0.283d0
    !real(8),parameter:: theta0=0.940d0 !0.917d0
    real(8) :: theta0
    real(8),parameter:: theta1=0.183d0 !0.137d0
    real(8),parameter:: tau_max = 0.396d0
    real(8) :: yb_cutoff
    
    real(8), parameter :: tau_estate = 0.1d0
    real(8), parameter :: a_estate = 0.5d0


    real(prec)::maxa0,maxa1,mina0,mina1
    real(prec)::maxtauk,mintauk

    integer:: a0c,a1c	
    integer:: taukc			

    !real(prec):: govcons
    real(8), parameter :: Govcons = 0.5d0 !25.5490651400000d0
    real(8) :: GovconsN

    real(prec)::socwelf2,optr,optw,optN,optK
    real(prec),dimension(5):: opttax
    real(prec),dimension(nty,ns,na,J,n_ofsh)::optvfun,optcfun,optlfun,optafun
    !real(prec),dimension(nty,ns,na,J,n_ofsh)::optafun


    !! Distribution over state space
    !real(prec),dimension(nty,ns,na,J,n_ofsh)::phi,phitot,optphi
    !real(8) :: Phi_work(nty,ns,na,Jr-1,n_ofsh)


    ! Household value and policy functions
    !real(prec),dimension(nty,ns,na,J,n_ofsh)::afun
    !real(prec),dimension(nty,ns,na,J,n_ofsh)::Vfun,cfun,lfun,vpfun,afun
    


    ! Asset distribution
    real(prec),dimension(na):: Adis


    ! Average variables by age
    !real(prec),dimension(J):: abar,lbar,labar,cbar,astartbar,meanearn,logmean,varlogearn


    ! Average variables by age and type
    real(prec),dimension(J,nty,n_ofsh):: abartype,astartbartype,lbartype,labartype
    real(prec),dimension(J,nty,n_ofsh)::meanearntype,logmeantype,varlogearntype,cbartype


    ! Prices of capital and labor; labor supply and capital stock, and other 
    real(prec):: r,w,N,LabS,K,As,Astart,Y,C,Tr,exdem,Totinctax,hours,Transagg,stdle,stdleini
    real(8) :: Rs, TaxS, RetS, Rs_aux, YauxS, AftTaxauxS, TaxCS, TaxaboveybS, TaxE
    real(8) :: rbar


    ! Social security taxes and benefits
    real(prec),parameter:: taup=0.124,maxSSrat=87000.0/37748.0

    real(prec):: maxSS
    real(prec):: SS,SSn,TotSStax


    ! Bequest 
    real(prec):: TrB,TrBn,Trstart


    ! Welfare measures
    real(prec),dimension(nty,ns,na,J):: equivar,equivarss
    real(prec):: equiss
    

    

    integer :: jj_glob


    !===========================================================================
    ! FUNCTIONS
    !===========================================================================

CONTAINS
    
    
    subroutine SetParams()
        !===========================================================================
        ! This subroutine sets the parameters for the model
        !===========================================================================
        call GRID
        call PREFERENCE
        call DEMOGRAPHICS
        call LABOR
        call RETURNS()
        
        ! Offshoring
        psi_vals = [0.1d0, 0.5d0, 1.5d0] 
        psi_prob = [1d0/3d0, 1d0/3d0, 1d0/3d0]
        
    end subroutine SetParams
    
    
    subroutine RETURNS()
        !use params
        real(8), parameter :: pi_theta_hl = 0.5d0
        real(8), parameter :: pi_theta_lh = 0.2d0
        real(8) :: p0(1, ntheta), p1(1, ntheta)
        real(8) :: dist
        integer :: i
        
        pi_theta(1,:) = [1.0d0 - pi_theta_hl, pi_theta_hl]
        pi_theta(2,:) = [pi_theta_lh, 1.0d0 - pi_theta_lh]
        
        p0(1,:) = 1d0/dble(ntheta) ![1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns)]

        dist = 1d0
        do i = 1, 500
            p1 = matmul(p0, pi_theta)
            dist = sum( (p1 - p0)**2d0 )
            p0 = p1
            if (dist < 1d-12) exit
        end do             
        pi_theta_stat = p0(1,:)
        !print *, sum(pi_theta_stat)
        
        !pi_theta = 1d0
        thetas = [1d0, 0.00d0] !1d0 ![1d0, 0d0]
        
        Kappas = [-sig_kappa, 0.0d0, +sig_kappa] ![0d0] 
        Kappas = exp(Kappas)
        pi_kappa = [0.25d0, 0.5d0, 0.25d0] !1d0 !
        rbar = 0d0
        
    end subroutine RETURNS   
    
    
    subroutine LABOR
        ! THIS SUBROUTINE DEFINES THE STOCHASTIC PROCESS FOR LABOR PRODUCTIVITY
        !use params
        use TAUCHEN_mod, only: tauchen_pareto

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
        call tauchen_pareto(sig_z, rho_z, nz, m_tauch, pareto_cutoff, alpha_pareto, eta, pi)
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
        real(prec),parameter::scale=150.0,curv=1.0 !2
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
        
        bbeta = 0.989d0!0.959d0
        chi = 17.4d0        

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

        ! open(unit=32,file='measpop.txt')
        ! rewind(32)
        ! write(32,fmt=*) mu
        ! rewind(32)
        ! close(32)

        topop=sum(Nu)

    end subroutine DEMOGRAPHICS !    
    
    
    function tax_income(y) result(res)
        real(8), intent(in) :: y
        real(8) :: res
        
        res = y - after_tax_income(y)
    end function tax_income
    
    
    function D_after_tax_income(y) result(res)
        real(8), intent(in) :: y
        real(8) :: res
        
        if (y >= yb_cutoff) then
            res = 1d0-tau_max    
        else
            res = theta0*(1d0-theta1)*y**(-theta1)
        end if
        
    end function D_after_tax_income
    
    
    elemental function after_tax_income(y) result(res)
        real(8), intent(in) :: y
        real(8) :: res
        
        res = theta0*min(yb_cutoff,y)**(1d0-theta1) + (1d0-tau_max)*max(0d0, y-yb_cutoff)
    
    end function after_tax_income
    
    
    function after_tax_income_aux(y) result(res)
        real(8), intent(in) :: y
        real(8) :: res
        
        if (y >= yb_cutoff) then
            res = yb_cutoff**(1d0-theta1)
            !res = theta0*yb_cutoff**(1d0-theta1)
        else
            res = y**(1d0-theta1)
            !res = theta0*y**(1d0-theta1)
        end if
    
    end function after_tax_income_aux
    
    function rfunc(a, theta, kappa)
        real(8), intent(in) :: a
        real(8), intent(in) :: theta
        real(8), intent(in) :: kappa
        real(8) :: om
        real(8) :: rfunc
        
        om = Omega(a,theta)
        rfunc = rbar + rF*(1.0d0-om) + rR*kappa*om
        
    end function rfunc
    
    function Da_rfunc(a, theta, kappa)
        real(8), intent(in) :: a
        real(8), intent(in) :: theta
        real(8), intent(in) :: kappa
        real(8) :: Da_om
        real(8) :: Da_rfunc
        
        Da_om = DaOmega(a, theta)
        Da_rfunc = (rR*kappa - rF)*Da_om
    end function Da_rfunc
    
    function Omega(a, theta)
        real(8), intent(in) :: a
        real(8), intent(in) :: theta
        real(8) :: Omega
        
        Omega = theta*(omegabar + min( omega1*( max(a - abar_omega, 0.0d0) )**gamma_omega, omega2 ) )
    end function Omega
    
    function DaOmega(a, theta)
        real(8), intent(in) :: a
        real(8), intent(in) :: theta
        real(8) :: DaOmega
        real(8) :: tmp
        
        tmp = omega1*( max(a - abar_omega, 0.0d0) )**gamma_omega
        if (a <= abar_omega .or. tmp >= omega2) then
            DaOmega = 0d0
        else
            DaOmega = theta*omega1*gamma_omega*(a - abar_omega)**(gamma_omega-1.0d0)
        end if
        
    end function DaOmega

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


    !==========================================================================
    function Tax2(b0,b1,b2,ear)
        !==========================================================================

        ! Tax function in Gouveia & Strauss (National Tax Journal, 1994)

        implicit none

        real(prec):: Tax2
        real(prec),intent(in):: b0,b1,b2,ear

        Tax2=0.0

        if ( b1 >= 0.001 ) then
            if ( ear > 0.0 ) then
                Tax2=b0*( ear - ( ear**(-b1)+b2 )**(-1.0/b1) )
            endif
        endif

        if (b1 < 0.001 ) then
            if ( ear > 0.0 ) then
                Tax2=b0*ear +b2
            endif
        endif

    end function Tax2
    !==========================================================================


    !==========================================================================
    function Tax3(l0,l1,l2,ca0,lear,caear)
        !==========================================================================

        real(prec):: Tax3
        real(prec),intent(in):: l0,l1,l2,ca0,lear,caear


        Tax3=0.0

        ! Labor Income Tax 

        if ( l1 >= 0.000001 ) then
            if ( lear >= 0.00000001 ) then
                Tax3 = l0*( lear - ( lear**(-l1)+l2 )**(-1.0/l1) )
            endif
        endif
        if (l1 < 0.000001 ) then
            if ( lear >= 0.00000001 ) then
                Tax3 = l0*lear+l2
            endif
        endif


        ! Capital Income Tax

        Tax3 = Tax3 + ca0*caear

    end function Tax3
    !==========================================================================


    !==========================================================================
    function Tax4(earncap,lear,jj)
        !==========================================================================
        real(8):: Tax4
        real(8):: ap, lear, earncap !, caear
        integer:: jj

        Tax4 = lear - theta0*min(yb_cutoff,lear)**(1.0d0-theta1) - (1.0d0-tau_max)*max(0d0,lear-yb_cutoff)
        Tax4 = Tax4 + tk*earncap

    end function Tax4
    !==========================================================================

    !==========================================================================
    function MarTax2(b0,b1,b2,capear)
        !==========================================================================

        ! Marginal Taxes for capital
        ! Tax function in Gouveia & Strauss (National Tax Journal, 1994)

        implicit none

        real(prec):: MarTax2
        real(prec),intent(in):: b0,b1,b2,capear

        MarTax2=0.0

        if ( b1 >= 0.001 ) then
            if ( capear > 0.0 ) then
                MarTax2=b0*( 1.0 - (1.0+b2*capear**b1)**(-1.0/b1 - 1.0) )
            endif
        endif
        if (b1 < 0.001 ) then
            if ( capear > 0.0 ) then
                MarTax2=b0
            endif
        endif

    end function MarTax2
    !==========================================================================


    !==========================================================================
    function MarTax3(l0,l1,l2,ca0,lear,caear)
        !==========================================================================

        real(prec):: MarTax3
        real(prec),intent(in):: l0,l1,l2,ca0,lear,caear

        ! Marginal Capital Income Tax
        MarTax3=ca0


    end function MarTax3
    !==========================================================================



    !==========================================================================
    subroutine basefun (grid_x,npx,x,vals,inds) 
        !==========================================================================

        implicit none
        ! this subroutine returns the values and the indices of the two basis
        ! functions that are positive on a given x in the grid_x

        real(prec),intent(in) :: x
        integer , intent(in):: npx
        real(prec), intent(in) :: grid_x (npx)
        real(prec), intent(out) ::vals(2)
        integer ,intent(out) ::inds(2)
        integer :: i,ju,jl,jm

        jl=1     
        ju=npx   	

        do

            if (ju-jl<=1)  exit
            jm=(ju+jl)/2
            if (x>=grid_x(jm)) then
                jl=jm
            else
                ju=jm
            endif

        end do

        i=jl+1
        vals(2)=( x-grid_x(i-1) )/(grid_x(i)-grid_x(i-1))
        vals(1)=( grid_x(i)-x )/(grid_x(i)-grid_x(i-1))
        inds(2)=i
        inds(1)=i-1

    end subroutine basefun
    !==========================================================================

ENDMODULE PARAMS
