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
