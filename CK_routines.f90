module CK_routines
    
    implicit none
    real(8):: newton_res(3)
    
contains
    
    subroutine GRID

        ! THIS SUBROUTINE DEFINES GRID FOR INDIVIDUAL ASSET HOLDING

        use params
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

        use params

        implicit none

        ! SELECT UTILITY FUNCTION 
        ind_pref = 2
        ! = 1 NON-SEPARABLE
        ! = 2 SEPARABLE

        ! SELECT TAX FUNCTION 
        indext = 2
        ! = 1 BENCHMARK: Gouveia-Strauss TAX FUNCTION
        ! = 2 SEPARABLE IN CAPITAL AND LABOR INCOME (defined in Gridtax.f90/Tax3/Martax3)


        IF (ind_pref.eq.1) THEN 
            bbeta	= beta_NS
            delta	= delta_NS
            govcons = govconsNS
        ELSE
            bbeta	= beta_S
            delta	= delta_S
            govcons = govconsS
        ENDIF

    end subroutine PREFERENCE   
    
    subroutine DEMOGRAPHICS

        ! THIS SUBROUTINE BUILDS GRIDS FOR LABOR EFFICIENCY UNITS, SURVIVAL RATES
        ! AND AGE-DISTRIBUTION

        use params
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

        open(unit=32,file='measpop.txt')
        rewind(32)
        write(32,fmt=*) mu
        rewind(32)
        close(32)

        topop=sum(Nu)

    end subroutine DEMOGRAPHICS    
    
    subroutine LABOR

        ! THIS SUBROUTINE DEFINES THE STOCHASTIC PROCESS FOR LABOR PRODUCTIVITY

        use params
        ! use toolbox, only: rouwenhorst
        ! use toolbox_kk

        !use MSIMSL
        use TAUCHEN_mod, only: tauchen_pareto

        implicit none


        real(prec)::pilab(ns,jr-1)
        real(prec)::stdlel(jr-1)
        real(prec)::iden(ns,ns)

        ! endowment process parameters

        real(prec):: rho
        real(prec):: sigmaeps
        real(prec):: sigma2alpha

        real(8)::pi_small(5,5)
        real(8)::pi_small_test(5,5)
        real(8)::eta_small(5)
        real(8)::eta_small_test(5)
        real(8)::pistat_small(5)
        real(8)::p_i6
        real(8)::p_66
        real(8)::p_67
        real(8)::p_77
        real(8)::e6,e7
        real(8)::pi_test(ns,ns)
        real(8)::eta_test(ns)
        integer::i
        real(8)::tmp, sigmaeps2_test
        real(8)::p0(1,ns), ptmp(ns,ns)
        real(8)::dist
        
        
        real(8) :: m_tauch, pareto_cutoff
        real(8) :: alpha_pareto        


        yb_cutoff = ( theta0*(1.0d0-theta1)/(1.0d0-tau_max) )**(1d0/theta1)

        !if (ns==1) then
        eta=1.0d0
        pi=1.0d0
        pini=1.0d0
        !go to 200
        !else
        !end if


        nvar=1
        nval(1)=ns
        nsm=ns

        rho = 0.9850d0
        !sigma2alpha = 0.2061d0*0.59d0 + 0.1517d0*0.41d0
        sigma2alpha = 0d0
        tmp = 0.0346d0*0.59d0 + 0.0180d0*0.41d0
        sigmaeps = sqrt(tmp)/sqrt(1.0d0-rho**2.0d0)

        allocate(mstae(ns**nvar))	
        allocate(mstates(ns**nvar))
        allocate(mprobs(ns**nvar,ns**nvar)) 

        theta(1)=0.0
        theta(2)=rho
        theta(3)=(sigmaeps*(1.0-rho**2)**(1/2))**2

        !call markov
        !pistat=mstae(1:ns)
        !eta=mstates(1:ns)
        !pi=mprobs(1:ns,1:ns)
        !
        !!print *, sum(pistat)
        !!print *, sum(pi(1,:))
        !
        !
        !eta=exp(eta)
        !eta=eta/(sum(pistat*eta))    ! normalization of e
        !						     ! e are ordered e1<e2...<en

        ! Kindermann & Krueger parameters for unskilled: =======================================
        !pi(1,:) = [0.969909d0, 0.029317d0, 0.000332d0, 0.00002d0, 0.0d0, 0.000440d0, 0d0]
        !pi(1,:) = pi(1,:)/sum(pi(1,:))
        !pi(2,:) = [0.007329d0, 0.970075d0, 0.021989d0, 0.000166d0, 0.0d0, 0.000440d0, 0d0]
        !pi(2,:) = pi(2,:)/sum(pi(2,:))
        !pi(3,:) = [0.000055d0, 0.014659d0, 0.970130d0, 0.014659d0, 0.000055d0, 0.000440d0, 0d0]
        !pi(3,:) = pi(3,:)/sum(pi(3,:))
        !pi(4,:) = [0d0, 0.000166d0, 0.021989d0, 0.970075d0, 0.007329d0, 0.000440d0, 0d0]
        !pi(4,:) = pi(4,:)/sum(pi(4,:))
        !pi(5,:) = [0d0, 0.000002d0, 0.000332d0, 0.029317d0, 0.989909d0, 0.00044d0, 0d0]
        !pi(5,:) = pi(5,:)/sum(pi(5,:))
        !pi(6,:) = [0d0, 0d0, 0.002266d0, 0d0, 0d0, 0.97d0, 0.027734d0]
        !pi(6,:) = pi(6,:)/sum(pi(6,:))
        !pi(7,:) = [0d0, 0d0, 0d0, 0d0, 0d0, 0.288746d0, 0.711254d0]
        !pi(7,:) = pi(7,:)/sum(pi(7,:))
        !!print *, sum(pi(1,:))
        !!print *, sum(pi(2,:))
        !!print *, sum(pi(3,:))
        !!print *, sum(pi(4,:))
        !!print *, sum(pi(5,:))
        !!print *, sum(pi(6,:))
        !!print *, sum(pi(7,:))
        !
        !!print *, 'test'
        !eta = [0.1354d0, 0.3680d0, 1.0d0, 2.7176d0, 7.3853d0, 19.7204d0, 654.0124d0]
        
        !pi(1,:) = [0.9729d0,    0.0265d0,     0.0002d0,        0.0d0,         0.0d0,         0.0d0,       0.0004d0,      0.0d0]
        !pi(2,:) = [0.0134d0,    0.9729d0,     0.0134d0,        0.0d0,         0.0d0,         0.0d0,       0.0004d0,      0.0d0]
        !pi(3,:) = [0.0002d0,    0.0265d0,     0.9729d0,        0.0d0,         0.0d0,         0.0d0,       0.0004d0,      0.0d0]
        !pi(4,:) = [0.0d0,       0.0d0,        0.0d0,           0.9729d0,      0.0265d0,      0.0002d0,    0.0004d0,      0.0d0]
        !pi(5,:) = [0.0d0,       0.0d0,        0.0d0,           0.0134d0,      0.9729d0,      0.0134d0,    0.0004d0,      0.0d0]
        !pi(6,:) = [0.0d0,       0.0d0,        0.0d0,           0.0002d0,      0.0265d0,      0.9729d0,    0.0004d0,      0.0d0]
        !pi(7,:) = [0.0042d0,    0.0042d0,     0.0042d0,        0.0042d0,      0.0042d0,      0.0042d0,    0.9690d0,      0.0057d0]
        !pi(8,:) = [0.0d0,       0.0d0,        0.0d0,           0.0d0,         0.0d0,         0.0d0,       0.0576d0,      0.9424d0]

          
        !do i = 1, ns
        !    pi(i,:) = pi(i,:)/sum(pi(i,:))
        !    !print *, sum(pi(i,:))
        !end do
        
        pareto_cutoff = 0.9d0
        m_tauch = 2.7d0  
        alpha_pareto = 1.9d0
        call tauchen_pareto(sig_z, rho_z, nz, m_tauch, pareto_cutoff, alpha_pareto, eta, pi)
        xi = exp(0.0d0) !exp([-sig_xi, 0.0d0, sig_xi])
        pi_xi = 1d0 ![0.25d0, 0.5d0, 0.25d0]

        ! Determine the initial Distribution of Labor Productivities

        ! A: Fixed Effects: sigma**2_alpha=0.247

        ep(1,1:J)=exp(-sqrt(sigma2alpha))*ephansen(1:J)
        !ep(2,1:J)=exp(sqrt(sigma2alpha))*ephansen(1:J)



        ! B: Initialization of the Persistent/Transitory Component
        ! From STY (2000), Figure 2
        ! The variance of log income at 25 is 0.2735 and at age 60 is about 0.9
        ! From the permanent component we get a variance of 0.242
        ! Use the variable "weight" to make the process as persistent as the consinuous time
        ! process (since Tauchen does not get the presistence OK). Since variance increase almost
        ! linearly with age, need persitence close to 1
        ! also Tauchen does not get the variability correct, so use adj to get the process sufficiently
        ! dispersed
        ! Look at Figure 2 of STY (2000) and at the variable stdlel below for the variance of labor productivity
        ! fed into the model; it provides a "reasonable" approximation. Of course the problem remains that we
        ! assume that all of the varaince in the data is due to productivity, rather than hours; this may
        ! indicate that our productivity process now is too dispersed


        !pini(1)=0.0
        !pini(2)=0.0
        !pini(6)=0.0
        !pini(7)=0.0
        !pini(3)=0.0
        !pini(5)=0.0
        !pini(4)=1.0
        !pini = [0.044d0, 0.412d0, 0.044d0, 0.044d0, 0.412d0, 0.044d0, 0d0, 0d0]

        !pini = pstat(1,:)
        !pini(7:8) = 0d0
        !pini = pini/sum(pini)
        !
        !!print *, sum(pini)
        !open(21,file='stat_pi.txt')
        !do i = 1, ns
        !    write(21,'(f11.6)') pini(i)    
        !end do
        !

        ! Do Markov Mixing with Identity Matrix to increase persistence

        iden=0.0
        do sc=1,ns
            iden(sc,sc)=1.0
        end do


        stdle=(sum( pistat* (log(eta) - sum(pistat*log(eta)))**2 ))
        stdleini=(sum( pini* (log(eta) - sum(pini*log(eta)))**2 ))


        ! Compute cross-sectional std dev of labor productivity over time

        pilab(1:ns,1)=pini
        do tc=2,jr-1
            do sc=1,ns
                pilab(sc,tc)=sum(pi(1:ns,sc)*pilab(1:ns,tc-1))
            end do
        end do

        do tc=1,jr-1
            stdlel(tc)=(sum( pilab(1:ns,tc)* (log(eta) - sum(pistat*log(eta)))**2 ))
        end do
        
        b_ret = 0.4d0*eta

    end subroutine LABOR 
    
    subroutine returns_heterogeneity()
        use params
        real(8), parameter :: pi_theta_hl = 0.5d0
        real(8), parameter :: pi_theta_lh = 0.2d0
        
        !pi_theta(1,:) = [1.0d0 - pi_theta_hl, pi_theta_hl]
        !pi_theta(2,:) = [pi_theta_lh, 1.0d0 - pi_theta_lh]
        pi_theta = 1d0
        thetas = 1d0 ![1d0, 0d0]
        
        Kappas = [0d0] ![-sig_kappa, 0.0d0, +sig_kappa]
        Kappas = exp(Kappas)
        pi_kappa = 1d0 ![0.25d0, 0.5d0, 0.25d0]
        rbar = 0d0
        
    end subroutine returns_heterogeneity
    
    subroutine initialize()
        use moments

        call GRID
        call PREFERENCE
        call DEMOGRAPHICS
        call LABOR
        call returns_heterogeneity()
        call set_moments()
        
    end subroutine initialize    
    
    subroutine resid(x1,x2,x3,x4,x5,fv1,fv2,fv3,fv4,fv5)
        use params
        use Mod_Household
        !use Mod_Distribution
        use int_tictoc

        implicit none

        real(prec),intent(in):: x1,x2,x3,x4,x5
        real(prec),intent(out):: fv1,fv2,fv3,fv4,fv5
        
        r   = x1
        N   = x2
        !a2  = x3
        Govcons = x3
        TrB = x4 
        SS  = x5 

        K = N*( (alpha*TFP) / (r+delta) )**(1.0/(1.0-alpha))        ! Capital Stock
        Y = TFP*(K**alpha)*(N**(1.0-alpha))					        ! Aggregate Output
        w = (1.0-alpha)*Y/N											! Wages


        maxSS = maxSSrat*Y/sum(Nu)  
        
        ! Solve the Household Problem

        !CALL HOUSEHOLD

        !qd = 1.0d0/(1.0d0 + r*(1.d0-tk))
        !qw = 1.0d0/(1.0d0 + r)
        qd = 1.0d0/(1.0d0 + r*(1.d0-tk))
        qw = 1.0d0/(1.0d0 + r)
        qwd = qw*frac_ofsh + qd*(1d0-frac_ofsh)

        call tic()
        call SolveHH(save_res=.true.)
        print *, 'Solution took: '
        call toc()


        ! Find stationary Distribution
        !call tic()
        !CALL DISTRIBUTION
        !call toc()
        call tic()
        !call Distribution_alt2(save_res=.false.)
        print *, 'Simulation took: '
        call toc()


        ! Compute residuals of functions we want to set to zero
        fv1=As-K*(1.0+nn)
        fv2=LabS-N
        fv3=Govcons-tauc*C-Totinctax
        fv4=TrB-TrBn
        fv5=SS-SSn

        GovconsN = tauc*C + Totinctax

        ! Goods Market Clearing
        Y	 = TFP * (K**alpha) * (LabS**(1.0-alpha))	
        exdem=( C + As - (1.0-delta)*K + Govcons - Y )
           
    end subroutine resid
    
    subroutine newton(fun,gues1,gues2,gues3,gues4,gues5)

        ! This subroutine computes the steady state interest rate and labor supply and intercept for tax system using the 
        ! the classical newton method; inputs are the guesses for r, N and Tint and the subroutine calls resid that
        ! delivers the residual from markets clearing in the asset market, labor market and the gvernment budget constraint

        use params

        implicit none

        integer::i
        real(prec)::tol=0.0025 
        real(prec):: gues1,gues2,gues3,gues4,gues5
        real(prec):: ngues1,ngues2,ngues3,ngues4,ngues5
        real(prec):: fval1,fval2,fval3,fval4,fval5
        real(prec),parameter:: adj=0.2
        real(prec)::errel,errabs
        real(prec),dimension(3,3)::deltamat,deltainv
        real(prec)::x,xguess,fnorm,epss,etas,high,low
        integer::nroot,info,itmax
        
        interface
            subroutine fun(x1,x2,x3,x4,x5,fv1,fv2,fv3,fv4,fv5)
                real(8),intent(in):: x1,x2,x3,x4,x5
                real(8),intent(out):: fv1,fv2,fv3,fv4,fv5
            end subroutine fun
        end interface

        !external taxfn

        do i=1,maxit

        print*,'____________________________________________________________________'
        print*,"Newton iteration ",i


        call fun(gues1,gues2,gues3,gues4,gues5,fval1,fval2,fval3,fval4,fval5)


        ngues1=TFP*alpha*( As/(LabS*(1.0+nn)) )**(alpha-1.0)-delta
        ngues2=LabS
        ngues4=TrBn 
        ngues5=SSn  

        ! With Gouveia-Strauss: Updating a2

        errel	=0.000000001
        errabs	=0.0000000001
        epss	=0.1
        etas	=1.0
        nroot	=1
        itmax	=1000
        xguess	=gues3
        low		=0.00000001
        high	=10.0**20

        !call dzbren(taxfn,errabs,errel,low,high,itmax)
        !ngues3=high
        ngues3=GovconsN

        if ( (abs(fval1)/Y <tol) .and. (abs(fval2)/ngues2<tol) .and. ( abs(fval3)/Y < tol ) .and. ( abs(fval4) < tol ) .and. ( abs(fval5) < tol ) ) then
            print*,'Convergence Achieved'
            exit
        endif


        PRINT '(a)',' ' 
        print '(a)',               "      <variable>       <old guess>    <new guess>   <error>" 
        print '(a,2f15.6,f15.11)', "  (1) interest rate ",gues1,ngues1,fval1/Y
        print '(a,2f15.6,f15.11)', "  (2) labor supply  ",gues2,ngues2,fval2/ngues2
        !print '(a,2f15.6,f15.11)', "  (3) parameter a2  ",gues3,ngues3,fval3/Y
        print '(a,2f15.6,f15.11)', "  (3) Govcons       ",gues3,ngues3,fval3/Y
        print '(a,2f15.6,f15.11)', "  (4) bequest TrB   ",gues4,ngues4,fval4
        print '(a,2f15.6,f15.11)', "  (5) SS benefit SS ",gues5,ngues5,fval5


        PRINT '(a)',' ' 
        !if (indext.eq.1) then 
        !	print '(a,3f10.4,f15.5)', 'Tax System [a0,a1,a2]',a0(a0c),a1(a1c),ngues3
        !else
        !	print '(a,3f10.4,f15.5)', 'Tax System [tauk,a0,a1,a2]',tauk(taukc),a0(a0c),a1(a1c),ngues3
        !endif
        print '(a,f15.10)', 'Excess Dem. Goods Market ',exdem/Y
        print '(a,f15.10)', 'Total Shares in GDP      ',(C+Govcons+As*(delta+nn)/(1.0+nn) )/Y

        print*,' ' 
        !if (indext.eq.1) then
        print*, "CALIBRATION TARGETS ARE: K/Y=2.7 I/Y=0.255 G/Y=0.17 h=1/3"
        print '(a,f10.7,a,f10.7)', "  K/Y=", As/((1.0+nn)*Y),"  I/Y=", (delta+nn)*As/((1.0+nn)*Y)
        print '(a,f10.7,a,f10.7)', "  G/Y=",Govcons/Y,       "  C/Y=", C/Y
        print '(a,f10.7)', "  Avg hrs wrkd= ", hours
        print '(a,f10.7)', "  Avg tax rate= ", Totinctax/(Y-delta*As/(1.0+nn))
        !endif

        gues1=(1.0-adj)*gues1+adj*ngues1
        gues2=(1.0-adj)*gues2+adj*ngues2
        gues3=(1.0-adj)*gues3+adj*ngues3
        gues4=(1.0-adj)*gues4+adj*ngues4
        gues5=(1.0-adj)*gues5+adj*ngues5

        open(11, file='equilibrium_tmp.txt')
        write(11, '(f20.8)') gues1
        write(11, '(f20.8)') gues2
        write(11, '(f20.8)') gues3
        write(11, '(f20.8)') gues4
        write(11, '(f20.8)') gues5
        close(11)

        end do

        print *, ' ' 
        print *, "Convergence achieved in ",i," Iterations"


        PRINT '(a)',' ' 
        print '(a)',               "      <variable>       <old guess>    <new guess>   <error>" 
        print '(a,2f15.6,f15.11)', "  (1) interest rate ",gues1,ngues1,fval1/Y
        print '(a,2f15.6,f15.11)', "  (2) labor supply  ",gues2,ngues2,fval2/Y
        print '(a,2f15.6,f15.11)', "  (3) parameter a2  ",gues3,ngues3,fval3/Y
        print '(a,2f15.6,f15.11)', "  (4) bequest TrB   ",gues4,ngues4,fval4
        print '(a,2f15.6,f15.11)', "  (5) SS benefit SS ",gues5,ngues5,fval5

        !	print*, 'Excess Dem. Goods Mark.',exdem/Y
        !	print*,' ' 
        !
        !	if (indext==1) then
        !		print*, "CALIBRATION OF BENCHMARK ECONOMY"
        !		if (ind_pref==1) then
        !			print*, "with Non-Separable utility"
        !			print '(a,3f10.7)',"(betaNS,deltaNS,gamma)=", beta_NS,delta_NS,gamma
        !		else
        !			print*, "with Separable utility"
        !			print '(a,3f10.7)',"(betaS,deltaS,chi)=", beta_S,delta_S,chi
        !		end if
        !
        !		print '(a,f10.7)', "  K/Y=", As/((1.0+nn)*Y)
        !		print '(a,f10.7)', "  I/Y=", (delta+nn)*As/((1.0+nn)*Y)
        !		print '(a,f10.7)', "  G/Y=", Govcons/Y
        !		print '(a,f10.7)', "  C/Y=", C/Y
        !		print '(a,f10.7)', "  Avg hrs wrkd= ", hours
        !		print '(a,f10.7)', "  Avg tax rate= ", Totinctax/(Y-delta*As/(1.0+nn))
        !		print '(a,3f12.7)',"  Total bequest=",Tr,Y,Tr/Y
        !	else
        !	print '(a,f10.7)', "  K/Y=", As/((1.0+nn)*Y)
        !	print '(a,f10.7)', "  I/Y=", (delta+nn)*As/((1.0+nn)*Y)
        !	print '(a,f10.7)', "  G/Y=", Govcons/Y
        !	print '(a,f10.7)', "  C/Y=", C/Y
        !	print '(a,f10.7)', "  Avg hrs wrkd= ", hours
        !	print '(a,f10.7)', "  Avg tax rate= ", Totinctax/(Y-delta*As/(1.0+nn))
        !	print '(a,3f12.7)',"  Total bequest=",Tr,Y,Tr/Y
        !	end if


        r	= ngues1
        N	= ngues2
        a2	= ngues3
        TrB	= ngues4
        SS	= ngues5

        open(11, file='equilibrium_tmp_backup.txt')
        write(11, '(f20.16)') ngues1
        write(11, '(f20.16)') ngues2
        write(11, '(f20.16)') ngues3
        write(11, '(f20.16)') ngues4
        write(11, '(f20.16)') ngues5
        write(11, '(f20.16)') chi
        write(11, '(f20.16)') beta_S
        close(11)
        
        newton_res(1) = r
        newton_res(2) = As/((1.0+nn)*Y)
        newton_res(3) = hours
        
    end subroutine newton    

    subroutine klp(x, F, nx, m, sim_moms_2save)
        use params
        use moments
        implicit none
        integer, intent(in) :: nx
        integer, intent(in) :: m
        real(8), intent(in) :: x(nx)
        real(8), intent(out) :: F(m)
        real(8), optional, intent(out) :: sim_moms_2save(6) ! Optional output
        !real(8), optional, intent(out) :: pi_save(8, 8) ! Optional output
        !real(8), optional, intent(inout):: eta_save(:) ! Optional output
        real(8) :: guesr,guesN,guesB,guesS
        real(8) :: guesGovcons
        real(8) :: pi_small(3,3)
        real(8) :: p_in, p_out, p_ll, p_hh, p_lh, p_hl
        integer :: i
        real(8) :: p0(1,ns)
        real(8) :: dist
        integer :: i_closest
        real(8) :: tmp
        
        psi_vals = 0.5d0 !100000d0 ![5.00000d0/5.0d0, 20d0/5.0d0, 50d0/5.0d0]
        !psi_vals = [110.00000d0/5.0d0, 173.33333d0/5.0d0, 250.33333d0/5.0d0]!*100
        !psi_vals = [0d0, 110.00000d0/5.0d0, 173.33333d0/5.0d0]!*100
        bbeta = 0.989d0!0.959d0
        chi = 17.4d0

        p0(1,:) = 1d0/dble(ns) ![1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns), 1d0/dble(ns)]

        dist = 1d0
        do i = 1, 500
            pstat = matmul(p0, pi)
            dist = sum( (pstat - p0)**2d0 )
            p0 = pstat
        end do        
        !pini = pstat(1,:)
        !pini(7:8) = 0d0
        !pini = pini/sum(pini)    
        
        open(11, file='equilibrium_tmp.txt')
        read(11, '(f20.8)') guesr
        read(11, '(f20.8)') guesN
        read(11, '(f20.8)') guesGovcons
        read(11, '(f20.8)') guesB
        read(11, '(f20.8)') guesS
        close(11)   
        psi_prob = 1d0 ![1d0/3d0, 1d0/3d0, 1d0/3d0]
        
        call newton(resid,guesr,guesN,guesGovcons,guesB,guesS)
        
        ![wealth_obs_top_01_001, wealth_obs_top_001, gini_wealth_obs, inc_top_01_001, inc_top_001, gini_inc, wealth_obs_top_001]
        ! wealth top1, top01, top001; income top1, top01, top001 
        if (present(sim_moms_2save)) then
            sim_moms_2save = [sim_moms_klp(1:2), sim_moms(2), sim_moms_klp(4:5), sim_moms(5)]
        end if
        
        F(1:6) = (klp_data_moms - sim_moms_klp)/max(klp_data_moms,0.001d0)
        F(7) = (wealth_p_99_99_p_100 - sim_moms(2))/wealth_p_99_99_p_100
        F(8) = (inc_p_99_99_p_100 - sim_moms(5))/inc_p_99_9_p_99_99
        !mom_dist = sum( ((klp_data_moms - sim_moms_klp)/max(klp_data_moms,0.001d0))**2.0d0 ) + &
        !            ((wealth_p_99_99_p_100 - sim_moms(2))/wealth_p_99_99_p_100)**2.0d0 + &
        !            ((inc_p_99_99_p_100 - sim_moms(5))/inc_p_99_9_p_99_99)**2.0d0        
    end subroutine klp
    
    
    subroutine offshoring_test()
        use params
        integer, parameter :: nparam = 6
        integer, parameter :: nx = nparam
        integer, parameter :: nm = 8
        real(8) :: x0(nx)    
        real(8) :: fmom(nm)
        real(8), dimension(6) :: sim_moms_2save
        character(len=50) :: filename
        integer :: i, j_value, unit

        x0 = [0.0004d0, 0.0042d0, 0.9690d0, 0.9424d0, 137.36d0, 1349.46d0]
        !p_in = min( max(x(1), 0d0), 1d0) !0.002d0 ! x(1)
        !p_out = min( max(x(2), 0d0), 1d0) !0.005d0 ! x(2)
        !p_ll = min( max(x(3), 0d0), 1d0) !0.968018785002481d0 ! x(3)
        !plh  = min( max(x(3), 0d0), 1d0)
        !p_hh = min( max(x(4), 0d0), 1d0) !0.946093025947746d0 ! x(4)
        
        x0(5) = x0(5)*1d0
        x0(6) = x0(6)*6d0
        x0(1) = x0(1)/1d0
        x0(2) = x0(2)*5d0
        x0(3) = 0.0210d0/20d0! /4d0


        x0(4) = 0.7424d0
        
        call initialize()

        call klp(x0, fmom, nx, nm, sim_moms_2save)  

        ! Get the filename from the user
        print*, "Enter the filename for saving the CSV file:"
        read*, filename

        ! Open the file for writing
        open(newunit=unit, file=trim(filename), status='replace')

        ! Write column headers
        write(unit, '(A)') 'wealthtop1,top01,top001,incometop1,top01,top001'

        ! Write data to CSV file
        write(unit, '(6(F12.6, ", "))') (sim_moms_2save(i)*100, i=1,6)
        write(unit, '(8(F12.6, ", "))') (pini(i), i=1,8)
        write(unit, '(8(F12.6, ", "))') (DistXW(i)*100, i=1,8)
        write(unit, '(8(F12.6, ", "))') (DistX(i)*100, i=1,8)
        write(unit, '(8(F12.6, ", "))') (eta(i), i=1,8)
        do j_value = 1, 8
            write(unit, '(8(F12.6, ", "))') (pi(j_value, i), i = 1, 8)
        end do
        write(unit, '(A)') 'r,K2Y,hours'
        write(unit, '(6(F12.6, ", "))') (newton_res(i), i=1,3)
        write(unit, '(A)') 'offshoring costs'
        write(unit, '(6(F12.6, ", "))') (psi_vals(i), i=1,3)
        write(unit, '(A)') ',P0-90,P90-99,P99-99.9,P99.9-99.99,P99.99-100, GINI'
        write(unit, '(A, 5(F12.6, ", "))') 'offshoring,', (share_off_dist(i), i=1,5)
        write(unit, '(A, 6(F12.6, ", "))') 'Income,', (income_dist(i), i=1,6)
        write(unit, '(A, 6(F12.6, ", "))') 'Wealth obs,', (wealth_obs_dist(i), i=1,6)
        write(unit, '(A, 6(F12.6, ", "))') 'Wealth tot,', (wealth_tot_dist(i), i=1,6)
        write(unit, '(A, 6(F12.6, ", "))') 'Labor Income,', (lab_income_dist(i), i=1,6)
        ! Close the file
        close(unit)
      
        !! Prompt user for filename
        !write(*, '(A)', advance='no') 'Enter the filename for saving the CSV file: '
        !read(*, '(A)') filename
        !
        !! Open file for writing
        !open(newunit=unit, file=trim(filename), status='replace', action='write')
        !
        !! Write column headers
        !write(unit, '(A,6(A,F12.6,", "))') 'sim_moms_2save:', &
        !'wealthtop1', 'top01', 'top001', 'incometop1', 'top01', 'top001'
        !
        !! Write data to CSV file
        !write(unit, '(6(F12.6, ", "))') (sim_moms_2save(i), i=1,6)
        !
        !! Close the file
        !close(unit)
        !
        !write(*, '(A)', advance='yes') 'CSV file has been saved successfully.'
        
    end subroutine offshoring_test    
    
end module CK_routines