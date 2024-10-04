module CK_routines
    
    implicit none
    real(8):: newton_res(3)
    
contains
    
    subroutine initialize()
        !use moments
        use PARAMS, only: SetParams
        use Mod_Distribution, only: init_distr

        call SetParams()
        !call set_moments()
        call init_distr()
        
        !call allocate_policy_fns()
        
    end subroutine initialize    
    
    subroutine resid(x1,x2,x3,x4,fv1,fv2,fv3,fv4)
        use params
        use Mod_Household
        use Mod_Distribution
        use int_tictoc

        implicit none

        real(prec),intent(in):: x1,x2,x3,x4 !,x5
        real(prec),intent(out):: fv1,fv2,fv3,fv4 !,fv5
        CHARACTER (LEN=*), PARAMETER :: outDir = "tmp/"
        INTEGER :: iunit
        integer :: get_new_soln
        real(8) :: KN, w_test
        real(8) :: Y2, exdem2
        
        r   = x1
        !N   = x2
        rbar = x2
        !a2  = x3
        !Govcons = x3
        theta0 = x3
        TrB = x4 
        !SS  = x5 

        !K = N*( (alpha*TFP) / (r+delta) )**(1.0/(1.0-alpha))        ! Capital Stock
        !Y = TFP*(K**alpha)*(N**(1.0-alpha))					        ! Aggregate Output
        !w = (1.0-alpha)*Y/N											! Wages
        
        KN = ( alpha*TFP/(r+delta) )**(1.0d0/(1.0d0-alpha))
        w = TFP*(1.0d0-alpha)*KN**alpha


        !maxSS = maxSSrat*Y/sum(Nu)  
        
        ! Solve the Household Problem       
        get_new_soln = 1
        if (get_new_soln == 1) then
            call tic()
            call SolveHH(save_res=.false.)
            print *, 'Solution took: '
            call toc()
            OPEN(NEWUNIT=iunit, FILE=outDir // "tmp.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
            WRITE (iunit) afun
            write (iunit) lfun
            write (iunit) cfun
            write (iunit) offshoring
            write (iunit) afun_ret
            write (iunit) cfun_ret
            write (iunit) offshoring_ret
            CLOSE(iunit)
        else
            OPEN(NEWUNIT=iunit, FILE=outdir // "tmp.bin", FORM="unformatted", ACCESS="stream", STATUS="old")
            READ (iunit) afun
            READ (iunit) lfun
            read (iunit) cfun
            read (iunit) offshoring
            read (iunit) afun_ret
            read (iunit) cfun_ret  
            read (iunit) offshoring_ret
            CLOSE(iunit)    
        end if
        call tic()
        call Distribution(save_res=.false.)
        print *, 'Simulation took: '
        call toc()


        ! Compute residuals of functions we want to set to zero
        !fv1=As-K*(1.0+nn)
        !fv2=LabS-N
        !fv3=Govcons-tauc*C-Totinctax
        !fv4=TrB-TrBn
        !fv5=SS-SSn
        K   = KN*LAgg
        fv1 = AAgg/LAgg - KN !K/N
        fv2 = r*K - RAgg
        fv3 = Govcons + RetAgg - TaxTot 
        fv4 = TrB - TrBn
        
        !print '(4(a7,f12.7))', 'fv1 = ', fv1, ' fv2 = ', fv2, ' fv3 = ', fv3, ' fv4 = ', fv4

        GovconsN = tauc*C + Totinctax

        ! Goods Market Clearing
        Y	 = TFP * (K**alpha) * (LAgg**(1.0-alpha))	
        exdem=( CAgg + AAgg + TotOffshCost + Govcons - (1.0d0-delta)*K  - Y )
        
        Y2 = TFP * (AAgg**alpha) * (LAgg**(1.0-alpha))
        exdem2 = ( CAgg + AAgg + TotOffshCost + Govcons - (1.0d0-delta)*AAgg  - Y2 )
        print *, 'K - AAgg = ', K - AAgg
        print *, 'exdem = ', exdem
        print *, 'exdem2 = ', exdem2
    end subroutine resid
    
    subroutine newton(fun,gues1,gues2,gues3,gues4)

        ! This subroutine computes the steady state interest rate and labor supply and intercept for tax system using the 
        ! the classical newton method; inputs are the guesses for r, N and Tint and the subroutine calls resid that
        ! delivers the residual from markets clearing in the asset market, labor market and the gvernment budget constraint

        use params

        implicit none

        integer::i
        real(prec)::tol=0.000001d0 !0.0025 
        real(prec):: gues1,gues2,gues3,gues4 !,gues5
        real(prec):: ngues1,ngues2,ngues3,ngues4 !,ngues5
        real(prec):: fval1,fval2,fval3,fval4 !,fval5
        real(prec),parameter:: adj=0.2
        real(prec)::errel,errabs
        real(prec)::x,xguess,fnorm,epss,etas,high,low
        integer::nroot,info,itmax
        !real(8) :: d1, d2, d3, d4, d5
        
        interface
            subroutine fun(x1,x2,x3,x4,fv1,fv2,fv3,fv4)
                real(8),intent(in):: x1,x2,x3,x4 
                real(8),intent(out):: fv1,fv2,fv3,fv4 
            end subroutine fun
        end interface

        !external taxfn

        do i=1,maxit

            print*,'____________________________________________________________________'
            print*,"Newton iteration ",i

            call fun(gues1,gues2,gues3,gues4,fval1,fval2,fval3,fval4)

            ngues1=TFP*alpha*( AAgg/(LAgg*(1.0+nn)) )**(alpha-1.0)-delta
            !ngues2=LabS
            ngues2 = gues1 - RauxAgg/AAgg 
            !ngues4=gues4
            !Govcons = 25.5490651400000d0
            ngues3 = (YfBelowYb + TaxC + TaxIncAboveYb + TaxE - Govcons - RetAgg)/DBelowYb
            ngues4=TrBn 
            !ngues5=SSn  
            !ngues5=gues5

            !d1 = abs(ngues1-gues1)
            !d2 = abs(ngues2-gues2)
            !d3 = abs(ngues3-gues3)
            !d4 = abs(ngues4-gues4)
            !print '(5(a7,f12.7))', 'd1 = ', d1, ' d2 = ', d2, ' d3 = ', d3, ' d4 = ', d4

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
            !ngues3=GovconsN

            if ( (abs(fval1)/Y <tol) .and. (abs(fval2)/ngues2<tol) .and. ( abs(fval3)/Y < tol ) .and. ( abs(fval4) < tol ) ) then
                print*,'Convergence Achieved'
                exit
            endif


            print '(a)',' ' 
            print '(a)', "      <variable>       <old guess>      <new guess>     <error>" 
            print '(a20,2f15.6,f15.9)', " (1) interest rate  ", gues1, ngues1, fval1 !/Y
            print '(a20,2f15.6,f15.9)', " (2) rbar           ", gues2, ngues2, fval2
            print '(a20,2f15.6,f15.9)', " (3) theta0         ", gues3, ngues3, fval3
            print '(a20,2f15.6,f15.9)', " (4) Tr             ", gues4, ngues4, fval4
            !print '(a,2f15.6,f15.11)', "  (2) labor supply  ",gues2,ngues2,fval2/ngues2
            !!print '(a,2f15.6,f15.11)', "  (3) parameter a2  ",gues3,ngues3,fval3/Y
            !print '(a,2f15.6,f15.11)', "  (3) Govcons       ",gues3,ngues3,fval3/Y
            !print '(a,2f15.6,f15.11)', "  (4) bequest TrB   ",gues4,ngues4,fval4
            !print '(a,2f15.6,f15.11)', "  (5) SS benefit SS ",gues5,ngues5,fval5


            !PRINT '(a)',' ' 
            !print '(a,f15.10)', 'Excess Dem. Goods Market ',exdem/Y
            !print '(a,f15.10)', 'Total Shares in GDP      ',(C+Govcons+As*(delta+nn)/(1.0+nn) )/Y

            !print*,' ' 
            !!if (indext.eq.1) then
            !print*, "CALIBRATION TARGETS ARE: K/Y=2.7 I/Y=0.255 G/Y=0.17 h=1/3"
            !print '(a,f10.7,a,f10.7)', "  K/Y=", As/((1.0+nn)*Y),"  I/Y=", (delta+nn)*As/((1.0+nn)*Y)
            !print '(a,f10.7,a,f10.7)', "  G/Y=",Govcons/Y,       "  C/Y=", C/Y
            !print '(a,f10.7)', "  Avg hrs wrkd= ", hours
            !print '(a,f10.7)', "  Avg tax rate= ", Totinctax/(Y-delta*As/(1.0+nn))
            !endif

            gues1=(1.0-adj)*gues1+adj*ngues1
            gues2=(1.0-adj)*gues2+adj*ngues2
            gues3=(1.0-adj)*gues3+adj*ngues3
            gues4=(1.0-adj)*gues4+adj*ngues4
            !gues5=(1.0-adj)*gues5+adj*ngues5

            ! open(11, file='equilibrium_tmp.txt')
            ! write(11, '(f20.8)') gues1
            ! write(11, '(f20.8)') gues2
            ! write(11, '(f20.8)') gues3
            ! write(11, '(f20.8)') gues4
            ! write(11, '(f20.8)') gues5
            ! close(11)

        end do

        print *, ' ' 
        print *, "Convergence achieved in ",i," Iterations"


        PRINT '(a)',' ' 
        print '(a)',               "      <variable>       <old guess>    <new guess>   <error>" 
        print '(a,2f15.6,f15.11)', "  (1) interest rate ",gues1,ngues1,fval1/Y
        print '(a,2f15.6,f15.11)', "  (2) labor supply  ",gues2,ngues2,fval2/Y
        print '(a,2f15.6,f15.11)', "  (3) parameter a2  ",gues3,ngues3,fval3/Y
        print '(a,2f15.6,f15.11)', "  (4) bequest TrB   ",gues4,ngues4,fval4
        !print '(a,2f15.6,f15.11)', "  (5) SS benefit SS ",gues5,ngues5,fval5

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
        !SS	= ngues5

        open(11, file='equilibrium_tmp_backup.txt')
        write(11, '(f20.16)') ngues1
        write(11, '(f20.16)') ngues2
        write(11, '(f20.16)') ngues3
        write(11, '(f20.16)') ngues4
        !write(11, '(f20.16)') ngues5
        write(11, '(f20.16)') chi
        write(11, '(f20.16)') beta_S
        close(11)
        
        newton_res(1) = r
        newton_res(2) = AAgg/((1.0+nn)*Y)
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
        real(8) :: guesr,guesN,guesB,guesS, guesrb, guestheta0
        real(8) :: guesGovcons
        real(8) :: pi_small(3,3)
        real(8) :: p_in, p_out, p_ll, p_hh, p_lh, p_hl
        integer :: i
        real(8) :: p0(1,ns)
        real(8) :: dist
        integer :: i_closest
        real(8) :: tmp
        
        !psi_vals = [0.1d0, 0.5d0, 1.5d0] 
        !bbeta = 0.989d0!0.959d0
        !chi = 17.4d0
        
        !open(11, file='equilibrium_tmp.txt')
        !read(11, '(f20.8)') guesr
        !read(11, '(f20.8)') guesN
        !read(11, '(f20.8)') guesGovcons
        !read(11, '(f20.8)') guesB
        !read(11, '(f20.8)') guesS
        !close(11) 
        
        psi_prob = [1d0/3d0, 1d0/3d0, 1d0/3d0] !1d0 !
        
        guesr = 0.04d0
        guesrb = 0d0
        guestheta0 = 0.940d0
        guesB = 0d0
        guesS = 0d0
        call newton(resid,guesr,guesrb,guestheta0,guesB)
        
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
