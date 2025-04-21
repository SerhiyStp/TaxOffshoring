module CK_routines
    
    implicit none
    !real(8):: newton_res(3)
    
contains
    
  
    
    subroutine resid(x1,x2,x3,x4,fv1,fv2,fv3,fv4)
        use params
        use Taxes
        use Het_returns, only: rF, rR, sig_kappa, rbar
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
        integer :: iu
        
        r   = x1
        !N   = x2
        rbar = x2
        !a2  = x3
        !Govcons = x3
        theta0 = x3
        TrB = x4 
        !SS  = x5 
        
        open(newunit=iu, file='last_x.txt')
        write(iu, '(f20.16)') x1
        write(iu, '(f20.16)') x2
        write(iu, '(f20.16)') x3
        write(iu, '(f20.16)') x4
        write(iu, '(a25,f20.16)') 'rF: ', rF 
        write(iu, '(a25,f20.16)') 'rR: ', rR
        write(iu, '(a25,f20.16)') 'sig_kappa: ', sig_kappa
        write(iu, '(a25,f20.16)') 'sig_z: ', sig_z
        write(iu, '(a25,f20.16)') 'rho_z: ', rho_z
        write(iu, '(a25,f20.16)') 'sig_xi: ',sig_xi
        write(iu, '(a25,f20.16)') 'frisch: ', frisch
        write(iu, '(a25,f20.16)') 'sigma (risk aversion): ', sig1
        write(iu, '(a25,f20.16)') 'delta: ', delta_S
        write(iu, '(a25,f20.16)') 'alpha (production): ', alpha
        write(iu, '(a25,f20.16)') 'TFP: ', TFP
        write(iu, '(a25,f20.16)') 'chi: ', chi
        write(iu, '(a25,f20.16)') 'beta: ', bbeta
        close(iu)

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
            !OPEN(NEWUNIT=iunit, FILE=outDir // "tmp.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
            !WRITE (iunit) afun
            !write (iunit) lfun
            !write (iunit) cfun
            !write (iunit) offshoring
            !write (iunit) afun_ret
            !write (iunit) cfun_ret
            !write (iunit) offshoring_ret
            !CLOSE(iunit)
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
        call GetDistribution(save_res=.false.)
        call SummarizeDistribution()
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
        print *, 'K/Y = ', K/Y
        print *, 'I/Y = ', delta*K/Y
    end subroutine resid
    
    subroutine newton(fun,gues1,gues2,gues3,gues4)

        ! This subroutine computes the steady state interest rate and labor supply and intercept for tax system using the 
        ! the classical newton method; inputs are the guesses for r, N and Tint and the subroutine calls resid that
        ! delivers the residual from markets clearing in the asset market, labor market and the gvernment budget constraint

        use params
        use Taxes

        implicit none

        integer::i
        real(prec)::tol=0.0001d0 !0.000001d0 !0.0025 
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

            ngues1= max( TFP*alpha*( AAgg/(LAgg*(1.0+nn)) )**(alpha-1.0)-delta, 0.01d0)
            !ngues2=LabS
            ngues2 = ngues1 - RauxAgg/AAgg 
            !ngues4=gues4
            !Govcons = 25.5490651400000d0
            ngues3 = max( (YfBelowYb + TaxC + TaxIncAboveYb + TaxE - Govcons - RetAgg)/DBelowYb, 0.01d0)
            if (ngues3 < 0d0) then
                print *,'WARING: Negative guess for theta0'
            end if
            ngues4=TrBn 
            !ngues5=SSn  
            !ngues5=gues5

            !d1 = abs(ngues1-gues1)
            !d2 = abs(ngues2-gues2)
            !d3 = abs(ngues3-gues3)
            !d4 = abs(ngues4-gues4)
            !print '(5(a7,f12.7))', 'd1 = ', d1, ' d2 = ', d2, ' d3 = ', d3, ' d4 = ', d4

            ! With Gouveia-Strauss: Updating a2

            !errel	=0.000000001
            !errabs	=0.0000000001
            !epss	=0.1
            !etas	=1.0
            !nroot	=1
            !itmax	=1000
            !xguess	=gues3
            !low		=0.00000001
            !high	=10.0**20

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
        !r	= ngues1
        !N	= ngues2
        !!a2	= ngues3
        !TrB	= ngues4
        !!SS	= ngues5
        !
        !open(11, file='equilibrium_tmp_backup.txt')
        !write(11, '(f20.16)') ngues1
        !write(11, '(f20.16)') ngues2
        !write(11, '(f20.16)') ngues3
        !write(11, '(f20.16)') ngues4
        !!write(11, '(f20.16)') ngues5
        !write(11, '(f20.16)') chi
        !write(11, '(f20.16)') beta_S
        !close(11)
        !
        !newton_res(1) = r
        !newton_res(2) = AAgg/((1.0+nn)*Y)
        !newton_res(3) = hours
    end subroutine newton    
    
end module CK_routines
