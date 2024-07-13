MODULE TAUCHEN_mod

    IMPLICIT NONE

Contains
    
    subroutine tauchen_pareto(s_e, rho, nz, m, pareto_cutoff, alpha, zvals, prob)
        real(8) :: s_e
        real(8) :: s_z
        real(8) :: rho
        integer :: nz
        real(8) :: m
        real(8) :: pareto_cutoff
        real(8) :: x_norm_cutoff
        real(8) :: z_norm_cutoff        
        real(8) :: alpha
        real(8) :: zvals(nz), logzvals(nz)
        real(8) :: prob(nz,nz)
        real(8) :: pstat(nz), cdfstat(nz), cdfstat_test(nz)      
        integer, parameter :: ntmp=1
        real(8) :: atmp(ntmp)
        real(8) :: ytmp(ntmp)     
        integer :: i
        integer :: npareto
        
        atmp = [pareto_cutoff]
        call vdcdfnorminv(ntmp, atmp, ytmp)
        x_norm_cutoff = ytmp(1)
        s_z = s_e/sqrt(1d0-rho**2d0)
        z_norm_cutoff = x_norm_cutoff*s_z        
        
        call tauchen(s_e, rho, nz, m, logzvals, prob, pstat)
    
        !do i = 1, nz
        !    print *, sum(prob(i,:))
        !end do
        call normal_01_cdf(logzvals(1)/s_z, cdfstat(1))
        do i = 2, nz
            call normal_01_cdf(logzvals(i)/s_z, cdfstat(i))
        end do
    
        !alpha = 1.9d0
        npareto = 0
        do i = 1, nz
            if (cdfstat(i) < pareto_cutoff) then
            !if (logzvals(i) < z_norm_cutoff) then
                zvals(i) = exp(logzvals(i))
            else
                zvals(i) = inv_pareto( (cdfstat(i)-pareto_cutoff)/(1.0d0-pareto_cutoff), exp(z_norm_cutoff), alpha)
                npareto = npareto + 1
            end if
            !zvals_tmp(i) = exp(logzvals(i))
        end do 
        
        print *, npareto, ' values transformed'
    
    end subroutine tauchen_pareto
    

    Subroutine tauchen(STDE,RHO,ZZ,M,SZ,PZEE,pi)
        !=============================================================
        !DISCRETIZE A NORMAL DISTRIBUTION BY USING THE TAUCHEN(1986)
        !==============================================================
        !	Input:  STDW(standard deviation)
        !		    Log(z_t+1)=RHO*Log(z_t)+e_t+1; e~N(0,STDE^2).
        !			M is the number of st. dev. within the distribution is approximated
        !			ZZ is the number of shocks
        !	Ouput:  SZ is the set of shock values,
        !	        PZEE is the transition matrix of z
        !		    pi is a stationary distribution of z

        Implicit NONE

        Real(8), Intent(In) :: STDE,RHO
        Integer, Intent(In) :: ZZ
        Real(8), Intent(Out) :: SZ(ZZ),PZEE(ZZ,ZZ),pi(ZZ)
        INTEGER i,j,k
        Integer, parameter :: draws_number=1000000
        Real(8)  STDW, step, M, temp1,dum
        Real(8), Dimension(draws_number) :: normal_draws
        real(8) :: cdf1, cdf2, tmp

        !M=1.5d0

        !==================================
        !USE PZW TO SET THE VALUE OF SHOCKS
        !==================================

        STDW=STDE/sqrt(1-RHO*RHO)
        SZ(ZZ)=M*STDW

        SZ(1)=-SZ(ZZ)
        step=2*STDW*M/FLOAT(ZZ-1)
        DO i=2,ZZ-1
            SZ(i)=SZ(i-1)+step
        Enddo

        dum=(5D-1)*step

        do i=1,ZZ
            !PZEE(i,1)=D_ANORDF((SZ(1)-RHO*SZ(i)+dum)/STDE)
            call normal_01_cdf((SZ(1)-RHO*SZ(i)+dum)/STDE, cdf2)
            PZEE(i,1) = cdf2
        end do

        do i=2,ZZ-1
            do j=1,ZZ
                !PZEE(j,i)=D_ANORDF((SZ(i)-RHO*SZ(j)+dum)/STDE)-D_ANORDF((SZ(i)-RHO*SZ(j)-dum)/STDE)
                call normal_01_cdf((SZ(i)-RHO*SZ(j)+dum)/STDE, cdf2)
                call normal_01_cdf((SZ(i)-RHO*SZ(j)-dum)/STDE, cdf1)
                PZEE(j,i) = cdf2 - cdf1
            end do
        end do

        do i=1,ZZ
            !PZEE(i,ZZ)=1D0-D_ANORDF((SZ(ZZ)-RHO*SZ(i)-dum)/STDE)
            call normal_01_cdf((SZ(ZZ)-RHO*SZ(i)-dum)/STDE, cdf1)
            PZEE(i,ZZ) = 1.0d0 - cdf1
        end do

        !===========================================================
        !FIND TRANSITION MATRIX BY TAUCHEN(1986):
        !Use Monte Carlo Integration
        !===========================================================

        !Draw a vector normal_draws from the normal distribution.

        !DO  i = 1, draws_number
        !	normal_draws(i) = D_RNNOF()
        !END DO
        !
        !
        !normal_draws=normal_draws*stde

        !Compute a share of draws in each interval.
        !PZEE=0
        !Do j=1,ZZ
        !	temp1=rho*SZ(j)
        !	Do i=1,draws_number
        !		If (normal_draws(i)<=SZ(1)+step/2-temp1) PZEE(j,1)=PZEE(j,1)+1
        !		If (normal_draws(i)> SZ(ZZ)-step/2-temp1) PZEE(j,ZZ)=PZEE(j,ZZ)+1
        !	Enddo
        !Enddo
        !
        !Do j=1,ZZ
        !	temp1=rho*SZ(j)
        !	Do i=1,draws_number
        !		Do k=2, ZZ-1
        !			If ((normal_draws(i)> SZ(k)-step/2-temp1) .and. (normal_draws(i)<=SZ(k)+step/2-temp1)) &
        !				 PZEE(j,k)=PZEE(j,k)+1
        !		Enddo
        !	Enddo
        !Enddo
        !
        !
        !PZEE=PZEE/draws_number
        !
        !DO i=1,ZZ
        !   SZ(i)=EXP(SZ(i))
        !Enddo



        !=============
        !Find the stationary distribution
        !=============
        pi=1d0/ZZ
        Do i=1,200
            pi=Matmul(pi,PZEE)
        Enddo

        RETURN

    END Subroutine tauchen
    
    
     Subroutine tauchen_low(STDE,RHO,ZZ,M1,Mm,SZ,PZEE,pi)
        !=============================================================
        !DISCRETIZE A NORMAL DISTRIBUTION BY USING THE TAUCHEN(1986)
        !==============================================================
        !	Input:  STDW(standard deviation)
        !		    Log(z_t+1)=RHO*Log(z_t)+e_t+1; e~N(0,STDE^2).
        !			Mm is the maximum number of st. dev. within the distribution is approximated
        !           M1 is the minimum number of st. dev. within the distribution is approximated (!negative)
        !			ZZ is the number of shocks
        !	Ouput:  SZ is the set of shock values,
        !	        PZEE is the transition matrix of z
        !		    pi is a stationary distribution of z

        Implicit NONE

        Real(8), Intent(In) :: STDE,RHO
        Integer, Intent(In) :: ZZ
        Real(8), Intent(Out) :: SZ(ZZ),PZEE(ZZ,ZZ),pi(ZZ)
        INTEGER i,j,k
        Integer, parameter :: draws_number=1000000
        Real(8)  STDW, step, M1, Mm, temp1,dum
        Real(8), Dimension(draws_number) :: normal_draws
        real(8) :: cdf1, cdf2, tmp

        !M=1.5d0

        !==================================
        !USE PZW TO SET THE VALUE OF SHOCKS
        !==================================

        STDW=STDE/sqrt(1-RHO*RHO)
        SZ(ZZ)=Mm*STDW

        SZ(1)=M1*STDW! - SZ(ZZ)
        step=STDW*(Mm-M1)/FLOAT(ZZ-1)
        DO i=2,ZZ-1
            SZ(i)=SZ(i-1)+step
        Enddo

        dum=(5D-1)*step

        do i=1,ZZ
            !PZEE(i,1)=D_ANORDF((SZ(1)-RHO*SZ(i)+dum)/STDE)
            call normal_01_cdf((SZ(1)-RHO*SZ(i)+dum)/STDE, cdf2)
            PZEE(i,1) = cdf2
        end do

        do i=2,ZZ-1
            do j=1,ZZ
                !PZEE(j,i)=D_ANORDF((SZ(i)-RHO*SZ(j)+dum)/STDE)-D_ANORDF((SZ(i)-RHO*SZ(j)-dum)/STDE)
                call normal_01_cdf((SZ(i)-RHO*SZ(j)+dum)/STDE, cdf2)
                call normal_01_cdf((SZ(i)-RHO*SZ(j)-dum)/STDE, cdf1)
                PZEE(j,i) = cdf2 - cdf1
            end do
        end do

        do i=1,ZZ
            !PZEE(i,ZZ)=1D0-D_ANORDF((SZ(ZZ)-RHO*SZ(i)-dum)/STDE)
            call normal_01_cdf((SZ(ZZ)-RHO*SZ(i)-dum)/STDE, cdf1)
            PZEE(i,ZZ) = 1.0d0 - cdf1
        end do

        !===========================================================
        !FIND TRANSITION MATRIX BY TAUCHEN(1986):
        !Use Monte Carlo Integration
        !===========================================================

        !Draw a vector normal_draws from the normal distribution.

        !DO  i = 1, draws_number
        !	normal_draws(i) = D_RNNOF()
        !END DO
        !
        !
        !normal_draws=normal_draws*stde

        !Compute a share of draws in each interval.
        !PZEE=0
        !Do j=1,ZZ
        !	temp1=rho*SZ(j)
        !	Do i=1,draws_number
        !		If (normal_draws(i)<=SZ(1)+step/2-temp1) PZEE(j,1)=PZEE(j,1)+1
        !		If (normal_draws(i)> SZ(ZZ)-step/2-temp1) PZEE(j,ZZ)=PZEE(j,ZZ)+1
        !	Enddo
        !Enddo
        !
        !Do j=1,ZZ
        !	temp1=rho*SZ(j)
        !	Do i=1,draws_number
        !		Do k=2, ZZ-1
        !			If ((normal_draws(i)> SZ(k)-step/2-temp1) .and. (normal_draws(i)<=SZ(k)+step/2-temp1)) &
        !				 PZEE(j,k)=PZEE(j,k)+1
        !		Enddo
        !	Enddo
        !Enddo
        !
        !
        !PZEE=PZEE/draws_number
        !
        !DO i=1,ZZ
        !   SZ(i)=EXP(SZ(i))
        !Enddo



        !=============
        !Find the stationary distribution
        !=============
        pi=1d0/ZZ
        Do i=1,200
            pi=Matmul(pi,PZEE)
        Enddo

        RETURN

    END Subroutine tauchen_low
    

    subroutine normal_01_cdf ( x, cdf )

        !*****************************************************************************80
        !
        !! NORMAL_01_CDF evaluates the Normal 01 CDF.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    10 February 1999
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Reference:
        !
        !    AG Adams,
        !    Algorithm 39,
        !    Areas Under the Normal Curve,
        !    Computer Journal,
        !    Volume 12, pages 197-198, 1969.
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) X, the argument of the CDF.
        !
        !    Output, real ( kind = 8 ) CDF, the value of the CDF.
        !
        implicit none

        real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
        real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
        real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
        real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
        real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
        real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
        real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
        real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
        real ( kind = 8 ), parameter :: b1 = 3.8052D-08
        real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
        real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
        real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
        real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
        real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
        real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
        real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
        real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
        real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
        real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
        real ( kind = 8 ) cdf
        real ( kind = 8 ) q
        real ( kind = 8 ) x
        real ( kind = 8 ) y
        !
        !  |X| <= 1.28.
        !
        if ( abs ( x ) <= 1.28D+00 ) then

            y = 0.5D+00 * x * x

            q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
                + a6 / ( y + a7 ) ) ) )
            !
            !  1.28 < |X| <= 12.7
            !
        else if ( abs ( x ) <= 12.7D+00 ) then

            y = 0.5D+00 * x * x

            q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
                + b2 / ( abs ( x ) + b3 &
                + b4 / ( abs ( x ) - b5 &
                + b6 / ( abs ( x ) + b7 &
                - b8 / ( abs ( x ) + b9 &
                + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
            !
            !  12.7 < |X|
            !
        else

            q = 0.0D+00

        end if
        !
        !  Take account of negative X.
        !
        if ( x < 0.0D+00 ) then
            cdf = q
        else
            cdf = 1.0D+00 - q
        end if

        return
    end subroutine normal_01_cdf
    
    
    function inv_pareto(y, xm, alpha) result(f)
        real(8) :: y
        real(8) :: xm
        real(8) :: alpha
        real(8) :: f
    
        f = xm/(1.0d0-y)**(1.0d0/alpha)
    end function inv_pareto 
    

END MODULE TAUCHEN_mod
