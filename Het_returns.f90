module Het_returns
    use nrtype, only: dp
    
    implicit none
    
    real(dp),parameter:: omega1 = 0.072d0
    real(dp),parameter:: omega2 = 0.20d0    
    real(dp),parameter:: gamma_omega = 0.30d0
    real(dp),parameter:: abar_omega = 6.15000016987324d0 !10.0d0 !0.1d0 !0.0d0
    real(dp),parameter:: omegabar = 0.4d0
    ! risk-free rate of return
    real(dp),parameter:: rF = 0.01d0 !0.00d0 !0.0d0 !0.03d0 !0.01d0
    ! risky rate of return
    real(dp),parameter:: rR = 0.15d0 !0.14d0 !0.13d0 !0.12d0 !0.11d0 !0.15d0 !0.09d0 !0.06d0
    ! adjustment factor
    real(dp):: rbar
    ! temporary shocks to risky rate of return 
    integer, parameter :: nkappa = 3 !1 !3 
    real(dp),parameter:: sig_kappa = 0.05d0
    real(dp):: Kappas(nkappa)
    real(dp):: pi_kappa(nkappa)    
    ! persistent types in risky rates of return
    integer, parameter :: ntheta=2 !1 !2     
    real(8) :: thetas(ntheta)  
    real(8) :: pi_theta(ntheta,ntheta)    
    real(8) :: pi_theta_stat(ntheta)
    
contains
    
    subroutine SET_RETURNS()
    
        implicit none
        
        real(dp), parameter :: pi_theta_hl = 0.5d0
        real(dp), parameter :: pi_theta_lh = 0.2d0
        real(dp) :: p0(1, ntheta), p1(1, ntheta)
        real(dp) :: dist
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
        
    end subroutine SET_RETURNS       
    
    
    function rfunc(a, theta, kappa)
        
        implicit none
        
        real(dp), intent(in) :: a
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: kappa
        real(dp) :: om
        real(dp) :: rfunc
        
        om = Omega(a,theta)
        rfunc = max(0d0, rbar + rF*(1.0d0-om) + rR*kappa*om)
        
    end function rfunc
    
    function Da_rfunc(a, theta, kappa)
        
        implicit none
        
        real(8), intent(in) :: a
        real(8), intent(in) :: theta
        real(8), intent(in) :: kappa
        real(8) :: Da_om
        real(8) :: Da_rfunc
        
        Da_om = DaOmega(a, theta)
        Da_rfunc = (rR*kappa - rF)*Da_om
    end function Da_rfunc      
    
    function rfunc_PE(a, theta, kappa)
        
        implicit none
        
        real(8), intent(in) :: a
        real(8), intent(in) :: theta
        real(8), intent(in) :: kappa
        real(8) :: om
        real(8) :: rfunc_PE

        om = Omega(a,theta)
        rfunc_PE = rR*kappa*om

    end function rfunc_PE    
    
    function Omega(a, theta)
        
        implicit none
        
        real(8), intent(in) :: a
        real(8), intent(in) :: theta
        real(8) :: Omega
        
        Omega = theta*(omegabar + min( omega1*( max(a - abar_omega, 0.0d0) )**gamma_omega, omega2 ) )
    end function Omega  
    
    function DaOmega(a, theta)
        
        implicit none
        
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
    
  
    
end module Het_returns
