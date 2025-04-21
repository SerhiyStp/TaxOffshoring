module Taxes
    ! Tax code and government spending related parameters and functions
    use nrtype, only: dp
    
    implicit none
    
    ! tax level 
    real(dp) :: theta0
    ! tax progressivity
    real(dp),parameter:: theta1=0.183d0 !0.137d0
    ! top marginal tax rate
    real(dp),parameter:: tau_max = 0.396d0
    ! income cutoff for the top marginal tax rate
    real(dp) :: yb_cutoff
    ! consumption tax rate
    real(dp),parameter:: tauc=0.05d0 	    
    ! bequests tax rate
    real(dp), parameter :: tau_estate = 0.1d0
    real(dp), parameter :: a_estate = 0.5d0		

    real(dp), parameter :: Govcons = 0.5d0 !25.5490651400000d0
    real(dp) :: GovconsN    
    
    
contains

    function tax_income(y) result(res)
        real(8), intent(in) :: y
        real(8) :: res
        
        res = y - after_tax_income(y)
    end function tax_income
    
    
    function D_after_tax_income(y) result(res)
        
        implicit none
        
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
    
end module Taxes