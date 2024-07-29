subroutine static_focs(x, f, n)
    use Mod_Household, only: lambda_mod, w_mod, theta_mod, kappa_mod, c_mod, ap_mod, psi_mod
    use PARAMS, only: rfunc, d_after_tax_income, after_tax_income, sig2, chi, tauc, TrB, frac_ofsh
    implicit none
    integer :: n
    real(8) :: x(n), f(n)
    real(8) :: a, h
    real(8) :: y, dyh, y_aft
    real(8) :: rcur
    real(8) :: net_inc_ofsh, net_inc_nofsh, net_inc
    
    a = x(1)
    h = x(2)
    
    rcur = rfunc(a, theta_mod, kappa_mod)
    y = rcur*a + w_mod*h
    net_inc_ofsh = frac_ofsh*y + after_tax_income((1d0-frac_ofsh)*y) - psi_mod
    net_inc_nofsh = after_tax_income(y)  
    net_inc = max(net_inc_ofsh, net_inc_nofsh)     
    if (net_inc_ofsh >= net_inc_nofsh) then
        dyh = (D_after_tax_income((1d0-frac_ofsh)*y)*(1d0-frac_ofsh) + frac_ofsh)*w_mod
    else
        dyh = D_after_tax_income(y)*w_mod
    end if       
    f(1) = lambda_mod*dyh - chi*h**sig2
    f(2) = net_inc + a + TrB - (1.0d0+tauc)*c_mod - ap_mod 
end subroutine static_focs
    
subroutine static_focs_nofsh(x, f, n)
    use Mod_Household, only: lambda_mod, w_mod, theta_mod, kappa_mod, c_mod, ap_mod, psi_mod
    use PARAMS, only: rfunc, d_after_tax_income, after_tax_income, sig2, chi, tauc, TrB, frac_ofsh
    implicit none
    integer :: n
    real(8) :: x(n), f(n)
    real(8) :: a, h
    real(8) :: y, dyh, y_aft
    real(8) :: rcur
    real(8) :: net_inc
    
    a = x(1)
    h = x(2)
    rcur = rfunc(a, theta_mod, kappa_mod)
    y = rcur*a + w_mod*h
    net_inc = after_tax_income(y)       
    dyh = D_after_tax_income(y)*w_mod      
    f(1) = lambda_mod*dyh - chi*h**sig2
    f(2) = net_inc + a + TrB - (1.0d0+tauc)*c_mod - ap_mod 
end subroutine static_focs_nofsh    
    
    
subroutine static_focs_ofsh(x, f, n)
    use Mod_Household, only: lambda_mod, w_mod, theta_mod, kappa_mod, c_mod, ap_mod, psi_mod
    use PARAMS, only: rfunc, d_after_tax_income, after_tax_income, sig2, chi, tauc, TrB, frac_ofsh
    implicit none
    integer :: n
    real(8) :: x(n), f(n)
    real(8) :: a, h
    real(8) :: y, dyh, y_aft
    real(8) :: rcur
    real(8) :: net_inc
    
    a = x(1)
    h = x(2)
    rcur = rfunc(a, theta_mod, kappa_mod)
    y = rcur*a + w_mod*h
    net_inc = frac_ofsh*y + after_tax_income((1d0-frac_ofsh)*y) - psi_mod    
    dyh = (D_after_tax_income((1d0-frac_ofsh)*y)*(1d0-frac_ofsh) + frac_ofsh)*w_mod      
    f(1) = lambda_mod*dyh - chi*h**sig2
    f(2) = net_inc + a + TrB - (1.0d0+tauc)*c_mod - ap_mod 
end subroutine static_focs_ofsh    
    

  