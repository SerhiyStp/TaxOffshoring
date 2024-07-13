module focs_mod
    
    implicit none
    
contains

    function static_focs_BL_h(hrs) 
        use params, only: chi, sig2, sig1, d_after_tax_income, after_tax_income, TrB, tauc, frac_ofsh
        !use Mod_Household, only: rcur_mod, w_mod, theta_mod, kappa_mod, c_mod, a_mod, ap_mod
        use glob_vars_mod, only: r_glob, a_glob, w_glob, ap_glob, c_glob, psi_glob
        !
        real(8), intent(in) :: hrs
        real(8) :: y, dyh, c, y_aft
        real(8) :: net_inc_ofsh, net_inc_nofsh, net_inc
        real(8) :: static_focs_BL_h
        !
        y = r_glob*a_glob + w_glob*hrs
        !y_aft = after_tax_income(y)
        !dyh = D_after_tax_income(y)
        
        net_inc_ofsh = frac_ofsh*y + after_tax_income((1d0-frac_ofsh)*y) - psi_glob
        net_inc_nofsh = after_tax_income(y)  
        net_inc = max(net_inc_ofsh, net_inc_nofsh)     
    
        !dyh = D_after_tax_income(y)
        if (net_inc_ofsh >= net_inc_nofsh) then
            !dD_ = 1.0d0 + (frac_ofsh + D_after_tax_income((1d0-frac_ofsh)*gross_inc)*(1d0-frac_ofsh))*(rcur + acur*drcur)
            dyh = (D_after_tax_income((1d0-frac_ofsh)*y)*(1d0-frac_ofsh) + frac_ofsh)*w_glob
        else
            !dD_ = (1.0d0 + D_after_tax_income(gross_inc)*(rcur + acur*drcur))
            dyh = D_after_tax_income(y)*w_glob
        end if             
        
        
        c = (dyh/(1.0d0+tauc)/chi/hrs**sig2)**(1.0d0/sig1)
        c_glob = c
        
        static_focs_BL_h = (1.0d0+tauc)*c + ap_glob - y_aft - a_glob - TrB
        !static_focs_BL_h = 0d0
    end function static_focs_BL_h    
    
end module focs_mod