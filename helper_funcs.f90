function fff_static(hrs)
    use Mod_Household, only: lambda_mod
    use params, only: chi, sig2, ep, eta, jc, sc, tyc, w, theta0, theta1, taup, tau_max, yb_cutoff
        
    real(8) :: hrs    
    real(8) :: wtot
    real(8) :: labinc        
    real(8) :: fff_static
        
    wtot = w*eta(sc)*ep(tyc,jc)
    labinc = wtot*hrs
    if (labinc < yb_cutoff) then
        fff_static = chi*hrs**sig2 - lambda_mod*( theta0*(1d0-theta1)*(wtot)**(1d0-theta1)*hrs**(-theta1) - taup*wtot)
    else
        fff_static = chi*hrs**sig2 - lambda_mod*( (1.0d0-tau_max-taup)*wtot )
    end if        
end function fff_static   
    
    
function fff_helper_func1(hrs)
    use Mod_Household, only: nonlabinc_mod
    use params, only: chi, sig2, ep, eta, jc, sc, tyc, w, theta0, theta1, taup, sig1, tauc, tau_max, yb_cutoff
    
    real(8) :: hrs
    real(8) :: fff_helper_func1
    real(8) :: cons
    real(8) :: labinc
        
    labinc = w*eta(sc)*ep(tyc,jc)*hrs
    if (labinc < yb_cutoff) then
        cons = ( theta0*(1d0-theta1)*(w*eta(sc)*ep(tyc,jc))**(1d0-theta1)*hrs**(-theta1) - taup*w*eta(sc)*ep(tyc,jc) )/(1d0+tauc)/chi/(hrs**sig2)
    else
        cons = ( (1-tau_max)*(w*eta(sc)*ep(tyc,jc)) - taup*w*eta(sc)*ep(tyc,jc) )/(1d0+tauc)/chi/(hrs**sig2)    
    end if
        
    cons = cons**(1d0/sig1)
    fff_helper_func1= nonlabinc_mod + theta0*min(labinc,yb_cutoff)**(1d0-theta1) + (1d0-tau_max)*max(labinc-yb_cutoff,0d0) - cons*(1d0+tauc) - taup*labinc       
end function fff_helper_func1 
    
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
    !y_aft = after_tax_income(y) !- psi_mod
    
    net_inc_ofsh = frac_ofsh*y + after_tax_income((1d0-frac_ofsh)*y) - psi_mod
    net_inc_nofsh = after_tax_income(y)  
    net_inc = max(net_inc_ofsh, net_inc_nofsh)     
    
    !dyh = D_after_tax_income(y)
    if (net_inc_ofsh >= net_inc_nofsh) then
        !dD_ = 1.0d0 + (frac_ofsh + D_after_tax_income((1d0-frac_ofsh)*gross_inc)*(1d0-frac_ofsh))*(rcur + acur*drcur)
        dyh = (D_after_tax_income((1d0-frac_ofsh)*y)*(1d0-frac_ofsh) + frac_ofsh)*w_mod
    else
        !dD_ = (1.0d0 + D_after_tax_income(gross_inc)*(rcur + acur*drcur))
        dyh = D_after_tax_income(y)*w_mod
    end if       
    
    f(1) = chi*h**sig2 - lambda_mod*dyh
    f(2) = (1.0d0+tauc)*c_mod + ap_mod - net_inc - a - TrB
    
end subroutine static_focs
    

  