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