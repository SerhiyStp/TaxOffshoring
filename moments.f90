module moments
    
    use params
    
    implicit none
    
    
    real(8), parameter :: wealth_p_0_p_90 = 0.27d0
    real(8), parameter :: wealth_p_90_p_99 = 0.38d0
    real(8), parameter :: wealth_p_99_p_99_9 = 0.17d0
    real(8), parameter :: wealth_p_99_9_p_99_99 = 0.09d0
    real(8), parameter :: wealth_p_99_99_p_100 = 0.08d0
    real(8), parameter :: wealth_gini = 0.83d0
    real(8), parameter :: inc_p_0_p_90 = 0.55d0
    real(8), parameter :: inc_p_90_p_99 = 0.27d0
    real(8), parameter :: inc_p_99_p_99_9 = 0.10d0
    real(8), parameter :: inc_p_99_9_p_99_99 = 0.04d0
    real(8), parameter :: inc_p_99_99_p_100 = 0.03d0
    real(8), parameter :: inc_gini = 0.57d0
    real(8), parameter :: share_off_p_0_p_90 = 0.00d0
    real(8), parameter :: share_off_p_90_p_99 = 0.02d0
    real(8), parameter :: share_off_p_99_p_99_9 = 0.13d0
    real(8), parameter :: share_off_p_99_9_p_99_99 = 0.35d0
    real(8), parameter :: share_off_p_99_99_p_100 = 0.53d0
    real(8), parameter :: K_Y_mom = 2.9d0
    real(8), parameter :: hrs_mom = 0.333d0
    
    real(8), parameter :: klp_wealth_p_99_p_100 = 0.37d0
    real(8), parameter :: klp_wealth_p_99_9_p_100 = 0.14d0
    real(8), parameter :: klp_wealth_gini = 0.85d0
    real(8), parameter :: klp_inc_p_99_p_100 = 0.23d0
    real(8), parameter :: klp_inc_p_99_9_p_100 = 0.08d0
    real(8), parameter :: klp_iinc_gini = 0.67d0
    
    real(8) :: data_moms(n_moms)
    real(8) :: sim_moms(n_moms)
    real(8) :: sim_moms_aux(7)
    real(8) :: data_moms_aux(7)
    real(8) :: KLP_data_moms(klp_n_moms)
    real(8) :: sim_moms_klp(6)
    
contains

    subroutine set_moments()
    
        data_moms = [wealth_p_99_9_p_99_99, &
                     wealth_p_99_99_p_100, &
                     wealth_gini, &
                     inc_p_99_9_p_99_99, &
                     inc_p_99_99_p_100, &
                     inc_gini]  
        
        data_moms_aux = [share_off_p_0_p_90, &
                         share_off_p_90_p_99, &
                         share_off_p_99_p_99_9, &
                         share_off_p_99_9_p_99_99, &
                         share_off_p_99_99_p_100, &
                         K_Y_mom, &
                         hrs_mom] 
        
        klp_data_moms = [klp_wealth_p_99_p_100, &
                         klp_wealth_p_99_9_p_100, &
                         klp_wealth_gini, &
                         klp_inc_p_99_p_100, &
                         klp_inc_p_99_9_p_100, &
                         klp_iinc_gini]
        
    end subroutine set_moments


    
end module moments