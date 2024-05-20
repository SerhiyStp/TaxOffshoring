module int_tictoc

! this module provides a set of variables and two subroutines
! "call tic" saves the current time (in seconds) in tic_time
! "call toc" saves the time elapsed (in seconds) since the last "call tic" in toc_time
! the calls can be nested, and the routines keep track of the level using tictoc_counter and tic_time_vec
!
! the variabels tictoc_number and tictoc_sum are meant to keep track of the total and average time
! it takes to complete a set of instructions, as follows
!
!   call tic
!   instructions
!   call toc
!   tictoc_number=tictoc_number+1
!   tictoc_sum=tictoc_sum+toc_time


integer, dimension(8)     :: tictoc_values
real(8)                   :: tic_time, toc_time
logical                   :: tictoc_print=.TRUE.
integer                   :: tictoc_counter=10
real(8), dimension(100)   :: tic_time_vec
integer, dimension(100)   :: tictoc_number=0
real(8), dimension(100)   :: tictoc_sum=0.0_8

contains
    subroutine tic

        tictoc_counter=tictoc_counter+1
        call date_and_time(VALUES=tictoc_values)
        tic_time=86400.0_8*tictoc_values(3)+3600.0_8*tictoc_values(5)+60.0_8*tictoc_values(6)+tictoc_values(7)+0.001_8*tictoc_values(8)    
        tic_time_vec(tictoc_counter)=tic_time
    end subroutine tic
    
    subroutine toc
    
        call date_and_time(VALUES=tictoc_values)        
        toc_time=86400.0_8*tictoc_values(3)+3600.0_8*tictoc_values(5)+60.0_8*tictoc_values(6)+tictoc_values(7)+0.001_8*tictoc_values(8)-tic_time_vec(tictoc_counter)
        tictoc_counter=tictoc_counter-1
        if (tictoc_print) print *, 'time elapsed in seconds: ', toc_time
        
    end subroutine toc
    
    subroutine tuc
    
        call date_and_time(VALUES=tictoc_values)        
        toc_time=86400.0_8*tictoc_values(3)+3600.0_8*tictoc_values(5)+60.0_8*tictoc_values(6)+tictoc_values(7)+0.001_8*tictoc_values(8)-tic_time_vec(tictoc_counter)
        if (tictoc_print) print ('(f6.2)'), 'time elapsed in seconds: ', toc_time
        
    end subroutine tuc
    
end module int_tictoc