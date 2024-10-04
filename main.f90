program main
    
    use CK_routines
    
    implicit none
    !call offshoring_test()
    call test_run()
    
    print *, 'End of program'
    
contains
    
    subroutine test_run()
        real(8) :: guesr,guesN,guesB,guesS, guesrb, guestheta0
        
        call initialize()
        
        guesr = 0.038275d0 !0.04d0  !
        guesrb = -0.002272d0 !0d0    !
        guestheta0 = 0.833632d0 !0.940d0  ! 
        guesB = 0.031670d0 !0d0     ! 
        !guesS = 0d0
        call newton(resid,guesr,guesrb,guestheta0,guesB)    
        
    end subroutine test_run
    
end program main