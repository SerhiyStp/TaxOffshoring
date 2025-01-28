program main
    
    use CK_routines
    
    implicit none
    !call offshoring_test()
    call test_run()
    
    print *, 'End of program'
    
contains
    
    subroutine test_run()
        real(8) :: guesr,guesN,guesB,guesS, guesrb, guestheta0
        integer :: iu
        
        call initialize()
        
        guesr = 0.036667d0 !0.036577d0 !0.036488d0 !0.036414d0 !0.036348d0  !0.036096d0 !0.036942d0 !0.038275d0 !0.04d0  !
        guesrb = -0.003176d0 !-0.001714d0 !-0.000257d0 !0.001214d0 !0.002687d0 !-0.002977d0 !0.002133d0 !-0.002272d0 !0d0    !
        guestheta0 = 0.865982d0 !0.865968d0 !0.865992d0 !0.865804d0 !0.865382d0 !0.865004d0 !0.835079d0 !0.833632d0 !0.940d0  ! 
        guesB = 0.032077d0 !0.032079d0 !0.032080d0 !0.032085d0 !0.032102d0 !0.032104d0 !0.032256d0 !0.031670d0 !0d0     ! 
        !guesS = 0d0
        
        open(newunit=iu, file='last_x.txt')
        read(iu, '(f20.16)') guesr
        read(iu, '(f20.16)') guesrb
        read(iu, '(f20.16)') guestheta0
        read(iu, '(f20.16)') guesB
        close(iu)        
        
        call newton(resid,guesr,guesrb,guestheta0,guesB)    
        
    end subroutine test_run
    
end program main