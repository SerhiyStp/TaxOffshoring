module io
    
    implicit none
    
    interface save_array
        module procedure save_array_d6, save_array_d7
    end interface
    
    interface read_array
        module procedure read_array_d6, read_array_d7
    end interface    
    
contains
    
    subroutine save_array_d6(array, path)
        real(8), intent(in) :: array(:,:,:,:,:,:)
        character(len=*), intent(in) :: path
        integer :: iunit
        OPEN(NEWUNIT=iunit, FILE=path, FORM="unformatted", ACCESS="stream", STATUS="unknown")
        WRITE (iunit) array
        CLOSE(iunit)        
    end subroutine save_array_d6    
    
    subroutine save_array_d7(array, path)
        real(8), intent(in) :: array(:,:,:,:,:,:,:)
        character(len=*), intent(in) :: path
        integer :: iunit
        OPEN(NEWUNIT=iunit, FILE=path, FORM="unformatted", ACCESS="stream", STATUS="unknown")
        WRITE (iunit) array
        CLOSE(iunit)        
    end subroutine save_array_d7
    
    subroutine read_array_d6(array, path)
        real(8) :: array(:,:,:,:,:,:)
        character(len=*), intent(in) :: path
        integer :: iunit
        OPEN(NEWUNIT=iunit, FILE=path, FORM="unformatted", ACCESS="stream", STATUS="unknown")
        read (iunit) array
        CLOSE(iunit)        
    end subroutine read_array_d6    
    
    subroutine read_array_d7(array, path)
        real(8) :: array(:,:,:,:,:,:,:)
        character(len=*), intent(in) :: path
        integer :: iunit
        OPEN(NEWUNIT=iunit, FILE=path, FORM="unformatted", ACCESS="stream", STATUS="unknown")
        read (iunit) array
        CLOSE(iunit)        
    end subroutine read_array_d7    
         
end module io