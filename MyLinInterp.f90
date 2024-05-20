module MyLinInterp 
    
    implicit none
    
contains
    
    function LinInterp_1d(x,xGrid,fVals,nx)
        real(8) :: LinInterp_1d
        real(8) :: x
        integer :: nx
        real(8) :: xGrid(nx), fVals(nx)
        integer :: i
        real(8) :: b
        
!        i = 2
!        do while (x > xGrid(i) .and. i < nx)
!            i = i + 1
!        end do
        do i = 1, nx-1
            if (xGrid(i+1) >= x) exit
        end do
        i = min(i+1,nx)
        
        b = (x - xGrid(i-1))/(xGrid(i) - xGrid(i-1))
        
        LinInterp_1d = (1d0-b)*fVals(i-1) + b*fVals(i)
        
        !if(b>1.001d0) Then
        !    LinInterp_1d=fVals(nx)+((x - xGrid(nx))/(xGrid(nx) - xGrid(nx-1)))*(fVals(nx)-fVals(nx-1))
        !end if
        
    end function LinInterp_1d
    
    function locate(xx,x,n) result(loc)
        integer :: loc
        integer :: n
        real(8) :: xx(n)
        real(8) :: x
        integer :: jl, ju, jm
        
        if (x <= xx(1)) then
            loc=1
        else if (x >= xx(n)) then
            loc=n-1
        else
            jl = 1
            ju = n   
            do
                if (ju-jl <= 1) exit
                jm=(ju+jl)/2
                if (x >= xx(jm)) then
                    jl=jm
                else
                    ju=jm
                end if
            end do            
            loc=jl
        end if        
        
    end function locate
    
    function lininterp_2d(xx,yy,V,nx,ny,arg) result(val)
        real(8) :: val
        integer :: nx, ny
        !real(8), dimension (:), allocatable :: xx, yy
        !real(8), dimension (:,:), allocatable :: V
        real(8) :: xx(nx), yy(ny)
        real(8) :: V(nx,ny)
        real(8) :: arg(2)
        real(8) :: x, y
        integer :: jx, jy
        real(8) :: dx, dy
        
        x = arg(1)
        y = arg(2)
        
        jx = locate(xx,x,nx)
        jy = locate(yy,y,ny)   
        dx = (x - xx(jx))/(xx(jx+1) - xx(jx))
        dy = (y - yy(jy))/(yy(jy+1) - yy(jy))        
    
        val = V(jx,jy)*(1d0-dx)*(1d0-dy) + V(min(jx+1,nx),jy)*dx*(1d0-dy) &
            + V(jx,min(jy+1,ny))*(1d0-dx)*dy + V(min(jx+1,nx),min(jy+1,ny))*dx*dy
    
    end function lininterp_2d
    
    function lininterp_3d(xx,yy,zz,V,nx,ny,nz,arg) result(val)
        real(8) :: val
        integer :: nx, ny, nz
        real(8) :: xx(nx), yy(ny), zz(nz)
        real(8) :: V(nx,ny,nz)
        real(8) :: arg(3)
        real(8) :: x, y, z
        integer :: jx, jy, jz
        real(8) :: dx, dy, dz
        real(8) :: c00, c10, c01, c11
        real(8) :: c0, c1
        
        x = arg(1)
        y = arg(2)
        z = arg(3)
        
        jx = locate(xx,x,nx)
        jy = locate(yy,y,ny)
        jz = locate(zz,z,nz)
        
        dx = (x - xx(jx))/(xx(jx+1) - xx(jx))
        dy = (y - yy(jy))/(yy(jy+1) - yy(jy))
        dz = (z - zz(jz))/(zz(jz+1) - zz(jz))
        
        c00 = V(jx,jy,jz)*(1-dx)                     + V(min(jx+1,nx),jy,jz)*dx
        c10 = V(jx,min(jy+1,ny),jz)*(1-dx)           + V(min(jx+1,nx),min(jy+1,ny),jz)*dx
        c01 = V(jx,jy,min(jz+1,nz))*(1-dx)           + V(min(jx+1,nx),jy,min(jz+1,nz))*dx
        c11 = V(jx,min(jy+1,ny),min(jz+1,nz))*(1-dx) + V(min(jx+1,nx),min(jy+1,ny),min(jz+1,nz))*dx
        
        c0 = c00*(1-dy) + c10*dy
        c1 = c01*(1-dy) + c11*dy
        
        val = c0*(1-dz) + c1*dz
    
    end function lininterp_3d

    function lininterp_4d(xx,yy,zz,ww,V,nx,ny,nz,nw,arg) result(val)
        real(8) :: val
        integer :: nx, ny, nz,nw
        real(8) :: xx(nx), yy(ny), zz(nz), ww(nw)
        real(8) :: V(nx,ny,nz,nw)
        real(8) :: arg(4)
        real(8) :: x, y, z, w
        integer,dimension(2) :: jx, jy, jz, jw
        real(8) :: dx, dy, dz, dw
        real(8) :: px, py, pz, pw
        real(8) :: fx111, fx211, fx121, fx112, fx221, fx212, fx122, fx222
        real(8) :: fxy11, fxy12, fxy21, fxy22
        real(8) :: fxyz1, fxyz2
        
        
        x = arg(1)
        y = arg(2)
        z = arg(3)
        w = arg(4)
        
        jx(1) = locate(xx,x,nx)
        jy(1) = locate(yy,y,ny)
        jz(1) = locate(zz,z,nz)
        jw(1) = locate(ww,w,nw)

        jx(2) = jx(1)+1
        jy(2) = jy(1)+1
        jz(2) = jz(1)+1
        jw(2) = jw(1)+1
        
        dx = (x - xx(jx(1)))/(xx(jx(2)) - xx(jx(1)))
        dy = (y - yy(jy(1)))/(yy(jy(2)) - yy(jy(1)))
        dz = (z - zz(jz(1)))/(zz(jz(2)) - zz(jz(1)))
        dw = (w - ww(jw(1)))/(ww(jw(2)) - ww(jw(1)))

        px = 1d0-dx
        py = 1d0-dy
        pz = 1d0-dz
        pw = 1d0-dw
        
        fx111 = V(jx(1),jy(1),jz(1),jw(1))*px + V(jx(2),jy(1),jz(1),jw(1))*dx
        fx211 = V(jx(1),jy(2),jz(1),jw(1))*px + V(jx(2),jy(2),jz(1),jw(1))*dx
        fx121 = V(jx(1),jy(1),jz(2),jw(1))*px + V(jx(2),jy(1),jz(2),jw(1))*dx
        fx221 = V(jx(1),jy(2),jz(2),jw(1))*px + V(jx(2),jy(2),jz(2),jw(1))*dx
        fx112 = V(jx(1),jy(1),jz(1),jw(2))*px + V(jx(2),jy(1),jz(1),jw(2))*dx
        fx212 = V(jx(1),jy(2),jz(1),jw(2))*px + V(jx(2),jy(2),jz(1),jw(2))*dx
        fx122 = V(jx(1),jy(1),jz(2),jw(2))*px + V(jx(2),jy(1),jz(2),jw(2))*dx
        fx222 = V(jx(1),jy(2),jz(2),jw(2))*px + V(jx(2),jy(2),jz(2),jw(2))*dx
        
        fxy11 = fx111*py + fx211*dy
        fxy21 = fx121*py + fx221*dy
        fxy12 = fx112*py + fx212*dy
        fxy22 = fx122*py + fx222*dy        
        
        fxyz1 = fxy11*pz + fxy21*dz
        fxyz2 = fxy12*pz + fxy22*dz
        
        val = fxyz1*pw + fxyz2*dw
    
    end function lininterp_4d    

end module MyLinInterp
