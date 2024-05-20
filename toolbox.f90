MODULE toolbox
	!USE mkl95_precision, ONLY: WP => DP
    USE f95_precision, ONLY: WP => DP
CONTAINS
SUBROUTINE lorenz(f,x,fx,gini)
! Compute Lorenz curve and Gini coefficient by sorting vector x and using density f
! Note: to link dlasrt2 add mkl_scalapack_core.lib to project's additional dependencies
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: f
	REAL(WP), DIMENSION(:), INTENT(INOUT) :: x
	REAL(WP), DIMENSION(:), INTENT(OUT) :: fx
	REAL(WP), INTENT(OUT) :: gini
	!INTEGER, DIMENSION(size(x)) :: key
    INTEGER, ALLOCATABLE :: key(:)
	INTEGER :: n,i,info
    
	n=size(x)
    allocate(key(n))
	IF (size(f)/=n) THEN
		PRINT '(a,i3)', 'lorenz: x and f must be of the same size ',n
		STOP 'program terminated by lorenz'
	END IF
	IF (size(fx)/=n) THEN
		PRINT '(a,i3)', 'lorenz: x and fx must be of the same size ',n
		STOP 'program terminated by lorenz'
	END IF
	key=(/ (i,i=1,n) /)
	CALL dlasrt2('I',n,x,key,info)
	CALL check('dlasrt2',info)
	fx=f(key)
	x=x*fx
	gini=x(1)*fx(1)
	DO i=2,n
		x(i)=x(i)+x(i-1)
		gini=gini+(x(i)+x(i-1))*fx(i)
		fx(i)=fx(i)+fx(i-1)
	END DO
	gini=1-gini/x(n)
	x=x/x(n)
    deallocate(key)
END SUBROUTINE lorenz

SUBROUTINE lorenz_mod(f,x,fx,gini,x_sorted,f_sorted)
! Compute Lorenz curve and Gini coefficient by sorting vector x and using density f
! Note: to link dlasrt2 add mkl_scalapack_core.lib to project's additional dependencies
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: f
	REAL(WP), DIMENSION(:), INTENT(INOUT) :: x
	REAL(WP), DIMENSION(:), INTENT(OUT) :: fx
    real(wp), dimension(:), intent(out) :: x_sorted
    real(wp), dimension(:), intent(out) :: f_sorted
	REAL(WP), INTENT(OUT) :: gini
	!INTEGER, DIMENSION(size(x)) :: key
    INTEGER, ALLOCATABLE :: key(:)
	INTEGER :: n,i,info
    
	n=size(x)
    allocate(key(n))
	IF (size(f)/=n) THEN
		PRINT '(a,i3)', 'lorenz: x and f must be of the same size ',n
		STOP 'program terminated by lorenz'
	END IF
	IF (size(fx)/=n) THEN
		PRINT '(a,i3)', 'lorenz: x and fx must be of the same size ',n
		STOP 'program terminated by lorenz'
	END IF
	key=(/ (i,i=1,n) /)
	CALL dlasrt2('I',n,x,key,info)
	CALL check('dlasrt2',info)
    x_sorted = x
	fx=f(key)
    f_sorted=f(key)
	x=x*fx
	gini=x(1)*fx(1)
	DO i=2,n
		x(i)=x(i)+x(i-1)
		gini=gini+(x(i)+x(i-1))*fx(i)
		fx(i)=fx(i)+fx(i-1)
	END DO
	gini=1-gini/x(n)
	x=x/x(n)
    deallocate(key)
END SUBROUTINE lorenz_mod

SUBROUTINE ucase(string)
! Purpose: To shift a character string to upper case on any processor, regardless of collating sequence.
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: string  
INTEGER :: i,length
length = len(string)
DO i=1,length
   IF (lge(string(i:i),'a') .and. lle(string(i:i),'z')) THEN
      string(i:i)=achar(iachar(string(i:i))-32)
   END IF
END DO
END SUBROUTINE ucase

SUBROUTINE tauchen(rho,sigma,p,y,s,spread)
! Tauchen method to approximate univariate AR(1) process by Markov chain
!      y(t) = rho y(t-1)+ sigma sqrt(1-rho^2) e(t),   e(t)~N(0,1)
! INPUTS: rho - serial correlation coefficient,
!         sigma - coefficient of variation
! OUTPUT: P is an n-by-n matrix of Markov transition probabilities
!         y is an n-by-1 vector of symmetric and evenly-spaced Markov state space
!         s is an n-by-1 vector of stationary distribution
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: rho, sigma
	REAL(WP), DIMENSION(:,:), INTENT(OUT) :: p
	REAL(WP), DIMENSION(:), INTENT(OUT) :: y,s
	REAL(WP), INTENT(IN), OPTIONAL :: spread
	INTEGER :: n
	REAL(WP) :: smin,emin
	n=size(y)
	IF (size(p,dim=1)/=n .or. size(p,dim=2)/=n) THEN
		PRINT '(a,i3,a,i3)', 'tauchen: p must be a square matrix of size ',n,' x ',n
		STOP 'program terminated by tauchen'
	END IF
	IF (size(s)/=n) THEN
		PRINT '(a,i3)', 'tauchen: y and s must be vectors of the same size ',n
		STOP 'program terminated by tauchen'
	END IF
	IF (present(spread)) THEN
		CALL tauch(spread)
	ELSE
		emin=brent(1.0_WP,2.5_WP,4.0_WP,err,smin)
!		PRINT '(a,f,4x,a,e)', 'smin=',smin,'emin=',emin
		CALL tauch(smin)
	END IF
CONTAINS
SUBROUTINE tauch(m)
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: m
	REAL(WP) :: ybar, dy, sr
	REAL(WP), DIMENSION(size(y)) :: yd
	INTEGER :: i
	ybar=m*sigma
	CALL grid(y,-ybar,ybar,1.0_WP)
	dy=ybar/(n-1)
	sr=sigma*sqrt(1-rho**2)
	DO i=1,n
		yd=(y-rho*y(i)+dy)/sr
		CALL vdcdfnorm(n,yd,p(i,:))
	END DO
	p(:,n)=1
	DO i=n,2,-1
		p(:,i)=p(:,i)-p(:,i-1)
	END DO
	CALL ergodic(p,s)
END SUBROUTINE tauch
FUNCTION err(m)
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: m
	REAL(WP) :: err, rho_, sigma_
	
	CALL tauch(m)
	CALL markovtest(p,y,s,rho_,sigma_)
	err=(log(1-rho_)/log(1-rho)-1)**2+(sigma_/sigma-1)**2;
END FUNCTION err
END SUBROUTINE tauchen

SUBROUTINE markovtest(p,y,s,rho,sigma)
	!USE mkl95_blas, ONLY: gemv
	USE blas95, ONLY: gemv
	IMPLICIT NONE
	REAL(WP), DIMENSION(:,:), INTENT(IN) :: p
	REAL(WP), DIMENSION(:), INTENT(IN) :: y,s
	REAL(WP), INTENT(OUT) :: rho, sigma
	REAL(WP), DIMENSION(size(y)) :: py
	REAL(WP) :: Eyy, Ey2, E2y
	INTEGER :: n
	n=size(y)
	IF (size(p,dim=1)/=n .or. size(p,dim=2)/=n) THEN
		PRINT '(a,i3,a,i3)', 'markovtest: p must be a square matrix of size ',n,' x ',n
		STOP 'program terminated by markovtest'
	END IF
	IF (size(s)/=n) THEN
		PRINT '(a,i3)', 'markovtest: y and s must be vectors of the same size ',n
		STOP 'program terminated by markovtest'
	END IF
	CALL gemv(p,y,py)
	Eyy=sum(s*y*py)
	Ey2=sum(s*y**2)
	E2y=sum(s*y)**2;
	rho=(Eyy-E2y)/(Ey2-E2y);
	sigma=sqrt(Ey2-E2y);
END SUBROUTINE markovtest

SUBROUTINE rouwenhorst(rho,sigma,p,y,s)
! Rowenhurst method to approximate univariate AR(1) process by Markov chain
!      y(t) = rho y(t-1)+ sigma sqrt(1-rho^2) e(t),   e(t)~N(0,1)
! INPUTS: rho - serial correlation coefficient,
!         sigma - coefficient of variation
! OUTPUT: P is an n-by-n matrix of Markov transition probabilities
!         y is an n-by-1 vector of symmetric and evenly-spaced Markov state space
!         s is an n-by-1 vector of stationary distribution (binomial)
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: rho, sigma
	REAL(WP), DIMENSION(:,:), INTENT(OUT) :: p
	REAL(WP), DIMENSION(:), INTENT(OUT) :: y,s
	REAL(WP) :: ybar, q
	INTEGER :: n
	n=size(y)
	IF (size(p,dim=1)/=n .or. size(p,dim=2)/=n) THEN
		PRINT '(a,i3,a,i3)', 'rouwenhorst: p must be a square matrix of size ',n,' x ',n
		STOP 'program terminated by rouwenhorst'
	END IF
	IF (size(s)/=n) THEN
		PRINT '(a,i3)', 'rouwenhorst: y and s must be vectors of the same size ',n
		STOP 'program terminated by rouwenhorst'
	END IF
	ybar=sigma*sqrt(real(n-1,WP))
	q=(1+rho)/2
	CALL rhmat(p)
	CALL grid(y,-ybar,ybar,1.0_WP)
	CALL binom(s)
CONTAINS
RECURSIVE SUBROUTINE rhmat(p)
	IMPLICIT NONE
	REAL(WP), DIMENSION(:,:), INTENT(OUT) :: p
	REAL(WP), DIMENSION(size(p,dim=1)-1,size(p,dim=2)-1) :: p1
	INTEGER :: h
	h=size(p,dim=1)
	IF (size(p,dim=2)/=h) STOP 'P must be a square matrix'
	IF (h<2) STOP 'P must be at least 2-by-2 matrix'
	IF (h==2) THEN
		p=reshape((/q,1-q,1-q,q/),(/2,2/))
	ELSE
		CALL rhmat(p1)
		p=0
		p(1:h-1,1:h-1)=q*p1
		p(1:h-1,2:h)=(1-q)*p1+p(1:h-1,2:h)
		p(2:h,1:h-1)=(1-q)*p1+p(2:h,1:h-1)
		p(2:h,2:h)=q*p1+p(2:h,2:h)
		p(2:h-1,:)=p(2:h-1,:)/2
	END IF
END SUBROUTINE rhmat
END SUBROUTINE rouwenhorst

SUBROUTINE binom(f)
! Binomial probability mass function with p=1/2
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(OUT) :: f
	INTEGER :: n,k
	n=size(f)
	f(1)=2.0_WP**(1-n)
	DO k=1,n-1
		f(k+1)=f(k)*(n-k)/k
	END DO
END SUBROUTINE binom

SUBROUTINE ergodic(p,s)
! Purpose: Compute ergodic distribution s of Markov transition matrix p
	!USE mkl95_lapack, ONLY: geev
	USE lapack95, ONLY: geev
	IMPLICIT NONE
	REAL(WP), DIMENSION(:,:), INTENT(IN) :: p
	REAL(WP), DIMENSION(:), INTENT(OUT) :: s
	REAL(WP), DIMENSION(size(s),size(s)) :: ip,vl
	REAL(WP), DIMENSION(size(s)) :: wr,wi,dw
	INTEGER :: m,info,uw,w1(1)
	REAL(WP) :: ds
	
	m=size(s)
	IF (size(p,dim=1)/=m .or. size(p,dim=2)/=m) THEN
		PRINT '(a,i3,a,i3)', 'sd: p must be a square matrix of size ',m,' x ',m
		STOP 'program terminated by sd'
	END IF
	ip=p
	CALL geev(ip,wr,wi,vl=vl,info=info)
	CALL check('geev',info)
	dw=abs(sqrt(wr*wr+wi*wi)-1)
	w1=minloc(dw)
	uw=count(dw<1000*epsilon(dw))
	IF (uw<1) PRINT '(a)', 'Warning: No unitary eigenvalue is found. Stationary distribution of Markov chain does not exist.'
	IF (uw>1) PRINT '(a)', 'Warning: More than one unitary eigenvalue is found. Stationary distribution of Markov chain is not unique.'
	IF (uw<1 .or. uw>1) PRINT '(a,f20.15,a,f20.15)', 'Using eigenvalue ',wr(w1(1)),'+i',wi(w1(1))
	s=vl(:,w1(1))/sum(vl(:,w1(1)))
	IF (any(s<0)) THEN
		PRINT '(a)', 'The stationary distribution of Markov chain has negative values. Rebalancing...'
		ds=sum(s,mask=s<0)/count(s>=0)
		WHERE(s<0)
			s=0
		ELSEWHERE
			s=s+ds
		ENDWHERE
	END IF
END SUBROUTINE ergodic

SUBROUTINE grid(x,xmin,xmax,s)
! Purpose: Generate grid x on [xmin,xmax] using spacing parameter s set as follows:
! s=1		linear spacing
! s>1		left skewed grid spacing with power s
! 0<s<1		right skewed grid spacing with power s
! s<0		geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
! s=-1		logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
! s=0		logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(OUT) :: x
	REAL(WP), INTENT(IN) :: xmin,xmax,s
	REAL(WP) :: c ! growth rate of grid subintervals for logarithmic spacing
	INTEGER :: n,i
	n=size(x)
	FORALL(i=1:n) x(i)=(i-1)/real(n-1,WP)
	IF (s>0.0_WP) THEN
		x=x**s*(xmax-xmin)+xmin
		IF (s==1.0_WP) THEN
!			PRINT '(a,i8,a,f6.3,a,f6.3,a)', 'Using ',n,' equally spaced grid points over domain [',xmin,',',xmax,']'
		ELSE
!			PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' skewed spaced grid points with power ',s,' over domain [',xmin,',',xmax,']'
		END IF
	ELSE
		IF (s==-1.0_WP) THEN
			c=xmax-xmin+1
!		ELSEIF (s==0.0_WP) THEN
!			IF (xmin>0.0_WP) THEN
!				c=xmax/xmin
!			ELSE
!				STOP 'grid: can not use logarithmic spacing for nonpositive values'
!			END IF
		ELSE
			c=-s
		END IF
!		PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' logarithmically spaced grid points with growth rate ',c,' over domain [',xmin,',',xmax,']'
		x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
	END IF
END SUBROUTINE grid

SUBROUTINE bsearch(func,imin,imax,is,ys)
! Binary Search
! Purpose: Find the maximum ys=func(is) of a strictly concave function func defined over integer grid [imin,imax]
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: imin,imax
	INTEGER, INTENT(OUT) :: is
	REAL(WP), INTENT(OUT) :: ys
	INTERFACE
		FUNCTION func(x)
		!USE mkl95_precision, ONLY: WP => DP
        USE f95_precision, ONLY: WP => DP
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: x
		REAL(WP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER :: il,im,iu,i3(3),i3s(1)
	REAL(WP),DIMENSION(3) :: y3

	IF (imin>imax) STOP 'imin>imax'
	IF (imin==imax) THEN
		is=imax
		ys=func(imax)
	ELSE
		il=imin
		iu=imax
		DO
			IF (iu-il <= 2) EXIT
			im=(iu+il)/2
			IF (func(im+1)>func(im)) THEN
				il=im
			ELSE
				iu=im+1
			END IF
		END DO
		i3=(/ il,il+1,iu /)
		y3=(/ func(il), func(il+1), func(iu) /)
		i3s=maxloc(y3)
		is=i3(i3s(1))
		ys=y3(i3s(1))
	END IF
END SUBROUTINE bsearch

SUBROUTINE maxrhs(func,dfunc,x,i,imax,xs,ys)
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x
	INTEGER, INTENT(IN) :: i,imax
	REAL(WP), INTENT(OUT) :: xs,ys
	INTERFACE
		FUNCTION func(x)
		!USE mkl95_precision, ONLY: WP => DP
        USE f95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP) :: func
		END FUNCTION func
		FUNCTION dfunc(x)
		!USE mkl95_precision, ONLY: WP => DP
        USE f95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	REAL(WP), PARAMETER :: E=0.0001_WP
	REAL(WP) :: xl,xu,xi,yi,xe,ye

	IF (i>1) THEN
		xl=x(i-1)
	END IF
	xi=x(i)
	yi=func(xi)
	IF (i<imax) THEN
		xu=x(i+1)
	END IF
	IF (i>1 .and. i<imax) THEN
			ys=dbrent(xl,xi,xu,func,dfunc,xs)
!			ys=golden(xl,xu,func,xs)
	ELSE
		IF (i==1) THEN
			xe=xi+(xu-xi)*E
			ye=func(xe)
			IF (ye<yi) THEN
				xs=xi
				ys=yi
			ELSE
					ys=dbrent(xi,(xi+xu)/2,xu,func,dfunc,xs)
!					ys=golden(xi,xu,func,xs)
			END IF
		ELSE
			IF (i==imax) THEN
				xe=xi-(xi-xl)*E
				ye=func(xe)
				IF (ye<yi) THEN
					xs=xi
					ys=yi
				ELSE
					ys=dbrent(xl,(xl+xi)/2,xi,func,dfunc,xs)
!					ys=golden(xl,xi,func,xs)
				END IF
			ELSE 
				STOP 'i<1 or i>imax - check binary search algorithm'
			END IF
		END IF
	END IF
END SUBROUTINE maxrhs

SUBROUTINE spline(x,y,y2,yp1,ypn)
	!USE mkl95_lapack, ONLY: pttrf, pttrs
	USE lapack95, ONLY: pttrf, pttrs
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(WP), DIMENSION(:), INTENT(OUT) :: y2
	REAL(WP), INTENT(IN), OPTIONAL :: yp1,ypn
	REAL(WP), DIMENSION(size(x)) :: b
	REAL(WP), DIMENSION(size(x)-1) :: a
	INTEGER :: n,info
	n=size(x)
	IF (size(y)/=n .or. size(y2)/=n) THEN
		PRINT *, 'spline: x,y and y2 must be of the same size'
		STOP 'program terminated by spline'
	END IF
	a=x(2:n)-x(1:n-1)
	b(1)=a(1)
	b(n)=a(n-1)
	b(2:n-1)=x(3:n)-x(1:n-2)
	b=2.0_WP*b
	y2(1:n-1)=(y(2:n)-y(1:n-1))/a
	IF (present(ypn)) THEN
		y2(n)=ypn
	ELSE
		y2(n)=y2(n-1)
		a(n-1)=0.0_WP
	END IF
	y2(2:n)=y2(2:n)-y2(1:n-1)
	IF (present(yp1)) THEN
		y2(1)=y2(1)-yp1
	ELSE
		y2(1)=0.0_WP
		a(1)=0.0_WP
	END IF
	y2=6.0_WP*y2
	CALL pttrf(b,a,info)
	CALL check('pttrf',info)
	CALL pttrs(b,a,y2,info)
	CALL check('pttrs',info)
END SUBROUTINE spline

FUNCTION splint(x,y,y2,xi)
! cubic interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y,y2
	REAL(WP), INTENT(IN) :: xi
	REAL(WP) :: splint
	REAL(WP) :: a,b,d
	INTEGER :: n,i
	n=size(x)
	IF (size(y)/=n .or. size(y2)/=n) THEN
		PRINT *, 'splint: x,y and y2 must be of the same size'
		STOP 'program terminated by splint'
	END IF
	i=max(min(locate(x,xi),n-1),1)
	d=x(i+1)-x(i)
	if (d == 0.0) STOP 'bad x input in splint'
	a=(x(i+1)-xi)/d
	b=(xi-x(i))/d
	splint=a*y(i)+b*y(i+1)+((a**3-a)*y2(i)+(b**3-b)*y2(i+1))*(d**2)/6.0_WP
END FUNCTION splint

FUNCTION dsplint(x,y,y2,xi)
! cubic interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y,y2
	REAL(WP), INTENT(IN) :: xi
	REAL(WP) :: dsplint
	REAL(WP) :: a,b,d
	INTEGER :: n,i
	n=size(x)
	IF (size(y)/=n .or. size(y2)/=n) THEN
		PRINT *, 'splint: x,y and y2 must be of the same size'
		STOP 'program terminated by splint'
	END IF
	i=max(min(locate(x,xi),n-1),1)
	d=x(i+1)-x(i)
	if (d == 0.0) STOP 'bad x input in dsplint'
	a=(x(i+1)-xi)/d
	b=(xi-x(i))/d
	dsplint=(y(i+1)-y(i))/d+((3*b**2-1)*y2(i+1)-(3*a**2-1)*y2(i))*d/6.0_WP
END FUNCTION dsplint

FUNCTION linint(x,y,xi)
! linear interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(WP), INTENT(IN) :: xi
	REAL(WP) :: linint
	REAL(WP) :: a,b,d
	INTEGER :: n,i
	n=size(x)
	IF (size(y)/=n) THEN
		PRINT *, 'linint: x and y must be of the same size'
		STOP 'program terminated by linint'
	END IF
	i=max(min(locate(x,xi),n-1),1)
	d=x(i+1)-x(i)
	IF (d == 0.0) STOP 'bad x input in splint'
	a=(x(i+1)-xi)/d
	b=(xi-x(i))/d
	linint=a*y(i)+b*y(i+1)
END FUNCTION linint

SUBROUTINE linintv(x,y,xi,yi)
! linear interpolation of function y on grid x at interpolation vector xi
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN)  :: x,y,xi
	REAL(WP), DIMENSION(:), INTENT(OUT) :: yi
	REAL(WP) :: a,b,d
	INTEGER :: m,n,i,j
	n=size(x)
	IF (size(y)/=n) THEN
		PRINT *, 'linintv: x and y must be of the same size'
		STOP 'program terminated by linintv'
	END IF
	m=size(xi)
	IF (size(yi)/=m) THEN
		PRINT *, 'linintv: xi and yi must be of the same size'
		STOP 'program terminated by linintv'
	END IF
	DO j=1,m
		i=max(min(locate(x,xi(j)),n-1),1)
		d=x(i+1)-x(i)
		IF (d == 0.0) THEN
			STOP 'bad x input in linintv'
		END IF
		a=(x(i+1)-xi(j))/d
		b=(xi(j)-x(i))/d
		yi(j)=a*y(i)+b*y(i+1)
	END DO
END SUBROUTINE linintv

PURE FUNCTION locate(xx,x)
	IMPLICIT NONE
	REAL(WP), DIMENSION(:), INTENT(IN) :: xx
	REAL(WP), INTENT(IN) :: x
	INTEGER :: locate
	INTEGER :: n,il,im,iu
	n=size(xx)
	il=0
	iu=n+1
	do
		if (iu-il <= 1) exit
		im=(iu+il)/2
		if (x >= xx(im)) then
			il=im
		else
			iu=im
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=il
	end if
END FUNCTION locate

SUBROUTINE check(func,info)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(IN) :: func
	INTEGER, INTENT(IN) :: info
	IF (info<0) THEN
		PRINT '(a,a,a,i2,a)', 'Error in ',func,': the ',-info,'-th parameter had an illegal value'
	END IF
	IF (info>0) THEN
		SELECT CASE(func)
		CASE('geev')
			PRINT '(a,a,a)', 'Error in ',func,': the QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed.'
		CASE('getrf')
			PRINT '(a,a,a)', 'Error in ',func,': the factorization has been completed, but U is exactly singular.'
		CASE DEFAULT
			PRINT '(a,a,a,i2)', 'Error in ',func,': info=',info
		END SELECT
	END IF
	IF (info/=0) STOP 'program terminated by check'
END SUBROUTINE check

FUNCTION brent(ax,bx,cx,func,xmin)
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: ax,bx,cx
	REAL(WP), INTENT(OUT) :: xmin
	REAL(WP) :: brent
	INTERFACE
		FUNCTION func(x)
		!USE mkl95_precision, ONLY: WP => DP
        USE f95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=100
	REAL(WP), PARAMETER :: TOL=sqrt(epsilon(ax)),CGOLD=0.381966011250105_WP,ZEPS=1.0e-3_WP*epsilon(ax)
	INTEGER :: iter
	REAL(WP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	fx=func(x)
	fv=fx
	fw=fx
	do iter=1,ITMAX
		xm=0.5_WP*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_WP*tol1
		if (abs(x-xm) <= (tol2-0.5_WP*(b-a))) then
			xmin=x
			brent=fx
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0_WP*(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5_WP*q*etemp) .or. &
				p <= q*(a-x) .or. p >= q*(b-x)) then
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			else
				d=p/q
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
			end if
		else
			e=merge(a-x,b-x, x >= xm )
			d=CGOLD*e
		end if
		u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
		fu=func(u)
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			v=w
			fv=fw
			w=x
			fw=fx
			x=u
			fx=fu
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				v=w
				fv=fw
				w=u
				fw=fu
			else if (fu <= fv .or. v == x .or. v == w) then
				v=u
				fv=fu
			end if
		end if
	end do
	STOP 'brent: exceed maximum iterations'
END FUNCTION brent

FUNCTION dbrent(ax,bx,cx,func,dfunc,xmax)
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: ax,bx,cx
	REAL(WP), INTENT(OUT) :: xmax
	REAL(WP) :: dbrent
	INTERFACE
		FUNCTION func(x)
		!USE mkl95_precision, ONLY: WP => DP
        USE f95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP) :: func
		END FUNCTION func
		FUNCTION dfunc(x)
		!USE mkl95_precision, ONLY: WP => DP
        USE f95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=100
	REAL(WP), PARAMETER :: TOL=sqrt(epsilon(ax)), ZEPS=1.0e-3_WP*epsilon(ax)
	INTEGER :: iter
	REAL(WP) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
	LOGICAL :: ok1,ok2
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0_WP
	fx=func(x)
	fv=fx
	fw=fx
	dx=dfunc(x)
	dv=dx
	dw=dx
	do iter=1,ITMAX
		xm=0.5_WP*(a+b)
		tol1=TOL*abs(x)+ZEPS
		tol2=2.0_WP*tol1
		if (abs(x-xm) <= (tol2-0.5_WP*(b-a))) exit
		if (abs(e) > tol1) then
			d1=2.0_WP*(b-a)
			d2=d1
			if (dw /= dx) d1=(w-x)*dx/(dx-dw)
			if (dv /= dx) d2=(v-x)*dx/(dx-dv)
			u1=x+d1
			u2=x+d2
			ok1=((a-u1)*(u1-b) > 0.0_WP) .and. (dx*d1 >= 0.0_WP)
			ok2=((a-u2)*(u2-b) > 0.0_WP) .and. (dx*d2 >= 0.0_WP)
			olde=e
			e=d
			if (ok1 .or. ok2) then
				if (ok1 .and. ok2) then
					d=merge(d1,d2, abs(d1) < abs(d2))
				else
					d=merge(d1,d2,ok1)
				end if
				if (abs(d) <= abs(0.5_WP*olde)) then
					u=x+d
					if (u-a < tol2 .or. b-u < tol2) &
						d=sign(tol1,xm-x)
				else
					e=merge(a,b, dx <= 0.0_WP)-x
					d=0.5_WP*e
				end if
			else
				e=merge(a,b, dx <= 0.0_WP)-x
				d=0.5_WP*e
			end if
		else
			e=merge(a,b, dx <= 0.0_WP)-x
			d=0.5_WP*e
		end if
		if (abs(d) >= tol1) then
			u=x+d
			fu=func(u)
		else
			u=x+sign(tol1,d)
			fu=func(u)
			if (fu < fx) exit
		end if
		du=dfunc(u)
		if (fu >= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			v=w
			fv=fw
			dv=dw
			w=x
			fw=fx
			dw=dx
			x=u
			fx=fu
			dx=du
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu >= fw .or. w == x) then
				v=w
				fv=fw
				dv=dw
				w=u
				fw=fu
				dw=du
			else if (fu >= fv .or. v == x .or. v == w) then
				v=u
				fv=fu
				dv=du
			end if
		end if
	end do
	if (iter > ITMAX) STOP 'dbrent: exceeded maximum iterations'
	xmax=x
	dbrent=fx
END FUNCTION dbrent

FUNCTION golden(ax,bx,func,xmax)
! Golden search
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: ax,bx
	REAL(WP), INTENT(OUT) :: xmax
	REAL(WP), PARAMETER :: TOL=sqrt(epsilon(1.0_WP))
	REAL(WP) :: golden
	INTERFACE
		FUNCTION func(x)
		!USE mkl95_precision, ONLY: WP => DP
        USE f95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(WP), PARAMETER :: R=0.618033988749895_WP,C=1.0_WP-R
	REAL(WP) :: f1,f2,x0,x1,x2,x3

	x0=ax
	x3=bx
	x1=R*ax+C*bx
	x2=C*ax+R*bx
	f1=func(x1)
	f2=func(x2)
	DO
		IF (abs(x3-x0) <= TOL*(abs(x1)+abs(x2))) EXIT
		IF (f2 > f1) THEN
			x0=x1
			x1=x2
			x2=R*x2+C*x3
			f1=f2
			f2=func(x2)
		ELSE
			x3=x2
			x2=x1
			x1=R*x1+C*x0
			f2=f1
			f1=func(x1)
		END IF
	END DO
	IF (f1 > f2) THEN
		golden=f1
		xmax=x1
	ELSE
		golden=f2
		xmax=x2
	END IF
END FUNCTION golden

FUNCTION zbrent(func,x1,x2,tol)
	IMPLICIT NONE
	REAL(WP), INTENT(IN) :: x1,x2,tol
	REAL(WP) :: zbrent
	INTERFACE
		FUNCTION func(x)
		!USE mkl95_precision, ONLY: WP => DP
        USE f95_precision, ONLY: WP => DP
		IMPLICIT NONE
		REAL(WP), INTENT(IN) :: x
		REAL(WP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=100
	REAL(WP), PARAMETER :: EPS=epsilon(x1)
	INTEGER :: iter
	REAL(WP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)
	if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
		STOP 'root must be bracketed for zbrent'
	c=b
	fc=fb
	do iter=1,ITMAX
		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
			fc=fa
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol1=2.0_WP*EPS*abs(b)+0.5_WP*tol
		xm=0.5_WP*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			zbrent=b
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_WP*xm*s
				q=1.0_WP-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_WP*xm*q*(q-r)-(b-a)*(r-1.0_WP))
				q=(q-1.0_WP)*(r-1.0_WP)*(s-1.0_WP)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_WP*p  <  min(3.0_WP*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
		fb=func(b)
	end do
	STOP 'zbrent: exceeded maximum iterations'
	zbrent=b
END FUNCTION zbrent

END MODULE toolbox
