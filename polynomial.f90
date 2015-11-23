module polynomial
  
  implicit none
  
  ! Generate quadrature coefficients and points for various types of numerical integration
  double precision, parameter :: PI = 3.14159265358979d0

  type :: poly_coef
     double precision, allocatable :: coefs(:)
     integer :: n
  end type poly_coef

  type, extends(poly_coef) :: poly_gauss
     double precision, allocatable :: roots(:)
     double precision, allocatable :: weights(:)     
  end type poly_gauss

  type(poly_gauss) :: legendre_data, laguerre_data
  
contains
  
  function factorial(n) result(fn)
    integer :: n, fn, i
    if (n<0) then
       write(*,*) 'Error, attempting factorial of negative number!'
       stop
    else if (n==0) then
       fn = 1
    else
       fn = 1
       do i = n, 1, -1
          fn = fn*i
       end do
    end if
    return
  end function factorial

  function double_factorial(n) result(fn)
    integer :: n, fn, i
    if (n<-1) then
       write(*,*) 'Error, attempting factorial of negative number less than -1!'
       stop
    else if (n==0 .or. n==-1) then
       fn = 1
    else if (mod(n,2)==0) then
       fn = 1
       do i = n, 2, -2
          fn = fn*i
       end do
    else if (mod(n,2)==1) then
       fn = 1
       do i = n, 1, -2
          fn = fn*i
       end do
    end if
    return
  end function double_factorial

  function half_integer_factorial(n) result(fn)
    integer :: n ! Return the factorial of n/2 for n > 0
    double precision :: fn
    if (mod(n,2)==0) then
       ! n is even
       fn = factorial(n/2)
    else
       ! n is odd
       fn = dble(double_factorial(n-2))*sqrt(PI)/dble(2**((n-1)/2))
    end if
    return
  end function half_integer_factorial

  function binomial_n_pick_k(n,k) result(bnk)
    double precision :: bnk
    integer :: n, k
    if (n < k) then
       write(*,*) 'Error, attempting binomial n pick k with n < k!'
       stop
    else
       bnk = dble(factorial(n))/dble(factorial(k)*factorial(n-k))
    end if
    return
  end function binomial_n_pick_k

  function binomial_r_pick_k(r,k) result(brk)
    double precision :: r, brk
    integer :: k, i
    brk = 1.0d0
    do i = 0, k-1
       brk = brk * (r - dble(i))
    end do
    brk = brk/dble(factorial(k))
    return
  end function binomial_r_pick_k

  function legendre_coef(n,k) result(lcoef)
    double precision :: lcoef
    integer :: n, k
    if (k < 0) then
       write(*,*) 'Error: attempting to calculate coefficient k < 0!'
    else if (k > n) then
       write(*,*) 'Error: attempting to calculate coefficient k > n!'
    else
       lcoef = dble(2**n) * binomial_n_pick_k(n, k) * binomial_r_pick_k(dble(n+k-1)/2.0d0, n)
    end if
    return
  end function legendre_coef

  function laguerre_coef(n,k) result(lcoef)
    double precision :: lcoef
    integer :: n, k, m
    if (k < 0) then
       write(*,*) 'Error: attempting to calculate coefficient k < 0!'
    else if (k > n) then
       write(*,*) 'Error: attempting to calculate coefficient k > n!'
    else
       if (mod(k,2)==0) then
          m = 1
       else
          m = -1
       end if
       lcoef = dble(m) * binomial_n_pick_k(n, k) / dble(factorial(k))
    end if
    return
  end function laguerre_coef

  subroutine poly_gauss_clr(poly_r)
    type(poly_gauss), intent(inout) :: poly_r
    poly_r%n = 0
    deallocate(poly_r%coefs)
    deallocate(poly_r%roots)
    return
  end subroutine poly_gauss_clr

  subroutine legendre_form(n)
    integer, intent(in) :: n
    integer :: i
    double precision :: lcoef_n

    if (n < 0) then
       write(*,*) 'Error: attempting to form legendre polynomial of degree < 0!'
    end if
    
    legendre_data%n = n
    allocate(legendre_data%coefs(legendre_data%n+1))
    allocate(legendre_data%roots(legendre_data%n))
    
    lcoef_n = legendre_coef(n,n)
    do i = 0, n
       if (i == n) then
          legendre_data%coefs(i+1) = 1.0d0
       else
          legendre_data%coefs(i+1) = legendre_coef(n,i)/lcoef_n
       end if
    end do
  end subroutine legendre_form

  subroutine laguerre_form(n)
    integer, intent(in) :: n
    integer :: i
    double precision :: lcoef_n

    if (n < 0) then
       write(*,*) 'Error: attempting to form laguerre polynomial of degree < 0!'
    end if
    
    laguerre_data%n = n
    allocate(laguerre_data%coefs(laguerre_data%n+1))
    allocate(laguerre_data%roots(laguerre_data%n))    
    
    lcoef_n = laguerre_coef(n,n)
    do i = 0, n
       if (i == n) then
          laguerre_data%coefs(i+1) = 1.0d0
       else
          laguerre_data%coefs(i+1) = laguerre_coef(n,i)/lcoef_n
       end if
    end do
  end subroutine laguerre_form
  
  subroutine poly_get_coefs(m, k, roots, coefs_m_k)
    double precision, dimension(:), intent(in) :: roots
    double precision, dimension(:,:), intent(inout) :: coefs_m_k
    double precision :: cmk
    integer :: k, m
    integer :: i
    
    ! Fills in the k-th coefficient for the m-th construction polynomial
    cmk = 0.0d0 ! Ignore the case m < k and don't do anything to save evaluation
    if ( k == m ) then ! m = k
       cmk = 1.0d0
    else if ( k < m ) then ! m > k
       if ( k == 0 ) then
          cmk = 1.0d0
          do i = 1, m
             cmk = cmk * (-roots(i))
          end do
       else if (k == m - 1) then
          cmk = 0.0d0
          do i = 1, m
             cmk = cmk - roots(i)
          end do
       else
          cmk = -roots(m) * coefs_m_k(m, k+1) + coefs_m_k(m, k)
       end if
    end if
    coefs_m_k(m+1, k+1) = cmk
!    write(*,*) 'Got m = ',m,', k = ',k,', and cmk = ',cmk
    return
  end subroutine poly_get_coefs
  
  subroutine poly_root_to_coef(roots, coefs)
    double precision, dimension(:), intent(in)  :: roots
    double precision, dimension(:), intent(out) :: coefs
    double precision, allocatable :: coefs_m_k(:,:)
    integer :: mi, ki, nr, nc, n

    nr = size(roots)
    nc = size(coefs)
    n = nc - 1

    allocate( coefs_m_k(n+1, n+1) )
    do mi = 0, n
       do ki = 0, mi
          call poly_get_coefs(mi, ki, roots, coefs_m_k)
       end do
    end do

    coefs(:) = coefs_m_k(n+1,:)
    deallocate( coefs_m_k )
    return
  end subroutine poly_root_to_coef
  
  function poly_eval(x,poly_d) result(f)
    ! Returns the polynomial poly_d at x
    type(poly_coef) :: poly_d
    double precision :: x, f
    integer :: i

    f = 0.0d0
    do i = 0, poly_d%n
       f = f + (x**i) * poly_d%coefs(i+1)
    end do
    return
  end function poly_eval
  
  function poly_deval(x, m, poly_d) result(f)
    ! Returns the m-th order derivative of poly_d at x
    type(poly_coef) :: poly_d
    double precision :: x, f
    integer :: i, j, m, fi
    
    f = 0.0d0
    if (poly_d%n == 0 .or. m >= poly_d%n+1) then
       return
    else if (m == 0) then
       f = poly_eval(x, poly_d)
       return
    else
       do i = m, poly_d%n
          fi = 1
          do j = 0, m-1
             fi = fi * (i - j)
          end do
          f = f + dble(fi) * (x**(i-m)) * poly_d%coefs(i+1)
       end do
    end if
    return
  end function poly_deval

  subroutine gauss_int_print(poly_g, save_file_name)
    ! Output Profile Misc
    character(len=*), intent(in) :: save_file_name
    type(poly_gauss), intent(in) :: poly_g
    
    character(len=50) :: rfmt
    character(len=50) :: hfmt
    integer, parameter :: ncols = 2
    integer :: k

    ! Set format for writing column entries
    write(hfmt,'(A,I5,A)') '(', ncols, '(1x,A25))'
    write(rfmt,'(A,I5,A)') '(', ncols, '(1x,ES25.14))'
    open(unit=2, file=save_file_name, action='write', &
         status='replace', recl=(25*ncols+10))
    write(2, fmt=hfmt) 'Abscissa', 'Weights'
    do k = 1, poly_g%n+1
       write(2, fmt=rfmt) poly_g%roots(k), poly_g%weights(k)
    end do
    close(unit=2)
  end subroutine gauss_int_print

  ! subroutine legendre_roots(n)
  !   integer, intent(in) :: n
  !   integer :: i, j, npts
  !   double precision :: xo, f, dfdx, xn, xxo, xxn
  !   double precision, parameter :: a = -1.0d0
    
  !   if (n == 0) then
  !      write(*,*) 'Error: attempted to get nonexistent abscissa for Legendre degree 0!'
  !      stop
  !   end if
  !   call legendre_form(n)

  !   npts = 1000*n
  !   xo = a
  !   legendre_data%na_found = 0
    
  !   do i = 1, npts
  !      f = legendre_eval(xo)
  !      dfdx = legendre_deval(xo)
  !      xn = xo - 0.1d0*f/dfdx ! Use a tiny step to avoid overshooting a root
  !      if (abs((xn-xo)/xn) < RTOL) then
  !         write(*,*) 'Legendre root found! x = ', xn
  !         legendre_data%na_found = legendre_data%na_found + 1
  !         legendre_data%roots(legendre_data%na_found) = xn
  !         if (legendre_data%na_found == legendre_data%n) then
  !            write(*,*) 'All Legendre roots found!'
  !            return
  !         else
  !            ! Get over the next hump
  !            xxo = xo
  !            do j = 1, npts
  !               f = legendre_deval(xxo)
  !               dfdx = legendre_ddeval(xxo)
  !               xxn = xxo + abs(1.1d0*f/dfdx)
  !               if (legendre_deval(xxn)/f < 0.0d0) then
  !                  ! sign of derivative changed, seek another root
  !                  exit
  !               else
  !                  cycle
  !               end if
  !               if (j == npts) then
  !                  write(*,*) 'Error: sign of Legendre derivative could not be changed!'
  !               end if
  !            end do
  !         end if
  !      else
  !         xo = xn
  !      end if
  !      if (i == npts) then
  !         write(*,*) 'Error: Not all Legendre roots could be found!'
  !      end if
  !   end do
    
  !   return
  ! end subroutine legendre_roots  

  ! subroutine legendre_clr()
  !   legendre_data%n = 0
  !   deallocate(legendre_data%coefs)
  !   deallocate(legendre_data%roots)    
  ! end subroutine legendre_clr

  ! subroutine laguerre_roots(n)
  !   integer, intent(in) :: n
  !   integer :: i, j, npts
  !   double precision :: xo, f, dfdx, xn, xxo, xxn
  !   double precision, parameter :: a = 1.0d0
    
  !   if (n == 0) then
  !      write(*,*) 'Error: attempted to get nonexistent abscissa for Laguerre degree 0!'
  !      stop
  !   end if
  !   call laguerre_form(n)

  !   npts = 1000*n
  !   xo = a
  !   laguerre_data%na_found = 0
    
  !   do i = 1, npts
  !      f = laguerre_eval(xo)
  !      dfdx = laguerre_deval(xo)
  !      xn = xo - 0.1d0*f/dfdx ! Use a tiny step to avoid overshooting a root
  !      if (abs((xn-xo)/xn) < RTOL) then
  !         write(*,*) 'Laguerre root found! x = ', xn
  !         laguerre_data%na_found = laguerre_data%na_found + 1
  !         laguerre_data%roots(laguerre_data%na_found) = xn
  !         if (laguerre_data%na_found == laguerre_data%n) then
  !            write(*,*) 'All Laguerre roots found!'
  !            return
  !         else
  !            ! Get over the next hump
  !            xxo = xo
  !            do j = 1, npts
  !               f = laguerre_deval(xxo)
  !               dfdx = laguerre_ddeval(xxo)
  !               xxn = xxo + abs(1.1d0*f/dfdx)
  !               if (laguerre_deval(xxn)/f < 0.0d0) then
  !                  ! sign of derivative changed, seek another root
  !                  exit
  !               else
  !                  cycle
  !               end if
  !               if (j == npts) then
  !                  write(*,*) 'Error: sign of Laguerre derivative could not be changed!'
  !               end if
  !            end do
  !         end if
  !      else
  !         xo = xn
  !      end if
  !      if (i == npts) then
  !         write(*,*) 'Error: Not all Laguerre roots could be found!'
  !      end if
  !   end do
    
  !   return
  ! end subroutine laguerre_roots
    
  ! subroutine laguerre_clr()
  !   laguerre_data%n = 0
  !   deallocate(laguerre_data%coefs)
  !   deallocate(laguerre_data%roots)    
  ! end subroutine laguerre_clr
  
end module polynomial
