program quadtest
  use polynomial
  
  implicit none

  integer :: j, i
  type(poly_gauss) :: poly_data

  ! Test polynomial functions/subroutines
  do j = 1, 10
     poly_data%n = j
     allocate(poly_data%coefs(poly_data%n+1))
     allocate(poly_data%roots(poly_data%n))
     do i = 1, poly_data%n
        poly_data%roots(i) = dble(i)
     end do
     
     ! See if I get the right coefficients
     call poly_root_to_coef(poly_data%roots, poly_data%coefs)

     write(*,*) '----------------------------------------'
     write(*,*) 'Polynomial order: ', poly_data%n
     write(*,*) 'The roots are: '
     do i = 1, poly_data%n
        write(*,*) poly_data%roots(i)
     end do
     write(*,*) 'The coefficients c_k for x^k are: '
     do i = 1, poly_data%n+1
        write(*,*) 'k, ck = ', i-1, ', ', poly_data%coefs(i)
     end do

     deallocate(poly_data%coefs)
     deallocate(poly_data%roots)
  end do
end program quadtest
