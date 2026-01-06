program bandfs_demo
  implicit none

  integer, parameter :: wp = kind(1.0d0) ! Double precision

  ! BANDFS procedure interfaces
  include "bandfs.fi"

  integer, parameter :: n = 5, m = 3  ! Dimension and bandwidth
  integer, parameter :: h = (m-1)/2   ! Half-bandwidth
  
  ! Array storage
  integer, parameter :: lda = n
  real(wp) :: a(lda,-h:h+h), b(n)
  integer  :: ipiv(n), info


  ! 1. Initialize array and RHS
  ! Gilbert Strang's favorite matrix (second order differences)
  a(:,-1) = [ 0,-1,-1,-1,-1]  ! Lower diagonal
  a(:, 0) = [ 2, 2, 2, 2, 2]  ! Main diagonal
  a(:, 1) = [-1,-1,-1,-1, 0]  ! Upper diagonal

  ! Right-hand side for a target solution of x = [1, 2, 3, 4, 5]
  b = [ 0, 0, 0, 0, 6]

  ! 2. Perform LU Factorization
  ! 'a' is modified in-place; 'ipiv' stores pivoting sequence
  call bandf(n, m, a, lda, ipiv, info)

  if (info /= 0) then
     write(*,'(A)') "Error: Matrix is singular or factorization failed."
     stop 1
  end if

  ! 3. Solve the system
  ! 'b' is overwritten with the solution vector 'x'
  call bands(n, m, a, lda, ipiv, b)

  write(*,'(A,/,*(2X,F6.3,:,/))') "Solution x:", b

end program
