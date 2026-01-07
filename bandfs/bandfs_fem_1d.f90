program bandfs_fem_1d
  implicit none

  integer, parameter :: wp = kind(1.0d0)

  ! BANDFS procedure interfaces
  include "bandfs.fi"

  ! Simulation parameters
  integer :: n = 50, lda, nel, info, i, argc
  integer, parameter :: m = 3           ! Bandwidth (tridiagonal)
  integer, parameter :: h = (m-1)/2     ! Half-bandwidth
  real(wp) :: phi = 1.0_wp, L, dx
  character(len=32) :: arg

  ! Arrays (Allocatable for dynamic sizing)
  real(wp), allocatable :: a(:,:), b(:)
  integer,  allocatable :: ipiv(:)

  real(wp) :: ke(2,2), me(2,2)
  integer  :: e

  ! 1. Simple Command Line Parser
  argc = command_argument_count()
  if (argc >= 1) then
      call get_command_argument(1, arg)
      read(arg, *) n
  end if
  if (argc >= 2) then
      call get_command_argument(2, arg)
      read(arg, *) phi
  end if

  ! 2. Initialization
  nel = n - 1
  L   = 1.0_wp
  dx  = L / real(nel, wp)

  lda = n
  allocate(a(lda, -h:h+h), b(n), ipiv(n))
  a = 0.0_wp
  b = 0.0_wp

  ! Local Element Matrices
  ! Stiffness: integral(u' * v')
  ke = reshape([1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp], [2,2]) / dx
  ! Mass: integral(u * v) * phi^2
  me = reshape([2.0_wp, 1.0_wp, 1.0_wp, 2.0_wp], [2,2]) * (dx / 6.0_wp) * (phi**2)

  ! 3. Global Assembly Loop
  do e = 1, nel
     ! Update Banded Matrix Diagonals using row slicing
     a(e:e+1,  0) = a(e:e+1,  0) + [ (ke(1,1) + me(1,1)), (ke(2,2) + me(2,2)) ]
     a(e,      1) = a(e,      1) + (ke(1,2) + me(1,2))
     a(e+1,   -1) = a(e+1,   -1) + (ke(2,1) + me(2,1))
  end do

  ! 4. Apply Boundary Conditions

  ! Node 1: Natural (Neumann) Boundary Condition (Symmetry)
  ! We "do nothing" here. The integration by parts naturally enforces du/dx = 0.

  ! Node n: Fixed Concentration - Dirichlet Boundary Condition.
  ! Set the n-th row to represent u_n = 1.0
  a(n, :-1)  = 0.0_wp  ! Zero the lower diagonal entry
  a(n,  0)  = 1.0_wp  ! Identity on main diagonal
  b(n)      = 1.0_wp  ! Target value

  ! 5. Solver Stage
  call bandf(n, m, a, lda, ipiv, info)
  if (info /= 0) stop "Factorization error: Matrix is singular."

  call bands(n, m, a, lda, ipiv, b, info)

  ! 6. Results - Using implied-DO loop and format repetition
  write(*,'("# n=", I0, ", phi=", F8.3)') n, phi
  write(*,'("#", A9, A12)') "x", "u(x)"
  write(*,'(F10.4, F12.6)') ((i-1)*dx, b(i), i = 1, n)

end program bandfs_fem_1d
