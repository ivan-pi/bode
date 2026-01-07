program bandfs_fem_1d_quadratic
  implicit none

  integer, parameter :: wp = kind(1.0d0)

  ! BANDFS procedure interfaces
  include "bandfs.fi"

  ! Simulation parameters
  integer :: n, lda, nel = 10, info, i, argc, node_idx
  integer, parameter :: m = 5       ! Bandwidth (Pentadiagonal for quadratic)
  integer, parameter :: h = (m-1)/2   ! Half-bandwidth (h=2)
  real(wp) :: phi = 1.0_wp, L, dx, Le
  character(len=32) :: arg

  ! Arrays
  real(wp), allocatable :: a(:,:), b(:)
  integer,  allocatable :: ipiv(:)

  real(wp) :: ke(3,3), me(3,3)
  integer  :: e

  ! 1. Simple Command Line Parser
  argc = command_argument_count()
  if (argc >= 1) then
      call get_command_argument(1, arg)
      read(arg, *) nel  ! Input now refers to number of elements
  end if
  if (argc >= 2) then
      call get_command_argument(2, arg)
      read(arg, *) phi
  end if

  ! 2. Initialization for Quadratic Elements
  n   = 2 * nel + 1        ! Total nodes for nel quadratic elements
  L   = 1.0_wp
  Le  = L / real(nel, wp)  ! Length of one quadratic element
  dx  = Le / 2.0_wp        ! Distance between adjacent nodes

  lda = n
  ! Banded storage: a(lda, -h:h+h)
  allocate(a(lda, -h:h+h), b(n), ipiv(n))
  a = 0.0_wp
  b = 0.0_wp

  ! 3. Local Element Matrices (Quadratic Lagrange Basis)
  ! Stiffness: integral(N_i' * N_j') over [0, Le]
  ke = reshape([ 7.0_wp, -8.0_wp,  1.0_wp, &
                -8.0_wp, 16.0_wp, -8.0_wp, &
                 1.0_wp, -8.0_wp,  7.0_wp  ], [3,3], order=[2,1]) / (3.0_wp * Le)

  ! Mass: integral(N_i * N_j) * phi^2 over [0, Le]
  me = reshape([ 4.0_wp,  2.0_wp, -1.0_wp, &
                 2.0_wp, 16.0_wp,  2.0_wp, &
                -1.0_wp,  2.0_wp,  4.0_wp  ], [3,3], order=[2,1]) * (Le / 30.0_wp) * (phi**2)

  ! 4. Global Assembly Loop
  do e = 1, nel
      node_idx = 2*e - 1  ! Left-most node of element e
      
      ! Add contributions to the banded matrix
      ! Row node_idx (i)
      a(node_idx,   0) = a(node_idx,   0) + (ke(1,1) + me(1,1))
      a(node_idx,   1) = a(node_idx,   1) + (ke(1,2) + me(1,2))
      a(node_idx,   2) = a(node_idx,   2) + (ke(1,3) + me(1,3))
      
      ! Row node_idx + 1 (i+1, Midpoint)
      a(node_idx+1, -1) = a(node_idx+1, -1) + (ke(2,1) + me(2,1))
      a(node_idx+1,  0) = a(node_idx+1,  0) + (ke(2,2) + me(2,2))
      a(node_idx+1,  1) = a(node_idx+1,  1) + (ke(2,3) + me(2,3))
      
      ! Row node_idx + 2 (i+2)
      a(node_idx+2, -2) = a(node_idx+2, -2) + (ke(3,1) + me(3,1))
      a(node_idx+2, -1) = a(node_idx+2, -1) + (ke(3,2) + me(3,2))
      a(node_idx+2,  0) = a(node_idx+2,  0) + (ke(3,3) + me(3,3))
  end do

  ! 5. Apply Boundary Conditions
  ! Node 1: Natural (Neumann), du/dx = 0; (Symmetry) is handled naturally.

  ! Node n: Dirichlet, u(L) = 1.0
  a(n, -2:-1) = 0.0_wp  ! Zero out lower bands for the last row
  a(n,  0)    = 1.0_wp  ! Diagonal to 1
  b(n)        = 1.0_wp  ! RHS to target value

  ! 6. Solver Stage
  call bandf(n, m, a, lda, ipiv, info)
  if (info /= 0) stop "Factorization error: Matrix is singular."

  call bands(n, m, a, lda, ipiv, b, info)

  ! 7. Results
  write(*,'("# Quadratic FEM: nel=", I0, ", n=", I0, ", phi=", F8.3)') nel, n, phi
  write(*,'("#", A9, A12)') "x", "u(x)"
  write(*,'(F10.4, F12.6)') ((i-1)*dx, b(i), i=1,n)

end program bandfs_fem_1d_quadratic
