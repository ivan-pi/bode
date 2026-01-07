program bandfs_reaction_diffusion
    implicit none
    integer, parameter :: wp = kind(1.0d0)

    ! BANDFS procedure interfaces
    include "bandfs.fi"

    ! Matrix storage: m=3 (bandwidth), h=1 (half-bandwidth)
    integer, parameter :: m = 3, h = (m-1)/2
    real(wp), allocatable :: a(:,:), b(:)
    integer, allocatable  :: ipiv(:)

    ! Simulation parameters
    integer :: n = 50, lda, info, i, argc
    real(wp) :: phi = 1.0_wp, dx
    character(len=32) :: arg

    ! 1. Simple Command Line Parser
    ! Usage: ./solve <n_points> <thiele_modulus_phi>
    argc = command_argument_count()
    if (argc >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) n
    end if
    if (argc >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) phi
    end if

    ! 2. Allocation
    ! a(n, -1:2) provides space for tridiagonal + pivoting fill-in
    allocate(a(n, -h:h+h), b(n), ipiv(n))
    lda = n
    
    dx = 1.0_wp / real(n-1, wp)
    a = 0.0_wp
    b = 0.0_wp

    ! 3. Setup SPD Matrix using array slicing
    ! Discretization: -[1, -2, 1] + (dx*phi)**2 * [0, 1, 0]
    a(:, -1) = -1.0_wp                           ! Lower diagonal
    a(:,  0) =  2.0_wp + (dx * phi)**2           ! Main diagonal
    a(:,  1) = -1.0_wp                           ! Upper diagonal

    ! 4. Apply Boundary Conditions
    ! Left Boundary (index 1): du/dx = 0 (Symmetry at center)
    ! Ghost node approach: u_0 = u_2. 
    ! Standard row: -u_0 + (2 + dx^2*phi^2)u_1 - u_2 = 0
    ! Becomes: (2 + dx^2*phi^2)u_1 - 2*u_2 = 0
    a(1,  0) =  2.0_wp + (dx * phi)**2
    a(1,  1) = -2.0_wp
    b(1)     =  0.0_wp

    ! Right Boundary (index n): u(1) = 1.0 (Dirichlet at edge)
    a(n, -1) =  0.0_wp
    a(n,  0) =  1.0_wp
    b(n)     =  1.0_wp

    ! 5. LU Factorization and Solve
    call bandf(n, m, a, lda, ipiv, info)
    if (info /= 0) stop "Error: Factorization failed."

    ! Solution vector x is returned in b
    call bands(n, m, a, lda, ipiv, b, info)

    ! 6. Output Results
    ! Uses an implied-DO list to print all pairs. 
    ! The last format repeats for each pair until the data is exhausted.
    write(*,'("# n=", I0, ", phi=", F8.3)') n, phi
    write(*,'("#", A9, A12)') "x", "u(x)"
    write(*,'(F10.4, F12.6)') ((i-1)*dx, b(i), i=1,n)

end program
