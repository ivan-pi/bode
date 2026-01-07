module diff_mod
use bode_mod, only: wp, ltri, tsol
implicit none
contains
    subroutine solve_banded(n,k,a,lda,ipiv,b,ifail)
        implicit none
        integer, intent(in) :: n, k, lda
        real(wp), intent(inout) :: a(lda,*), b(n)
        integer, intent(out) :: ipiv(n), ifail

        real(wp) :: x(n)
        real(wp) :: tl(n,k)
        integer :: m1, m21, m3
        integer :: i

        m1 = k       ! half-bandwidth
        m21 = 2*m1+1  ! total bandwidth
        m3 = m1    ! extra needed for factorization

!        print *, "m1  = ", m1
!        print *, "m21 = ", m21
!        print *, "m3  = ", m3
!        tl = -66

        call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
        if (ifail /= 0) return

!        print *, "FACTORIZATION"
!        do i = 1, n
!            print *, "tl(i,:) = ", tl(i,:), "|", ipiv(i)
!        end do

        ! Copy input values
        x(1:n) = b(1:n)
        call tsol(a,tl,lda, &
            m1=m1,m3=m3,m21=m21,n=n,&
            ipiv=ipiv,xin=x,xout=b)

    end subroutine

    ! Reaction-Diffusion Example 1
    ! Uses backward difference (first order) for the Neumann BC
    subroutine setup1(n, ab, ldab, b)
        integer, intent(in) :: n, ldab
        real(wp), intent(out) :: ab(ldab,3)
        real(wp), intent(out) :: b(n)

        real(wp) :: h, diff, reac
        integer :: i

        b = 0.0_wp
        b(1) = 1.0_wp ! Dirichlet boundary condition
        ! b(n) = 0.0_wp ! zero-flux boundary condition satisfied during initialization

        ! stepsize, diffusivity and reaction coefficient
        h = 1.0_wp/real(n-1,wp)
        diff = 1.0_wp
        reac = 1.0_wp

        ab = 0.0_wp ! set matrix to zero
        
        ! D [1 -2 1] u = h^2 k u
        do i = 2, n - 1
            ab(i,1:3) = [diff, -2._wp*diff - h**2*reac, diff]
        end do
        
        ! column 3 is the diagonal

        ab(1,2) = 1.0 ! Dirichlet boundary condition

        ! Zero-flux boundary condition on right side
        ab(n,1) = 1.0_wp
        ab(n,2) = -1.0_wp

    end subroutine

    ! Reaction-Diffusion Example 2
    ! Uses central difference (via ghost node) for the Neumann BC
    subroutine setup2(n, ab, ldab, b)
        integer, intent(in) :: n, ldab
        real(wp), intent(out) :: ab(ldab,3)
        real(wp), intent(out) :: b(n)

        real(wp) :: h, diff, reac
        integer :: i

        b = 0.0_wp
        b(1) = 1.0_wp ! Dirichlet boundary condition
        ! b(n) = 0.0_wp ! zero-flux boundary condition satisfied during initialization

        ! stepsize, diffusivity and reaction coefficient
        h = 1.0_wp/real(n-1,wp)
        diff = 1.0_wp
        reac = 1.0_wp

        ab = 0.0_wp ! set matrix to zero
        
        ! D [1 -2 1] u = h^2 k u
        do i = 2, n - 1
            ab(i,1:3) = [diff, -2._wp*diff - h**2*reac, diff]
        end do

        ab(1,2) = 1.0 ! Dirichlet boundary condition

        ! Zero-flux boundary condition on right side
        ! using ghost-node approach
        ab(n,1) =  2.0_wp*diff
        ab(n,2) = -2.0_wp*diff - h**2*reac

    end subroutine

    ! Reaction-Diffusion Example 3
    ! Uses backward difference (second order) for the Neumann BC
    subroutine setup3(n, ab, ldab, b)
        integer, intent(in) :: n, ldab
        real(wp), intent(out) :: ab(ldab,-2:2)
        real(wp), intent(out) :: b(n)

        real(wp) :: h, diff, reac

        ab = 0.0_wp
        b = 0.0_wp

        ! stepsize, diffusivity and reaction coefficient
        h = 1.0_wp/real(n-1,wp)
        diff = 1.0_wp
        reac = 1.0_wp

        ! Fill tridiagonal parts for all internal nodes
        ab(1:n, -1) = diff
        ab(1:n,  0) = -2.0_wp * diff - h**2 * reac
        ab(1:n,  1) = diff

        ! Note: the diagonals -2 and 2 are set above

        ! Row 1: Dirichlet boundary
        ab(1,0:2) = [1.0_wp, 0.0_wp, 0.0_wp]
        b(1)     = 1.0_wp

        ! Row n: Neumann boundary, backward differences (second order)
        ab(n, -2:0) = [1.0_wp, -4.0_wp, 3.0_wp]

    end subroutine

    ! Analytical solution
    elemental function sol(x,Th) result(c)
        real(wp), intent(in) :: x, Th
        real(wp) :: c
        c = cosh((1.0 - x)*Th)/cosh(Th)
    end function

end module

program test_banded2
    use bode_mod, only: wp
    use diff_mod, only: setup1, setup2, setup3, solve_banded, sol
    implicit none

    ! External procedure inclusion
    include 'bandfs.fi'

    ! --- Configuration ---
    integer :: n, icase
    logical, parameter :: verbose = .false.

    real(wp), allocatable :: ab(:,:), b(:), x(:), c(:), res(:)
    real(wp), allocatable :: ab_orig(:,:), b_orig(:)
    integer, allocatable  :: ipiv(:)

    real(wp) :: h, l2_err
    integer  :: kh, ldab, info, i
    character(len=32) :: arg

    ! --- Command Line Parsing ---
    if (command_argument_count() /= 2) then
        print *, "Usage: ./test_banded2 <icase> <n>"
        print *, "  icase : 1 (1st order), 2 (Central/Ghost), 3 (2nd order backward)"
        print *, "  n     : Number of nodes"
        stop 1
    end if

    ! Get Case (Argument 1)
    call get_command_argument(1, arg)
    read(arg, *) icase

    ! Get N (Argument 2)
    call get_command_argument(2, arg)
    read(arg, *) n

    ! --- Status Message ---
    write(*, '(">>> Running Case ", I0, " with ", I0, " nodes.")') icase, n
    select case(icase)
    case(1)
        write(*,'(4X,A)') "Mode: 1st Order Backward Difference Neumann BC"
    case(2)
        write(*,'(4X,A)') "Mode: 2nd Order Central Difference (Ghost Node) Neumann BC"
    case(3)
        write(*,'(4X,A)') "Mode: 2nd Order Backward Difference Neumann BC"
    case default
        error stop "Invalid case selection. Choose 1, 2, or 3."
    end select

    ldab = n
    allocate(b(n), x(n), c(n), res(n), ipiv(n), b_orig(n))

    ! --- Setup Phase ---
    select case(icase)
    case(1)
        kh = 1
        allocate(ab(ldab,-kh:kh))
        call setup1(n, ab, ldab, b)
    case(2)
        kh = 1
        allocate(ab(ldab,-kh:kh))
        call setup2(n, ab, ldab, b)
    case(3)
        kh = 2
        allocate(ab(ldab,-kh:kh))
        call setup3(n, ab, ldab, b)
    end select

    ! Backup for validation (residuals)
    allocate(ab_orig, source=ab)
    b_orig = b

    if (verbose) call print_debug()

    ! --- Linear Solve ---
    call solve_banded(n, kh, ab, ldab, ipiv, b, info)

    ! Early exit on solver failure
    if (info /= 0) then
        write(*, '(A, I0)') "Error: solve_banded failed with info = ", info
        stop 1
    end if

    ! --- Post-Processing ---

    ! 1. Calculate Residuals: res = b_orig - 1.0 * A * b_sol
    res = b_orig
    call bndmv(n, kh, kh, -1.0_wp, ab_orig, ldab, b, 1.0_wp, res)

    ! 2. Calculate L2 Error Norm
    h = 1.0_wp / real(n-1, wp)
    x = [((i-1)*h, i=1, n)]
    c = sol(x, Th=1.0_wp)

    ! L2 Error: ||numerical - analytical||
    l2_err = norm2(b - c)

    print *, "-----------------------------------------------------------"
    print *, "FINAL SOLUTION AND RESIDUALS"
    print *, "   x      Analytical      Numerical        Residual"
    do i = 1, n
        print "(F6.3, 2F16.10, E16.6)", x(i), c(i), b(i), res(i)
    end do
    print *, "-----------------------------------------------------------"
    print "(A, E12.6)", " L2 Error Norm: ", l2_err
    print "(A, E12.6)", " Residual Norm: ", norm2(res)

contains

    subroutine print_debug()
        integer :: j
        print *, "--- Matrix Setup (Case", icase, ") ---"
        do j = 1, n
            print "(I3, ' |', 5F9.4, ' |', F9.4)", j, ab(j, :), b(j)
        end do
    end subroutine

end program test_banded2
