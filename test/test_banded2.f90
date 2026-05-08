module diff_mod
    use bode_mod, only: wp, ltri, tsol
    implicit none
contains
    subroutine solve_banded(n,k,a,lda,ipiv,b,ifail)
        integer, intent(in) :: n, k, lda
        real(wp), intent(inout) :: a(lda,*), b(n)
        integer, intent(out) :: ipiv(n), ifail

        real(wp) :: x(n)
        real(wp) :: tl(n,k)
        integer :: m1, m21, m3

        m1 = k
        m21 = 2*m1+1
        m3 = m1

        call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
        if (ifail /= 0) return

        x(1:n) = b(1:n)
        call tsol(a,tl,lda,m1=m1,m3=m3,m21=m21,n=n,ipiv=ipiv,xin=x,xout=b)
    end subroutine

    subroutine setup1(n, ab, ldab, b)
        integer, intent(in) :: n, ldab
        real(wp), intent(out) :: ab(ldab,3), b(n)
        real(wp) :: h, diff, reac
        integer :: i

        b = 0.0_wp
        b(1) = 1.0_wp
        h = 1.0_wp/real(n-1,wp)
        diff = 1.0_wp
        reac = 1.0_wp

        ab = 0.0_wp
        do i = 2, n - 1
            ab(i,1:3) = [diff, -2._wp*diff - h**2*reac, diff]
        end do
        ab(1,2) = 1.0_wp
        ab(n,1) = 1.0_wp
        ab(n,2) = -1.0_wp
    end subroutine

    subroutine setup2(n, ab, ldab, b)
        integer, intent(in) :: n, ldab
        real(wp), intent(out) :: ab(ldab,3), b(n)
        real(wp) :: h, diff, reac
        integer :: i

        b = 0.0_wp
        b(1) = 1.0_wp
        h = 1.0_wp/real(n-1,wp)
        diff = 1.0_wp
        reac = 1.0_wp

        ab = 0.0_wp
        do i = 2, n - 1
            ab(i,1:3) = [diff, -2._wp*diff - h**2*reac, diff]
        end do
        ab(1,2) = 1.0_wp
        ab(n,1) =  2.0_wp*diff
        ab(n,2) = -2.0_wp*diff - h**2*reac
    end subroutine

    subroutine setup3(n, ab, ldab, b)
        integer, intent(in) :: n, ldab
        real(wp), intent(out) :: ab(ldab,-2:2), b(n)
        real(wp) :: h, diff, reac

        ab = 0.0_wp
        b = 0.0_wp

        h = 1.0_wp/real(n-1,wp)
        diff = 1.0_wp
        reac = 1.0_wp

        ab(1:n, -1) = diff
        ab(1:n,  0) = -2.0_wp * diff - h**2 * reac
        ab(1:n,  1) = diff

        ab(1,0:2) = [1.0_wp, 0.0_wp, 0.0_wp]
        b(1)     = 1.0_wp

        ab(n, -2:0) = [1.0_wp, -4.0_wp, 3.0_wp]
    end subroutine

    elemental function sol(x,theta) result(c)
        real(wp), intent(in) :: x, theta
        real(wp) :: c
        c = cosh((1.0_wp - x)*theta)/cosh(theta)
    end function

end module

program test_banded2
    use bode_mod, only: wp
    use diff_mod, only: setup1, setup2, setup3, solve_banded, sol
    use test_support, only: assert_true, assert_close_scalar, finalize_tests
    implicit none

    include 'bandfs.fi'

    character(len=256) :: ref_file

    if (command_argument_count() /= 1) then
        error stop "Usage: test_banded2 <reference_file>"
    end if

    call get_command_argument(1, ref_file)

    call run_case(1, 21, ref_file)
    call run_case(2, 21, ref_file)
    call run_case(3, 21, ref_file)

    call finalize_tests()

contains

    subroutine run_case(icase, n, ref_file)
        integer, intent(in) :: icase, n
        character(len=*), intent(in) :: ref_file

        real(wp), allocatable :: ab(:,:), b(:), x(:), c(:), res(:)
        real(wp), allocatable :: ab_orig(:,:), b_orig(:)
        integer, allocatable  :: ipiv(:)

        real(wp) :: h, l2_err, res_norm, ref_l2, ref_res
        integer  :: kh, ldab, info, i

        ldab = n
        allocate(b(n), x(n), c(n), res(n), ipiv(n), b_orig(n))

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
        case default
            error stop "invalid case"
        end select

        allocate(ab_orig, source=ab)
        b_orig = b

        call solve_banded(n, kh, ab, ldab, ipiv, b, info)
        call assert_true("solve_banded_case"//itoa(icase), info == 0)

        res = b_orig
        call bndmv(n, kh, kh, -1.0_wp, ab_orig, ldab, b, 1.0_wp, res)

        h = 1.0_wp / real(n-1, wp)
        x = [((i-1)*h, i=1, n)]
        c = sol(x, theta=1.0_wp)

        l2_err = norm2(b - c)
        res_norm = norm2(res)

        call read_reference(ref_file, icase, n, ref_l2, ref_res)

        call assert_close_scalar("l2_err_case"//itoa(icase), l2_err, ref_l2, 2.5e-4_wp)
        call assert_close_scalar("res_norm_case"//itoa(icase), res_norm, ref_res, 1.0e-12_wp)

        deallocate(ab, ab_orig, b, b_orig, x, c, res, ipiv)
    end subroutine

    subroutine read_reference(ref_file, icase, n, ref_l2, ref_res)
        character(len=*), intent(in) :: ref_file
        integer, intent(in) :: icase, n
        real(wp), intent(out) :: ref_l2, ref_res

        integer :: unit, ios, line_case, line_n
        real(wp) :: l2, res
        character(len=256) :: line
        logical :: found

        found = .false.
        open(newunit=unit, file=trim(ref_file), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "Could not open reference file"

        read(unit,'(A)',iostat=ios) line
        if (ios /= 0) then
            close(unit)
            error stop "Invalid reference header"
        end if
        do
            read(unit,*,iostat=ios) line_case, line_n, l2, res
            if (ios /= 0) exit
            if (line_case == icase .and. line_n == n) then
                ref_l2 = l2
                ref_res = res
                found = .true.
                exit
            end if
        end do
        close(unit)

        if (.not. found) then
            error stop "Reference row not found"
        end if
    end subroutine

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=16) :: s
        write(s,'(I0)') i
    end function

end program
