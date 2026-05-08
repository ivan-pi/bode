program test_banded1
    use bode_mod, only: wp, ltri, tsol
    use test_support, only: assert_true, assert_close_scalar, assert_close_vec, finalize_tests
    implicit none

    call test_1x1()
    call test_2x2()
    call test_4x4()
    call test_singular()

    call finalize_tests()

contains

    subroutine test_1x1()
        integer, parameter :: n = 1, lda = 1, m1 = 0, m21 = 1, m3 = 1
        real(wp) :: a(lda,m21), tl(lda,m3), b(n), x(n)
        integer :: ipiv(n), ifail

        a(1,1) = 5.0_wp
        b(1) = 1.0_wp

        call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
        call assert_true("ltri_1x1_success", ifail == 0)

        call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)
        call assert_close_scalar("tsol_1x1", x(1), 0.2_wp, 1.0e-12_wp)
    end subroutine

    subroutine test_2x2()
        integer, parameter :: n = 2, lda = 2, m1 = 1, m21 = 3, m3 = 2
        real(wp) :: a(lda,m21), tl(lda,m3), b(n), x(n)
        integer :: ipiv(n), ifail

        a(1,:) = [0.0_wp, 4.0_wp, 1.0_wp]
        a(2,:) = [1.0_wp, 3.0_wp, 0.0_wp]
        b = [1.0_wp, 2.0_wp]

        call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
        call assert_true("ltri_2x2_success", ifail == 0)

        call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)
        call assert_close_vec("tsol_2x2", x, [1.0_wp/11.0_wp, 7.0_wp/11.0_wp], 1.0e-12_wp)
    end subroutine

    subroutine test_4x4()
        integer, parameter :: n = 4, lda = 4, m1 = 1, m21 = 3, m3 = 2
        real(wp) :: a(lda,m21), tl(lda,m3), b(n), x(n)
        integer :: ipiv(n), ifail

        a(1,:) = [0.0_wp, -2.0_wp, 1.0_wp]
        a(2,:) = [1.0_wp, -2.0_wp, 1.0_wp]
        a(3,:) = [1.0_wp, -2.0_wp, 1.0_wp]
        a(4,:) = [1.0_wp, -2.0_wp, 0.0_wp]

        b = [2.0_wp, 2.0_wp, 4.0_wp - 2.0_wp*9.0_wp + 16.0_wp, 9.0_wp - 2.0_wp*16.0_wp]

        call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
        call assert_true("ltri_4x4_success", ifail == 0)

        call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)
        call assert_close_vec("tsol_4x4", x, [1.0_wp, 4.0_wp, 9.0_wp, 16.0_wp], 1.0e-11_wp)
    end subroutine

    subroutine test_singular()
        integer, parameter :: n = 3, lda = 3, m1 = 1, m21 = 3, m3 = 2
        real(wp) :: a(lda,m21), tl(lda,m3)
        integer :: ipiv(n), ifail

        a(1,:) = [0.0_wp, 1.0_wp, 2.0_wp]
        a(2,:) = [2.0_wp, 4.0_wp, 0.0_wp]
        a(3,:) = [0.0_wp, 1.0_wp, 0.0_wp]

        call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
        call assert_true("ltri_singular_detected", ifail /= 0)
    end subroutine

end program
