program test_bandfs_api
    use bode_mod, only: wp
    use test_support, only: assert_true, assert_close_vec, finalize_tests
    implicit none

    include 'bandfs.fi'

    call test_bandf_bands()
    call test_bndmv()
    call test_bndsl_multi_rhs()

    call finalize_tests()

contains

    subroutine test_bandf_bands()
        integer, parameter :: n = 5, k = 3, kh = 1, lda = n
        real(wp) :: a(lda,-kh:2*kh), rhs(n), x_expected(n)
        integer :: ipiv(n), info

        a(:, -1) = [0.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp]
        a(:,  0) = [2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp]
        a(:,  1) = [-1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, 0.0_wp]
        a(:,  2) = 0.0_wp

        rhs = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp]
        x_expected = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]

        call bandf(n, k, a, lda, ipiv, info)
        call assert_true("bandf_success", info == 0)

        call bands(n, k, a, lda, ipiv, rhs, info)
        call assert_true("bands_success", info == 0)
        call assert_close_vec("bands_solution", rhs, x_expected, 1.0e-11_wp)
    end subroutine

    subroutine test_bndmv()
        integer, parameter :: n = 5, kl = 1, ku = 1, lda = n
        real(wp) :: a(lda,-kl:ku), x(n), y(n), expected(n)

        a(:, -1) = [0.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp]
        a(:,  0) = [2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp]
        a(:,  1) = [-1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, 0.0_wp]

        x = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
        y = 0.0_wp

        call bndmv(n, kl, ku, 1.0_wp, a, lda, x, 0.0_wp, y)

        expected = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp]
        call assert_close_vec("bndmv_result", y, expected, 1.0e-12_wp)
    end subroutine

    subroutine test_bndsl_multi_rhs()
        integer, parameter :: n = 5, k = 3, kh = 1, lda = n, ldb = n, nrhs = 2
        real(wp) :: a(lda,-kh:2*kh), b(ldb,nrhs), expected(ldb,nrhs)
        integer :: ipiv(n), info

        a(:, -1) = [0.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp]
        a(:,  0) = [2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp]
        a(:,  1) = [-1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, 0.0_wp]
        a(:,  2) = 0.0_wp

        b(:,1) = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 6.0_wp]
        b(:,2) = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

        expected(:,1) = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
        expected(:,2) = [5.0_wp/6.0_wp, 2.0_wp/3.0_wp, 1.0_wp/2.0_wp, 1.0_wp/3.0_wp, 1.0_wp/6.0_wp]

        call dbndsl(n, k, nrhs, a, lda, ipiv, b, ldb, info)
        call assert_true("bndsl_success", info == 0)
        call assert_close_vec("bndsl_rhs1", b(:,1), expected(:,1), 1.0e-11_wp)
        call assert_close_vec("bndsl_rhs2", b(:,2), expected(:,2), 1.0e-11_wp)
    end subroutine

end program
