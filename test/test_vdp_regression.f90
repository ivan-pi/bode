module vdp_regression_mod
    use bode_mod, only: wp
    implicit none
    real(wp), parameter :: mu = 10.0_wp, alpha = 1.0_wp
end module

subroutine pmult(p,n,q)
    use vdp_regression_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: p(n)
    real(wp), intent(out) :: q(n)
    q = p
end subroutine

subroutine deriv(y,n,q,x,h)
    use vdp_regression_mod, only: wp, mu, alpha
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(out) :: q(n)
    real(wp) :: hh

    hh = alpha*h
    q(1) = y(2)
    q(2) = mu*(1.0_wp - y(1)**2)*y(2) - y(1)
    q = hh*q
end subroutine

subroutine pjacb(x,y,n,jac,ld,m1,h)
    use vdp_regression_mod, only: wp, mu, alpha
    integer, intent(in) :: n, ld, m1
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(inout) :: jac(ld,*)
    real(wp) :: j(2,2), hh

    hh = alpha*h
    j(1,1) = 0.0_wp
    j(1,2) = 1.0_wp
    j(2,1) = -2.0_wp*mu*y(1)*y(2) - 1.0_wp
    j(2,2) = mu*(1.0_wp - y(1)**2)

    jac(1:n,1:n) = 0.0_wp
    jac(1,1) = 1.0_wp + hh*j(1,1)
    jac(1,2) = hh*j(1,2)
    jac(2,1) = hh*j(2,1)
    jac(2,2) = 1.0_wp + hh*j(2,2)
end subroutine

program test_vdp_regression
    use bode_mod, only: wp, bode
    use test_support, only: assert_close_scalar, assert_true, finalize_tests
    implicit none

    integer, parameter :: neq = 2
    real(wp) :: y(neq), ymin(neq)
    real(wp) :: t, tout, emax, xstep
    integer :: iopt(3), ifail
    real(wp) :: ref_t, ref_y1, ref_y2
    character(len=256) :: ref_file

    if (command_argument_count() /= 1) then
        error stop "Usage: test_vdp_regression <reference_file>"
    end if
    call get_command_argument(1, ref_file)

    ymin = 1.0e-6_wp
    emax = 1.0e-6_wp
    xstep = 0.001_wp

    iopt(1) = 1
    iopt(2) = 1
    iopt(3) = 0

    t = 0.0_wp
    tout = 10.0_wp
    y = 1.0_wp

    call bode(t, tout, neq, y, ymin, emax, xstep, monit, 0, iopt, ifail)
    call assert_true("vdp_ifail_zero", ifail == 0)

    call read_reference(ref_file, ref_t, ref_y1, ref_y2)

    call assert_close_scalar("vdp_tout", tout, ref_t, 1.0e-12_wp)
    call assert_close_scalar("vdp_y1", y(1), ref_y1, 2.0e-4_wp)
    call assert_close_scalar("vdp_y2", y(2), ref_y2, 2.0e-3_wp)

    call finalize_tests()

contains

    subroutine monit(y,n,x,iha,qa)
        integer, intent(in) :: n, iha
        real(wp), intent(in) :: y(n), x, qa
    end subroutine

    subroutine read_reference(path, t, y1, y2)
        character(len=*), intent(in) :: path
        real(wp), intent(out) :: t, y1, y2
        integer :: u, ios
        open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "could not open vdp reference"
        read(u,*) t, y1, y2
        close(u)
    end subroutine

end program
