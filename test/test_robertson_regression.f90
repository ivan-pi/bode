module robertson_regression_mod
    use bode_mod, only: wp
    implicit none
end module

subroutine pmult(p,n,q)
    use robertson_regression_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: p(n)
    real(wp), intent(out) :: q(n)
    q = p
end subroutine

subroutine deriv(y,n,q,x,h)
    use robertson_regression_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(out) :: q(n)

    q(1) = -0.04_wp*y(1) + 1.0e4_wp*y(2)*y(3)
    q(2) =  0.04_wp*y(1) - 1.0e4_wp*y(2)*y(3) - 3.0e7_wp*y(2)**2
    q(3) =                                      3.0e7_wp*y(2)**2
    q = h*q
end subroutine

subroutine pjacb(x,y,n,jac,ld,m1,h)
    use robertson_regression_mod, only: wp
    integer, intent(in) :: n, ld, m1
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(inout) :: jac(ld,*)
    real(wp) :: j(3,3)

    j(1,1) = -0.04_wp
    j(2,1) =  0.04_wp
    j(3,1) =  0.0_wp

    j(1,2) =  1.0e4_wp*y(3)
    j(2,2) = -1.0e4_wp*y(3) - 2*3.0e7_wp*y(2)
    j(3,2) =                  2*3.0e7_wp*y(2)

    j(1,3) =  1.0e4_wp*y(2)
    j(2,3) = -1.0e4_wp*y(2)
    j(3,3) =  0.0_wp

    j = -h*j
    j(1,1) = 1.0_wp + j(1,1)
    j(2,2) = 1.0_wp + j(2,2)
    j(3,3) = 1.0_wp + j(3,3)

    jac(1:3,1:3) = j
end subroutine

program test_robertson_regression
    use bode_mod, only: wp, bode
    use test_support, only: assert_close_scalar, assert_true, finalize_tests
    implicit none

    integer, parameter :: neq = 3
    real(wp) :: y(neq), ymin(neq)
    real(wp) :: t, tout, emax, xstep
    integer :: iopt(3), ifail
    real(wp) :: ref_t, ref_y1, ref_y2, ref_y3
    character(len=256) :: ref_file

    if (command_argument_count() /= 1) then
        error stop "Usage: test_robertson_regression <reference_file>"
    end if
    call get_command_argument(1, ref_file)

    ymin = 1.0e-7_wp
    emax = 1.0e-5_wp
    xstep = 0.0001_wp

    iopt(1) = 2
    iopt(2) = 1
    iopt(3) = 0

    t = 0.0_wp
    tout = 1000.0_wp
    y = 0.0_wp
    y(1) = 1.0_wp

    call bode(t, tout, neq, y, ymin, emax, xstep, monit, 0, iopt, ifail)
    call assert_true("robertson_ifail_zero", ifail == 0)

    call read_reference(ref_file, ref_t, ref_y1, ref_y2, ref_y3)

    call assert_close_scalar("robertson_tout", tout, ref_t, 1.0e-12_wp)
    call assert_close_scalar("robertson_y1", y(1), ref_y1, 2.0e-4_wp)
    call assert_close_scalar("robertson_y2", y(2), ref_y2, 5.0e-8_wp)
    call assert_close_scalar("robertson_y3", y(3), ref_y3, 2.0e-4_wp)

    call finalize_tests()

contains

    subroutine monit(y,n,x,iha,qa)
        integer, intent(in) :: n, iha
        real(wp), intent(in) :: y(n), x, qa
    end subroutine

    subroutine read_reference(path, t, y1, y2, y3)
        character(len=*), intent(in) :: path
        real(wp), intent(out) :: t, y1, y2, y3
        integer :: u, ios
        open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "could not open robertson reference"
        read(u,*) t, y1, y2, y3
        close(u)
    end subroutine

end program
