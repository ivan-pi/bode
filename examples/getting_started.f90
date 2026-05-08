module getting_started_mod
    use bode_mod, only: wp
    implicit none
end module

subroutine pmult(p,n,q)
    use getting_started_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: p(n)
    real(wp), intent(out) :: q(n)
    q = p
end subroutine

subroutine deriv(y,n,q,x,h)
    use getting_started_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(out) :: q(n)
    q(1) = -h * y(1)
end subroutine

subroutine pjacb(x,y,n,jac,ld,m1,h)
    use getting_started_mod, only: wp
    integer, intent(in) :: n, ld, m1
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(inout) :: jac(ld,2*m1+1)
    jac = 0.0_wp
    jac(1,1) = 1.0_wp + h
end subroutine

program getting_started
    use bode_mod, only: wp, bode
    implicit none

    integer, parameter :: neq = 1
    real(wp) :: y(neq), ymin(neq), xin, xout, emax, xstep
    integer :: iopt(3), ifail

    xin = 0.0_wp
    xout = 1.0_wp
    y(1) = 1.0_wp

    ymin(1) = 1.0e-12_wp
    emax = 1.0e-7_wp
    xstep = 1.0e-3_wp

    iopt(1) = 0
    iopt(2) = 1
    iopt(3) = 0

    call bode(xin, xout, neq, y, ymin, emax, xstep, monit, 0, iopt, ifail)
    if (ifail /= 0) error stop "BODE integration failed"

    print '(A,F10.6)', 't = ', xout
    print '(A,F10.6)', 'y_numeric = ', y(1)
    print '(A,F10.6)', 'y_exact   = ', exp(-xout)

contains

    subroutine monit(y,n,x,iha,qa)
        integer, intent(in) :: n, iha
        real(wp), intent(in) :: y(n), x, qa
    end subroutine

end program
