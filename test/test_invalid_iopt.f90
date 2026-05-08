module invalid_iopt_mod
    use bode_mod, only: wp
    implicit none
end module

subroutine pmult(p,n,q)
    use invalid_iopt_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: p(n)
    real(wp), intent(out) :: q(n)
    q = p
end subroutine

subroutine deriv(y,n,q,x,h)
    use invalid_iopt_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(out) :: q(n)
    q(1) = h * y(2)
    q(2) = -h * y(1)
end subroutine

subroutine pjacb(x,y,n,jac,ld,m1,h)
    use invalid_iopt_mod, only: wp
    integer, intent(in) :: n, ld, m1
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(inout) :: jac(ld,*)
    jac(1:n,1:2*m1+1) = 0.0_wp
end subroutine

program test_invalid_iopt
    use bode_mod, only: wp, bode
    implicit none

    integer, parameter :: neq = 2
    real(wp) :: y(neq), ymin(neq), xin, xout, emax, xstep
    integer :: iopt(3), ifail

    xin = 0.0_wp
    xout = 1.0_wp
    y = [1.0_wp, 0.0_wp]
    ymin = 1.0e-8_wp
    emax = 1.0e-6_wp
    xstep = 0.01_wp

    iopt(1) = 1
    iopt(2) = 0
    iopt(3) = 0

    call bode(xin, xout, neq, y, ymin, emax, xstep, monit, 0, iopt, ifail)

contains

    subroutine monit(y,n,x,iha,qa)
        integer, intent(in) :: n, iha
        real(wp), intent(in) :: y(n), x, qa
    end subroutine

end program
