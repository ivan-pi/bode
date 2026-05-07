module dense_vs_banded_mod
    use bode_mod, only: wp
    implicit none
    logical :: use_dense_jac = .true.
end module

subroutine pmult(p,n,q)
    use dense_vs_banded_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: p(n)
    real(wp), intent(out) :: q(n)
    q = p
end subroutine

subroutine deriv(y,n,q,x,h)
    use dense_vs_banded_mod, only: wp
    integer, intent(in) :: n
    real(wp), intent(in) :: y(n), x, h
    real(wp), intent(out) :: q(n)

    q(1) = -0.04_wp*y(1) + 1.0e4_wp*y(2)*y(3)
    q(2) =  0.04_wp*y(1) - 1.0e4_wp*y(2)*y(3) - 3.0e7_wp*y(2)**2
    q(3) =                                      3.0e7_wp*y(2)**2
    q = h*q
end subroutine

subroutine pjacb(x,y,n,jac,ld,m1,h)
    use dense_vs_banded_mod, only: wp, use_dense_jac
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

    jac(1:n,1:2*m1+1) = 0.0_wp
    if (use_dense_jac) then
        jac(1:3,1:3) = j
    else
        jac(3,1) = j(3,1)
        jac(2,2) = j(2,1)
        jac(3,2) = j(3,2)
        jac(1,3) = j(1,1)
        jac(2,3) = j(2,2)
        jac(3,3) = j(3,3)
        jac(1,4) = j(1,2)
        jac(2,4) = j(2,3)
        jac(1,5) = j(1,3)
    end if
end subroutine

program dense_vs_banded_robertson
    use bode_mod, only: wp, bode
    use dense_vs_banded_mod, only: use_dense_jac
    implicit none

    real(wp) :: yd(3), yb(3), diff_norm

    use_dense_jac = .true.
    call solve_robertson(yd)

    use_dense_jac = .false.
    call solve_robertson(yb)

    diff_norm = norm2(yd - yb)

    print '(A,3(ES14.6,1X))', 'dense  :', yd
    print '(A,3(ES14.6,1X))', 'banded :', yb
    print '(A,ES14.6)', 'difference norm: ', diff_norm

contains

    subroutine solve_robertson(y)
        real(wp), intent(out) :: y(3)
        real(wp) :: ymin(3), t, tout, emax, xstep
        integer :: iopt(3), ifail

        ymin = 1.0e-7_wp
        emax = 1.0e-5_wp
        xstep = 1.0e-4_wp

        y = 0.0_wp
        y(1) = 1.0_wp

        t = 0.0_wp
        tout = 1000.0_wp

        iopt(1) = 2
        if (use_dense_jac) then
            iopt(2) = 1
            iopt(3) = 0
        else
            iopt(2) = 0
            iopt(3) = 1
        end if

        call bode(t, tout, 3, y, ymin, emax, xstep, monit, 0, iopt, ifail)
        if (ifail /= 0) error stop "BODE failed in dense_vs_banded example"
    end subroutine

    subroutine monit(y,n,x,iha,qa)
        integer, intent(in) :: n, iha
        real(wp), intent(in) :: y(n), x, qa
    end subroutine

end program
