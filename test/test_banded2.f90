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

        print *, "m1  = ", m1
        print *, "m21 = ", m21
        print *, "m3  = ", m3
        
        tl = -66
        call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
        if (ifail /= 0) return

        print *, "FACTORIZATION"
        do i = 1, n
            print *, "tl(i,:) = ", tl(i,:), "|", ipiv(i)
        end do

        ! Copy input values
        x(1:n) = b(1:n)
        call tsol(a,tl,lda, &
            m1=m1,m3=m3,m21=m21,n=n,&
            ipiv=ipiv,xin=x,xout=b)

    end subroutine

    ! Reaction-Diffusion Example 1
    ! Uses backward difference (first order) for the Neumann BC
    subroutine setup1(n, ab, b)
        integer, intent(in) :: n
        real(wp), intent(out) :: ab(n,3)
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
    subroutine setup2(n, ab, b)
        integer, intent(in) :: n
        real(wp), intent(out) :: ab(n,3)
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
    subroutine setup3
        error stop "setup3: Not Implemented"
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
use diff_mod, only: setup1, setup2, solve_banded, sol

implicit none

integer, parameter :: kl = 1, ku = 1
integer, parameter :: n = 11

integer, parameter :: bw = kl + ku + 1

real(wp) :: ab(n, 5), b(n), h, x(n), c(n)
integer :: ipiv(n), info, k, i

! Half-bandwidth
k = 1

! Setup reaction-diffusion equation
call setup2(n, ab, b)

print *, "ab = "
do i = 1, n
    print *, ab(i,:), "|", b(i)
end do

! General banded matrix driver
call solve_banded(n,k,ab,lda=n,ipiv=ipiv,b=b,ifail=info)
print *, "# factor and solve info = ", info

do i = 1, n
    print *, ab(i,:), "|", b(i) , "|", ipiv(i)
end do


h = 1.0_wp / real(n-1,wp)

x = [((i-1)*h,i=1,n)]
c = sol(x,Th=1.0_wp)

print *, "SOLUTION"
do i = 1, n
    print *, x(i), c(i), b(i)
end do

end program