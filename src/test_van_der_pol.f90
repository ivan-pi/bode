
!
! Van der Pol example
!

module vdp_mod
   use bode_mod, only: wp
   real(wp), parameter :: mu = 10.0_wp

   !real(wp) :: alpha = 2.333333_wp
   !real(wp), parameter :: alpha = 1.67_wp
   real(wp), parameter :: alpha = 1.0_wp

end module

! q := B p
subroutine pmult(p,n,q)
   use bode_mod, only: wp
   integer, intent(in) :: n
   real(wp), intent(in) :: p(n)
   real(wp), intent(out) :: q(n)

   ! Action of mass matrix

   ! Unit mass matrix, B = I
   q = p

end subroutine

! Evaluate the right-hand side, q := h*F(x,y)
subroutine deriv(y,n,q,x,h)
   use bode_mod, only: wp
   use vdp_mod, only: mu, alpha
   integer, intent(in) :: n
   real(wp), intent(in) :: y(n), x, h
   real(wp), intent(out) :: q(n)

   real(wp) :: hh

   hh = alpha*h

   q(1) = y(2)
   q(2) = mu*(1.0_wp - y(1)**2)*y(2) - y(1)
   !q(3) = y(3)

   q = hh * q

end subroutine

! Evaluate the Jacobian of F(x,y), dF/dy, banded version
subroutine pjacb(x,y,n,jac,ld,m1,h)
   use bode_mod, only: wp
   use vdp_mod, only: mu, alpha
   integer, intent(in) :: n
   real(wp), intent(in) :: y(n), x
   real(wp), intent(inout) :: jac(ld,2*m1+1)
   real(wp), intent(in) :: h

   real(wp) :: j(2,2), hh

   hh = alpha*h

   j(1,1) = 0.0_wp
   j(1,2) = 1.0_wp
   j(2,1) = mu*y(2)*(-2*y(1)) - 1.0_wp
   j(2,2) = mu*(1.0_wp - y(1)**2)


   jac(1:n,:) = 0
   !
   ! We must evaluate B + h J
   !

   ! J21
   jac(2,1) = hh*j(2,1)

   ! J12
   jac(1,3) = hh*j(1,2)

   ! Diagonal
   jac(1,2) = 1.0_wp + hh*j(1,1)
   jac(2,2) = 1.0_wp + hh*j(2,2)

end subroutine

program test

   use, intrinsic :: iso_fortran_env, only: error_unit

   use bode_mod, only: wp, bode

   implicit none
   integer, parameter :: neq = 2
   real(wp) :: yn(neq), ymin(neq)

   real(wp) :: t, tout, emax, xstep
   integer :: imn, m1, ifail
   integer :: info(2)

   ymin = 1.0e-6_wp
   emax = 1.0e-6_wp
   xstep = 0.001_wp     ! First step
   imn = 1  ! Monitor frequency
   m1 = 1   ! Bandwidth
   tout = 100.0

   info(1) = m1
   info(2) = 1   ! FD (0) or user-provided Jacobian (1)

   ! Initial condition
   t = 0
   yn = 1.0_wp
   !yn(3) = 0.0_wp

   call monit(yn,neq,t)
   call bode(t,tout,neq,yn,ymin,emax,xstep,monit,imn,info,ifail)
   if (ifail /= 0) then
      write(error_unit,'("BODE failed at t = ",G0," with IFAIL = ",I0)') tout, ifail
      error stop
   end if
   call monit(yn,neq,tout)

contains

   subroutine monit(y,n,x)
      integer, intent(in) :: n
      real(wp), intent(in) :: y(n), x

      print '(3(ES12.5,:,1X))', x, y
   end subroutine

end program