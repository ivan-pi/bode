
!
! This is the test of the Robertson Rate equations
!
! Some trivial equations are added to make the system banded.
! Further details can be found here:
! https://github.com/scipy/scipy/issues/10793
!
module robertson_mod
   public
   logical, parameter :: dense_jac = .true.
end module


! q := B p
subroutine pmult(p,n,q)
   use bode_mod, only: wp
   integer, intent(in) :: n
   real(wp), intent(in) :: p(n)
   real(wp), intent(out) :: q(n)

   ! Action of the mass matrix

   ! Unit mass matrix, B = I
   q = p

end subroutine

! Evaluate the right-hand side, q := h*F(x,y)
subroutine deriv(y,n,q,x,h)
   use bode_mod, only: wp
   integer, intent(in) :: n
   real(wp), intent(in) :: y(n), x, h
   real(wp), intent(out) :: q(n)

   q(1) = -0.04_wp*y(1) + 1.0e4_wp*y(2)*y(3)
   q(2) =  0.04_wp*y(1) - 1.0e4_wp*y(2)*y(3) - 3.0e7_wp*y(2)**2
   q(3) =                                      3.0e7_wp*y(2)**2

   q = h * q

end subroutine

! Evaluate the Jacobian of F(x,y), dF/dy, banded version
subroutine pjacb(x,y,n,jac,ld,m1,h)
   use bode_mod, only: wp
   use robertson_mod, only: dense_jac
   integer, intent(in) :: n
   real(wp), intent(in) :: y(n), x
   real(wp), intent(inout) :: jac(ld,*)
   real(wp), intent(in) :: h

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

   if (dense_jac) then
      jac(1:3,1:3) = j
   else
      !
      ! Pack into banded storage
      !

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

program test

   use, intrinsic :: iso_fortran_env, only: error_unit

   use bode_mod, only: wp, bode
   use robertson_mod, only: dense_jac

   implicit none
   integer, parameter :: neq = 3
   real(wp) :: yn(neq), ymin(neq)

   real(wp) :: t, tout, emax, xstep
   integer :: imn, m1, ifail
   integer :: iopt(3)

   ymin = 1.0e-7_wp
   emax = 1.0e-5_wp
   xstep = 0.0001_wp     ! First step
   imn = 1  ! Monitor frequency
   m1 = 2   ! Bandwidth

   iopt(1) = m1
   iopt(2) = 1   ! FD (0) or symbolic (1) Jacobian

   iopt(3) = merge(0,1,dense_jac)

   print *, "iopt = ", iopt

   t = 0

   ! Initial condition
   yn = 0
   yn(1) = 1

   tout = 1000.0

   call bode(t,tout,neq,yn,ymin,emax,xstep,monit,imn,iopt,ifail)
   if (ifail /= 0) then
      write(error_unit,'("BODE failed at t = ",G0," with IFAIL = ",I0)') tout, ifail
      error stop
   end if

   call print_stats

contains

#ifdef ORIG
   subroutine monit(y,n,x)
      integer, intent(in) :: n
      real(wp), intent(in) :: y(n), x
      print '(4(ES12.5,1X))', x, y
   end subroutine
#else
   subroutine monit(y,n,x,iha,qa)
      integer, intent(in) :: n, iha
      real(wp), intent(in) :: y(n), x, qa
      print '(4(ES12.5,1X),I0,1X,ES12.5)', x, y, iha, qa
   end subroutine
#endif

   subroutine print_stats
      use bode_mod, only: nfev, njev, nlu, nbsol
      print *, " nfev: ", nfev
      print *, " njev: ", njev
      print *, "  nlu: ", nlu
      print *, "nbsol: ", nbsol
   end subroutine

end program