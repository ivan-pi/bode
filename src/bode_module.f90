module bode_mod

implicit none
public

! Work precision
integer, parameter :: wp = kind(0.0)

! Leading dimension of factor arrays
integer, parameter :: ld = 75

abstract interface
   subroutine pmult(p,n,q)
     integer, intent(in) :: n
     real, intent(in) :: p(n)
     real, intent(out) :: q(n)
   end subroutine
   subroutine deriv(y,n,q,x,h)
      integer, intent(in) :: n 
      real, intent(in) :: y(n), x, h 
      real, intent(out) :: q(n)
   end subroutine
end interface

interface
   subroutine ltri(a,tl,n,ipiv,m21,m1,m3,ifail)
    import wp, ld
    implicit none
    integer, intent(in) :: m1, m3, m21, n
    integer, intent(out) :: ipiv(n), ifail
    real(wp), intent(inout) :: a(ld,m21), tl(ld,m3)
   end subroutine
   subroutine tsol(a,tl,m1,m3,m21,n,ipiv,xin,xout)
    import wp, ld
    implicit none
    integer, intent(in) :: m1, m3, m21, n, ipiv(n)
    real(wp), intent(in) :: a(ld,m21), tl(ld,m3), xin(n)
    real(wp), intent(out) :: xout(n)
   end subroutine
end interface

end module