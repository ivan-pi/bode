program banded_example

use, intrinsic :: iso_fortran_env, only: error_unit
implicit none

integer, parameter :: wp = kind(1.0d0)

! Constants
integer, parameter :: n = 5
integer, parameter :: ml = 1, mu = 2
integer, parameter :: m1 = max(ml,mu)
integer, parameter :: bw = 2*m1 + 1

! Dense matrix and vectors
real(wp) :: ad(n,n), b(n), xla(n), xsl(n)
real(wp) :: wrk(n)
integer :: ipiv(n)

! SLATEC format (extra columns needed for pivoting)
real(wp) :: as(n,bw+m1)
integer :: itask, ind

! LAPACK format (additional ml rows needed by factorization)
real(wp) :: ab(bw+m1,n)
integer :: info

! BANDFS procedure interfaces
include "bandfs.fi"

! Dense matrix
ad(1,:) = [1, 2, 3, 0, 0]
ad(2,:) = [4, 5, 6, 7, 0]
ad(3,:) = [0, 8, 9, 1, 2]
ad(4,:) = [0, 0, 3, 4, 5]
ad(5,:) = [0, 0, 0, 6, 7]

! Right-hand side
! A * x = [6, 22, 20, 12, 13]^T, x = [1]_n
b = sum(ad,dim=2)

! Convert dense to banded (extra storage is provided for the factorization)
ab = dense2bandedLA(ad,m1,m1)
as = dense2bandedSLATEC(ad,m1,m1)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
! Solve banded system A x = b with LAPACK
!

xla = b
call dgbsv(n,m1,m1,1,ab,size(ab,1),ipiv,xla,n,info)
if (info /= 0) then
   write(error_unit,*) "[dgbsv] an error occured"
   error stop info
end if

write(*,'(A,*(/,G0))') "LAPACK: x = ", xla

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
! Solve banded system A x = b with BANDF/BANDS
! (only structurally symmetric banded matrices are allowed)
!

xsl = b
call bandf(as,bw,n,ipiv,info)
if (info /= 0) then
   write(*,'("bandf failed with info = ",I0)') info
   error stop
end if
call bands(as,xsl,bw,n,ipiv)
write(*,'(A,*(/,G0))') "BANDFS: x = ", xsl

contains

!
! Convert to banded format using LAPACK storage
!
! Description found here:
! - https://www.netlib.org/lapack/lug/node124.html
function dense2bandedLA(A, ml, mu) result(AB)
   integer, intent(in) :: ml, mu
   real(wp), intent(in) :: A(:,:)
   real(wp) :: AB(2*ml+mu+1, size(A,1))
   integer :: j, i
   ab = 0
   do j = 1, size(A,1)
      do i = max(1, j - mu), min(size(A,1), j+ml)
         AB(ml+mu+1+i-j,j) = A(i,j)
      end do
   end do
end function

!
! Convert to banded using SLATEC storage
! (diagonals are stored as columns)
!
function dense2bandedSLATEC(A, ml, mu) result(ab)
   integer, intent(in) :: ml, mu
   real(wp), intent(in) :: A(:,:)

   real(wp) :: ab(size(A,1),2*ml+mu+1)

   integer :: i, j1, j2, j, k
   ab = 0
   do i = 1, n
      j1 = max(1,i-ml)
      j2 = min(n,i+mu)
      do j = j1, j2
         k = j - i + ml + 1
         ab(i,k) = a(i,j)
      end do
   end do
end function

end program