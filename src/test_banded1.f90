program test_banded1

use bode_mod, only: wp, ltri, tsol
implicit none

! lda ... maximum system size we may want to solve
!   m ... number of diagonals in band

integer, parameter :: lda = 4, m = 3
real(wp) :: a(lda,m), tl(lda,m), b(lda), x(lda)
integer :: ipiv(lda), m1, m21, m3, ifail, n

call test1
call test2
call test2_diag
call test3
call test4
call test4_2

contains

subroutine test1
print *, "test1"

! Solving a 1-by-1 system:
!
! | 5 | |x1|   |1|
!
! The solution is given by:  x1 = 1/5

n = 1
a(1,1) = 5
b(1) = 1

m1 = 0
m21 = 2*m1 +1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1)
print *, "expected = ", 1.0_wp/5.0_wp

end subroutine

! 2-by-2 system, bandwidth = 1
subroutine test2
print *, "test2"

! Solving a 2-by-2 system:
!
! | 4 1 | |x1|   |1|
! | 1 3 | |x2| = |2|
!
! The solution is given by:  x^T = [1/11, 7/11]

n = 2
a(1,:) = [0, 4, 1]
a(2,:) = [1, 3, 0]

b(1:2) = [1, 2]

m1 = 1
m21 = 2*m1 +1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1:2)
print *, "expected = ", [1.0_wp/11.0_wp, 7.0_wp/11.0_wp]

end subroutine

! 2-by-2 system, bandwidth = 0
subroutine test2_diag
print *, "test2_diag"

!
! | 4   | |x1|   |1|
! |   3 | |x2| = |1|
!
! The solution is given by:  x1 = 1/4, x2 = 1/3

n = 2
a(1:2,1) = [4, 3]
b(1:2) = [1, 1]

m1 = 0
m21 = 1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1:2)
print *, "expected = ", [1.0_wp/4.0_wp, 1.0_wp/3.0_wp]

end subroutine


subroutine test3
print *, "Test 3"

!---------------------------------------------------
!
! Let's solve G. Strang's favorite system:
!
! | -2  1    | |x1|   |  2|
! |  1 -2  1 | |x2| = |  2|
! |     1 -2 | |x3|   |-14|
!
! The solution is given by:  x^T = [1, 4, 9]
!

n = 3
a(1,:) = [0, -2, 1]
a(2,:) = [1, -2, 1]
a(3,:) = [1, -2, 0]

b(1:3) = [2, 2, -14]

m1 = 1
m21 = 2*m1 +1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1:3)
print *, "expected = ", [1.0_wp, 4.0_wp, 9.0_wp]

end subroutine

! 4-by-4 system, bandwidth=1
subroutine test4
print *, "Test 4"
! -----------------------------------------
!
! | -2  1       | |x1|   |  2|
! |  1 -2  1    | |x2| = |  2|
! |     1 -2  1 | |x3| = |  2|
! |        1 -2 | |x4|   |-23|
!
! The solution is given by:  x^T = [1, 4, 9, 16]
!

n = 4
a(1,:) = [0, -2, 1]
a(2,:) = [1, -2, 1]
a(3,:) = [1, -2, 1]
a(4,:) = [1, -2, 0]

b(1:4) = [2, 2, 4 - 2*9 + 16, 9 - 2*16]

call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1:4)
print *, "expected = ", [1.0_wp, 4.0_wp, 9.0_wp, 16.0_wp]

end subroutine

! 4-by-4 system, bandwidth=2
subroutine test4_2
print *, "Test 4"
! -----------------------------------------
!
! | -2  1       | |x1|   |  2|
! |  1 -2  1    | |x2| = |  2|
! |     1 -2  1 | |x3| = |  2|
! |        1 -2 | |x4|   |-23|
!
! The solution is given by:  x^T = [1, 4, 9, 16]
!

n = 4
a(1,:) = [0, -2, 1]
a(2,:) = [1, -2, 1]
a(3,:) = [1, -2, 1]
a(4,:) = [1, -2, 0]

b(1:4) = [2, 2, 4 - 2*9 + 16, 9 - 2*16]

call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1:4)
print *, "expected = ", [1.0_wp, 4.0_wp, 9.0_wp, 16.0_wp]

end subroutine

end program