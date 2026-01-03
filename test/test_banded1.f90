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

! Example from Martin & Wilkinson (1967)
call test5

! Examples from SlideShare
! https://www.slideshare.net/slideshow/lect07-249491248/249491248 (slide 2)
call test11
call test12
call test13

! Examples from Klein & Strzodka (2023)
call test14
call test15
call test16

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

! 7-by-7 system, bandwidth=2
subroutine test5

integer, parameter :: lda = 7, m = 5
real(wp) :: a(lda,m), tl(lda,m), b(lda), x(lda), exact(lda)
integer :: ipiv(lda), m1, m21, m3, ifail, n

print *, "Test 5"

! This test is taken from
!
! Martin, R. S., & Wilkinson, J. H. (1967)
! Solution of symmetric and unsymmetric band equations and the calculation
! of eigenvectors of band matrices. Numerische Mathematik, 9(4), 279-301.
! https://doi.org/10.1007/BF02162421
!
! The matrix of order seven is given by
!
! |  5 -4  1             | |x1|   |  0|
! | -4  6 -4  1          | |x2| = |  0|
! |  1 -4  6 -4  1       | |x3| = |  0|
! |     1 -4  6 -4  1    | |x4| = |  1|
! |        1 -4  6 -4  1 | |x5| = |  0|
! |           1 -4  6 -4 | |x6| = |  0|
! |              1 -4  5 | |x7| = |  0|
!
! The exact solution is given by x^T = (4, 15/2, 10, 11, 10, 15/2, 4)
!

n = 7
a(1,:) = [0,  0, 5, -4, 1]
a(2,:) = [0, -4, 6, -4, 1]
a(3,:) = [1, -4, 6, -4, 1]
a(4,:) = [1, -4, 6, -4, 1]
a(5,:) = [1, -4, 6, -4, 1]
a(6,:) = [1, -4, 6, -4, 0]
a(7,:) = [1, -4, 5,  0, 0]


b = [real(wp) :: 0, 0, 0, 1, 0, 0, 0]
exact = [4.0_wp, 7.5_wp, 10.0_wp, 11.0_wp, 10.0_wp, 7.5_wp, 4.0_wp]

m1 = 2
m21 = 2*m1 +1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1:n)
print *, "expected = ", exact

end subroutine


! subroutine test6

! print *, "Test 6"

! ! Matrix test taken from
! ! https://www.centerspace.net/doc/NMath/user/matrix-types-79829.htm
! !
! ! This is 5 x 5 symmetric banded matrix with half bandwidth of 1:
! !
! ! |  1 -2          |
! ! | -2 -1  3       |
! ! |     3  0 -1    |
! ! |       -1  1  4 |
! ! |           4  3 |

! n = 5
! A(1,:) = [ 0, 1,-2]
! A(2,:) = [-2,-1, 3]
! A(3,:) = [ 3, 0,-1]
! A(4,:) = [-1, 1, 4]
! A(5,:) = [ 4, 3, 0]

! end subroutine

! subroutine test7

! print *, "Test 7"

! ! Determinant of tridiagonal matrix example taken from
! ! https://math.stackexchange.com/questions/2522089/determinant-of-tridiagonal-banded-matrix

! n = 6
! A(1,:) = [0, 1, 1]
! A(2,:) = [1, 2, 2]
! A(3,:) = [2, 3, 3]
! A(4,:) = [3, 4, 4]
! A(5,:) = [4, 5, 5]
! A(6,:) = [0, 5, 6]

! end subroutine


! subroutine test8

! print *, "Test 8"
! !
! ! Matrix taken taken from the MIT news
! ! https://news.mit.edu/2010/faster-fourier-0729
! !

! A(:,4) = [9, 1, 5, 4, 8, 2, 0, 0]
! A(:,3) = [6, 4, 9, 6, 5, 6, 1, 0]
! A(:,2) = [5, 2, 7, 6, 1, 9, 4, 3]   ! DIAGONAL
! A(:,1) = [0, 3, 1, 3, 1, 7, 5, 8]

! end subroutine

! subroutine test9

! print *, "Test 9"
! ! Matrix example taken from
! ! https://math.stackexchange.com/questions/4777900/lu-decomposition-of-banded-matrix-with-partial-pivoting

! n = 5
! A(1,:) = [ 0, 5, 8]
! A(2,:) = [35,64, 4]
! A(3,:) = [40,22, 9]
! A(4,:) = [18,85, 7]
! A(5,:) = [ 0, 8,19]

! end subroutine


! subroutine test10

! ! Example taken from
! ! https://www.chegg.com/homework-help/questions-and-answers/3-10-points-8-8-banded-matrix-given-follows-6-0-0-0-0-0-1-0-0-0-0-1-5-0-0-0-6-7-3-0-0-3-5--q99214846


! n = 8
! A(:,6) = [ 3, 1, 5,-3, 8, 3, 0, 0]
! A(:,5) = [-6, 9,-1, 7, 2, 4, 1, 0]
! A(:,4) = [ 5, 7, 4, 6, 5, 1,-5, 6] ! DIAGONAL
! A(:,3) = [ 0, 2, 6, 1, 3, 4, 3,-8]
! A(:,2) = [ 0, 0,10,-2,12,11, 7, 9]
! A(:,1) = [ 0, 0, 0, 8,-9, 6, 8, 4]

! mu = 2
! ml = 3

! end subroutine

! 3-by-3 system, bandwidth = 2
subroutine test11
implicit none
integer, parameter :: n = 3, bw = 2
real(wp) :: a(n,2*bw+1), tl(n,2*bw+1), x(n), b(n)
integer :: ipiv(n), m1, m21, m3, ifail, lda

print *, "Test 11"

! Singular matrix
!
! | 1 2 3 |
! | 4 5 6 |
! | 7 8 9 |
!

lda = n

!a(1,:) = [0, 0, 1, 2, 3]
!a(2,:) = [0, 4, 5, 6, 0]
!a(3,:) = [7, 8, 9, 0, 0]

a(:,5) = [3,0,0]
a(:,4) = [2,6,0]
a(:,3) = [1,5,9]
a(:,2) = [0,4,8]
a(:,1) = [0,0,5]

call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail == 0) then
    write(*,*) "Expected singular matrix; LTRI returned with IFAIL = ", ifail
    error stop 1
end if
print *, "Pass"

end subroutine

! 3-by-3 system, bandwidth = 1 (tridiagonal matrix)
subroutine test12
implicit none
integer, parameter :: n = 3, bw = 1
real(wp) :: a(n,2*bw+1), tl(n,bw+1), x(n), b(n)
integer :: ipiv(n), m1, m21, m3, ifail, lda

print *, "Test 12"

!
! | 1 2 0 | |x1|   |  5 |
! | 4 5 6 | |x2| = | 32 |
! | 0 8 9 | |x3|   | 43 |
!
! The solution is given by:  x^T = [1, 2, 3]

! 1 + 4 = 5
! 4 + 10 + 18 = 32
! 0 + 16 + 27 = 43

lda = n

a(:,3) = [2,6,0]
a(:,2) = [1,5,9]
a(:,1) = [0,4,8]

b(1:3) = [5, 32, 43]

m1 = 1
m21 = 2*m1 + 1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1:3)
print *, "expected = ", [1.0_wp, 2.0_wp, 3.0_wp]

end subroutine


! 3-by-3 system, bandwidth = 0 (diagonal matrix)
subroutine test13
implicit none
integer, parameter :: n = 3, bw = 0
real(wp) :: a(n,2*bw+1), tl(n,bw+1), x(n), b(n)
integer :: ipiv(n), m1, m21, m3, ifail, lda

print *, "Test 13"

!
! | 1 0 0 | |x1|   | 1 |
! | 0 5 0 | |x2| = | 10 |
! | 0 0 9 | |x3|   | 27 |
!
! The solution is given by:  x^T = [1, 2, 3]

lda = n

a(:,1) = [1, 5, 9]

b(1:3) = [1,10,27]

m1 = 0
m21 = 2*m1 + 1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x(1:3)
print *, "expected = ", [1.0_wp, 2.0_wp, 3.0_wp]

end subroutine

! 6-by-6 system, bandwidth = 1 (tri-diagonal matrix)
subroutine test14
implicit none
integer, parameter :: n = 6, bw = 1
real(wp) :: a(n,2*bw+1), tl(n,bw+1), x(n), b(n)
integer :: ipiv(n), m1, m21, m3, ifail, lda

print *, "Test 14"

! Example of Problem Class Scalar in
! https://dl.acm.org/doi/epdf/10.1145/3580373

lda = n

a(:,1) = [ 0, 1, 2, 3, 4, 5]
a(:,2) = [ 6, 7, 8, 9,10,11]
a(:,3) = [12,13,14,15,16,17]

b(1:6) = [1,2,3,4,5,6]

m1 = bw
m21 = 2*m1 + 1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x

end subroutine

! 8-by-8 system, bandwidth = 2 (penta-diagonal matrix)
subroutine test15
implicit none
integer, parameter :: n = 8, bw = 3
real(wp) :: a(n,2*bw+1), tl(n,bw+1), x(n), b(n)
integer :: ipiv(n), m1, m21, m3, ifail, lda

print *, "Test 15"

! Example of Problem Class Block in
! https://dl.acm.org/doi/epdf/10.1145/3580373

lda = n

a(:,1) = [ 0, 0, 0, 2, 0, 6, 0,10]
a(:,2) = [ 0, 0, 1, 4, 5, 8, 9,12]
a(:,3) = [ 0,14, 3,18, 7,22,11,26]
a(:,4) = [13,16,17,20,21,24,25,28]
a(:,5) = [15,30,19,34,23,38,27, 0]
a(:,6) = [29,32,33,36,37,40, 0, 0]
a(:,7) = [31, 0,35, 0,39, 0, 0, 0]

b(1:8) = [ 1, 2, 3, 4, 5, 6, 7, 8]

m1 = bw
m21 = 2*m1 + 1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x

end subroutine

! 8-by-8 system, bandwidth = 2 (penta-diagonal matrix)
subroutine test16
implicit none
integer, parameter :: n = 8, bw = 2
real(wp) :: a(n,2*bw+1), tl(n,bw+1), x(n), b(n)
integer :: ipiv(n), m1, m21, m3, ifail, lda

print *, "Test 16"

! Example of Problem Class DIA in
! https://dl.acm.org/doi/epdf/10.1145/3580373

lda = n

a(:,1) = [ 0, 0, 1, 4, 5, 8, 9,12]
a(:,2) = [ 0,14, 3,18, 7,22,11,26]
a(:,3) = [13,16,17,20,21,24,25,28]
a(:,4) = [15,30,19,34,23,38,27, 0]
a(:,5) = [29,32,33,36,37,40, 0, 0]

b(1:8) = [1,2,3,4,5,6,7,8]

m1 = bw
m21 = 2*m1 + 1
m3 = m1 + 1
call ltri(a,tl,lda,n,ipiv,m21,m1,m3,ifail)
if (ifail /= 0) then
    write(*,*) "LTRI failed with error: ", ifail
    error stop 1
end if

call tsol(a,tl,lda,m1,m3,m21,n,ipiv,xin=b,xout=x)

print *, "x        = ", x

end subroutine


end program