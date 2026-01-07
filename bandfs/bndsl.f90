! sbndslv computes the solution to the system of linear equations
!   A*X = B for structurally symmetric banded matrices (simple driver)
subroutine sbndsl(n,k,nrhs,a,lda,ipiv,b,ldb,info)
implicit none

integer, parameter :: wp = kind(1.0e0)
integer, intent(in) :: n, k, nrhs, lda, ldb
real(wp), intent(inout) :: a(lda,*)
integer, intent(out) :: ipiv(n)
real(wp), intent(inout) :: b(ldb,nrhs)
integer, intent(out) :: info

external :: sbandf, sbands
integer :: i

! Test the input parameters
info = 0
if (n < 0) then
    info = -1
else if (k < 0) then
    info = -2
else if (nrhs < 0) then
    info = -3
else if (lda < n) then
    info = -4
else if (ldb < max(n,1)) then
    info = -8
end if

if (info /= 0) return

! Compute the LU Factorization of the structurally symmetric band matrix A
call sbandf(n, k, a, lda, ipiv, info)
if (info == 0) then
    ! Solve the system A*X = B, overwriting B with X
    do i = 1, nrhs
        call sbands(n, k, a, lda, ipiv, b(1,i), info)
        if (info /= 0) return
    end do
end if

end subroutine
!
! dbndslv computes the solution to the system of linear equations
!   A*X = B for structurally symmetric banded matrices (simple driver)
subroutine dbndsl(n,k,nrhs,a,lda,ipiv,b,ldb,info)
implicit none

integer, parameter :: wp = kind(1.0d0)
integer, intent(in) :: n, k, nrhs, lda, ldb
real(wp), intent(inout) :: a(lda,*)
integer, intent(out) :: ipiv(n)
real(wp), intent(inout) :: b(ldb,nrhs)
integer, intent(out) :: info

external :: dbandf, dbands
integer :: i

! Test the input parameters
info = 0
if (n < 0) then
    info = -1
else if (k < 0) then
    info = -2
else if (nrhs < 0) then
    info = -3
else if (lda < n) then
    info = -4
else if (ldb < max(n,1)) then
    info = -8
end if

if (info /= 0) return

! Compute the LU Factorization of the structurally symmetric band matrix A
call dbandf(n, k, a, lda, ipiv, info)
if (info == 0) then
    ! Solve the system A*X = B, overwriting B with X
    do i = 1, nrhs
        call dbands(n, k, a, lda, ipiv, b(1,i), info)
        if (info /= 0) return
    end do
end if

end subroutine
!
