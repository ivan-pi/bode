

!> Fortran version of the routine BANDET1
!< Martin & Wilkinson, Num. Math. Vol 9, p279-301 (1967)
!
! Factorization of a banded matrix A = LU
!
! a ... input matrix, modified on output
! tl ... factors
! n ... size of matrix
! ipiv ... pivot array
! m21,m1,m3 ... sizes ??
! ifail ... error flag (0 on success, 1 if singular)
! 
subroutine ltri(a,tl,ld,n,ipiv,m21,m1,m3,ifail)
    use bode_mod, only: wp
    implicit none
    integer, intent(in) :: ld, m1, m3, m21, n
    integer, intent(out) :: ipiv(n), ifail
    real(wp), intent(inout) :: a(ld,m21), tl(ld,m3)

    integer :: i, ik, j, jl, k, k1, l, m21l, m4
    real(wp) :: x

    ifail = 0

    l = m1
    do i = 1, m1
        m4 = m1 + i
        do j = 1, m4
            jl = j + l
            a(i,j) = a(i,jl)
        end do
        l = l - 1
        m21l = m21 - l
        do j = m21l, m21
            a(i,j) = 0.0_wp
        end do
    end do

    l = m1
    do k = 1, n
        x = a(k,1)
        i = k 
        if (l < n) l = l + 1
        k1 = k + 1
        do j = k1, l
            if (abs(a(j,1)) < abs(x)) cycle
            x = a(j,1)
            i = j
        end do
        ipiv(k) = i
        if (x == 0.0_wp) then
            ifail = 1    ! zero pivot
            return
        end if
        if (i /= k) then
            ! swap rows
            do j = 1, m21
                x = a(k,j)
                a(k,j) = a(i,j)
                a(i,j) = x
            end do
        end if
        k1 = k + 1
        do i = k1, l
            x = a(i,1)/a(k,1)
            tl(k,i-k) = x
            do j = 2, m21
                a(i,j-1) = a(i,j) - x*a(k,j)
            end do
            a(i,m21) = 0.0_wp
        end do
    end do 

end subroutine

!> Fortran version of the routine BANSOL1
!< Martin & Wilkinson, Num. Math. Vol 9, p279-301 (1967)
!
! Solves banded system A x = b
!
! A, TL ... matrix factors
! n ... size of linear system
! m1, m3, m21 ... sizes of ???
! ipiv ... vector of pivots
! xin ... b
! xout ... x
!
! TODO: modify routine to over-write vector b,
!       users should make a copy externally if needed
! 
subroutine tsol(a,tl,ld,m1,m3,m21,n,ipiv,xin,xout)
    use bode_mod, only: wp
    implicit none
    integer, intent(in) :: ld, m1, m3, m21, n, ipiv(n)
    real(wp), intent(in) :: a(ld,m21), tl(ld,m3), xin(n)
    real(wp), intent(out) :: xout(n)

    integer :: i, l, k, k1, ik, ii
    real(wp) :: x

    ! Copy input vector
    do i = 1, n
        xout(i) = xin(i)
    end do

    ! Solve Lf = P^(-1)b for f
    l = m1
    do k = 1, n
        i = ipiv(k)
        if (i /= k) then
            x = xout(k)
            xout(k) = xout(i)
            xout(i) = x
        end if
        if (l < n) l = l + 1
        k1 = k + 1
        do i = k1, l
          ik = i - k
          x = tl(k,ik)
          xout(i) = xout(i) - x*xout(k)
        end do
    end do

    ! Solve Ux = f for x
    l = 1
    do ii = 1, n
        i = n - ii + 1
        x = xout(i)
        do k = 2, l
            k1 = k + i - 1
            x = x - a(i,k) * xout(k1)
        end do
        xout(i) = x/a(i,1)
        if (l < m21) l = l + 1
    end do

end subroutine