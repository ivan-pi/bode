subroutine dbndmv(n, kl, ku, alpha, ab, x, beta, y)

    integer, parameter :: wp = kind(1.0d0)

    integer,  intent(in)  :: n             ! Matrix dimension
    integer,  intent(in)  :: kl, ku        ! Lower and upper bandwidth
    real(wp), intent(in)  :: alpha         ! Scalar multiplier for Ax
    real(wp), intent(in)  :: ab(n, -kl:ku) ! Banded storage with offset indexing
    real(wp), intent(in)  :: x(n)          ! Input vector
    real(wp), intent(in)  :: beta          ! Scalar multiplier for y
    real(wp), intent(inout) :: y(n)        ! Output vector (modified in place)
    
    integer :: i, j
    real(wp) :: row_sum

    ! 1. Handle the beta scaling first
    ! If beta is 0, we overwrite y. If beta is 1, we leave it.
    if (beta /= 1.0_wp) then
        if (beta == 0.0_wp) then
            y = 0.0_wp
        else
            y = y * beta
        end if
    end if

    ! 2. Compute alpha * A * x and add to y
    if (alpha == 0.0_wp) return ! Quick return if no matrix contribution

    do i = 1, n
        row_sum = 0.0_wp
        ! Inner loop over the band columns
        do j = max(1, i - kl), min(n, i + ku)
            ! k is the relative offset (j - i)
            row_sum = row_sum + ab(i, j - i) * x(j)
        end do
        y(i) = y(i) + alpha * row_sum
    end do
end subroutine