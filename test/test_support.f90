module test_support
    use, intrinsic :: iso_fortran_env, only: error_unit
    use bode_mod, only: wp
    implicit none
    private

    integer :: failures = 0

    public :: assert_true, assert_close_scalar, assert_close_vec, finalize_tests

contains

    subroutine assert_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition
        if (.not. condition) then
            failures = failures + 1
            write(error_unit,'("ASSERT_TRUE failed: ",A)') trim(name)
        end if
    end subroutine

    subroutine assert_close_scalar(name, actual, expected, tol)
        character(len=*), intent(in) :: name
        real(wp), intent(in) :: actual, expected, tol
        if (abs(actual - expected) > tol) then
            failures = failures + 1
            write(error_unit,'("ASSERT_CLOSE scalar failed: ",A,2X,"actual=",ES12.5,2X,"expected=",ES12.5,2X,"tol=",ES12.5)') &
                trim(name), actual, expected, tol
        end if
    end subroutine

    subroutine assert_close_vec(name, actual, expected, tol)
        character(len=*), intent(in) :: name
        real(wp), intent(in) :: actual(:), expected(:), tol
        real(wp) :: err
        err = maxval(abs(actual - expected))
        if (err > tol) then
            failures = failures + 1
            write(error_unit,'("ASSERT_CLOSE vector failed: ",A,2X,"maxerr=",ES12.5,2X,"tol=",ES12.5)') trim(name), err, tol
        end if
    end subroutine

    subroutine finalize_tests()
        if (failures > 0) then
            write(error_unit,'("Total assertion failures: ",I0)') failures
            error stop 1
        else
            write(*,'("All assertions passed")')
        end if
    end subroutine

end module
