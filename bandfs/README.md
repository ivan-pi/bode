The procedures in this folder are adapted from the work

> Jeff Thorson, Gaussian elimination on a banded matrix
> https://sep.stanford.edu/data/media/public/oldreports/sep20/20_11_abs.html

The author Jeffrey R. Thorson was a member of the
Stanford Exploration Project and student (?) of Jon Claerbout, geophysics
professor at Stanford University. Thorson's Dissertation (1984) is catalogued
at: https://search.worldcat.org/en/title/38632165

The routines use column storage of the matrix, i.e. the diagonals are stored
as columns, similar to SLATEC.

The factorization is performed using the procedure,
```fortran
! Banded factorization
subroutine bandf(a,m,n,p,ifail)
  integer, intent(in) :: m, n
  real(kind=rk), intent(inout) :: a(n,*)
  integer, intent(out) :: p(n)
  integer, intent(out) :: ifail
```
where
* `a`, the values of the banded matrix; the matrix should have an additional `(m-1)/2` columns as storage for the factorization
* `m`, the bandwidth of the matrix (i.e. `m = ml + mu + 1`)
* `n`, the dimension of the system of linear equations
* `p`, the vector of pivots
* `ifail`, status flag, zero on success

Once factorized, the the solution of the system can be obtained by calling `bands`:
```fortran
! Banded back-substitution (aka solve)
subroutine bands(a,b,m,n,p)
   integer, intent(in) :: m, n
   real(kind=rk), intent(in) :: a(n,*)
   real(kind=rk), intent(inout) :: b(n)
   integer, intent(in) :: p(n)
```
Upon entry, `b` should contain the values of the right-hand side of the system `A x = b`. Upon exit, `b` contains the solution values `x`.

The routines can also be referred to by their specific names `sbandf/dbandf` and `sbands/dbands` for single and double precision, respectively.

