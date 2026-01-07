# bandfs

This folder contains procedures for banded LU-factorization and back-substitution,
templated for single and double precision `real` kinds.

The procedures are adapted from the complex-type subroutines originally
developed by:

> Jeff Thorson, Gaussian elimination on a banded matrix, Stanford Exploration
> Project: Report No. 20, October 1979,
> https://sep.stanford.edu/data/media/public/oldreports/sep20/20_11_abs.html

The original complex routines can be found in the [`complex/`](./complex)
sub-folder.

---

**Contents:**
* [Usage](./README.md#usage)
* [Procedure Reference](./README.md#procedure-reference)
  - [Banded Storage](./README.md#banded-storage)
  - [bandf](./README.md#bandf)
  - [bands](./README.md#bands)
  - [bndmv](./README.md#bndmv)
  - [bndsl](./README.md#bndsl)
* [Historical Note](./README.md#historical-note)

# Usage

The following program demonstrates how to solve a 5 × 5 tridiagonal system
using the procedured `bandf` and `bands`.
Note how the custom arrays bounds are used to center the matrix in column 0.
The array is padded with `kh = (k-1)/2` columns necessary for pivoting.

```fortran
program bandfs_demo
implicit none

integer, parameter :: wp = kind(1.0d0) ! Double precision

! BANDFS procedure interfaces
include "bandfs.fi"

integer, parameter :: n = 5, k = 3  ! Dimension and bandwidth
integer, parameter :: kh = (k-1)/2   ! Half-bandwidth

! Array storage
integer, parameter :: lda = n
real(wp) :: a(lda,-kh:kh+kh), b(n)
integer  :: ipiv(n), info


! 1. Initialize array and RHS
! Gilbert Strang's favorite matrix (second order differences)
a(:,-1) = [ 0,-1,-1,-1,-1]  ! Lower diagonal
a(:, 0) = [ 2, 2, 2, 2, 2]  ! Main diagonal
a(:, 1) = [-1,-1,-1,-1, 0]  ! Upper diagonal

! Right-hand side for a target solution of x = [1, 2, 3, 4, 5]
b = [ 0, 0, 0, 0, 6]

! 2. Perform LU Factorization
! 'a' is modified in-place; 'ipiv' stores pivoting sequence
call bandf(n, k, a, lda, ipiv, info)

if (info /= 0) then
   write(*,'(A)') "Error: Matrix is singular or factorization failed."
   stop 1
end if

! 3. Solve the system
! 'b' is overwritten with the solution vector 'x'
call bands(n, k, a, lda, ipiv, b, info)

write(*,'(A,/,*(2X,F6.3,:,/))') "Solution x:", b

end program
```

# Procedure Reference

## Banded Storage

The diagonals of matrix $A$ are stored in array `A(1:n,1:k)` where `n`
is the dimension of the system and `k` is the bandwidth (`k` must be odd).

The matrix is stored in row-wise manner, with the diagonals of matrix $A$
loaded into columns of array `A`. Thus, element $A_{ij}$ is to be loaded
into element `A(i,j-i+(k-1)/2)`. **Note:** this row-wise storage order differs
from the column-wise format used in the LINPACK and LAPACK libraries.

The procedures `bandf` and `bands` are limited to _structurally_ symmetric
banded matrices,where the number of lower and upper diagonals are equal:
`kl = ku = (k-1)/2`.

Due to the fill-in during factorization, the second dimension of `A` must
be at least `k + (k-1)/2`. The extra `(k-1)/2` columns of padding provide
workspace for the interchanged rows during pivoting.

The same storage scheme is used by the SLATEC library procedures [`dnbfa`](https://netlib.org/slatec/src/dnbfa.f)/[`dnbsl`](https://netlib.org/slatec/src/dnbsl.f),
however the SLATEC procedures do not impose any symmetry restrictions.
(In theory, the procedures introduced below could just wrap the SLATEC procedures.)

### Examples

**Example 1 (tridiagonal):** If the original matrix is
```
11 22  0  0  0
21 22 23  0  0
 0 32 33 34  0
 0  0 43 44 45
 0  0  0 54 55
```
then `n = 5`, `k = 3`, and the array `A` should contain
```
 * 11 12 +     , * = not used
21 22 23 +     , + = used for pivoting
32 33 34 +
43 44 45 +
54 55  * +
```

**Example 2 (pentadiagonal):** If the original matrix is,
```
11 12 13  0  0  0
21 22 23 24  0  0
 0 32 33 34 35  0
 0  0 43 44 45 46
 0  0  0 54 55 56
 0  0  0 64 65 66
```
then `n = 6`, `k = 5`, and the array `A` should contain
```
 *  * 11 12 13  +  +      , * = not used
 * 21 22 23 24  +  +      , + = used for pivoting
 0 32 33 34 35  +  +
 0 43 44 45 45  +  +
 0 54 55 56  *  +  +
64 65 66  *  *  +  +
```

## `bandf()`

Banded LU factorization with partial pivoting

```fortran
subroutine bandf(n,k,a,lda,ipiv,info)
  integer, intent(in) :: n, k, lda
  real(kind=[sp,dp]), intent(inout) :: a(lda,*)
  integer, intent(out) :: p(n)
  integer, intent(out) :: info
```

Arguments:
* `n`: The dimension of the system (number of rows/columns).
* `k`: The full bandwidth of the matrix (total number of diagonals).
* `a`: The values of the banded matrix. The second dimension must be at least `k + (k-1)/2`. On output contains the factorized matrix.
* `lda`: Leading dimension of array `a`. `lda >= n`.
* `ipiv`: The pivot indices.
* `info`: Error flag (0 on success; if positive, the matrix is singular; if negative with `info = -i`, the `i`-th argument had an illegal value).

Available as generic `bandf` or specific names `sbandf`/`dbandf`.

## `bands()`

Solve the banded system $Ax = b$ using the factorization generated by `bandf`.

```fortran
subroutine bands(n,k,a,lda,ipiv,info)
  integer, intent(in) :: n, k, lda
  real(kind=[sp,dp]), intent(in) :: a(lda,*)
  integer, intent(in) :: ipiv(n)
  real(kind=[sp,dp]), intent(inout) :: b(n)
  integer, intent(out) :: info
```

Arguments:
- `n`: The dimension of the system (number of rows).
- `k`: The full bandwidth of the matrix (total number of diagonals).
- `a`: The factorized matrix returned by `bandf`.
- `lda`: Leading dimension of the array `a`. `lda >= n`.
- `ipiv`: Array of pivot indices from `bandf`.
- `b`: Upon entry, the right-hand side. Upon exit, the solution vector $x$.
- `info`: Error flag (0 on success; if negative with `info = -i`, the `i`-th argument had an illegal value).

Available as generic `bands` or specific names `sbands`/`dbands`.

## `bndmv()`

Performs the matrix-vector operation, `y := alpha*A*x + beta*y`

```fortran
subroutine bndmv(n, kl, ku, alpha, a, lda, x, beta, y)
  integer,  intent(in)    :: n, kl, ku, lda
  real(kind=[sp,dp]), intent(in)    :: alpha, beta
  real(kind=[sp,dp]), intent(in)    :: a(lda,-kl:ku), x(n)
  real(kind=[sp,dp]), intent(inout) :: y(n)
```

Arguments:
- `n`: The number of rows/columns of the matrix `A`. `N >= 0`.
- `kl`: The number of sub-diagonals of the matrix `A`. `KL >= 0`.
- `ku`: The number of super-diagonals of the matrix `A`. `KU >= 0`.
- `alpha`: The scalar alpha.
- `a`: The banded matrix in column storage.
- `lda`: The first dimension of `a` as declaring in the calling (sub)program. `LDA >= N`.
- `x`: The vector x.
- `beta`: The scalar beta.
- `y`: The vector y.

Available as generic `bands` or specific names `sbndmv`/`dbndmv`.

## `bndsl()`

Compute the solution to the system of linear equations `A * X = B` for structurally symmetric banded matrices.

```fortran
subroutine bndsl(n,k,nrhs,a,lda,ipiv,b,ldb,info)
  integer, intent(in) :: n, k, nrhs, lda, ldb
  real(kind=[sp,dp]), intent(inout) :: a(lda,*)
  integer, intent(out) :: ipiv(n)
  real(kind=[sp,dp]), intent(inout) :: b(ldb,nrhs)
  integer, intent(out) :: info
```

Arguments:
- `n`: Dimension of the system (number of rows/columns).
- `k`: The full bandwidth of the matrix (total number of diagonals).
- `nrhs`: Number of right-hand sides.
- `a`: The banded matrix in column storage.
- `lda`: Leading dimension of the array `a`. `lda >= n`.
- `ipiv`: The pivot indices that define the permutation matrix.
- `b`: On entry, the `n`-by-`nrhs` right hand side matrix `b`. On exit, if `info = 0`, the `n`-by-`nrhs` solution matrix $x$.
- `ldb`: The leading dimension of the array `b`. `ldb >= max(1,n)`.
- `info`: Error flag (0 on success; if positive, the matrix is singular; if negative with `info = -i`, the `i`-th argument had an illegal value).

Available as generic `bndsl` or specific names `sbndsl`/`dbndsl`.


# Historical Note

According to various internet sources, the author Jeffrey R. Thorson was a
member of the Stanford Exploration Project (SEP) led by [Jon Claerbout](https://sep.stanford.edu/sep/jon/) (geophysics professor at Stanford University).
Thorson's dissertation (1984) appears in the [WorldCat](https://search.worldcat.org/en/title/38632165) catalogue.

The work on Gaussian elimination is part on an [old report](https://sep.stanford.edu/data/media/public/oldreports/sep20/).