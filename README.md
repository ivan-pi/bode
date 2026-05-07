# bode

A low-order adaptive solver for implicit ODE systems,

\[
B\,\frac{dy}{dx} = F(x,y)
\]

with support for structurally-symmetric banded Jacobians.

## Behavior freeze milestone

Implementation code under `/src` is currently under a **behavior freeze** until the new test suite is fully established and green in CI.
During this milestone, only the following areas should change:

- tests (`/test`)
- build and packaging (`CMakeLists.txt`, `/cmake`)
- CI/CD (`.github/workflows`)
- documentation (`README`, `/docs`, `ford.md`)
- examples (`/examples`)

## Supported toolchain and dependencies

- Compilers: `gfortran` (primary; tested in CI)
- Build system: CMake (canonical)
- Required numerical libraries: BLAS + LAPACK
- Optional docs generator: [FORD](https://github.com/Fortran-FOSS-Programmers/ford)

## Quick start

```bash
cmake -S . -B build -DBUILD_TESTING=ON -DBODE_BUILD_EXAMPLES=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## API overview

Public API is exposed by `bode_mod`:

- `bode(...)`: adaptive implicit ODE integrator
- `ltri(...)`, `tsol(...)`: internal banded factor/solve helpers used by the solver
- Runtime counters: `nfev`, `njev`, `nlu`, `nbsol`

Callback contracts used by `bode`:

- `pmult(p,n,q)`: applies the mass matrix action `q = B p`
- `deriv(y,n,q,x,h)`: computes scaled RHS `q = h * F(x,y)`
- `pjacb(x,y,n,jac,ld,m1,h)`: computes Jacobian form expected by selected mode
- `monit(...)`: monitoring callback called during integration

For generated API pages, see **API documentation** below.

## Tests

CTest is the standard entrypoint:

```bash
ctest --test-dir build --output-on-failure
```

Test categories:

- `unit`: linear algebra and helper routines
- `regression`: fixed benchmark trajectories/tolerances
- `failure`: negative/invalid-path behavior

Reference data used by regression tests lives under `/test/data`.

## Examples

Examples are documented in `/examples/README.md`.

Build and run one example:

```bash
cmake --build build --target example_getting_started
./build/example_getting_started
```

## API documentation

Generate with FORD:

```bash
ford ford.md
```

Output is written to `build/docs/ford`.

## CI/CD

Workflows are defined in `.github/workflows`:

- `ci.yml`: build + test for push/PR
- `docs.yml`: FORD docs build and GitHub Pages deploy (main branch)
- `release.yml`: tagged source archive packaging

## Historical context

The solver has been extracted from:

> Hopkins, T. R., & Wait, R. (1978). A comparison of Galerkin, collocation and the method of lines for partial differential equations. International Journal for Numerical Methods in Engineering, 12(7), 1081-1107. https://doi.org/10.1002/nme.1620120703

The method used appeared originally in:

> Prothero, A., & Robinson, A. (1974). On the stability and accuracy of one-step methods for solving stiff systems of ordinary differential equations. Mathematics of Computation, 28(125), 145-162. https://doi.org/10.1090/S0025-5718-1974-0331793-2
