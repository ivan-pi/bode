# Examples

This directory contains executable examples for the `bode` solver.

## Getting started: scalar linear decay

- **Problem statement:** solve `y' = -y`, `y(0)=1` on `[0,1]` with `B = I`.
- **Expected output characteristics:** numerical value near `exp(-1)` with small integration error.
- **Run command:** `cmake --build build --target example_getting_started && ./build/example_getting_started`
- **Interpretation notes:** demonstrates minimum callback set (`pmult`, `deriv`, `pjacb`, `monit`) and `iopt` usage.

## Van der Pol oscillator (`test_van_der_pol.f90`)

- **Problem statement:** stiff Van der Pol system with `mu=10`.
- **Expected output characteristics:** bounded oscillatory trajectory with adaptive step control.
- **Run command:** `cmake --build build --target example_vdp && ./build/example_vdp`
- **Interpretation notes:** demonstrates banded Jacobian callback for a 2x2 system.

## Robertson kinetics (`robertson.F90`)

- **Problem statement:** stiff chemical kinetics benchmark with 3 species.
- **Expected output characteristics:** positivity-preserving concentration evolution and long-time mass balance.
- **Run command:** `cmake --build build --target example_robertson && ./build/example_robertson`
- **Interpretation notes:** includes dense Jacobian path and integration statistics.

## Advanced: dense vs banded Jacobian (`dense_vs_banded_robertson.f90`)

- **Problem statement:** compare dense and banded Jacobian paths on the same Robertson setup.
- **Expected output characteristics:** both paths converge to nearly identical final states.
- **Run command:** `cmake --build build --target example_dense_vs_banded && ./build/example_dense_vs_banded`
- **Interpretation notes:** use the reported difference norm to validate consistency between Jacobian modes.
