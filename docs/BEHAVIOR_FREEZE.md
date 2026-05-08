# Behavior Freeze Policy

## Temporary milestone objective

Before any implementation changes in `/src`, the repository must maintain a strong safety net:

1. Deterministic unit + regression + failure-path tests
2. CI that blocks merges on test failures
3. Baseline reference outputs tracked under `test/data`
4. Updated user/developer documentation

## Allowed changes during freeze

- tests, test data, and harnesses
- CMake/build/packaging
- CI/CD workflows
- docs and examples

## Explicitly disallowed during freeze

- algorithmic/behavior changes in solver implementation files under `/src`

## Exit criteria

The freeze can be lifted only when all gates are satisfied:

- Gate A: comprehensive tests are green
- Gate B: CI enforces build + tests
- Gate C: API docs generated from source comments
- Gate D: examples are documented and runnable
