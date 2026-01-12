# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.7.0] - 2026-01-12

### Added
- **Native Sparse Matrix Support**: Optimized several quantum functionals for large sparse matrices using iterative solvers from `Arpack.jl`. This includes:
  - `norm_trace`
  - `fidelity_sqrt`
  - `vonneumann_entropy`
  - `renyi_entropy`
  - `relative_entropy`
  - `negativity`
  - `ppt`
  - `concurrence` (via dense conversion as per qubit definition)
- **Dependency**: Added `Arpack.jl` for sparse eigenvalue and singular value decomposition.
- **Constructors**: Added typeless `SuperOperator(::Function, ::Int, ::Int)` constructor.

### Fixed
- **Restored `𝕀` function**: Brought back the identity gate function in `src/gates.jl`.
- **Type Stability & JET Warnings**: Resolved multiple type inference issues, particularly with `Hermitian` wrappers on sparse matrices.
- **Constructor Validation**: Fixed typeless constructors for `KrausOperators` and `POVMMeasurement` to correctly pass and validate dimensions.
- **Sparse `ptrace` and `permutesystems`**: Improved efficiency and fixed edge cases for sparse matrix operations.

### Changed
- **Test Suite Refactoring**: New tests distributed across module-specific files for better organization.
- **Increased Test Coverage**: Expanded the test suite with comprehensive error handling, predicate edge cases, and large matrix scenarios (461+ tests passing).
- **Default Iterative Solver Threshold**: Set to 32x32 for most sparse functional implementations for optimal performance.
