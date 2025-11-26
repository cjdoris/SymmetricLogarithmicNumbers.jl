# Agent Instructions

- The project name is **SymmetricLogarithmicNumbers.jl**; the exported number type is **SymLog{T}** (with float aliases `SymLogF16`, `SymLogF32`, `SymLogF64`).
- Run tests with `julia --project -e 'using Pkg; Pkg.test()'`.
- Julia compatibility is targeted at **1.10** (current LTS), and CI runs tests on Julia `1` and `1.10` on Linux via GitHub Actions.
- Keep this file updated with any new findings or corrections for future agents. If you discover an incorrect assumption, add the corrected information here.
