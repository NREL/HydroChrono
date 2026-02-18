---
layout: page
title: Unit Tests
parent_section: Developer Documentation
---

# Unit Tests

Unit tests verify internal correctness of HydroChrono and its integration with Project Chrono at the C++ level. They are self-contained (no external data files required) and typically run in under a second.

## Running unit tests

```powershell
ctest -C Release -L unit --test-dir build
```

Useful flags:

| Flag | Description |
|------|-------------|
| `-V` | Verbose — show full test output (matrix dumps, signatures) |
| `--output-on-failure` | Only show output when a test fails |
| `-R <pattern>` | Run only tests matching a name pattern |

Example with verbose output:
```powershell
ctest -C Release -L unit --test-dir build -V
```

## Test: added-mass assembly (`test_added_mass_determinism`)

Verifies that Chrono's `ChLoadHydrodynamics` (the default added-mass path) correctly assembles a cross-coupled added-mass matrix for a multi-body system.

### What it does

1. **Constructs** a fully dense, symmetric positive-definite 18×18 added-mass matrix for 3 bodies, with non-trivial off-diagonal coupling between all body pairs and all 6 DOFs
2. **Passes** the per-body 6×18 row blocks to `ChLoadHydrodynamics` (the same path HydroChrono uses at runtime)
3. **Extracts** the assembled system mass matrix from Chrono, subtracts body inertia, and asserts it matches the original input
4. **Runs 12 independent trials** with fresh system instances and heap perturbation between them, asserting bit-identical velocity results across all trials
5. **Confirms** the added mass is actively influencing the dynamics (non-trivial velocity change after one timestep)
6. **Sweeps across Chrono solver types** and reports a summary table showing which solvers correctly preserve the full added-mass matrix and produce consistent dynamics

### What it catches

- Incorrect scatter of added-mass blocks into the system matrix (e.g., wrong row/column placement)
- Ordering sensitivity from internal container changes (e.g., `unordered_map` vs `vector`)
- Non-deterministic behaviour across independent system instances
- Silent corruption of off-diagonal cross-coupling terms
- Solver-specific issues with added-mass handling (e.g., solvers that cannot handle stiffness/damping matrices, or iterative solvers with convergence differences)

### Solver sweep

The test runs the same added-mass assembly and single-step dynamics with every available Chrono solver type and prints a summary table:

| Column | Meaning |
|--------|---------|
| Assembly | Whether the extracted added-mass matrix matches the input |
| NaN/Inf | Whether velocities contain invalid values |
| Dynamics | Whether post-step velocities match the SPARSE_QR reference |

The solver sweep is **informational only** — it does not affect the test PASS/FAIL status. Some solvers may not support this problem type (e.g., PSOR rejects systems with stiffness/damping matrices) or may converge to slightly different results (e.g., APGD). The table makes these differences visible without false positives.

### Reference results

The table below summarises the results observed when the test was last run. Run the test with `-V` to regenerate it against your current Chrono build.

| Solver | Assembly | NaN/Inf | Dynamics | Notes |
|--------|----------|---------|----------|-------|
| SPARSE_QR | MATCH | clean | MATCH | Reference solver used by HydroChrono |
| SPARSE_LU | MATCH | clean | MATCH | |
| MINRES | MATCH | clean | MATCH | |
| GMRES | MATCH | clean | MATCH | |
| BICGSTAB | MATCH | clean | MATCH | |
| PSOR | ERROR | — | — | Cannot handle stiffness/damping matrices |
| BARZILAIBORWEIN | MATCH | clean | MATCH | |
| APGD | MATCH | clean | diff ~6.7e-3 | Assembly correct; iterative convergence differs slightly |

**Key takeaways:**

- **Assembly is solver-independent.** Every solver that can run produces an added-mass matrix that matches the input exactly, including all off-diagonal cross-coupling terms. Choosing a different solver does not lose structural information from the mass matrix.
- **PSOR is incompatible.** `ChLoadHydrodynamics` contributes stiffness/damping matrices that PSOR cannot handle. This is a Chrono limitation, not an assembly issue. HydroChrono should never be configured to use PSOR.
- **APGD assembles correctly but converges differently.** The ~6.7e-3 velocity difference is a convergence characteristic of the iterative APGD solver, not a mass-matrix problem. If tighter agreement is needed, use a direct solver.
- **All other solvers (direct and iterative) produce bit-identical dynamics.**

### Verbose output

When run with `-V`, the test prints the input blocks, the assembled matrix extracted from Chrono, eigenvalues, post-step velocity signatures, and the solver sweep summary table — making it straightforward to visually inspect what went into the solver and what came out.
