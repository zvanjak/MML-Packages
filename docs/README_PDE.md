# PDE Package

**Partial Differential Equation Solvers**

The PDE package provides a comprehensive framework for solving partial differential equations using finite difference methods, including sparse matrix infrastructure, iterative solvers, and specialized PDE solvers for elliptic, parabolic, and hyperbolic problems.

The ten Poisson (1D/2D/3D), Heat (1D/2D/3D), Wave (1D/2D), and Advection (1D/2D) solver classes are compiled into `MML_PDE`. Include the public `<pde/...>` headers and link `MML::PDE`; prebuilt solver instantiations are available for `float`, `double`, and `long double`, not arbitrary scalar types. The public grid, boundary-condition, result, sparse-matrix, and generic callable helper templates remain in headers, along with the MML core headers they depend on. Solver method definitions live in `src/pde/*.cpp`, which is not a consumer include directory.

The default build is static; configure with `-DBUILD_SHARED_LIBS=ON -DMML_STATIC=OFF` for a shared PDE library and make its DLL available to consumers on Windows. A prebuilt C++ library requires a compatible compiler, standard library, runtime, architecture, and ABI. Moving solver methods into the library does not hide inline helpers, public template code, or exported binary symbols, and does not provide complete IP protection. Packaging and `find_package` integration are separate work.

## Features

### Sparse Matrix Formats
- **COO (Coordinate)** - Construction format with efficient random insertion
- **CSR (Compressed Sparse Row)** - Efficient for sparse matrix-vector products
- **CSC (Compressed Sparse Column)** - Efficient for column operations

### Iterative Solvers
- **Conjugate Gradient (CG)** - For symmetric positive-definite systems
- **BiCGSTAB** - For general non-symmetric systems
- **GMRES** - Generalized Minimal Residual for general matrices

### Preconditioners
- **Jacobi** - Diagonal scaling (simple, parallel-friendly)
- **SSOR** - Symmetric Successive Over-Relaxation
- **ILU(0)** - Incomplete LU factorization

### Grid Infrastructure
- **1D/2D/3D Grids** - Uniform structured grids
- **Grid Functions** - Solutions defined on grids
- **Boundary Conditions** - Dirichlet, Neumann, Periodic, Robin
- **Stencils** - Finite difference stencil definitions

### PDE Solvers
- **Elliptic** - Poisson, Laplace equations (1D, 2D, 3D)
- **Parabolic** - Heat/diffusion equation with explicit, implicit, Crank-Nicolson methods
- **Hyperbolic** - Wave equation, advection (upwind, Lax-Friedrichs, Lax-Wendroff)

## Quick Start

See the [PDE quick-start guide](../../docs/pde/Quick_Start_Guide.md) and its runnable Poisson and Heat counterparts in `docs_demos/pde/example01_1d_poisson.cpp` and `docs_demos/pde/example02_1d_heat.cpp`. The standalone consumer in `tests/pde/consumer_link.cpp` exercises all ten compiled solver classes for each of the three supported scalar types.

## API Reference

### Sparse Matrix Classes

```cpp
// COO - Coordinate format (construction)
SparseMatrixCOO<T> coo(rows, cols);
coo.addEntry(i, j, value);
auto csr = coo.toCSR();
auto csc = coo.toCSC();

// CSR - Compressed Sparse Row (SpMV)
SparseMatrixCSR<T> A(coo);
std::vector<T> y = A * x;  // Matrix-vector product
A.apply(x, y);             // In-place: y = A*x

// CSC - Compressed Sparse Column
SparseMatrixCSC<T> B(coo);
```

### Iterative Solvers

| Solver | Use Case | Convergence |
|--------|----------|-------------|
| `ConjugateGradient<T>` | SPD matrices (Poisson, etc.) | O(√κ) |
| `BiCGSTAB<T>` | Non-symmetric systems | O(κ) |
| `GMRES<T>` | General matrices | O(n) in restart |

### Solver Configuration

```cpp
SolverConfig config;
config.maxIterations = 1000;    // Maximum iterations
config.relTolerance = 1e-8;     // Relative tolerance
config.absTolerance = 1e-12;    // Absolute tolerance
config.verbose = true;          // Print progress

solver.setConfig(config);
```

### Solver Results

```cpp
struct SolverResult<T> {
    SolverStatus status;      // Success, MaxIterations, Diverged
    int iterations;           // Iterations used
    T residualNorm;           // Final ||b - Ax||
    T relativeResidual;       // ||b - Ax|| / ||b||
    std::string message;      // Status message
};
```

### Preconditioners

| Preconditioner | Description | Cost per iteration |
|----------------|-------------|-------------------|
| `Jacobi` | Diagonal scaling | O(n) |
| `SSOR` | Symmetric SOR | O(nnz) |
| `ILU0` | Incomplete LU | O(nnz) |

## File Structure

```
pde/
├── README.md              # This file
├── CMakeLists.txt         # Build configuration
├── PDE.h                  # Main include header
│
├── include/
│   ├── sparse/            # Sparse matrix formats
│   │   ├── SparseMatrix.h      # Base class
│   │   ├── SparseMatrixCOO.h   # Coordinate format
│   │   ├── SparseMatrixCSR.h   # Compressed Row
│   │   └── SparseMatrixCSC.h   # Compressed Column
│   │
│   ├── solvers/           # Iterative linear solvers
│   │   ├── IterativeSolverBase.h  # Base class
│   │   ├── ConjugateGradient.h    # CG
│   │   ├── BiCGSTAB.h             # BiCGSTAB
│   │   ├── GMRES.h                # GMRES(m)
│   │   ├── Preconditioners.h      # Jacobi, SSOR, ILU
│   │   └── IterativeSolvers.h     # Convenience header
│   │
│   ├── grid/              # Grid infrastructure
│   │   ├── Grid.h               # Base grid concepts
│   │   ├── Grid1D.h             # 1D uniform grid
│   │   ├── Grid2D.h             # 2D uniform grid
│   │   ├── Grid3D.h             # 3D uniform grid
│   │   ├── GridFunction.h       # Functions on grids
│   │   ├── BoundaryConditions.h # BC types
│   │   └── Stencils.h           # FD stencils
│   │
│   ├── elliptic/          # Elliptic PDE solvers
│   │   └── PoissonSolver.h  # Poisson (1D, 2D, 3D)
│   │
│   ├── parabolic/         # Parabolic PDE solvers
│   │   └── HeatSolver.h     # Heat equation
│   │
│   ├── hyperbolic/        # Hyperbolic PDE solvers
│   │   ├── WaveSolver.h       # Wave equation
│   │   └── AdvectionSolver.h  # Advection equation
│   │
│   └── testbeds/          # Validation tools
│       ├── ManufacturedSolutions.h
│       └── BenchmarkProblems.h
│
└── tests/
    ├── sparse_tests.cpp
    ├── solver_tests.cpp
    └── pde_tests.cpp
```

## Mathematical Background

### Poisson Equation

Solves: $-\nabla^2 u = f$ in domain $\Omega$ with boundary conditions on $\partial\Omega$

**2D 5-point stencil:**
$$\nabla^2 u \approx \frac{u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1} - 4u_{i,j}}{h^2}$$

### Heat Equation

Solves: $\frac{\partial u}{\partial t} = \alpha \nabla^2 u$

**Time schemes:**
- Explicit (Forward Euler): Conditionally stable, $\Delta t < h^2/(2\alpha d)$
- Implicit (Backward Euler): Unconditionally stable, 1st order
- Crank-Nicolson: Unconditionally stable, 2nd order

### Wave Equation

Solves: $\frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u$

**Leapfrog scheme (CFL stable for $c\Delta t < h$):**
$$u^{n+1}_i = 2u^n_i - u^{n-1}_i + \frac{c^2\Delta t^2}{h^2}(u^n_{i+1} - 2u^n_i + u^n_{i-1})$$

## Performance Tips

1. **Matrix format selection:**
   - Use COO for construction, CSR for solving
   - CSR is optimal for row-major SpMV

2. **Preconditioner selection:**
   - Jacobi: Simple, good for well-conditioned
   - SSOR: Good for Laplacian-like problems
   - ILU(0): Best convergence, higher setup cost

3. **Grid resolution:**
   - Error ~ O(h²) for 2nd order FD
   - Doubling resolution increases unknowns by 4× (2D) or 8× (3D)

4. **Solver tuning:**
   - Set appropriate tolerance (usually 1e-8 to 1e-12)
   - Monitor residual for convergence issues

## References

- Saad, Y. (2003). "Iterative Methods for Sparse Linear Systems" (2nd ed.)
- Strikwerda, J.C. (2004). "Finite Difference Schemes and Partial Differential Equations"
- LeVeque, R.J. (2007). "Finite Difference Methods for Ordinary and Partial Differential Equations"
- Numerical Recipes in C++, 3rd Edition, Chapter 20: Less-Numerical Algorithms

## See Also

- [mml_packages/README.md](../README.md) - Package overview
- [Optimization Package](../optimization/README.md) - Optimization of PDE parameters
- [Systems Package](../systems/README.md) - ODE integration for time-dependent problems
