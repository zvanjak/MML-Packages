///////////////////////////////////////////////////////////////////////////////////////////
// PDE.h - Main Header for PDE Solver Module
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary
//
// Provides a single include for all PDE solver functionality:
//   - Sparse matrix formats (COO, CSR, CSC)
//   - Iterative linear solvers (CG, BiCGSTAB, GMRES)
//   - Preconditioners (Jacobi, SSOR, ILU)
//   - Grid infrastructure (1D, 2D, 3D grids, boundary conditions)
//   - Elliptic PDE solvers (Poisson, Laplace)
//   - Parabolic PDE solvers (Heat/Diffusion equation)
//   - Hyperbolic PDE solvers (Wave, Advection equations)
//
// Usage:
//   #include <pde/PDE.h>
//   using namespace MML::PDE;
//
// For selective includes, use individual headers:
//   #include <mml/base/SparseMatrix/SparseMatrix.h>
//   #include <mml/core/SparseSolvers/ConjugateGradient.h>
//   #include <pde/elliptic/PoissonSolver.h>
//   etc.
//
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_H
#define MML_PDE_H

// ============================================================================
// Sparse Matrix Infrastructure
// ============================================================================
// Efficient storage formats for large sparse matrices arising from PDE discretization

#include <mml/base/SparseMatrix/SparseMatrix.h>
#include <mml/base/SparseMatrix/SparseMatrixCOO.h>
#include <mml/base/SparseMatrix/SparseMatrixCSR.h>
#include <mml/base/SparseMatrix/SparseMatrixCSC.h>

// ============================================================================
// Iterative Linear Solvers
// ============================================================================
// Krylov subspace methods for solving Ax = b

#include <mml/core/SparseSolvers/IterativeSolverBase.h>
#include <mml/core/SparseSolvers/ConjugateGradient.h>
#include <mml/core/SparseSolvers/BiCGSTAB.h>
#include <mml/core/SparseSolvers/GMRES.h>
#include <mml/core/SparseSolvers/Preconditioners.h>
#include <mml/core/SparseSolvers/IterativeSolvers.h>

// ============================================================================
// Grid Infrastructure
// ============================================================================
// Structured grids, boundary conditions, and grid functions

#include "grid/Grid.h"                // Grid base concepts
#include "grid/Grid1D.h"              // 1D uniform grid
#include "grid/Grid2D.h"              // 2D uniform grid
#include "grid/Grid3D.h"              // 3D uniform grid
#include "grid/GridFunction.h"        // Functions defined on grids
#include "grid/BoundaryConditions.h"  // Dirichlet, Neumann, Periodic BCs
#include "grid/Stencils.h"            // Finite difference stencils

// ============================================================================
// Elliptic PDE Solvers
// ============================================================================
// Steady-state boundary value problems: -∇²u = f

#include "elliptic/PoissonSolver.h"   // Poisson equation solver (1D, 2D, 3D)

// ============================================================================
// Parabolic PDE Solvers
// ============================================================================
// Time-dependent diffusion: ∂u/∂t = α∇²u

#include "parabolic/HeatSolver.h"     // Heat/diffusion equation (1D, 2D, 3D)

// ============================================================================
// Hyperbolic PDE Solvers
// ============================================================================
// Wave propagation and transport: ∂²u/∂t² = c²∇²u, ∂u/∂t + c·∇u = 0

#include "hyperbolic/WaveSolver.h"      // Wave equation (1D, 2D)
#include "hyperbolic/AdvectionSolver.h" // Advection equation (1D, 2D)

// ============================================================================
// Testing and Validation
// ============================================================================
// Manufactured solutions and benchmark problems for solver validation

#include "testbeds/ManufacturedSolutions.h"  // Analytical solutions for testing
#include "testbeds/BenchmarkProblems.h"      // Standard benchmark problems

// ============================================================================
// Namespace Convenience
// ============================================================================

namespace MML {
namespace PDE {

/**
 * @brief PDE Module version information
 */
struct PDEVersion {
    static constexpr int major = 1;
    static constexpr int minor = 0;
    static constexpr int patch = 0;
    static constexpr const char* string = "1.0.0";
};

/**
 * @brief Print PDE module capabilities
 */
void printCapabilities();

} // namespace PDE
} // namespace MML

#endif // MML_PDE_H
