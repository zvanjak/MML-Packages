///////////////////////////////////////////////////////////////////////////////////////////
// PoissonSolver.h - Finite Difference Solver for Poisson and Laplace Equations
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Solves elliptic PDEs of the form:
//   ∇²u = f(x)     in Ω       (Poisson equation)
//   ∇²u = 0        in Ω       (Laplace equation)
//
// with Dirichlet, Neumann, or Robin boundary conditions.
//
// Features:
// - 1D, 2D, and 3D Poisson solvers
// - Standard 5-point (2D) and 7-point (3D) stencils
// - Integration with iterative solvers (CG, BiCGSTAB, GMRES)
// - Manufactured solution validation
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_POISSON_SOLVER_H
#define MML_PDE_POISSON_SOLVER_H

#include <mml/mml_export.h>
#include <functional>
#include <cmath>
#include <iostream>

#include <mml/base/SparseMatrix/SparseMatrixCOO.h>
#include <mml/base/SparseMatrix/SparseMatrixCSR.h>
#include <mml/core/SparseSolvers/IterativeSolvers.h>
#include "../grid/Grid.h"
#include "../grid/Grid1D.h"
#include "../grid/Grid2D.h"
#include "../grid/Grid3D.h"
#include "../grid/GridFunction.h"
#include "../grid/BoundaryConditions.h"

namespace MML::PDE
{
    using MML::SparseMatrix::SparseMatrixCOO;
    using MML::SparseMatrix::SparseMatrixCSR;
    using namespace MML::SparseSolvers;

    ///////////////////////////////////////////////////////////////////////
    ///                      1D Poisson Solver                          ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 1D Poisson equation: -u'' = f(x)
     * 
     * Uses second-order central differences:
     *   -u''(x_i) ≈ (-u_{i-1} + 2u_i - u_{i+1}) / h²
     * 
     * The negative sign convention gives a positive-definite matrix.
     * 
     * @tparam T Numeric type (float, double)
     */
    template<typename T>
    class PoissonSolver1D {
    private:
        Grid1D<T> grid_;
        BoundaryConditions1D<T> bc_;
        std::function<T(T)> source_;
        
    public:
        /**
         * @brief Construct solver for given grid and boundary conditions
         */
        PoissonSolver1D(const Grid1D<T>& grid, const BoundaryConditions1D<T>& bc);
        
        /**
         * @brief Set the source term f(x) in -∇²u = f
         */
        void setSource(std::function<T(T)> f);
        
        /**
         * @brief Assemble the discrete Laplacian matrix and RHS vector
         * 
         * Creates the matrix A and vector b for the system Au = b
         * where A is the negative discrete Laplacian with BCs applied.
         */
        std::pair<SparseMatrixCSR<T>, std::vector<T>> assemble() const;
        
        /**
         * @brief Solve the Poisson equation
         * 
         * Uses CG for symmetric systems (pure Dirichlet BCs) or BiCGSTAB
         * for non-symmetric systems (Neumann/Robin BCs make the matrix asymmetric).
         * 
         * @param tol Convergence tolerance for iterative solver
         * @param maxIter Maximum iterations
         * @return Solution as a GridFunction1D
         */
        GridFunction1D<T> solve(T tol = T(1e-10), int maxIter = 10000) const;
        
        /**
         * @brief Get the discretization matrix (for analysis/debugging)
         */
        SparseMatrixCSR<T> getMatrix() const;
        
        const Grid1D<T>& grid() const;
    };

    extern template class MML_PDE_TEMPLATE_API PoissonSolver1D<float>;
    extern template class MML_PDE_TEMPLATE_API PoissonSolver1D<double>;
    extern template class MML_PDE_TEMPLATE_API PoissonSolver1D<long double>;

    ///////////////////////////////////////////////////////////////////////
    ///                      2D Poisson Solver                          ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 2D Poisson equation: -∇²u = f(x,y)
     * 
     * Uses the standard 5-point stencil:
     *   -∇²u(x_i, y_j) ≈ (4u_{i,j} - u_{i-1,j} - u_{i+1,j} - u_{i,j-1} - u_{i,j+1}) / h²
     * 
     * For non-square grids (dx ≠ dy):
     *   -∇²u ≈ -(u_{i-1,j} - 2u_{i,j} + u_{i+1,j})/dx² 
     *          -(u_{i,j-1} - 2u_{i,j} + u_{i,j+1})/dy²
     * 
     * @tparam T Numeric type
     */
    template<typename T>
    class PoissonSolver2D {
    private:
        Grid2D<T> grid_;
        BoundaryConditions2D<T> bc_;
        std::function<T(T, T)> source_;
        
    public:
        PoissonSolver2D(const Grid2D<T>& grid, const BoundaryConditions2D<T>& bc);
        
        void setSource(std::function<T(T, T)> f);
        
        /**
         * @brief Assemble the discrete Laplacian matrix and RHS
         * 
         * Matrix ordering: lexicographic (row-major)
         *   idx = i + j * nx  where i is x-index, j is y-index
         */
        std::pair<SparseMatrixCSR<T>, std::vector<T>> assemble() const;
        
        /**
         * @brief Solve the Poisson equation
         * 
         * Uses CG for symmetric systems (pure Dirichlet BCs) or BiCGSTAB
         * for non-symmetric systems (Neumann/Robin BCs make the matrix asymmetric).
         */
        GridFunction2D<T> solve(T tol = T(1e-10), int maxIter = 10000) const;
        
        /**
         * @brief Solve and return detailed results
         */
        struct SolveResult {
            GridFunction2D<T> solution;
            bool converged;
            int iterations;
            T residual;
        };
        
        SolveResult solveDetailed(T tol = T(1e-10), int maxIter = 10000) const;
        
        SparseMatrixCSR<T> getMatrix() const;
        
        const Grid2D<T>& grid() const;
    };

    extern template class MML_PDE_TEMPLATE_API PoissonSolver2D<float>;
    extern template class MML_PDE_TEMPLATE_API PoissonSolver2D<double>;
    extern template class MML_PDE_TEMPLATE_API PoissonSolver2D<long double>;

    ///////////////////////////////////////////////////////////////////////
    ///                      3D Poisson Solver                          ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 3D Poisson equation: -∇²u = f(x,y,z)
     * 
     * Uses the standard 7-point stencil.
     */
    template<typename T>
    class PoissonSolver3D {
    private:
        Grid3D<T> grid_;
        BoundaryConditions3D<T> bc_;
        std::function<T(T, T, T)> source_;
        
    public:
        PoissonSolver3D(const Grid3D<T>& grid, const BoundaryConditions3D<T>& bc);
        
        void setSource(std::function<T(T, T, T)> f);
        
        std::pair<SparseMatrixCSR<T>, std::vector<T>> assemble() const;
        
        GridFunction3D<T> solve(T tol = T(1e-10), int maxIter = 10000) const;
        
        const Grid3D<T>& grid() const;
    };

    extern template class MML_PDE_TEMPLATE_API PoissonSolver3D<float>;
    extern template class MML_PDE_TEMPLATE_API PoissonSolver3D<double>;
    extern template class MML_PDE_TEMPLATE_API PoissonSolver3D<long double>;

    ///////////////////////////////////////////////////////////////////////
    ///                   Manufactured Solutions                        ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Collection of manufactured solutions for testing
     * 
     * These are exact solutions to Poisson's equation with known
     * source terms, used for validating numerical accuracy.
     */
    namespace ManufacturedSolutions {
        
        // 1D: u(x) = sin(πx), -u'' = π²sin(πx), u(0)=u(1)=0
        template<typename T>
        struct SinPi1D {
            static T solution(T x) { 
                return std::sin(M_PI * x); 
            }
            static T source(T x) { 
                return M_PI * M_PI * std::sin(M_PI * x); 
            }
            static T leftBC() { return T(0); }
            static T rightBC() { return T(0); }
        };
        
        // 1D: u(x) = x(1-x), -u'' = 2, u(0)=u(1)=0
        template<typename T>
        struct Parabola1D {
            static T solution(T x) { 
                return x * (T(1) - x); 
            }
            static T source(T) { 
                return T(2); 
            }
            static T leftBC() { return T(0); }
            static T rightBC() { return T(0); }
        };
        
        // 2D: u(x,y) = sin(πx)sin(πy), -∇²u = 2π²sin(πx)sin(πy)
        template<typename T>
        struct SinPi2D {
            static T solution(T x, T y) { 
                return std::sin(M_PI * x) * std::sin(M_PI * y); 
            }
            static T source(T x, T y) { 
                return T(2) * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y); 
            }
            static T boundaryValue(T, T) { return T(0); }
        };
        
        // 2D: u(x,y) = x(1-x)y(1-y), -∇²u = 2[y(1-y) + x(1-x)]
        template<typename T>
        struct Parabola2D {
            static T solution(T x, T y) {
                return x * (T(1) - x) * y * (T(1) - y);
            }
            static T source(T x, T y) {
                return T(2) * (y * (T(1) - y) + x * (T(1) - x));
            }
            static T boundaryValue(T, T) { return T(0); }
        };
        
        // 2D: u(x,y) = exp(x)sin(πy), tests non-homogeneous BCs
        template<typename T>
        struct ExpSin2D {
            static T solution(T x, T y) {
                return std::exp(x) * std::sin(M_PI * y);
            }
            static T source(T x, T y) {
                // -∇²u = -(∂²u/∂x² + ∂²u/∂y²)
                // ∂²u/∂x² = exp(x)sin(πy)
                // ∂²u/∂y² = -π²exp(x)sin(πy)
                // -∇²u = (π² - 1)exp(x)sin(πy)
                return (M_PI * M_PI - T(1)) * std::exp(x) * std::sin(M_PI * y);
            }
            static T boundaryValue(T x, T y) {
                return solution(x, y);
            }
        };
        
        // 3D: u(x,y,z) = sin(πx)sin(πy)sin(πz)
        template<typename T>
        struct SinPi3D {
            static T solution(T x, T y, T z) {
                return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
            }
            static T source(T x, T y, T z) {
                return T(3) * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
            }
            static T boundaryValue(T, T, T) { return T(0); }
        };
    }

    ///////////////////////////////////////////////////////////////////////
    ///                      Error Computation                          ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Compute L∞ (max) error between numerical and exact solution
     */
    template<typename T, typename F>
    T errorLInf(const GridFunction1D<T>& numerical, F exact) {
        T maxErr = T(0);
        const auto& grid = numerical.grid();
        for (int i = 0; i < grid.numNodes(); ++i) {
            T err = std::abs(numerical(i) - exact(grid.x(i)));
            maxErr = std::max(maxErr, err);
        }
        return maxErr;
    }
    
    template<typename T, typename F>
    T errorLInf(const GridFunction2D<T>& numerical, F exact) {
        T maxErr = T(0);
        const auto& grid = numerical.grid();
        for (int j = 0; j < grid.numNodesY(); ++j) {
            for (int i = 0; i < grid.numNodesX(); ++i) {
                T err = std::abs(numerical(i, j) - exact(grid.x(i), grid.y(j)));
                maxErr = std::max(maxErr, err);
            }
        }
        return maxErr;
    }
    
    /**
     * @brief Compute L² (RMS) error between numerical and exact solution
     */
    template<typename T, typename F>
    T errorL2(const GridFunction1D<T>& numerical, F exact) {
        T sumSq = T(0);
        const auto& grid = numerical.grid();
        for (int i = 0; i < grid.numNodes(); ++i) {
            T err = numerical(i) - exact(grid.x(i));
            sumSq += err * err;
        }
        return std::sqrt(sumSq * grid.dx());
    }
    
    template<typename T, typename F>
    T errorL2(const GridFunction2D<T>& numerical, F exact) {
        T sumSq = T(0);
        const auto& grid = numerical.grid();
        for (int j = 0; j < grid.numNodesY(); ++j) {
            for (int i = 0; i < grid.numNodesX(); ++i) {
                T err = numerical(i, j) - exact(grid.x(i), grid.y(j));
                sumSq += err * err;
            }
        }
        return std::sqrt(sumSq * grid.dx() * grid.dy());
    }

} // namespace MML::PDE

#endif // MML_PDE_POISSON_SOLVER_H
