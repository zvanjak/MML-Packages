///////////////////////////////////////////////////////////////////////////////////////////
// HeatSolver.h - Finite Difference Solver for Heat and Diffusion Equations
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Solves parabolic PDEs of the form:
//   ∂u/∂t = α∇²u + f(x,t)     in Ω × (0,T]     (Heat equation)
//
// with initial conditions and Dirichlet/Neumann boundary conditions.
//
// Time Integration Schemes:
// - Forward Euler (FTCS): Explicit, O(Δt), conditionally stable
// - Backward Euler (BTCS): Implicit, O(Δt), unconditionally stable
// - Crank-Nicolson: Implicit, O(Δt²), unconditionally stable (recommended)
// - θ-method: Generalized scheme, θ∈[0,1]
//
// Features:
// - 1D, 2D, and 3D heat equation solvers
// - Multiple time stepping schemes
// - Automatic CFL stability check for explicit methods
// - Observer pattern for solution snapshots
// - Energy/conservation monitoring
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_HEAT_SOLVER_H
#define MML_PDE_HEAT_SOLVER_H

#include <mml/mml_export.h>
#include <functional>
#include <cmath>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
    ///                     Time Stepping Schemes                       ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Time integration schemes for parabolic PDEs
     */
    enum class TimeScheme {
        ForwardEuler,    ///< Explicit, O(Δt), conditionally stable: Δt ≤ h²/(2α·dim)
        BackwardEuler,   ///< Implicit, O(Δt), unconditionally stable
        CrankNicolson,   ///< Implicit, O(Δt²), unconditionally stable (recommended)
        Theta            ///< θ-method: θ=0 (FE), θ=0.5 (CN), θ=1 (BE)
    };

    /**
     * @brief Result of a heat equation solve
     */
    template<typename T>
    struct HeatSolveResult {
        bool stable;              ///< Whether solution remained stable
        int steps;                ///< Number of time steps taken
        T finalTime;              ///< Final simulation time reached
        T maxTemperature;         ///< Maximum temperature at final time
        T minTemperature;         ///< Minimum temperature at final time
        T totalEnergy;            ///< Integral of u over domain (conservation check)
        std::string message;      ///< Status message
    };

    ///////////////////////////////////////////////////////////////////////
    ///                      1D Heat Solver                             ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 1D heat equation: ∂u/∂t = α·∂²u/∂x²
     * 
     * Spatial discretization: Central differences (2nd order)
     *   ∂²u/∂x² ≈ (u_{i-1} - 2u_i + u_{i+1}) / h²
     * 
     * Time discretization: θ-method
     *   (u^{n+1} - u^n)/Δt = α·[(1-θ)·L·u^n + θ·L·u^{n+1}]
     * 
     * where L is the discrete Laplacian operator.
     * 
     * @tparam T Numeric type (float, double)
     */
    template<typename T>
    class HeatSolver1D {
    private:
        Grid1D<T> grid_;
        BoundaryConditions1D<T> bc_;
        T alpha_;                              // Thermal diffusivity
        T currentTime_;
        std::vector<T> solution_;              // Current solution values
        std::function<T(T, T)> source_;        // Source term f(x, t)
        std::function<void(T, const std::vector<T>&)> observer_;
        T theta_;                              // θ parameter for θ-method
        
        // Precomputed matrices for implicit schemes
        SparseMatrixCSR<T> implicitMatrix_;
        bool matricesAssembled_;
        T assembledDt_;
        
    public:
        /**
         * @brief Construct heat solver
         * @param grid Spatial grid
         * @param bc Boundary conditions
         * @param alpha Thermal diffusivity (α > 0)
         */
                HeatSolver1D(const Grid1D<T>& grid, const BoundaryConditions1D<T>& bc, T alpha);
        
        /**
         * @brief Set initial condition u(x, 0) = u0(x)
         */
        void setInitialCondition(std::function<T(T)> u0);
        
        /**
         * @brief Set source term f(x, t) in ∂u/∂t = α∇²u + f
         */
        void setSource(std::function<T(T, T)> f);
        
        /**
         * @brief Set observer for solution snapshots
         * 
         * Observer is called after each time step with (time, solution).
         * Useful for animation, data collection, or convergence monitoring.
         */
        void setObserver(std::function<void(T, const std::vector<T>&)> obs);
        
        /**
         * @brief Set θ parameter for θ-method
         * @param theta 0 = Forward Euler, 0.5 = Crank-Nicolson, 1 = Backward Euler
         */
        void setTheta(T theta);
        
        /**
         * @brief Get current simulation time
         */
        T getTime() const;
        
        /**
         * @brief Get current solution
         */
        const std::vector<T>& getSolution() const;
        
        /**
         * @brief Get grid
         */
        const Grid1D<T>& getGrid() const;
        
        /**
         * @brief Compute CFL stability limit for explicit method
         * @return Maximum stable time step for Forward Euler
         */
        T computeCFLLimit() const;
        
        /**
         * @brief Check if given time step is stable for explicit method
         */
        bool isStable(T dt) const;
        
        /**
         * @brief Perform one time step using specified scheme
         * @param dt Time step size
         * @param scheme Time integration scheme
         * @return true if step was stable
         */
        bool step(T dt, TimeScheme scheme = TimeScheme::CrankNicolson);
        
        /**
         * @brief Solve until final time
         * @param finalTime Target simulation time
         * @param dt Time step size
         * @param scheme Time integration scheme
         * @return Solve result with statistics
         */
        HeatSolveResult<T> solve(T finalTime, T dt, TimeScheme scheme = TimeScheme::CrankNicolson);
        
    private:
        /**
         * @brief Forward Euler (FTCS) explicit step
         * 
         * u^{n+1} = u^n + Δt·α·L·u^n
         * 
         * Stability: Δt ≤ h²/(2α)
         */
        bool stepForwardEuler(T dt);
        
        /**
         * @brief θ-method implicit step
         * 
         * (I - θ·Δt·α·L)·u^{n+1} = (I + (1-θ)·Δt·α·L)·u^n + Δt·f^{n+θ}
         * 
         * θ=0: Forward Euler, θ=0.5: Crank-Nicolson, θ=1: Backward Euler
         */
        bool stepTheta(T dt);
        
        /**
         * @brief Assemble the implicit matrix (I - θ·Δt·α·L)
         * 
         * For Neumann BCs, uses ghost point method:
         * - du/dx = g => ghost point u_{-1} = u_1 - 2h·g
         * - This modifies the stencil at the boundary
         */
        void assembleImplicitMatrix(T dt);
        
        /**
         * @brief Apply boundary conditions to solution vector
         */
        void applyBoundaryConditions(std::vector<T>& u, T t);
        
        /**
         * @brief Apply boundary conditions to RHS for implicit solve
         * 
         * For Neumann BCs with ghost point method, we need to add the flux contribution
         */
        void applyBoundaryConditionsToRHS(std::vector<T>& rhs, T t);
    };

    extern template class MML_PDE_TEMPLATE_API HeatSolver1D<float>;
    extern template class MML_PDE_TEMPLATE_API HeatSolver1D<double>;
    extern template class MML_PDE_TEMPLATE_API HeatSolver1D<long double>;

    ///////////////////////////////////////////////////////////////////////
    ///                      2D Heat Solver                             ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 2D heat equation: ∂u/∂t = α(∂²u/∂x² + ∂²u/∂y²)
     * 
     * Uses the standard 5-point stencil for the Laplacian.
     * 
     * @tparam T Numeric type (float, double)
     */
    template<typename T>
    class HeatSolver2D {
    private:
        Grid2D<T> grid_;
        BoundaryConditions2D<T> bc_;
        T alpha_;
        T currentTime_;
        std::vector<T> solution_;
        std::function<T(T, T, T)> source_;  // f(x, y, t)
        std::function<void(T, const std::vector<T>&)> observer_;
        T theta_;
        
        SparseMatrixCSR<T> implicitMatrix_;
        bool matricesAssembled_;
        T assembledDt_;
        
    public:
        HeatSolver2D(const Grid2D<T>& grid, const BoundaryConditions2D<T>& bc, T alpha);
        
        void setInitialCondition(std::function<T(T, T)> u0);
        
        void setSource(std::function<T(T, T, T)> f);
        
        void setObserver(std::function<void(T, const std::vector<T>&)> obs);
        
        void setTheta(T theta);
        
        T getTime() const;
        const std::vector<T>& getSolution() const;
        const Grid2D<T>& getGrid() const;
        
        /**
         * @brief CFL limit for 2D: Δt ≤ h²/(4α) for equal spacing
         */
        T computeCFLLimit() const;
        
        bool isStable(T dt) const;
        
        bool step(T dt, TimeScheme scheme = TimeScheme::CrankNicolson);
        
        HeatSolveResult<T> solve(T finalTime, T dt, TimeScheme scheme = TimeScheme::CrankNicolson);
        
    private:
        bool stepForwardEuler(T dt);
        
        bool stepTheta(T dt);
        
        void assembleImplicitMatrix(T dt);
        
        /**
         * @brief Get the BC type for a boundary node, handling corners
         */
        BCType getBoundaryType(int i, int j, int nx, int ny) const;
        
        /**
         * @brief Assemble matrix row for Neumann boundary using ghost point method
         * 
         * For Neumann BC (du/dn = g), we use the ghost point method:
         * - At left (i=0): u_{-1} = u_1 - 2*hx*g, so d²u/dx² ≈ 2(u_1 - u_0)/hx²
         * - At right (i=nx-1): u_{nx} = u_{nx-2} + 2*hx*g, so d²u/dx² ≈ 2(u_{nx-2} - u_{nx-1})/hx²
         * - Similar for y-direction
         */
        void assembleNeumannBoundaryRow(SparseMatrixCOO<T>& coo, int i, int j, int idx,
                                        T rx, T ry, int nx, int ny) const;
        
        /**
         * @brief Compute RHS for Neumann boundary node using ghost point method
         * 
         * For the explicit part of the theta-method, we need to compute
         * the Laplacian at boundary nodes using the ghost point formula.
         */
        T computeNeumannBoundaryRHS(int i, int j, int idx, T explicitRx, T explicitRy,
                                    T dt, T tMid, int nx, int ny) const;
        
        /**
         * @brief Apply boundary conditions to solution vector
         * 
         * For Dirichlet: set u = g
         * For Neumann: use ghost point formula u_boundary = u_interior ± h*g
         */
        void applyBoundaryConditions(std::vector<T>& u, [[maybe_unused]] T t);
        
        /**
         * @brief Apply boundary conditions to RHS for implicit solve
         * 
         * For Dirichlet: set RHS to boundary value
         * For Neumann: include flux contribution from ghost point elimination
         */
        void applyBoundaryConditionsToRHS(std::vector<T>& rhs, [[maybe_unused]] T t);
    };

    extern template class MML_PDE_TEMPLATE_API HeatSolver2D<float>;
    extern template class MML_PDE_TEMPLATE_API HeatSolver2D<double>;
    extern template class MML_PDE_TEMPLATE_API HeatSolver2D<long double>;

    ///////////////////////////////////////////////////////////////////////
    ///                      3D Heat Solver                             ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 3D heat equation: ∂u/∂t = α(∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²)
     * 
     * Uses the standard 7-point stencil for the Laplacian.
     * 
     * @tparam T Numeric type (float, double)
     */
    template<typename T>
    class HeatSolver3D {
    private:
        Grid3D<T> grid_;
        BoundaryConditions3D<T> bc_;
        T alpha_;
        T currentTime_;
        std::vector<T> solution_;
        std::function<T(T, T, T, T)> source_;  // f(x, y, z, t)
        std::function<void(T, const std::vector<T>&)> observer_;
        T theta_;
        
        SparseMatrixCSR<T> implicitMatrix_;
        bool matricesAssembled_;
        T assembledDt_;
        
    public:
        HeatSolver3D(const Grid3D<T>& grid, const BoundaryConditions3D<T>& bc, T alpha);
        
        void setInitialCondition(std::function<T(T, T, T)> u0);
        
        void setSource(std::function<T(T, T, T, T)> f);
        
        void setObserver(std::function<void(T, const std::vector<T>&)> obs);
        
        void setTheta(T theta);
        
        T getTime() const;
        const std::vector<T>& getSolution() const;
        const Grid3D<T>& getGrid() const;
        
        /**
         * @brief CFL limit for 3D: Δt ≤ 1/(2α·(1/hx² + 1/hy² + 1/hz²))
         */
        T computeCFLLimit() const;
        
        bool isStable(T dt) const;
        
        bool step(T dt, TimeScheme scheme = TimeScheme::CrankNicolson);
        
        HeatSolveResult<T> solve(T finalTime, T dt, TimeScheme scheme = TimeScheme::CrankNicolson);
        
    private:
        bool stepForwardEuler(T dt);
        
        bool stepTheta(T dt);
        
        /**
         * @brief Get the BC type for a 3D boundary node
         * 
         * Priority for corners/edges: Left/Right > Bottom/Top > Front/Back
         */
        BCType getBoundaryType(int i, int j, int k, int nx, int ny, int nz) const;
        
        /**
         * @brief Get the boundary side for a 3D boundary node
         * For 3D: Left/Right (X), Bottom/Top (Y), Front/Back (Z)
         */
        BoundarySide getBoundarySide(int i, int j, int k, int nx, int ny, int nz) const;
        
        /**
         * @brief Compute RHS for Neumann boundary node using ghost point method (3D)
         */
        T computeNeumannBoundaryRHS3D(int i, int j, int k, int idx,
                                      T explicitRx, T explicitRy, T explicitRz,
                                      T dt, T tMid, int nx, int ny, int nz) const;
        
        void assembleImplicitMatrix(T dt);
        
        /**
         * @brief Assemble matrix row for Neumann boundary using ghost point method (3D)
         */
        void assembleNeumannBoundaryRow3D(SparseMatrixCOO<T>& coo, int i, int j, int k, int idx,
                                          T rx, T ry, T rz, int nx, int ny, int nz) const;
        
        void applyBoundaryConditions(std::vector<T>& u, [[maybe_unused]] T t);
        
        void applyBoundaryConditionsToRHS(std::vector<T>& rhs, [[maybe_unused]] T t);
    };

    extern template class MML_PDE_TEMPLATE_API HeatSolver3D<float>;
    extern template class MML_PDE_TEMPLATE_API HeatSolver3D<double>;
    extern template class MML_PDE_TEMPLATE_API HeatSolver3D<long double>;

    ///////////////////////////////////////////////////////////////////////
    ///                   Manufactured Solutions                        ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Analytical solutions for testing heat equation solvers
     */
    namespace HeatManufacturedSolutions {
        
        /**
         * @brief 1D Exponentially decaying sinusoid
         * 
         * u(x, t) = exp(-α·π²·t) · sin(πx)
         * 
         * Satisfies ∂u/∂t = α·∂²u/∂x² on [0,1] with u(0,t) = u(1,t) = 0
         */
        template<typename T>
        struct ExpDecaySin1D {
            T alpha;
            
            ExpDecaySin1D(T a = T(1)) : alpha(a) {}
            
            T solution(T x, T t) const {
                return std::exp(-alpha * M_PI * M_PI * t) * std::sin(M_PI * x);
            }
            
            T initial(T x) const {
                return std::sin(M_PI * x);
            }
            
            T source(T x, T t) const {
                return T(0);  // Homogeneous heat equation
            }
            
            T boundary(T x, T t) const {
                return T(0);  // Dirichlet zero at x=0 and x=1
            }
        };
        
        /**
         * @brief 2D Exponentially decaying sinusoid
         * 
         * u(x, y, t) = exp(-α·2π²·t) · sin(πx) · sin(πy)
         * 
         * Satisfies ∂u/∂t = α·∇²u on [0,1]² with zero Dirichlet BCs
         */
        template<typename T>
        struct ExpDecaySin2D {
            T alpha;
            
            ExpDecaySin2D(T a = T(1)) : alpha(a) {}
            
            T solution(T x, T y, T t) const {
                return std::exp(-alpha * T(2) * M_PI * M_PI * t) * std::sin(M_PI * x) * std::sin(M_PI * y);
            }
            
            T initial(T x, T y) const {
                return std::sin(M_PI * x) * std::sin(M_PI * y);
            }
            
            T source(T x, T y, T t) const {
                return T(0);
            }
            
            T boundary(T x, T y, T t) const {
                return T(0);
            }
        };
        
        /**
         * @brief 1D Heat with source term
         * 
         * Choose u(x,t) = sin(πx)·(1 + t), then
         * ∂u/∂t = sin(πx)
         * ∂²u/∂x² = -π²·sin(πx)·(1+t)
         * 
         * Source: f = sin(πx) + α·π²·sin(πx)·(1+t)
         */
        template<typename T>
        struct WithSource1D {
            T alpha;
            
            WithSource1D(T a = T(1)) : alpha(a) {}
            
            T solution(T x, T t) const {
                return std::sin(M_PI * x) * (T(1) + t);
            }
            
            T initial(T x) const {
                return std::sin(M_PI * x);
            }
            
            T source(T x, T t) const {
                return std::sin(M_PI * x) + alpha * M_PI * M_PI * std::sin(M_PI * x) * (T(1) + t);
            }
            
            T boundary(T x, T t) const {
                return T(0);
            }
        };
        
        /**
         * @brief Room cooling model (simplified)
         * 
         * Initial: Hot room at temperature T_hot
         * Boundary: One side at T_cold (window), others at T_hot (walls)
         */
        template<typename T>
        struct RoomCooling2D {
            T T_hot;
            T T_cold;
            
            RoomCooling2D(T hot = T(20), T cold = T(0)) : T_hot(hot), T_cold(cold) {}
            
            T initial(T x, T y) const {
                return T_hot;  // Room starts at hot temperature
            }
            
            T boundary(T x, T y, T t) const {
                // Cold at x=0 (window), hot everywhere else
                if (std::abs(x) < T(1e-10)) {
                    return T_cold;  // Window - cold air
                }
                return T_hot;  // Other walls - insulated/warm
            }
        };
    }

    ///////////////////////////////////////////////////////////////////////
    ///                      Error Computation                          ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Compute L∞ error for 1D heat solver against exact solution
     */
    template<typename T, typename F>
    T heatErrorLInf1D(const HeatSolver1D<T>& solver, F&& exact) {
        const auto& sol = solver.getSolution();
        const auto& grid = solver.getGrid();
        T t = solver.getTime();
        
        T maxErr = T(0);
        for (int i = 0; i < grid.numNodes(); ++i) {
            T err = std::abs(sol[i] - exact(grid.x(i), t));
            maxErr = std::max(maxErr, err);
        }
        return maxErr;
    }
    
    /**
     * @brief Compute L² error for 1D heat solver
     */
    template<typename T, typename F>
    T heatErrorL2_1D(const HeatSolver1D<T>& solver, F&& exact) {
        const auto& sol = solver.getSolution();
        const auto& grid = solver.getGrid();
        T t = solver.getTime();
        T h = grid.dx();
        
        T sumSq = T(0);
        for (int i = 0; i < grid.numNodes(); ++i) {
            T err = sol[i] - exact(grid.x(i), t);
            T weight = (i == 0 || i == grid.numNodes() - 1) ? h / T(2) : h;
            sumSq += weight * err * err;
        }
        return std::sqrt(sumSq);
    }
    
    /**
     * @brief Compute L∞ error for 2D heat solver
     */
    template<typename T, typename F>
    T heatErrorLInf2D(const HeatSolver2D<T>& solver, F&& exact) {
        const auto& sol = solver.getSolution();
        const auto& grid = solver.getGrid();
        T t = solver.getTime();
        
        T maxErr = T(0);
        grid.forEachNode([&](int i, int j, T x, T y) {
            T err = std::abs(sol[grid.index(i, j)] - exact(x, y, t));
            maxErr = std::max(maxErr, err);
        });
        return maxErr;
    }

} // namespace MML::PDE

#endif // MML_PDE_HEAT_SOLVER_H
