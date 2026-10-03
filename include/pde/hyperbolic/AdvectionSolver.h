///////////////////////////////////////////////////////////////////////////////////////////
// AdvectionSolver.h - Finite Difference Solver for Advection Equations
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Solves first-order hyperbolic PDEs of the form:
//   ∂u/∂t + c·∂u/∂x = 0        in Ω × (0,T]     (Linear advection)
//   ∂u/∂t + a·∇u = 0           (Multi-dimensional advection)
//
// with initial condition u(x, 0) = u₀(x) and appropriate boundary conditions.
//
// The advection equation models pure transport/convection without diffusion.
// The solution translates: u(x, t) = u₀(x - ct)
//
// Numerical Schemes:
// - Upwind (First-order): Stable, but diffusive
// - Lax-Friedrichs: Simple, stable, diffusive
// - Lax-Wendroff: Second-order, some dispersion
// - Beam-Warming: Second-order upwind
// - FTCS: Centered differences - UNCONDITIONALLY UNSTABLE (for reference only)
//
// CFL Stability Condition: |c|·Δt/h ≤ 1
//
// Features:
// - 1D and 2D advection solvers
// - Multiple numerical schemes
// - Periodic and inflow/outflow boundary conditions
// - Mass conservation monitoring
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_ADVECTION_SOLVER_H
#define MML_PDE_ADVECTION_SOLVER_H

#include <mml/mml_export.h>
#include <functional>
#include <cmath>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <numeric>

#include "../grid/Grid.h"
#include "../grid/Grid1D.h"
#include "../grid/Grid2D.h"
#include "../grid/BoundaryConditions.h"

namespace MML::PDE
{
    ///////////////////////////////////////////////////////////////////////
    ///                   Advection Equation Schemes                    ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Numerical schemes for advection equations
     */
    enum class AdvectionScheme {
        Upwind,          ///< First-order upwind - stable, O(h) diffusive
        LaxFriedrichs,   ///< First-order centered with averaging - simple, diffusive
        LaxWendroff,     ///< Second-order - less diffusion, some dispersion (wiggles)
        BeamWarming,     ///< Second-order upwind - alternative to Lax-Wendroff
        FTCS             ///< Forward Time, Centered Space - UNSTABLE! (for demonstration)
    };

    /**
     * @brief Result of an advection equation solve
     */
    template<typename T>
    struct AdvectionSolveResult {
        bool stable;              ///< Whether solution remained stable
        int steps;                ///< Number of time steps taken
        T finalTime;              ///< Final simulation time reached
        T maxValue;               ///< Maximum value at final time
        T minValue;               ///< Minimum value at final time
        T totalMass;              ///< Integral of u over domain (conservation check)
        T initialMass;            ///< Initial total mass
        T massError;              ///< Relative mass error |M - M₀|/M₀
        T maxAbsValue;            ///< Maximum absolute value (for oscillation detection)
        std::string message;      ///< Status message
    };

    ///////////////////////////////////////////////////////////////////////
    ///                    1D Advection Solver                          ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 1D linear advection equation: ∂u/∂t + c·∂u/∂x = 0
     * 
     * The advection equation describes transport of a quantity u at constant
     * velocity c. The exact solution is u(x, t) = u₀(x - ct), i.e., the 
     * initial profile translates without changing shape.
     * 
     * Key numerical challenges:
     * - Standard centered differences are unconditionally unstable
     * - Upwinding is needed for stability
     * - Numerical diffusion smears sharp features
     * - Numerical dispersion causes spurious oscillations
     * 
     * @tparam T Numeric type (float, double)
     */
    template<typename T>
    class AdvectionSolver1D {
    private:
        Grid1D<T> grid_;
        T velocity_;                           // Advection velocity c
        T currentTime_;
        std::vector<T> solution_;
        std::function<void(T, const std::vector<T>&)> observer_;
        T initialMass_;
        bool periodic_;                        // Use periodic boundary conditions
        
    public:
        /**
         * @brief Construct advection solver
         * @param grid Spatial grid
         * @param velocity Advection velocity c (can be positive or negative)
         * @param periodic Use periodic boundary conditions (default: true)
         */
            AdvectionSolver1D(const Grid1D<T>& grid, T velocity, bool periodic = true);
        
        /**
         * @brief Set initial condition u(x, 0) = u0(x)
         */
        void setInitialCondition(std::function<T(T)> u0);
        
        /**
         * @brief Set observer for solution snapshots
         */
        void setObserver(std::function<void(T, const std::vector<T>&)> obs);
        
        /**
         * @brief Enable/disable periodic boundary conditions
         */
        void setPeriodic(bool periodic);
        
        T getTime() const;
        const std::vector<T>& getSolution() const;
        const Grid1D<T>& getGrid() const;
        T getVelocity() const;
        
        /**
         * @brief Compute CFL stability limit
         * @return Maximum stable time step: Δt_max = h/|c|
         */
        T computeCFLLimit() const;
        
        /**
         * @brief Compute Courant number |c|·Δt/h
         */
        T computeCourantNumber(T dt) const;
        
        /**
         * @brief Check if time step satisfies CFL condition
         */
        bool isStable(T dt) const;
        
        /**
         * @brief Compute total mass (integral of u)
         */
        T computeMass() const;
        
        /**
         * @brief Perform one time step
         * @param dt Time step size
         * @param scheme Numerical scheme to use
         * @return true if step was stable
         */
        bool step(T dt, AdvectionScheme scheme = AdvectionScheme::Upwind);
        
        /**
         * @brief Solve until final time
         * @param finalTime Target simulation time
         * @param dt Time step size
         * @param scheme Numerical scheme
         * @return Solve result with statistics
         */
        AdvectionSolveResult<T> solve(T finalTime, T dt, AdvectionScheme scheme = AdvectionScheme::Upwind);
        
    private:
        /**
         * @brief Get index with periodic wrapping
         */
        int periodicIndex(int i) const;
        
        /**
         * @brief First-order upwind scheme
         * 
         * For c > 0 (rightward flow): u_i^{n+1} = u_i^n - r·(u_i^n - u_{i-1}^n)
         * For c < 0 (leftward flow):  u_i^{n+1} = u_i^n - r·(u_{i+1}^n - u_i^n)
         * 
         * Uses information from upwind direction (where flow comes from).
         * Stable for |r| ≤ 1. First-order accurate, introduces numerical diffusion.
         */
        bool stepUpwind(T dt);
        
        /**
         * @brief Lax-Friedrichs scheme
         * 
         * u_i^{n+1} = (u_{i+1}^n + u_{i-1}^n)/2 - (r/2)·(u_{i+1}^n - u_{i-1}^n)
         * 
         * Averages neighbors before applying centered difference.
         * Stable for |r| ≤ 1. Very diffusive.
         */
        bool stepLaxFriedrichs(T dt);
        
        /**
         * @brief Lax-Wendroff scheme
         * 
         * u_i^{n+1} = u_i^n - (r/2)·(u_{i+1}^n - u_{i-1}^n) 
         *                   + (r²/2)·(u_{i+1}^n - 2u_i^n + u_{i-1}^n)
         * 
         * Second-order accurate in both space and time.
         * Stable for |r| ≤ 1. Less diffusive than upwind but can produce
         * spurious oscillations (dispersion) near discontinuities.
         */
        bool stepLaxWendroff(T dt);
        
        /**
         * @brief Beam-Warming scheme (second-order upwind)
         * 
         * For c > 0:
         * u_i^{n+1} = u_i^n - (r/2)·(3u_i^n - 4u_{i-1}^n + u_{i-2}^n)
         *                   + (r²/2)·(u_i^n - 2u_{i-1}^n + u_{i-2}^n)
         * 
         * Second-order accurate with upwind bias.
         * Less oscillatory than Lax-Wendroff for smooth solutions.
         */
        bool stepBeamWarming(T dt);
        
        /**
         * @brief Forward Time, Centered Space (FTCS) - UNSTABLE!
         * 
         * u_i^{n+1} = u_i^n - (r/2)·(u_{i+1}^n - u_{i-1}^n)
         * 
         * This scheme is UNCONDITIONALLY UNSTABLE for the advection equation!
         * Included only for educational purposes to demonstrate instability.
         */
        bool stepFTCS(T dt);
    };

    extern template class MML_PDE_TEMPLATE_API AdvectionSolver1D<float>;
    extern template class MML_PDE_TEMPLATE_API AdvectionSolver1D<double>;
    extern template class MML_PDE_TEMPLATE_API AdvectionSolver1D<long double>;

    ///////////////////////////////////////////////////////////////////////
    ///                    2D Advection Solver                          ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 2D linear advection equation: ∂u/∂t + a·∂u/∂x + b·∂u/∂y = 0
     * 
     * Models transport in 2D at constant velocity (a, b).
     * 
     * CFL Stability Condition: |a|·Δt/hx + |b|·Δt/hy ≤ 1
     * 
     * @tparam T Numeric type (float, double)
     */
    template<typename T>
    class AdvectionSolver2D {
    private:
        Grid2D<T> grid_;
        T velocityX_;                          // x-component of velocity
        T velocityY_;                          // y-component of velocity
        T currentTime_;
        std::vector<T> solution_;
        std::function<void(T, const std::vector<T>&)> observer_;
        T initialMass_;
        bool periodicX_;
        bool periodicY_;
        
    public:
        /**
         * @brief Construct 2D advection solver
         * @param grid 2D spatial grid
         * @param vx x-component of advection velocity
         * @param vy y-component of advection velocity
         * @param periodicX Use periodic BCs in x-direction
         * @param periodicY Use periodic BCs in y-direction
         */
        AdvectionSolver2D(const Grid2D<T>& grid, T vx, T vy, 
                                                    bool periodicX = true, bool periodicY = true);
        
        /**
         * @brief Set initial condition u(x, y, 0) = u0(x, y)
         */
        void setInitialCondition(std::function<T(T, T)> u0);
        
        void setObserver(std::function<void(T, const std::vector<T>&)> obs);
        
        T getTime() const;
        const std::vector<T>& getSolution() const;
        const Grid2D<T>& getGrid() const;
        
        /**
         * @brief Compute CFL stability limit
         */
        T computeCFLLimit() const;
        
        /**
         * @brief Compute Courant number
         */
        T computeCourantNumber(T dt) const;
        
        bool isStable(T dt) const;
        
        /**
         * @brief Compute total mass
         */
        T computeMass() const;
        
        /**
         * @brief Perform one time step using dimensional splitting (upwind)
         */
        bool step(T dt, AdvectionScheme scheme = AdvectionScheme::Upwind);
        
        /**
         * @brief Solve until final time
         */
        AdvectionSolveResult<T> solve(T finalTime, T dt, AdvectionScheme scheme = AdvectionScheme::Upwind);
        
    private:
        int periodicIndexX(int i) const;
        
        int periodicIndexY(int j) const;
        
        /**
         * @brief Dimensional splitting with upwind
         */
        bool stepDimensionalSplit(T dt, AdvectionScheme scheme);
    };

    extern template class MML_PDE_TEMPLATE_API AdvectionSolver2D<float>;
    extern template class MML_PDE_TEMPLATE_API AdvectionSolver2D<double>;
    extern template class MML_PDE_TEMPLATE_API AdvectionSolver2D<long double>;

} // namespace MML::PDE

#endif // MML_PDE_ADVECTION_SOLVER_H
