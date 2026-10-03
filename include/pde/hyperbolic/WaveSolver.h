///////////////////////////////////////////////////////////////////////////////////////////
// WaveSolver.h - Finite Difference Solver for Wave Equations
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Solves hyperbolic PDEs of the form:
//   ∂²u/∂t² = c²∇²u + f(x,t)     in Ω × (0,T]     (Wave equation)
//
// with initial conditions:
//   u(x, 0) = u₀(x)        initial displacement
//   ∂u/∂t(x, 0) = v₀(x)    initial velocity
//
// and Dirichlet/Neumann/Periodic boundary conditions.
//
// Time Integration Schemes:
// - Leapfrog: Explicit, O(Δt²), conditionally stable, non-dissipative
// - Lax-Wendroff: Explicit, O(Δt²), some numerical dissipation
//
// Features:
// - 1D, 2D, and 3D wave equation solvers
// - CFL stability condition enforcement
// - Energy conservation monitoring
// - Observer pattern for solution snapshots
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_WAVE_SOLVER_H
#define MML_PDE_WAVE_SOLVER_H

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

#include "../grid/Grid.h"
#include "../grid/Grid1D.h"
#include "../grid/Grid2D.h"
#include "../grid/Grid3D.h"
#include "../grid/GridFunction.h"
#include "../grid/BoundaryConditions.h"

namespace MML::PDE
{
    ///////////////////////////////////////////////////////////////////////
    ///                     Wave Equation Schemes                       ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Time integration schemes for hyperbolic PDEs
     */
    enum class WaveScheme {
        Leapfrog,        ///< Centered in time and space, O(Δt²), non-dissipative
        LaxWendroff      ///< Second-order, some dissipation for smoothing
    };

    /**
     * @brief Result of a wave equation solve
     */
    template<typename T>
    struct WaveSolveResult {
        bool stable;              ///< Whether solution remained stable
        int steps;                ///< Number of time steps taken
        T finalTime;              ///< Final simulation time reached
        T maxDisplacement;        ///< Maximum displacement at final time
        T minDisplacement;        ///< Minimum displacement at final time
        T totalEnergy;            ///< Total energy (kinetic + potential)
        T initialEnergy;          ///< Initial total energy (for conservation check)
        T energyError;            ///< Relative energy error |E - E₀|/E₀
        std::string message;      ///< Status message
    };

    ///////////////////////////////////////////////////////////////////////
    ///                      1D Wave Solver                             ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 1D wave equation: ∂²u/∂t² = c²·∂²u/∂x²
     * 
     * The wave equation is a second-order hyperbolic PDE that models:
     * - Vibrating strings (guitar, violin)
     * - Acoustic waves in tubes
     * - Electromagnetic waves in 1D
     * 
     * Spatial discretization: Central differences (2nd order)
     *   ∂²u/∂x² ≈ (u_{i-1} - 2u_i + u_{i+1}) / h²
     * 
     * Time discretization: Leapfrog (2nd order)
     *   (u^{n+1} - 2u^n + u^{n-1})/Δt² = c²·(u_{i-1} - 2u_i + u_{i+1})/h²
     * 
     * Rearranged:
     *   u^{n+1}_i = 2u^n_i - u^{n-1}_i + r²·(u^n_{i-1} - 2u^n_i + u^n_{i+1})
     * 
     * where r = c·Δt/h is the Courant number.
     * 
     * CFL Stability Condition: r ≤ 1 (i.e., Δt ≤ h/c)
     * 
     * @tparam T Numeric type (float, double)
     */
    template<typename T>
    class WaveSolver1D {
    private:
        Grid1D<T> grid_;
        BoundaryConditions1D<T> bc_;
        T waveSpeed_;                          // Wave speed c
        T currentTime_;
        std::vector<T> u_curr_;                // Current displacement u^n
        std::vector<T> u_prev_;                // Previous displacement u^{n-1}
        std::vector<T> velocity_;              // Current velocity ∂u/∂t (for energy calculation)
        std::function<T(T, T)> source_;        // Source term f(x, t)
        std::function<void(T, const std::vector<T>&, const std::vector<T>&)> observer_;
        T initialEnergy_;                      // Store initial energy for conservation check
        bool energyComputed_;
        
    public:
        /**
         * @brief Construct wave solver
         * @param grid Spatial grid
         * @param bc Boundary conditions
         * @param waveSpeed Wave propagation speed (c > 0)
         */
                WaveSolver1D(const Grid1D<T>& grid, const BoundaryConditions1D<T>& bc, T waveSpeed);
        
        /**
         * @brief Set initial displacement u(x, 0) = u0(x)
         */
        void setInitialDisplacement(std::function<T(T)> u0);
        
        /**
         * @brief Set initial velocity ∂u/∂t(x, 0) = v0(x)
         * 
         * Must be called after setInitialDisplacement.
         * Uses the initial velocity to compute u^{-1} for leapfrog startup.
         */
        void setInitialVelocity(std::function<T(T)> v0);
        
        /**
         * @brief Initialize the solver with first time step
         * 
         * Must be called after setting initial conditions and before stepping.
         * Computes u^{-1} from initial displacement and velocity using Taylor expansion:
         *   u^{-1} ≈ u^0 - Δt·v^0 + (Δt²/2)·c²·∇²u^0
         * 
         * @param dt Time step size (used for first step setup)
         */
        void initialize(T dt);
        
        /**
         * @brief Set source term f(x, t) in ∂²u/∂t² = c²∇²u + f
         */
        void setSource(std::function<T(T, T)> f);
        
        /**
         * @brief Set observer for solution snapshots
         * 
         * Observer is called after each time step with (time, displacement, velocity).
         */
        void setObserver(std::function<void(T, const std::vector<T>&, const std::vector<T>&)> obs);
        
        /**
         * @brief Get current simulation time
         */
        T getTime() const;
        
        /**
         * @brief Get current displacement
         */
        const std::vector<T>& getDisplacement() const;
        
        /**
         * @brief Get current velocity
         */
        const std::vector<T>& getVelocity() const;
        
        /**
         * @brief Get grid
         */
        const Grid1D<T>& getGrid() const;
        
        /**
         * @brief Compute CFL stability limit
         * @return Maximum stable time step: Δt_max = h/c
         */
        T computeCFLLimit() const;
        
        /**
         * @brief Compute Courant number r = c·Δt/h
         */
        T computeCourantNumber(T dt) const;
        
        /**
         * @brief Check if given time step satisfies CFL condition
         */
        bool isStable(T dt) const;
        
        /**
         * @brief Compute total energy E = (1/2)∫[(∂u/∂t)² + c²(∂u/∂x)²]dx
         * 
         * For the wave equation, energy should be conserved (for f=0 and fixed BCs).
         */
        T computeEnergy(T dt) const;
        
        /**
         * @brief Perform one time step using Leapfrog scheme
         * @param dt Time step size
         * @return true if step was stable
         */
        bool step(T dt, WaveScheme scheme = WaveScheme::Leapfrog);
        
        /**
         * @brief Solve until final time
         * @param finalTime Target simulation time
         * @param dt Time step size
         * @param scheme Time integration scheme
         * @return Solve result with statistics
         */
        WaveSolveResult<T> solve(T finalTime, T dt, WaveScheme scheme = WaveScheme::Leapfrog);
        
    private:
        /**
         * @brief Apply boundary conditions to solution vector
         */
        void applyBoundaryConditions(std::vector<T>& u);
        
        /**
         * @brief Leapfrog (centered) explicit step
         * 
         * u^{n+1}_i = 2u^n_i - u^{n-1}_i + r²·(u^n_{i-1} - 2u^n_i + u^n_{i+1}) + Δt²·f^n_i
         */
        bool stepLeapfrog(T dt);
        
        /**
         * @brief Lax-Wendroff step (more dissipative than Leapfrog)
         * 
         * Two-step method that is second-order accurate in both time and space.
         * Has some numerical dissipation which can help smooth oscillations.
         */
        bool stepLaxWendroff(T dt);
    };

    extern template class MML_PDE_TEMPLATE_API WaveSolver1D<float>;
    extern template class MML_PDE_TEMPLATE_API WaveSolver1D<double>;
    extern template class MML_PDE_TEMPLATE_API WaveSolver1D<long double>;

    ///////////////////////////////////////////////////////////////////////
    ///                      2D Wave Solver                             ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Solves the 2D wave equation: ∂²u/∂t² = c²(∂²u/∂x² + ∂²u/∂y²)
     * 
     * Models vibrating membranes (drums), water surface waves, etc.
     * 
     * CFL Stability Condition: c·Δt/h ≤ 1/√2 ≈ 0.707
     * 
     * @tparam T Numeric type (float, double)
     */
    template<typename T>
    class WaveSolver2D {
    private:
        Grid2D<T> grid_;
        BoundaryConditions2D<T> bc_;
        T waveSpeed_;
        T currentTime_;
        std::vector<T> u_curr_;
        std::vector<T> u_prev_;
        std::vector<T> velocity_;
        std::function<T(T, T, T)> source_;
        std::function<void(T, const std::vector<T>&)> observer_;
        T initialEnergy_;
        bool energyComputed_;
        
    public:
        /**
         * @brief Construct 2D wave solver
         * @param grid 2D spatial grid
         * @param bc Boundary conditions
         * @param waveSpeed Wave propagation speed (c > 0)
         */
                WaveSolver2D(const Grid2D<T>& grid, const BoundaryConditions2D<T>& bc, T waveSpeed);
        
        /**
         * @brief Set initial displacement u(x, y, 0) = u0(x, y)
         */
        void setInitialDisplacement(std::function<T(T, T)> u0);
        
        /**
         * @brief Set initial velocity ∂u/∂t(x, y, 0) = v0(x, y)
         */
        void setInitialVelocity(std::function<T(T, T)> v0);
        
        /**
         * @brief Initialize solver for time stepping
         */
        void initialize(T dt);
        
        /**
         * @brief Set source term f(x, y, t)
         */
        void setSource(std::function<T(T, T, T)> f);
        
        /**
         * @brief Set observer for solution snapshots
         */
        void setObserver(std::function<void(T, const std::vector<T>&)> obs);
        
        T getTime() const;
        const std::vector<T>& getDisplacement() const;
        const std::vector<T>& getVelocity() const;
        const Grid2D<T>& getGrid() const;
        
        /**
         * @brief Compute CFL stability limit for 2D: Δt_max = h/(c·√2)
         */
        T computeCFLLimit() const;
        
        /**
         * @brief Compute Courant number
         */
        T computeCourantNumber(T dt) const;

        bool isStable(T dt) const;

        /**
         * @brief Compute total energy
         */
        T computeEnergy(T dt) const;
        
        /**
         * @brief Perform one time step
         */
        bool step(T dt, WaveScheme scheme = WaveScheme::Leapfrog);
        
        /**
         * @brief Solve until final time
         */
        WaveSolveResult<T> solve(T finalTime, T dt, WaveScheme scheme = WaveScheme::Leapfrog);
        
    private:
        void applyBoundaryConditions(std::vector<T>& u);

        bool stepLeapfrog(T dt);
    };

    extern template class MML_PDE_TEMPLATE_API WaveSolver2D<float>;
    extern template class MML_PDE_TEMPLATE_API WaveSolver2D<double>;
    extern template class MML_PDE_TEMPLATE_API WaveSolver2D<long double>;

} // namespace MML::PDE

#endif // MML_PDE_WAVE_SOLVER_H
