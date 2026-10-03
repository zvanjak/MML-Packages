///////////////////////////////////////////////////////////////////////////////////////////
// ManufacturedSolutions.h - Analytical Solutions for PDE Solver Validation
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Manufactured solutions are exact analytical solutions used to validate numerical
// solvers. The process:
//   1. Choose an analytical solution u(x, t)
//   2. Compute the required source term f and boundary conditions
//   3. Solve numerically and compare with the known solution
//
// This file consolidates all manufactured solutions for different PDE types.
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_TESTBEDS_MANUFACTURED_SOLUTIONS_H
#define MML_PDE_TESTBEDS_MANUFACTURED_SOLUTIONS_H

#include <cmath>
#include <functional>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace MML::PDE::Testbeds {

    ///////////////////////////////////////////////////////////////////////
    ///          Elliptic PDE Manufactured Solutions (Poisson)          ///
    ///////////////////////////////////////////////////////////////////////
    
    namespace Elliptic {
        
        /**
         * @brief 1D: u(x) = sin(πx), -u'' = π²sin(πx)
         * Domain: [0, 1], Dirichlet BCs: u(0) = u(1) = 0
         */
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
            
            // Error norms
            template<typename GridFunc>
            static T errorLInf(const GridFunc& u, int n) {
                T maxErr = T(0);
                for (int i = 0; i <= n; ++i) {
                    T x = T(i) / T(n);
                    T err = std::abs(u(i) - solution(x));
                    maxErr = std::max(maxErr, err);
                }
                return maxErr;
            }
        };
        
        /**
         * @brief 1D: u(x) = x(1-x), -u'' = 2
         * Polynomial solution - exactly recovered by 2nd order scheme
         */
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
        
        /**
         * @brief 1D: u(x) = x³ - x, tests cubic source
         * -u'' = -6x, u(0) = 0, u(1) = 0
         */
        template<typename T>
        struct Cubic1D {
            static T solution(T x) {
                return x * x * x - x;
            }
            static T source(T x) {
                return -T(6) * x;
            }
            static T leftBC() { return T(0); }
            static T rightBC() { return T(0); }
        };
        
        /**
         * @brief 2D: u(x,y) = sin(πx)sin(πy), -∇²u = 2π²sin(πx)sin(πy)
         * Domain: [0,1]², Dirichlet BCs: u = 0 on boundary
         */
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
        
        /**
         * @brief 2D: u(x,y) = x(1-x)y(1-y)
         * Polynomial - good for verifying O(h²) accuracy
         */
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
        
        /**
         * @brief 2D: u(x,y) = exp(x)sin(πy), non-homogeneous BCs
         */
        template<typename T>
        struct ExpSin2D {
            static T solution(T x, T y) {
                return std::exp(x) * std::sin(M_PI * y);
            }
            static T source(T x, T y) {
                return (M_PI * M_PI - T(1)) * std::exp(x) * std::sin(M_PI * y);
            }
            static T boundaryValue(T x, T y) {
                return solution(x, y);
            }
        };
        
        /**
         * @brief 3D: u(x,y,z) = sin(πx)sin(πy)sin(πz)
         * -∇²u = 3π²sin(πx)sin(πy)sin(πz)
         */
        template<typename T>
        struct SinPi3D {
            static T solution(T x, T y, T z) {
                return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
            }
            static T source(T x, T y, T z) {
                return T(3) * M_PI * M_PI * solution(x, y, z);
            }
            static T boundaryValue(T, T, T) { return T(0); }
        };
        
    } // namespace Elliptic
    
    ///////////////////////////////////////////////////////////////////////
    ///         Parabolic PDE Manufactured Solutions (Heat)             ///
    ///////////////////////////////////////////////////////////////////////
    
    namespace Parabolic {
        
        /**
         * @brief 1D: u(x,t) = exp(-απ²t)·sin(πx)
         * 
         * Exponentially decaying fundamental mode.
         * ∂u/∂t = α·∂²u/∂x² (homogeneous heat equation)
         */
        template<typename T>
        struct ExpDecaySin1D {
            T alpha;
            
            ExpDecaySin1D(T a = T(1)) : alpha(a) {}
            
            T solution(T x, T t) const {
                return std::exp(-alpha * M_PI * M_PI * t) * std::sin(M_PI * x);
            }
            
            T initial(T x) const { return std::sin(M_PI * x); }
            T source(T, T) const { return T(0); }
            T boundary(T, T) const { return T(0); }
            
            // Decay constant
            T decayRate() const { return alpha * M_PI * M_PI; }
        };
        
        /**
         * @brief 1D with time-dependent source
         * 
         * u(x,t) = sin(πx)·(1 + t)
         * f(x,t) = sin(πx) + απ²·sin(πx)·(1+t)
         */
        template<typename T>
        struct WithSource1D {
            T alpha;
            
            WithSource1D(T a = T(1)) : alpha(a) {}
            
            T solution(T x, T t) const {
                return std::sin(M_PI * x) * (T(1) + t);
            }
            
            T initial(T x) const { return std::sin(M_PI * x); }
            
            T source(T x, T t) const {
                return std::sin(M_PI * x) + alpha * M_PI * M_PI * std::sin(M_PI * x) * (T(1) + t);
            }
            
            T boundary(T, T) const { return T(0); }
        };
        
        /**
         * @brief 2D: u(x,y,t) = exp(-2απ²t)·sin(πx)·sin(πy)
         */
        template<typename T>
        struct ExpDecaySin2D {
            T alpha;
            
            ExpDecaySin2D(T a = T(1)) : alpha(a) {}
            
            T solution(T x, T y, T t) const {
                return std::exp(-alpha * T(2) * M_PI * M_PI * t) * 
                       std::sin(M_PI * x) * std::sin(M_PI * y);
            }
            
            T initial(T x, T y) const { 
                return std::sin(M_PI * x) * std::sin(M_PI * y); 
            }
            T source(T, T, T) const { return T(0); }
            T boundary(T, T, T) const { return T(0); }
        };
        
        /**
         * @brief Gaussian blob diffusion (no closed form, for qualitative tests)
         */
        template<typename T>
        struct GaussianBlob2D {
            T x0, y0;      // Initial center
            T sigma0;      // Initial width
            T alpha;       // Diffusion coefficient
            
            GaussianBlob2D(T cx = T(0.5), T cy = T(0.5), T s = T(0.1), T a = T(1))
                : x0(cx), y0(cy), sigma0(s), alpha(a) {}
            
            // Exact solution for free diffusion (infinite domain approximation)
            T solution(T x, T y, T t) const {
                T sigma2 = sigma0 * sigma0 + T(4) * alpha * t;
                T dx = x - x0;
                T dy = y - y0;
                T r2 = dx * dx + dy * dy;
                return (sigma0 * sigma0 / sigma2) * std::exp(-r2 / sigma2);
            }
            
            T initial(T x, T y) const { return solution(x, y, T(0)); }
            
            // Width grows as sqrt(t)
            T currentWidth(T t) const {
                return std::sqrt(sigma0 * sigma0 + T(4) * alpha * t);
            }
        };
        
    } // namespace Parabolic
    
    ///////////////////////////////////////////////////////////////////////
    ///         Hyperbolic PDE Manufactured Solutions (Wave)            ///
    ///////////////////////////////////////////////////////////////////////
    
    namespace Hyperbolic {
        
        /**
         * @brief 1D Standing wave: u(x,t) = sin(πx)·cos(πct)
         * 
         * Fundamental mode of vibrating string with fixed ends.
         */
        template<typename T>
        struct StandingWave1D {
            T waveSpeed;
            
            StandingWave1D(T c = T(1)) : waveSpeed(c) {}
            
            T solution(T x, T t) const {
                return std::sin(M_PI * x) * std::cos(M_PI * waveSpeed * t);
            }
            
            T initialDisplacement(T x) const { return std::sin(M_PI * x); }
            T initialVelocity(T) const { return T(0); }
            T boundary(T) const { return T(0); }
            
            // Period of oscillation
            T period() const { return T(2) / waveSpeed; }
            
            // Natural frequency
            T frequency() const { return waveSpeed / T(2); }
        };
        
        /**
         * @brief 1D Traveling wave: u(x,t) = f(x - ct)
         * 
         * Right-traveling Gaussian pulse.
         */
        template<typename T>
        struct TravelingGaussian1D {
            T waveSpeed;
            T x0;       // Initial center
            T sigma;    // Width
            
            TravelingGaussian1D(T c = T(1), T center = T(0.3), T width = T(0.05))
                : waveSpeed(c), x0(center), sigma(width) {}
            
            T solution(T x, T t) const {
                T xi = x - waveSpeed * t - x0;
                return std::exp(-xi * xi / (sigma * sigma));
            }
            
            T initialDisplacement(T x) const { return solution(x, T(0)); }
            
            // For pure right-traveling wave: v₀ = c·∂u₀/∂x
            T initialVelocity(T x) const {
                T xi = x - x0;
                T dudx = -T(2) * xi / (sigma * sigma) * std::exp(-xi * xi / (sigma * sigma));
                return waveSpeed * dudx;
            }
            
            // Current center position
            T centerPosition(T t) const { return x0 + waveSpeed * t; }
        };
        
        /**
         * @brief 2D Standing wave on square membrane
         * 
         * u(x,y,t) = sin(mπx)·sin(nπy)·cos(ωt)
         * where ω = πc·√(m² + n²)
         */
        template<typename T>
        struct MembraneMode2D {
            T waveSpeed;
            int m, n;    // Mode numbers
            
            MembraneMode2D(T c = T(1), int mode_x = 1, int mode_y = 1)
                : waveSpeed(c), m(mode_x), n(mode_y) {}
            
            T frequency() const {
                return M_PI * waveSpeed * std::sqrt(T(m * m + n * n));
            }
            
            T solution(T x, T y, T t) const {
                return std::sin(m * M_PI * x) * std::sin(n * M_PI * y) * 
                       std::cos(frequency() * t);
            }
            
            T initialDisplacement(T x, T y) const {
                return std::sin(m * M_PI * x) * std::sin(n * M_PI * y);
            }
            
            T initialVelocity(T, T) const { return T(0); }
            T boundary(T, T) const { return T(0); }
        };
        
        /**
         * @brief 1D Advection: u(x,t) = u₀(x - ct)
         * 
         * Translation of initial profile at constant speed.
         */
        template<typename T>
        struct AdvectedGaussian1D {
            T velocity;
            T x0;
            T sigma;
            
            AdvectedGaussian1D(T v = T(1), T center = T(0.3), T width = T(0.1))
                : velocity(v), x0(center), sigma(width) {}
            
            T solution(T x, T t) const {
                T xi = x - velocity * t - x0;
                return std::exp(-xi * xi / (sigma * sigma));
            }
            
            T initial(T x) const { return solution(x, T(0)); }
            
            T centerPosition(T t) const { return x0 + velocity * t; }
            
            // For testing mass conservation
            T totalMass() const {
                return sigma * std::sqrt(M_PI);  // Integral of Gaussian
            }
        };
        
    } // namespace Hyperbolic
    
    ///////////////////////////////////////////////////////////////////////
    ///                   Convergence Test Utilities                    ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Compute observed convergence rate from error sequence
     * 
     * For second-order scheme: rate ≈ 2 (error halves when h halves → 4x reduction)
     * 
     * @param errors Vector of errors for decreasing h (h, h/2, h/4, ...)
     * @return Observed convergence rate (should be ~2 for O(h²) scheme)
     */
    template<typename T>
    T computeConvergenceRate(const std::vector<T>& errors) {
        if (errors.size() < 2) return T(0);
        
        // Average the log ratio
        T totalRate = T(0);
        int count = 0;
        for (size_t i = 1; i < errors.size(); ++i) {
            if (errors[i] > T(1e-15) && errors[i-1] > T(1e-15)) {
                totalRate += std::log(errors[i-1] / errors[i]) / std::log(T(2));
                count++;
            }
        }
        
        return count > 0 ? totalRate / T(count) : T(0);
    }
    
    /**
     * @brief Verify second-order convergence
     * 
     * @param errors Error sequence from grid refinement
     * @param tolerance How close to 2.0 the rate should be (default 0.2)
     * @return true if observed rate is within tolerance of 2.0
     */
    template<typename T>
    bool verifySecondOrderConvergence(const std::vector<T>& errors, T tolerance = T(0.2)) {
        T rate = computeConvergenceRate(errors);
        return std::abs(rate - T(2)) < tolerance;
    }
    
} // namespace MML::PDE::Testbeds

#endif // MML_PDE_TESTBEDS_MANUFACTURED_SOLUTIONS_H
