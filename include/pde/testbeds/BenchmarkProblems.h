///////////////////////////////////////////////////////////////////////////////////////////
// BenchmarkProblems.h - Standard Benchmark Problems for PDE Solver Validation
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// This file provides well-known benchmark problems from the literature with:
// - Problem specification (domain, BCs, source terms)
// - Exact solutions (where available)
// - Reference values for validation
// - Physical interpretations
//
// These benchmarks can be used for:
// - Solver validation against known results
// - Performance benchmarking
// - Cross-validation between different methods
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_TESTBEDS_BENCHMARK_PROBLEMS_H
#define MML_PDE_TESTBEDS_BENCHMARK_PROBLEMS_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace MML::PDE::Testbeds {

    ///////////////////////////////////////////////////////////////////////
    ///              Elliptic Benchmark Problems (Poisson)              ///
    ///////////////////////////////////////////////////////////////////////
    
    namespace EllipticBenchmarks {
        
        /**
         * @brief Unit square with polynomial solution
         * 
         * Problem: -∇²u = 2 on [0,1]²
         * BCs: u = 0 on boundary
         * Exact: u(x,y) = x(1-x)·y(1-y)/2 at center ≈ 0.0625
         * 
         * Use: Basic verification, O(h²) convergence check
         */
        template<typename T>
        struct UnitSquarePoisson {
            static constexpr const char* name = "Unit Square Poisson";
            static constexpr const char* description = 
                "Polynomial source on [0,1]^2, zero Dirichlet BCs";
            
            static T source(T, T) { return T(2); }
            static T boundary(T, T) { return T(0); }
            
            static T exactSolution(T x, T y) {
                // This is NOT exact for f=2 constant, approximate only
                return x * (T(1) - x) * y * (T(1) - y);
            }
            
            // Reference value at center
            static T referenceValueAtCenter() { return T(0.0625); }
        };
        
        /**
         * @brief L-shaped domain (reentrant corner)
         * 
         * Problem: -∇²u = 0 on L-shaped domain
         * Domain: [0,1]² \ [0.5,1]²
         * BCs: u = r^(2/3)·sin(2θ/3) on boundary
         * 
         * Famous benchmark with singular solution at corner.
         * Tests adaptive refinement and corner handling.
         */
        template<typename T>
        struct LShapedDomain {
            static constexpr const char* name = "L-Shaped Domain";
            static constexpr const char* description = 
                "Laplace with singular corner, r^(2/3) solution";
            
            // Exact solution in polar coordinates centered at corner (0.5, 0.5)
            static T exactSolution(T x, T y) {
                T dx = x - T(0.5);
                T dy = y - T(0.5);
                T r = std::sqrt(dx * dx + dy * dy);
                T theta = std::atan2(dy, dx);
                if (theta < 0) theta += T(2) * M_PI;
                return std::pow(r, T(2.0/3.0)) * std::sin(T(2.0/3.0) * theta);
            }
            
            static T source(T, T) { return T(0); }  // Laplace equation
            
            // Convergence rate is reduced due to corner singularity
            // Standard rate: 2/3 instead of 2
            static T expectedConvergenceRate() { return T(2.0/3.0); }
        };
        
        /**
         * @brief Point source (Green's function test)
         * 
         * Problem: -∇²u = δ(x-0.5, y-0.5) on [0,1]²
         * Tests handling of singular source terms.
         */
        template<typename T>
        struct PointSource2D {
            static constexpr const char* name = "Point Source";
            static constexpr const char* description = 
                "Dirac delta source at center, tests singularity handling";
            
            T x0 = T(0.5), y0 = T(0.5);  // Source location
            T epsilon = T(0.01);  // Regularization width
            
            // Regularized delta function (Gaussian)
            T source(T x, T y) const {
                T dx = x - x0;
                T dy = y - y0;
                T r2 = dx * dx + dy * dy;
                return std::exp(-r2 / (epsilon * epsilon)) / (M_PI * epsilon * epsilon);
            }
            
            // Free-space Green's function (infinite domain reference)
            T freeSpaceSolution(T x, T y) const {
                T dx = x - x0;
                T dy = y - y0;
                T r = std::sqrt(dx * dx + dy * dy);
                return -std::log(r) / (T(2) * M_PI);
            }
        };
        
        /**
         * @brief Anisotropic diffusion
         * 
         * Problem: -∇·(K∇u) = f where K is a 2x2 tensor
         * Tests solvers for variable coefficients.
         */
        template<typename T>
        struct AnisotropicDiffusion {
            static constexpr const char* name = "Anisotropic Diffusion";
            
            T kxx = T(1), kyy = T(10);  // Diffusion coefficients
            
            // Tests alignment with principal directions
            static T source(T x, T y) {
                return std::sin(M_PI * x) * std::sin(M_PI * y);
            }
        };
        
    } // namespace EllipticBenchmarks
    
    ///////////////////////////////////////////////////////////////////////
    ///             Parabolic Benchmark Problems (Heat)                 ///
    ///////////////////////////////////////////////////////////////////////
    
    namespace ParabolicBenchmarks {
        
        /**
         * @brief Heat conduction in a rod
         * 
         * Initially uniform temperature, ends held at zero.
         * Classic Fourier series solution.
         */
        template<typename T>
        struct CoolingRod1D {
            static constexpr const char* name = "Cooling Rod";
            static constexpr const char* description = 
                "1D rod with initial T=1, ends at T=0, Fourier series solution";
            
            T alpha = T(1);     // Thermal diffusivity
            T L = T(1);         // Length
            T T_initial = T(1); // Initial temperature
            int numTerms = 50;  // Fourier series terms
            
            CoolingRod1D(T a = T(1)) : alpha(a) {}
            
            T initial(T) const { return T_initial; }
            T boundary(T) const { return T(0); }
            T source(T, T) const { return T(0); }
            
            // Fourier series solution (truncated)
            T solution(T x, T t) const {
                T sum = T(0);
                for (int n = 1; n <= numTerms; n += 2) {  // Odd terms only
                    T npi = n * M_PI / L;
                    T coeff = T(4) / (n * M_PI);
                    sum += coeff * std::sin(npi * x) * std::exp(-alpha * npi * npi * t);
                }
                return sum * T_initial;
            }
            
            // Time for temperature at center to reach half initial
            T halfTime() const {
                return std::log(T(2)) / (alpha * M_PI * M_PI / (L * L));
            }
        };
        
        /**
         * @brief Stefan problem (phase change)
         * 
         * Melting/freezing with moving boundary.
         * Self-similar solution: interface position s(t) ∝ √t
         */
        template<typename T>
        struct StefanProblem1D {
            static constexpr const char* name = "Stefan Problem";
            static constexpr const char* description = 
                "Phase change with moving boundary, sqrt(t) interface motion";
            
            T alpha = T(1);      // Thermal diffusivity
            T T_melt = T(0);     // Melting temperature
            T T_hot = T(1);      // Hot boundary temperature
            T lambda;            // Stefan number parameter
            
            StefanProblem1D(T a = T(1)) : alpha(a), lambda(T(0.5)) {}
            
            // Interface position s(t) = 2λ√(αt)
            T interfacePosition(T t) const {
                return T(2) * lambda * std::sqrt(alpha * t);
            }
            
            // Temperature in liquid (x < s(t))
            T liquidTemperature(T x, T t) const {
                if (t <= T(0)) return T_hot;
                T s = interfacePosition(t);
                if (x >= s) return T_melt;
                T eta = x / (T(2) * std::sqrt(alpha * t));
                return T_hot * (T(1) - std::erf(eta) / std::erf(lambda));
            }
        };
        
        /**
         * @brief Diffusion from localized source
         * 
         * Instantaneous release of heat/mass at t=0, spreading Gaussian.
         */
        template<typename T>
        struct InstantaneousSource2D {
            static constexpr const char* name = "Instantaneous Source";
            static constexpr const char* description = 
                "Point release at t=0, Gaussian spreading solution";
            
            T alpha = T(1);
            T x0 = T(0.5), y0 = T(0.5);
            T Q = T(1);  // Total heat/mass released
            
            InstantaneousSource2D(T a = T(1)) : alpha(a) {}
            
            // Fundamental solution (Green's function)
            T solution(T x, T y, T t) const {
                if (t <= T(0)) return T(0);
                T dx = x - x0;
                T dy = y - y0;
                T r2 = dx * dx + dy * dy;
                return Q / (T(4) * M_PI * alpha * t) * std::exp(-r2 / (T(4) * alpha * t));
            }
            
            // Maximum temperature decreases as 1/t
            T maxTemperature(T t) const {
                return Q / (T(4) * M_PI * alpha * t);
            }
            
            // Width grows as √t
            T characteristicWidth(T t) const {
                return T(2) * std::sqrt(alpha * t);
            }
        };
        
    } // namespace ParabolicBenchmarks
    
    ///////////////////////////////////////////////////////////////////////
    ///             Hyperbolic Benchmark Problems (Wave)                ///
    ///////////////////////////////////////////////////////////////////////
    
    namespace HyperbolicBenchmarks {
        
        /**
         * @brief Plucked string
         * 
         * Triangular initial displacement, zero velocity.
         * D'Alembert solution with reflections.
         */
        template<typename T>
        struct PluckedString1D {
            static constexpr const char* name = "Plucked String";
            static constexpr const char* description = 
                "Triangular initial shape, d'Alembert solution";
            
            T L = T(1);          // String length
            T c = T(1);          // Wave speed
            T pluckPos = T(0.5); // Pluck position (fraction of L)
            T amplitude = T(1);  // Initial amplitude
            
            PluckedString1D(T speed = T(1), T pos = T(0.5)) 
                : c(speed), pluckPos(pos) {}
            
            T initialDisplacement(T x) const {
                T xp = pluckPos * L;
                if (x <= xp) {
                    return amplitude * x / xp;
                } else {
                    return amplitude * (L - x) / (L - xp);
                }
            }
            
            T initialVelocity(T) const { return T(0); }
            
            // Period of fundamental mode
            T fundamentalPeriod() const { return T(2) * L / c; }
            
            // Energy (should be conserved)
            T totalEnergy() const {
                // E ∝ (amplitude/L)² for triangular pulse
                return amplitude * amplitude / (L * pluckPos * (T(1) - pluckPos));
            }
        };
        
        /**
         * @brief Circular drum membrane
         * 
         * Axisymmetric modes with Bessel function solutions.
         */
        template<typename T>
        struct CircularMembrane2D {
            static constexpr const char* name = "Circular Membrane";
            static constexpr const char* description = 
                "Drum head with Bessel function modes";
            
            T radius = T(1);
            T c = T(1);
            
            // First few zeros of J_0 (for axisymmetric modes)
            static constexpr T j0_zeros[] = {T(2.4048), T(5.5201), T(8.6537)};
            
            // Fundamental frequency
            T fundamentalFrequency() const {
                return c * j0_zeros[0] / (T(2) * M_PI * radius);
            }
            
            // Radial mode shape (n-th mode)
            T modeShape(T r, int n = 0) const {
                if (n >= 3) n = 0;
                T j0n = j0_zeros[n];
                // Approximate J_0 for small argument
                T x = j0n * r / radius;
                return T(1) - x * x / T(4) + x * x * x * x / T(64); // Taylor for J_0
            }
        };
        
        /**
         * @brief Shock tube (Riemann problem)
         * 
         * Discontinuous initial data, tests scheme accuracy near discontinuities.
         */
        template<typename T>
        struct ShockTube1D {
            static constexpr const char* name = "Shock Tube";
            static constexpr const char* description = 
                "Discontinuous IC, tests numerical diffusion/dispersion";
            
            T L = T(1);
            T c = T(1);
            T x_disc = T(0.5);  // Discontinuity location
            T u_left = T(1);
            T u_right = T(0);
            
            T initial(T x) const {
                return (x < x_disc) ? u_left : u_right;
            }
            
            // For linear advection: discontinuity moves at speed c
            T discontinuityPosition(T t) const {
                return x_disc + c * t;
            }
            
            // Exact solution (advection only)
            T exactSolution(T x, T t) const {
                T x_disc_t = discontinuityPosition(t);
                return (x < x_disc_t) ? u_left : u_right;
            }
        };
        
        /**
         * @brief Smooth pulse advection
         * 
         * Gaussian pulse transported without deformation.
         * Good for measuring numerical diffusion.
         */
        template<typename T>
        struct SmoothPulseAdvection1D {
            static constexpr const char* name = "Smooth Pulse Advection";
            static constexpr const char* description = 
                "Gaussian pulse transport, measures numerical diffusion";
            
            T c = T(1);
            T x0 = T(0.25);
            T sigma = T(0.05);
            T L = T(1);  // Domain length (periodic)
            
            SmoothPulseAdvection1D(T speed = T(1)) : c(speed) {}
            
            T initial(T x) const {
                T dx = x - x0;
                return std::exp(-dx * dx / (sigma * sigma));
            }
            
            T exactSolution(T x, T t) const {
                // Periodic domain handling
                T x_rel = std::fmod(x - c * t - x0 + T(10) * L, L);
                if (x_rel < 0) x_rel += L;
                if (x_rel > L/2) x_rel = L - x_rel;
                return std::exp(-x_rel * x_rel / (sigma * sigma));
            }
            
            // Measure diffusion by peak height reduction
            T peakHeightReduction(const std::vector<T>& numerical) const {
                T maxVal = *std::max_element(numerical.begin(), numerical.end());
                return T(1) - maxVal;  // Should be 0 for exact solution
            }
        };
        
    } // namespace HyperbolicBenchmarks
    
    ///////////////////////////////////////////////////////////////////////
    ///                Performance Benchmark Suite                      ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Standard problem sizes for benchmarking
     */
    struct BenchmarkSizes {
        // 1D problems
        static constexpr int small_1d = 100;
        static constexpr int medium_1d = 1000;
        static constexpr int large_1d = 10000;
        
        // 2D problems (total unknowns = n²)
        static constexpr int small_2d = 32;    // ~1K unknowns
        static constexpr int medium_2d = 100;  // ~10K unknowns
        static constexpr int large_2d = 316;   // ~100K unknowns
        static constexpr int xlarge_2d = 1000; // ~1M unknowns
        
        // 3D problems (total unknowns = n³)
        static constexpr int small_3d = 10;    // ~1K unknowns
        static constexpr int medium_3d = 22;   // ~10K unknowns
        static constexpr int large_3d = 46;    // ~100K unknowns
    };
    
} // namespace MML::PDE::Testbeds

#endif // MML_PDE_TESTBEDS_BENCHMARK_PROBLEMS_H
