///////////////////////////////////////////////////////////////////////////////////////////
// SingleObjective_TestProblems.h
//
// Standard unconstrained benchmark problems for single-objective optimization.
//
// Classic test functions used for validating and comparing optimization algorithms.
// All problems have known global optima for verification.
//
// Categories:
// - Unimodal: Sphere, Rosenbrock (single global minimum, good for convergence tests)
// - Multimodal: Rastrigin, Ackley, Griewank, Schwefel (many local minima)
// - Other: Beale, Booth, Matyas, Himmelblau (various difficulties)
//
// References:
// - Jamil, M., & Yang, X. S. (2013). A literature survey of benchmark functions for 
//   global optimisation problems. arXiv preprint arXiv:1308.4008.
// - https://www.sfu.ca/~ssurjano/optimization.html
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SINGLE_OBJECTIVE_TEST_PROBLEMS_H
#define MML_SINGLE_OBJECTIVE_TEST_PROBLEMS_H

#include "MMLBase.h"
#include "OptimizationProblem.h"

#include <cmath>

namespace MML
{
namespace TestProblems
{
    using namespace Optimization;  // For SingleObjectiveProblem, etc.

///////////////////////////////////////////////////////////////////////////////////////////
//                              UNIMODAL FUNCTIONS                                       //
///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////
// Sphere Function
//
// f(x) = Σ x_i²
//
// Properties:
// - Unimodal, convex, separable
// - Global minimum: f(0,...,0) = 0
// - Typical bounds: [-5.12, 5.12]^N
// - Easiest test function, good for sanity checks
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Sphere : public SingleObjectiveProblem
{
public:
    Sphere() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            for (int i = 0; i < N; ++i)
                sum += x[i] * x[i];
            return sum;
        })
    {
        this->SetBounds(-5.12, 5.12);
    }

    /// @brief Get the known global optimum location
    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 0;
        return opt;
    }

    /// @brief Get the known global minimum value
    static Real GetOptimalValue() { return 0; }

    /// @brief Get problem name
    static std::string GetName() { return "Sphere"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Rosenbrock Function (Banana Function)
//
// f(x) = Σ [100*(x_{i+1} - x_i²)² + (1 - x_i)²]
//
// Properties:
// - Unimodal for N ≤ 3, multimodal for N > 3
// - Global minimum: f(1,...,1) = 0
// - Typical bounds: [-5, 10]^N or [-2.048, 2.048]^N
// - Difficult due to narrow curved valley
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Rosenbrock : public SingleObjectiveProblem
{
    static_assert(N >= 2, "Rosenbrock requires at least 2 variables");

public:
    Rosenbrock() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            for (int i = 0; i < N - 1; ++i) {
                Real t1 = x[i + 1] - x[i] * x[i];
                Real t2 = 1 - x[i];
                sum += 100 * t1 * t1 + t2 * t2;
            }
            return sum;
        })
    {
        this->SetBounds(-5, 10);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 1;
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Rosenbrock"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Sum of Different Powers Function
//
// f(x) = Σ |x_i|^(i+1)
//
// Properties:
// - Unimodal, non-separable
// - Global minimum: f(0,...,0) = 0
// - Typical bounds: [-1, 1]^N
// - Sensitivity varies by dimension
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class SumOfDifferentPowers : public SingleObjectiveProblem
{
public:
    SumOfDifferentPowers() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            for (int i = 0; i < N; ++i)
                sum += std::pow(std::abs(x[i]), i + 2);
            return sum;
        })
    {
        this->SetBounds(-1, 1);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 0;
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "SumOfDifferentPowers"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
//                              MULTIMODAL FUNCTIONS                                     //
///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////
// Rastrigin Function
//
// f(x) = 10*N + Σ [x_i² - 10*cos(2πx_i)]
//
// Properties:
// - Highly multimodal with regular local minima
// - Global minimum: f(0,...,0) = 0
// - Typical bounds: [-5.12, 5.12]^N
// - Difficult due to massive number of local minima (approximately 11^N)
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Rastrigin : public SingleObjectiveProblem
{
public:
    Rastrigin() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            const Real A = 10;
            Real sum = A * N;
            for (int i = 0; i < N; ++i) {
                sum += x[i] * x[i] - A * std::cos(2 * Constants::PI * x[i]);
            }
            return sum;
        })
    {
        this->SetBounds(-5.12, 5.12);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 0;
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Rastrigin"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Ackley Function
//
// f(x) = -20*exp(-0.2*sqrt(Σx_i²/N)) - exp(Σcos(2πx_i)/N) + 20 + e
//
// Properties:
// - Multimodal with nearly flat outer region
// - Global minimum: f(0,...,0) = 0
// - Typical bounds: [-32.768, 32.768]^N
// - Central funnel surrounded by nearly flat region
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Ackley : public SingleObjectiveProblem
{
public:
    Ackley() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum1 = 0, sum2 = 0;
            for (int i = 0; i < N; ++i) {
                sum1 += x[i] * x[i];
                sum2 += std::cos(2 * Constants::PI * x[i]);
            }
            return -20 * std::exp(-0.2 * std::sqrt(sum1 / N)) 
                   - std::exp(sum2 / N) + 20 + Constants::E;
        })
    {
        this->SetBounds(-32.768, 32.768);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 0;
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Ackley"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Griewank Function
//
// f(x) = 1 + Σ(x_i²/4000) - Π cos(x_i/√i)
//
// Properties:
// - Multimodal with product term creating correlations
// - Global minimum: f(0,...,0) = 0
// - Typical bounds: [-600, 600]^N
// - Many widespread local minima
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Griewank : public SingleObjectiveProblem
{
public:
    Griewank() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            Real prod = 1;
            for (int i = 0; i < N; ++i) {
                sum += x[i] * x[i] / 4000;
                prod *= std::cos(x[i] / std::sqrt(i + 1));
            }
            return 1 + sum - prod;
        })
    {
        this->SetBounds(-600, 600);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 0;
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Griewank"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Schwefel Function
//
// f(x) = 418.9829*N - Σ x_i*sin(√|x_i|)
//
// Properties:
// - Multimodal with global minimum far from next best local minima
// - Global minimum: f(420.9687,...,420.9687) ≈ 0
// - Typical bounds: [-500, 500]^N
// - Deceptive: best local minima are distant from global optimum
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Schwefel : public SingleObjectiveProblem
{
public:
    Schwefel() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            for (int i = 0; i < N; ++i) {
                sum += x[i] * std::sin(std::sqrt(std::abs(x[i])));
            }
            return 418.9829 * N - sum;
        })
    {
        this->SetBounds(-500, 500);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 420.9687;
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Schwefel"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Levy Function
//
// f(x) = sin²(πw₁) + Σ(wᵢ-1)²[1+10sin²(πwᵢ+1)] + (wₙ-1)²[1+sin²(2πwₙ)]
// where wᵢ = 1 + (xᵢ-1)/4
//
// Properties:
// - Multimodal
// - Global minimum: f(1,...,1) = 0
// - Typical bounds: [-10, 10]^N
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Levy : public SingleObjectiveProblem
{
public:
    Levy() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            auto w = [](Real xi) { return 1 + (xi - 1) / 4; };
            
            Real w1 = w(x[0]);
            Real wn = w(x[N-1]);
            
            Real term1 = std::pow(std::sin(Constants::PI * w1), 2);
            
            Real sum = 0;
            for (int i = 0; i < N - 1; ++i) {
                Real wi = w(x[i]);
                sum += std::pow(wi - 1, 2) * (1 + 10 * std::pow(std::sin(Constants::PI * wi + 1), 2));
            }
            
            Real term3 = std::pow(wn - 1, 2) * (1 + std::pow(std::sin(2 * Constants::PI * wn), 2));
            
            return term1 + sum + term3;
        })
    {
        this->SetBounds(-10, 10);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 1;
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Levy"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
//                           TWO-DIMENSIONAL FUNCTIONS                                   //
///////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////
// Beale Function (2D only)
//
// f(x,y) = (1.5-x+xy)² + (2.25-x+xy²)² + (2.625-x+xy³)²
//
// Properties:
// - Multimodal
// - Global minimum: f(3, 0.5) = 0
// - Typical bounds: [-4.5, 4.5]²
///////////////////////////////////////////////////////////////////////////////////////////

class Beale : public SingleObjectiveProblem
{
public:
    Beale() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            Real t1 = 1.5 - x[0] + x[0] * x[1];
            Real t2 = 2.25 - x[0] + x[0] * x[1] * x[1];
            Real t3 = 2.625 - x[0] + x[0] * x[1] * x[1] * x[1];
            return t1 * t1 + t2 * t2 + t3 * t3;
        })
    {
        this->SetBounds(-4.5, 4.5);
    }

    static Vector<Real> GetOptimum() {
        return Vector<Real>{3, 0.5};
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Beale"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Booth Function (2D only)
//
// f(x,y) = (x + 2y - 7)² + (2x + y - 5)²
//
// Properties:
// - Unimodal
// - Global minimum: f(1, 3) = 0
// - Typical bounds: [-10, 10]²
///////////////////////////////////////////////////////////////////////////////////////////

class Booth : public SingleObjectiveProblem
{
public:
    Booth() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            Real t1 = x[0] + 2 * x[1] - 7;
            Real t2 = 2 * x[0] + x[1] - 5;
            return t1 * t1 + t2 * t2;
        })
    {
        this->SetBounds(-10, 10);
    }

    static Vector<Real> GetOptimum() {
        return Vector<Real>{1, 3};
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Booth"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Matyas Function (2D only)
//
// f(x,y) = 0.26*(x² + y²) - 0.48*x*y
//
// Properties:
// - Unimodal
// - Global minimum: f(0, 0) = 0
// - Typical bounds: [-10, 10]²
///////////////////////////////////////////////////////////////////////////////////////////

class Matyas : public SingleObjectiveProblem
{
public:
    Matyas() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            return 0.26 * (x[0] * x[0] + x[1] * x[1]) - 0.48 * x[0] * x[1];
        })
    {
        this->SetBounds(-10, 10);
    }

    static Vector<Real> GetOptimum() {
        return Vector<Real>{0, 0};
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Matyas"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Himmelblau Function (2D only)
//
// f(x,y) = (x² + y - 11)² + (x + y² - 7)²
//
// Properties:
// - Multimodal with 4 identical global minima
// - Global minima:
//   f(3, 2) = 0
//   f(-2.805118, 3.131312) = 0
//   f(-3.779310, -3.283186) = 0
//   f(3.584428, -1.848126) = 0
// - Typical bounds: [-5, 5]²
///////////////////////////////////////////////////////////////////////////////////////////

class Himmelblau : public SingleObjectiveProblem
{
public:
    Himmelblau() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            Real t1 = x[0] * x[0] + x[1] - 11;
            Real t2 = x[0] + x[1] * x[1] - 7;
            return t1 * t1 + t2 * t2;
        })
    {
        this->SetBounds(-5, 5);
    }

    /// @brief Get one of the 4 global optima (returns the first one)
    static Vector<Real> GetOptimum() {
        return Vector<Real>{3, 2};
    }

    /// @brief Get all 4 global optima
    static std::array<Vector<Real>, 4> GetAllOptima() {
        return {{
            Vector<Real>{3.0, 2.0},
            Vector<Real>{-2.805118, 3.131312},
            Vector<Real>{-3.779310, -3.283186},
            Vector<Real>{3.584428, -1.848126}
        }};
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Himmelblau"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Goldstein-Price Function (2D only)
//
// f(x,y) = [1 + (x+y+1)²(19-14x+3x²-14y+6xy+3y²)] *
//          [30 + (2x-3y)²(18-32x+12x²+48y-36xy+27y²)]
//
// Properties:
// - Multimodal
// - Global minimum: f(0, -1) = 3
// - Typical bounds: [-2, 2]²
///////////////////////////////////////////////////////////////////////////////////////////

class GoldsteinPrice : public SingleObjectiveProblem
{
public:
    GoldsteinPrice() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            Real x1 = x[0], x2 = x[1];
            
            Real t1 = x1 + x2 + 1;
            Real t2 = 19 - 14*x1 + 3*x1*x1 - 14*x2 + 6*x1*x2 + 3*x2*x2;
            Real factor1 = 1 + t1*t1 * t2;
            
            Real t3 = 2*x1 - 3*x2;
            Real t4 = 18 - 32*x1 + 12*x1*x1 + 48*x2 - 36*x1*x2 + 27*x2*x2;
            Real factor2 = 30 + t3*t3 * t4;
            
            return factor1 * factor2;
        })
    {
        this->SetBounds(-2, 2);
    }

    static Vector<Real> GetOptimum() {
        return Vector<Real>{0, -1};
    }

    static Real GetOptimalValue() { return 3; }
    static std::string GetName() { return "GoldsteinPrice"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Six-Hump Camel Function (2D only)
//
// f(x,y) = (4 - 2.1x² + x⁴/3)x² + xy + (-4 + 4y²)y²
//
// Properties:
// - Multimodal with 6 local minima (2 global)
// - Global minima: f(±0.0898, ∓0.7126) ≈ -1.0316
// - Typical bounds: x ∈ [-3, 3], y ∈ [-2, 2]
///////////////////////////////////////////////////////////////////////////////////////////

class SixHumpCamel : public SingleObjectiveProblem
{
public:
    SixHumpCamel() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            Real x1 = x[0], x2 = x[1];
            return (4 - 2.1*x1*x1 + std::pow(x1, 4)/3) * x1*x1 
                   + x1*x2 
                   + (-4 + 4*x2*x2) * x2*x2;
        })
    {
        // Different bounds for x and y
        Vector<Real> lower{-3, -2};
        Vector<Real> upper{3, 2};
        this->SetBounds(lower, upper);
    }

    static Vector<Real> GetOptimum() {
        return Vector<Real>{0.0898, -0.7126};
    }

    /// @brief Get both global optima
    static std::array<Vector<Real>, 2> GetAllOptima() {
        return {{
            Vector<Real>{0.0898, -0.7126},
            Vector<Real>{-0.0898, 0.7126}
        }};
    }

    static Real GetOptimalValue() { return -1.0316284534898; }
    static std::string GetName() { return "SixHumpCamel"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Three-Hump Camel Function (2D only)
//
// f(x,y) = 2x² - 1.05x⁴ + x⁶/6 + xy + y²
//
// Properties:
// - Multimodal
// - Global minimum: f(0, 0) = 0
// - Typical bounds: [-5, 5]²
///////////////////////////////////////////////////////////////////////////////////////////

class ThreeHumpCamel : public SingleObjectiveProblem
{
public:
    ThreeHumpCamel() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            Real x1 = x[0], x2 = x[1];
            return 2*x1*x1 - 1.05*std::pow(x1, 4) + std::pow(x1, 6)/6 + x1*x2 + x2*x2;
        })
    {
        this->SetBounds(-5, 5);
    }

    static Vector<Real> GetOptimum() {
        return Vector<Real>{0, 0};
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "ThreeHumpCamel"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Easom Function (2D only)
//
// f(x,y) = -cos(x)*cos(y)*exp(-((x-π)² + (y-π)²))
//
// Properties:
// - Unimodal with large nearly flat region
// - Global minimum: f(π, π) = -1
// - Typical bounds: [-100, 100]²
// - Difficult due to small attraction basin
///////////////////////////////////////////////////////////////////////////////////////////

class Easom : public SingleObjectiveProblem
{
public:
    Easom() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            Real t1 = (x[0] - Constants::PI) * (x[0] - Constants::PI);
            Real t2 = (x[1] - Constants::PI) * (x[1] - Constants::PI);
            return -std::cos(x[0]) * std::cos(x[1]) * std::exp(-(t1 + t2));
        })
    {
        this->SetBounds(-100, 100);
    }

    static Vector<Real> GetOptimum() {
        return Vector<Real>{Constants::PI, Constants::PI};
    }

    static Real GetOptimalValue() { return -1; }
    static std::string GetName() { return "Easom"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// McCormick Function (2D only)
//
// f(x,y) = sin(x+y) + (x-y)² - 1.5x + 2.5y + 1
//
// Properties:
// - Multimodal
// - Global minimum: f(-0.54719, -1.54719) ≈ -1.9133
// - Typical bounds: x ∈ [-1.5, 4], y ∈ [-3, 4]
///////////////////////////////////////////////////////////////////////////////////////////

class McCormick : public SingleObjectiveProblem
{
public:
    McCormick() : SingleObjectiveProblem(2,
        [](const Vector<Real>& x) -> Real {
            return std::sin(x[0] + x[1]) + std::pow(x[0] - x[1], 2) 
                   - 1.5*x[0] + 2.5*x[1] + 1;
        })
    {
        Vector<Real> lower{-1.5, -3};
        Vector<Real> upper{4, 4};
        this->SetBounds(lower, upper);
    }

    static Vector<Real> GetOptimum() {
        return Vector<Real>{-0.54719, -1.54719};
    }

    static Real GetOptimalValue() { return -1.9133; }
    static std::string GetName() { return "McCormick"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Styblinski-Tang Function
//
// f(x) = 0.5 * Σ (x_i⁴ - 16x_i² + 5x_i)
//
// Properties:
// - Multimodal
// - Global minimum: f(-2.903534,...,-2.903534) ≈ -39.16599*N
// - Typical bounds: [-5, 5]^N
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class StyblinskiTang : public SingleObjectiveProblem
{
public:
    StyblinskiTang() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            for (int i = 0; i < N; ++i) {
                Real xi = x[i];
                sum += std::pow(xi, 4) - 16*xi*xi + 5*xi;
            }
            return 0.5 * sum;
        })
    {
        this->SetBounds(-5, 5);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = -2.903534;
        return opt;
    }

    static Real GetOptimalValue() { return -39.16599 * N; }
    static std::string GetName() { return "StyblinskiTang"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Dixon-Price Function
//
// f(x) = (x₁-1)² + Σ i*(2x_i² - x_{i-1})²
//
// Properties:
// - Unimodal
// - Global minimum: x_i = 2^(-(2^i-2)/2^i) with f* = 0
// - Typical bounds: [-10, 10]^N
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class DixonPrice : public SingleObjectiveProblem
{
public:
    DixonPrice() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum = (x[0] - 1) * (x[0] - 1);
            for (int i = 1; i < N; ++i) {
                Real t = 2 * x[i] * x[i] - x[i-1];
                sum += (i + 1) * t * t;
            }
            return sum;
        })
    {
        this->SetBounds(-10, 10);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) {
            opt[i] = std::pow(2.0, -((std::pow(2.0, i+1) - 2) / std::pow(2.0, i+1)));
        }
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "DixonPrice"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Zakharov Function
//
// f(x) = Σx_i² + (0.5*Σi*x_i)² + (0.5*Σi*x_i)⁴
//
// Properties:
// - Unimodal
// - Global minimum: f(0,...,0) = 0
// - Typical bounds: [-5, 10]^N
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Zakharov : public SingleObjectiveProblem
{
public:
    Zakharov() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum1 = 0;  // Σx_i²
            Real sum2 = 0;  // 0.5*Σi*x_i
            for (int i = 0; i < N; ++i) {
                sum1 += x[i] * x[i];
                sum2 += 0.5 * (i + 1) * x[i];
            }
            return sum1 + sum2*sum2 + sum2*sum2*sum2*sum2;
        })
    {
        this->SetBounds(-5, 10);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i) opt[i] = 0;
        return opt;
    }

    static Real GetOptimalValue() { return 0; }
    static std::string GetName() { return "Zakharov"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Trid Function
//
// f(x) = Σ(x_i - 1)² - Σx_i*x_{i-1}
//
// Properties:
// - Unimodal
// - Global minimum: x_i = i(N+1-i) with f* = -N(N+4)(N-1)/6
// - Typical bounds: [-N², N²]^N
///////////////////////////////////////////////////////////////////////////////////////////

template<int N>
class Trid : public SingleObjectiveProblem
{
public:
    Trid() : SingleObjectiveProblem(N,
        [](const Vector<Real>& x) -> Real {
            Real sum1 = 0;  // Σ(x_i - 1)²
            Real sum2 = 0;  // Σx_i*x_{i-1}
            for (int i = 0; i < N; ++i) {
                sum1 += (x[i] - 1) * (x[i] - 1);
                if (i > 0)
                    sum2 += x[i] * x[i-1];
            }
            return sum1 - sum2;
        })
    {
        Real bound = static_cast<Real>(N * N);
        this->SetBounds(-bound, bound);
    }

    static Vector<Real> GetOptimum() {
        Vector<Real> opt(N);
        for (int i = 0; i < N; ++i)
            opt[i] = (i + 1) * (N + 1 - (i + 1));
        return opt;
    }

    static Real GetOptimalValue() { 
        return -static_cast<Real>(N * (N + 4) * (N - 1)) / 6; 
    }
    static std::string GetName() { return "Trid"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
//               DYNAMIC VECTOR VERSIONS FOR HEURISTIC OPTIMIZERS                        //
//                                                                                       //
// These free functions return std::function objects compatible with Vector<Real>        //
// (dynamic size vectors) for use with SimulatedAnnealing, ParticleSwarm, etc.           //
///////////////////////////////////////////////////////////////////////////////////////////

namespace Functions
{
    /// @brief Sphere function: f(x) = Σ x_i²
    /// Global minimum at origin with value 0
    inline std::function<Real(const Vector<Real>&)> Sphere() {
        return [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            for (int i = 0; i < x.size(); ++i)
                sum += x[i] * x[i];
            return sum;
        };
    }

    /// @brief Rosenbrock function: f(x) = Σ [100*(x_{i+1} - x_i²)² + (1 - x_i)²]
    /// Global minimum at (1,...,1) with value 0
    inline std::function<Real(const Vector<Real>&)> Rosenbrock() {
        return [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            for (int i = 0; i < x.size() - 1; ++i) {
                Real t1 = x[i + 1] - x[i] * x[i];
                Real t2 = 1 - x[i];
                sum += 100 * t1 * t1 + t2 * t2;
            }
            return sum;
        };
    }

    /// @brief Rastrigin function: f(x) = 10*N + Σ [x_i² - 10*cos(2πx_i)]
    /// Global minimum at origin with value 0
    inline std::function<Real(const Vector<Real>&)> Rastrigin() {
        return [](const Vector<Real>& x) -> Real {
            const Real A = 10;
            int n = x.size();
            Real sum = A * n;
            for (int i = 0; i < n; ++i) {
                sum += x[i] * x[i] - A * std::cos(2 * Constants::PI * x[i]);
            }
            return sum;
        };
    }

    /// @brief Ackley function
    /// Global minimum at origin with value 0
    inline std::function<Real(const Vector<Real>&)> Ackley() {
        return [](const Vector<Real>& x) -> Real {
            int n = x.size();
            Real sum1 = 0, sum2 = 0;
            for (int i = 0; i < n; ++i) {
                sum1 += x[i] * x[i];
                sum2 += std::cos(2 * Constants::PI * x[i]);
            }
            return -20 * std::exp(-0.2 * std::sqrt(sum1 / n)) 
                   - std::exp(sum2 / n) + 20 + Constants::E;
        };
    }

    /// @brief Griewank function
    /// Global minimum at origin with value 0
    inline std::function<Real(const Vector<Real>&)> Griewank() {
        return [](const Vector<Real>& x) -> Real {
            int n = x.size();
            Real sum = 0;
            Real prod = 1;
            for (int i = 0; i < n; ++i) {
                sum += x[i] * x[i] / 4000;
                prod *= std::cos(x[i] / std::sqrt(i + 1));
            }
            return 1 + sum - prod;
        };
    }

    /// @brief Schwefel function
    /// Global minimum at (420.9687,...,420.9687) with value ≈ 0
    inline std::function<Real(const Vector<Real>&)> Schwefel() {
        return [](const Vector<Real>& x) -> Real {
            int n = x.size();
            Real sum = 0;
            for (int i = 0; i < n; ++i) {
                sum += x[i] * std::sin(std::sqrt(std::abs(x[i])));
            }
            return 418.9829 * n - sum;
        };
    }

    /// @brief Levy function
    /// Global minimum at (1,...,1) with value 0
    inline std::function<Real(const Vector<Real>&)> Levy() {
        return [](const Vector<Real>& x) -> Real {
            auto w = [](Real xi) { return 1 + (xi - 1) / 4; };
            int n = x.size();
            
            Real w1 = w(x[0]);
            Real wn = w(x[n-1]);
            
            Real term1 = std::pow(std::sin(Constants::PI * w1), 2);
            
            Real sum = 0;
            for (int i = 0; i < n - 1; ++i) {
                Real wi = w(x[i]);
                sum += std::pow(wi - 1, 2) * (1 + 10 * std::pow(std::sin(Constants::PI * wi + 1), 2));
            }
            
            Real term3 = std::pow(wn - 1, 2) * (1 + std::pow(std::sin(2 * Constants::PI * wn), 2));
            
            return term1 + sum + term3;
        };
    }

    /// @brief Styblinski-Tang function
    /// Global minimum at (-2.903534,...,-2.903534) with value ≈ -39.16599*N
    inline std::function<Real(const Vector<Real>&)> StyblinskiTang() {
        return [](const Vector<Real>& x) -> Real {
            Real sum = 0;
            for (int i = 0; i < x.size(); ++i) {
                Real xi = x[i];
                sum += std::pow(xi, 4) - 16*xi*xi + 5*xi;
            }
            return 0.5 * sum;
        };
    }

    /// @brief Dixon-Price function
    /// Global minimum with value 0
    inline std::function<Real(const Vector<Real>&)> DixonPrice() {
        return [](const Vector<Real>& x) -> Real {
            Real sum = (x[0] - 1) * (x[0] - 1);
            for (int i = 1; i < x.size(); ++i) {
                Real t = 2 * x[i] * x[i] - x[i-1];
                sum += (i + 1) * t * t;
            }
            return sum;
        };
    }

    /// @brief Zakharov function
    /// Global minimum at origin with value 0
    inline std::function<Real(const Vector<Real>&)> Zakharov() {
        return [](const Vector<Real>& x) -> Real {
            Real sum1 = 0;
            Real sum2 = 0;
            for (int i = 0; i < x.size(); ++i) {
                sum1 += x[i] * x[i];
                sum2 += 0.5 * (i + 1) * x[i];
            }
            return sum1 + sum2*sum2 + sum2*sum2*sum2*sum2;
        };
    }

    //-------------------------------------------------------------------------------------
    // 2D-only functions
    //-------------------------------------------------------------------------------------

    /// @brief Beale function (2D)
    /// Global minimum at (3, 0.5) with value 0
    inline std::function<Real(const Vector<Real>&)> Beale() {
        return [](const Vector<Real>& x) -> Real {
            Real t1 = 1.5 - x[0] + x[0] * x[1];
            Real t2 = 2.25 - x[0] + x[0] * x[1] * x[1];
            Real t3 = 2.625 - x[0] + x[0] * x[1] * x[1] * x[1];
            return t1 * t1 + t2 * t2 + t3 * t3;
        };
    }

    /// @brief Booth function (2D)
    /// Global minimum at (1, 3) with value 0
    inline std::function<Real(const Vector<Real>&)> Booth() {
        return [](const Vector<Real>& x) -> Real {
            Real t1 = x[0] + 2 * x[1] - 7;
            Real t2 = 2 * x[0] + x[1] - 5;
            return t1 * t1 + t2 * t2;
        };
    }

    /// @brief Matyas function (2D)
    /// Global minimum at (0, 0) with value 0
    inline std::function<Real(const Vector<Real>&)> Matyas() {
        return [](const Vector<Real>& x) -> Real {
            return 0.26 * (x[0] * x[0] + x[1] * x[1]) - 0.48 * x[0] * x[1];
        };
    }

    /// @brief Himmelblau function (2D)
    /// Has 4 global minima, all with value 0
    inline std::function<Real(const Vector<Real>&)> Himmelblau() {
        return [](const Vector<Real>& x) -> Real {
            Real t1 = x[0] * x[0] + x[1] - 11;
            Real t2 = x[0] + x[1] * x[1] - 7;
            return t1 * t1 + t2 * t2;
        };
    }

    /// @brief Goldstein-Price function (2D)
    /// Global minimum at (0, -1) with value 3
    inline std::function<Real(const Vector<Real>&)> GoldsteinPrice() {
        return [](const Vector<Real>& x) -> Real {
            Real x1 = x[0], x2 = x[1];
            
            Real t1 = x1 + x2 + 1;
            Real t2 = 19 - 14*x1 + 3*x1*x1 - 14*x2 + 6*x1*x2 + 3*x2*x2;
            Real factor1 = 1 + t1*t1 * t2;
            
            Real t3 = 2*x1 - 3*x2;
            Real t4 = 18 - 32*x1 + 12*x1*x1 + 48*x2 - 36*x1*x2 + 27*x2*x2;
            Real factor2 = 30 + t3*t3 * t4;
            
            return factor1 * factor2;
        };
    }

    /// @brief Six-Hump Camel function (2D)
    /// Global minimum at (0.0898, -0.7126) with value ≈ -1.0316
    inline std::function<Real(const Vector<Real>&)> SixHumpCamel() {
        return [](const Vector<Real>& x) -> Real {
            Real x1 = x[0], x2 = x[1];
            return (4 - 2.1*x1*x1 + std::pow(x1, 4)/3) * x1*x1 
                   + x1*x2 
                   + (-4 + 4*x2*x2) * x2*x2;
        };
    }

    /// @brief Three-Hump Camel function (2D)
    /// Global minimum at (0, 0) with value 0
    inline std::function<Real(const Vector<Real>&)> ThreeHumpCamel() {
        return [](const Vector<Real>& x) -> Real {
            Real x1 = x[0], x2 = x[1];
            return 2*x1*x1 - 1.05*std::pow(x1, 4) + std::pow(x1, 6)/6 + x1*x2 + x2*x2;
        };
    }

    /// @brief Easom function (2D)
    /// Global minimum at (π, π) with value -1
    inline std::function<Real(const Vector<Real>&)> Easom() {
        return [](const Vector<Real>& x) -> Real {
            Real t1 = (x[0] - Constants::PI) * (x[0] - Constants::PI);
            Real t2 = (x[1] - Constants::PI) * (x[1] - Constants::PI);
            return -std::cos(x[0]) * std::cos(x[1]) * std::exp(-(t1 + t2));
        };
    }

    /// @brief McCormick function (2D)
    /// Global minimum at (-0.54719, -1.54719) with value ≈ -1.9133
    inline std::function<Real(const Vector<Real>&)> McCormick() {
        return [](const Vector<Real>& x) -> Real {
            return std::sin(x[0] + x[1]) + std::pow(x[0] - x[1], 2) 
                   - 1.5*x[0] + 2.5*x[1] + 1;
        };
    }

} // namespace Functions

} // namespace TestProblems
} // namespace MML

#endif // MML_SINGLE_OBJECTIVE_TEST_PROBLEMS_H
