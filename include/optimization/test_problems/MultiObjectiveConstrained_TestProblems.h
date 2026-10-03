///////////////////////////////////////////////////////////////////////////////////////////
// MultiObjectiveConstrained_TestProblems.h
//
// Standard benchmark problems for multi-objective optimization.
//
// ZDT Series: Classic bi-objective test problems from:
//   Zitzler, E., Deb, K., & Thiele, L. (2000). Comparison of multiobjective evolutionary
//   algorithms: Empirical results. Evolutionary computation, 8(2), 173-195.
//
// All ZDT problems have:
// - 2 objectives to minimize
// - Known Pareto front for validation
// - Various difficulty features (convex/concave front, multimodality, deception)
//
// All test problems inherit from MultiObjectiveProblem<N, M> (defined in OptimizationProblem.h)
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_MULTI_OBJECTIVE_TEST_PROBLEMS_H
#define MML_MULTI_OBJECTIVE_TEST_PROBLEMS_H

#include "MMLBase.h"
#include "OptimizationProblem.h"

#include <cmath>

namespace MML
{
namespace TestProblems
{
    using namespace Optimization;  // For MultiObjectiveProblem, etc.


///////////////////////////////////////////////////////////////////////////////////////////
// ZDT1 Problem
//
// Convex Pareto front, easy problem
//
// f1(x) = x1
// f2(x) = g(x) * h(f1, g)
// g(x) = 1 + 9/(n-1) * sum(x2..xn)
// h(f1, g) = 1 - sqrt(f1/g)
//
// Bounds: xi ∈ [0, 1] for all i
// Pareto front: x1 ∈ [0,1], x2..xn = 0
// True front: f2 = 1 - sqrt(f1)
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 30>
class ZDT1 : public MultiObjectiveProblem
{
    static_assert(N >= 2, "ZDT1 requires at least 2 variables");

public:
    ZDT1() : MultiObjectiveProblem(N, 2) {}

    Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override
    {
        Real f1 = x[0];
        
        Real g = 0;
        for (int i = 1; i < N; ++i) {
            g += x[i];
        }
        g = 1 + 9 * g / (N - 1);
        
        Real h = 1 - std::sqrt(f1 / g);
        Real f2 = g * h;
        
        return Vector<Real>{f1, f2};
    }

    std::pair<Real, Real> GetVariableBounds(int /*i*/) const override
    {
        return {0.0, 1.0};
    }

    std::string GetProblemName() const override { return "ZDT1"; }

    bool IsOnParetoFront(const Vector<Real>& x, Real tolerance = 1e-6) const override
    {
        // On Pareto front: x2..xn should be 0
        for (int i = 1; i < N; ++i) {
            if (std::abs(x[i]) > tolerance)
                return false;
        }
        return x[0] >= 0 && x[0] <= 1;
    }

    /// @brief Get a point on the true Pareto front for f1 value
    static Vector<Real> GetParetoFrontPoint(Real f1)
    {
        return Vector<Real>{f1, 1 - std::sqrt(f1)};
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// ZDT2 Problem
//
// Non-convex (concave) Pareto front
//
// f1(x) = x1
// f2(x) = g(x) * h(f1, g)
// g(x) = 1 + 9/(n-1) * sum(x2..xn)
// h(f1, g) = 1 - (f1/g)^2
//
// Bounds: xi ∈ [0, 1] for all i
// Pareto front: x1 ∈ [0,1], x2..xn = 0
// True front: f2 = 1 - f1^2
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 30>
class ZDT2 : public MultiObjectiveProblem
{
    static_assert(N >= 2, "ZDT2 requires at least 2 variables");

public:
    ZDT2() : MultiObjectiveProblem(N, 2) {}

    Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override
    {
        Real f1 = x[0];
        
        Real g = 0;
        for (int i = 1; i < N; ++i) {
            g += x[i];
        }
        g = 1 + 9 * g / (N - 1);
        
        Real h = 1 - (f1 / g) * (f1 / g);
        Real f2 = g * h;
        
        return Vector<Real>{f1, f2};
    }

    std::pair<Real, Real> GetVariableBounds(int /*i*/) const override
    {
        return {0.0, 1.0};
    }

    std::string GetProblemName() const override { return "ZDT2"; }

    bool IsOnParetoFront(const Vector<Real>& x, Real tolerance = 1e-6) const override
    {
        for (int i = 1; i < N; ++i) {
            if (std::abs(x[i]) > tolerance)
                return false;
        }
        return x[0] >= 0 && x[0] <= 1;
    }

    static Vector<Real> GetParetoFrontPoint(Real f1)
    {
        return Vector<Real>{f1, 1 - f1 * f1};
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// ZDT3 Problem
//
// Discontinuous Pareto front (multiple disconnected segments)
//
// f1(x) = x1
// f2(x) = g(x) * h(f1, g)
// g(x) = 1 + 9/(n-1) * sum(x2..xn)
// h(f1, g) = 1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1)
//
// Bounds: xi ∈ [0, 1] for all i
// Pareto front: Disconnected segments
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 30>
class ZDT3 : public MultiObjectiveProblem
{
    static_assert(N >= 2, "ZDT3 requires at least 2 variables");

public:
    ZDT3() : MultiObjectiveProblem(N, 2) {}

    Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override
    {
        Real f1 = x[0];
        
        Real g = 0;
        for (int i = 1; i < N; ++i) {
            g += x[i];
        }
        g = 1 + 9 * g / (N - 1);
        
        Real ratio = f1 / g;
        Real h = 1 - std::sqrt(ratio) - ratio * std::sin(10 * Constants::PI * f1);
        Real f2 = g * h;
        
        return Vector<Real>{f1, f2};
    }

    std::pair<Real, Real> GetVariableBounds(int /*i*/) const override
    {
        return {0.0, 1.0};
    }

    std::string GetProblemName() const override { return "ZDT3"; }

    bool IsOnParetoFront(const Vector<Real>& x, Real tolerance = 1e-6) const override
    {
        for (int i = 1; i < N; ++i) {
            if (std::abs(x[i]) > tolerance)
                return false;
        }
        return x[0] >= 0 && x[0] <= 1;
    }

    static Vector<Real> GetParetoFrontPoint(Real f1)
    {
        return Vector<Real>{f1, 1 - std::sqrt(f1) - f1 * std::sin(10 * Constants::PI * f1)};
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// ZDT4 Problem
//
// Multimodal - many local Pareto fronts, one global
//
// f1(x) = x1
// f2(x) = g(x) * h(f1, g)
// g(x) = 1 + 10(n-1) + sum(xi^2 - 10*cos(4*pi*xi)) for i=2..n
// h(f1, g) = 1 - sqrt(f1/g)
//
// Bounds: x1 ∈ [0, 1], xi ∈ [-5, 5] for i=2..n
// Global Pareto front: x1 ∈ [0,1], x2..xn = 0
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 10>
class ZDT4 : public MultiObjectiveProblem
{
    static_assert(N >= 2, "ZDT4 requires at least 2 variables");

public:
    ZDT4() : MultiObjectiveProblem(N, 2) {}

    Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override
    {
        Real f1 = x[0];
        
        Real g = 1 + 10 * (N - 1);
        for (int i = 1; i < N; ++i) {
            g += x[i] * x[i] - 10 * std::cos(4 * Constants::PI * x[i]);
        }
        
        Real h = 1 - std::sqrt(f1 / g);
        Real f2 = g * h;
        
        return Vector<Real>{f1, f2};
    }

    std::pair<Real, Real> GetVariableBounds(int i) const override
    {
        if (i == 0)
            return {0.0, 1.0};
        return {-5.0, 5.0};
    }

    std::string GetProblemName() const override { return "ZDT4"; }

    bool IsOnParetoFront(const Vector<Real>& x, Real tolerance = 1e-6) const override
    {
        for (int i = 1; i < N; ++i) {
            if (std::abs(x[i]) > tolerance)
                return false;
        }
        return x[0] >= 0 && x[0] <= 1;
    }

    static Vector<Real> GetParetoFrontPoint(Real f1)
    {
        return Vector<Real>{f1, 1 - std::sqrt(f1)};
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// ZDT6 Problem
//
// Non-uniform distribution along Pareto front
//
// f1(x) = 1 - exp(-4*x1) * sin^6(6*pi*x1)
// f2(x) = g(x) * h(f1, g)
// g(x) = 1 + 9 * (sum(x2..xn)/(n-1))^0.25
// h(f1, g) = 1 - (f1/g)^2
//
// Bounds: xi ∈ [0, 1] for all i
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 10>
class ZDT6 : public MultiObjectiveProblem
{
    static_assert(N >= 2, "ZDT6 requires at least 2 variables");

public:
    ZDT6() : MultiObjectiveProblem(N, 2) {}

    Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override
    {
        Real sinVal = std::sin(6 * Constants::PI * x[0]);
        Real f1 = 1 - std::exp(-4 * x[0]) * std::pow(sinVal, 6);
        
        Real sumX = 0;
        for (int i = 1; i < N; ++i) {
            sumX += x[i];
        }
        Real g = 1 + 9 * std::pow(sumX / (N - 1), 0.25);
        
        Real h = 1 - (f1 / g) * (f1 / g);
        Real f2 = g * h;
        
        return Vector<Real>{f1, f2};
    }

    std::pair<Real, Real> GetVariableBounds(int /*i*/) const override
    {
        return {0.0, 1.0};
    }

    std::string GetProblemName() const override { return "ZDT6"; }

    bool IsOnParetoFront(const Vector<Real>& x, Real tolerance = 1e-6) const override
    {
        for (int i = 1; i < N; ++i) {
            if (std::abs(x[i]) > tolerance)
                return false;
        }
        return x[0] >= 0 && x[0] <= 1;
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Schaffer N.1 (SCH1) - Simple 1D bi-objective problem
//
// f1(x) = x^2
// f2(x) = (x-2)^2
//
// Bounds: x ∈ [-10^5, 10^5]
// Pareto front: x ∈ [0, 2]
///////////////////////////////////////////////////////////////////////////////////////////

class SCH1 : public MultiObjectiveProblem
{
public:
    SCH1() : MultiObjectiveProblem(1, 2) {}

    Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override
    {
        Real f1 = x[0] * x[0];
        Real f2 = (x[0] - 2) * (x[0] - 2);
        return Vector<Real>{f1, f2};
    }

    std::pair<Real, Real> GetVariableBounds(int /*i*/) const override
    {
        return {-1e5, 1e5};
    }

    std::string GetProblemName() const override { return "SCH1"; }

    bool IsOnParetoFront(const Vector<Real>& x, Real /*tolerance*/ = 1e-6) const override
    {
        return x[0] >= 0 && x[0] <= 2;
    }

    static Vector<Real> GetParetoFrontPoint(Real x)
    {
        return Vector<Real>{x * x, (x - 2) * (x - 2)};
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Fonseca-Fleming (FON) - Concave Pareto front
//
// f1(x) = 1 - exp(-sum((xi - 1/sqrt(n))^2))
// f2(x) = 1 - exp(-sum((xi + 1/sqrt(n))^2))
//
// Bounds: xi ∈ [-4, 4] for all i
// Pareto front: xi = xi for all i, where x ∈ [-1/sqrt(n), 1/sqrt(n)]
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 3>
class FON : public MultiObjectiveProblem
{
public:
    FON() : MultiObjectiveProblem(N, 2) {}

    Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override
    {
        Real invSqrtN = 1.0 / std::sqrt(static_cast<Real>(N));
        
        Real sum1 = 0, sum2 = 0;
        for (int i = 0; i < N; ++i) {
            Real t1 = x[i] - invSqrtN;
            Real t2 = x[i] + invSqrtN;
            sum1 += t1 * t1;
            sum2 += t2 * t2;
        }
        
        Real f1 = 1 - std::exp(-sum1);
        Real f2 = 1 - std::exp(-sum2);
        
        return Vector<Real>{f1, f2};
    }

    std::pair<Real, Real> GetVariableBounds(int /*i*/) const override
    {
        return {-4.0, 4.0};
    }

    std::string GetProblemName() const override { return "FON"; }

    bool IsOnParetoFront(const Vector<Real>& x, Real tolerance = 1e-6) const override
    {
        // On Pareto front: all xi should be equal
        Real invSqrtN = 1.0 / std::sqrt(static_cast<Real>(N));
        for (int i = 1; i < N; ++i) {
            if (std::abs(x[i] - x[0]) > tolerance)
                return false;
        }
        return x[0] >= -invSqrtN && x[0] <= invSqrtN;
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Kursawe (KUR) - Non-convex, discontinuous Pareto front
//
// f1(x) = sum(-10 * exp(-0.2 * sqrt(xi^2 + x(i+1)^2))) for i=1..n-1
// f2(x) = sum(|xi|^0.8 + 5*sin(xi^3))
//
// Bounds: xi ∈ [-5, 5] for all i
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 3>
class KUR : public MultiObjectiveProblem
{
    static_assert(N >= 2, "KUR requires at least 2 variables");

public:
    KUR() : MultiObjectiveProblem(N, 2) {}

    Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override
    {
        Real f1 = 0, f2 = 0;
        
        for (int i = 0; i < N - 1; ++i) {
            f1 += -10 * std::exp(-0.2 * std::sqrt(x[i]*x[i] + x[i+1]*x[i+1]));
        }
        
        for (int i = 0; i < N; ++i) {
            f2 += std::pow(std::abs(x[i]), 0.8) + 5 * std::sin(x[i]*x[i]*x[i]);
        }
        
        return Vector<Real>{f1, f2};
    }

    std::pair<Real, Real> GetVariableBounds(int /*i*/) const override
    {
        return {-5.0, 5.0};
    }

    std::string GetProblemName() const override { return "KUR"; }

    bool IsOnParetoFront(const Vector<Real>& /*x*/, Real /*tolerance*/ = 1e-6) const override
    {
        // KUR doesn't have a simple analytical Pareto front
        return false;
    }
};

} // namespace TestProblems
} // namespace MML

#endif // MML_MULTI_OBJECTIVE_TEST_PROBLEMS_H
