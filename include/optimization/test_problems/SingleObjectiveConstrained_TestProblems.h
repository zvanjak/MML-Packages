///////////////////////////////////////////////////////////////////////////////////////////
// OptimizationTestProblems.h
//
// Standard benchmark problems for testing and validating optimization algorithms.
//
// G-Series: Classic constrained optimization benchmarks from the CEC competitions.
// These problems are widely used for evaluating constrained optimization algorithms.
//
// References:
// - Michalewicz, Z., & Schoenauer, M. (1996). Evolutionary algorithms for constrained
//   parameter optimization problems. Evolutionary computation, 4(1), 1-32.
// - Koziel, S., & Michalewicz, Z. (1999). Evolutionary algorithms, homomorphous mappings,
//   and constrained parameter optimization. Evolutionary computation, 7(1), 19-44.
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_OPTIMIZATION_TEST_PROBLEMS_H
#define MML_OPTIMIZATION_TEST_PROBLEMS_H

#include "MMLBase.h"
#include "OptimizationProblem.h"

namespace MML
{
namespace TestProblems
{
    using namespace Optimization;

///////////////////////////////////////////////////////////////////////////////////////////
// G1 Problem
//
// Minimize: f(x) = 5*sum(x1..x4) - 5*sum(x1..x4)^2 - sum(x5..x13)
// Subject to:
//   g1: 2*x1 + 2*x2 + x10 + x11 <= 10
//   g2: 2*x1 + 2*x3 + x10 + x12 <= 10
//   g3: 2*x2 + 2*x3 + x11 + x12 <= 10
//   g4: -8*x1 + x10 <= 0
//   g5: -8*x2 + x11 <= 0
//   g6: -8*x3 + x12 <= 0
//   g7: -2*x4 - x5 + x10 <= 0
//   g8: -2*x6 - x7 + x11 <= 0
//   g9: -2*x8 - x9 + x12 <= 0
//
// Bounds: 0 <= xi <= 1 (i=1..9), 0 <= xi <= 100 (i=10..12), 0 <= x13 <= 1
// Global minimum: f* = -15 at x* = (1,1,1,1,1,1,1,1,1,3,3,3,1)
///////////////////////////////////////////////////////////////////////////////////////////


template<int N = 13>
class G1Problem : public ConstrainedSingleObjectiveProblem
{
    static_assert(N == 13, "G1 problem requires exactly 13 variables");

private:
    static ProblemSpec CreateSpec()
    {
        ProblemSpec spec(N);
        for (int i = 0; i < 9; ++i)
            spec.SetVariable(i, OptVariableSpec::Continuous(0, 1));
        for (int i = 9; i < 12; ++i)
            spec.SetVariable(i, OptVariableSpec::Continuous(0, 100));
        spec.SetVariable(12, OptVariableSpec::Continuous(0, 1));
        return spec;
    }

    static std::shared_ptr<SingleObjectiveProblem> CreateObjective()
    {
        return std::make_shared<SingleObjectiveProblem>(
            [](const Vector<Real>& x) {
                Real sum1 = 0, sum2 = 0, sum3 = 0;
                for (int i = 0; i < 4; ++i) {
                    sum1 += x[i];
                    sum2 += x[i] * x[i];
                }
                for (int i = 4; i < 13; ++i) {
                    sum3 += x[i];
                }
                return 5*sum1 - 5*sum2 - sum3;
            },
            CreateSpec()
        );
    }

public:
    G1Problem()
        : ConstrainedSingleObjectiveProblem(CreateObjective())
    {
        // Add inequality constraints (g <= 0)
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 2*x[0] + 2*x[1] + x[9] + x[10] - 10;  // g1
        });
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 2*x[0] + 2*x[2] + x[9] + x[11] - 10;  // g2
        });
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 2*x[1] + 2*x[2] + x[10] + x[11] - 10;  // g3
        });
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -8*x[0] + x[9];  // g4
        });
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -8*x[1] + x[10];  // g5
        });
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -8*x[2] + x[11];  // g6
        });
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -2*x[3] - x[4] + x[9];  // g7
        });
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -2*x[5] - x[6] + x[10];  // g8
        });
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -2*x[7] - x[8] + x[11];  // g9
        });
    }

    /// @brief Check if problem is constrained
    bool IsConstrained() const { return this->GetNumConstraints() > 0; }

    /// @brief Get problem dimension
    static constexpr int GetDimension() { return N; }

    /// @brief Get the known global optimum value
    static Real GetOptimalValue() { return -15.0; }

    /// @brief Get the known optimal solution
    static Vector<Real> GetOptimalSolution()
    {
        return Vector<Real>{1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1};
    }

    /// @brief Get problem name
    static std::string GetName() { return "G1"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// G6 Problem
//
// Minimize: f(x) = (x1 - 10)^3 + (x2 - 20)^3
// Subject to:
//   g1: -(x1 - 5)^2 - (x2 - 5)^2 + 100 <= 0
//   g2: (x1 - 6)^2 + (x2 - 5)^2 - 82.81 <= 0
//
// Bounds: 13 <= x1 <= 100, 0 <= x2 <= 100
// Global minimum: f* = -6961.81381 at x* = (14.095, 0.84296)
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 2>
class G6Problem : public ConstrainedSingleObjectiveProblem
{
    static_assert(N == 2, "G6 problem requires exactly 2 variables");

private:
    static ProblemSpec CreateSpec()
    {
        ProblemSpec spec(N);
        spec.SetVariable(0, OptVariableSpec::Continuous(13, 100));
        spec.SetVariable(1, OptVariableSpec::Continuous(0, 100));
        return spec;
    }

    static std::shared_ptr<SingleObjectiveProblem> CreateObjective()
    {
        return std::make_shared<SingleObjectiveProblem>(
            [](const Vector<Real>& x) {
                Real t1 = x[0] - 10;
                Real t2 = x[1] - 20;
                return t1*t1*t1 + t2*t2*t2;
            },
            CreateSpec()
        );
    }

public:
    G6Problem()
        : ConstrainedSingleObjectiveProblem(CreateObjective())
    {
        // g1: -(x1 - 5)^2 - (x2 - 5)^2 + 100 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -(x[0] - 5)*(x[0] - 5) - (x[1] - 5)*(x[1] - 5) + 100;
        });

        // g2: (x1 - 6)^2 + (x2 - 5)^2 - 82.81 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return (x[0] - 6)*(x[0] - 6) + (x[1] - 5)*(x[1] - 5) - 82.81;
        });
    }

    bool IsConstrained() const { return this->GetNumConstraints() > 0; }
    static constexpr int GetDimension() { return N; }
    static Real GetOptimalValue() { return -6961.81381; }
    static Vector<Real> GetOptimalSolution() { return Vector<Real>{14.095, 0.84296}; }
    static std::string GetName() { return "G6"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// G7 Problem
//
// Minimize: f(x) = x1^2 + x2^2 + x1*x2 - 14*x1 - 16*x2 + (x3-10)^2 + 4*(x4-5)^2
//                  + (x5-3)^2 + 2*(x6-1)^2 + 5*x7^2 + 7*(x8-11)^2 + 2*(x9-10)^2
//                  + (x10-7)^2 + 45
// Subject to: 8 inequality constraints
//
// Bounds: -10 <= xi <= 10 for all i
// Global minimum: f* = 24.3062 at a specific point
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 10>
class G7Problem : public ConstrainedSingleObjectiveProblem
{
    static_assert(N == 10, "G7 problem requires exactly 10 variables");

private:
    static ProblemSpec CreateSpec()
    {
        ProblemSpec spec(N);
        spec.SetAllContinuous(-10, 10);
        return spec;
    }

    static std::shared_ptr<SingleObjectiveProblem> CreateObjective()
    {
        return std::make_shared<SingleObjectiveProblem>(
            [](const Vector<Real>& x) {
                return x[0]*x[0] + x[1]*x[1] + x[0]*x[1] - 14*x[0] - 16*x[1]
                     + (x[2]-10)*(x[2]-10) + 4*(x[3]-5)*(x[3]-5)
                     + (x[4]-3)*(x[4]-3) + 2*(x[5]-1)*(x[5]-1)
                     + 5*x[6]*x[6] + 7*(x[7]-11)*(x[7]-11)
                     + 2*(x[8]-10)*(x[8]-10) + (x[9]-7)*(x[9]-7) + 45;
            },
            CreateSpec()
        );
    }

public:
    G7Problem()
        : ConstrainedSingleObjectiveProblem(CreateObjective())
    {
        // g1: -105 + 4*x1 + 5*x2 - 3*x7 + 9*x8 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -105 + 4*x[0] + 5*x[1] - 3*x[6] + 9*x[7];
        });
        // g2: 10*x1 - 8*x2 - 17*x7 + 2*x8 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 10*x[0] - 8*x[1] - 17*x[6] + 2*x[7];
        });
        // g3: -8*x1 + 2*x2 + 5*x9 - 2*x10 - 12 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -8*x[0] + 2*x[1] + 5*x[8] - 2*x[9] - 12;
        });
        // g4: 3*(x1-2)^2 + 4*(x2-3)^2 + 2*x3^2 - 7*x4 - 120 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 3*(x[0]-2)*(x[0]-2) + 4*(x[1]-3)*(x[1]-3) + 2*x[2]*x[2] - 7*x[3] - 120;
        });
        // g5: 5*x1^2 + 8*x2 + (x3-6)^2 - 2*x4 - 40 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 5*x[0]*x[0] + 8*x[1] + (x[2]-6)*(x[2]-6) - 2*x[3] - 40;
        });
        // g6: x1^2 + 2*(x2-2)^2 - 2*x1*x2 + 14*x5 - 6*x6 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return x[0]*x[0] + 2*(x[1]-2)*(x[1]-2) - 2*x[0]*x[1] + 14*x[4] - 6*x[5];
        });
        // g7: 0.5*(x1-8)^2 + 2*(x2-4)^2 + 3*x5^2 - x6 - 30 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 0.5*(x[0]-8)*(x[0]-8) + 2*(x[1]-4)*(x[1]-4) + 3*x[4]*x[4] - x[5] - 30;
        });
        // g8: -3*x1 + 6*x2 + 12*(x9-8)^2 - 7*x10 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -3*x[0] + 6*x[1] + 12*(x[8]-8)*(x[8]-8) - 7*x[9];
        });
    }

    bool IsConstrained() const { return this->GetNumConstraints() > 0; }
    static constexpr int GetDimension() { return N; }
    static Real GetOptimalValue() { return 24.3062; }
    static Vector<Real> GetOptimalSolution() {
        return Vector<Real>{2.171996, 2.363683, 8.773926, 5.095984, 0.9906548,
                                1.430574, 1.321644, 9.828726, 8.280092, 8.375927};
    }
    static std::string GetName() { return "G7"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// G9 Problem
//
// Minimize: f(x) = (x1-10)^2 + 5*(x2-12)^2 + x3^4 + 3*(x4-11)^2 + 10*x5^6
//                  + 7*x6^2 + x7^4 - 4*x6*x7 - 10*x6 - 8*x7
// Subject to: 4 inequality constraints
//
// Bounds: -10 <= xi <= 10 for all i
// Global minimum: f* = 680.6300573 at a specific point
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 7>
class G9Problem : public ConstrainedSingleObjectiveProblem
{
    static_assert(N == 7, "G9 problem requires exactly 7 variables");

private:
    static ProblemSpec CreateSpec()
    {
        ProblemSpec spec(N);
        spec.SetAllContinuous(-10, 10);
        return spec;
    }

    static std::shared_ptr<SingleObjectiveProblem> CreateObjective()
    {
        return std::make_shared<SingleObjectiveProblem>(
            [](const Vector<Real>& x) {
                Real x6_2 = x[5]*x[5];
                return (x[0]-10)*(x[0]-10) + 5*(x[1]-12)*(x[1]-12)
                     + x[2]*x[2]*x[2]*x[2] + 3*(x[3]-11)*(x[3]-11)
                     + 10*x[4]*x[4]*x[4]*x[4]*x[4]*x[4]
                     + 7*x6_2 + x[6]*x[6]*x[6]*x[6]
                     - 4*x[5]*x[6] - 10*x[5] - 8*x[6];
            },
            CreateSpec()
        );
    }

public:
    G9Problem()
        : ConstrainedSingleObjectiveProblem(CreateObjective())
    {
        // g1: -127 + 2*x1^2 + 3*x2^4 + x3 + 4*x4^2 + 5*x5 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -127 + 2*x[0]*x[0] + 3*x[1]*x[1]*x[1]*x[1] + x[2] + 4*x[3]*x[3] + 5*x[4];
        });
        // g2: -282 + 7*x1 + 3*x2 + 10*x3^2 + x4 - x5 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -282 + 7*x[0] + 3*x[1] + 10*x[2]*x[2] + x[3] - x[4];
        });
        // g3: -196 + 23*x1 + x2^2 + 6*x6^2 - 8*x7 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -196 + 23*x[0] + x[1]*x[1] + 6*x[5]*x[5] - 8*x[6];
        });
        // g4: 4*x1^2 + x2^2 - 3*x1*x2 + 2*x3^2 + 5*x6 - 11*x7 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 4*x[0]*x[0] + x[1]*x[1] - 3*x[0]*x[1] + 2*x[2]*x[2] + 5*x[5] - 11*x[6];
        });
    }

    bool IsConstrained() const { return this->GetNumConstraints() > 0; }
    static constexpr int GetDimension() { return N; }
    static Real GetOptimalValue() { return 680.6300573; }
    static Vector<Real> GetOptimalSolution() {
        return Vector<Real>{2.330499, 1.951372, -0.4775414, 4.365726, -0.6244870,
                                1.038131, 1.594227};
    }
    static std::string GetName() { return "G9"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// G10 Problem (Pressure Vessel Design)
//
// Minimize: f(x) = 0.6224*x1*x3*x4 + 1.7781*x2*x3^2 + 3.1661*x1^2*x4 + 19.84*x1^2*x3
// Subject to: 4 inequality constraints
//
// Bounds: 0 <= x1, x2 <= 99, 10 <= x3, x4 <= 200
// Global minimum: f* = 7049.248 at x* approximately (0.8125, 0.4375, 42.0984, 176.6366)
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 4>
class G10Problem : public ConstrainedSingleObjectiveProblem
{
    static_assert(N == 4, "G10 problem requires exactly 4 variables");

private:
    static ProblemSpec CreateSpec()
    {
        ProblemSpec spec(N);
        spec.SetVariable(0, OptVariableSpec::Continuous(0, 99));
        spec.SetVariable(1, OptVariableSpec::Continuous(0, 99));
        spec.SetVariable(2, OptVariableSpec::Continuous(10, 200));
        spec.SetVariable(3, OptVariableSpec::Continuous(10, 200));
        return spec;
    }

    static std::shared_ptr<SingleObjectiveProblem> CreateObjective()
    {
        return std::make_shared<SingleObjectiveProblem>(
            [](const Vector<Real>& x) {
                Real x1_2 = x[0]*x[0];
                Real x3_2 = x[2]*x[2];
                return 0.6224*x[0]*x[2]*x[3] + 1.7781*x[1]*x3_2
                     + 3.1661*x1_2*x[3] + 19.84*x1_2*x[2];
            },
            CreateSpec()
        );
    }

public:
    G10Problem()
        : ConstrainedSingleObjectiveProblem(CreateObjective())
    {
        // g1: -x1 + 0.0193*x3 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -x[0] + 0.0193*x[2];
        });
        // g2: -x2 + 0.00954*x3 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return -x[1] + 0.00954*x[2];
        });
        // g3: -pi*x3^2*x4 - (4/3)*pi*x3^3 + 1296000 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            Real x3_2 = x[2]*x[2];
            Real x3_3 = x3_2*x[2];
            return -Constants::PI*x3_2*x[3] - (4.0/3.0)*Constants::PI*x3_3 + 1296000;
        });
        // g4: x4 - 240 <= 0
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return x[3] - 240;
        });
    }

    bool IsConstrained() const { return this->GetNumConstraints() > 0; }
    static constexpr int GetDimension() { return N; }
    static Real GetOptimalValue() { return 7049.248; }
    
    /// @brief Get the known optimal solution (rounded values from literature)
    /// @note The standard literature values are slightly infeasible due to rounding.
    ///       Use GetFeasibleNearOptimalSolution() for a strictly feasible point.
    static Vector<Real> GetOptimalSolution() { return Vector<Real>{0.8125, 0.4375, 42.0984, 176.6366}; }
    
    /// @brief Get a feasible solution that is near-optimal
    /// @note Adjusted from literature values to satisfy g1 constraint strictly
    static Vector<Real> GetFeasibleNearOptimalSolution() { 
        // g1 constraint: -x1 + 0.0193*x3 <= 0  =>  x1 >= 0.0193*x3
        // With x3 = 42.0984, we need x1 >= 0.8125, so we use x1 = 0.8126
        return Vector<Real>{0.8126, 0.4375, 42.0984, 176.6366}; 
    }
    static std::string GetName() { return "G10"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Himmelblau's Problem (classic 2D test)
//
// Minimize: f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
// Subject to:
//   g1: 26 - (x-5)^2 - y^2 >= 0  (inside circle)
//
// Bounds: -5 <= x, y <= 5
// Has multiple local minima
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 2>
class HimmelblauConstrained : public ConstrainedSingleObjectiveProblem
{
    static_assert(N == 2, "Himmelblau problem requires exactly 2 variables");

private:
    static ProblemSpec CreateSpec()
    {
        ProblemSpec spec(N);
        spec.SetAllContinuous(-5, 5);
        return spec;
    }

    static std::shared_ptr<SingleObjectiveProblem> CreateObjective()
    {
        return std::make_shared<SingleObjectiveProblem>(
            [](const Vector<Real>& x) {
                Real t1 = x[0]*x[0] + x[1] - 11;
                Real t2 = x[0] + x[1]*x[1] - 7;
                return t1*t1 + t2*t2;
            },
            CreateSpec()
        );
    }

public:
    HimmelblauConstrained()
        : ConstrainedSingleObjectiveProblem(CreateObjective())
    {
        // g1: (x-5)^2 + y^2 - 26 <= 0  (inside circle centered at (5,0))
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return (x[0]-5)*(x[0]-5) + x[1]*x[1] - 26;
        });
    }

    bool IsConstrained() const { return this->GetNumConstraints() > 0; }
    static constexpr int GetDimension() { return N; }
    static std::string GetName() { return "Himmelblau (Constrained)"; }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Simple Quadratic with Linear Constraint (for basic testing)
//
// Minimize: f(x) = x1^2 + x2^2
// Subject to: x1 + x2 >= 1
//
// Bounds: -10 <= xi <= 10
// Global minimum: f* = 0.5 at x* = (0.5, 0.5)
///////////////////////////////////////////////////////////////////////////////////////////

template<int N = 2>
class SimpleConstrainedQuadratic : public ConstrainedSingleObjectiveProblem
{
    static_assert(N == 2, "SimpleConstrainedQuadratic requires exactly 2 variables");

private:
    static ProblemSpec CreateSpec()
    {
        ProblemSpec spec(N);
        spec.SetAllContinuous(-10, 10);
        return spec;
    }

    static std::shared_ptr<SingleObjectiveProblem> CreateObjective()
    {
        return std::make_shared<SingleObjectiveProblem>(
            [](const Vector<Real>& x) {
                return x[0]*x[0] + x[1]*x[1];
            },
            CreateSpec()
        );
    }

public:
    SimpleConstrainedQuadratic()
        : ConstrainedSingleObjectiveProblem(CreateObjective())
    {
        // g1: 1 - x1 - x2 <= 0  (equivalent to x1 + x2 >= 1)
        this->AddInequalityConstraint([](const Vector<Real>& x) {
            return 1 - x[0] - x[1];
        });
    }

    bool IsConstrained() const { return this->GetNumConstraints() > 0; }
    static constexpr int GetDimension() { return N; }
    static Real GetOptimalValue() { return 0.5; }
    static Vector<Real> GetOptimalSolution() { return Vector<Real>{0.5, 0.5}; }
    static std::string GetName() { return "Simple Constrained Quadratic"; }
};

} // namespace TestProblems
} // namespace MML

#endif // MML_OPTIMIZATION_TEST_PROBLEMS_H
