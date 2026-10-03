///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        PenaltyMethod.h                                                     ///
///  Description: Penalty method for constrained-to-unconstrained conversion +        ///
///               factory functions (split from OptimizationProblem.h)               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_OPTIMIZATION_PENALTYMETHOD_H
#define MML_OPTIMIZATION_PENALTYMETHOD_H

#include "MMLBase.h"
#include "base/Vector/VectorN.h"
#include "interfaces/IFunction.h"
#include "mml/mml_export.h"

#include <memory>
#include <vector>
#include <array>
#include <string>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <random>
#include <functional>

#include "ProblemTypes.h"
#include "IConstraint.h"

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////
	///                        PENALTY METHOD                               ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Penalty types for constraint handling
	 */
	enum class PenaltyType {
		Quadratic,      ///< P = ρ * violation²  (smooth, differentiable)
		Linear,         ///< P = ρ * violation   (exact for inequality, non-smooth)
		Exponential     ///< P = ρ * (exp(violation) - 1)  (strongly penalizes large violations)
	};

	/**
	 * @brief Converts constrained problems to unconstrained using penalty methods
	 * 
	 * The penalty method transforms:
	 *   minimize f(x) subject to g_i(x) ≤ 0, h_j(x) = 0
	 * into:
	 *   minimize f(x) + ρ * Σ P(g_i(x)) + ρ * Σ P(h_j(x))
	 * 
	 * where P is the penalty function and ρ is the penalty parameter.
	 * 
	 * Usage with existing optimizers:
	 * @code
	 * // Create constrained problem
	 * SingleObjectiveProblem base(2, objective);
	 * ConstrainedSingleObjectiveProblem constrained(base);
	 * constrained.AddInequalityConstraint(g1);
	 * constrained.AddEqualityConstraint(h1);
	 * 
	 * // Wrap with penalty method
	 * PenaltyMethod penalized(constrained, 1000.0);
	 * 
	 * // Use with any unconstrained optimizer
	 * auto result = NelderMead::Minimize(penalized, x0);
	 * @endcode
	 * 
	 * @note For best results, start with small ρ and increase progressively.
	 *       Very large ρ can cause numerical issues (ill-conditioning).
	 */
	class PenaltyMethod : public SingleObjectiveProblem {
	private:
		std::shared_ptr<ConstrainedSingleObjectiveProblem> _ownedProblem;
		const ConstrainedSingleObjectiveProblem* _constrainedProblem;
		Real _penaltyParameter;
		PenaltyType _penaltyType;
		Real _equalityWeight;  // Additional weight for equality constraints

		static const ConstrainedSingleObjectiveProblem& RequireProblem(
			const std::shared_ptr<ConstrainedSingleObjectiveProblem>& problem) {
			if (!problem)
				throw OptimizationProblemError("PenaltyMethod requires a constrained problem");
			return *problem;
		}

		/// @brief Compute penalty for a single constraint violation
		Real ComputePenalty(Real violation) const {
			if (violation <= 0) return 0;  // Feasible, no penalty

			switch (_penaltyType) {
				case PenaltyType::Quadratic:
					return violation * violation;
				case PenaltyType::Linear:
					return violation;
				case PenaltyType::Exponential:
					return std::exp(violation) - 1.0;
				default:
					return violation * violation;
			}
		}

	public:
		//---------------------------------------------------------------------
		// Constructors
		//---------------------------------------------------------------------

		/// @brief Construct with non-owning reference to constrained problem
		/// @param problem The constrained problem to convert
		/// @param rho Penalty parameter (larger = stricter constraint enforcement)
		/// @param type Penalty function type (default: Quadratic)
		explicit PenaltyMethod(
			const ConstrainedSingleObjectiveProblem& problem,
			Real rho = 1000.0,
			PenaltyType type = PenaltyType::Quadratic)
			: _constrainedProblem(&problem)
			, _penaltyParameter(rho)
			, _penaltyType(type)
			, _equalityWeight(1.0)
		{
			if (problem.GetNumObjectives() != 1)
				throw OptimizationProblemError("PenaltyMethod requires exactly 1 objective");
			this->GetProblemSpecMutable() = problem.GetProblemSpec();
			SetupPenalizedObjective();
		}

		/// @brief Construct with owning pointer to constrained problem
		explicit PenaltyMethod(
			std::shared_ptr<ConstrainedSingleObjectiveProblem> problem,
			Real rho = 1000.0,
			PenaltyType type = PenaltyType::Quadratic)
			: _ownedProblem(std::move(problem))
			, _constrainedProblem(&RequireProblem(_ownedProblem))
			, _penaltyParameter(rho)
			, _penaltyType(type)
			, _equalityWeight(1.0)
		{
			this->GetProblemSpecMutable() = _constrainedProblem->GetProblemSpec();
			SetupPenalizedObjective();
		}

		//---------------------------------------------------------------------
		// Configuration
		//---------------------------------------------------------------------

		/// @brief Set the penalty parameter ρ
		void SetPenaltyParameter(Real rho) {
			if (rho <= 0)
				throw OptimizationProblemError("Penalty parameter must be positive");
			_penaltyParameter = rho;
		}

		/// @brief Get current penalty parameter
		Real GetPenaltyParameter() const { return _penaltyParameter; }

		/// @brief Set the penalty type
		void SetPenaltyType(PenaltyType type) { _penaltyType = type; }

		/// @brief Get current penalty type
		PenaltyType GetPenaltyType() const { return _penaltyType; }

		/// @brief Set additional weight for equality constraints
		/// @note Equality constraints are often harder to satisfy, so extra weight helps
		void SetEqualityWeight(Real weight) {
			if (weight <= 0)
				throw OptimizationProblemError("Equality weight must be positive");
			_equalityWeight = weight;
		}

		/// @brief Get equality constraint weight
		Real GetEqualityWeight() const { return _equalityWeight; }

		//---------------------------------------------------------------------
		// Problem access
		//---------------------------------------------------------------------

		/// @brief Get the underlying constrained problem
		const ConstrainedSingleObjectiveProblem& GetConstrainedProblem() const {
			return *_constrainedProblem;
		}

		/// @brief Get the base problem (unwrapped from constraints)
		const IOptimizationProblem& GetBaseProblem() const {
			return _constrainedProblem->GetBaseProblem();
		}

		//---------------------------------------------------------------------
		// Evaluation
		//---------------------------------------------------------------------

		/// @brief Evaluate penalized objective at x
		/// @return f(x) + ρ * Σ penalties
		Real EvaluatePenalized(const Vector<Real>& x) const {
			Real f = _constrainedProblem->Evaluate(x);

			// Add penalty for each constraint
			Real totalPenalty = 0;
			int nc = _constrainedProblem->GetNumConstraints();
			
			for (int i = 0; i < nc; ++i) {
				const auto& constraint = _constrainedProblem->GetConstraint(i);
				Real violation = constraint.GetViolation(x);
				Real penalty = ComputePenalty(violation);
				
				// Apply extra weight for equality constraints
				if (constraint.GetType() == ConstraintType::Equality) {
					penalty *= _equalityWeight;
				}
				
				totalPenalty += penalty;
			}

			return f + _penaltyParameter * totalPenalty;
		}

		//---------------------------------------------------------------------
		// Analysis utilities
		//---------------------------------------------------------------------

		/// @brief Get penalty contribution breakdown for analysis
		struct PenaltyBreakdown {
			Real objectiveValue;      ///< Original f(x)
			Real totalPenalty;        ///< Sum of all penalties (before ρ scaling)
			Real penalizedValue;      ///< f(x) + ρ * totalPenalty
			Real maxViolation;        ///< Largest constraint violation
			int numViolated;          ///< Number of violated constraints
			bool isFeasible;          ///< All constraints satisfied?
		};

		/// @brief Analyze penalty contributions at a point
		PenaltyBreakdown Analyze(const Vector<Real>& x, Real tolerance = 1e-8) const {
			PenaltyBreakdown result;
			
			// Get original objective
			result.objectiveValue = _constrainedProblem->Evaluate(x);
			
			// Analyze constraints
			result.totalPenalty = 0;
			result.maxViolation = 0;
			result.numViolated = 0;
			result.isFeasible = true;
			
			int nc = _constrainedProblem->GetNumConstraints();
			for (int i = 0; i < nc; ++i) {
				const auto& constraint = _constrainedProblem->GetConstraint(i);
				Real violation = constraint.GetViolation(x);
				
				if (violation > tolerance) {
					result.numViolated++;
					result.isFeasible = false;
					result.maxViolation = std::max(result.maxViolation, violation);
				}
				
				Real penalty = ComputePenalty(violation);
				if (constraint.GetType() == ConstraintType::Equality) {
					penalty *= _equalityWeight;
				}
				result.totalPenalty += penalty;
			}
			
			result.penalizedValue = result.objectiveValue + 
				_penaltyParameter * result.totalPenalty;
			
			return result;
		}

		/// @brief Check if current solution is feasible
		bool IsSolutionFeasible(const Vector<Real>& x, Real tolerance = 1e-8) const {
			return _constrainedProblem->IsFeasible(x, tolerance);
		}

	private:
		/// @brief Set up the objective function for this SingleObjectiveProblem
		void SetupPenalizedObjective() {
			// Capture 'this' to create the penalized objective
			this->SetObjective([this](const Vector<Real>& x) {
				return this->EvaluatePenalized(x);
			});
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                        FACTORY FUNCTIONS                            ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Create a simple continuous optimization problem
	 * 
	 * @param dim Number of decision variables
	 * @param objective Objective function
	 * @param lower Lower bound for all variables
	 * @param upper Upper bound for all variables
	 */
	inline SingleObjectiveProblem CreateContinuousProblem(
		int dim,
		std::function<Real(const Vector<Real>&)> objective,
		Real lower,
		Real upper) {
		SingleObjectiveProblem problem(dim, std::move(objective));
		problem.GetProblemSpecMutable().SetAllContinuous(lower, upper);
		return problem;
	}

	/**
	 * @brief Create a continuous optimization problem with variable-specific bounds
	 * 
	 * @param objective Objective function
	 * @param lower Lower bounds vector (its size defines the dimension)
	 * @param upper Upper bounds vector
	 */
	inline SingleObjectiveProblem CreateContinuousProblem(
		std::function<Real(const Vector<Real>&)> objective,
		const Vector<Real>& lower,
		const Vector<Real>& upper) {
		SingleObjectiveProblem problem(static_cast<int>(lower.size()), std::move(objective));
		problem.SetBounds(lower, upper);
		return problem;
	}

	/**
	 * @brief Create a penalized problem from a constrained problem
	 * 
	 * @param problem The constrained problem
	 * @param rho Penalty parameter
	 * @param type Penalty function type
	 */
	inline PenaltyMethod CreatePenalizedProblem(
		const ConstrainedSingleObjectiveProblem& problem,
		Real rho = 1000.0,
		PenaltyType type = PenaltyType::Quadratic) {
		return PenaltyMethod(problem, rho, type);
	}

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_PENALTYMETHOD_H
