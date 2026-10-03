///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IOptimizationProblem.h                                              ///
///  Description: Abstract base interface for all optimization problems               ///
///               (split from OptimizationProblem.h)                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_OPTIMIZATION_IPROBLEM_H
#define MML_OPTIMIZATION_IPROBLEM_H

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

#include "Variables.h"
#include "IConstraint.h"

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////
	///                  IOPTIMIZATION PROBLEM INTERFACE                    ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Abstract base interface for all optimization problems
	 * 
	 * This is the unified interface that all optimization problem types
	 * (single-objective, multi-objective, constrained) derive from.
	 * 
	 * Provides:
	 * - Variable specification (bounds, types)
	 * - Constraint management
	 * - Feasibility checking
	 * 
	 * Problem dimension is a runtime property, obtained from the ProblemSpec.
	 */
	class IOptimizationProblem {
	public:
		virtual ~IOptimizationProblem() = default;

		//---------------------------------------------------------------------
		// Dimension and variable specification
		//---------------------------------------------------------------------

		/// @brief Get problem dimension (runtime, derived from the problem spec)
		virtual int Dimension() const { return GetProblemSpec().Dimension(); }

		/// @brief Get variable specification
		virtual const ProblemSpec& GetProblemSpec() const = 0;

		/// @brief Get lower bounds
		Vector<Real> GetLowerBounds() const { return GetProblemSpec().GetLowerBounds(); }

		/// @brief Get upper bounds
		Vector<Real> GetUpperBounds() const { return GetProblemSpec().GetUpperBounds(); }

		//---------------------------------------------------------------------
		// Constraint access
		//---------------------------------------------------------------------

		/// @brief Get number of constraints
		virtual int GetNumConstraints() const = 0;

		/// @brief Get constraint by index
		virtual const IConstraint& GetConstraint(int index) const = 0;

		/// @brief Check if problem has any constraints
		bool HasConstraints() const { return GetNumConstraints() > 0; }

		//---------------------------------------------------------------------
		// Feasibility checking
		//---------------------------------------------------------------------

		/// @brief Check if solution is feasible (satisfies all constraints and bounds)
		virtual bool IsFeasible(const Vector<Real>& x, Real tolerance = 1e-8) const {
			// Check bounds
			if (!GetProblemSpec().IsValid(x))
				return false;

			// Check constraints
			for (int i = 0; i < GetNumConstraints(); ++i) {
				if (!GetConstraint(i).IsFeasible(x, tolerance))
					return false;
			}
			return true;
		}

		/// @brief Get total constraint violation
		virtual Real GetTotalViolation(const Vector<Real>& x) const {
			Real total = 0;
			for (int i = 0; i < GetNumConstraints(); ++i) {
				total += GetConstraint(i).GetViolation(x);
			}
			return total;
		}

		/// @brief Get violation for each constraint
		virtual Vector<Real> GetViolationVector(const Vector<Real>& x) const {
			int nc = GetNumConstraints();
			Vector<Real> violations(nc);
			for (int i = 0; i < nc; ++i) {
				violations[i] = GetConstraint(i).GetViolation(x);
			}
			return violations;
		}

		//---------------------------------------------------------------------
		// Problem type queries
		//---------------------------------------------------------------------

		/// @brief Get number of objectives (1 for single-objective, M for multi-objective)
		virtual int GetNumObjectives() const = 0;

		/// @brief Check if this is a single-objective problem
		bool IsSingleObjective() const { return GetNumObjectives() == 1; }

		/// @brief Check if this is a multi-objective problem
		bool IsMultiObjective() const { return GetNumObjectives() > 1; }
	};

	class ISingleObjectiveProblem : public IOptimizationProblem {
	public:
		virtual Real Evaluate(const Vector<Real>& x) const = 0;
	};

	class IMultiObjectiveProblem : public IOptimizationProblem {
	public:
		virtual Vector<Real> EvaluateObjectives(const Vector<Real>& x) const = 0;

		Vector<Real> EvaluateValidated(const Vector<Real>& x) const {
			if (Dimension() <= 0 || x.size() != static_cast<size_t>(Dimension()))
				throw OptimizationProblemError("Objective input dimension mismatch");
			if (GetNumObjectives() < 2)
				throw OptimizationProblemError("Multi-objective problem requires at least 2 objectives");
			auto values = EvaluateObjectives(x);
			if (values.size() != static_cast<size_t>(GetNumObjectives()))
				throw OptimizationProblemError("Objective result count mismatch");
			return values;
		}
	};

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_IPROBLEM_H
