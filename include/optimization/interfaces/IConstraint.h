///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IConstraint.h                                                       ///
///  Description: Constraint type enumeration and constraint interface                ///
///               (split from OptimizationProblem.h)                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_OPTIMIZATION_ICONSTRAINT_H
#define MML_OPTIMIZATION_ICONSTRAINT_H

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

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////
	///                    CONSTRAINT TYPE ENUMERATION                      ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Types of constraints in optimization problems
	 */
	enum class ConstraintType {
		Equality,     ///< g(x) = 0
		Inequality    ///< g(x) ≤ 0 (feasible when non-positive)
	};

	///////////////////////////////////////////////////////////////////////////
	///                      CONSTRAINT INTERFACE                           ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Interface for a single constraint function
	 * 
	 * Constraints are expressed in standard form:
	 * - Equality: g(x) = 0
	 * - Inequality: g(x) ≤ 0 (feasible when g(x) ≤ 0)
	 * 
	 * Dimension of the decision space is a runtime property of the problem.
	 */
	class IConstraint {
	public:
		virtual ~IConstraint() = default;

		/// @brief Get constraint type (equality or inequality)
		virtual ConstraintType GetType() const = 0;

		/// @brief Evaluate constraint function at x
		/// @return g(x) where feasibility is g(x) ≤ 0 for inequality, g(x) = 0 for equality
		virtual Real Evaluate(const Vector<Real>& x) const = 0;

		/// @brief Get constraint violation magnitude
		/// @return 0 if feasible, positive value indicating violation magnitude
		virtual Real GetViolation(const Vector<Real>& x) const {
			Real g = Evaluate(x);
			switch (GetType()) {
				case ConstraintType::Inequality:
					return (g > 0) ? g : 0;
				case ConstraintType::Equality:
					return std::abs(g);
			}
			return 0;
		}

		/// @brief Check if constraint is satisfied at x
		/// @param tolerance Tolerance for equality constraints
		virtual bool IsFeasible(const Vector<Real>& x, Real tolerance = 1e-8) const {
			return GetViolation(x) <= tolerance;
		}

		/// @brief Get optional name for this constraint
		virtual std::string GetName() const { return ""; }

		//---------------------------------------------------------------------
		// Optional gradient support (for gradient-based methods)
		//---------------------------------------------------------------------

		/// @brief Check if gradient is available
		virtual bool HasGradient() const { return false; }

		/// @brief Compute gradient of constraint function
		/// @note Only valid if HasGradient() returns true
		virtual Vector<Real> Gradient(const Vector<Real>& x) const {
			return Vector<Real>(x.size());
		}
	};

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_ICONSTRAINT_H
