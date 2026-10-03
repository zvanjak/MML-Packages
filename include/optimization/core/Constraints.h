///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Constraints.h                                                       ///
///  Description: Function-based and integrality constraint implementations           ///
///               (split from OptimizationProblem.h)                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_OPTIMIZATION_CONSTRAINTS_H
#define MML_OPTIMIZATION_CONSTRAINTS_H

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

#include "IConstraint.h"
#include "Variables.h"

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////
	///                  FUNCTION-BASED CONSTRAINT                          ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Constraint implementation using a function object
	 * 
	 * Wraps any callable (function, lambda, functor) as a constraint.
	 */
	class FunctionConstraint : public IConstraint {
	private:
		std::function<Real(const Vector<Real>&)> _function;
		ConstraintType _type;
		std::string _name;

	public:
		/// @brief Construct from function and constraint type
		FunctionConstraint(
			std::function<Real(const Vector<Real>&)> func,
			ConstraintType type,
			const std::string& name = "")
			: _function(std::move(func))
			, _type(type)
			, _name(name) {}

		ConstraintType GetType() const override { return _type; }
		
		Real Evaluate(const Vector<Real>& x) const override {
			return _function(x);
		}

		std::string GetName() const override { return _name; }
	};

	///////////////////////////////////////////////////////////////////////////
	///                   INTEGRALITY CONSTRAINT                            ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Constraint that checks if discrete variables have integer values
	 * 
	 * This is a meta-constraint used to enforce integrality in mixed-integer
	 * optimization when using continuous relaxation approaches.
	 * 
	 * The violation is the sum of distances from each discrete variable
	 * to its nearest integer value.
	 */
	class IntegralityConstraint : public IConstraint {
	private:
		const ProblemSpec* _spec;
		Real _tolerance;
		std::string _name;

	public:
		/// @brief Construct with problem spec
		explicit IntegralityConstraint(
			const ProblemSpec& spec,
			Real tolerance = 1e-9,
			const std::string& name = "integrality")
			: _spec(&spec)
			, _tolerance(tolerance)
			, _name(name) {}

		ConstraintType GetType() const override {
			return ConstraintType::Equality;  // Integrality is an equality constraint (violation = 0)
		}

		/// @brief Evaluate integrality violation
		/// @return Sum of distances from discrete variables to nearest integers
		Real Evaluate(const Vector<Real>& x) const override {
			return _spec->GetIntegralityViolation(x);
		}

		/// @brief Get violation (same as Evaluate for this constraint)
		Real GetViolation(const Vector<Real>& x) const override {
			return Evaluate(x);
		}

		/// @brief Check if all discrete variables are integer (within tolerance)
		bool IsFeasible(const Vector<Real>& x, Real /*tolerance*/) const override {
			return Evaluate(x) <= _tolerance;
		}

		std::string GetName() const override { return _name; }
	};

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_CONSTRAINTS_H
