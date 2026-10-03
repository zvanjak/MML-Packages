///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Variables.h                                                         ///
///  Description: Optimization variable specification and decision vectors            ///
///               (split from OptimizationProblem.h)                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_OPTIMIZATION_VARIABLES_H
#define MML_OPTIMIZATION_VARIABLES_H

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
	///                    OPTIMIZATION PROBLEM ERROR                       ///
	///////////////////////////////////////////////////////////////////////////

	class OptimizationProblemError : public std::runtime_error {
	public:
		explicit OptimizationProblemError(const std::string& message)
			: std::runtime_error("OptimizationProblemError: " + message) {}
	};

	///////////////////////////////////////////////////////////////////////////
	///                         VARIABLE TYPES                              ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Types of variables supported in optimization problems
	 * 
	 * Supports standard mathematical programming variable types:
	 * - Continuous: Real-valued variables (standard NLP)
	 * - Integer: Integer-valued variables (MILP, MINLP)
	 * - Binary: 0/1 variables (feature selection, yes/no decisions)
	 * - Categorical: Unordered discrete choices (encoded as integers)
	 * 
	 * @note Permutation/assignment problems use specialized interfaces
	 *       (IPermutationProblem) since they require different algorithms.
	 */
	enum class OptVariableType {
		Continuous,   ///< Real values in [lower, upper]
		Integer,      ///< Integers in [lower, upper]
		Binary,       ///< {0, 1} - special case of Integer(0,1)
		Categorical   ///< One of N unordered categories (encoded as 0, 1, ..., K-1)
	};

	///////////////////////////////////////////////////////////////////////////
	///                        VARIABLE SPECIFICATION                       ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Specification for a single optimization variable
	 * 
	 * Defines the type, bounds, and properties for one decision variable.
	 * Used by ProblemSpec to describe the complete variable space.
	 */
	struct MML_OPTIMIZATION_API OptVariableSpec {
		std::string name;           ///< Optional name for debugging/reporting
		OptVariableType type;       ///< Variable type
		Real lowerBound;            ///< Minimum value (for Continuous/Integer)
		Real upperBound;            ///< Maximum value (for Continuous/Integer)
		int numCategories;          ///< For Categorical: number of options

		/// @brief Default constructor creates a continuous variable in [0, 1]
		OptVariableSpec();

		//---------------------------------------------------------------------
		// Factory methods for creating variable specifications
		//---------------------------------------------------------------------

		/// @brief Create a continuous variable specification
		static OptVariableSpec Continuous(Real lower, Real upper, const std::string& varName = "");

		/// @brief Create an integer variable specification
		static OptVariableSpec Integer(int lower, int upper, const std::string& varName = "");

		/// @brief Create a binary (0/1) variable specification
		static OptVariableSpec Binary(const std::string& varName = "");

		/// @brief Create a categorical variable specification
		static OptVariableSpec Categorical(int numOptions, const std::string& varName = "");

		//---------------------------------------------------------------------
		// Validation and repair
		//---------------------------------------------------------------------

		/// @brief Check if a value is within bounds for this variable
		bool IsInBounds(Real value) const;

		/// @brief Repair a value to be valid for this variable type
		Real Repair(Real value) const;

		/// @brief Check if this is a discrete (non-continuous) variable type
		bool IsDiscrete() const;
	};

	///////////////////////////////////////////////////////////////////////////
	///                         PROBLEM SPECIFICATION                       ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Complete specification of optimization problem variables
	 * 
	 * Defines all variables with their types, bounds, and properties.
	 * This is the unified problem specification used by all optimizers.
	 * The number of variables is a runtime property.
	 */
	class ProblemSpec {
	private:
		std::vector<OptVariableSpec> _variables;
		static int CheckedDimension(int dimension) {
			if (dimension < 0)
				throw OptimizationProblemError("ProblemSpec dimension must be non-negative");
			return dimension;
		}

	public:
		/// @brief Default constructor - empty (no variables)
		ProblemSpec() = default;

		/// @brief Construct with a given number of variables (all continuous in [0, 1])
		explicit ProblemSpec(int n) : _variables(CheckedDimension(n)) {}

		//---------------------------------------------------------------------
		// Dimension
		//---------------------------------------------------------------------

		/// @brief Get number of variables (runtime)
		int Dimension() const { return static_cast<int>(_variables.size()); }

		//---------------------------------------------------------------------
		// Variable access
		//---------------------------------------------------------------------

		/// @brief Access variable specification by index (const)
		const OptVariableSpec& operator[](int i) const { 
			if (i < 0 || i >= Dimension())
				throw OptimizationProblemError("Variable index out of range");
			return _variables[i]; 
		}

		/// @brief Access variable specification by index (non-const)
		OptVariableSpec& operator[](int i) { 
			if (i < 0 || i >= Dimension())
				throw OptimizationProblemError("Variable index out of range");
			return _variables[i]; 
		}

		/// @brief Get all variables
		const std::vector<OptVariableSpec>& Variables() const { return _variables; }

		//---------------------------------------------------------------------
		// Fluent API for building problem spec
		//---------------------------------------------------------------------

		/// @brief Append a variable specification
		ProblemSpec& AddVariable(const OptVariableSpec& spec) {
			_variables.push_back(spec);
			return *this;
		}

		/// @brief Create a homogeneous continuous problem
		static ProblemSpec Continuous(int dimension, Real lower, Real upper) {
			if (dimension < 0)
				throw OptimizationProblemError("ProblemSpec dimension must be non-negative");
			ProblemSpec spec;
			for (int i = 0; i < dimension; ++i)
				spec.AddVariable(OptVariableSpec::Continuous(lower, upper));
			return spec;
		}

		/// @brief Create a continuous problem from per-variable bounds
		static ProblemSpec Continuous(const Vector<Real>& lower, const Vector<Real>& upper) {
			if (lower.size() != upper.size())
				throw OptimizationProblemError("Bound vectors must have the same dimension");
			ProblemSpec spec;
			for (int i = 0; i < lower.size(); ++i)
				spec.AddVariable(OptVariableSpec::Continuous(lower[i], upper[i]));
			return spec;
		}

		/// @brief Set a single variable specification
		ProblemSpec& SetVariable(int index, const OptVariableSpec& spec) {
			if (index < 0 || index >= Dimension())
				throw OptimizationProblemError("Variable index out of range");
			_variables[index] = spec;
			return *this;
		}

		/// @brief Set all variables to continuous with same bounds
		ProblemSpec& SetAllContinuous(Real lower, Real upper) {
			for (int i = 0; i < Dimension(); ++i)
				_variables[i] = OptVariableSpec::Continuous(lower, upper);
			return *this;
		}

		/// @brief Set all variables to integer with same bounds
		ProblemSpec& SetAllInteger(int lower, int upper) {
			for (int i = 0; i < Dimension(); ++i)
				_variables[i] = OptVariableSpec::Integer(lower, upper);
			return *this;
		}

		/// @brief Set all variables to binary
		ProblemSpec& SetAllBinary() {
			for (int i = 0; i < Dimension(); ++i)
				_variables[i] = OptVariableSpec::Binary();
			return *this;
		}

		//---------------------------------------------------------------------
		// Query variable properties
		//---------------------------------------------------------------------

		/// @brief Get variable type for index i
		OptVariableType GetType(int i) const { return (*this)[i].type; }

		/// @brief Check if problem has any integer variables
		bool HasIntegerVariables() const {
			for (int i = 0; i < Dimension(); ++i)
				if (_variables[i].type == OptVariableType::Integer)
					return true;
			return false;
		}

		/// @brief Check if problem has any binary variables
		bool HasBinaryVariables() const {
			for (int i = 0; i < Dimension(); ++i)
				if (_variables[i].type == OptVariableType::Binary)
					return true;
			return false;
		}

		/// @brief Check if problem has any categorical variables
		bool HasCategoricalVariables() const {
			for (int i = 0; i < Dimension(); ++i)
				if (_variables[i].type == OptVariableType::Categorical)
					return true;
			return false;
		}

		/// @brief Check if problem is purely continuous (no discrete variables)
		bool IsPurelyContinuous() const {
			for (int i = 0; i < Dimension(); ++i)
				if (_variables[i].IsDiscrete())
					return false;
			return true;
		}

		/// @brief Check if problem has mixed integer and continuous variables
		bool IsMixedInteger() const {
			bool hasContinuous = false;
			bool hasDiscrete = false;
			for (int i = 0; i < Dimension(); ++i) {
				if (_variables[i].type == OptVariableType::Continuous)
					hasContinuous = true;
				else
					hasDiscrete = true;
			}
			return hasContinuous && hasDiscrete;
		}

		//---------------------------------------------------------------------
		// Bounds access
		//---------------------------------------------------------------------

		/// @brief Get lower bounds as a vector
		Vector<Real> GetLowerBounds() const {
			Vector<Real> lb(Dimension());
			for (int i = 0; i < Dimension(); ++i)
				lb[i] = _variables[i].lowerBound;
			return lb;
		}

		/// @brief Get upper bounds as a vector
		Vector<Real> GetUpperBounds() const {
			Vector<Real> ub(Dimension());
			for (int i = 0; i < Dimension(); ++i)
				ub[i] = _variables[i].upperBound;
			return ub;
		}

		//---------------------------------------------------------------------
		// Validation and repair
		//---------------------------------------------------------------------

		/// @brief Check if a solution is valid (within bounds and satisfies type constraints)
		bool IsValid(const Vector<Real>& x) const {
			if (static_cast<int>(x.size()) != Dimension())
				return false;
			for (int i = 0; i < Dimension(); ++i) {
				if (!_variables[i].IsInBounds(x[i]))
					return false;

				// Check discrete variables have integer values
				if (_variables[i].IsDiscrete()) {
					if (std::abs(x[i] - std::round(x[i])) > 1e-9)
						return false;
				}
			}
			return true;
		}

		/// @brief Repair a solution to be valid
		Vector<Real> Repair(const Vector<Real>& x) const {
			if (static_cast<int>(x.size()) != Dimension())
				throw OptimizationProblemError("Decision vector dimension mismatch");
			Vector<Real> repaired(Dimension());
			for (int i = 0; i < Dimension(); ++i) {
				repaired[i] = _variables[i].Repair(x[i]);
			}
			return repaired;
		}

		/// @brief Get integrality violation (sum of distances to nearest integers for discrete vars)
		Real GetIntegralityViolation(const Vector<Real>& x) const {
			if (static_cast<int>(x.size()) != Dimension())
				throw OptimizationProblemError("Decision vector dimension mismatch");
			Real violation = 0;
			for (int i = 0; i < Dimension(); ++i) {
				if (_variables[i].IsDiscrete()) {
					Real frac = x[i] - std::floor(x[i]);
					violation += std::min<Real>(frac, Real(1.0) - frac);
				}
			}
			return violation;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                    INTEGER REPAIR STRATEGIES                        ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Strategies for rounding continuous values to integers
	 * 
	 * When continuous optimizers find solutions for mixed-integer problems,
	 * we need strategies to repair floating-point values to valid integers.
	 */
	enum class IntegerRepairStrategy {
		Nearest,        ///< Round to nearest integer (std::round)
		Floor,          ///< Round down (std::floor)
		Ceil,           ///< Round up (std::ceil)
		Probabilistic   ///< Round probabilistically based on fractional part
	};

	///////////////////////////////////////////////////////////////////////////
	///                       DECISION VECTOR                               ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Type-safe wrapper for optimization decision vectors
	 * 
	 * Provides type-safe access to mixed-integer decision variables.
	 * Knows the variable types from ProblemSpec and offers appropriate
	 * accessors for continuous, integer, binary, and categorical values.
	 * 
	 * Usage:
	 * @code
	 * ProblemSpec spec(3);
	 * spec.SetVariable(0, OptVariableSpec::Continuous(-10, 10, "x"));
	 * spec.SetVariable(1, OptVariableSpec::Integer(0, 100, "quantity"));
	 * spec.SetVariable(2, OptVariableSpec::Binary("use_feature"));
	 * 
	 * DecisionVector dv(spec);
	 * dv.SetReal(0, 3.14);
	 * dv.SetInteger(1, 42);
	 * dv.SetBinary(2, true);
	 * 
	 * // Type-safe access
	 * Real x = dv.GetReal(0);
	 * int qty = dv.GetInteger(1);
	 * bool useFeature = dv.GetBinary(2);
	 * @endcode
	 */
	class DecisionVector {
	private:
		Vector<Real> _values;
		const ProblemSpec* _spec;

	public:
		//---------------------------------------------------------------------
		// Constructors
		//---------------------------------------------------------------------

		/// @brief Default constructor (no spec attached)
		DecisionVector() : _spec(nullptr) {}

		/// @brief Construct with problem specification (values sized to the spec)
		explicit DecisionVector(const ProblemSpec& spec)
			: _values(spec.Dimension()), _spec(&spec) {}

		/// @brief Construct from raw values and spec
		DecisionVector(const Vector<Real>& values, const ProblemSpec& spec)
			: _values(values), _spec(&spec) {}

		/// @brief Construct from raw values without spec
		explicit DecisionVector(const Vector<Real>& values)
			: _values(values), _spec(nullptr) {}

		//---------------------------------------------------------------------
		// Raw access
		//---------------------------------------------------------------------

		/// @brief Get underlying vector (const)
		const Vector<Real>& Values() const { return _values; }

		/// @brief Get underlying vector (mutable)
		Vector<Real>& Values() { return _values; }

		/// @brief Implicit conversion to Vector
		operator const Vector<Real>&() const { return _values; }

		/// @brief Raw element access
		Real& operator[](int i) { return _values[i]; }
		const Real& operator[](int i) const { return _values[i]; }

		//---------------------------------------------------------------------
		// Problem spec access
		//---------------------------------------------------------------------

		/// @brief Check if problem spec is attached
		bool HasSpec() const { return _spec != nullptr; }

		/// @brief Get attached problem spec
		const ProblemSpec& GetSpec() const {
			if (!_spec)
				throw OptimizationProblemError("DecisionVector has no ProblemSpec attached");
			return *_spec;
		}

		/// @brief Attach a problem spec
		void SetSpec(const ProblemSpec& spec) { _spec = &spec; }

		//---------------------------------------------------------------------
		// Type-safe getters
		//---------------------------------------------------------------------

		/// @brief Get value as Real (for continuous variables)
		Real GetReal(int i) const {
			return _values[i];
		}

		/// @brief Get value as integer (for Integer variables)
		int GetInteger(int i) const {
			if (_spec && _spec->GetType(i) != OptVariableType::Integer &&
			            _spec->GetType(i) != OptVariableType::Categorical)
				throw OptimizationProblemError("Variable " + std::to_string(i) + " is not Integer/Categorical");
			return static_cast<int>(std::round(_values[i]));
		}

		/// @brief Get value as boolean (for Binary variables)
		bool GetBinary(int i) const {
			if (_spec && _spec->GetType(i) != OptVariableType::Binary)
				throw OptimizationProblemError("Variable " + std::to_string(i) + " is not Binary");
			return _values[i] >= 0.5;
		}

		/// @brief Get value as category index (for Categorical variables)
		int GetCategorical(int i) const {
			if (_spec && _spec->GetType(i) != OptVariableType::Categorical)
				throw OptimizationProblemError("Variable " + std::to_string(i) + " is not Categorical");
			return static_cast<int>(std::round(_values[i]));
		}

		//---------------------------------------------------------------------
		// Type-safe setters
		//---------------------------------------------------------------------

		/// @brief Set Real value (for continuous variables)
		void SetReal(int i, Real value) {
			_values[i] = value;
		}

		/// @brief Set integer value
		void SetInteger(int i, int value) {
			_values[i] = static_cast<Real>(value);
		}

		/// @brief Set binary value
		void SetBinary(int i, bool value) {
			_values[i] = value ? 1.0 : 0.0;
		}

		/// @brief Set categorical value (0-indexed category)
		void SetCategorical(int i, int category) {
			_values[i] = static_cast<Real>(category);
		}

		//---------------------------------------------------------------------
		// Validation and repair
		//---------------------------------------------------------------------

		/// @brief Check if all values are valid according to spec
		bool IsValid() const {
			if (!_spec) return true;  // No spec = assume valid
			return _spec->IsValid(_values);
		}

		/// @brief Repair values to be valid according to spec
		DecisionVector& Repair() {
			if (_spec) {
				_values = _spec->Repair(_values);
			}
			return *this;
		}

		/// @brief Repair with specific integer repair strategy
		DecisionVector& RepairWithStrategy(IntegerRepairStrategy strategy) {
			if (!_spec) return *this;

			for (int i = 0; i < static_cast<int>(_values.size()); ++i) {
				const auto& varSpec = (*_spec)[i];
				
				// Clamp to bounds first
				_values[i] = std::max(varSpec.lowerBound, 
				             std::min(varSpec.upperBound, _values[i]));

				// Apply repair for discrete variables
				if (varSpec.IsDiscrete()) {
					switch (strategy) {
						case IntegerRepairStrategy::Nearest:
							_values[i] = std::round(_values[i]);
							break;
						case IntegerRepairStrategy::Floor:
							_values[i] = std::floor(_values[i]);
							break;
						case IntegerRepairStrategy::Ceil:
							_values[i] = std::ceil(_values[i]);
							break;
						case IntegerRepairStrategy::Probabilistic: {
							Real frac = _values[i] - std::floor(_values[i]);
							// Use simple deterministic pseudo-random based on value
							// For true randomness, use external RNG
							_values[i] = (frac > 0.5) ? std::ceil(_values[i]) : std::floor(_values[i]);
							break;
						}
					}
					// Final bounds clamp after rounding
					_values[i] = std::max(varSpec.lowerBound, 
					             std::min(varSpec.upperBound, _values[i]));
				}
			}
			return *this;
		}

		/// @brief Get integrality violation (sum of distances to nearest integers)
		Real GetIntegralityViolation() const {
			if (!_spec) return 0;
			return _spec->GetIntegralityViolation(_values);
		}

		/// @brief Check if all discrete variables have integer values
		bool IsIntegral(Real tolerance = 1e-9) const {
			if (!_spec) return true;
			
			for (int i = 0; i < static_cast<int>(_values.size()); ++i) {
				if ((*_spec)[i].IsDiscrete()) {
					Real frac = std::abs(_values[i] - std::round(_values[i]));
					if (frac > tolerance)
						return false;
				}
			}
			return true;
		}

		//---------------------------------------------------------------------
		// Factory methods
		//---------------------------------------------------------------------

		/// @brief Create from a Vector with attached spec
		static DecisionVector FromVector(const Vector<Real>& v, const ProblemSpec& spec) {
			return DecisionVector(v, spec);
		}

		/// @brief Create with random values within bounds
		static DecisionVector Random(const ProblemSpec& spec, unsigned int seed = 0) {
			DecisionVector dv(spec);
			std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
			
			for (int i = 0; i < spec.Dimension(); ++i) {
				const auto& varSpec = spec[i];
				std::uniform_real_distribution<Real> dist(varSpec.lowerBound, varSpec.upperBound);
				dv._values[i] = dist(rng);
			}
			
			// Repair to ensure valid discrete values
			dv.Repair();
			return dv;
		}
	};

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_VARIABLES_H
