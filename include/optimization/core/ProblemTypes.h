///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ProblemTypes.h                                                      ///
///  Description: Single-objective, multi-objective, and constrained problem types    ///
///               (split from OptimizationProblem.h)                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_OPTIMIZATION_PROBLEMTYPES_H
#define MML_OPTIMIZATION_PROBLEMTYPES_H

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
#include <type_traits>

#include "IOptimizationProblem.h"
#include "Constraints.h"
#include "Variables.h"

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////
	///                  SINGLE-OBJECTIVE PROBLEM                           ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Unconstrained single-objective optimization problem
	 * 
	 * Wraps a scalar objective function f: R^n → R.
	 * Supports continuous, integer, binary, and categorical variables.
	 * The number of decision variables is a runtime property.
	 * 
	 * For constrained problems, wrap with ConstrainedSingleObjectiveProblem.
	 * 
	 * Usage:
	 * @code
	 * // Define objective function
	 * auto objective = [](const Vector<Real>& x) {
	 *     return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	 * };
	 *
	 * SingleObjectiveProblem problem(3, objective);
	 * problem.SetBounds(-10, 10);
	 * ConstrainedSingleObjectiveProblem constrained(problem);
	 * constrained.AddInequalityConstraint([](const auto& x) { return x[0] + x[1] - 5; });
	 * @endcode
	 */
	class SingleObjectiveProblem : public ISingleObjectiveProblem {
	private:
		std::function<Real(const Vector<Real>&)> _objective;
		ProblemSpec _problemSpec;

	public:
		//---------------------------------------------------------------------
		// Constructors
		//---------------------------------------------------------------------

		/// @brief Default constructor (objective and dimension must be set later)
		SingleObjectiveProblem() = default;

		/// @brief Construct with objective function only (dimension set later via spec)
		explicit SingleObjectiveProblem(std::function<Real(const Vector<Real>&)> objective)
			: _objective(std::move(objective)) {}

		/// @brief Construct with dimension and objective function
		SingleObjectiveProblem(int dim, std::function<Real(const Vector<Real>&)> objective)
			: _objective(std::move(objective))
			, _problemSpec(dim) {}

		/// @brief Construct with objective and problem spec
		SingleObjectiveProblem(
			std::function<Real(const Vector<Real>&)> objective,
			const ProblemSpec& spec)
			: _objective(std::move(objective))
			, _problemSpec(spec) {}

		//---------------------------------------------------------------------
		// Objective function
		//---------------------------------------------------------------------

		/// @brief Set the objective function
		void SetObjective(std::function<Real(const Vector<Real>&)> objective) {
			_objective = std::move(objective);
		}

		/// @brief Evaluate objective function at x
		Real Evaluate(const Vector<Real>& x) const override {
			if (Dimension() <= 0 || x.size() != static_cast<size_t>(Dimension()))
				throw OptimizationProblemError("Objective input dimension mismatch");
			if (!_objective)
				throw OptimizationProblemError("Objective function not set");
			return _objective(x);
		}

		/// @brief Evaluate objective (operator() syntax)
		Real operator()(const Vector<Real>& x) const {
			return Evaluate(x);
		}

		//---------------------------------------------------------------------
		// IOptimizationProblem interface implementation
		//---------------------------------------------------------------------

		const ProblemSpec& GetProblemSpec() const override {
			return _problemSpec;
		}

		/// @brief Get mutable access to problem spec for configuration
		ProblemSpec& GetProblemSpecMutable() {
			return _problemSpec;
		}

		/// @brief Unconstrained problem has no constraints
		int GetNumConstraints() const override { return 0; }

		/// @brief Throws - no constraints in unconstrained problem
		const IConstraint& GetConstraint(int /*index*/) const override {
			throw OptimizationProblemError("Unconstrained problem has no constraints");
		}

		int GetNumObjectives() const override { return 1; }

		//---------------------------------------------------------------------
		// Fluent API for problem configuration
		//---------------------------------------------------------------------

		/// @brief Set the problem dimension (resizes the spec to n continuous [0,1] variables)
		SingleObjectiveProblem& SetDimension(int n) {
			_problemSpec = ProblemSpec(n);
			return *this;
		}

		/// @brief Set bounds for all variables (continuous)
		SingleObjectiveProblem& SetBounds(Real lower, Real upper) {
			_problemSpec.SetAllContinuous(lower, upper);
			return *this;
		}

		/// @brief Set bounds using vectors
		SingleObjectiveProblem& SetBounds(
			const Vector<Real>& lower, 
			const Vector<Real>& upper) {
			if (_problemSpec.Dimension() != static_cast<int>(lower.size()))
				_problemSpec = ProblemSpec(static_cast<int>(lower.size()));
			for (int i = 0; i < _problemSpec.Dimension(); ++i) {
				_problemSpec.SetVariable(i, OptVariableSpec::Continuous(lower[i], upper[i]));
			}
			return *this;
		}

		/// @brief Set a single variable specification
		SingleObjectiveProblem& SetVariable(int index, const OptVariableSpec& spec) {
			_problemSpec.SetVariable(index, spec);
			return *this;
		}
	};


	///////////////////////////////////////////////////////////////////////////
	///                  MULTI-OBJECTIVE PROBLEM BASE                       ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Abstract base class for multi-objective optimization problems
	 * 
	 * Extends IOptimizationProblem to provide a standard interface for
	 * problems with M objectives. This is the base class for:
	 * - Custom multi-objective problems
	 * - Benchmark test problems (ZDT, DTLZ, etc.)
	 * 
	 * Both the number of variables and the number of objectives are runtime
	 * properties, supplied by derived classes through the protected constructor.
	 * 
	 * For single-objective problems, use SingleObjectiveProblem instead.
	 * 
	 * Usage:
	 * @code
	 * class MyBiObjectiveProblem : public MultiObjectiveProblem {
	 * public:
	 *     MyBiObjectiveProblem() : MultiObjectiveProblem(5, 2) {}
	 *     Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override {
	 *         Real f1 = x[0];  // First objective
	 *         Real f2 = 1 - sqrt(f1);  // Second objective
	 *         return Vector<Real>{f1, f2};
	 *     }
	 *     std::pair<Real, Real> GetVariableBounds(int i) const override {
	 *         return {0.0, 1.0};  // All variables in [0, 1]
	 *     }
	 *     std::string GetProblemName() const override { return "MyProblem"; }
	 * };
	 * @endcode
	 */
	class MultiObjectiveProblem : public IMultiObjectiveProblem {
	private:
		int _numVars;
		int _numObjectives;
		mutable ProblemSpec _cachedSpec;
		mutable bool _specInitialized = false;

	protected:
		/// @brief Construct with runtime variable and objective counts
		MultiObjectiveProblem(int numVars, int numObjectives)
			: _numVars(numVars), _numObjectives(numObjectives) {
			if (numVars <= 0)
				throw OptimizationProblemError("MultiObjectiveProblem requires a positive dimension");
			if (numObjectives < 2)
				throw OptimizationProblemError("MultiObjectiveProblem requires at least 2 objectives");
		}

	public:
		virtual ~MultiObjectiveProblem() = default;

		//---------------------------------------------------------------------
		// Core multi-objective interface (pure virtual - must implement)
		//---------------------------------------------------------------------

		/// @brief Evaluate all objectives at decision point x
		/// @param x Decision variable vector
		/// @return Vector of objective values (size = NumObjectives())
		virtual Vector<Real> EvaluateObjectives(const Vector<Real>& x) const = 0;

		/// @brief Get bounds for decision variable i
		/// @param i Variable index (0 to NumVariables()-1)
		/// @return Pair of (lower, upper) bounds
		virtual std::pair<Real, Real> GetVariableBounds(int i) const = 0;

		/// @brief Get human-readable problem name
		virtual std::string GetProblemName() const = 0;

		//---------------------------------------------------------------------
		// Optional: Pareto front validation (for benchmark problems)
		//---------------------------------------------------------------------

		/// @brief Check if a decision point lies on the true Pareto front
		/// Default implementation returns false (unknown).
		virtual bool IsOnParetoFront(const Vector<Real>& x, 
									 Real tolerance = 1e-6) const {
			(void)x; (void)tolerance;
			return false;  // Unknown by default
		}

		/// @brief Check if objective values are on the true Pareto front
		/// Default implementation returns false (unknown).
		virtual bool ObjectivesOnParetoFront(const Vector<Real>& objectives,
											 Real tolerance = 1e-6) const {
			(void)objectives; (void)tolerance;
			return false;  // Unknown by default
		}

		//---------------------------------------------------------------------
		// Dimension accessors
		//---------------------------------------------------------------------

		/// @brief Get number of decision variables (runtime)
		int NumVariables() const { return _numVars; }

		/// @brief Get number of objectives (runtime)
		int NumObjectives() const { return _numObjectives; }

		//---------------------------------------------------------------------
		// IOptimizationProblem interface implementation
		//---------------------------------------------------------------------

		/// @brief Get problem specification (bounds)
		const ProblemSpec& GetProblemSpec() const override {
			if (!_specInitialized) {
				_cachedSpec = ProblemSpec(_numVars);
				for (int i = 0; i < _numVars; ++i) {
					auto [lo, hi] = GetVariableBounds(i);
					_cachedSpec.SetVariable(i, OptVariableSpec::Continuous(lo, hi));
				}
				_specInitialized = true;
			}
			return _cachedSpec;
		}

		/// @brief Multi-objective problems have no constraints by default
		/// Override in derived classes to add constraints
		int GetNumConstraints() const override { return 0; }

		/// @brief Throws - base class has no constraints
		const IConstraint& GetConstraint(int /*index*/) const override {
			throw OptimizationProblemError("MultiObjectiveProblem base has no constraints");
		}

		/// @brief Returns number of objectives
		int GetNumObjectives() const override { return _numObjectives; }

		//---------------------------------------------------------------------
		// Utility methods
		//---------------------------------------------------------------------

		/// @brief Get all variable bounds as vectors
		void GetAllBounds(Vector<Real>& lower, Vector<Real>& upper) const {
			lower = Vector<Real>(_numVars);
			upper = Vector<Real>(_numVars);
			for (int i = 0; i < _numVars; ++i) {
				auto [lo, hi] = GetVariableBounds(i);
				lower[i] = lo;
				upper[i] = hi;
			}
		}

		/// @brief Evaluate single objective by index
		Real EvaluateSingleObjective(const Vector<Real>& x, int objectiveIndex) const {
			if (objectiveIndex < 0 || objectiveIndex >= _numObjectives)
				throw OptimizationProblemError("Objective index out of range");
			return EvaluateValidated(x)[objectiveIndex];
		}

		//---------------------------------------------------------------------
		// Aliases for compatibility with algorithms (NSGA2, MOEAD)
		//---------------------------------------------------------------------

		/// @brief Alias for EvaluateObjectives (backward compatibility)
		Vector<Real> Evaluate(const Vector<Real>& x) const {
			return EvaluateValidated(x);
		}

		/// @brief Alias for GetVariableBounds (backward compatibility)
		std::pair<Real, Real> GetBounds(int i) const {
			return GetVariableBounds(i);
		}

		/// @brief Alias for GetProblemName (backward compatibility)
		std::string GetName() const {
			return GetProblemName();
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                    CONSTRAINED PROBLEM DECORATOR                    ///
	///////////////////////////////////////////////////////////////////////////

	namespace detail {
	class ConstraintCollection {
	private:
		const IOptimizationProblem& _baseProblem;
		std::vector<std::shared_ptr<IConstraint>> _constraints;

	public:
		explicit ConstraintCollection(const IOptimizationProblem& baseProblem)
			: _baseProblem(baseProblem) {}

		//---------------------------------------------------------------------
		// Constraint management
		//---------------------------------------------------------------------

		int GetNumConstraints() const {
			return _baseProblem.GetNumConstraints() + static_cast<int>(_constraints.size());
		}

		const IConstraint& GetConstraint(int index) const {
			if (index < 0 || index >= GetNumConstraints())
				throw OptimizationProblemError("Constraint index out of range");
			int inherited = _baseProblem.GetNumConstraints();
			return index < inherited ? _baseProblem.GetConstraint(index) : *_constraints[index - inherited];
		}

		/// @brief Add a constraint
		void AddConstraint(std::shared_ptr<IConstraint> constraint) {
			if (!constraint)
				throw OptimizationProblemError("Constraint must not be null");
			_constraints.push_back(std::move(constraint));
		}

		/// @brief Add inequality constraint g(x) <= 0
		void AddInequalityConstraint(
			std::function<Real(const Vector<Real>&)> func,
			const std::string& name = "") {
			_constraints.push_back(std::make_shared<FunctionConstraint>(
				std::move(func), ConstraintType::Inequality, name));
		}

		/// @brief Add equality constraint g(x) = 0
		void AddEqualityConstraint(
			std::function<Real(const Vector<Real>&)> func,
			const std::string& name = "") {
			_constraints.push_back(std::make_shared<FunctionConstraint>(
				std::move(func), ConstraintType::Equality, name));
		}

		/// @brief Clear all constraints
		void ClearConstraints() {
			_constraints.clear();
		}

		//---------------------------------------------------------------------
		// Constraint statistics
		//---------------------------------------------------------------------

		/// @brief Get number of inequality constraints
		int GetNumInequalityConstraints() const {
			int count = 0;
			for (int i = 0; i < GetNumConstraints(); ++i)
				if (GetConstraint(i).GetType() == ConstraintType::Inequality) ++count;
			return count;
		}

		/// @brief Get number of equality constraints
		int GetNumEqualityConstraints() const {
			int count = 0;
			for (int i = 0; i < GetNumConstraints(); ++i)
				if (GetConstraint(i).GetType() == ConstraintType::Equality) ++count;
			return count;
		}
	};
	}

	template<typename Interface>
	class TypedConstrainedProblemBase : public Interface {
	private:
		std::shared_ptr<Interface> _ownedProblem;
		const Interface* _baseProblem;
		detail::ConstraintCollection _constraints;

		static const Interface& RequireBase(const Interface* problem) {
			if (!problem)
				throw OptimizationProblemError("Typed constrained problem requires a problem");
			return *problem;
		}

	public:
		explicit TypedConstrainedProblemBase(const Interface& problem)
			: _baseProblem(&problem), _constraints(problem) {}

		explicit TypedConstrainedProblemBase(std::shared_ptr<Interface> problem)
			: _ownedProblem(std::move(problem)), _baseProblem(_ownedProblem.get()),
			  _constraints(RequireBase(_baseProblem)) {}

		const Interface& GetBaseProblem() const { return *_baseProblem; }
		const ProblemSpec& GetProblemSpec() const override { return _baseProblem->GetProblemSpec(); }
		int GetNumObjectives() const override { return _baseProblem->GetNumObjectives(); }
		int GetNumConstraints() const override { return _constraints.GetNumConstraints(); }
		const IConstraint& GetConstraint(int index) const override { return _constraints.GetConstraint(index); }

		void AddConstraint(std::shared_ptr<IConstraint> constraint) { _constraints.AddConstraint(std::move(constraint)); }
		void AddInequalityConstraint(std::function<Real(const Vector<Real>&)> func, const std::string& name = "") {
			_constraints.AddInequalityConstraint(std::move(func), name);
		}
		void AddEqualityConstraint(std::function<Real(const Vector<Real>&)> func, const std::string& name = "") {
			_constraints.AddEqualityConstraint(std::move(func), name);
		}
		void ClearConstraints() { _constraints.ClearConstraints(); }
		int GetNumInequalityConstraints() const {
			int count = 0;
			for (int i = 0; i < GetNumConstraints(); ++i)
				if (GetConstraint(i).GetType() == ConstraintType::Inequality) ++count;
			return count;
		}
		int GetNumEqualityConstraints() const {
			int count = 0;
			for (int i = 0; i < GetNumConstraints(); ++i)
				if (GetConstraint(i).GetType() == ConstraintType::Equality) ++count;
			return count;
		}
	};

	class ConstrainedSingleObjectiveProblem : public TypedConstrainedProblemBase<ISingleObjectiveProblem> {
	public:
		using TypedConstrainedProblemBase::TypedConstrainedProblemBase;

		Real Evaluate(const Vector<Real>& x) const override {
			if (GetNumObjectives() != 1)
				throw OptimizationProblemError("Single-objective problem requires exactly 1 objective");
			if (Dimension() <= 0 || x.size() != static_cast<size_t>(Dimension()))
				throw OptimizationProblemError("Objective input dimension mismatch");
			return GetBaseProblem().Evaluate(x);
		}
	};

	class ConstrainedMultiObjectiveProblem : public TypedConstrainedProblemBase<IMultiObjectiveProblem> {
	public:
		using TypedConstrainedProblemBase::TypedConstrainedProblemBase;

		Vector<Real> EvaluateObjectives(const Vector<Real>& x) const override {
			return GetBaseProblem().EvaluateValidated(x);
		}
	};

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_PROBLEMTYPES_H
