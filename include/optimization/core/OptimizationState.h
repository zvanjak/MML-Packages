#if !defined MML_OPTIMIZATION_STATE_H
#define MML_OPTIMIZATION_STATE_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"

#include <vector>
#include <string>
#include <limits>
#include <memory>
#include <functional>
#include <optional>

///////////////////////////////////////////////////////////////////////////////////////////
// OptimizationState.h - Core optimization state, result, and diagnostics
//
// Provides:
//   - Algorithm-specific diagnostics (GradientDiagnostics, SimplexDiagnostics, etc.)
//   - OptimizationState: Current state during optimization (core + optional diagnostics)
//   - OptimizationResult: Final result with trajectory and diagnostics
//   - CreateResult: build an OptimizationResult from a final OptimizationState
//
// The pluggable interfaces (ITerminationCriterion, IOptimizationObserver) live in
// optimization/interfaces/; OptimizationCommon.h is an umbrella over this file + those.
///////////////////////////////////////////////////////////////////////////////////////////

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////////////////////
	// ALGORITHM-SPECIFIC DIAGNOSTIC STRUCTS
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Diagnostics for gradient-based optimizers (BFGS, CG, Newton, etc.)
	struct GradientDiagnostics {
		int gradEvals = 0;          ///< Total gradient evaluations
		Real gradNorm = 0;        ///< Current/final gradient norm
		Real stepSize = 0;        ///< Last step size taken (line search)
		Real directionChange = 0; ///< Direction change angle (for CG, Powell)

		GradientDiagnostics() = default;
		GradientDiagnostics(int evals, Real norm, Real step = 0, Real dirChange = 0)
			: gradEvals(evals), gradNorm(norm), stepSize(step), directionChange(dirChange) {}
	};

	/// @brief Diagnostics for simplex-based optimizers (Nelder-Mead)
	struct SimplexDiagnostics {
		Real simplexSize = 0; ///< Current/final simplex diameter

		SimplexDiagnostics() = default;
		explicit SimplexDiagnostics(Real size) : simplexSize(size) {}
	};

	/// @brief Diagnostics for stochastic optimizers (SA, GA, MCMC)
	struct StochasticDiagnostics {
		Real temperature = 0; ///< Current/final temperature (SA)
		int acceptedMoves = 0;  ///< Number of accepted moves
		int rejectedMoves = 0;  ///< Number of rejected moves

		StochasticDiagnostics() = default;
		StochasticDiagnostics(Real temp, int accepted, int rejected)
			: temperature(temp), acceptedMoves(accepted), rejectedMoves(rejected) {}

		/// @brief Get acceptance rate
		Real GetAcceptanceRate() const {
			int total = acceptedMoves + rejectedMoves;
			return total > 0 ? static_cast<Real>(acceptedMoves) / total : Real(0);
		}
	};

	/// @brief Diagnostics for constrained optimization
	struct ConstraintDiagnostics {
		Real totalViolation = 0;  ///< Sum of all constraint violations
		Real maxViolation = 0;    ///< Maximum single constraint violation
		int numViolated = 0;        ///< Number of violated constraints
		bool isFeasible = true;     ///< Whether solution satisfies all constraints

		// Best feasible solution tracking
		Vector<Real> xBestFeasible; ///< Best feasible solution found (empty if none)
		Real fBestFeasible = std::numeric_limits<Real>::max(); ///< Objective of best feasible

		ConstraintDiagnostics() = default;

		/// @brief Check if a feasible solution has been found
		bool HasFeasibleSolution() const { return xBestFeasible.size() > 0; }

		/// @brief Update from current evaluation
		void Update(Real violation, Real maxViol, int numViol, bool feasible,
		            const Vector<Real>& x, Real f) {
			totalViolation = violation;
			maxViolation = maxViol;
			numViolated = numViol;
			isFeasible = feasible;

			// Track best feasible solution
			if (feasible && f < fBestFeasible) {
				xBestFeasible = x;
				fBestFeasible = f;
			}
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// OPTIMIZATION STATE - Current state during optimization
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Encapsulates the current state of an optimization algorithm
	/// @details Core fields are universal; algorithm-specific diagnostics are optional
	struct OptimizationState {
		//----------------------------------------------------------------------
		// CORE STATE (universal for all optimizers)
		//----------------------------------------------------------------------

		// Iteration counters
		int iteration = 0;  ///< Current iteration number (0-based)
		int funcEvals = 0;  ///< Total function evaluations so far

		// Current solution
		Vector<Real> xCurrent; ///< Current solution vector
		Real fCurrent = std::numeric_limits<Real>::max(); ///< Current function value

		// Best solution so far
		Vector<Real> xBest; ///< Best solution found so far
		Real fBest = std::numeric_limits<Real>::max(); ///< Best function value

		// Stagnation tracking
		int iterSinceImprovement = 0; ///< Iterations without improvement
		Real fBestPrevious = std::numeric_limits<Real>::max(); ///< Previous best (for tracking)

		// Timing
		double elapsedTime = 0.0; ///< Wall-clock time since start (seconds)

		//----------------------------------------------------------------------
		// ALGORITHM-SPECIFIC DIAGNOSTICS (optional)
		//----------------------------------------------------------------------

		std::optional<GradientDiagnostics> gradient;     ///< For gradient-based methods
		std::optional<SimplexDiagnostics> simplex;       ///< For Nelder-Mead
		std::optional<StochasticDiagnostics> stochastic;  ///< For SA, GA
		std::optional<ConstraintDiagnostics> constraint; ///< For constrained optimization

		//----------------------------------------------------------------------
		// CONSTRUCTORS
		//----------------------------------------------------------------------

		OptimizationState() = default;

		//----------------------------------------------------------------------
		// HELPER METHODS
		//----------------------------------------------------------------------

		/// @brief Check if there was improvement in this iteration
		bool HasImproved() const {
			return fBest < fBestPrevious - std::numeric_limits<Real>::epsilon();
		}

		/// @brief Check if gradient diagnostics are available
		bool HasGradientInfo() const { return gradient.has_value(); }

		/// @brief Check if simplex diagnostics are available
		bool HasSimplexInfo() const { return simplex.has_value(); }

		/// @brief Check if stochastic diagnostics are available
		bool HasStochasticInfo() const { return stochastic.has_value(); }

		/// @brief Check if constraint diagnostics are available
		bool HasConstraintInfo() const { return constraint.has_value(); }

		/// @brief Check if any feasible solution has been found (constrained only)
		bool HasFeasibleSolution() const {
			return constraint && constraint->HasFeasibleSolution();
		}

		/// @brief Enable gradient diagnostics
		GradientDiagnostics& EnableGradient() {
			if (!gradient) gradient.emplace();
			return *gradient;
		}

		/// @brief Enable simplex diagnostics
		SimplexDiagnostics& EnableSimplex() {
			if (!simplex) simplex.emplace();
			return *simplex;
		}

		/// @brief Enable stochastic diagnostics
		StochasticDiagnostics& EnableStochastic() {
			if (!stochastic) stochastic.emplace();
			return *stochastic;
		}

		/// @brief Enable constraint diagnostics
		ConstraintDiagnostics& EnableConstraint() {
			if (!constraint) constraint.emplace();
			return *constraint;
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// OPTIMIZATION RESULT - Final result with diagnostics
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Enhanced optimization result with trajectory and diagnostics
	struct OptimizationResult {
		//----------------------------------------------------------------------
		// CORE RESULTS (universal)
		//----------------------------------------------------------------------

		// Solution
		Vector<Real> xBest; ///< Best solution found
		Real fBest = std::numeric_limits<Real>::max(); ///< Best function value

		// Convergence information
		bool converged = false;                      ///< Whether optimization converged
		std::string terminationReason = "Not started"; ///< Human-readable reason

		// Statistics
		int iterations = 0;    ///< Total iterations performed
		int funcEvals = 0;     ///< Total function evaluations
		double elapsedTime = 0.0; ///< Total wall-clock time (seconds)

		//----------------------------------------------------------------------
		// ALGORITHM-SPECIFIC DIAGNOSTICS (optional)
		//----------------------------------------------------------------------

		std::optional<GradientDiagnostics> gradient;     ///< For gradient-based methods
		std::optional<SimplexDiagnostics> simplex;       ///< For Nelder-Mead
		std::optional<StochasticDiagnostics> stochastic;  ///< For SA, GA
		std::optional<ConstraintDiagnostics> constraint; ///< For constrained optimization

		//----------------------------------------------------------------------
		// TRAJECTORY (optional)
		//----------------------------------------------------------------------

		std::vector<OptimizationState> trajectory; ///< Recorded states (if observer used)

		//----------------------------------------------------------------------
		// CONSTRUCTORS
		//----------------------------------------------------------------------

		OptimizationResult() = default;

		//----------------------------------------------------------------------
		// HELPER METHODS
		//----------------------------------------------------------------------

		/// @brief Check if trajectory was recorded
		bool HasTrajectory() const { return !trajectory.empty(); }

		/// @brief Check if gradient diagnostics are available
		bool HasGradientInfo() const { return gradient.has_value(); }

		/// @brief Check if simplex diagnostics are available
		bool HasSimplexInfo() const { return simplex.has_value(); }

		/// @brief Check if stochastic diagnostics are available
		bool HasStochasticInfo() const { return stochastic.has_value(); }

		/// @brief Check if constraint diagnostics are available
		bool HasConstraintInfo() const { return constraint.has_value(); }

		/// @brief Check if a feasible solution was found (constrained only)
		bool HasFeasibleSolution() const {
			return constraint && constraint->HasFeasibleSolution();
		}

		/// @brief Check if the best solution is feasible (constrained only)
		bool IsBestFeasible() const {
			return !constraint || constraint->isFeasible;
		}

		/// @brief Get best feasible objective value (max if none found)
		Real GetBestFeasibleValue() const {
			if (HasFeasibleSolution()) {
				return constraint->fBestFeasible;
			}
			return std::numeric_limits<Real>::max();
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// HELPER FUNCTIONS
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Create an OptimizationResult from final OptimizationState
	inline OptimizationResult CreateResult(const OptimizationState& state,
	                                         const std::string& reason, bool converged) {
		OptimizationResult result;

		// Core results
		result.xBest = state.xBest;
		result.fBest = state.fBest;
		result.converged = converged;
		result.terminationReason = reason;
		result.iterations = state.iteration;
		result.funcEvals = state.funcEvals;
		result.elapsedTime = state.elapsedTime;

		// Copy algorithm-specific diagnostics (if present)
		result.gradient = state.gradient;
		result.simplex = state.simplex;
		result.stochastic = state.stochastic;
		result.constraint = state.constraint;

		return result;
	}

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_STATE_H
