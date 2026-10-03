#if !defined MML_OPTIMIZATION_IOPTIMIZATIONOBSERVER_H
#define MML_OPTIMIZATION_IOPTIMIZATIONOBSERVER_H

#include "OptimizationState.h"

#include <string>

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////////////////////
	// OPTIMIZATION OBSERVER INTERFACE
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Observer interface for monitoring optimization progress
	class IOptimizationObserver {
	public:
		virtual ~IOptimizationObserver() = default;

		/// @brief Called at the start of optimization
		virtual void OnStart(const OptimizationState& state) = 0;

		/// @brief Called after each iteration
		/// @return true to continue, false to abort
		virtual bool OnIteration(const OptimizationState& state) = 0;

		/// @brief Called when optimization completes
		virtual void OnComplete(const OptimizationState& state, const std::string& reason) = 0;
	};

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_IOPTIMIZATIONOBSERVER_H
