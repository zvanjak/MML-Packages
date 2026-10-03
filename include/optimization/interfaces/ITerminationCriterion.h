#if !defined MML_OPTIMIZATION_ITERMINATIONCRITERION_H
#define MML_OPTIMIZATION_ITERMINATIONCRITERION_H

#include "OptimizationState.h"
#include "IIterativeAlgorithm.h"

#include <string>

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////////////////////
	// TERMINATION CRITERION INTERFACE
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Interface for termination criteria
	class ITerminationCriterion {
	public:
		virtual ~ITerminationCriterion() = default;

		/// @brief Check if optimization should terminate
		virtual bool ShouldTerminate(const OptimizationState& state) const = 0;

		/// @brief Get human-readable reason for termination
		virtual std::string GetReason() const = 0;
		virtual IterativeRunStatus GetStopStatus() const { return IterativeRunStatus::Stopped; }

		/// @brief Reset criterion for new optimization run
		virtual void Reset() = 0;
	};

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_ITERMINATIONCRITERION_H
