#if !defined MML_OPTIMIZATION_COMMON_H
#define MML_OPTIMIZATION_COMMON_H

///////////////////////////////////////////////////////////////////////////////////////////
// OptimizationCommon.h - Umbrella for the optimization framework's core infrastructure
//
// Aggregates:
//   - core/OptimizationState.h          : diagnostics, OptimizationState, OptimizationResult, CreateResult
//   - interfaces/ITerminationCriterion.h: pluggable termination-criterion interface
//   - interfaces/IOptimizationObserver.h: progress-monitoring observer interface
//
// Existing code that includes "OptimizationCommon.h" keeps compiling unchanged.
///////////////////////////////////////////////////////////////////////////////////////////

#include "OptimizationState.h"
#include "ITerminationCriterion.h"
#include "IOptimizationObserver.h"
#include "IIterativeAlgorithm.h"

#endif // MML_OPTIMIZATION_COMMON_H
