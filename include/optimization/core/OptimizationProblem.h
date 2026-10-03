///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        OptimizationProblem.h                                               ///
///  Description: Unified optimization problem interfaces                              ///
///               IOptimizationProblem<N> - Abstract base for all problem types       ///
///               ProblemSpec<N> - Variable specification (bounds, types)             ///
///               Constraint types and handling                                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_OPTIMIZATION_PROBLEM_H
#define MML_OPTIMIZATION_PROBLEM_H

#include "Variables.h"
#include "IConstraint.h"
#include "IOptimizationProblem.h"
#include "Constraints.h"
#include "ProblemTypes.h"
#include "PenaltyMethod.h"
#include "FunctionAdapters.h"

#endif // MML_OPTIMIZATION_PROBLEM_H
