#if !defined MML_OPTIMIZATION_IITERATIVE_ALGORITHM_H
#define MML_OPTIMIZATION_IITERATIVE_ALGORITHM_H

#include "OptimizationState.h"

#include <span>
#include <string>

namespace MML::Optimization {

class ObjectivePoint;

enum class IterativeRunStatus {
	Uninitialized,
	Running,
	Converged,
	LimitReached,
	Stopped,
	Aborted,
	Failed
};

struct IterationProgress {
	int iterations = 0;
	int funcEvals = 0;
	double elapsedTime = 0.0;
	IterativeRunStatus status = IterativeRunStatus::Uninitialized;
	std::string terminationReason;
};

struct ScalarIterationProgress : IterationProgress {
	const Vector<Real>* xCurrent = nullptr;
	Real fCurrent = std::numeric_limits<Real>::max();
	const Vector<Real>* xBest = nullptr;
	Real fBest = std::numeric_limits<Real>::max();
	int iterSinceImprovement = 0;
	std::optional<StochasticDiagnostics> stochastic;
};

struct ParetoIterationProgress : IterationProgress {
	int numObjectives = 0;
	int paretoFrontSize = 0;
	std::span<const ObjectivePoint> paretoFront;
	std::span<const ObjectivePoint> population;
	const Vector<Real>* idealPoint = nullptr;
};

class IIterativeAlgorithm {
public:
	virtual ~IIterativeAlgorithm() = default;
	virtual void Initialize() = 0;
	virtual void Step() = 0;
	virtual const IterationProgress& Progress() const = 0;
	virtual void Run() = 0;
};

} // namespace MML::Optimization

#endif // MML_OPTIMIZATION_IITERATIVE_ALGORITHM_H