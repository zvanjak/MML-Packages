#ifndef MML_OPTIMIZATION_ITERATIVE_RUNNER_H
#define MML_OPTIMIZATION_ITERATIVE_RUNNER_H

#include "IIterativeAlgorithm.h"
#include "OptimizationConfig.h"

#include <chrono>
#include <cmath>
#include <functional>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace MML::Optimization {

struct IterationStop {
	IterativeRunStatus status;
	std::string reason;
};

template<typename ProgressType>
class IterativeRunner {
public:
	struct Observer {
		std::function<void(const ProgressType&)> onStart;
		std::function<bool(const ProgressType&)> onIteration;
		std::function<void(const ProgressType&, const std::string&)> onComplete;
	};

	IterativeRunner(ProgressType& progress, int maxIterations, int maxEvaluations = 0, double maxSeconds = 0.0)
		: _progress(progress), _maxIterations(maxIterations), _maxEvaluations(maxEvaluations), _maxSeconds(maxSeconds) {
		if (maxIterations <= 0 || maxEvaluations < 0 || !std::isfinite(maxSeconds) || maxSeconds < 0.0)
			throw std::invalid_argument("Invalid iterative run limits");
	}

	void SetCriterion(std::function<std::optional<IterationStop>(const ProgressType&)> criterion,
	                  std::function<void()> reset = {}) {
		if (_progress.status == IterativeRunStatus::Running)
			throw std::logic_error("Cannot change an active run");
		_criterion = std::move(criterion);
		_resetCriterion = std::move(reset);
	}

	void AddObserver(Observer observer) {
		if (_progress.status == IterativeRunStatus::Running)
			throw std::logic_error("Cannot change an active run");
		_observers.push_back(std::move(observer));
	}

	void ConfigureScalar(OptimizationConfig& config, OptimizationState& state)
		requires std::is_same_v<ProgressType, ScalarIterationProgress> {
		if (_progress.status == IterativeRunStatus::Running)
			throw std::logic_error("Cannot change an active run");
		auto* criterion = config.GetCriterion();
		if (!criterion)
			throw std::invalid_argument("Scalar run requires a criterion");
		SetCriterion([criterion, &state](const ProgressType&) -> std::optional<IterationStop> {
			if (criterion->ShouldTerminate(state))
				return IterationStop{criterion->GetStopStatus(), criterion->GetReason()};
			return std::nullopt;
		}, [criterion] { criterion->Reset(); });
		_syncState = [&state](const ProgressType& progress) {
			state.iteration = progress.iterations;
			state.funcEvals = progress.funcEvals;
			state.elapsedTime = progress.elapsedTime;
		};
		_observers.clear();
		for (const auto& observer : config.GetObservers()) {
			AddObserver({[observer, &state](const auto&) { observer->OnStart(state); },
			             [observer, &state](const auto&) { return observer->OnIteration(state); },
			             [observer, &state](const auto&, const auto& reason) { observer->OnComplete(state, reason); }});
		}
	}

	void Initialize() {
		if (_progress.status == IterativeRunStatus::Running)
			throw std::logic_error("Cannot initialize an active run");
		if (_progress.funcEvals < 0 || (_maxEvaluations && _progress.funcEvals > _maxEvaluations)) {
			_progress.status = IterativeRunStatus::Failed;
			throw std::invalid_argument("Initial evaluations exceed the run limit");
		}
		_start = std::chrono::steady_clock::now();
		_progress.iterations = 0;
		_progress.elapsedTime = 0.0;
		_progress.terminationReason.clear();
		_progress.status = IterativeRunStatus::Running;
		try {
			if (_resetCriterion)
				_resetCriterion();
			if (_syncState)
				_syncState(_progress);
			for (auto& observer : _observers)
				if (observer.onStart)
					observer.onStart(_progress);
			CheckStop();
		} catch (...) {
			_progress.status = IterativeRunStatus::Failed;
			throw;
		}
	}

	template<typename StepAction>
	void Step(int nextEvaluations, StepAction&& action) {
		if (_progress.status != IterativeRunStatus::Running)
			throw std::logic_error("No active iterative run");
		if (nextEvaluations < 0)
			throw std::invalid_argument("Step evaluations cannot be negative");
		try {
			if (_progress.iterations >= _maxIterations) {
				Finish(IterativeRunStatus::LimitReached, "Maximum iterations reached");
				return;
			}
			if (_maxEvaluations && nextEvaluations > _maxEvaluations - _progress.funcEvals) {
				Finish(IterativeRunStatus::LimitReached, "Maximum evaluations reached");
				return;
			}
			if (_maxSeconds && std::chrono::duration<double>(std::chrono::steady_clock::now() - _start).count() >= _maxSeconds) {
				Finish(IterativeRunStatus::LimitReached, "Time limit reached");
				return;
			}
			std::forward<StepAction>(action)();
			++_progress.iterations;
			_progress.funcEvals += nextEvaluations;
			_progress.elapsedTime = std::chrono::duration<double>(std::chrono::steady_clock::now() - _start).count();
			if (_syncState)
				_syncState(_progress);
			bool continueRun = true;
			for (auto& observer : _observers)
				if (observer.onIteration && !observer.onIteration(_progress))
					continueRun = false;
			if (!continueRun)
				Finish(IterativeRunStatus::Aborted, "Aborted by observer");
			else
				CheckStop();
		} catch (...) {
			_progress.status = IterativeRunStatus::Failed;
			throw;
		}
	}

	template<typename StepAction>
	void Run(int nextEvaluations, StepAction&& action) {
		if (_progress.status == IterativeRunStatus::Uninitialized)
			throw std::logic_error("No initialized iterative run");
		while (_progress.status == IterativeRunStatus::Running)
			Step(nextEvaluations, action);
	}

	void Finish(IterativeRunStatus status, std::string reason) {
		if (_progress.status != IterativeRunStatus::Running)
			throw std::logic_error("No active iterative run");
		if (status == IterativeRunStatus::Running || status == IterativeRunStatus::Uninitialized)
			throw std::invalid_argument("Completion requires a terminal status");
		_progress.status = status;
		_progress.terminationReason = std::move(reason);
		try {
			for (auto& observer : _observers)
				if (observer.onComplete)
					observer.onComplete(_progress, _progress.terminationReason);
		} catch (...) {
			_progress.status = IterativeRunStatus::Failed;
			throw;
		}
	}

private:
	void CheckStop() {
		if (_criterion) {
			auto stop = _criterion(_progress);
			if (stop) {
				Finish(stop->status, std::move(stop->reason));
				return;
			}
		}
		if (_progress.iterations >= _maxIterations)
			Finish(IterativeRunStatus::LimitReached, "Maximum iterations reached");
		else if (_maxEvaluations && _progress.funcEvals >= _maxEvaluations)
			Finish(IterativeRunStatus::LimitReached, "Maximum evaluations reached");
		else if (_maxSeconds && _progress.elapsedTime >= _maxSeconds)
			Finish(IterativeRunStatus::LimitReached, "Time limit reached");
	}

	ProgressType& _progress;
	int _maxIterations;
	int _maxEvaluations;
	double _maxSeconds;
	std::chrono::steady_clock::time_point _start;
	std::function<std::optional<IterationStop>(const ProgressType&)> _criterion;
	std::function<void()> _resetCriterion;
	std::function<void(const ProgressType&)> _syncState;
	std::vector<Observer> _observers;
};

} // namespace MML::Optimization

#endif // MML_OPTIMIZATION_ITERATIVE_RUNNER_H