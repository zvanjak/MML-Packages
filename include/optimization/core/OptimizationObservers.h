#if !defined MML_OPTIMIZATION_OBSERVERS_H
#define MML_OPTIMIZATION_OBSERVERS_H

#include "OptimizationCommon.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

///////////////////////////////////////////////////////////////////////////////////////////
// OptimizationObservers.h - Standard observers for optimization progress monitoring
//
// Provides concrete implementations of IOptimizationObserver:
//   - ConsoleObserver: Prints progress to console
//   - TrajectoryObserver: Saves OptimizationState snapshots with CSV export
//   - CallbackObserver: Lambda-friendly custom callbacks
//   - NullObserver: No-op observer (for disabling observation)
//
// Observers enable:
//   - Progress monitoring and logging
//   - Intermediate result saving (checkpointing)
//   - Custom early stopping logic
//   - GUI updates, plotting, analysis
//
// Usage patterns:
//   - Single observer: direct use
//   - Multiple observers: wrap in MultiObserver
//   - Custom behavior: use CallbackObserver with lambdas
///////////////////////////////////////////////////////////////////////////////////////////

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////////////////////
	// NULL OBSERVER (for disabling observation)
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief No-op observer that does nothing
	/// @details Useful as a default when no observation needed
	class NullObserver : public IOptimizationObserver {
	public:
		void OnStart(const OptimizationState&) override {}
		bool OnIteration(const OptimizationState&) override { return true; }
		void OnComplete(const OptimizationState&, const std::string&) override {}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// CONSOLE OBSERVER
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Prints optimization progress to console
	/// @details Configurable output frequency and verbosity
	class ConsoleObserver : public IOptimizationObserver {
	private:
		int _printEvery;	///< Print every N iterations (1 = every iteration)
		bool _verbose;		///< Verbose output (more details)
		std::ostream& _out; ///< Output stream (default: std::cout)

	public:
		/// @brief Constructor
		/// @param printEvery Print every N iterations (default: 1)
		/// @param verbose Verbose output with more details (default: false)
		/// @param out Output stream (default: std::cout)
		explicit ConsoleObserver(int printEvery = 1, bool verbose = false, std::ostream& out = std::cout)
			: _printEvery(printEvery)
			, _verbose(verbose)
			, _out(out) {
			if (printEvery <= 0) {
				throw std::invalid_argument("ConsoleObserver: printEvery must be positive");
			}
		}

		void OnStart(const OptimizationState& state) override {
			_out << "\n=== Optimization Started ===" << std::endl;
			_out << "Initial f = " << state.fCurrent << std::endl;
			if (_verbose) {
				_out << "Dimension = " << state.xCurrent.size() << std::endl;
			}
			_out << std::endl;

			// Print header
			_out << std::setw(8) << "Iter" << std::setw(15) << "fBest" << std::setw(12) << "fEvals";

			if (_verbose) {
				_out << std::setw(12) << "Accepted" << std::setw(12) << "AcceptRate" << std::setw(12) << "Time(s)";
			}
			_out << std::endl;
			_out << std::string(80, '-') << std::endl;
		}

		bool OnIteration(const OptimizationState& state) override {
			if (state.iteration % _printEvery == 0) {
				_out << std::setw(8) << state.iteration << std::setw(15) << std::scientific << std::setprecision(6) << state.fBest
					 << std::setw(12) << state.funcEvals;

				if (_verbose && state.stochastic) {
					_out << std::setw(12) << state.stochastic->acceptedMoves << std::setw(12) << std::fixed << std::setprecision(3)
						 << state.stochastic->GetAcceptanceRate() << std::setw(12) << std::fixed << std::setprecision(2) << state.elapsedTime;
				} else if (_verbose) {
					_out << std::setw(12) << "-" << std::setw(12) << "-" << std::setw(12) << std::fixed << std::setprecision(2) << state.elapsedTime;
				}
				_out << std::endl;
			}
			return true; // Continue optimization
		}

		void OnComplete(const OptimizationState& state, const std::string& reason) override {
			_out << std::string(80, '-') << std::endl;
			_out << "\n=== Optimization Complete ===" << std::endl;
			_out << "Reason: " << reason << std::endl;
			_out << "Final fBest = " << state.fBest << std::endl;
			_out << "Iterations: " << state.iteration << std::endl;
			_out << "Function evaluations: " << state.funcEvals << std::endl;

			if (_verbose) {
				_out << "Elapsed time: " << state.elapsedTime << " seconds" << std::endl;
				if (state.stochastic) {
					_out << "Acceptance rate: " << state.stochastic->GetAcceptanceRate() << std::endl;
				}
				if (state.gradient && state.gradient->gradNorm > 0) {
					_out << "Final gradient norm: " << state.gradient->gradNorm << std::endl;
				}
			}
			_out << std::endl;
		}

		void SetPrintEvery(int printEvery) {
			if (printEvery <= 0) {
				throw std::invalid_argument("printEvery must be positive");
			}
			_printEvery = printEvery;
		}
		void SetVerbose(bool verbose) { _verbose = verbose; }
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// TRAJECTORY OBSERVER
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Saves OptimizationState snapshots for post-analysis
	/// @details Allows trajectory visualization, convergence analysis, etc.
	class TrajectoryObserver : public IOptimizationObserver {
	private:
		std::vector<OptimizationState> _trajectory;
		int _saveEvery;	 ///< Save every N iterations (1 = every iteration)
		size_t _maxSize; ///< Maximum trajectory size (0 = unlimited)

	public:
		/// @brief Constructor
		/// @param saveEvery Save state every N iterations (default: 1)
		/// @param maxSize Maximum trajectory size, 0 for unlimited (default: 0)
		explicit TrajectoryObserver(int saveEvery = 1, size_t maxSize = 0)
			: _saveEvery(saveEvery)
			, _maxSize(maxSize) {
			if (saveEvery <= 0) {
				throw std::invalid_argument("TrajectoryObserver: saveEvery must be positive");
			}
		}

		void OnStart(const OptimizationState& state) override {
			_trajectory.clear();
			_trajectory.push_back(state); // Save initial state
		}

		bool OnIteration(const OptimizationState& state) override {
			if (state.iteration % _saveEvery == 0) {
				// Check size limit
				if (_maxSize > 0 && _trajectory.size() >= _maxSize) {
					// Keep first and last, downsample middle
					DownsampleTrajectory();
				}
				_trajectory.push_back(state);
			}
			return true; // Continue
		}

		void OnComplete(const OptimizationState& state, const std::string&) override {
			// Always save final state if not already saved
			if (_trajectory.empty() || _trajectory.back().iteration != state.iteration) {
				_trajectory.push_back(state);
			}
		}

		/// @brief Get the recorded trajectory
		const std::vector<OptimizationState>& GetTrajectory() const { return _trajectory; }

		/// @brief Get trajectory size
		size_t GetSize() const { return _trajectory.size(); }

		/// @brief Clear trajectory (free memory)
		void Clear() { _trajectory.clear(); }

		/// @brief Export trajectory to CSV file
		/// @param filename Output file path
		/// @param includeX Include solution vector components (default: false, can be large)
		void ExportCSV(const std::string& filename, bool includeX = false) const {
			std::ofstream file(filename);
			if (!file.is_open()) {
				throw std::runtime_error("TrajectoryObserver: Cannot open file " + filename);
			}

			// Header - core fields always, optional diagnostics if any state has them
			file << "iteration,funcEvals,fCurrent,fBest,iterSinceImprovement,elapsedTime";
			
			// Check if any state has optional diagnostics
			bool hasGradient = false, hasSimplex = false, hasStochastic = false;
			for (const auto& s : _trajectory) {
				if (s.gradient) hasGradient = true;
				if (s.simplex) hasSimplex = true;
				if (s.stochastic) hasStochastic = true;
			}
			
			if (hasGradient) file << ",gradNorm,stepSize";
			if (hasSimplex) file << ",simplexSize";
			if (hasStochastic) file << ",temperature,acceptedMoves,rejectedMoves,acceptRate";

			if (includeX && !_trajectory.empty()) {
				size_t dim = _trajectory[0].xBest.size();
				for (size_t i = 0; i < dim; ++i) {
					file << ",x" << i;
				}
			}
			file << std::endl;

			// Data
			for (const auto& state : _trajectory) {
				file << state.iteration << "," << state.funcEvals << "," << state.fCurrent << "," << state.fBest
				     << "," << state.iterSinceImprovement << "," << state.elapsedTime;

				if (hasGradient) {
					if (state.gradient) {
						file << "," << state.gradient->gradNorm << "," << state.gradient->stepSize;
					} else {
						file << ",,";
					}
				}
				if (hasSimplex) {
					if (state.simplex) {
						file << "," << state.simplex->simplexSize;
					} else {
						file << ",";
					}
				}
				if (hasStochastic) {
					if (state.stochastic) {
						file << "," << state.stochastic->temperature << "," << state.stochastic->acceptedMoves
						     << "," << state.stochastic->rejectedMoves << "," << state.stochastic->GetAcceptanceRate();
					} else {
						file << ",,,,";
					}
				}

				if (includeX) {
					for (int i = 0; i < state.xBest.size(); ++i) {
						file << "," << state.xBest[i];
					}
				}
				file << std::endl;
			}

			file.close();
		}

		/// @brief Export convergence plot data (iteration vs fBest)
		void ExportConvergencePlot(const std::string& filename) const {
			std::ofstream file(filename);
			if (!file.is_open()) {
				throw std::runtime_error("TrajectoryObserver: Cannot open file " + filename);
			}

			file << "iteration,fBest" << std::endl;
			for (const auto& state : _trajectory) {
				file << state.iteration << "," << state.fBest << std::endl;
			}

			file.close();
		}

	private:
		/// @brief Downsample trajectory when size limit reached
		void DownsampleTrajectory() {
			if (_trajectory.size() <= 2)
				return;

			// Keep every other element (simple decimation)
			std::vector<OptimizationState> downsampled;
			downsampled.reserve(_trajectory.size() / 2 + 1);

			for (size_t i = 0; i < _trajectory.size(); i += 2) {
				downsampled.push_back(_trajectory[i]);
			}

			_trajectory = std::move(downsampled);
			_saveEvery *= 2; // Adjust save frequency
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// CALLBACK OBSERVER
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Lambda-friendly observer for custom callbacks
	/// @details Allows inline custom logic without creating new observer classes
	class CallbackObserver : public IOptimizationObserver {
	private:
		std::function<void(const OptimizationState&)> _onStartCallback;
		std::function<bool(const OptimizationState&)> _onIterCallback;
		std::function<void(const OptimizationState&, const std::string&)> _onCompleteCallback;

	public:
		CallbackObserver() = default;

		/// @brief Set callback for optimization start
		void SetOnStart(std::function<void(const OptimizationState&)> callback) { _onStartCallback = callback; }

		/// @brief Set callback for each iteration
		/// @note Callback should return true to continue, false to abort
		void SetOnIteration(std::function<bool(const OptimizationState&)> callback) { _onIterCallback = callback; }

		/// @brief Set callback for optimization completion
		void SetOnComplete(std::function<void(const OptimizationState&, const std::string&)> callback) {
			_onCompleteCallback = callback;
		}

		void OnStart(const OptimizationState& state) override {
			if (_onStartCallback) {
				_onStartCallback(state);
			}
		}

		bool OnIteration(const OptimizationState& state) override {
			if (_onIterCallback) {
				return _onIterCallback(state);
			}
			return true; // Continue by default
		}

		void OnComplete(const OptimizationState& state, const std::string& reason) override {
			if (_onCompleteCallback) {
				_onCompleteCallback(state, reason);
			}
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// MULTI OBSERVER (Composite)
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Combines multiple observers into one
	/// @details Calls all observers in sequence
	class MultiObserver : public IOptimizationObserver {
	private:
		std::vector<std::shared_ptr<IOptimizationObserver>> _observers;

	public:
		MultiObserver() = default;

		/// @brief Add an observer to the collection
		void AddObserver(std::shared_ptr<IOptimizationObserver> observer) {
			if (!observer) {
				throw std::invalid_argument("MultiObserver: Cannot add null observer");
			}
			_observers.push_back(observer);
		}

		/// @brief Add observer by constructing in-place
		template<typename ObserverType, typename... Args>
		void AddObserver(Args&&... args) {
			_observers.push_back(std::make_shared<ObserverType>(std::forward<Args>(args)...));
		}

		void OnStart(const OptimizationState& state) override {
			for (auto& obs : _observers) {
				obs->OnStart(state);
			}
		}

		bool OnIteration(const OptimizationState& state) override {
			// If any observer returns false, abort
			for (auto& obs : _observers) {
				if (!obs->OnIteration(state)) {
					return false;
				}
			}
			return true;
		}

		void OnComplete(const OptimizationState& state, const std::string& reason) override {
			for (auto& obs : _observers) {
				obs->OnComplete(state, reason);
			}
		}

		size_t GetObserverCount() const { return _observers.size(); }
		void Clear() { _observers.clear(); }
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// EXAMPLE CALLBACK PATTERNS
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Example: Checkpointing callback (saves state every N iterations)
	inline std::function<bool(const OptimizationState&)> MakeCheckpointCallback(int saveEvery, const std::string& baseFilename) {
		return [saveEvery, baseFilename](const OptimizationState& state) {
			if (state.iteration % saveEvery == 0) {
				std::string filename = baseFilename + "_iter" + std::to_string(state.iteration) + ".mml";
				std::ofstream file(filename);
				if (file.is_open()) {
					file << "iteration: " << state.iteration << std::endl;
					file << "fBest: " << state.fBest << std::endl;
					file << "xBest:";
					for (int i = 0; i < state.xBest.size(); ++i) {
						file << " " << state.xBest[i];
					}
					file << std::endl;
					file.close();
				}
			}
			return true; // Continue
		};
	}

	/// @brief Example: Early stopping if function value too good (sanity check)
	inline std::function<bool(const OptimizationState&)> MakeEarlyStoppingCallback(Real stopIfBelow) {
		return [stopIfBelow](const OptimizationState& state) {
			return state.fBest >= stopIfBelow; // Stop if too good
		};
	}

	/// @brief Example: Alert if stagnation detected
	inline std::function<bool(const OptimizationState&)> MakeStagnationAlertCallback(int alertAfter) {
		return [alertAfter](const OptimizationState& state) {
			if (state.iterSinceImprovement == alertAfter) {
				std::cout << "⚠ ALERT: No improvement for " << alertAfter << " iterations!" << std::endl;
			}
			return true; // Continue
		};
	}

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_OBSERVERS_H

