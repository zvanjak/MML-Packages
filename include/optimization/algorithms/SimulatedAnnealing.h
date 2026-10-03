///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SimulatedAnnealing.h                                                ///
///  Description: Simulated Annealing optimization algorithm                          ///
///               Global optimization for multimodal and non-convex problems          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PKG_SIMULATED_ANNEALING_H
#define MML_PKG_SIMULATED_ANNEALING_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Vector/Vector.h"

#include "IOptimizationProblem.h"
#include "OptimizationCommon.h"
#include "OptimizationConfig.h"
#include "IterativeRunner.h"
#include "mml/mml_export.h"

#include <random>
#include <functional>
#include <cmath>
#include <limits>
#include <chrono>

namespace MML::Optimization {
	/////////////////////////////////////////////////////////////////////
	///                   HEURISTIC OPTIMIZATION                      ///
	/////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	///                   HeuristicOptimizationError                        ///
	///////////////////////////////////////////////////////////////////////////
	class HeuristicOptimizationError : public std::runtime_error {
	public:
		explicit HeuristicOptimizationError(const std::string& message)
			: std::runtime_error("HeuristicOptimizationError: " + message) {}
	};

	///////////////////////////////////////////////////////////////////////////
	///                  HeuristicOptimizationResult                        ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Result structure for heuristic optimization methods
     */
	struct MML_OPTIMIZATION_API HeuristicOptimizationResult {
		Vector<Real> xbest; ///< Best solution found
		Real fbest;			///< Function value at best solution
		int iterations;		///< Number of iterations performed
		int funcEvals;		///< Total function evaluations
		int acceptedMoves;	///< Number of accepted moves (for SA)
		bool converged;		///< True if convergence criterion met

		HeuristicOptimizationResult();

		HeuristicOptimizationResult(const Vector<Real>& x, Real f, int iter, int fEvals, int accepted, bool conv);
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Cooling Schedule Interface                      ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Interface for temperature cooling schedules
     */
	class ICoolingSchedule {
	public:
		virtual ~ICoolingSchedule() = default;

		/**
         * @brief Get the temperature at a given iteration
         * @param iteration Current iteration number (0-based)
         * @param T0 Initial temperature
         * @return Temperature at this iteration
         */
		virtual Real Temperature(int iteration, Real T0) const = 0;
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Exponential Cooling Schedule                    ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Exponential (geometric) cooling: T(k) = T0 * alpha^k
     * 
     * This is the most common cooling schedule. Temperature decreases
     * geometrically with each iteration.
     */
	class MML_OPTIMIZATION_API ExponentialCooling : public ICoolingSchedule {
	private:
		Real _alpha; ///< Cooling rate (0 < alpha < 1, typically 0.9-0.99)

	public:
		/**
         * @brief Construct exponential cooling schedule
         * @param alpha Cooling rate (default 0.95)
         */
		explicit ExponentialCooling(Real alpha = 0.95);

		Real Temperature(int iteration, Real T0) const override;

		Real getAlpha() const;
		void setAlpha(Real alpha);
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Linear Cooling Schedule                         ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Linear cooling: T(k) = T0 * (1 - k/maxIter)
     * 
     * Temperature decreases linearly to zero at maxIter.
     */
	class MML_OPTIMIZATION_API LinearCooling : public ICoolingSchedule {
	private:
		int _maxIter; ///< Maximum iterations (temperature reaches 0)

	public:
		/**
         * @brief Construct linear cooling schedule
         * @param maxIter Maximum iterations
         */
		explicit LinearCooling(int maxIter);

		Real Temperature(int iteration, Real T0) const override;

		int getMaxIter() const;
		void setMaxIter(int maxIter);
	};

	///////////////////////////////////////////////////////////////////////////
	///                   Logarithmic Cooling Schedule                      ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Logarithmic (slow) cooling: T(k) = T0 / (1 + c * ln(1 + k))
     * 
     * Very slow cooling that theoretically guarantees convergence to
     * the global minimum, but may be impractically slow.L
     */
	class MML_OPTIMIZATION_API LogarithmicCooling : public ICoolingSchedule {
	private:
		Real _c; ///< Cooling constant

	public:
		/**
         * @brief Construct logarithmic cooling schedule
         * @param c Cooling constant (default 1.0)
         */
		explicit LogarithmicCooling(Real c = 1.0);

		Real Temperature(int iteration, Real T0) const override;

		Real getC() const;
		void setC(Real c);
	};

	///////////////////////////////////////////////////////////////////////////
	///                      Adaptive Cooling Schedule                      ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Adaptive cooling that adjusts based on acceptance rate
     * 
     * Slows cooling when acceptance rate is low (stuck),
     * speeds up when acceptance rate is high (too hot).
     */
	class MML_OPTIMIZATION_API AdaptiveCooling : public ICoolingSchedule {
	private:
		Real _alphaFast;		///< Fast cooling rate
		Real _alphaSlow;		///< Slow cooling rate
		Real _targetAcceptRate; ///< Target acceptance rate

		mutable Real _currentAlpha;
		mutable Real _currentTemp;
		mutable int _lastIter;

	public:
		/**
         * @brief Construct adaptive cooling schedule
         * @param alphaFast Fast cooling rate (default 0.99)
         * @param alphaSlow Slow cooling rate (default 0.8)
         * @param targetAcceptRate Target acceptance rate (default 0.3)
         */
		AdaptiveCooling(Real alphaFast = 0.99, Real alphaSlow = 0.8, Real targetAcceptRate = 0.3);

		Real Temperature(int iteration, Real T0) const override;

		/**
         * @brief Update cooling rate based on actual acceptance rate
         * @param acceptRate Current acceptance rate
         */
		void UpdateRate(Real acceptRate);

		void reset();
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Neighbor Generator Interface                    ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Interface for generating neighbor solutions
     */
	class INeighborGenerator {
	public:
		virtual ~INeighborGenerator() = default;

		/**
         * @brief Generate a neighbor of the current solution
         * @param current Current solution
         * @param temperature Current temperature (may affect step size)
         * @return Neighbor solution
         */
		virtual Vector<Real> Generate(const Vector<Real>& current, Real temperature) = 0;
	};

	///////////////////////////////////////////////////////////////////////////
	///                    Uniform Random Neighbor Generator                ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Generate neighbors by uniform random perturbation
     * 
     * Each component is perturbed by a uniform random value in [-delta, delta]
     * where delta can optionally scale with temperature.
     */
	class MML_OPTIMIZATION_API UniformNeighborGenerator : public INeighborGenerator {
	private:
		Real _delta;		 ///< Maximum perturbation per component
		bool _scaleWithTemp; ///< Whether to scale delta with temperature
		Real _T0;			 ///< Reference temperature for scaling

		std::mt19937 _rng;
		std::uniform_real_distribution<Real> _dist;

		Vector<Real> _lowerBounds; ///< Lower bounds (optional)
		Vector<Real> _upperBounds; ///< Upper bounds (optional)
		bool _hasBounds;

	public:
		/**
         * @brief Construct uniform neighbor generator
         * @param delta Maximum perturbation magnitude
         * @param scaleWithTemp Scale delta proportionally to temperature
         * @param T0 Reference temperature (initial temperature)
         * @param seed Random seed (0 for random)
         */
		UniformNeighborGenerator(Real delta = 1.0, bool scaleWithTemp = false, Real T0 = 1.0, unsigned int seed = 0);

		/**
         * @brief Set bounds for the search space
         */
		void SetBounds(const Vector<Real>& lower, const Vector<Real>& upper);

		void ClearBounds();

		Vector<Real> Generate(const Vector<Real>& current, Real temperature) override;

		Real getDelta() const;
		void setDelta(Real delta);
		void setT0(Real T0);
	};

	///////////////////////////////////////////////////////////////////////////
	///                  Gaussian Neighbor Generator                        ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Generate neighbors by Gaussian perturbation
     * 
     * Each component is perturbed by a Gaussian random value with 
     * standard deviation sigma.
     */
	class MML_OPTIMIZATION_API GaussianNeighborGenerator : public INeighborGenerator {
	private:
		Real _sigma;		 ///< Standard deviation of perturbation
		bool _scaleWithTemp; ///< Whether to scale sigma with temperature
		Real _T0;			 ///< Reference temperature for scaling

		std::mt19937 _rng;
		std::normal_distribution<Real> _dist;

		Vector<Real> _lowerBounds;
		Vector<Real> _upperBounds;
		bool _hasBounds;

	public:
		/**
         * @brief Construct Gaussian neighbor generator
         * @param sigma Standard deviation of perturbation
         * @param scaleWithTemp Scale sigma proportionally to sqrt(temperature)
         * @param T0 Reference temperature
         * @param seed Random seed (0 for random)
         */
		GaussianNeighborGenerator(Real sigma = 1.0, bool scaleWithTemp = false, Real T0 = 1.0, unsigned int seed = 0);

		void SetBounds(const Vector<Real>& lower, const Vector<Real>& upper);

		void ClearBounds();

		Vector<Real> Generate(const Vector<Real>& current, Real temperature) override;

		Real getSigma() const;
		void setSigma(Real sigma);
		void setT0(Real T0);
	};

	///////////////////////////////////////////////////////////////////////////
	///                      Simulated Annealing                            ///
	///////////////////////////////////////////////////////////////////////////
	/**
     * @brief Simulated Annealing optimization algorithm
     * 
     * Simulated Annealing is a probabilistic metaheuristic for global optimization.
     * It mimics the physical process of heating and slowly cooling a material
     * to decrease defects (minimizing energy).
     * 
     * Algorithm:
     * 1. Start with initial solution x and temperature T
     * 2. Generate neighbor x' of current solution
     * 3. If f(x') < f(x), accept x' (downhill move)
     * 4. If f(x') >= f(x), accept with probability exp(-(f(x')-f(x))/T)
     * 5. Reduce temperature according to cooling schedule
     * 6. Repeat until stopping criterion met
     * 
     * The key insight is that at high temperatures, uphill moves are likely
     * accepted (allowing escape from local minima), while at low temperatures
     * the algorithm behaves like hill-climbing.
     * 
     * Reference: Kirkpatrick, Gelatt, Vecchi (1983)
     */
	class MML_OPTIMIZATION_API SimulatedAnnealing : public IIterativeAlgorithm {
	public:
		/// Stopping criteria
		enum class StopCriteria {
			MaxIterations,	///< Stop after maxIter iterations
			MinTemperature, ///< Stop when temperature drops below threshold
			NoImprovement,	///< Stop after N iterations without improvement
			Combined		///< Use all criteria
		};

	private:
		Real _T0;				///< Initial temperature
		Real _Tmin;				///< Minimum temperature (stopping criterion)
		int _maxIter;			///< Maximum iterations
		int _stagnationLimit;	///< Iterations without improvement to stop
		StopCriteria _stopCrit; ///< Stopping criterion to use

		// Components
		std::unique_ptr<ICoolingSchedule> _coolingSchedule;
		std::unique_ptr<INeighborGenerator> _neighborGen;

		// Random number generation for acceptance
		mutable std::mt19937 _rng;
		mutable std::uniform_real_distribution<Real> _acceptDist;

		struct RunState {
			const ISingleObjectiveProblem* problem = nullptr;
			OptimizationConfig* config = nullptr;
			Vector<Real> start;
			OptimizationState state;
			ScalarIterationProgress progress;
			std::unique_ptr<IterativeRunner<ScalarIterationProgress>> runner;
		};
		std::unique_ptr<RunState> _run;

		void ValidateProblem(const ISingleObjectiveProblem& problem, const Vector<Real>& x0) const;

	public:
		/**
         * @brief Construct Simulated Annealing optimizer with default components
         * @param T0 Initial temperature (default 100.0)
         * @param Tmin Minimum temperature (default 1e-8)
         * @param maxIter Maximum iterations (default 10000)
         * @param stagnationLimit Iterations without improvement (default 1000)
         * @param stopCrit Stopping criterion (default Combined)
         * @param seed Random seed (0 for random)
         */
		SimulatedAnnealing(Real T0 = 100.0, Real Tmin = 1e-8, int maxIter = 10000, int stagnationLimit = 1000,
						   StopCriteria stopCrit = StopCriteria::Combined, unsigned int seed = 0);

		/**
         * @brief Construct with custom cooling schedule and neighbor generator
         */
		SimulatedAnnealing(std::unique_ptr<ICoolingSchedule> cooling, std::unique_ptr<INeighborGenerator> neighbor, Real T0 = 100.0,
						   Real Tmin = 1e-8, int maxIter = 10000, int stagnationLimit = 1000,
						   StopCriteria stopCrit = StopCriteria::Combined, unsigned int seed = 0);

		void Start(const ISingleObjectiveProblem& problem, const Vector<Real>& x0, OptimizationConfig* config = nullptr);
		void Initialize() override;
		void Step() override;
		const ScalarIterationProgress& Progress() const override;
		void Run() override;

		/**
		 * @brief Minimize a single-objective problem using simulated annealing
		 * @param problem Problem to minimize
         * @param x0 Initial solution
         * @return Optimization result
         */
		HeuristicOptimizationResult Minimize(const ISingleObjectiveProblem& problem, const Vector<Real>& x0) {
			Start(problem, x0);
			Run();
			const auto& progress = Progress();
			return HeuristicOptimizationResult(*progress.xBest, progress.fBest, progress.iterations, progress.funcEvals,
			                                   progress.stochastic->acceptedMoves, progress.status == IterativeRunStatus::Converged);
		}

		/**
		 * @brief Minimize using optimization configuration
		 * @param problem Problem to minimize
         * @param x0 Initial solution
         * @param config Optimization configuration (criteria, observers, etc.)
         * @return Enhanced OptimizationResult with trajectory if configured
         */
		OptimizationResult Minimize(const ISingleObjectiveProblem& problem, const Vector<Real>& x0, OptimizationConfig& config) {
			Start(problem, x0, &config);
			Run();
			const auto& progress = Progress();
			auto result = CreateResult(_run->state, progress.terminationReason, progress.status == IterativeRunStatus::Converged);
			if (auto trajectory = config.GetTrajectoryObserver())
				result.trajectory = trajectory->GetTrajectory();
			return result;
		}

		// Accessors
		Real getT0() const;
		void setT0(Real T0);

		Real getTmin() const;
		void setTmin(Real Tmin);

		int getMaxIter() const;
		void setMaxIter(int maxIter);

		int getStagnationLimit() const;
		void setStagnationLimit(int limit);

		StopCriteria getStopCriteria() const;
		void setStopCriteria(StopCriteria crit);

		/**
         * @brief Set custom cooling schedule
         */
		void setCoolingSchedule(std::unique_ptr<ICoolingSchedule> schedule);

		/**
         * @brief Set custom neighbor generator
         */
		void setNeighborGenerator(std::unique_ptr<INeighborGenerator> gen);
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Convenience Functions                           ///
	///////////////////////////////////////////////////////////////////////////

	/**
	* @brief Simple simulated annealing minimization
	* @param problem Problem to minimize
     * @param x0 Initial solution
     * @param T0 Initial temperature (default 100)
     * @param maxIter Maximum iterations (default 10000)
     * @return Optimization result
     */
	inline HeuristicOptimizationResult SimulatedAnnealingMinimize(const ISingleObjectiveProblem& problem, const Vector<Real>& x0, Real T0 = 100.0, int maxIter = 10000) {
		SimulatedAnnealing sa(T0, 1e-8, maxIter);
		return sa.Minimize(problem, x0);
	}
} // namespace MML::Optimization
#endif // MML_PKG_SIMULATED_ANNEALING_H