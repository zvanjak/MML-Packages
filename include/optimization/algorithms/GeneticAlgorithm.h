///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GeneticAlgorithm.h                                                  ///
///  Description: Genetic Algorithm for real-coded mixed-integer optimization         ///
///               Single-objective optimization with continuous, integer, and         ///
///               categorical variables                                               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GENETIC_ALGORITHM_H
#define MML_GENETIC_ALGORITHM_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Vector/Vector.h"

// Optimization framework integration
#include "OptimizationCommon.h"
#include "OptimizationConfig.h"
#include "IterativeRunner.h"
#include "IOptimizationProblem.h"
#include "mml/mml_export.h"

#include <concepts>
#include <random>
#include <functional>
#include <cmath>
#include <limits>
#include <chrono>
#include <algorithm>
#include <memory>
#include <string>
#include <type_traits>

#include "ga/Population.h"
#include "ga/GAConfig.h"
#include "ga/Selection.h"
#include "ga/Crossover.h"
#include "ga/Mutation.h"

namespace MML::Optimization 
{
	///////////////////////////////////////////////////////////////////////////
	///                     Genetic Algorithm                               ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Real-coded Genetic Algorithm for single-objective optimization
     * 
     * Features:
     * - Mixed-integer variables (continuous, integer, categorical)
     * - Pluggable selection strategies
     * - Pluggable crossover operators
     * - Pluggable mutation operators
     * - Elitism
     * - Multiple termination criteria
     * 
     * Usage:
     * @code
	* ProblemSpec spec = ProblemSpec::Continuous(10, -5.12, 5.12);
     * GeneticAlgorithm ga;
     * ga.SetProblem(spec);
     * GAResult result = ga.Minimize(objectiveFunction);
     * @endcode
     */
	class MML_OPTIMIZATION_API GeneticAlgorithm : public IIterativeAlgorithm {
	private:
		GAConfig _config;
		ProblemSpec _problemSpec;

		// Operators (owned)
		std::unique_ptr<ISelectionStrategy> _selection;
		std::unique_ptr<ICrossoverOperator> _crossover;
		std::unique_ptr<IMutationOperator> _mutation;

		// State
		Population _population;
		int _funcEvals;
		Real _bestFitness;
		bool _minimize = true;
		Vector<Real> _bestSolution;

		// RNG
		std::mt19937 _rng;
		std::function<Real(const Vector<Real>&)> _objective;
		OptimizationConfig* _runConfig = nullptr;
		OptimizationState _state;
		ScalarIterationProgress _progress;
		std::unique_ptr<IterativeRunner<ScalarIterationProgress>> _runner;
		GAResult _result;
		int _stagnationCount = 0;
		Real _previousBest = std::numeric_limits<Real>::max();

	public:
		/// @brief Default constructor with sensible defaults
		GeneticAlgorithm();

		/// @brief Construct with configuration
		explicit GeneticAlgorithm(const GAConfig& config);

		/// @brief Set the problem specification
		void SetProblem(const ProblemSpec& spec);

		/// @brief Set configuration
		void SetConfig(const GAConfig& config);

		/// @brief Get current configuration
		const GAConfig& GetConfig() const;

		/// @brief Set selection strategy
		void SetSelection(std::unique_ptr<ISelectionStrategy> selection);

		/// @brief Set crossover operator
		void SetCrossover(std::unique_ptr<ICrossoverOperator> crossover);

		/// @brief Set mutation operator
		void SetMutation(std::unique_ptr<IMutationOperator> mutation);

		/// @brief Access current population
		const Population& GetPopulation() const;

		/// @brief Get best fitness found
		Real GetBestFitness() const;

		/// @brief Get best solution found
		const Vector<Real>& GetBestSolution() const;

		/// @brief Minimize a typed scalar problem using its variable specification
		GAResult Minimize(const ISingleObjectiveProblem& problem);

		/// @brief Maximize a typed scalar problem using its variable specification
		GAResult Maximize(const ISingleObjectiveProblem& problem);

		/// @brief Optimize a typed scalar problem in the requested direction
		GAResult Optimize(const ISingleObjectiveProblem& problem, bool minimize = true);

		void Start(const ISingleObjectiveProblem& problem, bool minimize = true, OptimizationConfig* config = nullptr);
		template<typename Func> requires (std::invocable<Func&, const Vector<Real>&> &&
			!std::derived_from<std::remove_cvref_t<Func>, ISingleObjectiveProblem>)
		void Start(Func& func, bool minimize = true, OptimizationConfig* config = nullptr) {
			if (_progress.status == IterativeRunStatus::Running)
				throw std::logic_error("Cannot rebind an active GA run");
			if (_problemSpec.Dimension() == 0)
				throw GeneticAlgorithmError("Problem specification not set");
			_objective = [&func](const Vector<Real>& genes) { return func(genes); };
			_minimize = minimize;
			_runConfig = config;
			Initialize();
		}
		void Initialize() override;
		void Step() override;
		const ScalarIterationProgress& Progress() const override;
		void Run() override;

		/**
         * @brief Minimize an objective function
         * @tparam Func Callable with signature Real(const Vector<Real>&)
         * @param func Objective function to minimize
         * @return Optimization result
         */
		template<typename Func> requires (std::invocable<Func&, const Vector<Real>&> &&
			!std::derived_from<std::remove_cvref_t<Func>, ISingleObjectiveProblem>)
		GAResult Minimize(Func& func) {
			return Optimize(func, true);
		}

		/**
         * @brief Maximize an objective function
         * @tparam Func Callable with signature Real(const Vector<Real>&)
         * @param func Objective function to maximize
         * @return Optimization result
         */
		template<typename Func> requires (std::invocable<Func&, const Vector<Real>&> &&
			!std::derived_from<std::remove_cvref_t<Func>, ISingleObjectiveProblem>)
		GAResult Maximize(Func& func) {
			return Optimize(func, false);
		}

		/**
         * @brief Run optimization
         * @tparam Func Callable with signature Real(const Vector<Real>&)
         * @param func Objective function
         * @param minimize True for minimization, false for maximization
         * @return Optimization result
         */
		template<typename Func> requires (std::invocable<Func&, const Vector<Real>&> &&
			!std::derived_from<std::remove_cvref_t<Func>, ISingleObjectiveProblem>)
		GAResult Optimize(Func& func, bool minimize = true) {
			Start(func, minimize);
			Run();
			return Result();
		}

	private:
		GAResult Result() const;
		/// @brief Initialize random population
		void InitializePopulation();

		/// @brief Evaluate population fitness
		template<typename Func>
		void EvaluatePopulation(Func& func, bool minimize, int skipFirst = 0) {
			for (int i = skipFirst; i < _population.Size(); ++i) {
				if (!_population[i].evaluated) {
					Real fitness = func(_population[i].genes);
					if (!minimize)
						fitness = -fitness; // Convert to minimization
					_population[i].fitness = fitness;
					_population[i].evaluated = true;
					_funcEvals++;
				}
			}
		}

		/// @brief Update best solution
		void UpdateBest(bool minimize);
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Convenience Functions                           ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Simple GA minimization with automatic configuration
     */
	template<typename Func>
	GAResult GAMinimize(Func& func, const ProblemSpec& spec, int populationSize = 100, int maxGenerations = 1000) {
		GAConfig config;
		config.populationSize = populationSize;
		config.maxGenerations = maxGenerations;

		GeneticAlgorithm ga(config);
		ga.SetProblem(spec);
		return ga.Minimize(func);
	}

	/**
     * @brief GA minimization with bounds (all continuous variables)
     */
	template<typename Func>
	GAResult GAMinimize(Func& func, const Vector<Real>& lowerBounds, const Vector<Real>& upperBounds, int populationSize = 100,
											int maxGenerations = 1000) {
		auto spec = ProblemSpec::Continuous(lowerBounds, upperBounds);
		return GAMinimize(func, spec, populationSize, maxGenerations);
	}
} // namespace MML::Optimization
#endif // MML_GENETIC_ALGORITHM_H
