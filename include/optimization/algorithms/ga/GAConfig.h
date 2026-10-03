#if !defined MML_OPTIMIZATION_GA_CONFIG_H
#define MML_OPTIMIZATION_GA_CONFIG_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Vector/Vector.h"

// Optimization framework integration
#include "OptimizationCommon.h"
#include "OptimizationConfig.h"
#include "Variables.h"
#include "mml/mml_export.h"

#include <random>
#include <functional>
#include <cmath>
#include <limits>
#include <chrono>
#include <algorithm>
#include <memory>
#include <string>

namespace MML::Optimization 
{
	class GeneticAlgorithmError : public std::runtime_error {
	public:
		explicit GeneticAlgorithmError(const std::string& message)
			: std::runtime_error("GeneticAlgorithmError: " + message) {}
	};

	///////////////////////////////////////////////////////////////////////////
	///                         GAConfig                                    ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Configuration for Genetic Algorithm
     */
	struct GAConfig {
		// Population settings
		int populationSize = 100; ///< Number of individuals in population
		int numElites = 2;				///< Best individuals preserved unchanged

		// Operator probabilities
		Real crossoverRate = 0.9; ///< Probability of crossover (per pair)
		Real mutationRate = 0.1;	///< Per-gene mutation probability

		// Termination criteria
		int maxGenerations = 1000;																	 ///< Maximum generations
		int maxFuncEvals = 100000;																	 ///< Maximum function evaluations (0 = unlimited)
		int stagnationLimit = 100;																	 ///< Generations without improvement before stop
		Real targetFitness = -std::numeric_limits<Real>::infinity(); ///< Stop if reached
		Real fitnessThreshold = 1e-10;															 ///< Convergence threshold

		// Behavior
		bool minimize = true;			 ///< Minimize (true) or maximize (false)
		unsigned int seed = 0;		 ///< Random seed (0 = use time-based seed)
		bool verbose = false;			 ///< Print progress to stdout
		bool trackHistory = false; ///< Track fitness history per generation

		/// @brief Get actual seed (use time if seed == 0)
		MML_OPTIMIZATION_API unsigned int GetSeed() const;
	};

	///////////////////////////////////////////////////////////////////////////
	///                         GAResult                                    ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Result of genetic algorithm optimization
     */
	struct GAResult {
		Vector<Real> bestSolution;		 ///< Best solution found
		Real bestFitness;							 ///< Best fitness value
		int generations;							 ///< Number of generations completed
		int funcEvals;								 ///< Total function evaluations
		int stagnationCount;					 ///< Consecutive generations without improvement
		bool converged;								 ///< True if converged to target or stagnated
		std::string terminationReason; ///< Human-readable reason

		// Optional history
		std::vector<Real> fitnessHistory;		 ///< Best fitness per generation
		std::vector<Real> avgFitnessHistory; ///< Average fitness per generation

		/// @brief Default constructor
		MML_OPTIMIZATION_API GAResult();
	};
} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_GA_CONFIG_H
