#if !defined MML_OPTIMIZATION_GA_SELECTION_H
#define MML_OPTIMIZATION_GA_SELECTION_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/Vector/Vector.h"

// Optimization framework integration
#include "OptimizationCommon.h"
#include "OptimizationConfig.h"
#include "mml/mml_export.h"

#include <random>
#include <functional>
#include <cmath>
#include <limits>
#include <chrono>
#include <algorithm>
#include <memory>
#include <string>

#include "GAConfig.h"
#include "Population.h"

namespace MML::Optimization 
{
	///////////////////////////////////////////////////////////////////////////
	///                   Selection Strategy Interface                      ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Interface for parent selection strategies
     */
	class ISelectionStrategy {
	public:
		virtual ~ISelectionStrategy() = default;

		/**
         * @brief Select one individual from population
         * @param pop Population to select from
         * @return Index of selected individual
         */
		virtual int Select(const Population& pop) = 0;

		/**
         * @brief Select two parents for crossover
         * @param pop Population to select from
         * @return Pair of indices (parent1, parent2)
         */
		virtual std::pair<int, int> SelectParents(const Population& pop) {
			int p1 = Select(pop);
			int p2 = Select(pop);
			// Ensure different parents
			while (p2 == p1 && pop.Size() > 1)
				p2 = Select(pop);
			return {p1, p2};
		}

		/**
         * @brief Set random seed for reproducibility
         */
		virtual void SetSeed(unsigned int seed) = 0;
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Tournament Selection                            ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Tournament Selection: Pick k random individuals, return best
     * 
     * - k=2 (binary tournament): Moderate selection pressure
     * - k=3-5: Higher selection pressure
     * - Most commonly used in practice
     */
	class TournamentSelection : public ISelectionStrategy {
	private:
		int _tournamentSize;
		std::mt19937 _rng;

	public:
		/**
         * @brief Construct tournament selection
         * @param tournamentSize Number of individuals in each tournament (default 3)
         * @param seed Random seed (0 = time-based)
         */
		explicit TournamentSelection(int tournamentSize = 3, unsigned int seed = 0)
				: _tournamentSize(tournamentSize) {
			if (tournamentSize < 2)
				throw GeneticAlgorithmError("Tournament size must be at least 2");
			SetSeed(seed);
		}

		void SetSeed(unsigned int seed) override {
			if (seed == 0)
				seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
			_rng.seed(seed);
		}

		int Select(const Population& pop) override {
			if (pop.Size() == 0)
				throw GeneticAlgorithmError("Cannot select from empty population");

			std::uniform_int_distribution<int> dist(0, pop.Size() - 1);

			int bestIdx = dist(_rng);
			Real bestFit = pop[bestIdx].fitness;

			for (int i = 1; i < _tournamentSize; ++i) {
				int idx = dist(_rng);
				if (pop[idx].fitness < bestFit) // Minimization
				{
					bestIdx = idx;
					bestFit = pop[idx].fitness;
				}
			}

			return bestIdx;
		}

		int GetTournamentSize() const { return _tournamentSize; }
		void SetTournamentSize(int size) {
			if (size < 2)
				throw GeneticAlgorithmError("Tournament size must be at least 2");
			_tournamentSize = size;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///               Feasibility-First Tournament Selection                ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Tournament selection using Deb's feasibility rules for constrained optimization
     * 
     * Implements the constraint-handling technique from:
     * Deb, K. (2000). "An efficient constraint handling method for genetic algorithms"
     * 
     * Selection rules (for minimization):
     * 1. Feasible individual beats infeasible individual
     * 2. Among two feasible individuals, better fitness wins
     * 3. Among two infeasible individuals, lower constraint violation wins
     * 
     * This approach doesn't require penalty parameters and naturally guides
     * the population toward the feasible region while optimizing.
     */
	class FeasibilityTournamentSelection : public ISelectionStrategy {
	private:
		int _tournamentSize;
		std::mt19937 _rng;

	public:
		/**
         * @brief Construct feasibility-first tournament selection
         * @param tournamentSize Number of individuals in each tournament (default 2)
         * @param seed Random seed (0 = time-based)
         */
		explicit FeasibilityTournamentSelection(int tournamentSize = 2, unsigned int seed = 0)
				: _tournamentSize(tournamentSize) {
			if (tournamentSize < 2)
				throw GeneticAlgorithmError("Tournament size must be at least 2");
			SetSeed(seed);
		}

		void SetSeed(unsigned int seed) override {
			if (seed == 0)
				seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
			_rng.seed(seed);
		}

		/**
         * @brief Compare two individuals using Deb's feasibility rules
         * @return true if individual a is better than individual b
         */
		bool IsBetter(const Individual& a, const Individual& b) const {
			bool aFeasible = a.IsFeasible();
			bool bFeasible = b.IsFeasible();

			// Rule 1: Feasible beats infeasible
			if (aFeasible && !bFeasible)
				return true;
			if (!aFeasible && bFeasible)
				return false;

			// Rule 2: Both feasible - compare by fitness (minimization)
			if (aFeasible && bFeasible) {
				return a.fitness < b.fitness;
			}

			// Rule 3: Both infeasible - compare by constraint violation
			return a.constraintViolation < b.constraintViolation;
		}

		int Select(const Population& pop) override {
			if (pop.Size() == 0)
				throw GeneticAlgorithmError("Cannot select from empty population");

			std::uniform_int_distribution<int> dist(0, pop.Size() - 1);

			int bestIdx = dist(_rng);

			for (int i = 1; i < _tournamentSize; ++i) {
				int idx = dist(_rng);
				if (IsBetter(pop[idx], pop[bestIdx])) {
					bestIdx = idx;
				}
			}

			return bestIdx;
		}

		int GetTournamentSize() const { return _tournamentSize; }
		void SetTournamentSize(int size) {
			if (size < 2)
				throw GeneticAlgorithmError("Tournament size must be at least 2");
			_tournamentSize = size;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Rank-Based Selection                            ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Selection probability based on rank, not raw fitness
     * 
     * More stable than roulette when fitness values vary widely.
     * Linear ranking with configurable selective pressure.
     */
	class RankSelection : public ISelectionStrategy {
	private:
		Real _selectivePressure; ///< Typically 1.5 to 2.0
		std::mt19937 _rng;
		mutable std::vector<int> _sortedIndices;
		mutable bool _needsSort;

	public:
		/**
         * @brief Construct rank selection
         * @param selectivePressure Pressure parameter (1.0 = uniform, 2.0 = high pressure)
         * @param seed Random seed
         */
		explicit RankSelection(Real selectivePressure = 1.5, unsigned int seed = 0)
				: _selectivePressure(selectivePressure)
				, _needsSort(true) {
			if (selectivePressure < 1.0 || selectivePressure > 2.0)
				throw GeneticAlgorithmError("Selective pressure must be in [1.0, 2.0]");
			SetSeed(seed);
		}

		void SetSeed(unsigned int seed) override {
			if (seed == 0)
				seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
			_rng.seed(seed);
		}

		int Select(const Population& pop) override {
			if (pop.Size() == 0)
				throw GeneticAlgorithmError("Cannot select from empty population");

			int n = pop.Size();

			// Build sorted index array
			_sortedIndices.resize(n);
			for (int i = 0; i < n; ++i)
				_sortedIndices[i] = i;

			std::sort(_sortedIndices.begin(), _sortedIndices.end(), [&pop](int a, int b) { return pop[a].fitness < pop[b].fitness; });

			// Linear ranking probabilities
			// P(rank i) = (2 - s)/n + 2*i*(s-1)/(n*(n-1))
			// where s = selective pressure, i = rank (0 = worst, n-1 = best)
			std::uniform_real_distribution<Real> dist(0.0, 1.0);
			Real r = dist(_rng);

			Real cumProb = 0.0;
			Real s = _selectivePressure;
			for (int i = 0; i < n; ++i) {
				// Rank i (0 = best for minimization)
				Real prob = (2.0 - s) / n + 2.0 * (n - 1 - i) * (s - 1.0) / (n * (n - 1.0));
				cumProb += prob;
				if (r <= cumProb)
					return _sortedIndices[i];
			}

			return _sortedIndices[n - 1]; // Last resort
		}
	};
} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_GA_SELECTION_H
