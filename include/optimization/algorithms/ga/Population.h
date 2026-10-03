#if !defined MML_OPTIMIZATION_GA_POPULATION_H
#define MML_OPTIMIZATION_GA_POPULATION_H

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

namespace MML::Optimization 
{
	///////////////////////////////////////////////////////////////////////////
	///                         Individual                                  ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief A single individual (chromosome) in the GA population
     * 
     * Contains the solution vector (genes), fitness value, and metadata.
     */
	struct Individual {
		Vector<Real> genes;				///< Solution vector
		Real fitness;			///< Objective value (lower = better for minimization)
		Real constraintViolation; ///< Total constraint violation (0 = feasible)
		int age;									///< Generations since creation
		bool evaluated;						///< Has fitness been computed?

		/// @brief Default constructor
		Individual()
				: fitness(std::numeric_limits<Real>::max())
				, constraintViolation(0)
				, age(0)
				, evaluated(false) {}

		/// @brief Construct with dimension
		explicit Individual(int n)
				: genes(n)
				, fitness(std::numeric_limits<Real>::max())
				, constraintViolation(0)
				, age(0)
				, evaluated(false) {}

		/// @brief Construct from gene vector
		explicit Individual(const Vector<Real>& x)
				: genes(x)
				, fitness(std::numeric_limits<Real>::max())
				, constraintViolation(0)
				, age(0)
				, evaluated(false) {}

		/// @brief Check if this individual is feasible
		bool IsFeasible() const { return constraintViolation <= 0; }

		/// @brief Comparison by fitness (for sorting, minimization)
		bool operator<(const Individual& other) const { return fitness < other.fitness; }

		/// @brief Get dimension
		int Dimension() const { return genes.size(); }
	};

	///////////////////////////////////////////////////////////////////////////
	///                         Population                                  ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief A population of individuals for genetic algorithm
     * 
     * Manages a collection of individuals with statistics and operations.
     */
	class Population {
	private:
		std::vector<Individual> _individuals;
		int _generation;

	public:
		/// @brief Default constructor
		Population()
				: _generation(0) {}

		/// @brief Construct with initial capacity
		explicit Population(int size)
				: _generation(0) {
			_individuals.reserve(size);
		}

		/// @brief Get population size
		int Size() const { return static_cast<int>(_individuals.size()); }

		/// @brief Get current generation
		int Generation() const { return _generation; }

		/// @brief Increment generation counter
		void NextGeneration() { ++_generation; }

		/// @brief Access individual by index
		Individual& operator[](int i) { return _individuals[i]; }
		const Individual& operator[](int i) const { return _individuals[i]; }

		/// @brief Get best individual (lowest fitness for minimization)
		const Individual& Best() const {
			if (_individuals.empty())
				throw GeneticAlgorithmError("Cannot get best from empty population");

			return *std::min_element(_individuals.begin(), _individuals.end());
		}

		/// @brief Get worst individual (highest fitness for minimization)
		const Individual& Worst() const {
			if (_individuals.empty())
				throw GeneticAlgorithmError("Cannot get worst from empty population");

			return *std::max_element(_individuals.begin(), _individuals.end());
		}

		/// @brief Get best fitness value
		Real BestFitness() const { return Best().fitness; }

		/// @brief Get worst fitness value
		Real WorstFitness() const { return Worst().fitness; }

		/// @brief Calculate average fitness
		Real AverageFitness() const {
			if (_individuals.empty())
				return Real(0);

			Real sum = 0;
			for (const auto& ind : _individuals)
				sum += ind.fitness;
			return sum / static_cast<Real>(_individuals.size());
		}

		/// @brief Calculate fitness standard deviation
		Real FitnessStdDev() const {
			if (_individuals.size() < 2)
				return Real(0);

			Real avg = AverageFitness();
			Real sumSq = 0;
			for (const auto& ind : _individuals) {
				Real diff = ind.fitness - avg;
				sumSq += diff * diff;
			}
			return std::sqrt(sumSq / static_cast<Real>(_individuals.size() - 1));
		}

		/// @brief Sort population by fitness (ascending)
		void Sort() { std::sort(_individuals.begin(), _individuals.end()); }

		/// @brief Add an individual to the population
		void Add(const Individual& ind) { _individuals.push_back(ind); }

		/// @brief Add an individual (move)
		void Add(Individual&& ind) { _individuals.push_back(std::move(ind)); }

		/// @brief Clear the population
		void Clear() { _individuals.clear(); }

		/// @brief Reserve capacity
		void Reserve(int capacity) { _individuals.reserve(capacity); }

		/// @brief Resize population
		void Resize(int newSize) { _individuals.resize(newSize); }

		/// @brief Replace population with new individuals
		void Replace(std::vector<Individual>&& newIndividuals) { _individuals = std::move(newIndividuals); }

		/// @brief Iterator access
		auto begin() { return _individuals.begin(); }
		auto end() { return _individuals.end(); }
		auto begin() const { return _individuals.begin(); }
		auto end() const { return _individuals.end(); }

		/// @brief Get underlying vector
		std::vector<Individual>& Individuals() { return _individuals; }
		const std::vector<Individual>& Individuals() const { return _individuals; }
	};
} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_GA_POPULATION_H
