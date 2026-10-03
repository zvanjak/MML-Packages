#if !defined MML_OPTIMIZATION_GA_MUTATION_H
#define MML_OPTIMIZATION_GA_MUTATION_H

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
	///                     Mutation Operator Interface                     ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Interface for mutation operators
     */
	class IMutationOperator {
	public:
		virtual ~IMutationOperator() = default;

		/**
         * @brief Mutate a chromosome in place
         * @param chromosome Solution to mutate
         * @param spec Problem specification
         * @param mutationRate Probability of mutating each gene
         */
		virtual void Mutate(Vector<Real>& chromosome, const ProblemSpec& spec, Real mutationRate) = 0;

		/**
         * @brief Set random seed
         */
		virtual void SetSeed(unsigned int seed) = 0;
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Gaussian Mutation                               ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Add Gaussian noise to continuous variables
     * 
     * For each gene i (if mutating):
     *   x[i] += N(0, σ * range)
     * 
     * σ is proportional to variable range by default.
     */
	class MML_OPTIMIZATION_API GaussianMutation : public IMutationOperator {
	private:
		Real _sigma; ///< Standard deviation as fraction of range
		std::mt19937 _rng;

	public:
		/**
         * @brief Construct Gaussian mutation
         * @param sigma Std dev as fraction of variable range (default 0.1 = 10%)
         * @param seed Random seed
         */
		explicit GaussianMutation(Real sigma = 0.1, unsigned int seed = 0);

		void SetSeed(unsigned int seed) override;

		void Mutate(Vector<Real>& chromosome, const ProblemSpec& spec, Real mutationRate) override;

		Real GetSigma() const;
		void SetSigma(Real sigma);
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Polynomial Mutation                             ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Polynomial mutation (used in NSGA-II)
     * 
     * Uses polynomial distribution controlled by η_m:
     * - η_m = 20: Small perturbations (common)
     * - η_m = 100: Very small perturbations
     * 
     * Bounded: automatically respects variable bounds
     */
	class MML_OPTIMIZATION_API PolynomialMutation : public IMutationOperator {
	private:
		Real _eta; ///< Distribution index
		std::mt19937 _rng;

	public:
		/**
         * @brief Construct polynomial mutation
         * @param eta Distribution index (higher = smaller perturbations)
         * @param seed Random seed
         */
		explicit PolynomialMutation(Real eta = 20.0, unsigned int seed = 0);

		void SetSeed(unsigned int seed) override;

		void Mutate(Vector<Real>& chromosome, const ProblemSpec& spec, Real mutationRate) override;

		Real GetEta() const;
		void SetEta(Real eta);
	};
} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_GA_MUTATION_H
