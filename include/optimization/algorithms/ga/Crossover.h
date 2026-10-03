#if !defined MML_OPTIMIZATION_GA_CROSSOVER_H
#define MML_OPTIMIZATION_GA_CROSSOVER_H

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
	///                     Crossover Operator Interface                    ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Interface for crossover operators
     */
	class ICrossoverOperator {
	public:
		virtual ~ICrossoverOperator() = default;

		/**
         * @brief Create offspring from two parents
         * @param parent1 First parent chromosome
         * @param parent2 Second parent chromosome
         * @param spec Problem specification (for type-aware crossover)
         * @return Pair of offspring (child1, child2)
         */
		virtual std::pair<Vector<Real>, Vector<Real>> Crossover(const Vector<Real>& parent1, const Vector<Real>& parent2,
																														const ProblemSpec& spec) = 0;

		/**
         * @brief Set random seed
         */
		virtual void SetSeed(unsigned int seed) = 0;
	};

	///////////////////////////////////////////////////////////////////////////
	///                     BLX-α Crossover (Blend)                         ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Blend Crossover: offspring gene sampled from extended range
     * 
     * For gene i with parents p1[i], p2[i]:
     *   range = |p1[i] - p2[i]|
     *   child[i] ~ Uniform[min - α*range, max + α*range]
     * 
     * α = 0.5 is most common (BLX-0.5)
     * Excellent for continuous optimization
     */
	class MML_OPTIMIZATION_API BLXCrossover : public ICrossoverOperator {
	private:
		Real _alpha;
		std::mt19937 _rng;

	public:
		/**
         * @brief Construct BLX-α crossover
         * @param alpha Extension factor (0.5 is typical)
         * @param seed Random seed
         */
		explicit BLXCrossover(Real alpha = 0.5, unsigned int seed = 0);

		void SetSeed(unsigned int seed) override;

		std::pair<Vector<Real>, Vector<Real>> Crossover(const Vector<Real>& parent1, const Vector<Real>& parent2,
																		const ProblemSpec& spec) override;

		Real GetAlpha() const;
		void SetAlpha(Real alpha);
	};

	///////////////////////////////////////////////////////////////////////////
	///                     SBX Crossover (Simulated Binary)                ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Simulated Binary Crossover: mimics single-point crossover behavior
     * 
     * Uses polynomial probability distribution controlled by η_c:
     * - η_c = 2: High spread (more exploration)
     * - η_c = 5: Moderate spread
     * - η_c = 20: Low spread (children close to parents)
     * 
     * Widely used in NSGA-II and other modern GAs
     */
	class MML_OPTIMIZATION_API SBXCrossover : public ICrossoverOperator {
	private:
		Real _eta; ///< Distribution index (typically 2-20)
		std::mt19937 _rng;

	public:
		/**
         * @brief Construct SBX crossover
         * @param eta Distribution index (higher = children closer to parents)
         * @param seed Random seed
         */
		explicit SBXCrossover(Real eta = 15.0, unsigned int seed = 0);

		void SetSeed(unsigned int seed) override;

		std::pair<Vector<Real>, Vector<Real>> Crossover(const Vector<Real>& parent1, const Vector<Real>& parent2,
																		const ProblemSpec& spec) override;

		Real GetEta() const;
		void SetEta(Real eta);
	};

	///////////////////////////////////////////////////////////////////////////
	///                     Uniform Crossover                               ///
	///////////////////////////////////////////////////////////////////////////

	/**
     * @brief Uniform Crossover: each gene independently from either parent
     * 
     * For each gene i:
     *   child1[i] = (random < 0.5) ? parent1[i] : parent2[i]
     *   child2[i] = (random < 0.5) ? parent2[i] : parent1[i]
     * 
     * Good for problems where genes are independent
     */
	class MML_OPTIMIZATION_API UniformCrossover : public ICrossoverOperator {
	private:
		Real _swapProb;
		std::mt19937 _rng;

	public:
		/**
         * @brief Construct uniform crossover
         * @param swapProb Probability of swapping each gene (default 0.5)
         * @param seed Random seed
         */
		explicit UniformCrossover(Real swapProb = 0.5, unsigned int seed = 0);

		void SetSeed(unsigned int seed) override;

		std::pair<Vector<Real>, Vector<Real>> Crossover(const Vector<Real>& parent1, const Vector<Real>& parent2,
																		const ProblemSpec& spec) override;
	};
} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_GA_CROSSOVER_H
