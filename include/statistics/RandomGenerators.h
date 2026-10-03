///////////////////////////////////////////////////////////////////////////////////////////
// RandomGenerators.h
// 
// Random number generators (deviates) for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// All generators use Mersenne Twister (std::mt19937_64) for high-quality randomness
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_RANDOM_GENERATORS_H
#define MML_RANDOM_GENERATORS_H

#include <mml/algorithms/Statistics.h>

#include <random>

#include "mml/mml_export.h"

namespace MML
{
	namespace Statistics
	{
		/// @brief Exponential random deviate generator
		/// Generates random numbers from an exponential distribution using the
		/// inverse transform method: X = -ln(U)/λ where U ~ Uniform(0,1)
		/// Uses C++20 std::mt19937_64 for random number generation.
		class MML_STATISTICS_API ExponentialDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _lambda;

		public:
			ExponentialDeviate(Real rate, uint64_t seed = std::random_device{}());
			Real generate();
			Real operator()() { return generate(); }
		};

		/// @brief Logistic random deviate generator
		///
		/// Generates random numbers from a logistic distribution using the
		/// inverse CDF method: X = μ + σ·ln(U/(1-U))
		///
		/// Uses C++20 std::mt19937_64 for random number generation.
		class MML_STATISTICS_API LogisticDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _mu, _sigma;

		public:
			LogisticDeviate(Real location, Real scale, uint64_t seed = std::random_device{}());
			Real generate();
			Real operator()() { return generate(); }
		};

		/// @brief Normal (Gaussian) random deviate generator using Box-Muller transform
		///
		/// Box-Muller method generates pairs of independent normal variates from
		/// uniform random numbers. This implementation caches one value for efficiency.
		///
		/// Algorithm: If U1, U2 ~ Uniform(0,1), then
		/// X = √(-2ln(R²)) · V1,  Y = √(-2ln(R²)) · V2
		///
		/// where V1, V2 are uniform on unit circle, R² = V1² + V2²
		class MML_STATISTICS_API NormalDeviateBoxMuller
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _mu, _sigma;
			Real _storedValue;
			bool _hasStored;

		public:
			NormalDeviateBoxMuller(Real mean, Real stddev, uint64_t seed = std::random_device{}());
			Real generate();
			Real operator()() { return generate(); }
		};

		/// @brief Cauchy random deviate generator
		///
		/// Generates Cauchy variates using the ratio of two uniform random numbers.
		/// The Cauchy distribution has no defined mean or variance.
		///
		/// Algorithm: If V1, V2 uniform on unit circle, then X = V1/V2 ~ Cauchy(0,1)
		class MML_STATISTICS_API CauchyDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _mu, _sigma;

		public:
			CauchyDeviate(Real location, Real scale, uint64_t seed = std::random_device{}());
			Real generate();
			Real operator()() { return generate(); }
		};

		/// @brief Normal random deviate generator using Leva's ratio-of-uniforms method
		///
		/// Fast and accurate method for generating normal deviates. More efficient
		/// than Box-Muller for single values (no caching needed).
		///
		/// Uses a ratio-of-uniforms method with quick acceptance tests.
		class MML_STATISTICS_API NormalDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _mu, _sigma;

		public:
			NormalDeviate(Real mean, Real stddev, uint64_t seed = std::random_device{}());
			Real generate();
			Real operator()() { return generate(); }
		};

		/// @brief Gamma random deviate generator
		///
		/// Generates gamma distributed random numbers using Marsaglia and Tsang's method.
		/// For shape parameter α < 1, uses a transformation.
		///
		/// The gamma distribution is used in Bayesian statistics, queuing theory,
		/// and reliability analysis.
		class MML_STATISTICS_API GammaDeviate
		{
		private:
			NormalDeviate _normalGen;
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _alpha;
			Real _originalAlpha;
			Real _beta;
			Real _a1, _a2;

		public:
			GammaDeviate(Real shape, Real scale, uint64_t seed = std::random_device{}());
			Real generate();
			Real operator()() { return generate(); }
		};

		/// @brief Poisson random deviate generator
		///
		/// Generates Poisson distributed random integers. Uses different algorithms
		/// depending on mean λ:
		///
		/// - λ < 5: Direct method (multiplicative algorithm)
		/// - λ ≥ 5: Ratio-of-uniforms method (faster for large λ)
		///
		/// The Poisson distribution models the number of events in a fixed interval.
		class MML_STATISTICS_API PoissonDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _lambda;
			Real _sqrtLambda, _logLambda, _lambdaExp;
			Real _previousLambda;
			std::vector<Real> _logFactorial;

		public:
			PoissonDeviate(Real mean, uint64_t seed = std::random_device{}());
			int generate();
			int generate(Real mean);
			int operator()() { return generate(); }
		};

		/// @brief Binomial random deviate generator
		///
		/// Generates binomial distributed random integers B(n,p).
		/// Uses different algorithms based on n and p:
		///
		/// - Small n: Bit comparison method
		/// - Medium n·p: Inverse CDF method
		///
		/// - Large n·p: Ratio-of-uniforms rejection method
		/// Models the number of successes in n independent Bernoulli trials.
		class MML_STATISTICS_API BinomialDeviate
		{
		private:
			std::mt19937_64 _gen;
			std::uniform_real_distribution<Real> _uniform;
			Real _pp, _p, _pb;
			Real _expNP, _np, _glnp, _plog, _pclog, _sq;
			int _n;
			int _method;
			
			uint64_t _uz, _uo, _unfin, _diff, _rltp;
			int _pbits[5];
			
			Real _cdf[64];
			
			Real _logFactorial[1024];

		public:
			BinomialDeviate(int trials, Real probability, uint64_t seed = std::random_device{}());
			int generate();
			int operator()() { return generate(); }
		};

	} // namespace Statistics
}  // namespace MML

#endif // MML_RANDOM_GENERATORS_H
