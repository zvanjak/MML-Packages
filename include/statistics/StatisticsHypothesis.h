///////////////////////////////////////////////////////////////////////////////////////////
// StatisticsHypothesis.h  (MML-Packages mml_ext - moved from MML core for Release 2.0)
//
// Statistical hypothesis testing functions for MML (MinimalMathLibrary)
//
// Requires: MML core (Statistics.h, StatisticsBase.h, Distributions.h)
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_HYPOTHESIS_TESTING_H
#define MML_HYPOTHESIS_TESTING_H

#include <mml/algorithms/Statistics.h>
#include <mml/algorithms/Statistics/StatisticsBase.h>
#include <mml/algorithms/Statistics/Distributions.h>

#include "mml/mml_export.h"

namespace MML
{
	namespace Statistics
	{
		/// @brief Result structure for hypothesis tests
		/// Contains all relevant information from a hypothesis test:
		/// - Test statistic and its distribution
		/// - P-value and decision at given confidence level
		/// - Degrees of freedom if applicable
		struct MML_STATISTICS_API HypothesisTestResult
		{
			Real testStatistic;      ///< Computed test statistic (t, z, chi-square, F, etc.)
			Real pValue;             ///< P-value for the test
			Real criticalValue;      ///< Critical value at the given significance level
			bool rejectNull;         ///< Whether to reject the null hypothesis
			Real confidenceLevel;    ///< Confidence level used (e.g., 0.95 for 95%)
			int degreesOfFreedom;    ///< Degrees of freedom (if applicable, -1 otherwise)
			std::string testName;    ///< Descriptive name of the test
			///
			/// @brief Construct a hypothesis test result
			HypothesisTestResult(
				Real statistic = 0.0,
				Real p = 0.0,
				Real critical = 0.0,
				bool reject = false,
				Real confidence = 0.95,
				int df = -1,
				const std::string& name = "Hypothesis Test"
			) : testStatistic(statistic),
			    pValue(p),
			    criticalValue(critical),
			    rejectNull(reject),
			    confidenceLevel(confidence),
			    degreesOfFreedom(df),
			    testName(name)
			{}
		};

		/// @brief One-sample t-test
		///
		/// Tests the null hypothesis H0: μ = μ0 against H1: μ ≠ μ0
		/// where μ is the population mean.
		///
		/// Test statistic: t = (x̄ - μ0) / (s / √n)
		/// where x̄ is sample mean, s is sample standard deviation, n is sample size
		///
		/// @param sample Vector of sample values
		/// @param mu0 Hypothesized population mean
		///
		/// @param alpha Significance level (default: 0.05 for 95% confidence)
		/// @return HypothesisTestResult with test details
		MML_STATISTICS_API HypothesisTestResult OneSampleTTest(
			const Vector<Real>& sample,
			Real mu0,
			Real alpha = 0.05
		);

		/// @brief Two-sample t-test with equal variances (pooled)
		///
		/// Tests H0: μ1 = μ2 against H1: μ1 ≠ μ2
		/// Assumes equal population variances (use Welch's t-test if unequal).
		///
		/// Test statistic: t = (x̄1 - x̄2) / (sp · √(1/n1 + 1/n2))
		/// where sp² = ((n1-1)s1² + (n2-1)s2²) / (n1 + n2 - 2) is pooled variance
		///
		/// @param sample1 First sample
		/// @param sample2 Second sample
		///
		/// @param alpha Significance level (default: 0.05)
		/// @return HypothesisTestResult with test details
		MML_STATISTICS_API HypothesisTestResult TwoSampleTTest(
			const Vector<Real>& sample1,
			const Vector<Real>& sample2,
			Real alpha = 0.05
		);

		/// @brief Welch's t-test (unequal variances)
		///
		/// Tests H0: μ1 = μ2 against H1: μ1 ≠ μ2
		/// Does NOT assume equal variances (more robust than pooled t-test).
		///
		/// Test statistic: t = (x̄1 - x̄2) / √(s1²/n1 + s2²/n2)
		/// Degrees of freedom use Welch-Satterthwaite equation
		///
		/// @param sample1 First sample
		/// @param sample2 Second sample
		///
		/// @param alpha Significance level (default: 0.05)
		/// @return HypothesisTestResult with test details
		MML_STATISTICS_API HypothesisTestResult WelchTTest(
			const Vector<Real>& sample1,
			const Vector<Real>& sample2,
			Real alpha = 0.05
		);

		/// @brief Paired t-test
		///
		/// Tests H0: μd = 0 against H1: μd ≠ 0
		/// where μd is the mean difference between paired observations.
		///
		/// Equivalent to one-sample t-test on differences.
		/// @param before First measurement (e.g., before treatment)
		///
		/// @param after Second measurement (e.g., after treatment)
		/// @param alpha Significance level (default: 0.05)
		///
		/// @return HypothesisTestResult with test details
		MML_STATISTICS_API HypothesisTestResult PairedTTest(
			const Vector<Real>& before,
			const Vector<Real>& after,
			Real alpha = 0.05
		);

		/// @brief Chi-Square Goodness-of-Fit Test
		///
		/// Tests whether observed frequencies match expected frequencies.
		/// H0: Data follows the expected distribution
		///
		/// H1: Data does not follow the expected distribution
		/// Test statistic: χ² = Σ((O_i - E_i)² / E_i)
		///
		/// @param observed Vector of observed frequencies
		/// @param expected Vector of expected frequencies (must sum to same as observed)
		///
		/// @param alpha Significance level (default: 0.05)
		/// @return HypothesisTestResult with test details
		///
		/// @throws StatisticsError if sizes don't match, frequencies invalid, or expected has zeros
		MML_STATISTICS_API HypothesisTestResult ChiSquareGoodnessOfFit(
			const Vector<Real>& observed,
			const Vector<Real>& expected,
			Real alpha = 0.05
		);

		/// @brief Chi-Square Test of Independence
		///
		/// Tests whether two categorical variables are independent.
		/// H0: Variables are independent
		///
		/// H1: Variables are dependent (associated)
		/// @param contingencyTable Matrix of observed frequencies (rows × columns)
		///
		/// @param alpha Significance level (default: 0.05)
		/// @return HypothesisTestResult with test details
		///
		/// @throws StatisticsError if table too small or has invalid frequencies
		MML_STATISTICS_API HypothesisTestResult ChiSquareIndependence(
			const Matrix<Real>& contingencyTable,
			Real alpha = 0.05
		);

		/// @brief One-Way ANOVA (Analysis of Variance)
		///
		/// Tests whether the means of three or more groups are equal.
		/// H0: μ₁ = μ₂ = μ₃ = ... = μₖ (all group means equal)
		///
		/// H1: At least one mean differs
		/// F-statistic = (Between-group variance) / (Within-group variance)
		///
		/// = (SSB/df_between) / (SSW/df_within)
		/// @param groups Vector of groups, where each group is a vector of observations
		///
		/// @param alpha Significance level (default: 0.05)
		/// @return HypothesisTestResult with F-statistic and p-value
		///
		/// @throws StatisticsError if fewer than 2 groups, any group too small, or empty groups
		MML_STATISTICS_API HypothesisTestResult OneWayANOVA(
			const std::vector<Vector<Real>>& groups,
			Real alpha = 0.05
		);

		/******************************************************************************/
		/*****             Hypothesis Testing Detailed API                       *****/
		/******************************************************************************/

		/// Detailed variant of OneSampleTTest
		MML_STATISTICS_API HypothesisTestDetailedResult OneSampleTTestDetailed(
			const Vector<Real>& sample,
			Real mu0,
			Real alpha = 0.05,
			const StatisticsConfig& config = {});

		/// Detailed variant of TwoSampleTTest
		MML_STATISTICS_API HypothesisTestDetailedResult TwoSampleTTestDetailed(
			const Vector<Real>& sample1,
			const Vector<Real>& sample2,
			Real alpha = 0.05,
			const StatisticsConfig& config = {});

		/// Detailed variant of WelchTTest
		MML_STATISTICS_API HypothesisTestDetailedResult WelchTTestDetailed(
			const Vector<Real>& sample1,
			const Vector<Real>& sample2,
			Real alpha = 0.05,
			const StatisticsConfig& config = {});

		/// Detailed variant of PairedTTest
		MML_STATISTICS_API HypothesisTestDetailedResult PairedTTestDetailed(
			const Vector<Real>& before,
			const Vector<Real>& after,
			Real alpha = 0.05,
			const StatisticsConfig& config = {});

		/// Detailed variant of ChiSquareGoodnessOfFit
		MML_STATISTICS_API HypothesisTestDetailedResult ChiSquareGoodnessOfFitDetailed(
			const Vector<Real>& observed,
			const Vector<Real>& expected,
			Real alpha = 0.05,
			const StatisticsConfig& config = {});

		/// Detailed variant of ChiSquareIndependence
		MML_STATISTICS_API HypothesisTestDetailedResult ChiSquareIndependenceDetailed(
			const Matrix<Real>& contingencyTable,
			Real alpha = 0.05,
			const StatisticsConfig& config = {});

		/// Detailed variant of OneWayANOVA
		MML_STATISTICS_API HypothesisTestDetailedResult OneWayANOVADetailed(
			const std::vector<Vector<Real>>& groups,
			Real alpha = 0.05,
			const StatisticsConfig& config = {});

	} // namespace Statistics
}  // namespace MML

#endif // MML_HYPOTHESIS_TESTING_H
