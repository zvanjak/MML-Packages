///////////////////////////////////////////////////////////////////////////////////////////
// ConfidenceIntervals.h
//
// Confidence interval functions for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// Requires: Statistics.h (for basic stats functions and distributions)
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_CONFIDENCE_INTERVALS_H
#define MML_CONFIDENCE_INTERVALS_H

#include <mml/algorithms/Statistics.h>
#include <mml/algorithms/Statistics/StatisticsBase.h>
#include <mml/algorithms/Statistics/Distributions.h>

#include "mml/mml_export.h"

namespace MML {
	namespace Statistics {
		/// @brief Structure to hold confidence interval results
		struct MML_STATISTICS_API ConfidenceInterval {
			Real estimate;		   // Point estimate
			Real lowerBound;	   // Lower confidence limit
			Real upperBound;	   // Upper confidence limit
			Real marginOfError;	   // Half-width of interval
			Real confidenceLevel;  // Confidence level (e.g., 0.95 for 95%)
			std::string parameter; // What we're estimating (e.g., "Mean", "Proportion")

			ConfidenceInterval(Real est, Real lower, Real upper, Real margin, Real conf, const std::string& param)
				: estimate(est)
				, lowerBound(lower)
				, upperBound(upper)
				, marginOfError(margin)
				, confidenceLevel(conf)
				, parameter(param) {}
		};

		/// @brief Confidence Interval for Population Mean (single sample)
		///
		/// Uses t-distribution for unknown population variance.
		/// CI = x̄ ± t(α/2, df) × (s/√n)
		///
		/// @param sample Sample data
		/// @param confidenceLevel Confidence level (default: 0.95 for 95% CI)
		///
		/// @return ConfidenceInterval with bounds
		MML_STATISTICS_API ConfidenceInterval ConfidenceIntervalMean(const Vector<Real>& sample, Real confidenceLevel = 0.95);

		/// @brief Confidence Interval for Difference of Two Means (independent samples)
		///
		/// Uses pooled variance assuming equal population variances.
		/// CI = (x̄₁ - x̄₂) ± t(α/2, df) × SE
		///
		/// @param sample1 First sample
		/// @param sample2 Second sample
		///
		/// @param confidenceLevel Confidence level (default: 0.95)
		/// @return ConfidenceInterval for μ₁ - μ₂
		MML_STATISTICS_API ConfidenceInterval ConfidenceIntervalMeanDifference(const Vector<Real>& sample1, const Vector<Real>& sample2,
																			   Real confidenceLevel = 0.95);

		/// @brief Confidence Interval for Population Proportion
		///
		/// Uses normal approximation (valid for large samples).
		/// CI = p̂ ± z(α/2) × √(p̂(1-p̂)/n)
		///
		/// @param successes Number of successes
		/// @param trials Total number of trials
		///
		/// @param confidenceLevel Confidence level (default: 0.95)
		/// @return ConfidenceInterval for proportion
		MML_STATISTICS_API ConfidenceInterval ConfidenceIntervalProportion(int successes, int trials, Real confidenceLevel = 0.95);

		/// @brief Confidence Interval for Difference of Two Proportions
		///
		/// Uses normal approximation.
		/// CI = (p̂₁ - p̂₂) ± z(α/2) × √(p̂₁(1-p̂₁)/n₁ + p̂₂(1-p̂₂)/n₂)
		///
		/// @param successes1 Successes in first sample
		/// @param trials1 Trials in first sample
		///
		/// @param successes2 Successes in second sample
		/// @param trials2 Trials in second sample
		///
		/// @param confidenceLevel Confidence level (default: 0.95)
		/// @return ConfidenceInterval for p₁ - p₂
		MML_STATISTICS_API ConfidenceInterval ConfidenceIntervalProportionDifference(int successes1, int trials1, int successes2, int trials2,
																					 Real confidenceLevel = 0.95);

		/// @brief Confidence Interval for Mean of Paired Differences
		///
		/// For paired data (before-after, matched pairs).
		/// CI = d̄ ± t(α/2, df) × (s_d/√n)
		///
		/// @param before First measurements
		/// @param after Second measurements
		///
		/// @param confidenceLevel Confidence level (default: 0.95)
		/// @return ConfidenceInterval for mean difference
		MML_STATISTICS_API ConfidenceInterval ConfidenceIntervalPairedDifference(const Vector<Real>& before, const Vector<Real>& after,
																				 Real confidenceLevel = 0.95);

		/******************************************************************************/
		/*****             Confidence Interval Detailed API                      *****/
		/******************************************************************************/

		/// Detailed variant of ConfidenceIntervalMean
		MML_STATISTICS_API ConfidenceIntervalDetailedResult ConfidenceIntervalMeanDetailed(
			const Vector<Real>& sample,
			Real confidenceLevel = 0.95,
			const StatisticsConfig& config = {});

		/// Detailed variant of ConfidenceIntervalMeanDifference
		MML_STATISTICS_API ConfidenceIntervalDetailedResult ConfidenceIntervalMeanDifferenceDetailed(
			const Vector<Real>& sample1,
			const Vector<Real>& sample2,
			Real confidenceLevel = 0.95,
			const StatisticsConfig& config = {});

		/// Detailed variant of ConfidenceIntervalProportion
		MML_STATISTICS_API ConfidenceIntervalDetailedResult ConfidenceIntervalProportionDetailed(
			int successes,
			int trials,
			Real confidenceLevel = 0.95,
			const StatisticsConfig& config = {});

		/// Detailed variant of ConfidenceIntervalProportionDifference
		MML_STATISTICS_API ConfidenceIntervalDetailedResult ConfidenceIntervalProportionDifferenceDetailed(
			int successes1, int trials1,
			int successes2, int trials2,
			Real confidenceLevel = 0.95,
			const StatisticsConfig& config = {});

		/// Detailed variant of ConfidenceIntervalPairedDifference
		MML_STATISTICS_API ConfidenceIntervalDetailedResult ConfidenceIntervalPairedDifferenceDetailed(
			const Vector<Real>& before,
			const Vector<Real>& after,
			Real confidenceLevel = 0.95,
			const StatisticsConfig& config = {});

	} // namespace Statistics
} // namespace MML

#endif // MML_CONFIDENCE_INTERVALS_H
