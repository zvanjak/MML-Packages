///////////////////////////////////////////////////////////////////////////////////////////
// RankCorrelation.h
// 
// Rank-based correlation methods for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// Contains: Spearman's rho, Kendall's tau, and helper functions
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_RANK_CORRELATION_H
#define MML_RANK_CORRELATION_H

#include <mml/algorithms/Statistics.h>
#include <mml/algorithms/Statistics/StatisticsBase.h>

#include "mml/mml_export.h"

namespace MML
{
	namespace Statistics
	{
		/// @brief Compute ranks for a dataset
		/// Converts raw data values to ranks (1-based). Tied values receive
		/// the average of the ranks they would have received.
		/// Example: [5, 2, 2, 8] → [3, 1.5, 1.5, 4]
		/// @param data Vector of values to rank
		/// @return Vector of ranks (1-based, with average ranks for ties)
		/// @note Uses average ranking for ties: if values at positions 2,3,4
		/// are tied, they all receive rank (2+3+4)/3 = 3
		/// Complexity: O(n log n) for sorting
		MML_STATISTICS_API Vector<Real> ComputeRanks(const Vector<Real>& data);

		/// @brief Compute Spearman's rank correlation coefficient
		///
		/// ρ (rho) = 1 - 6·Σd²ᵢ / (n(n²-1))  [no ties formula]
		/// For data with ties, computes Pearson correlation on ranks.
		///
		/// Spearman correlation measures monotonic relationship (not just linear).
		/// Values range from -1 (perfect negative monotonic) to +1 (perfect positive).
		///
		/// @param x First variable (vector of values)
		/// @param y Second variable (must have same size as x)
		///
		/// @return Spearman correlation coefficient ρ in [-1, 1]
		/// @throws StatisticsError if vectors are empty, have different sizes,
		///
		/// or have fewer than 2 elements
		/// Complexity: O(n log n) due to ranking
		MML_STATISTICS_API Real SpearmanCorrelation(const Vector<Real>& x, const Vector<Real>& y);

		/// @brief Result structure for rank correlation with significance test
		struct MML_STATISTICS_API RankCorrelationResult
		{
			Real rho;              ///< Correlation coefficient (Spearman or Kendall)
			Real zScore;           ///< z-score for large sample approximation
			int n;                 ///< Sample size
			///
			/// @brief Check if correlation is significant at given alpha level
			/// Uses normal approximation (valid for n > 10)
			///
			/// Uses Abramowitz & Stegun rational approximation for inverse normal CDF
			/// to compute critical z-value for any alpha in (0, 1).
			bool IsSignificant(Real alpha = 0.05) const;
		};

		/// @brief Compute Spearman correlation with significance test
		///
		/// Returns correlation coefficient with z-score for significance testing.
		/// Uses large-sample approximation: z = ρ·√(n-1)
		///
		/// @param x First variable
		/// @param y Second variable
		///
		/// @return RankCorrelationResult with rho, zScore, and n
		/// @throws StatisticsError if vectors are invalid
		MML_STATISTICS_API RankCorrelationResult SpearmanCorrelationWithTest(const Vector<Real>& x, const Vector<Real>& y);

		/// @brief Compute Kendall's tau-b rank correlation coefficient
		///
		/// τ_b = (C - D) / √((C + D + Tx)(C + D + Ty))
		/// where:
		///
		/// - C = number of concordant pairs
		/// - D = number of discordant pairs
		///
		/// - Tx = pairs tied only in x
		/// - Ty = pairs tied only in y
		///
		/// Kendall's tau is more robust than Spearman for small samples
		/// and has a more intuitive interpretation (probability difference).
		///
		/// @param x First variable (vector of values)
		/// @param y Second variable (must have same size as x)
		///
		/// @return Kendall tau-b coefficient in [-1, 1]
		/// @throws StatisticsError if vectors are empty, have different sizes,
		///
		/// or have fewer than 2 elements
		/// Complexity: O(n²) - naive implementation
		MML_STATISTICS_API Real KendallCorrelation(const Vector<Real>& x, const Vector<Real>& y);

		/// @brief Compute Kendall correlation with significance test
		///
		/// Returns correlation coefficient with z-score for significance testing.
		/// Uses asymptotic normal approximation:
		///
		/// z = τ / √(2(2n+5) / 9n(n-1))
		/// @param x First variable
		///
		/// @param y Second variable
		/// @return RankCorrelationResult with tau (as rho), zScore, and n
		///
		/// @throws StatisticsError if vectors are invalid
		MML_STATISTICS_API RankCorrelationResult KendallCorrelationWithTest(const Vector<Real>& x, const Vector<Real>& y);

		/******************************************************************************/
		/*****             Rank Correlation Detailed API                         *****/
		/******************************************************************************/

		/// Detailed variant of SpearmanCorrelation
		MML_STATISTICS_API RankCorrelationDetailedResult SpearmanCorrelationDetailed(
			const Vector<Real>& x,
			const Vector<Real>& y,
			const StatisticsConfig& config = {});

		/// Detailed variant of SpearmanCorrelationWithTest
		MML_STATISTICS_API RankCorrelationDetailedResult SpearmanCorrelationWithTestDetailed(
			const Vector<Real>& x,
			const Vector<Real>& y,
			const StatisticsConfig& config = {});

		/// Detailed variant of KendallCorrelation
		MML_STATISTICS_API RankCorrelationDetailedResult KendallCorrelationDetailed(
			const Vector<Real>& x,
			const Vector<Real>& y,
			const StatisticsConfig& config = {});

		/// Detailed variant of KendallCorrelationWithTest
		MML_STATISTICS_API RankCorrelationDetailedResult KendallCorrelationWithTestDetailed(
			const Vector<Real>& x,
			const Vector<Real>& y,
			const StatisticsConfig& config = {});
	}
}

#endif // MML_RANK_CORRELATION_H
