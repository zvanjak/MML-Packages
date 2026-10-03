///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        TimeSeries.h                                                        ///
///  Description: Time series analysis functions including moving averages,           ///
///               rolling statistics, autocorrelation, and differencing operations    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_TIME_SERIES_H
#define MML_TIME_SERIES_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "mml/mml_export.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

#include <deque>
#include <algorithm>
#include <cmath>

namespace MML {
	namespace Statistics {
		namespace TimeSeries {
			/////////////////////////////////////////////////////////////////////////////////////
			///                           MOVING AVERAGES                                     ///
			/////////////////////////////////////////////////////////////////////////////////////
			///
			/// @brief Simple Moving Average (SMA)
			/// @details Computes the unweighted mean of the previous windowSize data points.
			///          Formula: SMA_t = (1/windowSize) * sum(x_{t-windowSize+1} to x_t)
			/// @param data Input time series data
			/// @param windowSize Number of points in the moving window
			/// @return Vector of size (n - windowSize + 1) containing smoothed values
			/// @throws StatisticsError if windowSize > data.size() or windowSize < 1 or data is empty
			MML_STATISTICS_API Vector<Real> SimpleMovingAverage(const Vector<Real>& data, int windowSize);

			/// @brief Exponential Moving Average (EMA) with smoothing factor alpha
			/// @details Applies exponential weighting where recent values have more influence.
			///          Formula: EMA_t = alpha * x_t + (1 - alpha) * EMA_{t-1}
			///          Initial EMA_0 = x_0
			/// @param data Input time series data
			/// @param alpha Smoothing factor in range (0, 1). Higher alpha = more weight on recent values.
			/// @return Vector of same size as input containing EMA values
			/// @throws StatisticsError if alpha not in (0,1) or data is empty
			MML_STATISTICS_API Vector<Real> ExponentialMovingAverage(const Vector<Real>& data, Real alpha);

			/// @brief Exponential Moving Average (EMA) with span parameter
			/// @details Computes alpha from span: alpha = 2 / (span + 1)
			///          This is the common convention used in pandas.
			/// @param data Input time series data
			/// @param span The span parameter (must be >= 1)
			/// @return Vector of same size as input containing EMA values
			/// @throws StatisticsError if span < 1 or data is empty
			MML_STATISTICS_API Vector<Real> ExponentialMovingAverageSpan(const Vector<Real>& data, int span);

			/// @brief Weighted Moving Average (WMA) with custom weights
			/// @details Computes weighted average over sliding window using provided weights.
			///          Formula: WMA_t = sum(w_i * x_{t-n+1+i}) / sum(w_i)
			/// @param data Input time series data
			/// @param weights Weight vector (window size = weights.size())
			/// @return Vector of size (n - windowSize + 1) containing weighted averages
			/// @throws StatisticsError if weights empty, window > data length, or all weights zero
			MML_STATISTICS_API Vector<Real> WeightedMovingAverage(const Vector<Real>& data, const Vector<Real>& weights);

			/// @brief Weighted Moving Average (WMA) with linear weights
			/// @details Uses linearly increasing weights: [1, 2, 3, ..., windowSize]
			///          Most recent value gets highest weight.
			/// @param data Input time series data
			/// @param windowSize Number of points in the moving window
			/// @return Vector of size (n - windowSize + 1) containing weighted averages
			MML_STATISTICS_API Vector<Real> WeightedMovingAverageLinear(const Vector<Real>& data, int windowSize);


			/////////////////////////////////////////////////////////////////////////////////////
			///                          ROLLING STATISTICS                                   ///
			/////////////////////////////////////////////////////////////////////////////////////
			///
			/// @brief Rolling Mean (equivalent to SMA)
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window
			/// @return Vector of size (n - windowSize + 1) containing rolling means
			inline Vector<Real> RollingMean(const Vector<Real>& data, int windowSize) { return SimpleMovingAverage(data, windowSize); }

			/// @brief Rolling Sum over a sliding window
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window
			/// @return Vector of size (n - windowSize + 1) containing rolling sums
			MML_STATISTICS_API Vector<Real> RollingSum(const Vector<Real>& data, int windowSize);

			/// @brief Rolling Variance using Welford's online algorithm for numerical stability
			/// @details Computes sample variance (n-1 denominator) over sliding window.
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window (must be >= 2)
			/// @return Vector of size (n - windowSize + 1) containing rolling variances
			MML_STATISTICS_API Vector<Real> RollingVariance(const Vector<Real>& data, int windowSize);

			/// @brief Rolling Standard Deviation
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window (must be >= 2)
			/// @return Vector of size (n - windowSize + 1) containing rolling standard deviations
			MML_STATISTICS_API Vector<Real> RollingStdDev(const Vector<Real>& data, int windowSize);

			/// @brief Rolling Minimum using deque-based O(n) algorithm
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window
			/// @return Vector of size (n - windowSize + 1) containing rolling minimums
			MML_STATISTICS_API Vector<Real> RollingMin(const Vector<Real>& data, int windowSize);

			/// @brief Rolling Maximum using deque-based O(n) algorithm
			/// @param data Input time series data
			/// @param windowSize Number of points in the rolling window
			/// @return Vector of size (n - windowSize + 1) containing rolling maximums
			MML_STATISTICS_API Vector<Real> RollingMax(const Vector<Real>& data, int windowSize);


			/////////////////////////////////////////////////////////////////////////////////////
			///                         AUTOCORRELATION                                       ///
			/////////////////////////////////////////////////////////////////////////////////////
			///
			/// @brief Single lag autocorrelation coefficient
			/// @details Measures correlation between x_t and x_{t-lag}
			///          Formula: ACF(k) = Cov(x_t, x_{t-k}) / Var(x_t)
			/// @param data Input time series data
			/// @param lag The lag value (must be >= 0 and < data.size())
			/// @return Autocorrelation coefficient in [-1, 1], ACF(0) = 1.0
			MML_STATISTICS_API Real Autocorrelation(const Vector<Real>& data, int lag);

			/// @brief Autocorrelation Function (ACF) for multiple lags
			/// @details Computes ACF for lags 0, 1, 2, ..., maxLag
			/// @param data Input time series data
			/// @param maxLag Maximum lag to compute (inclusive)
			/// @return Vector of size (maxLag + 1) with ACF values
			MML_STATISTICS_API Vector<Real> AutocorrelationFunction(const Vector<Real>& data, int maxLag);

			/// @brief Partial Autocorrelation using Durbin-Levinson recursion
			/// @details Measures direct correlation between x_t and x_{t-k}, removing
			///          indirect effects through intermediate lags.
			/// @param data Input time series data
			/// @param lag The lag value (must be >= 0 and < data.size())
			/// @return Partial autocorrelation coefficient at given lag
			MML_STATISTICS_API Real PartialAutocorrelation(const Vector<Real>& data, int lag);

			/// @brief Partial Autocorrelation Function (PACF) for multiple lags
			/// @details Computes PACF for lags 0, 1, 2, ..., maxLag using Durbin-Levinson
			/// @param data Input time series data
			/// @param maxLag Maximum lag to compute (inclusive)
			/// @return Vector of size (maxLag + 1) with PACF values
			MML_STATISTICS_API Vector<Real> PartialAutocorrelationFunction(const Vector<Real>& data, int maxLag);


			/////////////////////////////////////////////////////////////////////////////////////
			///                      DIFFERENCING AND LAG OPERATIONS                          ///
			/////////////////////////////////////////////////////////////////////////////////////
			///
			/// @brief Create lagged version of time series
			/// @details Returns data shifted by k positions: {x_k, x_{k+1}, ..., x_{n-1}}
			/// @param data Input time series data
			/// @param k Lag amount (number of positions to shift)
			/// @return Vector of size (n - k) with lagged values
			MML_STATISTICS_API Vector<Real> Lag(const Vector<Real>& data, int k);

			/// @brief Create lag matrix with columns for lags 0, 1, ..., maxLag
			/// @details Each column i contains the series lagged by i positions.
			///          Row j, Column i = data[j + i] for valid indices.
			/// @param data Input time series data
			/// @param maxLag Maximum lag (number of columns - 1)
			/// @return Matrix of size (n - maxLag) x (maxLag + 1)
			MML_STATISTICS_API Matrix<Real> LagMatrix(const Vector<Real>& data, int maxLag);

			/// @brief First difference of time series
			/// @details Computes d_i = x_{i+1} - x_i
			/// @param data Input time series data
			/// @return Vector of size (n - 1) containing first differences
			MML_STATISTICS_API Vector<Real> FirstDifference(const Vector<Real>& data);

			/// @brief Apply differencing multiple times
			/// @details order=1 is first difference, order=2 is difference of differences, etc.
			/// @param data Input time series data
			/// @param order Number of times to difference (must be >= 1)
			/// @return Vector of size (n - order) containing differenced values
			MML_STATISTICS_API Vector<Real> Difference(const Vector<Real>& data, int order);

			/// @brief Seasonal difference
			/// @details Computes d_i = x_i - x_{i-period}
			/// @param data Input time series data
			/// @param period Seasonal period (e.g., 12 for monthly, 4 for quarterly)
			/// @return Vector of size (n - period) containing seasonal differences
			MML_STATISTICS_API Vector<Real> SeasonalDifference(const Vector<Real>& data, int period);

			/// @brief Cumulative sum (undoes differencing)
			/// @details cumsum[0] = initial, cumsum[i] = cumsum[i-1] + differences[i-1]
			/// @param differences The differenced series
			/// @param initial Starting value (typically the first value of original series)
			/// @return Vector of size (differences.size() + 1) with reconstructed values
			MML_STATISTICS_API Vector<Real> CumulativeSum(const Vector<Real>& differences, Real initial);

			/// @brief Cumulative sum without initial value (starts at 0)
			/// @param differences The differenced series
			/// @return Vector of size differences.size() with cumulative sums
			MML_STATISTICS_API Vector<Real> CumulativeSum(const Vector<Real>& differences);

		} // namespace TimeSeries
	} // namespace Statistics
} // namespace MML

#endif // MML_TIME_SERIES_H
