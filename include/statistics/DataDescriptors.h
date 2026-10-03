///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataDescriptors.h                                                   ///
///  Description: Typed data analysis classes beyond Real numbers                     ///
///               BoolDataStats, TimeDataStats, DateDataStats, GeoDataStats,         ///
///               CategoryStats                                                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_DESCRIPTORS_H
#define MML_DATA_DESCRIPTORS_H

#include "MMLBase.h"
#include "MMLExceptions.h"
#include "base/Vector/Vector.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <sstream>
#include <iomanip>

#include "mml/mml_export.h"

namespace MML
{
namespace Statistics
{

/////////////////////////////////////////////////////////////////////////////////////
///                          BOOLEAN DATA STATISTICS                               ///
/////////////////////////////////////////////////////////////////////////////////////
///
/// @brief Statistics for boolean data
/// Computes proportion-based statistics including confidence intervals
///
/// and chi-square test for expected proportion.
struct MML_STATISTICS_API BoolDataStats
{
	size_t count;           ///< Total count
	size_t countTrue;       ///< Count of true values
	size_t countFalse;      ///< Count of false values
	Real proportionTrue;    ///< Proportion of true values (p = countTrue/count)
	Real standardError;     ///< Standard error of proportion: sqrt(p(1-p)/n)
	Real ciLower95;         ///< 95% confidence interval lower bound
	Real ciUpper95;         ///< 95% confidence interval upper bound

	static BoolDataStats Compute(const std::vector<bool>& data);
	Real ChiSquare(Real expectedProportion) const;
	void PrintSummary(std::ostream& os = std::cout) const;
};

/////////////////////////////////////////////////////////////////////////////////////
///                           TIME DATA STATISTICS                                 ///
/////////////////////////////////////////////////////////////////////////////////////
///
/// @brief Statistics for time/duration data (seconds since epoch or duration in seconds)
/// Supports both timestamps and durations. For timestamps, provides
///
/// distributions by hour of day and day of week.
struct MML_STATISTICS_API TimeDataStats
{
	size_t count;           ///< Number of time values
	Real minSeconds;        ///< Minimum time in seconds
	Real maxSeconds;        ///< Maximum time in seconds
	Real avgSeconds;        ///< Average time in seconds
	Real medianSeconds;     ///< Median time in seconds
	Real stdDevSeconds;     ///< Standard deviation in seconds
	Real rangeSeconds;      ///< Range (max - min) in seconds

	std::vector<size_t> hourDistribution;    ///< Count by hour (0-23)
	std::vector<size_t> dayOfWeekDistribution; ///< Count by day (0=Sun, 6=Sat)

	static TimeDataStats Compute(const Vector<Real>& data, bool isTimestamp = true);
	static std::string SecondsToHMS(Real seconds);
	static Real SecondsToDays(Real seconds) { return seconds / 86400.0; }
	void PrintSummary(std::ostream& os = std::cout) const;
};

/////////////////////////////////////////////////////////////////////////////////////
///                           DATE DATA STATISTICS                                 ///
/////////////////////////////////////////////////////////////////////////////////////
///
/// @brief Statistics for date data (days since epoch or structured dates)
/// Provides year/month/day-of-week frequency distributions and
///
/// weekday vs weekend analysis.
struct MML_STATISTICS_API DateDataStats
{
	size_t count;           ///< Number of date values
	int minYear;            ///< Minimum year
	int maxYear;            ///< Maximum year
	int spanDays;           ///< Span in days (max - min)
	
	std::map<int, size_t> yearDistribution;      ///< Count by year
	std::vector<size_t> monthDistribution;       ///< Count by month (0=Jan, 11=Dec)
	std::vector<size_t> dayOfWeekDistribution;   ///< Count by day (0=Sun, 6=Sat)
	size_t weekdayCount;    ///< Mon-Fri count
	size_t weekendCount;    ///< Sat-Sun count

	static DateDataStats Compute(const Vector<Real>& epochDays, bool isSeconds = false);

	int GetModeMonth() const
	{
		return static_cast<int>(std::max_element(monthDistribution.begin(), 
		                                          monthDistribution.end()) - monthDistribution.begin());
	}

	Real GetWeekdayProportion() const
	{
		return static_cast<Real>(weekdayCount) / count;
	}

	void PrintSummary(std::ostream& os = std::cout) const;
};

/////////////////////////////////////////////////////////////////////////////////////
///                      GEOGRAPHIC DATA STATISTICS                                ///
/////////////////////////////////////////////////////////////////////////////////////
///
/// @brief Statistics for geographic (lat/lon) data
/// Computes centroid, bounding box, and average distance from centroid.
///
/// Uses WGS84 coordinates (standard GPS).
struct MML_STATISTICS_API GeoDataStats
{
	size_t count;           ///< Number of points
	Real centroidLat;       ///< Centroid latitude
	Real centroidLon;       ///< Centroid longitude
	Real minLat, maxLat;    ///< Latitude bounding box
	Real minLon, maxLon;    ///< Longitude bounding box
	Real avgDistanceFromCentroid;  ///< Average distance from centroid (km)
	Real maxDistanceFromCentroid;  ///< Maximum distance from centroid (km)

	static constexpr Real EARTH_RADIUS_KM = 6371.0;

	static Real HaversineDistance(Real lat1, Real lon1, Real lat2, Real lon2);
	static GeoDataStats Compute(const Vector<Real>& latitudes, const Vector<Real>& longitudes);

	Real GetBoundingBoxDiagonal() const
	{
		return HaversineDistance(minLat, minLon, maxLat, maxLon);
	}

	void PrintSummary(std::ostream& os = std::cout) const;
};

/////////////////////////////////////////////////////////////////////////////////////
///                        CATEGORICAL DATA STATISTICS                             ///
/////////////////////////////////////////////////////////////////////////////////////
///
/// @brief Statistics for categorical (string) data
/// Computes frequency distribution, unique count, mode, and entropy.
struct MML_STATISTICS_API CategoryStats
{
	size_t count;           ///< Total number of values
	size_t uniqueCount;     ///< Number of unique categories
	std::string mode;       ///< Most frequent category
	size_t modeCount;       ///< Count of most frequent category
	Real entropy;           ///< Shannon entropy (bits)
	
	std::map<std::string, size_t> frequencies;  ///< Category -> count

	static CategoryStats Compute(const std::vector<std::string>& data);
	std::vector<std::pair<std::string, size_t>> GetTopN(size_t n) const;

	Real GetModeProportion() const
	{
		return static_cast<Real>(modeCount) / count;
	}

	void PrintSummary(std::ostream& os = std::cout) const;
};

/////////////////////////////////////////////////////////////////////////////////////
///                         COMBINED DATA DESCRIPTOR                               ///
/////////////////////////////////////////////////////////////////////////////////////
///
/// @brief Unified interface for computing type-appropriate statistics
/// Provides static methods to compute stats based on detected or specified data type.
struct MML_STATISTICS_API DataDescriptor
{
	enum class Type { REAL, BOOL, TIME, DATE, GEO, CATEGORY };
	static Type DetectType(const std::vector<std::string>& samples);
};

} // namespace Statistics
} // namespace MML

#endif // MML_DATA_DESCRIPTORS_H
