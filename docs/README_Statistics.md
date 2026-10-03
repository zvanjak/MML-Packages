# Statistics Package

**Statistical Computing and Analysis**

The Statistics package provides comprehensive tools for statistical analysis, including probability distributions, hypothesis testing, descriptive statistics, and time series analysis.

## Features

### Probability Distributions

**Continuous Distributions (13 total):**
- **Normal (Gaussian)** - Central limit theorem, general-purpose
- **Student's t** - Small sample inference, heavy tails
- **Chi-Square (χ²)** - Variance testing, goodness of fit
- **F-Distribution** - ANOVA, variance ratio tests
- **Cauchy** - Heavy-tailed, undefined mean/variance
- **Exponential** - Memoryless property, survival analysis
- **Logistic** - Logistic regression, neural networks
- **Gamma** - Generalizes exponential and chi-square
- **Beta** - Models proportions and probabilities on [0,1]
- **Weibull** - Reliability analysis, failure modeling
- **Pareto** - Power-law distribution (wealth, file sizes)
- **Log-Normal** - Multiplicative processes
- **Uniform** - Constant probability over interval

**Discrete Distributions (6 total):**
- **Bernoulli** - Single trial success/failure
- **Binomial** - Number of successes in n trials
- **Poisson** - Rare events, count data
- **Geometric** - Trials until first success
- **Negative Binomial** - Trials until r-th success
- **Hypergeometric** - Sampling without replacement

### Descriptive Statistics
- **Central Tendency** - Mean, median, mode, trimmed mean
- **Dispersion** - Explicit sample/population and frequency-weighted variance and standard deviation, range, IQR
- **Shape** - Skewness, kurtosis, moments
- **Percentiles** - Quartiles, deciles, arbitrary quantiles
- **Histograms** - Frequency distributions, binning

### Hypothesis Testing
- **t-Tests** - One-sample, pooled two-sample, Welch, paired
- **Chi-Square Tests** - Independence, goodness of fit
- **ANOVA** - One-way comparison of group means

### Correlation Analysis
- **Pearson** - Linear correlation coefficient
- **Weighted Pearson** - Frequency-weighted covariance and correlation
- **Matrix Helpers** - Dense sample covariance and correlation matrices
- **Spearman** - Rank correlation
- **Kendall's Tau** - Concordance measure

### Confidence Intervals
- **Means** - Single-sample, independent difference, paired difference
- **Proportions** - Normal-approximation single and difference intervals

### Specialized Data Types
- **Boolean Data** - Proportion statistics, chi-square tests
- **Time Data** - Duration analysis
- **Categorical Data** - Frequency tables, mode

### Time Series (Basic)
- **Autocorrelation** - ACF computation
- **Moving Averages** - Simple and exponential
- **Trend Detection** - Linear trend fitting

## Quick Start

### Descriptive Statistics

```cpp
#include "DataDescriptors.h"
#include "algorithms/Statistics.h"

using namespace MML;
using namespace MML::Statistics;

// Sample data
Vector<Real> data({2.3, 4.5, 3.1, 5.2, 4.8, 3.9, 4.2});

// Basic statistics
Real mean = Mean(data);
Real sampleVariance = SampleVariance(data);         // denominator n-1
Real populationVariance = PopulationVariance(data); // denominator n
Real sampleStdDev = SampleStdDev(data);
Real populationStdDev = PopulationStdDev(data);
Real median = Median(data);

// Higher moments
Real averageDeviation, skewness, kurtosis;
Moments(data, mean, averageDeviation, sampleStdDev, sampleVariance, skewness, kurtosis);

std::cout << "Mean: " << mean << "\n";
std::cout << "Sample Std Dev: " << sampleStdDev << "\n";
std::cout << "Skewness: " << skewness << "\n";
```

### Probability Distributions

```cpp
#include "algorithms/Statistics/Distributions.h"
#include "DiscreteDistributions.h"

using namespace MML;
using namespace MML::Statistics;

// Normal distribution
NormalDistribution normal(0.0, 1.0);  // μ=0, σ=1
double pdf = normal.pdf(1.5);         // P(X=1.5)
double cdf = normal.cdf(1.96);        // P(X ≤ 1.96) ≈ 0.975
double quantile = normal.inverseCdf(0.95);  // 95th percentile

// Student's t distribution
TDistribution t(10);  // 10 degrees of freedom
double tCritical = t.inverseCdf(0.975);  // Two-tailed 95% critical value

// Exponential distribution (λ = 0.5)
ExponentialDistribution exp(0.5);
double expMean = exp.mean();        // = 1/λ = 2.0
double expVar = exp.variance();     // = 1/λ² = 4.0

// Chi-square distribution
ChiSquareDistribution chi2(5);  // 5 degrees of freedom
double chiCritical = chi2.criticalValue(0.05); // 95th percentile
```

### Hypothesis Testing

```cpp
#include "algorithms/Statistics/StatisticsHypothesis.h"

using namespace MML;
using namespace MML::Statistics;

// Sample data
Vector<Real> sample1 = {23.1, 25.4, 22.8, 24.5, 26.1, 23.9};
Vector<Real> sample2 = {21.2, 22.5, 20.8, 21.9, 22.3, 21.7};

// One-sample t-test: H0: μ = 24.0
auto result1 = OneSampleTTest(sample1, 24.0);
std::cout << "t-statistic: " << result1.testStatistic << "\n";
std::cout << "p-value: " << result1.pValue << "\n";
std::cout << "Reject null: " << (result1.rejectNull ? "Yes" : "No") << "\n";

// Two-sample t-test: H0: μ1 = μ2
auto result2 = TwoSampleTTest(sample1, sample2);

// Paired t-test (for matched samples)
auto result3 = PairedTTest(before, after);

// Chi-square test for independence
auto chiResult = ChiSquareIndependence(contingencyTable);
```

### Confidence Intervals

```cpp
#include "algorithms/Statistics/StatisticsConfidence.h"

using namespace MML;
using namespace MML::Statistics;

Vector<Real> data = {...};

// 95% CI for mean
auto ci = ConfidenceIntervalMean(data, 0.95);
std::cout << "95% CI: [" << ci.lowerBound << ", " << ci.upperBound << "]\n";

// CI for proportion
int successes = 45;
int trials = 100;
auto ciProp = ConfidenceIntervalProportion(successes, trials, 0.95);
```

### Correlation Analysis

```cpp
#include "algorithms/Statistics/StatisticsRank.h"

using namespace MML;
using namespace MML::Statistics;

Vector<Real> x({1.2, 2.3, 3.1, 4.0, 5.2});
Vector<Real> y({1.5, 2.1, 3.5, 4.2, 5.0});

// Pearson correlation (linear)
Real r = PearsonCorrelation(x, y);

// Spearman rank correlation (monotonic)
Real rho = SpearmanCorrelation(x, y);

// Kendall's tau (concordance)
Real tau = KendallCorrelation(x, y);

std::cout << "Pearson: " << r << "\n";
std::cout << "Spearman: " << rho << "\n";
std::cout << "Kendall: " << tau << "\n";
```

### Random Number Generation

```cpp
#include "RandomGenerators.h"

using namespace MML;
using namespace MML::Statistics;

NormalDeviate normal(0.0, 1.0, 42); // Fixed seed for reproducibility
Real z = normal.generate();

ExponentialDeviate exponential(0.5, 42);
Real waitingTime = exponential.generate();

PoissonDeviate poisson(4.0, 42);
int eventCount = poisson.generate();
```

### Histogram Analysis

```cpp
#include "Histogram.h"

using namespace MML;
using namespace MML::Statistics;

Vector<Real> data({0.2, 0.4, 0.7, 1.1, 1.8});

// Create histogram with 20 bins
Histogram::HistogramResult hist = Histogram::ComputeHistogram(data, 20);

// Or with specified bin edges
Vector<Real> binEdges({0.0, 0.5, 1.0, 2.0, 5.0});
Histogram::HistogramResult hist2 = Histogram::ComputeHistogram(data, binEdges);

// Get bin counts
const auto& counts = hist.counts;
const auto& edges = hist.binEdges;

// Normalized histogram (probability density)
const auto& density = hist.density;
```

## API Reference

### Distribution Interface

Continuous distributions provide `pdf(x)` and `cdf(x)`. Quantile and moment
operations depend on the distribution:

| Method | Description |
|--------|-------------|
| `pdf(x)` | Probability density function |
| `cdf(x)` | Cumulative distribution function P(X ≤ x) |
| `inverseCdf(p)` | Quantile function where available |
| `criticalValue(alpha)` | Upper-tail critical value for chi-square and F |
| `mean()`, `variance()` | Moments where mathematically defined |

### Hypothesis Test Result

```cpp
struct HypothesisTestResult {
    Real testStatistic;   // Computed test statistic
    Real pValue;          // P-value
    Real criticalValue;   // Critical value at significance level
    bool rejectNull;      // Decision at given α
    Real confidenceLevel; // 1 - α
    int degreesOfFreedom; // df if applicable
    std::string testName; // Test description
};
```

### Descriptive Functions

| Function | Description |
|----------|-------------|
| `Mean(data)` | Arithmetic mean |
| `Median(data)` | 50th percentile |
| `Variance(data)` | Sample variance |
| `StdDev(data)` | Sample standard deviation |
| `Moments(data, ...)` | Mean, deviations, skewness, and excess kurtosis |
| `Percentile(data, p)` | p-th percentile (0-100) |
| `IQR(data)` | Interquartile range |

## File Structure

```
statistics/
├── README.md              # This file
├── CMakeLists.txt         # Build configuration
│
├── include/
│   ├── DiscreteDistributions.h # Bernoulli, binomial, Poisson, and related distributions
│   ├── DataDescriptors.h      # Typed data analysis
│   │
│   ├── HypothesisTesting.h    # t-tests, chi-square tests
│   ├── ConfidenceIntervals.h  # CI computation
│   ├── RankCorrelation.h      # Spearman, Kendall
│   │
│   ├── RandomGenerators.h     # MT19937-based RNG
│   ├── Histogram.h            # Histogram analysis
│   ├── TimeSeries.h           # Basic time series
│   │
│   └── (algorithms/Statistics.h) # Core functions in main MML
│
└── tests/
    ├── datadescriptors_tests.cpp
    ├── histogram_tests.cpp
    └── hypothesis_tests.cpp
```

## Mathematical Background

### Normal Distribution

$$f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}$$

### Student's t Distribution

$$f(t) = \frac{\Gamma(\frac{\nu+1}{2})}{\sqrt{\nu\pi}\Gamma(\frac{\nu}{2})} \left(1+\frac{t^2}{\nu}\right)^{-\frac{\nu+1}{2}}$$

where $\nu$ is degrees of freedom.

### One-Sample t-Statistic

$$t = \frac{\bar{x} - \mu_0}{s / \sqrt{n}}$$

where $\bar{x}$ is sample mean, $s$ is sample standard deviation, $n$ is sample size.

### Confidence Interval for Mean

$$\bar{x} \pm t_{\alpha/2,n-1} \cdot \frac{s}{\sqrt{n}}$$

## Common Use Cases

| Task | Function/Class |
|------|---------------|
| Summarize data | `Mean()`, `StdDev()`, `Median()` |
| Test if mean differs from value | `OneSampleTTest()` |
| Compare two group means | `TwoSampleTTest()` |
| Test independence | `ChiSquareIndependence()` |
| Measure linear relationship | `PearsonCorrelation()` |
| Measure monotonic relationship | `SpearmanCorrelation()` |
| Generate normal samples | `NormalDeviate` |
| Create confidence interval | `ConfidenceIntervalMean()` |

## References

- Rice, J.A. (2006). "Mathematical Statistics and Data Analysis" (3rd ed.)
- Casella, G., and Berger, R.L. (2002). "Statistical Inference" (2nd ed.)
- Press, W.H., et al. (2007). "Numerical Recipes" (3rd ed.), Chapters 14-15
- NIST/SEMATECH e-Handbook of Statistical Methods

## See Also

- [mml_packages/README.md](../README.md) - Package overview
- [Optimization Package](../optimization/README.md) - Parameter estimation, fitting
- [MML algorithms/Statistics.h](../../../mml/algorithms/) - Core statistical functions
