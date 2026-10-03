///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Correlation.h                                                       ///
///  Description: Cross-correlation and auto-correlation via FFT                      ///
///               Signal similarity analysis and lag detection                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_CORRELATION_H
#define MML_CORRELATION_H

#include "mml/mml_export.h"

#include <mml/algorithms/Fourier/Fourier.h>

#include "base/Vector/Vector.h"

#include <complex>
#include <algorithm>
#include <stdexcept>

namespace MML::Fourier 
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// Correlation - Fast correlation using FFT
	//
	// CORRELATION THEOREM:
	//   corr(x,y)[k] = Σ x[i] * conj(y[i+k])  (discrete cross-correlation)
	//   corr(x,y) = IFFT(FFT(x) ⊙ conj(FFT(y)))
	//
	// COMPLEXITY: O(n log n) using FFT vs O(n²) naive
	//
	// APPLICATIONS:
	//   - Signal processing (template matching, echo detection)
	//   - Time series analysis (lag detection, similarity)
	//   - Image processing (pattern recognition)
	//   - Neuroscience (spike train analysis)
	///////////////////////////////////////////////////////////////////////////////////////////
	class MML_FOURIER_API Correlation {
	public:
		///////////////////////////////////////////////////////////////////////////////////////////
		// Cross - Cross-correlation of two real signals
		//
		// OUTPUT LENGTH: n + m - 1 (same as linear convolution)
		//   - Negative lags: [0, m-2]
		//   - Zero lag: [m-1]
		//   - Positive lags: [m, n+m-2]
		//
		// ALGORITHM:
		//   Correlation = Convolution with time-reversed signal
		//   1. Reverse y → y_rev
		//   2. Convolve(x, y_rev) using FFT
		//
		// NOTE: Result is NOT symmetric unless x == y
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> Cross(const Vector<Real>& x, const Vector<Real>& y);

		///////////////////////////////////////////////////////////////////////////////////////////
		// Cross (complex version)
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Complex> Cross(const Vector<Complex>& x, const Vector<Complex>& y);

		///////////////////////////////////////////////////////////////////////////////////////////
		// Auto - Auto-correlation (correlation of signal with itself)
		//
		// PROPERTIES:
		//   - Symmetric around zero lag
		//   - Maximum at zero lag (lag = n-1 in output)
		//   - Decays as lag increases (for non-periodic signals)
		//
		// APPLICATIONS:
		//   - Periodicity detection
		//   - Signal energy estimation
		//   - Power spectral density (Wiener-Khinchin theorem)
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> Auto(const Vector<Real>& x) { return Cross(x, x); }

		///////////////////////////////////////////////////////////////////////////////////////////
		// CrossNormalized - Normalized cross-correlation coefficient
		//
		// RANGE: [-1, 1]
		//   -1: Perfect anti-correlation
		//    0: No correlation
		//   +1: Perfect correlation
		//
		// FORMULA:
		//   ρ[k] = corr(x,y)[k] / sqrt(E_x * E_y)
		//   where E_x = Σ x[i]², E_y = Σ y[i]²
		//
		// APPLICATIONS:
		//   - Template matching (find best lag for similarity)
		//   - Signal comparison (measure similarity independent of scale)
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> CrossNormalized(const Vector<Real>& x, const Vector<Real>& y);

		///////////////////////////////////////////////////////////////////////////////////////////
		// FindPeakLag - Find lag with maximum correlation
		//
		// RETURNS: Lag index where correlation is maximum
		//
		// NOTE: For output of length n+m-1:
		//   - Lag = 0 is at index n-1
		//   - Negative lags: indices [0, n-2]
		//   - Positive lags: indices [n, n+m-2]
		///////////////////////////////////////////////////////////////////////////////////////////
		static int FindPeakLag(const Vector<Real>& correlation) {
			FourierValidation::ValidateNonEmpty(correlation.size(), "Correlation::FindPeakLag");

			int max_idx = 0;
			Real max_val = correlation[0];
			for (int i = 1; i < correlation.size(); i++) {
				if (correlation[i] > max_val) {
					max_val = correlation[i];
					max_idx = i;
				}
			}

			return max_idx;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// ConvertLagToIndex - Convert lag value to index in correlation output
		//
		// For signals x[n] and y[m]:
		//   - Output length: n + m - 1
		//   - Zero lag at: index n - 1
		//   - Negative lag k: index (n-1) - k
		//   - Positive lag k: index (n-1) + k
		///////////////////////////////////////////////////////////////////////////////////////////
		static int ConvertLagToIndex(int lag, int x_size) { return (x_size - 1) + lag; }
	};

} // namespace MML::Fourier
#endif // MML_CORRELATION_H
