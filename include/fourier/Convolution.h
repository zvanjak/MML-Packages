///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Convolution.h                                                       ///
///  Description: Fast convolution and deconvolution via FFT                          ///
///               Circular and linear convolution, Wiener deconvolution               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CONVOLUTION_H
#define MML_CONVOLUTION_H

#include "mml/mml_export.h"

#include "MMLBase.h"
#include <mml/algorithms/Fourier/Fourier.h>
#include <mml/algorithms/Fourier/FourierConvolution.h>

namespace MML::Fourier 
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// Convolution - Fast convolution and deconvolution via FFT
	//
	// CONVOLUTION THEOREM:
	//   x * y = IFFT(FFT(x) ⊙ FFT(y))
	//
	// PERFORMANCE:
	//   - Direct convolution: O(n²)
	//   - FFT convolution: O(n log n)
	//   - Speedup: ~100x for n=1024
	//
	// REFERENCE: Numerical Recipes convlv.h
	///////////////////////////////////////////////////////////////////////////////////////////

	class MML_FOURIER_API Convolution 
  {
	public:
		///////////////////////////////////////////////////////////////////////////////////////////
		// Linear - Linear convolution (standard discrete convolution)
		//
		// OUTPUT LENGTH: n + m - 1
		//   - For signal of length n and kernel of length m
		//   - Produces complete convolution with no wraparound
		//
		// ALGORITHM:
		//   1. Zero-pad both to next power of 2 >= n+m-1
		//   2. FFT both
		//   3. Element-wise multiply
		//   4. IFFT
		//   5. Extract result (first n+m-1 samples)
		//
		// APPLICATIONS:
		//   - Filtering (FIR filters)
		//   - Smoothing
		//   - Edge detection
		//   - Pattern matching
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> Linear(const Vector<Real>& signal, const Vector<Real>& kernel);

		///////////////////////////////////////////////////////////////////////////////////////////
		// Linear (complex version)
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Complex> Linear(const Vector<Complex>& x, const Vector<Complex>& y);

		///////////////////////////////////////////////////////////////////////////////////////////
		// Circular - Circular convolution (periodic convolution)
		//
		// OUTPUT LENGTH: max(n, m)
		//   - Same length as longer input
		//   - Wraparound at boundaries (periodic)
		//
		// ALGORITHM:
		//   1. Pad both to same length (next power of 2)
		//   2. FFT both
		//   3. Element-wise multiply
		//   4. IFFT
		//
		// APPLICATIONS:
		//   - Periodic signals
		//   - Spectral analysis
		//   - Fast polynomial multiplication (mod x^n - 1)
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> Circular(const Vector<Real>& x, const Vector<Real>& y);

		///////////////////////////////////////////////////////////////////////////////////////////
		// Deconvolve - Inverse filtering (deconvolution)
		//
		// PROBLEM: Given y = x * h, recover x given y and h
		// SOLUTION: x = IFFT(FFT(y) / FFT(h))
		//
		// WARNING:
		//   - Unstable if kernel has small frequency components
		//   - May amplify noise
		//   - Consider regularization for real applications
		//
		// PARAMETERS:
		//   signal - Convolved signal y
		//   kernel - Known kernel h
		//   epsilon - Regularization parameter (prevents division by zero)
		//
		// APPLICATIONS:
		//   - Image deblurring
		//   - Channel equalization
		//   - System identification
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> Deconvolve(const Vector<Real>& signal, const Vector<Real>& kernel, Real epsilon = 1e-10);

		///////////////////////////////////////////////////////////////////////////////////////////
		// Deconvolve (complex version)
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Complex> Deconvolve(const Vector<Complex>& signal, const Vector<Complex>& kernel, Real epsilon = 1e-10);
	};

} // namespace MML::Fourier
#endif // MML_CONVOLUTION_H
