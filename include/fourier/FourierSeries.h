///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FourierSeries.h                                                     ///
///  Description: Fourier series expansion and coefficient computation                ///
///               Periodic function approximation via trigonometric series            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FOURIER_SERIES_H
#define MML_FOURIER_SERIES_H

#include "mml/mml_export.h"

#include "MMLBase.h"
#include <mml/algorithms/Fourier/Fourier.h>

#include "interfaces/IFunction.h"
#include "base/Vector/Vector.h"
#include <mml/core/Integration/Integration1D.h>

#include <cmath>

namespace MML::Fourier 
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// FourierSeries - Continuous Fourier series representation
	//
	// Represents a periodic function as:
	//   f(x) = a₀/2 + Σ[aₙ·cos(nπx/L) + bₙ·sin(nπx/L)]  for n = 1, 2, ..., N
	//
	// where:
	//   - L is the half-period (function is 2L-periodic)
	//   - aₙ are cosine coefficients
	//   - bₙ are sine coefficients
	//
	// Key features:
	//   - Approximates periodic functions with spectral convergence (for smooth functions)
	//   - Closed-form derivatives and integrals
	//   - Parseval's theorem for energy computation
	//   - Conversion to/from complex exponential form
	//
	// Usage:
	//   auto sin_func = [](Real x) { return std::sin(x); };
	//   FourierSeries fs(sin_func, Constants::PI, 10);  // Approximate sin(x) on [-π, π]
	//   Real value = fs(0.5);  // Evaluate at x=0.5
	///////////////////////////////////////////////////////////////////////////////////////////
	class MML_FOURIER_API FourierSeries : public IRealFunction {
	private:
		Vector<Real> _a; // Cosine coefficients [a₀, a₁, a₂, ..., aₙ]
		Vector<Real> _b; // Sine coefficients [b₁, b₂, ..., bₙ]  (b₀ = 0 by definition)
		Real _L;				 // Half-period (function has period 2L)
		int _N;					 // Number of terms

	public:
		///////////////////////////////////////////////////////////////////////////////////////////
		// Constructors
		///////////////////////////////////////////////////////////////////////////////////////////

		// Default constructor - zero series
		FourierSeries()
				: _L(Constants::PI)
				, _N(0) {}

		// Construct from coefficients
		FourierSeries(const Vector<Real>& a, const Vector<Real>& b, Real L)
				: _a(a)
				, _b(b)
				, _L(L) {
			FourierValidation::ValidateNonEmpty(static_cast<int>(_a.size()), "FourierSeries (coefficient a)");
			if (_b.size() != _a.size() - 1)
				throw std::invalid_argument("FourierSeries: _b must have size = _a.size() - 1");

			_N = static_cast<int>(_a.size()) - 1; // N terms (excluding a₀)
		}

		// Construct from function by computing Fourier coefficients
		// Uses adaptive quadrature for coefficient computation
		FourierSeries(const IRealFunction& func, Real L, int N, Real eps = 1e-10);

		///////////////////////////////////////////////////////////////////////////////////////////
		// Evaluation
		///////////////////////////////////////////////////////////////////////////////////////////

		// Evaluate Fourier series at x
		Real operator()(Real x) const override {
			if (_N == 0)
				return _a[0] / 2.0;

			Real sum = _a[0] / 2.0;

			for (int n = 1; n <= _N; n++) {
				Real arg = n * Constants::PI * x / _L;
				sum += _a[n] * std::cos(arg) + _b[n - 1] * std::sin(arg);
			}

			return sum;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Calculus Operations (closed-form for Fourier series!)
		///////////////////////////////////////////////////////////////////////////////////////////

		// Derivative: d/dx[aₙcos(nπx/L) + bₙsin(nπx/L)]
		//           = -aₙ(nπ/L)sin(nπx/L) + bₙ(nπ/L)cos(nπx/L)
		// Maps to: a'ₙ = bₙ(nπ/L), b'ₙ = -aₙ(nπ/L)
		FourierSeries Derivative() const;

		// Integral: ∫[aₙcos(nπx/L) + bₙsin(nπx/L)]dx
		//         = aₙ(L/nπ)sin(nπx/L) - bₙ(L/nπ)cos(nπx/L)
		// Maps to: A'ₙ = -bₙ(L/nπ), B'ₙ = aₙ(L/nπ)
		// Note: a₀ term integrates to (a₀/2)x, but we drop it (indefinite integral constant)
		FourierSeries Integral() const;

		///////////////////////////////////////////////////////////////////////////////////////////
		// Coefficient Access
		///////////////////////////////////////////////////////////////////////////////////////////

		const Vector<Real>& CosineCoefficients() const { return _a; }
		const Vector<Real>& SineCoefficients() const { return _b; }
		Real HalfPeriod() const { return _L; }
		Real Period() const { return 2.0 * _L; }
		int NumTerms() const { return _N; }

		///////////////////////////////////////////////////////////////////////////////////////////
		// Complex Exponential Form
		///////////////////////////////////////////////////////////////////////////////////////////

		// Convert to complex exponential form: f(x) = Σ cₙ·e^(inπx/L)
		// Relationship: cₙ = (aₙ - i·bₙ)/2 for n > 0
		//               c₀ = a₀/2
		//               c₋ₙ = (aₙ + i·bₙ)/2 = c̄ₙ (complex conjugate)
		Vector<Complex> ComplexCoefficients() const;

		// Construct from complex coefficients
		// Input: c = [c₋ₙ, ..., c₋₁, c₀, c₁, ..., cₙ] (length 2N+1)
		static FourierSeries FromComplexCoefficients(const Vector<Complex>& c, Real L);

		///////////////////////////////////////////////////////////////////////////////////////////
		// Energy and Norms (Parseval's Theorem)
		///////////////////////////////////////////////////////////////////////////////////////////

		// Parseval's theorem: ∫[-L,L] |f(x)|² dx = L·(a₀²/2 + Σ(aₙ² + bₙ²))
		Real Energy() const {
			Real energy = (_a[0] * _a[0]) / 2.0; // a₀²/2

			for (int n = 1; n <= _N; n++) {
				energy += _a[n] * _a[n] + _b[n - 1] * _b[n - 1];
			}

			return _L * energy;
		}

		// L² norm: ||f||₂ = √(Energy/2L)
		Real L2Norm() const { return std::sqrt(Energy() / (2.0 * _L)); }

		///////////////////////////////////////////////////////////////////////////////////////////
		// Series Operations
		///////////////////////////////////////////////////////////////////////////////////////////

		// Truncate to M terms (M < N)
		FourierSeries Truncate(int M) const {
			FourierValidation::ValidateNonNegative(M, "M", "FourierSeries::Truncate");
			if (M > _N)
				return *this;

			Vector<Real> new_a(M + 1);
			Vector<Real> new_b(M);

			for (int i = 0; i <= M; i++)
				new_a[i] = _a[i];

			for (int i = 0; i < M; i++)
				new_b[i] = _b[i];

			return FourierSeries(new_a, new_b, _L);
		}

		// Addition
		FourierSeries operator+(const FourierSeries& other) const {
			if (std::abs(_L - other._L) > 1e-10)
				throw std::invalid_argument("FourierSeries::operator+: periods must match");

			int max_N = std::max(_N, other._N);
			Vector<Real> new_a(max_N + 1, 0.0);
			Vector<Real> new_b(max_N, 0.0);

			for (int i = 0; i <= _N; i++)
				new_a[i] = _a[i];
			for (int i = 0; i < _N; i++)
				new_b[i] = _b[i];

			for (int i = 0; i <= other._N; i++)
				new_a[i] += other._a[i];
			for (int i = 0; i < other._N; i++)
				new_b[i] += other._b[i];

			return FourierSeries(new_a, new_b, _L);
		}

		// Subtraction
		FourierSeries operator-(const FourierSeries& other) const {
			if (std::abs(_L - other._L) > 1e-10)
				throw std::invalid_argument("FourierSeries::operator-: periods must match");

			int max_N = std::max(_N, other._N);
			Vector<Real> new_a(max_N + 1, 0.0);
			Vector<Real> new_b(max_N, 0.0);

			for (int i = 0; i <= _N; i++)
				new_a[i] = _a[i];
			for (int i = 0; i < _N; i++)
				new_b[i] = _b[i];

			for (int i = 0; i <= other._N; i++)
				new_a[i] -= other._a[i];
			for (int i = 0; i < other._N; i++)
				new_b[i] -= other._b[i];

			return FourierSeries(new_a, new_b, _L);
		}

		// Scalar multiplication
		FourierSeries operator*(Real scalar) const {
			Vector<Real> new_a(_N + 1);
			Vector<Real> new_b(_N);

			for (int i = 0; i <= _N; i++)
				new_a[i] = _a[i] * scalar;
			for (int i = 0; i < _N; i++)
				new_b[i] = _b[i] * scalar;

			return FourierSeries(new_a, new_b, _L);
		}

		friend FourierSeries operator*(Real scalar, const FourierSeries& fs) { return fs * scalar; }

		///////////////////////////////////////////////////////////////////////////////////////////
		// Convergence Analysis
		///////////////////////////////////////////////////////////////////////////////////////////

		// Get magnitude of nth coefficient (for convergence analysis)
		Real CoefficientMagnitude(int n) const {
			if (n == 0)
				return std::abs(_a[0]);
			if (n > _N)
				return 0.0;

			return std::sqrt(_a[n] * _a[n] + _b[n - 1] * _b[n - 1]);
		}

		// Maximum coefficient magnitude (indicator of convergence)
		Real MaxCoefficientMagnitude() const {
			Real max_mag = std::abs(_a[0]);

			for (int n = 1; n <= _N; n++) {
				Real mag = std::sqrt(_a[n] * _a[n] + _b[n - 1] * _b[n - 1]);
				if (mag > max_mag)
					max_mag = mag;
			}

			return max_mag;
		}
	};

} // namespace MML::Fourier
#endif // MML_FOURIER_SERIES_H
