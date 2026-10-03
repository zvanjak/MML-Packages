# Fourier Package

**Signal Processing and Spectral Analysis**

The Fourier package provides comprehensive tools for frequency-domain analysis, including fast transforms, spectral estimation, window functions, and convolution operations.

## Features

### Transform Algorithms
- **FFT** - Cooley-Tukey radix-2 algorithm, O(n log n), power-of-2 sizes
- **DFT** - Reference O(n²) implementation for arbitrary sizes
- **Real FFT** - Optimized transform for real-valued signals (half storage)
- **DCT** - Discrete Cosine Transform (types I-IV)
- **Fourier Series** - Classical Fourier series expansion

### Spectral Analysis
- **Power Spectrum** - Direct periodogram estimation
- **Welch's Method** - Averaged periodogram with overlapping segments
- **Magnitude/Phase** - Spectrum decomposition
- **Frequency bins** - Automatic bin calculation for sample rates

### Window Functions
- **Rectangular** - No windowing (narrowest main lobe, -13 dB sidelobes)
- **Hann (Hanning)** - General-purpose (-31 dB sidelobes)
- **Hamming** - Audio processing (-43 dB sidelobes)
- **Blackman** - Low sidelobe applications (-58 dB)
- **Blackman-Harris** - Very low sidelobe (-92 dB)
- **Kaiser** - Adjustable β parameter
- **Gaussian** - Smooth frequency response
- **Flat-Top** - Amplitude accuracy

### Convolution & Correlation
- **Linear Convolution** - Standard discrete convolution via FFT
- **Circular Convolution** - Wraparound convolution
- **Cross-Correlation** - Signal similarity analysis
- **Auto-Correlation** - Signal self-similarity
- **Wiener Deconvolution** - Noise-robust inverse filtering

## Quick Start

### Basic FFT

```cpp
#include "mml_packages/fourier/include/FFT.h"

using namespace MML::Fourier;

// Create a signal
Vector<Complex> signal(1024);
for (int i = 0; i < 1024; i++) {
    double t = i / 1000.0;  // 1 kHz sample rate
    signal[i] = Complex(std::sin(2 * M_PI * 50 * t), 0.0);  // 50 Hz sine
}

// Forward FFT (time → frequency)
auto spectrum = FFT::Forward(signal);

// Inverse FFT (frequency → time)
auto reconstructed = FFT::Inverse(spectrum);
```

### Power Spectrum Analysis

```cpp
#include "mml_packages/fourier/include/Spectrum.h"
#include "mml_packages/fourier/include/Windowing.h"

using namespace MML::Fourier;

// Real-valued signal
Vector<Real> data = {...};

// Apply window
auto windowed = Windows::Hamming(data.size());

// Compute power spectrum
auto power = PowerSpectrum::Compute(data, windowed);

// Welch's method for reduced variance
auto psd = PowerSpectrum::Welch(data, 256, 128);  // segment=256, overlap=128
```

### Convolution and Filtering

```cpp
#include "mml_packages/fourier/include/Convolution.h"

using namespace MML::Fourier;

Vector<Real> signal = {...};    // Input signal
Vector<Real> kernel = {...};    // Filter kernel (e.g., low-pass)

// Linear convolution (output length = n + m - 1)
auto filtered = Convolution::Linear(signal, kernel);

// Cross-correlation (for pattern matching)
auto corr = Correlation::CrossCorrelation(signal, template_signal);
```

### Window Functions

```cpp
#include "mml_packages/fourier/include/Windowing.h"

using namespace MML::Fourier;

int N = 1024;

// Generate windows
auto hann = Windows::Hann(N);
auto hamming = Windows::Hamming(N);
auto blackman = Windows::Blackman(N);
auto kaiser = Windows::Kaiser(N, 5.0);  // β = 5.0

// Apply window to signal
for (int i = 0; i < N; i++) {
    windowed[i] = signal[i] * hamming[i];
}
```

### DCT (Discrete Cosine Transform)

```cpp
#include "mml_packages/fourier/include/DCT.h"

using namespace MML::Fourier;

Vector<Real> data = {...};

// DCT-II (most common, used in JPEG/MP3)
auto coefficients = DCT::Forward(data);

// Inverse DCT (DCT-III)
auto reconstructed = DCT::Inverse(coefficients);
```

## API Reference

### FFT Class

| Method | Description |
|--------|-------------|
| `Forward(data)` | Forward FFT, returns frequency-domain signal |
| `Inverse(spectrum)` | Inverse FFT, returns time-domain signal |
| `Transform(data, isign)` | In-place FFT (isign: +1=forward, -1=inverse) |
| `IsPowerOfTwo(n)` | Check if n is a power of 2 |
| `NextPowerOfTwo(n)` | Find next power of 2 ≥ n |
| `ZeroPad(data)` | Zero-pad to next power of 2 |

### DFT Class

| Method | Description |
|--------|-------------|
| `Forward(data)` | Forward DFT (any size), O(n²) |
| `Inverse(spectrum)` | Inverse DFT, O(n²) |

### PowerSpectrum Class

| Method | Description |
|--------|-------------|
| `Compute(data, window)` | Periodogram (direct method) |
| `Welch(data, segment, overlap)` | Welch's method (reduced variance) |
| `Magnitude(spectrum)` | Magnitude |X[k]| |
| `Phase(spectrum)` | Phase arg(X[k]) |

### Windows Namespace

| Function | Side Lobe Level | Main Lobe Width | Use Case |
|----------|-----------------|-----------------|----------|
| `Rectangular(n)` | -13 dB | Narrowest | Exact periodic signals |
| `Hann(n)` | -31 dB | Moderate | General purpose |
| `Hamming(n)` | -43 dB | Moderate | Audio processing |
| `Blackman(n)` | -58 dB | Wide | Low leakage |
| `BlackmanHarris(n)` | -92 dB | Very wide | Precision spectral analysis |
| `Kaiser(n, beta)` | Adjustable | Adjustable | Flexible trade-off |
| `Gaussian(n, sigma)` | Smooth | Gaussian | Smooth transitions |
| `FlatTop(n)` | -44 dB | Very wide | Amplitude accuracy |

### Convolution Class

| Method | Description |
|--------|-------------|
| `Linear(signal, kernel)` | Linear convolution (output: n+m-1) |
| `Circular(signal, kernel)` | Circular convolution (output: max(n,m)) |
| `Wiener(blurred, kernel, snr)` | Wiener deconvolution |

### Correlation Class

| Method | Description |
|--------|-------------|
| `Auto(signal)` | Autocorrelation |
| `Cross(signal1, signal2)` | Cross-correlation |
| `Normalized(signal1, signal2)` | Normalized cross-correlation |

## File Structure

```
fourier/
├── README.md              # This file
├── CMakeLists.txt         # Build configuration
│
├── include/
│   ├── FFT.h              # Fast Fourier Transform (radix-2)
│   ├── DFT.h              # Discrete Fourier Transform (reference)
│   ├── RealFFT.h          # Optimized FFT for real signals
│   ├── DCT.h              # Discrete Cosine Transform
│   ├── FourierSeries.h    # Classical Fourier series
│   ├── FourierBasis.h     # Basis function utilities
│   │
│   ├── Spectrum.h         # Power spectrum estimation
│   ├── Windowing.h        # Window functions
│   ├── Convolution.h      # Convolution operations
│   ├── Correlation.h      # Correlation analysis
│   │
│   └── FourierValidation.h # Input validation utilities
│
└── tests/
    └── fourier_tests.cpp  # Comprehensive test suite
```

## Mathematical Background

### Discrete Fourier Transform

**Forward DFT:**
$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-2\pi i \cdot kn/N}$$

**Inverse DFT:**
$$x[n] = \frac{1}{N} \sum_{k=0}^{N-1} X[k] \cdot e^{2\pi i \cdot kn/N}$$

### Convolution Theorem

The DFT converts convolution in time domain to multiplication in frequency domain:

$$\mathcal{F}\{x * y\} = \mathcal{F}\{x\} \cdot \mathcal{F}\{y\}$$

This allows O(n log n) convolution via:
1. Compute FFT of both signals
2. Element-wise multiplication
3. Compute inverse FFT

### Window Functions

Windows reduce spectral leakage by tapering signal edges:

$$w[n] = \text{window}(n, N) \quad \text{for } n = 0, 1, \ldots, N-1$$

Applied signal: $x_w[n] = x[n] \cdot w[n]$

## Performance Characteristics

| Operation | Complexity | Notes |
|-----------|------------|-------|
| FFT | O(n log n) | Radix-2, requires power-of-2 size |
| DFT | O(n²) | Any size, reference implementation |
| Real FFT | O(n log n) | Half storage for real signals |
| Convolution (FFT) | O(n log n) | ~100x faster than direct for n=1024 |
| Convolution (Direct) | O(nm) | Simple but slow |

## References

- Cooley, J.W., and Tukey, J.W. (1965). "An Algorithm for the Machine Calculation of Complex Fourier Series"
- Harris, F.J. (1978). "On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform"
- Welch, P.D. (1967). "The Use of Fast Fourier Transform for the Estimation of Power Spectra"
- Numerical Recipes in C++, 3rd Edition, Chapter 13: Fourier and Spectral Applications

## See Also

- [mml_packages/README.md](../README.md) - Package overview
- [Optimization Package](../optimization/README.md) - Uses Fourier for frequency-domain optimization
- [PDE Package](../pde/README.md) - Spectral methods for PDE solving
