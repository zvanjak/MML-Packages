# MML Packages

**Extended Mathematical Computing Libraries**

The MML packages are specialized domain libraries built on top of the MinimalMathLibrary core. Each package provides production-ready implementations for specific mathematical computing domains, from signal processing to partial differential equations.

## Overview

| Package | Description | Key Features |
|---------|-------------|--------------|
| [**fourier**](fourier/README.md) | Signal Processing & Spectral Analysis | FFT, DFT, DCT, windowing, convolution, filtering |
| [**optimization**](optimization/README.md) | Optimization Algorithms | Nelder-Mead, BFGS, genetic algorithms, multi-objective |
| [**pde**](pde/README.md) | Partial Differential Equations | Elliptic, parabolic, hyperbolic solvers; sparse matrices |
| [**statistics**](statistics/README.md) | Statistical Computing | Distributions, hypothesis testing, time series |
| [**symbolic**](symbolic/README.md) | Symbolic Mathematics | Automatic differentiation, expression trees, code generation |

## Quick Start

Each package is a self-contained CMake project that integrates seamlessly with MML:

```cpp
// Include the package you need
#include "mml_packages/fourier/include/FFT.h"
#include "mml_packages/optimization/include/OptimizationMultidim.h"
#include "mml_packages/pde/include/PDE.h"
#include "mml_packages/statistics/include/Distributions.h"
#include "mml_packages/symbolic/include/Symbolic.h"

// Systems package is now in core MML
#include "systems/DynamicalSystem.h"
#include "systems/LinearSystem.h"
```

## Package Summaries

### Fourier - Signal Processing

Comprehensive Fourier analysis toolkit:
- **FFT/IFFT** - Cooley-Tukey radix-2 algorithm O(n log n)
- **DFT** - Reference O(n²) implementation for any size
- **DCT** - Discrete Cosine Transform (types I-IV)
- **Windowing** - Hamming, Hanning, Blackman, Kaiser, Gaussian
- **Spectrum Analysis** - Power spectrum, phase, magnitude
- **Convolution & Correlation** - Linear and circular

```cpp
using namespace MML::Fourier;

// Compute FFT
Vector<Complex> signal = {...};
auto spectrum = FFT::Forward(signal);

// Apply window and compute power spectrum
auto windowed = Windowing::Hamming(signal);
auto power = Spectrum::PowerSpectrum(windowed);
```

### Optimization - Numerical Optimization

Gradient-free and gradient-based optimization:
- **Nelder-Mead** - Simplex method for unconstrained minimization
- **Powell** - Direction set method
- **Conjugate Gradient** - For large-scale problems
- **BFGS** - Quasi-Newton method
- **Levenberg-Marquardt** - Nonlinear least squares
- **Genetic Algorithms** - Global optimization
- **NSGA-II/MOEA-D** - Multi-objective optimization

```cpp
using namespace MML::Optimization;

// Find minimum of Rosenbrock function
auto result = NelderMead<2>::Minimize(
    [](const VectorN<Real, 2>& x) {
        return 100*POW2(x[1] - x[0]*x[0]) + POW2(1 - x[0]);
    },
    {-2.0, 2.0},  // initial guess
    1e-8          // tolerance
);
```

### PDE - Partial Differential Equations

Complete PDE solving framework:
- **Sparse Matrices** - COO, CSR, CSC formats
- **Iterative Solvers** - CG, BiCGSTAB, GMRES
- **Preconditioners** - Jacobi, SSOR, ILU(0)
- **Grid Infrastructure** - 1D/2D/3D uniform grids
- **Elliptic PDEs** - Poisson, Laplace solvers
- **Parabolic PDEs** - Heat/diffusion equation
- **Hyperbolic PDEs** - Wave, advection equations

```cpp
using namespace MML::PDE;

// Solve 2D Poisson equation: -∇²u = f
PoissonSolver2D solver(grid);
solver.setBoundaryConditions(DirichletBC::Zero());
solver.setSource(source_function);
auto solution = solver.solve();
```

### Statistics - Statistical Computing

Probability and statistical analysis:
- **Distributions** - Normal, t, χ², F, Cauchy, Exponential, Logistic
- **Descriptive Statistics** - Mean, variance, skewness, kurtosis
- **Hypothesis Testing** - t-tests, χ² tests, ANOVA
- **Confidence Intervals** - For means, proportions, regression
- **Correlation** - Pearson, Spearman, Kendall
- **Time Series** - Autocorrelation, moving averages, ARIMA basics
- **Random Generation** - MT19937-based generators

```cpp
using namespace MML::Statistics;

// Statistical analysis
std::vector<double> data = {...};
auto stats = DataDescriptors::Compute(data);
std::cout << "Mean: " << stats.mean << ", SD: " << stats.stddev << "\n";

// Hypothesis testing
auto ttest = HypothesisTesting::OneSampleT(data, 100.0);  // H₀: μ = 100
```

### Symbolic - Symbolic Computing & Automatic Differentiation

Three approaches to computing derivatives:

1. **Symbolic Differentiation** - Build and manipulate expression trees
2. **Forward-Mode AD** - Dual numbers for few inputs → many outputs
3. **Reverse-Mode AD** - Tape-based backpropagation (ML-style gradients)

```cpp
using namespace MML::Symbolic;

// Symbolic differentiation
auto x = symbol("x");
auto f = sin_(x * x);              // f = sin(x²)
auto df = f->diff("x");            // df/dx = 2x·cos(x²)

// Forward-mode AD
Dual<double> x(2.0, 1.0);
Dual<double> y = sin(x * x);       // value and derivative in one pass

// Reverse-mode AD (gradient computation)
Tape tape;
ADVar a(&tape, 2.0), b(&tape, 3.0);
ADVar z = sin(a * b) + exp(a);
z.backward();                      // Compute ∂z/∂a and ∂z/∂b
```

### Systems - Dynamical Systems Analysis

Tools for nonlinear dynamics:
- **Fixed Point Analysis** - Location, stability, classification
- **Lyapunov Exponents** - Full spectrum computation
- **Bifurcation Diagrams** - Parameter sweeps
- **Phase Space Tools** - Poincaré sections, return maps
- **Stability Classification** - Nodes, foci, saddles, centers

```cpp
using namespace MML::Systems;

// Analyze the Lorenz system
LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);  // σ, ρ, β

// Find and classify fixed points
auto fixedPoints = DynamicalSystemAnalyzer::FindFixedPoints(lorenz);
for (auto& fp : fixedPoints) {
    std::cout << ToString(fp.type) << " at " << fp.location << "\n";
}

// Compute Lyapunov exponents
auto lyapunov = DynamicalSystemAnalyzer::ComputeLyapunovSpectrum(lorenz, x0, T);
std::cout << "λmax = " << lyapunov.maxExponent << "\n";
```

## Directory Structure

```
mml_packages/
├── README.md              # This file
├── CMakeLists.txt         # Package integration
├── tests/                 # Cross-package integration tests
│
├── fourier/               # Signal processing
│   ├── include/           # Headers (FFT, DFT, DCT, etc.)
│   ├── tests/             # Fourier-specific tests
│   └── CMakeLists.txt
│
├── optimization/          # Optimization algorithms
│   ├── include/           # Headers (multidim, LP, GA, etc.)
│   │   └── test_problems/ # Standard test functions
│   ├── tests/             # Optimization tests
│   └── CMakeLists.txt
│
├── pde/                   # PDE solvers
│   ├── include/           # Headers organized by type
│   │   ├── sparse/        # Sparse matrix formats
│   │   ├── solvers/       # Iterative solvers
│   │   ├── grid/          # Grid infrastructure
│   │   ├── elliptic/      # Poisson, Laplace
│   │   ├── parabolic/     # Heat equation
│   │   ├── hyperbolic/    # Wave, advection
│   │   └── testbeds/      # Validation tools
│   ├── tests/             # PDE solver tests
│   └── CMakeLists.txt
│
├── statistics/            # Statistical computing
│   ├── include/           # Headers (distributions, tests, etc.)
│   ├── tests/             # Statistics tests
│   └── CMakeLists.txt
│
├── symbolic/              # Symbolic computing
│   ├── include/           # Headers (Expr, AD, code gen)
│   ├── tests/             # Symbolic tests
│   └── CMakeLists.txt
│
└── systems/               # Dynamical systems
    ├── include/           # Headers (analysis, Lyapunov, etc.)
    ├── tests/             # Systems tests
    └── CMakeLists.txt
```

## Building

The packages are built as part of the main MML build:

```bash
cmake -B build
cmake --build build
ctest --test-dir build  # Run all tests including package tests
```

Individual packages can also be built independently if their dependencies are met.

## Dependencies

All packages depend on the MML core library:
- `mml/base/` - Vector, Matrix, Complex
- `mml/core/` - Linear algebra solvers
- `mml/algorithms/` - ODE integrators, root finding, etc.

Some packages have inter-dependencies:
- `pde` uses sparse matrix solvers internally
- `symbolic` adapters integrate with MML function interfaces
- `systems` uses ODE integrators from `algorithms`

## Testing

Each package includes comprehensive tests:

```bash
# Run all package tests
ctest --test-dir build -R "mml_packages"

# Run specific package tests
ctest --test-dir build -R "fourier"
ctest --test-dir build -R "optimization"
ctest --test-dir build -R "pde"
ctest --test-dir build -R "statistics"
ctest --test-dir build -R "symbolic"
ctest --test-dir build -R "systems"
```

## License

Licensed under the [MIT License](../LICENSE.md).

## See Also

- [MML Documentation](../docs/)
- [Examples](../examples/)
- [Main README](../README.md)