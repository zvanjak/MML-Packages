<div align="center">

# 📦 MML Packages

### **Domain-Specific Numerical Libraries for the Minimal Math Library**

*Symbolic math • Optimization • PDEs • Fourier analysis • Inferential statistics • Graph analytics*

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/std/the-standard)
[![Packages](https://img.shields.io/badge/packages-6-orange.svg)](#-the-packages)
[![Built on MML](https://img.shields.io/badge/built%20on-MML%202.0-9cf.svg)](https://github.com/zvanjak/MML)
[![License](https://img.shields.io/badge/license-free%20personal%2Feducation%20%7C%20commercial%20paid-blue.svg)](LICENSE.md)

[![Windows released package](https://github.com/zvanjak/MML-Packages/actions/workflows/release-windows.yml/badge.svg)](https://github.com/zvanjak/MML-Packages/actions/workflows/release-windows.yml)
[![Linux released package](https://github.com/zvanjak/MML-Packages/actions/workflows/release-linux.yml/badge.svg)](https://github.com/zvanjak/MML-Packages/actions/workflows/release-linux.yml)
[![macOS released package](https://github.com/zvanjak/MML-Packages/actions/workflows/release-macos.yml/badge.svg)](https://github.com/zvanjak/MML-Packages/actions/workflows/release-macos.yml)

**Extend [MML Core](https://github.com/zvanjak/MML) with focused C++20 packages** — solve multi-objective optimization problems, march PDEs on structured grids, run spectral analysis, test hypotheses, build symbolic expressions, and analyze graphs.

[The Packages](#-the-packages) • [Quick Start](#-quick-start) • [Building](#-building--linking) • [Docs](#-documentation)

</div>

---

## 🎯 What is MML Packages?

**MML Packages** extends MML Core's vectors, matrices, calculus, and solvers with five domain libraries plus the `mml_ext` extension header tree. Use all packages through `MML::Packages`, or link only the target you need (`MML::Optimization`, `MML::PDE`, `MML::Fourier`, `MML::Statistics`, `MML::Symbolic`, or `MML::Ext`).

The public [releases](https://github.com/zvanjak/MML-Packages/releases) provide prebuilt static libraries, matching C++20 headers, a vendored MML Core snapshot, documentation demos, and relocatable CMake targets for Windows x64, Linux x86-64, and macOS arm64. This public repository is a **release companion**: it contains headers, documentation, examples, release metadata, and vendored MML Core headers, but it does **not** contain the private package implementation sources or root source-build project.

The platform badges report public consumer checks of a published release: each runner verifies its archive, builds the documentation apps against it, and runs both apps. They do not represent a source build of the private test suite. A badge appears once its workflow has run.

This is an early 0.1 series: APIs and supported environments may change before 1.0.

---

## 🧩 The Packages

| Package | What it does | Highlights |
|---------|--------------|------------|
| 🎯 [**Optimization**](docs/README_Optimization.md) | Derivative-free, stochastic, constrained, and multi-objective optimization | Genetic Algorithms, Simulated Annealing, **NSGA-II** & **MOEA/D**, penalty methods, Revised Simplex LP, lifecycle/observer APIs, benchmark problem suite |
| 🌊 [**PDE**](docs/README_PDE.md) | Partial differential equation solvers | Structured 1D/2D/3D grids, boundary conditions, stencils; Poisson, Heat, Wave & Advection; manufactured-solution test beds |
| 📈 [**Fourier**](docs/README_Fourier.md) | Spectral and signal analysis | Fourier series, spectrum analysis, convolution, correlation, orthogonal Fourier bases |
| 📊 [**Statistics**](docs/README_Statistics.md) | Data and inferential statistics | Data descriptors, random generators, time series, hypothesis tests, confidence intervals, rank correlation |
| 🧮 [**Symbolic**](docs/README_Symbolic.md) | Symbolic computing and automatic differentiation | Expression trees, simplification, automatic differentiation, expression parsing, code generation |
| 🔗 [**mml_ext**](include/mml_ext/) | Extensions to MML Core | Spectral graph analytics (PageRank, centralities) and field-line tracing; headers, not a sixth package library |

---

## 🚀 Quick Start

1. Download your platform archive and `SHA256SUMS` from [Releases](https://github.com/zvanjak/MML-Packages/releases). The current first public release is [`v0.1.0`](https://github.com/zvanjak/MML-Packages/releases/tag/v0.1.0).
2. Verify the archive's checksum, extract it, and use the extracted directory as your installation prefix. Keep its headers and libraries together.
3. Point CMake at that prefix and link `MML::Packages` for all five libraries, or an individual target such as `MML::Optimization`. The [binary installation guide](docs/Binary_Release.md) covers the exact toolchains, ABI requirements, and checksum commands.

For a taste of the API, this runs NSGA-II on a small ZDT1-style bi-objective optimization problem:

```cpp
#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <optimization/algorithms/NSGA2.h>
#include <optimization/core/Variables.h>

#include <cmath>
#include <iostream>

using namespace MML;
using namespace MML::Optimization;

int main() {
    struct ZDT1 {
        Vector<Real> Evaluate(const Vector<Real>& x) {
            Real f1 = x[0];
            Real sum = 0.0;
            for (int i = 1; i < x.size(); ++i)
                sum += x[i];

            Real g = 1.0 + 9.0 * sum / (x.size() - 1);
            Real f2 = g * (1.0 - std::sqrt(f1 / g));
            return Vector<Real>{ f1, f2 };
        }
    } problem;

    ProblemSpec spec = ProblemSpec::Continuous(10, 0.0, 1.0);

    NSGA2Config config;
    config.populationSize = 80;
    config.maxGenerations = 150;
    config.seed = 42;

    NSGA2 nsga(config);
    nsga.SetProblem(spec);

    auto result = nsga.Optimize(problem);
    std::cout << "Pareto solutions: " << result.paretoFront.size() << "\n";
}
```

Link the example with `MML::Optimization`. The [optimization guides](docs/README_Optimization.md) and [runnable documentation demos](docs_demos/optimization/) cover NSGA-II, MOEA/D, simulated annealing, genetic algorithms, penalty methods, and iterative lifecycle control in more detail.

---

## 🎯 Optimization Visualization Example

NSGA-II returns a set of non-dominated objective points. For two-objective problems, those points form an approximate Pareto front in the `(f1, f2)` objective plane. The ZDT1 quick-start problem above produces the characteristic smooth trade-off curve: improving `f1` forces `f2` upward, and no point on the front dominates another.

The screenshot below exports the sorted Pareto front as an MML Visualizers 2D parametric curve, with objective values scaled for readability.

| NSGA-II ZDT1 Pareto front |
|:-------------------------:|
| ![NSGA-II ZDT1 Pareto front](<docs/images/example01_Pareto_front/Screenshot 2026-10-08 100630.png>) |

---

## 🌊 PDE Visualization Example

The PDE package also ships finite-difference solvers that can feed MML Visualizers. The 2D Laplace demo solves a steady-state temperature distribution on a square plate with a sinusoidally heated top edge and fixed zero temperature on the remaining edges:

```text
∇²u = 0 on [0,1] × [0,1]
u(x,0) = 0
u(x,1) = sin(πx)
u(0,y) = 0
u(1,y) = 0
```

Source: [docs_demos/pde/example03_2d_laplace.cpp](docs_demos/pde/example03_2d_laplace.cpp). The same solution can be exported as an MML Visualizers scalar-field file, making the temperature surface easy to inspect interactively.

| Surface view | Rotated view | Height-field detail |
|:------------:|:------------:|:-------------------:|
| ![2D Laplace temperature surface](<docs/images/example02_temp_distribution/Screenshot 2026-10-08 094835.png>) | ![2D Laplace rotated temperature surface](<docs/images/example02_temp_distribution/Screenshot 2026-10-08 094848.png>) | ![2D Laplace height-field detail](<docs/images/example02_temp_distribution/Screenshot 2026-10-08 094903.png>) |

---

## 🛠️ Building & Linking

The public repository is a **release companion**, not a standalone source build: it contains headers, documentation, examples, and a vendored MML snapshot, but not the private library build sources or root CMake project. Use a matching [release archive](https://github.com/zvanjak/MML-Packages/releases) for compiled libraries. The installed `MMLPackages` config exports `MML::Packages`, five individual package targets, `MML::Core`, and `MML::Ext`.

The [documentation demo application](docs_demos/) can be built against an extracted release prefix:

```sh
cmake -S docs_demos -B build-docs -DMML_INSTALLED_PREFIX=/path/to/extracted/MML-Packages-prefix
cmake --build build-docs --config Release
```

Pass the actual extracted directory as `MML_INSTALLED_PREFIX`. On multi-configuration generators, run the executable from `build-docs/Release/`; on single-configuration generators, from `build-docs/`. Tested platform toolchains and runtime requirements are listed in the [installation guide](docs/Binary_Release.md), which currently describes the first prerelease.

## Relationship to MML Core

[MML Core](https://github.com/zvanjak/MML) supplies the general numerical foundation; these packages add domain-specific algorithms and solvers on top. Each release bundles a matching MML Core header snapshot, so consumers do not need to mix it with a separate MML checkout. Functionality that proves broadly useful may eventually move into MML Core.

## Documentation

- [Package documentation](docs/) covers each domain and its APIs.
- [Runnable documentation demos](docs_demos/) back the examples; start with the [PDE gallery](docs/pde/Examples_Gallery.md) or the [optimization guides](docs/README_Optimization.md).
- [Binary installation](docs/Binary_Release.md) covers toolchains, CMake targets, checksums, and the release manifest.
- Read the release notes and `share/mml-packages/release-manifest.json` that accompany each archive rather than assuming another release's toolchain details apply.

## License

MML Packages and the bundled MML Core are MIT licensed. Binary archives contain both license notices under `licenses/`; keep them with redistributed copies. The [MML Core notice](licenses/MML-Core-LICENSE.md) is also available in this checkout.