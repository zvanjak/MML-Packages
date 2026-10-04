<div align="center">

# MML Packages

### Domain-specific numerical libraries for the Minimal Math Library

*Symbolic math / Optimization / PDEs / Fourier analysis / Statistics / MML extensions*

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/std/the-standard)
[![Libraries](https://img.shields.io/badge/libraries-5-orange.svg)](#the-packages)
[![Platforms](https://img.shields.io/badge/platforms-Windows%20%7C%20Linux%20%7C%20macOS-brightgreen.svg)](#quick-start)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](#license)

**Build on [MinimalMathLibrary (MML)](https://github.com/zvanjak/MinimalMathLibrary) with focused C++20 packages** for numerical, symbolic, and scientific computing.

[The packages](#the-packages) | [Quick start](#quick-start) | [Building & linking](#building--linking) | [Documentation](#documentation)

</div>

---

## What is MML Packages?

MML Packages extends MML Core's vectors, matrices, calculus, and solvers with five domain libraries. Use the packages together or link only the ones you need. The separate `mml_ext` header tree contains MML-shaped extensions such as spectral graph analytics and field-line tracing.

The public [releases](https://github.com/zvanjak/MML-Packages/releases) provide prebuilt static libraries, matching C++20 headers, a vendored MML Core snapshot, and relocatable CMake targets for Windows x64, Linux x86-64, and macOS arm64. This is an early 0.1 series: APIs and supported environments may change before 1.0.

## The packages

| Package | What it does | Highlights |
| --- | --- | --- |
| [Optimization](docs/README_Optimization.md) | Search and multi-objective optimization | Genetic algorithms, simulated annealing, NSGA-II, MOEA/D, constraints, solver lifecycle |
| [PDE](docs/README_PDE.md) | Partial differential equation solvers | Structured grids, boundary conditions, Poisson, heat, wave, and advection |
| [Fourier](docs/README_Fourier.md) | Spectral and signal analysis | Fourier series, power spectra, convolution, correlation, Fourier bases |
| [Statistics](docs/README_Statistics.md) | Data and inferential statistics | Descriptors, random generators, time series, hypothesis tests |
| [Symbolic](docs/README_Symbolic.md) | Symbolic computing and differentiation | Expression trees, simplification, automatic differentiation, code generation |
| [mml_ext](include/mml_ext/) | Extensions to MML Core | Spectral graph analysis and field-line tracing (headers, not a sixth package library) |

## Quick start

1. Download your platform archive and `SHA256SUMS` from [Releases](https://github.com/zvanjak/MML-Packages/releases). The [first published prerelease](https://github.com/zvanjak/MML-Packages/releases/tag/v0.1.0-rc.5) is `v0.1.0-rc.5`.
2. Verify the archive's checksum, extract it, and use the extracted directory as your installation prefix. Keep its headers and libraries together.
3. Point CMake at that prefix and link `MML::Packages` for all five libraries, or an individual target such as `MML::PDE`. The [binary installation guide](docs/Binary_Release.md) covers the exact toolchains, ABI requirements, and checksum commands.

For a taste of the API, this solves a one-dimensional Poisson problem on a structured grid:

```cpp
Interval<double> domain(0.0, 1.0);
Grid1D<double> grid(domain, 100);
auto bc = homogeneousDirichlet1D<double>();
PoissonSolver1D<double> poisson(grid, bc);
poisson.setSource(source_function);
auto solution = poisson.solve();
```

The [full demo source](docs_demos/pde/example01_1d_poisson.cpp) defines `source_function`, includes the required headers, and compares the result with the analytical solution. It is compiled into `MML_DocsApp`; enable its call in [the demo main](docs_demos/docs_app_main.cpp) to run this particular example.

## Building & linking

The public repository is a **release companion**, not a standalone source build: it contains headers, documentation, examples, and a vendored MML snapshot, but not the private library build sources or root CMake project. Use a matching [release archive](https://github.com/zvanjak/MML-Packages/releases) for compiled libraries. The installed `MMLPackages` config exports `MML::Packages`, five individual package targets, `MML::Core`, and `MML::Ext`.

The [documentation demo application](docs_demos/) can be built against an extracted release prefix:

```sh
cmake -S docs_demos -B build-docs -DMML_INSTALLED_PREFIX=/path/to/extracted/MML-Packages-prefix
cmake --build build-docs --config Release
```

Pass the actual extracted directory as `MML_INSTALLED_PREFIX`. On multi-configuration generators, run the executable from `build-docs/Release/`; on single-configuration generators, from `build-docs/`. Tested platform toolchains and runtime requirements are listed in the [installation guide](docs/Binary_Release.md), which currently describes the first prerelease.

## Relationship to MML Core

[MinimalMathLibrary](https://github.com/zvanjak/MinimalMathLibrary) supplies the general numerical foundation; these packages add domain-specific algorithms and solvers on top. Each release bundles a matching MML Core header snapshot, so consumers do not need to mix it with a separate MML checkout. Functionality that proves broadly useful may eventually move into MML Core.

## Documentation

- [Package documentation](docs/) covers each domain and its APIs.
- [Runnable documentation demos](docs_demos/) back the examples; start with the [PDE gallery](docs/pde/Examples_Gallery.md) or the [optimization guides](docs/README_Optimization.md).
- [Binary installation](docs/Binary_Release.md) covers toolchains, CMake targets, checksums, and the release manifest.
- [First prerelease notes](docs/Release_Notes_v0.1.0-rc.5.md) record validated assets and known limitations. For later releases, read their accompanying notes and `share/mml-packages/release-manifest.json` rather than assuming the first release's toolchain details apply.

## License

MML Packages and the bundled MML Core are MIT licensed. Binary archives contain both license notices under `licenses/`; keep them with redistributed copies. The [MML Core notice](licenses/MML-Core-LICENSE.md) is also available in this checkout.