# mml_ext — MML Extension Tree

Headers moved out of the MML core with the **Release 2.0 cull** (2026-09-26, epic
MinimalMathLibrary-ya0v), plus incubation code. The tree mirrors the MML directory
layout, and **all namespaces/guards are unchanged** — migrating from pre-2.0 MML is
just an include-path switch (`<mml/...>` → `<mml_ext/...>`).

## Contents

| Header | Provides |
|---|---|
| `algorithms/Graphs/GraphSpectral.h` | Laplacian eigenvalues, algebraic connectivity, degree/closeness/betweenness/eigenvector centrality, PageRank, density, clustering, diameter |
| `algorithms/Analyzers/FieldLineTracer.h` | Vector field line tracing for visualization (incubation) |

**Relocated (2026-09-27):** the statistics headers (`StatisticsHypothesis`, `StatisticsConfidence`,
`StatisticsRank`) and the `RevisedSimplexSolver` LP engine were initially placed here by the cull,
but belong to their domain packages and have moved to `include/statistics/` and
`include/optimization/` respectively — each split into separate `.h`/`.cpp` per the package
convention. Only genuinely extension/incubation algorithms remain in `mml_ext`.

Note: MML core's `SimulatedAnnealing` was **retired** in the same cull — the canonical SA
lives in this repo at `include/optimization/SimulatedAnnealing.h`.

## Usage

Link the `MML_Ext` target (header-only INTERFACE, depends on `MML_Core`):

```cpp
#include <mml_ext/algorithms/Graphs/GraphSpectral.h>

// Spectral graph metrics: algebraic connectivity, centralities, PageRank, ...
```

Tests live in `tests/mml_ext/` (label `mml_ext`).
