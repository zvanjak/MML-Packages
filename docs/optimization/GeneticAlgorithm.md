# Genetic Algorithm (GA)

## Overview

The Genetic Algorithm implementation provides a **real-coded, mixed-integer** optimization framework for single-objective problems. Variables can be continuous, integer, or categorical; binary variables require operators that are not provided by this implementation.

**Location:** `include/optimization/algorithms/GeneticAlgorithm.h`

## Features

- **Mixed Variable Types**: Continuous, Integer, and Categorical
- **Modular Operators**: Pluggable selection, crossover, and mutation strategies
- **Elitism**: Preserve best solutions across generations
- **Stagnation Detection**: Early termination when no improvement
- **Type-Aware Operations**: Different operators for different variable types

## Quick Start

### Basic Continuous Optimization

```cpp
#include "GeneticAlgorithm.h"
#include "ProblemTypes.h"

using namespace MML;
using namespace MML::Optimization;

// Define objective function
auto quadratic = [](const Vector<Real>& x) -> Real {
    Real dx = x[0] - 2.0;
    Real dy = x[1] - 3.0;
    return dx * dx + dy * dy;
};

ProblemSpec spec = ProblemSpec::Continuous(2, -10.0, 10.0);

// Run GA
GeneticAlgorithm ga;
ga.SetProblem(spec);
GAResult result = ga.Minimize(quadratic);

std::cout << "Best fitness: " << result.bestFitness << "\n";
std::cout << "Generations: " << result.generations << "\n";
```

For a typed single-objective problem, GA reads the variable spec and scalar evaluation directly:

```cpp
SingleObjectiveProblem typedProblem(quadratic, spec);
GeneticAlgorithm typedGA;
GAResult typedResult = typedGA.Minimize(typedProblem);
```

`ProblemSpec` describes decision variables and bounds, not GA chromosomes or operators.
Both typed and callable GA entry points use the same algorithm-owned real-coded
selection, crossover, mutation, and repair. The typed overloads take an
`ISingleObjectiveProblem` directly; callable overloads require `SetProblem(spec)`.
Continuous, integer, and categorical coordinates are supported, with discrete
coordinates repaired to valid values after variation. Binary variables and
constrained typed problems are rejected before evaluation. `Minimize` and `Maximize`
select the direction at the algorithm call; `Maximize` returns the original
objective value, not the internal fitness used for selection.

### Mixed-Integer Optimization

```cpp
ProblemSpec spec;
spec.AddVariable(OptVariableSpec::Continuous(0.0, 5.0, "x"));
spec.AddVariable(OptVariableSpec::Integer(1, 5, "n"));
spec.AddVariable(OptVariableSpec::Categorical(4, "category"));

auto evaluate = [](const Vector<Real>& x) -> Real {
    Real continuous = (x[0] - 2.5) * (x[0] - 2.5);
    Real integer = (x[1] - 3.0) * (x[1] - 3.0);
    int category = static_cast<int>(std::round(x[2]));
    return continuous + integer + (category == 2 ? 0.0 : 1.0);
};

GAConfig config;
config.populationSize = 50;
config.maxGenerations = 100;

GeneticAlgorithm ga(config);
ga.SetProblem(spec);
GAResult result = ga.Minimize(evaluate);
```

## Variable Types

### Continuous

Real-valued variables with any value in the specified range.

```cpp
// Variable x in [-10, 10]
spec.AddVariable(OptVariableSpec::Continuous(-10.0, 10.0, "x"));

// Or using factory
auto varSpec = OptVariableSpec::Continuous(-10.0, 10.0, "x");
```

### Integer

Integer-valued variables. Values are rounded during crossover and mutation.

```cpp
// Integer n in [1, 100]
spec.AddVariable(OptVariableSpec::Integer(1, 100, "n"));
```

### Categorical

Discrete set of options, encoded as integers 0, 1, 2, ...

```cpp
// 5 categories (values 0-4)
spec.AddVariable(OptVariableSpec::Categorical(5, "category"));
```

## Configuration

### GAConfig Structure

| Parameter | Default | Description |
|-----------|---------|-------------|
| `populationSize` | 100 | Number of individuals |
| `numElites` | 2 | Best individuals preserved unchanged |
| `crossoverRate` | 0.9 | Probability of crossover per pair |
| `mutationRate` | 0.1 | Per-gene mutation probability |
| `maxGenerations` | 1000 | Maximum generations |
| `maxFuncEvals` | 100000 | Max function evaluations (0=unlimited) |
| `stagnationLimit` | 100 | Generations without improvement before stop |
| `targetFitness` | -∞ | Stop if this fitness is reached |
| `fitnessThreshold` | 1e-10 | Convergence threshold |
| `minimize` | true | Minimize (true) or maximize (false) |
| `seed` | 0 | Random seed (0=time-based) |
| `verbose` | false | Print progress to stdout |
| `trackHistory` | false | Track fitness history per generation |

### Example Configuration

```cpp
GAConfig config;
config.populationSize = 200;
config.maxGenerations = 500;
config.crossoverRate = 0.8;
config.mutationRate = 0.05;
config.numElites = 5;
config.stagnationLimit = 50;
config.seed = 42;  // For reproducibility
config.trackHistory = true;

GeneticAlgorithm ga(config);
```

## Selection Strategies

### Tournament Selection (Default)

Picks k random individuals and selects the best.

```cpp
// Binary tournament (k=2) - moderate selection pressure
ga.SetSelection(std::make_unique<TournamentSelection>(2));

// Larger tournament (k=5) - higher selection pressure
ga.SetSelection(std::make_unique<TournamentSelection>(5));
```

### Rank Selection

Selection probability based on rank, not raw fitness. More stable when fitness values vary widely.

```cpp
// Selective pressure 1.5 (moderate)
ga.SetSelection(std::make_unique<RankSelection>(1.5));

// Selective pressure 2.0 (high)
ga.SetSelection(std::make_unique<RankSelection>(2.0));
```

## Crossover Operators

### BLX-α Crossover (Default)

Blend crossover samples offspring from an extended range around parents.

```cpp
// BLX-0.5 (most common)
ga.SetCrossover(std::make_unique<BLXCrossover>(0.5));

// BLX-0.0 (children strictly between parents)
ga.SetCrossover(std::make_unique<BLXCrossover>(0.0));
```

**How it works:**
- For gene i: `range = |parent1[i] - parent2[i]|`
- Child sampled from `[min - α*range, max + α*range]`

### SBX Crossover (Simulated Binary)

Mimics single-point binary crossover behavior. Widely used in NSGA-II.

```cpp
// η=15 (moderate spread)
ga.SetCrossover(std::make_unique<SBXCrossover>(15.0));

// η=2 (high spread, more exploration)
ga.SetCrossover(std::make_unique<SBXCrossover>(2.0));

// η=20 (low spread, children close to parents)
ga.SetCrossover(std::make_unique<SBXCrossover>(20.0));
```

### Uniform Crossover

Each gene independently comes from either parent.

```cpp
// 50% swap probability (default)
ga.SetCrossover(std::make_unique<UniformCrossover>(0.5));
```

## Mutation Operators

### Gaussian Mutation (Default)

Adds Gaussian noise proportional to variable range.

```cpp
// σ = 10% of variable range (default)
ga.SetMutation(std::make_unique<GaussianMutation>(0.1));

// σ = 5% of variable range (smaller perturbations)
ga.SetMutation(std::make_unique<GaussianMutation>(0.05));
```

**Type-aware behavior:**
- Continuous: Gaussian perturbation
- Integer/Categorical: Step mutation (±1) or resample

### Polynomial Mutation

Bounded polynomial distribution (NSGA-II style).

```cpp
// η=20 (small perturbations, common)
ga.SetMutation(std::make_unique<PolynomialMutation>(20.0));

// η=100 (very small perturbations)
ga.SetMutation(std::make_unique<PolynomialMutation>(100.0));
```

## Result Structure

### GAResult Fields

| Field | Type | Description |
|-------|------|-------------|
| `bestSolution` | `Vector<Real>` | Best solution found |
| `bestFitness` | `Real` | Best fitness value |
| `generations` | `int` | Generations completed |
| `funcEvals` | `int` | Total function evaluations |
| `stagnationCount` | `int` | Consecutive generations without improvement |
| `converged` | `bool` | True if converged/stagnated |
| `terminationReason` | `string` | Human-readable reason |
| `fitnessHistory` | `vector<Real>` | Best fitness per generation (if tracked) |
| `avgFitnessHistory` | `vector<Real>` | Average fitness per generation (if tracked) |

### Example Result Processing

```cpp
GAResult result = ga.Minimize(func);

std::cout << "=== GA Results ===\n";
std::cout << "Best fitness: " << result.bestFitness << "\n";
std::cout << "Generations: " << result.generations << "\n";
std::cout << "Function evals: " << result.funcEvals << "\n";
std::cout << "Termination: " << result.terminationReason << "\n";

std::cout << "Best solution: [";
for (int i = 0; i < result.bestSolution.size(); ++i) {
    if (i > 0) std::cout << ", ";
    std::cout << result.bestSolution[i];
}
std::cout << "]\n";

// If history was tracked
if (!result.fitnessHistory.empty()) {
    std::cout << "Fitness progression:\n";
    for (size_t i = 0; i < result.fitnessHistory.size(); i += 10) {
        std::cout << "  Gen " << i << ": " << result.fitnessHistory[i] << "\n";
    }
}
```

## Convenience Functions

### Quick Minimization

```cpp
// With problem spec
auto result = GAMinimize(func, spec, 100, 500);  // popSize=100, maxGen=500

// With bounds (all continuous)
Vector<Real> lb(5), ub(5);
// ... set bounds ...
auto result = GAMinimize(func, lb, ub, 100, 500);
```

### Problem Spec Factories

```cpp
// Uniform bounds
auto spec = ProblemSpec::Continuous(2, -10.0, 10.0);

// Different bounds per variable
Vector<Real> lb{-10.0, -5.0};
Vector<Real> ub{10.0, 5.0};
auto spec = ProblemSpec::Continuous(lb, ub);
```

## Benchmark Functions

The test suite includes several standard benchmark functions:

| Function | Dimension | Global Minimum | Characteristics |
|----------|-----------|----------------|-----------------|
| Sphere | n | 0 at origin | Unimodal, easy |
| Rastrigin | n | 0 at origin | Highly multimodal |
| Rosenbrock | 2 | 0 at (1,1) | Valley-shaped |
| Ackley | n | 0 at origin | Many local minima |

## Best Practices

### Population Size
- Small problems (n < 10): 50-100 individuals
- Medium problems (n < 50): 100-200 individuals
- Large problems (n > 50): 200-500 individuals

### Selection Pressure
- Start with tournament size k=3 (balanced)
- Increase to k=5 if premature convergence
- Use rank selection if fitness scaling is extreme

### Crossover Rate
- Keep high (0.8-0.95) for most problems
- Lower if population diversity is lost too quickly

### Mutation Rate
- Per-gene rate: 1/n to 0.1 (where n = dimension)
- Higher for multimodal problems
- Lower for fine-tuning near optimum

### Elitism
- Always use some elitism (2-5 individuals)
- Prevents loss of best solutions
- Ensures monotonic improvement in best fitness

## Comparison with Simulated Annealing

| Aspect | Genetic Algorithm | Simulated Annealing |
|--------|-------------------|---------------------|
| Population | Yes (parallel search) | No (single solution) |
| Memory | Higher (stores population) | Lower |
| Exploration | Crossover + mutation | Random neighbors |
| Mixed-integer | Integer and categorical supported | Not supported, even with a custom neighbor |
| Parallelizable | Yes | Limited |
| Best for | Multimodal, mixed-type | Smooth, continuous |

## See Also

- [Simulated annealing demo](../../docs_demos/optimization/docs_demo_simulated_annealing.cpp) - Continuous-variable alternative
- [Optimization package README](../../src/optimization/README.md) - Shared problem and algorithm policies
