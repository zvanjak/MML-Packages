# Optimization Package

## MML core plus package-specific optimization algorithms

Link to `MML::Optimization` to use both the package and its public `MML_Core` dependency.
Core algorithms live in the MML optimization headers; this package adds evolutionary
methods, dynamic problem types, penalties, observers, and termination policies. MML's
fixed-dimension solvers use `VectorN<Real, N>` and `IScalarFunction<N>` (or
`IDifferentiableScalarFunction<N>`). Package algorithms generally use runtime-sized
`Vector<Real>` and their own problem interfaces. Both sets of types use the
`MML::Optimization` namespace, except the 1D routines in `MML::Minimization`.

## MML Core Facilities

These headers come from MML (vendored here under `libs/mml/`), not from this package:

| Area | Facilities | Public header |
| --- | --- | --- |
| One-dimensional minimization | Minimum bracketing, golden-section search, Brent minimization (with or without derivatives) | `mml/algorithms/Optimization/Optimization.h` |
| Derivative-free multidimensional | Nelder-Mead simplex, Powell direction set, line minimization | `mml/algorithms/Optimization/OptimizationMultidim.h` |
| Gradient-based multidimensional | Conjugate gradient (Fletcher-Reeves / Polak-Ribiere), BFGS, limited-memory BFGS | `mml/algorithms/Optimization/OptimizationMultidim.h` |
| Box constraints | Bounds, projected gradient, box-constrained Nelder-Mead and Powell | `mml/algorithms/Optimization/OptimizationMultidim.h` |
| General nonlinear constraints | Augmented Lagrangian for equality and inequality constraints | `mml/algorithms/Optimization/Constraints/AugmentedLagrangian.h` |
| Dense continuous linear programming | LP model, simplex tableau and solver, LP result types and solve helpers | `mml/algorithms/Optimization/LinearProgramming.h` |

`OptimizationMultidim.h` aggregates the multidimensional solvers and box-constraint
headers; Augmented Lagrangian is included separately. The MML 1D and LP headers are
also separate. The package does not ship copies of these MML headers or solvers.

## Package Facilities

| Area | Facilities | Public header |
| --- | --- | --- |
| Evolutionary search | Genetic algorithm with continuous and discrete variable operators | `optimization/algorithms/GeneticAlgorithm.h` |
| Stochastic search | Simulated annealing | `optimization/algorithms/SimulatedAnnealing.h` |
| Multi-objective search | NSGA-II and MOEA/D, Pareto sorting and variation | `optimization/algorithms/NSGA2.h`, `optimization/algorithms/MOEAD.h`, `optimization/core/MultiObjective.h` |
| Problem model | Typed single- and multi-objective problems, variable specifications, constraints and penalty method | `optimization/core/OptimizationProblem.h` |
| Run control | Shared iterative lifecycle, optimization configuration, observers and termination criteria | `optimization/interfaces/IIterativeAlgorithm.h`, `optimization/core/IterativeRunner.h`, `optimization/core/OptimizationConfig.h`, `optimization/core/OptimizationObservers.h`, `optimization/core/TerminationCriteria.h` |
| LP extension | Revised simplex implementation alongside MML's LP model | `optimization/algorithms/LP/RevisedSimplexSolver.h` |
| Benchmarks | Single-objective, constrained and multi-objective test problems | `test_problems/` |

The `AsDynamic` adapter in `optimization/interfaces/FunctionAdapters.h` bridges an MML
`IScalarFunction<N>` into a runtime-sized scalar objective; it checks the input
dimension. A package `PenaltyMethod` takes a typed
`ConstrainedSingleObjectiveProblem`, whereas MML's `AugmentedLagrangian<N>` works
with fixed-dimension objectives and constraints.

### Shared problems, algorithm-specific policies

`IOptimizationProblem` describes dimension, `ProblemSpec` variables, and constraints;
it does not evaluate an objective. `ISingleObjectiveProblem::Evaluate` returns one
scalar, while `IMultiObjectiveProblem::EvaluateObjectives` returns a vector of the
declared objective count. Decision vectors use `Vector<Real>` for every variable type:
`ProblemSpec` defines the bounds and whether each coordinate is continuous, integer,
binary, or categorical. Problem metadata does not select a chromosome encoding,
repair strategy, neighbor generator, or constraint-handling rule.

| Algorithm | Typed problem | Variables | Constraints | Other entry points |
| --- | --- | --- | --- | --- |
| GA | `ISingleObjectiveProblem` | Continuous, integer, categorical; algorithm-owned chromosome and operators | Rejected | Callable `Minimize`/`Maximize` with `SetProblem` |
| SA | `ISingleObjectiveProblem` | Continuous; bounds from `ProblemSpec`, algorithm-owned cooling and neighbors | Rejected | No callable or explicit-bounds `Minimize` |
| NSGA-II | `IMultiObjectiveProblem` | Continuous, integer, categorical; repair after real-coded variation | Rejected | Callable `Optimize` |
| MOEA/D | `IMultiObjectiveProblem` | Continuous, integer, categorical; repair after real-coded variation | Rejected | Callable `Optimize` |

Binary variables are not supported by these algorithm entry points. Typed acceptance
does not imply support for every variable or constraint type. GA chooses minimization
or maximization at the call site; SA minimizes its scalar objective; NSGA-II and
MOEA/D minimize all objectives. Use the typed constrained problem and `PenaltyMethod`
for supported scalar penalty workflows, not as implicit constraint support in the
algorithms above.

The typed interfaces changed source signatures in the 0.1 package; there are no
compatibility aliases for the former GA-specific variable model or SA callable
overloads. C++ binary compatibility is not guaranteed across these changes: rebuild
clients and the package together with a compatible compiler and configuration.

### Real-only run control and GA types

Optimization state, results, diagnostics, criteria, observers, config presets,
GA individuals, populations and selection strategies now use MML's `Real`
without a numeric template argument. For example, replace
`OptimizationState<Real>` with `OptimizationState`, `OptimizationResult<Real>`
with `OptimizationResult`, `OptimizationConfig<Real>` with `OptimizationConfig`,
`MaxIterationsCriterion<Real>` with `MaxIterationsCriterion`, and
`CreateDefaultCriterion<Real>()` with `CreateDefaultCriterion()`. Likewise,
derive custom criteria from `ITerminationCriterion`, observers from
`IOptimizationObserver`, and GA selectors from `ISelectionStrategy`; their
state and population arguments are `const OptimizationState&` and
`const Population&`. Use `Individual` and `Population` without `<Real>`.
Config presets such as `QuickConfig()` no longer take `<Real>` either.

`Real` defaults to `double` in MML and is chosen consistently for the entire
build, not separately per optimizer. This API does not support mixing
`Vector<float>` and `Vector<double>` in one default-`Real` optimization run;
convert at the problem boundary when necessary. Numeric template removal
changes source and C++ ABI compatibility: rebuild clients, including custom
criterion, observer and selection implementations, and the package together.
Semantic templates such as `IterativeRunner<ScalarIterationProgress>`,
`AsDynamic<N>` and callable entry points remain templated.

### Iterative runs

GA, SA, NSGA-II and MOEA/D implement `IIterativeAlgorithm`. Bind inputs with the
algorithm's `Start` overload, then use `Progress`, `Step` and `Run` through the shared
interface. `Start` evaluates the initial point/population; each `Step` completes one
SA move, GA/NSGA-II generation, or MOEA/D subproblem sweep. `Run` continues from the
same state until it stops. Existing `Minimize`/`Maximize`/`Optimize` entry points bind
and run the same lifecycle to completion.

The following excerpt is exercised by
[docs_demo_iterative_lifecycle.cpp](../../docs_demos/optimization/docs_demo_iterative_lifecycle.cpp):

```cpp
nsga.SetProblem(spec);
moead.SetProblem(spec);
OptimizationConfig runConfig;
runConfig.WithMaxIterations(2).WithTrajectory();
ga.Start(scalar, true, &runConfig);
sa.Start(scalar, Vector<Real>{0.5, 0.5});
nsga.Start(biObjective);
moead.Start(biObjective);
std::array<IIterativeAlgorithm*, 4> algorithms{&ga, &sa, &nsga, &moead};
for (auto* algorithm : algorithms) {
    algorithm->Step();
    algorithm->Run();
    std::cout << "Completed " << algorithm->Progress().iterations
              << " steps and " << algorithm->Progress().funcEvals << " evaluations\n";
}
```

Progress starts at zero completed steps and counts initial evaluations. It includes
steady-clock elapsed time, terminal status and reason. GA/SA return
`ScalarIterationProgress` with current/best values; NSGA-II/MOEA/D return
`ParetoIterationProgress` with objective count and a non-owning population view.
NSGA-II reports its ranked front size; MOEA/D exposes its ideal point. Call
`ObserveFront()` before `Start` to request a per-step front snapshot on either
multi-objective algorithm; MOEA/D only computes its front size when enabled.
Neither Pareto algorithm invents a scalar best value. Progress views remain valid
until the next run/step; do not retain spans or pointers across those transitions.

An algorithm instance owns at most one active run. Rebinding configuration or a
problem, or reinitializing, during an active run is an error. Once stopped or
aborted, `Start` can begin a fresh run. A bound problem, callable and any borrowed
`OptimizationConfig` must outlive that run. Iteration, evaluation and time
limits are checked between complete steps; an evaluation cap too small for the
initial population is rejected, and a sweep/generation is skipped if it cannot
fit the remaining budget. Scalar SA/GA use `OptimizationConfig` criteria and
observers; NSGA-II/MOEA/D use typed Pareto `SetCriterion` and `AddObserver` hooks.
Observers see the initial state, then each completed step, then one completion;
an observer veto yields `Aborted`, while an exhausted limit yields `LimitReached`.

`SimulatedAnnealing::Minimize` is now non-`const` because the run state is owned by
the algorithm. Rebuild C++ clients against these headers: adding the interface and
changing SA's method qualification are source/API and ABI changes.

For executable package examples, see the
[genetic algorithm](../../docs_demos/optimization/docs_demo_genetic_algorithm.cpp),
[iterative lifecycle](../../docs_demos/optimization/docs_demo_iterative_lifecycle.cpp),
[simulated annealing](../../docs_demos/optimization/docs_demo_simulated_annealing.cpp),
[NSGA-II](../../docs_demos/optimization/docs_demo_nsga2.cpp), and
[linear programming](../../docs_demos/optimization/docs_demo_linear_programming.cpp)
demos. They are built into `MML_DocsApp`.

## File Structure

```text
MML-Packages-Private/
├── include/optimization/          # Package public headers
│   ├── algorithms/               # GA, SA, NSGA-II, MOEA/D, package LP solver
│   ├── core/                     # Problem types, constraints, penalty method
│   ├── interfaces/               # Typed optimization interfaces
│   └── test_problems/            # Benchmark problems
├── src/optimization/             # Package implementations, CMakeLists.txt, this README
├── tests/optimization/           # Package tests
└── docs_demos/optimization/      # Runnable package examples

MML core (vendored as libs/mml/):
└── algorithms/Optimization/
    ├── Optimization.h
    ├── OptimizationMultidim.h   # Nelder-Mead, Powell, BFGS, CG, L-BFGS
    ├── LinearProgramming.h      # MML linear programming
    ├── Multidim/
    ├── Constraints/
    └── LP/
```

## Algorithm Selection Guide

| Problem type | Start with | Provider |
| --- | --- | --- |
| Smooth, differentiable | BFGS or conjugate gradient | MML |
| Derivative-free local search | Nelder-Mead or Powell | MML |
| Large-scale differentiable | L-BFGS | MML |
| Box constraints | Projected gradient (with derivatives), box-constrained Nelder-Mead or Powell | MML |
| Nonlinear equality/inequality constraints | Augmented Lagrangian or the package penalty method | MML / package |
| Dense continuous linear program | LP simplex | MML |
| Discrete or mixed-integer variables | Genetic algorithm | Package |
| Stochastic search | Simulated annealing or genetic algorithm | Package |
| Multiple objectives | NSGA-II or MOEA/D | Package |
| One-dimensional scalar function | Brent or golden-section search | MML |

## Performance Tips

1. **Good initial guess** - Local solvers often converge faster from better starts
2. **Scaling** - Normalize variables to similar ranges for better conditioning
3. **Gradient accuracy** - Analytical gradients are faster and more accurate than numerical
4. **Population sizing** - GA/NSGA-II: roughly 10× number of variables
5. **Termination** - Set appropriate tolerance based on problem requirements

## References

- Nelder, J.A., and Mead, R. (1965). "A Simplex Method for Function Minimization"
- Powell, M.J.D. (1964). "An Efficient Method for Finding the Minimum of a Function of Several Variables Without Calculating Derivatives"
- Nocedal, J., and Wright, S.J. (2006). "Numerical Optimization" (2nd ed.)
- Deb, K., et al. (2002). "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II"
- Numerical Recipes in C++, 3rd Edition, Chapter 10: Minimization or Maximization of Functions

## See Also

- [Package overview](../README.md)
