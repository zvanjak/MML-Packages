# Symbolic Package

**Symbolic Computing and Automatic Differentiation**

The Symbolic package provides three complementary approaches to computing derivatives: symbolic differentiation using expression trees, forward-mode automatic differentiation with dual numbers, and reverse-mode automatic differentiation with tape-based backpropagation.

To build `MML::Symbolic` as a shared library, configure with both `-DBUILD_SHARED_LIBS=ON` and `-DMML_STATIC=OFF`. The default configuration remains static. In shared builds, the standalone consumer target copies the Symbolic DLL beside its executable for testing.

## Features

### Symbolic Differentiation
- **Expression Trees** - Build mathematical expressions programmatically
- **Parsing** - Build expressions from strings with `parse()`
- **Symbolic Derivatives** - Compute derivative formulas exactly
- **Simplification** - Algebraic simplification of expressions with `simplify()`

### Forward-Mode AD (Dual Numbers)
- **Dual Numbers** - Value and derivative in one object
- **Operator Overloading** - Natural mathematical syntax
- **Efficient for f: ℝ → ℝᵐ** - Few inputs, many outputs

### Reverse-Mode AD (Backpropagation)
- **Tape Recording** - Automatic operation recording
- **Backward Pass** - Efficient gradient computation
- **Efficient for f: ℝⁿ → ℝ** - Many inputs, one output (ML-style)

### Jacobian & Hessian Computation
- **Dense Jacobians** - Full Jacobian matrix computation
- **Hessian** - Second-order derivatives
- **Integration** - Adapters for MML function interfaces

## When to Use Which

| Approach | Best For | Complexity |
|----------|----------|------------|
| **Symbolic** | Formula derivation, education | O(expression size) |
| **Forward AD** | f: ℝ → ℝᵐ (e.g., curve tangent, directional derivative) | O(n × cost(f)) |
| **Reverse AD** | f: ℝⁿ → ℝ (e.g., ML loss, optimization objective) | O(cost(f)) |

## Quick Start

### Symbolic Differentiation

```cpp
#include "mml_packages/symbolic/include/Symbolic.h"

using namespace MML::Symbolic;

// Create symbols
auto x = symbol("x");
auto y = symbol("y");

// Build expression: f = sin(x² + y)
auto f = sin_(x * x + y);

// Symbolic differentiation
auto df_dx = f->diff("x");  // 2x·cos(x² + y)
auto df_dy = f->diff("y");  // cos(x² + y)

// Print formulas
std::cout << "f = " << f->toString() << "\n";
std::cout << "∂f/∂x = " << df_dx->toString() << "\n";

// Evaluate at a point
VarBindings bindings = {{"x", 1.0}, {"y", 0.0}};
double value = f->eval(bindings);
double derivative = df_dx->eval(bindings);

std::cout << "f(1, 0) = " << value << "\n";
std::cout << "∂f/∂x(1, 0) = " << derivative << "\n";
```

### Parsing & Simplification

Build expressions from strings with `parse()`, then differentiate, simplify, and evaluate:

```cpp
#include <Symbolic.h>

using namespace MML::Symbolic;

// Parse an expression from a string (scalar-only engine)
auto f = parse("sin(x^2) + 3*x");

// Symbolic differentiation produces a structurally noisy tree...
auto df = f->diff("x");
// ...which simplify() cleans up to: (cos((x^2)) * (2 * x)) + 3
auto dfSimplified = simplify(df);

std::cout << "f'(x)   = " << dfSimplified->toString() << "\n";
std::cout << "f'(1.3) = " << dfSimplified->eval({{"x", 1.3}}) << "\n";

// Constant folding + identity elimination:  (x*1 + 0) + (2+3)  ->  x + 5
auto g = simplify(parse("x*1 + 0 + (2 + 3)"));
```

**Supported syntax:** `+ - * /`, `^` (right-associative), unary `+`/`-`, parentheses, the
constants `pi` and `e`, and function calls `sin cos tan asin acos atan sinh cosh tanh exp`
`log`/`ln` `log10 sqrt abs` (one argument) and `pow(base, exponent)`. Any other identifier is a
symbol. Vector/matrix literals and index access are intentionally **not** supported - this is a
scalar engine. Use `tryParse()` for a non-throwing variant that returns `nullptr` on failure.

### Forward-Mode AD (Dual Numbers)

```cpp
#include "mml_packages/symbolic/include/ForwardAD.h"

using namespace MML::Symbolic;

// Dual<T>: value and derivative together
// Seed derivative = 1.0 for the variable you're differentiating

// Compute f(x) = x³ - 2x + 1 and f'(x) at x = 2
Dual<double> x(2.0, 1.0);  // x = 2, dx/dx = 1
Dual<double> f = x*x*x - 2.0*x + 1.0;

std::cout << "f(2) = " << f.value() << "\n";      // 5.0
std::cout << "f'(2) = " << f.deriv() << "\n";     // 10.0 = 3x² - 2

// Multi-variable: compute partial derivative
// ∂/∂x[sin(xy)] at (π/4, 2)
Dual<double> x2(M_PI/4, 1.0);  // Seed x
double y_val = 2.0;             // y is constant
Dual<double> result = sin(x2 * y_val);

std::cout << "∂sin(xy)/∂x = " << result.deriv() << "\n";  // y·cos(xy)
```

### Reverse-Mode AD (Gradient Computation)

```cpp
#include "mml_packages/symbolic/include/ReverseAD.h"

using namespace MML::Symbolic;

// Tape records all operations
Tape tape;

// Create input variables
ADVar x(&tape, 2.0);
ADVar y(&tape, 3.0);
ADVar z(&tape, 1.5);

// Forward pass: compute f = x²y + exp(z) + sin(xyz)
ADVar f = x*x*y + exp(z) + sin(x*y*z);

// Backward pass: compute all gradients in one sweep
f.backward();

// Read gradients
std::cout << "f = " << f.value() << "\n";
std::cout << "∂f/∂x = " << x.adjoint() << "\n";
std::cout << "∂f/∂y = " << y.adjoint() << "\n";
std::cout << "∂f/∂z = " << z.adjoint() << "\n";
```

### Jacobian Matrix Computation

```cpp
#include "mml_packages/symbolic/include/Jacobian.h"

using namespace MML::Symbolic;

// Vector function f: ℝ² → ℝ³
// f(x,y) = [x² + y, sin(xy), exp(x) - y²]

auto f = [](const std::vector<Dual<double>>& vars) {
    Dual<double> x = vars[0];
    Dual<double> y = vars[1];
    return std::vector<Dual<double>>{
        x*x + y,
        sin(x * y),
        exp(x) - y*y
    };
};

std::vector<double> point = {1.0, 2.0};

// Compute full Jacobian matrix
DenseMatrix J = jacobian(f, point, 2, 3);

// J is 3×2 matrix:
// | ∂f₁/∂x  ∂f₁/∂y |   | 2x   1     |   | 2    1  |
// | ∂f₂/∂x  ∂f₂/∂y | = | y·cos(xy)  x·cos(xy) | = | 2cos2  cos2 |
// | ∂f₃/∂x  ∂f₃/∂y |   | exp(x)   -2y  |   | e    -4 |

for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j)
        std::cout << J(i, j) << " ";
    std::cout << "\n";
}
```

### MML Integration

```cpp
#include "mml_packages/symbolic/include/MMLAdapters.h"

using namespace MML::Symbolic;

// Create symbolic function
auto x = symbol("x");
auto expr = sin_(x * x);

// Wrap as MML IRealFunction
SymbolicRealFunction f(expr, "x");

// Use with MML algorithms
double val = f(1.5);            // Evaluate
double deriv = f.derivative(1.5);  // Automatic derivative

// Use with root finding, integration, etc.
Real root = RootFinding::Bisection(f, 0.1, 2.0);
```

## API Reference

### Expression Building

| Factory | Creates | Example |
|---------|---------|---------|
| `constant(v)` | Numeric constant | `constant(3.14)` |
| `symbol(name)` | Variable symbol | `symbol("x")` |
| `x()`, `y()`, `z()`, `t()` | Common symbols | `x() + y()` |

### Expression Operations

| Operation | Syntax | Result |
|-----------|--------|--------|
| Addition | `a + b` | Sum expression |
| Subtraction | `a - b` | Difference expression |
| Multiplication | `a * b` | Product expression |
| Division | `a / b` | Quotient expression |
| Power | `pow_(a, b)` | Power expression |
| Negation | `-a` | Negated expression |

### Mathematical Functions

| Function | Call | Derivative |
|----------|------|------------|
| sin | `sin_(x)` | cos(x) |
| cos | `cos_(x)` | -sin(x) |
| tan | `tan_(x)` | sec²(x) |
| exp | `exp_(x)` | exp(x) |
| log | `log_(x)` | 1/x |
| sqrt | `sqrt_(x)` | 1/(2√x) |
| abs | `abs_(x)` | sign(x) |
| asin | `asin_(x)` | 1/√(1-x²) |
| acos | `acos_(x)` | -1/√(1-x²) |
| atan | `atan_(x)` | 1/(1+x²) |
| sinh | `sinh_(x)` | cosh(x) |
| cosh | `cosh_(x)` | sinh(x) |
| tanh | `tanh_(x)` | sech²(x) |

### Dual Number Operations

```cpp
Dual<T> x(value, derivative);  // Constructor
x.value();                      // Get value
x.deriv();                      // Get derivative

// All arithmetic and math functions are overloaded
Dual<T> y = sin(x) + exp(x);
```

### Tape Operations

```cpp
Tape tape;
tape.clear();                   // Reset tape
tape.size();                    // Number of operations

ADVar x(&tape, value);          // Create tracked variable
result.backward();              // Compute all gradients
x.adjoint();                    // Get ∂output/∂x
```

## File Structure

```
symbolic/
├── README.md              # This file
├── CMakeLists.txt         # Build configuration
│
├── include/
│   ├── Symbolic.h         # Main include header
│   │
│   ├── Expr.h             # Expression tree base class
│   ├── Operations.h       # Arithmetic operations
│   ├── Functions.h        # Mathematical functions
│   │
│   ├── ForwardAD.h        # Dual numbers (forward mode)
│   ├── ReverseAD.h        # Tape-based AD (reverse mode)
│   │
│   ├── Jacobian.h         # Jacobian/Hessian computation
│   ├── MMLAdapters.h      # Integration with MML interfaces
│   │
│   └── ...
│
└── tests/
    └── symbolic_tests.cpp
```

## Mathematical Background

### Forward-Mode AD (Dual Numbers)

Dual numbers: $a + b\epsilon$ where $\epsilon^2 = 0$

For function $f(x)$:
$$f(a + \epsilon) = f(a) + f'(a)\epsilon$$

This propagates derivatives alongside values.

### Reverse-Mode AD

Record computation as directed graph, then backpropagate:
$$\frac{\partial L}{\partial x_i} = \sum_{j \in \text{children}(i)} \frac{\partial L}{\partial x_j} \cdot \frac{\partial x_j}{\partial x_i}$$

One backward pass computes all partial derivatives.

### Symbolic Differentiation

Apply differentiation rules to expression tree:
- $\frac{d}{dx}[f + g] = \frac{df}{dx} + \frac{dg}{dx}$
- $\frac{d}{dx}[f \cdot g] = f \cdot \frac{dg}{dx} + g \cdot \frac{df}{dx}$
- $\frac{d}{dx}[f(g)] = f'(g) \cdot \frac{dg}{dx}$ (chain rule)

## Performance Considerations

| Method | Memory | Time | Best Use Case |
|--------|--------|------|---------------|
| Forward AD | O(n) per eval | O(n × cost) | Gradients for n ≤ 10 |
| Reverse AD | O(ops) tape | O(cost) | Gradients for n > 10 |
| Symbolic | Expression tree | Depends on expr | Formula derivation |

**Tips:**
- Use reverse AD for optimization/ML (many inputs → one output)
- Use forward AD for Jacobians of small systems
- Use symbolic differentiation when the derivative formula itself is needed

## References

- Griewank, A., and Walther, A. (2008). "Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation"
- Baydin, A.G., et al. (2018). "Automatic Differentiation in Machine Learning: a Survey"
- Aho, Sethi, Ullman (2006). "Compilers: Principles, Techniques, and Tools" (expression trees)

## See Also

- [mml_packages/README.md](../README.md) - Package overview
- [Optimization Package](../optimization/README.md) - Uses AD for gradient-based optimization
- [Systems Package](../systems/README.md) - Jacobians for dynamical systems

## Design Notes (for contributors)

The symbolic engine follows two deliberate design decisions (see
`analysis/05_Symbolic_support_overview.md` in the MinimalMathLibrary repo). Please preserve them:

1. **`shared_ptr<Expr>` + virtual-method dispatch.** Operations intrinsic to a node (`eval`,
   `diff`, `clone`, `toString`) are virtual methods on the node. New whole-tree operations
   (`simplify`, and the output side of `parse`) are implemented as **free-function passes over the
   tree**, not as a visitor rewrite. This keeps adding *node types* cheap, keeps call sites simple
   (`f->diff("x")`), and lets sub-trees be shared cheaply during differentiation. Do not convert
   the engine to a `unique_ptr` + visitor model.

2. **Scalar (`double`) only.** Vectors, matrices, and polynomials are first-class *numeric* types
   elsewhere in MML, but they are intentionally **not** modelled as symbolic operands in the
   `Expr` tree. Downstream consumers (e.g. SigmaEngine) keep any typed handling in a thin shim
   over this scalar engine. Do not push a typed/variant value system into `Expr` without first
   revisiting that analysis.
