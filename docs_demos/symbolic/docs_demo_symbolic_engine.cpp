///////////////////////////////////////////////////////////////////////////////////////////
// docs_demo_symbolic_engine.cpp - Symbolic engine: parse -> diff -> simplify -> eval
//
// Demonstrates the MML::Symbolic engine end-to-end:
// 1. Parsing an expression from a string
// 2. Evaluating it at a point
// 3. Symbolic differentiation
// 4. Algebraic simplification of the derivative
// 5. Constant folding + identity elimination
///////////////////////////////////////////////////////////////////////////////////////////

#include <Symbolic.h>

#include <iostream>

using namespace MML::Symbolic;

void Docs_Demo_Symbolic_Engine()
{
    std::cout << "\n=== Symbolic engine: parse -> diff -> simplify -> eval ===\n";

    // 1. Parse an expression from a string (scalar-only engine)
    auto f = parse("sin(x^2) + 3*x");
    std::cout << "f(x)             = " << f->toString() << "\n";

    // 2. Evaluate at a point
    std::cout << "f(1.3)           = " << f->eval({{"x", 1.3}}) << "\n";

    // 3. Symbolic differentiation
    auto df = f->diff("x");
    std::cout << "f'(x) raw        = " << df->toString() << "\n";

    // 4. Simplify the (structurally noisy) derivative tree
    auto dfSimplified = simplify(df);
    std::cout << "f'(x) simplified = " << dfSimplified->toString() << "\n";
    std::cout << "f'(1.3)          = " << dfSimplified->eval({{"x", 1.3}}) << "\n";

    // 5. Constant folding + identity elimination in one pass
    auto g = parse("x*1 + 0 + (2 + 3)");
    std::cout << "g(x)             = " << g->toString()
              << "  ->  " << simplify(g)->toString() << "\n";
}
