///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Simplifier.h                                                        ///
///  Description: Algebraic simplification of symbolic expression trees               ///
///               Implemented as a free-function pass over the Expr tree (per the     ///
///               shared_ptr + virtual-method design; no visitor rewrite).            ///
///                                                                                   ///
///  Rules: constant folding, identity elimination (x+0, x*1, x^1), zero absorption   ///
///         (x*0, 0/x), power rules (x^0, x^1, (x^a)^b), negation cleanup, a few      ///
///         structural rules (x+x -> 2*x, x*x -> x^2, x/x -> 1, x-x -> 0), and        ///
///         function identities (sqrt(0), exp(log(x)) -> x, log(exp(x)) -> x).        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_SIMPLIFIER_H
#define MML_SYMBOLIC_SIMPLIFIER_H

#include <mml/mml_export.h>

#include "Expr.h"

namespace MML::Symbolic
{
    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Options controlling algebraic simplification.
    //////////////////////////////////////////////////////////////////////////////////////////
    struct SimplifyOptions
    {
        bool foldConstants    = true;   ///< Evaluate variable-free subexpressions to a constant
        bool eliminateIdentity = true;  ///< Remove +0, *1, x^1, etc.
        bool eliminateZero    = true;   ///< Simplify *0, 0/x, etc.
        bool simplifyPowers   = true;   ///< Simplify x^0, x^1, (x^a)^b, etc.
        bool simplifyNegation = true;   ///< Simplify --x, (-a)*(-b), etc.
        bool collectPolynomials = true; ///< Collect univariate polynomial like-terms (3x^2+4x-6)
        int  maxIterations    = 16;     ///< Max fixpoint passes
        double epsilon        = 1e-15;  ///< Tolerance for comparison to 0/1
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Simplify a symbolic expression, returning a new (simplified) tree.
    ///
    /// The input tree is not modified; a fresh tree is returned. Simplification is applied
    /// bottom-up and repeated until a fixpoint (or maxIterations) is reached.
    ///
    /// @code
    /// auto x = symbol("x");
    /// auto f = (x + constant(0.0)) * constant(1.0) + (constant(2.0) + constant(3.0));
    /// auto s = simplify(f);           // s == x + 5
    /// @endcode
    ///
    /// @param expr The expression to simplify (must be non-null)
    /// @param opts Simplification options
    /// @return A new, simplified expression tree
    //////////////////////////////////////////////////////////////////////////////////////////
    MML_SYMBOLIC_API ExprPtr simplify(const ExprPtr& expr, const SimplifyOptions& opts = {});

} // namespace MML::Symbolic

#endif // MML_SYMBOLIC_SIMPLIFIER_H
