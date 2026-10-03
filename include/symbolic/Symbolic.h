///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Symbolic.h                                                          ///
///  Description: Main include header for symbolic computing and automatic            ///
///               differentiation module                                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_H
#define MML_SYMBOLIC_H

///////////////////////////////////////////////////////////////////////////////////////////
/// @file Symbolic.h
/// @brief Symbolic Computing and Automatic Differentiation for MML
///
/// This module provides three approaches to computing derivatives:
///
/// 1. **Symbolic Differentiation** (Expr.h, Functions.h)
///    Build expression trees and derive new expressions symbolically.
///    @code
///    auto x = symbol("x");
///    auto f = sin_(x * x);           // f = sin(x²)
///    auto df = f->diff("x");         // df = 2x*cos(x²)
///    double val = df->eval({{"x", 2.0}});  // Evaluate at x=2
///    @endcode
///
/// 2. **Forward-Mode AD** (ForwardAD.h)
///    Use dual numbers for efficient derivative computation.
///    Best for f: R → R^m (few inputs, many outputs).
///    @code
///    Dual<double> x(2.0, 1.0);   // x=2, dx/dx=1
///    Dual<double> y = sin(x*x); // Automatically computes derivative
///    // y.value = sin(4), y.deriv = cos(4)*2*2
///    @endcode
///
/// 3. **Reverse-Mode AD** (ReverseAD.h)
///    Tape-based backpropagation for gradient computation.
///    Best for f: R^n → R (many inputs, one output).
///    @code
///    Tape tape;
///    ADVar x(&tape, 2.0), y(&tape, 3.0);
///    ADVar z = sin(x * y) + exp(x);
///    z.backward();
///    // x.adjoint() = cos(6)*3 + exp(2) = gradient w.r.t. x
///    // y.adjoint() = cos(6)*2 = gradient w.r.t. y
///    @endcode
///
/// **When to use which:**
/// - Symbolic: When you need the derivative formula, code generation, or simplification
/// - Forward AD: When computing directional derivatives or Jacobians with few inputs
/// - Reverse AD: When computing gradients of scalar functions (optimization, ML)
///
///////////////////////////////////////////////////////////////////////////////////////////

#include "Expr.h"
#include "Operations.h"
#include "Functions.h"
#include <mml/core/Derivation/ForwardAD.h>
#include <mml/core/Derivation/ReverseAD.h>
#include "Jacobian.h"
#include "MMLAdapters.h"
#include "Simplifier.h"
#include "Parser.h"
#include "Pretty.h"

namespace MML
{
    // Bring symbolic types into MML namespace for convenience
    using Symbolic::ExprPtr;
    using Symbolic::VarBindings;
    using Symbolic::Expr;
    using Symbolic::Const;
    using Symbolic::Symbol;
    
    // Forward AD
    using AD::Dual;
    using AD::DualD;
    using AD::DualF;
    
    // Reverse AD
    using AD::Tape;
    using AD::ADVar;
    
    // Factory functions
    using Symbolic::constant;
    using Symbolic::symbol;
    using Symbolic::x;
    using Symbolic::y;
    using Symbolic::z;
    using Symbolic::t;

    // Algebraic simplification
    using Symbolic::simplify;
    using Symbolic::SimplifyOptions;

    // String parsing
    using Symbolic::parse;
    using Symbolic::tryParse;
    using Symbolic::ParseError;

    // Pretty-printing (minimal-parenthesized infix)
    using Symbolic::toInfix;
    
    // Symbolic functions
    using Symbolic::pow_;
    using Symbolic::sqr;
    using Symbolic::cube;
    using Symbolic::sin_;
    using Symbolic::cos_;
    using Symbolic::tan_;
    using Symbolic::asin_;
    using Symbolic::acos_;
    using Symbolic::atan_;
    using Symbolic::sinh_;
    using Symbolic::cosh_;
    using Symbolic::tanh_;
    using Symbolic::exp_;
    using Symbolic::log_;
    using Symbolic::log10_;
    using Symbolic::sqrt_;
    using Symbolic::abs_;
    
    // Forward AD utility functions
    using AD::derivative;
    using AD::gradient;

    // Reverse AD utility functions
    using AD::gradientReverse;
    
    // Jacobian/Hessian computation
    using Symbolic::DenseMatrix;
    using Symbolic::JacobianComputer;
    using Symbolic::jacobian;
    using Symbolic::hessian;
    
    // MML interface adapters
    using Symbolic::SymbolicRealFunction;
    using Symbolic::makeRealFunction;
    // SymbolicScalarFunction<N> and SymbolicVectorFunction<N> are templates
    // and should be accessed via MML::Symbolic:: namespace

} // namespace MML

#endif // MML_SYMBOLIC_H
