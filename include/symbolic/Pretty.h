///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Pretty.h                                                            ///
///  Description: Precedence-aware infix pretty-printing for symbolic expressions.     ///
///               Produces minimal-parenthesized, re-parseable output such as          ///
///               "3 * x ^ 2 + 4 * x - 6" (contrast Expr::toString(), which fully      ///
///               parenthesizes).                                                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_PRETTY_H
#define MML_SYMBOLIC_PRETTY_H

#include <mml/mml_export.h>

#include "Expr.h"

#include <string>

namespace MML::Symbolic
{
    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Render an expression as a minimal-parenthesized infix string.
    ///
    /// Uses operator precedence (+,- = 1; *,/ = 2; ^ = 3) so only necessary parentheses are
    /// emitted, e.g. `3 * x ^ 2 + 4 * x - 6`. The result is re-parseable by parse().
    /// This differs from Expr::toString(), which fully parenthesizes every operation.
    //////////////////////////////////////////////////////////////////////////////////////////
    MML_SYMBOLIC_API std::string toInfix(const ExprPtr& expr);

} // namespace MML::Symbolic

#endif // MML_SYMBOLIC_PRETTY_H
