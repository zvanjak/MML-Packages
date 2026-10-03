///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Parser.h                                                            ///
///  Description: Recursive-descent parser: expression string -> symbolic Expr tree.  ///
///                                                                                   ///
///  Grammar:                                                                         ///
///    expression = term (('+' | '-') term)*                                          ///
///    term       = power (('*' | '/') power)*                                        ///
///    power      = unary ('^' power)?          // right-associative                  ///
///    unary      = ('-' | '+')? primary                                             ///
///    primary    = NUMBER | CONST | VARIABLE | call | '(' expression ')'             ///
///    call       = IDENTIFIER '(' arglist? ')'                                       ///
///    arglist    = expression (',' expression)*                                      ///
///                                                                                   ///
///  Scalar-only engine: 'pi' and 'e' parse as constants; any other identifier is a   ///
///  symbol. Supported calls: sin cos tan asin acos atan sinh cosh tanh exp log/ln    ///
///  log10 sqrt abs (1 arg) and pow(base, exponent) (2 args). Vector/matrix literals  ///
///  and index access are intentionally NOT supported.                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_PARSER_H
#define MML_SYMBOLIC_PARSER_H

#include <mml/mml_export.h>

#include "Expr.h"

#include <stdexcept>
#include <string>

namespace MML::Symbolic
{
    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Exception thrown when parsing fails.
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API ParseError : public std::runtime_error
    {
    public:
        size_t position;   ///< Index in the input string where the error occurred

        ParseError(const std::string& message, size_t pos = 0)
            : std::runtime_error(message), position(pos) {}
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Parse an expression string into a symbolic Expr tree.
    ///
    /// @code
    /// auto f  = parse("sin(x^2) + 3*x");
    /// auto df = f->diff("x");
    /// double v = f->eval({{"x", 1.0}});
    /// @endcode
    ///
    /// @param input The expression string (e.g. "sin(x^2) + 3*x")
    /// @return Root of the parsed expression tree
    /// @throws ParseError on empty input, invalid characters, or unknown functions
    //////////////////////////////////////////////////////////////////////////////////////////
    MML_SYMBOLIC_API ExprPtr parse(const std::string& input);

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Parse without throwing; returns nullptr on any parse failure.
    //////////////////////////////////////////////////////////////////////////////////////////
    MML_SYMBOLIC_API ExprPtr tryParse(const std::string& input) noexcept;

} // namespace MML::Symbolic

#endif // MML_SYMBOLIC_PARSER_H
