///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Operations.h                                                        ///
///  Description: Binary operation expression nodes (Pow) and utilities               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_OPERATIONS_H
#define MML_SYMBOLIC_OPERATIONS_H

#include "Expr.h"

namespace MML::Symbolic
{
    // Forward declarations for functions used in Pow::diff
    class Log;
    class Exp;

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Power operation (base ^ exponent)
    /// 
    /// Handles both integer powers (x^2) and general powers (x^y).
    /// Uses the generalized power rule for differentiation:
    /// d/dx(f^g) = f^g * (g' * ln(f) + g * f'/f)
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Pow : public BinaryOp
    {
    public:
        using BinaryOp::BinaryOp;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// Factory function for power expressions
    //////////////////////////////////////////////////////////////////////////////////////////
    inline ExprPtr pow_(ExprPtr base, ExprPtr exponent) {
        return std::make_shared<Pow>(std::move(base), std::move(exponent));
    }

    inline ExprPtr pow_(ExprPtr base, double exponent) {
        return std::make_shared<Pow>(std::move(base), constant(exponent));
    }

    inline ExprPtr pow_(double base, ExprPtr exponent) {
        return std::make_shared<Pow>(constant(base), std::move(exponent));
    }

    /// Square function (x^2)
    inline ExprPtr sqr(ExprPtr x) {
        return pow_(std::move(x), 2.0);
    }

    /// Cube function (x^3)
    inline ExprPtr cube(ExprPtr x) {
        return pow_(std::move(x), 3.0);
    }

} // namespace MML::Symbolic

#endif // MML_SYMBOLIC_OPERATIONS_H
