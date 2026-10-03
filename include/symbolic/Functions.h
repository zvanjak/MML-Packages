///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Functions.h                                                         ///
///  Description: Transcendental function expression nodes (sin, cos, exp, log, etc.) ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_FUNCTIONS_H
#define MML_SYMBOLIC_FUNCTIONS_H

#include "Expr.h"
#include "Operations.h"

namespace MML::Symbolic
{
    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Base class for unary mathematical functions
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API UnaryFunc : public Expr
    {
    protected:
        ExprPtr arg_;

    public:
        explicit UnaryFunc(ExprPtr arg);
        const ExprPtr& arg() const;
        void collectVariables(std::vector<std::string>& vars) const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // Trigonometric Functions
    //////////////////////////////////////////////////////////////////////////////////////////

    /// @brief Sine function: sin(x)
    class MML_SYMBOLIC_API Sin : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Cosine function: cos(x)
    class MML_SYMBOLIC_API Cos : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Tangent function: tan(x)
    class MML_SYMBOLIC_API Tan : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // Inverse Trigonometric Functions
    //////////////////////////////////////////////////////////////////////////////////////////

    /// @brief Arc sine function: asin(x)
    class MML_SYMBOLIC_API Asin : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Arc cosine function: acos(x)
    class MML_SYMBOLIC_API Acos : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Arc tangent function: atan(x)
    class MML_SYMBOLIC_API Atan : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // Hyperbolic Functions
    //////////////////////////////////////////////////////////////////////////////////////////

    /// @brief Hyperbolic sine: sinh(x)
    class MML_SYMBOLIC_API Sinh : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Hyperbolic cosine: cosh(x)
    class MML_SYMBOLIC_API Cosh : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Hyperbolic tangent: tanh(x)
    class MML_SYMBOLIC_API Tanh : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // Exponential and Logarithmic Functions
    //////////////////////////////////////////////////////////////////////////////////////////

    /// @brief Natural exponential: exp(x) = e^x
    class MML_SYMBOLIC_API Exp : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Natural logarithm: log(x) = ln(x)
    class MML_SYMBOLIC_API Log : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Base-10 logarithm: log10(x)
    class MML_SYMBOLIC_API Log10 : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // Algebraic Functions
    //////////////////////////////////////////////////////////////////////////////////////////

    /// @brief Square root: sqrt(x)
    class MML_SYMBOLIC_API Sqrt : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    /// @brief Absolute value: abs(x)
    class MML_SYMBOLIC_API Abs : public UnaryFunc
    {
    public:
        using UnaryFunc::UnaryFunc;

        double eval(const VarBindings& bindings) const override;
        ExprPtr diff(const std::string& var) const override;
        ExprPtr clone() const override;
        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // Deferred diff() implementations (to avoid circular dependencies)
    //////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////
    // Factory functions for clean syntax
    //////////////////////////////////////////////////////////////////////////////////////////

    inline ExprPtr sin_(ExprPtr x) { return std::make_shared<Sin>(std::move(x)); }
    inline ExprPtr cos_(ExprPtr x) { return std::make_shared<Cos>(std::move(x)); }
    inline ExprPtr tan_(ExprPtr x) { return std::make_shared<Tan>(std::move(x)); }

    inline ExprPtr asin_(ExprPtr x) { return std::make_shared<Asin>(std::move(x)); }
    inline ExprPtr acos_(ExprPtr x) { return std::make_shared<Acos>(std::move(x)); }
    inline ExprPtr atan_(ExprPtr x) { return std::make_shared<Atan>(std::move(x)); }

    inline ExprPtr sinh_(ExprPtr x) { return std::make_shared<Sinh>(std::move(x)); }
    inline ExprPtr cosh_(ExprPtr x) { return std::make_shared<Cosh>(std::move(x)); }
    inline ExprPtr tanh_(ExprPtr x) { return std::make_shared<Tanh>(std::move(x)); }

    inline ExprPtr exp_(ExprPtr x) { return std::make_shared<Exp>(std::move(x)); }
    inline ExprPtr log_(ExprPtr x) { return std::make_shared<Log>(std::move(x)); }
    inline ExprPtr log10_(ExprPtr x) { return std::make_shared<Log10>(std::move(x)); }

    inline ExprPtr sqrt_(ExprPtr x) { return std::make_shared<Sqrt>(std::move(x)); }
    inline ExprPtr abs_(ExprPtr x) { return std::make_shared<Abs>(std::move(x)); }

} // namespace MML::Symbolic

#endif // MML_SYMBOLIC_FUNCTIONS_H
