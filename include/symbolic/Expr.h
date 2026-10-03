///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Expr.h                                                              ///
///  Description: Core expression tree infrastructure for symbolic computing          ///
///               Defines Expr base class, Symbol, Const, and expression operators   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_EXPR_H
#define MML_SYMBOLIC_EXPR_H

#include <mml/mml_export.h>

#include <memory>
#include <string>
#include <map>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace MML::Symbolic
{
    // Forward declarations
    class Expr;
    class Add;
    class Sub;
    class Mul;
    class Div;
    class Neg;

    /// Type alias for variable bindings (variable name -> value)
    using VarBindings = std::map<std::string, double>;
    
    /// Shared pointer to expression (enables tree structure with shared ownership)
    using ExprPtr = std::shared_ptr<Expr>;

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Abstract base class for all expression nodes in the symbolic expression tree
    /// 
    /// The expression tree represents mathematical expressions symbolically, enabling:
    /// - Symbolic differentiation
    /// - Expression evaluation with variable bindings
    /// - Code generation
    /// - Algebraic simplification
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Expr : public std::enable_shared_from_this<Expr>
    {
    public:
        virtual ~Expr();

        /// Evaluate the expression with given variable bindings
        /// @param bindings Map of variable names to their values
        /// @return The computed value
        /// @throws std::runtime_error if a required variable is not bound
        virtual double eval(const VarBindings& bindings) const = 0;

        /// Compute the symbolic derivative with respect to a variable
        /// @param var The variable name to differentiate with respect to
        /// @return A new expression tree representing the derivative
        virtual ExprPtr diff(const std::string& var) const = 0;

        /// Create a deep copy of this expression tree
        virtual ExprPtr clone() const = 0;

        /// Convert expression to a human-readable string
        virtual std::string toString() const = 0;

        /// Check if this expression is a constant (number)
        virtual bool isConstant() const;

        /// Check if this expression is a symbol (variable)
        virtual bool isSymbol() const;

        /// Get constant value (only valid if isConstant() returns true)
        virtual double getConstantValue() const;

        /// Get symbol name (only valid if isSymbol() returns true)
        virtual std::string getSymbolName() const;

        /// Collect all variable names used in this expression
        virtual void collectVariables(std::vector<std::string>& vars) const = 0;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Represents a numeric constant in the expression tree
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Const : public Expr
    {
    private:
        double value_;

    public:
        explicit Const(double value);

        double eval(const VarBindings&) const override;

        ExprPtr diff(const std::string&) const override;

        ExprPtr clone() const override;

        std::string toString() const override;

        bool isConstant() const override;
        double getConstantValue() const override;

        void collectVariables(std::vector<std::string>&) const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Represents a named variable (symbol) in the expression tree
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Symbol : public Expr
    {
    private:
        std::string name_;

    public:
        explicit Symbol(const std::string& name);

        double eval(const VarBindings& bindings) const override;

        ExprPtr diff(const std::string& var) const override;

        ExprPtr clone() const override;

        std::string toString() const override;

        bool isSymbol() const override;
        std::string getSymbolName() const override;

        const std::string& name() const;

        void collectVariables(std::vector<std::string>& vars) const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Unary negation operation (-x)
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Neg : public Expr
    {
    private:
        ExprPtr operand_;

    public:
        explicit Neg(ExprPtr operand);

        double eval(const VarBindings& bindings) const override;

        ExprPtr diff(const std::string& var) const override;

        ExprPtr clone() const override;

        std::string toString() const override;

        const ExprPtr& operand() const;

        void collectVariables(std::vector<std::string>& vars) const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Base class for binary operations (Add, Sub, Mul, Div, etc.)
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API BinaryOp : public Expr
    {
    protected:
        ExprPtr left_;
        ExprPtr right_;

    public:
        BinaryOp(ExprPtr left, ExprPtr right);

        const ExprPtr& left() const;
        const ExprPtr& right() const;

        void collectVariables(std::vector<std::string>& vars) const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Addition operation (a + b)
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Add : public BinaryOp
    {
    public:
        using BinaryOp::BinaryOp;

        double eval(const VarBindings& bindings) const override;

        ExprPtr diff(const std::string& var) const override;

        ExprPtr clone() const override;

        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Subtraction operation (a - b)
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Sub : public BinaryOp
    {
    public:
        using BinaryOp::BinaryOp;

        double eval(const VarBindings& bindings) const override;

        ExprPtr diff(const std::string& var) const override;

        ExprPtr clone() const override;

        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Multiplication operation (a * b)
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Mul : public BinaryOp
    {
    public:
        using BinaryOp::BinaryOp;

        double eval(const VarBindings& bindings) const override;

        ExprPtr diff(const std::string& var) const override;

        ExprPtr clone() const override;

        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Division operation (a / b)
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API Div : public BinaryOp
    {
    public:
        using BinaryOp::BinaryOp;

        double eval(const VarBindings& bindings) const override;

        ExprPtr diff(const std::string& var) const override;

        ExprPtr clone() const override;

        std::string toString() const override;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // Factory functions for cleaner expression building
    //////////////////////////////////////////////////////////////////////////////////////////

    /// Create a constant expression
    inline ExprPtr constant(double value) {
        return std::make_shared<Const>(value);
    }

    /// Create a symbol (variable) expression
    inline ExprPtr symbol(const std::string& name) {
        return std::make_shared<Symbol>(name);
    }

    /// Shorthand: create symbol 'x'
    inline ExprPtr x() { return symbol("x"); }
    /// Shorthand: create symbol 'y'
    inline ExprPtr y() { return symbol("y"); }
    /// Shorthand: create symbol 'z'
    inline ExprPtr z() { return symbol("z"); }
    /// Shorthand: create symbol 't'
    inline ExprPtr t() { return symbol("t"); }

    //////////////////////////////////////////////////////////////////////////////////////////
    // Operator overloads for intuitive expression building
    //////////////////////////////////////////////////////////////////////////////////////////

    inline ExprPtr operator+(ExprPtr lhs, ExprPtr rhs) {
        return std::make_shared<Add>(std::move(lhs), std::move(rhs));
    }

    inline ExprPtr operator-(ExprPtr lhs, ExprPtr rhs) {
        return std::make_shared<Sub>(std::move(lhs), std::move(rhs));
    }

    inline ExprPtr operator*(ExprPtr lhs, ExprPtr rhs) {
        return std::make_shared<Mul>(std::move(lhs), std::move(rhs));
    }

    inline ExprPtr operator/(ExprPtr lhs, ExprPtr rhs) {
        return std::make_shared<Div>(std::move(lhs), std::move(rhs));
    }

    inline ExprPtr operator-(ExprPtr operand) {
        return std::make_shared<Neg>(std::move(operand));
    }

    // Mixed operations with double on left side
    inline ExprPtr operator+(double lhs, ExprPtr rhs) {
        return std::make_shared<Add>(constant(lhs), std::move(rhs));
    }

    inline ExprPtr operator-(double lhs, ExprPtr rhs) {
        return std::make_shared<Sub>(constant(lhs), std::move(rhs));
    }

    inline ExprPtr operator*(double lhs, ExprPtr rhs) {
        return std::make_shared<Mul>(constant(lhs), std::move(rhs));
    }

    inline ExprPtr operator/(double lhs, ExprPtr rhs) {
        return std::make_shared<Div>(constant(lhs), std::move(rhs));
    }

    // Mixed operations with double on right side
    inline ExprPtr operator+(ExprPtr lhs, double rhs) {
        return std::make_shared<Add>(std::move(lhs), constant(rhs));
    }

    inline ExprPtr operator-(ExprPtr lhs, double rhs) {
        return std::make_shared<Sub>(std::move(lhs), constant(rhs));
    }

    inline ExprPtr operator*(ExprPtr lhs, double rhs) {
        return std::make_shared<Mul>(std::move(lhs), constant(rhs));
    }

    inline ExprPtr operator/(ExprPtr lhs, double rhs) {
        return std::make_shared<Div>(std::move(lhs), constant(rhs));
    }

} // namespace MML::Symbolic

#endif // MML_SYMBOLIC_EXPR_H
