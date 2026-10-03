///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLAdapters.h                                                       ///
///  Description: Adapters to integrate symbolic expressions with MML function        ///
///               interfaces (IRealFunction, IScalarFunction, IVectorFunction)        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_ADAPTERS_H
#define MML_SYMBOLIC_ADAPTERS_H

#include "Expr.h"
#include "Functions.h"
#include <mml/core/Derivation/ForwardAD.h>
#include "Jacobian.h"
#include "interfaces/IFunction.h"

#include <array>
#include <stdexcept>

namespace MML::Symbolic
{
    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Adapts a symbolic expression to IRealFunction interface
    /// 
    /// Allows using symbolic expressions anywhere MML expects an IRealFunction.
    /// Key benefit: provides EXACT analytic derivatives instead of numerical approximation.
    /// 
    /// Example:
    /// @code
    /// auto x = symbol("x");
    /// auto expr = sin_(x) * exp_(-x);
    /// SymbolicRealFunction f(expr, "x");
    /// 
    /// double val = f(1.5);              // Evaluate at x=1.5
    /// double deriv = f.derivative(1.5); // EXACT derivative at x=1.5
    /// 
    /// // Use with MML algorithms
    /// RootFinding::NewtonRaphson(f, 0.0, 2.0);  // Works!
    /// @endcode
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API SymbolicRealFunction : public IRealFunction
    {
    private:
        ExprPtr expr_;
        std::string var_;
        mutable ExprPtr derivative_expr_;  // Cached derivative
        
    public:
        /// Constructor
        /// @param expr The symbolic expression
        /// @param var Variable name (default "x")
        SymbolicRealFunction(ExprPtr expr, const std::string& var = "x");
        
        /// Evaluate the function
        Real operator()(Real x) const override;
        
        /// Get the underlying expression
        ExprPtr expression() const;
        
        /// Get the variable name
        const std::string& variable() const;
        
        /// Get the symbolic derivative expression
        ExprPtr derivativeExpr() const;
        
        /// Evaluate the derivative at a point (EXACT, not numerical!)
        Real derivative(Real x) const;
        
        /// Get a new SymbolicRealFunction representing the derivative
        SymbolicRealFunction getDerivative() const;
        
        /// Get the second derivative expression
        ExprPtr secondDerivativeExpr() const;
        
        /// Evaluate the second derivative
        Real secondDerivative(Real x) const;
        
        /// String representation of the expression
        std::string toString() const;
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Adapts a symbolic expression to IScalarFunction<N> interface
    /// 
    /// For multi-variable scalar functions f: R^N → R
    /// Provides analytic gradient computation.
    /// 
    /// Example:
    /// @code
    /// auto x = symbol("x");
    /// auto y = symbol("y");
    /// auto expr = x*x + y*y;  // f(x,y) = x² + y²
    /// SymbolicScalarFunction<2> f(expr, {"x", "y"});
    /// 
    /// VectorN<Real,2> point = {3.0, 4.0};
    /// double val = f(point);                    // = 25
    /// VectorN<Real,2> grad = f.gradient(point); // = [6, 8]
    /// @endcode
    //////////////////////////////////////////////////////////////////////////////////////////
    template<int N>
    class SymbolicScalarFunction : public IScalarFunction<N>
    {
    private:
        ExprPtr expr_;
        std::array<std::string, N> vars_;
        mutable std::array<ExprPtr, N> gradient_exprs_;  // Cached gradient
        mutable bool gradient_computed_ = false;
        
    public:
        /// Constructor
        /// @param expr The symbolic expression
        /// @param vars Array of variable names
        SymbolicScalarFunction(ExprPtr expr, const std::array<std::string, N>& vars)
            : expr_(expr), vars_(vars) {}
        
        /// Constructor from initializer list
        SymbolicScalarFunction(ExprPtr expr, const std::vector<std::string>& vars)
            : expr_(expr) {
            if (vars.size() != N) {
                throw std::invalid_argument("Variable count must match template dimension N");
            }
            std::copy(vars.begin(), vars.end(), vars_.begin());
        }
        
        /// Evaluate the function
        Real operator()(const VectorN<Real, N>& x) const override {
            VarBindings bindings;
            for (int i = 0; i < N; ++i) {
                bindings[vars_[i]] = x[i];
            }
            return expr_->eval(bindings);
        }
        
        /// Compute gradient at a point (EXACT!)
        VectorN<Real, N> gradient(const VectorN<Real, N>& x) const {
            computeGradientExprs();
            
            VarBindings bindings;
            for (int i = 0; i < N; ++i) {
                bindings[vars_[i]] = x[i];
            }
            
            VectorN<Real, N> grad;
            for (int i = 0; i < N; ++i) {
                grad[i] = gradient_exprs_[i]->eval(bindings);
            }
            return grad;
        }
        
        /// Get partial derivative expression w.r.t. variable i
        ExprPtr partialDerivativeExpr(int i) const {
            computeGradientExprs();
            return gradient_exprs_[i];
        }
        
        /// Compute Hessian matrix at a point
        DenseMatrix<Real> hessian(const VectorN<Real, N>& x) const {
            std::vector<std::string> varVec(vars_.begin(), vars_.end());
            VarBindings bindings;
            for (int i = 0; i < N; ++i) {
                bindings[vars_[i]] = x[i];
            }
            return JacobianComputer::hessianSymbolic(expr_, varVec, bindings);
        }
        
        /// Get underlying expression
        ExprPtr expression() const { return expr_; }
        
        /// Get variable names
        const std::array<std::string, N>& variables() const { return vars_; }
        
    private:
        void computeGradientExprs() const {
            if (!gradient_computed_) {
                for (int i = 0; i < N; ++i) {
                    gradient_exprs_[i] = expr_->diff(vars_[i]);
                }
                gradient_computed_ = true;
            }
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Adapts a vector of symbolic expressions to IVectorFunction<N> interface
    /// 
    /// For vector functions f: R^N → R^N
    /// Provides analytic Jacobian computation.
    /// 
    /// Example:
    /// @code
    /// auto x = symbol("x");
    /// auto y = symbol("y");
    /// // f(x,y) = [x*y, x² + y²]
    /// SymbolicVectorFunction<2> f({x*y, x*x + y*y}, {"x", "y"});
    /// 
    /// VectorN<Real,2> point = {1.0, 2.0};
    /// VectorN<Real,2> val = f(point);           // = [2, 5]
    /// auto J = f.jacobian(point);               // 2x2 Jacobian
    /// @endcode
    //////////////////////////////////////////////////////////////////////////////////////////
    template<int N>
    class SymbolicVectorFunction : public IVectorFunction<N>
    {
    private:
        std::array<ExprPtr, N> exprs_;
        std::array<std::string, N> vars_;
        
    public:
        /// Constructor
        /// @param exprs Array of N expressions (output functions)
        /// @param vars Array of N variable names (input variables)
        SymbolicVectorFunction(const std::array<ExprPtr, N>& exprs, 
                               const std::array<std::string, N>& vars)
            : exprs_(exprs), vars_(vars) {}
        
        /// Constructor from vectors
        SymbolicVectorFunction(const std::vector<ExprPtr>& exprs,
                               const std::vector<std::string>& vars) {
            if (exprs.size() != N || vars.size() != N) {
                throw std::invalid_argument("Expression and variable counts must match N");
            }
            std::copy(exprs.begin(), exprs.end(), exprs_.begin());
            std::copy(vars.begin(), vars.end(), vars_.begin());
        }
        
        /// Evaluate the function
        VectorN<Real, N> operator()(const VectorN<Real, N>& x) const override {
            VarBindings bindings;
            for (int i = 0; i < N; ++i) {
                bindings[vars_[i]] = x[i];
            }
            
            VectorN<Real, N> result;
            for (int i = 0; i < N; ++i) {
                result[i] = exprs_[i]->eval(bindings);
            }
            return result;
        }
        
        /// Compute Jacobian matrix at a point (EXACT!)
        DenseMatrix<Real> jacobian(const VectorN<Real, N>& x) const {
            std::vector<ExprPtr> exprVec(exprs_.begin(), exprs_.end());
            std::vector<std::string> varVec(vars_.begin(), vars_.end());
            
            VarBindings bindings;
            for (int i = 0; i < N; ++i) {
                bindings[vars_[i]] = x[i];
            }
            
            // JacobianComputer works in double; convert to DenseMatrix<Real>
            auto jac = JacobianComputer::jacobianSymbolic(exprVec, varVec, bindings);
            DenseMatrix<Real> result(jac.rows(), jac.cols());
            for (size_t i = 0; i < jac.rows(); ++i)
                for (size_t j = 0; j < jac.cols(); ++j)
                    result(i, j) = static_cast<Real>(jac(i, j));
            return result;
        }
        
        /// Get expression for output component i
        ExprPtr expression(int i) const { return exprs_[i]; }
        
        /// Get all expressions
        const std::array<ExprPtr, N>& expressions() const { return exprs_; }
        
        /// Get variable names
        const std::array<std::string, N>& variables() const { return vars_; }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Adapts symbolic expressions to IVectorFunctionNM<N,M> interface
    /// 
    /// For vector functions f: R^N → R^M where N ≠ M
    /// 
    /// Example:
    /// @code
    /// auto x = symbol("x");
    /// auto y = symbol("y");
    /// auto z = symbol("z");
    /// // f: R³ → R² = [x+y+z, x*y*z]
    /// SymbolicVectorFunctionNM<3,2> f({x+y+z, x*y*z}, {"x","y","z"});
    /// @endcode
    //////////////////////////////////////////////////////////////////////////////////////////
    template<int N, int M>
    class SymbolicVectorFunctionNM : public IVectorFunctionNM<N, M>
    {
    private:
        std::array<ExprPtr, M> exprs_;
        std::array<std::string, N> vars_;
        
    public:
        /// Constructor
        SymbolicVectorFunctionNM(const std::array<ExprPtr, M>& exprs,
                                 const std::array<std::string, N>& vars)
            : exprs_(exprs), vars_(vars) {}
        
        /// Constructor from vectors
        SymbolicVectorFunctionNM(const std::vector<ExprPtr>& exprs,
                                 const std::vector<std::string>& vars) {
            if (exprs.size() != M || vars.size() != N) {
                throw std::invalid_argument("Expression count must be M, variable count must be N");
            }
            std::copy(exprs.begin(), exprs.end(), exprs_.begin());
            std::copy(vars.begin(), vars.end(), vars_.begin());
        }
        
        /// Evaluate the function
        VectorN<Real, M> operator()(const VectorN<Real, N>& x) const override {
            VarBindings bindings;
            for (int i = 0; i < N; ++i) {
                bindings[vars_[i]] = x[i];
            }
            
            VectorN<Real, M> result;
            for (int i = 0; i < M; ++i) {
                result[i] = exprs_[i]->eval(bindings);
            }
            return result;
        }
        
        /// Compute M×N Jacobian matrix
        DenseMatrix<Real> jacobian(const VectorN<Real, N>& x) const {
            std::vector<ExprPtr> exprVec(exprs_.begin(), exprs_.end());
            std::vector<std::string> varVec(vars_.begin(), vars_.end());
            
            VarBindings bindings;
            for (int i = 0; i < N; ++i) {
                bindings[vars_[i]] = x[i];
            }
            
            // JacobianComputer works in double; convert to DenseMatrix<Real>
            auto jac = JacobianComputer::jacobianSymbolic(exprVec, varVec, bindings);
            DenseMatrix<Real> result(jac.rows(), jac.cols());
            for (size_t i = 0; i < jac.rows(); ++i)
                for (size_t j = 0; j < jac.cols(); ++j)
                    result(i, j) = static_cast<Real>(jac(i, j));
            return result;
        }
        
        /// Get expressions
        const std::array<ExprPtr, M>& expressions() const { return exprs_; }
        
        /// Get variable names
        const std::array<std::string, N>& variables() const { return vars_; }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Utility functions for creating adapters
    //////////////////////////////////////////////////////////////////////////////////////////
    
    /// Create IRealFunction from expression
    inline SymbolicRealFunction makeRealFunction(ExprPtr expr, const std::string& var = "x") {
        return SymbolicRealFunction(expr, var);
    }
    
    /// Create IScalarFunction<2> from expression
    inline SymbolicScalarFunction<2> makeScalarFunction2D(
        ExprPtr expr, 
        const std::string& var1 = "x", 
        const std::string& var2 = "y") {
        std::array<std::string, 2> vars = {var1, var2};
        return SymbolicScalarFunction<2>(expr, vars);
    }
    
    /// Create IScalarFunction<3> from expression
    inline SymbolicScalarFunction<3> makeScalarFunction3D(
        ExprPtr expr,
        const std::string& var1 = "x",
        const std::string& var2 = "y",
        const std::string& var3 = "z") {
        std::array<std::string, 3> vars = {var1, var2, var3};
        return SymbolicScalarFunction<3>(expr, vars);
    }

} // namespace MML::Symbolic

#endif // MML_SYMBOLIC_ADAPTERS_H
