///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Jacobian.h                                                          ///
///  Description: Jacobian and Hessian matrix computation using symbolic and AD       ///
///               methods                                                             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_SYMBOLIC_JACOBIAN_H
#define MML_SYMBOLIC_JACOBIAN_H

#include "Expr.h"
#include <mml/core/Derivation/ForwardAD.h>
#include <mml/core/Derivation/ReverseAD.h>

#include <vector>
#include <functional>
#include <stdexcept>

namespace MML::Symbolic
{
    using AD::Dual;
    using AD::Tape;
    using AD::ADVar;

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Simple matrix class for Jacobian/Hessian storage
    /// @tparam T Element type (typically double)
    /// 
    /// Provides basic row-major matrix storage without pulling in full MML Matrix
    //////////////////////////////////////////////////////////////////////////////////////////
    template<typename T>
    class DenseMatrix
    {
    private:
        std::vector<T> data_;
        size_t rows_;
        size_t cols_;

    public:
        DenseMatrix() : rows_(0), cols_(0) {}
        
        DenseMatrix(size_t rows, size_t cols, T init = T(0))
            : data_(rows * cols, init), rows_(rows), cols_(cols) {}
        
        size_t rows() const { return rows_; }
        size_t cols() const { return cols_; }
        
        T& operator()(size_t i, size_t j) {
            return data_[i * cols_ + j];
        }
        
        const T& operator()(size_t i, size_t j) const {
            return data_[i * cols_ + j];
        }
        
        /// Get row as vector
        std::vector<T> row(size_t i) const {
            std::vector<T> r(cols_);
            for (size_t j = 0; j < cols_; ++j)
                r[j] = (*this)(i, j);
            return r;
        }
        
        /// Get column as vector
        std::vector<T> col(size_t j) const {
            std::vector<T> c(rows_);
            for (size_t i = 0; i < rows_; ++i)
                c[i] = (*this)(i, j);
            return c;
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Jacobian and Hessian computation utilities
    /// 
    /// Provides multiple approaches:
    /// - Symbolic differentiation of expression trees
    /// - Forward-mode AD (efficient when n < m for f: R^n → R^m)
    /// - Reverse-mode AD (efficient when n > m for f: R^n → R^m)
    /// - Automatic mode selection based on dimensions
    //////////////////////////////////////////////////////////////////////////////////////////
    class MML_SYMBOLIC_API JacobianComputer
    {
    public:
        //////////////////////////////////////////////////////////////////////////////////
        /// @brief Compute Jacobian matrix symbolically
        /// @param funcs Vector of m expression trees (output functions)
        /// @param vars Variable names (n inputs)
        /// @param bindings Variable values for evaluation
        /// @return m×n Jacobian matrix where J[i][j] = ∂f_i/∂x_j
        //////////////////////////////////////////////////////////////////////////////////
        static DenseMatrix<double> jacobianSymbolic(
            const std::vector<ExprPtr>& funcs,
            const std::vector<std::string>& vars,
            const VarBindings& bindings);

        //////////////////////////////////////////////////////////////////////////////////
        /// @brief Compute Jacobian using forward-mode AD
        /// @tparam F Function type taking vector<Dual<double>> and returning vector<Dual<double>>
        /// @param func The function to differentiate
        /// @param point Evaluation point
        /// @return m×n Jacobian matrix
        /// 
        /// Runs n forward passes (one per input variable)
        /// Efficient when n < m (few inputs, many outputs)
        //////////////////////////////////////////////////////////////////////////////////
        template<typename F>
        static DenseMatrix<double> jacobianForwardAD(
            F func,
            const std::vector<double>& point)
        {
            size_t n = point.size();
            
            // Determine output size by running with all-zero derivatives
            std::vector<Dual<double>> testInput(n);
            for (size_t i = 0; i < n; ++i)
                testInput[i] = Dual<double>(point[i], 0.0);
            auto testOutput = func(testInput);
            size_t m = testOutput.size();
            
            DenseMatrix<double> J(m, n);
            
            // For each input variable, set its derivative to 1 and propagate
            for (size_t j = 0; j < n; ++j) {
                std::vector<Dual<double>> input(n);
                for (size_t k = 0; k < n; ++k) {
                    input[k] = Dual<double>(point[k], (k == j) ? 1.0 : 0.0);
                }
                
                auto output = func(input);
                
                for (size_t i = 0; i < m; ++i) {
                    J(i, j) = output[i].deriv;
                }
            }
            
            return J;
        }

        //////////////////////////////////////////////////////////////////////////////////
        /// @brief Compute Jacobian using reverse-mode AD
        /// @tparam F Function type compatible with Tape-based AD
        /// @param func The function (takes Tape&, returns vector<ADVar>)
        /// @param point Evaluation point
        /// @return m×n Jacobian matrix
        /// 
        /// Runs m backward passes (one per output)
        /// Efficient when n > m (many inputs, few outputs)
        //////////////////////////////////////////////////////////////////////////////////
        static DenseMatrix<double> jacobianReverseAD(
            std::function<std::vector<ADVar>(Tape&, const std::vector<ADVar>&)> func,
            const std::vector<double>& point);

        //////////////////////////////////////////////////////////////////////////////////
        /// @brief Compute Hessian matrix symbolically
        /// @param func Scalar function expression
        /// @param vars Variable names
        /// @param bindings Variable values
        /// @return n×n symmetric Hessian matrix where H[i][j] = ∂²f/∂x_i∂x_j
        //////////////////////////////////////////////////////////////////////////////////
        static DenseMatrix<double> hessianSymbolic(
            const ExprPtr& func,
            const std::vector<std::string>& vars,
            const VarBindings& bindings);

        //////////////////////////////////////////////////////////////////////////////////
        /// @brief Compute Hessian using forward-over-forward AD (nested duals)
        /// @tparam F Function type taking Dual<Dual<double>> and returning same
        /// @param func Scalar function
        /// @param point Evaluation point
        /// @return n×n Hessian matrix
        /// 
        /// Uses Dual<Dual<double>> for second derivatives
        /// Requires n² forward passes
        //////////////////////////////////////////////////////////////////////////////////
        template<typename F>
        static DenseMatrix<double> hessianForwardAD(
            F func,
            const std::vector<double>& point)
        {
            using Dual2 = Dual<Dual<double>>;
            size_t n = point.size();
            DenseMatrix<double> H(n, n);
            
            // For each pair (i, j), compute ∂²f/∂x_i∂x_j
            // Exploit symmetry: only compute upper triangle
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i; j < n; ++j) {
                    std::vector<Dual2> input(n);
                    
                    for (size_t k = 0; k < n; ++k) {
                        // Inner dual: derivative w.r.t. x_j
                        // Outer dual: derivative w.r.t. x_i
                        double inner_val = point[k];
                        double inner_deriv = (k == j) ? 1.0 : 0.0;
                        double outer_deriv_val = (k == i) ? 1.0 : 0.0;
                        double outer_deriv_deriv = 0.0;
                        
                        input[k] = Dual2(
                            Dual<double>(inner_val, inner_deriv),
                            Dual<double>(outer_deriv_val, outer_deriv_deriv)
                        );
                    }
                    
                    Dual2 result = func(input);
                    
                    // The second derivative is in result.deriv.deriv
                    double val = result.deriv.deriv;
                    H(i, j) = val;
                    if (i != j) H(j, i) = val;  // Symmetry
                }
            }
            
            return H;
        }

        //////////////////////////////////////////////////////////////////////////////////
        /// @brief Compute gradient (special case: Jacobian of scalar function)
        /// @param func Scalar expression
        /// @param vars Variable names
        /// @param bindings Variable values
        /// @return Gradient vector [∂f/∂x_1, ∂f/∂x_2, ..., ∂f/∂x_n]
        //////////////////////////////////////////////////////////////////////////////////
        static std::vector<double> gradientSymbolic(
            const ExprPtr& func,
            const std::vector<std::string>& vars,
            const VarBindings& bindings);

        //////////////////////////////////////////////////////////////////////////////////
        /// @brief Compute gradient using reverse-mode AD
        /// @param func Function taking Tape& and vector<ADVar>, returning ADVar
        /// @param point Evaluation point
        /// @return Gradient vector (single backward pass)
        //////////////////////////////////////////////////////////////////////////////////
        static std::vector<double> gradientReverseAD(
            std::function<ADVar(Tape&, const std::vector<ADVar>&)> func,
            const std::vector<double>& point);
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    /// @brief Convenience functions for Jacobian/Hessian computation
    //////////////////////////////////////////////////////////////////////////////////////////
    
    /// Compute Jacobian of vector function symbolically
    inline DenseMatrix<double> jacobian(
        const std::vector<ExprPtr>& funcs,
        const std::vector<std::string>& vars,
        const VarBindings& bindings)
    {
        return JacobianComputer::jacobianSymbolic(funcs, vars, bindings);
    }
    
    /// Compute Hessian of scalar function symbolically
    inline DenseMatrix<double> hessian(
        const ExprPtr& func,
        const std::vector<std::string>& vars,
        const VarBindings& bindings)
    {
        return JacobianComputer::hessianSymbolic(func, vars, bindings);
    }

} // namespace MML::Symbolic

#endif // MML_SYMBOLIC_JACOBIAN_H
