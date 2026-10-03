///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Optimization/LP/RevisedSimplexSolver.h                                      ///
///  Description: Revised simplex solver using explicit basis inverse                      ///
///               PRIMARY production LP solver - SolveLP dispatches here                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LP_REVISED_SIMPLEX_SOLVER_H
#define MML_LP_REVISED_SIMPLEX_SOLVER_H

#include <mml/algorithms/Optimization/LP/LinearProgram.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "mml/mml_export.h"

namespace MML::Optimization {

class MML_OPTIMIZATION_API RevisedSimplexSolver {
private:
    LPConfig _config;
    
    // Problem data (standard form)
    Matrix<Real> _A;                 ///< Constraint matrix (m x n)
    Vector<Real> _b;                 ///< RHS vector (m)
    Vector<Real> _c;                 ///< Objective coefficients (n)
    int _m;                          ///< Number of constraints
    int _n;                          ///< Number of variables (including slack/artificial)
    int _numOriginalVars;            ///< Original problem variables
    int _numSlack;                   ///< Slack variables
    int _numArtificial;              ///< Artificial variables
    std::vector<int> _artificialIndices;
    
    // Basis representation
    Matrix<Real> _Binv;              ///< Basis inverse B^(-1) (m x m)
    std::vector<int> _basis;         ///< Indices of basic variables
    std::vector<int> _nonbasis;      ///< Indices of non-basic variables
    Vector<Real> _xB;                ///< Basic variable values: x_B = B^(-1) * b
    Vector<Real> _cB;                ///< Objective coefficients for basic variables
    
    // Statistics
    int _iterations = 0;
    int _phase = 0;
    
    // Eta file for product-form updates (for potential future optimization)
    // struct EtaMatrix { int col; Vector<Real> eta; };
    // std::vector<EtaMatrix> _etaFile;
    
public:
    RevisedSimplexSolver() = default;
    explicit RevisedSimplexSolver(const LPConfig& config);
    
    void SetConfig(const LPConfig& config);
    const LPConfig& GetConfig() const;
    int Iterations() const;
    
    /// @brief Solve a linear program using revised simplex method
    LPResult Solve(const LinearProgram& lp);
    
private:
    //---------------------------------------------------------------------------------
    // Initialization
    //---------------------------------------------------------------------------------
    
    /// @brief Initialize basis with identity (slack/artificial variables)
    void InitializeBasis();
    
    /// @brief Rebuild non-basis list excluding artificials (for Phase 2)
    void RebuildNonbasis();
    
    /// @brief Update c_B vector from current basis
    void UpdateCB();
    
    //---------------------------------------------------------------------------------
    // Core Revised Simplex Operations
    //---------------------------------------------------------------------------------
    
    /// @brief Compute basic solution: x_B = B^(-1) * b
    void ComputeBasicSolution();
    
    /// @brief Compute reduced cost for variable j: c_j - c_B * B^(-1) * A_j
    Real ComputeReducedCost(int j) const;
    
    /// @brief Compute column of B^(-1) * A for entering variable: d = B^(-1) * A_j
    Vector<Real> ComputeDirection(int j) const;
    
    /// @brief Compute objective value: c_B * x_B
    Real ComputeObjectiveValue() const;
    
    /// @brief Select entering variable (Dantzig's or Bland's rule)
    /// @return Column index, or -1 if optimal
    int SelectEnteringVariable() const;
    
    /// @brief Select leaving variable using minimum ratio test
    /// @param direction The direction vector d = B^(-1) * A_j
    /// @return Row index in basis, or -1 if unbounded
    int SelectLeavingVariable(const Vector<Real>& direction) const;
    
    /// @brief Update B^(-1) after a basis change using eta matrix multiplication
    /// 
    /// When column s leaves and column e enters, we update:
    /// B_new^(-1) = E * B_old^(-1)
    /// where E is an eta matrix (identity with column s replaced by -d/d_s with 1/d_s on diagonal)
    ///
    /// @param leavingRow Row index of leaving variable in basis
    /// @param direction Direction vector d = B^(-1) * A_entering
    void UpdateBasisInverse(int leavingRow, const Vector<Real>& direction);
    
    /// @brief Perform a pivot operation
    void Pivot(int enteringCol, int leavingRow, const Vector<Real>& direction);
    
    //---------------------------------------------------------------------------------
    // Main Loop
    //---------------------------------------------------------------------------------
    
    /// @brief Run revised simplex iterations
    LPStatus RunRevisedSimplex();
    
    //---------------------------------------------------------------------------------
    // Helpers
    //---------------------------------------------------------------------------------
    
    bool IsArtificial(int j) const;
    
    bool HasArtificialInBasis() const;
    
    /// @brief Try to pivot out artificial variables at zero level
    bool PivotOutArtificials();
    
    /// @brief Extract solution from current basis
    void ExtractSolution(const LinearProgram& lp, LPResult& result) const;
    
public:
    //---------------------------------------------------------------------------------
    // Debug Output
    //---------------------------------------------------------------------------------
    
    void PrintBasis(std::ostream& os = std::cout) const;
    
    void PrintBasisInverse(std::ostream& os = std::cout) const;
};


} // namespace MML::Optimization
#endif // MML_LP_REVISED_SIMPLEX_SOLVER_H
