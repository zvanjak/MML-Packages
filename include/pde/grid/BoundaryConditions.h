#ifndef MML_BOUNDARY_CONDITIONS_H
#define MML_BOUNDARY_CONDITIONS_H

#include <functional>
#include <map>
#include <vector>
#include <stdexcept>

#include "Grid.h"
#include "Grid1D.h"
#include "Grid2D.h"
#include "Grid3D.h"
#include "GridFunction.h"
#include <mml/base/SparseMatrix/SparseMatrixCSR.h>

namespace MML::PDE
{
    // Bring Sparse namespace types into scope
    using MML::SparseMatrix::SparseMatrixCSR;
    
    ///////////////////////////////////////////////////////////////////////
    ///                     Boundary Condition Types                    ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Types of boundary conditions supported
     * 
     * - Dirichlet: u = g(x) - fixed value
     * - Neumann: ∂u/∂n = g(x) - fixed flux/gradient
     * - Robin: αu + β(∂u/∂n) = g(x) - mixed condition
     * - Periodic: u(x_min) = u(x_max) - wrapping
     */
    enum class BCType { Dirichlet, Neumann, Robin, Periodic };
    
    /**
     * @brief String representation of boundary condition type
     */
    inline const char* BCTypeToString(BCType type) {
        switch (type) {
            case BCType::Dirichlet: return "Dirichlet";
            case BCType::Neumann:   return "Neumann";
            case BCType::Robin:     return "Robin";
            case BCType::Periodic:  return "Periodic";
            default:                return "Unknown";
        }
    }

    ///////////////////////////////////////////////////////////////////////
    ///                   1D Boundary Condition                         ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Boundary condition specification for 1D problems
     * 
     * @tparam T Numeric type (float, double)
     * 
     * For Dirichlet: u = value(x)
     * For Neumann: du/dx = value(x)  (or -du/dx on right boundary)
     * For Robin: alpha*u + beta*du/dn = value(x)
     * For Periodic: u(left) = u(right)
     */
    template<typename T>
    struct BoundaryCondition1D {
        BCType type = BCType::Dirichlet;
        std::function<T(T)> value;      // g(x) - value or flux function
        T alpha = T(1);                  // Robin coefficient for u
        T beta = T(0);                   // Robin coefficient for du/dn
        
        // Convenience constructors
        static BoundaryCondition1D Dirichlet(std::function<T(T)> g) {
            BoundaryCondition1D bc;
            bc.type = BCType::Dirichlet;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition1D Dirichlet(T constantValue) {
            return Dirichlet([constantValue](T) { return constantValue; });
        }
        
        static BoundaryCondition1D Neumann(std::function<T(T)> g) {
            BoundaryCondition1D bc;
            bc.type = BCType::Neumann;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition1D Neumann(T constantFlux) {
            return Neumann([constantFlux](T) { return constantFlux; });
        }
        
        static BoundaryCondition1D Robin(T a, T b, std::function<T(T)> g) {
            BoundaryCondition1D bc;
            bc.type = BCType::Robin;
            bc.alpha = a;
            bc.beta = b;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition1D Periodic() {
            BoundaryCondition1D bc;
            bc.type = BCType::Periodic;
            return bc;
        }
    };

    ///////////////////////////////////////////////////////////////////////
    ///                   2D Boundary Condition                         ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Boundary condition specification for 2D problems
     * 
     * @tparam T Numeric type (float, double)
     * 
     * value function takes (x, y) coordinates on the boundary
     */
    template<typename T>
    struct BoundaryCondition2D {
        BCType type = BCType::Dirichlet;
        std::function<T(T, T)> value;   // g(x, y)
        T alpha = T(1);
        T beta = T(0);
        
        static BoundaryCondition2D Dirichlet(std::function<T(T, T)> g) {
            BoundaryCondition2D bc;
            bc.type = BCType::Dirichlet;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition2D Dirichlet(T constantValue) {
            return Dirichlet([constantValue](T, T) { return constantValue; });
        }
        
        static BoundaryCondition2D Neumann(std::function<T(T, T)> g) {
            BoundaryCondition2D bc;
            bc.type = BCType::Neumann;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition2D Neumann(T constantFlux) {
            return Neumann([constantFlux](T, T) { return constantFlux; });
        }
        
        static BoundaryCondition2D Robin(T a, T b, std::function<T(T, T)> g) {
            BoundaryCondition2D bc;
            bc.type = BCType::Robin;
            bc.alpha = a;
            bc.beta = b;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition2D Periodic() {
            BoundaryCondition2D bc;
            bc.type = BCType::Periodic;
            return bc;
        }
    };

    ///////////////////////////////////////////////////////////////////////
    ///                   3D Boundary Condition                         ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Boundary condition specification for 3D problems
     */
    template<typename T>
    struct BoundaryCondition3D {
        BCType type = BCType::Dirichlet;
        std::function<T(T, T, T)> value;  // g(x, y, z)
        T alpha = T(1);
        T beta = T(0);
        
        static BoundaryCondition3D Dirichlet(std::function<T(T, T, T)> g) {
            BoundaryCondition3D bc;
            bc.type = BCType::Dirichlet;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition3D Dirichlet(T constantValue) {
            return Dirichlet([constantValue](T, T, T) { return constantValue; });
        }
        
        static BoundaryCondition3D Neumann(std::function<T(T, T, T)> g) {
            BoundaryCondition3D bc;
            bc.type = BCType::Neumann;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition3D Neumann(T constantFlux) {
            return Neumann([constantFlux](T, T, T) { return constantFlux; });
        }
        
        static BoundaryCondition3D Robin(T a, T b, std::function<T(T, T, T)> g) {
            BoundaryCondition3D bc;
            bc.type = BCType::Robin;
            bc.alpha = a;
            bc.beta = b;
            bc.value = g;
            return bc;
        }
        
        static BoundaryCondition3D Periodic() {
            BoundaryCondition3D bc;
            bc.type = BCType::Periodic;
            return bc;
        }
    };

    ///////////////////////////////////////////////////////////////////////
    ///               1D Boundary Conditions Container                  ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Container for 1D boundary conditions (left and right)
     * 
     * Manages boundary conditions for both ends of a 1D domain and provides
     * methods to apply them to matrices and grid functions.
     */
    template<typename T>
    class BoundaryConditions1D {
    private:
        BoundaryCondition1D<T> left_;
        BoundaryCondition1D<T> right_;
        
    public:
        BoundaryConditions1D() = default;
        
        BoundaryConditions1D(BoundaryCondition1D<T> left, BoundaryCondition1D<T> right)
            : left_(std::move(left)), right_(std::move(right)) {}
        
        // Setters
        void setLeft(BoundaryCondition1D<T> bc) { left_ = std::move(bc); }
        void setRight(BoundaryCondition1D<T> bc) { right_ = std::move(bc); }
        void set(BoundarySide side, BoundaryCondition1D<T> bc) {
            if (side == BoundarySide::Left) left_ = std::move(bc);
            else if (side == BoundarySide::Right) right_ = std::move(bc);
        }
        
        // Getters
        const BoundaryCondition1D<T>& left() const { return left_; }
        const BoundaryCondition1D<T>& right() const { return right_; }
        const BoundaryCondition1D<T>& get(BoundarySide side) const {
            if (side == BoundarySide::Left) return left_;
            return right_;
        }
        
        /**
         * @brief Apply boundary conditions to sparse matrix and RHS vector
         * 
         * Modifies the matrix A and RHS vector b to incorporate boundary conditions.
         * 
         * For Dirichlet: Row i becomes identity, b[i] = g(x_i)
         * For Neumann: Uses one-sided difference approximation
         * For Robin: Combines value and derivative terms
         * For Periodic: Links first and last nodes
         */
        void apply(SparseMatrixCSR<T>& A, std::vector<T>& b, const Grid1D<T>& grid) const {
            int n = grid.numNodes();
            T dx = grid.dx();
            
            // Left boundary (node 0)
            applyAtNode(A, b, grid, 0, left_, BoundarySide::Left, dx);
            
            // Right boundary (node n-1)
            applyAtNode(A, b, grid, n - 1, right_, BoundarySide::Right, dx);
            
            // Handle periodic BCs
            if (left_.type == BCType::Periodic || right_.type == BCType::Periodic) {
                applyPeriodic(A, b, grid);
            }
        }
        
        /**
         * @brief Apply Dirichlet values directly to a grid function
         * 
         * Sets boundary node values according to Dirichlet conditions.
         * Useful for setting initial conditions or post-processing.
         */
        void applyToFunction(GridFunction1D<T>& u) const {
            const auto& grid = u.grid();
            
            if (left_.type == BCType::Dirichlet && left_.value) {
                T x = grid.x(0);
                u(0) = left_.value(x);
            }
            
            if (right_.type == BCType::Dirichlet && right_.value) {
                T x = grid.x(grid.numNodes() - 1);
                u(grid.numNodes() - 1) = right_.value(x);
            }
        }
        
    private:
        void applyAtNode(SparseMatrixCSR<T>& A, std::vector<T>& b,
                        const Grid1D<T>& grid, int i,
                        const BoundaryCondition1D<T>& bc,
                        BoundarySide side, T dx) const {
            T x = grid.x(i);
            int n = grid.numNodes();
            
            switch (bc.type) {
                case BCType::Dirichlet:
                    // Replace row: u_i = g(x_i)
                    A.zeroRow(i);
                    A.set(i, i, T(1));
                    b[i] = bc.value ? bc.value(x) : T(0);
                    break;
                    
                case BCType::Neumann:
                    // du/dn = g => use one-sided difference
                    // Note: outward normal at left is -x, at right is +x
                    // Left:  du/dn = -du/dx = g => du/dx = -g => (u_1 - u_0)/dx = -g
                    // Right: du/dn = +du/dx = g => du/dx = g  => (u_n - u_{n-1})/dx = g
                    A.zeroRow(i);
                    if (side == BoundarySide::Left) {
                        A.set(i, 0, T(-1));
                        A.set(i, 1, T(1));
                        b[i] = -(bc.value ? bc.value(x) : T(0)) * dx;  // Negative for left
                    } else {
                        A.set(i, n - 2, T(-1));
                        A.set(i, n - 1, T(1));
                        b[i] = (bc.value ? bc.value(x) : T(0)) * dx;   // Positive for right
                    }
                    break;
                    
                case BCType::Robin:
                    // alpha*u + beta*du/dn = g
                    // Note: outward normal at left is -x, at right is +x
                    // Left:  du/dn = -du/dx, so alpha*u - beta*du/dx = g
                    //        alpha*u_0 - beta*(u_1 - u_0)/dx = g
                    //        => (alpha + beta/dx)*u_0 - (beta/dx)*u_1 = g
                    // Right: du/dn = +du/dx, so alpha*u + beta*du/dx = g
                    //        alpha*u_n + beta*(u_n - u_{n-1})/dx = g
                    //        => (-beta/dx)*u_{n-1} + (alpha + beta/dx)*u_n = g
                    A.zeroRow(i);
                    if (side == BoundarySide::Left) {
                        A.set(i, 0, bc.alpha + bc.beta / dx);
                        A.set(i, 1, -bc.beta / dx);
                    } else {
                        A.set(i, n - 2, -bc.beta / dx);
                        A.set(i, n - 1, bc.alpha + bc.beta / dx);
                    }
                    b[i] = bc.value ? bc.value(x) : T(0);
                    break;
                    
                case BCType::Periodic:
                    // Handled separately in applyPeriodic
                    break;
            }
        }
        
        void applyPeriodic(SparseMatrixCSR<T>& A, std::vector<T>& b, const Grid1D<T>& grid) const {
            int n = grid.numNodes();
            
            // u_0 = u_{n-1}: Replace row 0 with u_0 - u_{n-1} = 0
            A.zeroRow(0);
            A.set(0, 0, T(1));
            A.set(0, n - 1, T(-1));
            b[0] = T(0);
            
            // The interior equation at node n-1 wraps around to use node 0
            // This requires modifying the stencil, which is done in the PDE assembler
        }
    };

    ///////////////////////////////////////////////////////////////////////
    ///               2D Boundary Conditions Container                  ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Container for 2D boundary conditions (left, right, bottom, top)
     */
    template<typename T>
    class BoundaryConditions2D {
    private:
        BoundaryCondition2D<T> left_;
        BoundaryCondition2D<T> right_;
        BoundaryCondition2D<T> bottom_;
        BoundaryCondition2D<T> top_;
        
    public:
        BoundaryConditions2D() = default;
        
        // Setters
        void setLeft(BoundaryCondition2D<T> bc) { left_ = std::move(bc); }
        void setRight(BoundaryCondition2D<T> bc) { right_ = std::move(bc); }
        void setBottom(BoundaryCondition2D<T> bc) { bottom_ = std::move(bc); }
        void setTop(BoundaryCondition2D<T> bc) { top_ = std::move(bc); }
        
        void set(BoundarySide side, BoundaryCondition2D<T> bc) {
            switch (side) {
                case BoundarySide::Left:   left_ = std::move(bc); break;
                case BoundarySide::Right:  right_ = std::move(bc); break;
                case BoundarySide::Bottom: bottom_ = std::move(bc); break;
                case BoundarySide::Top:    top_ = std::move(bc); break;
                default: break;
            }
        }
        
        /**
         * @brief Set all boundaries to the same condition
         */
        void setAll(BoundaryCondition2D<T> bc) {
            left_ = bc;
            right_ = bc;
            bottom_ = bc;
            top_ = bc;
        }
        
        // Getters
        const BoundaryCondition2D<T>& left() const { return left_; }
        const BoundaryCondition2D<T>& right() const { return right_; }
        const BoundaryCondition2D<T>& bottom() const { return bottom_; }
        const BoundaryCondition2D<T>& top() const { return top_; }
        
        const BoundaryCondition2D<T>& get(BoundarySide side) const {
            switch (side) {
                case BoundarySide::Left:   return left_;
                case BoundarySide::Right:  return right_;
                case BoundarySide::Bottom: return bottom_;
                case BoundarySide::Top:    return top_;
                default: return left_;
            }
        }
        
        /**
         * @brief Apply boundary conditions to sparse matrix and RHS vector
         */
        void apply(SparseMatrixCSR<T>& A, std::vector<T>& b, const Grid2D<T>& grid) const {
            int nx = grid.nx();
            int ny = grid.ny();
            T dx = grid.dx();
            T dy = grid.dy();
            
            // Apply to all boundary nodes
            grid.forEachBoundary([&](int i, int j, T x, T y, BoundarySide side) {
                int idx = grid.index(i, j);
                applyAtNode(A, b, grid, idx, i, j, x, y, side, dx, dy);
            });
        }
        
        /**
         * @brief Apply Dirichlet values directly to a grid function
         */
        void applyToFunction(GridFunction2D<T>& u) const {
            const auto& grid = u.grid();
            
            grid.forEachBoundary([&](int i, int j, T x, T y, BoundarySide side) {
                const auto& bc = get(side);
                if (bc.type == BCType::Dirichlet && bc.value) {
                    u(i, j) = bc.value(x, y);
                }
            });
        }
        
    private:
        void applyAtNode(SparseMatrixCSR<T>& A, std::vector<T>& b,
                        const Grid2D<T>& grid, int idx, int i, int j,
                        T x, T y, BoundarySide side, T dx, T dy) const {
            const auto& bc = get(side);
            int nx = grid.nx();
            int ny = grid.ny();
            
            switch (bc.type) {
                case BCType::Dirichlet:
                    A.zeroRow(idx);
                    A.set(idx, idx, T(1));
                    b[idx] = bc.value ? bc.value(x, y) : T(0);
                    break;
                    
                case BCType::Neumann: {
                    // One-sided difference in the normal direction
                    // Note: outward normal points in negative direction for left/bottom
                    A.zeroRow(idx);
                    T g = bc.value ? bc.value(x, y) : T(0);
                    
                    switch (side) {
                        case BoundarySide::Left:
                            // Outward normal is -x, so du/dn = -du/dx = g => du/dx = -g
                            // (u_{i+1,j} - u_{i,j})/dx = -g => -u_i + u_{i+1} = -g*dx
                            A.set(idx, idx, T(-1));
                            A.set(idx, grid.index(i + 1, j), T(1));
                            b[idx] = -g * dx;
                            break;
                        case BoundarySide::Right:
                            // Outward normal is +x, so du/dn = +du/dx = g
                            // (u_{i,j} - u_{i-1,j})/dx = g => -u_{i-1} + u_i = g*dx
                            A.set(idx, grid.index(i - 1, j), T(-1));
                            A.set(idx, idx, T(1));
                            b[idx] = g * dx;
                            break;
                        case BoundarySide::Bottom:
                            // Outward normal is -y, so du/dn = -du/dy = g => du/dy = -g
                            // (u_{i,j+1} - u_{i,j})/dy = -g => -u_{i,j} + u_{i,j+1} = -g*dy
                            A.set(idx, idx, T(-1));
                            A.set(idx, grid.index(i, j + 1), T(1));
                            b[idx] = -g * dy;
                            break;
                        case BoundarySide::Top:
                            // Outward normal is +y, so du/dn = +du/dy = g
                            // (u_{i,j} - u_{i,j-1})/dy = g => -u_{i,j-1} + u_{i,j} = g*dy
                            A.set(idx, grid.index(i, j - 1), T(-1));
                            A.set(idx, idx, T(1));
                            b[idx] = g * dy;
                            break;
                        default:
                            break;
                    }
                    break;
                }
                    
                case BCType::Robin: {
                    // alpha*u + beta*du/dn = g
                    // Note: outward normal points in negative direction for left/bottom
                    A.zeroRow(idx);
                    T g = bc.value ? bc.value(x, y) : T(0);
                    
                    switch (side) {
                        case BoundarySide::Left:
                            // du/dn = -du/dx, so alpha*u - beta*du/dx = g
                            // alpha*u_i - beta*(u_{i+1} - u_i)/dx = g
                            // => (alpha + beta/dx)*u_i - (beta/dx)*u_{i+1} = g
                            A.set(idx, idx, bc.alpha + bc.beta / dx);
                            A.set(idx, grid.index(i + 1, j), -bc.beta / dx);
                            break;
                        case BoundarySide::Right:
                            // du/dn = +du/dx, so alpha*u + beta*du/dx = g
                            // alpha*u_i + beta*(u_i - u_{i-1})/dx = g
                            // => (-beta/dx)*u_{i-1} + (alpha + beta/dx)*u_i = g
                            A.set(idx, grid.index(i - 1, j), -bc.beta / dx);
                            A.set(idx, idx, bc.alpha + bc.beta / dx);
                            break;
                        case BoundarySide::Bottom:
                            // du/dn = -du/dy, so alpha*u - beta*du/dy = g
                            // alpha*u_{i,j} - beta*(u_{i,j+1} - u_{i,j})/dy = g
                            // => (alpha + beta/dy)*u_{i,j} - (beta/dy)*u_{i,j+1} = g
                            A.set(idx, idx, bc.alpha + bc.beta / dy);
                            A.set(idx, grid.index(i, j + 1), -bc.beta / dy);
                            break;
                        case BoundarySide::Top:
                            // du/dn = +du/dy, so alpha*u + beta*du/dy = g
                            // alpha*u_{i,j} + beta*(u_{i,j} - u_{i,j-1})/dy = g
                            // => (-beta/dy)*u_{i,j-1} + (alpha + beta/dy)*u_{i,j} = g
                            A.set(idx, grid.index(i, j - 1), -bc.beta / dy);
                            A.set(idx, idx, bc.alpha + bc.beta / dy);
                            break;
                        default:
                            break;
                    }
                    b[idx] = g;
                    break;
                }
                    
                case BCType::Periodic:
                    // Handled separately
                    break;
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////
    ///               3D Boundary Conditions Container                  ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Container for 3D boundary conditions (6 faces)
     */
    template<typename T>
    class BoundaryConditions3D {
    private:
        BoundaryCondition3D<T> left_;    // x = xmin
        BoundaryCondition3D<T> right_;   // x = xmax
        BoundaryCondition3D<T> bottom_;  // y = ymin
        BoundaryCondition3D<T> top_;     // y = ymax
        BoundaryCondition3D<T> front_;   // z = zmin
        BoundaryCondition3D<T> back_;    // z = zmax
        
    public:
        BoundaryConditions3D() = default;
        
        // Setters
        void setLeft(BoundaryCondition3D<T> bc) { left_ = std::move(bc); }
        void setRight(BoundaryCondition3D<T> bc) { right_ = std::move(bc); }
        void setBottom(BoundaryCondition3D<T> bc) { bottom_ = std::move(bc); }
        void setTop(BoundaryCondition3D<T> bc) { top_ = std::move(bc); }
        void setFront(BoundaryCondition3D<T> bc) { front_ = std::move(bc); }
        void setBack(BoundaryCondition3D<T> bc) { back_ = std::move(bc); }
        
        void set(BoundarySide side, BoundaryCondition3D<T> bc) {
            switch (side) {
                case BoundarySide::Left:   left_ = std::move(bc); break;
                case BoundarySide::Right:  right_ = std::move(bc); break;
                case BoundarySide::Bottom: bottom_ = std::move(bc); break;
                case BoundarySide::Top:    top_ = std::move(bc); break;
                case BoundarySide::Front:  front_ = std::move(bc); break;
                case BoundarySide::Back:   back_ = std::move(bc); break;
            }
        }
        
        void setAll(BoundaryCondition3D<T> bc) {
            left_ = bc; right_ = bc;
            bottom_ = bc; top_ = bc;
            front_ = bc; back_ = bc;
        }
        
        // Getters
        const BoundaryCondition3D<T>& left() const { return left_; }
        const BoundaryCondition3D<T>& right() const { return right_; }
        const BoundaryCondition3D<T>& bottom() const { return bottom_; }
        const BoundaryCondition3D<T>& top() const { return top_; }
        const BoundaryCondition3D<T>& front() const { return front_; }
        const BoundaryCondition3D<T>& back() const { return back_; }
        
        const BoundaryCondition3D<T>& get(BoundarySide side) const {
            switch (side) {
                case BoundarySide::Left:   return left_;
                case BoundarySide::Right:  return right_;
                case BoundarySide::Bottom: return bottom_;
                case BoundarySide::Top:    return top_;
                case BoundarySide::Front:  return front_;
                case BoundarySide::Back:   return back_;
                default: return left_;
            }
        }
        
        /**
         * @brief Apply boundary conditions to sparse matrix and RHS vector
         */
        void apply(SparseMatrixCSR<T>& A, std::vector<T>& b, const Grid3D<T>& grid) const {
            T dx = grid.dx();
            T dy = grid.dy();
            T dz = grid.dz();
            
            grid.forEachBoundary([&](int i, int j, int k, T x, T y, T z, BoundarySide side) {
                int idx = grid.index(i, j, k);
                applyAtNode(A, b, grid, idx, i, j, k, x, y, z, side, dx, dy, dz);
            });
        }
        
        /**
         * @brief Apply Dirichlet values directly to a grid function
         */
        void applyToFunction(GridFunction3D<T>& u) const {
            const auto& grid = u.grid();
            
            grid.forEachBoundary([&](int i, int j, int k, T x, T y, T z, BoundarySide side) {
                const auto& bc = get(side);
                if (bc.type == BCType::Dirichlet && bc.value) {
                    u(i, j, k) = bc.value(x, y, z);
                }
            });
        }
        
    private:
        void applyAtNode(SparseMatrixCSR<T>& A, std::vector<T>& b,
                        const Grid3D<T>& grid, int idx, int i, int j, int k,
                        T x, T y, T z, BoundarySide side,
                        T dx, T dy, T dz) const {
            const auto& bc = get(side);
            
            switch (bc.type) {
                case BCType::Dirichlet:
                    A.zeroRow(idx);
                    A.set(idx, idx, T(1));
                    b[idx] = bc.value ? bc.value(x, y, z) : T(0);
                    break;
                    
                case BCType::Neumann: {
                    // Neumann BC: ∂u/∂n = g, where n is the OUTWARD normal
                    // The outward normal points in opposite directions at opposite boundaries:
                    // - Left/Bottom/Front: outward normal points in NEGATIVE direction
                    // - Right/Top/Back: outward normal points in POSITIVE direction
                    A.zeroRow(idx);
                    T g = bc.value ? bc.value(x, y, z) : T(0);
                    
                    switch (side) {
                        case BoundarySide::Left:
                            // Outward normal n = (-1,0,0), so ∂u/∂n = -∂u/∂x = g
                            // Discretize: -(u[i+1] - u[i])/dx = g → u[i] - u[i+1] = g*dx
                            // Actually: (u[i+1] - u[i])/dx = -g → -u[i] + u[i+1] = -g*dx
                            A.set(idx, idx, T(-1));
                            A.set(idx, grid.index(i + 1, j, k), T(1));
                            b[idx] = -g * dx;  // Note: NEGATIVE for left boundary
                            break;
                        case BoundarySide::Right:
                            // Outward normal n = (+1,0,0), so ∂u/∂n = +∂u/∂x = g
                            // Discretize: (u[i] - u[i-1])/dx = g → -u[i-1] + u[i] = g*dx
                            A.set(idx, grid.index(i - 1, j, k), T(-1));
                            A.set(idx, idx, T(1));
                            b[idx] = g * dx;   // POSITIVE for right boundary
                            break;
                        case BoundarySide::Bottom:
                            // Outward normal n = (0,-1,0), so ∂u/∂n = -∂u/∂y = g
                            A.set(idx, idx, T(-1));
                            A.set(idx, grid.index(i, j + 1, k), T(1));
                            b[idx] = -g * dy;  // Note: NEGATIVE for bottom boundary
                            break;
                        case BoundarySide::Top:
                            // Outward normal n = (0,+1,0), so ∂u/∂n = +∂u/∂y = g
                            A.set(idx, grid.index(i, j - 1, k), T(-1));
                            A.set(idx, idx, T(1));
                            b[idx] = g * dy;   // POSITIVE for top boundary
                            break;
                        case BoundarySide::Front:
                            // Outward normal n = (0,0,-1), so ∂u/∂n = -∂u/∂z = g
                            A.set(idx, idx, T(-1));
                            A.set(idx, grid.index(i, j, k + 1), T(1));
                            b[idx] = -g * dz;  // Note: NEGATIVE for front boundary
                            break;
                        case BoundarySide::Back:
                            // Outward normal n = (0,0,+1), so ∂u/∂n = +∂u/∂z = g
                            A.set(idx, grid.index(i, j, k - 1), T(-1));
                            A.set(idx, idx, T(1));
                            b[idx] = g * dz;   // POSITIVE for back boundary
                            break;
                    }
                    break;
                }
                    
                case BCType::Robin: {
                    // Robin BC: α*u + β*(∂u/∂n) = g
                    // The outward normal direction matters for ∂u/∂n:
                    // - Left/Bottom/Front: ∂u/∂n = -∂u/∂x (or -∂u/∂y, -∂u/∂z)
                    // - Right/Top/Back: ∂u/∂n = +∂u/∂x (or +∂u/∂y, +∂u/∂z)
                    A.zeroRow(idx);
                    T g = bc.value ? bc.value(x, y, z) : T(0);
                    
                    switch (side) {
                        case BoundarySide::Left:
                            // α*u + β*(-∂u/∂x) = g → α*u - β*(u[i+1]-u[i])/dx = g
                            // → (α + β/dx)*u[i] - (β/dx)*u[i+1] = g
                            A.set(idx, idx, bc.alpha + bc.beta / dx);
                            A.set(idx, grid.index(i + 1, j, k), -bc.beta / dx);
                            break;
                        case BoundarySide::Right:
                            // α*u + β*(+∂u/∂x) = g → α*u + β*(u[i]-u[i-1])/dx = g
                            // → -(β/dx)*u[i-1] + (α + β/dx)*u[i] = g
                            A.set(idx, grid.index(i - 1, j, k), -bc.beta / dx);
                            A.set(idx, idx, bc.alpha + bc.beta / dx);
                            break;
                        case BoundarySide::Bottom:
                            // α*u + β*(-∂u/∂y) = g
                            // → (α + β/dy)*u[j] - (β/dy)*u[j+1] = g
                            A.set(idx, idx, bc.alpha + bc.beta / dy);
                            A.set(idx, grid.index(i, j + 1, k), -bc.beta / dy);
                            break;
                        case BoundarySide::Top:
                            // α*u + β*(+∂u/∂y) = g
                            // → -(β/dy)*u[j-1] + (α + β/dy)*u[j] = g
                            A.set(idx, grid.index(i, j - 1, k), -bc.beta / dy);
                            A.set(idx, idx, bc.alpha + bc.beta / dy);
                            break;
                        case BoundarySide::Front:
                            // α*u + β*(-∂u/∂z) = g
                            // → (α + β/dz)*u[k] - (β/dz)*u[k+1] = g
                            A.set(idx, idx, bc.alpha + bc.beta / dz);
                            A.set(idx, grid.index(i, j, k + 1), -bc.beta / dz);
                            break;
                        case BoundarySide::Back:
                            // α*u + β*(+∂u/∂z) = g
                            // → -(β/dz)*u[k-1] + (α + β/dz)*u[k] = g
                            A.set(idx, grid.index(i, j, k - 1), -bc.beta / dz);
                            A.set(idx, idx, bc.alpha + bc.beta / dz);
                            break;
                    }
                    b[idx] = g;
                    break;
                }
                    
                case BCType::Periodic:
                    break;
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////
    ///                      Convenience Functions                      ///
    ///////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Create homogeneous Dirichlet BCs (u = 0 on all boundaries)
     */
    template<typename T>
    BoundaryConditions1D<T> homogeneousDirichlet1D() {
        return BoundaryConditions1D<T>(
            BoundaryCondition1D<T>::Dirichlet(T(0)),
            BoundaryCondition1D<T>::Dirichlet(T(0))
        );
    }
    
    template<typename T>
    BoundaryConditions2D<T> homogeneousDirichlet2D() {
        BoundaryConditions2D<T> bc;
        bc.setAll(BoundaryCondition2D<T>::Dirichlet(T(0)));
        return bc;
    }
    
    template<typename T>
    BoundaryConditions3D<T> homogeneousDirichlet3D() {
        BoundaryConditions3D<T> bc;
        bc.setAll(BoundaryCondition3D<T>::Dirichlet(T(0)));
        return bc;
    }
    
    /**
     * @brief Create homogeneous Neumann BCs (du/dn = 0 on all boundaries)
     */
    template<typename T>
    BoundaryConditions1D<T> homogeneousNeumann1D() {
        return BoundaryConditions1D<T>(
            BoundaryCondition1D<T>::Neumann(T(0)),
            BoundaryCondition1D<T>::Neumann(T(0))
        );
    }
    
    template<typename T>
    BoundaryConditions2D<T> homogeneousNeumann2D() {
        BoundaryConditions2D<T> bc;
        bc.setAll(BoundaryCondition2D<T>::Neumann(T(0)));
        return bc;
    }
    
    /**
     * @brief Create periodic boundary conditions
     */
    template<typename T>
    BoundaryConditions1D<T> periodic1D() {
        return BoundaryConditions1D<T>(
            BoundaryCondition1D<T>::Periodic(),
            BoundaryCondition1D<T>::Periodic()
        );
    }

} // namespace MML::PDE

#endif // MML_BOUNDARY_CONDITIONS_H
