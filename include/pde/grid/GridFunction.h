///////////////////////////////////////////////////////////////////////////////////////////
// GridFunction.h - Grid functions for storing field values on structured grids
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// GridFunction stores values at grid nodes and provides convenient access
// patterns for finite difference computations.
//
// Features:
// - Efficient contiguous storage
// - Multi-index and linear index access
// - Initialization from functions
// - Norm computations
// - Import/export to vectors
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_GRID_FUNCTION_H
#define MML_PDE_GRID_FUNCTION_H

#include "Grid1D.h"
#include "Grid2D.h"
#include "Grid3D.h"
#include <functional>
#include <cmath>
#include <numeric>
#include <algorithm>

namespace MML::PDE {

//=============================================================================
// GridFunction1D
//=============================================================================

template<typename T>
class GridFunction1D {
    Grid1D<T> grid_;
    std::vector<T> values_;
    
public:
    //=========================================================================
    // Construction
    //=========================================================================
    
    explicit GridFunction1D(const Grid1D<T>& grid)
        : grid_(grid), values_(grid.numNodes(), T(0)) {}
    
    GridFunction1D(const Grid1D<T>& grid, T initialValue)
        : grid_(grid), values_(grid.numNodes(), initialValue) {}
    
    GridFunction1D(const Grid1D<T>& grid, const std::vector<T>& values)
        : grid_(grid), values_(values) 
    {
        if (values_.size() != static_cast<size_t>(grid_.numNodes())) {
            throw std::invalid_argument("GridFunction1D: values size mismatch");
        }
    }
    
    // Initialize from function f(x)
    GridFunction1D(const Grid1D<T>& grid, std::function<T(T)> f)
        : grid_(grid), values_(grid.numNodes())
    {
        for (int i = 0; i <= grid_.numCells(); ++i) {
            values_[i] = f(grid_.x(i));
        }
    }
    
    //=========================================================================
    // Grid Access
    //=========================================================================
    
    const Grid1D<T>& grid() const { return grid_; }
    int size() const { return static_cast<int>(values_.size()); }
    
    //=========================================================================
    // Value Access
    //=========================================================================
    
    T& operator()(int i) { return values_[i]; }
    const T& operator()(int i) const { return values_[i]; }
    
    T& operator[](int i) { return values_[i]; }
    const T& operator[](int i) const { return values_[i]; }
    
    // Access underlying vector
    std::vector<T>& values() { return values_; }
    const std::vector<T>& values() const { return values_; }
    T* data() { return values_.data(); }
    const T* data() const { return values_.data(); }
    
    //=========================================================================
    // Initialization
    //=========================================================================
    
    void fill(T value) {
        std::fill(values_.begin(), values_.end(), value);
    }
    
    void setFromFunction(std::function<T(T)> f) {
        for (int i = 0; i <= grid_.numCells(); ++i) {
            values_[i] = f(grid_.x(i));
        }
    }
    
    void setInterior(T value) {
        for (int i = 1; i < grid_.numCells(); ++i) {
            values_[i] = value;
        }
    }
    
    void setBoundary(T left, T right) {
        values_[0] = left;
        values_[grid_.numCells()] = right;
    }
    
    //=========================================================================
    // Norms
    //=========================================================================
    
    // L∞ norm (max absolute value)
    T normInf() const {
        T maxVal = T(0);
        for (const T& v : values_) {
            maxVal = std::max(maxVal, std::abs(v));
        }
        return maxVal;
    }
    
    // L2 norm (sqrt of sum of squares)
    T normL2() const {
        T sum = T(0);
        for (const T& v : values_) {
            sum += v * v;
        }
        return std::sqrt(sum * grid_.dx());
    }
    
    // L1 norm (integral of absolute value)
    T normL1() const {
        T sum = T(0);
        for (const T& v : values_) {
            sum += std::abs(v);
        }
        return sum * grid_.dx();
    }
    
    //=========================================================================
    // Arithmetic Operations
    //=========================================================================
    
    GridFunction1D& operator+=(const GridFunction1D& other) {
        for (size_t i = 0; i < values_.size(); ++i) {
            values_[i] += other.values_[i];
        }
        return *this;
    }
    
    GridFunction1D& operator-=(const GridFunction1D& other) {
        for (size_t i = 0; i < values_.size(); ++i) {
            values_[i] -= other.values_[i];
        }
        return *this;
    }
    
    GridFunction1D& operator*=(T scalar) {
        for (T& v : values_) v *= scalar;
        return *this;
    }
    
    GridFunction1D& operator/=(T scalar) {
        for (T& v : values_) v /= scalar;
        return *this;
    }
    
    //=========================================================================
    // Interior Values (for linear solver interface)
    //=========================================================================
    
    // Extract interior values to a vector
    std::vector<T> interiorValues() const {
        std::vector<T> interior(grid_.numInteriorNodes());
        for (int i = 1; i < grid_.numCells(); ++i) {
            interior[i - 1] = values_[i];
        }
        return interior;
    }
    
    // Set interior values from a vector
    void setInteriorValues(const std::vector<T>& interior) {
        if (interior.size() != static_cast<size_t>(grid_.numInteriorNodes())) {
            throw std::invalid_argument("Interior values size mismatch");
        }
        for (int i = 1; i < grid_.numCells(); ++i) {
            values_[i] = interior[i - 1];
        }
    }
    
    //=========================================================================
    // Interpolation
    //=========================================================================
    
    // Linear interpolation at point x
    T interpolate(T xp) const {
        int i;
        T alpha;
        grid_.locate(xp, i, alpha);
        return (T(1) - alpha) * values_[i] + alpha * values_[i + 1];
    }
};

//=============================================================================
// GridFunction2D
//=============================================================================

template<typename T>
class GridFunction2D {
    Grid2D<T> grid_;
    std::vector<T> values_;
    
public:
    //=========================================================================
    // Construction
    //=========================================================================
    
    explicit GridFunction2D(const Grid2D<T>& grid)
        : grid_(grid), values_(grid.numNodes(), T(0)) {}
    
    GridFunction2D(const Grid2D<T>& grid, T initialValue)
        : grid_(grid), values_(grid.numNodes(), initialValue) {}
    
    GridFunction2D(const Grid2D<T>& grid, const std::vector<T>& values)
        : grid_(grid), values_(values)
    {
        if (values_.size() != static_cast<size_t>(grid_.numNodes())) {
            throw std::invalid_argument("GridFunction2D: values size mismatch");
        }
    }
    
    // Initialize from function f(x, y)
    GridFunction2D(const Grid2D<T>& grid, std::function<T(T, T)> f)
        : grid_(grid), values_(grid.numNodes())
    {
        for (int j = 0; j <= grid_.ny(); ++j) {
            for (int i = 0; i <= grid_.nx(); ++i) {
                values_[grid_.index(i, j)] = f(grid_.x(i), grid_.y(j));
            }
        }
    }
    
    //=========================================================================
    // Grid Access
    //=========================================================================
    
    const Grid2D<T>& grid() const { return grid_; }
    int size() const { return static_cast<int>(values_.size()); }
    
    //=========================================================================
    // Value Access
    //=========================================================================
    
    // Multi-index access
    T& operator()(int i, int j) { return values_[grid_.index(i, j)]; }
    const T& operator()(int i, int j) const { return values_[grid_.index(i, j)]; }
    
    // Linear index access
    T& operator[](int idx) { return values_[idx]; }
    const T& operator[](int idx) const { return values_[idx]; }
    
    // Access underlying vector
    std::vector<T>& values() { return values_; }
    const std::vector<T>& values() const { return values_; }
    T* data() { return values_.data(); }
    const T* data() const { return values_.data(); }
    
    //=========================================================================
    // Initialization
    //=========================================================================
    
    void fill(T value) {
        std::fill(values_.begin(), values_.end(), value);
    }
    
    void setFromFunction(std::function<T(T, T)> f) {
        for (int j = 0; j <= grid_.ny(); ++j) {
            for (int i = 0; i <= grid_.nx(); ++i) {
                values_[grid_.index(i, j)] = f(grid_.x(i), grid_.y(j));
            }
        }
    }
    
    void setInterior(T value) {
        for (int j = 1; j < grid_.ny(); ++j) {
            for (int i = 1; i < grid_.nx(); ++i) {
                values_[grid_.index(i, j)] = value;
            }
        }
    }
    
    void setBoundary(T value) {
        // Bottom and top
        for (int i = 0; i <= grid_.nx(); ++i) {
            values_[grid_.index(i, 0)] = value;
            values_[grid_.index(i, grid_.ny())] = value;
        }
        // Left and right
        for (int j = 1; j < grid_.ny(); ++j) {
            values_[grid_.index(0, j)] = value;
            values_[grid_.index(grid_.nx(), j)] = value;
        }
    }
    
    void setBoundary(BoundarySide side, T value) {
        auto nodes = grid_.boundaryNodes(side);
        for (int idx : nodes) {
            values_[idx] = value;
        }
    }
    
    void setBoundary(BoundarySide side, std::function<T(T, T)> f) {
        auto nodes = grid_.boundaryNodes(side);
        for (int idx : nodes) {
            T xc, yc;
            grid_.coords(idx, xc, yc);
            values_[idx] = f(xc, yc);
        }
    }
    
    //=========================================================================
    // Norms
    //=========================================================================
    
    T normInf() const {
        T maxVal = T(0);
        for (const T& v : values_) {
            maxVal = std::max(maxVal, std::abs(v));
        }
        return maxVal;
    }
    
    T normL1() const {
        T sum = T(0);
        for (const T& v : values_) {
            sum += std::abs(v);
        }
        return sum;
    }
    
    T normL2() const {
        T sum = T(0);
        for (const T& v : values_) {
            sum += v * v;
        }
        return std::sqrt(sum * grid_.dx() * grid_.dy());
    }
    
    //=========================================================================
    // Arithmetic Operations
    //=========================================================================
    
    GridFunction2D& operator+=(const GridFunction2D& other) {
        for (size_t i = 0; i < values_.size(); ++i) {
            values_[i] += other.values_[i];
        }
        return *this;
    }
    
    GridFunction2D& operator-=(const GridFunction2D& other) {
        for (size_t i = 0; i < values_.size(); ++i) {
            values_[i] -= other.values_[i];
        }
        return *this;
    }
    
    GridFunction2D& operator*=(T scalar) {
        for (T& v : values_) v *= scalar;
        return *this;
    }
    
    //=========================================================================
    // Interior Values
    //=========================================================================
    
    std::vector<T> interiorValues() const {
        auto interiorNodes = grid_.interiorNodes();
        std::vector<T> interior(interiorNodes.size());
        for (size_t k = 0; k < interiorNodes.size(); ++k) {
            interior[k] = values_[interiorNodes[k]];
        }
        return interior;
    }
    
    void setInteriorValues(const std::vector<T>& interior) {
        auto interiorNodes = grid_.interiorNodes();
        if (interior.size() != interiorNodes.size()) {
            throw std::invalid_argument("Interior values size mismatch");
        }
        for (size_t k = 0; k < interiorNodes.size(); ++k) {
            values_[interiorNodes[k]] = interior[k];
        }
    }
    
    //=========================================================================
    // Interpolation
    //=========================================================================
    
    // Bilinear interpolation at point (xp, yp)
    T interpolate(T xp, T yp) const {
        int i, j;
        T alpha, beta;
        grid_.locate(xp, yp, i, j, alpha, beta);
        
        T v00 = values_[grid_.index(i, j)];
        T v10 = values_[grid_.index(i + 1, j)];
        T v01 = values_[grid_.index(i, j + 1)];
        T v11 = values_[grid_.index(i + 1, j + 1)];
        
        return (T(1) - alpha) * (T(1) - beta) * v00 +
               alpha * (T(1) - beta) * v10 +
               (T(1) - alpha) * beta * v01 +
               alpha * beta * v11;
    }
};

//=============================================================================
// GridFunction3D
//=============================================================================

template<typename T>
class GridFunction3D {
    Grid3D<T> grid_;
    std::vector<T> values_;
    
public:
    //=========================================================================
    // Construction
    //=========================================================================
    
    explicit GridFunction3D(const Grid3D<T>& grid)
        : grid_(grid), values_(grid.numNodes(), T(0)) {}
    
    GridFunction3D(const Grid3D<T>& grid, T initialValue)
        : grid_(grid), values_(grid.numNodes(), initialValue) {}
    
    GridFunction3D(const Grid3D<T>& grid, const std::vector<T>& values)
        : grid_(grid), values_(values)
    {
        if (values_.size() != static_cast<size_t>(grid_.numNodes())) {
            throw std::invalid_argument("GridFunction3D: values size mismatch");
        }
    }
    
    // Initialize from function f(x, y, z)
    GridFunction3D(const Grid3D<T>& grid, std::function<T(T, T, T)> f)
        : grid_(grid), values_(grid.numNodes())
    {
        for (int k = 0; k <= grid_.nz(); ++k) {
            for (int j = 0; j <= grid_.ny(); ++j) {
                for (int i = 0; i <= grid_.nx(); ++i) {
                    values_[grid_.index(i, j, k)] = f(grid_.x(i), grid_.y(j), grid_.z(k));
                }
            }
        }
    }
    
    //=========================================================================
    // Grid Access
    //=========================================================================
    
    const Grid3D<T>& grid() const { return grid_; }
    int size() const { return static_cast<int>(values_.size()); }
    
    //=========================================================================
    // Value Access
    //=========================================================================
    
    T& operator()(int i, int j, int k) { return values_[grid_.index(i, j, k)]; }
    const T& operator()(int i, int j, int k) const { return values_[grid_.index(i, j, k)]; }
    
    T& operator[](int idx) { return values_[idx]; }
    const T& operator[](int idx) const { return values_[idx]; }
    
    std::vector<T>& values() { return values_; }
    const std::vector<T>& values() const { return values_; }
    T* data() { return values_.data(); }
    const T* data() const { return values_.data(); }
    
    //=========================================================================
    // Initialization
    //=========================================================================
    
    void fill(T value) {
        std::fill(values_.begin(), values_.end(), value);
    }
    
    void setFromFunction(std::function<T(T, T, T)> f) {
        for (int k = 0; k <= grid_.nz(); ++k) {
            for (int j = 0; j <= grid_.ny(); ++j) {
                for (int i = 0; i <= grid_.nx(); ++i) {
                    values_[grid_.index(i, j, k)] = f(grid_.x(i), grid_.y(j), grid_.z(k));
                }
            }
        }
    }
    
    void setInterior(T value) {
        for (int k = 1; k < grid_.nz(); ++k) {
            for (int j = 1; j < grid_.ny(); ++j) {
                for (int i = 1; i < grid_.nx(); ++i) {
                    values_[grid_.index(i, j, k)] = value;
                }
            }
        }
    }
    
    //=========================================================================
    // Norms
    //=========================================================================
    
    T normInf() const {
        T maxVal = T(0);
        for (const T& v : values_) {
            maxVal = std::max(maxVal, std::abs(v));
        }
        return maxVal;
    }
    
    T normL2() const {
        T sum = T(0);
        for (const T& v : values_) {
            sum += v * v;
        }
        return std::sqrt(sum * grid_.dx() * grid_.dy() * grid_.dz());
    }
    
    //=========================================================================
    // Arithmetic Operations
    //=========================================================================
    
    GridFunction3D& operator+=(const GridFunction3D& other) {
        for (size_t i = 0; i < values_.size(); ++i) {
            values_[i] += other.values_[i];
        }
        return *this;
    }
    
    GridFunction3D& operator-=(const GridFunction3D& other) {
        for (size_t i = 0; i < values_.size(); ++i) {
            values_[i] -= other.values_[i];
        }
        return *this;
    }
    
    GridFunction3D& operator*=(T scalar) {
        for (T& v : values_) v *= scalar;
        return *this;
    }
    
    //=========================================================================
    // Interior Values
    //=========================================================================
    
    std::vector<T> interiorValues() const {
        auto interiorNodes = grid_.interiorNodes();
        std::vector<T> interior(interiorNodes.size());
        for (size_t k = 0; k < interiorNodes.size(); ++k) {
            interior[k] = values_[interiorNodes[k]];
        }
        return interior;
    }
    
    void setInteriorValues(const std::vector<T>& interior) {
        auto interiorNodes = grid_.interiorNodes();
        if (interior.size() != interiorNodes.size()) {
            throw std::invalid_argument("Interior values size mismatch");
        }
        for (size_t k = 0; k < interiorNodes.size(); ++k) {
            values_[interiorNodes[k]] = interior[k];
        }
    }
};

//=============================================================================
// Free Functions
//=============================================================================

// Difference between two grid functions
template<typename T>
GridFunction1D<T> operator-(const GridFunction1D<T>& a, const GridFunction1D<T>& b) {
    GridFunction1D<T> result = a;
    result -= b;
    return result;
}

template<typename T>
GridFunction2D<T> operator-(const GridFunction2D<T>& a, const GridFunction2D<T>& b) {
    GridFunction2D<T> result = a;
    result -= b;
    return result;
}

template<typename T>
GridFunction3D<T> operator-(const GridFunction3D<T>& a, const GridFunction3D<T>& b) {
    GridFunction3D<T> result = a;
    result -= b;
    return result;
}

} // namespace MML::PDE

#endif // MML_PDE_GRID_FUNCTION_H
