///////////////////////////////////////////////////////////////////////////////////////////
// Grid.h - Base grid infrastructure for structured grids
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// This file provides the foundational infrastructure for structured grids
// used in finite difference PDE discretizations.
//
// Features:
// - Domain specification (bounding box)
// - Common grid properties and interfaces
// - Index conversion utilities
// - Boundary identification
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_GRID_H
#define MML_PDE_GRID_H

#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace MML::PDE {

//=============================================================================
// Domain Specification
//=============================================================================

// 1D interval [xmin, xmax]
template<typename T>
struct Interval {
    T xmin = T(0);
    T xmax = T(1);
    
    Interval() = default;
    Interval(T min, T max) : xmin(min), xmax(max) {}
    
    T min() const { return xmin; }
    T max() const { return xmax; }
    T length() const { return xmax - xmin; }
    T center() const { return (xmin + xmax) / T(2); }
    bool contains(T x) const { return x >= xmin && x <= xmax; }
};

// 2D rectangular domain [xmin, xmax] × [ymin, ymax]
template<typename T>
struct Rectangle {
    T xmin_ = T(0), xmax_ = T(1);
    T ymin_ = T(0), ymax_ = T(1);
    
    Rectangle() = default;
    Rectangle(T x0, T x1, T y0, T y1) : xmin_(x0), xmax_(x1), ymin_(y0), ymax_(y1) {}
    Rectangle(const Interval<T>& ix, const Interval<T>& iy) 
        : xmin_(ix.xmin), xmax_(ix.xmax), ymin_(iy.xmin), ymax_(iy.xmax) {}
    
    T xMin() const { return xmin_; }
    T xMax() const { return xmax_; }
    T yMin() const { return ymin_; }
    T yMax() const { return ymax_; }
    T width() const { return xmax_ - xmin_; }
    T height() const { return ymax_ - ymin_; }
    T area() const { return width() * height(); }
    bool contains(T x, T y) const { 
        return x >= xmin_ && x <= xmax_ && y >= ymin_ && y <= ymax_; 
    }
};

// 3D box domain [xmin, xmax] × [ymin, ymax] × [zmin, zmax]
template<typename T>
struct Box {
    T xmin_ = T(0), xmax_ = T(1);
    T ymin_ = T(0), ymax_ = T(1);
    T zmin_ = T(0), zmax_ = T(1);
    
    Box() = default;
    Box(T x0, T x1, T y0, T y1, T z0, T z1) 
        : xmin_(x0), xmax_(x1), ymin_(y0), ymax_(y1), zmin_(z0), zmax_(z1) {}
    Box(const Interval<T>& ix, const Interval<T>& iy, const Interval<T>& iz)
        : xmin_(ix.xmin), xmax_(ix.xmax), ymin_(iy.xmin), ymax_(iy.xmax),
          zmin_(iz.xmin), zmax_(iz.xmax) {}
    
    T xMin() const { return xmin_; }
    T xMax() const { return xmax_; }
    T yMin() const { return ymin_; }
    T yMax() const { return ymax_; }
    T zMin() const { return zmin_; }
    T zMax() const { return zmax_; }
    T xLength() const { return xmax_ - xmin_; }
    T yLength() const { return ymax_ - ymin_; }
    T zLength() const { return zmax_ - zmin_; }
    T volume() const { return xLength() * yLength() * zLength(); }
    bool contains(T x, T y, T z) const {
        return x >= xmin_ && x <= xmax_ && 
               y >= ymin_ && y <= ymax_ && 
               z >= zmin_ && z <= zmax_;
    }
};

//=============================================================================
// Boundary Types
//=============================================================================

enum class BoundarySide {
    // 1D
    Left, Right,
    // 2D (adds)
    Bottom, Top,
    // 3D (adds)
    Front, Back
};

// Node location classification
enum class NodeLocation {
    Interior,      // Not on any boundary
    Boundary,      // On domain boundary
    Corner         // On multiple boundaries (corner/edge)
};

//=============================================================================
// Index Utilities
//=============================================================================

// Convert multi-index to linear index (row-major order)
int toLinear(int i, int nx);
int toLinear(int i, int j, int nx, int ny);
int toLinear(int i, int j, int k, int nx, int ny, int nz);

// Convert linear index to multi-index
void fromLinear(int idx, int nx, int& i);
void fromLinear(int idx, int nx, int ny, int& i, int& j);
void fromLinear(int idx, int nx, int ny, int nz, int& i, int& j, int& k);

//=============================================================================
// Grid Base Class
//=============================================================================

template<typename T, int Dim>
class GridBase {
public:
    virtual ~GridBase() = default;
    
    // Dimension
    static constexpr int dimension() { return Dim; }
    
    // Total number of nodes
    virtual int totalNodes() const = 0;
    
    // Check if index is valid
    virtual bool isValid(int linearIdx) const {
        return linearIdx >= 0 && linearIdx < totalNodes();
    }
    
    // Node location (interior, boundary, corner)
    virtual NodeLocation nodeLocation(int linearIdx) const = 0;
    
    // Is node on boundary?
    bool isBoundary(int linearIdx) const {
        return nodeLocation(linearIdx) != NodeLocation::Interior;
    }
    
    // Is node in interior?
    bool isInterior(int linearIdx) const {
        return nodeLocation(linearIdx) == NodeLocation::Interior;
    }
};

} // namespace MML::PDE

#endif // MML_PDE_GRID_H
