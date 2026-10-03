///////////////////////////////////////////////////////////////////////////////////////////
// Stencils.h - Finite difference stencil operations on grid functions
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Provides finite difference operators for computing derivatives and
// Laplacians on structured grids.
//
// Stencil types:
// - First derivatives: forward, backward, central
// - Second derivatives: standard 3-point
// - Laplacian: 3-point (1D), 5-point (2D), 7-point (3D)
// - Higher-order stencils (4th order accurate)
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_STENCILS_H
#define MML_PDE_STENCILS_H

#include "GridFunction.h"

namespace MML::PDE {

//=============================================================================
// 1D Stencil Operations
//=============================================================================

// First derivative - forward difference: (u[i+1] - u[i]) / dx
// Order: O(dx)
template<typename T>
T D1Forward(const GridFunction1D<T>& u, int i) {
    return (u(i + 1) - u(i)) / u.grid().dx();
}

// First derivative - backward difference: (u[i] - u[i-1]) / dx
// Order: O(dx)
template<typename T>
T D1Backward(const GridFunction1D<T>& u, int i) {
    return (u(i) - u(i - 1)) / u.grid().dx();
}

// First derivative - central difference: (u[i+1] - u[i-1]) / (2*dx)
// Order: O(dx²)
template<typename T>
T D1Central(const GridFunction1D<T>& u, int i) {
    return (u(i + 1) - u(i - 1)) / (T(2) * u.grid().dx());
}

// Second derivative: (u[i+1] - 2*u[i] + u[i-1]) / dx²
// Order: O(dx²)
template<typename T>
T D2(const GridFunction1D<T>& u, int i) {
    T dx = u.grid().dx();
    return (u(i + 1) - T(2) * u(i) + u(i - 1)) / (dx * dx);
}

// 1D Laplacian (same as D2)
template<typename T>
T Laplacian1D(const GridFunction1D<T>& u, int i) {
    return D2(u, i);
}

// 4th-order accurate first derivative
// (-u[i+2] + 8*u[i+1] - 8*u[i-1] + u[i-2]) / (12*dx)
template<typename T>
T D1Central4(const GridFunction1D<T>& u, int i) {
    T dx = u.grid().dx();
    return (-u(i + 2) + T(8) * u(i + 1) - T(8) * u(i - 1) + u(i - 2)) / (T(12) * dx);
}

// 4th-order accurate second derivative
// (-u[i+2] + 16*u[i+1] - 30*u[i] + 16*u[i-1] - u[i-2]) / (12*dx²)
template<typename T>
T D2_4(const GridFunction1D<T>& u, int i) {
    T dx = u.grid().dx();
    return (-u(i + 2) + T(16) * u(i + 1) - T(30) * u(i) + T(16) * u(i - 1) - u(i - 2)) 
           / (T(12) * dx * dx);
}

//=============================================================================
// 2D Stencil Operations
//=============================================================================

// First derivative in x - central: (u[i+1,j] - u[i-1,j]) / (2*dx)
template<typename T>
T DxCentral(const GridFunction2D<T>& u, int i, int j) {
    return (u(i + 1, j) - u(i - 1, j)) / (T(2) * u.grid().dx());
}

// First derivative in y - central: (u[i,j+1] - u[i,j-1]) / (2*dy)
template<typename T>
T DyCentral(const GridFunction2D<T>& u, int i, int j) {
    return (u(i, j + 1) - u(i, j - 1)) / (T(2) * u.grid().dy());
}

// Second derivative in x: (u[i+1,j] - 2*u[i,j] + u[i-1,j]) / dx²
template<typename T>
T Dxx(const GridFunction2D<T>& u, int i, int j) {
    T dx = u.grid().dx();
    return (u(i + 1, j) - T(2) * u(i, j) + u(i - 1, j)) / (dx * dx);
}

// Second derivative in y: (u[i,j+1] - 2*u[i,j] + u[i,j-1]) / dy²
template<typename T>
T Dyy(const GridFunction2D<T>& u, int i, int j) {
    T dy = u.grid().dy();
    return (u(i, j + 1) - T(2) * u(i, j) + u(i, j - 1)) / (dy * dy);
}

// Mixed derivative: (u[i+1,j+1] - u[i+1,j-1] - u[i-1,j+1] + u[i-1,j-1]) / (4*dx*dy)
template<typename T>
T Dxy(const GridFunction2D<T>& u, int i, int j) {
    T dx = u.grid().dx();
    T dy = u.grid().dy();
    return (u(i + 1, j + 1) - u(i + 1, j - 1) - u(i - 1, j + 1) + u(i - 1, j - 1))
           / (T(4) * dx * dy);
}

// 2D Laplacian (5-point stencil): ∇²u = Dxx + Dyy
template<typename T>
T Laplacian2D(const GridFunction2D<T>& u, int i, int j) {
    return Dxx(u, i, j) + Dyy(u, i, j);
}

// Gradient magnitude: |∇u| = sqrt((∂u/∂x)² + (∂u/∂y)²)
template<typename T>
T GradientMagnitude2D(const GridFunction2D<T>& u, int i, int j) {
    T ux = DxCentral(u, i, j);
    T uy = DyCentral(u, i, j);
    return std::sqrt(ux * ux + uy * uy);
}

// 9-point Laplacian (more isotropic, 4th order on uniform grids)
// Uses [1 4 1; 4 -20 4; 1 4 1] / (6*h²) stencil
template<typename T>
T Laplacian2D_9pt(const GridFunction2D<T>& u, int i, int j) {
    if (u.grid().dx() != u.grid().dy()) {
        // Fall back to 5-point for non-square grids
        return Laplacian2D(u, i, j);
    }
    T h = u.grid().dx();
    T h2 = h * h;
    
    return (u(i - 1, j - 1) + u(i + 1, j - 1) + u(i - 1, j + 1) + u(i + 1, j + 1)
            + T(4) * (u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1))
            - T(20) * u(i, j)) / (T(6) * h2);
}

//=============================================================================
// 3D Stencil Operations
//=============================================================================

// First derivatives - central
template<typename T>
T DxCentral(const GridFunction3D<T>& u, int i, int j, int k) {
    return (u(i + 1, j, k) - u(i - 1, j, k)) / (T(2) * u.grid().dx());
}

template<typename T>
T DyCentral(const GridFunction3D<T>& u, int i, int j, int k) {
    return (u(i, j + 1, k) - u(i, j - 1, k)) / (T(2) * u.grid().dy());
}

template<typename T>
T DzCentral(const GridFunction3D<T>& u, int i, int j, int k) {
    return (u(i, j, k + 1) - u(i, j, k - 1)) / (T(2) * u.grid().dz());
}

// Second derivatives
template<typename T>
T Dxx(const GridFunction3D<T>& u, int i, int j, int k) {
    T dx = u.grid().dx();
    return (u(i + 1, j, k) - T(2) * u(i, j, k) + u(i - 1, j, k)) / (dx * dx);
}

template<typename T>
T Dyy(const GridFunction3D<T>& u, int i, int j, int k) {
    T dy = u.grid().dy();
    return (u(i, j + 1, k) - T(2) * u(i, j, k) + u(i, j - 1, k)) / (dy * dy);
}

template<typename T>
T Dzz(const GridFunction3D<T>& u, int i, int j, int k) {
    T dz = u.grid().dz();
    return (u(i, j, k + 1) - T(2) * u(i, j, k) + u(i, j, k - 1)) / (dz * dz);
}

// 3D Laplacian (7-point stencil): ∇²u = Dxx + Dyy + Dzz
template<typename T>
T Laplacian3D(const GridFunction3D<T>& u, int i, int j, int k) {
    return Dxx(u, i, j, k) + Dyy(u, i, j, k) + Dzz(u, i, j, k);
}

// Gradient magnitude 3D
template<typename T>
T GradientMagnitude3D(const GridFunction3D<T>& u, int i, int j, int k) {
    T ux = DxCentral(u, i, j, k);
    T uy = DyCentral(u, i, j, k);
    T uz = DzCentral(u, i, j, k);
    return std::sqrt(ux * ux + uy * uy + uz * uz);
}

//=============================================================================
// Apply Laplacian to Entire Grid Function
//=============================================================================

// Apply 1D Laplacian to all interior points
template<typename T>
GridFunction1D<T> applyLaplacian(const GridFunction1D<T>& u) {
    GridFunction1D<T> result(u.grid());
    
    const auto& grid = u.grid();
    for (int i = 1; i < grid.numCells(); ++i) {
        result(i) = Laplacian1D(u, i);
    }
    // Boundary values remain zero
    return result;
}

// Apply 2D Laplacian to all interior points
template<typename T>
GridFunction2D<T> applyLaplacian(const GridFunction2D<T>& u) {
    GridFunction2D<T> result(u.grid());
    
    const auto& grid = u.grid();
    for (int j = 1; j < grid.ny(); ++j) {
        for (int i = 1; i < grid.nx(); ++i) {
            result(i, j) = Laplacian2D(u, i, j);
        }
    }
    return result;
}

// Apply 3D Laplacian to all interior points
template<typename T>
GridFunction3D<T> applyLaplacian(const GridFunction3D<T>& u) {
    GridFunction3D<T> result(u.grid());
    
    const auto& grid = u.grid();
    for (int k = 1; k < grid.nz(); ++k) {
        for (int j = 1; j < grid.ny(); ++j) {
            for (int i = 1; i < grid.nx(); ++i) {
                result(i, j, k) = Laplacian3D(u, i, j, k);
            }
        }
    }
    return result;
}

//=============================================================================
// Upwind Schemes for Advection
//=============================================================================

// Upwind first derivative in x based on velocity sign
template<typename T>
T DxUpwind(const GridFunction2D<T>& u, int i, int j, T velocity) {
    T dx = u.grid().dx();
    if (velocity >= T(0)) {
        // Information coming from left
        return (u(i, j) - u(i - 1, j)) / dx;
    } else {
        // Information coming from right
        return (u(i + 1, j) - u(i, j)) / dx;
    }
}

template<typename T>
T DyUpwind(const GridFunction2D<T>& u, int i, int j, T velocity) {
    T dy = u.grid().dy();
    if (velocity >= T(0)) {
        return (u(i, j) - u(i, j - 1)) / dy;
    } else {
        return (u(i, j + 1) - u(i, j)) / dy;
    }
}

//=============================================================================
// Flux-Conservative Form
//=============================================================================

// For conservation laws: ∂u/∂t + ∂f(u)/∂x = 0
// Compute flux difference at cell faces

template<typename T>
T fluxDifference1D(const GridFunction1D<T>& flux, int i) {
    T dx = flux.grid().dx();
    // flux(i) represents flux at face i+1/2
    // ∂f/∂x ≈ (f_{i+1/2} - f_{i-1/2}) / dx
    return (flux(i) - flux(i - 1)) / dx;
}

//=============================================================================
// Stencil Coefficients
//=============================================================================

// Get Laplacian stencil coefficients for matrix assembly
// For 1D: returns [left, center, right] coefficients
template<typename T>
std::array<T, 3> laplacianCoeffs1D(T dx) {
    T dx2 = dx * dx;
    return {T(1) / dx2, T(-2) / dx2, T(1) / dx2};
}

// For 2D: returns coefficients for [center, left, right, bottom, top]
template<typename T>
std::array<T, 5> laplacianCoeffs2D(T dx, T dy) {
    T dx2 = dx * dx;
    T dy2 = dy * dy;
    return {
        T(-2) / dx2 + T(-2) / dy2,  // center
        T(1) / dx2,                  // left
        T(1) / dx2,                  // right
        T(1) / dy2,                  // bottom
        T(1) / dy2                   // top
    };
}

// For 3D: returns coefficients for [center, ±x, ±y, ±z]
template<typename T>
std::array<T, 7> laplacianCoeffs3D(T dx, T dy, T dz) {
    T dx2 = dx * dx;
    T dy2 = dy * dy;
    T dz2 = dz * dz;
    return {
        T(-2) / dx2 + T(-2) / dy2 + T(-2) / dz2,  // center
        T(1) / dx2,   // x-
        T(1) / dx2,   // x+
        T(1) / dy2,   // y-
        T(1) / dy2,   // y+
        T(1) / dz2,   // z-
        T(1) / dz2    // z+
    };
}

} // namespace MML::PDE

#endif // MML_PDE_STENCILS_H
