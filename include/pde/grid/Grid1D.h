///////////////////////////////////////////////////////////////////////////////////////////
// Grid1D.h - 1D uniform structured grid
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Implements a uniform 1D grid on interval [xmin, xmax] with n+1 nodes.
// Grid points: x_i = xmin + i * dx, for i = 0, 1, ..., n
// Grid spacing: dx = (xmax - xmin) / n
//
// Node 0 and node n are boundary nodes.
// Nodes 1 through n-1 are interior nodes.
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_GRID1D_H
#define MML_PDE_GRID1D_H

#include "Grid.h"

namespace MML::PDE {

template<typename T>
class Grid1D : public GridBase<T, 1> {
    Interval<T> domain_;
    int n_;        // Number of cells (nodes = n + 1)
    T dx_;         // Grid spacing
    
public:
    //=========================================================================
    // Construction
    //=========================================================================
    
    // Create grid on [0, 1] with n cells
    explicit Grid1D(int n) 
        : domain_(T(0), T(1)), n_(n), dx_(T(1) / T(n)) 
    {
        if (n < 1) throw std::invalid_argument("Grid1D: n must be >= 1");
    }
    
    // Create grid on [xmin, xmax] with n cells
    Grid1D(T xmin, T xmax, int n) 
        : domain_(xmin, xmax), n_(n), dx_((xmax - xmin) / T(n))
    {
        if (n < 1) throw std::invalid_argument("Grid1D: n must be >= 1");
        if (xmax <= xmin) throw std::invalid_argument("Grid1D: xmax must be > xmin");
    }
    
    // Create grid from interval with n cells
    Grid1D(const Interval<T>& interval, int n)
        : domain_(interval), n_(n), dx_(interval.length() / T(n))
    {
        if (n < 1) throw std::invalid_argument("Grid1D: n must be >= 1");
    }
    
    //=========================================================================
    // Grid Properties
    //=========================================================================
    
    // Domain
    const Interval<T>& domain() const { return domain_; }
    T xmin() const { return domain_.xmin; }
    T xmax() const { return domain_.xmax; }
    T length() const { return domain_.length(); }
    
    // Discretization
    int numCells() const { return n_; }
    int numNodes() const { return n_ + 1; }
    int numInteriorNodes() const { return n_ - 1; }
    T dx() const { return dx_; }
    
    // GridBase interface
    int totalNodes() const override { return n_ + 1; }
    
    //=========================================================================
    // Coordinate Access
    //=========================================================================
    
    // Get x coordinate of node i
    T x(int i) const {
        return domain_.xmin + i * dx_;
    }
    
    // Get node index containing point x (floor)
    int nodeIndex(T xp) const {
        int i = static_cast<int>((xp - domain_.xmin) / dx_);
        return std::max(0, std::min(i, n_));
    }
    
    // Get interpolation parameter: x = x(i) + alpha * dx, alpha in [0,1]
    void locate(T xp, int& i, T& alpha) const {
        T scaled = (xp - domain_.xmin) / dx_;
        i = static_cast<int>(scaled);
        i = std::max(0, std::min(i, n_ - 1));
        alpha = scaled - i;
    }
    
    //=========================================================================
    // Index Validation
    //=========================================================================
    
    bool isValidNode(int i) const {
        return i >= 0 && i <= n_;
    }
    
    //=========================================================================
    // Boundary Information
    //=========================================================================
    
    NodeLocation nodeLocation(int i) const override {
        if (i == 0 || i == n_) return NodeLocation::Boundary;
        if (i > 0 && i < n_) return NodeLocation::Interior;
        return NodeLocation::Boundary;  // Invalid treated as boundary
    }
    
    bool isLeftBoundary(int i) const { return i == 0; }
    bool isRightBoundary(int i) const { return i == n_; }
    
    // Get boundary node indices
    int leftBoundaryNode() const { return 0; }
    int rightBoundaryNode() const { return n_; }
    
    // Get all interior node indices
    std::vector<int> interiorNodes() const {
        std::vector<int> nodes;
        nodes.reserve(n_ - 1);
        for (int i = 1; i < n_; ++i) {
            nodes.push_back(i);
        }
        return nodes;
    }
    
    // Get all boundary node indices
    std::vector<int> boundaryNodes() const {
        return {0, n_};
    }
    
    //=========================================================================
    // Neighbor Access
    //=========================================================================
    
    // Get left neighbor (-1 if at boundary)
    int left(int i) const { return i > 0 ? i - 1 : -1; }
    
    // Get right neighbor (-1 if at boundary)
    int right(int i) const { return i < n_ ? i + 1 : -1; }
    
    // Check if node has neighbors in both directions
    bool hasFullStencil(int i) const {
        return i > 0 && i < n_;
    }
    
    //=========================================================================
    // Iteration Helpers
    //=========================================================================
    
    // Iterate over all nodes
    template<typename Func>
    void forEachNode(Func&& f) const {
        for (int i = 0; i <= n_; ++i) {
            f(i, x(i));
        }
    }
    
    // Iterate over interior nodes only
    template<typename Func>
    void forEachInterior(Func&& f) const {
        for (int i = 1; i < n_; ++i) {
            f(i, x(i));
        }
    }
    
    // Iterate over boundary nodes
    template<typename Func>
    void forEachBoundary(Func&& f) const {
        f(0, x(0), BoundarySide::Left);
        f(n_, x(n_), BoundarySide::Right);
    }
    
    //=========================================================================
    // Grid Refinement
    //=========================================================================
    
    // Create a refined grid with 2x resolution
    Grid1D refine() const {
        return Grid1D(domain_, 2 * n_);
    }
    
    // Create a coarsened grid with 0.5x resolution
    Grid1D coarsen() const {
        if (n_ < 2) throw std::runtime_error("Cannot coarsen grid with n < 2");
        return Grid1D(domain_, n_ / 2);
    }
    
    //=========================================================================
    // Output
    //=========================================================================
    
    void print(std::ostream& os = std::cout) const {
        os << "Grid1D: [" << domain_.xmin << ", " << domain_.xmax << "]\n";
        os << "  Cells: " << n_ << ", Nodes: " << numNodes() << "\n";
        os << "  dx = " << dx_ << "\n";
    }
};

// Factory functions
template<typename T>
Grid1D<T> makeGrid1D(int n) {
    return Grid1D<T>(n);
}

template<typename T>
Grid1D<T> makeGrid1D(T xmin, T xmax, int n) {
    return Grid1D<T>(xmin, xmax, n);
}

} // namespace MML::PDE

#endif // MML_PDE_GRID1D_H
