///////////////////////////////////////////////////////////////////////////////////////////
// Grid2D.h - 2D uniform structured grid
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Implements a uniform 2D grid on rectangle [xmin,xmax] × [ymin,ymax].
// Grid points: (x_i, y_j) = (xmin + i*dx, ymin + j*dy)
//   for i = 0,...,nx and j = 0,...,ny
//
// Linear indexing: idx = i + j * (nx + 1)  (row-major order)
//
// Boundary nodes: i=0, i=nx, j=0, j=ny
// Interior nodes: 0 < i < nx AND 0 < j < ny
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_GRID2D_H
#define MML_PDE_GRID2D_H

#include "Grid.h"

namespace MML::PDE {

template<typename T>
class Grid2D : public GridBase<T, 2> {
    Rectangle<T> domain_;
    int nx_, ny_;   // Number of cells in each direction
    T dx_, dy_;     // Grid spacing
    
public:
    //=========================================================================
    // Construction
    //=========================================================================
    
    // Create grid on [0,1]² with nx × ny cells
    Grid2D(int nx, int ny)
        : domain_(), nx_(nx), ny_(ny),
          dx_(T(1) / T(nx)), dy_(T(1) / T(ny))
    {
        validate();
    }
    
    // Create grid on [xmin,xmax] × [ymin,ymax] with nx × ny cells
    Grid2D(T xmin, T xmax, T ymin, T ymax, int nx, int ny)
        : domain_(xmin, xmax, ymin, ymax), nx_(nx), ny_(ny),
          dx_((xmax - xmin) / T(nx)), dy_((ymax - ymin) / T(ny))
    {
        validate();
    }
    
    // Create grid from rectangle with nx × ny cells
    Grid2D(const Rectangle<T>& rect, int nx, int ny)
        : domain_(rect), nx_(nx), ny_(ny),
          dx_(rect.width() / T(nx)), dy_(rect.height() / T(ny))
    {
        validate();
    }
    
    // Create square grid on [0,1]² with n × n cells
    explicit Grid2D(int n) : Grid2D(n, n) {}
    
private:
    void validate() {
        if (nx_ < 1 || ny_ < 1) 
            throw std::invalid_argument("Grid2D: nx and ny must be >= 1");
        if (domain_.width() <= 0 || domain_.height() <= 0)
            throw std::invalid_argument("Grid2D: domain must have positive size");
    }
    
public:
    //=========================================================================
    // Grid Properties
    //=========================================================================
    
    // Domain
    const Rectangle<T>& domain() const { return domain_; }
    T xmin() const { return domain_.xMin(); }
    T xmax() const { return domain_.xMax(); }
    T ymin() const { return domain_.yMin(); }
    T ymax() const { return domain_.yMax(); }
    
    // Discretization
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int numNodesX() const { return nx_ + 1; }
    int numNodesY() const { return ny_ + 1; }
    int numCells() const { return nx_ * ny_; }
    int numNodes() const { return (nx_ + 1) * (ny_ + 1); }
    int numInteriorNodes() const { return (nx_ - 1) * (ny_ - 1); }
    T dx() const { return dx_; }
    T dy() const { return dy_; }
    T h() const { return std::min(dx_, dy_); }  // Minimum spacing
    
    // GridBase interface
    int totalNodes() const override { return numNodes(); }
    
    //=========================================================================
    // Index Conversion
    //=========================================================================
    
    // Multi-index to linear
    int index(int i, int j) const {
        return toLinear(i, j, nx_ + 1, ny_ + 1);
    }
    
    // Linear to multi-index
    void index(int idx, int& i, int& j) const {
        fromLinear(idx, nx_ + 1, ny_ + 1, i, j);
    }
    
    // Get i from linear index
    int indexI(int idx) const { return idx % (nx_ + 1); }
    
    // Get j from linear index
    int indexJ(int idx) const { return idx / (nx_ + 1); }
    
    //=========================================================================
    // Coordinate Access
    //=========================================================================
    
    // Get x coordinate of node (i, j)
    T x(int i) const { return domain_.xMin() + i * dx_; }
    T x(int i, int /*j*/) const { return domain_.xMin() + i * dx_; }
    
    // Get y coordinate of node (i, j)
    T y(int j) const { return domain_.yMin() + j * dy_; }
    T y(int /*i*/, int j) const { return domain_.yMin() + j * dy_; }
    
    // Get coordinates from linear index
    void coords(int idx, T& xc, T& yc) const {
        int i, j;
        index(idx, i, j);
        xc = x(i);
        yc = y(j);
    }
    
    // Get coordinates as pair
    std::pair<T, T> coords(int i, int j) const {
        return {x(i), y(j)};
    }
    
    // Locate point in grid
    void locate(T xp, T yp, int& i, int& j, T& alpha, T& beta) const {
        T scaledX = (xp - domain_.xMin()) / dx_;
        T scaledY = (yp - domain_.yMin()) / dy_;
        i = static_cast<int>(scaledX);
        j = static_cast<int>(scaledY);
        i = std::max(0, std::min(i, nx_ - 1));
        j = std::max(0, std::min(j, ny_ - 1));
        alpha = scaledX - i;
        beta = scaledY - j;
    }
    
    //=========================================================================
    // Index Validation
    //=========================================================================
    
    bool isValidNode(int i, int j) const {
        return i >= 0 && i <= nx_ && j >= 0 && j <= ny_;
    }
    
    bool isValidIndex(int idx) const {
        return idx >= 0 && idx < numNodes();
    }
    
    //=========================================================================
    // Boundary Information
    //=========================================================================
    
    NodeLocation nodeLocation(int idx) const override {
        int i, j;
        index(idx, i, j);
        return nodeLocation(i, j);
    }
    
    NodeLocation nodeLocation(int i, int j) const {
        bool onXBoundary = (i == 0 || i == nx_);
        bool onYBoundary = (j == 0 || j == ny_);
        
        if (onXBoundary && onYBoundary) return NodeLocation::Corner;
        if (onXBoundary || onYBoundary) return NodeLocation::Boundary;
        return NodeLocation::Interior;
    }
    
    bool isInterior(int i, int j) const {
        return i > 0 && i < nx_ && j > 0 && j < ny_;
    }
    
    bool isBoundary(int i, int j) const {
        return !isInterior(i, j);
    }
    
    // Boundary side checks
    bool isLeftBoundary(int i, int /*j*/) const { return i == 0; }
    bool isRightBoundary(int i, int /*j*/) const { return i == nx_; }
    bool isBottomBoundary(int /*i*/, int j) const { return j == 0; }
    bool isTopBoundary(int /*i*/, int j) const { return j == ny_; }
    
    // Get boundary side for a node
    std::vector<BoundarySide> boundarySides(int i, int j) const {
        std::vector<BoundarySide> sides;
        if (i == 0) sides.push_back(BoundarySide::Left);
        if (i == nx_) sides.push_back(BoundarySide::Right);
        if (j == 0) sides.push_back(BoundarySide::Bottom);
        if (j == ny_) sides.push_back(BoundarySide::Top);
        return sides;
    }
    
    // Get all interior node indices
    std::vector<int> interiorNodes() const {
        std::vector<int> nodes;
        nodes.reserve(numInteriorNodes());
        for (int j = 1; j < ny_; ++j) {
            for (int i = 1; i < nx_; ++i) {
                nodes.push_back(index(i, j));
            }
        }
        return nodes;
    }
    
    // Get boundary nodes on a specific side
    std::vector<int> boundaryNodes(BoundarySide side) const {
        std::vector<int> nodes;
        switch (side) {
            case BoundarySide::Left:
                for (int j = 0; j <= ny_; ++j) nodes.push_back(index(0, j));
                break;
            case BoundarySide::Right:
                for (int j = 0; j <= ny_; ++j) nodes.push_back(index(nx_, j));
                break;
            case BoundarySide::Bottom:
                for (int i = 0; i <= nx_; ++i) nodes.push_back(index(i, 0));
                break;
            case BoundarySide::Top:
                for (int i = 0; i <= nx_; ++i) nodes.push_back(index(i, ny_));
                break;
            default:
                break;
        }
        return nodes;
    }
    
    // Get all boundary nodes
    std::vector<int> allBoundaryNodes() const {
        std::vector<int> nodes;
        // Bottom and top
        for (int i = 0; i <= nx_; ++i) {
            nodes.push_back(index(i, 0));
            nodes.push_back(index(i, ny_));
        }
        // Left and right (excluding corners already added)
        for (int j = 1; j < ny_; ++j) {
            nodes.push_back(index(0, j));
            nodes.push_back(index(nx_, j));
        }
        return nodes;
    }
    
    //=========================================================================
    // Neighbor Access
    //=========================================================================
    
    // Direct neighbors (4-connectivity)
    int left(int i, int j) const { return i > 0 ? index(i-1, j) : -1; }
    int right(int i, int j) const { return i < nx_ ? index(i+1, j) : -1; }
    int bottom(int i, int j) const { return j > 0 ? index(i, j-1) : -1; }
    int top(int i, int j) const { return j < ny_ ? index(i, j+1) : -1; }
    
    // From linear index
    int left(int idx) const { int i, j; index(idx, i, j); return left(i, j); }
    int right(int idx) const { int i, j; index(idx, i, j); return right(i, j); }
    int bottom(int idx) const { int i, j; index(idx, i, j); return bottom(i, j); }
    int top(int idx) const { int i, j; index(idx, i, j); return top(i, j); }
    
    // Diagonal neighbors (8-connectivity)
    int bottomLeft(int i, int j) const { 
        return (i > 0 && j > 0) ? index(i-1, j-1) : -1; 
    }
    int bottomRight(int i, int j) const { 
        return (i < nx_ && j > 0) ? index(i+1, j-1) : -1; 
    }
    int topLeft(int i, int j) const { 
        return (i > 0 && j < ny_) ? index(i-1, j+1) : -1; 
    }
    int topRight(int i, int j) const { 
        return (i < nx_ && j < ny_) ? index(i+1, j+1) : -1; 
    }
    
    // Check if node has all 4 direct neighbors (for 5-point stencil)
    bool hasFullStencil(int i, int j) const {
        return i > 0 && i < nx_ && j > 0 && j < ny_;
    }
    
    // Get all neighbors of a node
    std::vector<int> neighbors4(int i, int j) const {
        std::vector<int> n;
        if (i > 0) n.push_back(index(i-1, j));
        if (i < nx_) n.push_back(index(i+1, j));
        if (j > 0) n.push_back(index(i, j-1));
        if (j < ny_) n.push_back(index(i, j+1));
        return n;
    }
    
    std::vector<int> neighbors8(int i, int j) const {
        std::vector<int> n = neighbors4(i, j);
        if (i > 0 && j > 0) n.push_back(index(i-1, j-1));
        if (i < nx_ && j > 0) n.push_back(index(i+1, j-1));
        if (i > 0 && j < ny_) n.push_back(index(i-1, j+1));
        if (i < nx_ && j < ny_) n.push_back(index(i+1, j+1));
        return n;
    }
    
    //=========================================================================
    // Iteration Helpers
    //=========================================================================
    
    // Iterate over all nodes
    template<typename Func>
    void forEachNode(Func&& f) const {
        for (int j = 0; j <= ny_; ++j) {
            for (int i = 0; i <= nx_; ++i) {
                f(i, j, x(i), y(j));
            }
        }
    }
    
    // Iterate over interior nodes only
    template<typename Func>
    void forEachInterior(Func&& f) const {
        for (int j = 1; j < ny_; ++j) {
            for (int i = 1; i < nx_; ++i) {
                f(i, j, x(i), y(j));
            }
        }
    }
    
    // Iterate over boundary nodes
    template<typename Func>
    void forEachBoundary(Func&& f) const {
        // Bottom
        for (int i = 0; i <= nx_; ++i) f(i, 0, x(i), y(0), BoundarySide::Bottom);
        // Top
        for (int i = 0; i <= nx_; ++i) f(i, ny_, x(i), y(ny_), BoundarySide::Top);
        // Left (excluding corners)
        for (int j = 1; j < ny_; ++j) f(0, j, x(0), y(j), BoundarySide::Left);
        // Right (excluding corners)
        for (int j = 1; j < ny_; ++j) f(nx_, j, x(nx_), y(j), BoundarySide::Right);
    }
    
    //=========================================================================
    // Grid Refinement
    //=========================================================================
    
    Grid2D refine() const {
        return Grid2D(domain_, 2 * nx_, 2 * ny_);
    }
    
    Grid2D coarsen() const {
        if (nx_ < 2 || ny_ < 2) 
            throw std::runtime_error("Cannot coarsen grid with nx or ny < 2");
        return Grid2D(domain_, nx_ / 2, ny_ / 2);
    }
    
    //=========================================================================
    // Output
    //=========================================================================
    
    void print(std::ostream& os = std::cout) const {
        os << "Grid2D: [" << domain_.xMin() << ", " << domain_.xMax() << "] × ["
           << domain_.yMin() << ", " << domain_.yMax() << "]\n";
        os << "  Cells: " << nx_ << " × " << ny_ << " = " << numCells() << "\n";
        os << "  Nodes: " << numNodesX() << " × " << numNodesY() 
           << " = " << numNodes() << "\n";
        os << "  Interior: " << numInteriorNodes() << "\n";
        os << "  dx = " << dx_ << ", dy = " << dy_ << "\n";
    }
};

// Factory functions
template<typename T>
Grid2D<T> makeGrid2D(int n) {
    return Grid2D<T>(n);
}

template<typename T>
Grid2D<T> makeGrid2D(int nx, int ny) {
    return Grid2D<T>(nx, ny);
}

template<typename T>
Grid2D<T> makeGrid2D(T xmin, T xmax, T ymin, T ymax, int nx, int ny) {
    return Grid2D<T>(xmin, xmax, ymin, ymax, nx, ny);
}

} // namespace MML::PDE

#endif // MML_PDE_GRID2D_H
