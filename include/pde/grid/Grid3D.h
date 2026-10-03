///////////////////////////////////////////////////////////////////////////////////////////
// Grid3D.h - 3D uniform structured grid
///////////////////////////////////////////////////////////////////////////////////////////
// Part of MinimalMathLibrary - PDE Solver Module
//
// Implements a uniform 3D grid on box [xmin,xmax] × [ymin,ymax] × [zmin,zmax].
// Grid points: (x_i, y_j, z_k) = (xmin + i*dx, ymin + j*dy, zmin + k*dz)
//   for i = 0,...,nx, j = 0,...,ny, k = 0,...,nz
//
// Linear indexing: idx = i + j*(nx+1) + k*(nx+1)*(ny+1)  (row-major order)
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_PDE_GRID3D_H
#define MML_PDE_GRID3D_H

#include "Grid.h"

namespace MML::PDE {

template<typename T>
class Grid3D : public GridBase<T, 3> {
    Box<T> domain_;
    int nx_, ny_, nz_;   // Number of cells in each direction
    T dx_, dy_, dz_;     // Grid spacing
    
public:
    //=========================================================================
    // Construction
    //=========================================================================
    
    // Create grid on [0,1]³ with nx × ny × nz cells
    Grid3D(int nx, int ny, int nz)
        : domain_(), nx_(nx), ny_(ny), nz_(nz),
          dx_(T(1) / T(nx)), dy_(T(1) / T(ny)), dz_(T(1) / T(nz))
    {
        validate();
    }
    
    // Create cubic grid on [0,1]³ with n × n × n cells
    explicit Grid3D(int n) : Grid3D(n, n, n) {}
    
    // Create grid on specified box with nx × ny × nz cells
    Grid3D(T xmin, T xmax, T ymin, T ymax, T zmin, T zmax, int nx, int ny, int nz)
        : domain_(xmin, xmax, ymin, ymax, zmin, zmax), 
          nx_(nx), ny_(ny), nz_(nz),
          dx_((xmax - xmin) / T(nx)), 
          dy_((ymax - ymin) / T(ny)),
          dz_((zmax - zmin) / T(nz))
    {
        validate();
    }
    
    // Create grid from box with nx × ny × nz cells
    Grid3D(const Box<T>& box, int nx, int ny, int nz)
        : domain_(box), nx_(nx), ny_(ny), nz_(nz),
          dx_(box.xLength() / T(nx)),
          dy_(box.yLength() / T(ny)),
          dz_(box.zLength() / T(nz))
    {
        validate();
    }
    
private:
    void validate() {
        if (nx_ < 1 || ny_ < 1 || nz_ < 1)
            throw std::invalid_argument("Grid3D: nx, ny, nz must be >= 1");
        if (domain_.xLength() <= 0 || domain_.yLength() <= 0 || domain_.zLength() <= 0)
            throw std::invalid_argument("Grid3D: domain must have positive size");
    }
    
public:
    //=========================================================================
    // Grid Properties
    //=========================================================================
    
    // Domain
    const Box<T>& domain() const { return domain_; }
    T xmin() const { return domain_.xMin(); }
    T xmax() const { return domain_.xMax(); }
    T ymin() const { return domain_.yMin(); }
    T ymax() const { return domain_.yMax(); }
    T zmin() const { return domain_.zMin(); }
    T zmax() const { return domain_.zMax(); }
    
    // Discretization
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    int numNodesX() const { return nx_ + 1; }
    int numNodesY() const { return ny_ + 1; }
    int numNodesZ() const { return nz_ + 1; }
    int numCells() const { return nx_ * ny_ * nz_; }
    int numNodes() const { return (nx_ + 1) * (ny_ + 1) * (nz_ + 1); }
    int numInteriorNodes() const { return (nx_ - 1) * (ny_ - 1) * (nz_ - 1); }
    T dx() const { return dx_; }
    T dy() const { return dy_; }
    T dz() const { return dz_; }
    T h() const { return std::min({dx_, dy_, dz_}); }  // Minimum spacing
    
    // GridBase interface
    int totalNodes() const override { return numNodes(); }
    
    //=========================================================================
    // Index Conversion
    //=========================================================================
    
    // Multi-index to linear
    int index(int i, int j, int k) const {
        return toLinear(i, j, k, nx_ + 1, ny_ + 1, nz_ + 1);
    }
    
    // Linear to multi-index
    void index(int idx, int& i, int& j, int& k) const {
        fromLinear(idx, nx_ + 1, ny_ + 1, nz_ + 1, i, j, k);
    }
    
    int indexI(int idx) const { return idx % (nx_ + 1); }
    int indexJ(int idx) const { return (idx / (nx_ + 1)) % (ny_ + 1); }
    int indexK(int idx) const { return idx / ((nx_ + 1) * (ny_ + 1)); }
    
    //=========================================================================
    // Coordinate Access
    //=========================================================================
    
    T x(int i) const { return domain_.xMin() + i * dx_; }
    T y(int j) const { return domain_.yMin() + j * dy_; }
    T z(int k) const { return domain_.zMin() + k * dz_; }
    
    void coords(int idx, T& xc, T& yc, T& zc) const {
        int i, j, k;
        index(idx, i, j, k);
        xc = x(i);
        yc = y(j);
        zc = z(k);
    }
    
    std::tuple<T, T, T> coords(int i, int j, int k) const {
        return {x(i), y(j), z(k)};
    }
    
    //=========================================================================
    // Index Validation
    //=========================================================================
    
    bool isValidNode(int i, int j, int k) const {
        return i >= 0 && i <= nx_ && j >= 0 && j <= ny_ && k >= 0 && k <= nz_;
    }
    
    bool isValidIndex(int idx) const {
        return idx >= 0 && idx < numNodes();
    }
    
    //=========================================================================
    // Boundary Information
    //=========================================================================
    
    NodeLocation nodeLocation(int idx) const override {
        int i, j, k;
        index(idx, i, j, k);
        return nodeLocation(i, j, k);
    }
    
    NodeLocation nodeLocation(int i, int j, int k) const {
        int numBoundaries = 0;
        if (i == 0 || i == nx_) ++numBoundaries;
        if (j == 0 || j == ny_) ++numBoundaries;
        if (k == 0 || k == nz_) ++numBoundaries;
        
        if (numBoundaries >= 2) return NodeLocation::Corner;  // Edge or corner
        if (numBoundaries == 1) return NodeLocation::Boundary;
        return NodeLocation::Interior;
    }
    
    bool isInterior(int i, int j, int k) const {
        return i > 0 && i < nx_ && j > 0 && j < ny_ && k > 0 && k < nz_;
    }
    
    bool isBoundary(int i, int j, int k) const {
        return !isInterior(i, j, k);
    }
    
    // Boundary face checks
    bool isXMinBoundary(int i, int, int) const { return i == 0; }
    bool isXMaxBoundary(int i, int, int) const { return i == nx_; }
    bool isYMinBoundary(int, int j, int) const { return j == 0; }
    bool isYMaxBoundary(int, int j, int) const { return j == ny_; }
    bool isZMinBoundary(int, int, int k) const { return k == 0; }
    bool isZMaxBoundary(int, int, int k) const { return k == nz_; }
    
    // Get all interior node indices
    std::vector<int> interiorNodes() const {
        std::vector<int> nodes;
        nodes.reserve(numInteriorNodes());
        for (int k = 1; k < nz_; ++k) {
            for (int j = 1; j < ny_; ++j) {
                for (int i = 1; i < nx_; ++i) {
                    nodes.push_back(index(i, j, k));
                }
            }
        }
        return nodes;
    }
    
    // Get boundary nodes on a specific face
    std::vector<int> boundaryNodes(BoundarySide side) const {
        std::vector<int> nodes;
        switch (side) {
            case BoundarySide::Left:   // x = xmin
                for (int k = 0; k <= nz_; ++k)
                    for (int j = 0; j <= ny_; ++j)
                        nodes.push_back(index(0, j, k));
                break;
            case BoundarySide::Right:  // x = xmax
                for (int k = 0; k <= nz_; ++k)
                    for (int j = 0; j <= ny_; ++j)
                        nodes.push_back(index(nx_, j, k));
                break;
            case BoundarySide::Bottom: // y = ymin
                for (int k = 0; k <= nz_; ++k)
                    for (int i = 0; i <= nx_; ++i)
                        nodes.push_back(index(i, 0, k));
                break;
            case BoundarySide::Top:    // y = ymax
                for (int k = 0; k <= nz_; ++k)
                    for (int i = 0; i <= nx_; ++i)
                        nodes.push_back(index(i, ny_, k));
                break;
            case BoundarySide::Front:  // z = zmin
                for (int j = 0; j <= ny_; ++j)
                    for (int i = 0; i <= nx_; ++i)
                        nodes.push_back(index(i, j, 0));
                break;
            case BoundarySide::Back:   // z = zmax
                for (int j = 0; j <= ny_; ++j)
                    for (int i = 0; i <= nx_; ++i)
                        nodes.push_back(index(i, j, nz_));
                break;
        }
        return nodes;
    }
    
    //=========================================================================
    // Neighbor Access (6-connectivity)
    //=========================================================================
    
    int xMinus(int i, int j, int k) const { return i > 0 ? index(i-1, j, k) : -1; }
    int xPlus(int i, int j, int k) const { return i < nx_ ? index(i+1, j, k) : -1; }
    int yMinus(int i, int j, int k) const { return j > 0 ? index(i, j-1, k) : -1; }
    int yPlus(int i, int j, int k) const { return j < ny_ ? index(i, j+1, k) : -1; }
    int zMinus(int i, int j, int k) const { return k > 0 ? index(i, j, k-1) : -1; }
    int zPlus(int i, int j, int k) const { return k < nz_ ? index(i, j, k+1) : -1; }
    
    // Check if node has all 6 direct neighbors (for 7-point stencil)
    bool hasFullStencil(int i, int j, int k) const {
        return i > 0 && i < nx_ && j > 0 && j < ny_ && k > 0 && k < nz_;
    }
    
    // Get all 6 neighbors
    std::vector<int> neighbors6(int i, int j, int k) const {
        std::vector<int> n;
        if (i > 0) n.push_back(index(i-1, j, k));
        if (i < nx_) n.push_back(index(i+1, j, k));
        if (j > 0) n.push_back(index(i, j-1, k));
        if (j < ny_) n.push_back(index(i, j+1, k));
        if (k > 0) n.push_back(index(i, j, k-1));
        if (k < nz_) n.push_back(index(i, j, k+1));
        return n;
    }
    
    //=========================================================================
    // Iteration Helpers
    //=========================================================================
    
    // Iterate over all nodes
    template<typename Func>
    void forEachNode(Func&& f) const {
        for (int k = 0; k <= nz_; ++k) {
            for (int j = 0; j <= ny_; ++j) {
                for (int i = 0; i <= nx_; ++i) {
                    f(i, j, k, x(i), y(j), z(k));
                }
            }
        }
    }
    
    // Iterate over interior nodes only
    template<typename Func>
    void forEachInterior(Func&& f) const {
        for (int k = 1; k < nz_; ++k) {
            for (int j = 1; j < ny_; ++j) {
                for (int i = 1; i < nx_; ++i) {
                    f(i, j, k, x(i), y(j), z(k));
                }
            }
        }
    }
    
    // Iterate over boundary nodes
    // Callback: f(i, j, k, x, y, z, side)
    template<typename Func>
    void forEachBoundary(Func&& f) const {
        // Left face (i = 0)
        for (int k = 0; k <= nz_; ++k) {
            for (int j = 0; j <= ny_; ++j) {
                f(0, j, k, x(0), y(j), z(k), BoundarySide::Left);
            }
        }
        // Right face (i = nx)
        for (int k = 0; k <= nz_; ++k) {
            for (int j = 0; j <= ny_; ++j) {
                f(nx_, j, k, x(nx_), y(j), z(k), BoundarySide::Right);
            }
        }
        // Bottom face (j = 0), excluding edges already counted
        for (int k = 0; k <= nz_; ++k) {
            for (int i = 1; i < nx_; ++i) {
                f(i, 0, k, x(i), y(0), z(k), BoundarySide::Bottom);
            }
        }
        // Top face (j = ny), excluding edges already counted
        for (int k = 0; k <= nz_; ++k) {
            for (int i = 1; i < nx_; ++i) {
                f(i, ny_, k, x(i), y(ny_), z(k), BoundarySide::Top);
            }
        }
        // Front face (k = 0), excluding edges already counted
        for (int j = 1; j < ny_; ++j) {
            for (int i = 1; i < nx_; ++i) {
                f(i, j, 0, x(i), y(j), z(0), BoundarySide::Front);
            }
        }
        // Back face (k = nz), excluding edges already counted
        for (int j = 1; j < ny_; ++j) {
            for (int i = 1; i < nx_; ++i) {
                f(i, j, nz_, x(i), y(j), z(nz_), BoundarySide::Back);
            }
        }
    }
    
    //=========================================================================
    // Grid Refinement
    //=========================================================================
    
    Grid3D refine() const {
        return Grid3D(domain_, 2 * nx_, 2 * ny_, 2 * nz_);
    }
    
    Grid3D coarsen() const {
        if (nx_ < 2 || ny_ < 2 || nz_ < 2)
            throw std::runtime_error("Cannot coarsen grid with nx, ny, or nz < 2");
        return Grid3D(domain_, nx_ / 2, ny_ / 2, nz_ / 2);
    }
    
    //=========================================================================
    // Output
    //=========================================================================
    
    void print(std::ostream& os = std::cout) const {
        os << "Grid3D: [" << domain_.xmin << ", " << domain_.xmax << "] × ["
           << domain_.ymin << ", " << domain_.ymax << "] × ["
           << domain_.zmin << ", " << domain_.zmax << "]\n";
        os << "  Cells: " << nx_ << " × " << ny_ << " × " << nz_ 
           << " = " << numCells() << "\n";
        os << "  Nodes: " << numNodesX() << " × " << numNodesY() << " × " << numNodesZ()
           << " = " << numNodes() << "\n";
        os << "  Interior: " << numInteriorNodes() << "\n";
        os << "  dx = " << dx_ << ", dy = " << dy_ << ", dz = " << dz_ << "\n";
    }
};

// Factory functions
template<typename T>
Grid3D<T> makeGrid3D(int n) {
    return Grid3D<T>(n);
}

template<typename T>
Grid3D<T> makeGrid3D(int nx, int ny, int nz) {
    return Grid3D<T>(nx, ny, nz);
}

} // namespace MML::PDE

#endif // MML_PDE_GRID3D_H
