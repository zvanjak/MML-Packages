///////////////////////////////////////////////////////////////////////////////////////////
// MultiObjective.h
//
// Multi-objective optimization foundation for MML.
// Provides core data structures and algorithms for Pareto-based optimization.
//
// Key Components:
// - ObjectivePoint: Solution with decision variables and objective values (runtime sizes)
// - ParetoRelation: Dominance comparison results
// - Pareto dominance functions: Dominates, WeaklyDominates, IsIncomparable
// - ParetoArchive: Non-dominated solution storage with bounded size
// - FastNonDominatedSort: O(N²) population ranking algorithm
// - CrowdingDistance: Diversity preservation measure for NSGA-II
//
// References:
// - Deb, K. (2001). Multi-objective optimization using evolutionary algorithms.
// - Deb, K., et al. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II.
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_MULTI_OBJECTIVE_H
#define MML_MULTI_OBJECTIVE_H

#include "MMLBase.h"
#include "OptimizationCommon.h"
#include "OptimizationProblem.h"
#include "mml/mml_export.h"

#include <vector>
#include <algorithm>
#include <limits>
#include <numeric>
#include <span>
#include <utility>

namespace MML::Optimization {

///////////////////////////////////////////////////////////////////////////////////////////
// ParetoRelation - Result of dominance comparison between two solutions
///////////////////////////////////////////////////////////////////////////////////////////

enum class ParetoRelation
{
    Dominates,          ///< First solution dominates second (strictly better in all objectives)
    IsDominatedBy,      ///< First solution is dominated by second
    Incomparable,       ///< Neither dominates the other (Pareto-incomparable)
    Equal               ///< Solutions have identical objective values
};


///////////////////////////////////////////////////////////////////////////////////////////
// ObjectivePoint - A solution point in multi-objective space
//
// Stores:
// - Decision vector (runtime number of variables)
// - Objective values (runtime number of objectives)
// - Constraint violation (for constrained MO problems)
// - Metadata for NSGA-II: rank, crowding distance
///////////////////////////////////////////////////////////////////////////////////////////

class ObjectivePoint
{
public:
    using DecisionVec = Vector<Real>;
    using ObjectiveVec = Vector<Real>;

private:
    DecisionVec _x;                    ///< Decision variables
    ObjectiveVec _objectives;          ///< Objective function values
    Real _constraintViolation = 0;     ///< Total constraint violation (0 = feasible)
    
    // NSGA-II metadata
    int _rank = -1;                    ///< Pareto rank (0 = first front, etc.)
    Real _crowdingDistance = 0;        ///< Crowding distance for diversity

public:
    //---------------------------------------------------------------------------------
    // Constructors
    //---------------------------------------------------------------------------------

    /// @brief Default constructor
    ObjectivePoint() = default;

    /// @brief Construct with decision vector only (objectives computed later)
    explicit ObjectivePoint(const DecisionVec& x)
        : _x(x) {}

    /// @brief Construct with decision vector and objective values
    ObjectivePoint(const DecisionVec& x, const ObjectiveVec& objectives)
        : _x(x), _objectives(objectives) {}

    /// @brief Construct with all fields
    ObjectivePoint(const DecisionVec& x, const ObjectiveVec& objectives, 
                   Real constraintViolation)
        : _x(x), _objectives(objectives), _constraintViolation(constraintViolation) {}

    //---------------------------------------------------------------------------------
    // Accessors
    //---------------------------------------------------------------------------------

    /// @brief Get decision vector
    const DecisionVec& GetX() const { return _x; }
    DecisionVec& GetX() { return _x; }

    /// @brief Get objective values
    const ObjectiveVec& GetObjectives() const { return _objectives; }
    ObjectiveVec& GetObjectives() { return _objectives; }

    /// @brief Get specific objective value
    Real GetObjective(int index) const { return _objectives[index]; }

    /// @brief Set objective values
    void SetObjectives(const ObjectiveVec& objectives) { _objectives = objectives; }

    /// @brief Get constraint violation
    Real GetConstraintViolation() const { return _constraintViolation; }
    void SetConstraintViolation(Real cv) { _constraintViolation = cv; }

    /// @brief Check if feasible
    bool IsFeasible() const { return _constraintViolation <= 0; }

    //---------------------------------------------------------------------------------
    // NSGA-II Metadata
    //---------------------------------------------------------------------------------

    /// @brief Get/set Pareto rank
    int GetRank() const { return _rank; }
    void SetRank(int rank) { _rank = rank; }

    /// @brief Get/set crowding distance
    Real GetCrowdingDistance() const { return _crowdingDistance; }
    void SetCrowdingDistance(Real cd) { _crowdingDistance = cd; }

    //---------------------------------------------------------------------------------
    // Runtime dimension information
    //---------------------------------------------------------------------------------

    /// @brief Get number of decision variables
    int GetNumVariables() const { return static_cast<int>(_x.size()); }

    /// @brief Get number of objectives
    int GetNumObjectives() const { return static_cast<int>(_objectives.size()); }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Pareto Dominance Functions
//
// For MINIMIZATION problems:
// - a dominates b if: a[i] <= b[i] for all i, AND a[j] < b[j] for at least one j
// - a weakly dominates b if: a[i] <= b[i] for all i
// - a and b are incomparable if: neither dominates the other
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Check if objectives 'a' dominate objectives 'b' (strict Pareto dominance)
/// @param a First objective vector
/// @param b Second objective vector
/// @return true if a dominates b (a is better in all objectives, strictly better in at least one)
inline bool Dominates(const Vector<Real>& a, const Vector<Real>& b)
{
    bool strictlyBetterInOne = false;
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        if (a[i] > b[i])
            return false;  // a is worse in at least one objective
        if (a[i] < b[i])
            strictlyBetterInOne = true;
    }
    return strictlyBetterInOne;
}

/// @brief Check if objectives 'a' weakly dominate objectives 'b'
/// @param a First objective vector
/// @param b Second objective vector
/// @return true if a[i] <= b[i] for all i
inline bool WeaklyDominates(const Vector<Real>& a, const Vector<Real>& b)
{
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        if (a[i] > b[i])
            return false;
    }
    return true;
}

/// @brief Check if two objective vectors are incomparable (neither dominates the other)
inline bool IsIncomparable(const Vector<Real>& a, const Vector<Real>& b)
{
    return !Dominates(a, b) && !Dominates(b, a);
}

/// @brief Compare two objective vectors and return their Pareto relation
inline ParetoRelation CompareParetoRelation(const Vector<Real>& a, const Vector<Real>& b)
{
    bool aBetterInOne = false;
    bool bBetterInOne = false;
    
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        if (a[i] < b[i])
            aBetterInOne = true;
        else if (a[i] > b[i])
            bBetterInOne = true;
    }
    
    if (aBetterInOne && !bBetterInOne)
        return ParetoRelation::Dominates;
    if (bBetterInOne && !aBetterInOne)
        return ParetoRelation::IsDominatedBy;
    if (!aBetterInOne && !bBetterInOne)
        return ParetoRelation::Equal;
    return ParetoRelation::Incomparable;
}

/// @brief Check if ObjectivePoint 'a' dominates ObjectivePoint 'b'
/// @note For constrained problems, uses Deb's constraint handling rules
inline bool Dominates(const ObjectivePoint& a, const ObjectivePoint& b)
{
    // Deb's constraint handling for constrained multi-objective optimization:
    // 1. Feasible dominates infeasible
    // 2. Among feasible, use Pareto dominance
    // 3. Among infeasible, lower violation dominates
    
    bool aFeasible = a.IsFeasible();
    bool bFeasible = b.IsFeasible();
    
    if (aFeasible && !bFeasible)
        return true;   // Feasible dominates infeasible
    if (!aFeasible && bFeasible)
        return false;  // Infeasible doesn't dominate feasible
    if (!aFeasible && !bFeasible) {
        // Both infeasible: lower constraint violation dominates
        return a.GetConstraintViolation() < b.GetConstraintViolation();
    }
    
    // Both feasible: use standard Pareto dominance
    return Dominates(a.GetObjectives(), b.GetObjectives());
}


///////////////////////////////////////////////////////////////////////////////////////////
// Fast Non-Dominated Sorting
//
// Implements the O(MN²) algorithm from NSGA-II paper.
// Assigns Pareto ranks to each solution in the population.
//
// Returns: Vector of fronts, where fronts[0] is the Pareto-optimal front.
///////////////////////////////////////////////////////////////////////////////////////////

namespace detail {
struct ParetoView {
    std::span<const Real> objectives;
    Real constraintViolation;
};

MML_OPTIMIZATION_API std::vector<std::vector<int>> FastNonDominatedSort(
    const std::vector<ParetoView>& points, std::vector<int>& ranks);
MML_OPTIMIZATION_API void CalculateCrowdingDistance(const std::vector<ParetoView>& points,
                                                    const std::vector<int>& front,
                                                    std::vector<Real>& distances);
}

/// @brief Perform fast non-dominated sorting on a population
/// @param population Population to sort (will have ranks assigned)
/// @return Vector of fronts (each front is a vector of indices into population)
inline std::vector<std::vector<int>> FastNonDominatedSort(std::vector<ObjectivePoint>& population)
{
    std::vector<detail::ParetoView> points;
    points.reserve(population.size());
    for (const auto& point : population) {
        points.push_back({std::span<const Real>(&point.GetObjectives()[0], point.GetObjectives().size()),
                          point.GetConstraintViolation()});
    }

    std::vector<int> ranks(population.size(), -1);
    auto fronts = detail::FastNonDominatedSort(points, ranks);
    for (size_t index = 0; index < population.size(); ++index) {
        population[index].SetRank(ranks[index]);
    }
    return fronts;
}


///////////////////////////////////////////////////////////////////////////////////////////
// Crowding Distance Calculation
//
// Measures solution density in objective space.
// Higher crowding distance = more isolated = better for diversity.
// Boundary solutions get infinite distance.
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Calculate crowding distance for a single front
/// @param population Full population (crowding distances will be assigned)
/// @param front Indices of solutions in this front
inline void CalculateCrowdingDistance(std::vector<ObjectivePoint>& population,
                               const std::vector<int>& front)
{
    if (front.empty()) return;

    std::vector<detail::ParetoView> points;
    points.reserve(population.size());
    for (const auto& point : population) {
        points.push_back({std::span<const Real>(&point.GetObjectives()[0], point.GetObjectives().size()),
                          point.GetConstraintViolation()});
    }

    std::vector<Real> distances(population.size());
    detail::CalculateCrowdingDistance(points, front, distances);
    for (int idx : front) {
        population[idx].SetCrowdingDistance(distances[idx]);
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
// Crowded Comparison Operator
//
// For NSGA-II selection: prefers lower rank; if equal rank, prefers larger crowding distance.
///////////////////////////////////////////////////////////////////////////////////////////

/// @brief Compare two solutions using crowded comparison operator
/// @return true if 'a' is preferred over 'b'
inline bool CrowdedComparison(const ObjectivePoint& a, const ObjectivePoint& b)
{
    // Prefer lower rank (closer to Pareto front)
    if (a.GetRank() < b.GetRank())
        return true;
    if (a.GetRank() > b.GetRank())
        return false;
    
    // Same rank: prefer larger crowding distance (more isolated)
    return a.GetCrowdingDistance() > b.GetCrowdingDistance();
}


///////////////////////////////////////////////////////////////////////////////////////////
// ParetoArchive - Bounded archive of non-dominated solutions
//
// Maintains a set of non-dominated solutions with a maximum size.
// When archive is full and a new non-dominated solution is added,
// removes the most crowded solution to maintain diversity.
///////////////////////////////////////////////////////////////////////////////////////////

class ParetoArchive
{
public:
    using Point = ObjectivePoint;
    
private:
    std::vector<Point> _archive;
    int _maxSize;

public:
    //---------------------------------------------------------------------------------
    // Constructors
    //---------------------------------------------------------------------------------

    /// @brief Create archive with specified maximum size
    explicit ParetoArchive(int maxSize = 100)
        : _maxSize(maxSize) {}

    //---------------------------------------------------------------------------------
    // Archive operations
    //---------------------------------------------------------------------------------

    /// @brief Try to add a solution to the archive
    /// @param point Solution to add
    /// @return true if solution was added, false if dominated by existing solutions
    bool TryAdd(const Point& point)
    {
        // Check if point is dominated by any archive member
        for (const auto& member : _archive) {
            if (Dominates(member, point))
                return false;  // Point is dominated, reject
        }
        
        // Remove any archive members dominated by the new point
        _archive.erase(
            std::remove_if(_archive.begin(), _archive.end(),
                [&point](const Point& member) {
                    return Dominates(point, member);
                }),
            _archive.end());
        
        // Add the new point
        _archive.push_back(point);
        
        // If archive exceeds max size, remove most crowded solution
        if (static_cast<int>(_archive.size()) > _maxSize) {
            TruncateArchive();
        }
        
        return true;
    }

    /// @brief Clear the archive
    void Clear() { _archive.clear(); }

    /// @brief Get archive size
    int Size() const { return static_cast<int>(_archive.size()); }

    /// @brief Get maximum archive size
    int MaxSize() const { return _maxSize; }

    /// @brief Set maximum archive size
    void SetMaxSize(int maxSize) { _maxSize = maxSize; }

    /// @brief Check if archive is empty
    bool IsEmpty() const { return _archive.empty(); }

    /// @brief Get all solutions in archive
    const std::vector<Point>& GetSolutions() const { return _archive; }

    /// @brief Get solution at index
    const Point& operator[](int index) const { return _archive[index]; }

    /// @brief Iterator support
    auto begin() { return _archive.begin(); }
    auto end() { return _archive.end(); }
    auto begin() const { return _archive.begin(); }
    auto end() const { return _archive.end(); }

private:
    /// @brief Remove most crowded solution when archive exceeds max size
    void TruncateArchive()
    {
        if (_archive.size() <= 1) return;
        
        // Calculate crowding distances
        std::vector<int> indices(_archive.size());
        std::iota(indices.begin(), indices.end(), 0);
        CalculateCrowdingDistance(_archive, indices);
        
        // Find and remove solution with smallest crowding distance (most crowded)
        auto minIt = std::min_element(_archive.begin(), _archive.end(),
            [](const Point& a, const Point& b) {
                return a.GetCrowdingDistance() < b.GetCrowdingDistance();
            });
        
        _archive.erase(minIt);
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// MultiObjectiveResult - Result of multi-objective optimization
//
// Contains:
// - Final Pareto front (non-dominated solutions)
// - Optional: full archive history
// - Quality metrics (hypervolume, spread, etc.)
///////////////////////////////////////////////////////////////////////////////////////////

class MultiObjectiveResult
{
public:
    using Point = ObjectivePoint;
    
private:
    std::vector<Point> _paretoFront;      ///< Final non-dominated solutions
    int _functionEvaluations = 0;          ///< Total function evaluations
    int _generations = 0;                  ///< Number of generations/iterations
    Real _hypervolume = 0;                 ///< Hypervolume indicator (if computed)
    bool _converged = false;               ///< Whether algorithm converged

public:
    //---------------------------------------------------------------------------------
    // Constructors
    //---------------------------------------------------------------------------------

    MultiObjectiveResult() = default;

    explicit MultiObjectiveResult(const std::vector<Point>& paretoFront)
        : _paretoFront(paretoFront) {}

    MultiObjectiveResult(const ParetoArchive& archive)
        : _paretoFront(archive.GetSolutions()) {}

    //---------------------------------------------------------------------------------
    // Accessors
    //---------------------------------------------------------------------------------

    /// @brief Get Pareto front solutions
    const std::vector<Point>& GetParetoFront() const { return _paretoFront; }
    std::vector<Point>& GetParetoFront() { return _paretoFront; }

    /// @brief Get number of Pareto-optimal solutions
    int GetParetoFrontSize() const { return static_cast<int>(_paretoFront.size()); }

    /// @brief Get function evaluation count
    int GetFunctionEvaluations() const { return _functionEvaluations; }
    void SetFunctionEvaluations(int evals) { _functionEvaluations = evals; }

    /// @brief Get generation count
    int GetGenerations() const { return _generations; }
    void SetGenerations(int gens) { _generations = gens; }

    /// @brief Get hypervolume indicator
    Real GetHypervolume() const { return _hypervolume; }
    void SetHypervolume(Real hv) { _hypervolume = hv; }

    /// @brief Check if converged
    bool IsConverged() const { return _converged; }
    void SetConverged(bool converged) { _converged = converged; }

    //---------------------------------------------------------------------------------
    // Utility functions
    //---------------------------------------------------------------------------------

    /// @brief Get number of objectives (inferred from the front; 0 if empty)
    int GetNumObjectives() const
    {
        return _paretoFront.empty() ? 0 : _paretoFront.front().GetNumObjectives();
    }

    /// @brief Get ideal point (minimum value for each objective)
    Vector<Real> GetIdealPoint() const
    {
        int M = GetNumObjectives();
        Vector<Real> ideal(M);
        for (int m = 0; m < M; ++m) {
            ideal[m] = std::numeric_limits<Real>::infinity();
            for (const auto& pt : _paretoFront) {
                ideal[m] = std::min(ideal[m], pt.GetObjective(m));
            }
        }
        return ideal;
    }

    /// @brief Get nadir point (maximum value for each objective on Pareto front)
    Vector<Real> GetNadirPoint() const
    {
        int M = GetNumObjectives();
        Vector<Real> nadir(M);
        for (int m = 0; m < M; ++m) {
            nadir[m] = -std::numeric_limits<Real>::infinity();
            for (const auto& pt : _paretoFront) {
                nadir[m] = std::max(nadir[m], pt.GetObjective(m));
            }
        }
        return nadir;
    }

    /// @brief Get objective range for each objective
    Vector<Real> GetObjectiveRange() const
    {
        auto ideal = GetIdealPoint();
        auto nadir = GetNadirPoint();
        int M = GetNumObjectives();
        Vector<Real> range(M);
        for (int m = 0; m < M; ++m) {
            range[m] = nadir[m] - ideal[m];
        }
        return range;
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// Hypervolume Indicator (2D only for now)
//
// Measures the volume of objective space dominated by the Pareto front
// and bounded by a reference point. Higher is better.
//
// TODO: Implement general M-dimensional hypervolume (WFG algorithm)
///////////////////////////////////////////////////////////////////////////////////////////

namespace detail {
MML_OPTIMIZATION_API Real CalculateHypervolume2D(std::vector<std::pair<Real, Real>>& objectives,
                                                 const VectorN<Real, 2>& referencePoint);
MML_OPTIMIZATION_API Real CalculateMOEADHypervolume2D(std::vector<std::pair<Real, Real>>& objectives,
                                                      const VectorN<Real, 2>& referencePoint);
}

/// @brief Calculate hypervolume indicator for 2D problems
/// @param front Pareto front solutions
/// @param referencePoint Reference point (should dominate all front solutions)
/// @return Hypervolume value
inline Real CalculateHypervolume2D(const std::vector<ObjectivePoint>& front,
                            const VectorN<Real, 2>& referencePoint)
{
    std::vector<std::pair<Real, Real>> objectives;
    objectives.reserve(front.size());
    for (const auto& pt : front) {
        objectives.emplace_back(pt.GetObjective(0), pt.GetObjective(1));
    }
    return detail::CalculateHypervolume2D(objectives, referencePoint);
}

} // namespace MML::Optimization
#endif // MML_MULTI_OBJECTIVE_H
