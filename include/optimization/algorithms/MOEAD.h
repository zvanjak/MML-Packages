///////////////////////////////////////////////////////////////////////////////////////////
// MOEAD.h
//
// MOEA/D (Multi-Objective Evolutionary Algorithm based on Decomposition) implementation.
// Based on Zhang & Li (2007) "MOEA/D: A Multiobjective Evolutionary Algorithm Based on
// Decomposition".
//
// Key Features:
// - Decomposes multi-objective problem into scalar subproblems
// - Uses neighborhood structure for mating selection
// - Supports multiple aggregation approaches (Weighted Sum, Tchebycheff, PBI)
// - SBX crossover and polynomial mutation (shared with NSGA-II)
// - Efficient for problems with uniform Pareto fronts
//
// Usage:
// @code
//     MOEAD moead;
//     moead.SetProblem(problemSpec);          // ProblemSpec (runtime dimension)
//     auto result = moead.Optimize(zdt1Problem, 250);
// @endcode
//
// References:
// - Zhang, Q., & Li, H. (2007). "MOEA/D: A multiobjective evolutionary algorithm based
//   on decomposition". IEEE Transactions on Evolutionary Computation, 11(6), 712-731.
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_MOEAD_H
#define MML_MOEAD_H

#include "MMLBase.h"
#include "MultiObjective.h"
#include "GeneticAlgorithm.h"
#include "IOptimizationProblem.h"
#include "IterativeRunner.h"

#include <concepts>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <type_traits>

namespace MML::Optimization {

///////////////////////////////////////////////////////////////////////////////////////////
// AggregationType - Aggregation approach for scalar subproblems
///////////////////////////////////////////////////////////////////////////////////////////

enum class AggregationType
{
    WeightedSum,    ///< Simple weighted sum: g(x|w) = sum(w_i * f_i(x))
    Tchebycheff,    ///< Tchebycheff: g(x|w,z*) = max(w_i * |f_i(x) - z*_i|)
    PBI             ///< Penalty-based boundary intersection (for uniform fronts)
};

namespace detail {
MML_OPTIMIZATION_API Real MOEADScalarizedValue(std::span<const Real> objectives,
                                               std::span<const Real> weights,
                                               std::span<const Real> idealPoint,
                                               AggregationType aggregation, Real pbiTheta);
MML_OPTIMIZATION_API void MOEADSBXCrossover(std::span<Real> first, std::span<Real> second,
                                            const ProblemSpec& spec, std::mt19937& rng, Real eta);
MML_OPTIMIZATION_API void MOEADPolynomialMutation(std::span<Real> values, const ProblemSpec& spec,
                                                  std::mt19937& rng, Real eta, Real mutationProb);
MML_OPTIMIZATION_API std::vector<Real> MOEADGenerateWeights(int objectiveCount, int populationSize,
                                                           std::mt19937& rng);
MML_OPTIMIZATION_API std::vector<std::vector<int>> MOEADInitializeNeighborhoods(
    std::span<const Real> weights, int objectiveCount, int neighborhoodSize);
}


///////////////////////////////////////////////////////////////////////////////////////////
// MOEADConfig - Configuration for MOEA/D algorithm
///////////////////////////////////////////////////////////////////////////////////////////

struct MOEADConfig
{
    int populationSize = 100;        ///< Number of weight vectors / subproblems
    int maxGenerations = 250;        ///< Maximum generations
    int neighborhoodSize = 20;       ///< Number of neighbors for each subproblem
    Real neighborhoodSelectionProb = 0.9; ///< Probability of selecting from neighborhood
    AggregationType aggregation = AggregationType::Tchebycheff; ///< Aggregation function
    Real pbiTheta = 5.0;             ///< Penalty parameter for PBI approach
    Real crossoverProb = 1.0;        ///< Crossover probability
    Real mutationProb = -1.0;        ///< Mutation probability (-1 = 1/n)
    Real crossoverEta = 20.0;        ///< SBX distribution index
    Real mutationEta = 20.0;         ///< Polynomial mutation distribution index
    unsigned int seed = 0;           ///< Random seed (0 = time-based)
    bool verbose = false;            ///< Print progress
    
    /// @brief Get effective mutation probability
    Real GetMutationProb(int n) const {
        return mutationProb < 0 ? (1.0 / static_cast<Real>(n)) : mutationProb;
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// MOEADResult - Result container for MOEA/D optimization
///////////////////////////////////////////////////////////////////////////////////////////

struct MOEADResult
{
    using Point = ObjectivePoint;
    
    std::vector<Point> finalPopulation;   ///< Final population (one per subproblem)
    std::vector<Point> paretoFront;       ///< Extracted Pareto front
    int generations = 0;                   ///< Generations completed
    int functionEvaluations = 0;           ///< Total objective evaluations
    Vector<Real> idealPoint;               ///< Best found in each objective
    
    /// @brief Get ideal point (minimum in each objective found during optimization)
    const Vector<Real>& GetIdealPoint() const { return idealPoint; }
    
    /// @brief Get number of objectives (inferred from the front; 0 if empty)
    int GetNumObjectives() const {
        return paretoFront.empty() ? 0 : paretoFront.front().GetNumObjectives();
    }

    /// @brief Get nadir point (maximum in each objective among Pareto front)
    Vector<Real> GetNadirPoint() const {
        int M = GetNumObjectives();
        Vector<Real> nadir(M);
        for (int m = 0; m < M; ++m)
            nadir[m] = -std::numeric_limits<Real>::infinity();
        
        for (const auto& point : paretoFront) {
            for (int m = 0; m < M; ++m) {
                nadir[m] = std::max(nadir[m], point.GetObjective(m));
            }
        }
        return nadir;
    }
    
    /// @brief Compute hypervolume indicator for bi-objective problems
    Real GetHypervolume2D(const VectorN<Real, 2>& referencePoint) const {
        std::vector<std::pair<Real, Real>> objectives;
        objectives.reserve(paretoFront.size());
        for (const auto& pt : paretoFront) {
            objectives.emplace_back(pt.GetObjective(0), pt.GetObjective(1));
        }
        return detail::CalculateMOEADHypervolume2D(objectives, referencePoint);
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// MOEAD - Multi-Objective Evolutionary Algorithm based on Decomposition
//
// The number of decision variables is taken from the ProblemSpec (runtime) and the
// number of objectives is declared by typed problems or inferred from the first candidate.
///////////////////////////////////////////////////////////////////////////////////////////

class MOEAD : public IIterativeAlgorithm
{
public:
    using DecisionVec = Vector<Real>;
    using ObjectiveVec = Vector<Real>;
    using Point = ObjectivePoint;
    using Config = MOEADConfig;
    using Result = MOEADResult;

private:
    Config _config;
    ProblemSpec _problemSpec;
    std::mt19937 _rng;
    int _numObjectives = 0;      // Determined at runtime from the problem
    
    // Weight vectors and neighborhood
    std::vector<Vector<Real>> _weightVectors;
    std::vector<std::vector<int>> _neighborhoods;
    
    // Reference point (ideal point)
    ObjectiveVec _idealPoint;
    std::vector<Point> _population;
    std::vector<Point> _frontSnapshot;
    bool _observeFront = false;
    int _declaredObjectives = 0;
    int _generations = 0;
    bool _initialOnly = false;
    int _maxEvaluations = 0;
    double _maxSeconds = 0.0;
    std::function<ObjectiveVec(const DecisionVec&)> _evaluate;
    ParetoIterationProgress _progress;
    std::unique_ptr<IterativeRunner<ParetoIterationProgress>> _runner;
    std::function<std::optional<IterationStop>(const ParetoIterationProgress&)> _criterion;
    std::function<void()> _resetCriterion;
    std::vector<IterativeRunner<ParetoIterationProgress>::Observer> _observers;

public:
    //---------------------------------------------------------------------------------
    // Constructors
    //---------------------------------------------------------------------------------
    
    /// @brief Default constructor
    MOEAD() : _rng(std::random_device{}()) {}
    
    /// @brief Construct with configuration
    explicit MOEAD(const Config& config) : _config(config) {
        if (config.seed != 0) {
            _rng.seed(config.seed);
        } else {
            _rng.seed(static_cast<unsigned int>(
                std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        }
    }
    
    //---------------------------------------------------------------------------------
    // Configuration
    //---------------------------------------------------------------------------------
    
    /// @brief Set configuration
    void SetConfig(const Config& config) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active MOEA/D run");
        _config = config;
        if (config.seed != 0) {
            _rng.seed(config.seed);
        }
    }
    
    /// @brief Get configuration
    const Config& GetConfig() const { return _config; }
    
    /// @brief Set problem specification
    void SetProblem(const ProblemSpec& spec) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active MOEA/D run");
        if (spec.Dimension() == 0) {
            throw std::invalid_argument("Problem specification has zero dimension");
        }
        for (const auto& variable : spec.Variables()) {
            if (variable.type == OptVariableType::Binary)
                throw std::invalid_argument("MOEA/D does not support binary variables with real-coded operators");
        }
        _problemSpec = spec;
    }

    void ObserveFront(bool enabled = true) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active MOEA/D run");
        _observeFront = enabled;
    }

    void SetCriterion(std::function<std::optional<IterationStop>(const ParetoIterationProgress&)> criterion,
                      std::function<void()> reset = {}) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active MOEA/D run");
        _criterion = std::move(criterion);
        _resetCriterion = std::move(reset);
    }

    void AddObserver(IterativeRunner<ParetoIterationProgress>::Observer observer) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active MOEA/D run");
        _observers.push_back(std::move(observer));
    }

    void Start(const IMultiObjectiveProblem& problem, int generations = -1, int maxEvaluations = 0, double maxSeconds = 0.0) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot rebind an active MOEA/D run");
        const ProblemSpec& spec = problem.GetProblemSpec();
        if (problem.Dimension() <= 0 || problem.Dimension() != spec.Dimension())
            throw OptimizationProblemError("MOEA/D problem dimension mismatch");
        if (problem.GetNumObjectives() < 2)
            throw OptimizationProblemError("MOEA/D requires at least two objectives");
        if (problem.GetNumConstraints() != 0)
            throw OptimizationProblemError("MOEA/D does not support constrained problems");
        SetProblem(spec);
        _declaredObjectives = problem.GetNumObjectives();
        _evaluate = [&problem](const DecisionVec& genes) { return problem.EvaluateValidated(genes); };
        _initialOnly = generations == 0;
        _generations = _initialOnly ? 1 : generations < 0 ? _config.maxGenerations : generations;
        _maxEvaluations = maxEvaluations;
        _maxSeconds = maxSeconds;
        Initialize();
    }

    template<typename ProblemType>
        requires (!std::derived_from<std::remove_cvref_t<ProblemType>, IMultiObjectiveProblem>)
    void Start(ProblemType& problem, int generations = -1, int maxEvaluations = 0, double maxSeconds = 0.0) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot rebind an active MOEA/D run");
        if (_problemSpec.Dimension() == 0)
            throw std::runtime_error("Problem specification not set. Call SetProblem() first.");
        if constexpr (requires { problem.GetNumObjectives(); })
            _declaredObjectives = problem.GetNumObjectives();
        else
            _declaredObjectives = 0;
        if (_declaredObjectives != 0 && _declaredObjectives < 2)
            throw OptimizationProblemError("MOEA/D requires at least two objectives");
        _evaluate = [&problem](const DecisionVec& genes) { return problem.Evaluate(genes); };
        _initialOnly = generations == 0;
        _generations = _initialOnly ? 1 : generations < 0 ? _config.maxGenerations : generations;
        _maxEvaluations = maxEvaluations;
        _maxSeconds = maxSeconds;
        Initialize();
    }

    void Initialize() override {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot initialize an active MOEA/D run");
        if (!_evaluate || _problemSpec.Dimension() == 0)
            throw std::logic_error("MOEA/D has no bound problem");
        if (_config.populationSize <= 0 || (_maxEvaluations > 0 && _maxEvaluations < _config.populationSize))
            throw std::invalid_argument("MOEA/D initial population exceeds evaluation limit");
        if (_generations <= 0 || _maxEvaluations < 0 || !std::isfinite(_maxSeconds) || _maxSeconds < 0.0)
            throw std::invalid_argument("Invalid MOEA/D run limits");
        try {
            _runner.reset();
            _progress = {};
            _frontSnapshot.clear();
            _numObjectives = _declaredObjectives;
            _population = InitializePopulation();
            InitializeNeighborhoods();
            InitializeIdealPoint();
            for (const auto& point : _population)
                UpdateIdealPoint(point.GetObjectives());
            UpdateProgress();
            _progress.funcEvals = _config.populationSize;
            _runner = std::make_unique<IterativeRunner<ParetoIterationProgress>>(
                _progress, _generations, _maxEvaluations, _maxSeconds);
            _runner->SetCriterion(_criterion, _resetCriterion);
            for (const auto& observer : _observers)
                _runner->AddObserver(observer);
            _runner->Initialize();
            if (_initialOnly && _progress.status == IterativeRunStatus::Running)
                _runner->Finish(IterativeRunStatus::LimitReached, "Maximum iterations reached");
        } catch (...) {
            _progress.status = IterativeRunStatus::Failed;
            throw;
        }
    }

    void Step() override {
        if (!_runner)
            throw std::logic_error("MOEA/D has no active run");
        _runner->Step(_config.populationSize, [this] {
            for (int index = 0; index < _config.populationSize; ++index) {
                bool useNeighborhood = std::uniform_real_distribution<Real>(0, 1)(_rng)
                                       < _config.neighborhoodSelectionProb;
                const auto& matingPool = useNeighborhood ? _neighborhoods[index]
                                                         : GetWholePopulationIndices();
                int parent1 = matingPool[std::uniform_int_distribution<int>(0,
                                      static_cast<int>(matingPool.size()) - 1)(_rng)];
                int parent2 = matingPool[std::uniform_int_distribution<int>(0,
                                      static_cast<int>(matingPool.size()) - 1)(_rng)];
                DecisionVec child = _population[parent1].GetX();
                DecisionVec other = _population[parent2].GetX();
                if (std::uniform_real_distribution<Real>(0, 1)(_rng) < _config.crossoverProb)
                    SBXCrossover(child, other);
                PolynomialMutation(child, _config.GetMutationProb(_problemSpec.Dimension()));
                ObjectiveVec childObjectives = _evaluate(child);
                if (childObjectives.size() != static_cast<size_t>(_numObjectives))
                    throw OptimizationProblemError("MOEA/D objective result count mismatch");
                UpdateIdealPoint(childObjectives);
                Point childPoint(child, childObjectives);
                const auto& updatePool = useNeighborhood ? _neighborhoods[index]
                                                         : GetWholePopulationIndices();
                for (int neighbor : updatePool) {
                    Real childScalar = ComputeScalarizedValue(childObjectives, _weightVectors[neighbor]);
                    Real currentScalar = ComputeScalarizedValue(
                        _population[neighbor].GetObjectives(), _weightVectors[neighbor]);
                    if (childScalar < currentScalar)
                        _population[neighbor] = childPoint;
                }
            }
            UpdateProgress();
            if (_config.verbose && (_progress.iterations + 1) % 50 == 0)
                std::cout << "MOEA/D Generation " << (_progress.iterations + 1) << "/" << _generations << "\n";
        });
    }

    const ParetoIterationProgress& Progress() const override {
        if (_progress.status == IterativeRunStatus::Uninitialized)
            throw std::logic_error("MOEA/D has not been initialized");
        return _progress;
    }

    void Run() override {
        if (!_runner)
            throw std::logic_error("MOEA/D has no initialized run");
        while (_progress.status == IterativeRunStatus::Running)
            Step();
    }

    //---------------------------------------------------------------------------------
    // Main Optimization
    //---------------------------------------------------------------------------------

    Result Optimize(const IMultiObjectiveProblem& problem, int generations = -1)
    {
        Start(problem, generations);
        Run();
        return GetResult();
    }
    
    /// @brief Run MOEA/D optimization
    /// @tparam ProblemType Multi-objective problem with Evaluate(const Vector<Real>&) -> Vector<Real>
    /// @param problem The multi-objective problem
    /// @param generations Number of generations (-1 = use config)
    /// @return Optimization result with Pareto front
    template<typename ProblemType>
        requires (!std::derived_from<std::remove_cvref_t<ProblemType>, IMultiObjectiveProblem>)
    Result Optimize(ProblemType& problem, int generations = -1)
    {
        Start(problem, generations);
        Run();
        return GetResult();
    }

    Result GetResult() const {
        const auto& progress = Progress();
        Result result;
        result.generations = progress.iterations;
        result.functionEvaluations = progress.funcEvals;
        result.idealPoint = _idealPoint;
        result.finalPopulation = _population;
        result.paretoFront = ExtractParetoFront(_population);
        return result;
    }

private:
    void UpdateProgress() {
        _progress.numObjectives = _numObjectives;
        _progress.idealPoint = &_idealPoint;
        _progress.population = std::span<const Point>(_population.data(), _population.size());
        _frontSnapshot.clear();
        if (_observeFront)
            _frontSnapshot = ExtractParetoFront(_population);
        _progress.paretoFront = std::span<const Point>(_frontSnapshot.data(), _frontSnapshot.size());
        _progress.paretoFrontSize = static_cast<int>(_frontSnapshot.size());
    }
    //---------------------------------------------------------------------------------
    // Weight Vector Generation (Simplex Lattice Design)
    //---------------------------------------------------------------------------------
    
    void GenerateWeightVectors()
    {
        const int M = _numObjectives;
        _weightVectors.clear();
        auto weights = detail::MOEADGenerateWeights(M, _config.populationSize, _rng);
        _weightVectors.reserve(weights.size() / M);
        for (size_t index = 0; index < weights.size(); index += M) {
            Vector<Real> weight(M);
            for (int objective = 0; objective < M; ++objective)
                weight[objective] = weights[index + objective];
            _weightVectors.push_back(weight);
        }
    }
    
    //---------------------------------------------------------------------------------
    // Neighborhood Initialization
    //---------------------------------------------------------------------------------
    
    void InitializeNeighborhoods()
    {
        const int M = _numObjectives;
        std::vector<Real> weights;
        weights.reserve(_weightVectors.size() * M);
        for (const auto& weight : _weightVectors) {
            for (int objective = 0; objective < M; ++objective)
                weights.push_back(weight[objective]);
        }
        _neighborhoods = detail::MOEADInitializeNeighborhoods(weights, M, _config.neighborhoodSize);
    }
    
    //---------------------------------------------------------------------------------
    // Population Initialization
    //---------------------------------------------------------------------------------
    
    std::vector<Point> InitializePopulation()
    {
        std::vector<Point> population;
        population.reserve(_config.populationSize);
        
        const int dim = _problemSpec.Dimension();
        auto firstCandidateRng = _rng;
        for (int i = 0; i < _config.populationSize; ++i) {
            DecisionVec x(dim);
            for (int j = 0; j < dim; ++j) {
                const auto& spec = _problemSpec[j];
                Real lb = spec.lowerBound;
                Real ub = spec.upperBound;
                
                Real val = std::uniform_real_distribution<Real>(lb, ub)(i == 0 ? firstCandidateRng : _rng);
                
                if (spec.IsDiscrete()) {
					val = spec.Repair(val);
                }
                x[j] = val;
            }
            
            ObjectiveVec obj = _evaluate(x);
            if (_numObjectives == 0)
                _numObjectives = static_cast<int>(obj.size());
            if (_numObjectives < 2)
                throw OptimizationProblemError("MOEA/D requires at least two objectives");
            if (obj.size() != static_cast<size_t>(_numObjectives))
                throw OptimizationProblemError("MOEA/D objective result count mismatch");
            if (i == 0)
                GenerateWeightVectors();
            population.emplace_back(x, obj);
        }
        
        return population;
    }
    
    //---------------------------------------------------------------------------------
    // Ideal Point Management
    //---------------------------------------------------------------------------------
    
    void InitializeIdealPoint()
    {
        _idealPoint = Vector<Real>(_numObjectives);
        for (int m = 0; m < _numObjectives; ++m) {
            _idealPoint[m] = std::numeric_limits<Real>::infinity();
        }
    }
    
    void UpdateIdealPoint(const ObjectiveVec& obj)
    {
        for (int m = 0; m < _numObjectives; ++m) {
            _idealPoint[m] = std::min(_idealPoint[m], obj[m]);
        }
    }
    
    //---------------------------------------------------------------------------------
    // Scalarization Functions
    //---------------------------------------------------------------------------------
    
    Real ComputeScalarizedValue(const ObjectiveVec& obj, const Vector<Real>& weight) const
    {
        return detail::MOEADScalarizedValue(std::span<const Real>(&obj[0], obj.size()),
                                            std::span<const Real>(&weight[0], weight.size()),
                                            std::span<const Real>(&_idealPoint[0], _idealPoint.size()),
                                            _config.aggregation, _config.pbiTheta);
    }
    
    //---------------------------------------------------------------------------------
    // Variation Operators (SBX Crossover and Polynomial Mutation)
    //---------------------------------------------------------------------------------
    
    void SBXCrossover(DecisionVec& c1, DecisionVec& c2)
    {
        detail::MOEADSBXCrossover(std::span<Real>(&c1[0], c1.size()),
                                  std::span<Real>(&c2[0], c2.size()),
                                  _problemSpec, _rng, _config.crossoverEta);
    }
    
    void PolynomialMutation(DecisionVec& x, Real mutationProb)
    {
        detail::MOEADPolynomialMutation(std::span<Real>(&x[0], x.size()), _problemSpec, _rng,
                                        _config.mutationEta, mutationProb);
    }
    
    //---------------------------------------------------------------------------------
    // Helper Functions
    //---------------------------------------------------------------------------------
    
    std::vector<int> GetWholePopulationIndices() const
    {
        std::vector<int> indices(_config.populationSize);
        std::iota(indices.begin(), indices.end(), 0);
        return indices;
    }
    
    /// @brief Extract non-dominated solutions from population
    std::vector<Point> ExtractParetoFront(const std::vector<Point>& population) const
    {
        std::vector<Point> front;
        front.reserve(population.size());
        
        for (size_t i = 0; i < population.size(); ++i) {
            bool dominated = false;
            for (size_t j = 0; j < population.size(); ++j) {
                if (i != j && Dominates(population[j].GetObjectives(), 
                                        population[i].GetObjectives())) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) {
                front.push_back(population[i]);
            }
        }
        
        return front;
    }
};

} // namespace MML::Optimization
#endif // MML_MOEAD_H
