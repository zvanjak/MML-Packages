///////////////////////////////////////////////////////////////////////////////////////////
// NSGA2.h
//
// NSGA-II (Non-dominated Sorting Genetic Algorithm II) implementation for multi-objective
// optimization. Based on Deb et al. (2002) "A fast and elitist multiobjective genetic
// algorithm: NSGA-II".
//
// Key Features:
// - Fast non-dominated sorting for population ranking
// - Crowding distance for diversity preservation
// - Binary tournament selection with crowded comparison
// - SBX crossover and polynomial mutation
// - Elitism through population merging
//
// Usage:
// @code
//     NSGA2 nsga;
//     nsga.SetProblem(problemSpec);           // ProblemSpec (runtime dimension)
//     auto result = nsga.Optimize(zdt1Problem, 250);  // 250 generations
// @endcode
//
// References:
// - Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. A. M. T. (2002).
//   "A fast and elitist multiobjective genetic algorithm: NSGA-II"
//   IEEE Transactions on Evolutionary Computation, 6(2), 182-197.
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_NSGA2_H
#define MML_NSGA2_H

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

namespace detail {
MML_OPTIMIZATION_API void NSGA2SBXCrossover(std::span<Real> first, std::span<Real> second,
                                             const ProblemSpec& spec, std::mt19937& rng, Real eta);
MML_OPTIMIZATION_API void NSGA2PolynomialMutation(std::span<Real> values, const ProblemSpec& spec,
                                                  std::mt19937& rng, Real eta, Real mutationProb);
}

///////////////////////////////////////////////////////////////////////////////////////////
// NSGA2Config - Configuration for NSGA-II algorithm
///////////////////////////////////////////////////////////////////////////////////////////

struct NSGA2Config
{
    int populationSize = 100;        ///< Population size (should be even)
    int maxGenerations = 250;        ///< Maximum generations
    Real crossoverProb = 0.9;        ///< Crossover probability
    Real mutationProb = -1.0;        ///< Mutation probability (-1 = 1/n)
    Real crossoverEta = 15.0;        ///< SBX distribution index
    Real mutationEta = 20.0;         ///< Polynomial mutation distribution index
    unsigned int seed = 0;           ///< Random seed (0 = time-based)
    bool verbose = false;            ///< Print progress
    
    /// @brief Get effective mutation probability (1/n if -1)
    Real GetMutationProb(int n) const {
        return mutationProb < 0 ? (1.0 / static_cast<Real>(n)) : mutationProb;
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// NSGA2Result - Result container for NSGA-II optimization
///////////////////////////////////////////////////////////////////////////////////////////

struct NSGA2Result
{
    using Point = ObjectivePoint;
    
    std::vector<Point> paretoFront;      ///< Final Pareto front (rank 0)
    int generations = 0;                  ///< Generations completed
    int functionEvaluations = 0;          ///< Total objective evaluations
    
    /// @brief Get number of objectives (inferred from the front; 0 if empty)
    int GetNumObjectives() const {
        return paretoFront.empty() ? 0 : paretoFront.front().GetNumObjectives();
    }

    /// @brief Get ideal point (minimum in each objective)
    Vector<Real> GetIdealPoint() const {
        int M = GetNumObjectives();
        Vector<Real> ideal(M);
        for (int m = 0; m < M; ++m)
            ideal[m] = std::numeric_limits<Real>::infinity();
        
        for (const auto& point : paretoFront) {
            for (int m = 0; m < M; ++m) {
                ideal[m] = std::min(ideal[m], point.GetObjective(m));
            }
        }
        return ideal;
    }
    
    /// @brief Get nadir point (maximum in each objective among front)
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
    
    /// @brief Calculate hypervolume indicator for 2D problems
    /// @param referencePoint Reference point (should be dominated by all front points)
    Real GetHypervolume2D(const VectorN<Real, 2>& referencePoint) const {
        return CalculateHypervolume2D(paretoFront, referencePoint);
    }
};


///////////////////////////////////////////////////////////////////////////////////////////
// NSGA2 - NSGA-II Multi-Objective Optimizer
//
// The number of decision variables is taken from the ProblemSpec (runtime) and the
// number of objectives is declared by typed problems or inferred from the first candidate.
///////////////////////////////////////////////////////////////////////////////////////////

class NSGA2 : public IIterativeAlgorithm
{
public:
    using Point = ObjectivePoint;
    using DecisionVec = Vector<Real>;
    using ObjectiveVec = Vector<Real>;
    using Result = NSGA2Result;

private:
    NSGA2Config _config;
    ProblemSpec _problemSpec;
    std::mt19937 _rng;
    int _funcEvals = 0;
    int _numObjectives = 0;
    int _declaredObjectives = 0;
    int _generations = 0;
    int _maxEvaluations = 0;
    double _maxSeconds = 0.0;
    std::function<ObjectiveVec(const DecisionVec&)> _evaluate;
    std::vector<Point> _population;
    std::vector<std::vector<int>> _fronts;
    std::vector<Point> _frontSnapshot;
    bool _observeFront = false;
    ParetoIterationProgress _progress;
    std::unique_ptr<IterativeRunner<ParetoIterationProgress>> _runner;
    std::function<std::optional<IterationStop>(const ParetoIterationProgress&)> _criterion;
    std::function<void()> _resetCriterion;
    std::vector<IterativeRunner<ParetoIterationProgress>::Observer> _observers;

public:
    //-------------------------------------------------------------------------------------
    // Constructors
    //-------------------------------------------------------------------------------------

    /// @brief Default constructor
    NSGA2() = default;

    /// @brief Construct with configuration
    explicit NSGA2(const NSGA2Config& config)
        : _config(config) {}

    //-------------------------------------------------------------------------------------
    // Configuration
    //-------------------------------------------------------------------------------------

    /// @brief Set algorithm configuration
    void SetConfig(const NSGA2Config& config) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active NSGA-II run");
        _config = config;
    }

    /// @brief Get current configuration
    const NSGA2Config& GetConfig() const { return _config; }

    /// @brief Set problem specification (variable bounds)
    void SetProblem(const ProblemSpec& spec) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active NSGA-II run");
        if (spec.Dimension() == 0) {
            throw std::invalid_argument("Problem specification has zero dimension");
        }
        for (const auto& variable : spec.Variables()) {
            if (variable.type == OptVariableType::Binary)
                throw std::invalid_argument("NSGA-II does not support binary variables with real-coded operators");
        }
        _problemSpec = spec;
    }

    void ObserveFront(bool enabled = true) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active NSGA-II run");
        _observeFront = enabled;
    }

    void SetCriterion(std::function<std::optional<IterationStop>(const ParetoIterationProgress&)> criterion,
                      std::function<void()> reset = {}) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active NSGA-II run");
        _criterion = std::move(criterion);
        _resetCriterion = std::move(reset);
    }

    void AddObserver(IterativeRunner<ParetoIterationProgress>::Observer observer) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot change an active NSGA-II run");
        _observers.push_back(std::move(observer));
    }

    void Start(const IMultiObjectiveProblem& problem, int generations = -1, int maxEvaluations = 0, double maxSeconds = 0.0) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot rebind an active NSGA-II run");
        const ProblemSpec& spec = problem.GetProblemSpec();
        if (problem.Dimension() <= 0 || problem.Dimension() != spec.Dimension())
            throw OptimizationProblemError("NSGA-II problem dimension mismatch");
        if (problem.GetNumObjectives() < 2)
            throw OptimizationProblemError("NSGA-II requires at least two objectives");
        if (problem.GetNumConstraints() != 0)
            throw OptimizationProblemError("NSGA-II does not support constrained problems");
        SetProblem(spec);
        _declaredObjectives = problem.GetNumObjectives();
        _evaluate = [&problem](const DecisionVec& genes) { return problem.EvaluateValidated(genes); };
        _generations = generations > 0 ? generations : _config.maxGenerations;
        _maxEvaluations = maxEvaluations;
        _maxSeconds = maxSeconds;
        Initialize();
    }

    template<typename ProblemType>
        requires (!std::derived_from<std::remove_cvref_t<ProblemType>, IMultiObjectiveProblem>)
    void Start(ProblemType& problem, int generations = -1, int maxEvaluations = 0, double maxSeconds = 0.0) {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot rebind an active NSGA-II run");
        if (_problemSpec.Dimension() == 0)
            throw std::runtime_error("Problem specification not set. Call SetProblem() first.");
        if constexpr (requires { problem.GetNumObjectives(); })
            _declaredObjectives = problem.GetNumObjectives();
        else
            _declaredObjectives = 0;
        if (_declaredObjectives != 0 && _declaredObjectives < 2)
            throw OptimizationProblemError("NSGA-II requires at least two objectives");
        _evaluate = [&problem](const DecisionVec& genes) { return problem.Evaluate(genes); };
        _generations = generations > 0 ? generations : _config.maxGenerations;
        _maxEvaluations = maxEvaluations;
        _maxSeconds = maxSeconds;
        Initialize();
    }

    void Initialize() override {
        if (_progress.status == IterativeRunStatus::Running)
            throw std::logic_error("Cannot initialize an active NSGA-II run");
        if (!_evaluate || _problemSpec.Dimension() == 0)
            throw std::logic_error("NSGA-II has no bound problem");
        if (_config.populationSize <= 0 || (_maxEvaluations > 0 && _maxEvaluations < _config.populationSize))
            throw std::invalid_argument("NSGA-II initial population exceeds evaluation limit");
        if (_generations <= 0 || _maxEvaluations < 0 || !std::isfinite(_maxSeconds) || _maxSeconds < 0.0)
            throw std::invalid_argument("Invalid NSGA-II run limits");
        try {
            _runner.reset();
            _progress = {};
            _frontSnapshot.clear();
            _fronts.clear();
            unsigned int seed = _config.seed;
            if (seed == 0)
                seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
            _rng.seed(seed);
            _funcEvals = 0;
            _numObjectives = _declaredObjectives;
            _population = InitializePopulation(_config.populationSize);
            EvaluatePopulation(_population);
            _fronts = FastNonDominatedSort(_population);
            for (const auto& front : _fronts)
                CalculateCrowdingDistance(_population, front);
            UpdateProgress();
            _progress.funcEvals = _funcEvals;
            _runner = std::make_unique<IterativeRunner<ParetoIterationProgress>>(
                _progress, _generations, _maxEvaluations, _maxSeconds);
            _runner->SetCriterion(_criterion, _resetCriterion);
            for (const auto& observer : _observers)
                _runner->AddObserver(observer);
            _runner->Initialize();
        } catch (...) {
            _progress.status = IterativeRunStatus::Failed;
            throw;
        }
    }

    void Step() override {
        if (!_runner)
            throw std::logic_error("NSGA-II has no active run");
        _runner->Step(_config.populationSize, [this] {
            if (_config.verbose && _progress.iterations % 10 == 0)
                std::cout << "Generation " << _progress.iterations << ": " << _fronts[0].size() << " Pareto-optimal solutions\n";
            auto offspring = CreateOffspring(_population, _config.crossoverProb,
                                             _config.GetMutationProb(_problemSpec.Dimension()));
            EvaluatePopulation(offspring);
            std::vector<Point> combined;
            combined.reserve(2 * _config.populationSize);
            combined.insert(combined.end(), _population.begin(), _population.end());
            combined.insert(combined.end(), offspring.begin(), offspring.end());
            auto fronts = FastNonDominatedSort(combined);
            _population = SelectNextGeneration(combined, fronts, _config.populationSize);
            _fronts = FastNonDominatedSort(_population);
            for (const auto& front : _fronts)
                CalculateCrowdingDistance(_population, front);
            UpdateProgress();
        });
    }

    const ParetoIterationProgress& Progress() const override {
        if (_progress.status == IterativeRunStatus::Uninitialized)
            throw std::logic_error("NSGA-II has not been initialized");
        return _progress;
    }

    void Run() override {
        if (!_runner)
            throw std::logic_error("NSGA-II has no initialized run");
        while (_progress.status == IterativeRunStatus::Running)
            Step();
    }

    //-------------------------------------------------------------------------------------
    // Main Optimization
    //-------------------------------------------------------------------------------------

    Result Optimize(const IMultiObjectiveProblem& problem, int generations = -1) {
        Start(problem, generations);
        Run();
        return GetResult();
    }

    /// @brief Run NSGA-II optimization
    /// @tparam ProblemType Multi-objective problem with Evaluate(const Vector<Real>&) -> Vector<Real>
    /// @param problem The multi-objective problem to solve
    /// @param generations Number of generations (optional, overrides config)
    template<typename ProblemType>
        requires (!std::derived_from<std::remove_cvref_t<ProblemType>, IMultiObjectiveProblem>)
    Result Optimize(ProblemType& problem, int generations = -1) {
        Start(problem, generations);
        Run();
        return GetResult();
    }

    Result GetResult() const {
        const auto& progress = Progress();
        Result result;
        result.generations = progress.iterations;
        result.functionEvaluations = progress.funcEvals;
        for (const auto& point : _population)
            if (point.GetRank() == 0)
                result.paretoFront.push_back(point);
        return result;
    }

private:
    void UpdateProgress() {
        _progress.numObjectives = _numObjectives;
        _progress.population = std::span<const Point>(_population.data(), _population.size());
        _progress.paretoFrontSize = static_cast<int>(_fronts[0].size());
        _frontSnapshot.clear();
        if (_observeFront) {
            _frontSnapshot.reserve(_fronts[0].size());
            for (const auto& point : _population)
                if (point.GetRank() == 0)
                    _frontSnapshot.push_back(point);
        }
        _progress.paretoFront = std::span<const Point>(_frontSnapshot.data(), _frontSnapshot.size());
    }
    //-------------------------------------------------------------------------------------
    // Population Initialization
    //-------------------------------------------------------------------------------------

    std::vector<Point> InitializePopulation(int size) {
        std::vector<Point> population;
        population.reserve(size);
        
        const int dim = _problemSpec.Dimension();
        for (int i = 0; i < size; ++i) {
            DecisionVec x(dim);
            for (int j = 0; j < dim; ++j) {
                const auto& var = _problemSpec[j];
                std::uniform_real_distribution<Real> dist(var.lowerBound, var.upperBound);
                x[j] = dist(_rng);
                
                // Repair for integer/categorical variables
                if (var.IsDiscrete()) {
                    x[j] = var.Repair(x[j]);
                }
            }
            population.emplace_back(x);
        }
        
        return population;
    }

    //-------------------------------------------------------------------------------------
    // Population Evaluation
    //-------------------------------------------------------------------------------------

    void EvaluatePopulation(std::vector<Point>& population) {
        for (auto& point : population) {
            auto objectives = _evaluate(point.GetX());
            if (_numObjectives == 0)
                _numObjectives = static_cast<int>(objectives.size());
            if (_numObjectives < 2 || objectives.size() != static_cast<size_t>(_numObjectives))
                throw OptimizationProblemError("NSGA-II objective result count mismatch");
            point.SetObjectives(objectives);
            ++_funcEvals;
        }
    }

    //-------------------------------------------------------------------------------------
    // Binary Tournament Selection (with crowded comparison)
    //-------------------------------------------------------------------------------------

    int BinaryTournament(const std::vector<Point>& population) {
        std::uniform_int_distribution<int> dist(0, static_cast<int>(population.size()) - 1);
        int a = dist(_rng);
        int b = dist(_rng);
        
        // Ensure different individuals
        while (b == a && population.size() > 1) {
            b = dist(_rng);
        }
        
        // Use crowded comparison: prefer lower rank, then higher crowding distance
        if (CrowdedComparison(population[a], population[b])) {
            return a;
        }
        return b;
    }

    //-------------------------------------------------------------------------------------
    // Offspring Creation (Crossover + Mutation)
    //-------------------------------------------------------------------------------------

    std::vector<Point> CreateOffspring(const std::vector<Point>& population,
                                       Real crossoverProb, Real mutationProb) {
        const int popSize = static_cast<int>(population.size());
        std::vector<Point> offspring;
        offspring.reserve(popSize);
        
        std::uniform_real_distribution<Real> prob(0.0, 1.0);
        
        while (static_cast<int>(offspring.size()) < popSize) {
            // Select parents
            int p1 = BinaryTournament(population);
            int p2 = BinaryTournament(population);
            
            DecisionVec child1 = population[p1].GetX();
            DecisionVec child2 = population[p2].GetX();
            
            // Crossover
            if (prob(_rng) < crossoverProb) {
                SBXCrossover(child1, child2);
            }
            
            // Mutation
            PolynomialMutation(child1, mutationProb);
            PolynomialMutation(child2, mutationProb);
            
            offspring.emplace_back(child1);
            if (static_cast<int>(offspring.size()) < popSize) {
                offspring.emplace_back(child2);
            }
        }
        
        return offspring;
    }

    //-------------------------------------------------------------------------------------
    // SBX Crossover (Simulated Binary Crossover)
    //-------------------------------------------------------------------------------------

    void SBXCrossover(DecisionVec& c1, DecisionVec& c2) {
        detail::NSGA2SBXCrossover(std::span<Real>(&c1[0], c1.size()),
                                  std::span<Real>(&c2[0], c2.size()),
                                  _problemSpec, _rng, _config.crossoverEta);
    }

    //-------------------------------------------------------------------------------------
    // Polynomial Mutation
    //-------------------------------------------------------------------------------------

    void PolynomialMutation(DecisionVec& x, Real mutationProb) {
        detail::NSGA2PolynomialMutation(std::span<Real>(&x[0], x.size()), _problemSpec, _rng,
                                        _config.mutationEta, mutationProb);
    }

    //-------------------------------------------------------------------------------------
    // Environmental Selection (Next Generation)
    //-------------------------------------------------------------------------------------

    std::vector<Point> SelectNextGeneration(std::vector<Point>& combined,
                                            std::vector<std::vector<int>>& fronts,
                                            int targetSize) {
        std::vector<Point> newPopulation;
        newPopulation.reserve(targetSize);
        
        // Add fronts until we reach or exceed target size
        int frontIdx = 0;
        while (frontIdx < static_cast<int>(fronts.size())) {
            const auto& front = fronts[frontIdx];
            
            if (static_cast<int>(newPopulation.size()) + static_cast<int>(front.size()) <= targetSize) {
                // Add entire front
                for (int idx : front) {
                    newPopulation.push_back(combined[idx]);
                }
                ++frontIdx;
            } else {
                // Front doesn't fit entirely - use crowding distance selection
                CalculateCrowdingDistance(combined, front);
                
                // Sort this front by crowding distance (descending)
                std::vector<int> sortedFront = front;
                std::sort(sortedFront.begin(), sortedFront.end(),
                    [&combined](int a, int b) {
                        return combined[a].GetCrowdingDistance() > combined[b].GetCrowdingDistance();
                    });
                
                // Add the most crowded solutions
                int remaining = targetSize - static_cast<int>(newPopulation.size());
                for (int i = 0; i < remaining; ++i) {
                    newPopulation.push_back(combined[sortedFront[i]]);
                }
                break;
            }
        }
        
        return newPopulation;
    }
};


} // namespace MML::Optimization
#endif // MML_NSGA2_H
