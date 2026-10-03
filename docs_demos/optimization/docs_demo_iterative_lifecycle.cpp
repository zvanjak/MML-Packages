#include "GeneticAlgorithm.h"
#include "MOEAD.h"
#include "NSGA2.h"
#include "ProblemTypes.h"
#include "SimulatedAnnealing.h"

#include <array>
#include <iostream>

using namespace MML;
using namespace MML::Optimization;

void Docs_Demo_IterativeLifecycle()
{
    auto spec = ProblemSpec::Continuous(2, 0.0, 1.0);
    SingleObjectiveProblem scalar([](const Vector<Real>& x) { return x[0] * x[0] + x[1] * x[1]; }, spec);
    struct BiObjective {
        Vector<Real> Evaluate(const Vector<Real>& x) const { return Vector<Real>{x[0], x[1]}; }
    } biObjective;

    GAConfig gaConfig;
    gaConfig.populationSize = 6;
    gaConfig.maxGenerations = 2;
    gaConfig.seed = 42;
    GeneticAlgorithm ga(gaConfig);
    SimulatedAnnealing sa(10.0, 1e-8, 3, 3, SimulatedAnnealing::StopCriteria::Combined, 42);
    NSGA2Config nsgaConfig;
    nsgaConfig.populationSize = 6;
    nsgaConfig.maxGenerations = 2;
    nsgaConfig.seed = 42;
    NSGA2 nsga(nsgaConfig);
    MOEADConfig moeadConfig;
    moeadConfig.populationSize = 6;
    moeadConfig.neighborhoodSize = 3;
    moeadConfig.maxGenerations = 2;
    moeadConfig.seed = 42;
    MOEAD moead(moeadConfig);

    nsga.SetProblem(spec);
    moead.SetProblem(spec);
    OptimizationConfig runConfig;
    runConfig.WithMaxIterations(2).WithTrajectory();
    ga.Start(scalar, true, &runConfig);
    sa.Start(scalar, Vector<Real>{0.5, 0.5});
    nsga.Start(biObjective);
    moead.Start(biObjective);
    std::array<IIterativeAlgorithm*, 4> algorithms{&ga, &sa, &nsga, &moead};
    for (auto* algorithm : algorithms) {
        algorithm->Step();
        algorithm->Run();
        std::cout << "Completed " << algorithm->Progress().iterations
                  << " steps and " << algorithm->Progress().funcEvals << " evaluations\n";
    }
}