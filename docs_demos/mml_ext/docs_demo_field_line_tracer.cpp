///////////////////////////////////////////////////////////////////////////////////////////
// docs_demo_field_line_tracer.cpp - FieldLineTracer package demonstration
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml_ext/algorithms/Analyzers/FieldLineTracer.h>

#include <mml/tools/serializer/SerializerFieldLines.h>

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

class MovingChargeProjectionField2D : public MML::IVectorFunction<2> {
public:
    explicit MovingChargeProjectionField2D(Real beta)
        : _beta(beta)
    {
    }

    MML::VectorN<Real, 2> operator()(const MML::VectorN<Real, 2>& point) const override
    {
        const Real x = point[0];
        const Real y = point[1];
        const Real gamma = 1.0 / std::sqrt(1.0 - _beta * _beta);
        const Real scaledRadiusSquared = gamma * gamma * x * x + y * y;

        if (scaledRadiusSquared < 1e-12) {
            return MML::VectorN<Real, 2>{0.0, 0.0};
        }

        const Real radiusFactor = std::pow(scaledRadiusSquared, -1.5);
        return MML::VectorN<Real, 2>{x * radiusFactor, gamma * y * radiusFactor};
    }

private:
    Real _beta;
};

template<int N>
std::vector<std::vector<MML::VectorN<Real, N>>> ToPointLines(const std::vector<MML::FieldLine<N>>& fieldLines)
{
    std::vector<std::vector<MML::VectorN<Real, N>>> pointLines;
    pointLines.reserve(fieldLines.size());

    for (const auto& line : fieldLines) {
        if (!line.empty()) {
            pointLines.push_back(line.points);
        }
    }

    return pointLines;
}

} // namespace

void Docs_Demo_FieldLineTracer_MovingCharge()
{
    const std::string outputDir = "field_lines_output";
    std::filesystem::create_directories(outputDir);

    MML::FieldLineTracer<2>::Config config;
    config.stepSize = 0.05;
    config.maxLength = 5.0;
    config.maxPoints = 500;
    config.minFieldMagnitude = 1e-6;
    config.traceForward = true;
    config.traceBackward = false;

    const auto bounds = MML::BoundingBox<2>::Create2D(-5.0, 5.0, -5.0, 5.0);

    constexpr int numLines = 24;
    constexpr Real seedRadius = 0.1;
    std::vector<MML::VectorN<Real, 2>> seedPoints;
    seedPoints.reserve(numLines);

    for (int i = 0; i < numLines; ++i) {
        const Real angle = 2.0 * MML::Constants::PI * static_cast<Real>(i) / static_cast<Real>(numLines);
        seedPoints.push_back(MML::VectorN<Real, 2>{
            seedRadius * std::cos(angle),
            seedRadius * std::sin(angle)
        });
    }

    std::cout << "FieldLineTracer moving-charge projection demo\n";
    std::cout << "      beta   gamma   lines   output\n";

    for (const Real beta : {0.0, 0.5, 0.9, 0.99}) {
        MovingChargeProjectionField2D field(beta);
        MML::FieldLineTracer<2> tracer(config);
        const auto fieldLines = tracer.TraceFromSeedPoints(field, seedPoints, bounds);
        const auto pointLines = ToPointLines(fieldLines);

        const std::string fileName = outputDir + "/Pancake_beta_" + std::to_string(static_cast<int>(beta * 100)) + ".mml";
        const auto result = MML::Serializer::SaveFieldLines2D(
            pointLines,
            "FieldLineTracer moving-charge projection, beta = " + std::to_string(beta),
            fileName);

        const Real gamma = beta < 1e-10 ? 1.0 : 1.0 / std::sqrt(1.0 - beta * beta);
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(10) << beta << " "
                  << std::setw(7) << std::setprecision(3) << gamma << " "
                  << std::setw(7) << fieldLines.size() << "   "
                  << (result.success ? fileName : result.message) << "\n";
    }
}
