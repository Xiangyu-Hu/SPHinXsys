#ifndef NECKING_REGRESSION_H
#define NECKING_REGRESSION_H

#include "necking_case.h"
#include <memory>

namespace SPH
{
/** Standard SPHinXsys regression checks on five synchronized output snapshots. */
class NeckingRegression
{
    using EnergyRegression = RegressionTestEnsembleAverage<ReducedQuantityRecording<TotalKineticEnergy>>;
    using DamageRegression = RegressionTestEnsembleAverage<ReducedQuantityRecording<Average<QuantitySummation<Real>>>>;
    std::unique_ptr<EnergyRegression> energy_;
    std::unique_ptr<DamageRegression> damage_;

  public:
    NeckingRegression(SPHBody &body, bool enabled)
    {
        if (enabled)
        {
            energy_ = std::make_unique<EnergyRegression>(body);
            damage_ = std::make_unique<DamageRegression>(body, "JCDamage");
        }
    }

    /** The device caller synchronizes velocity and damage before recording. */
    void Record(size_t frame)
    {
        if (energy_)
        {
            energy_->writeToFile(frame);
            damage_->writeToFile(frame);
        }
    }

    void TestResult()
    {
        if (energy_)
        {
            energy_->testResult();
            damage_->testResult();
        }
    }
};
} // namespace SPH
#endif // NECKING_REGRESSION_H
