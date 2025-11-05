#include "all_simbody.h"

//------------------------------------------------------------------------------
//                               udf motion
//------------------------------------------------------------------------------
class UdfSinusoidMotionImplementation
    : public SimTK::Motion::Custom::Implementation
{
  public:
    UdfSinusoidMotionImplementation(double amplitude, double t0, double rate, double phase)
        : SimTK::Motion::Custom::Implementation(), defAmplitude(amplitude),
          defRate(rate), defPhase(phase), defT(t0) {};

    UdfSinusoidMotionImplementation *clone() const override
    {
        UdfSinusoidMotionImplementation *copy = new UdfSinusoidMotionImplementation(*this);
        return copy;
    };

    void calcPrescribedPosition(const SimTK::State &s, int nq, double *q) const
    {
        const double t = s.getTime();
        double out = defAmplitude * std::sin(defRate * t + defPhase);
        if (t < defT)
        {
            double q = t / defT;
            out *= pow(q, 3) * (10.0 - 15.0 * q + 6.0 * q * q);
        }
        for (int i = 0; i < nq; ++i)
            q[i] = out;
    };

    void calcPrescribedPositionDot(const SimTK::State &s, int nq, double *qdot) const
    {
        const double t = s.getTime();
        double outd = defAmplitude * defRate * std::cos(defRate * t + defPhase);
        if (t < defT)
        {
            double q = t / defT;
            outd *= pow(q, 3) * (10.0 - 15.0 * q + 6.0 * q * q);
            outd += defAmplitude * std::sin(defRate * t + defPhase) *
                    (pow(q, 3) * (12.0 * q - 15.0) + 3.0 * q * q * (6.0 * q * q - 15.0 * q + 10.0)) / defT;
        }
        for (int i = 0; i < nq; ++i)
            qdot[i] = outd;
    };

    void calcPrescribedPositionDotDot(const SimTK::State &s, int nq, double *qdotdot) const
    {
        const double t = s.getTime();
        double outdd = -defAmplitude * defRate * defRate * std::sin(defRate * t + defPhase);

        if (t < defT)
        {
            double q = t / defT;
            outdd *= pow(q, 3) * (10.0 - 15.0 * q + 6.0 * q * q);
            outdd += 2.0 * defAmplitude * defRate * std::cos(defRate * t + defPhase) *
                     (pow(q, 3) * (12.0 * q - 15.0) + 3.0 * q * q * (6.0 * q * q - 15.0 * q + 10.0)) / defT;

            outdd += defAmplitude * std::sin(defRate * t + defPhase) *
                     (12.0 * pow(q, 3) + 6.0 * q * q * (12.0 * q - 15.0) + 6.0 * q * (6.0 * q * q - 15.0 * q + 10.0)) / (defT * defT);
        }

        for (int i = 0; i < nq; ++i)
            qdotdot[i] = outdd;
    };

    SimTK::Motion::Level getLevel(const SimTK::State &) const
    {
        return SimTK::Motion::Level::Position;
    };

  private:
    double defAmplitude, defRate, defPhase;
    double defT;
};

class MyMotion : public SimTK::Motion::Custom
{
  public:
    MyMotion(SimTK::MobilizedBody &mobod, double A, double t0, double w, double phi)
        : SimTK::Motion::Custom(mobod, new UdfSinusoidMotionImplementation(A, t0, w, phi)) {}
};
