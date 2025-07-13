#include "all_simbody.h"

using namespace SimTK;

//------------------------------------------------------------------------------
//                               udf motion
//------------------------------------------------------------------------------
class UdfSinusoidMotionImplementation
	: public Motion::Custom::Implementation
{
public:
	UdfSinusoidMotionImplementation(Real amplitude, Real t0, Real rate, Real phase)
		: Motion::Custom::Implementation(), defAmplitude(amplitude), defRate(rate), defPhase(phase), defT(t0)
	{
	
	};

	UdfSinusoidMotionImplementation* clone() const override {
		UdfSinusoidMotionImplementation* copy = new UdfSinusoidMotionImplementation(*this);
		return copy;
	};

	void calcPrescribedPosition(const State &s, int nq, Real *q) const 
	{
		const Real t = s.getTime();
		Real out = defAmplitude*std::sin(defRate*t + defPhase);
		if (t < defT)
		{
			Real q = t / defT;
			out *= pow(q, 3) * (10.0 - 15.0 * q + 6.0 * q *q);
		} 
		for (int i = 0; i < nq; ++i)
			q[i] = out;
	};

	void calcPrescribedPositionDot(const State &s, int nq, Real *qdot) const 
	{
		const Real t = s.getTime();
		Real outd = defAmplitude * defRate * std::cos(defRate*t + defPhase);
		if (t < defT)
		{
			Real q = t / defT;
			outd *= pow(q, 3) * (10.0 - 15.0 * q + 6.0 * q *q);
			outd += defAmplitude * std::sin(defRate*t + defPhase) * ( pow(q, 3) * (12.0 * q - 15.0) 
				   + 3.0 * q * q * (6.0 * q * q - 15.0 * q + 10.0) ) /defT;
		}
		for (int i = 0; i < nq; ++i)
			qdot[i] = outd;
	};

	void calcPrescribedPositionDotDot(const State &s, int nq, Real *qdotdot) const 
	{
		const Real t = s.getTime();
		Real outdd = - defAmplitude * defRate * defRate * std::sin(defRate*t + defPhase);

		if (t < defT)
		{
			Real q = t / defT;
			outdd *= pow(q, 3) * (10.0 - 15.0 * q + 6.0 * q *q);
			outdd += 2.0 * defAmplitude * defRate * std::cos(defRate*t + defPhase) 
				* (pow(q, 3) * (12.0 * q - 15.0) + 3.0 * q * q * (6.0 * q * q - 15.0 * q + 10.0)) / defT;

			outdd += defAmplitude * std::sin(defRate*t + defPhase)
				* (12.0 * pow(q, 3) + 6.0 * q * q * (12.0 * q - 15.0) + 6.0 * q * (6.0 * q * q - 15.0 * q + 10.0) ) / (defT * defT);
		}

		for (int i = 0; i < nq; ++i)
			qdotdot[i] = outdd;
	};
	
	Motion::Level getLevel(const State &) const {
		return Motion::Level::Position;
	};

private:
	Real	defAmplitude, defRate, defPhase;
	Real	defT;
};


class MyMotion : public Motion::Custom {
public:
	MyMotion(MobilizedBody& mobod, Real A, Real t0, Real w, Real phi)
		: Motion::Custom(mobod, new UdfSinusoidMotionImplementation(A, t0, w, phi)) {}
};