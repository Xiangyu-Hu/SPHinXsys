#ifndef ENERGY_H
#define ENERGY_H

#include "complex_solid.h"
#include "elastic_dynamics.h"

namespace SPH
{
 /**
  * @class KineticEnergy
  * @brief Compute the kinematic energy
  */
	class KineticEnergy
		: public LocalDynamicsReduce<Real, ReduceSum<Real>>,
		public GeneralDataDelegateSimple
	{
	protected:
		StdLargeVec<Real>& mass_;
		StdLargeVec<Vecd>& vel_;

	public:
		KineticEnergy(SPHBody& sph_body);
		virtual ~KineticEnergy() {};

		Real reduce(size_t index_i, Real dt = 0.0);
	};

	/**
	* @class PotentialEnergy
	* @brief Compute the Potential energy
	*/
	class PotentialEnergy
		: public LocalDynamicsReduce<Real, ReduceSum<Real>>,
		public GeneralDataDelegateSimple
	{
	private:
		SharedPtrKeeper<Gravity> gravity_ptr_keeper_;

	protected:
		StdLargeVec<Real>& mass_;
		StdLargeVec<Vecd>& pos_;

	public:
		PotentialEnergy(SPHBody& sph_body);
		virtual ~PotentialEnergy() {};

		Real reduce(size_t index_i, Real dt = 0.0);
	};
} // namespace SPH
#endif // ENERGY_H
