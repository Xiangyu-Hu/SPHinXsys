#ifndef GRANULAR_BODY_H
#define GRANULAR_BODY_H


#include "base_body.h"
#include "base_body_part.h"
#include "fluid_body.h"

namespace SPH {
	class GranularBody : public FluidBody
	{
	public:
		explicit GranularBody(SPHSystem &system, SharedPtr<Shape> shape_ptr);
		virtual ~GranularBody() {};

		virtual GranularBody* ThisObjectPtr() override { return this; };

	};
}
#endif //GRANULAR_BODY_H