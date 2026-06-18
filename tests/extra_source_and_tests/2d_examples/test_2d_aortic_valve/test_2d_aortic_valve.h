#include "sphinxsys.h"
using namespace SPH;

constexpr Real mm_to_m = 1e-3;
constexpr Real L_to_m3 = 1e-3;
constexpr Real min_to_s = 60.0;
constexpr Real degree_to_rad = Pi / 180.0;

// geometry
const Real height = 20 * mm_to_m;
const Real length = 6 * height;
const Real width = 117 * mm_to_m;
const Real radius = 20 * mm_to_m;
const Real shell_length = 26 * mm_to_m;
const Real shell_thickness = 0.16 * mm_to_m;
const Real angle = 45 * degree_to_rad;

const Real dp_fluid = height / 20.0;
const Real dp_solid = dp_fluid;
const Real dp_shell = dp_fluid;

const Real buffer_length = dp_fluid * 3.0;
const Real wall_thickness = dp_fluid * 4.0;

const Vec2d circle_center(3.5 * height, height);

// material
// fluid
const Real rho0_f = 1000;                      /**< Density. */
const Real Q_max = 23.44 * L_to_m3 / min_to_s; /**< Maximum flow rate. */
const Real U_f = Q_max / (width * height);     /**< Characteristic velocity. */
const Real U_max = 1.5 * U_f;                  /**< Maximum velocity. */
const Real c_f = 10.0 * U_max;                 /**< Speed of sound. */
const Real mu_f = 4.3e-3;                      /**< Dynamics viscosity. */

// solid
const Real rho0_s = rho0_f;        /**< Reference density.*/
const Real poisson = 0.49;         /**< Poisson ratio.*/
const Real Youngs_modulus = 1.5e6; /**< Youngs modulus.*/

// Cycle
const Real time_flow_init = 0;
const Real time_cycle = 1.0;
const size_t num_cycles = 2;