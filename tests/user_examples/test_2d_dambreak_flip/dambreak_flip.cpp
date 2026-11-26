/**
 * @file dambreak_flip.cpp
 * @brief 2D dam break in a closed tank with the FLIP method on a staggered MAC grid (SPHinXsys CK API).
 *        First working version (work in progress).
 *
 * Method overview (one time step):
 *   1. Classify the grid cells as Air / Fluid / Solid from the particle positions.
 *   2. Particles -> grid (P2G): transfer the particle velocities to the MAC grid faces
 *      with bilinear weights, and keep a copy of the resulting face velocities.
 *   3. Add gravity on the grid and impose the wall (no-penetration) condition.
 *   4. Pressure projection: solve the pressure Poisson equation on the fluid cells with a
 *      matrix-free conjugate gradient (Neumann condition at walls, p = 0 in air),
 *      then subtract the pressure gradient from the face velocities.
 *   5. Grid -> particles (G2P): update each particle velocity with a FLIP/PIC blend,
 *      v = alpha * (v_particle + dv_grid) + (1 - alpha) * v_grid, with alpha = 0.95.
 *   6. Move the particles (explicit Euler) and keep them inside the tank.
 *
 * Layout of this file:
 *   - MacGrid        : data holder and serial host methods of the staggered grid.
 *   - *CK classes    : particle dynamics written in the CK style (kernel executed per particle).
 *   - main()         : bodies, dynamics registration and the FLIP time loop.
 *
 * Not included yet: velocity extrapolation into air, particle reseeding, RK2 advection,
 * a preconditioner for the CG solver, and a parallel (CK) version of the grid operations.
 */
#include "sphinxsys.h"

using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;                    /**< Water tank length. */
Real DH = 5.366;                    /**< Water tank height. */
Real LL = 2.0;                      /**< Water column length. */
Real LH = 1.0;                      /**< Water column height. */
Real particle_spacing_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Reference density of fluid. */
Real gravity_g = 1.0;                    /**< Gravity. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Artificial sound speed (only needed to define the fluid material). */
//----------------------------------------------------------------------
//	Parameters for the shapes.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;

//----------------------------------------------------------------------
//	Staggered MAC grid (Harlow & Welch, 1965).
//	- Pressure p lives at the cell centres.
//	- The horizontal velocity u lives on the vertical faces (between left and right cells).
//	- The vertical velocity v lives on the horizontal faces (between lower and upper cells).
//	The grid covers the tank plus one layer of solid cells around it.
//	Naming: nx(), ny() = number of cells inside the tank;
//	        Nx(), Ny() = total number of cells including the solid layer (nx + 3, ny + 3).
//----------------------------------------------------------------------
class MacGrid
{
  public:
    enum CellType : int { Air = 0, Fluid = 1, Solid = 2 };

    MacGrid(Real length, Real height, Real target_dx)
        : nx_(static_cast<int>(std::ceil(length / target_dx))),
          dx_(length / nx_),
          ny_(static_cast<int>(std::lround(height / dx_))),
          mesh_(Vecd(-dx_, -dx_), dx_, Arrayi(nx_ + 3, ny_ + 3)),
          p_("Pressure", mesh_.NumberOfCells()),
          u_("VelocityX", (mesh_.AllCells()[0] + 1) * mesh_.AllCells()[1]),
          v_("VelocityY", mesh_.AllCells()[0] * (mesh_.AllCells()[1] + 1)),
          cell_type_("CellType", mesh_.NumberOfCells(), static_cast<int>(Air)),
          u_weight_("WeightU", (mesh_.AllCells()[0] + 1) * mesh_.AllCells()[1]),
          v_weight_("WeightV", mesh_.AllCells()[0] * (mesh_.AllCells()[1] + 1)),
          u_old_("VelocityXOld", (mesh_.AllCells()[0] + 1) * mesh_.AllCells()[1]),
          v_old_("VelocityYOld", mesh_.AllCells()[0] * (mesh_.AllCells()[1] + 1))
    {
        // Mark the outer layer of cells as solid (the tank walls).
        for (int i = 0; i < Nx(); ++i)
            for (int j = 0; j < Ny(); ++j)
                if (i == 0 || j == 0 || i == Nx() - 1 || j == Ny() - 1)
                    cell_type_.Data()[cellIndex(i, j)] = Solid;
    }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int Nx() const { return mesh_.AllCells()[0]; }
    int Ny() const { return mesh_.AllCells()[1]; }
    Real dx() const { return dx_; }
    const Mesh &mesh() const { return mesh_; }

    // Linear indices: cells use the Mesh convention (i * Ny + j);
    // u has Nx + 1 columns of Ny values, v has Nx columns of Ny + 1 values.
    UnsignedInt cellIndex(int i, int j) const { return mesh_.LinearCellIndex(Arrayi(i, j)); }
    UnsignedInt uIndex(int i, int j) const { return i * Ny() + j; }
    UnsignedInt vIndex(int i, int j) const { return i * (Ny() + 1) + j; }
    UnsignedInt numberOfCells() const { return mesh_.NumberOfCells(); }
    UnsignedInt numberOfUFaces() const { return (Nx() + 1) * Ny(); }
    UnsignedInt numberOfVFaces() const { return Nx() * (Ny() + 1); }

    int *cellType() { return cell_type_.Data(); }
    DiscreteVariable<int> *dvCellType() { return &cell_type_; }

    // Set every non-solid cell to Air (runs once per step on the host, before particles mark Fluid).
    void resetCellTypes()
    {
        for (UnsignedInt c = 0; c < numberOfCells(); ++c)
            if (cell_type_.Data()[c] != Solid)
                cell_type_.Data()[c] = Air;
    }

    // Direct access to the face velocities and accumulated weights, and to the full variables.
    Real *u() { return u_.Data(); }
    Real *v() { return v_.Data(); }
    Real *uWeight() { return u_weight_.Data(); }
    Real *vWeight() { return v_weight_.Data(); }
    DiscreteVariable<Real> *dvU() { return &u_; }
    DiscreteVariable<Real> *dvV() { return &v_; }
    DiscreteVariable<Real> *dvUWeight() { return &u_weight_; }
    DiscreteVariable<Real> *dvVWeight() { return &v_weight_; }

    // Set to zero the velocities and weights of all faces (before accumulating the particles).
    void resetFaces()
    {
        for (UnsignedInt f = 0; f < numberOfUFaces(); ++f)
        {
            u_.Data()[f] = 0.0;
            u_weight_.Data()[f] = 0.0;
        }
        for (UnsignedInt f = 0; f < numberOfVFaces(); ++f)
        {
            v_.Data()[f] = 0.0;
            v_weight_.Data()[f] = 0.0;
        }
    }

    // Divide the accumulated velocity by the accumulated weight on every face that received weight.
    void normalizeFaces()
    {
        for (UnsignedInt f = 0; f < numberOfUFaces(); ++f)
            if (u_weight_.Data()[f] > 0.0)
                u_.Data()[f] /= u_weight_.Data()[f];
        for (UnsignedInt f = 0; f < numberOfVFaces(); ++f)
            if (v_weight_.Data()[f] > 0.0)
                v_.Data()[f] /= v_weight_.Data()[f];
    }

    // Access to the copies of the face velocities taken right after the P2G.
    Real *uOld() { return u_old_.Data(); }
    Real *vOld() { return v_old_.Data(); }
    DiscreteVariable<Real> *dvUOld() { return &u_old_; }
    DiscreteVariable<Real> *dvVOld() { return &v_old_; }

    // Save a copy of the current face velocities (needed to compute the grid velocity change in FLIP).
    void saveOldFaces()
    {
        for (UnsignedInt f = 0; f < numberOfUFaces(); ++f)
            u_old_.Data()[f] = u_.Data()[f];
        for (UnsignedInt f = 0; f < numberOfVFaces(); ++f)
            v_old_.Data()[f] = v_.Data()[f];
    }

    // Gravity: subtract g * dt from the vertical velocity of all v faces.
    void addGravity(Real gravity, Real dt)
    {
        for (UnsignedInt f = 0; f < numberOfVFaces(); ++f)
            v_.Data()[f] -= gravity * dt;
    }

    // Walls: set to zero the normal velocity of the faces that touch a solid cell or the mesh edge.
    void applySolidBoundary()
    {
        for (int i = 0; i <= Nx(); ++i)
            for (int j = 0; j < Ny(); ++j)
            {
                bool left_solid = (i == 0) || cell_type_.Data()[cellIndex(i - 1, j)] == Solid;
                bool right_solid = (i == Nx()) || cell_type_.Data()[cellIndex(i, j)] == Solid;
                if (left_solid || right_solid)
                    u_.Data()[uIndex(i, j)] = 0.0;
            }
        for (int i = 0; i < Nx(); ++i)
            for (int j = 0; j <= Ny(); ++j)
            {
                bool below_solid = (j == 0) || cell_type_.Data()[cellIndex(i, j - 1)] == Solid;
                bool above_solid = (j == Ny()) || cell_type_.Data()[cellIndex(i, j)] == Solid;
                if (below_solid || above_solid)
                    v_.Data()[vIndex(i, j)] = 0.0;
            }
    }

    // Direct access to the cell pressure array.
    Real *p() { return p_.Data(); }

    // Divergence of cell (i, j): (u_right - u_left + v_top - v_bottom) / dx.
    Real divergence(int i, int j)
    {
        return (u_.Data()[uIndex(i + 1, j)] - u_.Data()[uIndex(i, j)] +
                v_.Data()[vIndex(i, j + 1)] - v_.Data()[vIndex(i, j)]) / dx_;
    }

    // Largest absolute divergence among the fluid cells.
    Real maxDivergence()
    {
        Real max_div = 0.0;
        for (int i = 0; i < Nx(); ++i)
            for (int j = 0; j < Ny(); ++j)
                if (cell_type_.Data()[cellIndex(i, j)] == Fluid)
                    max_div = std::max(max_div, std::abs(divergence(i, j)));
        return max_div;
    }

    // Solve the pressure Poisson equation on the fluid cells with a matrix-free conjugate gradient (CG).
    // Wall: zero pressure gradient. Air: zero pressure. The discrete equation in a fluid cell is
    //   n * p_c - sum(p of fluid neighbours) = -rho * dx^2 * divergence / dt,
    // where n is the number of non-solid neighbours.
    // Returns the number of iterations used and, in relative_residual, the final residual
    // |b - A x| / |b| recomputed explicitly from scratch.
    int solvePressure(Real dt, Real rho, Real tolerance, int max_iterations, Real &relative_residual)
    {
        const UnsignedInt n = numberOfCells();
        int *type = cell_type_.Data();
        std::vector<Real> b(n, 0.0), x(n, 0.0), r(n, 0.0), d(n, 0.0), Ad(n, 0.0);

        // Matrix-vector product of the Poisson operator: out = A * in (zero outside the fluid cells).
        auto apply = [&](const std::vector<Real> &in, std::vector<Real> &out)
        {
            const int di[4] = {-1, 1, 0, 0};
            const int dj[4] = {0, 0, -1, 1};
            for (int i = 0; i < Nx(); ++i)
                for (int j = 0; j < Ny(); ++j)
                {
                    UnsignedInt c = cellIndex(i, j);
                    out[c] = 0.0;
                    if (type[c] != Fluid)
                        continue;
                    int count = 0;
                    Real sum = 0.0;
                    for (int k = 0; k < 4; ++k)
                    {
                        UnsignedInt nb = cellIndex(i + di[k], j + dj[k]);
                        if (type[nb] == Solid)
                            continue;
                        ++count;
                        if (type[nb] == Fluid)
                            sum += in[nb];
                    }
                    out[c] = count * in[c] - sum;
                }
        };
        // Dot product of two cell vectors.
        auto dot = [&](const std::vector<Real> &a, const std::vector<Real> &e)
        {
            Real s = 0.0;
            for (UnsignedInt c = 0; c < n; ++c)
                s += a[c] * e[c];
            return s;
        };

        // Right-hand side: -rho * dx^2 * divergence / dt in every fluid cell.
        for (int i = 0; i < Nx(); ++i)
            for (int j = 0; j < Ny(); ++j)
            {
                UnsignedInt c = cellIndex(i, j);
                if (type[c] == Fluid)
                    b[c] = -rho * dx_ * dx_ * divergence(i, j) / dt;
            }

        // Conjugate gradient: we start from x = 0, so the initial residual is b.
        r = b;
        d = r;
        Real b_norm = std::sqrt(dot(b, b));
        Real rr = dot(r, r);
        int iteration = 0;
        if (b_norm > 0.0)
        {
            while (iteration < max_iterations && std::sqrt(rr) > tolerance * b_norm)
            {
                apply(d, Ad);
                Real alpha = rr / dot(d, Ad);
                for (UnsignedInt c = 0; c < n; ++c)
                {
                    x[c] += alpha * d[c];
                    r[c] -= alpha * Ad[c];
                }
                Real rr_new = dot(r, r);
                Real beta = rr_new / rr;
                for (UnsignedInt c = 0; c < n; ++c)
                    d[c] = r[c] + beta * d[c];
                rr = rr_new;
                ++iteration;
            }
        }

        // Final residual recomputed explicitly: |b - A x| / |b|.
        apply(x, Ad);
        Real res = 0.0;
        for (UnsignedInt c = 0; c < n; ++c)
            res += (b[c] - Ad[c]) * (b[c] - Ad[c]);
        relative_residual = (b_norm > 0.0) ? std::sqrt(res) / b_norm : 0.0;

        for (UnsignedInt c = 0; c < n; ++c)
            p_.Data()[c] = x[c];
        return iteration;
    }

    // Projection: subtract the pressure gradient from the velocity of the faces that touch fluid
    // (and do not touch a wall): u -= dt / (rho * dx) * (p_right - p_left), and the same for v.
    void applyPressureGradient(Real dt, Real rho)
    {
        int *type = cell_type_.Data();
        Real scale = dt / (rho * dx_);
        for (int i = 1; i < Nx(); ++i)
            for (int j = 0; j < Ny(); ++j)
            {
                int left = type[cellIndex(i - 1, j)];
                int right = type[cellIndex(i, j)];
                if (left == Solid || right == Solid)
                    continue;
                if (left == Fluid || right == Fluid)
                    u_.Data()[uIndex(i, j)] -= scale * (p_.Data()[cellIndex(i, j)] - p_.Data()[cellIndex(i - 1, j)]);
            }
        for (int i = 0; i < Nx(); ++i)
            for (int j = 1; j < Ny(); ++j)
            {
                int below = type[cellIndex(i, j - 1)];
                int above = type[cellIndex(i, j)];
                if (below == Solid || above == Solid)
                    continue;
                if (below == Fluid || above == Fluid)
                    v_.Data()[vIndex(i, j)] -= scale * (p_.Data()[cellIndex(i, j)] - p_.Data()[cellIndex(i, j - 1)]);
            }
    }

  private:
    int nx_;
    Real dx_;
    int ny_;
    Mesh mesh_;
    DiscreteVariable<Real> p_, u_, v_;
    DiscreteVariable<int> cell_type_;
    // Weights accumulated on the u and v faces during the particle-to-grid transfer.
    DiscreteVariable<Real> u_weight_, v_weight_;
    // Copies of the face velocities after the P2G, used to compute the grid velocity change in the G2P.
    DiscreteVariable<Real> u_old_, v_old_;
};

// Classify the cells: a non-solid cell that contains at least one particle becomes Fluid.
class ClassifyCellsCK : public LocalDynamics
{
  public:
    ClassifyCellsCK(SPHBody &sph_body, MacGrid &mac_grid)
        : LocalDynamics(sph_body), mac_grid_(mac_grid),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_cell_type_(mac_grid.dvCellType()) {}
    virtual ~ClassifyCellsCK() {};

    // Before looping over the particles: reset every non-solid cell to Air.
    virtual void setupDynamics(Real dt = 0.0) override { mac_grid_.resetCellTypes(); }

    // Kernel: executed for each particle, marks the cell that contains it as Fluid.
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : mesh_(encloser.mac_grid_.mesh()),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              cell_type_(encloser.dv_cell_type_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            UnsignedInt c = mesh_.LinearCellIndexFromPosition(pos_[index_i]);
            if (cell_type_[c] != MacGrid::Solid)
                cell_type_[c] = MacGrid::Fluid;
        }

      protected:
        Mesh mesh_;
        Vecd *pos_;
        int *cell_type_;
    };

  protected:
    MacGrid &mac_grid_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<int> *dv_cell_type_;
};

// Particle-to-grid transfer (P2G): each particle adds its velocity, with bilinear weights,
// to the 4 nearest u faces and the 4 nearest v faces; at the end each face is divided by its weight.
class ParticlesToGridCK : public LocalDynamics
{
  public:
    ParticlesToGridCK(SPHBody &sph_body, MacGrid &mac_grid)
        : LocalDynamics(sph_body), mac_grid_(mac_grid),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_vel_(particles_->registerStateVariable<Vecd>("Velocity")),
          dv_u_(mac_grid.dvU()), dv_v_(mac_grid.dvV()),
          dv_u_weight_(mac_grid.dvUWeight()), dv_v_weight_(mac_grid.dvVWeight()) {}
    virtual ~ParticlesToGridCK() {};

    // Before looping over the particles: set the faces to zero (runs once, serially).
    virtual void setupDynamics(Real dt = 0.0) override { mac_grid_.resetFaces(); }

    // Kernel: executed for each particle, adds its velocity to the nearby faces.
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : lower_(encloser.mac_grid_.mesh().MeshLowerBound()),
              dx_(encloser.mac_grid_.dx()),
              num_u_i_(encloser.mac_grid_.Nx() + 1), num_u_j_(encloser.mac_grid_.Ny()),
              num_v_i_(encloser.mac_grid_.Nx()), num_v_j_(encloser.mac_grid_.Ny() + 1),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
              u_(encloser.dv_u_->DelegatedData(ex_policy)),
              v_(encloser.dv_v_->DelegatedData(ex_policy)),
              u_weight_(encloser.dv_u_weight_->DelegatedData(ex_policy)),
              v_weight_(encloser.dv_v_weight_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            // Particle position in cell units, relative to the lower corner of the mesh.
            Real x = (pos_[index_i][0] - lower_[0]) / dx_;
            Real y = (pos_[index_i][1] - lower_[1]) / dx_;
            // u faces sit at half-cell height, v faces at half-cell width (hence the 0.5 shifts).
            splat(x, y - 0.5, vel_[index_i][0], u_, u_weight_, num_u_i_, num_u_j_);
            splat(x - 0.5, y, vel_[index_i][1], v_, v_weight_, num_v_i_, num_v_j_);
        }

      protected:
        // Spread a value over the 4 faces around position (fx, fy), in cell units.
        // Atomic additions avoid race conditions when particles are processed in parallel.
        void splat(Real fx, Real fy, Real value, Real *data, Real *weight, int num_i, int num_j)
        {
            int i0 = static_cast<int>(std::floor(fx));
            int j0 = static_cast<int>(std::floor(fy));
            Real tx = fx - i0, ty = fy - j0;
            for (int di = 0; di < 2; ++di)
                for (int dj = 0; dj < 2; ++dj)
                {
                    int i = i0 + di, j = j0 + dj;
                    if (i < 0 || i >= num_i || j < 0 || j >= num_j)
                        continue;
                    Real w = (di == 0 ? 1.0 - tx : tx) * (dj == 0 ? 1.0 - ty : ty);
                    AtomicRef<Real> atomic_data(data[i * num_j + j]);
                    AtomicRef<Real> atomic_weight(weight[i * num_j + j]);
                    atomic_data.fetch_add(w * value);
                    atomic_weight.fetch_add(w);
                }
        }

        Vecd lower_;
        Real dx_;
        int num_u_i_, num_u_j_, num_v_i_, num_v_j_;
        Vecd *pos_, *vel_;
        Real *u_, *v_, *u_weight_, *v_weight_;
    };

    // After looping over the particles: divide each face by its accumulated weight
    // and save a copy of the result (needed by the G2P to compute the FLIP velocity change).
    class FinishDynamics
    {
      public:
        FinishDynamics(ParticlesToGridCK &encloser) : mac_grid_(encloser.mac_grid_) {}
        void operator()()
        {
            mac_grid_.normalizeFaces();
            mac_grid_.saveOldFaces();
        }

      protected:
        MacGrid &mac_grid_;
    };

  protected:
    MacGrid &mac_grid_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_vel_;
    DiscreteVariable<Real> *dv_u_, *dv_v_, *dv_u_weight_, *dv_v_weight_;
};

// Grid-to-particle transfer (G2P): each particle reads the velocities of the nearby faces
// and updates its own velocity, blending PIC and FLIP according to flip_ratio
// (1 = pure FLIP, 0 = pure PIC):
//   v_p = flip_ratio * (v_p + v_grid_new - v_grid_old) + (1 - flip_ratio) * v_grid_new
class GridToParticlesCK : public LocalDynamics
{
  public:
    GridToParticlesCK(SPHBody &sph_body, MacGrid &mac_grid, Real flip_ratio)
        : LocalDynamics(sph_body), mac_grid_(mac_grid), flip_ratio_(flip_ratio),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_vel_(particles_->registerStateVariable<Vecd>("Velocity")),
          dv_u_(mac_grid.dvU()), dv_v_(mac_grid.dvV()),
          dv_u_old_(mac_grid.dvUOld()), dv_v_old_(mac_grid.dvVOld()) {}
    virtual ~GridToParticlesCK() {};

    // Kernel: executed for each particle; it only reads the grid and writes its own velocity.
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : lower_(encloser.mac_grid_.mesh().MeshLowerBound()),
              dx_(encloser.mac_grid_.dx()),
              num_u_i_(encloser.mac_grid_.Nx() + 1), num_u_j_(encloser.mac_grid_.Ny()),
              num_v_i_(encloser.mac_grid_.Nx()), num_v_j_(encloser.mac_grid_.Ny() + 1),
              flip_ratio_(encloser.flip_ratio_),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
              u_(encloser.dv_u_->DelegatedData(ex_policy)),
              v_(encloser.dv_v_->DelegatedData(ex_policy)),
              u_old_(encloser.dv_u_old_->DelegatedData(ex_policy)),
              v_old_(encloser.dv_v_old_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            Real x = (pos_[index_i][0] - lower_[0]) / dx_;
            Real y = (pos_[index_i][1] - lower_[1]) / dx_;
            Real u_new = gather(x, y - 0.5, u_, num_u_i_, num_u_j_);
            Real u_prev = gather(x, y - 0.5, u_old_, num_u_i_, num_u_j_);
            Real v_new = gather(x - 0.5, y, v_, num_v_i_, num_v_j_);
            Real v_prev = gather(x - 0.5, y, v_old_, num_v_i_, num_v_j_);
            vel_[index_i][0] = flip_ratio_ * (vel_[index_i][0] + u_new - u_prev) + (1.0 - flip_ratio_) * u_new;
            vel_[index_i][1] = flip_ratio_ * (vel_[index_i][1] + v_new - v_prev) + (1.0 - flip_ratio_) * v_new;
        }

      protected:
        // Interpolate at position (fx, fy), in cell units, from the 4 neighbouring grid values.
        // Indices are clamped so that particles near the border never read outside the array.
        Real gather(Real fx, Real fy, const Real *data, int num_i, int num_j) const
        {
            int i0 = static_cast<int>(std::floor(fx));
            int j0 = static_cast<int>(std::floor(fy));
            Real tx = fx - i0, ty = fy - j0;
            Real result = 0.0;
            for (int di = 0; di < 2; ++di)
                for (int dj = 0; dj < 2; ++dj)
                {
                    int i = std::min(std::max(i0 + di, 0), num_i - 1);
                    int j = std::min(std::max(j0 + dj, 0), num_j - 1);
                    Real w = (di == 0 ? 1.0 - tx : tx) * (dj == 0 ? 1.0 - ty : ty);
                    result += w * data[i * num_j + j];
                }
            return result;
        }

        Vecd lower_;
        Real dx_;
        int num_u_i_, num_u_j_, num_v_i_, num_v_j_;
        Real flip_ratio_;
        Vecd *pos_, *vel_;
        Real *u_, *v_, *u_old_, *v_old_;
    };

  protected:
    MacGrid &mac_grid_;
    Real flip_ratio_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_vel_;
    DiscreteVariable<Real> *dv_u_, *dv_v_, *dv_u_old_, *dv_v_old_;
};

// Particle motion: position += velocity * dt (explicit Euler), keeping every particle inside the tank
// (between lower and upper, with a small margin so that it does not end up exactly on the wall).
class MoveParticlesCK : public LocalDynamics
{
  public:
    MoveParticlesCK(SPHBody &sph_body, Vecd lower, Vecd upper, Real margin)
        : LocalDynamics(sph_body), lower_(lower), upper_(upper), margin_(margin),
          dv_pos_(particles_->getVariableByName<Vecd>("Position")),
          dv_vel_(particles_->registerStateVariable<Vecd>("Velocity")) {}
    virtual ~MoveParticlesCK() {};

    // Kernel: executed for each particle, moves it and keeps it inside the tank.
    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : lower_(encloser.lower_), upper_(encloser.upper_), margin_(encloser.margin_),
              pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
              vel_(encloser.dv_vel_->DelegatedData(ex_policy)) {}

        void update(size_t index_i, Real dt = 0.0)
        {
            pos_[index_i] += vel_[index_i] * dt;
            for (int d = 0; d < 2; ++d)
                pos_[index_i][d] = std::min(std::max(pos_[index_i][d], lower_[d] + margin_), upper_[d] - margin_);
        }

      protected:
        Vecd lower_, upper_;
        Real margin_;
        Vecd *pos_, *vel_;
    };

  protected:
    Vecd lower_, upper_;
    Real margin_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_vel_;
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with initial shape, materials and particles.
    //	The wall body is only used to visualise the tank; the FLIP walls are
    //	the solid cells of the MAC grid.
    //----------------------------------------------------------------------
    GeometricShapeBox initial_water_block(Transform(water_block_translation), water_block_halfsize, "WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    ComplexShape wall_complex_shape("WallBoundary");
    wall_complex_shape.add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
    wall_complex_shape.subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    SolidBody wall_boundary(sph_system, wall_complex_shape);
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    // Define the solver container that registers the particle dynamics.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.getMainMethodContainer();
    //----------------------------------------------------------------------
    // Define the FLIP numerical methods.
    //----------------------------------------------------------------------
    // MAC grid with cell size about twice the particle spacing (about 4 particles per cell in 2D).
    MacGrid mac_grid(DL, DH, 2.0 * particle_spacing_ref);
    // Classification of the cells (Air / Fluid / Solid) from the particle positions.
    auto &classify_cells = main_methods.addStateDynamics<ClassifyCellsCK>(water_block, mac_grid);
    // Particle-to-grid transfer.
    auto &particles_to_grid = main_methods.addStateDynamics<ParticlesToGridCK>(water_block, mac_grid);
    // Grid-to-particle transfer with the FLIP/PIC blend: 95 % FLIP and 5 % PIC.
    auto &grid_to_particles = main_methods.addStateDynamics<GridToParticlesCK>(water_block, mac_grid, 0.95);
    // Particle motion inside the tank (from (0, 0) to (DL, DH)).
    auto &move_particles = main_methods.addStateDynamics<MoveParticlesCK>(water_block, Vecd(0.0, 0.0), Vecd(DL, DH), 0.001 * mac_grid.dx());
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Real>(water_block, "Density");
    body_state_recorder.addToWrite<Vecd>(water_block, "Velocity");
    //----------------------------------------------------------------------
    //	First output before the integration loop.
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    //----------------------------------------------------------------------
    //	FLIP time loop: in every step the method steps are repeated in order.
    //----------------------------------------------------------------------
    Vecd *water_vel = water_block.getBaseParticles().registerStateVariableData<Vecd>("Velocity");
    UnsignedInt n_water = water_block.getBaseParticles().TotalRealParticles();
    Real flip_end_time = 5.0;
    Real physical_time = 0.0;
    Real next_output_time = 0.1;
    size_t flip_step = 0;
    while (physical_time < flip_end_time)
    {
        // Time step from a CFL condition: no particle should cross more than half a cell per step.
        Real v_max = 0.0;
        for (UnsignedInt k = 0; k < n_water; ++k)
            v_max = std::max(v_max, water_vel[k].norm());
        Real dt = 0.5 * mac_grid.dx() / (v_max + std::sqrt(gravity_g * mac_grid.dx()));

        classify_cells.exec();                // 1. cells with water, from the particle positions
        particles_to_grid.exec();             // 2. particle velocities to the grid (keeps a copy)
        mac_grid.addGravity(gravity_g, dt);   // 3. gravity on the grid
        mac_grid.applySolidBoundary();        //    and walls
        Real residual = 0.0;
        int iterations = mac_grid.solvePressure(dt, rho0_f, 1.0e-8, 2000, residual); // 4. pressure solve
        mac_grid.applyPressureGradient(dt, rho0_f);
        mac_grid.applySolidBoundary();
        grid_to_particles.exec();             // 5. grid velocities to the particles (95 % FLIP)
        move_particles.exec(dt);              // 6. move the particles
        physical_time += dt;
        ++flip_step;

        if (flip_step % 20 == 0)
            std::cout << "step=" << flip_step << " t=" << physical_time << " dt=" << dt
                      << " poisson_iterations=" << iterations << " max|v|=" << v_max << std::endl;
        if (physical_time >= next_output_time)
        {
            // Update the internal system time so that the .vtp file name carries the physical time.
            sph_system.svPhysicalTime().setValue(physical_time);
            body_state_recorder.writeToFile();
            next_output_time += 0.1;
        }
    }
    std::cout << "FLIP finished: steps=" << flip_step << " t=" << physical_time << std::endl;

    return 0;
};