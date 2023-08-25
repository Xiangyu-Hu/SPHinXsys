/**
 * @file 	heart_volume_change.cpp
 * @brief 	This is the case studying the electromechanics on a biventricular heart model in 3D, including the volume change of the ventricles.
 * @author 	John Benjamin, Chi Zhang and Xiangyu Hu
 */
#include "TriangleMeshDistance.h" // for mesh operations
#include "sphinxsys.h"            // SPHinXsys Library.
#include <numeric>

using namespace SPH;

//----------------------------------------------------------------------
//	Identify Surfaces and Calculate volume change
//----------------------------------------------------------------------
void write_csv_files(
    const std::string &file_name,
    const std::string &parameter_1_name,
    const std::string &parameter_2_name,
    const std::string &parameter_3_name,
    const std::vector<Real> &parameter_1,
    const std::vector<Real> &parameter_2,
    const std::vector<Real> &parameter_3);

template <typename VectorType>
void write_particles_as_obj(const VectorType &pos, const IndexVector &ids, const std::string &file_path)
{
    std::ofstream stream;
    stream.open(file_path);
    for (auto id : ids)
        stream << "v " << pos[id][0] << " " << pos[id][1] << " " << pos[id][2] << "\n";
    stream.close();
}

class MeshData
{
  private:
    struct StlData
    {
        std::string name;
        std::uintptr_t ptr;
    };
    struct Mesh
    {
        std::vector<std::array<Real, 3>> vertices;
        std::vector<std::array<int, 3>> faces;
    };

    Mesh stl_mesh;
    tmd::TriangleMeshDistance mesh_sdf;

  public:
    void load(std::string path_to_mesh, Real scale);
    void initialize();
    void translate(const Vec3d &translation);
    IndexVector get_ids_close_to_surface(const StdLargeVec<Vec3d> &pos_0, const IndexVector &all_surface_ids, Real distance) const;
    IndexVector get_ids_close_to_surface(const StdLargeVec<Vec3d> &pos_0, Real distance) const;
};

// MyocardiumSurfaces
// Algorithm for finding lv, rv, and pericardium
// 1. get lv and rv particle ids and positions
// 2. remove the lv and rv particles from the full myocardium surface to get the pericardium
class MyocardiumSurfaces
{
  public:
    MyocardiumSurfaces(SolidBody &body, ElasticSolidParticles &particles) : particles_(particles){};
    ~MyocardiumSurfaces() = default;

    // mesh_offset: max distance between myocardium and ventricle mesh
    // important to use different offsets as the RV wall is much thinner usually
    // smoothing_length is used for dbscan
    void init_surfaces(
        const MeshData &myocardium_mesh, const MeshData &lv_mesh, const MeshData &rv_mesh,
        Real lv_mesh_offset, Real rv_mesh_offset, Real smoothing_length);

    // verbose
    void write_all_surfaces_as_obj(const std::string &output_path) const;

    // get functions
    inline const IndexVector &get_lv_ids() const { return lv_ids_; }
    inline const IndexVector &get_rv_ids() const { return rv_ids_; }
    inline const IndexVector &get_pericardium_ids() const { return pericardium_ids_; }

  private:
    ElasticSolidParticles &particles_;
    IndexVector myo_surface_ids_; // full myocardium surface
    IndexVector lv_ids_;
    IndexVector rv_ids_;
    IndexVector pericardium_ids_;
};

class SurfaceOperationsVentricle
{
  public:
    SurfaceOperationsVentricle(ElasticSolidParticles &particles, const IndexVector &ids, InnerRelation &inner_relation) : particles_(particles), ids_(ids),
                                                                                                                          srf_area_0_(ids_.size(), 0), srf_area_n_(ids_.size(), 0),
                                                                                                                          Q_current_(0), Q_prev_(0), dQ_dt_(0), delta_V_(0)
    {
        std::cout << "SurfaceOperationsVentricle number of particles: " << ids_.size() << std::endl;
        init_srf_area(inner_relation);
    };
    ~SurfaceOperationsVentricle() = default;

    // updates
    void update_srf_area();
    virtual void update_flow_rate(Real dt); // flow rate
    void update_flow_acc(Real dt);          // flow acceleration and delta volume

    // get functions
    inline Real get_Q_current() const { return Q_current_; }
    inline Real get_Q_prev() const { return Q_prev_; }
    inline Real get_dQ_dt() const { return dQ_dt_; }
    inline Real get_delta_V() const { return delta_V_; }

    // statistics
    inline size_t num_particles() const { return ids_.size(); }
    inline size_t particle_index(size_t i) const { return ids_[i]; }
    inline Real accumulate_srf_0() const { return std::accumulate(srf_area_0_.begin(), srf_area_0_.end(), 0.0); }
    inline Real accumulate_srf_n() const { return std::accumulate(srf_area_n_.begin(), srf_area_n_.end(), 0.0); }

  protected:
    // called in constructor
    void init_srf_area(InnerRelation &inner_relation);

    ElasticSolidParticles &particles_;
    // ids_, srf_area_0_, srf_area_n_ maintain particle correspondence
    const IndexVector &ids_;  // ids of srf particles of interest
    StdVec<Real> srf_area_0_; // initial surface area
    StdVec<Real> srf_area_n_; // current surface area
    // parameters for Windkessel models
    Real Q_current_; // current flow rate = flow rate
    Real Q_prev_;    // previous step flow rate = flow rate
    Real dQ_dt_;     // flow rate derivative
    // parameter for recording
    Real delta_V_; // volume difference in (only) the current step
};