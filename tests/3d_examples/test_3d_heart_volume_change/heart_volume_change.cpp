#include "heart_volume_change.h"
using namespace SPH;

void write_csv_files(
    const std::string &file_name,
    const std::string &parameter_1_name,
    const std::string &parameter_2_name,
    const std::string &parameter_3_name,
    const std::vector<Real> &parameter_1,
    const std::vector<Real> &parameter_2,
    const std::vector<Real> &parameter_3)
{
    std::ofstream myfile;
    myfile.open(file_name);
    myfile << parameter_1_name << ";" << parameter_2_name << ";" << parameter_3_name << "\n";
    for (size_t i = 0; i < parameter_1.size(); ++i)
    {
        myfile << parameter_1[i] << ";" << parameter_2[i] << ";" << parameter_3[i] << "\n";
    }
    myfile.close();
}

void MeshData::load(std::string path_to_mesh, Real scale)
{
    SimTK::PolygonalMesh inner_mesh;
    inner_mesh.loadStlFile(path_to_mesh);
    inner_mesh.scaleMesh(scale);
    stl_mesh.vertices.reserve(inner_mesh.getNumVertices());
    for (int i = 0; i < inner_mesh.getNumVertices(); i++)
    {
        const auto &p = inner_mesh.getVertexPosition(i);
        stl_mesh.vertices.push_back({Real(p[0]), Real(p[1]), Real(p[2])});
    }

    stl_mesh.faces.reserve(inner_mesh.getNumFaces());
    for (int i = 0; i < inner_mesh.getNumFaces(); i++)
    {
        auto f1 = inner_mesh.getFaceVertex(i, 0);
        auto f2 = inner_mesh.getFaceVertex(i, 1);
        auto f3 = inner_mesh.getFaceVertex(i, 2);
        stl_mesh.faces.push_back({f1, f2, f3});
    }
}

void MeshData::initialize()
{
    mesh_sdf.construct(stl_mesh.vertices, stl_mesh.faces);
}

void MeshData::translate(const Vec3d &translation)
{
    for (auto &vertex : stl_mesh.vertices)
    {
        vertex[0] += translation[0];
        vertex[1] += translation[1];
        vertex[2] += translation[2];
    }
}

IndexVector MeshData::get_ids_close_to_surface(const StdLargeVec<Vec3d> &pos_0, const IndexVector &all_surface_ids, Real distance) const
{
    // the ids are returned in ascending order naturally
    // very important for later algorithms
    IndexVector ids;
    for (auto index : all_surface_ids)
    {
        tmd::Result result = mesh_sdf.unsigned_distance(pos_0[index]);
        if (result.distance <= distance)
            ids.push_back(index);
    }
    return ids;
}

IndexVector MeshData::get_ids_close_to_surface(const StdLargeVec<Vec3d> &pos_0, Real distance) const
{
    // the ids are returned in ascending order naturally
    // very important for later algorithms
    IndexVector ids;
    for (size_t i = 0; i < pos_0.size(); ++i)
    {
        tmd::Result result = mesh_sdf.unsigned_distance(pos_0[i]);
        if (result.distance <= distance)
            ids.push_back(i);
    }
    return ids;
}

void MyocardiumSurfaces::init_surfaces(
    const MeshData &myocardium_mesh, const MeshData &lv_mesh, const MeshData &rv_mesh,
    Real lv_mesh_offset, Real rv_mesh_offset, Real smoothing_length)
{
    // step 0.
    // get all surface ids - not using BodySurface class as it uses dp instead of a fraction of smoothing_length
    // using 0.7 multiplier to smoothing length - surface particles should be 0.5*dp away from the surface
    // myo_surface_ids_ will be sorted by default
    // the particles were generated based on this mesh so no offset needed
    myo_surface_ids_ = myocardium_mesh.get_ids_close_to_surface(particles_.pos0_, smoothing_length * 0.7);

    // step 1.
    // first layer of particles are dp/2 away from the surface
    // mesh_offset accounts for potential surface mismatch between ventricle and myocardium mesh
    lv_ids_ = lv_mesh.get_ids_close_to_surface(particles_.pos0_, myo_surface_ids_, lv_mesh_offset + smoothing_length * 0.5);
    rv_ids_ = rv_mesh.get_ids_close_to_surface(particles_.pos0_, myo_surface_ids_, rv_mesh_offset + smoothing_length * 0.5);
    StdLargeVec<Vec3d> lv_pos; // lv particles positions
    StdLargeVec<Vec3d> rv_pos; // rv particles positions
    lv_pos.reserve(lv_ids_.size());
    rv_pos.reserve(rv_ids_.size());
    for (auto index : lv_ids_)
        lv_pos.push_back(particles_.pos0_[index]);
    for (auto index : rv_ids_)
        rv_pos.push_back(particles_.pos0_[index]);

    // step 2
    // first combine the lv and rv ids using std::merge to keep the ids sorted in the merged vector
    // lv and rv ids have to be sorted already
    IndexVector ventricles_ids;
    ventricles_ids.reserve(lv_ids_.size() + rv_ids_.size());
    std::merge(
        lv_ids_.begin(), lv_ids_.end(),
        rv_ids_.begin(), rv_ids_.end(),
        std::back_inserter(ventricles_ids));
    // set_difference strictly relies on sorted ranges!
    // myo_surface_ids_ is sorted by default, ventricles_ids was kept sorted
    std::set_difference(
        myo_surface_ids_.begin(), myo_surface_ids_.end(),
        ventricles_ids.begin(), ventricles_ids.end(),
        std::inserter(pericardium_ids_, pericardium_ids_.begin()));
}

void MyocardiumSurfaces::write_all_surfaces_as_obj(const std::string &output_path) const
{
    write_particles_as_obj(particles_.pos0_, lv_ids_, output_path + "lv_particles.obj");
    write_particles_as_obj(particles_.pos0_, rv_ids_, output_path + "rv_particles.obj");
    write_particles_as_obj(particles_.pos0_, pericardium_ids_, output_path + "pericardium_particles.obj");
    write_particles_as_obj(particles_.pos0_, myo_surface_ids_, output_path + "all_surface_particles.obj");
}

void SurfaceOperationsVentricle::init_srf_area(InnerRelation &inner_relation)
{
    // kernel_ptr
    Kernel *kernel_ptr = particles_.getSPHBody().sph_adaptation_->getKernel();
    // assuming uniform surface particles here
    // calculation based on integrating the area over the kernel
    Real smoothing_length = particles_.getSPHBody().sph_adaptation_->ReferenceSmoothingLength();
    std::cout << "smoothing_length: " << smoothing_length << std::endl;

    Real area_integral = 0.0;
    int num = 1e5; // discretization of radius
    Real radius = 2.0 * smoothing_length;
    Real dr = radius / Real(num);
    for (int i = 0; i < num; ++i)
    {
        Real r = (i + 0.5) * dr;
        Real kernel = kernel_ptr->W_3D(r / smoothing_length);
        Real dA = r * dr;
        area_integral += kernel * dA;
    }
    area_integral *= 2.0 * M_PI;
    std::cout << "non-normalized particle area: " << area_integral << std::endl;
    // fill the vector
    std::fill(srf_area_0_.begin(), srf_area_0_.end(), area_integral);

    // normalize the area by the kernel sum at the particle position
    // kernel sum is the sum of the neighbor kernel values at the particle position + the own kernel value (1 for KernelWendlandC2)
    for (size_t i = 0; i < ids_.size(); ++i)
    {
        size_t index_i = ids_[i];
        Neighborhood &inner_neighborhood = inner_relation.inner_configuration_[index_i];
        Real kernel_sum = kernel_ptr->W_3D(0);
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (std::find(ids_.begin(), ids_.end(), index_j) != ids_.end())
            {
                Real r_ij = inner_neighborhood.r_ij_[n];
                kernel_sum += kernel_ptr->W_3D(r_ij / smoothing_length);
            }
        }
        srf_area_0_[i] /= kernel_sum;
        srf_area_0_[i] *= 1.03; // empirical - accounting for not fully supported kernel - calibrated for sphere
    }

    // copy to current area
    srf_area_n_ = srf_area_0_;
    std::cout << "SurfaceOperationsVentricle init_srf_area: " << std::accumulate(srf_area_0_.begin(), srf_area_0_.end(), 0.0) << std::endl;
}

void SurfaceOperationsVentricle::update_srf_area()
{
    for (size_t i = 0; i < ids_.size(); ++i)
    {
        size_t index_i = ids_[i];
        srf_area_n_[i] = srf_area_0_[i] * (particles_.F_[index_i]).determinant() * ((particles_.F_[index_i]).inverse().transpose() * particles_.n0_[index_i]).norm();
    }
}

void SurfaceOperationsVentricle::update_flow_rate(Real dt)
{
    // increment previous flow rate
    Q_prev_ = Q_current_;
    // calculate current flow rate
    Q_current_ = 0.0;
    for (size_t i = 0; i < ids_.size(); ++i)
    {
        size_t index_i = ids_[i];
        // flow rate towards the normal is positive by definition
        // this will mean positive change in case of ventricular contraction
        Q_current_ += srf_area_n_[i] * (particles_.vel_[index_i]).dot(particles_.n_[index_i]);
    };
}

void SurfaceOperationsVentricle::update_flow_acc(Real dt)
{
    // first order backward difference formula
    dQ_dt_ = (Q_current_ - Q_prev_) / (dt + TinyReal);
    // record volume difference - trivial integration
    // the volume difference is opposite sign as the flow rate
    delta_V_ = -(Q_current_ + Q_prev_) / 2 * dt;
}