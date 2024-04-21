
#include "unstructured_mesh.h"
#include "mesh_helper.h"

namespace SPH
{
//=================================================================================================//
ANSYSMesh::ANSYSMesh(const std::string &full_path)
{
    getDataFromMeshFile(full_path);
    getElementCenterCoordinates();
    getMinimumDistanceBetweenNodes();
}
//=================================================================================================//
void ANSYSMesh::getDataFromMeshFile(const std::string &full_path)
{
    ifstream mesh_file; /*!< \brief File object for the Ansys ASCII mesh file. */
    mesh_file.open(full_path);
    if (mesh_file.fail())
    {
        cout << "Error:Check if the file exists." << endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    }
    string text_line;
    /*--- Read the dimension of the problem ---*/
    size_t dimension(0);
    MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);
    /*--- Check dimension ---*/
    if (dimension != Dimensions)
    {
        cout << "Error:the dimension of problem does not match input mesh." << endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    }
    /*--- Read the total number node points ---*/
    size_t number_of_points(0);
    MeshFileHelpers::numberofNodes(mesh_file, number_of_points, text_line);
   
    /*--- Read the node coordinates ---*/
    MeshFileHelpers::nodeCoordinates(mesh_file, node_coordinates_, text_line, dimension);
    /*--- Check number of node points ---*/
    if (node_coordinates_.size() != number_of_points)
    {
        cout << "Error:Total number of node points does not match data!" << endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    }
    size_t boundary_type(0);
    size_t number_of_elements(0);
    size_t mesh_type = 3;
    /*--- Read the total number of elements ---*/
    MeshFileHelpers::numberofElements(mesh_file, number_of_elements, text_line);
   
    /*Preparing and initializing the data structure of mesh topology and element node connection*/
    MeshFileHelpers::dataStruct(mesh_topology_, elements_nodes_connection_, number_of_elements, mesh_type, dimension);

    while (getline(mesh_file, text_line))
    {
        if (text_line.find("(13", 0) != string::npos && text_line.find(")(", 0) != string::npos)
        {
            /*--- find the type of boundary condition ---*/
            boundary_type = MeshFileHelpers::findBoundaryType(text_line, boundary_type);
            types_of_boundary_condition_.push_back(boundary_type);
            while (getline(mesh_file, text_line))
            {
                if (text_line.find(")", 0) == string::npos)
                {
                    Vecd nodes = MeshFileHelpers::nodeIndex(text_line);
                    Vec2d cells = MeshFileHelpers::cellIndex(text_line);
                    /*--- build up all topology---*/
                    bool check_neighbor_cell1 = 1;
                    bool check_neighbor_cell2 = 0;
                    for (int cell1_cell2 = 0; cell1_cell2 != cells.size(); ++cell1_cell2)
                    {
                            if (mesh_type == 3)
                            {
                                if (cells[check_neighbor_cell2] != 0)
                                {
                                    MeshFileHelpers::updateElementsNodesConnection(elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2);
                                    MeshFileHelpers::updateCellLists(mesh_topology_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2, boundary_type);
                                    MeshFileHelpers::updateBoundaryCellLists(mesh_topology_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2, boundary_type);    
                                }
                                if (cells[check_neighbor_cell2] == 0)
                                {
                                    break;
                                }
                            }
                    }
                }
                else
                    break;
            }
        }
        if (text_line.find("Zone Sections", 0) != string::npos)
            break;
        if (text_line.find(")") != string::npos)
            continue;
    }
    mesh_topology_.erase(mesh_topology_.begin());
}
//=================================================================================================//
void ANSYSMesh::getElementCenterCoordinates()
{
    elements_centroids_.resize(elements_nodes_connection_.size());
    elements_volumes_.resize(elements_nodes_connection_.size());
    for (std::size_t element = 1; element != elements_nodes_connection_.size(); ++element)
    {
        Vecd center_coordinate = Vecd::Zero();
        MeshFileHelpers::cellCenterCoordinates(elements_nodes_connection_, element, node_coordinates_, elements_centroids_, center_coordinate);
        MeshFileHelpers::elementVolume(elements_nodes_connection_, element, node_coordinates_, elements_volumes_);
    }
    elements_volumes_.erase(elements_volumes_.begin());
    elements_centroids_.erase(elements_centroids_.begin());
    elements_nodes_connection_.erase(elements_nodes_connection_.begin());
}
//=================================================================================================//
void ANSYSMesh::getMinimumDistanceBetweenNodes()
{
    vector<Real> all_data_of_distance_between_nodes;
    all_data_of_distance_between_nodes.resize(0);
    MeshFileHelpers::minimumdistance(all_data_of_distance_between_nodes, elements_volumes_, mesh_topology_, node_coordinates_);
    auto min_distance_iter = std::min_element(all_data_of_distance_between_nodes.begin(), all_data_of_distance_between_nodes.end());
    if (min_distance_iter != all_data_of_distance_between_nodes.end())
    {
        min_distance_between_nodes_ = *min_distance_iter;
    }
    else
    {
        cout << "The array of all distance between nodes is empty " << endl;
    }
}
//=================================================================================================//
void BaseInnerRelationInFVM::resetNeighborhoodCurrentSize()
{
    parallel_for(
        IndexRange(0, base_particles_.total_real_particles_),
        [&](const IndexRange &r)
        {
            for (size_t num = r.begin(); num != r.end(); ++num)
            {
                inner_configuration_[num].current_size_ = 0;
            }
        },
        ap);
}
//=================================================================================================//
BaseInnerRelationInFVM::BaseInnerRelationInFVM(RealBody &real_body, ANSYSMesh &ansys_mesh)
    : BaseInnerRelation(real_body), real_body_(&real_body),
      node_coordinates_(ansys_mesh.node_coordinates_),
      mesh_topology_(ansys_mesh.mesh_topology_)
{
    subscribeToBody();
    inner_configuration_.resize(base_particles_.real_particles_bound_, Neighborhood());
};
//=================================================================================================//
void NeighborBuilderInFVM::createRelation(Neighborhood &neighborhood, Real &distance,
                                          Real &dW_ij, Vecd &interface_normal_direction, size_t j_index) const
{
    neighborhood.j_.push_back(j_index);
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(interface_normal_direction);
    neighborhood.dW_ij_.push_back(dW_ij);
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborBuilderInFVM::initializeRelation(Neighborhood &neighborhood, Real &distance,
                                              Real &dW_ij, Vecd &interface_normal_direction, size_t j_index) const
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = j_index;
    neighborhood.dW_ij_[current_size] = dW_ij;
    neighborhood.r_ij_[current_size] = distance;
    neighborhood.e_ij_[current_size] = interface_normal_direction;
}
//=================================================================================================//
InnerRelationInFVM::InnerRelationInFVM(RealBody &real_body, ANSYSMesh &ansys_mesh)
    : BaseInnerRelationInFVM(real_body, ansys_mesh), get_inner_neighbor_(&real_body){};
//=================================================================================================//
template <typename GetParticleIndex, typename GetNeighborRelation>
void InnerRelationInFVM::searchNeighborsByParticles(size_t total_particles, BaseParticles &source_particles,
                                                    ParticleConfiguration &particle_configuration,
                                                    GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation)
{
    parallel_for(
        IndexRange(0, base_particles_.total_real_particles_),
        [&](const IndexRange &r)
        {
            StdLargeVec<Vecd> &pos_n = source_particles.pos_;
            StdLargeVec<Real> &Vol_n = source_particles.Vol_;
            for (size_t num = r.begin(); num != r.end(); ++num)
            {
                size_t index_i = get_particle_index(num);
                Vecd &particle_position = pos_n[index_i];
                Real &Vol_i = Vol_n[index_i];

                Neighborhood &neighborhood = particle_configuration[index_i];
                for (std::size_t neighbor = 0; neighbor != mesh_topology_[index_i].size(); ++neighbor)
                {
                    size_t index_j = mesh_topology_[index_i][neighbor][0] - 1;
                    size_t boundary_type = mesh_topology_[index_i][neighbor][1];
                    size_t interface_node1_index = mesh_topology_[index_i][neighbor][2];
                    size_t interface_node2_index = mesh_topology_[index_i][neighbor][3];
                    Vecd node1_position = Vecd(node_coordinates_[interface_node1_index][0], node_coordinates_[interface_node1_index][1]);
                    Vecd node2_position = Vecd(node_coordinates_[interface_node2_index][0], node_coordinates_[interface_node2_index][1]);
                    Vecd interface_area_vector = node1_position - node2_position;
                    Real interface_area_size = interface_area_vector.norm();
                    Vecd unit_vector = interface_area_vector / interface_area_size;
                    // normal unit vector
                    Vecd normal_vector = Vecd(unit_vector[1], -unit_vector[0]);
                    // judge the direction
                    Vecd node1_to_center_direction = particle_position - node1_position;
                    if (node1_to_center_direction.dot(normal_vector) < 0)
                    {
                        normal_vector = -normal_vector;
                    };
                    Real r_ij = 0; // we need r_ij to calculate the viscous force
                    // boundary_type == 2 means both of them are inside of fluid
                    if (boundary_type == 2)
                    {
                        r_ij = (particle_position - pos_n[index_j]).dot(normal_vector);
                    }
                    //this refer to the different types of wall boundary condtions
                    if ((boundary_type == 3) | (boundary_type == 4) | (boundary_type == 5) | (boundary_type == 7) | (boundary_type == 9) | (boundary_type == 10) | (boundary_type == 36))
                    {
                        r_ij = node1_to_center_direction.dot(normal_vector) * 2.0;
                    }
                    Real dW_ij = -interface_area_size  / (2.0 * Vol_i * Vol_n[index_j]);
                    get_neighbor_relation(neighborhood, r_ij, dW_ij, normal_vector, index_j);
                }
            }
        },
        ap);
}
//=================================================================================================//
void InnerRelationInFVM::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    searchNeighborsByParticles(base_particles_.total_real_particles_,
                               base_particles_, inner_configuration_,
                               get_particle_index_, get_inner_neighbor_);
}
//=================================================================================================//
GhostCreationFromMesh::GhostCreationFromMesh(RealBody &real_body, ANSYSMesh &ansys_mesh,
                                             Ghost<ReserveSizeFactor> &ghost_boundary)
    : GeneralDataDelegateSimple(real_body),
      ghost_boundary_(ghost_boundary),
      node_coordinates_(ansys_mesh.node_coordinates_),
      mesh_topology_(ansys_mesh.mesh_topology_),
      pos_(particles_->pos_), Vol_(particles_->Vol_),
      ghost_bound_(ghost_boundary.GhostBound())
{
    ghost_boundary.checkParticlesReserved();
    each_boundary_type_with_all_ghosts_index_.resize(50);
    each_boundary_type_with_all_ghosts_eij_.resize(50);
    each_boundary_type_contact_real_index_.resize(50);
    addGhostParticleAndSetInConfiguration();
}
//=================================================================================================//
void GhostCreationFromMesh::addGhostParticleAndSetInConfiguration()
{
    ghost_bound_.second = ghost_bound_.first;

    for (size_t index_i = 0; index_i != particles_->total_real_particles_; ++index_i)
    {
        for (size_t neighbor_index = 0; neighbor_index != mesh_topology_[index_i].size(); ++neighbor_index)
        {
            size_t boundary_type = mesh_topology_[index_i][neighbor_index][1];
            if (mesh_topology_[index_i][neighbor_index][1] != 2)
            {
                mutex_create_ghost_particle_.lock();
                size_t ghost_particle_index = ghost_bound_.second;
                ghost_bound_.second++;
                ghost_boundary_.checkWithinGhostSize(ghost_bound_);

                particles_->updateGhostParticle(ghost_particle_index, index_i);
                size_t node1_index = mesh_topology_[index_i][neighbor_index][2];
                size_t node2_index = mesh_topology_[index_i][neighbor_index][3];
                Vecd node1_position = node_coordinates_[node1_index];
                Vecd node2_position = node_coordinates_[node2_index];
                Vecd ghost_particle_position = 0.5 * (node1_position + node2_position);

                mesh_topology_[index_i][neighbor_index][0] = ghost_particle_index + 1;
                pos_[ghost_particle_index] = ghost_particle_position;
                mutex_create_ghost_particle_.unlock();

                std::vector<std::vector<size_t>> new_element;
                // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                std::vector<size_t> sub_element1 = {index_i + 1, boundary_type, node1_index, node2_index};
                new_element.push_back(sub_element1);
                // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                std::vector<size_t> sub_element2 = {index_i + 1, boundary_type, node1_index, node2_index};
                new_element.push_back(sub_element2);
                // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                std::vector<size_t> sub_element3 = {index_i + 1, boundary_type, node1_index, node2_index};
                new_element.push_back(sub_element3);
                // Add the new element to mesh_topology_
                mesh_topology_.push_back(new_element);
                // creating the boundary files with ghost particle index
                each_boundary_type_with_all_ghosts_index_[boundary_type].push_back(ghost_particle_index);
                // creating the boundary files with contact real particle index
                each_boundary_type_contact_real_index_[boundary_type].push_back(index_i);

                // creating the boundary files with ghost eij
                Vecd interface_area_vector = node1_position - node2_position;
                Real interface_area_size = interface_area_vector.norm();
                Vecd unit_vector = interface_area_vector / interface_area_size;
                // normal unit vector
                Vecd normal_vector = Vecd(unit_vector[1], -unit_vector[0]);
                // judge the direction
                Vecd particle_position = pos_[index_i];
                Vecd node1_to_center_direction = particle_position - node1_position;
                if (node1_to_center_direction.dot(normal_vector) < 0)
                {
                    normal_vector = -normal_vector;
                };
                each_boundary_type_with_all_ghosts_eij_[boundary_type].push_back(normal_vector);
            }
        }
    }
};
//=================================================================================================//
BodyStatesRecordingInMeshToVtp::BodyStatesRecordingInMeshToVtp(SPHBody &body, ANSYSMesh &ansys_mesh)
    : BodyStatesRecording(body), node_coordinates_(ansys_mesh.node_coordinates_),
      elements_nodes_connection_(ansys_mesh.elements_nodes_connection_){};
//=================================================================================================//
void BodyStatesRecordingInMeshToVtp::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated() && state_recording_)
        {
            // TODO: we can short the file name by without using SPHBody
            std::string filefullpath = io_environment_.output_folder_ + "/" + body->getName() + "_" + sequence + ".vtp";
            if (fs::exists(filefullpath))
            {
                fs::remove(filefullpath);
            }
            std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
            // begin of the XML file
            out_file << "<?xml version=\"1.0\"?>\n";
            out_file << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
            out_file << "<PolyData>\n";

            // Write point data
            out_file << "<Piece NumberOfPoints=\"" << node_coordinates_.size()
                     << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
                     << elements_nodes_connection_.size() << "\">\n";
            out_file << "<Points>\n";
            out_file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

            size_t total_nodes = node_coordinates_.size();
            for (size_t node = 0; node != total_nodes; ++node)
            {
                Vec3d particle_position = upgradeToVec3d(node_coordinates_[node]);
                out_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << "\n";
                
            }

            out_file << "</DataArray>\n";
            out_file << "</Points>\n";

            // Write face data
            out_file << "<Polys>\n";
            out_file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

            for (const auto &element : elements_nodes_connection_)
            {
                for (const auto &vertex : element)
                {
                    out_file << vertex << " ";
                }
                out_file << "\n";
            }

            out_file << "</DataArray>\n";
            out_file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";

            size_t offset = 0;
            for (const auto &face : elements_nodes_connection_)
            {
                offset += face.size();
                out_file << offset << " ";
            }

            out_file << "\n</DataArray>\n";
            out_file << "</Polys>\n";

            // Write face attribute data
            out_file << "<CellData>\n";
            body->writeParticlesToVtpFile(out_file);

            out_file << "</CellData>\n";

            // Write file footer
            out_file << "</Piece>\n";
            out_file << "</PolyData>\n";
            out_file << "</VTKFile>\n";

            out_file.close();
        }
        body->setNotNewlyUpdated();
    }
}
//=================================================================================================//
BoundaryConditionSetupInFVM::
    BoundaryConditionSetupInFVM(BaseInnerRelationInFVM &inner_relation, GhostCreationFromMesh &ghost_creation)
    : fluid_dynamics::FluidDataInner(inner_relation), rho_(particles_->rho_),
      Vol_(particles_->Vol_), mass_(particles_->mass_),
      p_(*particles_->getVariableByName<Real>("Pressure")),
      vel_(particles_->vel_), pos_(particles_->pos_),
      mom_(*particles_->getVariableByName<Vecd>("Momentum")),
      ghost_bound_(ghost_creation.ghost_bound_),
      each_boundary_type_with_all_ghosts_index_(ghost_creation.each_boundary_type_with_all_ghosts_index_),
      each_boundary_type_with_all_ghosts_eij_(ghost_creation.each_boundary_type_with_all_ghosts_eij_),
      each_boundary_type_contact_real_index_(ghost_creation.each_boundary_type_contact_real_index_){};
//=================================================================================================//
void BoundaryConditionSetupInFVM::resetBoundaryConditions()
{
    for (size_t boundary_type = 0; boundary_type < each_boundary_type_with_all_ghosts_index_.size(); ++boundary_type)
    {
        if (!each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
        {
            for (size_t ghost_number = 0; ghost_number != each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
            {
                size_t ghost_index = each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                size_t index_i = each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                Vecd e_ij = each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];

                // Dispatch the appropriate boundary condition
                switch (boundary_type)
                {
                case 3: // this refer to the different types of wall boundary conditions
                    applyNonSlipWallBoundary(ghost_index, index_i);
                    applyReflectiveWallBoundary(ghost_index, index_i, e_ij);
                    break;
                case 4:
                    applyTopBoundary(ghost_index, index_i);
                    break;
                case 5:
                    applyPressureOutletBC(ghost_index, index_i);
                    break; 
                case 7:
                    applySymmetryBoundary(ghost_index, index_i, e_ij);
                    break;
                case 9:
                    applyFarFieldBoundary(ghost_index);
                    break;
                case 10:
                    applyGivenValueInletFlow(ghost_index);
                    applyVelocityInletFlow(ghost_index, index_i);
                    break;
                case 36:
                    applyOutletBoundary(ghost_index, index_i);
                    break;
                }
            }
        }
    }
}
//=============================================================================================//
} // namespace SPH
