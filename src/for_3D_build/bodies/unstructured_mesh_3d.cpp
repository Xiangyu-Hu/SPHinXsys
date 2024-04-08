
#include "unstructured_mesh_3d.h"
#include "MeshHelper.h"
namespace SPH
{
    //=================================================================================================//
    ANSYSMesh_3d::ANSYSMesh_3d(const std::string &full_path)
    {
        getDataFromMeshFile3d(full_path);
        getElementCenterCoordinates();
        getMinimumDistanceBetweenNodes();
    }
//===
    //=================================================================================================//
    void ANSYSMesh_3d::getDataFromMeshFile3d(const std::string &full_path)
    {
        Real ICEM = 0;
        ifstream mesh_file; /*!< \brief File object for the Ansys ASCII mesh file. */
        mesh_file.open(full_path);
        if (mesh_file.fail())
        {
            cout << "Error:Check if the file exists." << endl;
        }
        string text_line;
        /*Check mesh file generation software*/
        (getline(mesh_file, text_line));
        text_line.erase(0, 4);
        istringstream value(text_line);
        if (text_line.find("Created by", 0) != string::npos)
        {
            ICEM = 1;
        }
       
        if (ICEM == 1)
        {

            /*--- Read the dimension of the mesh ---*/
            size_t dimension(0);
            MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);

            /*--- Read the node data (index is starting from zero) ---*/
            size_t number_of_points(0);
            MeshFileHelpers::numberofNodes(mesh_file, number_of_points, text_line);

            /* Prepare our data structure for the point coordinates. */
            node_coordinates_.resize(number_of_points);
            for (std::vector<std::vector<double>>::size_type point = 0; point != node_coordinates_.size(); ++point)
            {
                node_coordinates_[point].resize(dimension);
            }

            point_coordinates.resize(dimension);
            for (std::size_t k = 0; k < dimension; k++)
            {
                point_coordinates[k].reserve(number_of_points);
            }

            size_t node_index = 0;
            MeshFileHelpers::nodeCoordinates(mesh_file, node_index, node_coordinates_, text_line, dimension, point_coordinates);
            
            size_t boundary_type(0);
            size_t number_of_elements(0);
            size_t mesh_type = 4;
            MeshFileHelpers::numberofElements(mesh_file, number_of_elements, text_line);

            /*Preparing and initializing the data structure of mesh topology and element node connection*/
            MeshFileHelpers::dataStruct(mesh_topology_, elements_nodes_connection_, number_of_elements, mesh_type, dimension, elements_neighbors_connection_);

            /*--- find the elements lines ---*/
            while (getline(mesh_file, text_line))
            {
                if (text_line.find("(13", 0) != string::npos && text_line.find(")(", 0) != string::npos)
                {
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
                                if (mesh_type == 4)
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
        else /*This section is for mesh files created from fluent*/
        {
            size_t dimension(0);
            MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);

            /*--- Read the node data (index is starting from zero) ---*/
            size_t number_of_points(0);
            MeshFileHelpers::numberofNodesFluent(mesh_file, number_of_points, text_line);
            node_coordinates_.resize(number_of_points);
            for (std::vector<std::vector<double>>::size_type point = 0; point != node_coordinates_.size(); ++point)
            {
                node_coordinates_[point].resize(dimension);
            }
            point_coordinates.resize(dimension);
            for (std::size_t k = 0; k < dimension; k++)
            {
                point_coordinates[k].reserve(number_of_points);
            }
            
            size_t node_index = 0;
            size_t boundary_type(0);
            size_t number_of_elements(0);
            size_t mesh_type = 4;

            MeshFileHelpers::numberofElementsFluent(mesh_file, number_of_elements, text_line);
            MeshFileHelpers::dataStruct(mesh_topology_, elements_nodes_connection_, number_of_elements, mesh_type, dimension, elements_neighbors_connection_);
            MeshFileHelpers::nodeCoordinatesFluent( mesh_file, node_index, node_coordinates_, text_line, dimension, point_coordinates);

            while (getline(mesh_file, text_line))
            {

                if (text_line.find("(13", 0) != string::npos && text_line.find(") (", 0) != string::npos)
                {
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
                                if (mesh_type == 4)
                                {
                                    if (cells[check_neighbor_cell2] != 0)
                                    {
                                        MeshFileHelpers::updateElementsNodesConnection(elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2);
                                        MeshFileHelpers::updateCellLists(mesh_topology_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, 
                                                                        check_neighbor_cell2, boundary_type);
                                        MeshFileHelpers::updateBoundaryCellListsFluent(mesh_topology_, elements_nodes_connection_, nodes, cells,
                                                                        check_neighbor_cell1, check_neighbor_cell2, boundary_type);
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
                if (text_line.find("(12", 0) != string::npos && text_line.find("))", 0) != string::npos)
                    break;
                if (text_line.find(")") != string::npos)
                    continue;
            }
            mesh_topology_.erase(mesh_topology_.begin());
            
        }
        
    }
    //=================================================================================================//

        void ANSYSMesh_3d::getElementCenterCoordinates()
        {
            elements_centroids_.resize(elements_nodes_connection_.size());
            elements_volumes_.resize(elements_nodes_connection_.size());
            for (std::vector<std::vector<long unsigned int>>::size_type element = 1; element != elements_nodes_connection_.size(); ++element)
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

        void ANSYSMesh_3d::getMinimumDistanceBetweenNodes()
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
        
        GeneratingMethod<UnstructuredMesh_3d>::GeneratingMethod(ANSYSMesh_3d& ansys_mesh_3d)
            : elements_centroids_(ansys_mesh_3d.elements_centroids_), // Assume ANSYSMesh_3d has a similar interface
            elements_volumes_(ansys_mesh_3d.elements_volumes_) {}
        //=================================================================================================//
        ParticleGenerator<UnstructuredMesh_3d>::ParticleGenerator(SPHBody& sph_body, ANSYSMesh_3d& ansys_mesh_3d)
            : ParticleGenerator<Base>(sph_body), GeneratingMethod<UnstructuredMesh_3d>(ansys_mesh_3d) {}
        //=================================================================================================//
        void ParticleGenerator<UnstructuredMesh_3d>::initializeGeometricVariables()
        {
            for (size_t i = 0; i != elements_centroids_.size(); ++i)
            {
                initializePositionAndVolumetricMeasure(elements_centroids_[i], elements_volumes_[i]);
            }
        }
        //=================================================================================================//
    
        BaseInnerRelationInFVM_3d::BaseInnerRelationInFVM_3d(RealBody &real_body, ANSYSMesh_3d& ansys_mesh)
        : BaseInnerRelation(real_body), real_body_(&real_body), node_coordinates_(ansys_mesh.node_coordinates_),
        mesh_topology_(ansys_mesh.mesh_topology_)
        {
            subscribeToBody();
            inner_configuration_.resize(base_particles_.real_particles_bound_, Neighborhood());
        };
    //=================================================================================================//
    void BaseInnerRelationInFVM_3d::resetNeighborhoodCurrentSize()
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
    void NeighborBuilderInFVM_3d::createRelation(Neighborhood &neighborhood, Real &distance,
                                            Real &dW_ij, Vecd &interface_normal_direction, size_t j_index) const
    {
        neighborhood.j_.push_back(j_index);
        neighborhood.r_ij_.push_back(distance);
        neighborhood.e_ij_.push_back(interface_normal_direction);
        neighborhood.dW_ij_.push_back(dW_ij);
        neighborhood.allocated_size_++;
    }
    //=================================================================================================//
    void NeighborBuilderInFVM_3d::initializeRelation(Neighborhood &neighborhood, Real &distance,
                                                Real &dW_ij, Vecd &interface_normal_direction, size_t j_index) const
    {
        size_t current_size = neighborhood.current_size_;
        neighborhood.j_[current_size] = j_index;
        neighborhood.dW_ij_[current_size] = dW_ij;
        neighborhood.r_ij_[current_size] = distance;
        neighborhood.e_ij_[current_size] = interface_normal_direction;
    }

    //=================================================================================================//
    InnerRelationInFVM_3d::InnerRelationInFVM_3d(RealBody &real_body, ANSYSMesh_3d& ansys_mesh)
        : BaseInnerRelationInFVM_3d(real_body, ansys_mesh), get_inner_neighbor_(&real_body){};
    //=================================================================================================//
    template <typename GetParticleIndex, typename GetNeighborRelation>
    void InnerRelationInFVM_3d::searchNeighborsByParticles(size_t total_particles, BaseParticles &source_particles,
                                                        ParticleConfiguration &particle_configuration, GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation)
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
                    for (std::vector<std::vector<long unsigned int>>::size_type neighbor = 0; neighbor != mesh_topology_[index_i].size(); ++neighbor)
                    {
                        size_t index_j = mesh_topology_[index_i][neighbor][0] - 1;
                        size_t boundary_type = mesh_topology_[index_i][neighbor][1];
                        size_t interface_node1_index = mesh_topology_[index_i][neighbor][2];
                        size_t interface_node2_index = mesh_topology_[index_i][neighbor][3];
                        size_t interface_node3_index = mesh_topology_[index_i][neighbor][4];
                        Vecd node1_position = Vecd(node_coordinates_[interface_node1_index][0], node_coordinates_[interface_node1_index][1], node_coordinates_[interface_node1_index][2]);
                        Vecd node2_position = Vecd(node_coordinates_[interface_node2_index][0], node_coordinates_[interface_node2_index][1], node_coordinates_[interface_node2_index][2]);
                        Vecd node3_position = Vecd(node_coordinates_[interface_node3_index][0], node_coordinates_[interface_node3_index][1], node_coordinates_[interface_node3_index][2]);
                        Vecd interface_area_vector1 = node1_position - node2_position;
                        Vecd interface_area_vector2 = node1_position - node3_position;
                        Vecd normal_vector = interface_area_vector1.cross(interface_area_vector2);
                        Real magnitude = normal_vector.norm();
                        Vecd normalized_normal_vector = normal_vector / magnitude;
                        Vecd node1_to_center_direction = particle_position - node1_position; 
                        if (node1_to_center_direction.dot(normalized_normal_vector) < 0)
                        {
                            normalized_normal_vector = -normalized_normal_vector; // vector pointing towards i
                        };
                        Real r_ij = 0; // we need r_ij to calculate the viscous force
                        // boundary_type == 2 means both of them are inside of fluid
                        if (boundary_type == 2)
                        {
                            r_ij = (particle_position - pos_n[index_j]).dot(normalized_normal_vector);
                        }
                        // boundary_type == 3 means fulid particle with wall boundary
                        if ((boundary_type == 3) | (boundary_type == 4) | (boundary_type == 9) | (boundary_type == 10) | (boundary_type == 36) 
                            | (boundary_type == 5) | (boundary_type == 7))
                        {
                            r_ij = node1_to_center_direction.dot(normalized_normal_vector) * 2.0;
                        }
                        Real dW_ij = (-0.5 * magnitude) / (2.0 * Vol_i * Vol_n[index_j]);
                        get_neighbor_relation(neighborhood, r_ij, dW_ij, normalized_normal_vector, index_j);
                    }
                }
            },
            ap);
    }
    //=================================================================================================//
    void InnerRelationInFVM_3d::updateConfiguration()
    {
        resetNeighborhoodCurrentSize();
        searchNeighborsByParticles(base_particles_.total_real_particles_,
                                base_particles_, inner_configuration_,
                                get_particle_index_, get_inner_neighbor_);
    }
    //=================================================================================================//
    GhostCreationFromMesh_3d::GhostCreationFromMesh_3d(RealBody& real_body, ANSYSMesh_3d& ansys_mesh, Ghost<ReserveSizeFactor>& ghost_boundary)
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
    void GhostCreationFromMesh_3d::addGhostParticleAndSetInConfiguration()
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
                    size_t node3_index = mesh_topology_[index_i][neighbor_index][4];
                    Vecd node1_position = Vecd(node_coordinates_[node1_index][0], node_coordinates_[node1_index][1], node_coordinates_[node1_index][2]);
                    Vecd node2_position = Vecd(node_coordinates_[node2_index][0], node_coordinates_[node2_index][1], node_coordinates_[node2_index][2]);
                    Vecd node3_position = Vecd(node_coordinates_[node3_index][0], node_coordinates_[node3_index][1], node_coordinates_[node3_index][2]);
                    Vecd ghost_particle_position = (1.0 / 3.0) * (node1_position + node2_position + node3_position);

                    mesh_topology_[index_i][neighbor_index][0] = ghost_particle_index + 1;
                    pos_[ghost_particle_index] = ghost_particle_position;
                    mutex_create_ghost_particle_.unlock();

                    mesh_topology_.resize(ghost_particle_index);
                    std::vector<std::vector<size_t>> new_element;

                    // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                    std::vector<size_t> sub_element1 = { index_i + 1, boundary_type, node1_index, node2_index, node3_index };
                    new_element.push_back(sub_element1);

                    // Add (correspon ding_index_i,boundary_type,node1_index,node2_index) to the new element
                    std::vector<size_t> sub_element2 = { index_i + 1, boundary_type, node1_index, node2_index, node3_index };
                    new_element.push_back(sub_element2);

                    // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                    std::vector<size_t> sub_element3 = { index_i + 1, boundary_type, node1_index, node2_index, node3_index };
                    new_element.push_back(sub_element3);

                    // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                    std::vector<size_t> sub_element4 = { index_i + 1, boundary_type, node1_index, node2_index, node3_index };
                    new_element.push_back(sub_element4);

                    // Add the new element to all_needed_data_from_mesh_file_
                    mesh_topology_.push_back(new_element);
                    // all_needed_data_from_mesh_file_[ghost_particle_index][0][0].push_back(size_t(0);

                    // creating the boundary files with ghost particle index
                    each_boundary_type_with_all_ghosts_index_[boundary_type].push_back(ghost_particle_index);

                    // creating the boundary files with contact real particle index
                    each_boundary_type_contact_real_index_[boundary_type].push_back(index_i);

                    // creating the boundary files with ghost eij
                    Vecd interface_area_vector1 = node2_position - node1_position;
                    Vecd interface_area_vector2 = node3_position - node1_position;
                    Vecd normal_vector = interface_area_vector1.cross(interface_area_vector2);
                    Real magnitude = normal_vector.norm();
                    //Real interface_area_size = interface_area_vector.norm();
                    Vecd normalized_normal_vector = normal_vector / magnitude;
                    // judge the direction
                    Vecd particle_position = pos_[index_i];
                    Vecd node1_to_center_direction = particle_position - node1_position;
                    if (node1_to_center_direction.dot(normalized_normal_vector) < 0)
                    {
                        normalized_normal_vector = -normalized_normal_vector;
                    };
                    each_boundary_type_with_all_ghosts_eij_[boundary_type].push_back(normalized_normal_vector);
                }
            }
        }
    };
    //=================================================================================================//
    BodyStatesRecordingInMeshToVtu::BodyStatesRecordingInMeshToVtu(SPHBody& body, ANSYSMesh_3d& ansys_mesh)
        : BodyStatesRecording(body), node_coordinates_(ansys_mesh.node_coordinates_),
        elements_nodes_connection_(ansys_mesh.elements_nodes_connection_), bounds_(body) {};
    //=================================================================================================//
    void BodyStatesRecordingInMeshToVtu::writeWithFileName(const std::string & sequence)
    {
        for (SPHBody* body : bodies_)
        {
            if (body->checkNewlyUpdated() && state_recording_)
            {
                std::string filefullpath = io_environment_.output_folder_ + "/SPHBody_" + body->getName() + "_" + sequence + ".vtu";
                if (fs::exists(filefullpath))
                {
                    fs::remove(filefullpath);
                }
                    std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

                    MeshFileHelpers::vtuFileHeader(out_file);
                    Real rangemax = 0.0;
                    MeshFileHelpers::vtuFileNodeCoordinates(out_file, node_coordinates_, elements_nodes_connection_, bounds_, rangemax);

                    MeshFileHelpers::vtuFileInformationKey(out_file, rangemax);

                    MeshFileHelpers::vtuFileCellConnectivity(out_file, elements_nodes_connection_, node_coordinates_);

                    MeshFileHelpers::vtuFileOffsets(out_file, elements_nodes_connection_);

                    MeshFileHelpers::vtuFileTypeOfCell(out_file, elements_nodes_connection_);

                    //write Particle data to vtu file
                    out_file << "<CellData>\n";
                    body->writeParticlesToVtuFile(out_file);
                    out_file << "</CellData>\n";
                    // Write VTU file footer
                    out_file << "</Piece>\n";
                    out_file << "</UnstructuredGrid>\n";
                    out_file << "</VTKFile>\n";
                    out_file.close();
            }
            body->setNotNewlyUpdated();
        }
    } 
    //=================================================================================================//
    BoundaryConditionSetupInFVM_3d::BoundaryConditionSetupInFVM_3d(BaseInnerRelationInFVM_3d& inner_relation, GhostCreationFromMesh_3d& ghost_creation) 
        : fluid_dynamics::FluidDataInner(inner_relation), rho_(particles_->rho_), Vol_(particles_->Vol_), mass_(particles_->mass_),
        p_(*particles_->getVariableByName<Real>("Pressure")),
        vel_(particles_->vel_), pos_(particles_->pos_), mom_(*particles_->getVariableByName<Vecd>("Momentum")),
        ghost_bound_(ghost_creation.ghost_bound_),
        each_boundary_type_with_all_ghosts_index_(ghost_creation.each_boundary_type_with_all_ghosts_index_),
        each_boundary_type_with_all_ghosts_eij_(ghost_creation.each_boundary_type_with_all_ghosts_eij_),
        each_boundary_type_contact_real_index_(ghost_creation.each_boundary_type_contact_real_index_) {};
    //=================================================================================================//
    void BoundaryConditionSetupInFVM_3d::resetBoundaryConditions()
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
                    case 3: // this refer to the different types of wall boundary condtions
                        applyNonSlipWallBoundary(ghost_index, index_i);
                        applyReflectiveWallBoundary(ghost_index, index_i, e_ij);
                        break;
                    case 5:
                        applyPressureOutletBC(ghost_index, index_i);
                        break;
                    case 7:
                        applySymmetryBoundary(ghost_index, index_i, e_ij);
                        break;
                    case 10:
                        applyGivenValueInletFlow(ghost_index);
                        applyVelocityInletFlow(ghost_index, index_i);
                        break;
                    case 36:
                        applyOutletBoundary(ghost_index, index_i);
                        break;
                    case 4:
                        applyTopBoundary(ghost_index, index_i);
                        break;
                    case 9:
                        applyFarFieldBoundary(ghost_index);
                        break;
                    }
                }
            }
        }
    }

    
}// namespace SPH