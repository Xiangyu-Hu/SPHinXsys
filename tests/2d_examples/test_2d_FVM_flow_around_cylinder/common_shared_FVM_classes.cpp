
#include "common_shared_FVM_classes.h"
namespace SPH
{
//=================================================================================================//
ANSYSMesh::ANSYSMesh(const std::string &full_path)
{
    getDataFromMeshFile(full_path);
    getElementCenterCoordinates();
    gerMinimumDistanceBetweenNodes();
}
//=================================================================================================//
void ANSYSMesh::readNodeCoordinate(const std::string &text_line, StdLargeVec<Vec2d> &node_coordinates)
{
    size_t divide_position = text_line.find_first_of(" ");
    std::string x_part = text_line;
    std::string y_part = text_line;
    std::string x_coordinate_string = x_part.erase(divide_position);
    std::string y_coordinate_string = y_part.erase(0, divide_position);

    Vec2d coordinate = Vec2d::Zero();
    std::istringstream stream_x, stream_y;
    stream_x.str(x_coordinate_string);
    stream_y.str(y_coordinate_string);
    stream_x >> coordinate[0];
    stream_y >> coordinate[1];
    node_coordinates.push_back(coordinate);
}
//=================================================================================================//
void ANSYSMesh::readNodeCoordinate(const std::string &text_line, StdLargeVec<Vec3d> &node_coordinates)
{
    size_t first_divide_position = text_line.find_first_of(" ");
    size_t last_divide_position = text_line.find_last_of(" ");
    string x_part = text_line;
    string y_part = text_line;
    string z_part = text_line;
    string x_coordinate_string = x_part.erase(first_divide_position);
    string y_coordinate_string = y_part.erase(last_divide_position);
    y_coordinate_string = y_coordinate_string.erase(0, first_divide_position);
    string z_coordinate_string = z_part.erase(0, last_divide_position);

    Vec3d coordinate = Vec3d::Zero();
    istringstream stream_x, stream_y, stream_z;
    stream_x.str(x_coordinate_string);
    stream_y.str(y_coordinate_string);
    stream_z.str(z_coordinate_string);
    stream_x >> coordinate[0];
    stream_y >> coordinate[1];
    stream_z >> coordinate[2];
    node_coordinates.push_back(coordinate);
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
    int dimension(0);
    while (getline(mesh_file, text_line))
    {
        text_line.erase(0, 1);
        text_line.erase(1);
        istringstream value(text_line);
        if (text_line.find("2", 0) != string::npos)
        {
            dimension = atoi(text_line.c_str());
            break;
        }
    }
    /*--- Check dimension ---*/
    if (dimension != Dimensions)
    {
        cout << "Error:the dimension of problem does not match input mesh." << endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    }
    /*--- Read the total number node points ---*/
    size_t number_of_points(0);
    while (getline(mesh_file, text_line))
    {
        text_line.erase(0, 1);
        string text1(text_line);
        text_line.erase(3);
        string text2(text_line);
        text_line.erase(2);
        if (atoi(text_line.c_str()) == 10 && text1.find("))", 0) != string::npos)
        {
            text1.erase(0, 8);
            Real last_position = text1.find_last_of(")");
            text1.erase(last_position - 5);
            number_of_points = stoi(text1, nullptr, 16);
            break;
        }
    };
    /*--- Read the node coordinates ---*/
    while (getline(mesh_file, text_line))
    {
        if (text_line.find("(", 0) == string::npos && text_line.find("))", 0) == string::npos)
        {
            if (text_line.find(" ", 0) != string::npos)
            {
                readNodeCoordinate(text_line, node_coordinates_);
            }
        }
        if (text_line.find("))", 0) != string::npos)
        {
            break;
        }
    }
    /*--- Check number of node points ---*/
    if (node_coordinates_.size() != number_of_points)
    {
        cout << "Error:Total number of node points does not match data!" << endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
    }

    /*--- Read the elements of the problem ---*/
    /** boundary condition types
     * bc-type==2, interior boundary condition.
     * bc-type==3, wall boundary condition.
     * bc-type==9, pressure-far-field boundary condition.
     * Note that Cell0 means boundary condition.
     * mesh_type==3, unstructured mesh.
     * mesh_type==4, structured mesh.
     * mesh_topology_
     * {[(neighbor_cell_index, bc_type, node1_of_face, node2_of_face), (....), (.....)], []..... }.
     * {inner_neighbor1, inner_neighbor2, ..... }.
     */
    size_t boundary_type(0);
    size_t number_of_elements(0);
    size_t mesh_type = 3;
    while (getline(mesh_file, text_line))
    {
        text_line.erase(0, 1);
        string text1(text_line);
        text_line.erase(3);
        string text2(text_line);
        text_line.erase(2);
        if (atoi(text_line.c_str()) == 12)
        {
            text1.erase(0, 8);
            Real last_position = text1.find_last_of(")");
            text1.erase(last_position - 5);
            number_of_elements = stoi(text1, nullptr, 16);
            break;
        }
    };
    mesh_topology_.resize(number_of_elements + 1);
    for (std::size_t a = 0; a != number_of_elements + 1; ++a)
    {
        mesh_topology_[a].resize(mesh_type);
        for (std::vector<std::vector<long unsigned int>>::size_type b = 0; b != mesh_topology_[a].size(); ++b)
        {
            mesh_topology_[a][b].resize(dimension + 2);
            for (std::vector<long unsigned int>::size_type c = 0; c != mesh_topology_[a][b].size(); ++c)
            {
                mesh_topology_[a][b][c] = -1;
            }
        }
    }
    /*--- reinitialize the number of elements ---*/
    elements_nodes_connection_.resize(number_of_elements + 1);
    elements_neighbors_connection_.resize(number_of_elements + 1);
    mesh_topology_.resize(number_of_elements + 1);
    for (std::size_t element = 0; element != number_of_elements + 1; ++element)
    {
        elements_nodes_connection_[element].resize(3);
        for (std::vector<long unsigned int>::size_type node = 0; node != elements_nodes_connection_[element].size(); ++node)
        {
            elements_nodes_connection_[element][node] = -1;
        }
    }

    /*--- find the elements lines ---*/
    while (getline(mesh_file, text_line))
    {
        if (text_line.find("(13", 0) != string::npos && text_line.find(")(", 0) != string::npos)
        {
            /*--- find the type of boundary condition ---*/
            size_t position = text_line.find(")", 0);
            text_line = text_line.erase(0, position - 4);
            text_line = text_line.erase(2);
            boundary_type = stoi(text_line, nullptr, 16);
            types_of_boundary_condition_.push_back(boundary_type);
            while (getline(mesh_file, text_line))
            {
                if (text_line.find(")", 0) == string::npos)
                {
                    /*--- find the node1 between two cells ---*/
                    string node1_string_copy = text_line;
                    size_t first_divide_position = text_line.find_first_of(" ", 0);
                    string node1_index_string = node1_string_copy.erase(first_divide_position);
                    size_t node1_index_decimal = stoi(node1_index_string, nullptr, 16) - 1;

                    /*--- find the node2 between two cells---*/
                    string node2_string = text_line;
                    node2_string = node2_string.erase(0, first_divide_position + 1);
                    size_t second_divide_position = node2_string.find_first_of(" ", 0);
                    node2_string.erase(second_divide_position);
                    size_t node2_index_decimal = stoi(node2_string, nullptr, 16) - 1;
                    Vecd nodes = Vecd(node1_index_decimal, node2_index_decimal);

                    /*--- find the cell1---*/
                    string cell1_string = text_line;
                    cell1_string = cell1_string.erase(0, first_divide_position + 1);
                    cell1_string = cell1_string.erase(0, second_divide_position + 1);
                    size_t third_divide_position = cell1_string.find_first_of(" ", 0);
                    cell1_string.erase(third_divide_position);
                    size_t cell1_index_decimal = stoi(cell1_string, nullptr, 16);

                    /*--- find the cell2---*/
                    string cell2_string = text_line;
                    cell2_string = cell2_string.erase(0, first_divide_position + 1);
                    cell2_string = cell2_string.erase(0, second_divide_position + 1);
                    cell2_string.erase(0, third_divide_position + 1);
                    size_t cell2_index_decimal = stoi(cell2_string, nullptr, 16);
                    Vecd cells = Vecd(cell1_index_decimal, cell2_index_decimal);

                    /*--- build up all topology---*/
                    bool check_neighbor_cell1 = 1;
                    bool check_neighbor_cell2 = 0;
                    for (int cell1_cell2 = 0; cell1_cell2 != cells.size(); ++cell1_cell2)
                    {
                        while (true)
                        {
                            if (mesh_type == 3)
                            {
                                if (cells[check_neighbor_cell2] != 0)
                                {
                                    /*--- build up connection with element and nodes only---*/
                                    for (int node = 0; node != nodes.size(); ++node)
                                    {
                                        if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != nodes[node] && elements_nodes_connection_[cells[check_neighbor_cell2]][1] != nodes[node] && elements_nodes_connection_[cells[check_neighbor_cell2]][2] != nodes[node])
                                        {
                                            if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && (elements_nodes_connection_[cells[check_neighbor_cell2]][1] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1)) && elements_nodes_connection_[cells[check_neighbor_cell2]][2] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                            {
                                                elements_nodes_connection_[cells[check_neighbor_cell2]][2] = nodes[node];
                                            }
                                            if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && elements_nodes_connection_[cells[check_neighbor_cell2]][1] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                            {
                                                elements_nodes_connection_[cells[check_neighbor_cell2]][1] = nodes[node];
                                            }
                                            if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                            {
                                                elements_nodes_connection_[cells[check_neighbor_cell2]][0] = nodes[node];
                                            }
                                        }
                                        else
                                            continue;
                                    }
                                    /*--- build up all connection data with element and neighbor and nodes---*/
                                    if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] && mesh_topology_[cells[check_neighbor_cell2]][1][0] != cells[check_neighbor_cell1] && mesh_topology_[cells[check_neighbor_cell2]][2][0] != cells[check_neighbor_cell1])
                                    {
                                        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            /*--- inner neighbor index---*/
                                            mesh_topology_[cells[check_neighbor_cell2]][0][0] = cells[check_neighbor_cell1];
                                            /*--- boundary type---*/
                                            mesh_topology_[cells[check_neighbor_cell2]][0][1] = boundary_type;
                                            /*--- nodes of a face---*/
                                            mesh_topology_[cells[check_neighbor_cell2]][0][2] = nodes[0];
                                            mesh_topology_[cells[check_neighbor_cell2]][0][3] = nodes[1];

                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            mesh_topology_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
                                            mesh_topology_[cells[check_neighbor_cell2]][1][1] = boundary_type;
                                            mesh_topology_[cells[check_neighbor_cell2]][1][2] = nodes[0];
                                            mesh_topology_[cells[check_neighbor_cell2]][1][3] = nodes[1];
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            mesh_topology_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
                                            mesh_topology_[cells[check_neighbor_cell2]][2][1] = boundary_type;
                                            mesh_topology_[cells[check_neighbor_cell2]][2][2] = nodes[0];
                                            mesh_topology_[cells[check_neighbor_cell2]][2][3] = nodes[1];
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                    }
                                    if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] && mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                    {
                                        mesh_topology_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
                                        mesh_topology_[cells[check_neighbor_cell2]][2][1] = boundary_type;
                                        mesh_topology_[cells[check_neighbor_cell2]][2][2] = nodes[0];
                                        mesh_topology_[cells[check_neighbor_cell2]][2][3] = nodes[1];
                                        check_neighbor_cell1 = false;
                                        check_neighbor_cell2 = true;
                                        break;
                                    }
                                    else
                                    {
                                        check_neighbor_cell1 = false;
                                        check_neighbor_cell2 = true;
                                        break;
                                    }
                                }
                                if (cells[check_neighbor_cell2] == 0)
                                {
                                    break;
                                }
                            }
                            if (mesh_type == 4)
                            {
                                if (cells[check_neighbor_cell2] != 0)
                                {
                                    if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] && mesh_topology_[cells[check_neighbor_cell2]][1][0] != cells[check_neighbor_cell1] && mesh_topology_[cells[check_neighbor_cell2]][2][0] != cells[check_neighbor_cell1] && mesh_topology_[cells[check_neighbor_cell2]][3][0] != cells[check_neighbor_cell1])
                                    {
                                        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            /*--- inner neighbor index---*/
                                            mesh_topology_[cells[check_neighbor_cell2]][0][0] = cells[check_neighbor_cell1];
                                            /*--- boundary type---*/
                                            mesh_topology_[cells[check_neighbor_cell2]][0][1] = boundary_type;
                                            /*--- nodes of a face---*/
                                            mesh_topology_[cells[check_neighbor_cell2]][0][2] = nodes[0];
                                            mesh_topology_[cells[check_neighbor_cell2]][0][3] = nodes[1];

                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            mesh_topology_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
                                            mesh_topology_[cells[check_neighbor_cell2]][1][1] = boundary_type;
                                            mesh_topology_[cells[check_neighbor_cell2]][1][2] = nodes[0];
                                            mesh_topology_[cells[check_neighbor_cell2]][1][3] = nodes[1];
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            mesh_topology_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
                                            mesh_topology_[cells[check_neighbor_cell2]][2][1] = boundary_type;
                                            mesh_topology_[cells[check_neighbor_cell2]][2][2] = nodes[0];
                                            mesh_topology_[cells[check_neighbor_cell2]][2][3] = nodes[1];
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][2][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            mesh_topology_[cells[check_neighbor_cell2]][3][0] = cells[check_neighbor_cell1];
                                            mesh_topology_[cells[check_neighbor_cell2]][3][1] = boundary_type;
                                            mesh_topology_[cells[check_neighbor_cell2]][3][2] = nodes[0];
                                            mesh_topology_[cells[check_neighbor_cell2]][3][3] = nodes[1];
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                        // if it is three dimension
                                    }
                                    else
                                    {
                                        check_neighbor_cell1 = false;
                                        check_neighbor_cell2 = true;
                                        break;
                                    }
                                }
                                if (cells[check_neighbor_cell2] == 0)
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
    elements_centroids.resize(elements_nodes_connection_.size());
    elements_volumes_.resize(elements_nodes_connection_.size());
    for (std::vector<std::vector<long unsigned int>>::size_type element = 1; element != elements_nodes_connection_.size(); ++element)
    {
        Vecd center_coordinate = Vecd::Zero();
        for (std::vector<long unsigned int>::size_type node = 0; node != elements_nodes_connection_[element].size(); ++node)
        {
            center_coordinate += Vecd(node_coordinates_[elements_nodes_connection_[element][node]][0] / 3.0,
                                      node_coordinates_[elements_nodes_connection_[element][node]][1] / 3.0);
        }
        elements_centroids[element] = center_coordinate;

        // calculating each volume of element
        // get nodes position
        Vec3d nodes = Vec3d(elements_nodes_connection_[element][0], elements_nodes_connection_[element][1], elements_nodes_connection_[element][2]);
        Vecd node1_coordinate = Vecd(node_coordinates_[nodes[0]][0], node_coordinates_[nodes[0]][1]);
        Vecd node2_coordinate = Vecd(node_coordinates_[nodes[1]][0], node_coordinates_[nodes[1]][1]);
        Vecd node3_coordinate = Vecd(node_coordinates_[nodes[2]][0], node_coordinates_[nodes[2]][1]);
        // get each line length
        Real first_side_length = (node1_coordinate - node2_coordinate).norm();
        Real second_side_length = (node1_coordinate - node3_coordinate).norm();
        Real third_side_length = (node2_coordinate - node3_coordinate).norm();
        // half perimeter
        Real half_perimeter = (first_side_length + second_side_length + third_side_length) / 2.0;
        // get element volume
        Real element_volume =
            pow(half_perimeter * (half_perimeter - first_side_length) * (half_perimeter - second_side_length) * (half_perimeter - third_side_length), 0.5);
        elements_volumes_[element] = element_volume;
    }
    elements_volumes_.erase(elements_volumes_.begin());
    elements_centroids.erase(elements_centroids.begin());
    elements_nodes_connection_.erase(elements_nodes_connection_.begin());
}
//=================================================================================================//
void ANSYSMesh::gerMinimumDistanceBetweenNodes()
{
    vector<Real> all_data_of_distance_between_nodes;
    all_data_of_distance_between_nodes.resize(0);
    for (size_t element_index = 0; element_index != elements_volumes_.size(); ++element_index)
    {
        for (std::vector<std::vector<long unsigned int>>::size_type neighbor = 0; neighbor != mesh_topology_[element_index].size(); ++neighbor)
        {
            size_t interface_node1_index = mesh_topology_[element_index][neighbor][2];
            size_t interface_node2_index = mesh_topology_[element_index][neighbor][3];
            Vecd node1_position = Vecd(node_coordinates_[interface_node1_index][0], node_coordinates_[interface_node1_index][1]);
            Vecd node2_position = Vecd(node_coordinates_[interface_node2_index][0], node_coordinates_[interface_node2_index][1]);
            Vecd interface_area_vector = node1_position - node2_position;
            Real interface_area_size = interface_area_vector.norm();
            all_data_of_distance_between_nodes.push_back(interface_area_size);
        }
    }
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
        IndexRange(0, base_particles_.total_real_particles_ + base_particles_.total_ghost_particles_),
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
    resizeConfiguration();
};
//=================================================================================================//
void BaseInnerRelationInFVM::resizeConfiguration()
{
    size_t updated_size = base_particles_.real_particles_bound_ + base_particles_.total_ghost_particles_;
    inner_configuration_.resize(updated_size, Neighborhood());
}
//=================================================================================================//
ParticleGeneratorInFVM::ParticleGeneratorInFVM(SPHBody &sph_body, ANSYSMesh &ansys_mesh)
    : ParticleGenerator(sph_body), elements_centroids(ansys_mesh.elements_centroids),
      elements_volumes_(ansys_mesh.elements_volumes_) {}
//=================================================================================================//
void ParticleGeneratorInFVM::initializeGeometricVariables()
{
    for (size_t particle_index = 0; particle_index != elements_centroids.size(); ++particle_index)
    {
        initializePositionAndVolumetricMeasure(elements_centroids[particle_index], elements_volumes_[particle_index]);
    }
}
//=================================================================================================//
void NeighborBuilderInFVM::createRelation(Neighborhood &neighborhood, Real &distance,
                                          Real &dW_ijV_j, Vecd &interface_normal_direction, size_t j_index) const
{
    neighborhood.j_.push_back(j_index);
    neighborhood.r_ij_.push_back(distance);
    neighborhood.e_ij_.push_back(interface_normal_direction);
    neighborhood.dW_ijV_j_.push_back(dW_ijV_j);
    neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborBuilderInFVM::initializeRelation(Neighborhood &neighborhood, Real &distance,
                                              Real &dW_ijV_j, Vecd &interface_normal_direction, size_t j_index) const
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = j_index;
    neighborhood.dW_ijV_j_[current_size] = dW_ijV_j;
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
        IndexRange(0, base_particles_.total_real_particles_ + base_particles_.total_ghost_particles_),
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
                    // boundary_type == 3 means fluid particle with wall boundary
                    if ((boundary_type == 3) | (boundary_type == 4) | (boundary_type == 9) | (boundary_type == 10) | (boundary_type == 36))
                    {
                        r_ij = node1_to_center_direction.dot(normal_vector) * 2.0;
                    }
                    Real dW_ijV_j = -interface_area_size / (2.0 * Vol_i);
                    get_neighbor_relation(neighborhood, r_ij, dW_ijV_j, normal_vector, index_j);
                }
            }
        },
        ap);
}
//=================================================================================================//
void InnerRelationInFVM::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    searchNeighborsByParticles(base_particles_.total_real_particles_ + base_particles_.total_ghost_particles_,
                               base_particles_, inner_configuration_,
                               get_particle_index_, get_inner_neighbor_);
}
//=================================================================================================//
GhostCreationFromMesh::GhostCreationFromMesh(RealBody &real_body, ANSYSMesh &ansys_mesh)
    : GeneralDataDelegateSimple(real_body),
      node_coordinates_(ansys_mesh.node_coordinates_),
      mesh_topology_(ansys_mesh.mesh_topology_),
      pos_(particles_->pos_), Vol_(particles_->Vol_),
      total_ghost_particles_(particles_->total_ghost_particles_),
      real_particles_bound_(particles_->real_particles_bound_)
{
    each_boundary_type_with_all_ghosts_index_.resize(50);
    each_boundary_type_with_all_ghosts_eij_.resize(50);
    each_boundary_type_contact_real_index_.resize(50);
    ghost_particles_.resize(1);
    addGhostParticleAndSetInConfiguration();
}
//=================================================================================================//
void GhostCreationFromMesh::addGhostParticleAndSetInConfiguration()
{
    for (size_t i = 0; i != ghost_particles_.size(); ++i)
        ghost_particles_[i].clear();

    for (size_t index_i = 0; index_i != real_particles_bound_; ++index_i)
    {
        for (size_t neighbor_index = 0; neighbor_index != mesh_topology_[index_i].size(); ++neighbor_index)
        {
            size_t boundary_type = mesh_topology_[index_i][neighbor_index][1];
            if (mesh_topology_[index_i][neighbor_index][1] != 2)
            {
                mutex_create_ghost_particle_.lock();
                size_t ghost_particle_index = particles_->insertAGhostParticle(index_i);
                size_t node1_index = mesh_topology_[index_i][neighbor_index][2];
                size_t node2_index = mesh_topology_[index_i][neighbor_index][3];
                Vecd node1_position = node_coordinates_[node1_index];
                Vecd node2_position = node_coordinates_[node2_index];
                Vecd ghost_particle_position = 0.5 * (node1_position + node2_position);

                mesh_topology_[index_i][neighbor_index][0] = ghost_particle_index + 1;
                ghost_particles_[0].push_back(ghost_particle_index);
                pos_[ghost_particle_index] = ghost_particle_position;
                mutex_create_ghost_particle_.unlock();

                mesh_topology_.resize(ghost_particle_index);
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
                // mesh_topology_[ghost_particle_index][0][0].push_back(size_t(0);

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
BodyStatesRecordingInMeshToVtp::
    BodyStatesRecordingInMeshToVtp(IOEnvironment &io_environment, SPHBody &body, ANSYSMesh &ansys_mesh)
    : BodyStatesRecording(io_environment, body), node_coordinates_(ansys_mesh.node_coordinates_),
      elements_nodes_connection_(ansys_mesh.elements_nodes_connection_){};
//=================================================================================================//
void BodyStatesRecordingInMeshToVtp::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated() && state_recording_)
        {
            // TODO: we can short the file name by without using SPHBody
            std::string filefullpath = io_environment_.output_folder_ + "/SPHBody_" + body->getName() + "_" + sequence + ".vtp";
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
//=============================================================================================//
} // namespace SPH
  //=================================================================================================//