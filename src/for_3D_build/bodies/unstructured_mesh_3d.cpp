#include "mesh_helper.h"
#include "unstructured_mesh.h"

#include "base_particle_dynamics.h"
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
    Real ICEM = 0;
    std::ifstream mesh_file; /*!< \brief File object for the Ansys ASCII mesh file. */
    mesh_file.open(full_path);
    if (mesh_file.fail())
    {
        std::cout << "Error:Check if the file exists." << std::endl;
    }
    std::string text_line;
    /*Check mesh file generation software*/
    (getline(mesh_file, text_line));
    text_line.erase(0, 4);
    std::istringstream value(text_line);
    if (text_line.find("Created by", 0) != std::string::npos)
    {
        ICEM = 1;
        std::cout << "Reading 3D ICEM mesh." << std::endl;
    }

    if (ICEM == 1)
    {

        /*--- Read the dimension of the mesh ---*/
        size_t dimension(0);
        MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);

        /*--- Read the node data (index is starting from zero) ---*/
        size_t number_of_points(0);
        MeshFileHelpers::numberOfNodes(mesh_file, number_of_points, text_line);
        MeshFileHelpers::nodeCoordinates(mesh_file, node_coordinates_, text_line, dimension);

        size_t boundary_type(0);
        size_t number_of_elements(0);
        size_t mesh_type = 4;
        MeshFileHelpers::numberOfElements(mesh_file, number_of_elements, text_line);

        /*Preparing and initializing the data structure of mesh topology and element node connection*/
        MeshFileHelpers::dataStruct(mesh_topology_, elements_nodes_connection_, number_of_elements, mesh_type, dimension);

        /*--- find the elements lines ---*/
        while (getline(mesh_file, text_line))
        {
            if (text_line.find("(13", 0) != std::string::npos && text_line.find(")(", 0) != std::string::npos)
            {
                boundary_type = MeshFileHelpers::findBoundaryType(text_line, boundary_type);
                types_of_boundary_condition_.push_back(boundary_type);
                while (getline(mesh_file, text_line))
                {
                    if (text_line.find(")", 0) == std::string::npos)
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
            if (text_line.find("Zone Sections", 0) != std::string::npos)
                break;
            if (text_line.find(")") != std::string::npos)
                continue;
        }
        mesh_topology_.erase(mesh_topology_.begin());
    }
    else /*This section is for mesh files created from fluent*/
    {
        std::cout << "Reading 3D Fluent mesh." << std::endl;
        size_t dimension(0);
        MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);

        /*--- Read the node data (index is starting from zero) ---*/
        size_t number_of_points(0);
        MeshFileHelpers::numberOfNodesFluent(mesh_file, number_of_points, text_line);

        size_t boundary_type(0);
        size_t number_of_elements(0);
        size_t mesh_type = 4;

        MeshFileHelpers::numberOfElementsFluent(mesh_file, number_of_elements, text_line);
        MeshFileHelpers::dataStruct(mesh_topology_, elements_nodes_connection_, number_of_elements, mesh_type, dimension);
        MeshFileHelpers::nodeCoordinatesFluent(mesh_file, node_coordinates_, text_line, dimension);

        while (getline(mesh_file, text_line))
        {

            if (text_line.find("(13", 0) != std::string::npos && text_line.find(") (", 0) != std::string::npos)
            {
                boundary_type = MeshFileHelpers::findBoundaryType(text_line, boundary_type);
                types_of_boundary_condition_.push_back(boundary_type);
                while (getline(mesh_file, text_line))
                {

                    if (text_line.find(")", 0) == std::string::npos)
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
            if (text_line.find("(12", 0) != std::string::npos && text_line.find("))", 0) != std::string::npos)
                break;
            if (text_line.find(")") != std::string::npos)
                continue;
        }
        mesh_topology_.erase(mesh_topology_.begin());
    }
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
    StdVec<Real> all_data_of_distance_between_nodes;
    all_data_of_distance_between_nodes.resize(0);
    MeshFileHelpers::minimumDistance(all_data_of_distance_between_nodes, elements_volumes_, mesh_topology_, node_coordinates_);
    auto min_distance_iter = std::min_element(all_data_of_distance_between_nodes.begin(), all_data_of_distance_between_nodes.end());
    if (min_distance_iter != all_data_of_distance_between_nodes.end())
    {
        min_distance_between_nodes_ = *min_distance_iter;
    }
    else
    {
        std::cout << "The array of all distance between nodes is empty " << std::endl;
    }
}
//=================================================================================================//
void BaseInnerRelationInFVM::resetNeighborhoodCurrentSize()
{
    parallel_for(
        IndexRange(0, base_particles_.TotalRealParticles()),
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
        IndexRange(0, base_particles_.TotalRealParticles()),
        [&](const IndexRange &r)
        {
            for (size_t num = r.begin(); num != r.end(); ++num)
            {
                size_t index_i = get_particle_index(num);
                Vecd &particle_position = pos_[index_i];

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
                        normalized_normal_vector = -normalized_normal_vector;
                    };
                    Real r_ij = 0; // we need r_ij to calculate the viscous force
                    // boundary_type == 2 means both of them are inside of fluid
                    if (boundary_type == 2)
                    {
                        r_ij = (particle_position - pos_[index_j]).dot(normalized_normal_vector);
                    }
                    // this refer to the different types of wall boundary conditions
                    if ((boundary_type == 3) | (boundary_type == 4) | (boundary_type == 5) | (boundary_type == 7) | (boundary_type == 9) |
                        (boundary_type == 10) | (boundary_type == 36))
                    {
                        r_ij = node1_to_center_direction.dot(normalized_normal_vector) * 2.0;
                    }
                    Real dW_ij = (-0.5 * magnitude) / (2.0 * Vol_[index_i] * Vol_[index_j]);
                    get_neighbor_relation(neighborhood, r_ij, dW_ij, normalized_normal_vector, index_j);
                }
            }
        },
        ap);
}
//=================================================================================================//
void InnerRelationInFVM::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    searchNeighborsByParticles(base_particles_.TotalRealParticles(),
                               base_particles_, inner_configuration_,
                               get_particle_index_, get_inner_neighbor_);
}
//=================================================================================================//
} // namespace SPH
