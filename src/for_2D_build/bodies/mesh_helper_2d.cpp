#include "mesh_helper.h"

#include "base_particle_dynamics.h"
namespace SPH
{
//=================================================================================================//
void MeshFileHelpers::meshDimension(std::ifstream &mesh_file, size_t &dimension, std::string &text_line)
{
    while (getline(mesh_file, text_line))
    {
        text_line.erase(0, 1);
        text_line.erase(1);
        std::istringstream value(text_line);
        if (text_line.find("2", 0) != std::string::npos)
        {
            dimension = atoi(text_line.c_str());
            break;
        }
    }
}
//=================================================================================================//
void MeshFileHelpers::numberOfNodes(std::ifstream &mesh_file, size_t &number_of_points, std::string &text_line)
{
    while (getline(mesh_file, text_line))
    {
        text_line.erase(0, 1);
        std::string text1(text_line);
        text_line.erase(3);
        std::string text2(text_line);
        text_line.erase(2);
        if (atoi(text_line.c_str()) == 10 && text1.find("))", 0) != std::string::npos)
        {
            text1.erase(0, 8);
            Real last_position = text1.find_last_of(")");
            text1.erase(last_position - 5);
            number_of_points = stoi(text1, nullptr, 16);
            break;
        }
    }
}
//=================================================================================================//
void MeshFileHelpers::nodeCoordinates(std::ifstream &mesh_file, StdVec<Vecd> &node_coordinates_,
                                      std::string &text_line, size_t &dimension)
{
    while (getline(mesh_file, text_line))
    {
        if (text_line.find("(", 0) == std::string::npos && text_line.find("))", 0) == std::string::npos)
        {
            if (text_line.find(" ", 0) != std::string::npos)
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
                node_coordinates_.push_back(coordinate);
            }
        }
        if (text_line.find("))", 0) != std::string::npos)
        {
            break;
        }
    }
}
//=================================================================================================//
void MeshFileHelpers::numberOfElements(std::ifstream &mesh_file, size_t &number_of_elements, std::string &text_line)
{
    while (getline(mesh_file, text_line))
    {
        text_line.erase(0, 1);
        std::string text1(text_line);
        text_line.erase(3);
        std::string text2(text_line);
        text_line.erase(2);
        if (atoi(text_line.c_str()) == 12)
        {
            text1.erase(0, 8);
            Real last_position = text1.find_last_of(")");
            text1.erase(last_position - 5);
            number_of_elements = stoi(text1, nullptr, 16);
            break;
        }
    }
}
//=================================================================================================//
/*--- Initialize mesh topology ---*/
/** mesh_topology_
 * {[(neighbor_cell_index, bc_type, node1_of_face, node2_of_face), (....), (.....)], []..... }.
 * {inner_neighbor1, inner_neighbor2, ..... }.
 */
/*--- Initialize the number of elements ---*/
void MeshFileHelpers::dataStruct(StdVec<StdVec<StdVec<size_t>>> &mesh_topology_,
                                 StdVec<StdVec<size_t>> &elements_nodes_connection_,
                                 size_t number_of_elements, size_t mesh_type, size_t dimension)
{

    mesh_topology_.resize(number_of_elements + 1);
    for (std::size_t a = 0; a != number_of_elements + 1; ++a)
    {
        mesh_topology_[a].resize(mesh_type);
        for (std::size_t b = 0; b != mesh_topology_[a].size(); ++b)
        {
            mesh_topology_[a][b].resize(dimension + 2);
            for (std::size_t c = 0; c != mesh_topology_[a][b].size(); ++c)
            {
                mesh_topology_[a][b][c] = MaxUnsignedInt;
            }
        }
    }
    /*--- Initialize the number of elements ---*/
    elements_nodes_connection_.resize(number_of_elements + 1);
    mesh_topology_.resize(number_of_elements + 1);
    for (std::size_t element = 0; element != number_of_elements + 1; ++element)
    {
        elements_nodes_connection_[element].resize(3);
        for (std::size_t node = 0; node != elements_nodes_connection_[element].size(); ++node)
        {
            elements_nodes_connection_[element][node] = MaxUnsignedInt;
        }
    }
}
//=================================================================================================//
size_t MeshFileHelpers::findBoundaryType(std::string &text_line, size_t boundary_type)
{
    /*--- Read the elements of the problem (13 (id f1 f2 type 0) (---*/
    /*--- Read the elements of the problem ---*/
    /** boundary condition types
     * bc-type==2, interior boundary condition.
     * bc-type==3, wall boundary condition.
     * bc-type==9, pressure-far-field boundary condition.
     * Note that Cell 0 means boundary condition.
     * mesh_type==3, unstructured mesh.
     * mesh_type==4, structured mesh.
     * mesh_topology_
     * {[(neighbor_cell_index, bc_type, node1_of_face, node2_of_face), (....), (.....)], []..... }.
     * {inner_neighbor1, inner_neighbor2, ..... }.
     */
    /*--- find the type of boundary condition ---*/
    size_t position = text_line.find(")", 0);
    text_line = text_line.erase(0, position - 4);
    text_line = text_line.erase(2);
    boundary_type = std::stoi(text_line, nullptr, 16);
    return boundary_type;
}
//=================================================================================================//
Vecd MeshFileHelpers::nodeIndex(std::string &text_line)
{
    /*--- find the node1 between two cells ---*/
    std::string node1_string_copy = text_line;
    size_t first_divide_position = text_line.find_first_of(" ", 0);
    std::string node1_index_string = node1_string_copy.erase(first_divide_position);
    size_t node1_index_decimal = stoi(node1_index_string, nullptr, 16) - 1;

    /*--- find the node2 between two cells---*/
    std::string node2_string = text_line;
    node2_string = node2_string.erase(0, first_divide_position + 1);
    size_t second_divide_position = node2_string.find_first_of(" ", 0);
    node2_string.erase(second_divide_position);
    size_t node2_index_decimal = stoi(node2_string, nullptr, 16) - 1;
    Vecd nodes = Vecd(node1_index_decimal, node2_index_decimal);
    return nodes;
}
//=================================================================================================//
Vec2d MeshFileHelpers::cellIndex(std::string &text_line)
{

    /*--- find the cell1---*/
    std::string cell1_string = text_line;
    size_t first_divide_position = text_line.find_first_of(" ", 0);
    cell1_string = cell1_string.erase(0, first_divide_position + 1);
    size_t second_divide_position = cell1_string.find_first_of(" ", 0);
    cell1_string = cell1_string.erase(0, second_divide_position + 1);
    size_t third_divide_position = cell1_string.find_first_of(" ", 0);
    cell1_string.erase(third_divide_position);
    size_t cell1_index_decimal = std::stoi(cell1_string, nullptr, 16);

    /*--- find the cell2---*/
    std::string cell2_string = text_line;
    cell2_string = cell2_string.erase(0, first_divide_position + 1);
    cell2_string = cell2_string.erase(0, second_divide_position + 1);
    cell2_string.erase(0, third_divide_position + 1);
    size_t cell2_index_decimal = stoi(cell2_string, nullptr, 16);
    Vecd cells = Vecd(cell1_index_decimal, cell2_index_decimal);
    return cells;
}
//=================================================================================================//
void MeshFileHelpers::updateElementsNodesConnection(StdVec<StdVec<size_t>> &elements_nodes_connection_,
                                                    Vecd nodes, Vec2d cells,
                                                    bool &check_neighbor_cell1, bool &check_neighbor_cell2)
{

    /*--- build up connection with element and nodes only---*/
    for (int node = 0; node != nodes.size(); ++node)
    {
        if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != nodes[node] &&
            elements_nodes_connection_[cells[check_neighbor_cell2]][1] != nodes[node] &&
            elements_nodes_connection_[cells[check_neighbor_cell2]][2] != nodes[node])
        {
            if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) &&
                (elements_nodes_connection_[cells[check_neighbor_cell2]][1] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1)) &&
                elements_nodes_connection_[cells[check_neighbor_cell2]][2] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
            {
                elements_nodes_connection_[cells[check_neighbor_cell2]][2] = nodes[node];
            }
            if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) &&
                elements_nodes_connection_[cells[check_neighbor_cell2]][1] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
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
}
//=================================================================================================//
void MeshFileHelpers::updateCellLists(StdVec<StdVec<StdVec<size_t>>> &mesh_topology_,
                                      StdVec<StdVec<size_t>> &elements_nodes_connection_, Vecd nodes,
                                      Vec2d cells, bool &check_neighbor_cell1, bool &check_neighbor_cell2, size_t boundary_type)
{
    /*--- build up all connection data with element and neighbor and nodes---*/
    if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] &&
        mesh_topology_[cells[check_neighbor_cell2]][1][0] != cells[check_neighbor_cell1] &&
        mesh_topology_[cells[check_neighbor_cell2]][2][0] != cells[check_neighbor_cell1])
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
            return;
        }
        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) &&
            mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            mesh_topology_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
            mesh_topology_[cells[check_neighbor_cell2]][1][1] = boundary_type;
            mesh_topology_[cells[check_neighbor_cell2]][1][2] = nodes[0];
            mesh_topology_[cells[check_neighbor_cell2]][1][3] = nodes[1];
            check_neighbor_cell1 = false;
            check_neighbor_cell2 = true;
            return;
        }
        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) &&
            mesh_topology_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) &&
            mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            mesh_topology_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
            mesh_topology_[cells[check_neighbor_cell2]][2][1] = boundary_type;
            mesh_topology_[cells[check_neighbor_cell2]][2][2] = nodes[0];
            mesh_topology_[cells[check_neighbor_cell2]][2][3] = nodes[1];
            check_neighbor_cell1 = false;
            check_neighbor_cell2 = true;
            return;
        }
    }
}
//=================================================================================================//
void MeshFileHelpers::updateBoundaryCellLists(StdVec<StdVec<StdVec<size_t>>> &mesh_topology_,
                                              StdVec<StdVec<size_t>> &elements_nodes_connection_,
                                              Vecd nodes, Vec2d cells, bool &check_neighbor_cell1, bool &check_neighbor_cell2, size_t boundary_type)
{
    if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] &&
        mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) &&
        mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
    {
        mesh_topology_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
        mesh_topology_[cells[check_neighbor_cell2]][2][1] = boundary_type;
        mesh_topology_[cells[check_neighbor_cell2]][2][2] = nodes[0];
        mesh_topology_[cells[check_neighbor_cell2]][2][3] = nodes[1];
        check_neighbor_cell1 = false;
        check_neighbor_cell2 = true;
        return;
    }
    else
    {
        check_neighbor_cell1 = false;
        check_neighbor_cell2 = true;
        return;
    }
}
//=================================================================================================//
void MeshFileHelpers::cellCenterCoordinates(StdVec<StdVec<size_t>> &elements_nodes_connection_,
                                            std::size_t &element, StdVec<Vecd> &node_coordinates_,
                                            StdVec<Vecd> &elements_centroids_, Vecd &center_coordinate)
{

    for (std::size_t node = 0; node != elements_nodes_connection_[element].size(); ++node)
    {
        center_coordinate += node_coordinates_[elements_nodes_connection_[element][node]] / 3.0;
    }
    elements_centroids_[element] = center_coordinate;
}
//=================================================================================================//
void MeshFileHelpers::elementVolume(StdVec<StdVec<size_t>> &elements_nodes_connection_, std::size_t &element,
                                    StdVec<Vecd> &node_coordinates_, StdVec<Real> &elements_volumes_)
{
    // calculating each volume of element
    // get nodes position
    Vecd node1_coordinate = node_coordinates_[elements_nodes_connection_[element][0]];
    Vecd node2_coordinate = node_coordinates_[elements_nodes_connection_[element][1]];
    Vecd node3_coordinate = node_coordinates_[elements_nodes_connection_[element][2]];
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
//=================================================================================================//
void MeshFileHelpers::minimumDistance(StdVec<Real> &all_data_of_distance_between_nodes,
                                      StdVec<Real> &elements_volumes_, StdVec<StdVec<StdVec<size_t>>> &mesh_topology_,
                                      StdVec<Vecd> &node_coordinates_)
{

    for (size_t element_index = 0; element_index != elements_volumes_.size(); ++element_index)
    {

        for (std::size_t neighbor = 0; neighbor != mesh_topology_[element_index].size(); ++neighbor)
        {

            size_t interface_node1_index = mesh_topology_[element_index][neighbor][2];
            size_t interface_node2_index = mesh_topology_[element_index][neighbor][3];
            Vecd node1_position = node_coordinates_[interface_node1_index];
            Vecd node2_position = node_coordinates_[interface_node2_index];
            Vecd interface_area_vector = node1_position - node2_position;
            Real interface_area_size = interface_area_vector.norm();
            all_data_of_distance_between_nodes.push_back(interface_area_size);
        }
    }
}
//=================================================================================================//
} // namespace SPH
