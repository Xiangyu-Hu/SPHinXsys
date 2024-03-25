#include "MeshHelper.h"
namespace SPH
{

    void MeshFileHelpers::meshDimension(ifstream& mesh_file, size_t& dimension, string& text_line)
    {
         while (getline(mesh_file, text_line))
            {
                (getline(mesh_file, text_line));
                text_line.erase(0, 1);
                text_line.erase(0, 2);
                istringstream value(text_line);
                if (text_line.find("3", 0) != string::npos)
                {
                    dimension = atoi(text_line.c_str());
                    break;
                }
        
            }

    }

    void MeshFileHelpers::numberofNodes(ifstream& mesh_file, size_t& number_of_points, string& text_line)
    {

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
            }
    }

    void MeshFileHelpers::nodeCoordinates(ifstream& mesh_file, size_t& node_index,
                                          vector<vector<Real>>& node_coordinates_, string& text_line, size_t& dimension, vector<vector<Real>>& point_coordinates)
    {
        while (getline(mesh_file, text_line))
            {
                if (text_line.find("(", 0) == string::npos && text_line.find("))", 0) == string::npos)
                {
                    Real Coords[3] = { 0.0, 0.0, 0.0 };
                    if (text_line.find(" ", 0) != string::npos)
                    {
                        if (dimension == 3)
                        {
                            size_t first_devide_position = text_line.find_first_of(" ");
                            size_t last_devide_position = text_line.find_last_of(" ");
                            string x_part = text_line;
                            string y_part = text_line;
                            string z_part = text_line;
                            string x_coordinate_string = x_part.erase(first_devide_position);
                            string y_coordinate_string = y_part.erase(last_devide_position);
                            y_coordinate_string = y_coordinate_string.erase(0, first_devide_position);
                            string z_coordinate_string = z_part.erase(0, last_devide_position);
                            istringstream streamx, streamy, streamz;
                            streamx.str(x_coordinate_string);
                            streamy.str(y_coordinate_string);
                            streamz.str(z_coordinate_string);
                            streamx >> Coords[0];
                            streamy >> Coords[1];
                            streamz >> Coords[2];
                            
                            node_coordinates_[node_index][0] = Coords[0];
                            node_coordinates_[node_index][1] = Coords[1];
                            node_coordinates_[node_index][2] = Coords[2];
                            ++node_index;
                        }
                        for (std::size_t iDim = 0; iDim != dimension; ++iDim)
                        {
                            point_coordinates[iDim].push_back(Coords[iDim]);
                        }
                    }
                }
                if (text_line.find("))", 0) != string::npos)
                {
                    break;
                }
            }
    }

    void MeshFileHelpers::numberofElements(ifstream& mesh_file, size_t& number_of_elements, string& text_line)
    {

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
            }
    }

    void MeshFileHelpers::dataStruct(vector<vector<vector<size_t>>>& mesh_topology_, vector<vector<size_t>>& elements_nodes_connection_,
                                    size_t number_of_elements, size_t mesh_type, size_t dimension, StdLargeVec<Vecd>& elements_neighbors_connection_)
    {

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
        //cell_lists_.resize(number_of_elements + 1);
        for (std::size_t element = 0; element != number_of_elements + 1; ++element)
        {
            elements_nodes_connection_[element].resize(4);
            for (std::vector<long unsigned int>::size_type node = 0; node != elements_nodes_connection_[element].size(); ++node)
            {
                elements_nodes_connection_[element][node] = -1;
            }
        }
    }

    size_t MeshFileHelpers::findBoundaryType(string& text_line, size_t boundary_type)
    {
        /*--- Read the elements of the problem (13 (id f1 f2 type 0) (---*/
        /** differnet boundary conditions
         * bc-type==2, interior boundary condition.
         * bc-type==3, wall boundary condition.
         * bc-type==4, Pressure Inlet boundary condition
         * bc-type==5, Pressure Outlet boundary condition
         * bc-type==7, Symmetry boundary condition
         * bc-type==9, pressure-far-field boundary condition.
         * bc-type==a, Velocity Inlet boundary condition.
         * bc-type==c, Periodic boundary condition.
         * bc-type==e, porous jumps boundary condition.
         * bc-type==14,Mass Flow Inlet boundary condition.
         * Note that Cell0 means boundary condition.
         * mesh_type==4, unstructured mesh.
         * mesh_type==6, structured mesh.
         * cell_lists_
         * {[(neighbor_cell_index, bc_type, node1_of_face, node2_of_face, node3_of_face), (....), (.....)], []..... }.
         * {inner_neighbor1, inner_neighbor2,inner_neighbor3, inner_neighbor4 ..... }.
         */
        /*--- find the type of boundary condition ---*/
        size_t position = text_line.find(")", 0);
        text_line = text_line.erase(0, position - 4);
        text_line = text_line.erase(2);
        boundary_type = stoi(text_line, nullptr, 16);
        return boundary_type;
    }

    Vecd MeshFileHelpers::nodeIndex(string& text_line)
    
    {
            /*--- find the node1 between two cells ---*/
        string node1_string_copy = text_line;
        size_t first_devide_position = text_line.find_first_of(" ", 0);
        string node1_index_string = node1_string_copy.erase(first_devide_position);
        size_t node1_index_decimal = stoi(node1_index_string, nullptr, 16) - 1;

        /*--- find the node2 between two cells---*/
        string node2_string = text_line;
        string node3_string = text_line;
        node2_string = node2_string.erase(0, first_devide_position + 1);
        node3_string = node2_string;
        size_t second_devide_position = node2_string.find_first_of(" ", 0);
        node2_string.erase(second_devide_position);
        size_t node2_index_decimal = stoi(node2_string, nullptr, 16) - 1;

        /*--- find the node3 between two cells---*/

        node3_string = node3_string.erase(0, second_devide_position + 1);
        size_t third_devide_position = node3_string.find_first_of(" ", 0);
        node3_string.erase(third_devide_position);
        size_t node3_index_decimal = stoi(node3_string, nullptr, 16) - 1;

        Vecd nodes = Vecd(node1_index_decimal, node2_index_decimal, node3_index_decimal);
        return nodes;
    }

    Vec2d MeshFileHelpers::cellIndex(string& text_line)
    {

        /*--- find the cell1---*/
        string cell1_string = text_line;
        string cell2_string = text_line;
        string cell3_string = text_line;
        size_t first_devide_position = text_line.find_first_of(" ", 0);
        cell1_string = cell1_string.erase(0, first_devide_position + 1);
        size_t second_devide_position = cell1_string.find_first_of(" ", 0);
        cell2_string = cell1_string.erase(0, second_devide_position + 1);
        size_t third_devide_position = cell2_string.find_first_of(" ", 0);
        cell3_string = cell2_string.erase(0, third_devide_position + 1);
        size_t fourth_devide_position = cell3_string.find_first_of(" ", 0);
        cell3_string.erase(fourth_devide_position);
        size_t cell1_index_decimal = stoi(cell3_string, nullptr, 16);

        /*--- find the cell2---*/
        string cell4_string = text_line;
        cell4_string = cell4_string.erase(0, first_devide_position + 1);
        cell4_string = cell4_string.erase(0, second_devide_position + 1);
        cell4_string = cell4_string.erase(0, third_devide_position + 1);
        cell4_string.erase(0, fourth_devide_position + 1);
        size_t cell2_index_decimal = stoi(cell4_string, nullptr, 16);
        Vec2d cells = Vec2d(cell1_index_decimal, cell2_index_decimal);        
        return cells;
    }

    void MeshFileHelpers::updateElementsNodesConnection(vector<vector<size_t>>& elements_nodes_connection_, Vecd nodes, Vec2d cells,
                                                        bool& check_neighbor_cell1, bool& check_neighbor_cell2)
    {

    /*--- build up connection with element and nodes only---*/
    for (int node = 0; node != nodes.size(); ++node)
    {
    if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != nodes[node] && elements_nodes_connection_[cells[check_neighbor_cell2]][1] != nodes[node] && elements_nodes_connection_[cells[check_neighbor_cell2]][2] != nodes[node] && elements_nodes_connection_[cells[check_neighbor_cell2]][3] != nodes[node])
    {   
        if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && (elements_nodes_connection_[cells[check_neighbor_cell2]][1] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1)) && elements_nodes_connection_[cells[check_neighbor_cell2]][2] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && elements_nodes_connection_[cells[check_neighbor_cell2]][3] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            elements_nodes_connection_[cells[check_neighbor_cell2]][3] = nodes[node];
        }
        if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && elements_nodes_connection_[cells[check_neighbor_cell2]][1] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && elements_nodes_connection_[cells[check_neighbor_cell2]][2] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
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
    }

    void MeshFileHelpers::updateCellLists(vector<vector<vector<size_t>>>& mesh_topology_, vector<vector<size_t>>& elements_nodes_connection_, Vecd nodes, 
                                            Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type)
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
                mesh_topology_[cells[check_neighbor_cell2]][0][4] = nodes[2];

                check_neighbor_cell1 = false;
                check_neighbor_cell2 = true;
                return;
            }
            if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
            {
                mesh_topology_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
                mesh_topology_[cells[check_neighbor_cell2]][1][1] = boundary_type;
                mesh_topology_[cells[check_neighbor_cell2]][1][2] = nodes[0];
                mesh_topology_[cells[check_neighbor_cell2]][1][3] = nodes[1];
                mesh_topology_[cells[check_neighbor_cell2]][1][4] = nodes[2];
                check_neighbor_cell1 = false;
                check_neighbor_cell2 = true;
                return;
            }
            if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
            {
                mesh_topology_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
                mesh_topology_[cells[check_neighbor_cell2]][2][1] = boundary_type;
                mesh_topology_[cells[check_neighbor_cell2]][2][2] = nodes[0];
                mesh_topology_[cells[check_neighbor_cell2]][2][3] = nodes[1];
                mesh_topology_[cells[check_neighbor_cell2]][2][4] = nodes[2];
                check_neighbor_cell1 = false;
                check_neighbor_cell2 = true;
                return;
            }
            if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][2][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
            {
                mesh_topology_[cells[check_neighbor_cell2]][3][0] = cells[check_neighbor_cell1];
                mesh_topology_[cells[check_neighbor_cell2]][3][1] = boundary_type;
                mesh_topology_[cells[check_neighbor_cell2]][3][2] = nodes[0];
                mesh_topology_[cells[check_neighbor_cell2]][3][3] = nodes[1];
                mesh_topology_[cells[check_neighbor_cell2]][3][4] = nodes[2];
                check_neighbor_cell1 = false;
                check_neighbor_cell2 = true;
                return;
            }
        }

    
    }

    void MeshFileHelpers::updateBoundaryCellLists(vector<vector<vector<size_t>>>& mesh_topology_, vector<vector<size_t>>& elements_nodes_connection_,
                                                    Vecd nodes, Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type)
    {
        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            mesh_topology_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
            mesh_topology_[cells[check_neighbor_cell2]][2][1] = boundary_type;
            mesh_topology_[cells[check_neighbor_cell2]][2][2] = nodes[0];
            mesh_topology_[cells[check_neighbor_cell2]][2][3] = nodes[1];
            mesh_topology_[cells[check_neighbor_cell2]][2][4] = nodes[2];
            check_neighbor_cell1 = false;
            check_neighbor_cell2 = true;
            return;
        }
        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && mesh_topology_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            mesh_topology_[cells[check_neighbor_cell2]][3][0] = cells[check_neighbor_cell1];
            mesh_topology_[cells[check_neighbor_cell2]][3][1] = boundary_type;
            mesh_topology_[cells[check_neighbor_cell2]][3][2] = nodes[0];
            mesh_topology_[cells[check_neighbor_cell2]][3][3] = nodes[1];
            mesh_topology_[cells[check_neighbor_cell2]][3][4] = nodes[2];
            check_neighbor_cell1 = false;
            check_neighbor_cell2 = true;
            return;
        } 
        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && mesh_topology_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            mesh_topology_[cells[check_neighbor_cell2]][3][0] = cells[check_neighbor_cell1];
            mesh_topology_[cells[check_neighbor_cell2]][3][1] = boundary_type;
            mesh_topology_[cells[check_neighbor_cell2]][3][2] = nodes[0];
            mesh_topology_[cells[check_neighbor_cell2]][3][3] = nodes[1];
            mesh_topology_[cells[check_neighbor_cell2]][3][4] = nodes[2];
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




    void MeshFileHelpers::cellCenterCoordinates(vector<vector<size_t>>& elements_nodes_connection_, std::vector<std::vector<long unsigned int>>::size_type& element,
                                                vector<vector<Real>>& node_coordinates_, StdLargeVec<Vecd>& elements_centroids_, Vecd& center_coordinate)
    {
        
        for (std::vector<long unsigned int>::size_type node = 0; node != elements_nodes_connection_[element].size(); ++node)
        {
            center_coordinate += Vecd(node_coordinates_[elements_nodes_connection_[element][node]][0] / 4.0,
                node_coordinates_[elements_nodes_connection_[element][node]][1] / 4.0, node_coordinates_[elements_nodes_connection_[element][node]][2] / 4.0);
        }
        elements_centroids_[element] = center_coordinate;
    } 

    void MeshFileHelpers::elementVolume(vector<vector<size_t>>& elements_nodes_connection_, std::vector<std::vector<long unsigned int>>::size_type& element,
                                        vector<vector<Real>>& node_coordinates_, StdLargeVec<Real>& elements_volumes_)
    {
        // get nodes position
        Real Total_volume = 0;
        using Vec4d = Eigen::Matrix<Real, 4, 1>;
        Vec4d nodes = Vec4d(elements_nodes_connection_[element][0], elements_nodes_connection_[element][1], elements_nodes_connection_[element][2], elements_nodes_connection_[element][3]);
        Vecd node1_coordinate = Vecd(node_coordinates_[nodes[0]][0], node_coordinates_[nodes[0]][1], node_coordinates_[nodes[0]][2]);
        Vecd node2_coordinate = Vecd(node_coordinates_[nodes[1]][0], node_coordinates_[nodes[1]][1], node_coordinates_[nodes[1]][2]);
        Vecd node3_coordinate = Vecd(node_coordinates_[nodes[2]][0], node_coordinates_[nodes[2]][1], node_coordinates_[nodes[2]][2]);
        Vecd node4_coordinate = Vecd(node_coordinates_[nodes[3]][0], node_coordinates_[nodes[3]][1], node_coordinates_[nodes[3]][2]);
        Matd M;
        M << node2_coordinate - node1_coordinate, node3_coordinate - node1_coordinate, node4_coordinate - node1_coordinate;
        Real determinant = abs(M.determinant());
        Real element_volume = (static_cast<double>(1) / 6) * determinant;
        Total_volume += element_volume;
        elements_volumes_[element] = element_volume;
    }

    void MeshFileHelpers::minimumdistance(vector<Real>& all_data_of_distance_between_nodes, StdLargeVec<Real>& elements_volumes_, vector<vector<vector<size_t>>>& mesh_topology_,
                                            vector<vector<Real>>& node_coordinates_)
    {

        for (size_t element_index = 0; element_index != elements_volumes_.size(); ++element_index)
        {
           
        for (std::vector<std::vector<long unsigned int>>::size_type neighbor = 0; neighbor != mesh_topology_[element_index].size(); ++neighbor)
        {
           
            size_t interface_node1_index = mesh_topology_[element_index][neighbor][2];
            size_t interface_node2_index = mesh_topology_[element_index][neighbor][3];
            size_t interface_node3_index = mesh_topology_[element_index][neighbor][4];
            Vecd node1_position = Vecd(node_coordinates_[interface_node1_index][0], node_coordinates_[interface_node1_index][1], node_coordinates_[interface_node1_index][2]);
            Vecd node2_position = Vecd(node_coordinates_[interface_node2_index][0], node_coordinates_[interface_node2_index][1], node_coordinates_[interface_node2_index][2]);
            Vecd node3_position = Vecd(node_coordinates_[interface_node3_index][0], node_coordinates_[interface_node3_index][1], node_coordinates_[interface_node3_index][2]);
            Vecd interface_area_vector1 = node2_position - node1_position;
            Vecd interface_area_vector2 = node3_position - node1_position;
            Vecd area_vector = interface_area_vector1.cross(interface_area_vector2);
            Real triangle_area = 0.5 * area_vector.norm();
            Real distance = sqrt(triangle_area);
            all_data_of_distance_between_nodes.push_back(distance);
        }
    }
    }

    void MeshFileHelpers::vtuFileHeader(std::ofstream& out_file)
    {
        out_file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
        out_file << "<UnstructuredGrid>\n";
        out_file << "<FieldData>\n";
        out_file << "<DataArray type=\"Int32\" Name=\"ispatch\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n";
        out_file << "0\n";
        out_file << "</DataArray>\n";
        out_file << "</FieldData>\n";
    }


    void MeshFileHelpers::vtuFileNodeCoordinates(std::ofstream& out_file, vector<vector<Real>>& node_coordinates_, vector<vector<size_t>>& elements_nodes_connection_)
    {
        // Write point data
        out_file << "<Piece NumberOfPoints=\"" << node_coordinates_.size() << "\" NumberOfCells=\"" << elements_nodes_connection_.size() << "\">\n";
        out_file << "<PointData>\n";
        out_file << "</PointData>\n";
        out_file << "<CellData>\n";
        out_file << "</CellData>\n";
        out_file << "<Points>\n";
        out_file << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"4.124318125496386\">\n";

        size_t total_nodes = node_coordinates_.size();
        for (size_t node = 0; node != total_nodes; ++node)
        {
            out_file << node_coordinates_[node][0] << " " << node_coordinates_[node][1] << " " << node_coordinates_[node][2] << " \n";
        }
    
    }

    void MeshFileHelpers::vtuFileInformationKey(std::ofstream& out_file)
    {
        out_file << "<InformationKey name=\"L2_NORM_RANGE\" location=\"vtkDataArray\" length=\"2\">\n";
        out_file << "<Value index=\"0\">\n";
        out_file << "0\n";
        out_file << "</Value>\n";
        out_file << "<Value index=\"1\">\n";
        out_file << "4.1243181255\n";
        out_file << "</Value>\n";
        out_file << "</InformationKey>\n";
        out_file << "<InformationKey name=\"L2_NORM_FINITE_RANGE\" location=\"vtkDataArray\" length=\"2\">\n";
        out_file << "<Value index=\"0\">\n";
        out_file << "0\n";
        out_file << "</Value>\n";
        out_file << "<Value index=\"1\">\n";
        out_file << "4.1243181255\n";
        out_file << "</Value>\n";
        out_file << "</InformationKey>\n";
        out_file << "</DataArray>\n";
        out_file << "</Points>\n";
      
    }

    void MeshFileHelpers::vtuFileCellConnectivity(std::ofstream& out_file, vector<vector<size_t>>& elements_nodes_connection_)
    {
        out_file << "<Cells>\n";
        // Write Cell data
        out_file << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"11138\">\n";

        for (const auto& cell : elements_nodes_connection_)
        {
            for (const auto& vertex : cell)
            {
                out_file << vertex << " ";
            }
            out_file << "\n";
        }

        out_file << "</DataArray>\n";
    }


    void MeshFileHelpers::vtuFileOffsets(std::ofstream& out_file, vector<vector<size_t>>& elements_nodes_connection_)
    {
        out_file << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"4\" RangeMax=\"212604\">\n";

        size_t offset = 0;
        for (const auto& face : elements_nodes_connection_)
        {
            offset += face.size();
            out_file << offset << " ";
        }
        out_file << "\n</DataArray>\n";
    }

    void MeshFileHelpers::vtuFileTypeOfCell(std::ofstream& out_file, vector<vector<size_t>>& elements_nodes_connection_)
    {
        size_t type = 10;
        //Specifies type of cell 10 = tetrahedral
        out_file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"10\" RangeMax=\"10\">\n";
        for (const auto& types : elements_nodes_connection_)
        {
            for (const auto& vertex : types)
            {
                out_file << type << " ";
            }
            out_file << "\n";
        }
        // Write face attribute data
        out_file << "</DataArray>\n";
        out_file << "</Cells>\n";
    }

    void MeshFileHelpers::numberofNodesFluent(ifstream& mesh_file, size_t& number_of_points, string& text_line)
    {

        while (getline(mesh_file, text_line))

        {
            text_line.erase(0, 1);
            if (atoi(text_line.c_str()) == 10 && text_line.find("))", 0) != string::npos)
            {
                string text1(text_line);
                text1.erase(0, 8);
                Real last_position = text1.find_last_of(")");
                text1.erase(last_position - 3);
                number_of_points = stoi(text1, nullptr, 16);
                break;
            }
        }
    }

    void MeshFileHelpers::numberofElementsFluent(ifstream& mesh_file, size_t& number_of_elements, string& text_line)
    {

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
                text1.erase(last_position - 3);
                number_of_elements = stoi(text1, nullptr, 16);
                break;
            }
        }
    }

    void MeshFileHelpers::nodeCoordinatesFluent(ifstream& mesh_file, size_t& node_index, vector<vector<Real>>& node_coordinates_, string& text_line,
                                                size_t& dimension, vector<vector<Real>>& point_coordinates)
    {
        while (getline(mesh_file, text_line))
        {
            text_line.erase(0, 1);
            if (atoi(text_line.c_str()) == 10 && text_line.find(") (", 0) != string::npos)
            {
                while (getline(mesh_file, text_line))
                {
                    if (text_line.find("(", 0) == string::npos && text_line.find("))", 0) == string::npos && text_line.find(" ", 0) != string::npos)
                    {
 
                        Real Coords[3] = { 0.0, 0.0, 0.0 };
                        if (dimension == 3)
                        {
                            size_t first_devide_position = text_line.find_first_of(" ");
                            size_t last_devide_position = text_line.find_last_of(" ");
                            string x_part = text_line;
                            string y_part = text_line;
                            string z_part = text_line;
                            string x_coordinate_string = x_part.erase(first_devide_position);
                            string y_coordinate_string = y_part.erase(last_devide_position);
                            y_coordinate_string = y_coordinate_string.erase(0, first_devide_position);
                            string z_coordinate_string = z_part.erase(0, last_devide_position);
                            istringstream streamx, streamy, streamz;
                            streamx.str(x_coordinate_string);
                            streamy.str(y_coordinate_string);
                            streamz.str(z_coordinate_string);
                            streamx >> Coords[0];
                            streamy >> Coords[1];
                            streamz >> Coords[2];
                            node_coordinates_[node_index][0] = Coords[0];
                            node_coordinates_[node_index][1] = Coords[1];
                            node_coordinates_[node_index][2] = Coords[2];
                            ++node_index;
                        }
                        for (std::size_t iDim = 0; iDim != dimension; ++iDim)
                        {
                            point_coordinates[iDim].push_back(Coords[iDim]);
                        }
                    }
                    text_line.erase(0, 1);
                    if (atoi(text_line.c_str()) == 11 || atoi(text_line.c_str()) == 13 || atoi(text_line.c_str()) == 12)
                    {
                        break;
                    }
                }

            }

            if ((atoi(text_line.c_str()) == 11 && text_line.find(") (") != string::npos) || (atoi(text_line.c_str()) == 13 && text_line.find(") (") != string::npos) || (atoi(text_line.c_str()) == 12 && text_line.find(") (") != string::npos))
                break;
        }
    }

    void MeshFileHelpers::updateBoundaryCellListsFluent(vector<vector<vector<size_t>>>& mesh_topology_, vector<vector<size_t>>& elements_nodes_connection_, 
                                                        Vecd nodes, Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type)
    {

        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            mesh_topology_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
            mesh_topology_[cells[check_neighbor_cell2]][1][1] = boundary_type;
            mesh_topology_[cells[check_neighbor_cell2]][1][2] = nodes[0];
            mesh_topology_[cells[check_neighbor_cell2]][1][3] = nodes[1];
            mesh_topology_[cells[check_neighbor_cell2]][1][4] = nodes[2];
            check_neighbor_cell1 = false;
            check_neighbor_cell2 = true;
            return;
        }
        if (mesh_topology_[cells[check_neighbor_cell2]][0][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && mesh_topology_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && mesh_topology_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && mesh_topology_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            mesh_topology_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
            mesh_topology_[cells[check_neighbor_cell2]][2][1] = boundary_type;
            mesh_topology_[cells[check_neighbor_cell2]][2][2] = nodes[0];
            mesh_topology_[cells[check_neighbor_cell2]][2][3] = nodes[1];
            mesh_topology_[cells[check_neighbor_cell2]][2][4] = nodes[2];
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
}