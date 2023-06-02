
#include "common_shared_FVM_classes.h"

namespace SPH
{
	//=================================================================================================//
	void readMeshFile::getDataFromMeshFile()
	{
		ifstream mesh_file;  /*!< \brief File object for the Ansys ASCII mesh file. */
		mesh_file.open(full_path_);
		if (mesh_file.fail())
		{
			cout << "Error:Check that the file exists." << endl;
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
		/*--- Read the node data (index is starting from zero) ---*/
		int number_of_points(0);
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
		point_coordinates_2D_.resize(number_of_points);
		for (int point = 0; point != point_coordinates_2D_.size(); ++point)
		{
			point_coordinates_2D_[point].resize(dimension);
		}

		/* Prepare our data structure for the point coordinates. */
		point_coordinates.resize(dimension);
		for (int k = 0; k < dimension; k++)
		{
			point_coordinates[k].reserve(number_of_points);
		}
		int node_index = 0;
		while (getline(mesh_file, text_line))
		{
			if (text_line.find("(", 0) == string::npos && text_line.find("))", 0) == string::npos)
			{
				double Coords[3] = { 0.0, 0.0, 0.0 };
				if (text_line.find(" ", 0) != string::npos)
				{
					if (dimension == 2)
					{
						int devide_position = text_line.find_first_of(" ");
						string x_part = text_line;
						string y_part = text_line;
						string x_coordinate_string = x_part.erase(devide_position);
						string y_coordinate_string = y_part.erase(0, devide_position);

						istringstream streamx, streamy;
						streamx.str(x_coordinate_string);
						streamy.str(y_coordinate_string);
						streamx >> Coords[0];
						streamy >> Coords[1];
						point_coordinates_2D_[node_index][0] = Coords[0];
						point_coordinates_2D_[node_index][1] = Coords[1];
						++node_index;
					}
					if (dimension == 3)
					{
						int first_devide_position = text_line.find_first_of(" ");
						int last_devide_position = text_line.find_last_of(" ");
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
						point_coordinates_2D_[node_index][0] = Coords[0];
						point_coordinates_2D_[node_index][1] = Coords[1];
						point_coordinates_2D_[node_index][2] = Coords[2];
						++node_index;
					}
					for (int iDim = 0; iDim != dimension; ++iDim)
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

		/*--- Read the elements of the problem ---*/
		/** differnet boundary conditions
		* bc-type==2, interior boundary condition.
		* bc-type==3, wall boundary condition.
		* bc-type==9, pressure-far-field boundary condition.
		* Note that Cell0 means boundary condition.
		* mesh_type==3, unstructured mesh.
		* mesh_type==4, structured mesh.
		* cell_lists_
		* {[(neighbor_cell_index, bc_type, node1_of_face, node2_of_face), (....), (.....)], []..... }.
		* {inner_neighbor1, inner_neighbor2, ..... }.
		*/
		int boundary_type(0);
		int number_of_elements(0);
		int mesh_type = 3;
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
		cell_lists_.resize(number_of_elements + 1);
		for (int a = 0; a != number_of_elements + 1; ++a)
		{
			cell_lists_[a].resize(mesh_type);
			for (int b = 0; b != cell_lists_[a].size(); ++b)
			{
				cell_lists_[a][b].resize(dimension + 2);
				for (int c = 0; c != cell_lists_[a][b].size(); ++c)
				{
					cell_lists_[a][b][c] = -1;
				}

			}
		}
		/*--- reinitialize the number of elements ---*/
		elements_nodes_connection_.resize(number_of_elements + 1);
		elements_neighbors_connection_.resize(number_of_elements + 1);
		cell_lists_.resize(number_of_elements + 1);
		for (int element = 0; element != number_of_elements + 1; ++element)
		{
			elements_nodes_connection_[element].resize(3);
			for (int node = 0; node != elements_nodes_connection_[element].size(); ++node)
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
				int position = text_line.find(")", 0);
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
						int first_devide_position = text_line.find_first_of(" ", 0);
						string node1_index_string = node1_string_copy.erase(first_devide_position);
						int node1_index_decimal = stoi(node1_index_string, nullptr, 16) - 1;

						/*--- find the node2 between two cells---*/
						string node2_string = text_line;
						node2_string = node2_string.erase(0, first_devide_position + 1);
						int second_devide_position = node2_string.find_first_of(" ", 0);
						node2_string.erase(second_devide_position);
						int node2_index_decimal = stoi(node2_string, nullptr, 16) - 1;
						Vecd nodes = Vecd(node1_index_decimal, node2_index_decimal);

						/*--- find the cell1---*/
						string cell1_string = text_line;
						cell1_string = cell1_string.erase(0, first_devide_position + 1);
						cell1_string = cell1_string.erase(0, second_devide_position + 1);
						int third_devide_position = cell1_string.find_first_of(" ", 0);
						cell1_string.erase(third_devide_position);
						int cell1_index_decimal = stoi(cell1_string, nullptr, 16);

						/*--- find the cell2---*/
						string cell2_string = text_line;
						cell2_string = cell2_string.erase(0, first_devide_position + 1);
						cell2_string = cell2_string.erase(0, second_devide_position + 1);
						cell2_string.erase(0, third_devide_position + 1);
						int cell2_index_decimal = stoi(cell2_string, nullptr, 16);
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
										if (cells[check_neighbor_cell2] == 9431)
										{
											Real a = 1.0;
										}
										int num = cells[check_neighbor_cell2];
										/*--- build up connection with element and nodes only---*/
										for (int node = 0; node != nodes.size(); ++node)
										{
											if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != nodes[node]
												&& elements_nodes_connection_[cells[check_neighbor_cell2]][1] != nodes[node]
												&& elements_nodes_connection_[cells[check_neighbor_cell2]][2] != nodes[node])
											{
												if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != -1
													&& elements_nodes_connection_[cells[check_neighbor_cell2]][1] != -1
													&& elements_nodes_connection_[cells[check_neighbor_cell2]][2] == -1)
												{
													elements_nodes_connection_[cells[check_neighbor_cell2]][2] = nodes[node];
												}
												if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != -1
													&& elements_nodes_connection_[cells[check_neighbor_cell2]][1] == -1)
												{
													elements_nodes_connection_[cells[check_neighbor_cell2]][1] = nodes[node];
												}
												if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] == -1)
												{
													elements_nodes_connection_[cells[check_neighbor_cell2]][0] = nodes[node];
												}
											}
											else continue;
										}
										/*--- build up all connection data with element and neighbor and nodes---*/
										if (cell_lists_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1]
											&& cell_lists_[cells[check_neighbor_cell2]][1][0] != cells[check_neighbor_cell1]
											&& cell_lists_[cells[check_neighbor_cell2]][2][0] != cells[check_neighbor_cell1])
										{
											if (cell_lists_[cells[check_neighbor_cell2]][0][0] == -1)
											{
												/*--- inner neighbor index---*/
												cell_lists_[cells[check_neighbor_cell2]][0][0] = cells[check_neighbor_cell1];
												/*--- boundary type---*/
												cell_lists_[cells[check_neighbor_cell2]][0][1] = boundary_type;
												/*--- nodes of a face---*/
												cell_lists_[cells[check_neighbor_cell2]][0][2] = nodes[0];
												cell_lists_[cells[check_neighbor_cell2]][0][3] = nodes[1];

												check_neighbor_cell1 = false;
												check_neighbor_cell2 = true;
												break;
											}
											if (cell_lists_[cells[check_neighbor_cell2]][0][0] != -1 && cell_lists_[cells[check_neighbor_cell2]][1][0] == -1)
											{
												cell_lists_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
												cell_lists_[cells[check_neighbor_cell2]][1][1] = boundary_type;
												cell_lists_[cells[check_neighbor_cell2]][1][2] = nodes[0];
												cell_lists_[cells[check_neighbor_cell2]][1][3] = nodes[1];
												check_neighbor_cell1 = false;
												check_neighbor_cell2 = true;
												break;
											}
											if (cell_lists_[cells[check_neighbor_cell2]][0][0] != -1
												&& cell_lists_[cells[check_neighbor_cell2]][1][0] != -1
												&& cell_lists_[cells[check_neighbor_cell2]][2][0] == -1)
											{
												cell_lists_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
												cell_lists_[cells[check_neighbor_cell2]][2][1] = boundary_type;
												cell_lists_[cells[check_neighbor_cell2]][2][2] = nodes[0];
												cell_lists_[cells[check_neighbor_cell2]][2][3] = nodes[1];
												check_neighbor_cell1 = false;
												check_neighbor_cell2 = true;
												break;
											}
										}
										if (cell_lists_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1]
											&& cell_lists_[cells[check_neighbor_cell2]][1][0] == 0
											&& cell_lists_[cells[check_neighbor_cell2]][2][0] == -1)
										{
											cell_lists_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
											cell_lists_[cells[check_neighbor_cell2]][2][1] = boundary_type;
											cell_lists_[cells[check_neighbor_cell2]][2][2] = nodes[0];
											cell_lists_[cells[check_neighbor_cell2]][2][3] = nodes[1];
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
										int a = cells[check_neighbor_cell2];
										if (cell_lists_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1]
											&& cell_lists_[cells[check_neighbor_cell2]][1][0] != cells[check_neighbor_cell1]
											&& cell_lists_[cells[check_neighbor_cell2]][2][0] != cells[check_neighbor_cell1]
											&& cell_lists_[cells[check_neighbor_cell2]][3][0] != cells[check_neighbor_cell1])
										{
											if (cell_lists_[cells[check_neighbor_cell2]][0][0] == -1)
											{
												/*--- inner neighbor index---*/
												cell_lists_[cells[check_neighbor_cell2]][0][0] = cells[check_neighbor_cell1];
												/*--- boundary type---*/
												cell_lists_[cells[check_neighbor_cell2]][0][1] = boundary_type;
												/*--- nodes of a face---*/
												cell_lists_[cells[check_neighbor_cell2]][0][2] = nodes[0];
												cell_lists_[cells[check_neighbor_cell2]][0][3] = nodes[1];

												check_neighbor_cell1 = false;
												check_neighbor_cell2 = true;
												break;
											}
											if (cell_lists_[cells[check_neighbor_cell2]][0][0] != -1 && cell_lists_[cells[check_neighbor_cell2]][1][0] == -1)
											{
												cell_lists_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
												cell_lists_[cells[check_neighbor_cell2]][1][1] = boundary_type;
												cell_lists_[cells[check_neighbor_cell2]][1][2] = nodes[0];
												cell_lists_[cells[check_neighbor_cell2]][1][3] = nodes[1];
												check_neighbor_cell1 = false;
												check_neighbor_cell2 = true;
												break;
											}
											if (cell_lists_[cells[check_neighbor_cell2]][0][0] != -1
												&& cell_lists_[cells[check_neighbor_cell2]][1][0] != -1
												&& cell_lists_[cells[check_neighbor_cell2]][2][0] == -1)
											{
												cell_lists_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
												cell_lists_[cells[check_neighbor_cell2]][2][1] = boundary_type;
												cell_lists_[cells[check_neighbor_cell2]][2][2] = nodes[0];
												cell_lists_[cells[check_neighbor_cell2]][2][3] = nodes[1];
												check_neighbor_cell1 = false;
												check_neighbor_cell2 = true;
												break;
											}
											if (cell_lists_[cells[check_neighbor_cell2]][0][0] != -1
												&& cell_lists_[cells[check_neighbor_cell2]][1][0] != -1
												&& cell_lists_[cells[check_neighbor_cell2]][2][0] != -1
												&& cell_lists_[cells[check_neighbor_cell2]][3][0] == -1)
											{
												cell_lists_[cells[check_neighbor_cell2]][3][0] = cells[check_neighbor_cell1];
												cell_lists_[cells[check_neighbor_cell2]][3][1] = boundary_type;
												cell_lists_[cells[check_neighbor_cell2]][3][2] = nodes[0];
												cell_lists_[cells[check_neighbor_cell2]][3][3] = nodes[1];
												check_neighbor_cell1 = false;
												check_neighbor_cell2 = true;
												break;
											}
											//if it is three dimension
										}
										else
										{
											check_neighbor_cell1 = false;
											check_neighbor_cell2 = true;
											break;
										}
									}
									if (cells[check_neighbor_cell2] == 0)break;
								}
							}
						}
					}
					else break;
				}
			}
			if (text_line.find("Zone Sections", 0) != string::npos) break;
			if (text_line.find(")", 0) != string::npos != string::npos) continue;
		}
	}
    //=================================================================================================//
    void readMeshFile::getElementCenterCoordinates()
    {
        elements_center_coordinates_.resize(elements_nodes_connection_.size());
        elements_volumes_.resize(elements_nodes_connection_.size());
        for (int element = 1; element != elements_nodes_connection_.size(); ++element)
        {
            Vecd center_coordinate = Vecd::Zero();
            for (int node = 0; node != elements_nodes_connection_[element].size(); ++node)
            {
                int node_index = elements_nodes_connection_[element][node];
                center_coordinate += Vecd(point_coordinates_2D_[elements_nodes_connection_[element][node]][0] / 3.0,
                                            point_coordinates_2D_[elements_nodes_connection_[element][node]][1] / 3.0);
            }
            elements_center_coordinates_[element] = center_coordinate;

            // calculating each volume of element
            //get nodes position
            Vec3d nodes = Vec3d(elements_nodes_connection_[element][0], elements_nodes_connection_[element][1], elements_nodes_connection_[element][2]);
            Vecd node1_coordinate = Vecd(point_coordinates_2D_[nodes[0]][0], point_coordinates_2D_[nodes[0]][1]);
            Vecd node2_coordinate = Vecd(point_coordinates_2D_[nodes[1]][0], point_coordinates_2D_[nodes[1]][1]);
            Vecd node3_coordinate = Vecd(point_coordinates_2D_[nodes[2]][0], point_coordinates_2D_[nodes[2]][1]);
            //get each line length
            Real first_side_length = (node1_coordinate - node2_coordinate).norm();
            Real second_side_length = (node1_coordinate - node3_coordinate).norm();
            Real third_side_length = (node2_coordinate - node3_coordinate).norm();
            //half perimeter
            Real half_perimeter = (first_side_length + second_side_length + third_side_length) / 2.0;
            //get element volume
            Real element_volume =
                pow(half_perimeter * (half_perimeter - first_side_length) * (half_perimeter - second_side_length) * (half_perimeter - third_side_length), 0.5);
            elements_volumes_[element] = element_volume;
        }
    }
    //=================================================================================================//
    void NeighborhoodInFVM::removeANeighbor(size_t neighbor_n) 
	{
        current_size_--;
        j_[neighbor_n] = j_[current_size_];
        dW_ijV_j_[neighbor_n] = dW_ijV_j_[current_size_];
        r_ij_[neighbor_n] = r_ij_[current_size_];
        interface_size_[neighbor_n] = interface_size_[current_size_];
        boundary_type_[neighbor_n] = boundary_type_[current_size_];
	}
	//=================================================================================================//
	void BaseInnerRelationInFVM::resetNeighborhoodCurrentSize()
	{
		parallel_for
		(
            IndexRange(0, base_particles_.total_real_particles_),
            [&](const IndexRange &r)
            {
                for (size_t num = r.begin(); num != r.end(); ++num)
                {
                    inner_configuration_in_FVM_[num].current_size_ = 0;
                }
            },
        ap);
	}
	//=================================================================================================//
	 BaseInnerRelationInFVM::BaseInnerRelationInFVM(RealBody &real_body, vector<vector<vector<int>>> data_inpute, vector<vector<Real>> nodes_coordinates)
        : SPHRelation(real_body), real_body_(&real_body)
    {
        all_needed_data_from_mesh_file_ = data_inpute;
        nodes_coordinates_ = nodes_coordinates;
        subscribeToBody();
        resizeConfiguration();
    };
	//=================================================================================================//
    void BaseInnerRelationInFVM::resizeConfiguration() 
	{
        size_t updated_size = base_particles_.real_particles_bound_;
        inner_configuration_in_FVM_.resize(updated_size, NeighborhoodInFVM());
	}
	//=================================================================================================//
    ParticleGeneratorInFVM::ParticleGeneratorInFVM(SPHBody &sph_body, const StdLargeVec<Vecd> &positions, const StdLargeVec<Real> &elements_volumes)
        : ParticleGenerator(sph_body), elements_center_coordinates_(positions), elements_volumes_(elements_volumes) {}
	//=================================================================================================//
	void ParticleGeneratorInFVM::initializeGeometricVariables()
	{
        for (size_t particle_index = 1; particle_index != elements_center_coordinates_.size(); ++particle_index)
        {
            initializePositionAndVolumetricMeasure(elements_center_coordinates_[particle_index], elements_volumes_[particle_index]);
        }
	}
	//=================================================================================================//
	void NeighborBuilderInFVM::createRelation(NeighborhoodInFVM &neighborhood, Real &distance, Real &dW_ijV_j, Real &interface_size,
                            Vecd &interface_normal_direction, size_t bc_type, size_t j_index) const
	{
        neighborhood.j_.push_back(j_index);
        neighborhood.r_ij_.push_back(distance);
        neighborhood.interface_size_.push_back(interface_size);
        neighborhood.e_ij_.push_back(interface_normal_direction);
        neighborhood.dW_ijV_j_.push_back(dW_ijV_j);
        neighborhood.boundary_type_.push_back(bc_type);
        neighborhood.allocated_size_++;
	}
	//=================================================================================================//
	void NeighborBuilderInFVM::initializeRelation(NeighborhoodInFVM &neighborhood, Real &distance, Real &dW_ijV_j, Real &interface_size,
                                Vecd &interface_normal_direction, size_t bc_type, size_t j_index)const
	{
		size_t current_size = neighborhood.current_size_;
		neighborhood.j_[current_size] = j_index;
		neighborhood.dW_ijV_j_[current_size] = dW_ijV_j;
		neighborhood.r_ij_[current_size] = distance;
		neighborhood.interface_size_[current_size] = interface_size;
		neighborhood.e_ij_[current_size] = interface_normal_direction;
		neighborhood.boundary_type_[current_size] = bc_type;
	}
	//=================================================================================================//
	InnerRelationInFVM::InnerRelationInFVM(RealBody &real_body, vector<vector<vector<int>>> data_inpute, vector<vector<Real>> nodes_coordinates)
            : BaseInnerRelationInFVM(real_body, data_inpute, nodes_coordinates), get_inner_neighbor_(&real_body){};
	//=================================================================================================//
	template <typename GetParticleIndex, typename GetNeighborRelation>
	void InnerRelationInFVM::searchNeighborsByParticles(size_t total_real_particles, BaseParticles &source_particles,
        ParticleConfigurationInFVM &particle_configuration, GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation)
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

                    NeighborhoodInFVM &neighborhood = particle_configuration[index_i];
                    for (int neighbor = 0; neighbor != all_needed_data_from_mesh_file_[index_i + 1].size(); ++neighbor)
                    {
                        size_t index_j = all_needed_data_from_mesh_file_[index_i + 1][neighbor][0] - 1;
                        size_t boundary_type = all_needed_data_from_mesh_file_[index_i + 1][neighbor][1];
                        size_t interface_node1_index = all_needed_data_from_mesh_file_[index_i + 1][neighbor][2];
                        size_t interface_node2_index = all_needed_data_from_mesh_file_[index_i + 1][neighbor][3];
                        Vecd node1_position = Vecd(nodes_coordinates_[interface_node1_index][0], nodes_coordinates_[interface_node1_index][1]);
                        Vecd node2_position = Vecd(nodes_coordinates_[interface_node2_index][0], nodes_coordinates_[interface_node2_index][1]);
                        Vecd interface_area_vector = node1_position - node2_position;
                        Real interface_area_size = interface_area_vector.norm();
                        Vecd unit_vector = interface_area_vector / interface_area_size;
                        //normal unit vector
                        Vecd normal_vector = Vecd(unit_vector[1], -unit_vector[0]);
                        //judge the direction
                        Vecd node1_to_center_direction = particle_position - node1_position;
                        if (node1_to_center_direction.dot(normal_vector) < 0)
                        {
                            normal_vector = -normal_vector;
                        };
                        Real r_ij = 0; // we need r_ij to calculate the viscous force
                        //boundary_type == 2 means both of them are inside of fluid
                        if (boundary_type == 2)
                        {
                            r_ij = (particle_position - pos_n[index_j]).dot(normal_vector);
                        }
                        //boundary_type == 3 means fulid particle with wall boundary
                        if (boundary_type == 3 | boundary_type == 4 | boundary_type == 9 | boundary_type == 10 | boundary_type == 36)
                        {
                            r_ij = node1_to_center_direction.dot(normal_vector) * 2.0;
                        }
                        Real dW_ijV_j = -interface_area_size / (2.0 * Vol_i);
                        get_neighbor_relation(neighborhood, r_ij, dW_ijV_j, interface_area_size, normal_vector, boundary_type, index_j);
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
									base_particles_, inner_configuration_in_FVM_,
									get_particle_index_, get_inner_neighbor_);
	}
	//=================================================================================================//
}
//=================================================================================================//