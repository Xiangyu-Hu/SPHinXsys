#include "implicit_surface.h"
#include "base_data_type.h"

#include <sstream>
#include <iterator>

namespace SPH
{

    ImplicitSurface::~ImplicitSurface()
    {

    }

    void ImplicitSurface::initialize(BoundingBox bounding_box, Real grid_spacing)
	{
		grid_spacing_ = grid_spacing;
		bounding_box_ = bounding_box;

		mesh_lower_bound_ = bounding_box.first;
		mesh_upper_bound_ = bounding_box.second;
		for(auto i=0; i!=Dimensions; ++i)
			number_of_grid_points_[i] = 1 + static_cast<int>(ceil((mesh_upper_bound_[i] - mesh_lower_bound_[i]) / grid_spacing));

		number_of_cells_ = number_of_grid_points_ - Vecu(1);

		data_.resize(number_of_cells_[0] * number_of_cells_[1] * number_of_cells_[2]);

		dataOrigin_ = mesh_lower_bound_ + Vecd(0.5 * grid_spacing);
	}

	Vec3d ImplicitSurface::findClosestPoint(const Vec3d& input_pnt)
	{
		auto sdf = findSignedDistance(input_pnt);
		auto normal = findNormalDirection(input_pnt);
		Vec3d pt = input_pnt - sdf * normal;

		for(auto i=0; i<maxNumberOfIterations_; ++i)
		{
			sdf = findSignedDistance(pt);
			if(std::fabs(sdf) < epsilon_)
				return pt;
			else
			{
				pt = pt - sdf * normal;
			}
		}

		return pt;
	}
    
    BoundingBox ImplicitSurface::findBounds() 
    {
      	return bounding_box_;
    }

	Vecd ImplicitSurface::getDataOrigin() const
	{
		return dataOrigin_;
	}

	size_t ImplicitSurface::convertToOneDIndex(Vec3u index)
	{
		auto oneDIndex = index[0] + index[1] * number_of_cells_[0] + index[2] * number_of_cells_[0] * number_of_cells_[1];
		return oneDIndex;
	}

	Real ImplicitSurface::getCellCenterData(Vec3u index)
	{
		auto oneDIndex = convertToOneDIndex(index);
		return data_[oneDIndex];
	}

	Vec3d ImplicitSurface::getCellCenterGrad(Vec3u index)
	{
		auto oneDIndex = convertToOneDIndex(index);
		return Vec3d(gradDataX_[oneDIndex], gradDataY_[oneDIndex], gradDataZ_[oneDIndex]);
	}

	void ImplicitSurface::computeGradient()
	{
		gradDataX_.resize(number_of_cells_[0] * number_of_cells_[1] * number_of_cells_[2]);
		gradDataY_.resize(number_of_cells_[0] * number_of_cells_[1] * number_of_cells_[2]);
		gradDataZ_.resize(number_of_cells_[0] * number_of_cells_[1] * number_of_cells_[2]);

		for(size_t k=0; k< number_of_cells_[2]; ++k)
			for(size_t j=0; j<number_of_cells_[1]; ++j)
				for(size_t i=0; i<number_of_cells_[0]; ++i)
				{
					size_t ip1, im1, jp1, jm1, kp1, km1;
					ip1 = (i==number_of_cells_[0]-1) ? i : i + 1;
					im1 = (i==0) ? i : i - 1;
					jp1 = (j==number_of_cells_[1]-1) ? j : j + 1;
					jm1 = (j==0) ? j : j - 1;
					kp1 = (k==number_of_cells_[2]-1) ? k : k + 1;
					km1 = (k==0) ? k : k - 1;

					Vec3u left(im1, j, k);
					Vec3u right(ip1, j, k);
					Vec3u top(i, jp1, k);
					Vec3u down(i, jm1, k);
					Vec3u back(i, j, km1);
					Vec3u front(i, j, kp1);

					Real dx = 0.5 * (data_[convertToOneDIndex(right)] - data_[convertToOneDIndex(left)]) / grid_spacing_;
					Real dy = 0.5 * (data_[convertToOneDIndex(top)] - data_[convertToOneDIndex(down)]) / grid_spacing_;
					Real dz = 0.5 * (data_[convertToOneDIndex(front)] - data_[convertToOneDIndex(back)]) / grid_spacing_;
					auto index = convertToOneDIndex(Vec3u(i,j,k));
					gradDataX_[index] = dx;
					gradDataY_[index] = dy; 
					gradDataZ_[index] = dz;
			}
	}

	std::pair<std::array<Vec3u, 8>, std::array<Real, 8>> ImplicitSurface::computeIndexAndWeights(const Vec3d &point)
	{
		Vec3d index(0,0,0);
		index[0] = (point[0] - dataOrigin_[0]) / grid_spacing_;
		index[1] = (point[1] - dataOrigin_[1]) / grid_spacing_;
		index[2] = (point[2] - dataOrigin_[2]) / grid_spacing_;

		index[0] = clamp<Real>(index[0], 0., number_of_cells_[0]);
		index[1] = clamp<Real>(index[1], 0., number_of_cells_[1]);
		index[2] = clamp<Real>(index[2], 0., number_of_cells_[2]);

		auto lowerX = std::floor(index[0]);
		auto lowerY = std::floor(index[1]);
		auto lowerZ = std::floor(index[2]);

		auto weightX = index[0] - lowerX;
		auto weightY = index[1] - lowerY;
		auto weightZ = index[2] - lowerZ;

		size_t i = static_cast<size_t>(lowerX);
		size_t j = static_cast<size_t>(lowerY);
		size_t k = static_cast<size_t>(lowerZ);
		size_t ip1 = SMIN(i+1, number_of_cells_[0]-1);
		size_t jp1 = SMIN(j+1, number_of_cells_[1]-1);
		size_t kp1 = SMIN(k+1, number_of_cells_[2]-1);

		Vec3u cell000(i, j, k);
		Vec3u cell100(ip1, j, k);
		Vec3u cell010(i, jp1, k);
		Vec3u cell110(ip1, jp1, k);
		Vec3u cell001(i, j, kp1);
		Vec3u cell101(ip1, j, kp1);
		Vec3u cell011(i, jp1, kp1);
		Vec3u cell111(ip1, jp1, kp1);

		Real weight000 = (1 - weightX) * (1 - weightY) * (1 - weightZ);
		Real weight100 = weightX * (1 - weightY) * (1 - weightZ);
		Real weight010 = (1 - weightX) * weightY * (1 - weightZ);
		Real weight110 = weightX * weightY * (1- weightZ);
		Real weight001 = (1 - weightX) * (1 - weightY) * weightZ;
		Real weight101 = weightX * (1 - weightY) * weightZ;
		Real weight011 = (1 - weightX) * weightY * weightZ;
		Real weight111 = weightX * weightY * weightZ;

		std::array<Vec3u, 8> cellIndexs{cell000, cell100, cell010, cell110, cell001, cell101, cell011, cell111};
		std::array<Real, 8> weights{weight000, weight100, weight010, weight110, weight001, weight101, weight011, weight111};

		return std::make_pair(cellIndexs, weights);
	}

	void ImplicitSurface::addAImplicitSphere(Vec3d center, Real radius, BoundingBox bounding_box, Real grid_spacing)
	{
		initialize(bounding_box, grid_spacing);
		for(auto k=0; k<number_of_cells_[2]; ++k)
			for(auto j=0; j<number_of_cells_[1]; ++j)
				for(auto i=0; i<number_of_cells_[0]; ++i)
				{
					auto oneDIndex = convertToOneDIndex(Vec3u(i,j,k));
					auto position = mesh_lower_bound_;
					position[0] += (i + 0.5) * grid_spacing_;
					position[1] += (j + 0.5) * grid_spacing_;
					position[2] += (k + 0.5) * grid_spacing_;
					auto distance = (position - center).norm() - radius;
					data_[oneDIndex] = distance;
				}

		computeGradient();
	}

	void ImplicitSurface::loadFromMHDImge(const std::string &filename)
	{
		// A helper function to split a string
		auto split = [](const std::string &info)
		{
			std::istringstream ss(info);
			std::vector<std::string> vs;
			for(auto i=std::istream_iterator<std::string>(ss); i!=std::istream_iterator<std::string>(); ++i)
			{
				vs.push_back(*i);
			}

			return vs;
		};
		// Data type of the image data
		enum class DataType
		{
			UCHAR,
			SHORT,
			FLOAT
		};
		// The default value of the datatype is float.
		DataType dataType = DataType::FLOAT;

		std::ifstream infoFile(filename);
		if(!infoFile.is_open())
		{
			std::cout << "\n FAILURE: open information data file " << filename << " failed!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		std::string lineInfo;
		std::string rawDataFilePath;

		while(std::getline(infoFile, lineInfo))
		{
			if(lineInfo.find("DimSize = ")!=std::string::npos)
			{
				size_t pos = lineInfo.find("=");
				auto dimString = lineInfo.substr(pos+1);
				auto values = split(dimString);
				for(auto i=0; i<3; ++i)
					number_of_cells_[i] = std::atoi(values[i].c_str());
			}
			else if(lineInfo.find("ElementSpacing = ")!=std::string::npos || 
			        lineInfo.find("ElementSize = ")!=std::string::npos)
			{
				size_t pos = lineInfo.find("=");
				auto spacingString = lineInfo.substr(pos+1);
				auto values = split(spacingString);
				std::vector<double> spacingVec;
				for(auto i=0; i<3; ++i)
					spacingVec.push_back(std::atof(values[i].c_str()));

				if(std::fabs(spacingVec[0]-spacingVec[1])>1e-10 ||
				   std::fabs(spacingVec[0]-spacingVec[2])>1e-10 ||
				   std::fabs(spacingVec[2]-spacingVec[1])>1e-10)
				{
					std::cout << "\n FAILURE: the spacing must be all the same at all directions!" << std::endl;
					std::cout << __FILE__ << ':' << __LINE__ << std::endl;
					exit(1);
				}
				else
					grid_spacing_ = spacingVec[0];
			}
			else if(lineInfo.find("Position = ")!=std::string::npos || 
			        lineInfo.find("Offset = ")!=std::string::npos)
			{
				size_t pos = lineInfo.find("=");
				auto offsetString = lineInfo.substr(pos+1);
				auto values = split(offsetString);
				for(auto i=0; i<3; ++i)
					dataOrigin_[i] = std::atof(values[i].c_str());
			}
			else if(lineInfo.find("ElementType = ")!=std::string::npos)
			{
				size_t pos = lineInfo.find("=");
				auto typeString = lineInfo.substr(pos+1);
				if(typeString.compare(5,5,"UCHAR") == 0)
					dataType = DataType::UCHAR;
				else if(typeString.compare(5,5,"SHORT") == 0)
					dataType = DataType::SHORT;
			}
			else if(lineInfo.find("ElementDataFile = ")!=std::string::npos)
			{
				size_t pos = lineInfo.find("=");
				auto dataFileString = lineInfo.substr(pos+2);
				std::string prefixPath = filename.substr(0, filename.find_last_of("\\/"));
				rawDataFilePath = prefixPath + "/" + dataFileString;
			}

		}// end while
		infoFile.close();
		// initialize
		// Note: Here we should compute the lower bound by data origin.
		mesh_lower_bound_ = dataOrigin_ - Vecd(0.5*grid_spacing_,
		      								   0.5*grid_spacing_, 
											   0.5*grid_spacing_);

		mesh_upper_bound_ = mesh_lower_bound_ + Vecd(number_of_cells_[0]*grid_spacing_,
		      										 number_of_cells_[1]*grid_spacing_, 
													 number_of_cells_[2]*grid_spacing_);

		bounding_box_.first = mesh_lower_bound_;
		bounding_box_.second = mesh_upper_bound_;

		number_of_grid_points_ = number_of_cells_ + Vecu(1);
		data_.resize(number_of_cells_[0] * number_of_cells_[1] * number_of_cells_[2]);
			
		// read raw data
		std::ifstream rawDataFile(rawDataFilePath, std::ios::binary);
		if(!rawDataFile.is_open())
		{
			std::cout << "\n FAILURE: open raw data file " << rawDataFilePath << " failed!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		if(dataType==DataType::UCHAR)
		{
			std::vector<unsigned char> rawData(number_of_cells_[0]* number_of_cells_[1] * number_of_cells_[2],0);
			rawDataFile.read((char *)rawData.data(), rawData.size() * sizeof(unsigned char));
			for(auto i=0; i<rawData.size(); ++i)
				data_[i] = static_cast<Real>(rawData[i]);
		}
		else if(dataType==DataType::SHORT)
		{
			std::vector<short> rawData(number_of_cells_[0]* number_of_cells_[1] * number_of_cells_[2], 0);
			rawDataFile.read((char *)rawData.data(), rawData.size() * sizeof(short));
			for(auto i=0; i<rawData.size(); ++i)
				data_[i] = static_cast<Real>(rawData[i]);
		}
		else if(dataType==DataType::FLOAT)
		{
			std::vector<float> rawData(number_of_cells_[0]* number_of_cells_[1] * number_of_cells_[2], 0);
			rawDataFile.read((char *)rawData.data(), rawData.size() * sizeof(float));
			for(auto i=0; i<rawData.size(); ++i)
				data_[i] = static_cast<Real>(rawData[i]);
		}
		else
		{
			std::cout << "\n FAILURE: The image data type is not supported!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		rawDataFile.close();
		computeGradient();
	}
	
	Real ImplicitSurface::getValueAtArbitraryPoint(Vec3d point)
	{
		Vec3d index(0, 0, 0);
		index[0] = (point[0] - dataOrigin_[0]) / grid_spacing_;
		index[1] = (point[1] - dataOrigin_[1]) / grid_spacing_;
		index[2] = (point[2] - dataOrigin_[2]) / grid_spacing_;
		// if the point is at outside of the mesh, we just set the value of signed distance to infinity. Otherwise, 
		// When constructing the level set mesh, there will be some errors.
		if (index[0]<-0.5 || index[0] > number_of_cells_[0] - 0.5 || 
		    index[1]<-0.5 || index[1] > number_of_cells_[1] - 0.5 ||
			index[2]<-0.5 || index[2] > number_of_cells_[2] - 0.5)
		{
			return  INFINITY;
		}

		Real result = 0;
		auto indexAndWeights = computeIndexAndWeights(point);
		for(auto i=0; i<8; ++i)
			result += getCellCenterData(indexAndWeights.first[i]) * indexAndWeights.second[i];
		
		return result;
	}

	Vec3d ImplicitSurface::getGradientAtArbitraryPoint(Vec3d point)
	{
		Vec3d index(0, 0, 0);
		index[0] = (point[0] - dataOrigin_[0]) / grid_spacing_;
		index[1] = (point[1] - dataOrigin_[1]) / grid_spacing_;
		index[2] = (point[2] - dataOrigin_[2]) / grid_spacing_;

		if (index[0]<-0.5 || index[0] > number_of_cells_[0] - 0.5 || 
		    index[1]<-0.5 || index[1] > number_of_cells_[1] - 0.5 ||
			index[2]<-0.5 || index[2] > number_of_cells_[2] - 0.5)
			return Vec3d(1, 0, 0);


		Real dx = 0;
		Real dy = 0;
		Real dz = 0;
		auto indexAndWeights = computeIndexAndWeights(point);
		for(auto i=0; i<8; ++i)
		{
			auto cellGrad = getCellCenterGrad(indexAndWeights.first[i]);
			dx += cellGrad[0] * indexAndWeights.second[i];
			dy += cellGrad[1] * indexAndWeights.second[i];
			dz += cellGrad[2] * indexAndWeights.second[i];
		}
	

		return Vec3d(dx, dy, dz);
	}

	bool ImplicitSurface::checkContain(const Vec3d &input_pnt, bool BOUNDARY_INCLUDED)
	{
		auto distance = findSignedDistance(input_pnt);
		if(BOUNDARY_INCLUDED)
		{
			if(distance < 0.0 || fabs(distance)<Eps)
				return true;
			else
				return false;
		}
		else
		{
			if(distance<0.0)
				return true;
			else
				return false;
		}
	}

	bool ImplicitSurface::checkNotFar(const Vec3d &input_pnt, Real threshold)
	{
		if(checkContain(input_pnt)||checkNotFar(input_pnt, threshold))
			return true;
		else
			return false;
	}

	bool ImplicitSurface::checkNearSurface(const Vec3d &input_pnt, Real threshold)
	{
		auto distance = findSignedDistance(input_pnt);
		if(fabs(distance)<threshold)
			return true;
		else 
			return false;
	}

	Real ImplicitSurface::findSignedDistance(const Vec3d &input_pnt)
	{
		auto result = getValueAtArbitraryPoint(input_pnt);
		return result;
	}

	Vec3d ImplicitSurface::findNormalDirection(const Vec3d &input_pnt)
	{
		auto grad = getGradientAtArbitraryPoint(input_pnt);
		return grad.normalize();
	}

	void ImplicitSurface::writeMeshFieldToPlt(const std::string &file)
	{
		std::ofstream output_file(file.c_str());

		output_file << "\n";
		output_file << "title='View'" << "\n";
		output_file << "variables= " << "x, " << "y, " << "z, " << "phi ";
		output_file << "zone N=" << number_of_grid_points_[0] * number_of_grid_points_[1] * number_of_grid_points_[2] 
		  			<< "  E=" << number_of_cells_[0] * number_of_cells_[1] * number_of_cells_[2]
					<< "  DATAPACKING=BLOCK"<< "\n";
		output_file << "varlocation=([4]=cellcentered)" << std::endl;
		output_file << "zonetype=febrick" << std::endl;

		for (size_t k = 0; k != number_of_grid_points_[2]; ++k)
		{
			for (size_t j = 0; j != number_of_grid_points_[1]; ++j)
			{
				for (size_t i = 0; i != number_of_grid_points_[0]; ++i)
				{
					Vecd data_position = mesh_lower_bound_;
					data_position[0] += i * grid_spacing_;
					output_file << data_position[0] << " ";
				}
				output_file << " \n";
			}
			output_file << " \n";
		}

		for (size_t k = 0; k != number_of_grid_points_[2]; ++k)
		{
			for (size_t j = 0; j != number_of_grid_points_[1]; ++j)
			{
				for (size_t i = 0; i != number_of_grid_points_[0]; ++i)
				{
					Vecd data_position = mesh_lower_bound_;
					data_position[1] += j * grid_spacing_;
					output_file << data_position[1] << " ";
				}
				output_file << " \n";
			}
			output_file << " \n";
		}

		for (size_t k = 0; k != number_of_grid_points_[2]; ++k)
		{
			for (size_t j = 0; j != number_of_grid_points_[1]; ++j)
			{
				for (size_t i = 0; i != number_of_grid_points_[0]; ++i)
				{
					Vecd data_position = mesh_lower_bound_;
					data_position[2] += k * grid_spacing_;
					output_file << data_position[2] << " ";
				}
				output_file << " \n";
			}
			output_file << " \n";
		}

		for(size_t k = 0; k != number_of_cells_[2]; ++k)
		{
			for (size_t j = 0; j != number_of_cells_[1]; ++j)
			{
				for (size_t i = 0; i != number_of_cells_[0]; ++i)
				{
					// phi
					output_file << getCellCenterData(Vec3u(i,j,k)) << " ";
				}
				output_file << " \n";
			}
			output_file << " \n";
		}

		for(size_t k = 0; k != number_of_cells_[2]; ++k)
		{
			for (size_t j = 0; j != number_of_cells_[1]; ++j)
			{
				for (size_t i = 0; i != number_of_cells_[0]; ++i)
				{
					auto one = i + j * number_of_grid_points_[0] + k * number_of_grid_points_[0] * number_of_grid_points_[1];
					auto two = one + 1;
					auto three = two + number_of_grid_points_[0];
					auto four = one + number_of_grid_points_[0];
					auto five = one + number_of_grid_points_[0] * number_of_grid_points_[1];
					auto six = five + 1;
					auto seven = six + number_of_grid_points_[0];
					auto eight = five + number_of_grid_points_[0];
					output_file << one+1 << " " << two+1 << " " << three+1 << " " << four+1 << " "
								<< five+1 << " " << six+1 << " " << seven+1 << " " << eight+1 << " ";
					output_file << " \n";
				}
			}
		}


		output_file.close();

	}


	void ImplicitSurface::writeMeshFieldToVtu(const std::string &file)
	{
		std::ofstream output_file(file.c_str());
		// add quotes around the int number
		auto addQuotes = [](size_t number){
			std::stringstream sstream;
			sstream << "\"" << number << "\"";

			return sstream.str();
		};
		//begin of the XML file
		std::string content = R"(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
	<UnstructuredGrid>)";
		output_file << content << std::endl;
		output_file << "		<Piece NumberOfPoints=" << addQuotes(number_of_grid_points_[0] * number_of_grid_points_[1] * number_of_grid_points_[2]) 
					<< " NumberOfCells=" << addQuotes(number_of_cells_[0] * number_of_cells_[1] * number_of_cells_[2]) << ">" << std::endl;
		// points of the mesh
		content = R"(			<Points>
				<DataArray type="Float32" NumberOfComponents="3" format="ascii">)";
		output_file << content << std::endl;
		for(size_t k=0; k != number_of_grid_points_[2]; ++k)
			for (size_t j = 0; j != number_of_grid_points_[1]; ++j)
				for (size_t i = 0; i != number_of_grid_points_[0]; ++i)
				{
					Vecd data_position = mesh_lower_bound_;
					data_position[0] += i * grid_spacing_;
					data_position[1] += j * grid_spacing_;
					data_position[2] += k * grid_spacing_;
					output_file << "					" << data_position[0] << " " << data_position[1] << " " << data_position[2] << "\n";
				}
			
		content = R"(				</DataArray>
      		</Points>)";
		output_file << content << std::endl;
		// cells of the mesh
		content = R"(			<Cells>
				<DataArray type="Int32" Name="connectivity" format="ascii">)";
		output_file << content << std::endl;
		for(size_t k = 0; k != number_of_cells_[2]; ++k)
			for (size_t j = 0; j != number_of_cells_[1]; ++j)
				for (size_t i = 0; i != number_of_cells_[0]; ++i)
				{
					auto one = i + j * number_of_grid_points_[0] + k * number_of_grid_points_[0] * number_of_grid_points_[1];
					auto two = one + 1;
					auto three = two + number_of_grid_points_[0];
					auto four = one + number_of_grid_points_[0];
					auto five = one + number_of_grid_points_[0] * number_of_grid_points_[1];
					auto six = five + 1;
					auto seven = six + number_of_grid_points_[0];
					auto eight = five + number_of_grid_points_[0];
					output_file << "					" << one << " " << two << " " << three << " " << four << " "
								<< five << " " << six << " " << seven << " " << eight << " ";
					output_file << " \n";
				}
			
		content = R"(				</DataArray>)";
		output_file << content << std::endl;
		content = R"(				<DataArray type="Int32" Name="offsets" format="ascii">)";
		output_file << content << std::endl;
		for(size_t i=0; i< number_of_cells_[0]*number_of_cells_[1]*number_of_cells_[2]; ++i)
			output_file << "					" << (i+1) * 8 << "\n";
		content = R"(				</DataArray>)";
		output_file << content << std::endl;
		content = R"(				<DataArray type="UInt8" Name="types" format="ascii">)";
		output_file << content << std::endl;
		for(size_t i=0; i< number_of_cells_[0]*number_of_cells_[1]*number_of_cells_[2]; ++i)
			output_file << "					" << 12 << "\n";
		content = R"(				</DataArray>)";
		output_file << content << std::endl;
		output_file << "			</Cells>" << std::endl;
		// cell data
		content = R"(			<CellData Vectors="gradient_cell">
				<DataArray type="Float32" Name="phi_cell" format="ascii">)";
		output_file << content << std::endl;
		for(size_t k = 0; k != number_of_cells_[2]; ++k)
			for (size_t j = 0; j != number_of_cells_[1]; ++j)
				for (size_t i = 0; i != number_of_cells_[0]; ++i)
			{
				// phi
				output_file << "					" << getCellCenterData(Vec3u(i,j,k)) << "\n";
			}
		content = R"(				</DataArray>)";
		output_file << content << std::endl;
		content = R"(				<DataArray type="Float32" Name="gradient_cell" NumberOfComponents="3" format="ascii">)";
		output_file << content << std::endl;
		for(size_t k = 0; k != number_of_cells_[2]; ++k)
			for (size_t j = 0; j != number_of_cells_[1]; ++j)
				for (size_t i = 0; i != number_of_cells_[0]; ++i)
			{
				// grad
				const auto &grad = getCellCenterGrad(Vec3u(i,j,k));
				output_file << "					" << grad[0] << " " << grad[1] << " " << grad[2] << "\n";
			}
		content = R"(				</DataArray>)";
		output_file << content << std::endl;
		output_file << "			</CellData>" << std::endl;
		//============================== begin point data===============================================
		content = R"(			<PointData Vectors="gradient_point">
				<DataArray type="Float32" Name="phi_point" format="ascii">)";
		output_file << content << std::endl;
		for(size_t k = 0; k != number_of_grid_points_[2]; ++k)
			for (size_t j = 0; j != number_of_grid_points_[1]; ++j)
				for (size_t i = 0; i != number_of_grid_points_[0]; ++i)
				{
					// phi
					Vecd data_position = mesh_lower_bound_;
					data_position[0] += i * grid_spacing_;
					data_position[1] += j * grid_spacing_;
					data_position[2] += k * grid_spacing_;
					output_file << "					" << getValueAtArbitraryPoint(data_position) << "\n";
				}
		content = R"(				</DataArray>)";
		output_file << content << std::endl;
		// data array for graident at grid points
		content = R"(				<DataArray type="Float32" Name="gradient_point" NumberOfComponents="3" format="ascii">)";
		output_file << content << std::endl;
		for(size_t k = 0; k != number_of_grid_points_[2]; ++k)
			for (size_t j = 0; j != number_of_grid_points_[1]; ++j)
				for (size_t i = 0; i != number_of_grid_points_[0]; ++i)
				{
					// grad
					Vecd data_position = mesh_lower_bound_;
					data_position[0] += i * grid_spacing_;
					data_position[1] += j * grid_spacing_;
					data_position[2] += k * grid_spacing_;
					const auto &grad = getGradientAtArbitraryPoint(data_position);
					output_file << "					" << grad[0] << " " << grad[1] << " " << grad[2] << "\n";
				}
		content = R"(				</DataArray>)";
		output_file << content << std::endl;
		output_file << "			</PointData>" << std::endl;
		//============================== end point data===============================================
		// end of the file
		content = R"(		</Piece>
	</UnstructuredGrid>
</VTKFile>)";
		output_file << content << std::endl;

		output_file.close();
	}

} // end of namespace SPH