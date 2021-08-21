#include "image_mesh_shape.h"

namespace SPH
{
	//=================================================================================================//
	ImageMeshShape::ImageMeshShape(std::string file_path_name) :
	origin_(0.0, 0.0, 0.0),
	translation_(0.0, 0.0, 0.0),
	rotation_(1.0),
	spacing_(1.0, 1.0, 1.0),
	dimensions_(1, 1, 1),
	data_(nullptr),
    size_(1),
    width_(1),
    height_(1),
    depth_(1),
    max_distance_(-INFINITY),
    min_distance_(INFINITY)
	{
		std::fstream dataFile(file_path_name);
		if (dataFile.fail())
		{
			std::cout << "File can not open <<" << file_path_name << std::endl;
		}

		while (!dataFile.fail() && !dataFile.eof())
		{
		}
		dataFile.close();

	}
	//=================================================================================================//
	ImageMeshShape::ImageMeshShape(Vec3d halfsize, int resolution, Vec3d translation, Mat3d rotation) :
	origin_(0.0, 0.0, 0.0),
	translation_(translation),
	rotation_(1.0),
	spacing_(1.0, 1.0, 1.0),
	dimensions_(1, 1, 1),
	data_(nullptr),
    size_(1),
    width_(1),
    height_(1),
    depth_(1),
    max_distance_(-INFINITY),
    min_distance_(INFINITY)
	{

	}
	//=================================================================================================//
	ImageMeshShape::ImageMeshShape(Real radius, int resolution, Vec3d translation, Mat3d rotation) :
	origin_(radius, radius, radius),
	translation_(translation),
	rotation_(rotation),
	spacing_(1.0, 1.0, 1.0),
	dimensions_(1, 1, 1),
	data_(nullptr),
    size_(1),
    width_(1),
    height_(1),
    depth_(1),
    max_distance_(-INFINITY),
    min_distance_(INFINITY)
	{
        // origin_ = origin;
        rotation_ = rotation;
        translation_ = translation;
		int length = int(std::ceil(3.0*radius));
        dimensions_ = Vec3i(length, length, length);
        width_ = dimensions_[0];
        height_ = dimensions_[1];
        depth_ = dimensions_[2];
        size_ = width_*height_*depth_;
        data_ = new float[size_];

		std::ofstream output_file("sphere.dat",std::ofstream::out);

		output_file << "\n";
		output_file << "title='View'" << "\n";
		output_file << "variables= " << "x, " << "y, " << "z, " << "phi " << "\n";
		output_file << "zone i=" << width_ << "  j=" << height_ << "  k=" << depth_
			<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (int z = 0; z < depth_; z++)
		{
			for (int y = 0; y < height_; y++)
			{
				for (int x = 0; x < width_; x++)
				{
					//std::cout << data_[index] << std::endl;
					output_file << x << " ";
				}
				output_file << " \n";
			}
		}

		for (int z = 0; z < depth_; z++)
		{
			for (int y = 0; y < height_; y++)
			{
				for (int x = 0; x < width_; x++)
				{
					//std::cout << data_[index] << std::endl;
					output_file << y << " ";
				}
				output_file << " \n";
			}
		}
		for (int z = 0; z < depth_; z++)
		{
			for (int y = 0; y < height_; y++)
			{
				for (int x = 0; x < width_; x++)
				{
					//std::cout << data_[index] << std::endl;
					output_file << z << " ";
				}
				output_file << " \n";
			}
		}

        for(int z = 0; z < depth_; z++)
        {
            for (int y = 0; y < height_; y++)
            {
                for (int x = 0; x < width_; x++)
                {
                    int index = z*width_*height_ + y*width_ + x;
                    double distance = (Vec3d(x,y,z) - origin_).norm()*spacing_[0]-radius*spacing_[0];
                    if(distance < min_distance_) min_distance_ = distance;
                    if(distance > max_distance_) max_distance_ = distance;
					data_[index] = float(distance);
                    /*if(distance < radius)
                    {
                        data_[index] = -float(distance);
                    }
                    else
                    {
                        data_[index] = float(distance);
                    }*/
					//std::cout << data_[index] << std::endl;
					output_file << data_[index] << " ";
                }
				output_file << " \n";
            }
        }
		output_file.close();
	}
	//=================================================================================================//
	ImageMeshShape::ImageMeshShape(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution, Vec3d translation, Mat3d rotation) :
	origin_(0.0, 0.0, 0.0),
	translation_(translation),
	rotation_(rotation),
	spacing_(1.0, 1.0, 1.0),
	dimensions_(1, 1, 1),
	data_(nullptr),
    size_(1),
    width_(1),
    height_(1),
    depth_(1),
    max_distance_(-INFINITY),
    min_distance_(INFINITY)
	{

	}
	//=================================================================================================//
	ImageMeshShape::~ImageMeshShape()
	{
		if(data_)
		{
			delete data_;
			data_ = nullptr;
		}
	}
	//=================================================================================================//
	bool ImageMeshShape::checkContain(const Vec3d& input_pnt, bool BOUNDARY_INCLUDED)
	{
        Real value = findValueAtPoint(input_pnt);
        if (BOUNDARY_INCLUDED == true)
        {
            if (value > 0.0) return false;
            else return true;
        }
        else
        {
            if (value >= 0.0) return false;
            else return true;
        }
	}
	//=================================================================================================//
	Vec3d ImageMeshShape::findClosestPoint(const Vec3d& input_pnt)
	{
        Vec3i this_cell;
        std::vector<int> neighbors = findNeighbors(input_pnt, this_cell);
        Vec3d n_sum(0.0, 0.0, 0.0);
        double weight_sum = 0.0;
        double d_sum = 0.0;
        for (const int& i : neighbors)
        {
            // checkIndexBound(i);
            Vec3d nCj = computeNormalAtCell(i);
            double dCj = float(getValueAtCell(i));
            double weight_Cj = 1.0/(fabs(dCj)+Eps);
            n_sum = n_sum + weight_Cj*nCj;
            weight_sum = weight_sum + weight_Cj;
            d_sum = d_sum + dCj;
        }
        Vec3d n = n_sum/(weight_sum+Eps);
        double d = d_sum/(weight_sum+Eps);

        Vec3d p_image = Vec3d(this_cell[0], this_cell[1], this_cell[2]) + n.normalize()*d;
        Vec3d p = convertToPhysicalSpace(p_image);
        return p;

	}
	//=================================================================================================//
	BoundingBox ImageMeshShape::findBounds()
	{
        //initial reference values
		Vec3d lower_bound = Vec3d(Infinity);
		Vec3d upper_bound = Vec3d(-Infinity);
        
        for(int z = 0; z < depth_+1; z++)
        {
            for(int y = 0; y < height_+1; y++)
            {
                for(int x = 0; x < width_+1; x++)
                {
                    Vec3d p_image = Vec3d(x, y, z);
                    Vec3d vertex_position = convertToPhysicalSpace(p_image);
			        for (int j = 0; j != 3; ++j) {
				        lower_bound[j] = SMIN(lower_bound[j], vertex_position[j]);
				        upper_bound[j] = SMAX(upper_bound[j], vertex_position[j]);
			        }
                }
            }
        }

        return BoundingBox(lower_bound, upper_bound);

	}
	//=================================================================================================//
	Real ImageMeshShape::findValueAtPoint(const Vec3d& input_pnt)
	{
        Vec3i this_cell;
        std::vector<int> neighbors = findNeighbors(input_pnt, this_cell);
        double weight_sum = 0.0;
        double d_sum = 0.0;
        if(neighbors.size() > 0)
        {
            for (const int& i : neighbors)
            {
                // checkIndexBound(i);
                double dCj = float(getValueAtCell(i));
                double weight_Cj = 1.0/(fabs(dCj)+Eps);
                weight_sum = weight_sum + weight_Cj;
                d_sum = d_sum + dCj;
            }       
            return d_sum/(weight_sum+Eps); 
        }
        else
        {
            return max_distance_;
        }

	}
	//=================================================================================================//
	Vec3d ImageMeshShape::findNormalAtPoint(const Vec3d & input_pnt)
	{
        Vec3i this_cell;
        std::vector<int> neighbors = findNeighbors(input_pnt, this_cell);
        Vec3d n_sum(0.0, 0.0, 0.0);
        double weight_sum = 0.0;
        double d_sum = 0.0;
        if( neighbors.size() > 0)
        {
            for (const int& i : neighbors)
            {
                // checkIndexBound(i);
                Vec3d nCj = computeNormalAtCell(i);
                double dCj = float(getValueAtCell(i));
                double weight_Cj = 1.0/(fabs(dCj)+Eps);
                n_sum = n_sum + weight_Cj*nCj;
                weight_sum = weight_sum + weight_Cj;
                d_sum = d_sum + dCj;
            }
            Vec3d n = n_sum/(weight_sum+Eps);
            return n.normalize(); 
        }
        else
        {
            return Vec3d(1.0,1.0,1.0).normalize();
        }

	}
	//=================================================================================================//
	std::vector<int> ImageMeshShape::findNeighbors(const Vec3d& input_pnt, Vec3i& this_cell)
	{
        std::vector<int> neighbors;

        Vec3d image_coord = rotation_.invert()*(input_pnt - translation_);
        // std::cout <<"findNeighbor of " << input_pnt << " ........... " << image_coord << std::endl;

        int z = int(floor(image_coord[2]));
        int y = int(floor(image_coord[1]));
        int x = int(floor(image_coord[0]));

		//- cannot count cells in buffer zone
		if (x < 0 || x > width_ - 1 || y < 0 || y > height_ - 1 || z < 0 || z > depth_ - 1) return neighbors;

		for (int k = z - 1; k < z + 2; k = k + 2)
		{
			for (int j = y - 1; j < y + 2; j = j + 2)
			{
				for (int i = x - 1; i < x + 2; i = i + 2)
				{
					if (i<0 || i >width_ - 1 || j <0 || j >height_ - 1 || k<0 || k >depth_)
						continue;
					int index = z * width_*height_ + y * width_ + x;
					neighbors.push_back(index);
				}
			}
				
		}

        return neighbors;

	}
	//=================================================================================================//
	Vec3d ImageMeshShape::computeGradientAtCell(int i)
	{
        //- translate 1D index to 3D index
        int width = width_;
        int height = height_;
        int depth = depth_;
        int sliceSize = width*height;
        int z = i/sliceSize;
        int y = (i%sliceSize)/width;
        int x = (i%sliceSize)%width;

        double gradx = 0.0; double grady = 0.0; double gradz = 0.0;
        //- cds (if inner cell)
        //- otherwise back/forward scheme
        if(x == 0)
        {
            int indexHigh = z*sliceSize + y*width + (x+1);
            gradx = (getValueAtCell(indexHigh) - getValueAtCell(i));
        }
        else if(x == width - 1)
        {
            int indexLow = z*sliceSize + y*width + (x-1);
            gradx = -(getValueAtCell(indexLow) - getValueAtCell(i));
        }
		else if (x > 0 && x < width_ - 1)
        {
            int indexHigh = z*sliceSize + y*width + (x+1);
            int indexLow = z*sliceSize + y*width + (x-1);
            gradx = (getValueAtCell(indexHigh) - getValueAtCell(indexLow))/2.0;
        }

        if(y == 0)
        {
            int indexHigh = z*sliceSize + (y+1)*width + x;
            grady = (getValueAtCell(indexHigh) - getValueAtCell(i));
        }
        else if(y == height - 1)
        {
            int indexLow = z*sliceSize + (y-1)*width + x;
            grady = -(getValueAtCell(indexLow) - getValueAtCell(i));
        }
		else if (y > 0 && y < height_ - 1)
        {
            int indexHigh = z*sliceSize + (y+1)*width + x;
            int indexLow = z*sliceSize + (y-1)*width + x;
            grady = (getValueAtCell(indexHigh) - getValueAtCell(indexLow))/2.0;
        }

        if(z == 0)
        {
            int indexHigh = (z+1)*sliceSize + y*width + x;
            gradz = (getValueAtCell(indexHigh) - getValueAtCell(i));
        }
        else if(z == depth - 1)
        {
            int indexLow = (z-1)*sliceSize + y*width + x;
            gradz = -(getValueAtCell(indexLow) - getValueAtCell(i));
        }
		else if (z > 0 && z < depth_ - 1)
        {
            int indexHigh = (z+1)*sliceSize + y*width + x;
            int indexLow = (z-1)*sliceSize + y*width + x;
            gradz = (getValueAtCell(indexHigh) - getValueAtCell(indexLow))/2.0;
        }
		gradx = gradx / spacing_[0];
		grady = grady / spacing_[1];
		gradz = gradz / spacing_[2];
        return Vec3d(gradx, grady, gradz);

	}
	//=================================================================================================//
	Vec3d ImageMeshShape::computeNormalAtCell(int i)
	{
        Vec3d grad_phi = computeGradientAtCell(i);
        Vec3d n = grad_phi.normalize();
        return n;
    }

	float ImageMeshShape::getValueAtCell(int i)
    {
        if(i < 0 || i > size_)
        {
            return float(max_distance_);
        }
        else
        {
            return data_[i];            
        }

	}
	//=================================================================================================//
	Vec3d ImageMeshShape::convertToPhysicalSpace(Vec3d p)
	{
        Vec3d position = rotation_*p + translation_;
        for(int i = 0; i < position.size(); i++)
        {
            position[i] = position[i]*spacing_[i];
        }
        return position;

	}
	//=================================================================================================//

}