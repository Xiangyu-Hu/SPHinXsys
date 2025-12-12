//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

#include "base_data_type.h"
#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "large_data_containers.h"
#include "pressure_boundary.h"
#include "sphinxsys.h"
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <functional>
#include <gtest/gtest.h>
#include <numbers>
#include <ostream>
#include <string>
#include <vector>
using namespace SPH;
const Real t_ref = 0.025;
//----------------------------------------------------------------------
//	Circular buffer for checking convergence usage
//----------------------------------------------------------------------
template <typename T>
class CircularBuffer
{
  private:
    std::vector<T> buffer;
    size_t head = 0;
    size_t capacity;
    size_t count = 0;

  public:
    CircularBuffer(size_t size) : capacity(size)
    {
        buffer.resize(size);
    }

    void push(const T &item)
    {
        buffer[head] = item;
        head = (head + 1) % capacity;
        if (count < capacity)
            count++;
    }

    T &operator[](size_t index)
    {
        return buffer[(head - count + index) % capacity];
    }

    size_t size() const
    {
        return count;
    }

    bool is_full() const
    {
        return count == capacity;
    }

    T get_average() const
    {
        if (count < capacity)
            return T::Constant(std::numeric_limits<typename T::Scalar>::infinity()); // Return a vector with all elements as infinity if the buffer isn't full

        T sum = T::Zero(); // Initialize sum as a zero vector of appropriate dimensions
        for (size_t i = 0; i < count; ++i)
        {
            sum += buffer[i];
        }
        return sum / count; // Divide each element of sum by the count to get the average
    }
};
//----------------------------------------------------------------------
//	ConvergenceChecker
//----------------------------------------------------------------------
template <typename T>
class ConvergenceChecker
{
  private:
    CircularBuffer<T> buffer;
    T threshold;
    T percentage_difference_ = std::numeric_limits<T>::infinity();

    bool calculate_convergence()
    {
        if (!buffer.is_full())
        {
            return false; // Not enough data to determine convergence
        }

        T max_val = buffer[0];
        T min_val = buffer[0];
        for (size_t i = 1; i < buffer.size(); ++i)
        {
            if (buffer[i] > max_val)
                max_val = buffer[i];
            if (buffer[i] < min_val)
                min_val = buffer[i];
        }

        // Calculate the percentage difference
        T range = max_val - min_val;
        T average = (max_val + min_val) / 2;
        T percentage_difference = (range / average) * 100;
        percentage_difference_ = percentage_difference;
        std::cout << "converger percentage_difference_ :" << percentage_difference_ << "\n";
        return abs(percentage_difference) < threshold;
    }

  public:
    ConvergenceChecker(size_t size, T conv_threshold)
        : buffer(size), threshold(conv_threshold) {}

    bool update(T new_value)
    {
        buffer.push(new_value);
        return calculate_convergence();
    }

    T get_percentage_difference()
    {
        return percentage_difference_;
    }
};
struct AxialVelocityProfile
{
    double z;
    std::vector<std::pair<double, double>> velocity_data; // Pair of position and velocity

    // Default constructor
    AxialVelocityProfile() : z(0)
    {
    }

    // Constructor with initial z-value
    AxialVelocityProfile(double z) : z(z)
    {
    }

    void addData(double position, double velocity) { velocity_data.emplace_back(position, velocity); }
};

void readDataFromFile(
    const std::string &filename,
    std::string string_pattern,
    std::map<double, AxialVelocityProfile> &velocity_profiles)
{
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Failed to open file" << std::endl;
        return;
    }

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string token;

        // Look for the specific pattern indicating the start of the data we care about
        if (line.find(string_pattern) != std::string::npos)
        {
            double z_value;
            int count;

            // Read and discard the first token ("plot-profile-axial-velocity-at-z")
            iss >> token;   // Consumes "plot-profile-axial-velocity-at-z"
            iss >> z_value; // Next token should be the z value

            // Debug output
            std::cout << "Found profile for Z-value: " << z_value << std::endl;

            // Move to the next line to read the count of data points
            if (std::getline(file, line))
            {
                std::istringstream iss_count(line);
                iss_count >> count; // First integer on the next line should be the count

                std::cout << "Data count: " << count << std::endl;

                AxialVelocityProfile profile(z_value);

                // Read the actual data points
                for (int i = 0; i < count; ++i)
                {
                    double position, velocity;
                    if (file >> position >> velocity)
                    {
                        profile.addData(position, velocity);
                    }
                }

                // Use emplace to avoid unnecessary object creation
                auto result = velocity_profiles.emplace(z_value, std::move(profile));
                if (!result.second)
                {
                    std::cerr << "Profile for z=" << z_value << " already exists." << std::endl;
                }
            }
            else
            {
                std::cerr << "Expected data count line after z-value but got nothing." << std::endl;
            }
        }
    }
};

StdVec<Vecd> ObersverAxialGenerator(double x_min, double x_max, size_t number_of_axial_observer)
{
    StdVec<Vecd> observer_location;
    int ny = number_of_axial_observer + 1;
    double full_length = x_max - x_min;
    for (int i = 0; i < ny; ++i)
    {
        double x = full_length / (ny - 1) * i + x_min;
        observer_location.emplace_back(Vec3d(x, 0, 0));
    }
    return observer_location;
}

StdVec<Vecd> ObserverRadialGenerator(double x, double diameter)
{
    StdVec<Vecd> observer_location;

    int ny = 51;
    for (int i = 0; i < ny - 1; i++)
    {
        double z = diameter * i / (ny - 1);
        observer_location.emplace_back(Vec3d(x, 0, -diameter / 2.0 + z));
    }
    return observer_location;
}

double compute_pressure(double p)
{
    double run_time = GlobalStaticVariables::physical_time_;
    double pressure = run_time < 0.5 ? 0.5 * p * (1.0 - cos(M_PI * run_time / 0.5)) : p;

    // double pressure = p;
    return pressure;
}

struct FDA_nozzle_parameters
{ // using default blood parameters
    // resolution, particles per diameter
    double scale = 0.001;
    double inlet_diameter = 12 * scale; // Entering part diameter
    double thoat_diameter = 4 * scale;
    double outlet_diameter = 12 * scale;

    uint number_of_particles = 10;

    Vec3d inlet_normal = Vec3d::UnitX();
    Vec3d inlet_center = Vec3d(-120 * scale, 0, 0);

    Vec3d outlet_normal = Vec3d::UnitX() * -1;
    Vec3d outlet_center = Vec3d(120 * scale, 0, 0);

    // defualt inlet area
    double inlet_area = pow(inlet_diameter, 2) * M_PI * 0.25;
    double throat_area = pow(thoat_diameter, 2) * M_PI * 0.25;

    double Q_f = 5e-6; // for Re = 500

    // material paramaters, default: blood
    // Newtonian
    double rho0_f = 1056; //[kg/m3]
    double mu_f = 0.0035; //[N s / m2] from https://arxiv.org/pdf/2204.10566 3.5 mPa s

    // fluid flle path
    fs::path fluid_file_path = "./input/Fluid_expansion/Fluid_expansion.stl";
    // wall flle path
    fs::path wall_file_path;
    // total length of fluid
    double fluid_length = 200 * scale; // defined in mesh file
    // end time
    double end_time = 2.0;
};

void setup_directory(const std::string &path)
{
    // Check if the directory exists
    if (fs::exists(path))
    {
        // Try to remove the directory and its contents
        if (!fs::remove_all(path))
        {
            std::cerr << "Failed to remove existing directory: " << path << std::endl;
            return;
        }
    }

    // Create the directory
    if (!fs::create_directory(path))
    {
        std::cerr << "Failed to create directory: " << path << std::endl;
        return;
    }

    std::cout << "Directory set up at: " << path << std::endl;
}

void print_profiles_data(const std::map<double, AxialVelocityProfile> &profiles, const std::string &description)
{
    std::cout << description << std::endl;
    for (const auto &entry : profiles)
    {
        double z = entry.first;
        const AxialVelocityProfile &profile = entry.second;

        std::cout << "Profile for Z-value: " << z << std::endl;
        for (const auto &data : profile.velocity_data)
        {
            std::cout << "    Position: " << data.first << ", Value: " << data.second << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
    }
}

const AxialVelocityProfile get_profiles_data(const std::map<double, AxialVelocityProfile> &profiles, double z)
{
    auto it = profiles.find(z); // Use the find method to search for the z value directly
    return (it != profiles.end()) ? it->second
                                  : throw std::runtime_error("Profile not found for z-value: " + std::to_string(z));
}

/**
 * @brief 	Pressure boundary definition.
 */
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        return p_;
    }
};

struct RightInflowPressure
{
    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real &p_)
    {
        /*constant pressure*/
        Real pressure = 0.;
        return pressure;
    }
};

Real U_f = 0.;
Real inlet_diameter = 0;
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
// TODO: add x_max in emitter zone to determined effectiness
class TimeDependentAcceleration : public Gravity
{
    AlignedBoxShape &aligned_box_;
    Real t_ref_, u_ref_, du_ave_dt_;

  public:
    explicit TimeDependentAcceleration(BodyAlignedBoxByCell &aligned_box_part, int axis, Vecd gravity_vector, Real t_ref, Real u_ref)
        : Gravity(gravity_vector), aligned_box_(aligned_box_part.getAlignedBoxShape()), t_ref_(t_ref), u_ref_(u_ref), du_ave_dt_(0)
    {
    }

    virtual Vecd InducedAcceleration(const Vecd &position) override
    {
        // if (aligned_box_.checkUpperBound(position))
        // {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);

        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0, 0.0) : global_acceleration_;
        // }
        // else
        // {
        //     return global_acceleration_;
        // }
    }
};
//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_max_;
    Real radius_squared_; // radius at inlet
    Real t_ref_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_max_(2 * U_f), radius_squared_(0.0), t_ref_(t_ref)
    {
        std::cout << "u_max_ in parabolic flow is " << u_max_ << std::endl;
    }

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vec3d pos = position;
        pos[0] = 0; // set x direction to 0
        double radius_square = pos.squaredNorm();
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_max = run_time < t_ref_ ? 0.5 * u_max_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_max_;
        radius_squared_ = pow(inlet_diameter * 0.5, 2);
        double vel = u_max * (1 - radius_square / radius_squared_);
        Vec3d target_velocity = Vec3d::UnitX() * vel;

        return target_velocity;
    }

    void set_u_max(double u_max)
    {
        u_max_ = u_max;
    }

    void set_radius_squared_(double radius_squared)
    {
        radius_squared_ = radius_squared;
    }
};

class ObserverRadial : public ObserverBody
{
  public:
    ObserverRadial(SPHSystem &sph_system, std::string name, double x, double diameter) : ObserverBody(sph_system, name) // Assuming the observer has a predefined shape or characteristic name.
    {
        this->generateParticles<ObserverParticles>(ObserverRadialGenerator(x, diameter));
    }
};

// Function template for writing CSV data
template <typename T>
void write_data_to_csv(std::ofstream &csv_output, const StdLargeVec<Vecd> &positions, const StdLargeVec<T> &properties, const std::string &property_name)
{
    // Write the header for Vecd type
    if constexpr (std::is_same_v<T, Vecd>)
    {
        csv_output << "Position X,Position Y,Position Z," << property_name << " X," << property_name << " Y," << property_name << " Z\n";
        for (size_t i = 0; i < positions.size(); i++)
        {
            csv_output << positions[i][0] << "," << positions[i][1] << "," << positions[i][2] << ","
                       << properties[i][0] << "," << properties[i][1] << "," << properties[i][2] << "\n";
        }
    }
    // Write the header for Real type
    else if constexpr (std::is_same_v<T, Real>)
    {
        csv_output << "Position X,Position Y,Position Z," << property_name << "\n";
        for (size_t i = 0; i < positions.size(); i++)
        {
            csv_output << positions[i][0] << "," << positions[i][1] << "," << positions[i][2] << ","
                       << properties[i] << "\n";
        }
    }
    csv_output.flush();
}

void write_observer_properties_to_output(ObserverBody &observer_body, const std::string property_name, fs::path output_path)
{
    // Retrieve positions
    StdLargeVec<Vecd> &pos = observer_body.getBaseParticles().ParticlePositions();

    // Determine the type of the property based on the name
    if (property_name == "Velocity" || property_name == "avg_velocity")
    { // Example Vecd properties
        // TODO: need to do a precaution for invalid property_ptr
        auto &property_ptr = *observer_body.getBaseParticles().getVariableDataByName<Vecd>(property_name);

        std::ofstream csv_output(output_path / ("observer_" + observer_body.getName() + "_" + property_name + ".csv"));
        write_data_to_csv(csv_output, pos, property_ptr, property_name);
        csv_output.close();
    }
    else if (property_name == "Pressure")
    { // Example Real properties
        // TODO: need to do a precaution for invalid property_ptr
        auto &property_ptr = *observer_body.getBaseParticles().getVariableDataByName<Real>(property_name);

        std::ofstream csv_output(output_path / ("observer_" + observer_body.getName() + "_" + property_name + ".csv"));
        write_data_to_csv(csv_output, pos, property_ptr, property_name);
        csv_output.close();
    }
    else
    {
        throw std::runtime_error("Invalid property name: " + property_name);
    }
}

class RadialObserverContainer
{

  private:
    SPHSystem &sph_system_;
    std::vector<std::unique_ptr<ObserverRadial>> observers_;
    std::vector<std::unique_ptr<ContactRelation>> contact_relations_;
    std::vector<std::unique_ptr<ObservingAQuantity<Vec3d>>> observing_quantities_;
    std::vector<std::unique_ptr<BodyStatesRecordingToVtp>> body_state_recording_;
    std::vector<std::unique_ptr<CircularBuffer<Vecd>>> circular_buffer_avg_velocity_;
    std::vector<double> z_;

  public:
    RadialObserverContainer(SPHSystem &sph_system) : sph_system_(sph_system)
    {
    }
    const std::vector<double> get_position() { return z_; }

    void add_radial_observer(double x, double diameter, FluidBody &fluid_body, SolidBody &wall_body)
    {
        // record the position into vector
        z_.emplace_back(x);

        auto observer_radial =
            std::make_unique<ObserverRadial>(sph_system_, "ObserverRadial_" + std::to_string(x), x, diameter);
        observer_radial->getBaseParticles().registerSharedVariable<Vecd>("avg_velocity");
        observer_radial->getBaseParticles().registerSharedVariable<int>("circular_buffer_index");
        observers_.push_back(std::move(observer_radial));

        // Convert SPHBody pointers to RealBody pointers
        std::vector<RealBody *> real_body_ptrs;
        real_body_ptrs.push_back(static_cast<RealBody *>(&fluid_body));
        real_body_ptrs.push_back(static_cast<RealBody *>(&wall_body));

        // Now use std::make_unique with the correct type
        contact_relations_.emplace_back(std::make_unique<ContactRelation>(*observers_.back(), real_body_ptrs));

        auto observer_contact = contact_relations_.back().get(); // Assuming this is a correct type derivation or has access to BaseContactRelation
        auto observing_velocity = std::make_unique<ObservingAQuantity<Vec3d>>(*observer_contact, "Velocity");
        observing_quantities_.push_back(std::move(observing_velocity));
        auto body_state_recording = std::make_unique<BodyStatesRecordingToVtp>(*observers_.back());
        body_state_recording->addToWrite<Vec3d>(*observers_.back(), "Velocity");
        body_state_recording->addToWrite<Vec3d>(*observers_.back(), "avg_velocity");
        body_state_recording->addToWrite<int>(*observers_.back(), "circular_buffer_index");
        body_state_recording_.push_back(std::move(body_state_recording));
    }

    ObserverBody *get_observer(size_t index)
    {
        if (index < observers_.size())
            return observers_[index].get();
        return nullptr;
    }

    ContactRelation *get_contact_relation(size_t index)
    {
        if (index < contact_relations_.size())
            return contact_relations_[index].get();
        return nullptr;
    }

    void update_configureations()
    {
        for (auto &contact : contact_relations_)
        {
            contact->updateConfiguration(); // Call the update function for each contact relation
        }
    }

    void init_avg_velocity()
    {
        circular_buffer_avg_velocity_.clear();
        for (auto &observer : observers_)
        {
            {
                auto &circular_buffer_indices = *observer->getBaseParticles().getVariableDataByName<int>("circular_buffer_index");
                // Create a buffer for each particle and assign an index.
                for (size_t i = 0; i < circular_buffer_indices.size(); ++i)
                {
                    // Create a new buffer for this particle.
                    circular_buffer_avg_velocity_.emplace_back(std::make_unique<CircularBuffer<Vecd>>(10));

                    // Assign the index of this new buffer to the particle's index.
                    circular_buffer_indices[i] = circular_buffer_avg_velocity_.size() - 1;
                }
            }
        }
    }

    void update_avg_velocity()
    {
        for (auto &observer : observers_)
        {
            auto &velocity = *observer->getBaseParticles().getVariableDataByName<Vecd>("Velocity");
            auto &avg_velocity = *observer->getBaseParticles().getVariableDataByName<Vecd>("avg_velocity");
            auto &circular_buffer_indices = *observer->getBaseParticles().getVariableDataByName<int>("circular_buffer_index");
            for (size_t i = 0; i < circular_buffer_indices.size(); ++i)
            {
                auto circular_buffer_indice = circular_buffer_indices[i];
                auto vel = velocity[i];
                circular_buffer_avg_velocity_[circular_buffer_indice]->push(vel);
                avg_velocity[i] = circular_buffer_avg_velocity_[circular_buffer_indice]->get_average();
            }
        }
    }

    ObservingAQuantity<Vec3d> *get_observing_quantity(size_t index)
    {
        if (index < observing_quantities_.size())
            return observing_quantities_[index].get();
        return nullptr;
    }

    void update_observing_quantities()
    {
        for (auto &quantity : observing_quantities_)
        {
            quantity->exec(); // Execute the observing function for each quantity
        }
    }

    void writeToFile()
    {
        for (auto &body_state_recording : body_state_recording_)
        {
            body_state_recording->writeToFile(); // Write the body state output
        }
    }

    void write_velocity_to_csv(fs::path output_path)
    {
        for (auto &observers : observers_)
        {
            write_observer_properties_to_output(*observers, "Velocity", output_path);
            write_observer_properties_to_output(*observers, "avg_velocity", output_path);
        }
    }
};

int FDA_nozzle(int ac, char *av[], std::string additional_output_name, bool is_TVLCCC, bool is_wallrelaxation, bool is_linearVisc, bool is_refineWall, bool is_oldoutlet, FDA_nozzle_parameters &params, const double Re = 500, const size_t number_of_axial_observer = 50)
{
    // Map to hold each profile accessed by the z-value as a key (obtaining radial velocity)
    std::map<double, AxialVelocityProfile> radial_velocity_profiles;
    std::map<double, AxialVelocityProfile> axial_velocity_profiles;
    std::map<double, AxialVelocityProfile> axial_pressure_profiles;

    readDataFromFile(
        "./input/PIV_Sudden_Expansion_500_999.txt", "plot-profile-axial-velocity-at-z", radial_velocity_profiles);
    readDataFromFile(
        "./input/PIV_Sudden_Expansion_500_999.txt", "plot-z-distribution-axial-velocity", axial_velocity_profiles);
    readDataFromFile(
        "./input/PIV_Sudden_Expansion_500_999.txt", "plot-z-distribution-pressure", axial_pressure_profiles);
    GlobalStaticVariables::physical_time_ = 0;

    const double scale = params.scale;
    inlet_diameter = params.thoat_diameter;
    const double Q_f = params.Q_f;
    U_f = Q_f / params.throat_area;
    double U_max = 0.8; // From exp data
    double c_f = U_max * 10;
    // double U_max = Q_f / params.throat_area * 2.5;                                      // U max suppose to happend at throat
    // const double max_pressure_diff = SMAX(250., 450.);                                  // maximum pressure difference found from exp data is 250. low rest max pressure diff is about 450
    // double c_f = SMAX(10. * U_max, pow(max_pressure_diff / 0.01 / params.rho0_f, 0.5)); // Speed of sound
    // U_max = SMAX(c_f * 0.1, U_max);

    std::cout << "U_f (velocity at inlet): " << U_f << std::endl;
    std::cout << "U_max                  : " << U_max << std::endl;
    std::cout << "Re                     : " << Re << std::endl;

    // radial observer at x axis
    std::vector<double> radiao_observer_X;
    for (const auto &entry : radial_velocity_profiles)
    {
        double z = entry.first;
        // Print the Z-value
        std::cout << "Profile for Z-value: " << z << std::endl;
        radiao_observer_X.emplace_back(z);
    }

    // Inside your FDA_nozzle function or the main function after data has been loaded
    std::cout << "Printing Axial Velocity Profiles Data:" << std::endl;
    print_profiles_data(axial_velocity_profiles, "Axial Velocity Profiles:");

    std::cout << "Printing Axial Pressure Profiles Data:" << std::endl;
    print_profiles_data(axial_pressure_profiles, "Axial Pressure Profiles:");

    const uint number_of_particles = params.number_of_particles;
    const double resolution_ref = params.thoat_diameter / number_of_particles;
    const double wall_thickness = 4 * resolution_ref;
    std::cout << "resolution_ref: " << resolution_ref << std::endl;

    // const double resolution_wall = resolution_ref;
    // const uint simtk_resolution = 10;

    const double rho0_f = params.rho0_f;
    const double mu_f = params.mu_f;

    // GEOMETRY
    auto fluid_shape = makeShared<ComplexShape>("fluid");
    fluid_shape->add<TriangleMeshShapeSTL>(params.fluid_file_path, Vec3d::Zero(), scale);

    const double x_min_domain = fluid_shape->getBounds().first_.x();
    const double x_max_domain = fluid_shape->getBounds().second_.x();
    std::cout << "x_min_domain" << x_min_domain << std::endl;
    std::cout << "x_max_domain" << x_max_domain << std::endl;

    Real boundary_width = 3. * resolution_ref;
    Real half_boundary_width = boundary_width * 0.5;
    Vecd bidirectional_buffer_halfsize = Vec3d(half_boundary_width, params.inlet_diameter * 0.5, params.inlet_diameter * 0.5);
    Vec3d left_disposer_translation = Vec3d(x_min_domain + half_boundary_width, 0, 0);
    Vec3d right_disposer_translation = Vec3d(x_max_domain - half_boundary_width, 0, 0);

    Vecd left_bidirectional_translation = Vec3d(x_min_domain + half_boundary_width, 0, 0);
    Vecd right_bidirectional_translation = Vec3d(x_max_domain - half_boundary_width, 0, 0);

    Real old_emitter_boundary_width = 20. * resolution_ref;
    Real half_old_emitter_boundary_width = old_emitter_boundary_width * 0.5;
    Vecd old_emitter_halfsize = Vec3d(half_old_emitter_boundary_width, params.inlet_diameter * 0.5, params.inlet_diameter * 0.5);
    Vecd left_old_emitter_translation = Vec3d(x_min_domain + half_old_emitter_boundary_width, 0, 0);

    Real old_boundary_width = 20. * resolution_ref;
    Real half_old_boundary_width = old_boundary_width * 0.5;
    Vecd old_buffer_halfsize = Vec3d(half_old_boundary_width, params.inlet_diameter * 0.5, params.inlet_diameter * 0.5);
    Vecd left_old_buffer_translation = Vec3d(x_min_domain + half_old_boundary_width, 0, 0);

    auto wall_shape = makeShared<ComplexShape>("wall");
    wall_shape->add<TriangleMeshShapeSTL>(params.wall_file_path, Vec3d::Zero(), scale);

    // SYSTEM
    BoundingBox bbox = wall_shape->getBounds();

    bbox.first_[0] = x_min_domain;
    bbox.first_ -= Vec3d::Ones() * wall_thickness;
    bbox.second_ += Vec3d::Ones() * wall_thickness;
    SPHSystem sph_system(bbox, resolution_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(true);

    /** handle command line arguments. */
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();

    IOEnvironment in_output(sph_system);
    in_output.output_folder_ += "_ExpanOnly_" + additional_output_name;
    if (is_TVLCCC)
        in_output.output_folder_ += "_TVLCCC_";
    else
        in_output.output_folder_ += "_TVCCC_";
    if (is_linearVisc)
        in_output.output_folder_ += "_ViscCorrection_";
    else
        in_output.output_folder_ += "_NoViscCorrection_";
    if (is_refineWall)
        in_output.output_folder_ += "_refinewall_";
    else
        in_output.output_folder_ += "_Norefinewall_";
    if (is_wallrelaxation)
        in_output.output_folder_ += "_Wallrelaxation_";
    else
        in_output.output_folder_ += "_NoWallrelaxation_";
    in_output.output_folder_ += std::to_string(params.number_of_particles);
    setup_directory(in_output.output_folder_);

    // Check if the folder exists
    std::filesystem::path folderPath = in_output.output_folder_;
    if (std::filesystem::exists(folderPath) && std::filesystem::is_directory(folderPath))
    {
        // Delete the folder if it exists
        std::filesystem::remove_all(folderPath);
    }

    // Create the folder
    std::filesystem::create_directory(folderPath);
    // FLUID
    FluidBody water_block(sph_system, fluid_shape);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.75);
    water_block.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet(0);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, wall_shape);
    if (is_refineWall)
        wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet(0);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    InnerRelation wall_boundary_inner(wall_boundary);
    //----------------------------------------------------------------------
    //	SPH Particle relaxation section
    //----------------------------------------------------------------------
    /** check whether run particle relaxation for body fitted particle distribution. */
    if (is_wallrelaxation)
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(wall_boundary);
        /** Write the body state to Vtp file. */
        // BodyStatesRecordingToVtp write_wall_body_to_vtp(wall_boundary);
        // Write the particle reload files.
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        // A  Physics relaxation step.
        // RelaxationStepInner relaxation_step_inner(wall_boundary_inner);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);

        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        // write_wall_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Particle relaxation loop.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                // write_wall_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        // Output particles position for reload.
        // write_particle_reload_files.writeToFile(0);
    }
    //----------------------------------------------------------------------
    //	Topology
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, {&wall_boundary});
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    BodyAlignedBoxByCell left_time_dependent_region(
        water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vecd(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    TimeDependentAcceleration time_dependent_acceleration(left_time_dependent_region, xAxis, Vecd::Zero(), t_ref, U_f);
    SimpleDynamics<GravityForce> apply_gravity_force(water_block, time_dependent_acceleration);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    /** zeroth order consistency */
    InteractionDynamics<NablaWVComplex> kernel_summation(water_block_inner, water_block_contact);
    /** surface particle identification */
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        boundary_indicator(water_block_inner, water_block_contact);

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> kernel_correction_complex(water_block_inner, water_block_contact);

    /** momentum equation. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    /** mass equation. */
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall>
        viscous_acceleration(water_block_inner, water_block_contact);
    /** Computing viscous acceleration. */
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWallCorrection>
        viscous_acceleration_correction(water_block_inner, water_block_contact);
    /** Impose transport velocity. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityLimitedCorrectionCorrectedComplex<BulkParticles>>
        transport_velocity_correction_TVLCCC(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionCorrectedComplex<BulkParticles>>
        transport_velocity_correction_TVCCC(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);

    //----------------------------------------------------------------------
    //	Set up boundary condition
    //----------------------------------------------------------------------
    /** delete outflow particles */
    BodyAlignedBoxByCell left_disposer(
        water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(Vecd(-1, 0, 0))), Vecd(0, -1, 0)), Vecd(left_disposer_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> left_disposer_outflow_deletion(left_disposer);
    BodyAlignedBoxByCell right_disposer(
        water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vecd(right_disposer_translation)), bidirectional_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> right_disposer_outflow_deletion(right_disposer);
    /** bidirectional buffer */
    BodyAlignedBoxByCell left_emitter(
        water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vecd(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    // fluid_dynamics::BidirectionalBuffer<LeftInflowPressure, SequencedPolicy> left_emitter_inflow_injection(left_emitter, in_outlet_particle_buffer, xAxis);
    fluid_dynamics::NonPrescribedPressureBidirectionalBuffer left_emitter_inflow_injection(left_emitter, in_outlet_particle_buffer);
    BodyAlignedBoxByCell right_emitter(
        water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation3d(std::acos(Eigen::Vector3d::UnitX().dot(Vecd(-1, 0, 0))), Vecd(0, -1, 0)), Vecd(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_emitter_inflow_injection(right_emitter, in_outlet_particle_buffer); /** output parameters */

    BodyAlignedBoxByParticle old_emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec3d(left_old_emitter_translation)), old_emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> old_emitter_inflow_injection(old_emitter, in_outlet_particle_buffer);
    BodyAlignedBoxByCell old_inflow_buffer(
        water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vecd(left_old_buffer_translation)), old_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(old_inflow_buffer);

    /** pressure boundary condition. */
    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    /** density correction in pressure-driven flow */
    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_inner, water_block_contact);

    //----------------------------------------------------------------------
    //	Observer
    //----------------------------------------------------------------------
    // Axial
    ObserverBody observer_axial(sph_system, "fluid_observer_axial");
    observer_axial.generateParticles<ObserverParticles>(ObersverAxialGenerator(x_min_domain, x_max_domain, number_of_axial_observer));
    ContactRelation axial_velocity_observer_contact(observer_axial, {&water_block});
    ObservingAQuantity<Vecd>
        update_axial_observer_velocity(axial_velocity_observer_contact, "Velocity");
    ObservingAQuantity<Real>
        update_axial_observer_pressure(axial_velocity_observer_contact, "Pressure");
    // Radial
    RadialObserverContainer observer_radial_container(sph_system);
    for (auto &pos : radiao_observer_X)
    {
        observer_radial_container.add_radial_observer(
            pos, params.inlet_diameter, water_block, wall_boundary);
    }
    observer_radial_container.init_avg_velocity();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp observer_axial_recording(observer_axial);
    observer_axial_recording.addToWrite<Vecd>(observer_axial, "Velocity");
    observer_axial_recording.addToWrite<Real>(observer_axial, "Pressure");
    BodyStatesRecordingToVtp water_body_recording(water_block);
    water_body_recording.addToWrite<Real>(water_block, "Pressure");
    water_body_recording.addToWrite<int>(water_block, "Indicator");
    water_body_recording.addToWrite<Real>(water_block, "Density");
    water_body_recording.addToWrite<int>(water_block, "BufferParticleIndicator");
    water_body_recording.addToWrite<Vec3d>(water_block, "KernelSummation");
    BodyStatesRecordingToVtp wall_body_recording(wall_boundary);
    wall_body_recording.addToWrite<Vec3d>(wall_boundary, "NormalDirection");
    //----------------------------------------------------------------------
    //	Defined convergence checker
    //----------------------------------------------------------------------
    ConvergenceChecker<double> conv_checker(10, 0.5); // Buffer of 10 values, convergence threshold of 0.5 percent
    bool is_converged = false;
    size_t mid_index_of_observer = number_of_axial_observer / 2.0;
    size_t first_index_of_observer = 0; // the observer near by interface
    size_t last_index_of_observer = 0;  // the observer near by interface
    {
        Real init_dist_x = std::numeric_limits<Real>::max(); // Initialize with max value to ensure any first comparison is smaller
        Real last_dist_x = std::numeric_limits<Real>::max(); // Initialize with max value for correct comparison
        for (size_t index = 0; index < number_of_axial_observer; index++)
        {
            auto x = observer_axial.getBaseParticles().ParticlePositions()[index][0];
            Real dist_to_min = abs(x_min_domain + boundary_width - x);
            Real dist_to_max = abs(0.088 - x);

            if (init_dist_x > dist_to_min)
            {
                init_dist_x = dist_to_min;
                first_index_of_observer = index;
            }
            if (last_dist_x > dist_to_max)
            {
                last_dist_x = dist_to_max;
                last_index_of_observer = index;
            }
        }
        std::cout << "first_index_of_observer is " << first_index_of_observer << ", distance is " << init_dist_x
                  << ", x location is " << observer_axial.getBaseParticles().ParticlePositions()[first_index_of_observer][0]
                  << "\nlast_index_of_observer is " << last_index_of_observer << ", distance is " << last_dist_x
                  << ", x location is " << observer_axial.getBaseParticles().ParticlePositions()[last_index_of_observer][0] << std::endl;
    }
    // StdLargeVec<Vecd> &pos_radial = observer_axial.getBaseParticles().ParticlePositions();
    StdLargeVec<Vecd> &vel_radial = *observer_axial.getBaseParticles().getVariableDataByName<Vecd>("Velocity");

    auto &vel_of_mid_index_observer = vel_radial[mid_index_of_observer][0];
    auto &vel_of_first_index_observer = vel_radial[first_index_of_observer][0];
    auto &vel_of_last_index_observer = vel_radial[last_index_of_observer][0];
    // auto &pos_y_of_mid_index_observer = pos_radial[mid_index_of_observer][1];

    size_t convergence_checker_output_interval = 10;
    std::ofstream file_mid_observer(sph_system.getIOEnvironment().output_folder_ + "/output_velocity_of_mid_observer_NP" + std::to_string(number_of_particles) + ".csv");
    file_mid_observer << "Time,Velocity X,Velocity X at first index,Velocity X at last index,Convergence Rate\n";
    file_mid_observer.flush(); // Ensure data is written to the file immediately
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    if (!is_oldoutlet)
    {
        left_emitter_inflow_injection.tag_buffer_particles.exec();
        right_emitter_inflow_injection.tag_buffer_particles.exec();
    }
    wall_boundary_normal_direction.exec();
    GlobalStaticVariables::physical_time_ = 0.0;
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = params.end_time; /**< End time. */
    Real Output_Time = 0.00005;      /**< Time stamps for output of body states. */
    Real dt = 0.0;                   /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    update_axial_observer_velocity.exec();
    update_axial_observer_pressure.exec();
    observer_axial_recording.writeToFile();
    water_body_recording.writeToFile();
    wall_body_recording.writeToFile();
    kernel_correction_complex.exec();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            apply_gravity_force.exec();
            time_instance = TickCount::now();
            Real Dt = get_fluid_advection_time_step_size.exec();
            boundary_indicator.exec();
            if (is_oldoutlet)
                update_density_by_summation.exec();
            else
                update_fluid_density.exec();
            kernel_correction_complex.exec();
            if (is_linearVisc)
                viscous_acceleration_correction.exec();
            else
                viscous_acceleration.exec();

            if (is_oldoutlet)
                transport_velocity_correction.exec();
            else
            {
                if (is_TVLCCC)
                    transport_velocity_correction_TVLCCC.exec();
                else
                    transport_velocity_correction_TVCCC.exec();
            }
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                pressure_relaxation.exec(dt);
                if (!is_oldoutlet)
                {
                    kernel_summation.exec();
                    left_inflow_pressure_condition.exec(dt);
                    right_inflow_pressure_condition.exec(dt);
                    inflow_velocity_condition.exec();
                }
                density_relaxation.exec(dt);
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                if (is_oldoutlet)
                    parabolic_inflow.exec();
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "	convergence = " << conv_checker.get_percentage_difference() << "    is_converged = " << is_converged << "\n";
            }

            number_of_iterations++;

            time_instance = TickCount::now();

            if (!is_oldoutlet)
            {
                left_emitter_inflow_injection.injection.exec();
                right_emitter_inflow_injection.injection.exec();
            }
            else
            {
                old_emitter_inflow_injection.exec();
            }
            // left_disposer_outflow_deletion.exec();
            // right_disposer_outflow_deletion.exec();
            // water_block.updateCellLinkedListWithParticleSort(100);
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
            // boundary_indicator.exec();
            if (!is_oldoutlet)
            {
                left_emitter_inflow_injection.tag_buffer_particles.exec();
                right_emitter_inflow_injection.tag_buffer_particles.exec();
            }
            if (number_of_iterations % convergence_checker_output_interval == 0)
            {
                axial_velocity_observer_contact.updateConfiguration();
                update_axial_observer_velocity.exec();
                is_converged = conv_checker.update(vel_of_mid_index_observer);
                std::cout << "add value to analysis convergence :" << vel_of_mid_index_observer << std::endl;
                if (is_converged)
                {
                    std::cout << "Converged at iteration " << vel_of_mid_index_observer << std::endl;
                }
                if (std::isfinite(conv_checker.get_percentage_difference()))
                    file_mid_observer << GlobalStaticVariables::physical_time_ << "," << vel_of_mid_index_observer << "," << vel_of_first_index_observer << "," << vel_of_last_index_observer << "," << conv_checker.get_percentage_difference() << "\n";
                else
                    file_mid_observer << GlobalStaticVariables::physical_time_ << "," << vel_of_mid_index_observer << "," << vel_of_first_index_observer << "," << vel_of_last_index_observer << "," << 0 << "\n";

                file_mid_observer.flush(); // Ensure data is written to the file immediately
                // calculate_observer_avg_vel();
            }
        }
        TickCount t2 = TickCount::now();
        observer_radial_container.update_configureations();
        observer_radial_container.update_observing_quantities();
        observer_radial_container.update_avg_velocity();
        observer_radial_container.writeToFile();
        observer_radial_container.write_velocity_to_csv(in_output.output_folder_);
        axial_velocity_observer_contact.updateConfiguration();
        update_axial_observer_velocity.exec();
        update_axial_observer_pressure.exec();
        observer_axial_recording.writeToFile();
        write_observer_properties_to_output(observer_axial, "Velocity", in_output.output_folder_);
        write_observer_properties_to_output(observer_axial, "Pressure", in_output.output_folder_);
        water_body_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    file_mid_observer.close();
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
    return 0;
}

int main(int argc, char *argv[])
{
    const double Re = 500;
    FDA_nozzle_parameters params;

    // std::vector<int> number_of_particles = {10, 15, 20}; // Create a vector with desired particle counts
    std::vector<int> number_of_particles = {10, 15};

    for (int particles : number_of_particles) // Loop through each number of particles
    {
        params.number_of_particles = particles;
        {
            std::stringstream wall_file;
            size_t np = params.number_of_particles;
            if (np > 20)
                np = 20;
            wall_file << "./input/Fluid_expansion/Wall_expansion_N" << np << ".stl";
            params.wall_file_path = wall_file.str();
        }
        params.end_time = 2.0;

        bool is_TVLCCC = false;
        bool is_wallrelaxation = true;
        bool is_linearVisc = true;
        bool is_refineWall = false;
        bool is_oldoutlet = true;
        // Get current time and format it as a string
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::stringstream datetime_ss;
        datetime_ss << std::put_time(std::localtime(&now_time), "%Y%m%d_%H%M%S");
        std::string datetime = datetime_ss.str();

        std::string additional_output_tag = "_2ndhalfNoRieman_OldOutlet_" + datetime;

        FDA_nozzle(argc, argv, additional_output_tag, is_TVLCCC, is_wallrelaxation, is_linearVisc, is_refineWall, is_oldoutlet, params, Re); // Call the function with current settings
    }

    return 0;
}