# Algorithm for one-sided FSI problem
The basic idea is to retrieve kernel completeness by using dummy shell contact particles. The algorithm is summarized below:

1. Create a new class `NeighborBuilderOneSidedSurfaceContact`.
2. Protected variable: `StdLargeVec<Vecd>& n_`, `StdLargeVec<double>& k_`, `double particle_distance_`, which will be obtained from `contact_body` in the initializer. The first is the normal direction, the second is the curvature and the last is the reference resolution of shell particles.
3. Pseudo code of operator():

```
    size_t index_j = std::get<0>(list_data_j);
    Vecd n_j = n_[index_j]; // normal direction of shell particle in cell linked list with index j
    double radius_j = 1 / k_[index_j]; // curvature at shell particle j
    
    Vecd pos_j = std::get<1>(list_data_j);
    Vecd displacement = pos_i - pos_j;
    Real distance = displacement.norm();

    Real Vol_j = std::get<2>(list_data_j);

    // correct normal direction, make sure it points from fluid to shell
    Vecd n_j_corrected = sgn(displacement.dot(n_j)) * n_j; // sign of r_ij and n_j

    while (distance < kernel_->CutOffRadius())
    {
        // create new neighborhood
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j, Vol_j)
            : initializeNeighbor(neighborhood, distance, displacement, index_j, Vol_j);
        neighborhood.current_size_++; 
        
        // calculate the position and volume of the next dummy particle
        pos_j += n_j_corrected * particle_distance_;
        displacement = pos_i - pos_j;
        distance = displacement.norm();

        // index will not change, so that the other variables (e.g. vel_ave_, acc_ave_) of the dummy particles will be the same as index_j 
        
        // for 3D, the volume (surface) can be approximated by S2/S1 = (R+dp)^2/R^2
        // for 2D, the volume (diameter) can be approximated by D2/D1 = (R+dp)/R
        double Vol_j *= std::pow(1+particle_distance_/radius_j, Dimensions);

        // go to the next loop, repeat until distance is larger than cut off radius
    }
```
4. Write a new class `OneSidedSurfaceContactRelation` with protected variable `UniquePtrsKeeper<NeighborBuilderOneSidedSurfaceContact> neighbor_builder_contact_ptrs_keeper_`

5. write auxiliary classes to calculate shell curvature