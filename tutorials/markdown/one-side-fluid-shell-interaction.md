# Algorithm for one-sided FSI problem
The basic idea is to retrieve kernel completeness by using dummy shell contact particles. The algorithm is summarized below:

1. Create a new class `NeighborBuilderOneSidedSurfaceContact`.
2. Protected variable: `StdLargeVec<Vecd>& n_`, `StdLargeVec<double>& k_`, `double particle_distance_`, which will be obtained from `contact_body` in the initializer. The first is the normal direction, the second is the curvature and the last is the reference resolution of shell particles.
3. override `initializeNeighbor` and `createNeighbor`
```
Vecd get_dWij_Vj_e_ij_ttl(const Real &distance,
                        const Vecd &displacement, const Real &Vol_j, const Vecd &pos_i, const Vecd &pos_j, const Vecd &n_j, const Real &radius_j)
{
    // calculate summation of dWij_Vj*e_ij
    Vecd dWij_Vj_e_ij_ttl = kernel_->dW(distance, displacement) * Vol_j * displacement / (distance + TinyReal); 
    
    while (distance < kernel_->CutOffRadius())
    {
        // calculate the position and volume of the next dummy particle
        pos_j += n_j * particle_distance_;
        displacement = pos_i - pos_j;
        distance = displacement.norm();
        
        // for 3D, the volume (surface) can be approximated by S2/S1 = (R+dp)^2/R^2
        // for 2D, the volume (diameter) can be approximated by D2/D1 = (R+dp)/R
        double Vol_j *= std::pow(1+particle_distance_/radius_j, Dimensions);

        // calculate dWij and eij
        double dW_ij = kernel_->dW(distance, displacement);
        Vecd e_ij = displacement / (distance + TinyReal);

        // go to the next loop, repeat until distance is larger than cut off radius
        dWij_Vj_e_ij_ttl += dW_ij * e_ij * Vol_j;
    }
    
    return dWij_Vj_e_ij_ttl;
}
```

```
void NeighborBuilderOneSidedSurfaceContact::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                     const Vecd &displacement, size_t index_j, const Real &Vol_j, const Vecd &pos_i, const Vecd &n_j, const Real &radius_j)
{
    neighborhood.j_.push_back(index_j);
    neighborhood.W_ij_.push_back(kernel_->W(distance, displacement));
    neighborhood.r_ij_.push_back(distance);

    // calculate summation of dWij_Vj*e_ij
    Vecd dWij_Vj_e_ij_ttl = get_dWij_Vj_e_ij_ttl(distance, displacement, Vol_j, pos_i, pos_j, n_j, radius_j); 

    dWij_Vj_corrected = dWij_Vj_e_ij_ttl.norm();
    e_ij_corrected = dWij_Vj_e_ij_ttl / dWij_Vj_corrected;
    
    neighborhood.dW_ijV_j_.push_back(dWij_Vj_corrected);
    neighborhood.e_ij_.push_back(e_ij_corrected);

    neighborhood.allocated_size_++;
}
```

```
void NeighborBuilderOneSidedSurfaceContact::createNeighbor(Neighborhood &neighborhood, const Real &distance,
                                     const Vecd &displacement, size_t index_j, const Real &Vol_j, const Vecd &pos_i, const Vecd &n_j, const Real &radius_j)
{
    size_t current_size = neighborhood.current_size_;
    neighborhood.j_[current_size] = index_j;
    neighborhood.W_ij_[current_size] = kernel_->W(distance, displacement);
    neighborhood.r_ij_[current_size] = distance;
    
    // calculate summation of dWij_Vj*e_ij
    Vecd dWij_Vj_e_ij_ttl = get_dWij_Vj_e_ij_ttl(distance, displacement, Vol_j, pos_i, pos_j, n_j, radius_j); 

    dWij_Vj_corrected = dWij_Vj_e_ij_ttl.norm();
    e_ij_corrected = dWij_Vj_e_ij_ttl / dWij_Vj_corrected;
    neighborhood.dW_ijV_j_[current_size] = kernel_->dW(distance, displacement) * Vol_j;
    neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
}
```
3. Pseudo code of operator():

```
void NeighborBuilderOneSidedSurfaceContact::operator()(Neighborhood &neighborhood,
                                        const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = std::get<0>(list_data_j);
    Vecd displacement = pos_i - std::get<1>(list_data_j);
    Real distance = displacement.norm();

    Vecd n_j = n_[index_j]; // normal direction of shell particle in cell linked list with index j
    double radius_j = 1 / k_[index_j]; // curvature at shell particle j

    // correct normal direction, make sure it points from fluid to shell
    Vecd n_j_corrected = sgn(displacement.dot(n_j)) * n_j; // sign of r_ij and n_j

    if (distance < kernel_->CutOffRadius())
    {
        neighborhood.current_size_ >= neighborhood.allocated_size_
            ? createNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j), pos_i, n_j_corrected, radius_j)
            : initializeNeighbor(neighborhood, distance, displacement, index_j, std::get<2>(list_data_j));
        neighborhood.current_size_++;
    }
}
```
4. Write a new class `OneSidedSurfaceContactRelation` with protected variable `UniquePtrsKeeper<NeighborBuilderOneSidedSurfaceContact> neighbor_builder_contact_ptrs_keeper_`

5. write auxiliary classes to calculate shell curvature


# Shell curvature
We can compute H and K by:

$$H=\nabla \cdot \mathbf{n}_r$$

$$K=\frac{1}{2}(H^2-\sum_i \sum_j (\frac{\partial n_j}{\partial x_i})^2)$$

where $\mathbf{n}_r$ is the normal direction of the surface. 


# Reference
- Nitschke, Ingo, Axel Voigt, and Jörg Wensch. "A finite element approach to incompressible two-phase flow on manifolds." Journal of Fluid Mechanics 708 (2012): 418-438.