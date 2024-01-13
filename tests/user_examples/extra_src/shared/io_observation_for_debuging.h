/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	io_observation.h
 * @brief 	Classes for input and output with vtk (Paraview) files.
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#pragma once

#include "io_base.h"
#include "base_fluid_dynamics.h"
#include "io_plt.h"

namespace SPH
{
    //class FluidDataSimple;

     typedef DataDelegateSimple<BaseParticles> BaseDataSimple;
     
     template <class ReturnType, class OperationType>
     class ReducedQuantityRecordingForDebuging : public BaseIO, public BaseDataSimple
    {
       protected:
         PltEngine plt_engine_;
         std::string sph_body_name_;
         const std::string quantity_name_;
         std::string filefullpath_output_;
         StdLargeVec<ReturnType> &varriable_for_output_;
         OperationType operation_;
         ReturnType reference_;
         OperationType &getOperation() { return operation_; };

       public:
         ReducedQuantityRecordingForDebuging(IOEnvironment &io_environment, SPHBody &sph_body, ReturnType reference, std::string quantity_name)
             : BaseIO(io_environment), BaseDataSimple(sph_body), plt_engine_(),
               sph_body_name_(sph_body.getName()), varriable_for_output_(*particles_->getVariableByName<ReturnType>(quantity_name)),
               quantity_name_(quantity_name), reference_(reference)
         {
             /** output for .dat file. */
             filefullpath_output_ = io_environment_.output_folder_ + "/" + sph_body_name_ + "_" + quantity_name_ + ".dat";
             std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
             out_file << "\"run_time\""
                      << "   ";
             plt_engine_.writeAQuantityHeader(out_file, reference_, quantity_name_);
             out_file << "\n";
             out_file.close();
         };
         virtual ~ReducedQuantityRecordingForDebuging(){};

         virtual void writeToFile(size_t iteration_step = 0) override
         {
             std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
             out_file << GlobalStaticVariables::physical_time_ << "   ";
             plt_engine_.writeAQuantity(out_file, this->particle_reduce(reference_, this->getOperation(), varriable_for_output_));
             out_file << "\n";
             out_file.close();
         };

         template <class ReturnType, typename OperationType>
         inline ReturnType particle_reduce(ReturnType temp, OperationType &&operation,
                                           const StdLargeVec<ReturnType> &varriable_for_output)
         {
             for (size_t i = 0; i < varriable_for_output.size(); ++i)
             {
                 temp = operation(temp, varriable_for_output[i]);
             }
             return temp;
         }
     };

      template <class ReturnType>
     class GlobalQuantityRecordingForDebuging : public BaseIO, public BaseDataSimple
    {
       protected:
         std::string sph_body_name_;
         const std::string quantity_name_;
         StdLargeVec<ReturnType> &varriable_for_output_;
         ReturnType reference_;
         StdLargeVec<Vecd>& varriable_position_;

       public:
         GlobalQuantityRecordingForDebuging(IOEnvironment &io_environment, SPHBody &sph_body, ReturnType reference, std::string quantity_name)
             : BaseIO(io_environment), BaseDataSimple(sph_body),
               sph_body_name_(sph_body.getName()), varriable_for_output_(*particles_->getVariableByName<ReturnType>(quantity_name)),
               quantity_name_(quantity_name), reference_(reference), varriable_position_(*particles_->getVariableByName<Vecd>("Position"))
         {
            
         };
         virtual ~GlobalQuantityRecordingForDebuging(){};

         void writeToFile( size_t iteration_step = static_cast<size_t>(GlobalStaticVariables::physical_time_ * 100000) ) override
         {
              /** output for .dat file. */
             std::string filefullpath_output = io_environment_.output_folder_ + "/" + sph_body_name_ + "_" + quantity_name_ + "_"+ std::to_string(iteration_step) + ".dat";
             std::ofstream out_file(filefullpath_output.c_str(), std::ios::app);
             this->writeAQuantityWithPosition(out_file, reference_);
             out_file << "\n";
             out_file.close();
         };

       protected:    
        void writeAQuantityWithPosition(std::ofstream &out_file, Real reference)
        {
            for (size_t n = 0; n != varriable_for_output_.size(); ++n)
            {
                for (int i = 0; i < Dimensions; ++i)
                {
                    out_file << std::fixed << std::setprecision(9) << varriable_position_[n][i] << "   ";
                }
                out_file << n << "   ";
                out_file << std::fixed << std::setprecision(9) << varriable_for_output_[n]<< "   ";  
                out_file << "\n";
            } 
        }

        void writeAQuantityWithPosition(std::ofstream &out_file, Vecd reference)
        {
            for (size_t n = 0; n != varriable_for_output_.size(); ++n)
            {
                for (int i = 0; i < Dimensions; ++i)
                {
                    out_file << std::fixed << std::setprecision(9) << varriable_position_[n][i] << "   ";
                }
                out_file << n << "   ";
                for (int j = 0; j < Dimensions; ++j)
                {
                    out_file << std::fixed << std::setprecision(9) << varriable_for_output_[n][j] << "   ";
                }
                out_file << "\n";
            } 
        }
     };

    
} // namespace SPH
