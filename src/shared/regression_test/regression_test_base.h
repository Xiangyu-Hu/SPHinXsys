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
 * @file regression_test_base.h
 * @brief Base classes for comparisons between validated and tested results.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "all_physical_dynamics.h"
#include "io_all.h"
#include "xml_engine.h"

namespace SPH
{
/**
 * @class 	RegressionTestBase
 * @brief 	The base of regression test for various method (time-averaged, ensemble-averaged, dynamic time warping)
 * @details The results of current run is saved in a vector of vector.
 * 			The inner vector gives the values on the observations points. The outer vector gives the snap shots of the observations.
 * 			The results of all runs is saved in a triple vector, where the outermost vector gives the observations of the runs.
 */
template <class ObserveMethodType>
class RegressionTestBase : public ObserveMethodType
{
    /*identify the variable type from the parent class. */
    using VariableType = decltype(ObserveMethodType::type_indicator_);

  protected:
    std::string input_folder_path_;      /*< the folder path for the input folder. (folder) */
    std::string in_output_filefullpath_; /*< the file path for current result from xml memory to xml file. */
    std::string result_filefullpath_;    /*< the file path for all run results. (.xml)*/
    std::string runtimes_filefullpath_;  /*< the file path for run times information. (.dat)*/
    std::string converged;               /*< the tag for result converged, default false. */

    XmlMemoryIO xmlmemory_io_; /*< xml memory in_output operator, which has defined several
                                                                              methods to read and write data from and into xml memory,
                                                                                  including one by one, or all result in the same time. */

    XmlEngine observe_xml_engine_;    /*< xml engine for current result io_environment. */
    XmlEngine result_xml_engine_in_;  /*< xml engine for input result. */
    XmlEngine result_xml_engine_out_; /*< xml engine for output result. */
                                      /*< the XmlEngine can operate the node name and elements in xml memory. */

    StdVec<std::string> element_tag_;             /*< the container of the tag of current result. */
    BiVector<VariableType> current_result_;       /*< the container of current run result stored as snapshot * observation. */
    BiVector<VariableType> current_result_trans_; /*< the container of current run result with snapshot & observations transposed,
                                                                                               because this data structure is required in TA and DTW method. */
    BiVector<VariableType> result_in_;            /*< the temporary container of each result when reading from .xml memory, and
                                                                                                  it may be snapshot * observations (reading all result for averaged methods),
                                                                                                      or observations * snapshot (reading specified result for TA and DTW method.) */
    TriVector<VariableType> result_;              /*< the container of results in all runs (run * snapshot * observation) */

    int snapshot_, observation_; /*< the size of each layer of current result vector. */
    int difference_;             /*< the length difference of snapshot between old and new result,
                                                                                which is used to trim each new run result to be a unified length. */
    int number_of_run_;          /*< the times of run. */
    int label_for_repeat_;       /*< the label used stable convergence (several convergence). */
    int number_of_snapshot_old_; /*< the snapshot size of last trimmed result. */

  public:
    template <typename... ConstructorArgs>
    explicit RegressionTestBase(ConstructorArgs &&...args) : ObserveMethodType(std::forward<ConstructorArgs>(args)...), xmlmemory_io_(),
                                                             observe_xml_engine_("xml_observe_reduce", this->quantity_name_),
                                                             result_xml_engine_in_("result_xml_engine_in", "result"),
                                                             result_xml_engine_out_("result_xml_engine_out", "result")
    {
        input_folder_path_ = this->io_environment_.input_folder_;
        in_output_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + ".xml";
        result_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + "_result.xml";
        runtimes_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + "_runtimes.dat";

        if (!fs::exists(runtimes_filefullpath_))
        {
            converged = "false";
            number_of_run_ = 1;
            label_for_repeat_ = 0;
        }
        else
        {
            std::ifstream in_file(runtimes_filefullpath_.c_str());
            in_file >> converged;
            in_file >> number_of_run_;
            in_file >> label_for_repeat_;
            in_file.close();
        };
    };
    virtual ~RegressionTestBase();

    void writeToXml(ObservedQuantityRecording<VariableType> *observe_method, size_t iteration = 0);
    template <typename ReduceType>
    void writeToXml(ReducedQuantityRecording<ReduceType> *reduce_method, size_t iteration = 0);
    /* read current result from xml file into xml memory. */
    void readFromXml(ObservedQuantityRecording<VariableType> *observe_method);
    template <typename ReduceType>
    void readFromXml(ReducedQuantityRecording<ReduceType> *reduce_method);

    void transposeTheIndex();                  /** transpose the current result (from snapshot*observation to observation*snapshot). */
    void readResultFromXml();                  /** read the result from the .xml file. (all result) */
    void writeResultToXml();                   /** write the result to the .xml file. (all result) */
    void readResultFromXml(int index_of_run_); /* read the result from the .xml file with the specified index. (DTW method, TA method) */
    void writeResultToXml(int index_of_run_);  /* write the result to the .xml file with the specified index. (DTW method, TA method) */

    /** the interface to write observed quantity into xml memory. */
    void writeToFile(size_t iteration = 0) override
    {
        ObserveMethodType::writeToFile(iteration); /* used for visualization (.dat)*/
        writeToXml(this, iteration);               /* used for regression test. (.xml) */
    };

    /*The new run result is stored in the xml memory, and can't be copied and used directly
     *so they should be output as .xml file and reloaded from .xml file into the memory.*/
    /** the interface to write xml memory into Xml file. */
    void writeXmlToXmlFile()
    {
        observe_xml_engine_.writeToXmlFile(in_output_filefullpath_); /* two method are same.*/
    };

    /** read current result from Xml file into Xml memory. */
    void readXmlFromXmlFile()
    {
        readFromXml(this);
    };
};
}; // namespace SPH