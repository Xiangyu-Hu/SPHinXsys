/**
 * @file total_artificial_heart.cpp
 * @brief This is the example of total artificial heart implantation path simulation
 * @author John Benjamin, Bence Rochlitz - Virtonomy GmbH
 */
#include "bernoulli_beam.h"

#ifdef __EMSCRIPTEN__

#include <emscripten.h>
#include <emscripten/bind.h>

EMSCRIPTEN_BINDINGS(SPHINXSYS)
{
    emscripten::value_array<std::array<double, 3>>("ArrayDouble3")
        .element(emscripten::index<0>())
        .element(emscripten::index<1>())
        .element(emscripten::index<2>());

    emscripten::value_object<StlData>("StlData")   
        .field("name", &StlData::name)
        .field("ptr", &StlData::ptr);
    
    emscripten::register_vector<std::string>("StringVector");
    emscripten::register_vector<double>("DoubleVector");
    emscripten::register_vector<IndexVector::value_type>("UIntVector");
    emscripten::register_vector<StlData>("StlsList");
    emscripten::register_vector<IndexVector>("BodiesList");
    emscripten::register_map<std::string, std::string>("VtuData");

    emscripten::value_object<SimTotalArtificialHeartInput>("SimTotalArtificialHeartInput")
        .field("material_model_name", &SimTotalArtificialHeartInput::material_model_name)
        .field("scale_stl", &SimTotalArtificialHeartInput::scale_stl)
        .field("resolution", &SimTotalArtificialHeartInput::resolution)
        .field("rho_0", &SimTotalArtificialHeartInput::rho_0)
        .field("poisson", &SimTotalArtificialHeartInput::poisson)
        .field("Youngs_modulus", &SimTotalArtificialHeartInput::Youngs_modulus)
        .field("Youngs_modulus_tah", &SimTotalArtificialHeartInput::Youngs_modulus_tah)
        .field("physical_viscosity", &SimTotalArtificialHeartInput::physical_viscosity)
        .field("translation_tah", &SimTotalArtificialHeartInput::translation_tah)
        .field("stls", &SimTotalArtificialHeartInput::stls)
        .field("relative_input_path", &SimTotalArtificialHeartInput::relative_input_path)
        .field("contacting_bodies_list", &SimTotalArtificialHeartInput::contacting_bodies_list);

    emscripten::class_<SimTotalArtificialHeartJS>("SimTotalArtificialHeart")
        .constructor<SimTotalArtificialHeartInput>()
        .function("runSimulation", &SimTotalArtificialHeartJS::runSimulation)
        .function("onError", &SimTotalArtificialHeartJS::onError)
        .property("vtuData", &SimTotalArtificialHeartJS::getVtuData);
}

#else

int main()
{
    // use simulation collection id to get json file for simulation definittion
    // download JSON file and fill in the following struct

    SimTotalArtificialHeartInput input;
    input.scale_stl = 0.001;
    input.resolution = { 8.0, 8.0, 8.0, 8.0, 8.0, 8.0 };
    // in order
    // resolution_tah = 8.0;
    // resolution_aorta = 8.0;
    // resolution_diaphragm = 8.0;
    // resolution_latrium = 8.0;
    // resolution_partery = 8.0;
    // resolution_ratrium = 8.0;
    input.rho_0 = 6.45e3;
    input.poisson = 0.3;
    input.Youngs_modulus = 5e8;
    input.Youngs_modulus_tah = 1e6;
    input.physical_viscosity = 200.0;
    input.translation_tah = { 0, -200.0, 0 };

    string tah_stl = "TAH_basic2_pos.stl";
    string aorta_stl = "Aorta.stl";
    string diaphragm_stl = "Diaphragm.stl";
    string latrium_stl = "LA.stl";
    string partery_stl = "PA.stl";
    string ratrium_stl = "RA.stl";

    input.stls = { tah_stl, aorta_stl, diaphragm_stl, latrium_stl, partery_stl, ratrium_stl };
    input.relative_input_path = "./input/";
    input.contacting_bodies_list = {
        { 1, 2, 3, 4, 5 },
        { 0, 4 },
        { 0 },
        { 0, 4 },
        { 0, 3, 5 },
        { 0, 4 }
        };

    /* DOWNLOAD STLs files at this point */

    try
    {
        // set up the simulation
        SimTotalArtificialHeartJS simTotalArtificialHeart(input);
        int number_of_steps = 700;
        std::cout << "About to run the simulation" << std::endl;
        simTotalArtificialHeart.runSimulation(number_of_steps);
    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}

#endif