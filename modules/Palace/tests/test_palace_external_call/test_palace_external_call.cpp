// modules/Palace/tests/test_palace_external_call/test_palace_external_call.cpp

#include <cstdlib>
#include <iostream>
#include <string>

int main()
{
    std::cout << "[SPHinXsys-Palace] External Palace call test\n";


    // PALACE_EXECUTABLE is provided by CMake (via compile definitions).
    // It points to the *actual palace binary* (e.g., palace/build/bin/palace).
    const char *palace_exe = PALACE_EXECUTABLE;


    // Hardcoded relative path; copied into the build directory by CMake.
    std::string config_file = "data/rings.json";


    // Build shell command: "<palace_executable> data/rings.json"
    std::string cmd = PALACE_EXECUTABLE;  // palace_exe
    cmd += " ";
    cmd += config_file;

    std::cout << "Running command:\n  " << cmd << std::endl;

    // Execute external process.
    // This tests the case where SPHinXsys calls Palace "as a subprocess",
    // rather than linking via the static library interface.
    int ret = std::system(cmd.c_str());
    if (ret != 0)
    {
        std::cerr << "Palace process failed with code " << ret << std::endl;
        return ret;
    }

    std::cout << "Palace external run finished successfully.\n";
    return 0;
}
