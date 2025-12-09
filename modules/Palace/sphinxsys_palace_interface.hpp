#pragma once

#include <string>

namespace sphinxsys_palace
{

/**
 * @brief Palace main-style entry point embedded into SPHinXsys.
 *
 * This function is a nearly direct copy of Palace's original main routine:
 *   - It parses command line options (help, version, dry-run).
 *   - It initializes MPI, MFEM device, libCEED, Hypre, and SLEPc (if enabled).
 *   - It constructs the appropriate Palace solver based on the JSON configuration.
 *   - It runs mesh refinement and the full solve loop.
 *   - It prints timing information and writes Palace output files.
 *
 * @param argc  Argument count (same as in a regular main function).
 * @param argv  Argument vector; argv[1] must be the path to a JSON config file.
 * @return int  0 on success, non-zero on error.
 */
int palace_main_embed(int argc, char *argv[]);

/**
 * @brief Convenient wrapper that runs Palace with a single JSON config file.
 *
 * This helper builds a minimal (argc, argv) and forwards the call to
 * palace_main_embed, so that callers only need to provide the JSON file path.
 *
 * @param json_file  Path to the Palace JSON configuration file.
 * @return int       0 on success, non-zero on error.
 */
int RunPalaceFromJson(const std::string &json_file);

}  // namespace sphinxsys_palace
