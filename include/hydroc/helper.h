#pragma once

#include <string>

/**@brief Base namespace for HydroChrono library
 * 
*/
namespace hydroc {

/**@brief Set initial environment
 * 
 * Use command line argument or env variables to set s
 * ome static data at initialization.
 * 
 * Set the main data directory ..
 * 
 * @param argc number of argument (same as for main function)
 * @param argv arguments of main function 
 * @return 1 on error 0 else
*/
int setInitialEnvironment(int argc, char* argv[]) noexcept;

/**@brief Get base name of data directory
 * 
 * @return the string containing the path in standard format
*/
std::string getDataDir() noexcept;


} // end namespace hydroc