/*********************************************************************
 * @file  moordyn_config.h
 * @brief Configuration struct for MoorDyn mooring coupling.
 *********************************************************************/

#ifndef HYDROC_CONFIG_MOORDYN_CONFIG_H
#define HYDROC_CONFIG_MOORDYN_CONFIG_H

#include <string>
#include <vector>

struct MoorDynConfig {
    bool enabled = false;
    std::string input_file;
    /// 0-based indices into the HydroChrono body list that correspond to
    /// MoorDyn coupled bodies, in the order they appear in the MoorDyn
    /// input file's Body List.
    std::vector<int> coupled_body_indices;
};

#endif  // HYDROC_CONFIG_MOORDYN_CONFIG_H
