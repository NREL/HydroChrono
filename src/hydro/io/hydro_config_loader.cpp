/*********************************************************************
 * @file  hydro_config_loader.cpp
 *
 * @brief Implementation of the simple hydrodynamics config loader.
 *********************************************************************/

#include "hydro_config_loader.h"
#include "../../hydro_yaml_parser.h"

namespace hydrochrono {
namespace hydro {

YAMLHydroData LoadHydroConfigFromYaml(const std::string& yaml_path) {
    return ReadHydroYAML(yaml_path);
}

}  // namespace hydro
}  // namespace hydrochrono

