/*********************************************************************
 * @file  hydro_config_loader.h
 *
 * @brief Simple helper to load hydrodynamics config from YAML.
 *
 * This wrapper keeps the handwritten parser in one place and gives
 * callers a stable, easy-to-read entry point.
 *********************************************************************/

#ifndef HYDRO_IO_HYDRO_CONFIG_LOADER_H
#define HYDRO_IO_HYDRO_CONFIG_LOADER_H

#include "../../hydro_types.h"
#include <string>

namespace hydrochrono {
namespace hydro {

/**
 * @brief Load hydrodynamics configuration from a YAML file.
 *
 * This function keeps the existing parsing behaviour exactly the same.
 * It simply forwards to the legacy ReadHydroYAML helper.
 *
 * @param yaml_path Absolute or relative path to hydro.yaml.
 * @return Parsed YAMLHydroData structure.
 * @throws std::runtime_error on parsing errors (same as before).
 */
YAMLHydroData LoadHydroConfigFromYaml(const std::string& yaml_path);

}  // namespace hydro
}  // namespace hydrochrono

#endif  // HYDRO_IO_HYDRO_CONFIG_LOADER_H

