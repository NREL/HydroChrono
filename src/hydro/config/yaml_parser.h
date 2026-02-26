/*********************************************************************
 * @file  yaml_parser.h
 *
 * @brief YAML parser for hydro.yaml files.
 *
 * This is an internal header. Library consumers should use the public
 * LoadHydroConfigFromYaml() function from config_loader.h instead.
 *********************************************************************/

#ifndef HYDRO_CONFIG_YAML_PARSER_H
#define HYDRO_CONFIG_YAML_PARSER_H

#include <hydroc/config/hydro_config.h>
#include <string>

/**
 * @brief Read and parse a hydro.yaml file into YAMLHydroData structure.
 *
 * @param hydro_file_path Path to the hydro.yaml file to parse.
 * @return YAMLHydroData containing the parsed configuration.
 * @throws std::runtime_error if the file cannot be read or parsed.
 */
YAMLHydroData ReadHydroYAML(const std::string& hydro_file_path);

#endif  // HYDRO_CONFIG_YAML_PARSER_H

