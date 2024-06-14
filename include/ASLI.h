/* ==========================================================================
 *  This file is part of ASLI (A Simple Lattice Infiller)
 *  Copyright (C) KU Leuven, 2019-2024
 *
 *  ASLI is free software: you can redistribute it and/or modify it under the 
 *  terms of the GNU Affero General Public License as published by the Free 
 *  Software Foundation, either version 3 of the License, or (at your option) 
 *  any later version.
 *
 *  ASLI is distributed in the hope that it will be useful, but WITHOUT ANY 
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 *  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for 
 *  more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *  Please read the terms carefully and use this copy of ASLI only if you
 *  accept them.
 * ==========================================================================*/

#ifndef ASLI_H
#define ASLI_H

#include "BasicGeometries.h"
#include "Infill.h"
#include "Mesh.h"
#include "TrilinearInterpolation.h"
#include "Filter.h"
#include "IO_ASLI.h"

#include "ExceptionClass.h" // Custom exception class

/* alglib headers */
#include "stdafx.h"
#include "alglibmisc.h"

/* yaml-cpp headers */
#include "yaml-cpp/yaml.h"

/* Standard library headers */
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <string>//for std::stod
#include <set>

class ASLI {
public:
	/* Constructors|Destructors */
	ASLI(std::string configFile);
	~ASLI();

	/* Parameters */
	// Constants
	const double PI = 4.0*std::atan(1.0);

	// Outer shell
	polygonSoup shell;

	// Lattice
	latticeType lt_type;
	latticeTypeTree lt_typeTree;
	latticeSize lt_size;
	latticeFeature lt_feature;

	// Mesh
	meshSettings me_settings;

private:
	const std::set<std::string> MESHERS = {"CGAL", "MMG", "CGAL_OLD"};
	const std::set<std::string> SIDES = {"scaffold", "void", "all"};
	const std::set<std::string> FEATURES = {"wallSize", "poreSize"};
	const std::set<std::string> MODES = {"relative", "absolute"};

	void SetUp(std::string configFile);
	void SetUpLattice();
	void SetUpKdtree(const alglib::real_2d_array &alglib_data, 
	                 const alglib::integer_1d_array &tags,
	                 alglib::kdtree &kdt);
	void SetUpTypeInterpolator(const alglib::real_2d_array &coordinates,
	                           const std::vector<double> &weights, 
	                           modelData &interpolationModel);

	double checkFeatureValueValidity(const latticeType &lt_type,
	                                 const latticeSize &lt_size,
	                                 const latticeFeature &lt_feature);
	double checkFeatureFieldValidity(const Point &p, const double &featureValue,
	                                 const latticeType &lt_type,
	                                 const latticeSize &lt_size,
	                                 const latticeFeature &lt_feature);

	std::vector<std::string> split(const std::string &s, 
	                               const std::string &delimiter);

	void LoadInputFiles();

	/* Parameters */
	std::string tapFile, sapFile, fapFile;
};

// Error messages
namespace ASLI_ERRMSG {
	const std::string INVALID_NUMBER_OF_COMMAND_LINE_ARGUMENTS = "ERROR: Invalid number of command line arguments.";
	const std::string CONFIGURATION_FILE_WITHOUT_EXTENSION = std::string("YAML configuration file not specified or ") + 
		"specified without extension. Please provide the path to the configuration file including filename and extension.";

	const std::string UNDEFINED_NODE = "Undefined node";
	const std::string UNDEFINED_KEY = " not defined in configuration file";

	// Input files
	const std::string FAILED_TO_OPEN_CONFIG_FILE = "Unable to open configuration file";
	const std::string FAILED_TO_PARSE_CONFIG_FILE = "Failed to parse configuration file";

	const std::string FAILED_TO_OPEN_STL_FILE = "Unable to open .stl file";
	const std::string INTERNAL_GEOMETRY_DIMENSION_ZERO = "Dimension(s) of internal geometry must be larger than zero.";
	const std::string FAILED_TO_OPEN_TAP_FILE = "Unable to open .tap file";
	const std::string FAILED_TO_OPEN_SAP_FILE = "Unable to open .sap file";
	const std::string FAILED_TO_OPEN_FAP_FILE = "Unable to open .fap file";

	const std::string NO_INPUT_GEOMETRY_SPECIFIED = "Incorrect or missing input geometry information";
	const std::string FAILED_TO_READ_TAP_KEY = "Unable to read .tap file key";
	const std::string FAILED_TO_READ_SAP_KEY = "Unable to read .sap file key";
	const std::string FAILED_TO_READ_FAP_KEY = "Unable to read .fap file key";

	const std::string INPUT_GEOMETRY_IS_NOT_CLOSED = "Input geometry does not form a closed surface";

	const std::string INVALID_TYPE_FIELD = " invalid unit cell type(s) found in .tap file";
	const std::string INVALID_SIZE_FIELD = " invalid unit cell size entry/entries in .sap file";
	const std::string INVALID_FEATURE_FIELD = " invalid feature value entry/entries in .fap file";

	// Lattice parameters
	const std::string INVALID_TPMS = "Invalid lattice type specified";
	const std::string INVALID_SIDE = "Invalid scaffold side specified";
	const std::string INVALID_FILER_RADIUS = "Hybrid mode filter radius cannot be negative";
	const std::string INVALID_CORRECTION_FACTOR = "Hybrid mode correction factor cannot be negative";

	const std::string INVALID_UNIT_CELL_SIZE = "Invalid unit cell size specified";

	const std::string INVALID_FEATURE = "Invalid unit cell feature specified";
	const std::string INVALID_FEATURE_VALUE = "Unit cell feature value cannot be negative";
	const std::string INVALID_FEATURE_MODE = "Invalid unit cell feature mode specified";

	// Mesh settings
	const std::string INVALID_MESHER = "Invalid mesher specified";
	const std::string INVALID_THREAD_NUMBER = "Number of threads cannot be an integer smaller than one";
	const std::string INVALID_MESH_SIZE = "Mesh size must be larger than zero";
	const std::string INVALID_EDGE_PROTECTION_ANGLE = "Edge protection angle cannot be negative";

	const std::string INVALID_FACET_DISTANCE = "Facet distance must be positive";
	const std::string INVALID_CELL_RADIUS_EDGE_RATIO = "Cell radius edge ratio must be larger than two";
	const std::string INVALID_FACET_ANGLE = "Facet angle must be larger than zero and not larger than thirty degrees";
	const std::string INVALID_RELATIVE_ERROR_BOUND = "Relative error bound must be positive";
	const std::string INVALID_FACET_SIZE = "Facet size must be larger than zero";
	const std::string INVALID_EDGE_SIZE = "Edge size must be larger than zero";
	const std::string INVALID_OFFSET = "Poisson offset must be larger than 0.1";

	const std::string INVALID_HAUSD = "hausd parameter cannot be negative";
	const std::string INVALID_HGRAD = "hgrad must be larger than 1";
	const std::string INVALID_HINITIAL = "hinitial cannot be negative";
	const std::string INVALID_MAX_MEMORY = "Max. memory cannot be negative";

	const std::string INVALID_THRESHOLD = "Invalid threshold value";
	
}

template <typename T>
std::optional<T> getOptional(const YAML::Node& node, const std::string& key, const bool required, std::set<std::string>& accessed_keys) {
	try {
		if (node.IsDefined()) {
			if (node[key] && node[key].IsDefined()) {
				std::optional<T> data = node[key].as<T>();
				accessed_keys.insert(key);
				return data;
			} else if (required) {
				throw std::runtime_error("\"" + key + "\"" + ASLI_ERRMSG::UNDEFINED_KEY);
			}
		} else if (required) {
			throw std::runtime_error(ASLI_ERRMSG::UNDEFINED_NODE);
		}
	} catch (const std::exception& e) {
		throw std::runtime_error("Failed to parse \"" + key + "\" (" + e.what() + ")");
	}
	
	return std::nullopt;
}

#endif