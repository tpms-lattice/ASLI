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

#ifndef INFILL_H
#define INFILL_H

#include "TrilinearInterpolation.h"
#include "ExceptionClass.h" // Custom exception class

/* Standard library headers */
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <list>
#include <unordered_map>

struct latticeType {
	std::string type;                     // Unit cell type
	std::string side        = "scaffold"; // Unit cell side
	std::string filterType  = "gaussian"; // filter (if type is hybrid)
	double filterRadius     = 1.0;        // filter radius (if type is hybrid)
	double correctionFactor = 0.25;       // filter correction (if type is hybrid)

	std::vector<modelData> interpModel_linear; // Interpolation model
	std::vector<std::string> typeVector;       // Unit cell types described by the
	                                           // corresponding interpolation model
};

struct latticeSize {
	double size;             // Unit cell size
	double scaling;          // Unit cell scaling
	double minUnitCellSize;  // Smallest unit cell size in design
	double maxUnitCellSize;  // Largest unit cell size in design
	double meanUnitCellSize; // Mean unit cell size

	std::vector<std::vector<double>> data; // Unit cell scaling data
	modelData interpModel_linear;          // Interpolation model
};

struct latticeUDF { 
	std::string userDefinedFeature; // User defined feature label
	double A;                       // User defined constant 1
	double B;                       // User defined constant 2
	double C;                       // User defined constant 3
	double D;                       // User defined constant 4
	double E;                       // User defined constant 5
};
struct latticeFeature {
	std::string feature;           // Unit cell feature label
	double feature_val;            // Unit cell feature value
	std::string mode = "relative"; // Type of sizing specification
	latticeUDF udf;                // User defined feature parameters

	std::vector<std::vector<double>> data; // Unit cell feature data
	modelData interpModel_linear;          // Interpolation model
};

struct featureSize {
	double wallSize;
	double poreSize;
};


template<typename CoordinateType>
class ASLI_point {
public:
	ASLI_point() : x_val(0), y_val(0), z_val(0) {};

	// Store data
	ASLI_point(CoordinateType const &x, CoordinateType const &y, 
	           CoordinateType const &z) : x_val(x), y_val(y), z_val(z) {};

	// Retrive data
	CoordinateType const & x() const { return x_val; };
	CoordinateType const & y() const { return y_val; };
	CoordinateType const & z() const { return z_val; };

private:
	CoordinateType x_val, y_val, z_val;
};


typedef ASLI_point<double> Point;

namespace Infill {
	// Constants
	const double PI = 4.0*std::atan(1.0);
	
	const std::unordered_map<std::string, std::pair<double, double>> TPMS_av = {
		{"gyroid",          {0.05, 0.95}},
		{"sheet_gyroid",    {0.05, 0.95}},
		{"diamond",         {0.10, 0.90}},
		{"sheet_diamond",   {0.05, 0.80}},
		{"primitive",       {0.22, 0.78}},
		{"sheet_primitive", {0.05, 0.55}},
		{"IWP",             {0.05, 0.90}},
		{"sheet_IWP",       {0.05, 0.85}},
		{"cubic",           {0.05, 0.90}}
	};

	const std::list<std::string> feature_av = {
		{"volumeFraction"},
		{"isovalue"},
		{"wallSize"},
		{"poreSize"},
		{"userDefined"}
	};

	// Signed distance functions
	double TPMS_function(const Point &p, const std::string &type,
		const double &scaling, const double &t);

	double TPMS_function(const Point &p, const latticeType &lt_type,
		const latticeSize &lt_size, const latticeFeature &lt_feature);

	// Mesh sizing functions
	double sizing_function(const Point &p, const latticeSize &lt_size,
		const std::string &mode);

	featureSize featureSize_function(const Point &p, const latticeType &lt_type, 
		const latticeSize &lt_size, const latticeFeature &lt_feature);

	// Parameter conversions
	namespace internal {
		double input2level(const std::string &type, const double &scaling,
			const std::string &feature, const double &featureValue, 
			const latticeUDF &userDefinedParameters, const std::string &featureMode);

		// Normalizations
		double unnormalizeLevel(const double &t_normalized, const std::string &type);

		// Convertions to isovalue (level-set constant)
		double vFraction2level(double vFraction, const std::string &type);

		double wallSize2level(double wallSize, const double &scaling,
			const std::string &type, const std::string &mode);
		double poreSize2level(double poreSize, const double &scaling,
			const std::string &type, const std::string &mode);

		// Conversions from isovalue (level-set constant)
		double level2wallSize(double t, const double &scaling,
			const std::string &type);
		double level2poreSize(double t, const double &scaling,
			const std::string &type);

		// Others
		double userDefinedInput2vFraction(const double &userDefinedInput, 
			const latticeUDF &userDefinedParameters);
	}
};

// Error messages
namespace INFILL_ERRMSG {
	const std::string INVALID_TPMS = "Invalid TPMS";
	const std::string INVALID_FEATURE_REQUEST = "Invalid feature requested";
	const std::string INVALID_USER_DEFINED_FEATURE_REQUEST = "Invalid user defined feature requested";
	const std::string NO_STRESS_CONVERSION_DEFINED = "Stress conversion has not been defined";
}

#endif