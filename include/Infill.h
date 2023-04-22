/* ==========================================================================
 *  This file is part of ASLI (A Simple Lattice Infiller)
 *  Copyright (C) KU Leuven, 2019-2023
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

/* Standard library headers */
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>

struct latticeType {
	std::string type;        // Unit cell type
	std::string filterType;  // filter (if type is hybrid)
	double filterRadius;     // filter radius (if type is hybrid)
	double correctionFactor; // filter correction (if type is hybrid)

	std::vector<modelData> interpModel_linear;
	std::vector<std::string> typeVector; // Unit cell type described by the
	                                     // corresponding interpolation model
};

struct latticeSize {
	double size;             // Unit cell size
	double scaling;          // Unit cell scaling
	double minUnitCellSize;  // Smallest unit cell size in design
	double maxUnitCellSize;  // Largest unit cell size in design
	double meanUnitCellSize; // Mean unit cell size

	std::vector<std::vector<double>> data; // Unit cell scaling data
	modelData interpModel_linear;
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
	std::string feature;       // Unit cell feature label
	double feature_val;        // Unit cell feature value
	std::string mode;          // Type of sizing specification
	latticeUDF udf;            // User defined feature parameters

	std::vector<std::vector<double>> data; // Unit cell feature data
	modelData interpModel_linear;
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

//	void setPoints(CoordinateType const &x, CoordinateType const &y, 
//	               CoordinateType const &z) { x_val = x; y_val = y; z_val = z; };//Not working

//	void x(CoordinateType const &x) { x_val = x; };//Not working
//	void y(CoordinateType const &y) { y_val = y; };//Not working
//	void z(CoordinateType const &z) { z_val = z; };//Not working

	// Retrive data
	CoordinateType const & x() { return x_val; };
	CoordinateType const & y() { return y_val; };
	CoordinateType const & z() { return z_val; };

private:
	CoordinateType x_val, y_val, z_val;
};


typedef ASLI_point<double> Point;

namespace Infill {
	// Constants
	const double PI = 4.0*std::atan(1.0);

	// Signed distance functions
	double TPMS_function(Point p, latticeType *lt_type, latticeSize *lt_size, 
	                     latticeFeature *lt_feature);
	double TPMS_function(Point p, std::string type, double scaling, double t);

	// Mesh sizing functions
	double sizing_function(Point p, latticeSize *lt_size, std::string mode);
	featureSize featureSize_function(Point p, latticeType *lt_type, 
	                           latticeSize *lt_size, latticeFeature *lt_feature);

	// Parameter conversions
	namespace internal {
		double input2level(std::string type, double scaling, std::string feature, 
		                  double featureValue, latticeUDF userDefinedParameters,
		                  std::string featureMode);

		// Normalizations
		double unnormalizeLevel(double t_normalized, std::string type);

		// Convertions to isovalue (level-set constant)
		double vFraction2level(double vFraction, std::string type);
		double wallSize2level(double wallSize, double scaling, std::string type, 
		                      std::string mode);
		double poreSize2level(double poreSize, double scaling, std::string type, 
		                      std::string mode);

		// Conversions from isovalue (level-set constant)
		double level2wallSize(double t, double scaling, std::string type);
		double level2poreSize(double t, double scaling, std::string type);

		// Others
		double userDefinedInput2vFraction(double userDefinedInput, 
		                                  latticeUDF userDefinedParameters);
	}
};
#endif