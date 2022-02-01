/* ==========================================================================
 *  This file is part of ASLI (A Simple Lattice Infiller)
 *  Copyright (C) KU Leuven, 2019-2022
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

#include "Infill.h"
#include "Mesh.h"
#include "TrilinearInterpolation.h"
#include "Filter.h"
#include "IO_ASLI.h"

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
#include <filesystem>
#include <cstdlib>
#include <cstdio>

class ASLI {
public:
	/* Constructors|Destructors */
	ASLI(std::string configFile);
	ASLI(std::string _stlFile, std::string _tapFile, std::string _sapFile,
	     std::string _fapFile, latticeType _lt_type, latticeSize _lt_size,
	     latticeFeature _lt_feature, meshSettings _me_settings);
	~ASLI();

	/* Parameters */
	// Constants
	const double PI = 4.0*std::atan(1.0);

	// Outer shell
	outerShell shell;

	// Lattice
	latticeType lt_type;
	latticeTypeTree lt_typeTree;
	latticeSize lt_size;
	latticeFeature lt_feature;

	// Mesh
	meshSettings me_settings;

private:
	// Allocate variables and set values
	void SetUp(std::string configFile);
	void SetUpLattice();
	void SetUpKdtree(alglib::real_2d_array *alglib_data, 
	                 alglib::integer_1d_array *tags, alglib::kdtree *kdt);
	void SetUpTypeInterpolator(alglib::real_2d_array coordinates,
	                           std::vector<double> weights, 
	                           modelData *interpolationModel);

	void LoadInputFiles();

	// Input files
	std::string stlFile;
	std::string tapFile, sapFile, fapFile;
};

#endif