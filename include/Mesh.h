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

#ifndef MESH_H
#define MESH_H

/* CGAL headers */
#ifdef CGAL_MESH
	#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	#include <CGAL/IO/STL.h>
#endif

/* tetgen headers */
#ifdef MMG_MESH
	#include "tetgen.h"
#endif

/* Standard library headers */
#include <string>

#ifdef CGAL_MESH
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_3 Point_3_;
#endif

#if defined CGAL_MESH && MMG_MESH
	struct outerShell {
		tetgenio tetgenPoints;
		std::vector<Point_3_> points;
		std::vector<std::vector<std::size_t>> polygons;
	};
#elif defined CGAL_MESH
	struct outerShell {
		std::vector<Point_3_> points;
		std::vector<std::vector<std::size_t>> polygons;
	};
#elif defined MMG_MESH
	struct outerShell {
		tetgenio tetgenPoints;
	};
#endif

struct meshSettings {
	// General settings
	std::string mesher; 	// Mesher to use
	int n_threads;      	// Number of parallel processes
	std::string side;   	// Side of scaffolds to mesh
	bool volumeMesh; 
	std::string output; 	// Output location
	std::string STLFormat;// STL format: ASCII or Binary
	bool exportls; 				// If true only level-set sol data is computed and saved

	// CGAL settings
	double CGAL_facetAngle;
	double CGAL_facetSize;
	double CGAL_facetDistance;
	double CGAL_cellRadiusEdgeRatio;
	double CGAL_cellSize;
	bool CGAL_preserveEdges;
	double CGAL_poissonOffset;
	double CGAL_edgesProtectionAngle;

	// TETGEN settings
	double TETGEN_hvol;

	// MMG settings
	double MMG_hinitial;
	double MMG_hmin;
	double MMG_hmax;
	double MMG_hausd;
	double MMG_hgrad;
};

#endif