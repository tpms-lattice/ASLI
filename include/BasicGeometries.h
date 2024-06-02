/* ==========================================================================
 *  This file is part of ASLI (A Simple Lattice Infiller)
 *  Copyright (C) KU Leuven, 2024
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

#ifndef BASIC_GEOMETRIES_H
#define BASIC_GEOMETRIES_H

#include "ExceptionClass.h" // Custom exception class

/* CGAL headers */
#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

// Surface mesh
#include <CGAL/Surface_mesh/Surface_mesh.h>

// Volume mesh
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

// Volume mesh (implicit domain)
#include <CGAL/Labeled_mesh_domain_3.h>

/* Standard library headers */
#include <algorithm>
#include <cmath>

// Aliases
using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = K::FT;
using Point_3 = K::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;

using Implicit_domain = CGAL::Labeled_mesh_domain_3<K>;
using Implicit_domain_with_features = CGAL::Mesh_domain_with_polyline_features_3<Implicit_domain>;
using VTr_IF = CGAL::Mesh_triangulation_3<Implicit_domain_with_features, CGAL::Default, CGAL::Sequential_tag>::type;
using C3t3_IF = CGAL::Mesh_complex_3_in_triangulation_3<VTr_IF, Implicit_domain_with_features::Corner_index, Implicit_domain_with_features::Curve_index>;
using Mesh_criteria_IF = CGAL::Mesh_criteria_3<VTr_IF>;

namespace BasicGeometries {

	// Declarations
	void implicit2surface(std::string geometry, std::vector<double> dimensions,
		std::vector<Point_3> &points, std::vector<std::vector<std::size_t>> &polygons);

	double implicit_function(std::string geometry, std::vector<double> dimensions, Point_3 p);

	CGAL::Bbox_3 bounding_box(std::string geometry, std::vector<double> dimensions);

	std::vector<std::vector<Point_3>> protected_edges(std::string geometry, std::vector<double> dimensions);

	namespace internal {
		std::vector<Point_3> elipse(Point_3 center, std::vector<double> radius, int nPoints);
	}

	// Error messages
	const std::string INVALID_BASIC_GEOMETRY = "Unavailable internal geometry selected.";
	const std::string INVALID_NUMBER_OF_PARAMETERS = "Number of parameters pased to the internal geometry is not valid.";
	const std::string INVALID_PARAMETER_VALUE = "Dimension(s) of internal geometry must be larger than zero.";
	const std::string INVALID_DIMENSIONS = "External Dimensions of the internally generated geometry must lie between 1e-3 and 1e3.";
	const std::string LARGE_ORDER_OF_MAGNITUTE_DIFFERENCE = "Order of difference between the external dimensions of the internally generated geometry cannot be more than two.";
	const std::string MESH_WITH_BORDER_EDGES = "Polygon has border edges.";
}
#endif