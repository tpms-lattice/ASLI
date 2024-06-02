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

/* BasicGeometries contains the functions required to generate a few basic input 
 * geometries, i.e. cuboid, cylinder and spheroid.
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

#include "BasicGeometries.h"

void BasicGeometries::implicit2surface(std::string geometry, std::vector<double> dimensions,
	std::vector<Point_3> &points, std::vector<std::vector<std::size_t>> &polygons) {
	/* Discretizes an implicit surface.
	 * Inputs :
	 *  geometry   : Label of geometry to be discretized
	 *  dimensions : Dimensions of the basic geometry to discretize
	 * Outputs :
	 *  points     : Vertices
	 *  polygons   : Faces
	 */

	// User settings
	FT facetAngle = 30; // Set to at most 30 degrees
	FT facetDistance = 0.001;
	FT relativeErrorBound = 1e-6; 

	// 
	CGAL::Bbox_3 boundingBox = bounding_box(geometry, dimensions);

	FT minSize = std::min(std::min((boundingBox.xmax() - boundingBox.xmin()), 
	                               (boundingBox.ymax() - boundingBox.ymin())), 
	                               (boundingBox.zmax() - boundingBox.zmin()));

	FT maxSize = std::max(std::max((boundingBox.xmax() - boundingBox.xmin()), 
	                               (boundingBox.ymax() - boundingBox.ymin())), 
	                               (boundingBox.zmax() - boundingBox.zmin()));

	if ( minSize < 1e-3 || maxSize > 1e3)
		throw std::runtime_error(INVALID_DIMENSIONS);
	if ( static_cast<int>(std::floor(std::log10(maxSize))) -
	     static_cast<int>(std::floor(std::log10(minSize))) > 2)
		throw std::runtime_error(LARGE_ORDER_OF_MAGNITUTE_DIFFERENCE);

	// 
	FT minDiameter;
	if ( geometry == "cuboid" )
		minDiameter = 1; // Set to 1 since all faces are flat
	else if ( geometry == "cylinder" )
		minDiameter = 2*(*std::min_element(std::begin(dimensions) + 1, std::end(dimensions)));
	else
		minDiameter = 2*(*std::min_element(std::begin(dimensions), std::end(dimensions)));

	// Domain
	auto myImplicitFunction = [geometry, dimensions](const Point_3 &p) { return implicit_function(geometry, dimensions, p); }; // To pass function as a static member
	Implicit_domain_with_features domain = Implicit_domain::create_implicit_mesh_domain(myImplicitFunction,
	                                                                                    boundingBox,
	                                                                                    CGAL::parameters::relative_error_bound(relativeErrorBound * minDiameter));

	// Protected features
	std::vector<std::vector<Point_3>> protected_features = protected_edges(geometry, dimensions);
	if ( !protected_features.empty() )
		domain.add_features(protected_features.begin(), protected_features.end());

	Mesh_criteria_IF criteria(//CGAL::parameters::edge_distance = facetDistance, // Uncomment with update to CGAL 6.x
	                          CGAL::parameters::facet_distance = facetDistance * minDiameter,//
	                          CGAL::parameters::facet_angle = facetAngle);

	// Generate mesh
	C3t3_IF c3t3 = CGAL::make_mesh_3<C3t3_IF>(domain, criteria,
	                                          CGAL::parameters::manifold(),
	                                          CGAL::parameters::no_exude(),
	                                          CGAL::parameters::no_perturb(),
	                                          CGAL::parameters::no_lloyd(),
	                                          CGAL::parameters::no_odt());

	// Extract surface triangulation
	Surface_mesh surfaceMesh;
	CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, surfaceMesh);
	CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(surfaceMesh, 
	                                                            points, polygons);

	// Clean up vertices
	CGAL::Polygon_mesh_processing::remove_isolated_vertices(surfaceMesh);

	// Reorient all faces coherently
	if ( CGAL::is_closed(surfaceMesh) == true )
		CGAL::Polygon_mesh_processing::orient(surfaceMesh, CGAL::parameters::outward_orientation(true));
	else
		throw ExceptionError(MESH_WITH_BORDER_EDGES, nullptr);

	//std::cout << "\n  SURFACE TRIANGULATION: " << std::endl;
	//std::cout << "  Number of vertices: " << surfaceMesh.number_of_vertices() << std::endl;
	//std::cout << "  Number of facets: " << surfaceMesh.number_of_faces() << "\n" << std::endl;

	// Cleanup
	c3t3.clear();
}

double BasicGeometries::implicit_function(std::string geometry, std::vector<double> dimensions, Point_3 p) {
	/* Collection of level-set function of basic geometries
	*  Inputs:
	*    geometry   : Geometry label
	*    dimensions : Dimensions of the geometry requested
	*    p          : Point at which to evaluate the implicit function of the 
	*                 requested geometry
	*  Return:
	*    Output of the implicit function evaluation
	*/

	Point_3 point(p.x(), p.y(), p.z());

	// Check input dimensions for negative sizes
	for(const double& i : dimensions) 
		if ( i <= 0 )
			throw std::runtime_error(INVALID_PARAMETER_VALUE);

	// Internal geometries
	if ( geometry == "cuboid" ) {
		double length_x, length_y, length_z;
		if ( dimensions.size() == 1 ) {
			length_x = dimensions[0];
			length_y = dimensions[0];
			length_z = dimensions[0];
		} else if ( dimensions.size() == 3 ) {
			length_x = dimensions[0];
			length_y = dimensions[1];
			length_z = dimensions[2];
		} else {
			throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
		}

		const double xabs = std::abs(p.x()), yabs = std::abs(p.y()), zabs = std::abs(p.z());
		return std::max(std::max(xabs/(length_x/2.0), yabs/(length_y/2.0)), zabs/(length_z/2.0)) - 1;

	} else if ( geometry == "cylinder" ) {
		double height, r_x, r_y;
		if ( dimensions.size() == 2 ) {
			height = dimensions[0];
			r_x = dimensions[1];
			r_y = dimensions[1];
		}	else if ( dimensions.size() == 3 ) {
			height = dimensions[0];
			r_x = dimensions[1];
			r_y = dimensions[2];
		} else {
			throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
		}

		return std::max(std::sqrt(std::pow(p.x(), 2)/std::pow(r_x, 2) + 
		                          std::pow(p.y(), 2)/std::pow(r_y, 2)) - 1, 
		                std::abs(p.z()) - height/2.0);

	} else if ( geometry == "spheroid" ) {
		double r_x, r_y, r_z;
		if ( dimensions.size() == 1 ) {
			r_x = dimensions[0];
			r_y = dimensions[0];
			r_z = dimensions[0];
		} else if ( dimensions.size() == 3 ) {
			r_x = dimensions[0];
			r_y = dimensions[1];
			r_z = dimensions[2];
		} else {
			throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
		}

		return std::pow(p.x()/r_x, 2) + std::pow(p.y()/r_y, 2) + std::pow(p.z()/r_z, 2) - 1;

	} else {
		throw std::runtime_error(INVALID_BASIC_GEOMETRY);
	}
}

CGAL::Bbox_3 BasicGeometries::bounding_box(std::string geometry, std::vector<double> dimensions) {
	/* Bounding boxes of basic geometries
	*  Inputs:
	*    geometry   : Geometry label
	*    dimensions : Dimensions of the geometry requested
	*  Return:
	*    Bounding box of the selected geometry
	*/

	// Check input dimensions for negative sizes
	for(const double& i : dimensions) 
		if ( i <= 0 )
			throw std::runtime_error(INVALID_PARAMETER_VALUE);

		// Bounding boxes of internal geometries
		if ( geometry == "cuboid" ) {
			double length_x, length_y, length_z;
			if ( dimensions.size() == 1 ) {
				length_x = dimensions[0];
				length_y = dimensions[0];
				length_z = dimensions[0];
			} else if ( dimensions.size() == 3 ) {
				length_x = dimensions[0];
				length_y = dimensions[1];
				length_z = dimensions[2];
			} else {
				throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
			}

		CGAL::Bbox_3 boundingBox(-length_x/2.0, -length_y/2.0, -length_z/2.0,
		                          length_x/2.0,  length_y/2.0,  length_z/2.0);
		return boundingBox;

	} else if ( geometry == "cylinder" ) {
		double height, r_x, r_y;
		if ( dimensions.size() == 2 ) {
			height = dimensions[0];
			r_x = dimensions[1];
			r_y = dimensions[1];
		}	else if ( dimensions.size() == 3 ) {
			height = dimensions[0];
			r_x = dimensions[1];
			r_y = dimensions[2];
		} else {
			throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
		}

		CGAL::Bbox_3 boundingBox(-r_x, -r_y, -height/2.0, r_x, r_y, height/2.0);
		return boundingBox;
		
	} else if ( geometry == "spheroid" ) {
		double r_x, r_y, r_z;
		if ( dimensions.size() == 1 ) {
			r_x = dimensions[0];
			r_y = dimensions[0];
			r_z = dimensions[0];
		} else if ( dimensions.size() == 3 ) {
			r_x = dimensions[0];
			r_y = dimensions[1];
			r_z = dimensions[2];
		} else {
			throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
		}

		CGAL::Bbox_3 boundingBox(-r_x, -r_y, -r_z, r_x, r_y, r_z);
		return boundingBox;

	} else {
		throw std::runtime_error(INVALID_BASIC_GEOMETRY);
	}
}

std::vector<std::vector<Point_3>> BasicGeometries::protected_edges(std::string geometry, std::vector<double> dimensions) {
	/* Edges of basic geometries
	*  Inputs:
	*    geometry   : Geometry label
	*    dimensions : Dimensions of the geometry requested
	*  Return:
	*    Edges of selected geometry as polylines
	*/

	// Check input dimensions for negative sizes
	for(const double& i : dimensions) 
		if ( i <= 0 )
			throw std::runtime_error(INVALID_PARAMETER_VALUE); 

	// Protected edges of internal geometries
	if ( geometry == "cuboid" ) {
		double length_x, length_y, length_z;
		if ( dimensions.size() == 1 ) {
			length_x = dimensions[0];
			length_y = dimensions[0];
			length_z = dimensions[0];
		} else if ( dimensions.size() == 3 ) {
			length_x = dimensions[0];
			length_y = dimensions[1];
			length_z = dimensions[2];
		} else {
			throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
		}

		std::vector<std::vector<Point_3>> polylines = {
			{Point_3(-length_x/2.0, -length_y/2.0, -length_z/2.0), Point_3( length_x/2.0, -length_y/2.0, -length_z/2.0)},
			{Point_3(-length_x/2.0, -length_y/2.0, -length_z/2.0), Point_3(-length_x/2.0,  length_y/2.0, -length_z/2.0)},
			{Point_3(-length_x/2.0, -length_y/2.0, -length_z/2.0), Point_3(-length_x/2.0, -length_y/2.0,  length_z/2.0)},
			{Point_3( length_x/2.0, -length_y/2.0, -length_z/2.0), Point_3( length_x/2.0,  length_y/2.0, -length_z/2.0)},
			{Point_3( length_x/2.0, -length_y/2.0, -length_z/2.0), Point_3( length_x/2.0, -length_y/2.0,  length_z/2.0)},
			{Point_3(-length_x/2.0,  length_y/2.0, -length_z/2.0), Point_3( length_x/2.0,  length_y/2.0, -length_z/2.0)},
			{Point_3(-length_x/2.0,  length_y/2.0, -length_z/2.0), Point_3(-length_x/2.0,  length_y/2.0,  length_z/2.0)},
			{Point_3(-length_x/2.0, -length_y/2.0,  length_z/2.0), Point_3( length_x/2.0, -length_y/2.0,  length_z/2.0)},
			{Point_3(-length_x/2.0, -length_y/2.0,  length_z/2.0), Point_3(-length_x/2.0,  length_y/2.0,  length_z/2.0)},
			{Point_3( length_x/2.0,  length_y/2.0, -length_z/2.0), Point_3( length_x/2.0,  length_y/2.0,  length_z/2.0)},
			{Point_3(-length_x/2.0,  length_y/2.0,  length_z/2.0), Point_3( length_x/2.0,  length_y/2.0,  length_z/2.0)},
			{Point_3( length_x/2.0, -length_y/2.0,  length_z/2.0), Point_3( length_x/2.0,  length_y/2.0,  length_z/2.0)}
		};
		return polylines;

	} else if ( geometry == "cylinder" ) {
		double height;
		std::vector<double> r;
		if ( dimensions.size() == 2 ) {
			height = dimensions[0];
			r.push_back(dimensions[1]);
		}	else if ( dimensions.size() == 3 ) {
			height = dimensions[0];
			r.push_back(dimensions[1]);
			r.push_back(dimensions[2]);
		} else {
			throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
		}

		std::vector<std::vector<Point_3>> polylines;
		polylines.push_back(internal::elipse(Point_3(0,0,height/2), r, 360));
		polylines.push_back(internal::elipse(Point_3(0,0,-height/2), r, 360));
		return polylines;

	} else if ( geometry == "spheroid") {
		return std::vector<std::vector<Point_3>>();

	} else {
		throw std::runtime_error(INVALID_BASIC_GEOMETRY);
	}
}

std::vector<Point_3> BasicGeometries::internal::elipse(Point_3 center, std::vector<double> radius, int nPoints) {
	/* Edges of elipse
	*  Inputs:
	*    center : Position of the elipse center
	*    radius : Radius of the circe
	*    nPoints: Number of points in which the elipse's circumference is divided
	*  Return:
	*    The elipse's circumference as a polyline 
	*/

	double r_x, r_y;
	if ( radius.size() == 1 ) {
		r_x = radius[0];
		r_y = radius[0];
	} else if ( radius.size() == 2 ) {
		r_x = radius[0];
		r_y = radius[1];
	} else {
		throw std::runtime_error(INVALID_NUMBER_OF_PARAMETERS);
	}

	std::vector<Point_3> polyline;
	for ( int i = 0; i < nPoints; ++i ) {
		double angle = 2.0 * CGAL_PI * i / nPoints;
		double x = center.x() + r_x * cos(angle);
		double y = center.y() + r_y * sin(angle);
		polyline.push_back(Point_3(x, y, center.z()));
	}
	if ( polyline.front() != polyline.back() ) {
		polyline.push_back(polyline.front());  // Close the polyline
	}
	return polyline;
}