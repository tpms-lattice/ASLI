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

#ifndef MESH_CGAL_H
#define MESH_CGAL_H

#include "Infill.h"
#include "TicToc.h"

#include "BasicGeometries.h"
#include "ExceptionClass.h" // Custom exception class

/* CGAL headers */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// IO
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/STL.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/write_xyz_points.h>

// Polygon processing
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Polygon_mesh_processing/repair.h> // Used as temporary fix for a bug, i.e. remove isolated vertices from polygon mesh

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h> // to read stl in mmg surface implicit...???

// Boost graph tools 
#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/copy_face_graph.h>

// Poisson surface reconstruction
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
// Point Set Processing
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/remove_outliers.h>

//// Periodic meshing
//#ifndef CGAL_CONCURRENT_MESH_3
//	#include <CGAL/Periodic_3_mesh_3/config.h>
//	#include <CGAL/make_periodic_3_mesh_3.h>
//	#include <CGAL/optimize_periodic_3_mesh_3.h>
//	#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
//	#include "periodic_replicator.h"
//	#include <CGAL/Periodic_3_mesh_triangulation_3.h>
//	#include <CGAL/number_type_config.h> // CGAL_PI
//#endif

// 3D surface mesh
#include <CGAL/Surface_mesh/Surface_mesh.h>
//#include "OffsetCGAL.h"

// 3D volume mesh
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

// 3D volume mesh (implicit domain)
#include <CGAL/Labeled_mesh_domain_3.h>

// 3D volume mesh (polyhedral domain)
#include <CGAL/Polyhedral_mesh_domain_3.h>

// 3D volume mesh (polyhedral domain with features)
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

// 3D volume mesh (hybrid domain with features)
#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/tags.h>

// Other
#include <CGAL/Bbox_3.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/remove_far_points_in_mesh_3.h>

/* Standard library headers */
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <time.h>
#include <cstdlib>
#include <utility>      // std::move
//#include <functional>

/* CGAL types */
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Sphere_3 Sphere_3;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

// Poisson surface reconstruction
typedef Kernel::Vector_3 Vector_3;

typedef std::pair<Point_3, Vector_3> Point_with_normal;
typedef std::vector<Point_with_normal> PointList;

typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;

//// Polygon processing
//typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> PolygonMesh;
//typedef boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor_; // Edge property to indicate whether constrained

// Surface mesh
typedef CGAL::Surface_mesh<Point_3> SurfaceMesh;
//typedef Triangle_mesh::Property_map<SurfaceMesh::Edge_index,bool> Constrained_edge_map;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor edge_descriptor;  // Edge property to indicate whether constrained
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor; 

// 3D volume meshing (implicit domain)
typedef CGAL::Labeled_mesh_domain_3<Kernel> Implicit_domain;
typedef CGAL::Mesh_triangulation_3<Implicit_domain, CGAL::Default, Concurrency_tag>::type ITr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<ITr> C3t3;
typedef CGAL::Mesh_criteria_3<ITr> Mesh_criteria;

// 3D volume meshing (polyhedron domain)
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel> Polyhedron_domain;
typedef CGAL::Mesh_triangulation_3<Polyhedron_domain, CGAL::Default, Concurrency_tag>::type PTr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<PTr> P_C3t3;
typedef CGAL::Mesh_criteria_3<PTr> P_Mesh_criteria;

// 3D volume meshing (implicit domain with features)
typedef CGAL::Mesh_domain_with_polyline_features_3<Implicit_domain> F_Implicit_domain;
typedef CGAL::Mesh_triangulation_3<F_Implicit_domain, CGAL::Default, Concurrency_tag>::type F_ITr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<F_ITr> F_I_C3t3;
typedef CGAL::Mesh_criteria_3<F_ITr> F_I_Mesh_criteria;

// 3D volume meshing (polyhedron domain with features)
typedef CGAL::Mesh_polyhedron_3<Kernel>::type F_Polyhedron;

typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel> F_Polyhedron_domain;
typedef CGAL::Mesh_triangulation_3<F_Polyhedron_domain, CGAL::Default, Concurrency_tag>::type F_PTr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<F_PTr, F_Polyhedron_domain::Corner_index, F_Polyhedron_domain::Curve_index> F_C3t3;
typedef CGAL::Mesh_criteria_3<F_PTr> F_Mesh_criteria;

//// Periodic meshing
//#ifndef CGAL_CONCURRENT_MESH_3
//	typedef Kernel::Iso_cuboid_3 Iso_cuboid;
//
//	typedef FT (Function)(const Point_3&);
//
//	typedef CGAL::Periodic_3_mesh_triangulation_3<Implicit_domain>::type Periodic_Tr;
//	typedef CGAL::Mesh_complex_3_in_triangulation_3<Periodic_Tr> Periodic_C3t3;
//	typedef CGAL::Mesh_criteria_3<Periodic_Tr> Periodic_mesh_criteria;
//#endif



using Subdomain_index = int;
using Surface_patch_index = unsigned char;
using Curve_index = char;
using Corner_index = short;
using Cb = CGAL::Simplicial_mesh_cell_base_3<Subdomain_index, Surface_patch_index>;
using Vb = CGAL::Simplicial_mesh_vertex_base_3<Kernel, Subdomain_index, Surface_patch_index,
                                               Curve_index, Corner_index>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, Concurrency_tag>;
using Triangulation = CGAL::Triangulation_3<Kernel, Tds>;


class Hybrid_domain {
	const Implicit_domain& implicit_domain;
	const Polyhedron_domain& polyhedron_domain;
	std::string side_;

	public:
	Hybrid_domain(const Implicit_domain& implicit_domain, const Polyhedron_domain& polyhedron_domain, const std::string& side)
	              : implicit_domain(implicit_domain) , polyhedron_domain(polyhedron_domain) , side_(side) {}
	
	// types required by the 'MeshDomain_3' concept
	typedef int Surface_patch_index;
	typedef int Subdomain_index;
	typedef int Index;
	typedef Kernel R;
	typedef Kernel::Point_3 Point_3;
	typedef std::tuple<Point_3, Index, int> Intersection;
  
	CGAL::Bbox_3 bbox() const { return implicit_domain.bbox() + polyhedron_domain.bbox(); }

	struct Construct_initial_points {
		Construct_initial_points(const Hybrid_domain& domain) : r_domain_(domain) {}
		template<class OutputIterator> OutputIterator operator()(OutputIterator pts, const int n = 20) const {
			//construct initial points on implicit domain
			typedef Implicit_domain::Index Implicit_Index;
			std::vector<std::pair<Point_3, Implicit_Index> > implicit_points_vector;
			Implicit_domain::Construct_initial_points cstr_implicit_initial_points =
				r_domain_.implicit_domain.construct_initial_points_object();
			cstr_implicit_initial_points(std::back_inserter(implicit_points_vector), n/2);
			for(std::size_t i = 0, end = implicit_points_vector.size(); i<end; ++i) {
				*pts++ = std::make_pair(implicit_points_vector[i].first, 2);
			}
			//construct initial points on polyhedral domain
			typedef Polyhedron_domain::Index Polyhedron_Index;
			std::vector<std::pair<Point_3, Polyhedron_Index> > polyhedron_points_vector;
			Polyhedron_domain::Construct_initial_points cstr_polyhedron_initial_points =
				r_domain_.polyhedron_domain.construct_initial_points_object();
			cstr_polyhedron_initial_points(std::back_inserter(polyhedron_points_vector), n/2);
			for(std::size_t i = 0, end = polyhedron_points_vector.size(); i<end; ++i) {
				*pts++ = std::make_pair(polyhedron_points_vector[i].first, 1);
			}
			return pts;
		}

	private:
	const Hybrid_domain& r_domain_;
	}; // end Construct_initial_points_object

	Construct_initial_points construct_initial_points_object() const {
		return Construct_initial_points(*this);
	}

	struct Is_in_domain {
		Is_in_domain(const Hybrid_domain& domain) : r_domain_(domain) {}
		boost::optional<Subdomain_index> operator()(const Kernel::Point_3& p) const {
			boost::optional<Subdomain_index> implicit_subdomain_index = 
				r_domain_.implicit_domain.is_in_domain_object()(p);
			boost::optional<Subdomain_index> polyhedron_subdomain_index = 
				r_domain_.polyhedron_domain.is_in_domain_object()(p);

			if ( r_domain_.side_ == "scaffold" ) {
				if (!implicit_subdomain_index && polyhedron_subdomain_index)
					return 0;
				else
					return r_domain_.polyhedron_domain.is_in_domain_object()(p);
			} else if ( r_domain_.side_ == "void" ) {
				if(implicit_subdomain_index && polyhedron_subdomain_index)
					return 0;
				else
					return r_domain_.polyhedron_domain.is_in_domain_object()(p);
			} else {
				if (!implicit_subdomain_index && polyhedron_subdomain_index)
					return 2;
				else 
					return r_domain_.polyhedron_domain.is_in_domain_object()(p);
			}

		}

		private:
		const Hybrid_domain& r_domain_;
	}; // end Is_in_domain_object

	Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }
  
	struct Construct_intersection {
		Construct_intersection(const Hybrid_domain& domain) : r_domain_(domain) {}
		template <typename Query> Intersection operator()(const Query& query) const {
			
			using boost::get;
		
			// intersection with implicit domain (within the polyhedral domain)
			Implicit_domain::Intersection implicit_inter =
				r_domain_.implicit_domain.construct_intersection_object()(query);

			// intersection with polyhedral domain
			Polyhedron_domain::Intersection polyhedron_inter =
				r_domain_.polyhedron_domain.construct_intersection_object()(query);
				
			if ( r_domain_.side_ == "scaffold" ) {
				// if found intersection with implicit domain (within the polyhedral domain), return it
				if ( get<2>(implicit_inter) != 0 ) {
					const Point_3 inter_point = get<0>(implicit_inter);
					if ( r_domain_.polyhedron_domain.is_in_domain_object()(inter_point) ) // Inner implicit surface
						return Intersection(inter_point, 3, get<2>(implicit_inter));
				}

				// if found intersection with polyhedral domain, return it
				if ( get<2>(polyhedron_inter) != 0 ) {
					const Point_3 inter_point = get<0>(polyhedron_inter);
					if ( r_domain_.implicit_domain.is_in_domain_object()(inter_point) ) // Scaffold
						return Intersection(inter_point, 3, get<2>(polyhedron_inter));
					else // Void
						return Intersection(inter_point, 4, get<2>(polyhedron_inter));
				}

			} else if ( r_domain_.side_ == "void" ) {
				// if found intersection with implicit domain (within the polyhedral domain), return it
				if ( get<2>(implicit_inter) != 0 ) {
					const Point_3 inter_point = get<0>(implicit_inter);
					if ( r_domain_.polyhedron_domain.is_in_domain_object()(inter_point) ) // Inner implicit surface
						return Intersection(inter_point, 3, get<2>(implicit_inter));
				}

				// if found intersection with polyhedral domain, return it
				if ( get<2>(polyhedron_inter) != 0 ) {
					const Point_3 inter_point = get<0>(polyhedron_inter);
					if ( r_domain_.implicit_domain.is_in_domain_object()(inter_point) ) // Scaffold
						return Intersection(inter_point, 4, get<2>(polyhedron_inter));
					else // Void
						return Intersection(inter_point, 3, get<2>(polyhedron_inter));
				}

			} else {
				// if found intersection with implicit domain (within the polyhedral domain), return it
				if ( get<2>(implicit_inter) != 0 ) {
					const Point_3 inter_point = get<0>(implicit_inter);
					if ( r_domain_.polyhedron_domain.is_in_domain_object()(inter_point) ) // Inner implicit surface
						return Intersection(inter_point, 3, get<2>(implicit_inter));
				}

				// if found intersection with polyhedral domain, return it
				if ( get<2>(polyhedron_inter) != 0 ) {
					const Point_3 inter_point = get<0>(polyhedron_inter);
					if ( r_domain_.implicit_domain.is_in_domain_object()(inter_point) ) // Scaffold
						return Intersection(inter_point, 4, get<2>(polyhedron_inter));
					else // Void
						return Intersection(inter_point, 5, get<2>(polyhedron_inter));
				}
			}

			//no intersection found
			return Intersection();
		}

		private:
		const Hybrid_domain& r_domain_;
	}; // end Construct_intersection_object

	Construct_intersection construct_intersection_object() const {
		return Construct_intersection(*this);
	}

	//Index types converters
	Index index_from_surface_patch_index(const Surface_patch_index& index) const
		{ return index; }
	Index index_from_subdomain_index(const Subdomain_index& index) const
		{ return index; }
	Surface_patch_index surface_patch_index(const Index& index) const
		{ return index; }
	Subdomain_index subdomain_index(const Index& index) const
		{ return index; }
}; // end Hybrid_domain class

// 3D volume meshing (hybrid domain with features)
typedef CGAL::Mesh_domain_with_polyline_features_3<Hybrid_domain> H_domain;

typedef CGAL::Mesh_triangulation_3<H_domain, CGAL::Default, Concurrency_tag>::type H_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<H_Tr> H_C3t3;
typedef CGAL::Mesh_criteria_3<H_Tr>     H_Mesh_criteria;
typedef H_Mesh_criteria::Edge_criteria  H_Edge_criteria;
typedef H_Mesh_criteria::Facet_criteria H_Facet_criteria;
typedef H_Mesh_criteria::Cell_criteria  H_Cell_criteria;



/* Constrained maps templates */
template <typename ECMap>
struct Is_constrained {
	ECMap ecmap;
	template <typename Edge>
	bool operator()(const Edge& e) const { return get(ecmap, e); }
};
template <typename VPMap>
struct Is_feature {
	VPMap vpmap;
	template <typename Vertex>
	bool operator()(const Vertex& v) const { return get(vpmap, v); }
};

struct Polyline_visitor {
	std::vector<std::vector<Point_3>> &polylines;
	const SurfaceMesh &mesh;

	void start_new_polyline() { polylines.emplace_back(); }
	void add_node(vertex_descriptor vd) {	polylines.back().push_back(mesh.point(vd));	}
	void end_polyline(){}
};

// Error messages
namespace CGAL_ERRMSG {
	const std::string FAILED_TO_OPEN_FILE = "Error opening file: ";

	const std::string INVALID_WALL_SIZE = "Invalid wall size";
	const std::string INVALID_PORE_SIZE = "Invalid pore size";
}

#endif