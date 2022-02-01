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

#ifndef MESHCGAL_H
#define MESHCGAL_H

#include "Mesh.h"
#include "Infill.h"
#include "TicToc.h"

//#include "meshMMG.h" // TEST!!!

/* CGAL headers */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// IO
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/IO/facets_in_complex_3_to_triangle_mesh.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/Polyhedron_builder_from_STL.h>
#include <CGAL/IO/STL_writer.h>
//#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/IO/write_xyz_points.h>

// Polygon processing
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Polygon_mesh_processing/repair.h> // Used as temporary fix for a bug, i.e. remove isolated vertices from polygon mesh

// Poisson surface reconstruction
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
// Point Set Processing
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/remove_outliers.h>



// 
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include "OffsetCGAL.h"



// 3D volume mesh
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

// 3D volume mesh (with features)
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

// Other
#include <CGAL/Bbox_3.h>
//#include <CGAL/Mesh_3/Dump_c3t3.h>

/* TBB headers */
#ifdef CGAL_CONCURRENT_MESH_3
	//#include <tbb/task_scheduler_init.h> // Deprecated
	#define TBB_PREVIEW_GLOBAL_CONTROL 1
	#include <tbb/global_control.h>
#endif

/* Standard library headers */
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <time.h>
#include <cstdlib>
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

// Polygon processing
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> PolygonMesh;



// Surface mesh
typedef CGAL::Surface_mesh<Point_3> SurfaceMesh;



// 3D volume meshing (from implicit domain)
typedef CGAL::Labeled_mesh_domain_3<Kernel> Implicit_domain;

typedef CGAL::Mesh_triangulation_3<Implicit_domain, CGAL::Default, Concurrency_tag>::type VTr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<VTr> C3t3;

typedef CGAL::Mesh_criteria_3<VTr> Mesh_criteria;

// 3D volume meshing (from polyhedron domain)
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel> Polyhedron_domain;

typedef CGAL::Mesh_triangulation_3<Polyhedron_domain, CGAL::Default, Concurrency_tag>::type PTr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<PTr> P_C3t3;

typedef CGAL::Mesh_criteria_3<PTr> P_Mesh_criteria;

// 3D volume meshing (from polyhedron domain with features)
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel> F_Polyhedron_domain;
//typedef CGAL::Polyhedral_complex_mesh_domain_3<Kernel> F_Polyhedron_domain;

typedef CGAL::Mesh_polyhedron_3<Kernel>::type F_Polyhedron;

typedef CGAL::Mesh_triangulation_3<F_Polyhedron_domain, CGAL::Default, Concurrency_tag>::type F_VTr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<F_VTr, F_Polyhedron_domain::Corner_index, F_Polyhedron_domain::Curve_index> F_C3t3;

typedef CGAL::Mesh_criteria_3<F_VTr> F_Mesh_criteria;

/* Sizing fields */
// Sizing fields for implicit domain
struct TPMS_dependent_wallsize_field_I {
	typedef ::FT _FT;
	typedef Point_3 _Point_3;
	typedef Implicit_domain::Index _Index;
	
	latticeType *lt_type;
	latticeSize *lt_size;
	latticeFeature *lt_feature;
	double *parameter;

	_FT operator()(const _Point_3& p, const int, const _Index&) const {
		Point point(p.x(), p.y(), p.z());
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);
		return localSize.wallSize * (*parameter);
	}
};

struct TPMS_dependent_unitcellsize_field_I {
	typedef ::FT _FT;
	typedef Point_3 _Point_3;
	typedef Implicit_domain::Index _Index;

	latticeSize *lt_size;
	double *parameter;

	_FT operator()(const _Point_3& p, const int, const _Index&) const {
		Point point(p.x(), p.y(), p.z());
		double size = Infill::sizing_function(point, lt_size, "");
	  return size * (*parameter);
  }
};

// Sizing fields for polyhedron domain with features
struct TPMS_dependent_wallsize_field_P {
	typedef ::FT _FT;
	typedef Point_3 _Point_3;
	typedef F_Polyhedron_domain::Index _Index;

	latticeType *lt_type;
	latticeSize *lt_size;
	latticeFeature *lt_feature;
	double *parameter;

	_FT operator()(const _Point_3& p, const int, const _Index&) const {
		Point point(p.x(), p.y(), p.z());
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);
		return localSize.wallSize * (*parameter);
	}
};

struct TPMS_dependent_unitcellsize_field_P {
	typedef ::FT _FT;
	typedef Point_3 _Point_3;
	typedef F_Polyhedron_domain::Index _Index;

	latticeSize *lt_size;
	double *parameter;

	_FT operator()(const _Point_3& p, const int, const _Index&) const {
		Point point(p.x(), p.y(), p.z());
		double size = Infill::sizing_function(point, lt_size, "");
		return size * (*parameter);
	}
};

namespace MeshCGAL {
	// Implicit function meshing
	bool implicit2volume(outerShell &shell, latticeType lt_type, 
		latticeSize lt_size, latticeFeature lt_feature, meshSettings me_settings,
		std::filesystem::path &outputFile_string);

	// Polehedral meshing
	bool polehedral2volume(std::vector<Point_3> points,
		std::vector<std::vector<std::size_t>> polygons, latticeType lt_type,
		latticeSize lt_size, latticeFeature lt_feature, meshSettings me_settings,
		std::filesystem::path outputFile);

	// Implicit poisson surface reconstruction
	void implicitPoissonReconstruction(outerShell &shell,
		 double offset, std::vector<double> &boundingBox, 
		Poisson_reconstruction_function **poissonSurfaceReconstruction);

	namespace internal {
		double Infill(Point_3 p, latticeType *lt_type, latticeSize *lt_size, 
			latticeFeature *lt_feature);
		double signedDistance(Point_3 p, latticeType *lt_type, 
			latticeSize *lt_size, latticeFeature *lt_feature, std::string side,
			Poisson_reconstruction_function *poissonReconstruction);
	}
};
#endif