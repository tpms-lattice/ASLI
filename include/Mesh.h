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

#ifndef MESH_H
#define MESH_H

#include "MeshCGAL.h"
#include "MeshMMG.h"
#include "ExceptionClass.h" // Custom exception class

/* CGAL headers */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/STL.h>

/* Standard library headers */
#include <string>
#include <filesystem>

struct polygonSoup {
	std::vector<CGAL::Exact_predicates_inexact_constructions_kernel::Point_3> points;
	std::vector<std::vector<std::size_t>> polygons;
};

struct meshSettings {
	// General settings
	std::string mesher;                             // Mesher to use
	int n_threads                   = 1;            // Number of parallel processes
	bool isVolumeMesh               = true;
	double edgeProtectionAngle      = 70;
	std::filesystem::path inputFile = "";           // Input geometry
	std::string internalGeometry;                   // 
	std::vector<double> internalGeometryParameters; // 
	std::filesystem::path output    = "outputs/";   // Output path
	std::string STLFormat           = "Binary";     // STL format: ASCII or Binary
	double verbosity                = -1;           // Value between -1 and 6
	
	double elementSize;
	double threshold                = 1.0/30.0;     // as a fraction of the local unit cell size

	// CGAL settings
	double CGAL_facetAngle          = 30;
	double CGAL_facetSize;
	double CGAL_facetDistance       = 0.01;
	double CGAL_cellRadiusEdgeRatio = 3.0;
	double CGAL_relativeErrorBound  = 1e-3;
	double CGAL_cellSize;
	double CGAL_edgeSize;
	double CGAL_minEdgeSize         = 0;
	double CGAL_poissonOffset       = 0.5;

	// MMG settings
	double MMG_hinitial = 0.36;
	double MMG_hsiz;
	double MMG_hausd    = 0.3;
	double MMG_hgrad    = 1.3;
	bool MMG_exportLS   = false;        // If true only level-set sol data is computed and saved
	double MMG_memory   = 0;
};

/* Sizing fields */
template <typename Domain> 
struct TPMS_dependent_cellsize_field {
	typedef ::FT _FT;
	typedef Point_3 _Point_3;
	typedef typename Domain::Index _Index;
	
	latticeType lt_type;
	latticeSize lt_size;
	latticeFeature lt_feature;
	double threshold;
	double parameter;

	_FT operator()(const _Point_3& p, const int, const _Index&) const {
		Point point(p.x(), p.y(), p.z());
		double unitCellSize = Infill::sizing_function(point, lt_size, "");
		unitCellSize *= parameter;

		//if (unitCellSize > 0)
			return unitCellSize; // Unit cell must be positive!
		//else
		//	throw ExceptionError("Unit cell size is less or equal to zero", nullptr);
	}
};

template <typename Domain> 
struct TPMS_dependent_wallsize_field {
	typedef ::FT _FT;
	typedef Point_3 _Point_3;
	typedef typename Domain::Index _Index;
	
	latticeType lt_type;
	latticeSize lt_size;
	latticeFeature lt_feature;
	double threshold;
	double parameter;

	_FT operator()(const _Point_3& p, const int, const _Index&) const {
		Point point(p.x(), p.y(), p.z());
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);
		double minFeature = (Infill::TPMS_function(point,lt_type,lt_size,lt_feature) <= 0) ? localSize.wallSize : localSize.poreSize;

		double unitCellSize = Infill::sizing_function(point, lt_size, "");
		if (minFeature <= unitCellSize*threshold) { minFeature = unitCellSize; }

		minFeature *= parameter;
		//if (minFeature > 0)
			return minFeature; // Unit cell must be positive!
		//else
		//	throw ExceptionError("minFeature size is less or equal to zero", nullptr);
	}
};

template <typename Domain> 
struct TPMS_dependent_featuresize_field {
	typedef ::FT _FT;
	typedef Point_3 _Point_3;
	typedef typename Domain::Index _Index;
	
	latticeType lt_type;
	latticeSize lt_size;
	latticeFeature lt_feature;
	double threshold;
	double parameter;

	_FT operator()(const _Point_3& p, const int, const _Index&) const {
		Point point(p.x(), p.y(), p.z());
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);

		double feature_size;
		if (lt_type.side == "scaffold")
			feature_size = localSize.wallSize;
		else if (lt_type.side == "void")
			feature_size = localSize.poreSize;
		else
			feature_size = std::min(localSize.wallSize, localSize.wallSize);
		
		double unitCellSize = Infill::sizing_function(point, lt_size, "");
		if (feature_size <= unitCellSize*threshold) { feature_size = unitCellSize; }
		
		feature_size *= parameter;
		//if (feature_size > 0)
			return feature_size; // Unit cell must be positive!
		//else
		//	throw ExceptionError("minFeature size is less or equal to zero", nullptr);
	}
};

template <typename Domain> 
struct TPMS_dependent_minfeaturesize_field {
	typedef ::FT _FT;
	typedef Point_3 _Point_3;
	typedef typename Domain::Index _Index;
	
	latticeType lt_type;
	latticeSize lt_size;
	latticeFeature lt_feature;
	double threshold;
	double parameter;

	_FT operator()(const _Point_3& p, const int, const _Index&) const {
		Point point(p.x(), p.y(), p.z());
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);
		double minFeature = std::min(localSize.wallSize, localSize.poreSize);

		double unitCellSize = Infill::sizing_function(point, lt_size, "");
		if (minFeature <= unitCellSize*threshold) { minFeature = unitCellSize; }

		minFeature *= parameter;
		//if (minFeature > 0)
			return minFeature; // Unit cell must be positive!
		//else
		//	throw ExceptionError("minFeature size is less or equal to zero", nullptr);
	}
};

namespace Mesh {
	// MMG workflow functions that depend on the CGAL library
	void polehedral2volume(const polygonSoup &shell, const latticeType &lt_type, const latticeSize &lt_size,
		const latticeFeature &lt_feature, const meshSettings &me_settings, F_C3t3 &c3t3);
	void c3t3_to_MMG3D(const F_C3t3 &c3t3, MMG5_pMesh mmgMesh);

	// CGAL workflow functions that depend on the MMG library
	void extractEdges(const polygonSoup &shell, const latticeType &lt_type, const latticeSize &lt_size,
		const latticeFeature &lt_feature, const meshSettings &me_settings, SurfaceMesh &surface_mesh,
		SurfaceMesh::Property_map<edge_descriptor, bool> &ecmap);
};

namespace MeshCGAL {
	// Implicit function meshing
	void implicit2volume(const polygonSoup &shell, const latticeType &lt_type,
		const latticeSize &lt_size, const latticeFeature &lt_feature,
		const meshSettings &me_settings);

	// Implicit function meshing
	void implicit2volume_old(const polygonSoup &shell, const latticeType &lt_type, 
		const latticeSize &lt_size, const latticeFeature &lt_feature, const meshSettings &me_settings);

	// Polehedral meshing
	void polehedral2volume(const F_Polyhedron &polygonSurface,
		const std::vector<std::vector<Point_3>> &protected_features,
		const latticeType &lt_type, const latticeSize &lt_size, const latticeFeature &lt_feature,
		const meshSettings &me_settings, F_C3t3 &c3t3);

	// Polehedral surface meshing
	void surfaceRemesh(const polygonSoup &shell, const latticeType &lt_type,
		const latticeSize &lt_size, const latticeFeature &lt_feature,
		const meshSettings &me_settings, SurfaceMesh &surfaceMeshOut);

	// Implicit poisson surface reconstruction
	void implicitPoissonReconstruction(const polygonSoup &shell, const double &offset, 
		std::vector<double> &boundingBox, Poisson_reconstruction_function *&poissonReconstruction);

	namespace internal {
		double infill(const Point_3 &p, const latticeType &lt_type, const latticeSize &lt_size,
			const latticeFeature &lt_feature, const std::vector<double> &BBoxSize,
			const std::vector<double> &BBoxCenter);
		double signedDistance_old(const Point_3 &p, const latticeType &lt_type, 
			const latticeSize &lt_size, const latticeFeature &lt_feature,
			const Poisson_reconstruction_function *const &poissonReconstruction);

		void ecmap2polylines(const SurfaceMesh::Property_map<edge_descriptor, bool> &ecmap,
			SurfaceMesh &surfaceMesh, std::vector<std::vector<Point_3>> &polylines);

		void writePolylinesToFile(const std::string &filename, const std::vector<std::vector<Point_3>> &polylines);
		template <typename C3T3> void remove_from_complex(const typename C3T3::Subdomain_index &sd_index, C3T3 &c3t3);
	}
};

namespace MeshMMG { 
	void implicit2volume(const polygonSoup &shell, const latticeType &lt_type,
		const latticeSize &lt_size, const latticeFeature &lt_feature,
		const meshSettings &me_settings);

	namespace internal {
		void shell_to_MMGS(const polygonSoup &shell, MMG5_pMesh mmgMesh);
		bool MMG3D_saveSurfaceAsSTL(const MMG5_pMesh &mmgMesh,
			const std::string &format, const char *filename);
	}
};
#endif