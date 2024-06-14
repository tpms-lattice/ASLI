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

#include "Mesh.h"

/* MESH
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

void Mesh::polehedral2volume(const polygonSoup &shell, const latticeType &lt_type,
                             const latticeSize &lt_size, const latticeFeature &lt_feature,
                             const meshSettings &me_settings, F_C3t3 &c3t3) {
	/* 
	 * Inputs :
	 *  shell       : Structure containing the general geometry
	 *  lt_type     : Lattice type parameters
	 *  lt_size     : Lattice size parameters
	 *  lt_feature  : Lattice feature parameters
	 *  me_settings : Mesh settings
	 * Outputs :
	 *  c3t3        : 
	 */

	// Mesh user settings
	FT facetAngle = me_settings.CGAL_facetAngle;
	FT facetDistance = me_settings.CGAL_facetDistance;

	FT cellRadiusEdgeRatio = me_settings.CGAL_cellRadiusEdgeRatio;
	FT cellSize = me_settings.CGAL_cellSize;     

	FT edgeProtectionAngle = me_settings.edgeProtectionAngle;

	F_Polyhedron polygonSurface;
	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(shell.points, shell.polygons, polygonSurface);

	// Setup domain
	F_Polyhedron_domain domain(polygonSurface);
	if (edgeProtectionAngle > 0) { domain.detect_features(edgeProtectionAngle); }

	// Discretize volume
	if (lt_type.type != "hybrid" && lt_size.size > 0 && lt_feature.feature_val > 0) { // If lattice is uniform
		featureSize localSize = Infill::featureSize_function(Point(0, 0, 0), lt_type, lt_size, lt_feature);
		double minFeature = std::min(localSize.wallSize, localSize.poreSize);

		// Mesh criteria
		F_Mesh_criteria criteria(CGAL::parameters::edge_size = cellSize * minFeature,
		                         CGAL::parameters::facet_angle = facetAngle, // Min triangle angle (degrees)
		                         CGAL::parameters::facet_distance = facetDistance * cellSize * minFeature, // Surface approximation error
		                         CGAL::parameters::cell_radius_edge_ratio = cellRadiusEdgeRatio, // Mesh tetrahedra radius|edge ratio upper bound
		                         CGAL::parameters::cell_size = cellSize * minFeature); //

		// Mesh generation
		c3t3 = CGAL::make_mesh_3<F_C3t3>(domain, criteria, CGAL::parameters::manifold());
		
	} else {
		// Sizing fields
		TPMS_dependent_wallsize_field<F_Polyhedron_domain> cellSizeField;
		cellSizeField.lt_type = lt_type;
		cellSizeField.lt_size = lt_size;
		cellSizeField.lt_feature = lt_feature;
		cellSizeField.threshold = me_settings.threshold;
		cellSizeField.parameter = cellSize;

		TPMS_dependent_cellsize_field<F_Polyhedron_domain> facetDistanceField;
		facetDistanceField.lt_type = lt_type;
		facetDistanceField.lt_size = lt_size;
		facetDistanceField.lt_feature = lt_feature;
		facetDistanceField.threshold = me_settings.threshold;
		facetDistanceField.parameter = facetDistance;

		// Mesh criteria
		F_Mesh_criteria criteria(CGAL::parameters::edge_size = cellSizeField,
														CGAL::parameters::facet_angle = facetAngle, // Min facet angle (degrees)
														CGAL::parameters::facet_distance = facetDistanceField, // Surface approximation error
														CGAL::parameters::cell_radius_edge_ratio = cellRadiusEdgeRatio, // Mesh tetrahedra radius|edge ratio upper bound
														CGAL::parameters::cell_size = cellSizeField); // Mesh tetrahedra circumradii upper bound

		// Mesh generation
		c3t3 = CGAL::make_mesh_3<F_C3t3>(domain, criteria, CGAL::parameters::manifold());
	}
	
	// Clean up vertices
	c3t3.remove_isolated_vertices();
}

void Mesh::extractEdges(const polygonSoup &shell, const latticeType &lt_type, const latticeSize &lt_size,
                        const latticeFeature &lt_feature, const meshSettings &me_settings,
                        SurfaceMesh &surface_mesh, SurfaceMesh::Property_map<edge_descriptor, bool> &ecmap) {
	/* Splits a surface mesh according to the implicit function.???
	 * Inputs :
	 *  shell        : Structure containing the general geometry
	 *  lt_type      : Lattice type parameters
	 *  lt_size      : Lattice size parameters
	 *  lt_feature   : Lattice feature parameters
	 *  me_settings  : Mesh settings
	 * Outputs :
	 *  surface_mesh : 
	 *  ecmap        :
	 */

	MMG5_pMesh mmgMesh;
	MMG5_pSol mmgLs, mmgMet; // mmgLs = level-set solution, mmgMet = size metric
	int np, nt, na; // np = #vertices, nt = #triangles, na = #edges
	int ierr;

	// Initialize mesh and sol structures
	mmgMesh = NULL;  mmgLs = NULL; mmgMet = NULL;

	MMGS_Init_mesh(MMG5_ARG_start,
	               MMG5_ARG_ppMesh, &mmgMesh,
	               MMG5_ARG_ppLs,   &mmgLs,
	               MMG5_ARG_ppMet,  &mmgMet,
	               MMG5_ARG_end);

	/* ------------------------------- Step 1 ------------------------------- */
	// Convert shell to MMGS mesh.
	MeshMMG::internal::shell_to_MMGS(shell, mmgMesh);

	// Save inital MMGS mesh for DEBUG PURPOSES!!!
	if (me_settings.verbosity > 5) { 
		std::string tempFile = "MMGS_STEP_1.mesh";
		if ( MMGS_saveMesh(mmgMesh, tempFile.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_MESH, nullptr);
	}
	
	/* ------------------------------- Step 2 ------------------------------- */
	// Remesh the surface such that the elements have a size appropiate to capture
	// the level-set describing the infill geometry.

	// Get mesh info and set solution size
	if ( MMGS_Get_meshSize(mmgMesh, &np, &nt, &na) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_MESH_SIZE, nullptr);
	if ( MMGS_Set_solSize(mmgMesh,mmgLs,MMG5_Vertex,np,MMG5_Scalar) != 1 ) 
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL_SIZE, nullptr);
	if ( MMGS_Set_solSize(mmgMesh,mmgMet,MMG5_Vertex,np,MMG5_Scalar) != 1 ) 
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL_SIZE, nullptr);

	// Set met values
	TicToc::tic(); std::cout << "    Computing local size metric... " << std::flush;
	for (int i=1 ; i<=np ; i++) {
		// Get the current point
		double point[3];
		if ( MMGS_Get_vertex(mmgMesh, &(point[0]), &(point[1]), &(point[2]), NULL, NULL, NULL) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_VERTEX, nullptr);
		Point p(point[0], point[1], point[2]);

		// Compute met value for current point
		featureSize localFeatureSize = Infill::featureSize_function(p, lt_type, lt_size, lt_feature);
		double solValueMet = (localFeatureSize.wallSize < localFeatureSize.poreSize) 
			? localFeatureSize.wallSize : localFeatureSize.poreSize;
		solValueMet *= me_settings.MMG_hinitial;

		// Make non-positive element sizes (or bellow some set threshold) equal to the size
		// of the unit cell at that location since there is no infill to capture anyway.
		double unitCellSize = Infill::sizing_function(p, lt_size, "");
		if ( solValueMet <= unitCellSize * me_settings.threshold )
			solValueMet = me_settings.MMG_hinitial * unitCellSize;

		// Set size metric
		if (solValueMet > 0) {
			if ( MMGS_Set_scalarSol(mmgMet, solValueMet, i) != 1 )
				throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL, nullptr);
		} else {
			throw ExceptionError(MMG_ERRMSG::INVALID_SIZE_METRIC, nullptr);
		}
	}

	// Check if the number of given entities match with mesh size
	if ( MMGS_Chk_meshData(mmgMesh, mmgMet) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_MESH_DATA_CHECK, nullptr);

	TicToc::toc("completed in ");

	// Remesh to size adecuate to capture the level-set
	if ( MMGS_Set_iparameter(mmgMesh, mmgMet, MMGS_IPARAM_verbose, me_settings.verbosity) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_IPARAMETER + "(MMGS_IPARAM_verbose)", nullptr);

	TicToc::tic(); std::cout << "    Remeshing input surface based on local metric... " << std::flush;
	ierr = MMGS_mmgslib(mmgMesh,mmgMet);
	TicToc::toc("completed in ");
	
	if ( ierr == MMG5_STRONGFAILURE )
		throw ExceptionError(MMG_ERRMSG::BAD_ENDING_OF_MMGSLIB, nullptr);
	else if ( ierr == MMG5_LOWFAILURE )
		std::cerr << "BAD ENDING OF MMGSLIB" << std::endl;

	// Save remeshed mesh for DEBUG PURPOSES!!!
	if (me_settings.verbosity > 5) { 
		std::string tempFile = "MMGS_STEP_2.mesh";
		if ( MMGS_saveMesh(mmgMesh, tempFile.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_MESH, nullptr);
	}

	/* TESTING ANISOTROPIC MESH START!!! */
	// Optimize mesh sizing distribution by computing an anisotropic metric 
	// and remeshing the surface based on it.
	//MeshMMG::internal::get_adapted_meshS(mmgMesh, mmgMet, lt_type, lt_size, lt_feature, me_settings);
	/* TESTING ANISOTROPIC MESH END!!! */

	/* ------------------------------- Step 3 ------------------------------- */
	// Compute the level-set describing the lattice infill

	// Get mesh info and set solution size
	if ( MMGS_Get_meshSize(mmgMesh, &np, &nt, &na) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_MESH_SIZE, nullptr);
	if ( MMGS_Set_solSize(mmgMesh,mmgLs,MMG5_Vertex,np,MMG5_Scalar) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL_SIZE, nullptr);
	if ( MMGS_Set_solSize(mmgMesh,mmgMet,MMG5_Vertex,np,MMG5_Scalar) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL_SIZE, nullptr);

	// Compute the level-set values within the volume to be provided with an 
	// infill and set them as the sol values
	TicToc::tic(); std::cout << "    Computing the level-set solution values... " << std::flush;
	for (int i=1 ; i<=np ; i++) {
		// Get the coordinates of the current point
		double point[3];
		if ( MMGS_Get_vertex(mmgMesh, &(point[0]), &(point[1]), &(point[2]), NULL, NULL, NULL) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_VERTEX, nullptr);
		Point p(point[0], point[1], point[2]);

		// Compute the level-set solution value at the current coordinate
		double isovalue = Infill::TPMS_function(p, lt_type, lt_size, lt_feature);
		if ( MMGS_Set_scalarSol(mmgLs, isovalue, i) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL, nullptr);
	}

	// Check if the number of given entities match with mesh size
	if ( MMGS_Chk_meshData(mmgMesh, mmgLs) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_MESH_DATA_CHECK, nullptr);

	TicToc::toc("completed in ");

	/* ------------------------------- Step 4 ------------------------------- */
	// Discretize the level-set describing the lattice infill

	// Global meshing parameters
	if ( MMGS_Set_iparameter(mmgMesh, mmgLs, MMGS_IPARAM_iso, 1) != 1 ) // Set to level-set mode
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_IPARAMETER + "MMGS_IPARAM_iso", nullptr);

	if (me_settings.MMG_hausd > 0) {
		if ( MMGS_Set_dparameter(mmgMesh, mmgLs, MMGS_DPARAM_hausd, me_settings.MMG_hausd*lt_size.meanUnitCellSize) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMGS_DPARAM_hausd", nullptr);
	}

	if ( MMGS_Set_dparameter(mmgMesh, mmgLs, MMGS_DPARAM_angleDetection, me_settings.edgeProtectionAngle) != 1 ) 
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMGS_DPARAM_angleDetection", nullptr);

	if ( MMGS_Set_iparameter(mmgMesh, mmgLs, MMGS_IPARAM_verbose, me_settings.verbosity) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_IPARAMETER + "MMGS_IPARAM_verbose", nullptr);

	// Isovalue discretization
	TicToc::tic(); std::cout << "    Discretizing edges to protect... " << std::flush;
	ierr = MMGS_mmgsls(mmgMesh,mmgLs,NULL);
	TicToc::toc("completed in ");

	if ( ierr == MMG5_STRONGFAILURE )
		throw ExceptionError(MMG_ERRMSG::BAD_ENDING_OF_MMGSLIB, nullptr);
	else if ( ierr == MMG5_LOWFAILURE )
		std::cerr << "BAD ENDING OF MMGSLIB" << std::endl;

	// Save MMGS discretzed edges for DEBUG PURPOSES!!!
	if (me_settings.verbosity > 5) { 
		std::string tempFile = "MMGS_STEP_3.mesh";
		if ( MMGS_saveMesh(mmgMesh, tempFile.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_MESH, nullptr);
	}

	/* ------------------------------- Step 4 ------------------------------- */
	// Convert the MMGS output to a CGAL::Surface_mesh
	MMG5_int v0, v1, v2, ref, prop_id, isRidge;	
	int isRequired;

	if ( MMGS_Get_meshSize(mmgMesh, &np,  &nt, &na) != 1) // Get mesh info
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_MESH_SIZE, nullptr);

	// Vertices
	for (size_t i=1 ; i<=np ; i++) {
		double x,y,z;
		if (MMGS_Get_vertex(mmgMesh, &x, &y, &z, NULL, NULL, NULL) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_VERTEX, nullptr);
		surface_mesh.add_vertex(Point_3(x, y, z));
	}

	// Edges (and add their index to the ecmap so that they can later be constrained)
	ecmap = surface_mesh.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
	for (size_t i=0 ; i<na ; i++) {
		if (MMGS_Get_edge(mmgMesh, &v0, &v1, &prop_id, &isRidge, &isRequired) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_EDGE, nullptr);

		SurfaceMesh::Vertex_index v_0(v0-1); // Numbering shifted back by 1 to match CGAL counting!
		SurfaceMesh::Vertex_index v_1(v1-1); // Numbering shifted back by 1 to match CGAL counting!
		SurfaceMesh::Halfedge_index he_idx = surface_mesh.add_edge(v_0, v_1);

		ecmap[SurfaceMesh::Edge_index(he_idx)] = true;
	}

	// Triangles
	for (size_t i=1 ; i<=nt ; i++) {
		if (MMGS_Get_triangle(mmgMesh, &v0, &v1, &v2, &ref, &isRequired) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_TRIANGLE, nullptr);

		SurfaceMesh::Vertex_index v_0(v0-1); // Numbering shifted back by 1 to match CGAL counting!
		SurfaceMesh::Vertex_index v_1(v1-1); // Numbering shifted back by 1 to match CGAL counting!
		SurfaceMesh::Vertex_index v_2(v2-1); // Numbering shifted back by 1 to match CGAL counting!
		surface_mesh.add_face(v_0, v_1, v_2);
	}

	// Cleanup
	MMGS_Free_all(MMG5_ARG_start,
	              MMG5_ARG_ppMesh, &mmgMesh,
	              MMG5_ARG_ppLs,   &mmgLs,
	              MMG5_ARG_ppMet,  &mmgMet,
	              MMG5_ARG_end);
}

void Mesh::c3t3_to_MMG3D(const F_C3t3 &c3t3, MMG5_pMesh mmgMesh) {
	/* Convert c3t3 to MMG3D_mesh
	 * Inputs:
	 *   c3t3    : CGAL 3D-mesh structure.
	 * Output:
	 *   mmgMesh : MMG 3D-mesh structure.
	 */

	std::unordered_map<F_C3t3::Vertex_handle, int> V; // Vertex indices

	// Manually set the mesh size: mesh, vertices, tetra, prisms, triangles, quads, edges
	int np = c3t3.triangulation().number_of_vertices();
	int nt = std::distance(c3t3.facets_in_complex_begin(), c3t3.facets_in_complex_end());
	int ne = std::distance(c3t3.cells_in_complex_begin(), c3t3.cells_in_complex_end());

	if ( MMG3D_Set_meshSize(mmgMesh, np, ne, 0, nt, 0, 0) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_MESH_SIZE, nullptr);

	// Vertices
	int inum = 1;
	int vertexIDX = 1;
	for (typename F_PTr::Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin(); vit != c3t3.triangulation().finite_vertices_end(); ++vit) {
		const F_C3t3::Triangulation::Point &p = c3t3.triangulation().point(vit);
		if ( MMG3D_Set_vertex(mmgMesh, CGAL::to_double(p.x()), CGAL::to_double(p.y()), 
		                               CGAL::to_double(p.z()), 0, vertexIDX++) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_VERTEX, nullptr);

		V[vit] = inum++;
	}

	// Tetrahedra
	int cellIDX = 1;
	for (F_C3t3::Cell_handle cit : c3t3.cells_in_complex()) {
		if ( MMG3D_Set_tetrahedron(mmgMesh, V[cit->vertex(0)], V[cit->vertex(1)], 
		                                    V[cit->vertex(2)], V[cit->vertex(3)],
		                                    1, cellIDX++) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_TETRA, nullptr);
	}

	// Triangles
	int facetIDX = 1;
	for (F_C3t3::Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit) {
		if ( MMG3D_Set_triangle(mmgMesh, V[fit->first->vertex((fit->second+1)%4)], 
		                                 V[fit->first->vertex((fit->second+2)%4)],
                                     V[fit->first->vertex((fit->second+3)%4)],
																		 2, facetIDX++) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_TRIANGLE, nullptr);
	}
}