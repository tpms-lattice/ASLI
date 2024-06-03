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

#include "Mesh.h"

/* MESHCGAL discretizes the provided closed surface with the requested infill
 * using tools from the Computational Geometry Algorithms Library (CGAL).
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

void MeshCGAL::implicit2volume(const polygonSoup &shell, const latticeType &lt_type,
                               const latticeSize &lt_size, const latticeFeature &lt_feature,
                               const meshSettings &me_settings) {
	/* Discretizes the volume of a shell with an implicitly defined infill.
	 * Inputs :
	 *  shell             : Structure containing the general geometry
	 *  lt_type           : Lattice type parameters
	 *  lt_size           : Lattice size parameters
	 *  lt_feature        : Lattice feature parameters
	 *  me_settings       : Mesh settings
	 * Outputs :
	 *  .stl file containing the surface triangulation
	 *  .mesh file containing the volume mesh (optional)
	 */

	// Mesh user settings
	FT me_facetSize = me_settings.CGAL_facetSize;
	FT me_facetDistance = me_settings.CGAL_facetDistance;
	FT me_facetAngle = me_settings.CGAL_facetAngle;

	FT me_edgeSize = me_settings.CGAL_edgeSize;
	FT me_minEdgeSize = me_settings.CGAL_minEdgeSize;

	FT me_cellSize = me_settings.elementSize;   
	FT me_cellRadiusEdgeRatio = me_settings.CGAL_cellRadiusEdgeRatio;

	const bool isUniform = (lt_type.type != "hybrid" && lt_size.size != 0 && lt_feature.feature_val != 0 
		&& (lt_type.side == "scaffold" || lt_type.side == "void" ));
	std::filesystem::path outputPath = me_settings.output;

	// Determine bounding box
	std::vector<double> BB = {HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
	for (size_t i=0; i<shell.points.size(); ++i) {
		BB[0] = std::min(BB[0], shell.points[i].x());
		BB[1] = std::min(BB[1], shell.points[i].y());
		BB[2] = std::min(BB[2], shell.points[i].z());
		BB[3] = std::max(BB[3], shell.points[i].x());
		BB[4] = std::max(BB[4], shell.points[i].y());
		BB[5] = std::max(BB[5], shell.points[i].z());
	}
	CGAL::Bbox_3 boundingBox(BB[0], BB[1], BB[2], BB[3], BB[4], BB[5]);

	const std::vector<double> BBoxSize = {boundingBox.xmax()-boundingBox.xmin(),
	                                      boundingBox.ymax()-boundingBox.ymin(),
	                                      boundingBox.zmax()-boundingBox.zmin()};

	const std::vector<double> BBoxCenter = {(boundingBox.xmin()+boundingBox.xmax())/2.0,
	                                        (boundingBox.ymin()+boundingBox.ymax())/2.0,
	                                        (boundingBox.zmin()+boundingBox.zmax())/2.0};

	// Implicit domain
	auto myImplicitFunction = [lt_type, lt_size, lt_feature, BBoxSize, BBoxCenter](const Point_3 &p) 
		{ return internal::infill(p, lt_type, lt_size, lt_feature, BBoxSize, BBoxCenter); }; // To be able to pass funciton on as a static member
	Implicit_domain my_implicit_domain = Implicit_domain::create_implicit_mesh_domain(myImplicitFunction,
		boundingBox, CGAL::parameters::relative_error_bound(me_settings.CGAL_relativeErrorBound));

	// Polyhedra domain
	SurfaceMesh surface_mesh;

	Polyhedron polygonSurface;
	std::vector<std::vector<Point_3>> protected_features;

	SurfaceMesh::Property_map<edge_descriptor, bool> ecmap;
	if ( me_settings.edgeProtectionAngle > 0 ) {
		std::cout << "  Computing edges to protect..." << std::endl;

		Mesh::extractEdges(shell, lt_type, lt_size, lt_feature, me_settings, surface_mesh, ecmap);
		CGAL::copy_face_graph(surface_mesh, polygonSurface);

		// Construct polyline with edges to protect
		typedef SurfaceMesh::Property_map<edge_descriptor, bool> ECMap;
		typedef SurfaceMesh::Property_map<vertex_descriptor, bool> VPMap;
		typedef boost::filtered_graph<SurfaceMesh, Is_constrained<ECMap>, Is_feature<VPMap>> FG;

		auto v_pmap = surface_mesh.add_property_map<vertex_descriptor, bool>("v:on_feature_curve", false).first;
		for(auto e : edges(surface_mesh)) {
			auto v1 = source(e, surface_mesh);
			auto v2 = target(e, surface_mesh);
			if (get(ecmap, e)) {
				put(v_pmap, v1, true);
				put(v_pmap, v2, true);
			}
		}

		Polyline_visitor visitor{protected_features, surface_mesh};
		Is_constrained<ECMap> is_constrained{ecmap};
		Is_feature<VPMap> is_feature{v_pmap};
		
		FG constrained_edges{surface_mesh, is_constrained, is_feature};
		CGAL::split_graph_into_polylines(constrained_edges, visitor);

	} else {
		CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(shell.points, shell.polygons, polygonSurface);
	}
	Polyhedron_domain my_polyhedron_domain(polygonSurface);

	// Hybrid domain
	H_domain domain(my_implicit_domain, my_polyhedron_domain, lt_type.side);
	if ( me_settings.edgeProtectionAngle > 0 ) // Insert edges in domain
		domain.add_features(protected_features.begin(), protected_features.end());

	// FOR DEBUG PURPOSES
	if (me_settings.verbosity > 5) { 
		internal::writePolylinesToFile("CGAL_STEP_2_POLYLINES.txt", protected_features);
		std::cout << "  Polylines written to: CGAL_STEP_2_POLYLINES.txt" << std::endl;
	}

	// Meshing criteria
	double feature_size;
	TPMS_dependent_minfeaturesize_field<H_domain> facetDistanceField;
	TPMS_dependent_minfeaturesize_field<H_domain> facetSizeField;
	TPMS_dependent_minfeaturesize_field<H_domain> edgeSizeField;
	TPMS_dependent_minfeaturesize_field<H_domain> cellSizeField;

	//FT scaledFacetDistance = me_facetDistance * me_cellSize;
	if (isUniform) {
		Point point(0, 0, 0);
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);
		if (lt_type.side == "scaffold")
			feature_size = localSize.wallSize;
		else if (lt_type.side == "void")
			feature_size = localSize.poreSize;
		else
			feature_size = std::min(localSize.wallSize, localSize.wallSize);
		
		if ( feature_size <= 0 )
			throw ExceptionError("Feature size is less or equal to zero", nullptr);

	} else {
		// Sizing field for facet distance
		facetDistanceField.lt_type = lt_type;
		facetDistanceField.lt_size = lt_size;
		facetDistanceField.lt_feature = lt_feature;
		facetDistanceField.threshold = me_settings.threshold;
		facetDistanceField.parameter = &me_facetDistance;

		// Sizing field for feature dependent facet size
		facetSizeField.lt_type = lt_type;
		facetSizeField.lt_size = lt_size;
		facetSizeField.lt_feature = lt_feature;
		facetSizeField.threshold = me_settings.threshold;
		facetSizeField.parameter = &me_facetSize;

		// Sizing field for feature dependent edge size
		edgeSizeField.lt_type = lt_type;
		edgeSizeField.lt_size = lt_size;
		edgeSizeField.lt_feature = lt_feature;
		edgeSizeField.threshold = me_settings.threshold;
		edgeSizeField.parameter = &me_edgeSize;

		// Sizing field for feature dependent mesh size
		cellSizeField.lt_type = lt_type;
		cellSizeField.lt_size = lt_size;
		cellSizeField.lt_feature = lt_feature;
		cellSizeField.threshold = me_settings.threshold;
		cellSizeField.parameter = &me_cellSize;
	}

	// 
	TicToc::tic();
	if (me_settings.isVolumeMesh == true) std::cout << "  Meshing volume" << std::flush;
	else std::cout << "  Meshing surface" << std::flush;
	if (me_settings.n_threads > 1) std::cout << " (using up to " << me_settings.n_threads << " threads)" << std::flush;
	std::cout << "... " << std::flush;

	H_C3t3 c3t3;

	if (me_settings.isVolumeMesh) {
		if (isUniform) {
			H_Edge_criteria hedge_criteria(me_edgeSize*feature_size, me_minEdgeSize*feature_size); // Edge: size, min size
			H_Facet_criteria hfacet_criteria(me_facetAngle, 0, me_facetDistance*feature_size); // Facet: angle, size, approximation
			H_Cell_criteria hcell_criteria(me_cellRadiusEdgeRatio, me_cellSize*feature_size); // Cell: radius-edge ratio, size
			H_Mesh_criteria criteria(hedge_criteria, hfacet_criteria, hcell_criteria);

			c3t3 = CGAL::make_mesh_3<H_C3t3>(domain, criteria);
		} else {
			H_Edge_criteria hedge_criteria(edgeSizeField, me_minEdgeSize); // Edge: size, min size
			H_Facet_criteria hfacet_criteria(me_facetAngle, 0, facetDistanceField); // Facet: angle, size, approximation
			H_Cell_criteria hcell_criteria(me_cellRadiusEdgeRatio, cellSizeField); // Cell: radius-edge ratio, size
			H_Mesh_criteria criteria(hedge_criteria, hfacet_criteria, hcell_criteria);

			c3t3 = CGAL::make_mesh_3<H_C3t3>(domain, criteria);
		}
	} else { // If surface mesh
		if (isUniform) {
			H_Edge_criteria hedge_criteria(me_edgeSize*feature_size, me_minEdgeSize*feature_size); // Edge: size, min size
			H_Facet_criteria hfacet_criteria(me_facetAngle, me_facetSize*feature_size, me_facetDistance*feature_size); // Facet: angle, size, approximation
			H_Cell_criteria hcell_criteria(me_cellRadiusEdgeRatio,0); // Cell: radius-edge ratio, size
			H_Mesh_criteria criteria(hedge_criteria, hfacet_criteria, hcell_criteria);

			c3t3 = CGAL::make_mesh_3<H_C3t3>(domain, criteria, CGAL::parameters::no_perturb().no_exude().no_lloyd().no_odt());
		} else {
			H_Edge_criteria hedge_criteria(edgeSizeField, me_minEdgeSize); // Edge: size, min size
			H_Facet_criteria hfacet_criteria(me_facetAngle, facetSizeField, facetDistanceField); // Facet: angle, size, approximation
			H_Cell_criteria hcell_criteria(me_cellRadiusEdgeRatio,0); // Cell: radius-edge ratio, size
			H_Mesh_criteria criteria(hedge_criteria, hfacet_criteria, hcell_criteria);

			c3t3 = CGAL::make_mesh_3<H_C3t3>(domain, criteria, CGAL::parameters::no_perturb().no_exude().no_lloyd().no_odt());
		}
	}

	//std::cout << " Removing facets and vertices not in complex... " << std::flush;
	// Remove facets not beloging to the selected side
	if ( lt_type.side=="scaffold" || lt_type.side=="void" )
		MeshCGAL::internal::remove_from_complex(1, c3t3);

	c3t3.remove_isolated_vertices();

	TicToc::toc("completed in ");

	// Save output for DEBUG PURPOSES!!! //TEMP!!!
	if (me_settings.verbosity > 5) { dump_c3t3(c3t3, "CGAL_STEP_1"); }

	// Extract surface mesh of selected domain
	SurfaceMesh surface_scaffold;
	CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, surface_scaffold);

	if (lt_type.side == "scaffold")
		CGAL::Polygon_mesh_processing::keep_largest_connected_components(surface_scaffold,1);
	else if (lt_type.side == "void")
		CGAL::Polygon_mesh_processing::keep_largest_connected_components(surface_scaffold,2);

	if ( lt_type.side == "scaffold" || lt_type.side == "void" ) {
		if (CGAL::is_closed(surface_scaffold))
			CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(surface_scaffold);
		else
			std::cout << "WARNING: Surface is open" << std::endl;
	}
	// Mesh info
	std::cout << "\n    Surface triangulation: " << std::endl;
	std::cout << "      Number of vertices: " << surface_scaffold.number_of_vertices() << std::endl;
	std::cout << "      Number of facets: " << surface_scaffold.number_of_faces() << std::endl;
	if (me_settings.isVolumeMesh == TRUE) {
		std::cout << "\n    Volume mesh: " << std::endl;
		std::cout << "      Number of vertices: " << c3t3.triangulation().number_of_vertices() << std::endl;
		std::cout << "      Number of facets: " << 2 * c3t3.number_of_facets() << std::endl;
		std::cout << "      Number of tetrahedra: " << c3t3.number_of_cells() << std::endl;
	}

	/* ------------------------------- Step 3 ------------------------------- */
	// Store output and cleanup
	std::cout << "\nSAVING OUTPUT(S)... " << std::endl;
	std::cout << "  Output file(s): " << std::flush;

	// Save volume mesh
	if (me_settings.isVolumeMesh) {
		std::ofstream medit_file(outputPath.replace_extension(".mesh"));
		c3t3.output_to_medit(medit_file);
		medit_file.close();		

		std::cout << outputPath.filename() << " & "  << std::flush;
	}
	
	// Save surface mesh
	std::ofstream stl_fileOut(outputPath.replace_extension(".stl"), std::ios::out | std::ios::binary);
	if ( me_settings.STLFormat == "ASCII" ) {
		CGAL::set_mode(stl_fileOut, CGAL::IO::ASCII);
		stl_fileOut.precision(std::numeric_limits<double>::digits10 + 2);
	} else {
		CGAL::set_mode(stl_fileOut, CGAL::IO::BINARY);
	}
	CGAL::IO::write_STL(stl_fileOut, surface_scaffold);
	stl_fileOut.close();

	std::cout << outputPath.filename() << std::endl;

	c3t3.clear();
}

void MeshCGAL::implicit2volume_old(const polygonSoup &shell, const latticeType &lt_type,
                                   const latticeSize &lt_size, const latticeFeature &lt_feature,
                                   const meshSettings &me_settings) {
	/* Discretizes the volume of a shell with an implicitly defined infill.
	 * Inputs :
	 *  shell             : Structure containing the general geometry
	 *  lt_type           : Lattice type parameters
	 *  lt_size           : Lattice size parameters
	 *  lt_feature        : Lattice feature parameters
	 *  me_settings       : Mesh settings
	 * Outputs :
	 *  .stl file containing the surface triangulation
	 *  .mesh file containing the volume mesh (optional)
	 */

	// Mesh user settings
	FT me_facetSize = me_settings.CGAL_facetSize;
	FT me_facetDistance = me_settings.CGAL_facetDistance;
	FT me_facetAngle = me_settings.CGAL_facetAngle;

	FT me_edgeSize = me_settings.CGAL_edgeSize;

	FT me_cellSize = me_settings.elementSize;   
	FT me_cellRadiusEdgeRatio = me_settings.CGAL_cellRadiusEdgeRatio;

	std::filesystem::path outputPath = me_settings.output;

	double offset = 0.0;
	if ( me_settings.edgeProtectionAngle > 0 )
		offset = me_settings.CGAL_poissonOffset*lt_size.maxUnitCellSize;

	// Build a Poisson surface representation of the input geometry
	TicToc::tic(); std::cout << "  Creating implicit surface description of the outer shell... " << std::flush;
	Poisson_reconstruction_function *poissonSurfaceReconstruction = nullptr;
	std::vector<double> poissonBounds;

	implicitPoissonReconstruction(shell, offset, poissonBounds, poissonSurfaceReconstruction);
  TicToc::toc("completed in ");

	// Determine bounding box (and add some padding)
	CGAL::Bbox_3 boundingBox(poissonBounds[0] - 0.1*(poissonBounds[1]-poissonBounds[0]),
	                         poissonBounds[2] - 0.1*(poissonBounds[3]-poissonBounds[2]), 
	                         poissonBounds[4] - 0.1*(poissonBounds[5]-poissonBounds[4]), 
	                         poissonBounds[1] + 0.1*(poissonBounds[1]-poissonBounds[0]), 
	                         poissonBounds[3] + 0.1*(poissonBounds[3]-poissonBounds[2]), 
	                         poissonBounds[5] + 0.1*(poissonBounds[5]-poissonBounds[4]));

	// Create mesh domain
	auto myImplicitFunction = [poissonSurfaceReconstruction, lt_type, lt_size, lt_feature](const Point_3 &p) 
		{ return internal::signedDistance_old(p, lt_type, lt_size, lt_feature, poissonSurfaceReconstruction); }; // To be able to pass funciton on as a static member
	F_Implicit_domain domain = F_Implicit_domain::create_implicit_mesh_domain(myImplicitFunction,
		boundingBox, CGAL::parameters::relative_error_bound(me_settings.CGAL_relativeErrorBound));

	// Define meshing criteria
	double feature_size;
	TPMS_dependent_minfeaturesize_field<F_Implicit_domain> facetDistanceField;
	TPMS_dependent_minfeaturesize_field<F_Implicit_domain> facetSizeField;
	TPMS_dependent_minfeaturesize_field<F_Implicit_domain> edgeSizeField;
	TPMS_dependent_minfeaturesize_field<F_Implicit_domain> cellSizeField;

	FT scaledFacetDistance = me_facetDistance * me_cellSize;
	if ( lt_type.type == "hybrid" || lt_size.size == 0 || lt_feature.feature_val == 0) {
		// Sizing field for facet distance
		facetDistanceField.lt_type = lt_type;
		facetDistanceField.lt_size = lt_size;
		facetDistanceField.lt_feature = lt_feature;
		facetDistanceField.threshold = me_settings.threshold;
		facetDistanceField.parameter = &scaledFacetDistance;

		// Sizing field for feature dependent facet size
		facetSizeField.lt_type = lt_type;
		facetSizeField.lt_size = lt_size;
		facetSizeField.lt_feature = lt_feature;
		facetSizeField.threshold = me_settings.threshold;
		facetSizeField.parameter = &me_facetSize;

		// Sizing field for feature dependent edge size
		edgeSizeField.lt_type = lt_type;
		edgeSizeField.lt_size = lt_size;
		edgeSizeField.lt_feature = lt_feature;
		edgeSizeField.threshold = me_settings.threshold;
		edgeSizeField.parameter = &me_edgeSize;

		// Sizing field for feature dependent mesh size
		cellSizeField.lt_type = lt_type;
		cellSizeField.lt_size = lt_size;
		cellSizeField.lt_feature = lt_feature;
		cellSizeField.threshold = me_settings.threshold;
		cellSizeField.parameter = &me_cellSize;

	} else {
		Point point(0, 0, 0);
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);
		if ( lt_type.side == "void" )
			feature_size = localSize.poreSize;
		else
			feature_size = localSize.wallSize;
	}

	//
	TicToc::tic();
	if (me_settings.isVolumeMesh == true) std::cout << "  Meshing volume" << std::flush;
	else std::cout << "  Meshing surface" << std::flush;
	if (me_settings.n_threads > 1) std::cout << " (using up to " << me_settings.n_threads << " threads)" << std::flush;
	std::cout << "... " << std::flush;

	F_I_C3t3 c3t3;
	if (me_settings.isVolumeMesh == TRUE && me_settings.edgeProtectionAngle <= 0) {
		if ( lt_type.type == "hybrid" || lt_size.size == 0 || lt_feature.feature_val == 0) {
			F_I_Mesh_criteria criteria(CGAL::parameters::facet_angle = me_facetAngle, // Min triangle angle (degrees)
			                       CGAL::parameters::facet_distance = facetDistanceField, // Surface approximation error
			                       CGAL::parameters::cell_radius_edge_ratio = me_cellRadiusEdgeRatio, // Mesh tetrahedra radius|edge ratio upper bound
			                       CGAL::parameters::cell_size = cellSizeField);

			c3t3 = CGAL::make_mesh_3<F_I_C3t3>(domain, criteria);
		} else {
			F_I_Mesh_criteria criteria(CGAL::parameters::facet_angle = me_facetAngle, // Min triangle angle (degrees)
			                       CGAL::parameters::facet_distance = me_facetDistance*me_cellSize*feature_size, // Surface approximation error
			                       CGAL::parameters::cell_radius_edge_ratio = me_cellRadiusEdgeRatio, // Mesh tetrahedra radius|edge ratio upper bound
			                       CGAL::parameters::cell_size = me_cellSize*feature_size);

			c3t3 = CGAL::make_mesh_3<F_I_C3t3>(domain, criteria);
		}
		
	} else { // If surface mesh
		if ( lt_type.type == "hybrid" || lt_size.size == 0 || lt_feature.feature_val == 0) {
			F_I_Mesh_criteria criteria(CGAL::parameters::facet_angle = me_facetAngle, // Min triangle angle (degrees)
			                       CGAL::parameters::facet_size = facetSizeField, //me_facetSize*lt_size.maxUnitCellSize, // Max triangle size
			                       CGAL::parameters::facet_distance = facetDistanceField, // Surface approximation error
			                       CGAL::parameters::cell_radius_edge_ratio = me_cellRadiusEdgeRatio); // Mesh tetrahedra radius|edge ratio upper bound

			c3t3 = CGAL::make_mesh_3<F_I_C3t3>(domain, criteria, CGAL::parameters::no_perturb().no_exude().no_lloyd().no_odt());
		} else {
			F_I_Mesh_criteria criteria(CGAL::parameters::facet_angle = me_facetAngle, // Min triangle angle (degrees)
			                       CGAL::parameters::facet_size = me_facetSize*feature_size, //me_facetSize*lt_size.maxUnitCellSize, // Max triangle size
			                       CGAL::parameters::facet_distance = me_facetDistance*me_cellSize*feature_size, // Surface approximation error
			                       CGAL::parameters::cell_radius_edge_ratio = me_cellRadiusEdgeRatio); // Mesh tetrahedra radius|edge ratio upper bound

			c3t3 = CGAL::make_mesh_3<F_I_C3t3>(domain, criteria, CGAL::parameters::no_perturb().no_exude().no_lloyd().no_odt());
		}
	}
	c3t3.remove_isolated_vertices();

	TicToc::toc("completed in ");

	// Save initial scaffold mesh for DEBUG PURPOSES!!!
	if ( me_settings.verbosity > 5 ) {
		std::ofstream medit_file("CGAL_STEP_1.mesh");
		c3t3.output_to_medit(medit_file);
		medit_file.close();
	}

	// Postprocessing
	TicToc::tic(); std::cout << "  Post-processing mesh... " << std::endl;

	SurfaceMesh surface_scaffold;
	if ( me_settings.edgeProtectionAngle > 0 ) {
		// Extract surface triangulation
		SurfaceMesh scaffold;
		CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, scaffold);

		// Save extracted surface for DEBUG PURPOSES!!!
		if ( me_settings.verbosity > 5 )
			std::ofstream("CGAL_STEP_2.off") << std::setprecision(17) << scaffold;

		// Remove all but the largest component
		if ( lt_type.side == "void" )
			CGAL::Polygon_mesh_processing::keep_largest_connected_components(scaffold,2);
		else
			CGAL::Polygon_mesh_processing::keep_largest_connected_components(scaffold,1);


		// Check if scaffold surface is closed and try to fix if not
		if (!CGAL::is_closed(scaffold)) {
			CGAL::Polygon_mesh_processing::experimental::snap_borders(scaffold);
			CGAL::Polygon_mesh_processing::stitch_borders(scaffold);
			if (!CGAL::is_closed(scaffold)) {
				std::cerr << "WARNING: Open scaffold." << std::endl;
				dump_c3t3(c3t3, "dump_shell");
				goto skip_to_here;
			}
		}
		CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(scaffold);

		// Clip scaffold (using intersection) and collect constrained edges in a property map
		SurfaceMesh outer_shell;
		SurfaceMesh::Property_map<edge_descriptor, bool> ecmap_shell = outer_shell.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
		surfaceRemesh(shell, lt_type, lt_size, lt_feature, me_settings, outer_shell); // Remesh outer shell to avoid assertion violations due to overly large elements when computing the intersection
		if (me_settings.edgeProtectionAngle > 0) CGAL::Polygon_mesh_processing::detect_sharp_edges(outer_shell, me_settings.edgeProtectionAngle, ecmap_shell);

		// Check if remeshed outer shell is closed and try to fix if not
		if (!CGAL::is_closed(outer_shell)) {
			CGAL::Polygon_mesh_processing::experimental::snap_borders(outer_shell);
			CGAL::Polygon_mesh_processing::stitch_borders(outer_shell);
			if (!CGAL::is_closed(outer_shell)) {
				std::cerr << "WARNING: Open geometry." << std::endl;
				dump_c3t3(c3t3, "dump_shell");
				goto skip_to_here;
			}
		}
		CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(outer_shell);

		SurfaceMesh processedScaffold;
		SurfaceMesh::Property_map<edge_descriptor, bool> ecmap = processedScaffold.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
		try {
			if(CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(scaffold, outer_shell, processedScaffold, CGAL::parameters::default_values(), 
				CGAL::parameters::edge_is_constrained_map(ecmap_shell), CGAL::parameters::edge_is_constrained_map(ecmap)) == false) {
				std::cerr << "\nWARNING: corefine_and_compute_intersection output may have only been corefined.\n" << std::endl;
				if (!CGAL::is_closed(processedScaffold)) {
					if (me_settings.verbosity > 5) std::ofstream("CGAL_STEP_3.off") << std::setprecision(17) << processedScaffold;
					std::cerr << "ERROR: Surface is not closed. Remaing steps will be skipped." << std::endl;
					dump_c3t3(c3t3, "dump_processed_scaffold");
					goto skip_to_here;
				}
			}
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			dump_c3t3(c3t3, "dump_processed_scaffold");
			goto skip_to_here;
		}

		// Save cut surface for DEBUG PURPOSES!!!
		if (me_settings.verbosity > 5) std::ofstream("CGAL_STEP_3.off") << std::setprecision(17) << processedScaffold;

		// Compute/Recompute volume mesh (when surface was meshed and edges are to be preserved)
		if (me_settings.isVolumeMesh == true) {
			F_Polyhedron polygonSurface;
			CGAL::copy_face_graph(processedScaffold, polygonSurface);

			std::vector<std::vector<Point_3>> polylines;
			internal::ecmap2polylines(ecmap, processedScaffold, polylines);

			F_C3t3 c3t3_with_edges;
			polehedral2volume(polygonSurface, polylines, lt_type, lt_size, lt_feature, me_settings, c3t3_with_edges);

			// Write output to .mesh file
			std::ofstream medit_file(outputPath.replace_extension(".mesh"));
			c3t3_with_edges.output_to_medit(medit_file);
			medit_file.close();

			std::cout << outputPath.filename() << " & "  << std::flush;

			// Extract surface mesh of selected domain
			CGAL::facets_in_complex_3_to_triangle_mesh(c3t3_with_edges, surface_scaffold);
		}

	} else {
		// Write output to .mesh file
		if (me_settings.isVolumeMesh == true) {
			std::ofstream medit_file(outputPath.replace_extension(".mesh"));
			c3t3.output_to_medit(medit_file);
			medit_file.close();		

			std::cout << outputPath.filename() << " & "  << std::flush;
		}

		// Extract surface mesh of selected domain
		CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, surface_scaffold);
	}

	skip_to_here:
	// UNCOMMENT FOR DEBUGING...
	if ( me_settings.verbosity > 5 )
		std::ofstream("CGAL_STEP_4.off") << std::setprecision(17) << surface_scaffold;

	if ( lt_type.side=="void" )
		CGAL::Polygon_mesh_processing::keep_largest_connected_components(surface_scaffold,2);
	else
		CGAL::Polygon_mesh_processing::keep_largest_connected_components(surface_scaffold,1);

	// Reorient faces to bound a volume
	if (CGAL::is_closed(surface_scaffold) == true)
		CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(surface_scaffold);
	else
		std::cout << "WARNING: Surface is not closed!" << std::endl;

	// Clean up vertices
	TicToc::toc("completed in ");

	// Mesh info
	std::cout << "\n    Surface triangulation: " << std::endl;
	std::cout << "      Number of vertices: " << surface_scaffold.number_of_vertices() << std::endl;
	std::cout << "      Number of facets: " << surface_scaffold.number_of_faces() << std::endl;
	if (me_settings.isVolumeMesh == TRUE) {
		std::cout << "\n    Volume mesh: " << std::endl;
		std::cout << "      Number of vertices: " << c3t3.triangulation().number_of_vertices() << std::endl;
		std::cout << "      Number of facets: " << 2 * c3t3.number_of_facets() << std::endl;
		std::cout << "      Number of tetrahedra: " << c3t3.number_of_cells() << std::endl;
	}

	/* ------------------------------- Step 3 ------------------------------- */
	// Store output and cleanup
	std::cout << "\nSAVING OUTPUT(S)... " << std::endl;
	std::cout << "  Output file(s): " << std::flush;

	// Save surface mesh
	std::ofstream stl_fileOut(outputPath.replace_extension(".stl"), std::ios::out | std::ios::binary);
	if ( me_settings.STLFormat == "ASCII" ) {
		CGAL::set_mode(stl_fileOut, CGAL::IO::ASCII);
		stl_fileOut.precision(std::numeric_limits<double>::digits10 + 2);
	} else {
		CGAL::set_mode(stl_fileOut, CGAL::IO::BINARY);
	}
	CGAL::IO::write_STL(stl_fileOut, surface_scaffold);
	stl_fileOut.close();

	std::cout << outputPath.filename() << std::endl;

	// Cleanup
	delete poissonSurfaceReconstruction;
	poissonSurfaceReconstruction = nullptr;
	c3t3.clear();
}

void MeshCGAL::polehedral2volume (const F_Polyhedron &polygonSurface, 
                                  const std::vector<std::vector<Point_3>> &protected_features, 
                                  const latticeType &lt_type, const latticeSize &lt_size, 
                                  const latticeFeature &lt_feature, const meshSettings &me_settings,
                                  F_C3t3 &c3t3) {
	/* Discretizes the volume of a closed triangulation.
	 * Inputs:
	 *  surfaceMesh    : Surface mesh
	 *  ecmap          : Edge is constrained map
	 *  lt_type        : Lattice type data
	 *  lt_size        : Lattice size data
	 *  lt_feature     : Lattice feature data
	 *  me_settings    : Mesh settings
	 * Outputs:
	 *  c3t3           : F_C3t3 containing the triangulated volume
	 */

	// Mesh user settings
	FT me_facetSize = me_settings.CGAL_facetSize;
	FT me_facetDistance = me_settings.CGAL_facetDistance;
	FT me_facetAngle = me_settings.CGAL_facetAngle;

	FT me_edgeSize = me_settings.CGAL_edgeSize;

	FT me_cellSize = me_settings.elementSize;   
	FT me_cellRadiusEdgeRatio = me_settings.CGAL_cellRadiusEdgeRatio;

	TicToc::tic(); std::cout << "  Meshing volume (using up to " << me_settings.n_threads << " threads)... " << std::flush;

	// Create domain
	F_Polyhedron_domain domain(polygonSurface);
	domain.add_features(protected_features.begin(), protected_features.end());

	// Define meshing criteria
	double feature_size;
	TPMS_dependent_minfeaturesize_field<F_Polyhedron_domain> facetDistanceField;
	TPMS_dependent_minfeaturesize_field<F_Polyhedron_domain> facetSizeField;
	TPMS_dependent_minfeaturesize_field<F_Polyhedron_domain> edgeSizeField;
	TPMS_dependent_minfeaturesize_field<F_Polyhedron_domain> cellSizeField;

	FT scaledFacetDistance = me_facetDistance * me_cellSize;
	if ( lt_type.type == "hybrid" || lt_size.size == 0 || lt_feature.feature_val == 0) {
		// Sizing field for facet distance
		facetDistanceField.lt_type = lt_type;
		facetDistanceField.lt_size = lt_size;
		facetDistanceField.lt_feature = lt_feature;
		facetDistanceField.threshold = me_settings.threshold;
		facetDistanceField.parameter = &scaledFacetDistance;

		// Sizing field for feature dependent facet size
		facetSizeField.lt_type = lt_type;
		facetSizeField.lt_size = lt_size;
		facetSizeField.lt_feature = lt_feature;
		facetSizeField.threshold = me_settings.threshold;
		facetSizeField.parameter = &me_facetSize;

		// Sizing field for feature dependent edge size
		edgeSizeField.lt_type = lt_type;
		edgeSizeField.lt_size = lt_size;
		edgeSizeField.lt_feature = lt_feature;
		edgeSizeField.threshold = me_settings.threshold;
		edgeSizeField.parameter = &me_edgeSize;

		// Sizing field for feature dependent mesh size
		cellSizeField.lt_type = lt_type;
		cellSizeField.lt_size = lt_size;
		cellSizeField.lt_feature = lt_feature;
		cellSizeField.threshold = me_settings.threshold;
		cellSizeField.parameter = &me_cellSize;

	} else {
		Point point(0, 0, 0);
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);
		if ( lt_type.side == "scaffold" )
			feature_size = localSize.wallSize;
		else if ( lt_type.side == "void" )
			feature_size = localSize.poreSize;
		else
			feature_size = std::min(localSize.wallSize, localSize.wallSize);
	}


	if (lt_type.type == "hybrid" && lt_size.size == 0 && lt_feature.feature_val == 0) { // If lattice is uniform
		// Mesh criteria
		F_Mesh_criteria criteria(CGAL::parameters::edge_size = edgeSizeField,
		                         CGAL::parameters::facet_angle = me_facetAngle, // Min triangle angle (degrees)
		                         //CGAL::parameters::facet_size = facetSizeField, // Max triangle size
		                         CGAL::parameters::facet_distance = facetDistanceField, // Surface approximation error
		                         CGAL::parameters::cell_radius_edge_ratio = me_cellRadiusEdgeRatio, // Mesh tetrahedra radius|edge ratio upper bound
		                         CGAL::parameters::cell_size = cellSizeField); // Note: Sizing field not yet suported by CHAL so fixed limit for now ...

		// Mesh generation
		c3t3 = CGAL::make_mesh_3<F_C3t3>(domain, criteria);//, 
		                                 //CGAL::parameters::manifold()); // Causes assertion violation errors!
		
	} else {
		Point point(0, 0, 0);
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);

		// Mesh criteria
		F_Mesh_criteria criteria(CGAL::parameters::edge_size = me_edgeSize*feature_size,
		                         CGAL::parameters::facet_angle = me_facetAngle, // Min triangle angle (degrees)
		                         CGAL::parameters::facet_size = me_facetSize * feature_size, // Max triangle size
		                         CGAL::parameters::facet_distance = me_facetDistance*me_cellSize*feature_size, // Surface approximation error
		                         CGAL::parameters::cell_radius_edge_ratio = me_cellRadiusEdgeRatio, // Mesh tetrahedra radius|edge ratio upper bound
		                         CGAL::parameters::cell_size = me_cellSize*feature_size, // 
		                         CGAL::parameters::edge_min_size = me_cellSize * feature_size/5); //

		// Mesh generation
		c3t3 = CGAL::make_mesh_3<F_C3t3>(domain, criteria);//, 
		                                 //CGAL::parameters::manifold()); // Causes assertion violation errors!
	}
	
	// Clean up vertices
	c3t3.remove_isolated_vertices();

	std::cout << "Finished!" << std::endl; TicToc::toc();

	std::cout << "\n  VOLUME MESH: " << std::endl;
	std::cout << "  Number of vertices: " << c3t3.triangulation().number_of_vertices() << std::endl;
	std::cout << "  Number of facets: " << 2 * c3t3.number_of_facets() << std::endl;
	std::cout << "  Number of tetrahedra: " << c3t3.number_of_cells() << std::endl;
}

void MeshCGAL::surfaceRemesh(const polygonSoup &shell, const latticeType &lt_type, const latticeSize &lt_size,
	const latticeFeature &lt_feature, const meshSettings &me_settings, SurfaceMesh &surfaceMeshOut) {
	/* 
	 * Inputs :
	 *  
	 * Outputs :
	 *  
	 */

	// Mesh user settings
	FT me_facetSize = 0.42;
	FT me_facetAngle = me_settings.CGAL_facetAngle;
	FT me_facet_distance = 0.001;

	FT edgeProtectionAngle = me_settings.edgeProtectionAngle;

	#ifdef ASLI_VERBOSE
		TicToc::tic(); std::cout << "  Remeshing bounding geometry (using up to " << me_settings.n_threads << " threads)... " << std::flush;
	#endif

	F_Polyhedron polygonSurface;
	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(shell.points, shell.polygons, polygonSurface);

	// Create a polyhedral domain with one polyhedron (and no "bounding polyhedron") such that the volumetric part of the domain stays empty
	std::vector<F_Polyhedron*> poly_ptrs_vector(1, &polygonSurface);
	F_Polyhedron_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());

	// Get sharp features
	if ( edgeProtectionAngle > 0 )
		domain.detect_features(edgeProtectionAngle);

	F_C3t3 c3t3;
	if (lt_type.type != "hybrid" && lt_size.minUnitCellSize == lt_size.maxUnitCellSize &&  lt_feature.feature_val > 0) { // If lattice is uniform
		Point point(0, 0, 0);
		featureSize localSize = Infill::featureSize_function(point, lt_type, lt_size, lt_feature);

		// Mesh criteria
		F_Mesh_criteria criteria(CGAL::parameters::edge_size = me_facetSize * localSize.wallSize,
		                         CGAL::parameters::facet_angle = me_facetAngle,
		                         CGAL::parameters::facet_size = me_facetSize * localSize.wallSize,
		                         CGAL::parameters::facet_distance = me_facet_distance);

		// Mesh generation
		c3t3 = CGAL::make_mesh_3<F_C3t3>(domain, criteria,
		                                         CGAL::parameters::no_exude(),
		                                         CGAL::parameters::no_perturb(),
		                                         CGAL::parameters::no_lloyd(),
		                                         CGAL::parameters::no_odt(),
		                                         CGAL::parameters::manifold());
	} else {
	// Define sizing field
	TPMS_dependent_wallsize_field<F_Polyhedron_domain> facetSize;
	facetSize.lt_type = lt_type;
	facetSize.lt_size = lt_size;
	facetSize.lt_feature = lt_feature;
	facetSize.threshold = me_settings.threshold;
	facetSize.parameter = &me_facetSize;

	// Mesh criteria
	F_Mesh_criteria criteria(CGAL::parameters::edge_size = facetSize,
	                         CGAL::parameters::facet_angle = me_facetAngle,
	                         CGAL::parameters::facet_size = facetSize,
	                         CGAL::parameters::facet_distance = me_facet_distance);

	// Mesh generation
	c3t3 = CGAL::make_mesh_3<F_C3t3>(domain, criteria,
	                                 CGAL::parameters::no_exude(),
	                                 CGAL::parameters::no_perturb(),
	                                 CGAL::parameters::no_lloyd(),
	                                 CGAL::parameters::no_odt(),
	                                 CGAL::parameters::manifold());
	}

	// Extract surface triangulation
	CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, surfaceMeshOut);

	// Reorient all face normals to point outward
	if (CGAL::is_closed(surfaceMeshOut) == true) {
		CGAL::Polygon_mesh_processing::orient(surfaceMeshOut, CGAL::parameters::outward_orientation(true));
	} else {
		std::cout << "\n  WARNING: Polygon has border edges!" << std::endl;
	}

	#ifdef ASLI_VERBOSE
		std::cout << "Finished!" << std::endl; TicToc::toc();
	#endif
}

void MeshCGAL::implicitPoissonReconstruction(const polygonSoup &shell, const double &offset, 
	std::vector<double> &boundingBox, Poisson_reconstruction_function *&poissonReconstruction) {
	/* Creates an implicit surface reconstruction using the Poisson method.
	 * Inputs:
	 *   shell  : Polygon soup
	 *   offset : Offset to give to reconstructed surface
   * Output:
	 *   boundingBox           : Bounding box of poisson surface reconstruction
	 *   poissonReconstruction : Poisson surface reconstruction
	 */

	const int nNeighbors = 18; // Number of neigbours for normal estimation

	PointList scaledPoints;
	boundingBox = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
//	bool useFaster = true; // TEMP!! Note: The faster offset is much faster but not completelly correct (although it usually works just fine). The slower one is significantlly slower but correct.
//	if (useFaster == 1) {
		// Convert points to Point_with_normal and scale
		PointList pointsCGAL; // scaled points
		for (size_t i=0; i<shell.points.size(); ++i) {
			Point_with_normal p(Point_3(shell.points[i].x(), shell.points[i].y(), shell.points[i].z()), 
			                    Vector_3(0.0, 0.0, 0.0));
			pointsCGAL.push_back(p);
		}

	//	// Check if points bound a volume
	//	if (CGAL::is_closed(pointsCGAL) == false)
	//		throw "ERROR: Input outer surface is not closed";

		// Estimate and orient normals
		CGAL::jet_estimate_normals<Concurrency_tag>(pointsCGAL, nNeighbors,
			CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
		                    normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

		// Remove points with an unoriented normal
		std::vector<Point_with_normal>::iterator unoriented_points_begin = 
			CGAL::mst_orient_normals(pointsCGAL, nNeighbors,
			CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
			                  normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

		pointsCGAL.erase(unoriented_points_begin, pointsCGAL.end());

		// Scaling
		for (size_t i=0; i<pointsCGAL.size(); ++i) {
			Point_with_normal p(Point_3(pointsCGAL[i].first.x() + offset*pointsCGAL[i].second.x(), 
			                            pointsCGAL[i].first.y() + offset*pointsCGAL[i].second.y(),
			                            pointsCGAL[i].first.z() + offset*pointsCGAL[i].second.z()), 
													Vector_3(pointsCGAL[i].second.x(), 
													         pointsCGAL[i].second.y(), 
													         pointsCGAL[i].second.z()));
			scaledPoints.push_back(p);

			// Compute bounding box of scaled point cloud
			boundingBox[0] = std::min(boundingBox[0], scaledPoints[i].first.x());
			boundingBox[1] = std::max(boundingBox[1], scaledPoints[i].first.x());
			boundingBox[2] = std::min(boundingBox[2], scaledPoints[i].first.y());
			boundingBox[3] = std::max(boundingBox[3], scaledPoints[i].first.y());
			boundingBox[4] = std::min(boundingBox[4], scaledPoints[i].first.z());
			boundingBox[5] = std::max(boundingBox[5], scaledPoints[i].first.z());
		}

//	} else {
//		TriangleMesh triangle_mesh;
//		CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(shell.points, shell.polygons, triangle_mesh);
//		const double averageSpacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(shell.points, nNeighbors);
//
////	// Check if points bound a volume
////	if (CGAL::is_closed(triangle_mesh) == false)
////		throw "ERROR: Input outer surface is not closed";
////	if (CGAL::Polygon_mesh_processing::does_bound_a_volume(triangle_mesh) == false)
////		throw "ERROR: Input outer surface does not bound a volume";
//
//		// Compute offset
//		polygonSoup offset_shell;
//		if (offset > 0) {
//			double facetAngle = 30;
//			double facetSize = averageSpacing/2;
//			double facetDistance = 0.1;
//			offset_shell = cgal_off_meshing(triangle_mesh, offset, facetAngle, facetSize, facetDistance);
//
//		} else {
//			offset_shell = shell;
//		}
//
//		// Convert to Point_with_normal
//		for (size_t i=0; i<offset_shell.points.size(); ++i) {
//			Point_with_normal p(Point_3(offset_shell.points[i].x(), offset_shell.points[i].y(), offset_shell.points[i].z()), 
//													Vector_3(0.0, 0.0, 0.0));
//			scaledPoints.push_back(p);
//
//			// Compute bounding box of scaled point cloud
//			boundingBox[0] = std::min(boundingBox[0], scaledPoints[i].first.x());
//			boundingBox[1] = std::max(boundingBox[1], scaledPoints[i].first.x());
//			boundingBox[2] = std::min(boundingBox[2], scaledPoints[i].first.y());
//			boundingBox[3] = std::max(boundingBox[3], scaledPoints[i].first.y());
//			boundingBox[4] = std::min(boundingBox[4], scaledPoints[i].first.z());
//			boundingBox[5] = std::max(boundingBox[5], scaledPoints[i].first.z());
//		}
//
//		// Estimate and orient normals
//		CGAL::jet_estimate_normals<Concurrency_tag>(scaledPoints, nNeighbors,
//																								CGAL::parameters::
//																									point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
//																									normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));
//
//		std::vector<Point_with_normal>::iterator unoriented_points_begin = 
//			CGAL::mst_orient_normals(scaledPoints, nNeighbors,
//															CGAL::parameters::
//																point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
//																normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));
//
//		// Remove points with an unoriented normal
//		scaledPoints.erase(unoriented_points_begin, scaledPoints.end());
//
//		// Scott Meyer's "swap trick" to trim excess capacity (Seems to solve the occasional 'CGAL ERROR: assertion violation! .../Refine_facets_3.h Line: 862')
//		std::vector<Point_with_normal>(scaledPoints).swap(scaledPoints);
//
////	{
////			// Save surface to .stl file (FOR DEBUG PURPOSES!)
////			TriangleMesh triangle_mesh_temp;
////			CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(offset_shell.points, offset_shell.polygons, triangle_mesh_temp);
////
////			std::ofstream stl_fileOut("out_test.stl", std::ios::out | std::ios::binary);
////			if ( me_settings.STLFormat == "ASCII" ) {
////				CGAL::set_mode(stl_fileOut, CGAL::IO::ASCII);
////				stl_fileOut.precision(std::numeric_limits<double>::digits10 + 2);
////			} else {
////				CGAL::set_mode(stl_fileOut, CGAL::IO::BINARY);
////			}
////			CGAL::IO::write_STL(stl_fileOut, triangle_mesh_temp);
////			stl_fileOut.close();
////	}
//	}

	// Create implicit function
	poissonReconstruction = new Poisson_reconstruction_function(scaledPoints.begin(), scaledPoints.end(),
	                                                            CGAL::First_of_pair_property_map<Point_with_normal>(),
	                                                            CGAL::Second_of_pair_property_map<Point_with_normal>());

	// Compute the Poisson indicator function f() at each vertex of the triangulation.
	if ( !poissonReconstruction->compute_implicit_function() )
		throw ExceptionError("ERROR_POISSON_LINEAR_SOLVER_FAILED", nullptr);
}

double MeshCGAL::internal::infill(const Point_3 &p, const latticeType &lt_type,
	const latticeSize &lt_size, const latticeFeature &lt_feature,
	const std::vector<double> &BBoxSize, const std::vector<double> &BBoxCenter) {
	/* 
	 * Inputs :
	 *  
	 * Returns :
	 *  
	 */

	Point point(p.x(), p.y(), p.z());

	const double x_abs = std::abs(p.x()-BBoxCenter[0])/BBoxSize[0];
	const double y_abs = std::abs(p.y()-BBoxCenter[1])/BBoxSize[1];
	const double z_abs = std::abs(p.z()-BBoxCenter[2])/BBoxSize[2];

	double a = std::max(std::max(x_abs, y_abs), z_abs) - 1.0/2.0; // BBox
	double b = Infill::TPMS_function (point, lt_type, lt_size, lt_feature); // Infill

	return (a > b) ? a : b;
}

double MeshCGAL::internal::signedDistance_old(const Point_3 &p, const latticeType &lt_type,
	const latticeSize &lt_size, const latticeFeature &lt_feature,
	const Poisson_reconstruction_function *const &poissonReconstruction) {
	/* 
	 * Inputs :
	 *  
	 * Returns :
	 *  
	 */

	Point point(p.x(), p.y(), p.z());

//	int mode = 0;
//	switch (mode) {
//		case 1: // Signed distance function of outer shape (For debug purposes)
//			return poissonReconstruction->operator()(p);
//
//		case 2: { // Signed distance function of a sphere (For debug purposes)
//			return BasicGeometries::implicit_function("spheroid", {0.5}, p);
//		}
//		case 3: { // Signed distance function of a lattice filled sphere (For debug purposes)
//			double a = BasicGeometries::implicit_function("spheroid", {0.5}, p);
//			double b = Infill::TPMS_function (point, lt_type, lt_size, lt_feature); // Infill
//				
//			return (a > b) ? a : b;
//
//		} default: { // Signed distance function of lattice filled outer shape
			double a = poissonReconstruction->operator()(p); // Outer shape
			double b = Infill::TPMS_function (point, lt_type, lt_size, lt_feature); // Infill

			 // Select a sub-domain to keep
			if (lt_type.side == "void")
				return (a > -b) ? a : -b;
			else
				return (a > b) ? a : b;

//		}
//	}
}

void MeshCGAL::internal::ecmap2polylines(const SurfaceMesh::Property_map<edge_descriptor, bool> &ecmap, 
                                         SurfaceMesh &surfaceMesh, std::vector<std::vector<Point_3>> &polylines) {
	/* 
	 * Inputs:
	 *  ecmap       : Edge constrained map
	 *  surfaceMesh : 
	 * Outputs:
	 *  surfaceMesh : 
	 *  polylines   : 
	 */

	// Construct polyline with edges to protect
	typedef SurfaceMesh::Property_map<edge_descriptor, bool> ECMap;
	typedef SurfaceMesh::Property_map<vertex_descriptor, bool> VPMap;
	typedef boost::filtered_graph<SurfaceMesh, Is_constrained<ECMap>, Is_feature<VPMap>> FG;

	std::vector<std::vector<Point_3>> protected_features;

	auto v_pmap = surfaceMesh.add_property_map<vertex_descriptor, bool>("v:on_feature_curve", false).first;
	for(auto e : edges(surfaceMesh)) {
		auto v1 = source(e, surfaceMesh);
		auto v2 = target(e, surfaceMesh);
		if (get(ecmap, e)) {
			put(v_pmap, v1, true);
			put(v_pmap, v2, true);
		}
	}

	Polyline_visitor visitor{protected_features, surfaceMesh};
	Is_constrained<ECMap> is_constrained{ecmap};
	Is_feature<VPMap> is_feature{v_pmap};
	
	FG constrained_edges{surfaceMesh, is_constrained, is_feature};
	CGAL::split_graph_into_polylines(constrained_edges, visitor);
}

template <typename C3T3>
void MeshCGAL::internal::remove_from_complex(const typename C3T3::Subdomain_index &sd_index, C3T3 &c3t3) {
	/* Removes the selected facet subdomain from the complex
	 * Inputs :
	 *  sd_index : Subdomain_index
	 *  c3t3     : Mesh_complex_3_in_triangulation_3
	 * Outputs :
	 *  Mesh with selected domain removed
	 */

	std::set<typename C3T3::Facet> facets_to_remove;
	for(typename C3T3::Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit) {
		typename C3T3::Subdomain_index cell_sd=c3t3.subdomain_index(fit->first);
		typename C3T3::Subdomain_index opp_sd=c3t3.subdomain_index(fit->first->neighbor(fit->second));

		if (cell_sd==sd_index || opp_sd==sd_index)
			continue;

		facets_to_remove.insert(*fit);
	}

	for (typename std::set<typename C3T3::Facet>::iterator it = facets_to_remove.begin(); it != facets_to_remove.end(); ++it)
		c3t3.remove_from_complex(*it);
}

void MeshCGAL::internal::writePolylinesToFile(const std::string &filename,
	const std::vector<std::vector<Point_3>> &polylines) {
	/* Writes polylines to a file
	 * Inputs :
	 *  filename  : 
	 *  polylines : 
	 * Outputs :
	 *   File containing polylines
	 */

	std::ofstream outfile(filename);
	if (!outfile.is_open())
		throw ExceptionError(CGAL_ERRMSG::FAILED_TO_OPEN_FILE + filename, nullptr);

	for (const auto& polyline : polylines) {
		outfile << polyline.size() << " ";
		for (const auto& point : polyline)
			outfile << point.x() << " " << point.y() << " " << point.z() << " ";
		outfile << std::endl;
	}
	outfile.close();
}