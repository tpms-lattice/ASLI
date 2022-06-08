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

#include "MeshCGAL.h"

/* MESHCGAL discretizes the provided closed surface with the requested infill
 * using tools from the Computational Geometry Algorithms Library (CGAL).
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

bool MeshCGAL::implicit2volume(outerShell &shell, latticeType lt_type, 
                               latticeSize lt_size, latticeFeature lt_feature,
                               meshSettings me_settings,
															 std::filesystem::path &outputFile_string) {
	/* Discretizes the volume of a shell with an implicitly defined infill.
	 * Inputs :
	 *  shell             : Structure containing the outer shell geometry
	 *  lt_type           : Lattice type data
	 *  lt_size           : Lattice size data
	 *  lt_feature        : Lattice feature data
	 *  me_settings       : Mesh settings
	 *  outputFile_string : Output file path
	 * Outputs :
	 *  .mesh file containing the volume mesh
	 *  .stl file containing the surface triangulation
	 */

	// Mesh user settings
	FT me_facetAngle = me_settings.CGAL_facetAngle;
	FT me_facetSize = me_settings.CGAL_facetSize;
	FT me_facetDistance = me_settings.CGAL_facetDistance;

	FT me_cellRadiusEdgeRatio = me_settings.CGAL_cellRadiusEdgeRatio;
	FT me_cellSize = me_settings.CGAL_cellSize;     

	std::string me_side = me_settings.side;

	// Concurrency settings
	int n_threads = 1;
	#ifdef CGAL_CONCURRENT_MESH_3
		if (me_settings.n_threads < 1) {
			me_settings.n_threads = 1;
		} else if (me_settings.n_threads > tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism)) {
			me_settings.n_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism); // tbb::this_task_arena::max_concurrency()
			std::cout << "  Warning: Requested number of threads exceeds the maximum allowed " 
			          << "parallelism. Number of threads set to the maximum allowed parallelism." << std::endl;
		}
		//tbb::task_scheduler_init init(me_settings.n_threads); // Deprecated
		//n_threads = tbb::this_task_arena::max_concurrency();
		
		tbb::global_control control(tbb::global_control::max_allowed_parallelism, me_settings.n_threads);
		n_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
	#endif

	// Compute implicit function of outer surface
	double offset = 0.0;
	if (me_settings.CGAL_preserveEdges == TRUE)
		offset = me_settings.CGAL_poissonOffset*lt_size.maxUnitCellSize;

	Poisson_reconstruction_function *poissonSurfaceReconstruction;
	poissonSurfaceReconstruction = NULL;
	std::vector<double> poissonBounds;

	implicitPoissonReconstruction(shell, offset, poissonBounds, &poissonSurfaceReconstruction);

	std::cout << "  Meshing " << std::flush;
	if (me_settings.CGAL_preserveEdges == FALSE && me_settings.volumeMesh == TRUE)
		std::cout << "volume " << std::flush;
	else
		std::cout << "surface " << std::flush;
	std::cout << "(using " << n_threads << " threads)... " << std::flush; TicToc::tic();

	// Determine bounding box (and add some padding)
	CGAL::Bbox_3 boundingBox(poissonBounds[0] - 0.1*(poissonBounds[1]-poissonBounds[0]),
	                         poissonBounds[2] - 0.1*(poissonBounds[3]-poissonBounds[2]), 
	                         poissonBounds[4] - 0.1*(poissonBounds[5]-poissonBounds[4]), 
	                         poissonBounds[1] + 0.1*(poissonBounds[1]-poissonBounds[0]), 
	                         poissonBounds[3] + 0.1*(poissonBounds[3]-poissonBounds[2]), 
	                         poissonBounds[5] + 0.1*(poissonBounds[5]-poissonBounds[4]));

	// Create mesh domain
	auto myImplicitFunction = [poissonSurfaceReconstruction, &lt_type, &lt_size, &lt_feature, me_side](const Point_3 &p) { return internal::signedDistance(p, &lt_type, &lt_size, &lt_feature, me_side, poissonSurfaceReconstruction); };// To be able to pass funciton on as a static member
	Implicit_domain domain = Implicit_domain::create_implicit_mesh_domain(myImplicitFunction,
	                                                                      boundingBox);

	// Sizing fields for facet size and distance
	TPMS_dependent_unitcellsize_field_I facetSize;
	facetSize.lt_size = &lt_size;
	facetSize.parameter = &me_facetSize;
	
	TPMS_dependent_unitcellsize_field_I facetDistance;
	facetDistance.lt_size = &lt_size;
	facetDistance.parameter = &me_facetDistance;

  // Meshing
	C3t3 c3t3;
	if (me_settings.CGAL_preserveEdges == TRUE) {
		// Meshing criteria (surface mesh)
		Mesh_criteria criteria(CGAL::parameters::facet_angle = me_facetAngle, // Min triangle angle (degrees)
		                       CGAL::parameters::facet_size = facetSize, //me_facetSize*lt_size.maxUnitCellSize, // Max triangle size
		                       CGAL::parameters::facet_distance = facetDistance);//me_facetDistance * lt_size.meanUnitCellSize); // Approximation error

		// Generate mesh
		//c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
		//                               CGAL::parameters::manifold());
		c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
		//                               CGAL::parameters::manifold(), // Leads to an assertion violation occasionally
		                               CGAL::parameters::no_exude(),
		                               CGAL::parameters::no_perturb(),
				                           CGAL::parameters::no_lloyd(),
		                               CGAL::parameters::no_odt());

	} else {
		// Sizing field for cell size
		TPMS_dependent_wallsize_field_I cellSize;
		cellSize.lt_type = &lt_type;
		cellSize.lt_size = &lt_size;
		cellSize.lt_feature = &lt_feature;
		cellSize.parameter = &me_cellSize;

		// Meshing criteria (volume mesh)
		Mesh_criteria criteria(CGAL::parameters::facet_angle = me_facetAngle,     // Min triangle angle (degrees)
		                       CGAL::parameters::facet_size = facetSize, //me_facetSize*lt_size.maxUnitCellSize,        // Max triangle size
		                       CGAL::parameters::facet_distance = facetDistance, // Surface approximation error
		                       CGAL::parameters::cell_radius_edge_ratio = me_cellRadiusEdgeRatio, // Mesh tetrahedra radius|edge ratio upper bound
		                       CGAL::parameters::cell_size = cellSize);

		// Generate mesh
		c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);//,
		//                               CGAL::parameters::manifold()); // Leads to an assertion violation occasionally
	}

	std::cout << "Finished!" << std::endl; TicToc::toc();

	// Determine current time to append to filenames
  time_t t;
  char currentTime[50];
	std::time(&t);
  std::strftime(currentTime, sizeof(currentTime), "_%Y-%m-%d_%H%M", localtime(&t)); 
	outputFile_string = me_settings.output + currentTime + "_cgal.stl";

//	// Write output to .mesh file (prior to post-processing FOR DEBUGING)
//	std::ofstream medit_file(me_settings.output + currentTime + "_cgal_DEBUG.mesh");
//	c3t3.output_to_medit(medit_file);
//	medit_file.close();
//
//	// Write output to .off file (prior to post-processing FOR DEBUGING)
//	std::ofstream off_file(me_settings.output + currentTime + "_cgal_DEBUG.off");
//	c3t3.output_boundary_to_off(off_file);
//	off_file.close();

	// Postprocessing
	std::cout << "  Post-processing mesh... " << std::flush; TicToc::tic();

	PolygonMesh processedScaffold;
	if (me_settings.CGAL_preserveEdges == TRUE) {
		// Extract surface triangulation
		PolygonMesh scaffold;
		CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, scaffold);
//		std::ofstream("scaffold_dump.off") << std::setprecision(17) << scaffold; // UNCOMMENT FOR DEBUGING...

		// Remove all but the largest component
		CGAL::Polygon_mesh_processing::keep_largest_connected_components(scaffold, 1);

		// Clip scaffold (using intersection)
		PolygonMesh outer_shell;
		CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(shell.points, shell.polygons, outer_shell);

		if( CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(scaffold, outer_shell, processedScaffold) == false)
			std::cout << "\n  WARNING: corefine_and_compute_intersection output may have only been corefined!" << std::endl;

	} else { // If preserveEdges != TRUE
		// Save volume mesh
		if (me_settings.volumeMesh == TRUE) {
			// Write output to .mesh file
			std::ofstream medit_file(me_settings.output + currentTime + "_cgal.mesh");
			c3t3.output_to_medit(medit_file);
			medit_file.close();

			//// Write output to .vtu file
			//std::ofstream vtu_file(me_settings.output + currentTime + "_cgal.vtu");
			//CGAL::output_to_vtu(vtu_file, c3t3);
			//vtu_file.close();
		}	

		// Extract surface triangulation
		CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, processedScaffold);
//		std::ofstream("scaffold_dump.off") << std::setprecision(17) << processedScaffold; // UNCOMMENT FOR DEBUGING...
	}

	// Remove all but the largest component from the surface triangulation
	CGAL::Polygon_mesh_processing::keep_largest_connected_components(processedScaffold, 1);

	// Reorient all faces coherently
	if (CGAL::is_closed(processedScaffold) == true) {
		CGAL::Polygon_mesh_processing::orient(processedScaffold, CGAL::parameters::outward_orientation(true));
		//CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(processedScaffold);
	} else {
		std::cout << "\n  WARNING: Polygon has border edges!" << std::endl;
	}

	std::cout << "Finished!" << std::endl; TicToc::toc();

	// Save surface to .stl file
	std::ofstream stl_fileOut(me_settings.output + currentTime + "_cgal.stl", std::ios::out | std::ios::binary);
	if ( me_settings.STLFormat == "ASCII" ) {
		CGAL::set_mode(stl_fileOut, CGAL::IO::ASCII);
		stl_fileOut.precision(std::numeric_limits<double>::digits10 + 2);
	} else {
		CGAL::set_mode(stl_fileOut, CGAL::IO::BINARY);
	}
	CGAL::IO::write_STL(stl_fileOut, processedScaffold);
	stl_fileOut.close();

	std::cout << "\n  SURFACE TRIANGULATION: " << std::endl;
	std::cout << "  Number of vertices: " << processedScaffold.size_of_vertices() << std::endl;
	std::cout << "  Number of facets: " << processedScaffold.size_of_facets() << "\n" << std::endl;

	// Compute/Recompute volume mesh (when surface was meshed and edges are to be preserved)
	if (me_settings.CGAL_preserveEdges == TRUE && me_settings.volumeMesh == TRUE) {
//		std::string testFile = me_settings.output + currentTime + "_cgal.stl";
//		std::filesystem::path testFileOut = me_settings.output + currentTime + "_cgal_mmg_test.mesh";
//		MeshMMG::surface2volume(testFile, lt_type, lt_size, lt_feature, me_settings, testFileOut);

		std::vector<Point_3> points;
		std::vector<std::vector<std::size_t>> polygons;
		CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(processedScaffold, points, polygons);
		polehedral2volume(points, polygons, lt_type, lt_size, lt_feature, me_settings,
		                  me_settings.output + currentTime + "_cgal.mesh");
	}

	// ...
	//#ifdef CGAL_CONCURRENT_MESH_3
	//	//init.~task_scheduler_init(); // Deprecated
	//	//control.~global_control();
	//#endif
	
	// Cleanup
	delete poissonSurfaceReconstruction;

	return EXIT_SUCCESS;
}

bool MeshCGAL::polehedral2volume (std::vector<Point_3> points,
                                  std::vector<std::vector<std::size_t>> polygons,
                                  latticeType lt_type, latticeSize lt_size, 
	                                latticeFeature lt_feature, meshSettings me_settings,
                                  std::filesystem::path outputFile) {
	/* Discretizes the volume of a closed triangulation.
	 * Inputs:
	 *  points         : Structure containing the outer shell geometry
	 *  polygons       : Lattice type data
	 *  outputFile     : Output filename
	 *  surf           : Lattice size data
	 * Outputs:
	 *  .mesh file containing the volume mesh
	 */

	// Mesh user settings
	FT me_facetAngle = me_settings.CGAL_facetAngle;
	FT me_facetSize = me_settings.CGAL_facetSize;
	FT me_facetDistance = me_settings.CGAL_facetDistance;

	FT me_cellRadiusEdgeRatio = me_settings.CGAL_cellRadiusEdgeRatio;
	FT me_cellSize = me_settings.CGAL_cellSize;     

	FT edgesProtectionAngle = me_settings.CGAL_edgesProtectionAngle;

	// Concurrency settings
	int n_threads = 1;
	#ifdef CGAL_CONCURRENT_MESH_3
		if (me_settings.n_threads < 1) {
			me_settings.n_threads = 1;
		} else if (me_settings.n_threads > tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism)) {
			me_settings.n_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism); // tbb::this_task_arena::max_concurrency()
			std::cout << "  Warning: Requested number of threads exceeds the maximum allowed " 
			          << "parallelism. Number of threads set to the maximum allowed parallelism." << std::endl;
		}
		//tbb::task_scheduler_init init(me_settings.n_threads); // Deprecated
		//n_threads = tbb::this_task_arena::max_concurrency();
		
		tbb::global_control control(tbb::global_control::max_allowed_parallelism, me_settings.n_threads);
		n_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
	#endif

	std::cout << "  Meshing volume (using " << n_threads << " threads)... " << std::flush; TicToc::tic();

	// Convert to polygon mesh
	F_Polyhedron polygonSurface;
	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, polygonSurface);

	// Create domain
	F_Polyhedron_domain domain(polygonSurface);

  // Get sharp features
  domain.detect_features(edgesProtectionAngle);
	
	// Sizing fields
	TPMS_dependent_unitcellsize_field_P facetSize;
	facetSize.lt_size = &lt_size;
	facetSize.parameter = &me_facetSize;

	TPMS_dependent_unitcellsize_field_P facetDistance;
	facetDistance.lt_size = &lt_size;
	facetDistance.parameter = &me_facetDistance;

	TPMS_dependent_wallsize_field_P cellSize;
	cellSize.lt_type = &lt_type;
	cellSize.lt_size = &lt_size;
	cellSize.lt_feature = &lt_feature;
	cellSize.parameter = &me_cellSize;

	// Mesh criteria
	F_Mesh_criteria criteria(CGAL::parameters::edge_size = facetSize,
	                         CGAL::parameters::facet_angle = me_facetAngle,     // Min triangle angle (degrees)
	                         CGAL::parameters::facet_size = facetSize,         // Max triangle size
	                         CGAL::parameters::facet_distance = facetDistance, // Surface approximation error
	                         CGAL::parameters::cell_radius_edge_ratio = me_cellRadiusEdgeRatio,	// Mesh tetrahedra radius|edge ratio upper bound
	                         CGAL::parameters::cell_size = cellSize);

	// Mesh generation
	F_C3t3 c3t3 = CGAL::make_mesh_3<F_C3t3>(domain, criteria);//, 
	                                        //CGAL::parameters::manifold()); // Causes assertion violation error!

	// Write output to .mesh file
	outputFile.replace_extension(".mesh");
	std::ofstream medit_file(outputFile.string());
	c3t3.output_to_medit(medit_file);
	medit_file.close();

	// ...
	//#ifdef CGAL_CONCURRENT_MESH_3
	//	//init.~task_scheduler_init(); // Deprecated
	//	//control.~global_control();
	//#endif

	std::cout << "Finished!" << std::endl; TicToc::toc();

	std::cout << "\n  VOLUME MESH: " << std::endl;
	//std::cout << "  Number of vertices: " << c3t3.number_of_vertices_in_complex() << std::endl; // Returns incorrect value
	//std::cout << "  Number of edges: " << c3t3.number_of_edges() << std::endl; // Returns incorrect value
	//std::cout << "  Number of facets: " << 2 * c3t3.number_of_facets() << std::endl; // ???
	std::cout << "  Number of tetrahedra: " << c3t3.number_of_cells() << std::endl;

	return EXIT_SUCCESS;
}

void MeshCGAL::implicitPoissonReconstruction(outerShell &shell, double offset, 
	std::vector<double> &boundingBox, Poisson_reconstruction_function **poissonSurfaceReconstruction) {
	/* Creates an implicit surface reconstruction using the Poisson method.
	 * Inputs:
	 *   points : Point cloud to be reconstructed
	 *   offset : Offset to give to reconstructed surface
   * Output:
	 *   poissonSurfaceReconstruction : Poisson surface reconstruction
	 */

	const int nNeighbors = 18; // Number of neigbours for normal estimation

	std::cout << "  Creating implicit surface description of the outer shell... " << std::flush; TicToc::tic();

	PointList scaledPoints;
	boundingBox = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
	bool useFaster = true; // TEMP!! Note: The faster offset is much faster but not completelly correct (although it usually works just fine). The slower one is significantlly slower but correct.
	if (useFaster == 1) {
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
	//	if (CGAL::Polygon_mesh_processing::does_bound_a_volume(pointsCGAL) == false)
	//		throw "ERROR: Input outer surface does not bound a volume";

		// Estimate and orient normals
		CGAL::jet_estimate_normals<Concurrency_tag>(pointsCGAL,
																								nNeighbors,
																								CGAL::parameters::
																									point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
																									normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

		std::vector<Point_with_normal>::iterator unoriented_points_begin = 
			CGAL::mst_orient_normals(pointsCGAL, nNeighbors,
															CGAL::parameters::
																point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
																normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

		// Remove points with an unoriented normal
		pointsCGAL.erase(unoriented_points_begin, pointsCGAL.end());

		// Scaling
		for (size_t i=0; i<pointsCGAL.size(); ++i) {
//			pointsCGAL[i].first = Point_3(pointsCGAL[i].first.x() + offset*pointsCGAL[i].second.x(),
//																		pointsCGAL[i].first.y() + offset*pointsCGAL[i].second.y(),
//																		pointsCGAL[i].first.z() + offset*pointsCGAL[i].second.z());
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

	} else {
		TriangleMesh triangle_mesh;
		CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(shell.points, shell.polygons, triangle_mesh);
		const double averageSpacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(shell.points, nNeighbors);

//	// Check if points bound a volume
//	if (CGAL::is_closed(triangle_mesh) == false)
//		throw "ERROR: Input outer surface is not closed";
//	if (CGAL::Polygon_mesh_processing::does_bound_a_volume(triangle_mesh) == false)
//		throw "ERROR: Input outer surface does not bound a volume";

		// Compute offset
		outerShell offset_shell;
		if (offset > 0) {
			double facetAngle = 30;
			double facetSize = averageSpacing/2;
			double facetDistance = 0.1;
			offset_shell = cgal_off_meshing(triangle_mesh, offset, facetAngle, facetSize, facetDistance);

		} else {
			offset_shell = shell;
		}

		// Convert to Point_with_normal
		for (size_t i=0; i<offset_shell.points.size(); ++i) {
			Point_with_normal p(Point_3(offset_shell.points[i].x(), offset_shell.points[i].y(), offset_shell.points[i].z()), 
													Vector_3(0.0, 0.0, 0.0));
			scaledPoints.push_back(p);

			// Compute bounding box of scaled point cloud
			boundingBox[0] = std::min(boundingBox[0], scaledPoints[i].first.x());
			boundingBox[1] = std::max(boundingBox[1], scaledPoints[i].first.x());
			boundingBox[2] = std::min(boundingBox[2], scaledPoints[i].first.y());
			boundingBox[3] = std::max(boundingBox[3], scaledPoints[i].first.y());
			boundingBox[4] = std::min(boundingBox[4], scaledPoints[i].first.z());
			boundingBox[5] = std::max(boundingBox[5], scaledPoints[i].first.z());
		}

		// Estimate and orient normals
		CGAL::jet_estimate_normals<Concurrency_tag>(scaledPoints, nNeighbors,
																								CGAL::parameters::
																									point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
																									normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

		std::vector<Point_with_normal>::iterator unoriented_points_begin = 
			CGAL::mst_orient_normals(scaledPoints, nNeighbors,
															CGAL::parameters::
																point_map(CGAL::First_of_pair_property_map<Point_with_normal>()).
																normal_map(CGAL::Second_of_pair_property_map<Point_with_normal>()));

		// Remove points with an unoriented normal
		scaledPoints.erase(unoriented_points_begin, scaledPoints.end());

		// Scott Meyer's "swap trick" to trim excess capacity (Seems to solve the occasional 'CGAL ERROR: assertion violation! .../Refine_facets_3.h Line: 862')
		std::vector<Point_with_normal>(scaledPoints).swap(scaledPoints);

//	{
//			// Save surface to .stl file (FOR DEBUG PURPOSES!)
//			TriangleMesh triangle_mesh_temp;
//			CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(offset_shell.points, offset_shell.polygons, triangle_mesh_temp);
//
//			std::ofstream stl_fileOut("out_test.stl", std::ios::out | std::ios::binary);
//			if ( me_settings.STLFormat == "ASCII" ) {
//				CGAL::set_mode(stl_fileOut, CGAL::IO::ASCII);
//				stl_fileOut.precision(std::numeric_limits<double>::digits10 + 2);
//			} else {
//				CGAL::set_mode(stl_fileOut, CGAL::IO::BINARY);
//			}
//			CGAL::IO::write_STL(stl_fileOut, triangle_mesh_temp);
//			stl_fileOut.close();
//	}
	}

	// Create implicit function
	*poissonSurfaceReconstruction = new Poisson_reconstruction_function(scaledPoints.begin(), scaledPoints.end(),
	                                                                   CGAL::First_of_pair_property_map<Point_with_normal>(),
	                                                                   CGAL::Second_of_pair_property_map<Point_with_normal>());
	// Compute the Poisson indicator function f() at each vertex of the triangulation.
	if ( ! (*poissonSurfaceReconstruction)->compute_implicit_function() )
		throw "ERROR_POISSON_LINEAR_SOLVER_FAILED";
	
  std::cout << "Finished!" << std::endl; TicToc::toc();

	//return EXIT_SUCCESS;
}



double MeshCGAL::internal::signedDistance (Point_3 p, latticeType *lt_type, 
                      latticeSize *lt_size, latticeFeature *lt_feature, std::string side,
                      Poisson_reconstruction_function *poissonReconstruction) {
	Point point(p.x(), p.y(), p.z());

	int mode = 0;
	switch (mode) {
		case 1: // Signed distance function of outer shape (For debug purposes)
			return poissonReconstruction->operator()(p);

		case 2: { // Signed distance function of a sphere (For debug purposes)
			const double x2 = p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
			return x2+y2+z2-1;
		}
		case 3: { // Signed distance function of a lattice filled sphere (For debug purposes)
			const double x2 = p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
			double a = (x2+y2+z2-1); // Outer shape
			double b = Infill::TPMS_function (point, lt_type, lt_size, lt_feature); // Infill
				
			if (a > b) {return a;} else {return b;}
		}
		default: { // Signed distance function of lattice filled outer shape
			double a = poissonReconstruction->operator()(p); // Outer shape
			double b = Infill::TPMS_function (point, lt_type, lt_size, lt_feature); // Infill

			if (side == "void") // Select a sub-domain to keep
				if (a > -b) {return a;} else {return -b;}
			else {
				if (a > b) {return a;} else {return b;}
			}
		}
	}
}

double MeshCGAL::internal::Infill (Point_3 p, latticeType *lt_type, latticeSize *lt_size,
                                 latticeFeature *lt_feature) {
	Point point(p.x(), p.y(), p.z());
	return Infill::TPMS_function (point, lt_type, lt_size, lt_feature);
}