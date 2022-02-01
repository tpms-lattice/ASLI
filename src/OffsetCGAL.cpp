/* ==========================================================================
 *  Based on CGAL's offset meshing plugin. See CGAL demos.
 * ==========================================================================*/

#include "OffsetCGAL.h"

CGAL::Offset_function<TriangleMesh, EPICK>
offset_function(TriangleMesh* surface_mesh_ptr, double offset_value) {
	return { *surface_mesh_ptr, offset_value };
}

outerShell cgal_off_meshing(TriangleMesh &triangle_mesh, double offset,
                            double facetAngle, double facetSize,
                            double facetDistance) {
	typedef EPICK GT;
	typedef CGAL::Labeled_mesh_domain_3<GT, int, int> Mesh_domain;
	typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
	typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
	typedef GT::Sphere_3 Sphere_3;

	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(triangle_mesh);

	GT::Point_3 center((bbox.xmax()+bbox.xmin())/2,
	                   (bbox.ymax()+bbox.ymin())/2,
	                   (bbox.zmax()+bbox.zmin())/2);
	double sqrad = 0.6 * std::sqrt( CGAL::square(bbox.xmax()-bbox.xmin())+
	                                CGAL::square(bbox.ymax()-bbox.ymin())+
	                                CGAL::square(bbox.zmax()-bbox.zmin()) )
	              + offset;
	sqrad = CGAL::square(sqrad);

	Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain (
		offset_function(&triangle_mesh, offset),
		Sphere_3(center, sqrad),
		CGAL::parameters::relative_error_bound = 1e-7,
		CGAL::parameters::construct_surface_patch_index = [](int i, int j) { return (i * 1000 + j); }
	);

	CGAL::Mesh_facet_topology topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH;
	topology = CGAL::Mesh_facet_topology(topology | CGAL::MANIFOLD);
	Mesh_criteria criteria(CGAL::parameters::facet_angle = facetAngle,
	                       CGAL::parameters::facet_size = facetSize,
	                       CGAL::parameters::facet_distance = facetDistance,
	                       CGAL::parameters::facet_topology = topology);

	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
	                                    CGAL::parameters::no_perturb(),
	                                    CGAL::parameters::no_exude());

	// Convert to polygon soup
	TriangleMesh surfaceMesh;
	outerShell offsetShell;
	CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, surfaceMesh);
	CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(surfaceMesh, offsetShell.points, offsetShell.polygons);

	return offsetShell;
}