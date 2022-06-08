/* ==========================================================================
 *  Based on CGAL's offset meshing plugin. See CGAL demos.
 * ==========================================================================*/

#ifndef OFFSETCGAL_H
#define OFFSETCGAL_H

#include "Mesh.h"
#include "SMesh_type.h"

/* CGAL headers */
// IO
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

// Polygon processing
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Side_of_triangle_mesh.h>

// 3D Fast Intersection and Distance Computation
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

// 3D volume mesh
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_triangulation_3.h>

/* Boost headers */
#include <boost/array.hpp>

/* Standard library headers */
#include <iostream>

namespace CGAL {
	template <class TriangleMesh, class GeomTraits>
	class Offset_function {
		typedef AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
		typedef AABB_traits<GeomTraits, Primitive> Traits;
		typedef AABB_tree<Traits> Tree;
		typedef Side_of_triangle_mesh<TriangleMesh, GeomTraits> Side_of;

	public:
		Offset_function(TriangleMesh& tm, double offset_distance)
			: m_tree_ptr(new Tree(boost::begin(faces(tm)),
			                      boost::end(faces(tm)),
			                      tm) )
			, m_side_of_ptr( new Side_of(*m_tree_ptr) )
			, m_offset_distance(offset_distance)
			, m_is_closed( is_closed(tm) )
		{
			CGAL_assertion(!m_tree_ptr->empty());
		}

		double operator()(const typename GeomTraits::Point_3& p) const {
			using CGAL::sqrt;

			Bounded_side side = m_is_closed?m_side_of_ptr->operator()(p):ON_UNBOUNDED_SIDE;
			if (side==ON_BOUNDARY) return m_offset_distance;

			typename GeomTraits::Point_3 closest_point = m_tree_ptr->closest_point(p);
			double distance = std::sqrt(squared_distance(p, closest_point));

			return (side == ON_UNBOUNDED_SIDE ? -distance : distance) + m_offset_distance;
		}

	private:
		boost::shared_ptr<Tree> m_tree_ptr;
		boost::shared_ptr<Side_of> m_side_of_ptr;
		double m_offset_distance;
		bool m_is_closed;
	};
} //end of CGAL namespace

outerShell cgal_off_meshing(TriangleMesh &triangle_mesh, double offset,
                            double facetAngle, double facetSize,
                            double facetDistance);
#endif