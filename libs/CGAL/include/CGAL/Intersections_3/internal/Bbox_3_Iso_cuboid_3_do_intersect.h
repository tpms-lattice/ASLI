// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// Copyright (c) 2010, 2014 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.2.3/Intersections_3/include/CGAL/Intersections_3/internal/Bbox_3_Iso_cuboid_3_do_intersect.h $
// $Id: Bbox_3_Iso_cuboid_3_do_intersect.h 52164b1 2019-10-19T15:34:59+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_DO_INTERSECT_H

#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
bool do_intersect(const CGAL::Bbox_3& bb,
                  const typename K::Iso_cuboid_3& ic,
                  const K& /* k */)
{
  if (bb.xmax() < ic.xmin() || ic.xmax() < bb.xmin())
    return false;
  if (bb.ymax() < ic.ymin() || ic.ymax() < bb.ymin())
    return false;
  if (bb.zmax() < ic.zmin() || ic.zmax() < bb.zmin())
    return false;
  return true;
}

} // namespace internal
} // namespace Intersections
} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_DO_INTERSECT_H
