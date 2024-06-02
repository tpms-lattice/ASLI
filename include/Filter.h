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

#ifndef FILTER_H
#define FILTER_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

struct neighbourData {
	std::vector<int> IDX;
	std::vector<double> distance;
};

namespace Filter {
	double filter(const neighbourData &neighbours, 
	              const std::vector<double> &unfilteredValues,
	              const std::string &weightFunction,
	              const double &filterRadius);

	namespace internal {
		double weightFunction(const std::string &type,
		                      const double &distance,
		                      const double &radius);
	}
};

#endif