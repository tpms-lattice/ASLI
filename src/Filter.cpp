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

#include "Filter.h"
/* Filter
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

double Filter::filter(neighbourData neighbours, 
                      std::vector<double> unfilteredValues,
                      std::string filterType, double filterRadius){
	/* Filters the data using the chosen weight function
	*  Inputs:
	*    neighbours        : Structure containing IDX and distance information
	*                        of the values in the neighbourhood of the unfiltered
	*                        value to be filtered
	*    unfilteredValues  : Unfiltered values
  *    filterType        : Filter type
	*    filterRadius      : Filter radius
	*  Return:
	*    Filtered value
	*/

	double a = 0, b = 0;

	for (size_t i = 0; i < neighbours.distance.size(); i++) {
		double weight = internal::weightFunction(filterType, 
		                                         neighbours.distance[i],
		                                         filterRadius);
															
		a += weight * unfilteredValues[neighbours.IDX[i]];
		b += weight;
	}

	return a/b;
}

double Filter::internal::weightFunction(std::string type, double distance,
                                        double radius) {

	if (type == "linear") { // Linear weight function
		if (std::abs(distance) <= radius)
			return 1 - (1/radius) * std::abs(distance);
		else
			return 0;

	} else { // Default: Gaussian weight function
		double sigma = radius/2;
		return std::exp( -std::pow(distance/sigma, 2) );
	}
}