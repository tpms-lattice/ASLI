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

#ifndef TRILINEAR_INTERPOLATION_H
#define TRILINEAR_INTERPOLATION_H

/* Eigen headers */
#include <Eigen/Core>
#include <Eigen/Dense>

/* Standard library headers */
#include <iostream>

struct modelData {
	Eigen::VectorXd interpolationCoefficients;
	Eigen::Matrix3Xd x;
};

namespace TrilinearInterpolation {
	bool setup(const std::vector<std::vector<double>> &data,
	           modelData &interpolationModel);
	bool setup(const Eigen::Matrix3Xd &coordinates,
	           const Eigen::VectorXd &y,
	           modelData &interpolationModel);

	double evaluate(const Eigen::Vector3d p,
	                const modelData &interpolationModel);

	namespace internal {
		double kernel(const double &r);

		Eigen::MatrixXd distanceMatrix(const Eigen::Matrix3Xd &A,
		                               const Eigen::Matrix3Xd &B);
		Eigen::VectorXd distanceVector(const Eigen::Matrix3Xd &A,
		                               const Eigen::Vector3d &B);
	}
}
#endif