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

#include "TrilinearInterpolation.h"
/* RBF interpolation for 3D data using a linear kernel
 *
 * Author(s): F.P.B. (KU Leuven)
 *            M.B. (KU Leuven)
 */

bool TrilinearInterpolation::setup(const std::vector<std::vector<double>> &data,
                                   modelData &interpolationModel) {
	int nDataPoints = data.size();
	Eigen::Matrix3Xd coordinates(3, nDataPoints);
	Eigen::VectorXd y(nDataPoints);

	for (int i = 0; i < nDataPoints; i++) {
		coordinates.col(i) = Eigen::Vector3d(data[i][0],
		                                     data[i][1],
		                                     data[i][2]);
		y(i) = data[i][3];				
	}

	return setup(coordinates, y, interpolationModel);
}

bool TrilinearInterpolation::setup(const Eigen::Matrix3Xd &coordinates,
                                   const Eigen::VectorXd &y, 
                                   modelData &interpolationModel) {
	/* Computes the interpolation coefficients.
	 *  Inputs:
	 *   coordinates : Coordinates of the scattered data to be interpolated
	 *   y           : Outputs corresponding to the provided coordinates
	 *  Output:
	 *   interpolationModel : Structure containing the interpolation model data
	 */

	// Add data points to model data structutre
	interpolationModel.x = coordinates;

	// Compute the interpolation matrix
	Eigen::MatrixXd distanceMatrix = internal::distanceMatrix(coordinates, coordinates);
	Eigen::MatrixXd interpolationMatrix = distanceMatrix.unaryExpr(&internal::kernel);

//	Eigen::MatrixXd interpolationMatrix(coordinates.cols(), coordinates.cols());
//	for (size_t i = 0; i < coordinates.cols(); i++) {
//		for (size_t j = 0; j < coordinates.cols(); j++) {
//			interpolationMatrix(i,j) = TrilinearInterpolation::internal::kernel(distanceMatrix(i,j));
//		}
//	}

	// Augmented interpolation matrix (includes polynomial term)
	Eigen::VectorXd onesVec = Eigen::VectorXd::Ones(coordinates.cols());
	Eigen::MatrixXd P(y.size(), coordinates.rows() + 1);
	P << onesVec, coordinates.transpose();
	Eigen::MatrixXd zerosMat = Eigen::MatrixXd::Zero(coordinates.rows()+1, coordinates.rows()+1);
	Eigen::MatrixXd interpolationMatrix_augmented(interpolationMatrix.rows()+coordinates.rows()+1,
	                                              interpolationMatrix.cols()+coordinates.rows()+1);

	interpolationMatrix_augmented << interpolationMatrix, P, P.transpose(), zerosMat;
	
	// Augmented y (includes polynomial term)
	Eigen::VectorXd y_augmented(y.size() + coordinates.rows() + 1);
	Eigen::VectorXd zerosVec =  Eigen::VectorXd::Zero(coordinates.rows()+1);
	y_augmented << y, zerosVec;

	// Compute the interpolation coefficents. Using LU decomposition with partial
	// pivoting since distance matrix will always be invertible.
	interpolationModel.interpolationCoefficients = 
		interpolationMatrix_augmented.partialPivLu().solve(y_augmented);

	return EXIT_SUCCESS;
}

double TrilinearInterpolation::evaluate(const Eigen::Vector3d p,
                                        const modelData &interpolationModel) {
	/* Evaluates the RBF-interpolation model
	 * Inputs:
	 *   p                  : Point for which output is to be estimated
	 *   interpolationModel : Structure containing the interpolation model data
	 * Return:
	 *   Output estimation at point p
	 */

	// Compute the interpolation vector
	Eigen::VectorXd interpolationVector = internal::distanceVector(interpolationModel.x, p);
	for (size_t i = 0; i < interpolationVector.size(); i++) {
		interpolationVector(i) = TrilinearInterpolation::internal::kernel(interpolationVector(i));
	}

	// Augmented interpolation vector (includes polynomial term)
	Eigen::VectorXd interpolationVector_augmented(interpolationVector.size()+p.size()+1);
	interpolationVector_augmented << interpolationVector, 1, p;

	// Estimate approximation
	return interpolationVector_augmented.dot(interpolationModel.interpolationCoefficients);
}

Eigen::MatrixXd TrilinearInterpolation::internal::distanceMatrix (const Eigen::Matrix3Xd &A,
                                                                  const Eigen::Matrix3Xd &B) {
	/* Computes a matrix containing the pairwise distances between the elements
	 * of A and B.
	 */

	int nDataPointsA = A.cols();
	int nDataPointsB = B.cols();

	Eigen::MatrixXd distances(nDataPointsB, nDataPointsA);
	for (size_t i = 0; i < nDataPointsA; i++) {
		distances.col(i) = (B.colwise() - A.col(i)).matrix().colwise().norm();

		//for (size_t j = 0; j < nDataPointsB; j++) {
		//	distances(j, i) = std::sqrt(std::pow(B->col(j)[0] - A->col(i)[0], 2.0)+
		//	                            std::pow(B->col(j)[1] - A->col(i)[1], 2.0)+
		//	                            std::pow(B->col(j)[2] - A->col(i)[2], 2.0));
		//}
	}

	return distances;
}

Eigen::VectorXd TrilinearInterpolation::internal::distanceVector(const Eigen::Matrix3Xd &A,
                                                                 const Eigen::Vector3d &B) {
	/* Computes a vector containing the pairwise distances between the elements
	 * of A and point B.
	 */

	int nDataPointsA = A.cols();

	Eigen::VectorXd distances(nDataPointsA);
	for (size_t i = 0; i < nDataPointsA; i++) {
		distances(i) = (B.transpose().transpose() - A.col(i)).norm();

		//distances(i) = (B->colwise() - A->col(i)).norm();

		//distances(i) = std::sqrt(std::pow((*B)(0) - A->col(i)[0], 2.0)+
		//                         std::pow((*B)(1) - A->col(i)[1], 2.0)+
		//                         std::pow((*B)(2) - A->col(i)[2], 2.0));
	}

	return distances;
}

double TrilinearInterpolation::internal::kernel(const double &r) {
	/* Basis function */
	return r; // Linear kernel
}