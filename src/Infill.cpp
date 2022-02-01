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

#include "Infill.h"

/* Infill contains the functions required to create lattice infills based on 
 * triply periodic minimal surfaces. 
 * 
 * Note: Strictly speaking Infill can handle infills described by any signed 
 * distance function, not only triply periodic minimal surface based infills.
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

double Infill::internal::input2level(std::string type, double scaling, std::string feature, 
                                     double featureValue, latticeMaterial materialParameters,
                                     std::string featureMode) {
	/* Estimates the level-set constant that corresponds to the provided input
	 * parameter. Estimation avaliable for: gyroid, diamond, primitive and IWP
	 * in both strut and sheet variants.
	 * Inputs:
	 *   type               : Unit cell type
	 *   scaling            : Unit cell scaling
	 *   feature            : Feature to be converted
	 *   featureValue       : Feature value
	 *   materialParameters : Material model parameters (used if the feature 
	 *                        provided is the elastic mudulus)
	 *   featureMode        : 
	 * Return:
	 *   Level-set constant
	 */

	if (feature == "isovalue") {
		return unnormalizeLevel(featureValue, type);
	} else if (feature == "volumeFraction") {
		return vFraction2level(featureValue, type);
	} else if (feature == "wallSize") {
		return wallSize2level(featureValue, scaling, type, featureMode);
	} else if (feature == "poreSize") {
		return poreSize2level(featureValue, scaling, type, featureMode);
	} else if (feature == "elasticModulus") {
		return vFraction2level(eModulus2vFraction(featureValue, materialParameters), type);
	} else {
		std::cerr << "ERROR_INVALID_FEATURE_REQUEST" << std::endl;
		exit(EXIT_FAILURE);
	}
}

double Infill::internal::unnormalizeLevel(double t_normalized, std::string type) {
	/* Reverts the min-max normalization of the provided normalized (range 
	 * [0, 1]) level-set constant. Avaliable for: gyroid, diamond, primitive
	 * and IWP in both strut and sheet variants.
	 * Inputs:
	 *   t_normalized : Nomralized level-set constant (in the range [0, 1])
	 *   type   : Unit cell type
	 * Return:
	 *   Level-set constant
	 */

	if (type == "gyroid") {
		return (t_normalized * (1.5 + 1.5)) + -1.5;

	} else if (type == "sheet_gyroid") {
		return (t_normalized * (1.55 - 0)) + 0;

	} else if (type == "diamond") {
		return (t_normalized * (1.425 - -1.425)) + -1.425;
		
	} else if (type == "sheet_diamond") {
		return (t_normalized * (1.45 - 0)) + 0;

	} else if (type == "primitive") {
		return (t_normalized * (3 - -3)) + -3;
		
	} else if (type == "sheet_primitive") {
		return (t_normalized * (3.05 - 0)) + 0;

	} else if (type == "IWP") {
		return (t_normalized * (3.025 - -5.025)) + -5.025;

	} else if (type == "sheet_IWP") {
		return (t_normalized * (5.05 - 0)) + 0;

	} else {
		std::cerr << "ERROR_INVALID_TPMS" << std::endl;
		exit(EXIT_FAILURE);
	}
}

double Infill::internal::vFraction2level(double volFraction, std::string type) {
	/* Estimates the level-set constant that corresponds to the requested volume
	 * fraction. Estimation avaliable for: gyroid, diamond, primitive and IWP in 
	 * both strut and sheet variants.
	 * Inputs:
	 *   volFraction : Volume fraction to be converted
	 *   type        : Unit cell type
	 * Return:
	 *   Level-set constant
	 */

	// Polinomial coeficients
	double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;
	
	// Enforce volume fraction bounds [0, 1]
	if (volFraction < 0) { volFraction = 0; }
	else if (volFraction > 1) { volFraction = 1; }

	// Set the polynomial coefficients (for the requested unit cell type)
	if (type == "gyroid") {
		a = -1.5; b = 3.0;

	} else if (type == "sheet_gyroid") {
		if (volFraction > 0.999) { volFraction = 0.999; } // TEMP!! Check why it fails with value > 0.999...
		a = -0.00224e-3; b = 1.59	; c = -0.117;

	} else if (type == "diamond") {
		a = -1.36; b = 5.61; c = -19.59; d = 49.46; e = -54.60; f = 21.84;

	} else if (type == "sheet_diamond") {
		a = -0.0108; b = 1.82; c = -5.48; d = 17.74; e = -23.55; f = 10.88;

	} else if (type == "primitive") {
		a = -2.86; b = 18.57; c = -74.78; d = 170.67; e = -181.17; f = 72.43;

	} else if (type == "sheet_primitive") {
		a = -0.0250; b = 2.83; c = -8.89; d = 26.90; e = -34.93; f = 17.03;

	} else if (type == "IWP") {
		a = 3.08; b = -6.61; c = 54.92; d = -393.14; e = 1219.46; f = -1936.45;
		g = 1529.72; h = -475.82;

	} else if (type == "sheet_IWP") {
		a = -0.0343; b = 6.70; c = -48.91; d = 325.14; e = -1049.81; f = 1749.51;
		g = -1448.85; h = 471.14;

	} else {
		std::cerr << "ERROR_INVALID_TPMS" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Compute the level-set constant
	return a + b*volFraction + c*std::pow(volFraction,2) + 
	                           d*std::pow(volFraction,3) + 
	                           e*std::pow(volFraction,4) + 
	                           f*std::pow(volFraction,5) + 
	                           g*std::pow(volFraction,6) + 
	                           h*std::pow(volFraction,7);
}

double Infill::internal::wallSize2level(double wallSize, double scaling, std::string type, std::string mode) {
	/* Estimates the level-set constant that corresponds to the requested wall
	 * thickness. Estimation avaliable for: gyroid, diamond, primitive and IWP in
	 * both strut and sheet variants.
	 * Inputs:
	 *   wallSize : Wall thickness to be converted
	 *   scaling  : Unit cell scaling
	 *   type     : Unit cell type
	 *   mode     : 
	 *  Return:
	 *   Level-set constant
	 */

	// Polinomial coeficients
	double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;

	// Enforce lower wall thickness bound [0, __]
	if (wallSize < 0) {wallSize = 0;}

	// Compute 'unit' wall size
	if (mode == "absolute")
		wallSize = wallSize / (2*PI / scaling);

	// Set the polynomial coefficients and enforce the upper wall thickness bound
	// (for the requested unit cell type)
	if (type == "gyroid") {
		a = -1.43; b = 1.05; c = 2.80; d = 19.83; e = -37.50; f = 16.70;
		if (wallSize > 0.848) {wallSize = 0.848;}

	} else if (type == "sheet_gyroid") {
		a = -0.0350; b = 6.55; c = -6.99;
		if (wallSize > 0.425) {wallSize = 0.425;}

	} else if (type == "diamond") {
		a = -0.984; b = 0.904; c = 10.63; d = -10.19;
		if (wallSize > 0.71) {wallSize = 0.71;}

	} else if (type == "sheet_diamond") {
		a = -0.0396; b = 7.07; c = -13.42; d = 10.71;
		if (wallSize > 0.48) {wallSize = 0.48;}

	} else if (type == "primitive") {
		a = -0.999; b = 0.604; c = 3.55; d = 2.18; e = -4.21; f = 1.19;
		if (wallSize > 1.285) {wallSize = 1.285;}

	} else if (type == "sheet_primitive") {
		a = 7.80e-05; b = 5.44; c = 0.0212; d = -3.07; e = 0.170	; f = 0.355;
		if (wallSize > 0.86) {wallSize = 0.86;}

	} else if (type == "IWP") {
		a = 3.00; b = -0.156; c = -29.07; d = -3.51; e = 50.44; f = -24.81;
		if (wallSize > 0.79) {wallSize = 0.79;}

	} else if (type == "sheet_IWP") {
		a = -0.0134; b = 17.61; c = -77.60; d = 943.10; e = -5727.87; f = 15460.84;
		g = -18974.65; h = 8707.20;
		if (wallSize > 0.67) {wallSize = 0.67;}

	} else {
		std::cerr << "ERROR_INVALID_TPMS" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Compute the level-set constant
	return (a + b*wallSize + c*std::pow(wallSize,2) + 
	                         d*std::pow(wallSize,3) + 
	                         e*std::pow(wallSize,4) + 
	                         f*std::pow(wallSize,5) + 
	                         g*std::pow(wallSize,6) + 
	                         h*std::pow(wallSize,7));
}

double Infill::internal::poreSize2level(double poreSize, double scaling, std::string type, std::string mode) {
	/* Estimates the level-set constant that corresponds to the requested pore
	 * size. Estimation avaliable for: gyroid, diamond, primitive and IWP in 
	 * both strut and sheet variants.
	 * Inputs:
	 *   poreSize : Pore size to be converted
	 *   scaling  : Unit cell scaling
	 *   type     : Unit cell type
	 *   mode     : 
	 *  Return:
	 *   Level-set constant
	 */

	// Polinomial coeficients
	double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;

	// Enforce lower pore size bound [0, __]
	if (poreSize < 0) {poreSize = 0;} 

	// Compute 'unit' pore size
	if (mode == "absolute") 
		poreSize = poreSize / (2*PI / scaling);

	// Set the polynomial coefficients and enforce the upper pore size bound
	// (for the requested unit cell type)
	if (type == "gyroid") {
		a = 1.52; b = -0.486; c = -7.57; d = -6.28; e = 21.03; f = -9.61;
		if (poreSize > 0.848) {poreSize = 0.848;}

	} else if (type == "sheet_gyroid") {
		a = 1.47; b = -0.314; c = -6.84; d = -15.23; e = 50.37; f = -40.61;
		if (poreSize > 0.425) {poreSize = 0.425;}

	} else if (type == "diamond") {
		a = 1.40; b = -0.497; c = -11.42; d = 9.63;
		if (poreSize > 0.78) {poreSize = 0.78;}

	} else if (type == "sheet_diamond") {
		a = 1.41; b = -1.09; c = 3.81; d = -65.70; e = 153.92; f = -119.84;
		if (poreSize > 0.41) {poreSize = 0.41;}

	} else if (type == "primitive") {
		a = 3.00; b = -0.504; c = -3.83; d = -1.56; e = 2.68; f = -0.619;
		if (poreSize > 1.61) {poreSize = 1.61;}

	} else if (type == "sheet_primitive") {
		a = 3.05; b = -1.42; c = 1.23; d = -13.67; e = 15.75; f = -5.86;
		if (poreSize > 0.84) {poreSize = 0.84;}

	} else if (type == "IWP") {
		a = -4.94; b = -0.408; c = 44.62; d = -19.91; e = -64.79; f = 50.55;
		if (poreSize > 0.77) {poreSize = 0.77;}

	} else if (type == "sheet_IWP") {
		a = 4.96; b = -0.511; c = -33.51; d = -32.62; e = 173.24; f = -133.09;
		if (poreSize > 0.42) {poreSize = 0.42;}

	} else {
		std::cerr << "ERROR_INVALID_TPMS" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Compute the level-set constant
	return (a + b*poreSize + c*std::pow(poreSize,2) + 
	                         d*std::pow(poreSize,3) + 
	                         e*std::pow(poreSize,4) + 
	                         f*std::pow(poreSize,5) + 
	                         g*std::pow(poreSize,6) + 
	                         h*std::pow(poreSize,7));
}


double Infill::internal::level2wallSize(double t, double scaling, std::string type) {
	/* Estimates the wall thickness that corresponds to the requested level-set
	 * constant. Estimation avaliable for: gyroid, diamond, primitive and IWP in
	 * both strut and sheet variants.
	 * Inputs:
	 *   t       : Level-set constant
	 *   scaling : Unit cell scaling
	 *   type    : Unit cell type
	 *  Return:
	 *   Wall thickness
	 */

	// Polinomial coeficients
	double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;

	// Enforce the level-set constant bounds and set the polynomial coefficients
	// (for the requested unit cell type)
	if (type == "gyroid") {
		if (t < -1.5) {t = -1.5;}
		else if (t > 1.5) {t = 1.5;}

		a = 0.415; b = 0.197; c = -0.00343; d = -0.0202; e = -0.00393; f = 0.0273;

	} else if (type == "sheet_gyroid") {
		if (t < 0) {t = 0;}
		else if (t > 1.55) {t = 1.55;}

		a = -0.00584; b = 0.387; c = -1.24; d = 2.73; e = -2.43; f = 0.767;

	} else if (type == "diamond") {
		if (t < -1.425) {t = -1.425;}
		else if (t > 1.425) {t = 1.425;}

		a = 0.305; b = 0.222; c = 0.00838; d = 0.00314; e =-0.0462; f = 0.0470;

	} else if (type == "sheet_diamond") {
		if (t < 0) {t = 0;}
		else if (t > 1.45) {t = 1.45;}

		a = -0.00219; b = 0.226; c = -0.155; d = 0.160;

	} else if (type == "primitive") {
		if (t < -3) {t = -3;}
		else if (t > 3) {t = 3;}

		a = 0.442; b = 0.251; c = -0.0551; d = 0.0690; e = -0.0370; f = 0.00717;

	} else if (type == "sheet_primitive") {
		if (t < 0) {t = 0;}
		else if (t > 3.05) {t = 3.05;}

		a = -0.00488; b = 0.236; c = -0.140; d = 0.150	; e = -0.0652; f = 0.0106;

	} else if (type == "IWP") {
		if (t < -5.025) {t = -5.025;}
		else if (t > 3.025) {t = 3.025;}

		a = 0.339; b = -0.0668; c = 0.000523; d = -0.000707; e = -0.000660;
		f = -0.000132;

	} else if (type == "sheet_IWP") { // VERY ROUGH APPROXIMATION!!!
		if (t < 0) {t = 0;}
		else if (t > 5.05) {t = 5.05;}

		a = 0.0212; b = -0.149; c = 0.534; d = -0.518; e = 0.226; f = -0.0437;
		g = 0.00310;

	} else {
		std::cerr << "ERROR_INVALID_TPMS" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Compute the wall thickness
	double wallSize = (a + b*t + c*std::pow(t,2) + 
	                             d*std::pow(t,3) + 
	                             e*std::pow(t,4) + 
	                             f*std::pow(t,5) + 
	                             g*std::pow(t,6) + 
	                             h*std::pow(t,7)) * 
	                  (2*PI / scaling);

	if (wallSize < 0) { return 0; }
	else { return wallSize; }
}

double Infill::internal::level2poreSize(double t, double scaling, std::string type) {
	/* Estimates the pore size that corresponds to the requested level-set
	 * constant. Estimation avaliable for: gyroid, diamond, primitive and IWP in
	 * both strut and sheet variants.
	 * Inputs:
	 *   t       : Level-set constant
	 *   scaling : Unit cell scaling
	 *   type    : Unit cell type
	 *  Return:
	 *   Pore size
	 */
	
	// Polinomial coeficients
	double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;

	// Enforce the level-set constant bounds and set the polynomial coefficients
	// (for the requested unit cell type)
	if (type == "gyroid") {
		if (t < -1.5) {t = -1.5;}
		else if (t > 1.5) {t = 1.5;}
		
		a = 0.428; b = -0.191; c = 0.00166; d = 0.00822; e = -0.00133; f = -0.0169;

	} else if (type == "sheet_gyroid") {
		if (t < 0) {t = 0;}
		else if (t > 1.55) {t = 1.55;}

		a = 0.432; b = -0.264; c = 0.418; d = -0.893; e = 0.782; f = -0.253;

	} else if (type == "diamond") {
		if (t < -1.425) {t = -1.425;}
		else if (t > 1.425) {t = 1.425;}

		a = 0.395; b = -0.208; c = 0.0180; d = 0.0144; e = -0.00916; f = -0.0251;

	} else if (type == "sheet_diamond") {
		if (t < 0) {t = 0;}
		else if (t > 1.45) {t = 1.45;}

		a = 0.421; b = -0.288; c = 0.568; d = -1.26; e = 1.16; f = -0.397;

	} else if (type == "primitive") {
		if (t < -3) {t = -3;}
		else if (t > 3) {t = 3;}

		a = 0.843; b = -0.193; c = 0.00123; d = 0.00326; e = -0.000306;
		f = -0.00123;
		
	} else if (type == "sheet_primitive") {
		if (t < 0) {t = 0;}
		else if (t > 3.05) {t = 3.05;}

		a = 0.847; b = -0.225; c = 0.113; d = -0.125; e = 0.0558; f = -0.00943;

	} else if (type == "IWP") {
		if (t < -5.025) {t = -5.025;}
		else if (t > 3.025) {t = 3.025;}

		a = 0.426; b = 0.0662; c = 0.000860; d = 0.000956; e = 0.000764;
		f = 0.000146;

	} else if (type == "sheet_IWP") {
		if (t < 0) {t = 0;}
		else if (t > 5.05) {t = 5.05;}

		a = 0.428; b = -0.0906; c = 0.0434; d = -0.0258; e = 0.00663;
		f = -0.000627;

	} else {
		std::cerr << "ERROR_INVALID_TPMS" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Compute the pore size
	double poreSize = (a + b*t + c*std::pow(t,2) + 
	                             d*std::pow(t,3) + 
	                             e*std::pow(t,4) + 
	                             f*std::pow(t,5) + 
	                             g*std::pow(t,6) + 
	                             h*std::pow(t,7)) * 
	                  (2*PI / scaling);

	if (poreSize < 0) { return 0; }
	else { return poreSize; }
}

double Infill::internal::eModulus2vFraction(double eModulus, 
                                          latticeMaterial materialParameters) {
	/* Converts the elastic modulus to a volume fraction (relative density)
	 * according to the material model selected. The Gibson & Ashby equations are
	 * provided with the code. Users are strongly adviced to use their own 
	 * equations.
	 * Inputs:
	 *   eModulus           : Elastic modulus to be converted
	 *   materialParameters : Material model parameters
	 * Return:
	 *   Volume fraction
	 */

	double relativeModulus = eModulus/materialParameters.eModulusSolid; // Relative density
	if (materialParameters.materialModel == "gibson_ashby") { // User defined Gibson-Ashby based model
		if (relativeModulus < 0) {
			return 0;
		} else if (relativeModulus > 1) {
			return 1;
		} else {
			return std::pow( (eModulus/(materialParameters.C*materialParameters.eModulusSolid)), (1.0/materialParameters.n) );
		}

	} else if (materialParameters.materialModel == "GyroidYan2015") { // Yan2015|Gibson-Ashby model
		if (relativeModulus < 0) {
			return 0;
		} else if (relativeModulus <= 0.01) { // Yan2015 model for a Gyroid TPMS (80-95 porosity)
			return std::pow( (eModulus/(0.19*materialParameters.eModulusSolid)), (1.0/1.71) );
		} else if (relativeModulus > 1) {
			return 1;
		} else {
			return std::sqrt(eModulus/materialParameters.eModulusSolid); // Gibson-Ashby model
		}

	} else {
		std::cerr << "ERROR_INVALID_MATERIAL_MODEL_REQUEST" << std::endl;
		exit(EXIT_FAILURE);
	}
}


double Infill::TPMS_function(Point p, std::string type, double scaling, double t) {
	/* Returns the signed distance at point p for the requested unit cell type.
	 * Currently available for: gyroid, diamond, primitive and IWP in both strut
	 * and sheet variants.
	 * Inputs:
	 *   p       : Coordinates of the point to be evaluated
	 *   type    : Unit cell type (at the point)
	 *   scaling : Unit cell scaling (at the point)
	 *   t       : Level-set constant (at the point)
	 * Return:
	 *   Signed distance for the chosen unit cell at the point
	 */

	if (type == "gyroid") {
		return sin(scaling*p.x())*cos(scaling*p.y()) +
		       sin(scaling*p.y())*cos(scaling*p.z()) +
		       sin(scaling*p.z())*cos(scaling*p.x()) -
		       t;

	} else if (type == "sheet_gyroid") {
		return std::pow((sin(scaling*p.x())*cos(scaling*p.y()) +
		                 sin(scaling*p.y())*cos(scaling*p.z()) +
		                 sin(scaling*p.z())*cos(scaling*p.x())), 2) -
		       std::pow(t, 2); 

	} else if (type == "diamond") {
		return sin(scaling*p.x())*sin(scaling*p.y())*sin(scaling*p.z()) +
		       sin(scaling*p.x())*cos(scaling*p.y())*cos(scaling*p.z()) +
		       cos(scaling*p.x())*sin(scaling*p.y())*cos(scaling*p.z()) +
		       cos(scaling*p.x())*cos(scaling*p.y())*sin(scaling*p.z()) -
		       t;

	} else if (type == "sheet_diamond") {
		return std::pow(sin(scaling*p.x())*sin(scaling*p.y())*sin(scaling*p.z()) +
		                sin(scaling*p.x())*cos(scaling*p.y())*cos(scaling*p.z()) +
		                cos(scaling*p.x())*sin(scaling*p.y())*cos(scaling*p.z()) +
		                cos(scaling*p.x())*cos(scaling*p.y())*sin(scaling*p.z()), 2) -
		       std::pow(t, 2);
	
	} else if (type == "primitive") {
		return -(cos(scaling*p.x()) + cos(scaling*p.y()) + cos(scaling*p.z())) - t;

	} else if (type == "sheet_primitive") {
		return std::pow(-(cos(scaling*p.x()) +
		                  cos(scaling*p.y()) +
		                  cos(scaling*p.z())), 2) -
		       std::pow(t, 2);

	} else if (type == "IWP") {
		return (cos(2*scaling*p.x())+cos(2*scaling*p.y())+cos(2*scaling*p.z())) -
		       2*(cos(scaling*p.x())*cos(scaling*p.y()) +
		          cos(scaling*p.y())*cos(scaling*p.z()) +
		          cos(scaling*p.z())*cos(scaling*p.x())) +
		       t;

	} else if (type == "sheet_IWP") {
		return std::pow((-cos(2*scaling*p.x()) - 
		                  cos(2*scaling*p.y()) - 
		                  cos(2*scaling*p.z())) +
		       				  2*(cos(scaling*p.x())*cos(scaling*p.y()) +
		                   cos(scaling*p.y())*cos(scaling*p.z()) +
		                   cos(scaling*p.z())*cos(scaling*p.x())), 2) -
		       std::pow(t, 2);

	} else { 
		std::cerr << "ERROR_INVALID_TPMS" << std::endl;
		exit(EXIT_FAILURE);
	}
}

double Infill::TPMS_function(Point p, latticeType *lt_type, latticeSize *lt_size,
                           latticeFeature *lt_feature) {
	/* Returns the signed distance at point p for a lattice of requested 
	 * properties. Currently available for: gyroid, diamond, primitive and 
	 * IWP in both strut and sheet variants.
	 * Inputs:
	 *   p          : Coordinates of point to be evaluated
	 *   lt_type    : Lattice type data
	 *   lt_size    : Lattice size data
	 *   lt_feature : Lattice feature data
	 * Return:
	 *   Signed distance at point p
	 */

	Eigen::Vector3d pEigen;
	double scaling, featureVal, signedDistance;

	// Process point p
	pEigen = Eigen::Vector3d(p.x(), p.y(), p.z());

	// Determine unit cell scaling
	scaling = sizing_function(p, lt_size, "scaling");

	// Determine the feature value
	if (lt_feature->feature_val > 0) // If feature is provided as constant
		featureVal = lt_feature->feature_val;
	else // If feature is provided as variable
		featureVal = TrilinearInterpolation::evaluate(pEigen, 
		                                        &(lt_feature->interpModel_linear));

	// Detemine the typing and compute the signed distance at point p
	if (lt_type->type != "hybrid") { // If unit cell is fixed
		double t = internal::input2level(lt_type->type, scaling, lt_feature->feature, featureVal, 
		                                 lt_feature->mp, lt_feature->mode);
		signedDistance = TPMS_function(p, lt_type->type, scaling, t);

	} else { // If unit cell is variable
		// Compute weights
		std::vector<double> weights(lt_type->interpModel_linear.size());
		double norm = 0;
		for (size_t i = 0; i < lt_type->interpModel_linear.size(); i++) {
			weights[i] = TrilinearInterpolation::evaluate(pEigen, 
			                                      &(lt_type->interpModel_linear[i]));
			norm += weights[i] * weights[i];
		}
		norm = sqrt(norm);

		// Determine the weighted signed distance at point p
		signedDistance = 0;
		for (size_t i = 0; i < lt_type->interpModel_linear.size(); i++) {
			// Compute the level correction
			double correction = 1 + lt_type->correctionFactor * (1-std::pow(2*(weights[i]/norm)-1 ,2));

			// Compute the corrected feature value (Compensates thinning in hybrid region)
			featureVal = correction * featureVal;

			// Compute the weighted signed distance
			double t = internal::input2level(lt_type->typeVector[i], scaling, lt_feature->feature, featureVal, lt_feature->mp, lt_feature->mode);
			signedDistance += weights[i]/norm * TPMS_function(p, lt_type->typeVector[i], scaling, t);
		}
	}

	// (Weighted) signed distance at point p
	return signedDistance;
}

double Infill::sizing_function(Point p, latticeSize *lt_size, std::string mode) {
	/* Returns the size or scaling at point p.
	 * Inputs:
	 *   p       : Coordinates of point to be evaluated
	 *   lt_size : Lattice size data
	 *   mode    : Value to be returned: "scaling" or "size"
	 * Returns:
	 *   Size or scaling at point p
	 */

	double scaling;

	if (lt_size->size > 0) { // If size is provided as constant
		scaling = lt_size->scaling;

	} else { // If size is provided as variable
		Eigen::Vector3d pEigen = Eigen::Vector3d(p.x(), p.y(), p.z());
		scaling = TrilinearInterpolation::evaluate(pEigen, 
		                                           &(lt_size->interpModel_linear));
	}

	if (mode == "scaling") // Return scaling
		return scaling;
	else                   // Return size
		return 2*PI/scaling;
}

featureSize Infill::featureSize_function(Point p, latticeType *lt_type,
                                       latticeSize *lt_size, 
                                       latticeFeature *lt_feature) {
	/* Returns the feature size at point p for a lattice of requested 
	 * properties. Currently available for: gyroid, diamond, primitive and 
	 * IWP in both strut and sheet variants.
	 * Inputs:
	 *   p          : Coordinates of point to be evaluated
	 *   lt_type    : Lattice type data
	 *   lt_size    : Lattice size data
	 *   lt_feature : Lattice feature data
	 * Return:
	 *   Wall and pore size at point p
	 */

	Eigen::Vector3d pEigen;
	double scaling, featureVal;

	featureSize output;
		output.wallSize = HUGE_VAL;
		output.poreSize = HUGE_VAL;

	// Process point p
	pEigen = Eigen::Vector3d(p.x(), p.y(), p.z());

	// Determine unit cell scaling
	scaling = sizing_function(p, lt_size, "scaling");

	// Determine the feature value
	if (lt_feature->feature_val > 0) // If feature is provided as constant
		featureVal = lt_feature->feature_val;
	else // If feature is provided as variable
		featureVal = TrilinearInterpolation::evaluate(pEigen, 
		                                        &(lt_feature->interpModel_linear));

	// Determine wall and pore size
	if (lt_type->type != "hybrid") { // If unit cell is fixed
		// Compute the level-set constant
		double t = internal::input2level(lt_type->type, scaling, lt_feature->feature, featureVal, 
		                                 lt_feature->mp, lt_feature->mode);
	
		// Compute wall and pore size
		output.wallSize = internal::level2wallSize(t, scaling, lt_type->type);
		output.poreSize = internal::level2poreSize(t, scaling, lt_type->type);

	} else { // If unit cell is variable
		// Compute  weights
		std::vector<double> weights(lt_type->interpModel_linear.size());
		double norm = 0;
		for (size_t i = 0; i < lt_type->interpModel_linear.size(); i++) {
			weights[i] = TrilinearInterpolation::evaluate(pEigen, &(lt_type->interpModel_linear[i]));
			norm += weights[i] * weights[i];
		}
		norm = sqrt(norm);

		// Compute wall and pore size
		for (size_t i = 0; i < lt_type->interpModel_linear.size(); i++) {
			double t = internal::input2level(lt_type->typeVector[i], scaling, 
			                                lt_feature->feature, featureVal,
			                                lt_feature->mp, lt_feature->mode);
			
			if (weights[i]/norm > 0.1) {
				output.wallSize = std::min(output.wallSize, internal::level2wallSize(t, scaling, lt_type->typeVector[i]));
				output.poreSize = std::min(output.poreSize, internal::level2poreSize(t, scaling, lt_type->typeVector[i]));
			}
		}
	}

	return output;
}