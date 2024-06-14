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

#include "Infill.h"

/* Infill contains the functions required to create lattice infills based on 
 * triply periodic minimal surfaces. 
 * 
 * Note: Strictly speaking Infill can handle infills described by any signed 
 * distance function, not only triply periodic minimal surface based infills.
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

double Infill::internal::input2level(const std::string &type, const double &scaling, const std::string &feature, 
                                     const double &featureValue, const latticeUDF &userDefinedParameters,
                                     const std::string &featureMode) {
	/* Estimates the level-set constant that corresponds to the provided input
	 * parameter. Estimation avaliable for: gyroid, diamond, primitive and IWP
	 * in both strut and sheet variants.
	 * Inputs:
	 *   type                  : Unit cell type
	 *   scaling               : Unit cell scaling
	 *   feature               : Feature to be converted
	 *   featureValue          : Feature value
	 *   userDefinedParameters : User defined feature parameters (Used only 
	 *                           when making use of a user defined feature)
	 *   featureMode           : 
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
	} else if (feature == "userDefined") {
		return vFraction2level(userDefinedInput2vFraction(featureValue, userDefinedParameters), type);
	} else {
		throw std::runtime_error(INFILL_ERRMSG::INVALID_FEATURE_REQUEST);
	}
}

double Infill::internal::unnormalizeLevel(const double &t_normalized, const std::string &type) {
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
		return (t_normalized * (1.5 - 0)) + 0;

	} else if (type == "diamond") {
		return (t_normalized * (1.42 - -1.42)) + -1.42;
		
	} else if (type == "sheet_diamond") {
		return (t_normalized * (1.42 - 0)) + 0;

	} else if (type == "primitive") {
		return (t_normalized * (3.01 - -3)) + -3;
		
	} else if (type == "sheet_primitive") {
		return (t_normalized * (3.01 - 0)) + 0;

	} else if (type == "IWP") {
		return (t_normalized * (5 - -3)) + -3;

	} else if (type == "sheet_IWP") {
		return (t_normalized * (5 - 0)) + 0;

	} else if (type == "cubic") { // Without smoothing
		return (t_normalized * (0.71 - 0)) + 0;

	} else {
		throw std::runtime_error(INFILL_ERRMSG::INVALID_TPMS);
	}
}

double Infill::internal::vFraction2level(double volFraction, const std::string &type) {
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
		a = 0; b = 1.5486	; c = -0.0486;

	} else if (type == "diamond") {
		a = -1.42; b = 9.5121; c = -74.9721; d = 364.3329; e = -930.7729; f = 1291.0637;
		g = -919.0071; h = 262.6834;

	} else if (type == "sheet_diamond") {
		a = 0; b = 0.7912; c = 6.5983; d = -33.7214; e = 75.7986; f = -77.7199;
		g = 29.6732;

	} else if (type == "primitive") {
		a = -3; b = 21.9305; c = -96.7125; d = 227.6736; e = -244.8562; f = 97.9746;

	} else if (type == "sheet_primitive") {
		a = 0; b = 3.2223; c = -29.3794; d = 200.6042; e = -633.4513; f = 1009.5031;
		g = -788.9678; h = 241.4790;

	} else if (type == "IWP") {
		a = -3; b = 5.8392; c = -59.4099; d = 462.7837; e = -1495.7952; f = 2430.5479;
		g = -1945.4597; h = 609.4940;

	} else if (type == "sheet_IWP") {
		a = 0; b = 6.3703; c = -53.5909; d = 383.7502; e = -1279.0549; f = 2160.0684;
		g = -1795.0308; h = 582.4877;

	} else if (type == "cubic") { // Without smoothing
		a = 0; b = 2.8402; c = -30.1344; d = 173.9104; e = -509.3130; f = 782.2430;
		g = -600.4298; h = 181.5936;

	} else {
		throw std::runtime_error(INFILL_ERRMSG::INVALID_TPMS);
	}

	// Compute the level-set constant
	return a + b*volFraction + c*std::pow(volFraction,2) + 
	                           d*std::pow(volFraction,3) + 
	                           e*std::pow(volFraction,4) + 
	                           f*std::pow(volFraction,5) + 
	                           g*std::pow(volFraction,6) + 
	                           h*std::pow(volFraction,7);
}

double Infill::internal::wallSize2level(double wallSize, const double &scaling, const std::string &type, const std::string &mode) {
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
		a = -1.4292; b = 1.1123; c = 2.0586; d = 22.7985; e = -42.2134; f = 19.2938;
		if (wallSize > 0.7761) {wallSize = 0.7761;}

	} else if (type == "sheet_gyroid") {
		a = -0.0018; b = 5.7957; c = -7.9287; d = 47.1949; e = -172.5956; f = 175.2846;
		if (wallSize > 0.4769) {wallSize = 0.4769;}

	} else if (type == "diamond") {
		a = -1.0143; b = 1.2004; c = 9.8086; d = -9.5197;
		if (wallSize > 0.6864) {wallSize = 0.6864;}

	} else if (type == "sheet_diamond") {
		a = 0; b = 5.8008; c = -6.0142;
		if (wallSize > 0.4365) {wallSize = 0.4365;}

	} else if (type == "primitive") {
		a = -0.9970; b = 0.1082; c = 6.3135; d = -3.1200;
		if (wallSize > 1.2202) {wallSize = 1.2202;}

	} else if (type == "sheet_primitive") {
		a = -0.0007; b = 5.5064; c = -0.4280; d = -2.2617;
		if (wallSize > 0.7541) {wallSize = 0.7541;}

	} else if (type == "IWP") {
		a = -3; b = -0.9792; c = 40.0955; d = -33.0667;
		if (wallSize > 0.7346) {wallSize = 0.7346;}

	} else if (type == "sheet_IWP") {
		a = 0; b = 12.9334; c = 58.1260; d =-454.0008; e = 947.9448; f = -626.0384;
		if (wallSize > 0.6415) {wallSize = 0.6415;}

	} else if (type == "cubic") { // Without smoothing
		a = 0; b = 0.9585; c = -3.1702; d = 10.1052; e = -16.9656; f = 15.5225;
		g = -7.3270; h = 1.3964;
		if (wallSize > 1.3565) {wallSize = 1.3565;}

	} else {
		throw std::runtime_error(INFILL_ERRMSG::INVALID_TPMS);
	}

	// Compute the level-set constant
	return (a + b*wallSize + c*std::pow(wallSize,2) + 
	                         d*std::pow(wallSize,3) + 
	                         e*std::pow(wallSize,4) + 
	                         f*std::pow(wallSize,5) + 
	                         g*std::pow(wallSize,6) + 
	                         h*std::pow(wallSize,7));
}

double Infill::internal::poreSize2level(double poreSize, const double &scaling, const std::string &type, const std::string &mode) {
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
		a = 1.4961; b = -0.1128; c = -9.7480; d = -0.6083; e = 14.2529; f = -6.5775;
		if (poreSize > 0.7761) {poreSize = 0.7761;}

	} else if (type == "sheet_gyroid") {
		a = 1.4886; b = 0.1637; c = -12.3632; d = 8.9919;
		if (poreSize > 0.4212) {poreSize = 0.4212;}

	} else if (type == "diamond") {
		a = 1.3822; b = -0.3794; c = -11.5654; d = 9.6452;
		if (poreSize > 0.7601) {poreSize = 0.7601;}

	} else if (type == "sheet_diamond") {
		a = 1.3939; b = -0.4082; c = -5.9832; d = -9.7175; e = 15.3286; f = 4.1403;
		if (poreSize > 0.4097) {poreSize = 0.4097;}

	} else if (type == "primitive") {
		a = 2.9258; b = 0.3522; c = -6.7912; d = 2.6777;
		if (poreSize > 1.5391) {poreSize = 1.5391;}

	} else if (type == "sheet_primitive") {
		a = 2.9750; b = -0.0967; c = -5.7828; d = 2.0385;
		if (poreSize > 0.8309) {poreSize = 0.8309;}

	} else if (type == "IWP") {
		a = 5.0085; b = -1.1781; c = -40.9767; d = 38.2058;
		if (poreSize > 0.7098) {poreSize = 0.7098;}

	} else if (type == "sheet_IWP") {
		a = 4.9112; b = 1.1684; c = -51.6997; d = 51.0962;
		if (poreSize > 0.4164) {poreSize = 0.4164;}


	} else if (type == "cubic") { // Without smoothing
		a = 0.6945; b = -0.4706; c = -0.2586; d = 0.6325; e = -0.5959; f = 0.1914;
		if (poreSize > 1.7093) {poreSize = 1.7093;}

	} else {
		throw std::runtime_error(INFILL_ERRMSG::INVALID_TPMS);
	}

	// Compute the level-set constant
	return (a + b*poreSize + c*std::pow(poreSize,2) + 
	                         d*std::pow(poreSize,3) + 
	                         e*std::pow(poreSize,4) + 
	                         f*std::pow(poreSize,5) + 
	                         g*std::pow(poreSize,6) + 
	                         h*std::pow(poreSize,7));
}


double Infill::internal::level2wallSize(double t, const double &scaling, const std::string &type) {
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

		a = 0.4133; b = 0.1942; c = 0.0061; d = -0.0126; e = -0.0104; f = 0.0237;

	} else if (type == "sheet_gyroid") {
		if (t < 0) {t = 0;}
		else if (t > 1.5) {t = 1.5;}

		a = -0.0048; b = 0.3776; c = -2.0275; d = 8.8079; e = -18.8187; f = 20.9671;
 		g = -11.6445; h = 2.5463;

	} else if (type == "diamond") {
		if (t < -1.42) {t = -1.42;}
		else if (t > 1.42) {t = 1.42;}

		a = 0.3096; b = 0.2034; c = -0.0336; d = 0.0659;

	} else if (type == "sheet_diamond") {
		if (t < 0) {t = 0;}
		else if (t > 1.42) {t = 1.42;}

		a = 0; b = 0.2118; c = -0.1328; d = 0.1494;

	} else if (type == "primitive") {
		if (t < -3) {t = -3;}
		else if (t > 3.01) {t = 3.01;}

		a = 0.4445; b = 0.2437; c = -0.0617; d = 0.0858; e = -0.0454; f = 0.0084;

	} else if (type == "sheet_primitive") {
		if (t < 0) {t = 0;}
		else if (t > 3.01) {t = 3.01;}

		a = -0.0028; b = 0.2320; c = -0.1465; d = 0.1641; e = -0.0731; f = 0.0120;

	} else if (type == "IWP") {
		if (t < -3) {t = -3;}
		else if (t > 5) {t = 5;}

		a = 3.4176e-1; b = 7.1522e-2; c = -5.7481e-3; d = -1.5987e-3; e = 1.4963e-3; f = 2.0003e-4;
		g = -1.9120e-4; h = 2.3700e-5;

	} else if (type == "sheet_IWP") {
		if (t < 0) {t = 0;}
		else if (t > 5) {t = 5;}

		a = 0; b = -0.0218; c = 0.2852; d = -0.2908; e = 0.1171; f = -0.0157;
		g = -5.8286e-4; h=1.9297e-4;






	} else if (type == "cubic") { // Without smoothing
		if (t < 0) {t = 0;}
		else if (t > 0.71) {t = 0.71;}

		a = 0; b = 0.9949; c = 7.4779; d = -23.7118; e = 33.4551; f = -17.1937;

	} else {
		throw std::runtime_error(INFILL_ERRMSG::INVALID_TPMS);
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

double Infill::internal::level2poreSize(double t, const double &scaling, const std::string &type) {
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
		
		a = 0.4276; b = -0.1906; c = 0.0005; d = 0.0072; e = -0.0006; f = -0.0166;

	} else if (type == "sheet_gyroid") {
		if (t < 0) {t = 0;}
		else if (t > 1.5) {t = 1.5;}

		a = 0.4295; b = -0.2419; c = 0.3462; d = -0.7699; e = 0.6843; f = -0.2224;

	} else if (type == "diamond") {
		if (t < -1.42) {t = -1.42;}
		else if (t > 1.42) {t = 1.42;}

		a = 0.3955; b = -0.2086; c = 0.0175; d = 0.0164; e = -0.0086; f = -0.0265;

	} else if (type == "sheet_diamond") {
		if (t < 0) {t = 0;}
		else if (t > 1.42) {t = 1.42;}

		a = 0.4181; b = -0.2666; c = 0.5460; d = -1.3124; e = 1.2599; f = -0.4374;

	} else if (type == "primitive") {
		if (t < -3) {t = -3;}
		else if (t > 3.01) {t = 3.01;}

		a = 0.8432; b = -0.1902; c = 0.0003; d = 0.0016; e = -0.0002; f = -0.0010;
		
	} else if (type == "sheet_primitive") {
		if (t < 0) {t = 0;}
		else if (t > 3.01) {t = 3.01;}

		a = 0.8469; b = -0.2328; c = 0.1407; d = -0.1563; e = 0.0697; f = -0.0116;

	} else if (type == "IWP") {
		if (t < -3) {t = -3;}
		else if (t > 5) {t = 5;}

		a = 4.2556e-1; b = -6.6094e-2; c = 1.8489e-3; d = -1.1156e-3; e = 6.2127e-4; 
		f = -1.1685e-4;

	} else if (type == "sheet_IWP") {
		if (t < 0) {t = 0;}
		else if (t > 5) {t = 5;}

		a = 0.4264; b = -0.0862; c = 0.0395; d = -0.0245; e =0.0064; f = -6.1871e-4;


	} else if (type == "cubic") { // Without smoothing
		if (t < 0) {t = 0;}
		else if (t > 0.71) {t = 0.71;}
		
		a = 1.7093, b = -13.0874, c = 123.2556, d = -631.6172, e = 1690.8477, f = -2427.3976,
		g = 1753.0583, h = -492.3583;

	} else {
		throw std::runtime_error(INFILL_ERRMSG::INVALID_TPMS);
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

double Infill::internal::userDefinedInput2vFraction(const double &userDefinedInput, 
                                                    const latticeUDF &userDefinedParameters) {
	/* Converts a user defined input to a volume fraction (relative density)
	 * according to the user defined model selected. Users can implement 
	 * the equation(s) that convert from their desired inputs to volume fractions
	 * by adding an else if statement bellow that contains said conversion 
	 * equation(s).
	 * Inputs:
	 *   userDefinedInput      : User defined input to be converted, i.e. 
	 *                           interpolated FAP file value at point being
	 *                           evaluated
	 *   userDefinedParameters : Strucutre containing the user defined model 
	 *                           parameters
	 * Return:
	 *   Volume fraction
	 */

	if (userDefinedParameters.userDefinedFeature == "elasticModulus") { // Gibson-Ashby based model
		double relativeModulus = userDefinedInput/userDefinedParameters.A; // Relative density
		if (relativeModulus < 0) {
			return 0;
		} else if (relativeModulus > 1) {
			return 1;
		} else {
			return std::pow( (userDefinedInput/(userDefinedParameters.B*userDefinedParameters.A)), (1.0/userDefinedParameters.C) );
		}

	} else if (userDefinedParameters.userDefinedFeature == "stresses") {
		throw std::runtime_error(INFILL_ERRMSG::NO_STRESS_CONVERSION_DEFINED);
		
	} else {
		throw std::runtime_error(INFILL_ERRMSG::INVALID_USER_DEFINED_FEATURE_REQUEST);
	}
}


double Infill::TPMS_function(const Point &p, const std::string &type,
	const double &scaling, const double &t) {
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
		          cos(scaling*p.z())*cos(scaling*p.x())) -
		       t;

	} else if (type == "sheet_IWP") {
		return std::pow((-cos(2*scaling*p.x()) - 
		                  cos(2*scaling*p.y()) - 
		                  cos(2*scaling*p.z())) +
		       				  2*(cos(scaling*p.x())*cos(scaling*p.y()) +
		                   cos(scaling*p.y())*cos(scaling*p.z()) +
		                   cos(scaling*p.z())*cos(scaling*p.x())), 2) -
		       std::pow(t, 2);

	} else if (type == "cubic") { // Non-TPMS unit cell
		bool smooth = false; // If changed all other equations for the cubic unit cell need to be adjusted!
		double smoothingFactor = 2.5; // If changed all other equations for the cubic unit cell need to be adjusted!

		// Period (2PI is used to make definition similar to that of TPMS unit cells implemented)
		double period = 2*PI;

		// Construct triply periodic (TP) unit cell from cylinders
		double x = std::remainder(p.x(), (1/scaling)*period);
		double y = std::remainder(p.y(), (1/scaling)*period);
		double z = std::remainder(p.z(), (1/scaling)*period);

		double X1 = (std::sqrt(std::pow(x, 2) + std::pow(y, 2)) / t) / (2*PI/scaling);
		double X2 = (std::sqrt(std::pow(y, 2) + std::pow(z, 2)) / t) / (2*PI/scaling);
		double X3 = (std::sqrt(std::pow(x, 2) + std::pow(z, 2)) / t) / (2*PI/scaling);

		// Intersection  smoothing
		double smoothing = 0;
		if ( smooth == true )
			smoothing = std::pow( std::max({smoothingFactor-std::abs(X1-X2)-std::abs(X2-X3)-std::abs(X1-X3), 0.0}), 3 ) / (6*std::pow(smoothingFactor, 2));

		return std::min({X1, X2, X3}) - 1 - smoothing; // Boolean union of individual cylinders (with smoothing)

	} else { 
		throw std::runtime_error(INFILL_ERRMSG::INVALID_TPMS);
	}
}

double Infill::TPMS_function(const Point &p, const latticeType &lt_type,
	const latticeSize &lt_size, const latticeFeature &lt_feature) {
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
	if (lt_feature.feature_val > 0) // If feature is provided as constant
		featureVal = lt_feature.feature_val;
	else // If feature is provided as variable
		featureVal = TrilinearInterpolation::evaluate(pEigen, 
		                                        lt_feature.interpModel_linear);

	// Detemine the typing and compute the signed distance at point p
	if (lt_type.type != "hybrid") { // If unit cell is fixed
		double t = internal::input2level(lt_type.type, scaling, lt_feature.feature, featureVal, 
		                                 lt_feature.udf, lt_feature.mode);
		signedDistance = TPMS_function(p, lt_type.type, scaling, t);

	} else { // If unit cell is variable
		// Compute weights
		std::vector<double> weights(lt_type.interpModel_linear.size());
		double norm = 0;
		for (size_t i = 0; i < lt_type.interpModel_linear.size(); i++) {
			weights[i] = TrilinearInterpolation::evaluate(pEigen, 
			                                      lt_type.interpModel_linear[i]);
			norm += weights[i] * weights[i];
		}
		norm = sqrt(norm);

		// Determine the weighted signed distance at point p
		signedDistance = 0;
		for (size_t i = 0; i < lt_type.interpModel_linear.size(); i++) {
			// Compute the level correction
			double correction = 1 + lt_type.correctionFactor * (1-std::pow(2*(weights[i]/norm)-1 ,2));

			// Compute the corrected feature value (Compensates thinning in hybrid region)
			featureVal = correction * featureVal;

			// Compute the weighted signed distance
			double t = internal::input2level(lt_type.typeVector[i], scaling, lt_feature.feature, featureVal, lt_feature.udf, lt_feature.mode);
			signedDistance += weights[i]/norm * TPMS_function(p, lt_type.typeVector[i], scaling, t);
		}
	}

	// (Weighted) signed distance at point p
	return signedDistance;
}

double Infill::sizing_function(const Point &p, const latticeSize &lt_size, const std::string &mode) {
	/* Returns the size or scaling at point p.
	 * Inputs:
	 *   p       : Coordinates of point to be evaluated
	 *   lt_size : Lattice size data
	 *   mode    : Value to be returned: "scaling" or "size"
	 * Returns:
	 *   Size or scaling at point p
	 */

	double scaling;

	if (lt_size.size > 0) { // If size is provided as constant
		scaling = lt_size.scaling;

	} else { // If size is provided as variable
		Eigen::Vector3d pEigen = Eigen::Vector3d(p.x(), p.y(), p.z());
		scaling = TrilinearInterpolation::evaluate(pEigen, 
		                                           lt_size.interpModel_linear);
	}

	if (mode == "scaling") // Return scaling
		return scaling;
	else                   // Return size
		return 2*PI/scaling;
}

featureSize Infill::featureSize_function(const Point &p, const latticeType &lt_type,
	const latticeSize &lt_size, const latticeFeature &lt_feature) {
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
	if (lt_feature.feature_val > 0) // If feature is provided as constant
		featureVal = lt_feature.feature_val;
	else // If feature is provided as variable
		featureVal = TrilinearInterpolation::evaluate(pEigen, lt_feature.interpModel_linear);

	// Determine wall and pore size
	if (lt_type.type != "hybrid") { // If unit cell is fixed
		// Compute the level-set constant
		double t = internal::input2level(lt_type.type, scaling, lt_feature.feature, featureVal, 
		                                 lt_feature.udf, lt_feature.mode);
	
		// Compute wall and pore size
		output.wallSize = internal::level2wallSize(t, scaling, lt_type.type);
		output.poreSize = internal::level2poreSize(t, scaling, lt_type.type);

	} else { // If unit cell is variable
		// Compute  weights
		std::vector<double> weights(lt_type.interpModel_linear.size());
		double norm = 0;
		for (size_t i = 0; i < lt_type.interpModel_linear.size(); i++) {
			weights[i] = TrilinearInterpolation::evaluate(pEigen, lt_type.interpModel_linear[i]);
			norm += weights[i] * weights[i];
		}
		norm = sqrt(norm);

		// Compute wall and pore size
		for (size_t i = 0; i < lt_type.interpModel_linear.size(); i++) {
			double t = internal::input2level(lt_type.typeVector[i], scaling, 
			                                lt_feature.feature, featureVal,
			                                lt_feature.udf, lt_feature.mode);
			
			if (weights[i]/norm > 0.1) {
				output.wallSize = std::min(output.wallSize, internal::level2wallSize(t, scaling, lt_type.typeVector[i]));
				output.poreSize = std::min(output.poreSize, internal::level2poreSize(t, scaling, lt_type.typeVector[i]));
			}
		}
	}

	return output;
}