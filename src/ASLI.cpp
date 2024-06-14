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

#include "ASLI.h"

/* The ASLI class setups all the data required to provide the requested object
 * with a (heterogeneous) lattice infill.
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

ASLI::ASLI(std::string configFile) {
	SetUp(configFile);
}

ASLI::~ASLI() {}

void ASLI::SetUp(std::string configFile) {
	/* Reads the data from the configuration file and subsequently calls the 
	 * necessary setup functions.
	 */

	YAML::Node config;
	std::set<std::string> accessed_keys;

	// Load configuration file
	std::cout << "  Opening " << configFile << std::endl;
	config = YAML::LoadFile(configFile);
		
	// Lattice type
	auto type = getOptional<std::string>(config, "lt_type", true, accessed_keys);
	if (type.has_value()) lt_type.type = *type;
	else throw std::runtime_error(ASLI_ERRMSG::INVALID_TPMS);
	if (Infill::UNIT_CELLS.find(lt_type.type) == Infill::UNIT_CELLS.end()) throw std::runtime_error(ASLI_ERRMSG::INVALID_TPMS);

	auto side = getOptional<std::string>(config, "lt_type_side", false, accessed_keys);
	if (side.has_value()) lt_type.side = *side;
	if (SIDES.find(lt_type.side) == SIDES.end()) throw std::runtime_error(ASLI_ERRMSG::INVALID_SIDE);

	if (lt_type.type == "hybrid") {
		auto filterRadius = getOptional<double>(config, "lt_type_filterRadius", false, accessed_keys);
		if (filterRadius.has_value()) lt_type.filterRadius = *filterRadius;
		if (lt_type.filterRadius < 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_FILER_RADIUS);

		auto correctionFactor = getOptional<double>(config, "lt_type_correctionFactor", false, accessed_keys);
		if (correctionFactor.has_value()) lt_type.correctionFactor = *correctionFactor;
		if (lt_type.correctionFactor < 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_CORRECTION_FACTOR);
	}

	// Lattice size
	auto size = getOptional<double>(config, "lt_size", true, accessed_keys);
	if (size.has_value()) lt_size.size = *size;
	else throw std::runtime_error(ASLI_ERRMSG::INVALID_UNIT_CELL_SIZE);

	// Lattice features
	auto feature = getOptional<std::string>(config, "lt_feature", true, accessed_keys);
	if (feature.has_value()) lt_feature.feature = *feature;
	if (Infill::FEATURES.find(lt_feature.feature) == Infill::FEATURES.end()) throw std::runtime_error(ASLI_ERRMSG::INVALID_FEATURE);

	auto feature_val = getOptional<double>(config, "lt_feature_val", true, accessed_keys);
	if (feature_val.has_value()) lt_feature.feature_val = *feature_val;
	if (lt_feature.feature_val < 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_FEATURE_VALUE);

	if (FEATURES.find(lt_feature.feature) != FEATURES.end()) {
		auto mode = getOptional<std::string>(config, "lt_feature_mode", false, accessed_keys);
		if (mode.has_value()) lt_feature.mode = *mode;
		if (MODES.find(lt_feature.mode) == MODES.end()) throw std::runtime_error(ASLI_ERRMSG::INVALID_FEATURE_MODE);
	}

	// User defined feature parameters
	if ( lt_feature.feature == "userDefined" ) {
		auto userDefinedFeature = getOptional<std::string>(config, "udf_feature", false, accessed_keys);
		if (userDefinedFeature.has_value()) lt_feature.udf.userDefinedFeature = *userDefinedFeature;
		auto udf_A = getOptional<double>(config, "udf_A", false, accessed_keys);
		if (udf_A.has_value()) lt_feature.udf.A = *udf_A;
		auto udf_B = getOptional<double>(config, "udf_B", false, accessed_keys);
		if (udf_B.has_value()) lt_feature.udf.B = *udf_B;
		auto udf_C = getOptional<double>(config, "udf_C", false, accessed_keys);
		if (udf_C.has_value()) lt_feature.udf.C = *udf_C;
		auto udf_D = getOptional<double>(config, "udf_D", false, accessed_keys);
		if (udf_D.has_value()) lt_feature.udf.D = *udf_D;
		auto udf_E = getOptional<double>(config, "udf_E", false, accessed_keys);
		if (udf_E.has_value()) lt_feature.udf.E = *udf_E;
	}

	// Mesh settings
	auto mesher = getOptional<std::string>(config, "me_mesher", true, accessed_keys);
	if (mesher.has_value()) me_settings.mesher = *mesher;
	if (MESHERS.find(me_settings.mesher) == MESHERS.end()) throw std::runtime_error(ASLI_ERRMSG::INVALID_MESHER);

	auto isVolumeMesh = getOptional<bool>(config, "me_volumeMesh", false, accessed_keys);
	if (isVolumeMesh.has_value()) me_settings.isVolumeMesh = *isVolumeMesh;

	auto n_threads = getOptional<int>(config, "me_nThreads", false, accessed_keys);
	if (n_threads.has_value()) me_settings.n_threads = *n_threads;
	if (me_settings.n_threads < 1) throw std::runtime_error(ASLI_ERRMSG::INVALID_THREAD_NUMBER);

	auto elementSize = getOptional<double>(config, "me_elementSize", true, accessed_keys);
	if (elementSize.has_value()) me_settings.elementSize = *elementSize;
	if (me_settings.elementSize <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_MESH_SIZE);

	auto edgeProtectionAngle = getOptional<double>(config, "me_edgeProtectionAngle", false, accessed_keys);
	if (edgeProtectionAngle.has_value()) me_settings.edgeProtectionAngle = *edgeProtectionAngle;
	if (me_settings.edgeProtectionAngle < 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_EDGE_PROTECTION_ANGLE);
	
	auto STLFormat = getOptional<std::string>(config, "STLFormat", false, accessed_keys);
	if (STLFormat.has_value() && *STLFormat == "ASCII") me_settings.STLFormat = *STLFormat;

	auto threshold = getOptional<double>(config, "threshold", false, accessed_keys);
	if (threshold.has_value()) me_settings.threshold = *threshold;
	if (me_settings.threshold < 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_THRESHOLD);

	// CGAL specific settings
	if (me_settings.mesher == "CGAL" || me_settings.mesher == "CGAL_OLD") {
		me_settings.CGAL_cellSize = me_settings.elementSize;

		auto CGAL_facetDistance = getOptional<double>(config, "CGAL_facetDistance", false, accessed_keys);
		if (CGAL_facetDistance.has_value()) me_settings.CGAL_facetDistance = *CGAL_facetDistance;
		if (me_settings.CGAL_facetDistance <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_FACET_DISTANCE);

		auto CGAL_relativeErrorBound = getOptional<double>(config, "CGAL_relativeErrorBound", false, accessed_keys);
		if (CGAL_relativeErrorBound.has_value() ) me_settings.CGAL_relativeErrorBound = *CGAL_relativeErrorBound;
		if (me_settings.CGAL_relativeErrorBound <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_RELATIVE_ERROR_BOUND);

		auto CGAL_cellRadiusEdgeRatio = getOptional<double>(config, "CGAL_cellRadiusEdgeRatio", false, accessed_keys);
		if (CGAL_cellRadiusEdgeRatio.has_value() ) me_settings.CGAL_cellRadiusEdgeRatio = *CGAL_cellRadiusEdgeRatio;
		if (me_settings.CGAL_cellRadiusEdgeRatio <= 2) throw std::runtime_error(ASLI_ERRMSG::INVALID_CELL_RADIUS_EDGE_RATIO);

		// Aditional settings
		auto CGAL_facetAngle = getOptional<double>(config, "CGAL_facetAngle", false, accessed_keys);
		if (CGAL_facetAngle.has_value()) me_settings.CGAL_facetAngle = *CGAL_facetAngle;
		if (me_settings.CGAL_facetAngle <= 0 || me_settings.CGAL_facetAngle > 30) throw std::runtime_error(ASLI_ERRMSG::INVALID_FACET_ANGLE);

	} 
	if ( me_settings.mesher == "CGAL" || me_settings.mesher == "CGAL_OLD" ) {//Separated for now to keep segmentation fault bug caused by yaml-cpp from ocurring
		me_settings.CGAL_facetSize = me_settings.elementSize;
		auto CGAL_facetSize = getOptional<double>(config, "CGAL_facetSize", false, accessed_keys);
		if (CGAL_facetSize.has_value()) me_settings.CGAL_facetSize = *CGAL_facetSize;
		if (me_settings.CGAL_facetSize <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_FACET_SIZE);

		me_settings.CGAL_edgeSize = me_settings.elementSize;
		auto CGAL_edgeSize = getOptional<double>(config, "CGAL_edgeSize", false, accessed_keys);
		if (CGAL_edgeSize.has_value()) me_settings.CGAL_edgeSize = *CGAL_edgeSize;
		if (me_settings.CGAL_edgeSize <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_EDGE_SIZE);

		auto CGAL_minEdgeSize = getOptional<double>(config, "CGAL_minEdgeSize", false, accessed_keys);
		if (CGAL_minEdgeSize.has_value()) me_settings.CGAL_minEdgeSize = *CGAL_minEdgeSize;
		if (me_settings.CGAL_minEdgeSize < 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_EDGE_SIZE);

		auto CGAL_poissonOffset = getOptional<double>(config, "CGAL_poissonOffset", false, accessed_keys);
		if (CGAL_poissonOffset.has_value()) me_settings.CGAL_poissonOffset = *CGAL_poissonOffset;
		if (me_settings.CGAL_poissonOffset < 0.1) throw std::runtime_error(ASLI_ERRMSG::INVALID_OFFSET);

			// Set CGAL's workflow MMG settings
			me_settings.MMG_hinitial = 0.2;
			auto MMG_hinitial = getOptional<double>(config, "CGAL_hinitial", false, accessed_keys);
			if (MMG_hinitial.has_value()) me_settings.MMG_hinitial = *MMG_hinitial;
			if (me_settings.MMG_hinitial <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_HINITIAL);

			me_settings.MMG_hausd = 0.01;
			auto MMG_hausd = getOptional<double>(config, "CGAL_hausd", false, accessed_keys);
			if (MMG_hausd.has_value()) me_settings.MMG_hausd = *MMG_hausd;
			if (me_settings.MMG_hausd <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_HAUSD);
	}

	// MMG specific settings
	if (me_settings.mesher == "MMG") {
		me_settings.MMG_hsiz = me_settings.elementSize;

		auto MMG_hausd = getOptional<double>(config, "MMG_hausd", false, accessed_keys);
		if (MMG_hausd.has_value()) me_settings.MMG_hausd = *MMG_hausd;
		if (me_settings.MMG_hausd <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_HAUSD);
			
		auto MMG_hgrad = getOptional<double>(config, "MMG_hgrad", false, accessed_keys);
		if (MMG_hgrad.has_value()) me_settings.MMG_hgrad = *MMG_hgrad;
		if (me_settings.MMG_hgrad <= 1) throw std::runtime_error(ASLI_ERRMSG::INVALID_HGRAD);

		// Aditional settings
		auto MMG_hinitial = getOptional<double>(config, "MMG_hinitial", false, accessed_keys);
		if (MMG_hinitial.has_value()) me_settings.MMG_hinitial = *MMG_hinitial;
		if (me_settings.MMG_hinitial <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_HINITIAL);

		auto MMG_memory = getOptional<double>(config, "MMG_memory", false, accessed_keys);
		if (MMG_memory.has_value()) me_settings.MMG_memory = *MMG_memory;
		if (me_settings.MMG_memory < 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_MAX_MEMORY);

		auto MMG_exportLS = getOptional<bool>(config, "MMG_exportLS", false, accessed_keys);
		if (MMG_exportLS.has_value()) me_settings.MMG_exportLS = *MMG_exportLS;

			// Set MMG's workflow CGAL settings
			me_settings.CGAL_cellSize = me_settings.MMG_hinitial;
			me_settings.CGAL_facetDistance = 0.1; // Scales with local mesh size
			auto CGAL_facetDistance = getOptional<double>(config, "MMG_facetDistance", false, accessed_keys);
			if (CGAL_facetDistance.has_value()) me_settings.CGAL_facetDistance = *CGAL_facetDistance;
			if (me_settings.CGAL_facetDistance <= 0) throw std::runtime_error(ASLI_ERRMSG::INVALID_FACET_DISTANCE);
			me_settings.CGAL_cellRadiusEdgeRatio = 3.0;
			me_settings.CGAL_facetAngle = 30;
	}

	// Input and output files
	{	std::string fileKey = "files";
		YAML::Node model_files = config[fileKey]; 
		accessed_keys.insert(fileKey);

		// Input geometry
		auto inputFile = getOptional<std::string>(model_files, "stl", false, accessed_keys);
		if ( inputFile.has_value() ) 
			me_settings.inputFile = *inputFile;
		else {
			auto internalGeometry = getOptional<std::string>(model_files, "internalGeometry", false, accessed_keys);
			if (internalGeometry.has_value()) {
				me_settings.internalGeometry = *internalGeometry;

				std::vector<std::string> v = split (me_settings.internalGeometry, " ");
				me_settings.internalGeometry = v[0];
				for (int i = 1; i < v.size(); i++)
					me_settings.internalGeometryParameters.push_back(std::stod(v[i])); 

				if (std::any_of(std::begin(me_settings.internalGeometryParameters), std::end(me_settings.internalGeometryParameters), [](auto i) { return i <= 0; }))
					throw std::runtime_error(ASLI_ERRMSG::INTERNAL_GEOMETRY_DIMENSION_ZERO);

			} else {
				throw std::runtime_error(ASLI_ERRMSG::NO_INPUT_GEOMETRY_SPECIFIED);
			}
		}

		// Unit cell type data file
		if (lt_type.type == "hybrid") {
			auto tap = getOptional<std::string>(model_files, "tap", true, accessed_keys);
			if (tap.has_value()) tapFile = *tap;
			else throw std::runtime_error(ASLI_ERRMSG::FAILED_TO_READ_TAP_KEY);
		}

		// Unit cell size data file
		if (lt_size.size == 0) {
			auto sap = getOptional<std::string>(model_files, "sap", true, accessed_keys);
			if (sap.has_value()) sapFile = *sap;
			else throw std::runtime_error(ASLI_ERRMSG::FAILED_TO_READ_SAP_KEY);
		}

		// Unit cell feature data file
		if (lt_feature.feature_val == 0) {
			auto fap = getOptional<std::string>(model_files, "fap", true, accessed_keys);
			if (fap.has_value()) fapFile = *fap;
			else throw std::runtime_error(ASLI_ERRMSG::FAILED_TO_READ_FAP_KEY);
		}

		// Output directory
		std::filesystem::path output_;
		auto output = getOptional<std::string>(model_files, "output", false, accessed_keys);
		if (output.has_value()) output_ = std::filesystem::u8path(*output);
		if (std::filesystem::is_directory(output_.parent_path())) {
			me_settings.output = std::filesystem::u8path(*output);
		} else {
			if (!std::filesystem::is_directory(me_settings.output.parent_path()) || !std::filesystem::exists(me_settings.output.parent_path()))
				std::filesystem::create_directory(me_settings.output.parent_path());

			if (output.has_value())
				std::cout << "  WARNING: Requested output location " << *output << " was not found. Defaulting to " 
					<< me_settings.output.parent_path() << "." << std::endl;
		}

		// Output filename
		if (me_settings.output.has_filename()) {
			if (!me_settings.output.extension().compare("vtu"))
				me_settings.output.replace_extension(".mesh");
		} else { // Append current time to input filename and use as name for output file
			time_t t = std::time(nullptr);
			std::ostringstream currentTime;
			currentTime << std::put_time(std::localtime(&t), "_%Y-%m-%d_%H%M");

			me_settings.output = me_settings.output.string() + me_settings.inputFile.stem().string() + currentTime.str() + ".mesh";
		}

		// Check for unused input/output file keys
		for (YAML::const_iterator it = model_files.begin(); it != model_files.end(); ++it) {
				std::string key = it->first.as<std::string>();
				if (accessed_keys.find(key) == accessed_keys.end()) {
					std::cout << "  WARNING: Unused \"" << fileKey << " -> " << key << "\" key" << std::endl;
				}
		}
	}

	// Other settings
	auto verbosity = getOptional<double>(config, "verbosity", false, accessed_keys);
	if (verbosity.has_value() && *verbosity >= -1) me_settings.verbosity = *verbosity;

	// Check for remaining unused keys
	for (YAML::const_iterator it = config.begin(); it != config.end(); ++it) {
			std::string key = it->first.as<std::string>();
			if (accessed_keys.find(key) == accessed_keys.end()) {
				std::cout << "  WARNING: Unused \"" << key << "\" key" << std::endl;
			}
	}

	// Load data from the user input files
	LoadInputFiles();

	// Compute the required lattice variables and create the necessary 
	// interpolation models
	SetUpLattice();
}

void ASLI::SetUpLattice() {
	/* Setups the lattice
	 * 
	 */

	/* Setup the unit cell scaling (sizing) variables */
	// Determine the bounding box of the outer shell 
	std::vector<double> bounds = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
	for (size_t i=0; i<shell.points.size(); ++i) {
		bounds[0] = std::min(bounds[0], shell.points[i].x());
		bounds[1] = std::max(bounds[1], shell.points[i].x());
		bounds[2] = std::min(bounds[2], shell.points[i].y());
		bounds[3] = std::max(bounds[3], shell.points[i].y());
		bounds[4] = std::min(bounds[4], shell.points[i].z());
		bounds[5] = std::max(bounds[5], shell.points[i].z());
	}

	std::vector<double> boundingBoxDimensions;
	for (size_t i = 0; i < bounds.size(); i+=2)
		boundingBoxDimensions.push_back(bounds[i+1] - bounds[i]);

	// Setup the required unit cell scaling variables
	lt_size.minUnitCellSize = HUGE_VAL;
	lt_size.maxUnitCellSize = -HUGE_VAL;

	if (lt_size.size == 0 ) { // If a .sap file is used to specify the sizing
		lt_size.meanUnitCellSize = 0;
		for (size_t i = 0; i < lt_size.data.size(); i++) {
			lt_size.meanUnitCellSize += lt_size.data[i][3];

			lt_size.minUnitCellSize = std::min(lt_size.minUnitCellSize, lt_size.data[i][3]);
			lt_size.maxUnitCellSize = std::max(lt_size.maxUnitCellSize, lt_size.data[i][3]);

			lt_size.data[i][3] =  2*PI/lt_size.data[i][3]; // Convert the sizing data to scaling data
		}
		lt_size.meanUnitCellSize = lt_size.meanUnitCellSize/lt_size.data.size();

		// Setup interpolator
		std::cout << "  Setting up size data interpolator (Linear kernel)... " << std::flush;
		TrilinearInterpolation::setup(lt_size.data, lt_size.interpModel_linear);
		std::cout << "Finished!" << std::endl;

	} else {
		if (lt_size.size > 0) { // If the size of the unit cell is specified
			lt_size.scaling = 2*PI/lt_size.size;

		} else if (lt_size.size < 0) { // If instead the number of unit cells is specified
			double boundingBoxMinimumDimension = *std::min_element(std::begin(boundingBoxDimensions), std::end(boundingBoxDimensions));
			lt_size.scaling = std::abs(lt_size.size) * 2*PI / boundingBoxMinimumDimension;
			lt_size.size = 2*PI/lt_size.scaling;

			std::cout << "  WARNING: Using \"# unit cells in the minimum dimension\" mode!" << std::endl;
		}

		lt_size.minUnitCellSize = lt_size.size;
		lt_size.maxUnitCellSize = lt_size.size;
		lt_size.meanUnitCellSize = lt_size.size;
	}

	/* Setup the required unit cell type variables when a .tap file is used to specify the typing */
	if (lt_type.type == "hybrid") {
		// Construct kd-tree for nearest neighbour serch
		SetUpKdtree(lt_typeTree.alglib_data.alglib_type_kdtree, lt_typeTree.alglib_data.IDX, lt_typeTree.kdt);

		// Construct the interpolation models (one for each unit cell type in the design)
		lt_type.typeVector.resize(lt_typeTree.alglib_data.unitCells.size());
		lt_type.interpModel_linear.resize(lt_typeTree.alglib_data.unitCells.size());

		std::cout << "  Setting up type data interpolator (Linear kernel) " << std::flush;
		for (size_t i = 0; i < lt_typeTree.alglib_data.unitCells.size(); i++) {
			std::cout << i << "... " << std::flush;
			lt_type.typeVector[i] = lt_typeTree.alglib_data.unitCells[i].type;
			SetUpTypeInterpolator(lt_typeTree.alglib_data.alglib_type_kdtree, 
			                      lt_typeTree.alglib_data.unitCells[i].weights, 
			                      lt_type.interpModel_linear[i]);
		}
		std::cout << "Finished!" << std::endl;
	}

	/* Setup the required unit cell feature variables when a .fap file is used to specify the feature information */
	if (lt_feature.feature_val == 0) {
		// Validity check
		int counter = 0;
		double t;
		for(const auto& vec : lt_feature.data) {
			Point p;
			p = Point(vec[0], vec[1], vec[2]);
				
			if (ASLI::checkFeatureFieldValidity(p, vec[3], lt_type, lt_size, lt_feature) == 1)
				counter += 1;
		}
		if (counter > 0) {
			std::string TEXT = "  Feature value outside range adviced for " + std::to_string(counter) + " unit cell " + ((counter > 1) ? " types" : " type") + " in .fap file";
			throw ExceptionError(TEXT, nullptr);
		}

		// Interpolate data
		std::cout << "  Setting up feature data interpolator (Linear kernel)... " << std::flush;
		TrilinearInterpolation::setup(lt_feature.data, lt_feature.interpModel_linear);
		std::cout << "Finished!" << std::endl;
	} else {
		double counter = ASLI::checkFeatureValueValidity(lt_type, lt_size, lt_feature);
		if (counter > 0) {
			std::string TEXT = "  Feature value outside range adviced for " + std::to_string(counter) + " unit cell " + ((counter > 1) ? " types" : " type");
			throw ExceptionError(TEXT, nullptr);
		}
	}


	/* Print inputs used for setup to the command line */

	// Lattice parameters
	std::cout << "\n  Unit cell type: " << lt_type.type 
		        << " (" << lt_type.side << ")" << std::flush;
	if (lt_type.type == "hybrid") {
		std::cout << ", Filter radius: " << lt_type.filterRadius
		          << ", Correction factor: " << lt_type.correctionFactor << std::flush;
	}
	std::cout << std::endl;

	if ( lt_size.minUnitCellSize != lt_size.maxUnitCellSize) {
		std::cout << "  Range of unit cell sizes: " << lt_size.minUnitCellSize
		          << " to " << lt_size.maxUnitCellSize << std::endl;
	} else
		std::cout << "  Unit cell size: " << lt_size.size << std::endl;

	std::cout << "  " << lt_feature.feature << ": " << lt_feature.feature_val << std::endl;
	if ( lt_feature.feature == "wallSize" || lt_feature.feature == "poreSize")
		std::cout << ", Feature mode: " << lt_feature.mode << std::flush;

	if ( lt_feature.feature == "userDefined") {
		std::cout << "  User defined feature: " << lt_feature.udf.userDefinedFeature
		          << ", User defined constants: [" << lt_feature.udf.A
		          << ", " << lt_feature.udf.B
		          << ", " << lt_feature.udf.C
		          << ", " << lt_feature.udf.D
		          << ", " << lt_feature.udf.E << "]"<< std::endl;
	}

	// Mesh settings
	std::cout << "  Mesher: " << me_settings.mesher << std::flush;
	#ifdef CGAL_CONCURRENT_MESH_3
		if ( me_settings.mesher == "CGAL" || me_settings.mesher == "CGAL_OLD" )
			std::cout << " (# threads = " << me_settings.n_threads << ")" << std::flush;
	#endif
	if ( me_settings.mesher == "MMG" )
		if ( me_settings.MMG_memory > 0 )
			std::cout << " (Max. memory = " << me_settings.MMG_memory << " MB)" << std::flush;

	std::cout << ", Mesh size: " << me_settings.elementSize << std::flush;

	if ( me_settings.edgeProtectionAngle > 0 )
		std::cout << ", Edge protection angle: " << me_settings.edgeProtectionAngle << " deg" << std::flush;
	else
		std::cout << ", Edge protection angle: disabled" << std::flush;
	
	std::cout << ", Volume mesh: " << std::boolalpha << me_settings.isVolumeMesh << std::noboolalpha << std::endl;

	if (me_settings.mesher == "MMG") {
		std::cout << "    MMG settings -> hausd: " << me_settings.MMG_hausd
		          << ", hgrad: " << me_settings.MMG_hgrad
		          << ", hinitial: " << me_settings.MMG_hinitial << std::flush;

	} else if ( me_settings.mesher == "CGAL" || me_settings.mesher == "CGAL_OLD" ) {
		std::cout << "    CGAL settings -> Facet distance: " << me_settings.CGAL_facetDistance
		          << ", Facet angle: " << me_settings.CGAL_facetAngle
		          << ", Cell radius edge ratio: " << me_settings.CGAL_cellRadiusEdgeRatio
							<< ", Relative error bound: " << me_settings.CGAL_relativeErrorBound
							<< std::flush;

		if ( me_settings.CGAL_facetSize != me_settings.CGAL_cellSize )
			std::cout << ", Facet size: " << me_settings.CGAL_facetSize << std::flush;
				
		if ( me_settings.CGAL_edgeSize != me_settings.CGAL_cellSize )
				std::cout << ", Edge size: " << me_settings.CGAL_edgeSize << std::flush;

		if ( me_settings.CGAL_minEdgeSize != me_settings.CGAL_cellSize )
			std::cout << ", Facet size: " << me_settings.CGAL_minEdgeSize << std::flush;

		if ( me_settings.edgeProtectionAngle > 0 )
			std::cout << ", Poisson offset: " << me_settings.CGAL_poissonOffset << std::flush;
	}
	std::cout << std::endl;
}

void ASLI::SetUpKdtree(const alglib::real_2d_array &coordinates, const alglib::integer_1d_array &tags, alglib::kdtree &kdt) {
	/* Constructs a kd-tree ...
	 * Inputs:
	 *   coordinates : Coordinates (each row corresponds to one point)
	 *   tags        : Coordinates ID's
	 * Output:
	 *   kdt         : Kd-tree
	 */

	std::cout << "  Building kd-tree... " << std::flush;

	// Norm type (If the norm is changed the metricRadius passed on to the
	// kd-tree query function MUST be changed accordinglly!)
	alglib::ae_int_t normtype = 2; // 0 = infinity-norm, 1 = 1-norm, 2 = 2-norm

	// Built k-d tree
	try {
		alglib::kdtreebuildtagged(coordinates, tags, 3, 1, normtype, kdt);
	} catch (const alglib::ap_error& e) {
		throw ExceptionError(e.msg.c_str(), nullptr);
	}

	std::cout << "Finished!" << std::endl;
}

void ASLI::SetUpTypeInterpolator(const alglib::real_2d_array &coordinates, const std::vector<double> &weights, modelData &interpolationModel) {
	/* Creates a linear interpolation model of the provided typing data.
	 * Inputs:
	 *   coordinates        : Coordinates (each row corresponds to one point)
	 *   weights            : Weights corresponding to each coordinate
	 * Output:
	 *   interpolationModel : Structure containing interpolation model data
	 */

	size_t nWeights = weights.size();
	Eigen::Matrix3Xd dataPoints(3, nWeights);
	Eigen::VectorXd filteredWeights(nWeights);
	for (int i = 0; i < nWeights; i++) {
		dataPoints.col(i) = Eigen::Vector3d(coordinates[i][0],
		                                    coordinates[i][1],
		                                    coordinates[i][2]);

		// Compute local unit cell size
		Point p(coordinates[i][0], coordinates[i][1], coordinates[i][2]);
		double size = Infill::sizing_function(p, lt_size, "");
		double localFilterSize = lt_type.filterRadius * size;

		// Compute the weights
		if (lt_type.filterRadius > 0) { // If radius was specified
			try {
				alglib::real_1d_array x;
				alglib::ae_int_t nNeightboursFound;

				std::vector<double> xTemp{dataPoints(0,i), dataPoints(1,i), dataPoints(2,i)};
				x.setcontent(xTemp.size(), &(xTemp[0]));

				// Nearest neighbour search (Using k-d tree)
				nNeightboursFound = alglib::kdtreequeryrnnu(lt_typeTree.kdt, x, localFilterSize);// radius = size unit cell * userparmeter/2 // Radius must be in the same norm used to construct the k-d tree.
				if (nNeightboursFound == 0) { // If no radius is specified or no neighbours are found within requested radious
					nNeightboursFound = alglib::kdtreequeryknn(lt_typeTree.kdt, x, 1);
					if (nNeightboursFound == 0) {
						throw ExceptionError("Error: No neighbours found", nullptr);
					}
				}
				alglib::integer_1d_array IDX;
				alglib::kdtreequeryresultstags(lt_typeTree.kdt, IDX);

				alglib::real_1d_array distances;
				alglib::kdtreequeryresultsdistances(lt_typeTree.kdt, distances);

				// Filter weights
				neighbourData neighbours;
				for (int jj = 0; jj < IDX.length(); jj++) {
					neighbours.IDX.push_back(IDX[jj]);
					neighbours.distance.push_back(distances[jj]);
				}
				filteredWeights(i) = Filter::filter(neighbours, weights, lt_type.filterType, localFilterSize);

			} catch (const alglib::ap_error& e) {
				throw ExceptionError(e.msg.c_str(), nullptr);
			}

		} else { // Pass on the weights without filtereing
			filteredWeights(i) = weights[i];
		}
	}

	// Compute interpolation coefficients
	bool out = TrilinearInterpolation::setup(dataPoints, filteredWeights, interpolationModel);
}

void ASLI::LoadInputFiles() {
	/* Loads the required data from the provided input files.
	 * 
	 */

	// Load input geometry or generate internally
	if ( !me_settings.inputFile.empty() ) { // Load input file
		std::cout << "  Opening " << me_settings.inputFile.string() << std::endl;
		std::ifstream stl_file(me_settings.inputFile.string(), std::ios::binary);
		if ( CGAL::IO::read_STL(stl_file, shell.points, shell.polygons) == true ) {
			stl_file.close();
		} else {
			stl_file.close();
			throw std::runtime_error(ASLI_ERRMSG::FAILED_TO_OPEN_STL_FILE);
		}
	} else { // Generate geometry internally
		t_point t0 = TicToc::tic("");
		std::cout << "  Generating " << me_settings.internalGeometry << " geometry... " << std::flush; 
		BasicGeometries::implicit2surface(me_settings.internalGeometry, 
			me_settings.internalGeometryParameters, shell.points, shell.polygons);
		TicToc::toc(t0, "completed in ");
	}

	// Check if input geometry satisfies  basic requirements
	{	Polyhedron polygon_surface;
		CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(shell.points, shell.polygons, polygon_surface);
		
		// Save surface for DEBUG PURPOSES!!!
		if (me_settings.verbosity > 5) {
			std::ofstream stl_fileOut("initial_geometry.stl", std::ios::out | std::ios::binary);
			CGAL::set_mode(stl_fileOut, CGAL::IO::BINARY);
			CGAL::IO::write_STL(stl_fileOut, polygon_surface);
			stl_fileOut.close();
		}

		if (!CGAL::is_closed(polygon_surface)) { throw std::runtime_error(ASLI_ERRMSG::INPUT_GEOMETRY_IS_NOT_CLOSED); }
	}

	// Read the unit cell typing file (.tap file)
	if (lt_type.type == "hybrid") {
		ASLI_IO::read_tap(tapFile, ',', lt_typeTree.alglib_data);

		// Validity check
		int counter = 0;
		for (const auto& uc : lt_typeTree.alglib_data.unitCells) {
			if (Infill::UNIT_CELLS.find(uc.type) == Infill::UNIT_CELLS.end() || uc.type == "hybrid") { counter += 1; }
		}
		if (counter > 0) { throw std::runtime_error(std::to_string(counter) + ASLI_ERRMSG::INVALID_TYPE_FIELD); }
	}

	// Read the unit cell sizing file (.sap file)
	if (lt_size.size == 0) {
		ASLI_IO::read_csv(sapFile, ',', lt_size.data);

		// Validity check
		int counter = 0;
		for(const auto& vec : lt_size.data) {	if (vec[3] <= 0) { counter +=1; }	}
		if (counter > 0) { throw std::runtime_error(std::to_string(counter) + ASLI_ERRMSG::INVALID_SIZE_FIELD); }
	}

	// Read the unit cell feature data file (.fap file)
	if (lt_feature.feature_val == 0) {
		ASLI_IO::read_csv(fapFile, ',', lt_feature.data);

		// Initial validity check
		int counter = 0;
		for(const auto& vec : lt_feature.data) {	if (vec[3] < 0) { counter +=1; }	}
		if (counter > 0) { throw std::runtime_error(std::to_string(counter) + ASLI_ERRMSG::INVALID_FEATURE_FIELD); }
	}
}



std::vector<std::string> ASLI::split(const std::string &s, const std::string &delimiter) {
	// for string delimiter
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	std::string token;
	std::vector<std::string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
		token = s.substr (pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back (token);
	}

	res.push_back (s.substr (pos_start));
	return res;
}

double ASLI::checkFeatureValueValidity(const latticeType &lt_type, const latticeSize &lt_size, const latticeFeature &lt_feature) {
	/* Checks if feature value is within recommended bounds
	 * Inputs:
	 *   lt_type    : Lattice type data
	 *   lt_size    : Lattice size data
	 *   lt_feature : Lattice feature data
	 * Return:
	 *   Number of unit cell types for which feature value is ouside recomented bounds
	 */

	double scaling = Infill::sizing_function(Point(0,0,0), lt_size, "scaling");

	if (lt_type.type != "hybrid") { // If unit cell is fixed
		double t = Infill::internal::input2level(lt_type.type, scaling, lt_feature.feature, lt_feature.feature_val, lt_feature.udf, lt_feature.mode);
		if (lt_type.side == "scaffold") {
			return (t < std::get<1>(Infill::UNIT_CELLS.at(lt_type.type))) ? 1 : 0;
		}	else if (lt_type.side == "void") {
			return (t > std::get<2>(Infill::UNIT_CELLS.at(lt_type.type))) ? 1 : 0;
		}

	} else { // If unit cell is variable
		int counter = 0;
		for (size_t i = 0; i < lt_type.interpModel_linear.size(); i++) {
			double t = Infill::internal::input2level(lt_type.typeVector[i], scaling, lt_feature.feature, lt_feature.feature_val, lt_feature.udf, lt_feature.mode);
			if (lt_type.side == "scaffold") {
				return (t < std::get<1>(Infill::UNIT_CELLS.at(lt_type.typeVector[i]))) ? counter += 1 : 0;
			}	else if (lt_type.side == "void") {
				return (t > std::get<2>(Infill::UNIT_CELLS.at(lt_type.typeVector[i]))) ? counter += 1 : 0;
			}
		}
		if (counter > 0)
			return counter;
	}
	return 0;
}

double ASLI::checkFeatureFieldValidity(const Point &p, const double &featureValue, const latticeType &lt_type,
	const latticeSize &lt_size, const latticeFeature &lt_feature) {
	/* Checks if feature values are within recommended bounds
	 * Inputs:
	 *   p            : Coordinates of point to be evaluated
	 *   featureValue :
	 *   lt_type      : Lattice type data
	 *   lt_size      : Lattice size data
	 *   lt_feature   : Lattice feature data
	 * Return:
	 *   ???
	 */

	double scaling = Infill::sizing_function(p, lt_size, "scaling");

	if (lt_type.type != "hybrid") { // If unit cell is fixed
		double t = Infill::internal::input2level(lt_type.type, scaling, lt_feature.feature, featureValue, 
			lt_feature.udf, lt_feature.mode);
		if (lt_type.side == "scaffold") {
			return (t < std::get<1>(Infill::UNIT_CELLS.at(lt_type.type))) ? 1 : 0;
		}	else if (lt_type.side == "void") {
			return (t > std::get<2>(Infill::UNIT_CELLS.at(lt_type.type))) ? 1 : 0;
		}

	} else { // If unit cell is variable
		std::vector<double> weights(lt_type.interpModel_linear.size());
		for (size_t i = 0; i < lt_type.interpModel_linear.size(); i++) {
			Eigen::Vector3d pEigen = Eigen::Vector3d(p.x(), p.y(), p.z());
			weights[i] = TrilinearInterpolation::evaluate(pEigen, lt_type.interpModel_linear[i]);

			if (weights[i] > 0.5) { 
				double t = Infill::internal::input2level(lt_type.typeVector[i], scaling, lt_feature.feature, featureValue, lt_feature.udf, lt_feature.mode);
				if (lt_type.side == "scaffold") {
					return (t < std::get<1>(Infill::UNIT_CELLS.at(lt_type.typeVector[i]))) ? 1 : 0;
				}	else if (lt_type.side == "void") {
					return (t > std::get<2>(Infill::UNIT_CELLS.at(lt_type.typeVector[i]))) ? 1 : 0;
				}
			}
		}

	}
	return 0;
}