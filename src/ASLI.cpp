/* ==========================================================================
 *  This file is part of ASLI (A Simple Lattice Infiller)
 *  Copyright (C) KU Leuven, 2019-2023
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

ASLI::ASLI(std::string configFile, bool exportls) {
	me_settings.exportls = exportls;

	SetUp(configFile);
}

ASLI::ASLI(std::string _stlFile, std::string _tapFile, std::string _sapFile,
           std::string _fapFile, latticeType _lt_type, latticeSize _lt_size,
           latticeFeature _lt_feature, meshSettings _me_settings) {

	// Lattice settings
	lt_type = _lt_type;
	lt_size = _lt_size;
	lt_feature = _lt_feature;

	// Mesh settings
	me_settings = _me_settings;

	// User input files
	stlFile = _stlFile;   // Closed outer surface file
	if (lt_type.type == "hybrid")
		tapFile = _tapFile; // Unit cell typing file
	if (lt_size.size == 0)
		sapFile = _sapFile; // Unit cell sizing file
	if (lt_feature.feature_val == 0)
		fapFile = _fapFile; // Unit cell feature data file

	// Load the data from the user input files
	LoadInputFiles();

	// Compute the required lattice variables and create the necessary 
	// interpolation models
	SetUpLattice();
}

ASLI::~ASLI() {}

void ASLI::SetUp(std::string configFile) {
	/* Reads the data from the configuration file and subsequently calls the 
	 * necessary setup functions.
	 */

	// Load configuration file
	YAML::Node config;
	try {
		config = YAML::LoadFile(configFile);
		std::cout << "  Opening " << configFile << std::endl;
	} catch (const std::exception& e) {
		std::cerr << e.what() << "\nERROR: Unable to open " << configFile << std::endl;
		exit(EXIT_FAILURE);
	}

	// Set user lattice input parameters
	try {
  	lt_type.type = config["lt_type"].as<std::string>();
		if (lt_type.type == "hybrid") {
			lt_type.filterType = "gaussian";//config["lt_type_filterType"].as<std::string>();
			lt_type.filterRadius = config["lt_type_filterRadius"].as<double>();
			lt_type.correctionFactor = config["lt_type_correctionFactor"].as<double>();
			if (lt_type.filterRadius < 0 || lt_type.correctionFactor < 0) {
				std::cerr << "\nERROR: Filter radius and correction factor cannot be negative" << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		lt_size.size = config["lt_size"].as<double>();

		lt_feature.feature = config["lt_feature"].as<std::string>();
		lt_feature.feature_val = config["lt_feature_val"].as<double>();
		lt_feature.mode = config["lt_feature_mode"].as<std::string>();
		if (lt_feature.mode != "absolute")
			lt_feature.mode = "relative";
		if (lt_feature.feature_val < 0) {
			std::cerr << "\nERROR: Feature value (" << lt_feature.feature << ") cannot be negative" << std::endl;
			exit(EXIT_FAILURE);
		}

  } catch (const std::exception& e) {
    std::cerr << e.what() << "\nERROR: Incorrectly defined or missing lattice parameter(s) in " << configFile << std::endl;
		exit(EXIT_FAILURE);
  }

	// Set user defined feature parameters
	try {
		lt_feature.udf.userDefinedFeature = config["udf_userDefinedFeature"].as<std::string>();
		lt_feature.udf.A = config["udf_A"].as<double>();
		lt_feature.udf.B = config["udf_B"].as<double>();
		lt_feature.udf.C = config["udf_C"].as<double>();
		lt_feature.udf.D = config["udf_D"].as<double>();
		lt_feature.udf.E = config["udf_E"].as<double>();

  } catch (const std::exception& e) {
    std::cerr << e.what() << "\nERROR: Incorrectly defined or missing material parameter(s) in " << configFile << std::endl;
		exit(EXIT_FAILURE);
  }

	// Determine user input files
	try {
		YAML::Node model_files = config["files"];
		stlFile = model_files["stl"].as<std::string>();   // Closed outer surface file
		if (lt_type.type == "hybrid")
			tapFile = model_files["tap"].as<std::string>(); // Unit cell typing file
		if (lt_size.size == 0)
			sapFile = model_files["sap"].as<std::string>(); // Unit cell sizing file
		if (lt_feature.feature_val == 0)
			fapFile = model_files["fap"].as<std::string>(); // Unit cell feature data file

		me_settings.output = model_files["output"].as<std::string>(); // Output location/filename

		// Check if requested output folder exists, and if not, fall back on default output folder 
		if (!std::filesystem::is_directory(me_settings.output) || !std::filesystem::exists(me_settings.output)) {
			if (!std::filesystem::is_directory("outputs") || !std::filesystem::exists("outputs"))
				std::filesystem::create_directory("outputs");
			std::cout << "  WARNING: Requested output location \"" << me_settings.output << "\" was not found. Defaulting to the \"outputs\" folder." << std::endl;
			me_settings.output = "outputs/";
		}
		std::filesystem::path tempPath(me_settings.output);
		std::filesystem::path tempFile(stlFile);
		std::filesystem::path outputPath(tempPath / tempFile.filename().replace_extension(""));
		me_settings.output = outputPath.u8string();

  } catch (const std::exception& e) {
    std::cerr << e.what() << "\nERROR: Incorrectly defined or missing filename(s) in " << configFile << std::endl;
		exit(EXIT_FAILURE);
  }

	// Set user mesh settings
	try {
		me_settings.STLFormat = "Binary";//config["me_STLFormat"].as<std::string>();

		me_settings.mesher = config["me_mesher"].as<std::string>();
		me_settings.side = config["me_side"].as<std::string>();

		me_settings.volumeMesh = config["me_volumeMesh"].as<bool>();
		me_settings.n_threads = config["me_nThreads"].as<int>();

		me_settings.CGAL_facetAngle = config["me_facetAngle"].as<double>();
		if (me_settings.CGAL_facetAngle <= 0) me_settings.CGAL_facetAngle = 30.0;
		me_settings.CGAL_facetSize = config["me_facetSize"].as<double>();
		if (me_settings.CGAL_facetSize <= 0) me_settings.CGAL_facetSize = 1.0;
		me_settings.CGAL_facetDistance = config["me_facetDistance"].as<double>();
		if (me_settings.CGAL_facetDistance <= 0) {
			std::cerr << "\nERROR: Facet distance must be positive" << std::endl;
			exit(EXIT_FAILURE);
		}
		me_settings.CGAL_cellRadiusEdgeRatio = config["me_cellRadiusEdgeRatio"].as<double>();
		if (me_settings.CGAL_cellRadiusEdgeRatio <= 2) me_settings.CGAL_cellRadiusEdgeRatio = 3.0;
		me_settings.CGAL_cellSize = config["me_cellSize"].as<double>();
		if (me_settings.CGAL_cellSize <= 0) me_settings.CGAL_cellSize = 1.0;
		me_settings.CGAL_preserveEdges = true;//config["me_preserveEdges"].as<bool>();
		me_settings.CGAL_poissonOffset = config["me_poissonOffset"].as<double>();
		if (me_settings.CGAL_poissonOffset < 0.1) me_settings.CGAL_poissonOffset = 0.5;
		me_settings.CGAL_edgesProtectionAngle = 60; //config["me_edgesProtectionAngle"].as<double>();

		me_settings.TETGEN_hvol = config["me_hvol"].as<double>();
		if (me_settings.TETGEN_hvol <= 0) me_settings.TETGEN_hvol = 1.5;
		me_settings.MMG_hinitial = config["me_hinitial"].as<double>();
		if (me_settings.MMG_hinitial <= 0) me_settings.MMG_hinitial = 0.35;
		me_settings.MMG_hmin = config["me_hmin"].as<double>();
		if (me_settings.MMG_hmin < 0) me_settings.MMG_hmin = 0.0;
		me_settings.MMG_hmax = 	config["me_hmax"].as<double>();
		if (me_settings.MMG_hmax < 0) me_settings.MMG_hmax = 0.0;
		me_settings.MMG_hausd = config["me_hausd"].as<double>();
		if (me_settings.MMG_hausd <= 0) {std::cerr << "\nERROR: hausd must be positive" << std::endl; exit(EXIT_FAILURE);}
		me_settings.MMG_hgrad = config["me_hgrad"].as<double>();
		if (me_settings.MMG_hgrad <= 0) me_settings.MMG_hgrad = 1.3;

  } catch (const std::exception& e) {
    std::cerr << e.what() << "\nERROR: Incorrectly defined or missing mesh setting(s) in " << configFile << std::endl;
		exit(EXIT_FAILURE);
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
	
	// Determine the bounding box of the outer shell 
	std::vector<double> bounds = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
	#ifdef CGAL_MESH // Place somewhere else
		if( me_settings.mesher == "CGAL" ) {
			for (size_t i=0; i<shell.points.size(); ++i) {
				bounds[0] = std::min(bounds[0], shell.points[i].x());
				bounds[1] = std::max(bounds[1], shell.points[i].x());
				bounds[2] = std::min(bounds[2], shell.points[i].y());
				bounds[3] = std::max(bounds[3], shell.points[i].y());
				bounds[4] = std::min(bounds[4], shell.points[i].z());
				bounds[5] = std::max(bounds[5], shell.points[i].z());
			}
		}
	#endif
	#ifdef MMG_MESH // Place somewhere else
		if( me_settings.mesher == "MMG" ) {
			for (size_t i=0; i<shell.tetgenPoints.numberofpoints; i+=3) {
				bounds[0] = std::min(bounds[0], shell.tetgenPoints.pointlist[i]);
				bounds[1] = std::max(bounds[1], shell.tetgenPoints.pointlist[i]);
				bounds[2] = std::min(bounds[2], shell.tetgenPoints.pointlist[i+1]);
				bounds[3] = std::max(bounds[3], shell.tetgenPoints.pointlist[i+1]);
				bounds[4] = std::min(bounds[4], shell.tetgenPoints.pointlist[i+2]);
				bounds[5] = std::max(bounds[5], shell.tetgenPoints.pointlist[i+2]);
			}
		}
	#endif


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
		TrilinearInterpolation::setup(lt_size.data, &lt_size.interpModel_linear);
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

	// Setup the required unit cell feature variables when a .fap file is used to specify the feature information
	if (lt_feature.feature_val == 0) {
		std::cout << "  Setting up feature data interpolator (Linear kernel)... " << std::flush;
		TrilinearInterpolation::setup(lt_feature.data, &lt_feature.interpModel_linear);
		std::cout << "Finished!" << std::endl;
	}

	// Setup the required unit cell type variables when a .tap file is used to specify the typing
	if (lt_type.type == "hybrid") {
		// Construct kd-tree for nearest neighbour serch
		SetUpKdtree(&lt_typeTree.alglib_data.alglib_type_kdtree, &lt_typeTree.alglib_data.IDX, &(lt_typeTree.kdt));

		// Construct the interpolation models (one for each unit cell type in the design)
		lt_type.typeVector.resize(lt_typeTree.alglib_data.unitCells.size());
		lt_type.interpModel_linear.resize(lt_typeTree.alglib_data.unitCells.size());

		std::cout << "  Setting up type data interpolator (Linear kernel) " << std::flush;
		for (size_t i = 0; i < lt_typeTree.alglib_data.unitCells.size(); i++) {
			std::cout << i << "... " << std::flush;
			lt_type.typeVector[i] = lt_typeTree.alglib_data.unitCells[i].type;
			SetUpTypeInterpolator(lt_typeTree.alglib_data.alglib_type_kdtree, 
			                      lt_typeTree.alglib_data.unitCells[i].weights, 
			                      &lt_type.interpModel_linear[i]);
		}
		std::cout << "Finished!" << std::endl;
	}

	// Print inputs used for setup to the command line
	std::cout << "  Unit cell type: " << lt_type.type << std::flush;
	if (lt_type.type == "hybrid") {
		std::cout << ", Filter radius: " << lt_type.filterRadius
		          << ", Correction factor: " << lt_type.correctionFactor << std::flush;
	}
	std::cout << std::endl;

	std::cout << "  Min. unit cell size in design: " << lt_size.minUnitCellSize
	          << ", Max. unit cell size in design: " << lt_size.maxUnitCellSize << std::endl;

	std::cout << "  Unit cell feature: " << lt_feature.feature << std::flush;
	if ( lt_feature.feature == "wallSize" || lt_feature.feature == "poreSize")
		std::cout << ", Feature mode: " << lt_feature.mode << std::flush;
	std::cout << ", Feature value: " << lt_feature.feature_val << std::endl;

	if ( lt_feature.feature == "userDefined") {
		std::cout << "  User defined feature: " << lt_feature.udf.userDefinedFeature
		          << ", User defined constant 1: " << lt_feature.udf.A
		          << ", User defined constant 2: " << lt_feature.udf.B
		          << ", User defined constant 3: " << lt_feature.udf.C
		          << ", User defined constant 4: " << lt_feature.udf.D
		          << ", User defined constant 5: " << lt_feature.udf.E << std::endl;
	}

	std::cout << "  Mesher: " << me_settings.mesher 
		      << ", Side: " << me_settings.side
	          << ", Volume mesh: " << std::boolalpha << me_settings.volumeMesh << std::flush;
	if (me_settings.mesher == "MMG") {
		std::cout << ", hvol: " << me_settings.TETGEN_hvol
		          << ", hinitial: " << me_settings.MMG_hinitial
		          << ",\n               hmin: " << me_settings.MMG_hmin
		          << ", hmax: " << me_settings.MMG_hmax
		          << ", hausd: " << me_settings.MMG_hausd
		          << ", hgrad: " << me_settings.MMG_hgrad << std::flush;
	} else if (me_settings.mesher == "CGAL") {
		std::cout << ", Facet angle: " << me_settings.CGAL_facetAngle
		          << ",\n                Facet size: " << me_settings.CGAL_facetSize
		          << ", Facet distance: " << me_settings.CGAL_facetDistance
		          << ", Cell radius edge ratio: " << me_settings.CGAL_cellRadiusEdgeRatio
		          << ",\n                Cell size: " << me_settings.CGAL_cellSize
		          << ", Preserve edges: " << std::boolalpha << me_settings.CGAL_preserveEdges << std::flush;
		if (me_settings.CGAL_preserveEdges == true)
			std::cout << ", Poisson offset: " << me_settings.CGAL_poissonOffset << std::flush;
	}
	std::cout << std::endl;
}

void ASLI::SetUpKdtree(alglib::real_2d_array *coordinates, alglib::integer_1d_array *tags, alglib::kdtree *kdt) {
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
		alglib::kdtreebuildtagged(*coordinates, *tags, 3, 1, normtype, *kdt);
	} catch (const alglib::ap_error& e) {
		std::cerr << e.msg.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}

	std::cout << "Finished!" << std::endl;

	//return EXIT_SUCCESS
}

void ASLI::SetUpTypeInterpolator(alglib::real_2d_array coordinates, std::vector<double> weights, modelData *interpolationModel) {
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
		double size = Infill::sizing_function(p, &lt_size, "");
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
						std::cerr << "Error: No neighbours found" << std::endl;
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
				std::cerr << e.msg.c_str() << std::endl;
				exit(EXIT_FAILURE);
			}

		} else { // Pass on the weights without filtereing
			filteredWeights(i) = weights[i];
		}
	}

	// Compute interpolation coefficients
	bool out = TrilinearInterpolation::setup(dataPoints, filteredWeights, interpolationModel);

	//return EXIT_SUCCESS
}

void ASLI::LoadInputFiles() {
	/* Loads the required data from the provided input files.
	 * 
	 */

	// Read the stl file
	#ifdef CGAL_MESH
		if( me_settings.mesher == "CGAL" ) {
			std::cout << "  Opening " << stlFile << std::endl;
			std::ifstream stl_file(stlFile, std::ios::binary);
			CGAL::IO::read_STL(stl_file, shell.points, shell.polygons);
			stl_file.close();
		}
	#endif
	#ifdef MMG_MESH
		if( me_settings.mesher == "MMG" ) {
			// The 'pointlist' and 'facetlist' together return the polyhedron. Note: After
			// calling load_stl() or read_stl() there will be many duplicated points in 
			// 'pointlist'. These are unified during the Delaunay tetrahedralization process.
			char* stlFileChar = new char[stlFile.size() + 1];
			strcpy(stlFileChar, stlFile.c_str());
			std::cout << "  " << std::flush;
			if (shell.tetgenPoints.load_stl(stlFileChar) == 0) { // For 1.6 (if 1.4 use read_stl instead)
				delete stlFileChar;
				exit(EXIT_FAILURE);
			}
			delete stlFileChar;
		}
	#endif

	// Read the unit cell typing file (.tap file)
	if (lt_type.type == "hybrid") {
		if (ASLI_IO::read_tap(tapFile, ',', lt_typeTree.alglib_data) != EXIT_SUCCESS)
			exit(EXIT_FAILURE);
	}

	// Read the unit cell sizing file (.sap file)
	if (lt_size.size == 0) {
		if (ASLI_IO::read_csv(sapFile, ',', lt_size.data) != EXIT_SUCCESS)
			exit(EXIT_FAILURE);
	}

	// Read the unit cell feature data file (.fap file)
	if (lt_feature.feature_val == 0) {
		if (ASLI_IO::read_csv(fapFile, ',', lt_feature.data) != EXIT_SUCCESS)
			exit(EXIT_FAILURE);
	}
}