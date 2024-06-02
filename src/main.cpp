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

#include "version.h"
#include "ASLI.h"           // Class containing all necessary data
#include "Mesh.h"           // Implicit meshing functions
#include "ExceptionClass.h" // Custom exception class
#include "TicToc.h"         // A simple timer

/* TBB headers */
#ifdef CGAL_CONCURRENT_MESH_3
	#define TBB_PREVIEW_GLOBAL_CONTROL 1
	#include <tbb/global_control.h>
#endif

/* Standard library headers */
#include <iomanip>

/* ASLI (A Simple Lattice Infiller)
 * 
 * ASLI is a small program that allows to add a heterogeneous triply periodic
 * lattice infill to 3D objects.
 * 
 * Author(s): F.P.B. (KU Leuven)
 *            M.B. (KU Leuven) 
 */

void show_version() {
	/* THE VERSION SCREEN */
	#ifdef CGAL_CONCURRENT_MESH_3
		std::cout << "ASLI v" << ASLI_VERSION_MAJOR << "." << ASLI_VERSION_MINOR 
		<< " (P)" << std::endl; // P = Parallel
	#else
		std::cout << "ASLI v" << ASLI_VERSION_MAJOR << "." << ASLI_VERSION_MINOR 
		<< " (S)" << std::endl; // S = Serial
	#endif
	std::cout << "Copyright 2019-" << ASLI_RELEASE_YEAR << " KU Leuven" << std::endl;
	std::cout << "ASLI is free software; see the source for copying conditions. There is NO\n" <<
	"warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n" << std::endl;
}

void show_help() {
	/* THE HELP SCREEN */
	std::cout << "Usage" << std::endl;
	std::cout 
		<< "  ./ASLI\n"
		<< "  ./ASLI <path-to-configuration-file/name-of-configuration-file.yml>\n"
	<< std::endl;
	std::cout << "Options" << std::endl;
	std::cout
		<< "  -h,--help        : Display usage information\n"
		<< "  -v,--version     : Display version number   \n"
		<< "  \n"
		<< "Example of a (minimal) configuration file\n"
		<< "  +-------------------------------+      \n"
		<< "  |# --- CONFIGURATION FILE  --- #|      \n"
		<< "  +-------------------------------+      \n"
		<< "  |files:                         |      \n"
		<< "  |  internalGeometry: cuboid 1   |      \n"
		<< "  |lt_type: gyroid                |      \n"
		<< "  |lt_size: 1                     |      \n"
		<< "  |lt_feature: volumeFraction     |      \n"
		<< "  |lt_feature_val: 0.5            |      \n"
		<< "  |me_mesher: CGAL #CGAL or MMG   |      \n"
		<< "  |me_elementSize: 0.3            |      \n"
		<< "  +-------------------------------+      \n"
		<< "  \n"
		<< "List of configuration file keys\n"
		<< "  files -> internalGeometry   me_mesher                \n"
		<< "  files -> stl                me_nThreads              \n"
		<< "  files -> tap                me_elementSize           \n"
		<< "  files -> sap                me_elementSizeSolid      \n"
		<< "  files -> fap                me_edgeProtectionAngle   \n"
		<< "  files -> output             me_volumeMesh            \n"
		<< "  lt_type                     CGAL_facetDistance       \n"
		<< "  lt_type_side                CGAL_relativeErrorBound  \n"
		<< "  lt_type_filterRadius        CGAL_facetAngle          \n"
		<< "  lt_type_correctionFactor    CGAL_cellRadiusEdgeRatio \n"
		<< "  lt_size                     CGAL_poissonOffset       \n"
		<< "  lt_feature                  MMG_hausd                \n"
		<< "  lt_feature_val              MMG_hgrad                \n"
		<< "  lt_feature_mode             MMG_hinitial             \n"
		<< "  udf_feature                 MMG_exportLS             \n"
		<< "  udf_A                       verbosity                \n"
		<< "  udf_B                                                \n"
		<< "  udf_C                                                \n"
		<< "  udf_D                                                \n"
		<< "  udf_E                                                \n"
		<< " \n"
		<< "For more details regarding the use of ASLI users are\n"
		<< "referred to the user's manual.\n"
	<< std::endl;
	std::cout << std::string(49,' ') << "    __QQ" << std::endl;
	std::cout << std::string(49,' ') << "`~~(_)_\">\n\a" << std::endl;
}

void show_welcome_screen() {
	/* THE WELCOME SCREEN */
	int indent = 15;
	std::cout << "\n"
		<< std::string(indent,' ') << "_____/\\/\\________/\\/\\/\\/\\/\\__/\\/\\________/\\/\\/\\/\\_\n"
		<< std::string(indent,' ') << "___/\\/\\/\\/\\____/\\/\\__________/\\/\\__________/\\/\\___\n"
		<< std::string(indent,' ') << "_/\\/\\____/\\/\\____/\\/\\/\\/\\____/\\/\\__________/\\/\\___\n"
		<< std::string(indent,' ') << "_/\\/\\/\\/\\/\\/\\__________/\\/\\__/\\/\\__________/\\/\\___\n"
		<< std::string(indent,' ') << "_/\\/\\____/\\/\\__/\\/\\/\\/\\/\\____/\\/\\/\\/\\/\\__/\\/\\/\\/\\_\n"
		<< std::string(indent,' ') << "__________________________________________________" << std::endl;

	#ifdef CGAL_CONCURRENT_MESH_3
		std::cout << std::string(indent,' ') 
			<< "A Simple Lattice Infiller v" << ASLI_VERSION_MAJOR << "." 
			<< ASLI_VERSION_MINOR << " (P)       KU Leuven" << std::endl; // P = parallel
	#else
		std::cout << std::string(indent,' ')
			<< "A Simple Lattice Infiller v" << ASLI_VERSION_MAJOR << "." 
			<< ASLI_VERSION_MINOR << " (S)       KU Leuven" << std::endl; // S = Serial
	#endif
	std::cout << "\n" << std::endl;
}

int main(int argc, char *argv[]) {

	// Version or help info if requested
	for (int i = 1; i < argc; i++) {
		if( strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0 ) {
			show_version();
			return EXIT_SUCCESS;
		} else if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 ) {
			show_help();
			return EXIT_SUCCESS;
		}
	}

	if ( argc > 2 ) {
		std::cerr << ASLI_ERRMSG::INVALID_NUMBER_OF_COMMAND_LINE_ARGUMENTS << "\n" << std::endl;
		show_help();
		return(EXIT_FAILURE);
	}

	show_welcome_screen();
	std::cout << "Hello! You are running ASLI." << std::endl;

	/* STEP 0: PROCESS COMMAND LINE INPUT */
	std::filesystem::path configFile = "config.yml";
	if ( argc == 2 ) {
		configFile = argv[1];
		if (configFile.has_extension() == false) {
			std::cerr << ASLI_ERRMSG::CONFIGURATION_FILE_WITHOUT_EXTENSION << "\a" << std::endl;
			return(EXIT_FAILURE);
		}
	}

	t_point t0 = TicToc::tic(0); // START TIMER

	/* STEP 1: SETUP THE PROBLEM (INFILL AND OUTER SURFACE DATA) */
	std::cout << "\nSETTING UP PROBLEM..." << std::endl;
	//ASLI *data = new ASLI(configFile.string());
	
	ASLI *data = nullptr;
	try { 
		data = new ASLI(configFile.string());	

	} catch (YAML::BadFile e) {
		std::cerr << ASLI_ERRMSG::FAILED_TO_OPEN_CONFIG_FILE << " (" << e.what() << ")" << "\a" << std::endl;
		return(EXIT_FAILURE);

	} catch (YAML::ParserException e) {
		std::cerr << ASLI_ERRMSG::FAILED_TO_PARSE_CONFIG_FILE << " (" << e.what() << ")" << "\a" << std::endl;
		return(EXIT_FAILURE);

	} catch (const ExceptionError& e) {
		std::cerr << e.what() << "\n" << e.where() << "\a" << std::endl;
		return(EXIT_FAILURE);

	} catch (const std::exception& e) {
		std::cerr << e.what() << "\a" << std::endl;
		return(EXIT_FAILURE);

	} catch (...) {
		std::cerr << "Unknown error: ... \a" << std::endl;
		return(EXIT_FAILURE);
	}

	// Concurrency settings
	#ifdef CGAL_CONCURRENT_MESH_3
		if (data->me_settings.n_threads < 1) {
			data->me_settings.n_threads = 1;
		} else if (data->me_settings.n_threads > tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism)) {
			data->me_settings.n_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
			std::cout << "  Warning: Requested number of threads exceeds the maximum allowed. " 
			          << "Number of threads set to the maximum allowed, i.e. " << data->me_settings.n_threads 
			          << "." << std::endl;
		}

		static tbb::global_control control(tbb::global_control::max_allowed_parallelism, data->me_settings.n_threads);
		data->me_settings.n_threads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
	#endif

	/* STEP 2: MESH THE OUTER SURFACE WITH THE REQUESTED INFILL (AND SAVE
	THE OUTPUT) */
	std::cout << "\nDISCRETIZING SCAFFOLDED STRUCTURE..." << std::endl;

	try {
		if (data->me_settings.mesher == "CGAL_OLD")
			MeshCGAL::implicit2volume_old(data->shell, data->lt_type, data->lt_size,
				data->lt_feature, data->me_settings);
		else if (data->me_settings.mesher == "CGAL")
			MeshCGAL::implicit2volume(data->shell, data->lt_type, data->lt_size,
				data->lt_feature, data->me_settings);
		else if (data->me_settings.mesher == "MMG")
			MeshMMG::implicit2volume(data->shell, data->lt_type, data->lt_size,
				data->lt_feature, data->me_settings);
//		else if (data->me_settings.mesher == "TEST")
//			MeshNG::implicit2volume(data->shell, data->lt_type, data->lt_size, 
//				data->lt_feature, data->me_settings);
		else
			throw ExceptionError(ASLI_ERRMSG::INVALID_MESHER, nullptr);

	} catch (const ExceptionError& e) {
		std::cerr << e.what() << "\n" << e.where() << "\a" << std::endl;
		return EXIT_FAILURE;
	} catch (const std::exception& e) {
		std::cerr << e.what() << "\a" << std::endl;
		return EXIT_FAILURE;
	} catch (...) {
		std::cerr << "Unknown error: ... \a" << std::endl;
		return(EXIT_FAILURE);
	}

	TicToc::toc(t0, "\nTotal elapsed time: "); // STOP TIMER
	
	/* STEP 3: COPY CONFIG FILE USED TO OUTPUT FOLDER */
	try {
		std::filesystem::copy_file(configFile, data->me_settings.output.replace_extension(".yml"),
			std::filesystem::copy_options::overwrite_existing);
	} catch (std::filesystem::filesystem_error& e) {
		std::cout << "WARNING: Could not copy config file to the outputs folder." << std::endl;
		return(EXIT_FAILURE);
	}

	/* STEP 4: CLEAN UP */
	delete data;
	data = nullptr; 

	/* DONE! */
	std::cout << "\nDONE! (/^v^)/" << std::endl;
	std::cout << "\nCiao ciao!" << "\a" << std::endl;

	return EXIT_SUCCESS;
}