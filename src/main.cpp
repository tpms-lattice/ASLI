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

#include "version.h"
#include "ASLI.h"       // Class containing all the required data
#ifdef MMG_MESH
	#include "MeshMMG.h"  // Implicit surface meshing using MMG
#endif
#ifdef CGAL_MESH
	#include "MeshCGAL.h" // Implicit surface meshing using CGAL
#endif

#include "TicToc.h"     // A simple timer

/* Standard library headers */
#include <iomanip>
//#include <cstring>
#include <filesystem>

/* ASLI (A Simple Lattice Infiller)
 * 
 * ASLI is a small program that allows to add a heterogeneous TPMS lattice
 * infill to 3D objects.
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

int main(int argc, char *argv[]) {
	/* THE VERSION SCREEN */
	if(argc == 2 && (strcmp(argv[1], "-version")==0 || strcmp(argv[1], "--version")==0)) {
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

		return EXIT_SUCCESS;
	}

	/* THE HELP SCREEN */
	if(argc == 2 && (strcmp(argv[1], "-help")==0 || strcmp(argv[1], "--help")==0)) {
		std::cout << "Usage \n" << std::endl;
		std::cout 
			<< "  ./ASLI\n"
			<< "  ./ASLI <path-to-configuration-file> \n" 
		<< std::endl;
		std::cout << "Options" << std::endl;
		std::cout
			<< "  --help,-help           : Display usage information\n"
			<< "  --version,-version     : Display version number\n"
		<< std::endl;

		std::cout << std::string(55,' ') << "    __QQ" << std::endl;
		std::cout << std::string(55,' ') << "`~~(_)_\">\n" << std::endl;

		return EXIT_SUCCESS;
	}

	/* THE WELCOME SCREEN */
	std::cout << std::string(15,' ') 
		<< "_____/\\/\\________/\\/\\/\\/\\/\\__/\\/\\________/\\/\\/\\/\\_" 
		<< std::endl;
	std::cout << std::string(15,' ') 
		<< "___/\\/\\/\\/\\____/\\/\\__________/\\/\\__________/\\/\\___" 
		<< std::endl;
	std::cout << std::string(15,' ') 
		<< "_/\\/\\____/\\/\\____/\\/\\/\\/\\____/\\/\\__________/\\/\\___" 
		<< std::endl;
	std::cout << std::string(15,' ') 
		<< "_/\\/\\/\\/\\/\\/\\__________/\\/\\__/\\/\\__________/\\/\\___" 
		<< std::endl;
	std::cout << std::string(15,' ') 
		<< "_/\\/\\____/\\/\\__/\\/\\/\\/\\/\\____/\\/\\/\\/\\/\\__/\\/\\/\\/\\_" 
		<< std::endl;
	std::cout << std::string(15,' ') 
		<< "__________________________________________________" << std::endl; 

	#ifdef CGAL_CONCURRENT_MESH_3       // P = parallel
		std::cout << std::string(15,' ') 
			<< "A Simple Lattice Infiller v" << ASLI_VERSION_MAJOR << "." 
			<< ASLI_VERSION_MINOR << " (P)       KU Leuven" << std::endl;
	#else
		std::cout << std::string(15,' ')  // S = Serial
			<< "A Simple Lattice Infiller v" << ASLI_VERSION_MAJOR << "." 
			<< ASLI_VERSION_MINOR << " (S)       KU Leuven" << std::endl;
	#endif
	std::cout << "\n" << std::endl;

	std::cout << "Hello! You are running ASLI." << std::endl;

	/* STEP 0: PROCESS COMMAND LINE INPUT */
	std::filesystem::path configFile;
	if (argc >= 2) {
		configFile = argv[1];
		if (configFile.has_extension() == false)
			configFile.replace_extension(".yml");

	} else {
		configFile = "config.yml";
	}

	t_point t0 = TicToc::tic(""); // START TIMER

	/* STEP 1: SETUP THE PROBLEM (INFILL AND OUTER SURFACE DATA) */
	std::cout << "\nSETTING UP PROBLEM..." << std::endl;
	ASLI *data = new ASLI(configFile.u8string());	
	
	/* STEP 2: MESH THE OUTER SURFACE WITH THE REQUESTED INFILL (AND SAVE
	THE OUTPUT) */
	std::cout << "\nDISCRETIZING SCAFFOLDED STRUCTURE..." << std::endl;

	std::filesystem::path outputFile_string;
	if (data->me_settings.mesher == "CGAL") { 
		#ifdef CGAL_MESH
			MeshCGAL::implicit2volume(data->shell, data->lt_type, data->lt_size, 
			                          data->lt_feature, data->me_settings,
																outputFile_string);
		#else
			std::cerr << "ERROR: Compile with the CGAL_MESH tag to enable this option" 
				<< std::endl;
			return(EXIT_FAILURE);
		#endif

	} else if (data->me_settings.mesher == "MMG") {
		#ifdef MMG_MESH
			MeshMMG::implicit2volume(&data->shell.tetgenPoints, data->lt_type, 
			                         data->lt_size, data->lt_feature, 
			                         data->me_settings, outputFile_string);
		#else
			std::cerr << "ERROR: Compile with the MMG_MESH tag to enable this option" 
				<< std::endl;
			return(EXIT_FAILURE);
		#endif

//	} else if (data->me_settings.mesher == "TEST") {
//			Test::implicit2volume_H(data->shell, data->lt_type, data->lt_size,
//			                        data->lt_feature, data->me_settings, 
//			                        outputFile_string);

	} else {
		std::cout << data->me_settings.mesher << std::endl;
		std::cerr << "ERROR_INVALID_MESHER_REQUEST" << std::endl;
		return(EXIT_FAILURE);
	}

	/* STEP 3: CLEAN UP */
	delete data;

	TicToc::toc(t0, "\nTotal elapsed time: "); // STOP TIMER
	std::cout << "Output file: " <<  outputFile_string.replace_extension(".stl") << std::endl;

	/* DONE! */
	std::cout << "\nDONE! (/^v^)/" << std::endl;
	std::cout << "\nCiao ciao!" << std::endl;
	std::cout << "\a";

	return EXIT_SUCCESS;
}