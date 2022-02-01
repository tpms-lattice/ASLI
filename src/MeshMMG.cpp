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

#include "MeshMMG.h"
/* MESHMMG discretizes the provided closed surface with the requested infill
 * using the tetgen and mmg libraries.
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

int MeshMMG::implicit2volume(tetgenio *points, latticeType lt_type, latticeSize lt_size,
                             latticeFeature lt_feature, meshSettings me_settings,
                             std::filesystem::path &outputFile_string) {
	/* Discretizes a closed volume with an implicitly defined infill.
	 * Inputs :
	 *  points            : Point cloud data
	 *  lt_type           : Lattice type data
	 *  lt_size           : Lattice size data
	 *  lt_feature        : Lattice feature data
	 *  me_settings       : Mesh settings
	 *  outputFile_string : Output file path
	 * Outputs :
	 *  .stl file containing the surface triangulation
	 *  .mesh file containing the volume mesh
	 */

	MMG5_pMesh mmgMesh;
	MMG5_pSol  mmgLs, mmgMet;
	int np, ne, nt, na; // np = #vertices, ne = #tetra, nt = #triangles, na = #edges
	int ierr;
	char *outputFile, *tempFile;

	double minFeatureSize = HUGE_VAL;
	double maxFeatureSize = -HUGE_VAL;

	std::string temporaryFolder = "temp";
	std::string temporaryFile = temporaryFolder + "/mmg_temp.mesh";

	// Check if the folder to store temporary data exists (and create if missing)
	if (!std::filesystem::is_directory(temporaryFolder) || !std::filesystem::exists(temporaryFolder)) {
		std::filesystem::create_directory(temporaryFolder);
	}

	// Initialize mesh and sol structures
	mmgMesh = NULL; mmgLs  = NULL; mmgMet = NULL;

	MMG3D_Init_mesh(MMG5_ARG_start,
	                MMG5_ARG_ppMesh, &mmgMesh, 
	                MMG5_ARG_ppLs,   &mmgLs,
	                MMG5_ARG_ppMet,  &mmgMet,
	                MMG5_ARG_end);

	/* ------------------------------- Step 1 ------------------------------- */
	// Create a volume mesh of the closed manifold surface provided as input
	tempFile = (char *) calloc(strlen(temporaryFile.c_str()) + 1, sizeof(char));
	strcpy(tempFile,temporaryFile.c_str());

	fprintf(stdout,"\n  Creating initial volume mesh...\n");
	internal::surface2volume(points, tempFile, 
	                         me_settings.TETGEN_hvol*lt_size.minUnitCellSize);

	/* ------------------------------- Step 2 ------------------------------- */
	// Remesh the volume such that the elements have a size appropiate to capture
	// the level-set describing the infill geometry.
	fprintf(stdout,"  Remeshing inital volume mesh...\n");

	// Read .mesh(b) file
	if ( MMG3D_loadMesh(mmgMesh, tempFile) != 1 )  exit(EXIT_FAILURE);

	// Set met size
	MMG3D_Get_meshSize(mmgMesh, &np, &ne, NULL, &nt, NULL, &na);
	if ( MMG3D_Set_solSize(mmgMesh, mmgMet, MMG5_Vertex,np, MMG5_Scalar) != 1 )
		exit(EXIT_FAILURE);

	// Set met values
	fprintf(stdout,"    Computing local feature size metric...\n");
	for (int i=1 ; i<=np ; i++) {
		// Get the current point
		double point[3];
		MMG3D_Get_vertex(mmgMesh, &(point[0]), &(point[1]), &(point[2]), NULL, NULL, NULL);
		Point p(point[0], point[1], point[2]);

		// Compute met value for current point
		double solValueMet;
		featureSize localFeatureSize;
		localFeatureSize = Infill::featureSize_function(p, &lt_type, &lt_size, &lt_feature);
		double minFeature = std::min(localFeatureSize.wallSize, localFeatureSize.poreSize);
		double maxFeature = std::max(localFeatureSize.wallSize, localFeatureSize.poreSize);

		solValueMet = me_settings.MMG_hinitial * minFeature;
		if ( MMG3D_Set_scalarSol(mmgMet, solValueMet, i) != 1 ) exit(EXIT_FAILURE);

		// Determine minimum and maximum feature size in model
		minFeatureSize = std::min(minFeatureSize, minFeature);
		maxFeatureSize = std::max(maxFeatureSize, maxFeature);
	}

	// Set level of verbosity (-1 ... 10)
	if ( MMG3D_Set_iparameter(mmgMesh, mmgMet, MMG3D_IPARAM_verbose, -1) != 1 )  exit(EXIT_FAILURE);

	// Remesh volume
	fprintf(stdout,"    Refining volume mesh based on feature metric...\n");
	ierr = MMG3D_mmg3dlib(mmgMesh, mmgMet);
	if ( ierr == MMG5_STRONGFAILURE ) {
		fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
		return(ierr);
	} else if ( ierr == MMG5_LOWFAILURE )
		fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

	if ( MMG3D_saveMesh(mmgMesh, tempFile) != 1 ) {
		fprintf(stdout, "UNABLE TO SAVE VOLUME MESH\n");
		return(MMG5_STRONGFAILURE);
	}





/* TESTING MSHMET IMPROVEMENT START!!! */
	bool testImprovement = 1;
	if (testImprovement == 1) {
		/* ------------------------------ Step 2.1 ------------------------------ */
		// Optimize mesh sizing distribution by computing a signed distance based metric.
		// Get mesh info

		MMG3D_Get_meshSize(mmgMesh, &np, &ne, NULL, &nt, NULL, &na);

		// Set sol size
		if ( MMG3D_Set_solSize(mmgMesh, mmgMet, MMG5_Vertex,np, MMG5_Scalar) != 1 )
			exit(EXIT_FAILURE);

		//
		std::vector<double> solValueLs (np);
		double minLS = HUGE_VAL;
		double maxLS = -HUGE_VAL;
		
		fprintf(stdout,"    Computing level-set based metric... \n");
		for (size_t i=0 ; i<np ; i++) {
			// Get the current point value
			double point[3];
			MMG3D_Get_vertex(mmgMesh, &(point[0]), &(point[1]), &(point[2]), NULL, NULL, NULL);
			Point p(point[0], point[1], point[2]);

			// Compute the sol value for current point
			solValueLs[i] = Infill::TPMS_function(p, &lt_type, &lt_size, &lt_feature);
			solValueLs[i] = std::abs(solValueLs[i]); //solValueLs[i] = std::pow(std::abs(solValueLs[i]), 1.0/3.0);

			minLS = std::min(minLS, solValueLs[i]);
			maxLS = std::max(maxLS, solValueLs[i]);
		}

		double a = std::max(1e-6, minLS);//me_settings.MMG_hinitial
		double b = 1e2*maxLS;
		for (size_t i=0 ; i<np ; i++) {
			double solVals = (solValueLs[i] - minLS) / (maxLS - minLS);

			solVals = std::pow(solVals, 1.0/3.0);

			//solVals = a + ( ((solValueLs[i] - minLS)*(b-a)) / (maxLS - minLS) );

			// Set sol value
			if ( MMG3D_Set_scalarSol(mmgMet, solVals, i+1) != 1 )
				exit(EXIT_FAILURE);
		}
		
		// 
		if ( MMG3D_saveSol(mmgMesh, mmgMet, tempFile) != 1 ) {
			fprintf(stdout, "UNABLE TO SAVE SOL\n");
			return(MMG5_STRONGFAILURE);
		}

		// this will 'fill' the string command with the right stuff, assuming myFile and convertedFile are strings themselves
		char *mshmetName;
		std::string temporaryFileB = temporaryFolder + "/mmg_temp.new.mesh";

		mshmetName = (char *) calloc(strlen(temporaryFileB.c_str()) + 1, sizeof(char));
		strcpy(mshmetName,temporaryFileB.c_str());

		// Determine mshmet input parameters
		std::ostringstream hmin, hmax, mshmet_outFile;
		hmin << 0.05*minFeatureSize;
		hmax << 1*maxFeatureSize;
		mshmet_outFile << tempFile;

		std::vector<std::string> arguments = {"./mshmet", 
		                                      "-hmin", hmin.str(), 
		                                      "-hmax", hmax.str(), 
		                                      "-hgrad", "4.5",
		                                      "-eps", "0.05",
		                                      "-v", "0",
		                                      mshmet_outFile.str()};

		// Conctruct argv like input to pass on to mshmet
		std::vector<char*> argv;
		for (const auto& arg : arguments)
			argv.push_back((char*)arg.data());
		argv.push_back(nullptr);

		fprintf(stdout, "\n    ");
		for (size_t i=0; i < argv.size() - 1; i++)
			fprintf(stdout, "%s ", arguments[i].c_str());
		fprintf(stdout, "\n");

		// Construct signed distance based metric using mshmet
		main_mshmet(argv.size() - 1, argv.data());

		// Load computed metric
		if ( MMG3D_loadSol(mmgMesh, mmgMet, mshmetName) != 1 ) exit(EXIT_FAILURE);

		if ( MMG3D_Set_iparameter(mmgMesh, mmgMet, MMG3D_IPARAM_verbose, -1) != 1 )  exit(EXIT_FAILURE); // Set level of verbosity (-1 ... 10)
		if ( MMG3D_Set_dparameter(mmgMesh, mmgMet, MMG3D_DPARAM_hgrad, 4.5) != 1 )  exit(EXIT_FAILURE);

		// Remesh volume
		fprintf(stdout,"    Remeshing volume based on level-set metric...\n");
		ierr = MMG3D_mmg3dlib(mmgMesh, mmgMet);
		if ( ierr == MMG5_STRONGFAILURE ) {
			fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
			return(ierr);
		} else if ( ierr == MMG5_LOWFAILURE )
			fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

		if ( MMG3D_saveMesh(mmgMesh, tempFile) != 1 ) {
			fprintf(stdout, "UNABLE TO SAVE VOLUME MESH\n");
			return(MMG5_STRONGFAILURE);
		}
	}
/* TESTING MSHMET IMPROVEMENT END!!! */





	/* ------------------------------- Step 3 ------------------------------- */
	// Discretize the level-set describing the lattice infill
	fprintf(stdout,"\n  Discretizing level-set...\n");

	// Get mesh info
	MMG3D_Get_meshSize(mmgMesh, &np, &ne, NULL, &nt, NULL, &na);

	// Set sol size
	if ( MMG3D_Set_solSize(mmgMesh, mmgLs, MMG5_Vertex,np, MMG5_Scalar) != 1 )
		exit(EXIT_FAILURE);

	// Compute the level-set values within the volume to be provided with an 
	// infill and set them as the sol values
	fprintf(stdout,"    Computing the level-set solution values...\n");
	for (int i=0 ; i<np ; i++) {
		// Get the current point value
		double point[3];
		MMG3D_Get_vertex(mmgMesh, &(point[0]), &(point[1]), &(point[2]), NULL, NULL, NULL);
		Point p(point[0], point[1], point[2]);

		// Compute the sol value for current point
		double solValueLs;
		solValueLs = Infill::TPMS_function(p, &lt_type, &lt_size, &lt_feature);

		// Set sol value
		if ( MMG3D_Set_scalarSol(mmgLs, solValueLs, i+1) != 1 )
			exit(EXIT_FAILURE);
	}

	// Recompute local feature size metric
//	if (me_settings.volumeMesh == true) {
//		fprintf(stdout,"    Recomputing local feature size metric... ");
//		for (int i=1 ; i<=np ; i++) {
//			// Get the current point
//			double point[3];
//			MMG3D_Get_vertex(mmgMesh, &(point[0]), &(point[1]), &(point[2]), NULL, NULL, NULL);
//			Point p(point[0], point[1], point[2]);
//
//			// Compute met value for current point
//			double solValueMet;
//			featureSize localFeatureSize;
//			localFeatureSize = Infill::featureSize_function(p, &lt_type, &lt_size, &lt_feature);
////			double minFeature = std::min(localFeatureSize.wallSize, localFeatureSize.poreSize);
////			double maxFeature = std::max(localFeatureSize.wallSize, localFeatureSize.poreSize);
//
//			double hlocal;
//			if (me_settings.MMG_hmax==0)
//				hlocal = 1;
//			else
//				hlocal = me_settings.MMG_hmax;
//
//			solValueMet = hlocal * localFeatureSize.wallSize;
//			if ( MMG3D_Set_scalarSol(mmgMet, solValueMet, i) != 1 ) exit(EXIT_FAILURE);
//		}
//	}

	// SAVE MESH AND SOL FILE TO BE USED FOR LEVEL_SET DISCRETIZATION (FOR
	// DEBUG PURPOSES) COMMENT WHEN NOT REQUIRED!!!
//	if ( MMG3D_saveMesh(mmgMesh, tempFile) != 1 ) {
//		fprintf(stdout, "UNABLE TO SAVE VOLUME MESH\n");
//		return(MMG5_STRONGFAILURE);
//	}
////	if ( MMG3D_saveSol(mmgMesh, mmgMet, tempFile) != 1 ) {
////		fprintf(stdout, "UNABLE TO SAVE SOL\n");
////		return(MMG5_STRONGFAILURE);
////	}
//	if ( MMG3D_saveSol(mmgMesh, mmgLs, tempFile) != 1 ) {
//		fprintf(stdout, "UNABLE TO SAVE SOL\n");
//		return(MMG5_STRONGFAILURE);
//	}

	// Set global meshing parameters
	//if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_mem, 5000) != 1 )  exit(EXIT_FAILURE); // Set max. memory size in MB

	if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_iso, 1) != 1 ) exit(EXIT_FAILURE);     // Use level-set
	if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_rmc, 1e-1) != 1 )  exit(EXIT_FAILURE); // Remove small solid parasitic components
	//if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_rmcvoid, 1e-5) != 1 )  exit(EXIT_FAILURE); // Remove void parasitic components

	if (me_settings.MMG_hmin > 0)
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hmin, me_settings.MMG_hmin*minFeatureSize) != 1 )  exit(EXIT_FAILURE); // Minimal mesh size
	if (me_settings.MMG_hmax > 0)
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hmax, me_settings.MMG_hmax*maxFeatureSize) != 1 )  exit(EXIT_FAILURE); // Maximal mesh size
	if (me_settings.MMG_hausd > 0)
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hausd, me_settings.MMG_hausd*(minFeatureSize+maxFeatureSize)/2) != 1 )  exit(EXIT_FAILURE); // Control global Hausdorff distance (on all the boundary surfaces of the mesh) CHECK!!
	
	if (me_settings.volumeMesh == true) {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hgrad, me_settings.MMG_hgrad) != 1 )  exit(EXIT_FAILURE); // 
	} else {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hgrad, -1) != 1 )  exit(EXIT_FAILURE);
	}

	if (me_settings.side == "void") { // Select a sub-domain to keep
		if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_numsubdomain, 2) != 1 )  exit(EXIT_FAILURE);
	} else if (me_settings.side == "scaffold") {
		if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_numsubdomain, 3) != 1 )  exit(EXIT_FAILURE);
	}

	if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_verbose, 2) != 1 )  exit(EXIT_FAILURE); // Tune level of verbosity -1 to 10 

	// Dicretize the level-set
	if (me_settings.volumeMesh == true) {
		fprintf(stdout,"    Discretizing volume...\n");
 		ierr = MMG3D_mmg3dls(mmgMesh, mmgLs, NULL);
//		ierr = MMG3D_mmg3dls(mmgMesh, mmgLs, mmgMet);
	} else {
		fprintf(stdout,"    Discretizing surface...\n");
		ierr = MMG3D_mmg3dls(mmgMesh, mmgLs, NULL);
	}
	if ( ierr == MMG5_STRONGFAILURE ) {
		fprintf(stdout,"BAD ENDING OF MMG3DLS: UNABLE TO SAVE MESH\n");
		return(ierr);
	} else if ( ierr == MMG5_LOWFAILURE )
		fprintf(stdout,"BAD ENDING OF MMG3DLS\n");

	/* ------------------------------- Step 4 ------------------------------- */
	// Store output and cleanup

	//Determine current time to append to filenames
	time_t t;
	char currentTime[50];
	std::time(&t);
	std::strftime(currentTime, sizeof(currentTime), "_%Y-%m-%d_%H%M", localtime(&t)); 

	outputFile_string = me_settings.output + currentTime + "_mmg.mesh";
	outputFile = (char *) calloc(outputFile_string.string().size() + 1, sizeof(char));
	strcpy(outputFile, outputFile_string.string().c_str());

	// Save volume mesh
	if (me_settings.volumeMesh == true) {
		if ( MMG3D_saveMesh(mmgMesh, outputFile) != 1 ) {
			fprintf(stdout, "UNABLE TO SAVE VOLUME MESH\n");
			return(MMG5_STRONGFAILURE);
		}
	}

	// Save surface triangulation
	if ( internal::MMG3D_saveSurfaceAsSTL (mmgMesh, me_settings.STLFormat, outputFile) != EXIT_SUCCESS ) {
		fprintf(stdout, "UNABLE TO SAVE SURFACE TRIANGULATION\n");
		return(MMG5_STRONGFAILURE);
	}

	// Delete temporary folder and its contents
	std::filesystem::remove_all(temporaryFolder);

	// Free the MMG3D5 structures
	MMG3D_Free_all(MMG5_ARG_start,
	               MMG5_ARG_ppMesh, &mmgMesh,
	               MMG5_ARG_ppLs,   &mmgLs,
	               MMG5_ARG_ppMet,  &mmgMet,
	               MMG5_ARG_end);

	// Cleanup
	free(outputFile);
	outputFile = NULL;

	free(tempFile);
	tempFile = NULL;

	return(ierr);
}


bool MeshMMG::internal::surface2volume(tetgenio *points, char *fileout, double hvol) {
	/* Converts a closed surface triangulation into a volume mesh.
	 * Inputs:
	 * 	points  : Point cloud data
	 *  fileout : Output filename
	 *  hvol    : Max. volume constraint
	 * Outputs:
	 * 	.mesh file containing the volume mesh
	 */
	
	FILE *fileID;
	std::string inputString;
	tetgenio out;

	// Tetgen input string
	inputString = "p";        // Tetrahedralize a piecewise linear complex (PLC).
//	inputString = "p/0.001";  // Tetrahedralize a piecewise linear complex (PLC). // 1.6?
	inputString += "q1.414";  // Refines mesh (to improve mesh quality).
	inputString += "a" + std::to_string( std::pow(hvol, 3) ); // max. tetrahedron volume constraint
	inputString += "g";       // Outputs mesh to .mesh
	inputString += "Q";       // Quiet: No terminal output except errors.

	/* ------------------------------- Step 1 ------------------------------- */
	// Create volume mesh
//	char tetgenInputCommand[inputString.size() + 1];
//	strcpy(tetgenInputCommand, inputString.c_str());
//	tetrahedralize(tetgenInputCommand, points, &out);

	char* tetgenInputCommand = new char[inputString.size() + 1]; 
	strcpy(tetgenInputCommand, inputString.c_str());
	tetrahedralize(tetgenInputCommand, points, &out);
	delete tetgenInputCommand;

	/* ------------------------------- Step 2 ------------------------------- */
	// Store output
	if( !(fileID = fopen(fileout,"w")) ) {
		fprintf(stderr,"  ** UNABLE TO OPEN OUTPUT MESH FILE.\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fileID,"MeshVersionFormatted 2\n");
	fprintf(fileID,"\nDimension 3\n");

	// Loop over points
	fprintf(fileID,"\nVertices\n%d\n", out.numberofpoints);
	for(size_t i = 0; i < out.numberofpoints; i++) {
		fprintf(fileID, "%.15lg %.15lg %.15lg %d \n", out.pointlist[i*3+0], out.pointlist[i*3+1], out.pointlist[i*3+2], 0);
	}

	// Loop over tetrahedra
	fprintf(fileID, "\nTetrahedra\n%d\n", out.numberoftetrahedra);
	int T[4];
	for(size_t i = 0; i < out.numberoftetrahedra; i++) {
		for(size_t j = 0; j<4; j++) {
			int index = out.tetrahedronlist[i * out.numberofcorners + j];
			T[j] = index;
		}
		fprintf(fileID, "%d %d %d %d %d \n", T[0], T[1], T[2], T[3], 1);
	}

	fprintf(fileID, "\nEnd\n");
	fclose(fileID);

	return EXIT_SUCCESS;
}

int MeshMMG::surface2volume(std::string inputFile, latticeType lt_type, latticeSize lt_size,
                            latticeFeature lt_feature, meshSettings me_settings,
                            std::filesystem::path &outputFile_string) {
	/* Discretizes a closed volume with an implicitly defined infill.
	 * Inputs :
	 *  inputFile         : Surface .stl file
	 *  lt_type           : Lattice type data
	 *  lt_size           : Lattice size data
	 *  lt_feature        : Lattice feature data
	 *  me_settings       : Mesh settings
	 *  outputFile_string : Output file
	 * Outputs :
	 *  .stl file containing the surface triangulation
	 *  .mesh file containing the volume mesh
	 */

	tetgenio points;
	MMG5_pMesh mmgMesh;
	MMG5_pSol  mmgLs, mmgMet;
	int np, ne, nt, na; // np = #vertices, ne = #tetra, nt = #triangles, na = #edges
	int ierr;
	char *outputFile, *tempFile;

	double minFeatureSize = HUGE_VAL;
	double maxFeatureSize = -HUGE_VAL;

	// Read input file
	char* stlFileChar = new char[inputFile.size() + 1];
	strcpy(stlFileChar, inputFile.c_str());
	std::cout << "  " << std::flush;
	if (points.load_stl(stlFileChar) == 0) { //1.6
	//if (points.read_stl(stlFileChar) == 0) { //1.4
		delete stlFileChar;
		exit(EXIT_FAILURE);
	}
	delete stlFileChar;

	// Check if the folder to store temporary data exists (and create if missing)
	std::string temporaryFolder = "temp";
	std::string temporaryFile = temporaryFolder + "/mmg_temp.mesh";
	if (!std::filesystem::is_directory(temporaryFolder) || !std::filesystem::exists(temporaryFolder)) {
		std::filesystem::create_directory(temporaryFolder);
	}

	// Initialize mesh and sol structures
	mmgMesh = NULL; mmgLs  = NULL; mmgMet = NULL;

	MMG3D_Init_mesh(MMG5_ARG_start,
	                MMG5_ARG_ppMesh, &mmgMesh,
	                MMG5_ARG_ppMet,  &mmgMet,
	                MMG5_ARG_end);

	/* ------------------------------- Step 1 ------------------------------- */
	// Create a volume mesh of the closed manifold surface provided as input
	tempFile = (char *) calloc(strlen(temporaryFile.c_str()) + 1, sizeof(char));
	strcpy(tempFile,temporaryFile.c_str());

	fprintf(stdout,"\n  Creating volume mesh...\n");
	internal::surface2volume(&points, tempFile, 
	                         me_settings.TETGEN_hvol*lt_size.minUnitCellSize);

	/* ------------------------------- Step 2 ------------------------------- */
	// Remesh the volume to ensure elements have disires quality.
	fprintf(stdout,"  Remeshing volume...\n");

	// Read .mesh(b) file
	if ( MMG3D_loadMesh(mmgMesh, tempFile) != 1 )  exit(EXIT_FAILURE);

	// Set met size
	MMG3D_Get_meshSize(mmgMesh, &np, &ne, NULL, &nt, NULL, &na);
	if ( MMG3D_Set_solSize(mmgMesh, mmgMet, MMG5_Vertex,np, MMG5_Scalar) != 1 )
		exit(EXIT_FAILURE);

	// Set met values
	fprintf(stdout,"    Computing the metValues... ");
	for (int i=1 ; i<=np ; i++) {
		// Get the current point
		double point[3];
		MMG3D_Get_vertex(mmgMesh, &(point[0]), &(point[1]), &(point[2]), NULL, NULL, NULL);
		Point p(point[0], point[1], point[2]);

		// Compute met value for current point
		double solValueMet;
		featureSize localFeatureSize;
		localFeatureSize = Infill::featureSize_function(p, &lt_type, &lt_size, &lt_feature);
		double minFeature = std::min(localFeatureSize.wallSize, localFeatureSize.poreSize);
		double maxFeature = std::max(localFeatureSize.wallSize, localFeatureSize.poreSize);

		solValueMet = me_settings.MMG_hinitial * minFeature;
		if ( MMG3D_Set_scalarSol(mmgMet, solValueMet, i) != 1 ) exit(EXIT_FAILURE);

		// Determine minimum and maximum feature size in model
		minFeatureSize = std::min(minFeatureSize, minFeature);
		maxFeatureSize = std::max(maxFeatureSize, maxFeature);
	}
	fprintf(stdout,"Finished!\n");

if (me_settings.MMG_hmin > 0)
		if ( MMG3D_Set_dparameter(mmgMesh, mmgMet, MMG3D_DPARAM_hmin, me_settings.MMG_hmin*minFeatureSize) != 1 )  exit(EXIT_FAILURE);		// Minimal mesh size
	if (me_settings.MMG_hmax > 0)
		if ( MMG3D_Set_dparameter(mmgMesh, mmgMet, MMG3D_DPARAM_hmax, me_settings.MMG_hmax*maxFeatureSize) != 1 )  exit(EXIT_FAILURE);		// Maximal mesh size
	if (me_settings.MMG_hausd > 0)
		if ( MMG3D_Set_dparameter(mmgMesh, mmgMet, MMG3D_DPARAM_hausd, me_settings.MMG_hausd*(minFeatureSize+maxFeatureSize)/2) != 1 )  exit(EXIT_FAILURE); 	// Control global Hausdorff distance (on all the boundary surfaces of the mesh) CHECK!!
	if ( MMG3D_Set_dparameter(mmgMesh, mmgMet, MMG3D_DPARAM_hgrad, me_settings.MMG_hgrad) != 1 )  exit(EXIT_FAILURE);	// 

	// Set level of verbosity (-1 ... 10)
	if ( MMG3D_Set_iparameter(mmgMesh, mmgMet, MMG3D_IPARAM_verbose, 2) != 1 )  exit(EXIT_FAILURE);

	// Remesh volume
	fprintf(stdout,"    Remeshing... \n");
	ierr = MMG3D_mmg3dlib(mmgMesh, mmgMet);
	if ( ierr == MMG5_STRONGFAILURE ) {
		fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
		return(ierr);
	} else if ( ierr == MMG5_LOWFAILURE )
		fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");
	fprintf(stdout,"Finished!\n");

	/* ------------------------------- Step 4 ------------------------------- */
	// Store output and cleanup

	//Determine current time to append to filenames
  time_t t;
  char currentTime[50];
	std::time(&t);
  std::strftime(currentTime, sizeof(currentTime), "_%Y-%m-%d_%H%M", localtime(&t)); 

	outputFile = (char *) calloc(outputFile_string.string().size() + 1, sizeof(char));
	strcpy(outputFile, outputFile_string.string().c_str());

	// Save volume mesh
	if (me_settings.volumeMesh == true) {
		if ( MMG3D_saveMesh(mmgMesh, outputFile) != 1 ) {
			fprintf(stdout, "UNABLE TO SAVE VOLUME MESH\n");
			return(MMG5_STRONGFAILURE);
		}
	}

	// Save surface triangulation
	if ( internal::MMG3D_saveSurfaceAsSTL (mmgMesh, me_settings.STLFormat, outputFile) != EXIT_SUCCESS ) {
		fprintf(stdout, "UNABLE TO SAVE SURFACE TRIANGULATION\n");
		return(MMG5_STRONGFAILURE);
	}

	// Free the MMG3D5 structures
	MMG3D_Free_all(MMG5_ARG_start,
	               MMG5_ARG_ppMesh, &mmgMesh,
	               MMG5_ARG_ppMet,  &mmgMet,
	               MMG5_ARG_end);

	// Cleanup
	free(outputFile);
	outputFile = NULL;

	free(tempFile);
	tempFile = NULL;

	return(ierr);

}

bool MeshMMG::internal::MMG3D_saveSurfaceAsSTL(MMG5_pMesh mesh, std::string format, char *filename) {
	/* Extracts and saves the surface of an MMG5_pMesh to an .stl file.
	 * Inputs :
	 * 	 mesh     : Mesh whose surface is to be stored as an .stl
	 *   format   : Output file format: ASCII or Binary (Default)
	 *   filename : Output path including filename
	 * Output:
	 *   .stl file containing the surface triangulation
	 */

	FILE *fileID;
	char *ptr, *outputFile;

	// Open output file for writing
	MMG5_SAFE_CALLOC(outputFile, strlen(filename)+6, char, return 0);
  strcpy(outputFile, filename);
  ptr = strstr(outputFile,".mesh");
	if ( ptr ) { // if filename contains the mesh extension
    ptr = strstr(outputFile,".mesh");
    if ( ptr ) *ptr = '\0';
    strcat(outputFile,".stl");
  }

	std::ofstream myOutputFile(outputFile, std::ios::out | std::ios::binary);
	if(!myOutputFile) {
		fprintf(stderr,"  ** UNABLE TO OPEN OUTPUT MESH FILE.\n");
		exit(EXIT_FAILURE);
	}

	if ( mesh->info.imprim >= 0 ) 
		fprintf(stdout,"  %%%% %s OPENED\n",outputFile);

	// Start of file
	if ( format == "ASCII" ) {
		myOutputFile << std::setprecision(15);
		myOutputFile << "solid "<< outputFile << " (STL generated by ASLI)" << "\n";
	} else { // BINARY
		// 80 bytes ASCII header
		myOutputFile << "FileType: Binary" + std::string(64,' ');

		// Number of facets
		const uint32_t nFacets = static_cast<uint32_t>(mesh->nt);
		myOutputFile.write(reinterpret_cast<const char*>(&nFacets), sizeof(nFacets));
	}

	// Body of file
	if ( mesh->xp && mesh->xpoint ) {
		for (size_t k=1; k<=mesh->nt; k++) {
			MMG5_xPoint pxp;

			MMG5_pTria ptt = &mesh->tria[k];
			double faceNormal[3] = {0,0,0};
			for (size_t i=0; i<3; i++) {
				MMG5_Point point = mesh->point[ptt->v[i]];
				pxp = mesh->xpoint[point.xp];
				faceNormal[0] += pxp.n1[0]/3;
				faceNormal[1] += pxp.n1[1]/3;
				faceNormal[2] += pxp.n1[2]/3;
			}
			double myNorm = std::sqrt(faceNormal[0]*faceNormal[0] + 
			                          faceNormal[1]*faceNormal[1] +
			                          faceNormal[2]*faceNormal[2]);
			faceNormal[0] = faceNormal[0]/myNorm;
			faceNormal[1] = faceNormal[1]/myNorm;
			faceNormal[2] = faceNormal[2]/myNorm;

			if ( format == "ASCII" ) {
				myOutputFile << "  facet normal  " << faceNormal[0] << " " 
					<< faceNormal[1] << " " << faceNormal[2] <<"\n";
				myOutputFile << "    outer loop\n";
				for (size_t i=0; i<3; i++) {
					MMG5_Point point = mesh->point[ptt->v[i]];
					myOutputFile << "      vertex   " << point.c[0] << " " << point.c[1] 
						<< " " << point.c[2] <<"\n";
				}
				myOutputFile << "    endloop\n";
				myOutputFile << "  endfacet\n";

			} else { // BINARY
				// Normal
				for (size_t i = 0; i < 3; i++) {
					float f = static_cast<float>(faceNormal[i]); //int32_t
					myOutputFile.write(reinterpret_cast<const char*>(&f), sizeof(f));
				}

				// Vertex
				for (size_t i = 0; i < 3; i++) {
					MMG5_Point point = mesh->point[ptt->v[i]];
					for (size_t j = 0; j < 3; j++) {
						float p = static_cast<float>(point.c[j]);//int32_t
						myOutputFile.write(reinterpret_cast<const char*>(&p), sizeof(p));
					}
				}

				// Atribute
				uint16_t attribute = 0;
				myOutputFile.write(reinterpret_cast<const char*>(&attribute), 
					sizeof(attribute));
			}
		}
	}

	// End of file
	if ( format == "ASCII" ) {
		myOutputFile << "endsolid " << outputFile << std::flush;
	} else { // BINARY
		myOutputFile << std::flush;
	}

	myOutputFile.close();

	// Cleanup
	MMG5_SAFE_FREE(outputFile);
	outputFile = NULL;

	return EXIT_SUCCESS;
}