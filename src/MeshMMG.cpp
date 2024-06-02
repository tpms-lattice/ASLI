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

#include "Mesh.h"

/* MESHMMG discretizes the provided closed surface with the requested infill 
 * using the CGAL library to generate a volume mesh of the provided surface 
 * and the mmg library to discretize the geometry with the prescribed infill.
 * 
 * Author(s): F.P.B. (KU Leuven)
 */

void MeshMMG::implicit2volume(const polygonSoup &shell, const latticeType &lt_type,
                              const latticeSize &lt_size, const latticeFeature &lt_feature,
                              const meshSettings &me_settings) {
	/* Discretizes a closed volume with an implicitly defined infill.
	 * Inputs :
	 *  shell             : Structure containing the general geometry
	 *  lt_type           : Lattice type parameters
	 *  lt_size           : Lattice size parameters
	 *  lt_feature        : Lattice feature parameters
	 *  me_settings       : Mesh settings
	 * Outputs :
	 *  .stl file containing the surface triangulation
	 *  .mesh file containing the volume mesh (optional)
	 */

	MMG5_pMesh mmgMesh;
	MMG5_pSol  mmgLs, mmgMet; // mmgLs = level-set solution, mmgMet = size metric
	int np, ne, nt, na; // np = #vertices, ne = #tetra, nt = #triangles, na = #edges
	int ierr;

	// Initialize mesh and sol structures
	mmgMesh = NULL; mmgLs  = NULL; mmgMet = NULL;

	MMG3D_Init_mesh(MMG5_ARG_start,
	                MMG5_ARG_ppMesh, &mmgMesh, 
	                MMG5_ARG_ppLs,   &mmgLs,
	                MMG5_ARG_ppMet,  &mmgMet,
	                MMG5_ARG_end);

	const bool isUniform = (lt_type.type != "hybrid" && lt_size.size != 0 && lt_feature.feature_val != 0 
		&& (lt_type.side == "scaffold" || lt_type.side == "void" ));
	std::filesystem::path outputPath = me_settings.output;

	/* ------------------------------- Step 1 ------------------------------- */
	// Create a volume mesh of the closed manifold input surface with elements
	// of a size appropiate to capture the level-set describing the infill.
	TicToc::tic(); std::cout << "  Generating initial volume mesh... " << std::flush;

	F_C3t3 c3t3;
	Mesh::polehedral2volume(shell, lt_type, lt_size, lt_feature, me_settings, c3t3);
	Mesh::c3t3_to_MMG3D(c3t3, mmgMesh);
	c3t3.clear();

	// Save initial volume mesh for DEBUG PURPOSES!!!
	if (me_settings.verbosity > 5) {
		std::string tempFile = "MMG_Step_1.mesh";
		if ( MMG3D_saveMesh(mmgMesh, tempFile.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_MESH, nullptr);
	}

	TicToc::toc("completed in ");

	/* ------------------------------- Step 2 ------------------------------- */
	// Compute the level-set describing the lattice infill

	// Get mesh info and set solution size
	if ( MMG3D_Get_meshSize(mmgMesh, &np, &ne, NULL, &nt, NULL, &na) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_MESH_SIZE, nullptr);
	if ( MMG3D_Set_solSize(mmgMesh, mmgLs, MMG5_Vertex, np, MMG5_Scalar) != 1 ) 
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL_SIZE, nullptr);
	if (!isUniform) {
		if ( MMG3D_Set_solSize(mmgMesh, mmgMet, MMG5_Vertex, np, MMG5_Scalar) != 1 ) 
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL_SIZE, nullptr);
	}
	
	// Compute the level-set values within the volume to be provided with an 
	// infill and set them as the solution values
	double minFeatureSize = HUGE_VAL, maxFeatureSize = -HUGE_VAL;
	double wallSize, poreSize;

	TicToc::tic();
	if (isUniform)
		std::cout << "  Computing the level-set solution values... " << std::flush;
	else
		std::cout << "  Computing the level-set solution values and size metrics... " << std::flush;

	for (int i=1 ; i<=np ; i++) {
		// Get the coordinates of the current point
		double point[3];
		if ( MMG3D_Get_vertex(mmgMesh, &(point[0]), &(point[1]), &(point[2]), NULL, NULL, NULL) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_GET_VERTEX, nullptr);

		Point p(point[0], point[1], point[2]);

		// Compute the level-set solution value at the current coordinate
		double isovalue = Infill::TPMS_function(p, lt_type, lt_size, lt_feature);
		if ( MMG3D_Set_scalarSol(mmgLs, isovalue, i) != 1 ) 
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL, nullptr);

		// Compute the size metric at the current coordinate
		if (!isUniform) {
			// Compute the size of the walls and pores at the current coordinate
			featureSize localFeatureSize = Infill::featureSize_function(p, lt_type, lt_size, lt_feature);

			double solValueMet = (isovalue > 0) ? localFeatureSize.poreSize : localFeatureSize.wallSize;
			solValueMet *= me_settings.elementSize;

			// Make non-positive element sizes (or bellow some set threshold) equal to the size
			// of the unit cell at that location since there is no infill to capture anyway.
			double unitCellSize = Infill::sizing_function(p, lt_size, "");
			if (solValueMet <= unitCellSize*me_settings.threshold) {solValueMet = unitCellSize;}

			// Set size metric
			if (solValueMet > 0) {
				if ( MMG3D_Set_scalarSol(mmgMet, solValueMet, i) != 1 )
					throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_SOL_SIZE, nullptr);
			} else {
				throw ExceptionError(MMG_ERRMSG::INVALID_SIZE_METRIC, nullptr);
			}

			// Track the minimum and maximum feature size in the geometry
			minFeatureSize = std::min(minFeatureSize, std::min(localFeatureSize.wallSize, localFeatureSize.poreSize));
			maxFeatureSize = std::max(maxFeatureSize, std::max(localFeatureSize.wallSize, localFeatureSize.poreSize));
		}
	}

	// Check if the number of given entities match with mesh size
	if ( MMGS_Chk_meshData(mmgMesh, mmgLs) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_MESH_DATA_CHECK, nullptr);

	// Global size
	double hsiz;
	if (isUniform) {
		Point p(0, 0, 0);
		featureSize globalFeatureSize = Infill::featureSize_function(p, lt_type, lt_size, lt_feature);
		if (lt_type.side == "scaffold") {
			hsiz = me_settings.elementSize*globalFeatureSize.wallSize;
		} else {
			hsiz = me_settings.elementSize*globalFeatureSize.poreSize;
		}

		if ( hsiz <= 0 )
			throw ExceptionError(MMG_ERRMSG::INVALID_SIZE_METRIC, nullptr);

	} else { // Check validity of size metric
		if ( MMGS_Chk_meshData(mmgMesh, mmgMet) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_MESH_DATA_CHECK, nullptr);
	}

	TicToc::toc("completed in ");

	// Save mesh and metrics for DEBUG PURPOSES!!!
	if (me_settings.verbosity > 5) {
		std::string tempFile = "MMG_Step_2.mesh";
		if ( MMG3D_saveMesh(mmgMesh, tempFile.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_MESH, nullptr);

		std::string tempFileLS = "MMG_Step_2_LS.sol";
		if ( MMG3D_saveSol(mmgMesh, mmgLs, tempFileLS.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_SOL, nullptr);

		if (!isUniform) {
			std::string tempFileSZ = "MMG_Step_2_SZ.sol";
			if ( MMG3D_saveSol(mmgMesh, mmgMet, tempFileSZ.c_str()) != 1 )
				throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_SOL, nullptr);
		}
	}

	// Save MESH and level-set SOL file for external use (and skip to cleanup)
	if (me_settings.MMG_exportLS == true) {
		TicToc::tic(); std::cout << "  Exporting base mesh and level-set sol file... " << std::flush;

		outputPath.replace_extension(".mesh");
		std::string tempFile = outputPath.string(); //Extra step to ensure correct compilation in Windows

		if ( MMG3D_saveMesh(mmgMesh, tempFile.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_MESH, nullptr);

		if ( MMG3D_saveSol(mmgMesh, mmgLs, tempFile.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_SOL, nullptr);

		TicToc::toc("completed in ");

		goto cleanup;
	}

	/* ------------------------------- Step 3 ------------------------------- */
	// Discretize the level-set describing the lattice infill

	// Set global meshing parameters
	if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_iso, 1) != 1 ) // Set to level-set mode
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_IPARAMETER + "MMG3D_IPARAM_iso", nullptr);

	if (lt_type.side == "scaffold") {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_rmc, 1e-1) != 1 ) // Remove small solid parasitic components
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_rmc", nullptr);
	} else if (lt_type.side == "void") {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_rmcvoid, 1e-1) != 1 ) // Remove small void parasitic components
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_rmcvoid", nullptr);
	}

	if (isUniform) {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hsiz, hsiz) != 1 ) // Constant mesh size
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_hsiz", nullptr);
	} else {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hmin, (me_settings.elementSize-0.1*me_settings.elementSize)*minFeatureSize) != 1 ) // Set min. mesh size
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_hmin", nullptr);
		//if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hmax, (me_settings.elementSize+0.1*me_settings.elementSize)*maxFeatureSize) != 1 ) // Set max. mesh size
		//	throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_hmax", nullptr);
	} 

	if (me_settings.MMG_hausd > 0) {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hausd, me_settings.MMG_hausd*lt_size.meanUnitCellSize) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_hausd", nullptr);
	}

	if (me_settings.isVolumeMesh) {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hgrad, me_settings.MMG_hgrad) != 1 ) // 
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_hgrad", nullptr);
	} else {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_hgradreq, -1) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_hgradreq", nullptr);
	}

	if (me_settings.edgeProtectionAngle > 0) {
		if ( MMG3D_Set_dparameter(mmgMesh, mmgLs, MMG3D_DPARAM_angleDetection, me_settings.edgeProtectionAngle) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_DPARAM_angleDetection", nullptr);
	} else {
		if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_angle, 0) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_DPARAMETER + "MMG3D_IPARAM_angle", nullptr);
	}

	if (lt_type.side == "scaffold") { // Select sub-domain to keep
		if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_numsubdomain, 3) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_IPARAMETER + "MMG3D_IPARAM_numsubdomain", nullptr);
	} else if ( lt_type.side == "void" ) {
		if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_numsubdomain, 2) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_IPARAMETER + "MMG3D_IPARAM_numsubdomain", nullptr);
	}

	if (me_settings.MMG_memory > 0) {
		if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_mem, me_settings.MMG_memory) != 1 ) // Max. memory size in MB
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_IPARAMETER + "MMG3D_IPARAM_mem", nullptr);
	}

	if ( MMG3D_Set_iparameter(mmgMesh, mmgLs, MMG3D_IPARAM_verbose, me_settings.verbosity) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_IPARAMETER + "MMG3D_IPARAM_verbose", nullptr);

	// Dicretize the geometry with its lattice infill
	TicToc::tic();
	if (me_settings.isVolumeMesh) std::cout << "  Meshing volume... " << std::flush;
	else std::cout << "  Meshing surface... " << std::flush;
	if (isUniform)
		ierr = MMG3D_mmg3dls(mmgMesh, mmgLs, NULL);
	else
		ierr = MMG3D_mmg3dls(mmgMesh, mmgLs, mmgMet);
	TicToc::toc("completed in ");

	if (ierr == MMG5_STRONGFAILURE)
		throw ExceptionError(MMG_ERRMSG::BAD_ENDING_OF_MMG3DLS, nullptr);
	else if (ierr == MMG5_LOWFAILURE)
		std::cerr << "ERROR: BAD ENDING OF MMG3DLS!" << std::endl;

	// Mesh info
	MMG3D_Get_meshSize(mmgMesh, &np, &ne, NULL, &nt, NULL, &na);
	std::cout << "\n    Surface triangulation: " << std::endl;
	std::cout << "      Number of vertices: N/A" << std::endl;
	//std::cout << "      Number of vertices: ~" << double(2 + 0.5 * nt) << std::endl;
	std::cout << "      Number of facets: " << nt << std::endl;
	if (me_settings.isVolumeMesh) {
		std::cout << "\n    Volume mesh: " << std::endl;
		std::cout << "      Number of vertices: " << np << std::endl;
		std::cout << "      Number of facets: " << nt << std::endl;
		std::cout << "      Number of tetrahedra: " << ne << std::endl;
	}

	/* ------------------------------- Step 4 ------------------------------- */
	// Store output and cleanup
	std::cout << "\nSAVING OUTPUT(S)... " << std::endl;
	std::cout << "  Output file(s): " << std::flush;

	// Save volume mesh
	if (me_settings.isVolumeMesh) {
		outputPath.replace_extension(".mesh");
		std::string tempFile = outputPath.string();

		if ( MMG3D_saveMesh(mmgMesh, tempFile.c_str()) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_MESH, nullptr);

		std::cout << outputPath.filename() << " & "  << std::flush;
	}

	// Save surface triangulation
	{	outputPath.replace_extension(".stl");
		std::string tempFile = outputPath.string();
		if ( internal::MMG3D_saveSurfaceAsSTL (mmgMesh, me_settings.STLFormat, tempFile.c_str()) != EXIT_SUCCESS )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SAVE_STL, nullptr);

		std::cout << outputPath.filename() << std::endl;
	}
	cleanup:
	// Free the MMG3D5 structures
	MMG3D_Free_all(MMG5_ARG_start,
	               MMG5_ARG_ppMesh, &mmgMesh,
	               MMG5_ARG_ppLs,   &mmgLs,
	               MMG5_ARG_ppMet,  &mmgMet,
	               MMG5_ARG_end);
}

void MeshMMG::internal::shell_to_MMGS(const polygonSoup &shell, MMG5_pMesh mmgMesh) {
	/* Convert shell to MMGS_mesh
	 * Inputs:
	 *   shell   : Polygon soup
	 * Output:
	 *   mmgMesh : MMG surface mesh
	 */

	int np = shell.points.size();
	int nt = shell.polygons.size();

	// Set mesh size: mesh, # vertices, # triangles, # edges
	if ( MMGS_Set_meshSize(mmgMesh, np, nt, 0) != 1 )
		throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_MESH_SIZE, nullptr);

	// Vertices
	for (unsigned int i = 0; i < np; i++) {
		if ( MMGS_Set_vertex(mmgMesh, shell.points[i].x(), shell.points[i].y(), shell.points[i].z(), 0, i+1) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_VERTEX, nullptr);
	}

	// Triangles
	for (unsigned int i = 0; i < nt; i++) {
		if ( MMGS_Set_triangle(mmgMesh,  shell.polygons[i][0]+1,  shell.polygons[i][1]+1,  shell.polygons[i][2]+1, 3, i+1) != 1 )
			throw ExceptionError(MMG_ERRMSG::FAILED_TO_SET_TRIANGLE, nullptr);
	}
}

bool MeshMMG::internal::MMG3D_saveSurfaceAsSTL(const MMG5_pMesh &mmgMesh, const std::string &format, const char *filename) {
	/* Extracts and saves the surface of an MMG5_pMesh to an .stl file.
	 * Inputs :
	 *   mmgMesh  : MMG mesh whose surface is to be stored as an .stl
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

	if ( mmgMesh->info.imprim >= 0 ) 
		fprintf(stdout,"  %%%% %s OPENED\n",outputFile);

	// Start of file
	if ( format == "ASCII" ) {
		myOutputFile << std::setprecision(15);
		myOutputFile << "solid "<< outputFile << " (STL generated by ASLI)" << "\n";
	} else { // BINARY
		// 80 bytes ASCII header
		myOutputFile << "FileType: Binary" + std::string(64,' ');

		// Number of facets
		const uint32_t nFacets = static_cast<uint32_t>(mmgMesh->nt);
		myOutputFile.write(reinterpret_cast<const char*>(&nFacets), sizeof(nFacets));
	}

	// Body of file
	if ( mmgMesh->xp && mmgMesh->xpoint ) {
		for (size_t k=1; k<=mmgMesh->nt; k++) {
			MMG5_xPoint pxp;

			MMG5_pTria ptt = &mmgMesh->tria[k];
			double faceNormal[3] = {0,0,0};
			for (size_t i=0; i<3; i++) {
				MMG5_Point point = mmgMesh->point[ptt->v[i]];
				pxp = mmgMesh->xpoint[point.xp];
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
					MMG5_Point point = mmgMesh->point[ptt->v[i]];
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
					MMG5_Point point = mmgMesh->point[ptt->v[i]];
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