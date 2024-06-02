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

#ifndef MESH_MMG_H
#define MESH_MMG_H

#include "Infill.h"
#include "ExceptionClass.h" // Custom exception class

/* Mmg headers */
#include "mmg/mmg3d/libmmg3d.h"
#include "mmg/mmgs/libmmgs.h"

/* Standard library headers */
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <filesystem>

static inline // (Taken from MMG [mmgcommon.h])
size_t myfree(void *ptr) { 
	size_t s;
	char * ptr_c = (char*)ptr;

	if ( !ptr ) return 0;

	ptr_c = ptr_c-sizeof(size_t);
	s = *((size_t*)ptr_c);
	free(ptr_c);

	return s;
}
static inline // (Taken from MMG [mmgcommon.h])
void * mycalloc(size_t c, size_t s) {
	char *ptr;
	ptr = (char *)calloc(c*s+sizeof(size_t),1);
	if (ptr == NULL)
		return NULL;
	else {
		*((size_t*)ptr)=c*s;
		ptr+=sizeof(size_t);
		return (void*)ptr;
	}
}
// Safe deallocation (Taken from MMG [mmgcommon.h])
#define MMG5_SAFE_FREE(ptr) do                \
	{ myfree(ptr);                            \
		ptr = NULL;                           \
	} while(0)
// Safe allocation with calloc (Taken from MMG [mmgcommon.h])
#define MMG5_SAFE_CALLOC(ptr,size,type,law) do   \
	{ ptr = (type*)mycalloc(size,sizeof(type));  \
		if ( !ptr ) {                            \
		perror("  ## Memory problem: calloc");   \
		law;                                     \
		}                                        \
	} while(0)

// Error messages
namespace MMG_ERRMSG {
	const std::string FAILED_TO_SET_MESH_SIZE = "Failed to set mesh size";
	const std::string FAILED_TO_SET_VERTEX = "Failed to set vertex";
	const std::string FAILED_TO_SET_TRIANGLE = "Failed to set triangle";
	const std::string FAILED_TO_SET_TETRA = "Failed to set tretrahedral";
	const std::string FAILED_TO_SET_SOL_SIZE = "Failed to set sol size";
	const std::string FAILED_TO_SET_SOL = "Failed to set sol value";
	const std::string FAILED_MESH_DATA_CHECK = "Failed mesh data check";
	const std::string FAILED_TO_SET_IPARAMETER = "Failed to set ";
	const std::string FAILED_TO_SET_DPARAMETER = "Failed to set ";

	const std::string FAILED_TO_GET_MESH_SIZE = "Failed to get mesh size";
	const std::string FAILED_TO_GET_VERTEX = "Failed to get vertex";
	const std::string FAILED_TO_GET_EDGE = "Failed to get edge";
	const std::string FAILED_TO_GET_TRIANGLE = "Failed to get triangle";

	const std::string INVALID_SIZE_METRIC = "Invalid size metric";

	const std::string BAD_ENDING_OF_MMG3DLS = "BAD ENDING OF MMG3DLS: UNABLE TO SAVE MESH";
	const std::string BAD_ENDING_OF_MMGSLIB = "BAD ENDING OF MMGSLIB: UNABLE TO SAVE MESH";

	const std::string FAILED_TO_SAVE_MESH = "Unable to save volume mesh";
	const std::string FAILED_TO_SAVE_STL = "Unable to save surface triangulation";
	const std::string FAILED_TO_SAVE_SOL = "Unable to save sol data";
}

#endif