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

#ifndef MESHMMG_H
#define MESHMMG_H

#include "Mesh.h"
#include "Infill.h"
#include "MeshCGAL.h"

/* Mmg headers */
#include "mmg/mmg3d/libmmg3d.h"

/* MshMet headers */
extern "C" {
	#include "mshmet_.h"
}

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

namespace MeshMMG { 
	int implicit2volume(outerShell &shell, latticeType lt_type, latticeSize lt_size,
                     latticeFeature lt_feature, meshSettings me_settings,
                     std::filesystem::path &outputFile_string);

	namespace internal {
		bool MMG3D_saveSurfaceAsSTL(MMG5_pMesh mesh, std::string format, char *filename) ;
	}
};
#endif