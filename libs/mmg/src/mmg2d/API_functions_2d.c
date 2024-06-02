/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file mmg2d/API_functions_2d.c
 * \brief C API functions definitions for MMG2D library.
 * \author Algiane Froehly (Bx INP/Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmg2d/liblibmmg2d_private.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * C API for MMG2D library.
 *
 */

#include "libmmg2d.h"
#include "libmmg2d_private.h"

int MMG2D_Init_mesh(const int starter,...) {
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = MMG2D_Init_mesh_var(argptr);

  va_end(argptr);

  return ier;
}

void MMG2D_Init_fileNames(MMG5_pMesh mesh,MMG5_pSol sol
  ) {

  MMG5_Init_fileNames(mesh,sol);
  return;
}

int MMG2D_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin) {

  return MMG5_Set_inputMeshName(mesh,meshin);
}

int MMG2D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin) {
  return MMG5_Set_inputSolName(mesh,sol,solin);
}

int MMG2D_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout) {

  return MMG5_Set_outputMeshName(mesh,meshout);
}

int MMG2D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout) {
  return MMG5_Set_outputSolName(mesh,sol,solout);
}

void MMG2D_Init_parameters(MMG5_pMesh mesh) {

  /* Init common parameters for mmg2d, mmgs and mmg2d. */
  MMG5_Init_parameters(mesh);

  /* default values for integers */
  mesh->info.lag      = MMG5_LAG;
  mesh->info.setfem   = MMG5_FEM;
  mesh->info.optim    = MMG5_OFF;
  /* [0/1]   ,avoid/allow surface modifications */
  mesh->info.nosurf   =  MMG5_OFF;
  /* [0/1/2] , set 3D mode for 2D mesh  */
  mesh->info.renum    = MMG5_OFF;
  /* [0/1]   , set off/on normal regularization */
  mesh->info.nreg     = MMG5_OFF;
  /* [0/1]   , set off/on vertex regularization */
  mesh->info.xreg     = MMG5_OFF;
  /* default values for doubles */
  /* level set value */
  mesh->info.ls       = MMG5_LS;

  /* Ridge detection */
  mesh->info.dhd      = MMG5_ANGEDG;
}


int MMG2D_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam, MMG5_int val){
  int k;

  switch ( iparam ) {
    /* Integer parameters */
  case MMG2D_IPARAM_verbose :
    mesh->info.imprim   = val;
    break;
  case MMG2D_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf(stderr,"\n  ## Warning: %s: maximal memory authorized must"
              " be strictly positive.\n",__func__);
      fprintf(stderr,"  Reset to default value.\n");
    }
    else
      mesh->info.mem      = val;
    if ( !MMG2D_memOption(mesh) ) return 0;
    break;
  case MMG2D_IPARAM_debug :
    mesh->info.ddebug   = val;
    break;
  case MMG2D_IPARAM_angle :
    /* free table that may contains old ridges */
    if ( mesh->htab.geom )
      MMG5_DEL_MEM(mesh,mesh->htab.geom);
    if ( mesh->xpoint )
      MMG5_DEL_MEM(mesh,mesh->xpoint);
    if ( mesh->xtetra )
      MMG5_DEL_MEM(mesh,mesh->xtetra);
    if ( !val )
      mesh->info.dhd    = MMG5_NR;
    else {
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: angle detection parameter"
                " set to default value\n",__func__);
      mesh->info.dhd    = MMG5_ANGEDG;
    }
    break;
  case MMG2D_IPARAM_nofem :
    mesh->info.setfem = (val==1)? 0 : 1;
    break;
  case MMG2D_IPARAM_opnbdy :
    mesh->info.opnbdy = val;
    break;
  case MMG2D_IPARAM_iso :
    mesh->info.iso      = val;
    break;
  case MMG2D_IPARAM_isoref :
    mesh->info.isoref   = val;
    break;
  case MMG2D_IPARAM_isosurf :
    mesh->info.isosurf = val;
    break;
  case MMG2D_IPARAM_lag :
#ifdef USE_ELAS
    if ( val < 0 || val > 2 )
      return 0;
    mesh->info.lag = val;
    /* No connectivity changes unless lag >= 2 */
    if ( val < 2 ) {
      if ( !MMG2D_Set_iparameter(mesh,sol,MMG2D_IPARAM_noinsert,1) )
        return 0;
    }
#else
    fprintf(stderr,"\n  ## Error: %s"
            " \"lagrangian motion\" option unavailable (-lag):\n"
            " set the USE_ELAS CMake's flag to ON when compiling the mmg3d"
            " library to enable this feature.\n",__func__);
    return 0;
#endif
    break;
  case MMG2D_IPARAM_3dMedit :
    mesh->info.renum = val;
    break;
  case MMG2D_IPARAM_numsubdomain :
    mesh->info.nsd = val;
    break;
  case MMG2D_IPARAM_optim :
    mesh->info.optim = val;
    break;
  case MMG2D_IPARAM_noinsert :
    mesh->info.noinsert = val;
    break;
  case MMG2D_IPARAM_noswap :
    mesh->info.noswap   = val;
    break;
  case MMG2D_IPARAM_nomove :
    mesh->info.nomove   = val;
    break;
  case MMG2D_IPARAM_nosurf :
    mesh->info.nosurf   = val;
    break;
  case MMG2D_IPARAM_nreg :
    mesh->info.nreg     = val;
    break;
  case MMG2D_IPARAM_xreg :
    mesh->info.xreg     = val;
    break;
  case MMG2D_IPARAM_nosizreq :
    mesh->info.nosizreq = val;
    break;
  case MMG2D_IPARAM_numberOfLocalParam :
    if ( mesh->info.par ) {
      MMG5_DEL_MEM(mesh,mesh->info.par);
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: new local parameter values\n",__func__);
    }
    mesh->info.npar   = val;
    mesh->info.npari  = 0;
    mesh->info.parTyp = 0;

    MMG5_ADD_MEM(mesh,mesh->info.npar*sizeof(MMG5_Par),"parameters",
                  printf("  Exit program.\n");
                  return 0);
    MMG5_SAFE_CALLOC(mesh->info.par,mesh->info.npar,MMG5_Par,return 0);

    MMG5_int inival;
    if ( sizeof(MMG5_int) == 8 ) {
      inival = LONG_MAX;
    }
    else {
      inival = INT_MAX;
    }

    for (k=0; k<mesh->info.npar; k++) {
      mesh->info.par[k].elt   = MMG5_Noentity;
      mesh->info.par[k].ref   = inival;
      mesh->info.par[k].hausd = mesh->info.hausd;
      mesh->info.par[k].hmin  = mesh->info.hmin;
      mesh->info.par[k].hmax  = mesh->info.hmax;
    }
    break;
  case MMG2D_IPARAM_numberOfLSBaseReferences :
    if ( mesh->info.br ) {
      MMG5_DEL_MEM(mesh,mesh->info.br);
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: new level-set based references values\n",__func__);
    }
    mesh->info.nbr   = val;
    mesh->info.nbri  = 0;
    MMG5_ADD_MEM(mesh,mesh->info.nbr*sizeof(MMG5_int),"References",
                  printf("  Exit program.\n");
                  return 0);
    MMG5_SAFE_CALLOC(mesh->info.br,mesh->info.nbr,MMG5_int,return 0);

    for (k=0; k<mesh->info.nbr; k++)
      mesh->info.br[k] = 0;

    break;

  case MMG2D_IPARAM_numberOfMat :
    if ( mesh->info.mat ) {
      MMG5_DEL_MEM(mesh,mesh->info.mat);
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: new multi materials values\n",__func__);
    }
    mesh->info.nmat   = val;
    mesh->info.nmati  = 0;

    MMG5_ADD_MEM(mesh,(mesh->info.nmat)*sizeof(MMG5_Mat),"multi material",
                 printf("  Exit program.\n");
                 return 0);
    MMG5_SAFE_CALLOC(mesh->info.mat,mesh->info.nmat,MMG5_Mat,return 0);
    break;

  case MMG2D_IPARAM_anisosize :
    mesh->info.ani = val;
    break;
  default :
    fprintf(stderr,"\n  ## Error: %s: unknown type of parameter\n",__func__);
    return 0;
  }
  /* other options */

  return 1;
}

int MMG2D_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val){

  switch ( dparam ) {
    /* double parameters */
  case MMG2D_DPARAM_angleDetection :
    mesh->info.dhd = val;
    mesh->info.dhd = MG_MAX(0.0, MG_MIN(180.0,mesh->info.dhd));
    mesh->info.dhd = cos(mesh->info.dhd*M_PI/180.0);
    break;
  case MMG2D_DPARAM_hmin :
    mesh->info.sethmin  = 1;
    mesh->info.hmin     = val;
    if ( mesh->info.sethmax && ( mesh->info.hmin >=  mesh->info.hmax ) ) {
      fprintf(stderr,"\n  ## Error: hmin value must be strictly lower than hmax one"
              " (hmin = %lf  hmax = %lf ).\n",mesh->info.hmin, mesh->info.hmax);
      return 0;
    }
    if ( val <= 0. ) {
      fprintf(stderr,"\n  ## Error: hmin must be strictly positive "
              "(minimal edge length).\n");
      return 0;
    }

    break;
  case MMG2D_DPARAM_hmax :
    mesh->info.sethmax  = 1;
    mesh->info.hmax     = val;
    if ( mesh->info.sethmin && ( mesh->info.hmin >=  mesh->info.hmax ) ) {
      fprintf(stderr,"\n  ## Error: hmin value must be strictly lower than hmax one"
              " (hmin = %lf  hmax = %lf ).\n",mesh->info.hmin, mesh->info.hmax);
      return 0;
    }
    if ( val <= 0. ) {
      fprintf(stderr,"\n  ## Error: hmax must be strictly positive "
              "(maximal edge length).\n");
      return 0;
    }

    break;
  case MMG2D_DPARAM_hsiz :
    mesh->info.hsiz     = val;
    break;
  case MMG2D_DPARAM_hgrad :
    mesh->info.hgrad    = val;
    if ( mesh->info.hgrad <= 0.0 ) {
      mesh->info.hgrad = MMG5_NOHGRAD;
    }
    else {
      mesh->info.hgrad = log(mesh->info.hgrad);
    }
    break;
  case MMG2D_DPARAM_hgradreq :
    mesh->info.hgradreq    = val;
    if ( mesh->info.hgradreq <= 0.0 ) {
      mesh->info.hgradreq = MMG5_NOHGRAD;
    }
    else {
      mesh->info.hgradreq = log(mesh->info.hgradreq);
    }
    break;
  case MMG2D_DPARAM_hausd :
    if ( val <=0 ) {
      fprintf(stderr,"\n  ## Error: %s: hausdorff number must be"
              " strictly positive.\n",__func__);
      return 0;
    }
    else
      mesh->info.hausd    = val;
    break;
  case MMG2D_DPARAM_ls :
    mesh->info.ls       = val;
    break;
  case MMG2D_DPARAM_rmc :
    if ( !val ) {
      /* Default value */
      mesh->info.rmc      = MMG5_VOLFRAC;
    }
    else {
      /* User customized value */
      mesh->info.rmc      = val;
    }
    break;
  default :
    fprintf(stderr,"\n  ## Error: %s: unknown type of parameter\n",
            __func__);
    return 0;
  }
  return 1;
}

int MMG2D_Set_localParameter(MMG5_pMesh mesh,MMG5_pSol sol, int typ, MMG5_int ref,
                             double hmin,double hmax,double hausd){
  MMG5_pPar par;
  int k;

  if ( !mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of local"
            " parameters",__func__);
    fprintf(stderr," with the MMG2D_Set_iparameters function before setting");
    fprintf(stderr," values in local parameters structure. \n");
    return 0;
  }
  if ( mesh->info.npari >= mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new local parameter.\n",
            __func__);
    fprintf(stderr,"    max number of local parameters: %d\n",mesh->info.npar);
    return 0;
  }
  if ( typ != MMG5_Triangle && typ != MMG5_Edg ) {
    fprintf(stderr,"\n  ## Warning: %s: you must apply your local parameters",
            __func__);
    fprintf(stderr," on triangles (MMG5_Triangle or %d) or edges"
            " (MMG5_Edg or %d).\n",MMG5_Triangle,MMG5_Edg);
    fprintf(stderr,"\n  ## Unknown type of entity: ignored.\n");
    return 0;
  }
  if ( ref < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative references are not allowed.\n",
            __func__);
    return 0;
  }
  if ( hmin <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative hmin value is not allowed.\n",
            __func__);
    return 0;
  }
  if ( hmax <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative hmax value is not allowed.\n",
            __func__);
    return 0;
  }
  if ( hausd <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative hausd value is not allowed.\n",
            __func__);
    return 0;
  }

  for (k=0; k<mesh->info.npari; k++) {
    par = &mesh->info.par[k];

    if ( par->elt == typ && par->ref == ref ) {
      par->hausd = hausd;
      par->hmin  = hmin;
      par->hmax  = hmax;
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
        fprintf(stderr,"\n  ## Warning: %s: new parameters (hausd, hmin and hmax)",
                __func__);
        fprintf(stderr," for entities of type %d and of ref %" MMG5_PRId "\n",typ,ref);
      }
      return 1;
    }
  }

  mesh->info.par[mesh->info.npari].elt   = typ;
  mesh->info.par[mesh->info.npari].ref   = ref;
  mesh->info.par[mesh->info.npari].hmin  = hmin;
  mesh->info.par[mesh->info.npari].hmax  = hmax;
  mesh->info.par[mesh->info.npari].hausd = hausd;

  switch ( typ )
  {
  case ( MMG5_Triangle ):
    mesh->info.parTyp |= MG_Tria;
    break;
  case ( MMG5_Edg ):
    mesh->info.parTyp |= MG_Edge;
    break;
  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected entity type: %s.\n",
            __func__,MMG5_Get_entitiesName(typ));
    return 0;
  }

  mesh->info.npari++;

  return 1;
}

int MMG2D_Set_multiMat(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int ref,
                       int split,MMG5_int rin,MMG5_int rout){
  return MMG5_Set_multiMat(mesh,sol,ref,split,rin,rout);
}

int MMG2D_Set_lsBaseReference(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int br){
  return MMG5_Set_lsBaseReference(mesh,sol,br);
}


int MMG2D_Set_meshSize(MMG5_pMesh mesh, MMG5_int np, MMG5_int nt, MMG5_int nquad, MMG5_int na) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->tria || mesh->edge) )
    fprintf(stderr,"\n  ## Warning: %s: old mesh deletion.\n",__func__);

  if ( mesh->point )
    MMG5_DEL_MEM(mesh,mesh->point);
  if ( mesh->tria )
    MMG5_DEL_MEM(mesh,mesh->tria);
  if ( mesh->quadra )
    MMG5_DEL_MEM(mesh,mesh->quadra);
  if ( mesh->edge )
    MMG5_DEL_MEM(mesh,mesh->edge);

  mesh->np  = np;
  mesh->nt  = nt;
  mesh->nquad  = nquad;
  mesh->na  = na;
  mesh->npi = mesh->np;
  mesh->nti = mesh->nt;
  mesh->nai = mesh->na;

  /*tester si -m definie : renvoie 0 si pas ok et met la taille min dans info.mem */
  if( mesh->info.mem > 0) {
    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->namax < mesh->na) ) {
      if ( !MMG2D_memOption(mesh) )  return 0;
    } else if(mesh->info.mem < 39) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory (%d).\n",
              __func__,mesh->info.mem);
      return 0;
    }
  } else {
    if ( !MMG2D_memOption(mesh) )  return 0;
  }

  /* Mesh allocation and linkage */
  if ( !MMG2D_setMeshSize_alloc( mesh ) ) return 0;

  return 1;
}

int MMG2D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, MMG5_int np, int typSol) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && sol->m )
    fprintf(stderr,"\n  ## Warning: %s: old solution deletion.\n",__func__);

  if ( typEntity != MMG5_Vertex ) {
    fprintf(stderr,"\n  ## Error: %s: mmg2d need a solution imposed on vertices.\n",
            __func__);
    return 0;
  }

  sol->type = typSol;

  if ( typSol == MMG5_Scalar ) {
    sol->size = 1;
  }
  else if ( typSol == MMG5_Vector ) {
    sol->size = 2;
  }
  else if ( typSol == MMG5_Tensor ) {
    sol->size = 3;
  }
  else {
    fprintf(stderr,"\n  ## Error: %s: type of solution not yet implemented.\n",
            __func__);
    return 0;
  }

  sol->dim = 2;
  if ( np ) {
    sol->np  = np;
    sol->npi = np;
    if ( sol->m )
      MMG5_DEL_MEM(mesh,sol->m);

    sol->npmax = mesh->npmax;
    MMG5_ADD_MEM(mesh,(sol->size*(sol->npmax+1))*sizeof(double),"initial solution",
                  printf("  Exit program.\n");
                  return 0);
    MMG5_SAFE_CALLOC(sol->m,(sol->size*(sol->npmax+1)),double,return 0);
  }
  return 1;
}

int MMG2D_Set_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol *sol,int nsols,
                                 MMG5_int np, int *typSol) {
  MMG5_pSol psl;
  char      data[18];
  int       j;

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && mesh->nsols ) {
    if ( *sol ) {
      fprintf(stderr,"\n  ## Warning: %s: old solutions array deletion.\n",
              __func__);
      MMG5_DEL_MEM(mesh,*sol);
    }
  }

  /** Sol tab allocation */
  mesh->nsols = nsols;

  MMG5_ADD_MEM(mesh,nsols*sizeof(MMG5_Sol),"solutions array",
                return 0);
  MMG5_SAFE_CALLOC(*sol,nsols,MMG5_Sol,return 0);

  for ( j=0; j<nsols; ++j ) {
    psl = *sol + j;
    psl->ver = 2;

    /* Give an arbitrary name to the solution */
    sprintf(data,"sol_%d",j);
    if ( !MMG2D_Set_inputSolName(mesh,psl,data) ) {
      return 0;
    }
    /* Give an arbitrary name to the solution */
    sprintf(data,"sol_%d.o",j);
    if ( !MMG2D_Set_outputSolName(mesh,psl,data) ) {
      return 0;
    }

    if ( !MMG2D_Set_solSize(mesh,psl,MMG5_Vertex,mesh->np,typSol[j]) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to set the size of the"
              " solution num %d.\n",__func__,j);
      return 0;
    }
  }
  return 1;
}

int MMG2D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, MMG5_int* np,
                      int* typSol) {

  if ( typEntity != NULL )
    *typEntity = MMG5_Vertex;

  if ( typSol != NULL ) {
    if ( sol->size == 1 )
      *typSol    = MMG5_Scalar;
    else if ( sol->size == 2 )
      *typSol    = MMG5_Vector;
    else if ( sol->size == 3 )
      *typSol    = MMG5_Tensor;
    else
      *typSol    = MMG5_Notype;
  }

  assert( (!sol->np) || (sol->np == mesh->np));

  if ( np != NULL )
    *np = sol->np;

  return 1;
}

int MMG2D_Get_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol *sol, int *nsols,
                                 MMG5_int* np, int* typSol) {
  MMG5_pSol psl;
  MMG5_int  j;

  if ( !mesh ) {
    fprintf(stderr,"\n  ## Error: %s: your mesh structure must be allocated"
            " and filled\n",__func__);
    return 0;
  }

  if ( nsols != NULL )
    *nsols = mesh->nsols;

  for ( j=0; j<mesh->nsols; ++j ) {
    psl = *sol + j;

    if ( typSol != NULL ) {
      typSol[j]    = psl->type;
    }

    assert( (!psl->np) || (psl->np == mesh->np));
  }
  if ( np != NULL )
    *np = mesh->np;

  return 1;
}

int MMG2D_Get_meshSize(MMG5_pMesh mesh, MMG5_int* np, MMG5_int* nt, MMG5_int* nquad, MMG5_int* na) {
  MMG5_int k;

  if ( np != NULL )
    *np = mesh->np;
  if ( nt != NULL )
    *nt = mesh->nt;

  if ( nquad != NULL )
    *nquad = mesh->nquad;

  if ( na != NULL ) {
    // Edges are not packed, thus we must count it.
    *na = 0;
    if ( mesh->na ) {
      for (k=1; k<=mesh->na; k++) {
        if ( mesh->edge[k].a ) ++(*na);
      }
    }
  }

  return 1;
}

int MMG2D_Set_vertex(MMG5_pMesh mesh, double c0, double c1, MMG5_int ref, MMG5_int pos) {

  if ( !mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: you must set the number of points with the",
            __func__);
    fprintf(stderr," MMG2D_Set_meshSize function before setting vertices in mesh\n");
    return 0;
  }

  if ( pos > mesh->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new point.\n",
            __func__);
    fprintf(stderr,"    max number of points: %" MMG5_PRId "\n",mesh->npmax);
    MMG5_INCREASE_MEM_MESSAGE();
    return 0;
  }

  if ( pos > mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new vertex at position %" MMG5_PRId ".",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of vertices: %" MMG5_PRId "\n",mesh->np);
    fprintf(stderr,"  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the vertex.\n");
    return 0;
  }

  mesh->point[pos].c[0] = c0;
  mesh->point[pos].c[1] = c1;
  mesh->point[pos].ref  = ref;
  if ( mesh->nt )
    mesh->point[pos].tag  = MG_NUL;
  else
    mesh->point[pos].tag  &= ~MG_NUL;

  mesh->point[pos].flag = 0;
  mesh->point[pos].tmp = 0;

  return 1;
}

int MMG2D_Set_corner(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_CRN;
  return 1;
}

int MMG2D_Unset_corner(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag &= ~MG_CRN;
  return 1;
}

int MMG2D_Set_requiredVertex(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_REQ;
  mesh->point[k].tag &= ~MG_NUL;
  return 1;
}

int MMG2D_Unset_requiredVertex(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag &= ~MG_REQ;
  return 1;
}

int MMG2D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, MMG5_int* ref,
                     int* isCorner, int* isRequired) {

  if ( mesh->npi == mesh->np ) {
    mesh->npi = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG2D_Get_vertex function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",mesh->np);
    }
  }

  mesh->npi++;

  if ( mesh->npi > mesh->np ) {
    fprintf(stderr,"  ## Error: %s: unable to get point.\n",__func__);
    fprintf(stderr,"     The number of call of MMG2D_Get_vertex function");
    fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",mesh->np);
    return 0;
  }

  return MMG2D_GetByIdx_vertex( mesh,c0,c1,ref,isCorner,isRequired,mesh->npi);
}

int MMG2D_GetByIdx_vertex(MMG5_pMesh mesh, double* c0, double* c1, MMG5_int* ref,
                          int* isCorner, int* isRequired, MMG5_int idx) {

  if ( idx < 1 || idx > mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get point at position %" MMG5_PRId ".\n",
            __func__,idx);
    fprintf(stderr,"     Your vertices numbering goes from 1 to %" MMG5_PRId "\n",mesh->np);
    return 0;
  }

  *c0  = mesh->point[idx].c[0];
  *c1  = mesh->point[idx].c[1];
  if ( ref != NULL )
    *ref = mesh->point[idx].ref;

  if ( isCorner != NULL ) {
    if ( mesh->point[idx].tag & MG_CRN )
      *isCorner = 1;
    else
      *isCorner = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->point[idx].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}

int  MMG2D_Set_vertices(MMG5_pMesh mesh, double *vertices,MMG5_int *refs) {
  MMG5_pPoint ppt;
  MMG5_int i,j;

  /*coordinates vertices*/
  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*2;
    ppt->c[0]  = vertices[j];
    ppt->c[1]  = vertices[j+1];

    ppt->flag = 0;
    ppt->tmp = 0;

    if ( refs != NULL )
      ppt->ref   = refs[i-1];

    if ( mesh->nt )
      ppt->tag  = MG_NUL;
    else
      ppt->tag  &= ~MG_NUL;
  }

  return 1;
}

int  MMG2D_Get_vertices(MMG5_pMesh mesh, double* vertices, MMG5_int* refs,
                        int* areCorners, int* areRequired) {
  MMG5_pPoint ppt;
  MMG5_int i,j;

  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*2;
    vertices[j] = ppt->c[0];
    vertices[j+1] = ppt->c[1];

    j = i-1;
    if ( refs != NULL )
      refs[j] = ppt->ref;

    if ( areCorners !=NULL ) {
      if ( ppt->tag & MG_CRN )
        areCorners[j] = 1;
      else
        areCorners[j] = 0;
    }

    if ( areRequired != NULL ) {
      if ( ppt->tag & MG_REQ )
        areRequired[j] = 1;
      else
        areRequired[j] = 0;
    }
  }

  return 1;
}

int MMG2D_Set_triangle(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2, MMG5_int ref, MMG5_int pos) {
  MMG5_pPoint ppt;
  MMG5_pTria  pt;
  double      vol;
  int         i,j,ip;
  MMG5_int    tmp;

  if ( !mesh->nt ) {
    fprintf(stderr,"  ## Error: %s: You must set the number of elements with the",
            __func__);
    fprintf(stderr," MMG2D_Set_meshSize function before setting elements in mesh\n");
    return 0;
  }

  if ( pos > mesh->ntmax ) {
    fprintf(stderr,"  ## Error: %s: unable to allocate a new element.\n",
            __func__);
    fprintf(stderr,"    max number of element: %" MMG5_PRId "\n",mesh->ntmax);
    MMG5_INCREASE_MEM_MESSAGE();
    return 0;
  }

  if ( pos > mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new triangle at position %" MMG5_PRId ".",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of triangle: %" MMG5_PRId "\n",mesh->nt);
    fprintf(stderr,"  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the triangle.\n");
    return 0;
  }

  pt = &mesh->tria[pos];
  pt->v[0] = v0;
  pt->v[1] = v1;
  pt->v[2] = v2;
  pt->ref  = ref;

  mesh->point[pt->v[0]].tag &= ~MG_NUL;
  mesh->point[pt->v[1]].tag &= ~MG_NUL;
  mesh->point[pt->v[2]].tag &= ~MG_NUL;

  for(i=0 ; i<3 ; i++)
    pt->edg[i] = 0;

  vol = MMG2D_quickarea(mesh->point[pt->v[0]].c,mesh->point[pt->v[1]].c,
                       mesh->point[pt->v[2]].c);

  if ( vol == 0.0 ) {
    fprintf(stderr,"\n  ## Error: %s: triangle %" MMG5_PRId " has null area.\n",
            __func__,pos);
    for ( ip=0; ip<3; ip++ ) {
      ppt = &mesh->point[pt->v[ip]];
      for ( j=0; j<3; j++ ) {
        if ( fabs(ppt->c[j])>0. ) {
          fprintf(stderr," Check that you don't have a sliver triangle.\n");
          return 0;
        }
      }
    }
  }
  else if(vol < 0) {
    tmp = pt->v[2];
    pt->v[2] = pt->v[1];
    pt->v[1] = tmp;
    /* mesh->xt temporary used to count reoriented tria */
    mesh->xt++;
  }
  if ( mesh->info.ddebug && (mesh->nt == pos) && mesh->xt > 0 ) {
    fprintf(stderr,"\n  ## Warning: %s: %" MMG5_PRId " triangles reoriented\n",
            __func__,mesh->xt);
    mesh->xt = 0;
  }

  return 1;
}

int MMG2D_Set_requiredTriangle(MMG5_pMesh mesh, MMG5_int k) {
  MMG5_pTria pt;

  assert ( k <= mesh->nt );
  pt = &mesh->tria[k];

  pt->tag[0] |= MG_REQ;
  pt->tag[1] |= MG_REQ;
  pt->tag[2] |= MG_REQ;

  return 1;
}

int MMG2D_Unset_requiredTriangle(MMG5_pMesh mesh, MMG5_int k) {
  MMG5_pTria pt;

  assert ( k <= mesh->nt );
  pt = &mesh->tria[k];

  pt->tag[0] &= ~MG_REQ;
  pt->tag[1] &= ~MG_REQ;
  pt->tag[2] &= ~MG_REQ;

  return 1;
}

int MMG2D_Get_triangle(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* ref
                       ,int* isRequired) {
  MMG5_pTria  ptt;

  if ( mesh->nti == mesh->nt ) {
    mesh->nti = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of"
              " triangles.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG2D_Get_triangle function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of triangles: %" MMG5_PRId "\n ",mesh->nt);
    }
  }

  mesh->nti++;

  if ( mesh->nti > mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get triangle.\n",
            __func__);
    fprintf(stderr,"    The number of call of MMG2D_Get_triangle function");
    fprintf(stderr," can not exceed the number of triangles: %" MMG5_PRId "\n ",mesh->nt);
    return 0;
  }

  ptt = &mesh->tria[mesh->nti];

  *v0  = ptt->v[0];
  *v1  = ptt->v[1];
  *v2  = ptt->v[2];
  if ( ref != NULL )
    *ref = ptt->ref;

  if ( isRequired != NULL ) {
    if ( (ptt->tag[0] & MG_REQ) && (ptt->tag[1] & MG_REQ) &&
         (ptt->tag[2] & MG_REQ) )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}

int  MMG2D_Set_triangles(MMG5_pMesh mesh, MMG5_int *tria, MMG5_int *refs) {
  MMG5_pPoint ppt;
  MMG5_pTria ptt;
  double vol;
  int ip;
  MMG5_int  j,i,tmp;

  mesh->xt = 0;
  for (i=1;i<=mesh->nt;i++)
  {
    j = (i-1)*3;
    ptt = &mesh->tria[i];
    ptt->v[0] = tria[j]  ;
    ptt->v[1] = tria[j+2];
    ptt->v[2] = tria[j+1];
    if ( refs != NULL )
      ptt->ref  = refs[i-1];

    mesh->point[ptt->v[0]].tag &= ~MG_NUL;
    mesh->point[ptt->v[1]].tag &= ~MG_NUL;
    mesh->point[ptt->v[2]].tag &= ~MG_NUL;

    int ii;
    for( ii=0 ; ii<3 ; ii++)
      ptt->edg[ii] = 0;

    vol = MMG2D_quickarea(mesh->point[ptt->v[0]].c,mesh->point[ptt->v[1]].c,
                         mesh->point[ptt->v[2]].c);

    if ( vol == 0.0 ) {
      fprintf(stderr,"\n  ## Error: %s: triangle %" MMG5_PRId " has null area.\n",
              __func__,i);
      for ( ip=0; ip<3; ip++ ) {
        ppt = &mesh->point[ptt->v[ip]];
        int jj;
        for ( jj=0; jj<3; jj++ ) {
          if ( fabs(ppt->c[jj])>0. ) {
            fprintf(stderr," Check that you don't have a sliver triangle.\n");
            return 0;
          }
        }
      }
    }
    else if(vol < 0) {
      tmp = ptt->v[2];
      ptt->v[2] = ptt->v[1];
      ptt->v[1] = tmp;
      /* mesh->xt temporary used to count reoriented triangles */
      mesh->xt++;
    }
    if ( mesh->info.ddebug && mesh->xt > 0 ) {
      fprintf(stderr,"\n  ## Warning: %s: %" MMG5_PRId " triangles reoriented\n",
              __func__,mesh->xt);
    }
  }
  return 1;
}

int  MMG2D_Get_triangles(MMG5_pMesh mesh, MMG5_int* tria, MMG5_int* refs,
                         int* areRequired) {
  MMG5_pTria ptt;
  MMG5_int   i, j;

  for (i=1;i<=mesh->nt;i++)
  {
    j = (i-1)*3;
    ptt = &mesh->tria[i];
    tria[j]   = ptt->v[0];
    tria[j+1] = ptt->v[1];
    tria[j+2] = ptt->v[2];

    if ( refs!=NULL )
      refs[i-1]  = ptt->ref ;
    if ( areRequired != NULL ) {
      if ( (ptt->tag[0] & MG_REQ) && (ptt->tag[1] & MG_REQ) &&
           (ptt->tag[2] & MG_REQ) )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}

int MMG2D_Set_quadrilateral(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2, MMG5_int v3, MMG5_int ref, MMG5_int pos) {
  MMG5_pQuad  pq;

  if ( !mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of quadrilaterals with the",
            __func__);
    fprintf(stderr," MMG2D_Set_meshSize function before setting elements in mesh\n");
    return 0;
  }

  if ( pos > mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new quad at position %" MMG5_PRId ".",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of quads: %" MMG5_PRId "\n",mesh->nquad);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the quad.\n");
    return 0;
  }

  pq = &mesh->quadra[pos];
  pq->v[0] = v0;
  pq->v[1] = v1;
  pq->v[2] = v2;
  pq->v[3] = v3;
  pq->ref  = ref;

  mesh->point[pq->v[0]].tag &= ~MG_NUL;
  mesh->point[pq->v[1]].tag &= ~MG_NUL;
  mesh->point[pq->v[2]].tag &= ~MG_NUL;
  mesh->point[pq->v[3]].tag &= ~MG_NUL;

  return 1;
}

int MMG2D_Get_quadrilateral(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* v3,
                            MMG5_int* ref, int* isRequired) {
  static MMG5_int nqi = 0;

  if ( nqi == mesh->nquad ) {
    nqi = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of"
              " quadrilaterals.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG2D_Get_quadrilateral function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of quadrilaterals: %" MMG5_PRId "\n ",mesh->nquad);
    }
  }

  ++nqi;

  if ( nqi > mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get quadra.\n",__func__);
    fprintf(stderr,"    The number of call of MMG2D_Get_quadrilateral function");
    fprintf(stderr," can not exceed the number of quadra: %" MMG5_PRId "\n ",mesh->nquad);
    return 0;
  }

  *v0  = mesh->quadra[nqi].v[0];
  *v1  = mesh->quadra[nqi].v[1];
  *v2  = mesh->quadra[nqi].v[2];
  *v3  = mesh->quadra[nqi].v[3];

  if ( ref != NULL ) {
    *ref = mesh->quadra[nqi].ref;
  }

  if ( isRequired != NULL ) {
    if ( ( mesh->quadra[nqi].tag[0] & MG_REQ) && ( mesh->quadra[nqi].tag[1] & MG_REQ) &&
         ( mesh->quadra[nqi].tag[2] & MG_REQ) && ( mesh->quadra[nqi].tag[3] & MG_REQ) ) {
      *isRequired = 1;
    }
    else
      *isRequired = 0;
  }

  return 1;
}

int  MMG2D_Set_quadrilaterals(MMG5_pMesh mesh, MMG5_int *quadra, MMG5_int *refs) {
  MMG5_pQuad pq;
  MMG5_int   i,j;

  for (i=1;i<=mesh->nquad;i++)
  {
    j = (i-1)*4;
    pq = &mesh->quadra[i];
    pq->v[0]  = quadra[j];
    pq->v[1]  = quadra[j+1];
    pq->v[2]  = quadra[j+2];
    pq->v[3]  = quadra[j+3];

    if ( refs != NULL )
      pq->ref   = refs[i-1];

    mesh->point[pq->v[0]].tag &= ~MG_NUL;
    mesh->point[pq->v[1]].tag &= ~MG_NUL;
    mesh->point[pq->v[2]].tag &= ~MG_NUL;
    mesh->point[pq->v[3]].tag &= ~MG_NUL;
  }

  return 1;
}

int  MMG2D_Get_quadrilaterals(MMG5_pMesh mesh, MMG5_int *quadra, MMG5_int *refs, int * areRequired) {
  MMG5_pQuad pq;
  MMG5_int   i, j;

  for (i=1;i<=mesh->nquad;i++)
  {
    j = (i-1)*4;
    pq = &mesh->quadra[i];
    quadra[j]   = pq->v[0];
    quadra[j+1] = pq->v[1];
    quadra[j+2] = pq->v[2];
    quadra[j+3] = pq->v[3];

    if ( refs!=NULL )
      refs[i-1]  = pq->ref ;

    if ( areRequired != NULL ) {
     if ( (pq->tag[0] & MG_REQ) && (pq->tag[1] & MG_REQ) &&
          (pq->tag[2] & MG_REQ) && (pq->tag[3] & MG_REQ) ) {
       areRequired[i-1] = 1;
     }
      else
        areRequired[i-1] = 0;
    }

  }
  return 1;
}


int MMG2D_Set_edge(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int ref, MMG5_int pos) {
  MMG5_pEdge pt;

  if ( !mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: you must set the number of elements"
            " with the",__func__);
    fprintf(stderr," MMG2D_Set_meshSize function before setting elements in mesh\n");
    return 0;
  }

  if ( pos > mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new edge at position %" MMG5_PRId ".",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of edge: %" MMG5_PRId "\n",mesh->na);
    fprintf(stderr,"  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the edge.\n");
    return 0;
  }

  pt = &mesh->edge[pos];
  pt->a = v0;
  pt->b = v1;
  pt->ref  = ref;
  pt->tag &= MG_REF + MG_BDY;

  mesh->point[pt->a].tag &= ~MG_NUL;
  mesh->point[pt->b].tag &= ~MG_NUL;

  return 1;
}

int MMG2D_Set_requiredEdge(MMG5_pMesh mesh, MMG5_int k) {
  MMG5_pEdge  ped;

  assert ( k <= mesh->na );

  ped = &mesh->edge[k];

  ped->tag |= MG_REQ;

  return 1;
}

int MMG2D_Unset_requiredEdge(MMG5_pMesh mesh, MMG5_int k) {
  MMG5_pEdge  ped;

  assert ( k <= mesh->na );

  ped = &mesh->edge[k];

  ped->tag &= ~MG_REQ;

  return 1;
}

int MMG2D_Set_parallelEdge(MMG5_pMesh mesh, MMG5_int k) {
  MMG5_pPoint ppt;
  MMG5_pEdge  ped;

  assert ( k <= mesh->na );

  ped = &mesh->edge[k];

  ped->tag |= MG_PARBDY;

  ppt = &mesh->point[ped->a];
  ppt->tag |= MG_PARBDY;
  ppt = &mesh->point[ped->b];
  ppt->tag |= MG_PARBDY;

  return 1;
}

int MMG2D_Get_edge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref
                   ,int* isRidge, int* isRequired) {
  MMG5_pEdge        ped;

  if ( mesh->nai == mesh->na ) {
    mesh->nai = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of edges.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG2D_Get_edge function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of edges.\n ");
      fprintf(stderr,"     Please, call the MMG2D_Get_meshSize function to get"
              " this number.\n ");
    }
  }

  mesh->nai++;

  if ( mesh->nai > mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get edge.\n",__func__);
    fprintf(stderr,"    The number of call of MMG2D_Get_edge function");
    fprintf(stderr," can not exceed the number of edges: %" MMG5_PRId "\n ",mesh->na);
    return 0;
  }

  ped = &mesh->edge[mesh->nai];

  while ( !ped->a && ++mesh->nai <= mesh->na ) {
    ped = &mesh->edge[mesh->nai];
  }

  *e0  = ped->a;
  *e1  = ped->b;

  if ( ref != NULL )
    *ref = mesh->edge[mesh->nai].ref;

  if ( isRidge != NULL ) {
    if ( mesh->edge[mesh->nai].tag & MG_GEO )
      *isRidge = 1;
    else
      *isRidge = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->edge[mesh->nai].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}

int MMG2D_Set_edges(MMG5_pMesh mesh, MMG5_int *edges, MMG5_int *refs) {
  MMG5_int i,j;

  for (i=1;i<=mesh->na;i++)
  {
    j = (i-1)*2;

    mesh->edge[i].a    = edges[j];
    mesh->edge[i].b    = edges[j+1];
    if ( refs != NULL )
      mesh->edge[i].ref  = refs[i-1];
    mesh->edge[i].tag &= MG_REF+MG_BDY;

    mesh->point[mesh->edge[i].a].tag &= ~MG_NUL;
    mesh->point[mesh->edge[i].b].tag &= ~MG_NUL;
  }

  return 1;
}

int MMG2D_Get_edges(MMG5_pMesh mesh, MMG5_int* edges,MMG5_int *refs,int* areRidges,int* areRequired) {
  MMG5_int i,j;

  for (i=1;i<=mesh->na;i++)
  {
    j = (i-1)*2;
    edges[j]   = mesh->edge[i].a;
    edges[j+1] = mesh->edge[i].b;

    if ( refs!=NULL )
      refs[i-1] = mesh->edge[i].ref;

    if ( areRidges != NULL ) {
      if ( mesh->edge[i].tag & MG_GEO )
        areRidges[i-1] = 1;
      else
        areRidges[i-1] = 0;
    }

    if ( areRequired != NULL ) {
      if ( mesh->edge[i].tag & MG_REQ )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }

  return 1;
}

double MMG2D_Get_triangleQuality(MMG5_pMesh mesh,MMG5_pSol met, MMG5_int k) {
  double qual = 0.;
  MMG5_pTria pt;

  if ( k < 1 || k > mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: unable to access to triangle %" MMG5_PRId ".\n",
            __func__,k);
    fprintf(stderr,"     Tria numbering goes from 1 to %" MMG5_PRId "\n",mesh->nt);
    return 0.;
  }
  pt = &mesh->tria[k];
  assert ( MG_EOK(pt) );

  if ( (!met) || (!met->m) || met->size==1 ) {
    /* iso quality */
    qual =  MMG2D_ALPHAD * MMG2D_caltri_iso(mesh,NULL,pt);
  }
  else {
    qual = MMG2D_ALPHAD * MMG2D_caltri_ani(mesh,met,pt);
  }

  return qual;
}

int MMG2D_Set_scalarSol(MMG5_pSol met, double s, MMG5_int pos) {

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG2D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }

  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",
            __func__);
    fprintf(stderr,"    max number of solutions: %" MMG5_PRId "\n",met->npmax);
    return 0;
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution"
            " at position %" MMG5_PRId ".",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %" MMG5_PRId "\n",met->np);
    fprintf(stderr,"  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return 0;
  }

  met->m[pos] = s;
  return 1;
}

int  MMG2D_Get_scalarSol(MMG5_pSol met, double* s)
{
  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter"
              " of points.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG2D_Get_scalarSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",
            __func__);
    fprintf(stderr,"     The number of call of MMG2D_Get_scalarSol function");
    fprintf(stderr," can not exceed the number of points: %" MMG5_PRId "\n ",met->np);
    return 0;
  }

  *s  = met->m[met->npi];

  return 1;
}

int  MMG2D_Set_scalarSols(MMG5_pSol met, double *s) {
  MMG5_int k;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: you must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG2D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }

  for ( k=0; k<met->np; ++k )
    met->m[k+1] = s[k];

  return 1;
}

int  MMG2D_Get_scalarSols(MMG5_pSol met, double* s) {
  MMG5_int k;

  for ( k=0; k<met->np; ++k )
    s[k]  = met->m[k+1];

  return 1;
}

int MMG2D_Set_vectorSol(MMG5_pSol met, double vx,double vy, MMG5_int pos) {
  MMG5_int isol;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: you must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG2D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }
  if ( pos < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",
            __func__);
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return 0;
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",
            __func__);
    fprintf(stderr,"    max number of solutions: %" MMG5_PRId "\n",met->npmax);
    return 0;
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution"
            " at position %" MMG5_PRId ".",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %" MMG5_PRId "\n",met->np);
    fprintf(stderr,"\n  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return 0;
  }

  isol = (pos-1) * met->size + 1;

  met->m[isol]   = vx;
  met->m[isol+1] = vy;

  return 1;
}


int MMG2D_Get_vectorSol(MMG5_pSol met, double* vx, double* vy) {

  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG2D_Get_vectorSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",__func__);
    fprintf(stderr,"     The number of call of MMG2D_Get_vectorSol function");
    fprintf(stderr," can not exceed the number of points: %" MMG5_PRId "\n ",met->np);
    return 0;
  }

  *vx = met->m[met->size*(met->npi-1)+1];
  *vy = met->m[met->size*(met->npi-1)+2];

  return 1;
}

int MMG2D_Set_vectorSols(MMG5_pSol met, double *sols) {
  double *m;
  MMG5_int k,j;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: you must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG2D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }

  for ( k=0; k<met->np; ++k ) {
    j = 2*k;
    m = &met->m[j];
    m[1] = sols[j];
    m[2] = sols[j+1];
  }

  return 1;
}

int MMG2D_Get_vectorSols(MMG5_pSol met, double* sols) {
  double *m;
  MMG5_int k, j;

  for ( k=0; k<met->np; ++k ) {
    j = 2*k;
    m = &met->m[j];

    sols[j]   = m[1];
    sols[j+1] = m[2];
  }

  return 1;
}


int MMG2D_Set_tensorSol(MMG5_pSol met, double m11, double m12, double m22,
                        MMG5_int pos) {
  MMG5_int isol;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: you must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG2D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }

  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",
            __func__);
    fprintf(stderr,"    max number of solutions: %" MMG5_PRId "\n",met->npmax);
    return 0;
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution "
            "at position %" MMG5_PRId ".",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %" MMG5_PRId "\n",met->np);
    fprintf(stderr,"  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return 0;
  }
  isol = pos * met->size;
  met->m[isol    ] = m11;
  met->m[isol + 1] = m12;
  met->m[isol + 2] = m22;
  return 1;
}

int MMG2D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12,double *m22)
{
  int ddebug = 0;
  MMG5_int isol;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG2D_Get_tensorSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",__func__);
    fprintf(stderr,"     The number of call of MMG2D_Get_tensorSol function");
    fprintf(stderr," can not exceed the number of points: %" MMG5_PRId "\n ",met->np);
    return 0;
  }

  isol = met->size*(met->npi);
  *m11 = met->m[isol  ];
  *m12 = met->m[isol+1];
  *m22 = met->m[isol+2];

  return 1;
}

int MMG2D_Set_tensorSols(MMG5_pSol met, double *sols) {
  double *m;
  MMG5_int k,j;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number"
            " of solution with the",__func__);
    fprintf(stderr," MMG2D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j+3];

    m[0] = sols[j];
    m[1] = sols[j+1];
    m[2] = sols[j+2];
  }
  return 1;
}

int MMG2D_Get_tensorSols(MMG5_pSol met, double *sols) {
  double *m;
  MMG5_int k,j;

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j+3];

    sols[j]   = m[0];
    sols[j+1] = m[1];
    sols[j+2] = m[2];
  }

  return 1;
}

int  MMG2D_Set_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double *s) {
  MMG5_pSol psl;

  /* Warning: users give indices from 1 to nsols */
  psl = sol + (i-1);

  switch ( psl->type ) {
  case MMG5_Scalar:
    return MMG2D_Set_scalarSols(psl,s);
    break;

  case MMG5_Vector:
    MMG2D_Set_vectorSols(psl,s);
    break;

  case MMG5_Tensor:
    MMG2D_Set_tensorSols(psl,s);
    break;

  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected type of solution: %s.\n",
            __func__,MMG5_Get_typeName(psl->type));
    return 0;
  }

  return 1;
}

int  MMG2D_Get_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double *s) {
  MMG5_pSol psl;

  /* Warning: users give indices from 1 to nsols */
  psl = sol + (i-1);

  switch ( psl->type ) {
  case MMG5_Scalar:
    return MMG2D_Get_scalarSols(psl,s);
    break;

  case MMG5_Vector:
    MMG2D_Get_vectorSols(psl,s);
    break;

  case MMG5_Tensor:
    MMG2D_Get_tensorSols(psl,s);
    break;

  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected type of solution: %s\n",
            __func__,MMG5_Get_typeName(psl->type));
    return 0;
  }

  return 1;
}

int  MMG2D_Set_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos) {
  MMG5_pSol psl;

  /* Warning: users give indices from 1 to nsols */
  psl = sol + (i-1);

  switch ( psl->type ) {
  case MMG5_Scalar:
    return MMG2D_Set_scalarSol(psl,s[0],pos);
    break;

  case MMG5_Vector:
    MMG2D_Set_vectorSol(psl,s[0],s[1],pos);
    break;

  case MMG5_Tensor:
    MMG2D_Set_tensorSol(psl,s[0],s[1],s[2],pos);
    break;

  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected type of solution: %s.\n",
            __func__,MMG5_Get_typeName(psl->type));
    return 0;
  }
  return 1;
}

int  MMG2D_Get_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double *s,MMG5_int pos) {
  MMG5_pSol psl;

  /* Warning: users give indices from 1 to nsols */
  psl = sol + (i-1);

  psl->npi = pos-1;

  switch ( psl->type ) {
  case MMG5_Scalar:
    return MMG2D_Get_scalarSol(psl,&s[0]);
    break;

  case MMG5_Vector:
    MMG2D_Get_vectorSol(psl,&s[0],&s[1]);
    break;

  case MMG5_Tensor:
    MMG2D_Get_tensorSol(psl,&s[0],&s[1],&s[2]);
    break;

  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected type of solution: %s\n",
            __func__,MMG5_Get_typeName(psl->type));
    return 0;
  }

  return 1;
}

int MMG2D_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( (mesh->npi != mesh->np) || (mesh->nti != mesh->nt) ) {
    fprintf(stderr,"\n  ## Error: %s: if you don't use the MMG2D_loadMesh function,",
            __func__);
    fprintf(stderr," you must call the MMG2D_Set_meshSize function to have a");
    fprintf(stderr," valid mesh.\n");
    fprintf(stderr," Missing datas.\n");
    return 0;
  }

  if ( met->npi != met->np ) {
    fprintf(stderr,"\n  ## Error: %s: if you don't use the MMG2D_loadMet function,",
            __func__);
    fprintf(stderr," you must call the MMG2D_Set_solSize function to have a");
    fprintf(stderr," valid solution.\n");
    fprintf(stderr," Missing datas.\n");
    return 0;
  }

  /*  Check mesh data */
  if ( mesh->info.ddebug ) {
    if ( (!mesh->np) || (!mesh->point) ||
         (!mesh->nt)  ) {
      fprintf(stderr,"  ** MISSING DATA.\n");
      fprintf(stderr," Check that your mesh contains points.\n");
      fprintf(stderr," Exit program.\n");
      return 0;
    }
  }

  if ( mesh->dim != 2 ) {
    fprintf(stderr,"  ** 2 DIMENSIONAL MESH NEEDED. Exit program.\n");
    return 0;
  }
  if ( met->dim != 2 ) {
    fprintf(stderr,"  ** WRONG DIMENSION FOR METRIC. Exit program.\n");
    return 0;
  }
  if ( !mesh->ver )  mesh->ver = 2;
  if ( !met ->ver )  met ->ver = 2;

  return 1;
}

int MMG2D_Free_allSols(MMG5_pMesh mesh,MMG5_pSol *sol) {

  return MMG5_Free_allSols(mesh,sol);
}


int MMG2D_Free_all(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = MMG2D_Free_all_var(argptr);

  va_end(argptr);

  return ier;
}

int MMG2D_Free_structures(const int starter,...)
{
  int ier;

  va_list argptr;

  va_start(argptr, starter);

  ier = MMG2D_Free_structures_var(argptr);

  va_end(argptr);

  return ier;
}

int MMG2D_Free_names(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = MMG2D_Free_names_var(argptr);

  va_end(argptr);

  return ier;
}
