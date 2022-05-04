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
 * \file mmg3d/analys_3d.c
 * \brief Mesh analysis.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"

/**
 * \param mesh pointer towarad the mesh structure.
 *
 * Set all boundary triangles to required and add a tag to detect that they are
 * not realy required.
 *
 */
void MMG3D_set_reqBoundaries(MMG5_pMesh mesh) {
  MMG5_pTria     ptt;
  int            k;

  /* The MG_REQ+MG_NOSURF tag mark the boundary edges that we dont want to touch
   * but that are not really required (-nosurf option) */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];

    if ( mesh->info.nosurf  && (!(ptt->tag[0] & MG_REQ)) ) {
      ptt->tag[0] |= MG_REQ;
      ptt->tag[0] |= MG_NOSURF;
    }

    if ( ptt->tag[0] & MG_PARBDY ) {
      ptt->tag[0] |= MG_NOSURF;
      ptt->tag[0] |= MG_REQ;
    }

    if ( mesh->info.nosurf && (!(ptt->tag[1] & MG_REQ)) ) {
      ptt->tag[1] |= MG_REQ;
      ptt->tag[1] |= MG_NOSURF;
    }

    if ( ptt->tag[1] & MG_PARBDY ) {
      ptt->tag[1] |= MG_NOSURF;
      ptt->tag[1] |= MG_REQ;
    }

    if ( mesh->info.nosurf && (!(ptt->tag[2] & MG_REQ)) ) {
      ptt->tag[2] |= MG_REQ;
      ptt->tag[2] |= MG_NOSURF;
    }

    if ( ptt->tag[2] & MG_PARBDY ) {
      ptt->tag[2] |= MG_NOSURF;
      ptt->tag[2] |= MG_REQ;
    }
  }

  return;
}


/**
 * \param mesh pointer towarad the mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * topology: set tria adjacency, detect Moebius, flip faces, count connected
 * comp.
 *
 * \remark: as all triangles are mesh boundaries, we do not need to mark their
 * adges as MG_BDY so the MG_BDY tag may be used inside geometrical triangles
 * (external non-parallel, or internal parallel) to tag edges on the
 * intersection with purely parallel (non-geometrical) triangles.
 * The MG_PARBDYBDY tag is also added, as it does not have a supporting triangle
 * to inherit this tag from.
 *
 * \remark REQ, NOSURF, etc... tags are added only inside xtetra.
 *
 */
int MMG5_setadj(MMG5_pMesh mesh){
  MMG5_pTria   pt,pt1;
  int          *adja,*adjb,adji1,adji2,*pile,iad,ipil,ip1,ip2,gen;
  int          k,kk,iel,jel,nvf,nf,nr,nt,nre,nreq,ncc,ned,ref;
  int16_t      tag;
  int8_t       i,ii,i1,i2,ii1,ii2,voy;

  nvf = nf = ncc = ned = 0;

  MMG5_SAFE_MALLOC(pile,mesh->nt+1,int,return 0);

  pile[1] = 1;
  ipil    = 1;

  while ( ipil > 0 ) {
    ncc++;

    do {
      k  = pile[ipil--];
      pt = &mesh->tria[k];
      pt->flag = ncc;
      if ( !MG_EOK(pt) )  continue;

      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        if( ((pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY)) ||
            (pt->tag[i] & MG_BDY) ) continue;
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];

        if ( !mesh->point[ip1].tmp )  mesh->point[ip1].tmp = ++nvf;
        if ( !mesh->point[ip2].tmp )  mesh->point[ip2].tmp = ++nvf;

        if ( MG_EDG(pt->tag[i]) || pt->tag[i] & MG_REQ ) {
          tag = mesh->point[ip1].tag;
          mesh->point[ip1].tag |= pt->tag[i];
          // Remove the MG_NOSURF tag if the vertex is really required.
          if ( (tag & MG_REQ) && !(tag & MG_NOSURF) ) {
            mesh->point[ip1].tag &= ~MG_NOSURF;
          }
          tag = mesh->point[ip2].tag;
          mesh->point[ip2].tag |= pt->tag[i];
          // Remove the MG_NOSURF tag if the vertex is really required.
          if ( (tag & MG_REQ) && !(tag & MG_NOSURF) ) {
            mesh->point[ip2].tag &= ~MG_NOSURF;
          }
        }

        /* open boundary */
        tag = MG_GEO;
        if ( mesh->info.opnbdy ) tag += MG_OPNBDY;
        if ( !adja[i] ) {
          tag += MG_NOM;
          pt->tag[i] |= tag;
          mesh->point[ip1].tag |= tag;
          mesh->point[ip2].tag |= tag;
          ned++;
          continue;
        }

        kk = adja[i] / 3;
        ii = adja[i] % 3;
        if ( kk > k )  ned++;

        /* store adjacent */
        pt1 = &mesh->tria[kk];

        /* non manifold edge */
        if ( pt->tag[i] & MG_NOM ) {
          mesh->point[ip1].tag |= MG_NOM;
          mesh->point[ip2].tag |= MG_NOM;
          continue;
        }

        if ( abs(pt1->ref) != abs(pt->ref) ) {
          pt->tag[i]   |= MG_REF;
          pt1->tag[ii] |= MG_REF;
          mesh->point[ip1].tag |= MG_REF;
          mesh->point[ip2].tag |= MG_REF;
        }

        /* store adjacent */
        if ( !pt1->flag ) {
          pt1->flag    = ncc;
          pile[++ipil] = kk;
        }

        /* check orientation */
        ii1 = MMG5_inxt2[ii];
        ii2 = MMG5_iprv2[ii];
        if ( pt1->v[ii1] == ip1 ) {
          /* Moebius strip */
          if ( pt1->base < 0 ) {
            fprintf(stderr,"\n  ## Error: %s: Triangle orientation problem (1):"
                    " Moebius strip?\n",__func__);
            MMG5_SAFE_FREE(pile);
            return 0;
          }
          /* flip orientation */
          else {
            pt1->base   = -pt1->base;
            pt1->v[ii1] = ip2;
            pt1->v[ii2] = ip1;

            /* update adj */
            iad   = 3*(kk-1)+1;
            adjb  = &mesh->adjt[iad];
            adji1 = mesh->adjt[iad+ii1];
            adji2 = mesh->adjt[iad+ii2];
            adjb[ii1] = adji2;
            adjb[ii2] = adji1;

            /* modif tag + ref */
            tag = pt1->tag[ii1];
            pt1->tag[ii1] = pt1->tag[ii2];
            pt1->tag[ii2] = tag;
            ref = pt1->edg[ii1];
            pt1->edg[ii1] = pt1->edg[ii2];
            pt1->edg[ii2] = ref;

            /* modif voyeurs */
            if ( adjb[ii1] ) {
              iel = adjb[ii1] / 3;
              voy = adjb[ii1] % 3;
              mesh->adjt[3*(iel-1)+1+voy] = 3*kk + ii1;
            }
            if ( adjb[ii2] ) {
              iel = adjb[ii2] / 3;
              voy = adjb[ii2] % 3;
              mesh->adjt[3*(iel-1)+1+voy] = 3*kk + ii2;
            }
            nf++;
          }
        }
        else {
          /* Mark triangles that have a consistent orientation with their
           * neighbours */
          pt1->base =  -pt1->base;
        }
      }
    }
    while ( ipil > 0 );

    /* find next unmarked triangle */
    ipil = 0;
    for (kk=1; kk<=mesh->nt; kk++) {
      pt = &mesh->tria[kk];
      if ( MG_EOK(pt) && (pt->flag == 0) ) {
        ipil = 1;
        pile[ipil] = kk;
        pt->flag   = ncc+1;
        break;
      }
    }
  }

  /* bilan */
  nr = nre = nreq = nt = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    nt++;
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( ( !MG_EDG(pt->tag[i]) ) && ( !(pt->tag[i] & MG_REQ) ) )  continue;

      jel  = adja[i] / 3;
      if ( !jel || jel > k ) {
        if ( pt->tag[i] & MG_GEO )  nr++;
        if ( pt->tag[i] & MG_REF )  nre++;
        if ( pt->tag[i] & MG_REQ )  nreq++;
      }
    }
  }

  if ( mesh->info.ddebug ) {
    fprintf(stdout,"  a- ridges: %d found.\n",nr);
    fprintf(stdout,"  a- requir: %d found.\n",nreq);
    fprintf(stdout,"  a- connex: %d connected component(s)\n",ncc);
    fprintf(stdout,"  a- orient: %d flipped\n",nf);
  }
  else if ( abs(mesh->info.imprim) > 3 ) {
    gen = (2 - nvf + ned - nt) / 2;
    fprintf(stdout,"     Connected component: %d,  genus: %d,   reoriented: %d\n",ncc,gen,nf);
    fprintf(stdout,"     Edges: %d,  tagged: %d,  ridges: %d, required: %d, refs: %d\n",
            ned,nr+nre+nreq,nr,nreq,nre);
  }

  MMG5_SAFE_FREE(pile);
  return 1;
}

/** check for ridges: dihedral angle */
int MMG5_setdhd(MMG5_pMesh mesh) {
  MMG5_pTria    pt,pt1;
  double        n1[3],n2[3],dhd;
  int          *adja,k,kk,ne,nr;
  int8_t        i,ii,i1,i2;

  ne = nr = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    /* triangle normal */
    MMG5_nortri(mesh,pt,n1);
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( ((pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY)) ||
           (pt->tag[i] & MG_BDY) ) continue;

      kk  = adja[i] / 3;
      ii  = adja[i] % 3;
      if ( !kk ) {
        pt->tag[i] |= MG_GEO;
        i1 = MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];
        mesh->point[pt->v[i1]].tag |= MG_GEO;
        mesh->point[pt->v[i2]].tag |= MG_GEO;
        nr++;
      }
      else if ( k < kk ) {
        pt1 = &mesh->tria[kk];
        /* reference curve */
        if ( pt1->ref != pt->ref ) {
          pt->tag[i]   |= MG_REF;
          pt1->tag[ii] |= MG_REF;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].tag |= MG_REF;
          mesh->point[pt->v[i2]].tag |= MG_REF;
          ne++;
        }
        /* check angle w. neighbor */
        MMG5_nortri(mesh,pt1,n2);
        dhd = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
        if ( dhd <= mesh->info.dhd ) {
          pt->tag[i]   |= MG_GEO;
          pt1->tag[ii] |= MG_GEO;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].tag |= MG_GEO;
          mesh->point[pt->v[i2]].tag |= MG_GEO;
          nr++;
        }
      }
    }
  }
  if ( abs(mesh->info.imprim) > 3 && nr > 0 )
    fprintf(stdout,"     %d ridges, %d edges updated\n",nr,ne);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1.
 *
 * check subdomains connected by a vertex and mark these vertex as CRN and REQ.
 *
 */
int MMG5_chkVertexConnectedDomains(MMG5_pMesh mesh){
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   ppt;
  int           k,lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilists,ilistv,i0,ier;
  int8_t        i,j;
  static int8_t mmgWarn = 0;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->s = 0;
    ppt->flag = mesh->mark;
  }
  ++mesh->mark;

  /*count the number of tet around a point*/
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )   continue;

    for (i=0; i<4; i++) {
      mesh->point[pt->v[i]].s++;
    }
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )   continue;

    /* point j on face i */
    for (i=0; i<4; i++) {
      for (j=0; j<3; j++) {
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
        }
        else  pxt = 0;

        i0  = MMG5_idir[i][j];
        ppt = &mesh->point[pt->v[i0]];
        if ( !(ppt->tag & MG_BDY) ) continue;
        if ( ppt->flag == mesh->mark ) continue;

        /* Catch a boundary point by a boundary face */
        if ( (!pt->xt) || !(MG_BDY & pxt->ftag[i]) )  continue;
        if( ppt->tag & MG_NOM ){
          if ( mesh->adja[4*(k-1)+1+i] ) continue;
          ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,1);
        } else {
          ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
        }
        if ( ier != 1 && !mmgWarn ) {
          mmgWarn = 1;
          printf("  ## Warning: %s: unable to check that we don't have"
                 " non-connected domains.\n",__func__);
        }

        if(ilistv != ppt->s) {
          if(!(ppt->tag & MG_REQ) ) {
            ppt->tag |= MG_REQ;
            ppt->tag |= MG_CRN;
          }
        }
        ppt->flag = mesh->mark;
      }
    }
  }
  return 1;
}

/** check for singularities */
int MMG5_singul(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p1,p2;
  double         ux,uy,uz,vx,vy,vz,dd;
  int            list[MMG3D_LMAX+2],listref[MMG3D_LMAX+2],k,nc,xp,nr,ns,nre;
  int8_t         i;

  nre = nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) || ( ppt->tag & MG_NOM ) ||
          ( ppt->tag & MG_PARBDY ) ) continue;
      else if ( MG_EDG(ppt->tag) ) {
        /* Store the number of ridges passing through the point (xp) and the
         * number of ref edges (nr) */
        ns = MMG5_bouler(mesh,mesh->adjt,k,i,list,listref,&xp,&nr,MMG3D_LMAX);

        if ( !ns )  continue;
        if ( (xp+nr) > 2 ) {
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        else if ( (xp == 1) && (nr == 1) ) {
          ppt->tag |= MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
        }
        else if ( xp == 1 && !nr ){
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        else if ( nr == 1 && !xp ){
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        /* check ridge angle */
        else {
          p1 = &mesh->point[list[1]];
          p2 = &mesh->point[list[2]];
          ux = p1->c[0] - ppt->c[0];
          uy = p1->c[1] - ppt->c[1];
          uz = p1->c[2] - ppt->c[2];
          vx = p2->c[0] - ppt->c[0];
          vy = p2->c[1] - ppt->c[1];
          vz = p2->c[2] - ppt->c[2];
          dd = (ux*ux + uy*uy + uz*uz) * (vx*vx + vy*vy + vz*vz);
          if ( fabs(dd) > MMG5_EPSD ) {
            dd = (ux*vx + uy*vy + uz*vz) / sqrt(dd);
            if ( dd > -mesh->info.dhd ) {
              ppt->tag |= MG_CRN;
              nc++;
            }
          }
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 3 && nre > 0 )
    fprintf(stdout,"     %d corners, %d singular points detected\n",nc,nre);
  return 1;
}

/** compute normals at C1 vertices, for C0: tangents */
int MMG5_norver(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt;
  MMG5_xPoint    *pxp;
  double         n[3],dd;
  int            *adja,k,kk,ng,nn,nt,nf,nnr;
  int8_t         i,ii,i1;

  /* recomputation of normals only if mesh->xpoint has been freed */
  if ( mesh->xpoint ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: %s: no research of boundary points"
              " and normals of mesh. mesh->xpoint must be freed to enforce"
              " analysis.\n",__func__);
    }
    return 1;
  }

  /* identify boundary points */
  ++mesh->base;
  mesh->xp = 0;
  nnr      = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->flag == mesh->base )  continue;
      else {
        ++mesh->xp;
        ppt->flag = mesh->base;
        if ( mesh->nc1 ) {
          if ( ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2] > 0 ) {
            if ( ppt->tag & MG_PARBDY || ppt->tag & MG_CRN || ppt->tag & MG_NOM || MG_EDG(ppt->tag) ) {
              ++nnr;
              continue;
            }
            ppt->xp = -1;
          }
        }
      }
    }
  }

  /* memory to store normals for boundary points */
  mesh->xpmax  = MG_MAX( (long long)(1.5*mesh->xp),mesh->npmax);

  MMG5_ADD_MEM(mesh,(mesh->xpmax+1)*sizeof(MMG5_xPoint),"boundary points",return 0);
  MMG5_SAFE_CALLOC(mesh->xpoint,mesh->xpmax+1,MMG5_xPoint,return 0);

  /* compute normals + tangents */
  nn = ng = nt = nf = 0;
  mesh->xp = 0;
  ++mesh->base;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->tag & MG_PARBDY || ppt->tag & MG_CRN || ppt->tag & MG_NOM || ppt->flag == mesh->base )  continue;

      /* C1 point */
      if ( !MG_EDG(ppt->tag) ) {

        if ( (!mesh->nc1) ||
             ppt->n[0]*ppt->n[0]+ppt->n[1]*ppt->n[1]+ppt->n[2]*ppt->n[2]<=MMG5_EPSD2 ) {
          if ( !MMG5_boulen(mesh,mesh->adjt,k,i,ppt->n) ) {
            ++nf;
            continue;
          }
          else ++nn;
        }

        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                             "larger xpoint table",
                             mesh->xp--;return 0;);
        }
        ppt->xp = mesh->xp;
        pxp = &mesh->xpoint[ppt->xp];
        memcpy(pxp->n1,ppt->n,3*sizeof(double));
        ppt->n[0] = ppt->n[1] = ppt->n[2] = 0.;
        ppt->flag = mesh->base;

      }

      /* along ridge-curve */
      i1  = MMG5_inxt2[i];
      if ( !MG_EDG(pt->tag[i1]) )  continue;
      else if ( !MMG5_boulen(mesh,mesh->adjt,k,i,n) ) {
        ++nf;
        continue;
      }
      ++mesh->xp;
      if(mesh->xp > mesh->xpmax){
        MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                           "larger xpoint table",
                           mesh->xp--;return 0;);
      }
      ppt->xp = mesh->xp;
      pxp = &mesh->xpoint[ppt->xp];
      memcpy(pxp->n1,n,3*sizeof(double));

      if ( pt->tag[i1] & MG_GEO && adja[i1] > 0 ) {
        kk = adja[i1] / 3;
        ii = adja[i1] % 3;
        ii = MMG5_inxt2[ii];
        if ( !MMG5_boulen(mesh,mesh->adjt,kk,ii,n) ) {
          ++nf;
          continue;
        }
        memcpy(pxp->n2,n,3*sizeof(double));

        /* compute tangent as intersection of n1 + n2 */
        ppt->n[0] = pxp->n1[1]*pxp->n2[2] - pxp->n1[2]*pxp->n2[1];
        ppt->n[1] = pxp->n1[2]*pxp->n2[0] - pxp->n1[0]*pxp->n2[2];
        ppt->n[2] = pxp->n1[0]*pxp->n2[1] - pxp->n1[1]*pxp->n2[0];
        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
        if ( dd > MMG5_EPSD2 ) {
          dd = 1.0 / sqrt(dd);
          ppt->n[0] *= dd;
          ppt->n[1] *= dd;
          ppt->n[2] *= dd;
        }
        ppt->flag = mesh->base;
        ++nt;
        continue;
      }

      /* compute tgte */
      ppt->flag = mesh->base;
      ++nt;
      if ( !MMG5_boulec(mesh,mesh->adjt,k,i,ppt->n) ) {
        ++nf;
        continue;
      }
      dd = pxp->n1[0]*ppt->n[0] + pxp->n1[1]*ppt->n[1] + pxp->n1[2]*ppt->n[2];
      ppt->n[0] -= dd*pxp->n1[0];
      ppt->n[1] -= dd*pxp->n1[1];
      ppt->n[2] -= dd*pxp->n1[2];
      dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        ppt->n[0] *= dd;
        ppt->n[1] *= dd;
        ppt->n[2] *= dd;
      }
    }
  }
  mesh->nc1 = 0;

  if ( abs(mesh->info.imprim) > 3 && nn+nt > 0 ) {
    if ( nnr )
      fprintf(stdout,"     %d input normals ignored\n",nnr);
    fprintf(stdout,"     %d normals,  %d tangents updated  (%d failed)\n",nn,nt,nf);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 0 if fail, 1 otherwise
 *
 * Define continuous geometric support at non manifold vertices, using volume
 * information.
 *
 */
int MMG3D_nmgeom(MMG5_pMesh mesh){
  MMG5_pTetra     pt;
  MMG5_pPoint     p0;
  MMG5_pxPoint    pxp;
  int             k,base;
  int             *adja;
  double          n[3],t[3];
  int8_t          i,j,ip,ier;

  base = ++mesh->base;
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if( !MG_EOK(pt) ) continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      if ( adja[i] ) continue;
      for (j=0; j<3; j++) {
        ip = MMG5_idir[i][j];
        p0 = &mesh->point[pt->v[ip]];
        if ( p0->flag == base )  continue;
        else if ( !(p0->tag & MG_NOM) || (p0->tag & MG_PARBDY) )  continue;

        p0->flag = base;
        ier = MMG5_boulenm(mesh,k,ip,i,n,t);

        if ( ier < 0 )
          return 0;
        else if ( !ier ) {
          p0->tag |= MG_REQ;
          p0->tag &= ~MG_NOSURF;
        }
        else {
          if ( !p0->xp ) {
            ++mesh->xp;
            if(mesh->xp > mesh->xpmax){
              MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                                 "larger xpoint table",
                                 mesh->xp--;
                                 fprintf(stderr,"  Exit program.\n");return 0;);
            }
            p0->xp = mesh->xp;
          }
          pxp = &mesh->xpoint[p0->xp];
          memcpy(pxp->n1,n,3*sizeof(double));
          memcpy(p0->n,t,3*sizeof(double));
        }
      }
    }
  }
  /* Deal with the non-manifold points that do not belong to a surface
   * tetra (a tetra that has a face without adjacent)*/
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if( !MG_EOK(pt) ) continue;
    
    for (i=0; i<4; i++) {
      p0 = &mesh->point[pt->v[i]];
      if ( p0->tag & MG_REQ || !(p0->tag & MG_NOM) ||
           p0->xp || (p0->tag & MG_PARBDY) ) continue;
      ier = MMG5_boulenmInt(mesh,k,i,t);
      if ( ier ) {
        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                            "larger xpoint table",
                            mesh->xp--;
                            fprintf(stderr,"  Exit program.\n");return 0;);
        }
        p0->xp = mesh->xp;
        pxp = &mesh->xpoint[p0->xp];
        memcpy(p0->n,t,3*sizeof(double));
        pxp->nnor = 1;
      }
      else {
        p0->tag |= MG_REQ;
        p0->tag &= ~MG_NOSURF;
      }
    }
  }

  /*for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !(p0->tag & MG_NOM) || p0->xp ) continue;
    p0->tag |= MG_REQ;
    p0->tag &= ~MG_NOSURF;
  }*/
  
  return 1;
}

/**
 * \param mesh pointer toward mesh structure
 * \return 1 if success, 0 if fail
 *
 * Assign surface references demending on geometric features: portion of
 * surfaces shared by non-manifold or ridge edges will have different
 * references.
 */
int MMG3D_Set_faceColors(MMG5_pMesh mesh) {
  MMG5_pTria   pt,pt1;
  int          *adja,*pile,ipil;
  int          k,kk,ncc=0;
  int8_t       i;

  MMG5_SAFE_MALLOC(pile,mesh->nt+1,int,return 0);

  /* Reset flags */
  for (kk=1; kk<=mesh->nt; kk++) {
    mesh->tria[kk].flag = 0;
  }

  pile[1] = 1;
  ipil    = 1;
  while ( ipil > 0 ) {
    ncc++;

    /* Pile up triangles that can be seen from the initial one without crossing
     * a ridge or non manifold edge. Change surface reference of stacked
     * triangles */
    do {
      k  = pile[ipil--];
      pt = &mesh->tria[k];
      pt->flag = ncc;

      /* Assign new reference to portion of surface */
      pt->ref  = ncc;

      if ( !MG_EOK(pt) )  continue;

      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        /* Don't cross ridges or non-manifold edges */
        if ( (pt->tag[i] & MG_GEO) || (pt->tag[i] & MG_NOM) ) continue;

        /* store adjacent */
        kk = adja[i] / 3;
        pt1 = &mesh->tria[kk];
        if ( !pt1->flag ) {
          pt1->flag    = ncc;
          pile[++ipil] = kk;
        }
      }
    }
    while ( ipil > 0 );

    /* find next unmarked triangle */
    ipil = 0;
    for (kk=1; kk<=mesh->nt; kk++) {
      pt = &mesh->tria[kk];
      if ( MG_EOK(pt) && (pt->flag == 0) ) {
        ipil = 1;
        pile[ipil] = kk;
        pt->flag   = ncc+1;
        break;
      }
    }
  }

  /* bilan */
  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
    fprintf(stdout,"     %d surface area detected: assign new references\n",ncc);
  }

  MMG5_SAFE_FREE(pile);
  return 1;
}



/** preprocessing stage: mesh analysis */
int MMG3D_analys(MMG5_pMesh mesh) {
  MMG5_Hash hash;
  int       ier;

  /**--- stage 1: data structures for surface */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** SURFACE ANALYSIS\n");

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
    return 0;
  }

  if ( mesh->info.iso && mesh->info.opnbdy ) {
    ier = MMG3D_update_xtetra ( mesh );
    if ( !ier ) {
      fprintf(stderr,"\n  ## Problem when updating the xtetra data after ls discretization."
              " Exit program.\n");
      return 0;
    }
  }

  /* create prism adjacency */
  if ( !MMG3D_hashPrism(mesh) ) {
    fprintf(stderr,"\n  ## Prism hashing problem. Exit program.\n");
    return 0;
  }

  /* compatibility triangle orientation w/r tetras */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }
  
  /* identify surface mesh */
  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }

  MMG5_freeXTets(mesh);
  MMG5_freeXPrisms(mesh);

  /* Set surface triangles to required in nosurf mode or for parallel boundaries */
  MMG3D_set_reqBoundaries(mesh);

  /* create surface adjacency */
  memset ( &hash, 0x0, sizeof(MMG5_Hash));
  if ( !MMG3D_hashTria(mesh,&hash) ) {
    MMG5_DEL_MEM(mesh,hash.item);
    fprintf(stderr,"\n  ## Hashing problem (2). Exit program.\n");
    return 0;
  }

  /* build hash table for geometric edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /**--- stage 2: surface analysis */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  /* identify connexity and flip orientation of faces if needed */
  if ( !MMG5_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* check for ridges */
  if ( mesh->info.dhd > MMG5_ANGLIM && !MMG5_setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* identify singularities */
  if ( !MMG5_singul(mesh) ) {
    fprintf(stderr,"\n  ## MMG5_Singularity problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

  /* define (and regularize) normals */
  if ( !MMG5_norver(mesh) ) {
    fprintf(stderr,"\n  ## Normal problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }
  if ( mesh->info.nreg && !MMG5_regnor(mesh) ) {
    fprintf(stderr,"\n  ## Normal regularization problem. Exit program.\n");
    return 0;
  }

  if ( getenv("MMG_COLOR_SURFAREA") ) {
    if ( !MMG3D_Set_faceColors(mesh) )
      return 0;
  }

  /* set bdry entities to tetra and fill the orientation field */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* set non-manifold edges sharing non-intersecting multidomains as required */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** UPDATING TOPOLOGY AT NON-MANIFOLD POINTS\n");

  if ( !MMG5_setNmTag(mesh,&hash) ) {
    fprintf(stderr,"\n  ## Non-manifold topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* check subdomains connected by a vertex and mark these vertex as corner and
     required */
  MMG5_chkVertexConnectedDomains(mesh);

  /* build hash table for geometric edges */
  if ( !mesh->na && !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /* Update edges tags and references for xtetras */
  if ( !MMG5_bdryUpdate(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* define geometry for non manifold points */
  if ( !MMG3D_nmgeom(mesh) ) return 0;

#ifdef USE_POINTMAP
  /* Initialize source point with input index */
  int ip;
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].src = ip;
#endif

  /* release memory */
  MMG5_DEL_MEM(mesh,mesh->htab.geom);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  if ( mesh->nprism ) MMG5_DEL_MEM(mesh,mesh->adjapr);

  return 1;
}
