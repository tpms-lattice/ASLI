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
 * \file mmg3d/anisomovpt_3d.c
 * \brief Functions to move a point in the mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "inlined_functions_3d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param list pointer toward the volumic ball of the point.
 * \param ilist size of the volumic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.
 *
 * \return 0 if we can't move the point, 1 if we can.
 *
 * Move internal point whose volumic is passed.
 *
 * \remark the metric is not interpolated at the new position.
 * \remark we don't check if we break the hausdorff criterion.
 *
 */
int MMG5_movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met, MMG3D_pPROctree PROctree, int *list,int ilist,
                       int improve) {


  MMG5_pTetra          pt,pt0;
  MMG5_pPoint          p0,p1,p2,p3,ppt0;
  double               vol,totvol,m[6];
  double               calold,calnew,callist[MMG3D_LMAX+2],det;
  int                  k,iel,i0;

  assert ( ilist > 0 );
  if ( ilist <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s:"
            " volumic ball has null or negative size (%d)\n",
            __func__,ilist);
    return 0;
  }

  pt0    = &mesh->tetra[0];
  ppt0   = &mesh->point[0];
  memset(ppt0,0,sizeof(MMG5_Point));

  if ( met->m ) {
    iel = list[0] / 4;
    i0  = list[0] % 4;
    memcpy(&met->m[0],&met->m[met->size*mesh->tetra[iel].v[i0]],met->size*sizeof(double));
  }

  /* Coordinates of optimal point */
  calold = DBL_MAX;
  totvol = 0.0;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    pt = &mesh->tetra[iel];
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];
    vol= MMG5_det4pt(p0->c,p1->c,p2->c,p3->c);

    if ( !MMG5_moymet(mesh,met,pt,m) ) {
      // MMG5_moymet must succeed because we have at least 1 point of the tet
      // that is internal.
      return 0;
    }

    det = m[0] * ( m[3]*m[5] - m[4]*m[4]) - m[1] * ( m[1]*m[5] - m[2]*m[4])
      + m[2] * ( m[1]*m[4] - m[2]*m[3]);
    if ( det < MMG5_EPSD2 ) {
      return 0;
    }

    vol *= sqrt(det);

    totvol += vol;
    /* barycenter */
    ppt0->c[0] += 0.25 * vol*(p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]);
    ppt0->c[1] += 0.25 * vol*(p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]);
    ppt0->c[2] += 0.25 * vol*(p0->c[2] + p1->c[2] + p2->c[2] + p3->c[2]);
    calold = MG_MIN(calold, pt->qual);
  }
  if (totvol < MMG5_EPSD2) {
    return 0;
  }

  totvol = 1.0 / totvol;
  ppt0->c[0] *= totvol;
  ppt0->c[1] *= totvol;
  ppt0->c[2] *= totvol;

  /* Check new position validity */
  calnew = DBL_MAX;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    i0  = list[k] % 4;
    pt  = &mesh->tetra[iel];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    callist[k] = MMG5_orcal(mesh,met,0);
    if (callist[k] < MMG5_NULKAL) {
      return 0;
    }
    calnew = MG_MIN(calnew,callist[k]);
  }
  if (calold < MMG5_EPSOK && calnew <= calold) {
    return 0;
  }
  else if (calnew < MMG5_EPSOK) {
    return 0;
  }
  else if ( improve && calnew < 1.02* calold ) {
    return 0;
  }
  else if ( calnew < 0.3 * calold ) {
    return 0;
  }

  /* update position */
  if ( PROctree )
    MMG3D_movePROctree(mesh, PROctree, pt->v[i0], ppt0->c, p0->c);

  p0 = &mesh->point[pt->v[i0]];
  p0->c[0] = ppt0->c[0];
  p0->c[1] = ppt0->c[1];
  p0->c[2] = ppt0->c[2];
  for (k=0; k<ilist; k++) {
    (&mesh->tetra[list[k]/4])->qual=callist[k];
    (&mesh->tetra[list[k]/4])->mark=mesh->mark;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.

 * \return 0 if we can't move the point, 1 if we can, -1 if we fail.
 *
 * \remark we don't check if we break the hausdorff criterion.
 * \remark the metric is not interpolated at the new position.
 *
 * Move boundary regular point, whose volumic and surfacic balls are passed.
 *
 */
int MMG5_movbdyregpt_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG3D_pPROctree PROctree, int *listv,
                          int ilistv,int *lists,int ilists,
                          int improveSurf, int improveVol) {
  MMG5_pTetra       pt,pt0;
  MMG5_pxTetra      pxt;
  MMG5_pPoint       p0;
  MMG5_Tria         tt;
  MMG5_pxPoint      pxp;
  MMG5_Bezier       b;
  double            *n,r[3][3],lispoi[3*MMG3D_LMAX+1],det2d;
  double            detloc,gv[2],step,lambda[3];
  double            o[3],no[3],*m0,ncur[3],nprev[3],nneighi[3];
  double            calold,calnew,caltmp,callist[MMG3D_LMAX+2];
  int               k,kel,iel,l,ip0,nxp,ier;
  uint8_t           i0,iface,i;
  static int        warn = 0;
  static int8_t     mmgErr0=0;

  step = 0.1;
  if ( ilists < 2 )      return 0;

  k  = listv[0] / 4;
  i0 = listv[0] % 4;
  pt = &mesh->tetra[k];
  ip0 = pt->v[i0];
  p0 = &mesh->point[ip0];
  m0 = &met->m[6*ip0];
  assert( p0->xp && !MG_EDG(p0->tag) );

  n = &(mesh->xpoint[p0->xp].n1[0]);

  /** Step 1 : rotation matrix that sends normal n to the third coordinate vector of R^3 */
  if ( !MMG5_rotmatrix(n,r) ) {
    return 0;
  }

  /** Step 2 : rotation of the oriented surfacic ball with r : lispoi[k] is the common edge
      between faces lists[k-1] and lists[k] */
  if ( !MMG3D_rotate_surfacicBall(mesh,lists,ilists,ip0,r,lispoi) ) {
    return 0;
  }

  /** Step 3 :  Compute gradient towards optimal position = centre of mass of the
      ball, projected to tangent plane */
  gv[0] = 0.0;
  gv[1] = 0.0;

  for (k=0; k<ilists; k++) {
    iel    = lists[k] / 4;
    iface  = lists[k] % 4;
    pt     = &mesh->tetra[iel];
    pxt    = &mesh->xtetra[pt->xt];

    MMG5_tet2tri(mesh,iel,iface,&tt);

    if(!MMG5_bezierCP(mesh,&tt,&b,MG_GET(pxt->ori,iface))){
      if( !mmgErr0 ) {
        mmgErr0 = 1;
        fprintf(stderr,"\n  ## Error: %s: function MMG5_bezierCP return 0.\n",
                __func__);
      }
      return -1;
    }

    /* Compute integral of sqrt(T^J(xi)  M(P(xi)) J(xi)) * P(xi) over the triangle */
    if ( !MMG5_elementWeight(mesh,met,&tt,p0,&b,r,gv) ) {
      if ( !warn ) {
        ++warn;
        fprintf(stderr,"\n  ## Warning: %s:"
                " unable to compute optimal position for at least"
                " 1 point.\n",__func__);
      }
      return 0;
    }
  }

  /* At this point : gv = - gradient of V = direction to follow */
  /** Step 4 : locate new point in the ball, and compute its barycentric coordinates */
  det2d = lispoi[1]*gv[1] - lispoi[2]*gv[0];
  kel = 0;
  if ( det2d >= 0.0 ) {
    for (k=0; k<ilists; k++) {
      detloc = gv[0]*lispoi[3*(k+1)+2] - gv[1]*lispoi[3*(k+1)+1];
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == ilists ) {
      return 0;
    }
  }
  else {
    for (k=ilists-1; k>=0; k--) {
      detloc = lispoi[3*k+1]*gv[1] - lispoi[3*k+2]*gv[0];
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == -1 ) {
      return 0;
    }
  }

  /* Sizing of time step : make sure point does not go out corresponding triangle. */
  det2d = -gv[1]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]) + \
    gv[0]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]);
  if ( fabs(det2d) < MMG5_EPSD2 ) {
    return 0;
  }

  det2d = 1.0 / det2d;
  step *= det2d;

  det2d = lispoi[3*(kel)+1]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]) - \
    lispoi[3*(kel)+2 ]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]);
  step *= det2d;
  step  = fabs(step);
  gv[0] *= step;
  gv[1] *= step;

  /* Computation of the barycentric coordinates of the new point in the corresponding triangle. */
  det2d = lispoi[3*kel+1]*lispoi[3*(kel+1)+2] - lispoi[3*kel+2]*lispoi[3*(kel+1)+1];
  if ( det2d < MMG5_EPSD2 ) {
    return 0;
  }
  det2d = 1.0 / det2d;
  lambda[1] = lispoi[3*(kel+1)+2]*gv[0] - lispoi[3*(kel+1)+1]*gv[1];
  lambda[2] = -lispoi[3*(kel)+2]*gv[0] + lispoi[3*(kel)+1]*gv[1];
  lambda[1]*= (det2d);
  lambda[2]*= (det2d);
  lambda[0] = 1.0 - lambda[1] - lambda[2];

  /** Step 5 : come back to original problem, compute patch in triangle iel and
   * check that geometric approx has not been degraded too much */
  nxp = MMG3D_movbdyregpt_geom(mesh,lists,kel,ip0,n,lambda,o,no);
  if ( nxp < 0 ) {
    return -1;
  }
  else if ( !nxp ) {
    return 0;
  }
  pxp = &mesh->xpoint[nxp];

  /* parallel transport of metric at p0 to new point. */
  if ( !MMG5_paratmet(p0->c,n,m0,o,no,&met->m[0]) ) {
    return 0;
  }

  /* For each surfacic triangle build a virtual displaced triangle for check
   * purposes :
   *      - check the new triangle qualities;
   *      - check normal deviation with the adjacent through the edge facing ip0
   *        and the previous one */
  k           = lists[ilists-1] / 4;
  iface       = lists[ilists-1] % 4;

  MMG5_tet2tri(mesh,k,iface,&tt);
  for( i=0 ; i<3 ; i++ )
    if ( tt.v[i] == ip0 )      break;
  assert ( i<3 );
  if ( i>=3 ) return 0;
  tt.v[i] = 0;

  if ( !MMG5_nortri(mesh, &tt, nprev) ) return 0;

  calold = calnew = DBL_MAX;
  for (l=0; l<ilists; l++) {
    k     = lists[l] / 4;
    iface = lists[l] % 4;

    MMG5_tet2tri(mesh,k,iface,&tt);
    calold = MG_MIN(calold,MMG5_caltri(mesh,met,&tt));

    for( i=0 ; i<3 ; i++ )
      if ( tt.v[i] == ip0 )      break;

    assert ( i<3 );
    if ( i>=3 ) return 0;
    tt.v[i] = 0;

    caltmp = MMG5_caltri(mesh,met,&tt);
    if ( caltmp < MMG5_EPSD2 ) {
      /* We don't check the input triangle qualities, thus we may have a very
       * bad triangle in our mesh */
      return 0;
    }
    calnew = MG_MIN(calnew,caltmp);

    if ( !MMG5_nortri(mesh, &tt, ncur) ) return 0;

    if ( ( !(tt.tag[i] & MG_GEO) ) && ( !(tt.tag[i] & MG_NOM) ) ) {
      /* Check normal deviation between k and the triangle facing ip0 */
      ier = MMG3D_normalAdjaTri(mesh,k,iface, i,nneighi);
      if ( ier <= 0 ) {
        return 0;
      }
      ier =  MMG5_devangle( ncur, nneighi, mesh->info.dhd );
      if ( ier <= 0 ) {
        return 0;
      }
    }

    i = MMG5_iprv2[i];

    if ( ( !(tt.tag[i] & MG_GEO) ) && ( !(tt.tag[i] & MG_NOM) ) ) {
      /* Check normal deviation between k and the previous triangle */
      ier =  MMG5_devangle( ncur, nprev, mesh->info.dhd );
      if ( ier<=0 ) {
        return 0;
      }
    }
    memcpy(nprev, ncur, 3*sizeof(double));

  }
  if ( calold < MMG5_EPSOK && calnew <= calold ) {
    return 0;
  }  else if (calnew < MMG5_EPSOK) {
    return 0;
  }
  else if (improveSurf && calnew < 1.02*calold) {
    return 0;
  }
  else if ( calnew < 0.3*calold ) {
    return 0;
  }
  memset(pxp,0,sizeof(MMG5_xPoint));

  /* Test : check whether all volumes remain positive with new position of the point */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilistv; l++) {
    k    = listv[l] / 4;
    i0 = listv[l] % 4;
    pt = &mesh->tetra[k];
    pt0 = &mesh->tetra[0];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    calold = MG_MIN(calold, pt->qual);
    callist[l]=MMG5_orcal(mesh,met,0);
    if ( callist[l] < MMG5_NULKAL )  {
      return 0;
    }
    calnew = MG_MIN(calnew,callist[l]);
  }

  if ( calold < MMG5_EPSOK && calnew <= calold ) {
    return 0;
  }
  else if (calnew < MMG5_EPSOK) {
    return 0;
  }
  else if (improveVol && calnew < calold) {
    return 0;
  }
  else if ( calnew < 0.3*calold ) {
    return 0;
  }

  /* When all tests have been carried out, update coordinates, normals and metrics*/
  if ( PROctree )
    MMG3D_movePROctree(mesh, PROctree, ip0, o, p0->c);

  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  n[0] = no[0];
  n[1] = no[1];
  n[2] = no[2];

  memcpy(m0,&met->m[0],6*sizeof(double));

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual= callist[l];
    (&mesh->tetra[listv[l]/4])->mark= mesh->mark;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.
 * \param edgTag Type of edge on which we move (MG_REF, MG_NOM or MG_GEO).
 *
 * \return 0 if fail, 1 if success.
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 * Move boundary reference, ridge or non-manifold point, whose volumic and
 * surfacic balls are passed.
 *
 */
static inline
int MMG3D_movbdycurvept_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG3D_pPROctree PROctree, int *listv,
                            int ilistv, int *lists, int ilists,int improve,const int16_t edgTag){
  MMG5_pTetra           pt,pt0;
  MMG5_pPoint           p0,ppt0;
  MMG5_Tria             tt;
  MMG5_pxPoint          pxp;
  double                step,ll1old,ll2old,l1new,l2new;
  double                o[3],no[3],no2[3],to[3], ncur[3],nprev[3],nneighi[3];
  double                calold,calnew,caltmp,callist[MMG3D_LMAX+2];
  int                   l,iel,ip0,ipa,ipb,iptmpa,iptmpb,ip1,ip2,ip,nxp;
  int16_t               tag,ier;
  uint8_t               i,i0,ie,iface,iea,ieb,isrid;

  step      = 0.1;
  ip1 = ip2 = 0;
  pt        = &mesh->tetra[listv[0]/4];
  ip0       = pt->v[listv[0]%4];
  p0        = &mesh->point[ip0];

  /* Compute if the edge is a simple ridge to know if we have to compute a
   * second normal at point */
  isrid     = ((MG_GEO & edgTag) && !(MG_NOM & edgTag));

  assert ( edgTag & p0->tag );

  /* Travel surfacic ball and recover the two ending points of curve : two
     senses must be used */
  iel = lists[0]/4;
  iface = lists[0]%4;
  pt = &mesh->tetra[iel];
  ipa = ipb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[MMG5_idir[iface][i]] != ip0 ) {
      if ( !ipa )
        ipa = pt->v[MMG5_idir[iface][i]];
      else
        ipb = pt->v[MMG5_idir[iface][i]];
    }
  }
  assert(ipa && ipb);

  for (l=1; l<ilists; l++) {
    iel = lists[l]/4;
    iface = lists[l]%4;
    pt = &mesh->tetra[iel];
    iea = ieb = 0;
    for (i=0; i<3; i++) {
      ie = MMG5_iarf[iface][i]; //edge i on face iface
      if ( (pt->v[MMG5_iare[ie][0]] == ip0) || (pt->v[MMG5_iare[ie][1]] == ip0) ) {
        if ( !iea )
          iea = ie;
        else
          ieb = ie;
      }
    }
    if ( pt->v[MMG5_iare[iea][0]] != ip0 )
      iptmpa = pt->v[MMG5_iare[iea][0]];
    else {
      assert(pt->v[MMG5_iare[iea][1]] != ip0);
      iptmpa = pt->v[MMG5_iare[iea][1]];
    }
    if ( pt->v[MMG5_iare[ieb][0]] != ip0 )
      iptmpb = pt->v[MMG5_iare[ieb][0]];
    else {
      assert(pt->v[MMG5_iare[ieb][1]] != ip0);
      iptmpb = pt->v[MMG5_iare[ieb][1]];
    }
    if ( (iptmpa == ipa) || (iptmpa == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[iea];
      else  tag = 0;
      if ( edgTag & tag ) {
        ip1 = iptmpa;
        break;
      }
    }
    if ( (iptmpb == ipa) || (iptmpb == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[ieb];
      else  tag = 0;
      if ( edgTag & tag ) {
        ip1 = iptmpb;
        break;
      }
    }
    ipa = iptmpa;
    ipb = iptmpb;
  }

  /* Now travel surfacic list in the reverse sense so as to get the second curve */
  iel = lists[0]/4;
  iface = lists[0]%4;
  pt = &mesh->tetra[iel];
  ipa = ipb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[MMG5_idir[iface][i]] != ip0 ) {
      if ( !ipa )
        ipa = pt->v[MMG5_idir[iface][i]];
      else
        ipb = pt->v[MMG5_idir[iface][i]];
    }
  }
  assert(ipa && ipb);

  for (l=ilists-1; l>0; l--) {
    iel         = lists[l] / 4;
    iface = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    iea         = ieb = 0;
    for (i=0; i<3; i++) {
      ie = MMG5_iarf[iface][i]; //edge i on face iface
      if ( (pt->v[MMG5_iare[ie][0]] == ip0) || (pt->v[MMG5_iare[ie][1]] == ip0) ) {
        if ( !iea )
          iea = ie;
        else
          ieb = ie;
      }
    }
    if ( pt->v[MMG5_iare[iea][0]] != ip0 )
      iptmpa = pt->v[MMG5_iare[iea][0]];
    else {
      assert(pt->v[MMG5_iare[iea][1]] != ip0);
      iptmpa = pt->v[MMG5_iare[iea][1]];
    }
    if ( pt->v[MMG5_iare[ieb][0]] != ip0 )
      iptmpb = pt->v[MMG5_iare[ieb][0]];
    else {
      assert(pt->v[MMG5_iare[ieb][1]] != ip0);
      iptmpb = pt->v[MMG5_iare[ieb][1]];
    }
    if ( (iptmpa == ipa) || (iptmpa == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[iea];
      else  tag = 0;
      if ( edgTag & tag ) {
        ip2 = iptmpa;
        break;
      }
    }
    if ( (iptmpb == ipa) || (iptmpb == ipb) ) {
      if ( (MG_GEO & edgTag) && (!pt->xt) ) {
        tag = 0;
      }
      else {
        assert(pt->xt);
        tag = mesh->xtetra[pt->xt].tag[ieb];
      }
      if ( edgTag & tag ) {
        ip2 = iptmpb;
        break;
      }
    }
    ipa = iptmpa;
    ipb = iptmpb;
  }
  if ( !(ip1 && ip2 && (ip1 != ip2)) )  return 0;

  /* At this point, we get the point extremities of the curve passing through
     ip0 : ip1, ip2, along with support tets it1,it2, the surface faces
     iface1,iface2, and the associated edges ie1,ie2.*/

  /* Changes needed for choice of time step : see manuscript notes */
  ll1old = MMG5_lenSurfEdg(mesh,met,ip0,ip1,isrid);
  ll2old = MMG5_lenSurfEdg(mesh,met,ip0,ip2,isrid);

  if ( (!ll1old) || (!ll2old) ) return 0;

  if ( ll1old < ll2old ) { //move towards p2
    ip = ip2;
  }
  else {
    ip = ip1;
  }

  /* Compute support of the associated edge, and features of the new position */
  if ( MG_NOM & edgTag ) {
    if ( !(MMG5_BezierNom(mesh,ip0,ip,step,o,no,to)) ) {
      return 0;
    }
  }
  else if ( MG_GEO & edgTag ) {
    if ( !(MMG5_BezierRidge(mesh,ip0,ip,step,o,no,no2,to)) ) {
      return 0;
    }
  }
  else if ( MG_REF & edgTag ) {
    if ( !(MMG5_BezierRef(mesh,ip0,ip,step,o,no,to)) ) {
      return 0;
    }
  }
  else {
    assert ( 0 && "Unexpected edge tag in this function");
    return 0;
  }

  /* Test : make sure that geometric approximation has not been degraded too much */
  ppt0 = &mesh->point[0];
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];
  ppt0->tag  = p0->tag;
  ppt0->ref  = p0->ref;


  nxp = mesh->xp + 1;
  if ( nxp > mesh->xpmax ) {
    MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                       "larger xpoint table",return 0;);
  }
  ppt0->xp = nxp;
  pxp = &mesh->xpoint[nxp];
  memcpy(pxp,&(mesh->xpoint[p0->xp]),sizeof(MMG5_xPoint));

  ppt0->n[0] = to[0];
  ppt0->n[1] = to[1];
  ppt0->n[2] = to[2];

  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

  if ( isrid ) {
    /* Copy the second normal for ridge point */
    pxp->n2[0] = no2[0];
    pxp->n2[1] = no2[1];
    pxp->n2[2] = no2[2];
  }

  if ( (MG_GEO & edgTag) && !(MG_NOM & edgTag) ) {
    /* Interpolation of metric between ip0 and ip2 along ridge */
    if ( !MMG5_intridmet(mesh,met,ip0,ip,step,no,&met->m[0]) ) {
      return 0;
    }
  }
  else {
    /* Interpolation of metric between ip0 and ip2 along non manifold or ref edge */
    if ( !MMG5_paratmet(p0->c,mesh->xpoint[p0->xp].n1,&met->m[6*ip0],o,no,&met->m[0]) ) {
      return 0;
    }
  }

  /* Check whether proposed move is admissible under consideration of distances */
  l1new = MMG5_lenSurfEdg(mesh,met,0,ip1,isrid);
  l2new = MMG5_lenSurfEdg(mesh,met,0,ip2,isrid);

  if ( (!l1new) || (!l2new) ) return 0;

  if ( fabs(l2new -l1new) >= fabs(ll2old -ll1old) )
    return 0;

  /* For each surfacic triangle build a virtual displaced triangle for check
   * purposes :
   *      - check the new triangle qualities;
   *      - check normal deviation with the adjacent through the edge facing ip0
   *        and the previous one */
  iel         = lists[ilists-1] / 4;
  iface       = lists[ilists-1] % 4;

  MMG5_tet2tri(mesh,iel,iface,&tt);
  for( i=0 ; i<3 ; i++ ) {
    if ( tt.v[i] == ip0 ) {
      break;
    }
  }

  assert ( i<3 );
  if ( i>=3 ) return 0;
  tt.v[i] = 0;

  if ( !MMG5_nortri(mesh, &tt, nprev) ) return 0;

  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilists ; l++ ){
    iel         = lists[l] / 4;
    iface       = lists[l] % 4;

    MMG5_tet2tri(mesh,iel,iface,&tt);
    caltmp = MMG5_caltri(mesh,met,&tt);
    calold = MG_MIN(calold,caltmp);

    for( i=0 ; i<3 ; i++ ) {
      if ( tt.v[i] == ip0 ) {
        break;
      }
    }
    assert(i<3);
    if ( i==3 ) return 0;

    tt.v[i] = 0;

    caltmp = MMG5_caltri(mesh,met,&tt);
    if ( caltmp < MMG5_EPSD2 ) {
      /* We don't check the input triangle qualities, thus we may have a very
       * bad triangle in our mesh */
      return 0;
    }
    calnew = MG_MIN(calnew,caltmp);

    if ( !MMG5_nortri(mesh, &tt, ncur) ) return 0;

    if ( ( !(tt.tag[i] & MG_GEO) ) && ( !(tt.tag[i] & MG_NOM) ) ) {
      /* Check normal deviation between iel and the triangle facing ip0 */
      ier = MMG3D_normalAdjaTri(mesh,iel,iface, i,nneighi);
      if ( ier <=0 ) {
        return 0;
      }
      ier =  MMG5_devangle( ncur, nneighi, mesh->info.dhd );
      if ( ier <= 0 ) {
        return 0;
      }
    }

    i = MMG5_iprv2[i];

    if ( ( !(tt.tag[i] & MG_GEO) ) && ( !(tt.tag[i] & MG_NOM) ) ) {
      /* Check normal deviation between k and the previous triangle */
      ier =  MMG5_devangle( ncur, nprev, mesh->info.dhd );
      if ( ier<=0 ) {
        return 0;
      }
    }
    memcpy(nprev, ncur, 3*sizeof(double));
  }
  if ( calold < MMG5_EPSOK && calnew <= calold )    return 0;
  else if ( calnew < calold )    return 0;
  memset(pxp,0,sizeof(MMG5_xPoint));

  /* Test : check whether all volumes remain positive with new position of the point */
  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilistv ; l++ ){
    iel = listv[l] / 4;
    i0  = listv[l] % 4;
    pt  = &mesh->tetra[iel];
    pt0 = &mesh->tetra[0];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    calold = MG_MIN(calold, pt->qual);
    callist[l] = MMG5_orcal(mesh,met,0);
    if (callist[l] < MMG5_NULKAL) {
      return 0;
    }
    calnew = MG_MIN(calnew,callist[l]);
  }
  if ((calold < MMG5_EPSOK && calnew <= calold) ||
      (calnew < MMG5_EPSOK) || (calnew <= 0.3*calold)) {
    return 0;
  } else if (improve && calnew < calold) {
    return 0;
  }

  /* Update coordinates, normals, for new point */
  if ( PROctree )
    MMG3D_movePROctree(mesh, PROctree, ip0, o, p0->c);

  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  pxp = &mesh->xpoint[p0->xp];
  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

  p0->n[0] = to[0];
  p0->n[1] = to[1];
  p0->n[2] = to[2];

  if ( isrid ) {
    /* Copy the second normal for ridge point */
    pxp->n2[0] = no2[0];
    pxp->n2[1] = no2[1];
    pxp->n2[2] = no2[2];
  }

  memcpy(&met->m[6*ip0],&met->m[0],6*sizeof(double));

  for( l=0 ; l<ilistv ; l++ ){
    (&mesh->tetra[listv[l]/4])->qual = callist[l];
    (&mesh->tetra[listv[l]/4])->mark = mesh->mark;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.

 * \return 0 if fail, 1 if success.
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 * Move boundary reference point, whose volumic and surfacic balls are passed.
 *
 */
int MMG5_movbdyrefpt_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG3D_pPROctree PROctree, int *listv,
                          int ilistv, int *lists, int ilists,int improve){

  return MMG3D_movbdycurvept_ani(mesh,met,PROctree,listv,ilistv,lists,ilists,improve,MG_REF);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.
 * \return 0 if fail, 1 if success.
 *
 * Move boundary non manifold point, whose volumic and (exterior)
 * surfacic balls are passed
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 */
int MMG5_movbdynompt_ani(MMG5_pMesh mesh,MMG5_pSol met, MMG3D_pPROctree PROctree, int *listv,
                          int ilistv, int *lists, int ilists,
                          int improve){

  return MMG3D_movbdycurvept_ani(mesh,met,PROctree,listv,ilistv,lists,ilists,improve,MG_NOM);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * \return 0 if fail, 1 if success.
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 * Move boundary ridge point, whose volumic and surfacic balls are passed.
 *
 */
int MMG5_movbdyridpt_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG3D_pPROctree PROctree, int *listv,
                          int ilistv,int *lists,int ilists,
                          int improve) {

  return MMG3D_movbdycurvept_ani(mesh,met,PROctree,listv,ilistv,lists,ilists,improve,MG_GEO);
}
