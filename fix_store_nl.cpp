/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mike Parks (SNL), Ezwanur Rahman, J.T. Foster (UTSA)
------------------------------------------------------------------------- */

#include <cmath>
#include "fix_store_nl.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "lattice.h"
#include "memory.h"
#include "error.h"

#include <iostream>
using namespace std;

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixStoreNL::FixStoreNL(LAMMPS *lmp,int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  cout << "\n-----------------\n\n\n"
       << "inside fix store nl constructor!"
       << "\n\n\n-----------------\n";
  // we need these flags
  restart_global = 1;
  restart_peratom = 1;
  first = 1;

  // perform initial allocation of atom-based arrays
  // register with atom class
  // set maxpartner = 1 as placeholder

  maxpartner = 1;
  npartner = NULL;
  partner = NULL;

  // keep for now
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize npartner to 0 so atom migration is OK the 1st time
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;

  // set comm sizes needed by this fix
  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixStoreNL::~FixStoreNL()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->destroy(npartner);
  memory->destroy(partner);
}

/* ---------------------------------------------------------------------- */

int FixStoreNL::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStoreNL::init()
{
  if (!first) return;

  // need a full neighbor list once

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

}

/* ---------------------------------------------------------------------- */

void FixStoreNL::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   create initial list of neighbor partners via call to neighbor->build()
   must be done in setup (not init) since fix init comes before neigh init
------------------------------------------------------------------------- */

void FixStoreNL::setup(int /*vflag*/)
{
  cout << "\n-----------------\n\n\n"
       << "inside fix store nl setup!"
       << "\n\n\n-----------------\n";
  int i,j,ii,jj,itype,jtype,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh;
  int **firstneigh;


  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;


  // only build list of bonds on very first run
  if (!first) return;
  first = 0;


  // build full neighbor list, will copy or build as necessary
  neighbor->build_one(list);
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  // scan neighbor list to set maxpartner
  //Pair *anypair = force->pair_match("peri",0);
  //double **cutsq = anypair->cutsq;
  // cout << "storenl: first for loop\n";
  double cutsq = 2.25;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // if (rsq <= cutsq[itype][jtype]) npartner[i]++;
      if (rsq <= cutsq) npartner[i]++;
    }
    // cout << "atom " << i << " has npartners: "<< npartner[i]<<endl;
  }

  // cout << "storenl: END first for loop\n";


  // cout << "storenl: start mpi\n";
  // global max neighbor list size
  maxpartner = 0;
  for (i = 0; i < nlocal; i++) maxpartner = MAX(maxpartner,npartner[i]);
  int maxall;
  MPI_Allreduce(&maxpartner,&maxall,1,MPI_INT,MPI_MAX,world);
  maxpartner = maxall;
  // cout << "storenl: END mpi\n";


  // realloc arrays with correct value for maxpartner
  memory->destroy(partner);
  memory->destroy(npartner);
  npartner = NULL;
  partner = NULL;
  grow_arrays(atom->nmax);


  // save actual values of nl
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // cout << "starting for loop for " << i << endl;
    npartner[i] = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // if (rsq <= cutsq[itype][jtype]) {
      if (rsq <= cutsq){
        // store in array

        partner[i][npartner[i]] = atom->tag[j];

        // increment count
        npartner[i]++;
      }
    }
    // cout << "npartner=" << npartner[i] << endl;
  }

  // sanity check: does any atom appear twice in any neigborlist?
  // should only be possible if using pbc and domain < 2*delta

  // Need this in our code

  // if (domain->xperiodic || domain->yperiodic || domain->zperiodic) {
  //   for (i = 0; i < nlocal; i++) {
  //     jnum = npartner[i];
  //     for (jj = 0; jj < jnum; jj++) {
  //       for (int kk = jj+1; kk < jnum; kk++) {
  //         if (partner[i][jj] == partner[i][kk])
  //           error->one(FLERR,"Duplicate particle in PeriDynamic bond - "
  //                      "simulation box is too small");
  //       }
  //     }
  //   }
  // }
  cout << "finsished store nl"<<endl;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixStoreNL::memory_usage()
{
  int nmax = atom->nmax;
  int bytes = nmax * sizeof(int);
  bytes += nmax*maxpartner * sizeof(tagint);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixStoreNL::grow_arrays(int nmax)
{
   memory->grow(npartner,nmax,"peri_neigh:npartner");
   memory->grow(partner,nmax,maxpartner,"peri_neigh:partner");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixStoreNL::copy_arrays(int i, int j, int /*delflag*/)
{
  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[j]; m++) {
    partner[j][m] = partner[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixStoreNL::pack_exchange(int i, double *buf)
{
  // compact list by eliminating partner = 0 entries
  // set buf[0] after compaction

  int m = 1;
  for (int n = 0; n < npartner[i]; n++) {
    if (partner[i][n] == 0) continue;
    buf[m++] = partner[i][n];
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixStoreNL::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  // size of packed nl
  npartner[nlocal] = static_cast<int> (buf[m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<tagint> (buf[m++]);
  }

  return m;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixStoreNL::write_restart(FILE *fp)
{
  int n = 0;
  double list[2];
  list[n++] = first;
  list[n++] = maxpartner;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixStoreNL::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixStoreNL::unpack_restart(int nlocal, int nth)
{

  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<tagint> (extra[nlocal][m++]);
  }
}
