/* -*- c++ -*- ----------------------------------------------------------
  Compute CCF
------------------------------------------------------------------------- */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include "compute_ccf.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;


/* ---------------------------------------------------------------------- */

ComputeCCF::ComputeCCF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  qlist(NULL), distsq(NULL), nearest(NULL), rlist(NULL),
  qnarray(NULL), qnm_r(NULL), qnm_i(NULL)
{
  // validate and process args
    // cutoff,... etc

  // check if enough arguments
  // at least needs:
  // 1: compute ID
  // 2: group ID
  // 3: style name (ccf)
  if (narg < 3 ) error->all(FLERR,"Illegal compute ccf command");

  // set default values for optional args
  cutsq = 0.0;

  // process optional args
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      // need at least one more arg for value of cutoff
      if (iarg+1 >= narg)
        error->all(FLERR,"Illegal compute ccf command: more than one value for cutoff");
      double cutoff = force->numeric(FLERR,arg[iarg+1]);
      // no negative cutoff
      if (cutoff <= 0.0)
        error->all(FLERR,"Illegal compute ccf command: no negative cutoff values");
      cutsq = cutoff*cutoff;
      iarg += 2;
    } 
    else error->all(FLERR,"Illegal compute ccf command: keyword does not exist");
  }
}

/* ---------------------------------------------------------------------- */

ComputeCCF::~ComputeCCF()
{
  memory->destroy(qnarray);
  memory->destroy(distsq);
  memory->destroy(rlist);
  memory->destroy(nearest);
  memory->destroy(qlist);
  memory->destroy(qnm_r);
  memory->destroy(qnm_i);
}

// set parameters for neighbor list
void ComputeCCF::init()
{
  // neighbor list at t0 and tt
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

}

// initializes neighbor list
void ComputeCCF::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

// function that creates and uses neighbor list
void ComputeCCF::compute_ccf()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;


  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  double **x = atom->x;
  int *mask = atom->mask;

  // for every atom id 0 to inum
  for (ii = 0; ii < inum; ii++) {
    // get numneigh index for ii
    i = ilist[ii];
    // jum = length of numneigh for atom ii
    jnum = numneigh[i];


      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // rlist[] = distance vector to each
      // nearest[] = atom indices of neighbors

      // for every neighbor with "force cutoff"??
      // -ask jose
    int ncount = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      // looks like we are checking the distance of every atom for cutoffs
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // if distance less than cutoff store id
      if (rsq < cutsq) {
        distsq[ncount] = rsq;
        rlist[ncount][0] = delx;
        rlist[ncount][1] = dely;
        rlist[ncount][2] = delz;
        nearest[ncount++] = j;
      }
    }
  }
}


