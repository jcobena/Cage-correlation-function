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

// constructor
ComputeCCF::ComputeCCF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), nearest(NULL)
{
  // Validate and Process Arguments
  if (narg < 3 ) error->all(FLERR,"Illegal compute ccf command");
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      // need at least one more arg for value of cutoff
      if (iarg+1 >= narg) error->all(FLERR,"Illegal compute ccf command: more than one value for cutoff");
      // convert cutoff using force.numeric()
      double cutoff = force->numeric(FLERR,arg[iarg+1]);
      // no negative cutoff
      if (cutoff <= 0.0) error->all(FLERR,"Illegal compute ccf command: no negative cutoff values");
      // save cutoff value
      cutsq = cutoff*cutoff;
      iarg += 2;
    } 
    else error->all(FLERR,"Illegal compute ccf command: keyword does not exist");
  }

  // Initialize member variables
  cutsq = 0.0;
  maxneigh = 0;

  // Set flags for compute style methods that will be implemented/executed
  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

// deconstructor
ComputeCCF::~ComputeCCF()
{
  memory->destroy(nearest);
}

/* ---------------------------------------------------------------------- */

// set parameters for neighbor list
void ComputeCCF::init()
{
  // Validate cutoff
  // if cutsq == 0, set cutsq = (force cutoff)^2
  if (cutsq == 0.0) 
    cutsq = force->pair->cutforce * force->pair->cutforce;
  // if sqrt(cutsq) > force cut error!
  else if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute ccf cutoff is longer than pairwise cutoff");  
  
  // Initialize neighbor list parameters
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  // Make sure there is only 1 compute style "ccf"
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"ccf") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute ccf");
}

/* ---------------------------------------------------------------------- */

// initializes neighbor list
void ComputeCCF::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

// function that creates and uses neighbor list
double ComputeCCF::compute_scalar()
{
  // variables
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // update t?
  invoked_scalar = update->ntimestep;

  // build neighbor list
  neighbor->build_one(list);

  // retrieve neighbor list information
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // retrieve atom coordinates at time t
  double **x = atom->x;
  int *mask = atom->mask;

  // for every atom id 0 to inum
  for (ii = 0; ii < inum; ii++) {
    // get numneigh index for ii
    i = ilist[ii];

    // if atom is part of selected group?
    if (mask[i] & groupbit) {
      // get coordinates of atom i
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      // get the list and size of neighbors at atom i
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // make sure nearest array has enough space
      if (jnum > maxneigh) {
        memory->destroy(nearest);
        maxneigh = jnum;
        memory->create(nearest,maxneigh,"ccf:nearest");
      }


      // for all atom i neighbors within force cutoff
      int ncount = 0;
      for (jj = 0; jj < jnum; jj++) {
        // get unique id j
        j = jlist[jj];
        j &= NEIGHMASK;
        
        // get coordinates at atom j
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        
        // compute rsq
        rsq = delx*delx + dely*dely + delz*delz;

        // if atom j is within compute cutoff
        if (rsq < cutsq) {
          // save id
          nearest[ncount++] = j;
        }
      }

      // neighbor list of atom i should be withing nearest[] array
      // print?
    }
  }

  // return value
  scalar=66.6;
  return scalar;
}