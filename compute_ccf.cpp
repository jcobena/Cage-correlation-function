/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

/* ----------------------------------------------------------------------
   Contributing author:  Jose Cobena-Reyes
                         Spencer Ortega
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

// #include "fix_store_nl.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;



/* ---------------------------------------------------------------------- */

ComputeCCF::ComputeCCF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  cout << "\n-----------------\n\n\n"
       << "inside compute ccf constructor!"
       << "\n\n\n-----------------\n";


  // set default values for optional args

  //nnn = 12;
  cutsq = 0.0;
  // qlcompflag = 0;

  // specify which orders to request

  // nqlist = 5;
  // memory->create(qlist,nqlist,"ccf:qlist");
  // qlist[0] = 4;
  // qlist[1] = 6;
  // qlist[2] = 8;
  // qlist[3] = 10;
  // qlist[4] = 12;
  // qmax = 12;


  // process args
  if (narg < 5 ) error->all(FLERR,"Illegal compute ccf command");

  int iarg = 3;
  while (iarg < narg) {
      if (strcmp(arg[iarg],"cutoff") == 0) {
      double cutoff = force->numeric(FLERR,arg[iarg+1]);
      if (cutoff <= 0.0)
        error->all(FLERR,"Illegal compute ccf command");
      cutsq = cutoff*cutoff;
      iarg += 2;
    } else error->all(FLERR,"Illegal compute ccf command");
  }



  // if (qlcompflag) ncol = nqlist + 2*(2*qlcomp+1);
  // else ncol = nqlist;

  peratom_flag = 1;
  size_peratom_cols = 0; //0 means vector, same than in damage atom
  size_vector_variable = 1;

  //nmax = 1;
  // maxneigh = 0;

  // char **fixarg = new char*[3];
  // fixarg[0] = (char *) "store_nl";
  // fixarg[1] = (char *) "all";
  // fixarg[2] = (char *) "store_nl";
  // modify->add_fix(3,fixarg);
  // delete [] fixarg;
}

/* ---------------------------------------------------------------------- */

ComputeCCF::~ComputeCCF()
{
  // memory->destroy(qlist);
  // memory->destroy(qnarray);
  // memory->destroy(nearest);
}

/* ---------------------------------------------------------------------- */

void ComputeCCF::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute ccf requires a "
               "pair style be defined");
  if (cutsq == 0.0) cutsq = force->pair->cutforce * force->pair->cutforce;
  else if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute ccf cutoff is "
               "longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"ccf") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute ccf");

  // ifix_storenl = -1;
  // for (int i = 0; i < modify->nfix; i++)
  //   if (strcmp(modify->fix[i]->style,"store_nl") == 0) ifix_storenl = i;
  // if (ifix_storenl == -1)
  //   error->all(FLERR,"Compute ccf requires store nl fix");

}

/* ---------------------------------------------------------------------- */

void ComputeCCF::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCCF::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  // if (atom->nmax > nmax) {
  //   // memory->destroy(qnarray);
  //   nmax = atom->nmax;
  //   // memory->create(qnarray,nmax,ncol,"ccf:qnarray");
  //   // array_atom = qnarray;
  // }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute order parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *mask = atom->mask;

  // Retrieve NL at time = 0 from fix
  // tagint **partner = ((FixStoreNL *) modify->fix[ifix_storenl])->partner;
  // int *npartner = ((FixStoreNL *) modify->fix[ifix_storenl])->npartner;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // insure distsq and nearest arrays are long enough

      // if (jnum > maxneigh) {
      //   memory->destroy(nearest);
      //   maxneigh = jnum;
      //   memory->create(nearest,maxneigh,"ccf:nearest");
      // }

      // loop over list of all neighbors within force cutoff
      // nearest[] = atom indices of neighbors

      vector_atom = new double[jnum];
      // for(int iii = 0; iii < jnum; iii++){
      //   vector_atom[iii] = 0;
      // }
      int* nearest = new int[jnum];

      int ncount = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {

          nearest[ncount++] = j;
        }
      }

      // compare nearest with partner[i]

      // for(int iii = 0; )

      delete [] nearest;
      delete [] vector_atom;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCCF::memory_usage()
{
  // double bytes = ncol*nmax * sizeof(double);
  // bytes += (qmax*(2*qmax+1)+maxneigh*4) * sizeof(double);
  // bytes += (nqlist+maxneigh) * sizeof(int);
  // return bytes;
  double bytes = nmax * sizeof(double);
  return bytes;
}