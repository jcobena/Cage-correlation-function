

/* ----------------------------------------------------------------------
   Contributing author:  Jose Cobena-Reyes
                         Spencer Ortega
------------------------------------------------------------------------- */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include "compute_ccf.h"
#include "fix_store.h"
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
#include <string>
#include <sstream>
#include <iostream>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

/* ---------------------------------------------------------------------- */

// Constructor function invoked when:
//     'compute ccf_test all ccf cutoff ...'
// is called in the input file
ComputeCCF::ComputeCCF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  qlist(NULL), distsq(NULL), nearest(NULL), rlist(NULL),
  qnarray(NULL), qnm_r(NULL), qnm_i(NULL)
{
  error->warning(FLERR,"Entering Constructor");
  // Validate and Process Arguments
  if (narg < 3 ) error->all(FLERR,"Illegal compute ccf command");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      // no cutoff value
      if(iarg+1 == narg)
        error->all(FLERR,"Illegal compute ccf command: no value for cutoff");
      // need at least one more arg for value of cutoff
      else if (iarg+2 < narg) 
        error->all(FLERR,"Illegal compute ccf command: more than one value for cutoff");
      
      // convert cutoff using force.numeric()
      double cutoff = force->numeric(FLERR,arg[iarg+1]);

      // no negative cutoff
      if (cutoff < 0.0) 
        error->all(FLERR,"Illegal compute ccf command: no negative cutoff values");
      else{
        // compute cutoff sqrd value
        cutsq = cutoff*cutoff;
        iarg += 2;
      }

    } 
    else error->all(FLERR,"Illegal compute ccf command: keyword does not exist");
  }

  // Initialize member variables
  maxneigh = 0;

  // Set flags for compute style methods that will be implemented/executed
  peratom_flag = 1;

  //flag = FALSE;

}

/* ---------------------------------------------------------------------- */

ComputeCCF::~ComputeCCF()
{
  error->warning(FLERR,"Entering deconstructor");
  memory->destroy(qnarray);
  memory->destroy(distsq);
  memory->destroy(rlist);
  memory->destroy(nearest);
  memory->destroy(qlist);
  memory->destroy(qnm_r);
  memory->destroy(qnm_i);

}

/* ---------------------------------------------------------------------- */

// init function invoked when:
//     'fix 2nve all nve
//      run       20000'
// is called in the input file
void ComputeCCF::init()
{
  error->warning(FLERR,"Entering init");

  // Validate cutoff
  // if cutsq == 0, set cutsq = (force cutoff)^2
  if(cutsq == 0.0) {
    cutsq = force->pair->cutforce * force->pair->cutforce;
  }
  // if sqrt(cutsq) > force cut error!
  else if (sqrt(cutsq) > force->pair->cutforce){
    error->all(FLERR,"Compute ccf cutoff is longer than pairwise cutoff");
  }

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
// init function invoked when:
//     'fix 2nve all nve
//      run       20000'
// is called in the input file
// after init()
void ComputeCCF::init_list(int /*id*/, NeighList *ptr)
{
  error->warning(FLERR,"Entering init_list");
  list = ptr;
}

/* ---------------------------------------------------------------------- */

// init function invoked when:
//     'fix 2nve all nve
//      run       20000'
// is called in the input file
// after init_list()
void ComputeCCF::compute_peratom()
{
  error->warning(FLERR,"Entering compute_peratom");

  // variables
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // update t?
  invoked_peratom = update->ntimestep;

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

  // stringstream ss;
  // ss << inum; 
  // string message = "\ninum = " + ss.str() + "\n";
  cout << "cutsq=" << cutsq << endl;
  cout << "inum=" << inum << endl;
  
  // for every atom id 0 to inum
  for (ii = 0; ii < inum; ii++) {

    cout << "i=" << ii << " : ";

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

      

      cout << "jnum=" << jnum << " : " ;

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
      
      cout << "ncount=" << ncount << endl;

      // if (flag == FALSE){
      //  //store t0
      //  flag = TRUE;
      // }
      // else{}

    }

  }
  error->warning(FLERR,"finishing compute peratom");
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCCF::memory_usage()
{
  error->warning(FLERR,"Entering memory_usage");
  double bytes = ncol*nmax * sizeof(double);
  bytes += (qmax*(2*qmax+1)+maxneigh*4) * sizeof(double);
  bytes += (nqlist+maxneigh) * sizeof(int);
  return bytes;
}
   //JCO added the parameters to create the fix store, the value of nwewardg[5] should be the length of the vector array of the neighbots of atom i
	// create a new fix STORE style for reference positions
	// id = compute-ID + COMPUTE_STORE, fix group = compute group

	int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
	id_fix = new char[n];
	strcpy(id_fix, id);
	strcat(id_fix, "_COMPUTE_STORE");

	char **newarg = new char*[6];
	newarg[0] = id_fix;
	newarg[1] = group->names[igroup];
	newarg[2] = (char *) "STORE";
	newarg[3] = (char *) "peratom";
	newarg[4] = (char *) "1";
	newarg[5] = (char *) "3";
	modify->add_fix(6, newarg);
	fix = (FixStore *)modify->fix[modify->nfix - 1];
	delete[] newarg;
}

//JCO added to start testing the storing positions 
void ComputeCCF::set_arrays(int i)
{
	double **xoriginal = fix->astore;
	double **x = x->atom;
	xoriginal[i][0] = x[i][0];
	xoriginal[i][1] = x[i][1];
	xoriginal[i][2] = x[i][2];
	
}

