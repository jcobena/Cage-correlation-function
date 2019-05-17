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

// Add on
#include "fix_store_nl.h"

// DEBUGING PURPOSES
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>




using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;



/* ---------------------------------------------------------------------- */

ComputeCCF::ComputeCCF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  cout << "ComputeCCF-constructor" << endl;

  // cout << "\n-----------------\n\n\n"
  //      << "inside compute ccf constructor!"
  //      << "\n\n\n-----------------\n";


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

  char **fixarg = new char*[4];
  fixarg[0] = (char *) "store_nl";
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "store_nl";
  fixarg[3] = (char *) arg[4];
  modify->add_fix(4,fixarg);
  delete [] fixarg;
}

/* ---------------------------------------------------------------------- */

ComputeCCF::~ComputeCCF()
{
  cout << "ComputeCCF deconstructor" << endl;


  // memory->destroy(qlist);
  // memory->destroy(qnarray);
  // memory->destroy(nearest);
}

/* ---------------------------------------------------------------------- */

void ComputeCCF::init()
{

  cout << "ComputeCCF - init" << endl;

  // check pair style
  if (force->pair == NULL)
    error->all(FLERR,"Compute ccf requires a "
               "pair style be defined");

  // check cutoff
  if (cutsq == 0.0) cutsq = force->pair->cutforce * force->pair->cutforce;
  else if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute ccf cutoff is "
               "longer than pairwise cutoff");
  
  // set up neighbor list request
  // need an occasional full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;


  // verify only one compute ccf called
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"ccf") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute ccf");

  // verify that fix store nl is called with compute ccf
  ifix_storenl = -1;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"store_nl") == 0) ifix_storenl = i;
  if (ifix_storenl == -1)
    error->all(FLERR,"Compute ccf requires store nl fix");

}

/* ---------------------------------------------------------------------- */

void ComputeCCF::init_list(int /*id*/, NeighList *ptr)
{
  cout << "ComputeCCF-initlist" << endl;

  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCCF::compute_peratom()
{

  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  cout << "ComputeCCF-computeperatom-" << invoked_peratom << endl;


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


  //-------------------------------------------------------------------
  // DEBUGGING
  // create file stream object
  ofstream myfile;
  // stringify mpi rank number, time index
  stringstream ss, tt;
  // mpi rank
  ss << comm->me;
  string p;
  ss >> p;
  // time index
  tt << invoked_peratom;
  string time;
  tt >> time;
  // open file
  string filename = "./out_compute/t" + time + "-mpi"+ p +".txt";
  myfile.open(filename.c_str());

  // Retrieve NL from time = 0 via fix store nl
  myfile << "before fix store call" << endl;
  cout << "before fix store call" << comm->me << endl;
  tagint **partner = ((FixStoreNL *) modify->fix[ifix_storenl])->partner;
  int *npartner = ((FixStoreNL *) modify->fix[ifix_storenl])->npartner;
  myfile << "after fix store call" << endl;
  cout << "after fix store call" << comm->me << endl;


  // for every local atom in rank p
  for(int i = 0; i < atom->nlocal; i++){
    // print real id, number of neighbors at t=0, real ids of neighbors
    // assumption: 'npartner'/'partner' ordered index 'i' corresponds to local atom id
    // ex:
    //     local atom 0: number of neighbors = npartner[0]
    //                   neighbor list       = partner[0]
    //     local atom 9: number of neighbors = npartner[9]
    //                   neighbor list       = partner[9]
    
    // myfile << "real id = " << atom->tag[i] << "\nnum neighbors = " << npartner[i] << endl; 
    // myfile << atom->tag[i] << endl; 
    // for every neighbor (j) of i
    // for(int j = 0; j < npartner[i]; j++){
    //   myfile << partner[i][j] << endl;
    // }
    // myfile << "-----------------" <<endl;
  }


  //-------------------------------------------------------------------


  // vector_atom gets outputted to dump file
  vector_atom = new double[inum];
  for(int i = 0; i < inum; i++){
    vector_atom[i] = -1;
  }


  // for every local atom
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];

    myfile << i << " , " << atom->tag[i] << endl;
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


      // int* nearest = new int[jnum];

      int ncount = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        myfile << "    " << j << " , " << atom->tag[j] << endl;
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          ncount++;

        }
      }

      // store neighbor count for each local atom i
      vector_atom[ii] = ncount;
    }
  }

  myfile << "finished compute_peratom" << endl;
  cout << "finished compute_peratom" << comm->me << endl;
  t++;
  myfile.close();



}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCCF::memory_usage()
{

  cout << "ComputeCCF-memoryusage" << endl;

  // double bytes = ncol*nmax * sizeof(double);
  // bytes += (qmax*(2*qmax+1)+maxneigh*4) * sizeof(double);
  // bytes += (nqlist+maxneigh) * sizeof(int);
  // return bytes;
  if (comm->me == 0){
  cout << atom->nmax << endl;
  cout << nmax << endl;}
  double bytes = atom->nmax * sizeof(double);
  return bytes;
}