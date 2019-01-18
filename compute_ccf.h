/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories

------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(ccf,ComputeCCF)

#else

#ifndef LMP_COMPUTE_CCF_H
#define LMP_COMPUTE_CCF_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCCF : public Compute {
 public:
  ComputeCCF(class LAMMPS *, int, char **);
  ~ComputeCCF();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();
  double cutsq;
  int iqlcomp, qlcomp, qlcompflag;
  int *qlist;
  int nqlist;

 private:
  int nmax,maxneigh,ncol,nnn;
  class NeighList *list;
  double *distsq;
  int *nearest;
  double **rlist;
  int qmax;
  double **qnarray;
  double **qnm_r;
  double **qnm_i;
  int ifix_storenl;

  void select3(int, int, double *, int *, double **);
  void calc_boop(double **rlist, int numNeighbors,
                 double qn[], int nlist[], int nnlist);
  double dist(const double r[]);

  double polar_prefactor(int, int, double);
  double associated_legendre(int, int, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:
E: Illegal ... command
Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.
E: Compute orientorder/atom requires a pair style be defined
Self-explanatory.
E: Compute orientorder/atom cutoff is longer than pairwise cutoff
Cannot compute order parameter beyond cutoff.
W: More than one compute orientorder/atom
It is not efficient to use compute orientorder/atom more than once.
*/