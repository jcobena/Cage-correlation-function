/* -*- c++ -*- ----------------------------------------------------------
  Compute CCF
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
  // constructor/deconstructor
  // check user args, set flags, initialize variables
  ComputeCCF(class LAMMPS *, int, char **);
  // delete memory
  ~ComputeCCF();
  
  // virtual functions from compute class
  // set up neighborlist parameters
  void init();
  // initializes neighborlist
  void init_list(int, NeighList *ptr);
  // compute function
  void compute_peratom();


 private:
  class NeighList *list;
  int cutsq;

}

#endif
#endif