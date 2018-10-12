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
  ComputeCCF(class LAMMPS *, int, char **);
  ~ComputeCCF();
  void init();
  void init_list(int, NeighList *ptr);

 private:
  class NeighList *list;

}

#endif
#endif