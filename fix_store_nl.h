/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(STORE_NL,FixStoreNl)

#else

#ifndef LMP_FIX_STORE_NL_H
#define LMP_FIX_STORE_NL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStoreNl : public Fix {
  friend class ComputeCCF;

 public:
  // constructor
  FixStoreNl(class LAMMPS *,int, char **);
  // deconstructor
  virtual ~FixStoreNl();
  // required
  int setmask();
  // before run
  void init();
  // initialize neighbor list
  void init_list(int, class NeighList *);
  // function that stores nls
  void setup(int);

  // need for data migration
  double memory_usage();
  // allocate memory for atom-based arrays
  void grow_arrays(int);
  // copy atom info when an atom migrates to a new processor
  void copy_arrays(int, int, int);
  // store atom’s data in a buffer (optional)
  int pack_exchange(int, double *);
  // retrieve atom’s data from a buffer (optional)
  int unpack_exchange(int, double *);

  // checkpoint
  void write_restart(FILE *);
  // void restart(char *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  // int size_restart(int);
  // int maxsize_restart();


 protected:
  int first;                 // flag for first time initialization
  int maxpartner;            // max # of peridynamic neighs for any atom
  int *npartner;             // # of neighbors for each atom
  tagint **partner;          // neighs for each atom, stored as global IDs
  // double **deviatorextention; // Deviatoric extention
  // double **deviatorBackextention; // Deviatoric back extention
  // double **deviatorPlasticextension; // Deviatoric plastic extension
  // double *lambdaValue;
  // double **r0;               // initial distance to partners
  // double **r1;               // instanteneous distance to partners
  // double *thetaValue;        // dilatation
  // double *vinter;            // sum of vfrac for bonded neighbors
  // double *wvolume;           // weighted volume of particle
  // int isPMB,isLPS,isVES,isEPS;  // which flavor of PD

  class NeighList *list;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Duplicate particle in PeriDynamic bond - simulation box is too small

This is likely because your box length is shorter than 2 times
the bond length.

*/
