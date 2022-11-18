//
// Created by Ezra Cerpac on 18/11/2022.
//

#ifdef FIX_CLASS
// clang-format off
FixStyle(rdf,FixRDF);
// clang-format on
#else

#ifndef LAMMPS_FIX_RDF_H
#define LAMMPS_FIX_RDF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRDF : public Fix {
public:
    FixRDF(class LAMMPS *, int, char **);
    ~FixRDF() override;
    int setmask() override;
    void init() override;
    void setup(int) override;
    void post_force(int) override;
    double compute_vector(int) override;

private:
    int nbin;                // # of rdf bins
    bigint natom;            // # of atoms
    double Temp, Box, *Rxx, *Ryy, *Rzz, *Fxx, *Fyy, *Fzz, *hist;
};

}   // namespace LAMMPS_NS

#endif
#endif