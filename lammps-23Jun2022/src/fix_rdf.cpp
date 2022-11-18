//
// Created by Ezra Cerpac on 18/11/2022.
//
// Structure taken from fix_aveforce.cpp

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_rdf.h"

#include "atom.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

FixRDF::FixRDF(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    Rxx(nullptr), Ryy(nullptr), Rzz(nullptr), Fxx(nullptr), Fyy(nullptr), Fzz(nullptr),
    hist(nullptr)
{
    if (narg < 6) error->all(FLERR, "Illegal fix rdf command");

    nbin = utils::inumeric(FLERR, arg[3], false, lmp);

//    minTimestep = 20000;
    Temp = 299;             // Find method to get this
    Box = 10;               // Random number, Find method to get this
    natom = atom->natoms;

    memory->create(Rxx, natom, "rdf:Rxx");  // Find method to get this
    memory->create(Ryy, natom, "rdf:Ryy");
    memory->create(Rzz, natom, "rdf:Rzz");
    memory->create(Fxx, natom, "rdf:Fxx");
    memory->create(Fyy, natom, "rdf:Fyy");
    memory->create(Fzz, natom, "rdf:Fzz");
    memory->create(hist, nbin, "rdf:hist");

}

/* ---------------------------------------------------------------------- */

FixRDF::~FixRDF() {
    memory->destroy(Rxx);
    memory->destroy(Ryy);
    memory->destroy(Rzz);
    memory->destroy(Fxx);
    memory->destroy(Fyy);
    memory->destroy(Fzz);
    memory->destroy(hist);
}

/* ---------------------------------------------------------------------- */

int FixRDF::setmask() {
    int mask = 0;
    mask |= FixConst::POST_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixRDF::init() {
    // Check precondition
}

/* ---------------------------------------------------------------------- */

void FixRDF::setup(int vflag) {
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRDF::post_force(int vflag) {
    int i, j, bin;
    double const *x, *y, *z, *fx, *fy, *fz;
    double Dx, Dy, Dz;
    double Rij, Ff, R;
    double HalfBox = Box / 2.;
    const double Delta = 0.005; // Why?

    x = atom->x[0];
    y = atom->x[1];
    z = atom->x[2];
    fx = atom->f[0];
    fy = atom->f[1];
    fz = atom->f[2];

    for (i = 0; i < natom; i++) {
        Rxx[i] = x[i];
        Ryy[i] = y[i];
        Rzz[i] = z[i];
        Fxx[i] = fx[i];
        Fyy[i] = fy[i];
        Fzz[i] = fz[i];
        hist[i] = 0.;
    }

    for (i = 0; i < natom; i++) {
        for (j = i + 1; j < natom; j++) {
            // Calculate the distance between the two particles, R, accounting
            // for periodic boundary conditions
            Dx = Rxx[i] - Rxx[j];
            Dy = Ryy[i] - Ryy[j];
            Dz = Rzz[i] - Rzz[j];

            if (Dx > HalfBox) {
                Dx = Dx - Box;
            } else if (Dx < -HalfBox) {
                Dx = Dx + Box;
            }

            if (Dy > HalfBox) {
                Dy = Dy - Box;
            } else if (Dy< -HalfBox) {
                Dy = Dy + Box;
            }

            if (Dz > HalfBox) {
                Dz = Dz - Box;
            } else if (Dz < -HalfBox) {
                Dz = Dz + Box;
            }

            Rij = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);

            // RDF is only valid till half the box size
            if (Rij < HalfBox) {
                Ff = (Fxx[i]*Dx + Fyy[i]*Dy + Fzz[i]*Dz) / (Rij*Rij*Rij);

                for (bin = 0; bin < nbin; bin++) {
                    R = (bin + 0.5) * Delta;

                    if (R < HalfBox && R > Rij)
                        hist[bin] += Ff;
                }
            }
        }
    }

    const double c = Box * Box * Box / (2 * M_PI * Temp * natom * natom);

    for (bin = 0; bin < nbin; bin++) {
        hist[bin] *= c;
    }
}

/* ---------------------------------------------------------------------- */

double FixRDF::compute_vector(int n) {
    return hist[n];
}
