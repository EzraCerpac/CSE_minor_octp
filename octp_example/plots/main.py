#!/usr/local/bin/python3.10

import matplotlib.pyplot as plt
import numpy as np


def main():
    rdf_octp = np.loadtxt('../rdf.dat')
    rdf_lammps = np.loadtxt('../rdf_lammps.dat', skiprows=6)
    plt.plot(rdf_octp[:,0], rdf_octp[:,1], label='octp')
    plt.plot(rdf_lammps[:,1], rdf_lammps[:,2], label='lammps')
    # plt.plot(df[:,0], df[:,2])
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.legend()
    # plt.savefig('rdf.pdf')
    plt.show()


if __name__ == '__main__':
    main()
