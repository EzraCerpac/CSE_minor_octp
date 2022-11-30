#!/usr/local/bin/python3.10

import matplotlib.pyplot as plt
import numpy as np
import pickle


def main():
    path = "data.pkl"
    reload = True #if input("Reload data? [y/N]: ") == "y" else False
    if reload:
        rdf_octp, rdf_force = load_data(path)
    else:
        try:
            with open(path, 'rb') as f:
                rdf_octp, rdf_force = pickle.load(f)
        except FileNotFoundError:
            rdf_octp, rdf_force = load_data(path)
    plt.plot(rdf_octp[:, 0], rdf_octp[:, 1], label='octp')
    plt.plot(rdf_force[:, 1], rdf_force[:, 2], label='force')
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.legend()
    plt.grid()
    plt.minorticks_on()
    # plt.savefig('rdf.pdf')
    plt.show()


def load_data(path: str = 'data.pkl') -> tuple[np.ndarray, np.ndarray]:
    rdf_octp = np.loadtxt('../rdf.dat', max_rows=2000)
    rdf_force = np.loadtxt('../rdf_force.dat', skiprows=4, max_rows=2000)
    with open(path, 'wb') as f:
        pickle.dump((rdf_octp, rdf_force), f)
    return rdf_octp, rdf_force


if __name__ == '__main__':
    main()
