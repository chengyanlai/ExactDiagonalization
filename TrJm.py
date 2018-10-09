#!/usr/bin/env python
import numpy as np
import h5py
import argparse

Beta = np.logspace(1, -2, 121)

def main(filename, prefix):
    f = h5py.File(filename, "r")
    for i in f.keys():
        try:
            Jm = f[i]["Jmm"]["Elem"][:]
            Em = f[i]["Em"]["Elem"][:]
            Jm = Jm * Jm
            TrEm = []
            TrJm = []
            Ef = []
            for b in Beta:
                Z = np.exp(-b*Em)
                Ef.append(Z[-1])
                TrEm.append(np.sum(Z))
                TrJm.append( np.sum( Z * Jm ) )
            TrJm = np.array(TrJm)
            TrEm = np.array(TrEm)
            Ef = np.array(Ef)
            np.savetxt(prefix + i + ".txt", np.vstack([Beta, TrJm, TrEm, Ef]).T)
        except:
            pass
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="TrJm", usage='%(prog)s [options] -f file',
                                     description="Chenyen's personal command to make life easier.",
                                     epilog="All's well that ends well.")
    parser.add_argument('-i', '--input', default='Holstein.K.h5', help='Filename.')
    parser.add_argument('-o', '--output', default='TrJm', help='Filename Prefix.')
    args = vars(parser.parse_args())
    main(args["input"], args["output"])
