#!/usr/bin/env python
import numpy as np
import h5py
import argparse

L = 3

def main(filename, prefix):
    f = h5py.File(filename, "r")
    Nf1 = []
    Nf2 = []
    Nf3 = []
    Jf1 = []
    Jf2 = []
    Jf3 = []
    for i in f.keys():
        try:
            Fmn = f[i]["Fermion"]["Elem"][:]
            Nf1.append( Fmn[0,0].real )
            Nf2.append( Fmn[1,1].real )
            Nf3.append( Fmn[2,2].real )
            Jf1.append( Fmn[0,1].imag )
            Jf2.append( Fmn[1,2].imag )
            Jf3.append( Fmn[0,2].imag )
        except:
            pass
    Nf1 = np.array( Nf1 )
    Nf2 = np.array( Nf2 )
    Nf3 = np.array( Nf3 )
    Jf1 = np.array( Jf1 )
    Jf2 = np.array( Jf2 )
    Jf3 = np.array( Jf3 )
    np.savetxt(prefix + ".txt", np.vstack([Nf1, Nf1, Nf3, Jf1, Jf2, Jf3]).T )
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="TrJm", usage='%(prog)s [options] -f file',
                                     description="Chenyen's personal command to make life easier.",
                                     epilog="All's well that ends well.")
    parser.add_argument('-i', '--input', default='QuenchState-Z-0.h5', help='Filename.')
    parser.add_argument('-o', '--output', default='Jtot', help='Filename Prefix.')
    args = vars(parser.parse_args())
    main(args["input"], args["output"])
