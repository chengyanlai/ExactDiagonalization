#!/usr/bin/env python
import numpy as np
import h5py
import os
import argparse

L = 3

def main(InputFile, OutputFile, Prefix):
    f = h5py.File(os.path.join(Prefix, InputFile), "r")
    Nf1 = []
    Nf2 = []
    Nf3 = []
    Jf1 = []
    Jf2 = []
    Jf3 = []
    Np1 = []
    Np2 = []
    Np3 = []
    Xp1 = []
    Xp2 = []
    Xp3 = []
    cnt = 0
    gname = "Obs-" + str(cnt)
    e = gname in f.keys()
    while e:
        Fmn = f[gname]["Fermion"]["Elem"]["real"][:] + 1.0j * f[gname]["Fermion"]["Elem"]["imag"][:]
        Fmn = np.reshape(Fmn, (L,L))
        Nf1.append( Fmn[0,0].real )
        Nf2.append( Fmn[1,1].real )
        Nf3.append( Fmn[2,2].real )
        Jf1.append( Fmn[0,1].imag )
        Jf2.append( Fmn[1,2].imag )
        Jf3.append( Fmn[2,0].imag )
        Pm = f[gname]["Phonon"]["Elem"]["real"][:]
        Np1.append( Pm[0] )
        Np2.append( Pm[1] )
        Np3.append( Pm[2] )
        Xm = f[gname]["Ei"][:]
        Xp1.append( Xm[0] )
        Xp2.append( Xm[1] )
        Xp3.append( Xm[2] )
        cnt += 20
        gname = "Obs-" + str(cnt)
        e = gname in f.keys()
    Nf1 = np.array( Nf1 )
    Nf2 = np.array( Nf2 )
    Nf3 = np.array( Nf3 )
    Jf1 = np.array( Jf1 )
    Jf2 = np.array( Jf2 )
    Jf3 = np.array( Jf3 )
    Np1 = np.array( Np1 )
    Np2 = np.array( Np2 )
    Np3 = np.array( Np3 )
    Xp1 = np.array( Xp1 )
    Xp2 = np.array( Xp2 )
    Xp3 = np.array( Xp3 )
    np.savetxt(os.path.join(Prefix, OutputFile + ".txt"), np.vstack([Nf1, Nf2, Nf3, Jf1, Jf2, Jf3, Np1, Np2, Np3, Xp1, Xp2, Xp3]).T )
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="TrJm", usage='%(prog)s [options] -f file',
                                     description="Chenyen's personal command to make life easier.",
                                     epilog="All's well that ends well.")
    parser.add_argument('-i', '--input', default='QuenchState-Z-0.h5', help='Filename.')
    parser.add_argument('-o', '--output', default='Jtot', help='Filename Prefix.')
    parser.add_argument('--from', type=int, default=0, help='subfolder sequence.')
    parser.add_argument('--to', type=int, default=1000, help='subfolder sequence.')
    args = vars(parser.parse_args())
    for i in range(args["from"], args["to"], 1):
        prefix = "S"+str(i).zfill(4)
        main(args["input"], args["output"], prefix)
