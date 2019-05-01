#!/usr/bin/env python
import numpy as np
import h5py
import os
import argparse

def main(L, InputFile, OutputFile, Prefix, Every):
    f = h5py.File(os.path.join(Prefix, InputFile), "r")
    Nf = {}
    Jf = {}
    Np = {}
    Xp = {}
    for i in range(L):
        Nf[str(i)] = []
        Jf[str(i)] = []
        Np[str(i)] = []
        Xp[str(i)] = []
    cnt = 0
    gname = "Obs-" + str(cnt)
    e = gname in f.keys()
    while e:
        Fmn = f[gname]["Fermion"]["Elem"]["real"][:] + 1.0j * f[gname]["Fermion"]["Elem"]["imag"][:]
        Fmn = np.reshape(Fmn, (L,L))
        Pm = f[gname]["Phonon"]["Elem"]["real"][:]
        Xm = f[gname]["Ei"][:]
        for i in range(L):
            Nf[str(i)].append( Fmn[i,i].real )
            Np[str(i)].append( Pm[i] )
            Xp[str(i)].append( Xm[i] )
            if i == L - 1:
                Jf[str(i)].append( Fmn[i,0].imag )
            else:
                Jf[str(i)].append( Fmn[i,i+1].imag )
        cnt += Every
        gname = "Obs-" + str(cnt)
        e = gname in f.keys()
    for i in range(L):
        try:
            out = np.vstack( [out, np.array(Nf[str(i)])] )
        except UnboundLocalError:
            out = np.array(Nf[str(i)])
        except:
            raise
    for i in range(L):
        try:
            out = np.vstack( [out, np.array(Jf[str(i)])] )
        except:
            raise
    for i in range(L):
        try:
            out = np.vstack( [out, np.array(Np[str(i)])] )
        except:
            raise
    for i in range(L):
        try:
            out = np.vstack( [out, np.array(Xp[str(i)])] )
        except:
            raise
    np.savetxt(os.path.join(Prefix, OutputFile + ".txt"), out.T )
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="TrJm", usage='%(prog)s [options] -f file',
                                     description="Chenyen's personal command to make life easier.",
                                     epilog="All's well that ends well.")
    parser.add_argument('-l', '--length', type=int, default=3, help='L.')
    parser.add_argument('-i', '--input', default='QuenchState-Z-0.h5', help='Filename.')
    parser.add_argument('-o', '--output', default='Jtot', help='Filename Prefix.')
    parser.add_argument('--from', type=int, default=0, help='subfolder sequence.')
    parser.add_argument('--to', type=int, default=1000, help='subfolder sequence.')
    parser.add_argument('--every', type=int, default=20, help='load every.')
    args = vars(parser.parse_args())
    for i in range(args["from"], args["to"], 1):
        prefix = "S"+str(i).zfill(4)
        main(args["length"], args["input"], args["output"], prefix, args["every"])
