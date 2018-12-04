#!/usr/bin/env python
import numpy as np
import h5py
import argparse

L = 4

def main(filename, prefix):
    f = h5py.File(filename, "r")
    Nf1 = []
    Nf2 = []
    Nf3 = []
    Nf4 = []
    Jf1 = []
    Jf2 = []
    Jf3 = []
    Jf4 = []
    AvP = []
    cnt = 0
    gname = "Obs-" + str(cnt)
    e = gname in f.keys()
    while e:
        Fmn = f[gname]["Fermion"]["Elem"]["real"][:] + 1.0j * f[gname]["Fermion"]["Elem"]["imag"][:]
        Fmn = np.reshape(Fmn, (L,L))
        Nf1.append( Fmn[0,0].real )
        Nf2.append( Fmn[1,1].real )
        Nf3.append( Fmn[2,2].real )
        Nf4.append( Fmn[3,3].real )
        Jf1.append( Fmn[0,1].imag )
        Jf2.append( Fmn[1,2].imag )
        Jf3.append( Fmn[2,3].imag )
        Jf4.append( Fmn[3,0].imag )
        Pm = f[gname]["Phonon"]["Elem"]["real"][:]
        AvP.append( np.average(Pm) )
        cnt += 20
        gname = "Obs-" + str(cnt)
        e = gname in f.keys()
    Nf1 = np.array( Nf1 )
    Nf2 = np.array( Nf2 )
    Nf3 = np.array( Nf3 )
    Nf4 = np.array( Nf4 )
    Jf1 = np.array( Jf1 )
    Jf2 = np.array( Jf2 )
    Jf3 = np.array( Jf3 )
    Jf4 = np.array( Jf4 )
    AvP = np.array( AvP )
    np.savetxt(prefix + ".txt", np.vstack([Nf1, Nf2, Nf3, Nf4, Jf1, Jf2, Jf3, Jf4, AvP]).T )
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="TrJm", usage='%(prog)s [options] -f file',
                                     description="Chenyen's personal command to make life easier.",
                                     epilog="All's well that ends well.")
    parser.add_argument('-s', '--sites', default=3, help='How many sites?')
    parser.add_argument('-i', '--input', default='QuenchState-Z-0.h5', help='Filename.')
    parser.add_argument('-o', '--output', default='Jtot', help='Filename Prefix.')
    parser.add_argument('--from', type=int, default=0, help='subfolder sequence.')
    parser.add_argument('--to', type=int, default=1000, help='subfolder sequence.')
    args = vars(parser.parse_args())
    for i in range(args["from"], args["to"], 1):
        prefix = "S"+str(i).zfill(4)
        main(prefix+"/"+args["input"], prefix+"/"+args["output"])
