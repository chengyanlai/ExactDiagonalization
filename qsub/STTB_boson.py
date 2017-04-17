#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import subprocess
import platform
import socket
import numpy as np
import h5py
import Script_Helpers as shp
from Clusters import *

NumThreads = 12

QSUB = True

L = 8
N = 8
jaa = 1.0
jab = np.sqrt(2.)
TBloc = 0
# Uin = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0,]
Uin = [0.1, 5.0]
# GammaList = [0, 5., 8., 10., 20.]
GammaList = [1., 3., 5., 7., 9.]
# GammaList = [0,]

Vin = 0.0
# NOTE: Dynamics parameters
Tsteps = 6000# Tstep * dt is final time
dt = 0.005

for U in Uin:
    for Gamma in GammaList:
        APPs = []
        if Gamma:
            GSfolder = "-".join(["STP", "".join(["J", str(jaa)]), "".join(["L", str(L)]), str(N)])
            OldJobName =  "-".join(["".join(["U", str(U)]), "".join(["V", str(Vin)])])
            GSfile = os.path.join(EXEC_DIR, GSfolder, OldJobName, "GS.h5")
            LN1N2 = "-".join(["STTB", "".join(["J", str(jaa)]), "".join(["L", str(L)]), str(N)])
            Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["V", str(Vin)]),
                                  "".join(["G", str(Gamma)]), "".join(["TB", str(TBloc)])] )
            APPs.append(os.path.join(SRC_DIR, "build", "STTB.b"))
            Exac_program = "\n".join(APPs)
        else:
            LN1N2 = "-".join(["STP", "".join(["J", str(jaa)]), "".join(["L", str(L)]), str(N)])
            Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["V", str(Vin)])])
            APPs.append(os.path.join(SRC_DIR, "build", "ST.b"))
            Exac_program = "\n".join(APPs)
        DATADIR = os.path.join(EXEC_DIR, LN1N2)
        workdir = os.path.join(DATADIR, Job_Name)

        os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
        if Gamma:
            CopyGSCmd = " ".join(["cp", GSfile, workdir])
            subprocess.call(CopyGSCmd, shell=True)
        with shp.cd(workdir):
            if os.path.isfile('STTBb.h5'):
                pass
            else:
                f = h5py.File('conf.h5', 'w')
                para = f.create_group("Parameters")
                dset = para.create_dataset("L", data=L)
                dset = para.create_dataset("N", data=N)
                dset = para.create_dataset("jab", data=jab)
                dset = para.create_dataset("jaa", data=jaa)
                dset = para.create_dataset("U", data=U)
                dset = para.create_dataset("V", data=Vin)
                dset = para.create_dataset("Tsteps", data=Tsteps)
                dset = para.create_dataset("dt", data=dt)
                dset = para.create_dataset("TBloc", data=TBloc)
                dset = para.create_dataset("Gamma", data=Gamma)
                f.close()

                if socket.gethostname() == 'kagome.ucmerced.edu' or \
                  socket.gethostname() == 'edgestate.rcc.ucmerced.edu' or \
                  socket.gethostname() == 'atomtronics.ucmerced.edu':
                    shp.WriteQsubPBS("qsub.s", Job_Name, Exac_program, workdir,
                                     Nodes="1", NumCore=NumThreads)
                if platform.system() == "Darwin" and QSUB:
                    for i in APPs:
                        print(i)
                        f = open(Job_Name, "w")
                        subprocess.call(i, shell=True, stdout=f)
                        f.close()
                else:
                    qsub_script = " ".join([qsub_cmd, "qsub.s"])
                    if QSUB: subprocess.call(qsub_script, shell=True)
