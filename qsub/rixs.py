#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
from Clusters import *

# Get all platform settings
from Clusters import *

L = 14
N1 = 7
N2 = N1

Uinit = 0
# Uinit = 9
UVls = [(9,-3)]

Tsteps = np.arange(0, 2100, 100, dtype=np.int)
dt = 0.005

S2 = 0# no flip spin

APPs = []
APPs.append(os.path.join(SrcDir, "build", "rixs.mps"))
Exac_program = "\n".join(APPs)

LN1N2 = "-".join(["rixs", "".join(["L", str(L)]), str(N1), str(N2)])
DATADIR = os.path.join(ExecDir, "ED", LN1N2)

for U, V in UVls:
  Job_Name =  "-".join(["".join(["Ui", str(Uinit)]), "".join(["dU", str(U)]), "".join(["dV", str(V)])])

  for Tstep in Tsteps:
    TPrefix = "T" + str(Tstep)

    workdir = os.path.join(DATADIR, Job_Name, TPrefix)

    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
    f = h5py.File(os.path.join(workdir, 'conf.h5'), 'w')
    para = f.create_group("Parameters")
    dset = para.create_dataset("L", data=L)
    dset = para.create_dataset("N1", data=N1)
    dset = para.create_dataset("N2", data=N2)
    dset = para.create_dataset("Uinit", data=Uinit)
    dset = para.create_dataset("U", data=U)
    dset = para.create_dataset("V", data=V)
    dset = para.create_dataset("Tsteps", data=Tstep)
    dset = para.create_dataset("dt", data=dt)
    dset = para.create_dataset("S2", data=S2)
    f.close()
