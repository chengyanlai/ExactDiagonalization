#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

# Get all platform settings
from Clusters import *

L = 14
OBC = 1# 1:True
N1 = 6
N2 = N1
CHloc = np.int(L / 2)

Uinit = 0
UVls = [(9., -3.)]

dynamics = 1
Tsteps = 2000
dt = 0.005

APPs = []
APPs.append(os.path.join(SrcDir, "build", "xas.f"))

if OBC:
  LN1N2 = "-".join(["xasO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  LN1N2 = "-".join(["xasP", "".join(["L", str(L)]), str(N1), str(N2)])
DATADIR = os.path.join(ExecDir, "ED", LN1N2)

for U, V in UVls:
  JobName =  "-".join(["".join(["Ui", str(Uinit)]), "".join(["U", str(U)]), "".join(["V", str(V)])])

  workdir = os.path.join(DATADIR, JobName)

  os.makedirs(workdir, exist_ok=True)  # Python >= 3.2

  f = h5py.File(os.path.join(workdir, 'conf.h5'), 'w')
  para = f.create_group("Parameters")
  dset = para.create_dataset("L", data=L)
  dset = para.create_dataset("OBC", data=OBC)
  dset = para.create_dataset("N1", data=N1)
  dset = para.create_dataset("N2", data=N2)
  dset = para.create_dataset("CHloc", data=CHloc)
  Uls = np.ones(L, dtype=np.float64) * np.float64(Uinit)
  dset = para.create_dataset("Uinit", data=Uls)
  Uls[CHloc] += U
  dset = para.create_dataset("U", data=Uls)
  Vls = np.zeros(L, dtype=np.float64)
  Vls[CHloc] += V
  dset = para.create_dataset("V", data=Vls)
  dset = para.create_dataset("dynamics", data=dynamics)
  dset = para.create_dataset("Tsteps", data=Tsteps)
  dset = para.create_dataset("dt", data=dt)
  f.close()

  sg.GenerateScript("SLURM", os.path.join(workdir, 'job'), JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='16:00:00', Partition='standard', ProjectName='s17_cint')
