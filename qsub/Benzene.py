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

L = 6
OBC = 0# 1:True
N1 = 3
N2 = N1

J1 = 1.0
J2 = 0.5
Phases = np.linspace(-1, 1, 101)
Uin = 0

dynamics = 0
Tsteps = 0
dt = 0.005

if OBC:
  LN1N2 = "-".join(["FO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  LN1N2 = "-".join(["Benzene", "".join(["L", str(L)]), str(N1), str(N2)])
DATADIR = os.path.join(ExecDir, "ED", LN1N2)

for Phase in Phases:
  JobName =  "".join(["Phi", str(Phase)])

  workdir = os.path.join(DATADIR, "".join(["Uin", str(Uin)]), JobName )

  os.makedirs(workdir, exist_ok=True)  # Python >= 3.2

  f = h5py.File(os.path.join(workdir, 'conf.h5'), 'w')
  para = f.create_group("Parameters")
  dset = para.create_dataset("L", data=L)
  dset = para.create_dataset("OBC", data=OBC)
  dset = para.create_dataset("N1", data=N1)
  dset = para.create_dataset("N2", data=N2)
  dset = para.create_dataset("J1", data=J1)
  dset = para.create_dataset("J2", data=J2)
  dset = para.create_dataset("Phase", data=Phase)
  Uls = np.ones(L, dtype=np.float64) * np.float64(Uin)
  dset = para.create_dataset("U", data=Uls)
  Vls = np.zeros(L, dtype=np.float64)
  dset = para.create_dataset("V", data=Vls)
  dset = para.create_dataset("dynamics", data=dynamics)
  dset = para.create_dataset("Tsteps", data=Tsteps)
  dset = para.create_dataset("dt", data=dt)
  f.close()

  if Cluster == "Kagome":
    Filename = os.path.join(workdir, 'job')
    Executable = [" ".join([os.path.join(SrcDir, "build", "1D.f"), str(L),]), ]
    sg.GenerateScript("PBS", Filename, JobName, Executable, workdir, Nodes=1, NumCore=8, WallTime='1:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
  elif Cluster == "Merced":
    Filename = os.path.join(workdir, 'job')
    Executable = [" ".join([os.path.join(SrcDir, "build", "1D.f"), str(L),]), ]
    sg.GenerateScript("TORQUE", Filename, JobName, Executable, workdir, Nodes=1, NumCore=1, WallTime='1:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
