#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py

# Get all platform settings
from Clusters import *

# To generate scripts
import ScriptGenerator as sg

L = 12
N1 = 4
N2 = N1

Uinit = 0
# Uinit = 9
UVls = [(9,-3)]

Tsteps = np.arange(0, 2100, 100, dtype=np.int)
dt = 0.005

S2 = 0# no flip spin

LN1N2 = "-".join(["rixs", "".join(["L", str(L)]), str(N1), str(N2)])
DATADIR = os.path.join(ExecDir, "ED", LN1N2)

for U, V in UVls:
  UVPrefix =  "-".join(["".join(["Ui", str(Uinit)]), "".join(["dU", str(U)]), "".join(["dV", str(V)])])

  for Tstep in Tsteps:
    TPrefix = "T" + str(Tstep)
    Job_Name = TPrefix + UVPrefix.replace("-", "")
    workdir = os.path.join(DATADIR, UVPrefix, TPrefix)
    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
    if Cluster == "Kagome":
      Filename = os.path.join(workdir, 'job')
      Executable = [" ".join([os.path.join(SrcDir, "build", "rixs"), str(L), str(N1), str(S2), str(Tstep), str(Uinit), str(V), str(U)]),]
      sg.GenerateScript("PBS", Filename, Job_Name, Executable, workdir, Nodes=1, NumCore=8, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "Merced":
      Filename = os.path.join(workdir, 'job')
      Executable = [" ".join([os.path.join(SrcDir, "build", "rixs"), str(L), str(N1), str(S2), str(Tstep), str(Uinit), str(V), str(U)]),]
      sg.GenerateScript("TORQUE", Filename, Job_Name, Executable, workdir, Nodes=1, NumCore=10, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "LANL":
      Filename = os.path.join(workdir, 'job.mpi')
      Executable = [" ".join(["mpirun", "-n", str(L), "-ppn", "1", os.path.join(SrcDir, "build", "rixs.mpi"), str(L), str(N1), str(S2), str(Tstep), str(Uinit), str(V), str(U)]),]
      sg.GenerateScript("SLURM", Filename, Job_Name, Executable, workdir, Nodes=196, NumCore=20, WallTime='16:00:00', Partition='standard', ProjectName='s17_cint', MPI=1, PPN=1)
