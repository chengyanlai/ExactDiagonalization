#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

L = 12
OBC = 1# 1:True
N1 = 6
N2 = N1
Us = [0, 0.01, 4, 8, 12]

OpSites = [0, ]

if OBC:
  Prefix = "-".join(["BFKO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  Prefix = "-".join(["BFKP", "".join(["L", str(L)]), str(N1), str(N2)])
DataDir = os.path.join(ExecDir, "ED", Prefix)

APPs = []
APPs.append(os.path.join(SrcDir, "build", "fhm.1d 0"))
APPs.append(os.path.join(SrcDir, "build", "fhm.1d 2 1"))
APPs.append(os.path.join(SrcDir, "build", "fhm.1d 2 2"))
APPs.append(os.path.join(SrcDir, "build", "fhm.1d 4 1"))
APPs.append(os.path.join(SrcDir, "build", "fhm.1d 4 2"))
APPs.append("/bin/touch DONE")

for Uinit in Us:
  for OpSite in OpSites:
    Pump = "NP"
    JobName =  "-".join([ "".join(["Ui", str(Uinit)]), Pump, str(OpSite) ])

    workdir = os.path.join(DataDir, JobName)

    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2

    f = h5py.File(os.path.join(workdir, 'conf.h5'), 'w')
    para = f.create_group("Parameters")
    dset = para.create_dataset("L", data=L)
    dset = para.create_dataset("OBC", data=OBC)
    dset = para.create_dataset("N1", data=N1)
    dset = para.create_dataset("N2", data=N2)
    if OBC: Jls = np.ones(L-1, dtype=np.float64)
    else: Jls = np.ones(L, dtype=np.float64)
    dset = para.create_dataset("J", data=Jls)
    Uls = np.zeros(L, dtype=np.float64)
    Uls[OpSite] = Uinit
    dset = para.create_dataset("U", data=Uls)
    Vls = np.zeros(L, dtype=np.float64)
    Vls[OpSite] = -Uinit / 2.
    dset = para.create_dataset("V", data=Vls)
    Wls = np.zeros(L, dtype=np.float64)
    dset = para.create_dataset("W", data=Wls)
    # DOS
    dset = para.create_dataset("Uf", data=Uls)
    dset = para.create_dataset("Vf", data=Vls)
    dset = para.create_dataset("TStepsX", data=5000)
    dset = para.create_dataset("dtX", data=0.01)
    dset = para.create_dataset("CoreHole", data=OpSite)
    dset = para.create_dataset("Species", data=0)
    f.close()

    Filename = os.path.join(workdir, 'job')
    if Cluster == "Kagome":
      sg.GenerateScript("PBS", Filename, JobName, APPs, workdir, Nodes=1, NumCore=12, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "Merced":
      sg.GenerateScript("TORQUE", Filename, JobName, APPs, workdir, Nodes=1, NumCore=2, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "LANL":
      sg.GenerateScript("SLURM", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='16:00:00', Partition='standard', ProjectName='s17_cint')
