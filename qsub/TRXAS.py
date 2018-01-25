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
N1 = 7
N2 = N1
Uinit = 9

# XAS Parameters
CHloc = np.int(L / 2)
UVls = [(0., -3.)]
Tsteps = 4000
dt = 0.005

# For pumping pulse
A0 = 1# Amplitude
Tau = 2#
W0 = 3#
Td = np.int(Tau * np.rint(np.sqrt(2. * np.log(100 * A0))) )
tl = np.arange(0, 2*Td, dt)
def getAt(tlist, td, tau=12, W=3, A0=1):
  p = []
  for t in tlist:
      val = A0 * np.exp( -(t - td) * (t - td) / (2. * tau * tau) ) * np.cos(W * (t - td))
      p.append(val)
  return np.array(p)
At = getAt(tl, td=Td, tau=Tau, W=W0, A0=A0)

APPs = []
APPs.append(os.path.join(SrcDir, "build", "trxas.f"))

if OBC:
  Prefix = "-".join(["trxasO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  Prefix = "-".join(["trxasP", "".join(["L", str(L)]), str(N1), str(N2)])
DataDir = os.path.join(ExecDir, "ED", Prefix)

for U, V in UVls:
  JobName =  "-".join(["".join(["Ui", str(Uinit)]), "".join(["U", str(U), "V", str(V)]), "".join(["T",str(Tau), "W", str(W0), "A", str(A0)]) ])

  workdir = os.path.join(DataDir, JobName)

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
  dset = para.create_dataset("At", data=At)
  dset = para.create_dataset("Tsteps", data=Tsteps)
  dset = para.create_dataset("dt", data=dt)
  f.close()

  Filename = os.path.join(workdir, 'job')
  if Cluster == "Kagome":
    sg.GenerateScript("PBS", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
  elif Cluster == "Merced":
    sg.GenerateScript("TORQUE", Filename, JobName, APPs, workdir, Nodes=1, NumCore=20, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
  elif Cluster == "LANL":
    sg.GenerateScript("SLURM", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='16:00:00', Partition='standard', ProjectName='s17_cint')
