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

L = 12
OBC = 0# 1:True
N1 = 6
N2 = N1
Uinit = 6
Winit = 3

# ARPES Parameters
OpSites = range(L)
# Type = 1# absorption
Type = -1# emission
Tsteps = 4000
dt = 0.005

# For pumping pulse
def GetAt(tlist, td, tau=2, W=3, A0=1, phi=0):
    p = []
    for t in tlist:
        val = A0 * np.exp( -(t - td) * (t - td) / (2. * tau * tau) ) * np.cos(W * (t - td) + np.pi * phi)
        p.append(val)
    return np.array(p)

A0 = 0.1
Tau = 1
Omega = 4
Td = 4
Phi = 0
if A0:
  tl = np.arange(0, 2*Td+dt, dt)
  At = GetAt(tl, td=Td, tau=Tau, W=Omega, A0=A0, phi=Phi)
else:
  At = np.array([])

APPs = []
if A0:
  APPs.append(os.path.join(SrcDir, "build", "fhm.1d 2"))
  APPs.append(os.path.join(SrcDir, "build", "fhm.1d 3"))
else:
  APPs.append(os.path.join(SrcDir, "build", "fhm.1d 20"))
APPs.append("/bin/touch DONE")

if OBC:
  Prefix1 = "-".join(["trArpesO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  Prefix1 = "-".join(["trArpesP", "".join(["L", str(L)]), str(N1), str(N2)])

Prefix2 = "".join(["Ui", str(Uinit)])

DataDir = os.path.join(ExecDir, "ED", Prefix1, Prefix2)

for OpSite in OpSites:
  if A0:
    Pump = "".join(["T",str(Tau), "W", str(Omega), "A", str(A0)])
  else:
    Pump = "NP"
  OPSite = "".join(["OP", str(OpSite), "T", str(Type)])
  workdir = os.path.join(DataDir, Pump, OPSite)
  os.makedirs(workdir, exist_ok=True)  # Python >= 3.2

  JobName =  "-".join([Prefix2, OPSite, Pump ])
  f = h5py.File(os.path.join(workdir, 'conf.h5'), 'w')
  para = f.create_group("Parameters")
  dset = para.create_dataset("L", data=L)
  dset = para.create_dataset("OBC", data=OBC)
  dset = para.create_dataset("N1", data=N1)
  dset = para.create_dataset("N2", data=N2)
  if OBC: Jls = np.ones(L-1, dtype=np.float64)
  else: Jls = np.ones(L, dtype=np.float64)
  dset = para.create_dataset("J", data=Jls)
  Uls = np.ones(L, dtype=np.float64) * np.float64(Uinit)
  dset = para.create_dataset("U", data=Uls)
  Vls = np.zeros(L, dtype=np.float64)
  dset = para.create_dataset("V", data=Vls)
  Wls = np.ones(L, dtype=np.float64) * np.float64(Winit)
  dset = para.create_dataset("W", data=Wls)
  # Pump
  dset = para.create_dataset("At", data=At)
  dset = para.create_dataset("TStepsA", data=At.shape[0])
  dset = para.create_dataset("dtA", data=dt)
  # XAS
  dset = para.create_dataset("Uf", data=Uls)
  dset = para.create_dataset("Vf", data=Vls)
  dset = para.create_dataset("TStepsX", data=Tsteps)
  dset = para.create_dataset("dtX", data=dt)
  dset = para.create_dataset("CoreHole", data=OpSite)
  dset = para.create_dataset("Type", data=Type)
  dset = para.create_dataset("Species", data=0)
  f.close()

  Filename = os.path.join(workdir, 'job')
  if Cluster == "Kagome":
    sg.GenerateScript("PBS", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
  elif Cluster == "Merced":
    sg.GenerateScript("TORQUE", Filename, JobName, APPs, workdir, Nodes=1, NumCore=2, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
  elif Cluster == "LANL":
    sg.GenerateScript("SLURM", Filename, JobName, APPs, workdir, Nodes=1, NumCore=2, WallTime='16:00:00', Partition='standard', ProjectName='w18_xasrixs')
