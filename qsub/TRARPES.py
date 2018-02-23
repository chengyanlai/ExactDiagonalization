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

# ARPES Parameters
OpSites = range(L)
Tsteps = 4000
dt = 0.005

# For pumping pulse
A0 = 1
Tau = 2
W0 = 3
def getAt(tlist, td, tau=12, W=3, A0=1):
  p = []
  for t in tlist:
      val = A0 * np.exp( -(t - td) * (t - td) / (2. * tau * tau) ) * np.cos(W * (t - td))
      p.append(val)
  return np.array(p)
if A0:
  Td = np.int(Tau * np.rint(np.sqrt(2. * np.log(100 * A0))) )
  tl = np.arange(0, 2*Td, dt)
  At = getAt(tl, td=Td, tau=Tau, W=W0, A0=A0)
else:
  At = np.array([])

APPs = []
if A0: APPs.append(os.path.join(SrcDir, "build", "fhm.1d 210"))
else: APPs.append(os.path.join(SrcDir, "build", "fhm.1d 20"))
APPs.append("/bin/touch DONE")

if OBC:
  Prefix = "-".join(["trArpesO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  Prefix = "-".join(["trArpesP", "".join(["L", str(L)]), str(N1), str(N2)])
DataDir = os.path.join(ExecDir, "ED", Prefix)

for OpSite in OpSites:
  if A0:
    Pump = "".join(["T",str(Tau), "W", str(W0), "A", str(A0)])
  else:
    Pump = "NP"
  JobName =  "-".join(["".join(["Ui", str(Uinit)]), "".join(["OP", str(OpSite)]), Pump ])

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
  Uls = np.ones(L, dtype=np.float64) * np.float64(Uinit)
  dset = para.create_dataset("U", data=Uls)
  Vls = np.zeros(L, dtype=np.float64)
  dset = para.create_dataset("V", data=Vls)
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
  dset = para.create_dataset("Type", data=1)
  dset = para.create_dataset("Species", data=0)
  f.close()

  Filename = os.path.join(workdir, 'job')
  if Cluster == "Kagome":
    sg.GenerateScript("PBS", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
  elif Cluster == "Merced":
    sg.GenerateScript("TORQUE", Filename, JobName, APPs, workdir, Nodes=1, NumCore=20, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
  elif Cluster == "LANL":
    sg.GenerateScript("SLURM", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='16:00:00', Partition='standard', ProjectName='s17_cint')
