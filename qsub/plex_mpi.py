#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import subprocess
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

BL = 1
FL = 5
maxLocalB = 1
if BL == 1:
  Jbbs = [0.0,]
else:
  Jbbs = [0.0, 0.05, 0.10]
if FL == 1:
  Jffs = [0.0,]
  Uffs = [0.00]
else:
  Jffs = [0.0, 0.05, 0.10]
  Uffs = [0.0,-3.0,-4.0, 0.10, 0.50]
  # Uffs = [-4.1,]
Vbbs = [3.20, 3.50, 3.80]
Vffs = [3.20, 3.40, 3.60]
DeltaDCs = [-0.0080, -0.020, -0.050, -0.080]
CouplingForm = "uniform"# uniform, angle1
# CouplingForm = "angle1"
if CouplingForm == "angle1":
  if BL > FL:
    CouplingForm = "P2"
  if BL < FL:
    CouplingForm = "COS2"

# Dynamics
dT = 0.005
Tf = 3000
tlist = np.arange(0, 12, dT)
phi = 0
def pulse(tlist, td=6, tau=2, W=1, A0=1, Phase=0.0e0):
  p = []
  for t in tlist:
    val = A0 * np.exp( -(t - td) * (t - td) / (2. * tau * tau) ) * np.cos(W * (t - td) + Phase)
    p.append(val)
  return np.array(p)

NumCore = np.int(MaxNumThreads / 3)
WallTime = MaxWallTime

def DCcoupling(Delta, BL, FL, fi=0, form="uniform"):
  if form == "uniform":
    return Delta * np.ones(BL, dtype=np.float64)
  elif form == "P2" and BL > FL:# Bosons surround Fermion
    angles = np.linspace(0, 2.*np.pi, BL+1)
    work = -1.0 * Delta * np.polynomial.legendre.legval(np.cos(angles), [0,0,1])
    return work[:-1]
  elif form == "COS2" and BL < FL:# Fermion surround Bosons
    angles = np.linspace(0, 2.*np.pi, FL+1)
    work = -1.0 * Delta * np.cos(angles[fi]) * np.cos(angles[fi]) * np.ones(BL, dtype=np.float64)
    return work
  else:
    print("Not supported yet")
    sys.exit()

DataDir = os.path.join( ExecDir, "plex", "".join(["B", str(BL), "F", str(FL), "mB", str(maxLocalB)]), CouplingForm)
Paths = []
SetCount = 0
Job_Name = "".join(["PlExB", str(BL), "F", str(FL), "mB", str(maxLocalB), CouplingForm])
for Uff in Uffs:
  Prefix0 = "".join([ "Ue", str(Uff) ])
  for Jbb in Jbbs:
    for Jff in Jffs:
      Prefix1 = "".join([ "Jd", str(Jbb), "Je", str(Jff) ])
      for Vbb in Vbbs:
        for Vff in Vffs:
          for DeltaDC in DeltaDCs:
            Prefix2 = "".join(["Vd", str(Vbb), "Ve", str(Vff), "Ddc", str(DeltaDC)])
            workdir = os.path.join(DataDir, Prefix0, Prefix1, Prefix2)
            os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
            f = h5py.File(os.path.join(workdir, 'confs.h5'), 'w')
            para = f.create_group("Input")
            dset = para.create_dataset("BL", data=BL)
            dset = para.create_dataset("FL", data=FL)
            dset = para.create_dataset("maxLocalB", data=maxLocalB)
            dset = para.create_dataset("Jbb", data=Jbb)
            dset = para.create_dataset("Jff", data=Jff)
            dset = para.create_dataset("Vbb", data=Vbb)
            dset = para.create_dataset("Vff", data=Vff)
            dset = para.create_dataset("Uff", data=Uff)
            for i in range(FL):
              DeltaDCarr = DCcoupling(DeltaDC, BL, FL, fi=i, form=CouplingForm)
              gname = "DeltaDC-" + str(i)
              dset = para.create_dataset(gname, data=DeltaDCarr)
            SetCount += 1
            Paths.append(workdir)
            f.close()
            fp = h5py.File(os.path.join(workdir, 'pulse.h5'), 'w')
            para = fp.create_group("Input")
            Vfft = pulse(tlist, Phase=phi*np.pi)
            dset = para.create_dataset("Phase", data=phi)
            dset = para.create_dataset("Vfft", data=Vfft)
            dset = para.create_dataset("Tf", data=Tf+Vfft.shape[0])
            dset = para.create_dataset("dT", data=dT)
            fp.close()

MPIFolders = open(os.path.join(DataDir, "MPIFolders"), "w")
for item in Paths:
  MPIFolders.write("%s\n" % item)
APPs = []
APPs.append("mpirun -n " + str(SetCount) + " -ppn 3 " + os.path.join(SrcDir, "build", "plex.mpi 1"))
APPs.append("/bin/touch DONE")
Exac_program = "\n".join(APPs)
sg.GenerateScript("SLURM", os.path.join(DataDir, "job.mpi"), Job_Name, APPs, DataDir, np.int(SetCount/3), NumCore, WallTime, Partition, Project, MPI=SetCount, PPN=3)
