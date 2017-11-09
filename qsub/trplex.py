#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import subprocess
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

BL = 2
FL = 1
maxLocalB = 1
if BL == 1:
  Jds = [0.0,]
else:
  Jds = [0.0, 0.05, 0.10]
if FL == 1:
  Jcs = [0.0,]
  Vcs = [0.00]
else:
  Jcs = [0.0, 0.05, 0.10]
  Vcs = [0.0,]#-3.0,-4.0, 0.10, 0.50]
  # Vcs = [-4.1,]
Eds = [3.50, ]#3.20, 3.80]
Ecs = [3.40, 3.20, 3.60]
DeltaDCs = [-0.080, ]#-0.020, -0.050, -0.0080]
CouplingForm = "uniform"# uniform, angle1
# CouplingForm = "angle1"
if CouplingForm == "angle1":
  if BL > FL:
    CouplingForm = "P2"
  if BL < FL:
    CouplingForm = "COS2"

# Dynamics
dT = 0.005
Tf = 100000
# pulse
Phis = [0, 0.25, 0.5]
A0s = [0.5, 1.5, 3.0]
def pulse(dt, tau=2, W=1, A0=2.5, Phase=0):
  td = np.int(tau * np.rint(np.sqrt(2. * np.log(100 * A0))) )
  tlist = np.arange(0, 2*td, dt)
  p = []
  for t in tlist:
    val = A0 * np.exp( -(t - td) * (t - td) / (2. * tau * tau) ) * np.cos(W * (t - td) + np.pi * Phase)
    p.append(val)
  return np.array(p)
  # return np.zeros(len(tlist), dtype=np.float64)

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

DataDir = os.path.join( ExecDir, "trplex", "".join(["B", str(BL), "F", str(FL), "mB", str(maxLocalB)]), CouplingForm)
Paths = []
SetCount = 0
Job_Name = "".join(["PlExB", str(BL), "F", str(FL), "mB", str(maxLocalB), CouplingForm])
for Vc in Vcs:
  for Jd in Jds:
    for Jc in Jcs:
      Prefix1 = "".join([ "Jd", str(Jd), "Jc", str(Jc), "Vc", str(Vc) ])
      for Ed in Eds:
        for Ec in Ecs:
          for DeltaDC in DeltaDCs:
            Prefix2 = "".join([ "Ed", str(Ed), "Ec", str(Ec), "Ddc", str(DeltaDC) ])
            for A0 in A0s:
              for phi in Phis:
                Prefix3 = "".join(["A", str(A0), "P", str(phi)])
                workdir = os.path.join(DataDir, Prefix1, Prefix2, Prefix3)
                os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
                f = h5py.File(os.path.join(workdir, 'confs.h5'), 'w')
                para = f.create_group("Input")
                dset = para.create_dataset("BL", data=BL)
                dset = para.create_dataset("FL", data=FL)
                dset = para.create_dataset("maxLocalB", data=maxLocalB)
                dset = para.create_dataset("Jd", data=Jd)
                dset = para.create_dataset("Jc", data=Jc)
                dset = para.create_dataset("Ed", data=Ed)
                dset = para.create_dataset("Ec", data=Ec)
                dset = para.create_dataset("Vc", data=Vc)
                for i in range(FL):
                  DeltaDCarr = DCcoupling(DeltaDC, BL, FL, fi=i, form=CouplingForm)
                  gname = "DeltaDC-" + str(i)
                  dset = para.create_dataset(gname, data=DeltaDCarr)
                f.close()

                fp = h5py.File(os.path.join(workdir, 'pulse.h5'), 'w')
                para = fp.create_group("Input")
                Ect = pulse(dt=dT, A0=A0, Phase=phi) + Ec
                dset = para.create_dataset("Phase", data=phi)
                dset = para.create_dataset("Ect", data=Ect)
                dset = para.create_dataset("Tf", data=Tf+Ect.shape[0])
                dset = para.create_dataset("dT", data=dT)
                fp.close()

                SetCount += 1
                Paths.append(workdir)

MPIFolders = open(os.path.join(DataDir, "MPIFolders"), "w")
for item in Paths:
  MPIFolders.write("%s\n" % item)
APPs = []
APPs.append("mpirun -n " + str(SetCount) + " -ppn 3 " + os.path.join(SrcDir, "build", "plex.mpi 1"))
APPs.append("/bin/touch DONE")
Exac_program = "\n".join(APPs)
sg.GenerateScript("SLURM", os.path.join(DataDir, "job.mpi"), Job_Name, APPs, DataDir, np.int(SetCount/3), NumCore, WallTime, Partition, Project, MPI=SetCount, PPN=3)
