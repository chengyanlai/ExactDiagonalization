#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import Script_Helpers as shp

# Get all platform settings
from Clusters import *

def SetV(L, vtype="Box"):
  if vtype == "Box":
    V = np.zeros(L, dtype=np.float64)
  else:
    print("Vtype is not existed!")
    sys.exit()
  return V

NumThreads = 10
WallTime = '24:0:0'

L = 12
# J12ratio = np.linspace(0.1, 1.0, 10)
# J12ratio = np.linspace(0.92, 0.98,  4)
# J12ratio = [1.00, ]
J12ratio = [0.1, ]
# OBC = 1#1:True
OBC = 0# 0:False
if L % 2 == 1:
  N1 = np.int((L + 1) / 2)
  N2 = np.int((L - 1) / 2)
else:
  N1 = np.int(L / 2)
  N2 = np.int(L / 2)
# Uin = np.linspace(0.0, 10.0, 11)
# Uin = np.linspace(0.0, 10.0, 51)
Uin = [10.0,]
Phils = [0, ]
# Phils = np.linspace(0, L, 66)
Vtype = "Box"
Vin = SetV(L, vtype=Vtype)

dynamics = 1
Tsteps = 4000
dt = 0.005

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "SSH.f"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["FSSHO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  LN1N2 = "-".join(["FSSHP", "".join(["L", str(L)]), str(N1), str(N2)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

for nphi in Phils:
  phi = nphi * 2.0 * np.pi / np.float64(L)
  for J12 in J12ratio:
    for U in Uin:
      if OBC:
        Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["Jr", str(J12)]),
          Vtype])
      else:
        Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["Jr", str(J12)]),
          "".join(["n", str(nphi)]), Vtype])

      workdir = os.path.join(DATADIR, Job_Name)

      os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
      with shp.cd(workdir):
        if os.path.isfile('FSSH.h5'):
          print("".join([workdir, " already done!"]))
          pass
        elif os.path.isfile('conf.h5'):
          print("".join([workdir, " is schaduled!?"]))
          pass
        else:
          f = h5py.File('conf.h5', 'w')
          para = f.create_group("Parameters")
          dset = para.create_dataset("L", data=L)
          dset = para.create_dataset("J12", data=J12)
          dset = para.create_dataset("OBC", data=OBC)
          dset = para.create_dataset("N1", data=N1)
          dset = para.create_dataset("N2", data=N2)
          dset = para.create_dataset("U", data=U)
          dset = para.create_dataset("phi", data=phi)
          dset = para.create_dataset("V", data=Vin)
          dset = para.create_dataset("dynamics", data=dynamics)
          dset = para.create_dataset("Tsteps", data=Tsteps)
          dset = para.create_dataset("dt", data=dt)
          f.close()

          if Cluster == "Comet" or Cluster == "Stampede":
            shp.WriteQsubSBATCH("job", Job_Name, Exac_program, workdir, \
              NumCore=NumThreads, WallTime=WallTime, partition=Partition)
          elif Cluster == "Merced":
            shp.WriteQsubSGE("job", Job_Name, Exac_program, workdir, \
              NumCore=NumThreads, WallTime=WallTime)
          else:
            shp.WriteQsubPBS("job", Job_Name, Exac_program, workdir, \
              NumCore=NumThreads, WallTime=WallTime)

          if QSUB:
            qsub_script = " ".join([qsub_cmd, "job"])
            subprocess.call(qsub_script, shell=True)
