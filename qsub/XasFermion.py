#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import Script_Helpers as shp
from Clusters import *

# Get all platform settings
from Clusters import *

NumThreads = 12
WallTime = '240:0:0'
qsub_cmd = "qsub -q short.q"

L = 14
OBC = 1# 1:True
N1 = 7
N2 = N1
CHloc = np.int(L / 2)

Uinit = 0
Uinit = 9
UVls = [(40., -10.)]
# UVls = [
#         (5.0,  -1.), (10.0,  -1.), (15.0,  -1.), (20.0,  -1.), (25.0,  -1.), (30.0,  -1.), (40.0,  -1.),
#         (5.0,  -5.), (10.0,  -5.), (15.0,  -5.), (20.0,  -5.), (25.0,  -5.), (30.0,  -5.), (40.0,  -5.),
#         (5.0, -10.), (10.0, -10.), (15.0, -10.), (20.0, -10.), (25.0, -10.), (30.0, -10.), (40.0, -10.),
#         (5.0, -15.), (10.0, -15.), (15.0, -15.), (20.0, -15.), (25.0, -15.), (30.0, -15.), (40.0, -15.),
#         (5.0, -20.), (10.0, -20.), (15.0, -20.), (20.0, -20.), (25.0, -20.), (30.0, -20.), (40.0, -20.),
#        ]

dynamics = 1
Tsteps = 2000
dt = 0.005

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "xas.f"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["xasO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  LN1N2 = "-".join(["xasP", "".join(["L", str(L)]), str(N1), str(N2)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

for U, V in UVls:
  Job_Name =  "-".join(["".join(["Ui", str(Uinit)]), "".join(["U", str(U)]), "".join(["V", str(V)])])

  workdir = os.path.join(DATADIR, Job_Name)

  os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
  with shp.cd(workdir):
    if os.path.isfile('XASDYN.h5'):
      print("".join([workdir, " already done!"]))
      pass
    elif os.path.isfile('conf.h5'):
      print("".join([workdir, " is schaduled!?"]))
      pass
    else:
      f = h5py.File('conf.h5', 'w')
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
