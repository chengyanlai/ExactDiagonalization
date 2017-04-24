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

NumThreads = 10
WallTime = '24:0:0'

L = 14
OBC = 1# 1:True
N1 = 6
N2 = 6
CHloc = np.int(L / 2)

Uin = [0.5, 2.5, 5.0, 10.0]
Vin = [-2.5, -5.0, -10.0, -20.0, -40.0]

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

for V in Vin:
  for U in Uin:
    Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["V", str(V)])])

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
        Uls = np.zeros(L, dtype=np.float64)
        Uls[CHloc] = U
        dset = para.create_dataset("U", data=Uls)
        Vls = np.zeros(L, dtype=np.float64)
        Vls[CHloc] = V
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
