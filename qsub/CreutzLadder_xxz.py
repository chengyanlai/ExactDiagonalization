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

L = 6
Nup = 3# Ndn = L - Nup
# OBC = 1# 1:True
OBC = 0# 0:False (PBC)

JAA = 1.0
JBB = 1.0 * JAA
JvAB = 0.0 * JAA
JdAB = 1.0 * JAA
JdBA = 1.0 * JAA
Delta = 0.0 * JAA
Phi = 0.25# will multiply pi later on

dynamics = 0
Tsteps = 10000
dt = 0.005

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "CreutzLadder.xxz"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["CLxxzO", "".join(["L", str(L)]), str(Nup)])
else:
  LN1N2 = "-".join(["CLxxzP", "".join(["L", str(L)]), str(Nup)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

Job_Name = "".join(["PHI", str(Phi)])

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
    dset = para.create_dataset("JAA", data=JAA)
    dset = para.create_dataset("JBB", data=JBB)
    dset = para.create_dataset("JvAB", data=JvAB)
    dset = para.create_dataset("JdAB", data=JdAB)
    dset = para.create_dataset("JdBA", data=JdBA)
    dset = para.create_dataset("Delta", data=Delta)
    dset = para.create_dataset("Phi", data=Phi*np.pi)
    dset = para.create_dataset("OBC", data=OBC)
    dset = para.create_dataset("N", data=Nup)
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
