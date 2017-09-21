#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import subprocess
import platform
import socket
import numpy as np
import h5py
import Script_Helpers as shp
from Clusters import *

BL = 10
FL = 1
maxLocalB = 1
if BL == 1:
  Jbbs = [0.0,]
else:
  Jbbs = [0.01, 0.05, 0.10]
if FL == 1:
  Jffs = [0.0,]
  Uffs = [0.00]
else:
  Jffs = [0.01, 0.05, 0.10]
  Uffs = [0.00, 0.10, 0.50]
Vbbs = [3.20, 3.50, 3.80]
Vffs = [3.20, 3.40, 3.60]
DeltaDCs = [-0.0080, -0.020, -0.050, -0.080]

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "plex.mpi"))
APPs.append("/bin/touch DONE")
Exac_program = "\n".join(APPs)

QSUB = False
NumCores = 10
WallTime = MaxWallTime

DATADIR = os.path.join( EXEC_DIR, "plex", "".join(["B", str(BL), "F", str(FL), "mB", str(maxLocalB)]) )
for Uff in Uffs:
  for Jbb in Jbbs:
    for Jff in Jffs:
      Job_Name = "".join(["PlExJB", str(Jbb), "JF", str(Jff), "Uf", str(Uff)])
      workdir = os.path.join(DATADIR, Job_Name)
      os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
      with shp.cd(workdir):
        if os.path.isfile('DONE'):
          print(workdir, " is DONE!")
          continue
        SetCount = 0
        for Vbb in Vbbs:
          for Vff in Vffs:
            for DeltaDC in DeltaDCs:
              f = h5py.File('confs.h5', 'a')
              para = f.create_group("Input-" + str(SetCount))
              dset = para.create_dataset("BL", data=BL)
              dset = para.create_dataset("FL", data=FL)
              dset = para.create_dataset("maxLocalB", data=maxLocalB)
              dset = para.create_dataset("Jbb", data=Jbb)
              dset = para.create_dataset("Jff", data=Jff)
              dset = para.create_dataset("Vbb", data=Vbb)
              dset = para.create_dataset("Vff", data=Vff)
              dset = para.create_dataset("Uff", data=Uff)
              dset = para.create_dataset("DeltaDC", data=DeltaDC)
              f.close()
              SetCount += 1

        shp.WriteQsubSBATCH_MPI("job", Job_Name, Exac_program, workdir, Nodes=SetCount, \
                                NumCore=NumCores, WallTime=WallTime, partition=Partition)

