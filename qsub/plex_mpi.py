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
FL = 2
maxLocalB = 1
Jbbs = [0.10,]
Jffs = [0.10,]
Vbbs = [2.80, 3.20, 3.60]
Vffs = [3.00, 3.50, 4.00]
Uffs = [0.0, 0.10, 0.50,]
DeltaDCs = [-0.0010, -0.0050, -0.010, -0.050, -0.080]

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "plex.mpi"))
APPs.append("/bin/touch DONE")
Exac_program = "\n".join(APPs)

QSUB = False
NumCores = 10
WallTime = MaxWallTime

DATADIR = os.path.join( EXEC_DIR, "plex", "".join(["B", str(BL), "F", str(FL), "mB", str(maxLocalB)]) )
for Uff in Uffs:
  Job_Name = "".join(["PlExB", str(BL), "F", str(FL), "Uf", str(Uff)])
  workdir = os.path.join(DATADIR, Job_Name)
  os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
  with shp.cd(workdir):
    if os.path.isfile('DONE'):
      pass
    SetCount = 0
    for Jbb in Jbbs:
      for Jff in Jffs:
        for Vbb in Vbbs:
          for Vff in Vffs:
            for DeltaDC in DeltaDCs:
              f = h5py.File('confs.h5', 'a')
              para = f.create_group("Input-" + str(SetCount))
              print(SetCount)
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

    shp.WriteQsubSBATCH_MPI("job.mpi", Job_Name, Exac_program, workdir, Nodes=SetCount, \
                            NumCore=NumCores, WallTime=WallTime, partition=Partition)

