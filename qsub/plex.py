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
Jbb = 0.10
Jff = 0.10
Vbbs = [3.00, 3.20, 3.40]
Vffs = [3.00, 3.50, 4.00]
Uffs = [0.10, 1.00, 10.0]
DeltaDC = -0.080

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "plex"))
Exac_program = "\n".join(APPs)

QSUB = True
NumThreads = 1
WallTime = MaxWallTime

for Uff in Uffs:
  for Vbb in Vbbs:
    for Vff in Vffs:
        DATADIR = os.path.join( EXEC_DIR, "Plex", "".join(["B", str(BL), "F", str(FL), "mB", str(maxLocalB)]) )
        Job_Name =  "".join(["D", str(DeltaDC), "BJ", str(Jbb), "V", str(Vbb), "FJ", str(Jff), "V", str(Vff), "U", str(Uff)])
        workdir = os.path.join(DATADIR, Job_Name)

        os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
        with shp.cd(workdir):
          if os.path.isfile('plex.h5'):
            pass
          elif os.path.isfile('confs.h5'):
            print("".join([workdir, " is schaduled!?"]))
            pass
          else:
            f = h5py.File('confs.h5', 'w')
            para = f.create_group("Parameters")
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

            if Cluster == "Comet" or Cluster == "Stampede":
              shp.WriteQsubSBATCH("job", Job_Name, Exac_program, workdir, \
                NumCore=NumThreads, WallTime=WallTime, partition=Partition)
            elif Cluster == "Merced":
              shp.WriteQsubSGE("job", Job_Name, Exac_program, workdir, \
                NumCore=NumThreads, WallTime=WallTime)
            else:
              shp.WriteQsubPBS("job", Job_Name, Exac_program, workdir, \
                NumCore=NumThreads, WallTime=WallTime)

            if platform.system() == "Darwin":
              f = open(Job_Name, "w")
              for i in APPs:
                # print(i)
                subprocess.call(i, shell=True, stdout=f)
              f.close()
            elif QSUB:
              qsub_script = " ".join([qsub_cmd, "job;sleep 1"])
              if QSUB: subprocess.call(qsub_script, shell=True)

