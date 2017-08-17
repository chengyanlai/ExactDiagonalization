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

# t13List = np.linspace(0.9, 1.1, 5)
t13List = np.linspace(0.1, 1.9, 19)
t23 = 1.0
Vin = 0.0
Uin = [0.0, 1.0, 0.1, 5.0, 10.0]
# GammaList = np.logspace(-0.3, 1.8, num=50)
# GammaList = np.logspace(0., 1.2, num=20)
GammaList = [0.1,]# Small gamma search for j13 = 0.1

# NOTE: Dynamics parameters
Tsteps = 8000# Tstep * dt is final time
dt = 0.005

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "loop.f"))
Exac_program = "\n".join(APPs)

QSUB = True
NumThreads = 1
WallTime = MaxWallTime

for U in Uin:
  for t13 in t13List:
    DATADIR = os.path.join( EXEC_DIR, "TriPX", "".join(["TriU", str(U), "-t", str(t13)]) )
    cntG = 0
    for GammaL in GammaList:
      GammaR = GammaL
      # Job_Name =  "".join(["G", str(cntG),])
      Job_Name =  "".join(["Search", str(cntG),])
      workdir = os.path.join(DATADIR, Job_Name)

      os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
      with shp.cd(workdir):
        if os.path.isfile('TriPE.h5'):
          pass
        elif os.path.isfile('conf.h5'):
          print("".join([workdir, " is schaduled!?"]))
          pass
        else:
          f = h5py.File('conf.h5', 'w')
          para = f.create_group("Parameters")
          dset = para.create_dataset("t13", data=t13)
          dset = para.create_dataset("t23", data=t23)
          Uls = np.ones(3, dtype=np.float64) * U
          dset = para.create_dataset("U", data=Uls)
          Vls = np.zeros(3, dtype=np.float64)
          dset = para.create_dataset("V", data=Vls)
          dset = para.create_dataset("gammaL", data=GammaL)
          dset = para.create_dataset("gammaR", data=GammaR)
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

          if platform.system() == "Darwin":
            f = open(Job_Name, "w")
            for i in APPs:
              print(i)
              # subprocess.call(i, shell=True, stdout=f)
            f.close()
          elif QSUB:
            qsub_script = " ".join([qsub_cmd, "job"])
            if QSUB: subprocess.call(qsub_script, shell=True)
      cntG += 1
