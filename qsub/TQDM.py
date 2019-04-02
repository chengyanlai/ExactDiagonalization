#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import subprocess
import platform
import socket
import numpy as np
import h5py
import ScriptGenerator as sg

# Get all platform settings
from Clusters import *

# Remain const.
t23 = 1.0
Vin = 0.0

SearchJ = 0

if SearchJ:
  t13List = np.linspace(0.7, 1.0, 13)# For searching
  cntGStart = 0
  GammaList = [1.0,]# Not used
  UinTj = [(0.1, 0.18), (1.0, 0.16), (3.0, 0.14)]
else:
  # remember to set FIXJ to 0 in main program
  cntGStart = 0
  t13List = np.linspace(0.1, 1.9, 19)
  Uin = np.linspace(1, 6, 11)
  GammaList = np.logspace(0, 1.2, num=20)
  TargetJ = 0.0# Not used


# NOTE: Dynamics parameters
dt = 0.005

# for U, TargetJ in UinTj:
for U in Uin:
  for t13 in t13List:
    if SearchJ:
      DataDir = os.path.join( ExecDir, "ED", "TQDM", "SearchJ" + str(SearchJ), "".join(["U", str(U), "-t", str(t13)]) )
    else:
      DataDir = os.path.join( ExecDir, "ED", "TQDM", "".join(["U", str(U), "-t", str(np.round(t13,2))]) )
    cntG = cntGStart
    for GammaL in GammaList:
      GammaR = GammaL
      if GammaL > 800.: dt = 0.0001
      elif GammaL > 600.: dt = 0.0002
      elif GammaL > 400.: dt = 0.0005
      elif GammaL > 200.: dt = 0.001
      if SearchJ:
        JobName =  "".join(["U", str(U), "-t", str(np.round(t13,2)), "SJ", str(SearchJ)])
        workdir = DataDir
      else:
        JobName =  "".join(["U", str(U), "-t", str(np.round(t13,2)), "-G", str(cntG),])
        workdir = os.path.join(DataDir, "G"+str(cntG))

      os.makedirs(workdir, exist_ok=True)  # Python >= 3.2

      f = h5py.File( os.path.join(workdir, 'conf.h5'), 'w')
      para = f.create_group("Parameters")
      dset = para.create_dataset("t13", data=t13)
      dset = para.create_dataset("t23", data=t23)
      Uls = np.ones(3, dtype=np.float64) * U
      # Disorder U
      # Uls = np.zeros(3, dtype=np.float64)
      # Uls[1] = U
      dset = para.create_dataset("U", data=Uls)
      Vls = np.zeros(3, dtype=np.float64)
      dset = para.create_dataset("V", data=Vls)
      dset = para.create_dataset("GammaL", data=GammaL)
      dset = para.create_dataset("GammaR", data=GammaR)
      dset = para.create_dataset("dt", data=dt)
      dset = para.create_dataset("TargetJ", data=TargetJ)
      f.close()

      if Cluster == "Kagome":
        Filename = os.path.join(workdir, 'job')
        Executable = [" ".join([os.path.join(SrcDir, "build", "tqdm " + str(SearchJ)), ]), "/bin/touch DONE"]
        sg.GenerateScript("PBS", Filename, JobName, Executable, workdir, Nodes=1, NumCore=1, WallTime='4:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
      elif Cluster == "Merced":
        Filename = os.path.join(workdir, 'job')
        Executable = [" ".join([os.path.join(SrcDir, "build", "tqdm " + str(SearchJ)), ]), "/bin/touch DONE"]
        sg.GenerateScript("TORQUE", Filename, JobName, Executable, workdir, Nodes=1, NumCore=1, WallTime='4:00:00', Partition='', ProjectName='', MPI=0, PPN=1)

      cntG += 1
