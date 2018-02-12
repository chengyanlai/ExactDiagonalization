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
  t13List = np.linspace(0.9, 1.2, 31)# For searching
  cntGStart = 0
  GammaList = [1.0,]# Small gamma search
  # cntGStart = 1
  # GammaList = [9.0,]# Large gamma search
  Uin = [0.0,]
  TargetJ = 0.2# for searching U = 0.0
  # Uin = [0.1,]
  # TargetJ13 = 0.19# for searching U = 0.1
  # Uin = [1.0,]
  # TargetJ13 = 0.19# for searching U = 1.0
  # Uin = [5.0,]
  # TargetJ13 = 0.13# for searching U = 5.0
else:
  # remember to set FIXJ to 0 in main program
  cntGStart = 0
  t13List = np.array([0.6, 0.8, 1.0, 1.2, 1.4])
  Uin = [0.0, 0.1, 1.0, 5.0]
  GammaList = np.logspace(-2.0, 2.5, num=251)


# NOTE: Dynamics parameters
dt = 0.005

for U in Uin:
  for t13 in t13List:
    if SearchJ:
      DataDir = os.path.join( ExecDir, "ED", "TQDM", "SearchJ", "".join(["U", str(U), "-t", str(t13)]) )
    else:
      DataDir = os.path.join( ExecDir, "ED", "TQDM", "".join(["U", str(U), "-t", str(t13)]) )
    cntG = cntGStart
    for GammaL in GammaList:
      GammaR = GammaL
      if GammaL > 800.: dt = 0.0001
      elif GammaL > 600.: dt = 0.0002
      elif GammaL > 400.: dt = 0.0005
      elif GammaL > 200.: dt = 0.001
      JobName =  "".join(["G", str(cntG),])
      workdir = os.path.join(DataDir, JobName)

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
      if SearchJ: dset = para.create_dataset("TargetJ", data=TargetJ)
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
