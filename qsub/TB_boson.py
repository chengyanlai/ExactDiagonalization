#!/usr/bin/env python
# coding=utf-8

import os
import sys
import subprocess
import platform
import socket
import numpy as np
import h5py
import Script_Helpers as shp

if(platform.system() == "Linux"):
  QSUB = True
  if socket.gethostname() == 'stargate.phys.nthu.edu.tw':
    qsub_cmd = "qsub -q short.q"
    SRC_DIR = "/home/chenyen/GitRepo/ExactDiagonalization"
    EXEC_DIR = os.path.join(SRC_DIR, "data")
  elif socket.gethostname() == 'kagome.rcc.ucmerced.edu':
    SRC_DIR = "/condensate1/GitRepo/ExactDiagonalization"
    if sys.argv[1] == "k":
      NodeName = "kagome.rcc.ucmerced.edu"
      qsub_cmd = "qsub -q short.q"
      EXEC_DIR = "/home/chengyanlai/data/SSH-ED"
    elif sys.argv[1] == "c":
      NodeName = "condensate.rcc.ucmerced.edu"
      qsub_cmd = "qsub -q batch.q"
      EXEC_DIR = os.path.join(SRC_DIR, "data")
  elif socket.gethostname() == 'edgestate.rcc.ucmerced.edu':
    SRC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization"
    NodeName = "edgestate.rcc.ucmerced.edu"
    qsub_cmd = "qsub -q batch.q"
    EXEC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization/data"
  elif socket.gethostname() == 'atomtronics.ucmerced.edu':
    SRC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization"
    NodeName = "atomtronics.ucmerced.edu"
    qsub_cmd = "qsub -q batch.q"
    EXEC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization/data"
elif(platform.system() == "Darwin"):
  QSUB = False
  SRC_DIR = "/Volumes/Files/GitRepo/ExactDiagonalization"
  EXEC_DIR = os.path.join(SRC_DIR, "data")

NumThreads = 2

L = 4
OBC = 1# 1:True
# OBC = 0# 0:False
N = L - 1
Uin = [0.0, ]#1.0, 3.0, 5.0, 7.0, 9.0]
Vin = 0.0
# NOTE: Dynamics parameters
Tsteps = 2000# Tstep * dt is final time
dt = 0.01
TBloc = L - 1

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "TB.b"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["TBO", "".join(["L", str(L)]), str(N)])
else:
  LN1N2 = "-".join(["TBP", "".join(["L", str(L)]), str(N)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

for U in Uin:
  Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["V", str(Vin)])])
  workdir = os.path.join(DATADIR, Job_Name)

  os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
  with shp.cd(workdir):
    if os.path.isfile('TB.h5'):
      pass
    else:
      f = h5py.File('conf.h5', 'w')
      para = f.create_group("Parameters")
      dset = para.create_dataset("L", data=L)
      dset = para.create_dataset("OBC", data=OBC)
      dset = para.create_dataset("N", data=N)
      dset = para.create_dataset("U", data=U)
      dset = para.create_dataset("V", data=Vin)
      dset = para.create_dataset("Tsteps", data=Tsteps)
      dset = para.create_dataset("dt", data=dt)
      dset = para.create_dataset("TBloc", data=TBloc)
      f.close()

      if socket.gethostname() == 'kagome.rcc.ucmerced.edu' or \
         socket.gethostname() == 'edgestate.rcc.ucmerced.edu' or \
         socket.gethostname() == 'atomtronics.ucmerced.edu':
        shp.WriteQsubPBS("qsub.ucmerced", Job_Name, Exac_program, workdir,
                         NodeName=NodeName, NumCore=NumThreads)
      if not QSUB:
        f = open(Job_Name, "w")
        for i in APPs:
          print(i)
          subprocess.call(i, shell=True, stdout=f)
        f.close()
      else:
        qsub_script = " ".join([qsub_cmd, "qsub.ucmerced"])
        subprocess.call(qsub_script, shell=True)
