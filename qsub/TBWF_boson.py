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

def Copy_Initial_WF(target, destination):
  if os.path.isfile(target):
    cmd = " ".join(["cp", target, destination])
    subprocess.call(cmd, shell=True)
  else:
    print("Initial wavefunction NOT found!", target)
    sys.exit()

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
    qsub_cmd = "qsub -q batch"
    EXEC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization/data"
elif(platform.system() == "Darwin"):
  QSUB = False
  SRC_DIR = "/Volumes/Files/GitRepo/ExactDiagonalization"
  EXEC_DIR = os.path.join(SRC_DIR, "data")

NumThreads = 2

# NOTE: This depends on the algorithm!!
Algorithm = 2
L = 8
N = L - 1
if Algorithm == 2:
  L += 1
OBC = 1# 1:True
# OBC = 0# 0:False
Uin = [0.0, 0.5, 1.0, 3.0, 5.0, 7.0, 9.0]
Vin = 0.0
# NOTE: Dynamics parameters
Tsteps = 2000# Tstep * dt is final time
dt = 0.01
TBloc = L - 1
# TBloc = np.int(L / 2)

#/* 1.0 is HARD_CUT; 0.0 doesn't make sense */
FactorList = [0.01, 0.03, 0.05, 0.1, 0.2, 0.4, 0.6, 1.0]

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "TBWF.b"))
Exac_program = "\n".join(APPs)

if OBC:
  # NOTE: This depends on the algorithm!!
  if Algorithm == 1:
    IWF_DIR = "-".join(["BSSHO", "".join(["L", str(L)]), str(N)])
  elif Algorithm == 2:
    IWF_DIR = "-".join(["BSSHO", "".join(["L", str(L - 1)]), str(N)])
  LN1N2 = "-".join(["TBO", "".join(["L", str(L)]), str(N)])
else:
  # NOTE: This depends on the algorithm!!
  if Algoritm == 1:
    IWF_DIR = "-".join(["BSSHP", "".join(["L", str(L)]), str(N)])
  elif Algoritm == 2:
    IWF_DIR = "-".join(["BSSHP", "".join(["L", str(L - 1)]), str(N)])
  LN1N2 = "-".join(["TBP", "".join(["L", str(L)]), str(N)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

for Factor in FactorList:
  for U in Uin:
    Job_Name =  "-".join(["".join(["U", str(U)]), "Jr1.0-Box"])
    IWF_FILE = os.path.join(EXEC_DIR, IWF_DIR, Job_Name, "BSSH.h5")
    Job_Name =  "-".join([ "".join(["U", str(U)]), "".join(["V", str(Vin)]),
                           "".join(["F", str(Factor)]), "".join(["TB", str(TBloc)])] )
    workdir = os.path.join(DATADIR, Job_Name)

    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
    with shp.cd(workdir):
      if os.path.isfile('TBWF.h5'):
        pass
      else:
        Copy_Initial_WF(IWF_FILE, workdir)
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
        dset = para.create_dataset("FACTOR", data=Factor)
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
