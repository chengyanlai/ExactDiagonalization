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

def SetV(L, Val1=0.0, Val2=0.0, vtype="Uniform"):
  if vtype == "Uniform":
    V = np.ones(L, dtype=np.float64) * Val1
  elif vtype == "SinkCenter":
    V = np.ones(L, dtype=np.float64) * Val1
    if L % 2 == 1:
      V[np.int((L - 1) / 2)] = Val2
    else:
      V[np.int(L / 2)] = Val2
  elif vtype == "SinkEdgeLeft":
    V = np.ones(L, dtype=np.float64) * Val1
    V[0] = Val2
  elif vtype == "SinkEdgeRight":
    V = np.ones(L, dtype=np.float64) * Val1
    V[-1] = Val2
  else:
    print("Vtype is not existed!")
    sys.exit()
  return V

NumThreads = 2

L = 13
OBC = 1# 1:True
# OBC = 0# 0:False
N = L
Uin = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
VtypeEqm = "Uniform"
Val1 = 0.0
Vin = SetV(L, Val1=Val1, vtype=VtypeEqm)
# NOTE: Dynamics parameters
Tsteps = 2000# Tstep * dt is final time
dt = 0.01
VtypeDyn = "SinkCenter"
Val2List = [-4.0, ]

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "SourceSinkDyn.b"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["SSO", "".join(["L", str(L)]), str(N)])
else:
  LN1N2 = "-".join(["SSP", "".join(["L", str(L)]), str(N)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

for Val2 in Val2List:
  for U in Uin:
    Job_Name =  "-".join(["".join(["U", str(U)]),
      "".join([VtypeEqm, str(Val1)]), "".join([VtypeDyn, str(Val2)])])
    workdir = os.path.join(DATADIR, Job_Name)

    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
    with shp.cd(workdir):
      if os.path.isfile('SourceSink.h5'):
        pass
      else:
        f = h5py.File('Eqm.h5', 'w')
        para = f.create_group("Parameters")
        dset = para.create_dataset("L", data=L)
        dset = para.create_dataset("OBC", data=OBC)
        dset = para.create_dataset("N", data=N)
        dset = para.create_dataset("U", data=U)
        dset = para.create_dataset("V", data=Vin)
        f.close()
        f = h5py.File('Dyn.h5', 'w')
        para = f.create_group("Parameters")
        dset = para.create_dataset("Tsteps", data=Tsteps)
        dset = para.create_dataset("dt", data=dt)
        dset = para.create_dataset("U", data=U)
        Vt = SetV(L, Val2=Val2, vtype=VtypeDyn)
        dset = para.create_dataset("V", data=Vt)
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
