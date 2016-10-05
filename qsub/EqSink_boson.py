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
      EXEC_DIR = "/home/chengyanlai/data/ED"
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

L = 15
N = L - 2
OBC = 1# 1:True, 0:False
Uin = [1.0, ]
# Uin = [0.0, ]
# Uin = np.linspace(0.1, 5.0, 50)
# Uin = np.linspace(0.0, 5.0, 51)
Vtype = "SinkEdgeLeft"
V1 = 0.0
# V2List = [0.0, ]
V2List = np.arange(0.1, , 0.1)
# V2List = [-10.0, ]#-20.0, -30.0, -40.0]

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "SSH.b"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["EqmSink", "".join(["L", str(L)]), str(N)])
else:
  LN1N2 = "-".join(["EqmSinkP", "".join(["L", str(L)]), str(N)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

job_id = 0
RUN_NOW = True
QSUB = True
for V2 in V2List:
  Vin = SetV(L, Val1=V1, Val2=V2, vtype=Vtype)
  for U in Uin:
    Job_Name =  "-".join(["".join(["U", str(U)]), "".join([Vtype, str(V2)]) ])

    workdir = os.path.join(DATADIR, Job_Name)

    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
    with shp.cd(workdir):
      if os.path.isfile('BSSH.h5'):
        pass
      else:
        f = h5py.File('conf.h5', 'w')
        para = f.create_group("Parameters")
        dset = para.create_dataset("L", data=L)
        dset = para.create_dataset("J12", data=1.0)
        dset = para.create_dataset("OBC", data=OBC)
        dset = para.create_dataset("N", data=N)
        dset = para.create_dataset("U", data=U)
        dset = para.create_dataset("phi", data=0)
        dset = para.create_dataset("V", data=Vin)
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
            # subprocess.call(i, shell=True, stdout=f)
          f.close()
        else:
          if RUN_NOW:
            qsub_script = " ".join([qsub_cmd, "qsub.ucmerced"])
            if QSUB: subprocess.call(qsub_script, shell=True)
          elif job_id:
            after_id = "".join(["depend=afterany:", str(job_id), ".", socket.gethostname()])
            qsub_script = " ".join([qsub_cmd, "-W", after_id, "qsub.ucmerced"])
            print("Run Command - ", qsub_script)
            if QSUB: subprocess.call(qsub_script, shell=True)
            job_id += 1
            RUN_NOW = False
          else:
            print("Nothing Done.")
