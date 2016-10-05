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
    qsub_cmd = "qsub -q LM.q"
    EXEC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization/data"
elif(platform.system() == "Darwin"):
  QSUB = False
  SRC_DIR = "/Volumes/Files/GitRepo/ExactDiagonalization"
  EXEC_DIR = os.path.join(SRC_DIR, "test")
  # EXEC_DIR = os.path.join(SRC_DIR, "data")

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
N = L - 2
Uin = [0.0, 0.5, 1.0, 1.5, 2.0, ]
# Uin = np.linspace(0.0, 5.0, 51)
# Uin = [0.0, 0.5, 1.0, 1.5, 2.0,]# Source-Sink
# Uin = [0.5, 1.0, ]
VtypeEqm = "Uniform"
# VtypeEqm = "SinkEdgeLeft"
Val1 = 0.0
# NOTE: Dynamics parameters
Tsteps = 3000# Tstep * dt is final time
dt = 0.01
VtypeDyn = "SinkCenter"
# VtypeDyn = "SinkEdgeRight"
# Val2List = np.linspace(-2.0, -36.0, 35)
# Val2List = [-3.0, -9.0, -12.0, -15.0, -18.0, -27.0]
# Val2List = np.linspace(-0.2, -10.0, 50)
# Val2List = np.linspace(-9.0, -9.8, 5)
Val2List = [-3.0, -5.0]#-0.4, -0.6, -0.8, -10.0]

# T0stepList = [0, ]
T0stepList = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "SSWF.b"))
Exac_program = "\n".join(APPs)

for T0step in T0stepList:
  if OBC:
    LN1N2 = "-".join(["BSSO", "".join(["L", str(L)]), str(N), "".join(["T", str(T0step)])])
  else:
    LN1N2 = "-".join(["BSSP", "".join(["L", str(L)]), str(N), "".join(["T", str(T0step)])])
  DATADIR = os.path.join(EXEC_DIR, LN1N2)

  job_id = 0
  RUN_NOW = True

  for Val2 in Val2List:
    for U in Uin:
      Job_Name =  "-".join(["".join(["U", str(U)]),
        "".join([VtypeEqm, str(Val1)]), "".join([VtypeDyn, str(Val2)])])
      workdir = os.path.join(DATADIR, Job_Name)

      os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
      with shp.cd(workdir):
        if os.path.isfile('SourceSink.h5'):
          pass
        elif os.path.isfile('Eqm.h5') and os.path.isfile('Dyn.h5'):
          print("".join([workdir, " is schaduled!?"]))
          pass
        else:
          f = h5py.File('Eqm.h5', 'w')
          para = f.create_group("Parameters")
          dset = para.create_dataset("L", data=L)
          dset = para.create_dataset("OBC", data=OBC)
          dset = para.create_dataset("N", data=N)
          dset = para.create_dataset("U", data=U)
          Vin = SetV(L, Val1=Val1, Val2=Val2, vtype=VtypeEqm)
          dset = para.create_dataset("V", data=Vin)
          f.close()
          f = h5py.File('Dyn.h5', 'w')
          para = f.create_group("Parameters")
          dset = para.create_dataset("Tsteps", data=Tsteps)
          dset = para.create_dataset("dt", data=dt)
          dset = para.create_dataset("U", data=U)
          Vt = SetV(L, Val1=Val1, Val2=Val2, vtype=VtypeDyn)
          dset = para.create_dataset("V", data=Vt)
          dset = para.create_dataset("T0step", data=T0step)
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
            if RUN_NOW:
              qsub_script = " ".join([qsub_cmd, "qsub.ucmerced"])
              if QSUB: subprocess.call(qsub_script, shell=True)
            elif job_id:
              after_id = "".join(["depend=afterany:", str(job_id), ".", socket.gethostname()])
              qsub_script = " ".join([qsub_cmd, "-W", after_id, "qsub.ucmerced"])
              print("Run Command - ", qsub_script)
              subprocess.call(qsub_script, shell=True)
              job_id += 1
              RUN_NOW = False
