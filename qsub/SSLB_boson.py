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

if platform.system() == "Linux":
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
elif platform.system() == "Darwin":
  SRC_DIR = "/Volumes/Files/GitRepo/ExactDiagonalization"
  EXEC_DIR = os.path.join(SRC_DIR, "data")

def SetV(L, Val1=0.0, Val2=0.0, vtype="Uniform"):
  if vtype == "Uniform":
    V = np.ones(L, dtype=np.float64) * Val1
    S = []
  elif vtype == "SinkCenter":
    V = np.ones(L, dtype=np.float64) * Val1
    if L % 2 == 1:
      V[np.int((L - 1) / 2)] = Val2
      S = [np.int((L - 1) / 2), ]
    else:
      V[np.int(L / 2)] = Val2
      S = [np.int(L / 2), ]
  elif vtype == "SinkEdgeLeft":
    V = np.ones(L, dtype=np.float64) * Val1
    V[0] = Val2
    S = [0, ]
  elif vtype == "SinkEdgeRight":
    V = np.ones(L, dtype=np.float64) * Val1
    V[-1] = Val2
    S = [L - 1, ]
  else:
    print("Vtype is not existed!")
    sys.exit()
  return V, S

NumThreads = 2

L = 8
# Uin = [0.0, 0.5, 1.0, 3.0, 5.0, 7.0, 9.0]
# GammaList = [0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 7.0, 9.0, 10.0, 15.0, 20.0, 25.0]
# L = 9
Uin = [0.0, ]#
GammaList = [3.50, ]#

OBC = 1# 1:True
# OBC = 0# 0:False

# Particle number
N = L - 1

# Potential
VtypeEqm = "Uniform"
Val1 = 0.0
Vin, Sin = SetV(L, Val1=Val1, vtype=VtypeEqm)
VtypeDyn = "SinkCenter"
# Val2List = [-3.0, -6.0, -9.0, -12.0, -18.0, -27.0]
Val2List = [-27.0, ]

# Dissapation sites
SitesType = "All"
# SitesType = "Sink"

# NOTE: Dynamics parameters
Tsteps = 2000# Tstep * dt is final time
dt = 0.01

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "SSLB.b"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["SSLBO", "".join(["L", str(L)]), str(N)])
else:
  LN1N2 = "-".join(["SSLBP", "".join(["L", str(L)]), str(N)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

job_id = 0
RUN_NOW = True
QSUB = True
for Gamma in GammaList:
  for Val2 in Val2List:
    for U in Uin:
      Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["G", str(Gamma)]),
        "".join([VtypeEqm, str(Val1)]), "".join([VtypeDyn, str(Val2)]), SitesType])
      workdir = os.path.join(DATADIR, Job_Name)

      os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
      with shp.cd(workdir):
        if os.path.isfile('SSLB.h5'):
          pass
        elif os.path.isfile('conf.h5'):
          print("".join([workdir, " is schaduled!?"]))
          pass
        else:
          f = h5py.File('conf.h5', 'w')
          para = f.create_group("Parameters")
          dset = para.create_dataset("L", data=L)
          dset = para.create_dataset("OBC", data=OBC)
          dset = para.create_dataset("N", data=N)
          dset = para.create_dataset("U", data=U)
          dset = para.create_dataset("Tsteps", data=Tsteps)
          dset = para.create_dataset("dt", data=dt)
          dset = para.create_dataset("Gamma", data=Gamma)
          dset = para.create_dataset("Veqm", data=Vin)
          Vt, St = SetV(L, Val2=Val2, vtype=VtypeDyn)
          dset = para.create_dataset("Vdyn", data=Vt)
          if SitesType == "Sink":
            dset = para.create_dataset("Sites", data=St)
          elif SitesType == "All":
            dset = para.create_dataset("Sites", data=np.arange(L))
          f.close()

          if socket.gethostname() == 'kagome.rcc.ucmerced.edu' or \
             socket.gethostname() == 'edgestate.rcc.ucmerced.edu' or \
             socket.gethostname() == 'atomtronics.ucmerced.edu':
            shp.WriteQsubPBS("qsub.ucmerced", Job_Name, Exac_program, workdir,
                             NodeName=NodeName, NumCore=NumThreads)
          if platform.system() == "Darwin":
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
