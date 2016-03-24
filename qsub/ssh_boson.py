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
# J12ratio = [0.10, ]#0.20, 0.40, 0.60, 0.80, 1.00]
# NOTE: Prepare for Terminator Beam
# J12ratio = [1.00, ]
# OBC = 1# 1:True, 0:False
# N = L - 1
# Uin = [0.0, 0.5, 1.0, 3.0, 5.0, 7.0, 9.0]
# Vtype = "Uniform"
# Vin = SetV(L, vtype=Vtype)
# NOTE: Prepare for Terminator Beam
# NOTE: Eqm. Sink
J12ratio = [1.00, ]
OBC = 1# 1:True, 0:False
N = L - 2
Uin = [1.0, ]
# Uin = np.linspace(0.1, 5.0, 50)
# Uin = np.linspace(0.0, 5.0, 51)
# Vtype = "SinkCenter"
Vtype = "SinkEdgeLeft"
V1 = 0.0
V2 = -5.0
Vin = SetV(L, Val1=V1, Val2=V2, vtype=Vtype)
# NOTE: Eqm. Sink
if OBC:
  Phils = [0, ]
else:
  Phils = np.linspace(0, L, 66)

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "SSH.b"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["BSSHO", "".join(["L", str(L)]), str(N)])
else:
  LN1N2 = "-".join(["BSSHP", "".join(["L", str(L)]), str(N)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

job_id = 2597
RUN_NOW = False
QSUB = True
for nphi in Phils:
  phi = nphi * 2.0 * np.pi / np.float64(L)
  for J12 in J12ratio:
    for U in Uin:
      if OBC:
        Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["Jr", str(J12)]),
          Vtype])
      else:
        Job_Name =  "-".join(["".join(["U", str(U)]), "".join(["Jr", str(J12)]),
          "".join(["n", str(nphi)]), Vtype])

      workdir = os.path.join(DATADIR, Job_Name)

      os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
      with shp.cd(workdir):
        if os.path.isfile('BSSH.h5'):
          pass
        else:
          f = h5py.File('conf.h5', 'w')
          para = f.create_group("Parameters")
          dset = para.create_dataset("L", data=L)
          dset = para.create_dataset("J12", data=J12)
          dset = para.create_dataset("OBC", data=OBC)
          dset = para.create_dataset("N", data=N)
          dset = para.create_dataset("U", data=U)
          dset = para.create_dataset("phi", data=phi)
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
