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
    # print("Run Command - ", cmd)
    subprocess.call(cmd, shell=True)
  else:
    print("Initial wavefunction NOT found!", target)
    sys.exit()

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
    qsub_cmd = "qsub -q LM.q"
    EXEC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization/data"
elif platform.system() == "Darwin":
  SRC_DIR = "/Volumes/Files/GitRepo/ExactDiagonalization"
  EXEC_DIR = os.path.join(SRC_DIR, "data")

NumThreads = 2

L = 8
TBloc = L - 1
# TBloc = np.int(L / 2)
# Uin = [0.0, 0.5, 1.0, 3.0, 5.0, 7.0, 9.0]
# GammaList = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 10.0, 15.0, 20.0, 25.0]
# L = 9
Uin = [0.0, ]#
# TBloc = L - 1
GammaList = [20.0, ]#

OBC = 1# 1:True
# OBC = 0# 0:False
N = L - 1
Vin = 0.0
# NOTE: Dynamics parameters
Tsteps = 2000# Tstep * dt is final time
dt = 0.01

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "TBLB.b"))
Exac_program = "\n".join(APPs)

if OBC:
  IWF_DIR = "-".join(["BSSHO", "".join(["L", str(L)]), str(N)])
  LN1N2 = "-".join(["TBO", "".join(["L", str(L)]), str(N)])
  LN1N2 = "-".join(["TBO", "".join(["L", str(L)]), str(N),"test"])
else:
  IWF_DIR = "-".join(["BSSHP", "".join(["L", str(L)]), str(N)])
  LN1N2 = "-".join(["TBP", "".join(["L", str(L)]), str(N)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

job_id = 0
RUN_NOW = True
QSUB = True
for Gamma in GammaList:
  for U in Uin:
    Job_Name =  "-".join(["".join(["U", str(U)]), "Jr1.0-Box"])
    IWF_FILE = os.path.join(EXEC_DIR, IWF_DIR, Job_Name, "BSSH.h5")
    Job_Name =  "-".join( ["".join(["U", str(U)]), "".join(["V", str(Vin)]),
                           "".join(["G", str(Gamma)]), "".join(["TB", str(TBloc)])] )
    workdir = os.path.join(DATADIR, Job_Name)

    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
    with shp.cd(workdir):
      if os.path.isfile('TBLB.h5'):
        pass
      elif os.path.isfile('conf.h5'):
        print("".join([workdir, " is schaduled!?"]))
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
        dset = para.create_dataset("Gamma", data=Gamma)
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
