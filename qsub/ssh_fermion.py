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
    qsub_cmd = "qsub -q LM.q"
    EXEC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization/data"
elif(platform.system() == "Darwin"):
  QSUB = False
  SRC_DIR = "/Volumes/Files/GitRepo/ExactDiagonalization"
  EXEC_DIR = os.path.join(SRC_DIR, "data")

def SetV(L, vtype="Box"):
  if vtype == "Box":
    V = np.zeros(L, dtype=np.float64)
  else:
    print("Vtype is not existed!")
    sys.exit()
  return V

NumThreads = 1
L = 13
# J12ratio = np.linspace(0.0, 1.0, 51)
J12ratio = [1.00, ]
# OBC = 1#1:True
OBC = 1#0:False
if L % 2 == 1:
  # N1 = np.int((L + 1) / 2)
  # N2 = np.int((L - 1) / 2)
  N1 = 11
  N2 = 2
else:
  N1 = np.int(L / 2)
  N2 = np.int(L / 2)
Uin = [10.0, ]#np.linspace(0.0, 10.0, 51)
if OBC:
  Phils = [0, ]
else:
  Phils = np.linspace(0, L, 66)
Vtype = "Box"
Vin = SetV(L, vtype=Vtype)

APPs = []
APPs.append(os.path.join(SRC_DIR, "build", "SSH.f"))
Exac_program = "\n".join(APPs)

if OBC:
  LN1N2 = "-".join(["SSHO", "".join(["L", str(L)]), str(N1), str(N2)])
else:
  LN1N2 = "-".join(["SSHP", "".join(["L", str(L)]), str(N1), str(N2)])
DATADIR = os.path.join(EXEC_DIR, LN1N2)

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
        if os.path.isfile('SSH.h5'):
          print("".join([workdir, " already done!"]))
          pass
        elif os.path.isfile('SSHconf.h5'):
          print("".join([workdir, " is schaduled!?"]))
          pass
        else:
          f = h5py.File('SSHconf.h5', 'w')
          para = f.create_group("Parameters")
          dset = para.create_dataset("L", data=L)
          dset = para.create_dataset("J12", data=J12)
          dset = para.create_dataset("OBC", data=OBC)
          dset = para.create_dataset("N1", data=N1)
          dset = para.create_dataset("N2", data=N2)
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
              subprocess.call(i, shell=True, stdout=f)
            f.close()
          else:
            qsub_script = " ".join([qsub_cmd, "qsub.ucmerced"])
            subprocess.call(qsub_script, shell=True)