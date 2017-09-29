#!/usr/bin/env python3
# coding=utf-8

# install_name_tool -change "@rpath/libarpack.2.dylib" "/Users/chengyanlai/.bin/lib/libarpack.2.0.0.dylib" plex

import os
import sys
import subprocess
import numpy as np
import h5py
import Script_Helpers as shp
from Clusters import *

BL = 1
FL = 10
maxLocalB = 1
if BL == 1:
  Jbbs = [0.0,]
else:
  Jbbs = [0.0, 0.01, 0.05, 0.10]
if FL == 1:
  Jffs = [0.0,]
  Uffs = [0.00]
else:
  Jffs = [0.0, 0.01, 0.05, 0.10]
  Uffs = [0.00, 0.10, 0.50,-1.0,-5.0]
Vbbs = [3.20, 3.50, 3.80]
Vffs = [3.20, 3.40, 3.60]
DeltaDCs = [-0.0080, -0.020, -0.050, -0.080]
CouplingForm = "uniform"# uniform, angle1
# CouplingForm = "angle1"
if CouplingForm == "angle1":
  if BL > FL:
    CouplingForm = "P2"
  if BL < FL:
    CouplingForm = "COS2"

if BL + FL > 11:
  NumCores = 6
else:
  NumCores = 1
WallTime = MaxWallTime

def DCcoupling(Delta, BL, FL, fi=0, form="uniform"):
  if form == "uniform":
    return Delta * np.ones(BL, dtype=np.float64)
  elif form == "P2" and BL > FL:# Bosons surround Fermion
    angles = np.linspace(0, 2.*np.pi, BL+1)
    work = -1.0 * Delta * np.polynomial.legendre.legval(np.cos(angles), [0,0,1])
    return work[:-1]
  elif form == "COS2" and BL < FL:# Fermion surround Bosons
    angles = np.linspace(0, 2.*np.pi, FL+1)
    work = -1.0 * Delta * np.cos(angles[fi]) * np.cos(angles[fi]) * np.ones(BL, dtype=np.float64)
    return work
  else:
    print("Not supported yet")
    sys.exit()

DATADIR = os.path.join( EXEC_DIR, "plex", "".join(["B", str(BL), "F", str(FL), "mB", str(maxLocalB)]), CouplingForm)
for Uff in Uffs:
  for Jbb in Jbbs:
    for Jff in Jffs:
      Job_Name = "".join(["Pl", str(BL), "Ex", str(FL), "JB", str(Jbb), "JF", str(Jff), "Uf", str(Uff)])
      Folder = "".join(["PlExJB", str(Jbb), "JF", str(Jff), "Uf", str(Uff)])
      workdir = os.path.join(DATADIR, Folder)
      os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
      with shp.cd(workdir):
        if os.path.isfile('DONE'):
          print(workdir, " is DONE!")
          continue
        SetCount = 0
        f = h5py.File('confs.h5', 'w')
        for Vbb in Vbbs:
          for Vff in Vffs:
            for DeltaDC in DeltaDCs:
              para = f.create_group("Input-" + str(SetCount))
              dset = para.create_dataset("BL", data=BL)
              dset = para.create_dataset("FL", data=FL)
              dset = para.create_dataset("maxLocalB", data=maxLocalB)
              dset = para.create_dataset("Jbb", data=Jbb)
              dset = para.create_dataset("Jff", data=Jff)
              dset = para.create_dataset("Vbb", data=Vbb)
              dset = para.create_dataset("Vff", data=Vff)
              dset = para.create_dataset("Uff", data=Uff)
              for i in range(FL):
                DeltaDCarr = DCcoupling(DeltaDC, BL, FL, fi=i, form=CouplingForm)
                gname = "DeltaDC-" + str(i)
                dset = para.create_dataset(gname, data=DeltaDCarr)
              SetCount += 1
        f.close()
        if Cluster == "LANL":
          APPs = []
          APPs.append("/bin/cp confs.h5 confs.h5.backup")
          APPs.append("mpirun -n " + str(SetCount) + " -ppn " + str(NumCores) + " " + os.path.join(SRC_DIR, "build", "plex.mpi"))
          Exac_program = "\n".join(APPs)
          shp.WriteQsubSBATCH("job", Job_Name, Exac_program, workdir, \
            Nodes=SetCount, NumCore=NumCores, WallTime=WallTime, partition=Partition)
        elif Cluster == "Kagome":
          APPs = []
          APPs.append(os.path.join(SRC_DIR, "build", "plex " + str(SetCount) ))
          Exac_program = "\n".join(APPs)
          shp.WriteQsubPBS("job", Job_Name, Exac_program, workdir, \
            NumCore=NumCores, WallTime=WallTime)
