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

BL = 4
FL = 1
maxLocalB = 1
if BL == 1:
  Jbbs = [0.0,]
else:
  Jbbs = [0.01, 0.05, 0.10]
if FL == 1:
  Jffs = [0.0,]
  Uffs = [0.00]
else:
  Jffs = [0.01, 0.05, 0.10]
  Uffs = [0.00, 0.10, 0.50,-0.10,-0.50]
Vbbs = [3.20, 3.50, 3.80]
Vffs = [3.20, 3.40, 3.60]
DeltaDCs = [-0.0080, -0.020, -0.050, -0.080]
CouplingForm = "uniform"#

NumCores = 10
WallTime = MaxWallTime

def DCcoupling(Delta, FL, BL, form="uniform"):
  if form == "uniform":
    return np.array([Delta,] * BL)
  elif form == "dipole" and BL > FL:
    angles = np.linspace(0, 2.*np.pi, BL)
    work = np.zeros(BL, dtype=np.float64)
    return work
  elif form == "dipole" and BL < FL:
    angles = np.linspace(0, 2.*np.pi, BL)
    work = np.zeros(BL, dtype=np.float64)
    return work

def multiNodes():
  DATADIR = os.path.join( EXEC_DIR, "plex", "".join(["B", str(BL), "F", str(FL), "mB", str(maxLocalB)]), CouplingForm)
  for Uff in Uffs:
    for Jbb in Jbbs:
      for Jff in Jffs:
        Job_Name = "".join(["PlExJB", str(Jbb), "JF", str(Jff), "Uf", str(Uff)])
        workdir = os.path.join(DATADIR, Job_Name)
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
                DeltaDCarr = DCcoupling(DeltaDC, BL, FL, form=CouplingForm)
                for i in range(FL):
                  dset = para.create_dataset("DeltaDC-"+str(i), data=DeltaDCarr)
                SetCount += 1
          f.close()
          APPs = []
          APPs.append("/bin/cp confs.h5 confs.h5.backup")
          APPs.append("mpirun -n " + str(SetCount) + " -ppn " + str(NumCores) + " " + os.path.join(SRC_DIR, "build", "plex.mpi"))
          Exac_program = "\n".join(APPs)
          shp.WriteQsubSBATCH("job", Job_Name, Exac_program, workdir, \
            Nodes=SetCount, NumCore=NumCores, WallTime=WallTime, partition=Partition)

def singleNode():
  APPs = []
  APPs.append(os.path.join(SRC_DIR, "build", "plex"))
  Exac_program = "\n".join(APPs)
  DATADIR = os.path.join( EXEC_DIR, "Plex", "".join(["B", str(BL), "F", str(FL), "mB", str(maxLocalB), CouplingForm]) )
  for Uff in Uffs:
    for Jbb in Jbbs:
      for Jff in Jffs:
        for Vbb in Vbbs:
          for Vff in Vffs:
            for DeltaDC in DeltaDCs:
              Job_Name =  "".join(["D", str(DeltaDC), "BJ", str(Jbb), "V", str(Vbb), "FJ", str(Jff), "V", str(Vff), "U", str(Uff)])
              workdir = os.path.join(DATADIR, Job_Name)

              os.makedirs(workdir, exist_ok=True)  # Python >= 3.2
              with shp.cd(workdir):
                if os.path.isfile('DONE'):
                  print(workdir, " is DONE!")
                  continue
                else:
                  f = h5py.File('confs.h5', 'w')
                  para = f.create_group("Input-0")
                  dset = para.create_dataset("BL", data=BL)
                  dset = para.create_dataset("FL", data=FL)
                  dset = para.create_dataset("maxLocalB", data=maxLocalB)
                  dset = para.create_dataset("Jbb", data=Jbb)
                  dset = para.create_dataset("Jff", data=Jff)
                  dset = para.create_dataset("Vbb", data=Vbb)
                  dset = para.create_dataset("Vff", data=Vff)
                  dset = para.create_dataset("Uff", data=Uff)
                  DeltaDCarr = DCcoupling(DeltaDC, BL, FL, form=CouplingForm)
                  for i in range(FL):
                    dset = para.create_dataset("DeltaDC-"+str(i), data=DeltaDCarr)
                  f.close()

                  if Cluster == "Comet" or Cluster == "Stampede":
                    shp.WriteQsubSBATCH("job", Job_Name, Exac_program, workdir, \
                      NumCore=NumCores, WallTime=WallTime, partition=Partition)
                  elif Cluster == "Merced":
                    shp.WriteQsubSGE("job", Job_Name, Exac_program, workdir, \
                      NumCore=NumCores, WallTime=WallTime)
                  else:
                    shp.WriteQsubPBS("job", Job_Name, Exac_program, workdir, \
                      NumCore=NumCores, WallTime=WallTime)
                    f = open(Job_Name, "w")
                    for i in APPs:
                      subprocess.call(i, shell=True, stdout=f)
                    f.close()

if __name__ == '__main__':
  if sys.argv[1] == "mpi":
    multiNodes()
  else:
    singleNode()

