#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

dt = 0.005
TSteps = 30000
L = 3
Nh = 20 * L
Prefix1 = "".join([ "L", str(L), "N", str(Nh) ])

Omega, Lambda = 0.4, 1.0
Prefix2 = "-".join([ "".join(["W", str(Omega)]), "".join(["G", str(Lambda)]) ])

APPs = []
APPs.append(os.path.join(SrcDir, "build", "holstein.r 1 Z 0"))
APPs.append("/bin/touch DONE")

DataDir = os.path.join(ExecDir, "ED", "HTP-R", Prefix1, Prefix2)

for cnt in range(1, 40, 1):
    Prefix3 = "S"+str(cnt).zfill(4)
    workdir = os.path.join(DataDir, Prefix3)

    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2

    JobName = Prefix3 + Prefix2 + Prefix1
    f = h5py.File(os.path.join(workdir, 'conf.h5'), 'w')
    g = f.create_group("Parameters")
    dset = g.create_dataset("L", data=L)
    dset = g.create_dataset("N", data=Nh)
    dset = g.create_dataset("W", data=Omega)
    dset = g.create_dataset("G", data=Lambda)
    AlphaReal = np.random.normal(2., 0.2, L)
    AlphaPhase = np.random.uniform(-np.pi, np.pi, L)
    dset = g.create_dataset("AlphaReal", data=AlphaReal)
    dset = g.create_dataset("AlphaPhase", data=AlphaPhase)
    dset = g.create_dataset("TSteps", data=TSteps)
    dset = g.create_dataset("dt", data=dt)
    f.close()

    Filename = os.path.join(workdir, 'job')
    if Cluster == "Kagome":
        sg.GenerateScript("PBS", Filename, JobName, APPs, workdir, Nodes=1, NumCore=4, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "Merced":
        sg.GenerateScript("TORQUE", Filename, JobName, APPs, workdir, Nodes=1, NumCore=20, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "LANL":
        sg.GenerateScript("SLURM", Filename, JobName, APPs, workdir, Nodes=1, NumCore=18, WallTime='16:00:00', Partition='standard', ProjectName='w18_xasrixs', MPI=1, PPN=1)
