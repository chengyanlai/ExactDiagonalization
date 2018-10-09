#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

UseProject = "s17_cint"
NumCores = MaxNumThreads

# Space = "R"
Space = "K"

OLs = [(0.30, 1.00), (0.40, 1.00), (0.50, 1.00),
       (0.60, 0.90), (0.60, 1.00), (0.60, 1.10),]

dt = 0.005

APPs = []

if Space == "R":
    TSteps = 30000
    Arpack = 1000#* No use
    L = 3
    Nh = 20 * L
    Prefix1 = "".join([ "L", str(L), "N", str(Nh) ])
    AlphaReal = np.random.normal(2., 0.2, L)
    AlphaPhase = np.random.uniform(-np.pi, np.pi, L)
    # Prefix1 = "".join([ "L", str(L), "N", str(Nh), "D" ])
    # AlphaReal = np.array([3.0, 2.8, 3.2])
    # AlphaReal = np.ones(L) * 3.0
    # AlphaPhase = np.linspace(0., 2., L, endpoint=False) * np.pi
    APPs.append(os.path.join(SrcDir, "build", "holstein." + Space.lower() + " 1 Z 3"))
elif Space == "K":
    TSteps = 0
    Arpack = 0
    L = 4
    Nh = 16 * L
    Prefix1 = "".join([ "L", str(L), "N", str(Nh) ])
    APPs.append(os.path.join(SrcDir, "build", "holstein." + Space.lower() + " 0 0"))

APPs.append("/bin/touch DONE")
DataDir = os.path.join(ExecDir, "ED", "HTP-"+Space, Prefix1)

for (Omega, Lambda) in OLs:
    Prefix2 = "-".join([ "".join(["W", str(Omega)]), "".join(["G", str(Lambda)]) ])

    workdir = os.path.join(DataDir, Prefix2)

    os.makedirs(workdir, exist_ok=True)  # Python >= 3.2

    JobName = Prefix2 + Prefix1
    f = h5py.File(os.path.join(workdir, 'conf.h5'), 'w')
    g = f.create_group("Parameters")
    dset = g.create_dataset("L", data=L)
    dset = g.create_dataset("N", data=Nh)
    dset = g.create_dataset("W", data=Omega)
    dset = g.create_dataset("G", data=Lambda)
    dset = g.create_dataset("ARPACK", data=Arpack)
    if Space == "R":
        dset = g.create_dataset("AlphaReal", data=AlphaReal)
        dset = g.create_dataset("AlphaPhase", data=AlphaPhase)
    if TSteps:
        dset = g.create_dataset("TSteps", data=TSteps)
        dset = g.create_dataset("dt", data=dt)
    f.close()

    Filename = os.path.join(workdir, 'job')
    if Cluster == "Kagome":
        sg.GenerateScript("PBS", Filename, JobName, APPs, workdir, Nodes=1, NumCore=4, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "Merced":
        sg.GenerateScript("TORQUE", Filename, JobName, APPs, workdir, Nodes=1, NumCore=20, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "LANL":
        sg.GenerateScript("SLURM", Filename, JobName, APPs, workdir, Nodes=1, NumCore=NumCores, WallTime='16:00:00', Partition='standard', ProjectName=UseProject)
