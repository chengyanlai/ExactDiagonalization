#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

Space = "R"
L = 3
Nh = 15 * L

OLs = [(0.30, 1.0), (0.40, 1.0), (0.50, 1.0),
       (0.40, 3.0), (0.40, 5.0), (0.40,10.0),]

TSteps = 1000
dt = 0.005

# AlphaReal = np.random.randn(5., 0.2, L)
# AlphaPhase = np.random.uniform(-np.pi, np.pi, L)
AlphaReal = np.array([5., 5., 5.])
AlphaPhase = np.array([0., 2./3.,-2./3.]) * np.pi

APPs = []
Prefix1 = "".join([ "L", str(L), "N", str(Nh) ])
APPs.append(os.path.join(SrcDir, "build", "holstein." + Space.lower() + " 0 400"))
# APPs.append(os.path.join(SrcDir, "build", "holstein.k 1 E 178 179"))
DataDir = os.path.join(ExecDir, "ED", "HTP-"+Space, Prefix1)
APPs.append("/bin/touch DONE")

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
    dset = g.create_dataset("EShift", data=0)
    dset = g.create_dataset("Method", data=0)
    dset = g.create_dataset("AlphaReal", data=AlphaReal)
    dset = g.create_dataset("AlphaPhase", data=AlphaPhase)
    if TSteps:
        dset = g.create_dataset("TSteps", data=TSteps)
        dset = g.create_dataset("dt", data=dt)
    f.close()

    Filename = os.path.join(workdir, 'job')
    if Cluster == "Kagome":
        sg.GenerateScript("PBS", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "Merced":
        sg.GenerateScript("TORQUE", Filename, JobName, APPs, workdir, Nodes=1, NumCore=20, WallTime='336:00:00', Partition='', ProjectName='', MPI=0, PPN=1)
    elif Cluster == "LANL":
        sg.GenerateScript("SLURM", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='16:00:00', Partition='standard', ProjectName='s17_cint')
