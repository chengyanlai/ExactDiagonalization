#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

Nh = 20
# L = 2 * Nh# Gives max LFS, for thermodynamic limit.
L = 16# test iQ in PRB 94 014304 (2016)

## Will multiply pi int the c++ code
# Momentum = np.linspace(0., 1., 101)
Momentum = np.array([1.0,])# test iQ in PRB 94 014304 (2016)

# OLs = [(1.0, 1.0).]
OLs = [(2.0, np.sqrt(2.0)),]# test iQ in PRB 94 014304 (2016)

# Dynamics - for iQ
TSteps = 10000
dt = 0.005

APPs = []
Prefix1 = "".join([ "L", str(L), "Nh", str(Nh) ])
if TSteps:
    APPs.append(os.path.join(SrcDir, "build", "holstein.1dInfty 1"))
    DataDir = os.path.join(ExecDir, "ED", "HTP-iQ", Prefix1)
else:
    APPs.append(os.path.join(SrcDir, "build", "holstein.1dInfty 0"))
    DataDir = os.path.join(ExecDir, "ED", "HTP", Prefix1)
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
    dset = g.create_dataset("Momentum", data=Momentum)
    Wls = np.ones(L, dtype=np.float64) * Omega
    dset = g.create_dataset("W", data=Wls)
    Gls = np.ones(L, dtype=np.float64) * Lambda
    dset = g.create_dataset("G", data=Gls)
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
