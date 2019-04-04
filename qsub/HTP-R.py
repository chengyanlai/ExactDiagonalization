#!/usr/bin/env python
# coding=utf-8
import subprocess
import os
import numpy as np
import h5py
import ScriptGenerator as sg
from Clusters import *

dt = 0.005
#L, Nh = 3, 90
L, Nh = 4, 100
Prefix1 = "".join([ "L", str(L), "N", str(Nh) ])

#Omega, Lambda, Alpha = 0.1, 0.3, 5.0
#Omega, Lambda, Alpha = 0.3, 0.3, 5.0
#Omega, Lambda, Alpha = 0.5, 0.3, 5.0
#Omega, Lambda, Alpha = 0.7, 0.3, 5.0
#Omega, Lambda, Alpha = 0.9, 0.3, 5.0
Omega, Lambda, Alpha = 0.3, 0.4, 3.750
#Omega, Lambda, Alpha = 0.3, 0.5, 3.0
#Omega, Lambda, Alpha = 0.3, 0.6, 2.50
#Omega, Lambda, Alpha = 0.3, 0.8, 1.8750
#Omega, Lambda, Alpha = 0.3, 0.5, 1.0
#Omega, Lambda, Alpha = 0.3, 0.5, 4.0
#Omega, Lambda, Alpha = 0.3, 0.5, 5.50
#Omega, Lambda, Alpha = 0.3, 0.4, 0
TSteps = np.int(200 * 6.5 * 10 / Omega)# Roughly 10 Phonon Period as 6.5 > 2 * pi
if Alpha > 0.0:
    Prefix2 = "-".join([ "".join(["W", str(Omega)]), "".join(["G", str(Lambda)]), str(Alpha) ])
else:
    Prefix2 = "-".join([ "".join(["W", str(Omega)]), "".join(["G", str(Lambda)]) ])

APPs = []
APPs.append(os.path.join(SrcDir, "build", "holstein.r 1 Z 0"))
APPs.append("/bin/touch DONE")

DataDir = os.path.join(ExecDir, "ED", "HTP-R", Prefix1, Prefix2)
if Alpha > 0.0:
    From = 1001
    To = 1002
else:
    From = 0
    To = 1000
for cnt in range(From, To, 1):
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
    if Alpha > 0.0:
        AlphaReal = np.ones(L, dtype=np.float64) * Alpha
        AlphaPhase = np.linspace(0., 2.*np.pi, L, endpoint=0)
    else:
        AlphaReal = np.random.normal(3.5, 0.5, L)
        AlphaPhase = np.random.uniform(0, 2*np.pi, L)
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
        sg.GenerateScript("SLURM", Filename, JobName, APPs, workdir, Nodes=1, NumCore=16, WallTime='16:00:00', Partition='standard', ProjectName='s17_cint', MPI=0, PPN=1)
