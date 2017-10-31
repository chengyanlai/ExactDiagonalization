#!/usr/bin/env python3
# coding=utf-8
import platform
import socket
import os

if platform.system() == "Linux":
    if socket.gethostname() == 'kagome.ucmerced.edu':
        Cluster = "Kagome"
        Partition = ""
        Project = ""
        MaxNumThreads = 16
        MaxWallTime = '720:0:0'
        SrcDir = "/home/chengyanlai/GitRepo/ExactDiagonalization"
        ExecDir = "/home/chengyanlai/Data"
    elif socket.gethostname() == 'merced.cluster':
        Cluster = "Merced"
        Partition = ""
        Project = ""
        MaxNumThreads = 20
        MaxWallTime = '24:0:0'
        SrcDir = "/home/clai24/GitRepo/ExactDiagonalization"
        ExecDir = "/data/clai24"
    elif socket.gethostname()[2:5] == "-fe":
        Cluster = "LANL"
        Partition = "standard"
        Project = "s17_cint"
        MaxNumThreads = 16
        if socket.gethostname()[:2] == "gr":
            MaxNumThreads = 20
        MaxWallTime = '16:0:0'
        SrcDir = "/usr/projects/cint/cint_sces/ExactDiagonalization"
        ExecDir = "/net/scratch3/chengyanlai"
if platform.system() == "Darwin":
    if socket.gethostname() == 'pn1716764.lanl.gov':
        Cluster = "LANLMacPro"
        Partition = ""
        Project = ""
        MaxNumThreads = 12
        MaxWallTime = '48:0:0'
        SrcDir = "/Users/chengyanlai/GitRepo/ExactDiagonalization"
        ExecDir = "/Users/chengyanlai/data"
    else:
        Cluster = "MacBookPro"
        MaxNumThreads = 1
        MaxWallTime = '48:0:0'
        qsub_cmd = ""
        SrcDir = "/Volumes/Files/GitRepo/ExactDiagonalization"
        ExecDir = os.path.join(SrcDir, "data")
