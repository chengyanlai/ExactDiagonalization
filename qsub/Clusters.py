#!/usr/bin/env python3
# coding=utf-8
import platform
import socket
import sys
import os

if platform.system() == "Linux":
  QSUB = True
  if socket.gethostname() == 'kagome.ucmerced.edu':
    Cluster = "Kagome"
    MaxNumThreads = 16
    MaxWallTime = '720:0:0'
    qsub_cmd = "qsub"
    if len(sys.argv) > 1:
      qsub_cmd = ' '.join([ qsub_cmd, ' '.join(sys.argv[1:]) ])
      if any("fast" in s for s in sys.argv):
        MaxNumThreads = 12
    SRC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization"
    EXEC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization/data"
  elif socket.gethostname() == 'merced.cluster':
    Cluster = "Merced"
    MaxNumThreads = 20
    MaxWallTime = '24:0:0'
    qsub_cmd = "qsub"
    if len(sys.argv) > 1:
      qsub_cmd = ' '.join([ qsub_cmd, ' '.join(sys.argv[1:]) ])
    SRC_DIR = "/home/clai24/GitRepo/ExactDiagonalization"
    EXEC_DIR = "/home/clai24/GitRepo/ExactDiagonalization/data"
  elif socket.gethostname() == 'braid.cnsi.ucsb.edu':
    Cluster = "Braid"
    MaxNumThreads = 20
    MaxWallTime = '672:0:0'
    qsub_cmd = "qsub -q batch"
    if len(sys.argv) > 1:
      qsub_cmd = ' '.join([ qsub_cmd, ' '.join(sys.argv[1:]) ])
    SRC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization"
    EXEC_DIR = "/home/chengyanlai/GitRepo/ExactDiagonalization/data"
elif platform.system() == "Darwin":
  QSUB = False
  Cluster = "MacBookPro"
  MaxNumThreads = 1
  MaxWallTime = '48:0:0'
  qsub_cmd = ""
  SRC_DIR = "/Volumes/Files/GitRepo/ExactDiagonalization"
  EXEC_DIR = os.path.join(SRC_DIR, "data")
