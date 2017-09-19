#!/usr/bin/env python
# coding=utf-8
import os

class cd:

  """Context manager for changing the current working directory"""

  def __init__(self, newPath):
    self.newPath = newPath

  def __enter__(self):
    self.savedPath = os.getcwd()
    os.chdir(self.newPath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self.savedPath)

def WriteQsubPBS(
  fname,
  Job_Name,
  EXAC_Name,
  Folder_Name,
  Nodes='1',
  NumCore=1,
  WallTime='48:00:00'):
  job_string = """#!/bin/bash
#PBS -V
#PBS -l nodes=%s:ppn=%d
#PBS -l walltime=%s
#PBS -j oe
#PBS -M clai24@ucmerced.edu
#PBS -d %s
#PBS -N %s
export OMP_NUM_THREADS=%d
cd $PBS_O_WORKDIR
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
%s
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
exit 0
""" % (Nodes, NumCore, WallTime, Folder_Name, Job_Name, NumCore, EXAC_Name)
  f = open(fname, 'w')
  f.write(job_string)
  f.close()
  return

def WriteQsubSBATCH(
  fname,
  Job_Name,
  EXAC_Name,
  Folder_Name,
  Nodes='1',
  NumCore=1,
  WallTime='48:00:00',
  partition='compute',
  ProjectName='TG-DMR150075'):
  f = open(fname, 'w')
  job_string = """#!/bin/bash
#SBATCH -A %s
#SBATCH --job-name=%s
#SBATCH --partition=%s
#SBATCH --nodes=%s
#SBATCH --ntasks-per-node=%d
#SBATCH --export=ALL
#SBATCH -t %s
#SBATCH --mail-user=clai24@ucmerced.edu
#SBATCH --mail-type=end
export OMP_NUM_THREADS=%d
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
%s
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
exit 0
""" % (ProjectName, Job_Name, partition, Nodes, NumCore, WallTime, NumCore, EXAC_Name)
  f.write(job_string)
  f.close()
  return

def WriteQsubSBATCH_MPI(
  fname,
  Job_Name,
  EXAC_Name,
  Folder_Name,
  Nodes='1',
  NumCore=1,
  WallTime='16:00:00',
  partition='standard',
  ProjectName='s17_cint'):
  f = open(fname, 'w')
  job_string = """#!/bin/bash
#SBATCH -A %s
#SBATCH --job-name=%s
#SBATCH --partition=%s
#SBATCH --nodes=%d
#SBATCH --ntasks=%d
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=%d
#SBATCH --time=%s
#SBATCH --export=ALL
export OMP_NUM_THREADS=%d
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
mpirun -n %d %s
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
exit 0
""" % (ProjectName, Job_Name, partition, Nodes, Nodes, NumCore, WallTime, NumCore, Nodes, NumCore, EXAC_Name)
  f.write(job_string)
  f.close()
  return

def WriteQsubSGE(
  fname,
  Job_Name,
  EXAC_Name,
  Folder_Name,
  Nodes='1',
  NumCore=1,
  WallTime='24:00:00'):
  f = open(fname, 'w')
  job_string = """#!/bin/bash
# The job is located in the current working directory.
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -pe smp %d
#$ -l h_rt=%s
# Job name.
#$ -N %s
# Send email to user.
#$ -M clai24@ucmerced.com
#$ -cwd
cd %s
export OMP_NUM_THREADS=%d
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
%s
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
exit 0
""" % (NumCore, WallTime, Job_Name, Folder_Name, NumCore, EXAC_Name)
  f.write(job_string)
  f.close()
  return
