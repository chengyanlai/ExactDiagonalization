#!/usr/bin/env python
# coding=utf-8
<<<<<<< HEAD

=======
>>>>>>> dev
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

<<<<<<< HEAD

def WriteQsubSGE(
  fname,
  Job_Name,
  EXAC_Name,
  whentomail='n',
  mailto="chengyanlai@gmail.com",
  NumCore=1):
  f = open(fname, 'w')
  job_string = """#!/bin/bash
# The job is located in the current working directory.
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
# Job name.
#$ -N %s
# Send email to user.
#$ -M %s
# When will email be sent? begin, end, abort, suspend, none.
#$ -m %s
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
%s
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
exit 0
""" % (NumCore, Job_Name, mailto, whentomail, EXAC_Name)
  f.write(job_string)
  f.close()
  return


=======
>>>>>>> dev
def WriteQsubPBS(
  fname,
  Job_Name,
  EXAC_Name,
  Folder_Name,
<<<<<<< HEAD
  whentomail='n',
  mailto="chengyanlai@gmail.com",
  NodeName="kagome.rcc.ucmerced.edu",
  NumCore=1):
  f = open(fname, 'w')
  job_string = """#!/bin/bash
# Pass the current environment variables to the job.
#PBS -V
#PBS -l nodes=%s:ppn=%d
#PBS -l walltime=720:00:00
#PBS -j oe
#PBS -d %s
# Job name.
#PBS -N %s
# Change to current working directory (directory where qsub was executed)
# within PBS job (workaround for SGE option "-cwd")
=======
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
>>>>>>> dev
cd $PBS_O_WORKDIR
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
%s
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
exit 0
<<<<<<< HEAD
""" % (NodeName, NumCore, Folder_Name, Job_Name, EXAC_Name)
=======
""" % (Nodes, NumCore, WallTime, Folder_Name, Job_Name, NumCore, EXAC_Name)
  f = open(fname, 'w')
>>>>>>> dev
  f.write(job_string)
  f.close()
  return

def WriteQsubSBATCH(
  fname,
  Job_Name,
  EXAC_Name,
  Folder_Name,
<<<<<<< HEAD
  whentomail='n',
  mailto="chengyanlai@gmail.com",
=======
  Nodes='1',
>>>>>>> dev
  NumCore=1,
  WallTime='48:00:00',
  partition='compute',
  ProjectName='TG-DMR150075'):
  f = open(fname, 'w')
  job_string = """#!/bin/bash
#SBATCH -A %s
#SBATCH --job-name=%s
#SBATCH --partition=%s
<<<<<<< HEAD
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=%d
#SBATCH --export=ALL
#SBATCH -t %s

#SET the number of openmp threads
export OMP_NUM_THREADS=%d

=======
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
>>>>>>> dev
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
%s
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
exit 0
<<<<<<< HEAD
""" % (ProjectName, Job_Name, partition, NumCore, WallTime, NumCore, \
       EXAC_Name)
=======
""" % (NumCore, WallTime, Job_Name, Folder_Name, NumCore, EXAC_Name)
>>>>>>> dev
  f.write(job_string)
  f.close()
  return
