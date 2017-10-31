#!/usr/bin/env python
# coding=utf-8
import os

def GenerateScript(QueueingSystem, Filename, JobName, Executable, FolderName, Nodes=1, NumCore=1, WallTime='48:00:00', Partition='', ProjectName='', MPI=0, PPN=1):
    JobScript = ["#!/bin/bash",]
    if QueueingSystem == "PBS":
        JobScript.append("#PBS -V")
        JobScript.append("#PBS -j oe")
        JobScript.append("".join(["#PBS -l nodes=", str(Nodes), ":ppn=", str(NumCore)]))
        JobScript.append("".join(["#PBS -l walltime=", WallTime]))
        JobScript.append("".join(["#PBS -d ", FolderName]))
        JobScript.append("".join(["#PBS -N ", JobName]))
        JobScript.append("cd $PBS_O_WORKDIR")
    elif QueueingSystem == "TORQUE":
        JobScript.append("#$ -V ")
        JobScript.append("#$ -j y")
        JobScript.append("#$ -S /bin/bash")
        JobScript.append("".join(["#$ -pe smp", str(NumCore)]))
        JobScript.append("".join(["#$ -l h_rt=", WallTime]))
        JobScript.append("".join(["#$ -N ", JobName]))
        JobScript.append("".join(["#$ -cwd"]))
        JobScript.append("".join(["cd ", FolderName]))
    elif QueueingSystem == "SLURM":
        if len(ProjectName):
            JobScript.append("".join(["#SBATCH --account=", ProjectName]))
        if len(Partition):
            JobScript.append("".join(["#SBATCH --partition=", Partition]))
        JobScript.append("".join(["#SBATCH --nodes=", str(Nodes)]))
        if MPI:
            JobScript.append("".join(["#SBATCH --ntasks=", str(MPI)]))
            JobScript.append("".join(["#SBATCH --ntasks-per-node=", str(PPN)]))
        JobScript.append("".join(["#SBATCH --cpus-per-task=", str(NumCore)]))
        JobScript.append("".join(["#SBATCH --time=", WallTime]))
        JobScript.append("".join(["#SBATCH --job-name=", JobName]))
        JobScript.append("".join(["#SBATCH --export=all"]))
    JobScript.append("".join(['export OMP_NUM_THREADS=', str(NumCore)]))
    JobScript.append('echo ----------------------------------------------')
    JobScript.append('echo "Job started on" `date`')
    JobScript.append('echo ----------------------------------------------')
    JobScript.extend(Executable)
    JobScript.append('echo ----------------------------------------------')
    JobScript.append('echo "Job ended on" `date`')
    JobScript.append('echo ----------------------------------------------')
    JobScript.append('exit 0')
    JobScript.append('\n')
    f = open(Filename, 'w')
    f.write( "\n".join(JobScript) )
    f.close()
