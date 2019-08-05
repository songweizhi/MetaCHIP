#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=11:59:00
#PBS -j oe
#PBS -m ae

# load needed modules (optional, specify if needed)
module load blast+/2.6.0
