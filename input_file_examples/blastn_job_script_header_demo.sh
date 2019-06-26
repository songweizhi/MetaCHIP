#!/bin/bash
#PBS -l nodes=1:ppn=3
#PBS -l mem=10gb
#PBS -l walltime=11:59:00
#PBS -j oe
#PBS -M email_address@gmail.com
#PBS -m ae

# optional, if blastn executable is not in your PATH environment
module load blast+/2.6.0
