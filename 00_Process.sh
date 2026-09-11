#!/bin/bash
#2026-09-11
#Honggang Lyu
#bash 00_Process.sh

#prepare the working directories
mkdir -p data output

#=================================================
# MR analysis
#=================================================
echo "[INFO] TwoSampleMR analysis start!"
Rscript 01_MR.R
echo "[INFO] MR analysis Finished!!!"

#=================================================
# forest plot
#=================================================
echo "[INFO] Forest plot start!"
Rscript 02_Forest.R
echo "[INFO] Forest plot Finished!!!"
