#!/bin/bash
#$ -N comp_sim2
#$ -q LT
#$ -pe smp 2
#$ -wd /Users/ssrivastva/icd_biom/sim_class_2/code/
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/
#$ -j y

module load stack/2020.1
module load r-matrix/1.2-17_gcc-9.2.0
module load r-glmnet/2.0-18_gcc-9.2.0
module load r-boot/1.3-23_gcc-9.2.0
module load r-stringr/1.4.0_gcc-9.2.0

R CMD BATCH --no-save --no-restore "--args 1 $SGE_TASK_ID" submit.R samp/comp2_$SGE_TASK_ID.rout
