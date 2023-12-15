#!/bin/bash
#SBATCH -J 3.05-9
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-25
#SBATCH --mem=5G #[M|G|T]
#SBATCH -t 1000:00:00 #[days-hh:mm:ss]
#SBATCH -o ./temp/%x-%j.log #%x: job_name, %j: job_id
#SBATCH -e ./temp/%x-%j.err #leave out -e flag to join error/log files 
##SBATCH --mail-user=zl.lu@siat.ac.cn
##SBATCH --mail-type=BEGIN,END #common valid events: NONE, BEGIN, END, FAIL, REQUEUE, ALL

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
################################
##1st speed-limit step(1/2): 150 G data costs about 24 hours
#####################################################################################################

declare -a INNAME=('152_N_E' '152_N_I' '153_N_E' '153_N_I' '155_N_E' '155_N_I' '156_N_E' '156_N_I' '19_N_E' '19_T1_E' '19_T1_I' '19_T3_E' '19_T3_I' '19_T4_E' '19_T4_I' '19_T5_E' '19_T5_I' '25_N_E' '26_N_E' '34_N_E' '37_N5_E' '37_N5_I' '38_N5_E' '38_N5_I' '40_N_E')
SAMPLE=${INNAME[$SLURM_ARRAY_TASK_ID]}
INBAM="/data2/kantian/LineageTracing/SMALT/3.IBD/03.bam/"
INUC="/data2/kantian/LineageTracing/SMALT/3.IBD/04.uc0.99/"
OUTDIR="/data2/kantian/LineageTracing/SMALT/3.IBD/05.consensus0.99/"
mkdir -p ${OUTDIR}
echo ${SAMPLE}

python 05.generate_consensus.py $SAMPLE $INBAM $INUC $OUTDIR

#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
