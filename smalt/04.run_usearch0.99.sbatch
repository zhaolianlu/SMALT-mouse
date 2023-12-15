#!/bin/bash
#SBATCH -J 3.04_99
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-25
#SBATCH --mem=2G #[M|G|T]
#SBATCH -t 1000:00:00 #[days-hh:mm:ss]
#SBATCH -o ./temp/%x-%j.log #%x: job_name, %j: job_id
#SBATCH -e ./temp/%x-%j.err #leave out -e flag to join error/log files 
##SBATCH --mail-user=zl.lu@siat.ac.cn
##SBATCH --mail-type=BEGIN,END #common valid events: NONE, BEGIN, END, FAIL, REQUEUE, ALL

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
################################
##3rd speed-limit step(4/4): 150 G data costs about 6 hours
#####################################################################################################
REF="/data2/kantian/LineageTracing/SMALT/reference/3k_HMFonly.mmi"

UMIFASTQ="/data2/kantian/LineageTracing/SMALT/3.IBD/02.fastq/"
UCDIR="/data2/kantian/LineageTracing/SMALT/3.IBD/04.uc0.99/"

declare -a INNAME=('152_N_E' '152_N_I' '153_N_E' '153_N_I' '155_N_E' '155_N_I' '156_N_E' '156_N_I' '19_N_E' '19_T1_E' '19_T1_I' '19_T3_E' '19_T3_I' '19_T4_E' '19_T4_I' '19_T5_E' '19_T5_I' '25_N_E' '26_N_E' '34_N_E' '37_N5_E' '37_N5_I' '38_N5_E' '38_N5_I' '40_N_E')
SAMPLE_NAME=${INNAME[$SLURM_ARRAY_TASK_ID]}
echo ${SAMPLE_NAME}

INDIR=${UMIFASTQ}${SAMPLE_NAME}
OUTDIR=${UCDIR}${SAMPLE_NAME}

mkdir -p $OUTDIR

cd $INDIR

for sample in *.fastq  
    do  
        echo $sample

        # if more than two sequence
        numLine=`cat $sample | wc -l`
        if [ $numLine -gt 9 ];then  ##  && UMI reads: >= 3
            describer=$(echo ${sample} | sed 's/.fastq//')  
            echo $describer  
            /data/kantian/software/usearch -cluster_fast $sample -id 0.99  -gapopen 3.0I/2.0E -gapext 1.0I/0.5E -match +2.0 -mismatch -20 -sizeout  -uc $OUTDIR/${describer}.uc
        else
            echo "no more than 3 reads"
        fi
        
    done



#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
