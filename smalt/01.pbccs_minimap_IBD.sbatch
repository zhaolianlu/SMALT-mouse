#!/bin/bash
#SBATCH -J ccsIBD
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --array=0-1
##SBATCH --mem=100G #[M|G|T]
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
WORKDIR="/data/zhaolian/LineageTracing/SMALT/3.IBD/01.bam"
#REF="/data/zhaolian/LineageTracing/DSS/PacBio/reference/3k_BC_UMI_HMF_BC.mmi"
REF="/data/zhaolian/LineageTracing/Apc_YC/PacBio/reference/3kYC_BC_UMI_HMF_BC.mmi"

declare -a inbams=('/data2/zhaolian/LineageTracing/SMALT/IBD/1.rawdata/m64083_230609_233538.subreads.bam' '/data2/zhaolian/LineageTracing/SMALT/CRC19/1.rawdata/m64033_230612_032243.subreads.bam')
declare -a outname=('IBD' 'CRC19')
INBAM=${inbams[$SLURM_ARRAY_TASK_ID]}
SAMPLE=${outname[$SLURM_ARRAY_TASK_ID]} 
echo $SAMPLE
echo $INBAM

#####################################################################################################
# 1. ccs, by pbccs
mkdir -p $WORKDIR
cd $WORKDIR
ccs $INBAM ${SAMPLE}_stranded_ccs.bam --min-length=1000 --num-threads=10 --by-strand --min-passes=3 --min-rq=0.99 --report-file=${SAMPLE}_stranded_ccs.report
samtools fastq ${SAMPLE}_stranded_ccs.bam > ${SAMPLE}_stranded_ccs.fastq
rm ${SAMPLE}_stranded_ccs.bam
# 2. map, by minimap2
#minimap2 -d 3k_BC_UMI_HMF_BC.mmi 3k_BC_UMI_HMF_BC.fasta
#minimap2 -t ${threads} -O 2,12 -E 2,1 --score-N 0 --end-bonus 100 -a --MD -x map-pb ${ref} -
minimap2 -t 10 -A 4 -B 12 -O 10,15 -E 2,1 --score-N 0 --end-bonus 10 -a --MD -x map-pb $REF ${SAMPLE}_stranded_ccs.fastq> ${SAMPLE}_aln.sam
samtools view -bS ${SAMPLE}_aln.sam > ${SAMPLE}_unsorted.bam
samtools sort -@ 10 -o ${SAMPLE}.bam -O bam -T ${SAMPLE} ${SAMPLE}_unsorted.bam
samtools index ${SAMPLE}.bam
rm ${SAMPLE}_unsorted.bam

#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
