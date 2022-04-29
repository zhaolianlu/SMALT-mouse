##1st speed-limit step(1/2): 150 G data costs about 24 hours
#####################################################################################################
#BAMT="/data/zhaolian/LineageTracing/DSS/PacBio/Result-X101SC21114710-Z01-J002-B1-6-20220309/X101SC21114710-Z01_PacBio_Rawdata_XXXX/DSS-Tumor/all/FKDN220062840-1A/m64164_220304_192311.subreads.bam"
#BAMN="/data/zhaolian/LineageTracing/DSS/PacBio/Result-X101SC21114710-Z01-J002-B1-6-20220309/X101SC21114710-Z01_PacBio_Rawdata_XXXX/DSS-A/all/FKDN220062841-1A/m64164_220304_091036.subreads.bam"
WORKDIR="/data/zhaolian/LineageTracing/DSS/PacBio/1.bam"
REF="/data/zhaolian/LineageTracing/DSS/PacBio/reference/3k_BC_UMI_HMF_BC.mmi"
INBAM=$1
SAMPLE=$2
#####################################################################################################
# 1. ccs, by pbccs
mkdir -p $WORKDIR
cd $WORKDIR
ccs $INBAM ${SAMPLE}_stranded_ccs.bam --min-length=1000 --num-threads=10 --by-strand --min-passes=2 --min-rq=0.66 --report-file=${SAMPLE}_stranded_ccs.report
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
##usage nohup ./1.pbccs_minimap.sh /data/zhaolian/LineageTracing/DSS/PacBio/Result-X101SC21114710-Z01-J002-B1-6-20220309/X101SC21114710-Z01_PacBio_Rawdata_XXXX/DSS-A/all/FKDN220062841-1A/m64164_220304_091036.subreads.bam DSSN &
##usage nohup ./1.pbccs_minimap.sh /data/zhaolian/LineageTracing/DSS/PacBio/Result-X101SC21114710-Z01-J002-B1-6-20220309/X101SC21114710-Z01_PacBio_Rawdata_XXXX/DSS-Tumor/all/FKDN220062840-1A/m64164_220304_192311.subreads.bam DSST &
