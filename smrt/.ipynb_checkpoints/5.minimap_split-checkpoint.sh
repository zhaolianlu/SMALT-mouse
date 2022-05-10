!#/usr/bin/bash

INDIR="/data/zhaolian/LineageTracing/DSS/PacBio/CCS5passes/2.fastq/"
OUTDIR="/data/zhaolian/LineageTracing/DSS/PacBio/CCS5passes/5.bam_split/"
REF="/data/zhaolian/LineageTracing/DSS/PacBio/reference/3k_HMFonly.mmi"
# PIPDIR="/data/zhaolian/LineageTracing/DSS/PacBio/CCS5passes/scripts/"
mkdir -p $OUTDIR
cd $INDIR
echo "running 5.minimap_split.sh ............"


for sample in *.fastq
  do 
    echo $sample
    describer=$(echo ${sample} | sed 's/.fastq//')
    echo $describer
    minimap2 -t 28 -A 4 -B 12 -O 10,15 -E 2,1 --score-N 0 --end-bonus 10 -a --MD -x map-pb $REF $sample > $OUTDIR/${describer}.sam
    samtools view -@ 10 -bS $OUTDIR/${describer}.sam > $OUTDIR/${describer}_unsorted.bam
    samtools sort -@ 10 -o $OUTDIR/${describer}.bam -O bam -T ${describer} $OUTDIR/${describer}_unsorted.bam
    samtools index $OUTDIR/${describer}.bam
    rm $OUTDIR/${describer}_unsorted.bam


  done

##usage nohup ./5.minimap_split.sh &
