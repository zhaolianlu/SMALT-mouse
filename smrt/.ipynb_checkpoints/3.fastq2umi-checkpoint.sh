##!/usr/bin/bash

INDIR="/data/zhaolian/LineageTracing/DSS/PacBio/2.fastq/"
OUTDIR="/data/zhaolian/LineageTracing/DSS/PacBio/3.umi/"
PIPDIR="/data/zhaolian/LineageTracing/DSS/PacBio/scripts/"
mkdir -p $OUTDIR
cd $INDIR
echo "running 3.fastq2umi ............"

for sample in *.fastq
  do
    echo $sample
    describer=$(echo $sample | sed 's/.fastq//')
    cat ${describer}.fastq | awk 'NR%4==1' > $OUTDIR/${describer}.umi
    python $PIPDIR/3.01_umi2fastq.py $OUTDIR $describer
  done

##usage ./3.fastq2umi.sh
