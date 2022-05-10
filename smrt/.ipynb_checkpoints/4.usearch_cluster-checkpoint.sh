#!/usr/bin/bash

WORKDIR="/data/zhaolian/LineageTracing/DSS/PacBio/3.umi"
OUTDIR="/data/zhaolian/LineageTracing/DSS/PacBio/4.umi_clustering"
mkdir -p $OUTDIR
cd $WORKDIR
echo "running 4.usearch_cluster.sh ............"

for sample in *_umi.fq
  do 
    echo $sample
    describer=$(echo ${sample} | sed 's/_umi.fq//')
    echo ${describer}
    usearch -cluster_fast $sample -id 0.95 -gapopen 3.0I/2.0E -gapext 1.0I/0.5E -match +2.0 -mismatch -20.0 -sizeout -uc $OUTDIR/${describer}_umi_UsearchClusters.uc -consout $OUTDIR/${describer}_umi_UsearchConsensus.fa -threads 20
  done

##usage ./4.usearch_cluster.sh
