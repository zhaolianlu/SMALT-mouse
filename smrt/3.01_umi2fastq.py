import os
import sys

UMIDIR=sys.argv[1]
SAMPLE_NAME=sys.argv[2]

UMI = UMIDIR+SAMPLE_NAME+".umi"
FASTQ = UMIDIR+SAMPLE_NAME+"_umi.fq"

with open(UMI,"r") as in_f, open(FASTQ,"w") as out_f:
    for line in in_f.readlines():
        s = line.replace("\n", "").split("\t")
        if(len(s[2])>0):
            out_f.write(f"{s[0]}\n{s[2]}\n{str('+')}\n{s[3]}\n")
        
##usage python $PIPDIR/3.01_umi2fastq.py UMIDIR $describer