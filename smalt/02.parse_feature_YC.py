import multiprocessing as mp
import os
import sys
from dataclasses import dataclass
from typing import Tuple

import edlib
import pysam
import utils
from collections import defaultdict

# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
#########################################################################################
## 0. set class
#########################################################################################
@dataclass
class Seq:
    seq: str = ""
    qual: str = ""

@dataclass
class MatchedSeq:
    sam5: Seq
    umi: Seq
    bio: Seq
    sam3: Seq

#########################################################################################
## 1. 
#########################################################################################
def edit_barcode(barcode_seq: str,BC) -> Tuple[str, int]:
    """map the nearest sample barcode."""
    # correct barcode
    min_distance = len(barcode_seq)
    sample_matched = ""
    if min_distance == 0:
        return "", 0
    if BC == "P4":        
        for sample, barcode in BARCODE_DICT.items():
            result = edlib.align(barcode_seq, barcode[0], mode="NW", task="distance")
            dis = result["editDistance"]
            if dis <= min_distance:
                min_distance = dis
                sample_matched = sample
    if BC == "P5":
        for sample, barcode in BARCODE_DICT.items():
            result = edlib.align(barcode_seq, barcode[1], mode="NW", task="distance")
            dis = result["editDistance"]
            if dis <= min_distance:
                min_distance = dis
                sample_matched = sample        

    return sample_matched, min_distance

def check_sample_barcode(b1, b2, dis_cutoff=0.1):
    """check pair end barcode.
    dis_cutoff: proportion of mistach
    """
    # too strick filter
    if len(b1) == 0 or len(b2) == 0:
        return None
    org1, dis1 = edit_barcode(b1,"P4")
    org2, dis2 = edit_barcode(b2,"P5")
    if (
        org1 == org2 != ""
        and dis1 / len(b1) < dis_cutoff
        and dis2 / len(b2) < dis_cutoff
    ):
        return org1
    return None

def lookup_pos(pos, span_dict=None):
#     """
#     3K primers
#     [0:12] : sample barcode =>
#     [36:52] : UMI
#     [377:3381] : HMF =>
#     [3386:3398] : sample barcode <=
#     """
#     span_dict = {
#         range(0, 12): "sam5",
#         range(36 - 2, 52 + 2): "umi",
#         range(377, 3381): "bio",
#         range(3386, 3398): "sam3",
#     }
    """
    3KYC primers
    [6:12] : sample barcode =>
    [12:34] : P2 =>
    [36:52] : UMI  #UMI [X-2, X+3]
    [55:3059] : HMF =>
    [3059:3065] : sample barcode <=
    """

    span_dict = {
        range(6, 12): "sam5",
        range(36 - 2, 52 + 3): "umi",
        range(55, 3059): "bio",
        range(3059, 3065): "sam3",
    }
    for span in span_dict:
        if pos in span:
            return span_dict[span]
        

def parse_read(read):
    """parse bioseq to fastq."""
    #  dele = "-"
    dele = ""

    # check if there is MD tag
    # check whether query sequence is exist
    if not read.has_tag("MD") or (query_seq := read.query_sequence) is None:
        return None

    query_qual = read.query_qualities

    # init matched result object
    matched = MatchedSeq(
        sam5=Seq(seq="", qual=""),
        umi=Seq(seq="", qual=""),
        bio=Seq(seq="", qual=""),
        sam3=Seq(seq="", qual=""),
    )

    ref_pos_pointer = 0
    for query_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
        if ref_pos is not None:
            ref_pos_pointer = ref_pos
        feature_type = lookup_pos(ref_pos_pointer)
        if feature_type:
            getattr(matched, feature_type).seq += (
                query_seq[query_pos] if query_pos is not None else dele
            )
        if feature_type in ["bio", "umi"]:
            getattr(matched, feature_type).qual += (
                chr(query_qual[query_pos] + 33) ## 33??
                if query_pos is not None
                else dele
            )

    #  print(
    #  len(matched.bio.seq) > 2000,
    #  check_library_barcode(matched.lib5.seq, matched.lib3.seq),
    #  (sample := check_sample_barcode(matched.sam5.seq, matched.sam3.seq)) is not None
    #  )
    if (
        len(matched.bio.seq) > 2000
#         and check_library_barcode(matched.lib5.seq, matched.lib3.seq)
        and (sample := check_sample_barcode(matched.sam5.seq, matched.sam3.seq))
        is not None
    ):
        return (
            sample,
            matched.umi.seq,
            matched.umi.qual,
            matched.bio.seq,
            matched.bio.qual,
        )
    return None



#########################################################################################
## 1. run main scripts
#########################################################################################
if __name__ == "__main__":
    """
    2nd speed-limit step(2/2): 150 G data costs about 8 hours
    """
    SAMPLE_NAME = sys.argv[1]      # SAMPLE_NAME = "DSS2"
    INDIR=sys.argv[2]
    OUTDIR=sys.argv[3]
    SAMPLE_BARCODES = "/data/zhaolian/LineageTracing/SMALT/sampleBarcodes_"+SAMPLE_NAME+".txt"
    INPUT_BAM_FILE = INDIR+SAMPLE_NAME+".bam"
    BARCODE_DICT = {}
    with open(SAMPLE_BARCODES,"r") as f_barcode:
        for line in f_barcode.readlines():
            sample, P4, P5 = line.replace("\n", "").split("\t")
            BARCODE_DICT[sample] = P4, P5
    with pysam.AlignmentFile(
        INPUT_BAM_FILE, "rb", check_sq=False
    ) as infile:
        with open(SAMPLE_BARCODES,"r") as f_barcode:
            for line in f_barcode.readlines():
                sample = line.replace("\n", "").split("\t")[0]
#                 f[sample] = open(OUTDIR+str(sample)+".txt","a")
                if not os.path.exists(OUTDIR+str(sample)):
                    os.mkdir(OUTDIR+str(sample))
                
        res_write = defaultdict(str)
        for read in infile.fetch(until_eof=True):
            result = parse_read(read)
            if result: 
                if len(result[1]) >= 10 and len(result[1]) <= 30:
                    filename = OUTDIR+str(result[0])+"/"+str(result[1])+".fastq"
                    # if os.path.exists(filename):
                    # if filename in files
                    #     append_write = 'a' # append if already exists
                    # else:
                    #     append_write = 'w' # make a new file if not
                        
                    # outfile = open(filename, 'a')
                    # outfile.write(f"@{read.qname}\n{result[3]}\n+\n{result[4]}\n")
                    # outfile.close()
                    res_write[filename] = res_write[filename] + f"@{read.qname}\n{result[3]}\n+\n{result[4]}\n"
        for filename in res_write:
            with open(filename, 'w') as f:
                f.write(res_write[filename])