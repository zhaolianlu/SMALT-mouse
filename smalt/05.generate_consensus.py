import pandas as pd
import os
import sys
from Bio import SeqIO, AlignIO
import pysam
# from Bio.Align.Applications import PrankCommandline
# import subprocess
from collections import Counter
import statistics
import re
import glob

#########################################################################################
## 1. get dominant cluster from each UMI group
#########################################################################################
def domiCluster(ucfile):
    uc = pd.read_csv(ucfile,delimiter="\t",header=None,usecols=[0,1,8])
    uc.columns=["a","cluster","name"]
    uc = uc.loc[uc["a"].isin(["C","H"])]
    cluster1 = uc["cluster"].value_counts().idxmax()
    uc=uc.loc[uc["cluster"]==cluster1]
    return(uc["name"].tolist())

#########################################################################################
## 2. get mutation allele table 
#########################################################################################
def dictBam(bamfile,dqnames):
    """
    i
    """
    df=pd.DataFrame({'qname': [],'mut_type': [],'mut_qual':[],'indel_type': [],'indel_qual':[]})
    with pysam.AlignmentFile(bamfile, "rb", check_sq=False) as infile:
        for read in infile.fetch(until_eof=True):
            if read.qname in dqnames:
                if read.has_tag("MD") and (query_seq := read.query_sequence) is not None and not read.is_secondary:
                    readDict = callMutation(read)
                    temp=pd.DataFrame({'qname': read.qname,'mut_type': [readDict[0]],'mut_qual':[readDict[1]],'indel_type': [readDict[2]],'indel_qual':[readDict[3]]})
                    df=pd.concat([df,temp],ignore_index=True)
    return(df)

def callMutation(read):
    """
    get mutations and indels from a read which has MD tag in sam/bam
    """

    mut_type = []
    mut_qual= []
    indel_type = []
    indel_qual =[]

    ref_index = read.reference_start-1
    query_index = -1
    I_index = [0]
    I_qual = []
    D_index = []
    D_qual = []

    
    
    for query_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True): #matches_only=True, # if True, no None on either side

        if query_pos is None: # Deletion
            ref_index += 1
            D_index.append(ref_index)
            D_qual.append(read.query_qualities[query_index])
#                     if D_index[-1]-D_index[-2] >1:
#                         indel_type.append(str(ref_index)+"D")
#                         indel_qual.append(read.query_qualities[query_index])

        elif ref_pos is None: # Insertion
            query_index +=1
            I_index.append(ref_index)
            I_qual.append(read.query_qualities[query_index])
            if I_index[-1]-I_index[-2] >1:
                indel_type.append(str(ref_index)+"I")
                indel_qual.append(read.query_qualities[query_pos])

        elif ref_base.islower():
            ref_index += 1
            query_index +=1
            mut_type.append(str(ref_pos)+ref_base.upper()+read.query_sequence[query_pos])
            mut_qual.append(read.query_qualities[query_pos])

        else:
            ref_index += 1
            query_index +=1
    D = concatDel(D_index, D_qual,"D")
    indel_type = indel_type + D[0]
    indel_qual = indel_qual + D[1]
    return(mut_type, mut_qual, indel_type, indel_qual)


def concatDel(A,B,type="D"):
    A.append(-1)
    l = len(A)
    i = 0
    j = 0
    inds = []
    OUT= []
    QUAL = []
    while i < l-2:
        while j < l-1:
            if not (A[j+1] - A[j] - 1):
                j += 1
            else:
                i = j + 1
                j = i
                inds.append(j)
    inds.insert(0, 0)
    for i in range(len(inds)-1):
        if not inds[i+1] - inds[i] - 1:
            OUT.append(str(A[inds[i]])+type)
            QUAL.append(B[inds[i]])
        else:
            OUT.append('{}-{}{}'.format(A[inds[i]], A[inds[i+1]-1],type))
            qual = B[inds[i]:inds[i+1]]
            qual_avg = int(sum(qual)/len(qual))
            QUAL.append(qual_avg)
    return(OUT, QUAL)

#########################################################################################
## 3. filter mutations and get mutation matrix
#########################################################################################
def ccSeq(sub, frequency_threshold, quality_threshold):
    """
    find cluster consensus mutations in each uc group
    """
    num = len(sub)
    ccMuts = []
    
    muts =flatList(sub["mut_type"])
    muts = removeBracket(muts, qual = False)
    mut_quals =flatList(sub["mut_qual"])
    indels =flatList(sub["indel_type"])
    indel_quals =flatList(sub["indel_qual"])
      
    fM = filterMuts(muts, mut_quals, num, frequency_threshold, quality_threshold)
    fID = filterMuts(indels, indel_quals, num, frequency_threshold, quality_threshold)
    ccMuts.append(fM)
    ccMuts.append(fID)
    ccMuts = flatList(ccMuts)
    return(ccMuts)

def removeBracket(s,qual=True):
    """
    ori: ["['2048TC']" '[]' "['1283CT', '1529GC']"]
    transformed: ['2048TC','1283CT','1529GC']
    """
    su = []
    if s is not None:
        for s_sub in s:
            s_sub = s_sub.replace("[","")
            s_sub = s_sub.replace("]","")
            s_sub = s_sub.replace("'","")
            s_sub = s_sub.split(",")
            for s_sub_sub in s_sub:
                if s_sub_sub != '':
                    s_sub_sub = s_sub_sub.replace(" ","")
                    if qual:
                        su.append(int(s_sub_sub))
                    else:
                        su.append(str(s_sub_sub))
    return(su)
    
def filterMuts(muts, mut_quals, num, frequency_threshold = 0.5, quality_threshold = 50):
    """
    filter mutations/indels with frequency >= frequency_threshold
    """
    cs = [ele for ele, cnt in Counter(muts).items() if cnt >= num*frequency_threshold]
    csMuts = []
    if len(cs) >0:
        for mutType in cs:
            quals = []
            for b,q in zip(muts, mut_quals):
                if b == mutType:
                    quals.append(q)
            if statistics.mean(quals) >= quality_threshold:
                csMuts.append(mutType)
    return(csMuts)


def flatList(Alist):
    flat_list = []
    for sublist in Alist:
        for item in sublist:
            flat_list.append(item)        
    return(flat_list)

#########################################################################################
## 4. get Consensus, binary & seqs
#########################################################################################
def parseDel(D,maxLen=200):
    for s in D:
        pos = re.findall("(\d+)",s)
        if len(pos) ==2:
            if int(pos[1])-int(pos[0]) > maxLen:
                return("LongDels")
def parseSubs(subs, ref):
    pos=[re.findall(r'[0-9]+',s)[0] for s in subs]
    types=[re.findall(r'[A-Z]+',s)[0] for s in subs]
    pos=[re.findall(r'[0-9]+',s)[0] for s in subs]
    types=[re.findall(r'[A-Z]+',s)[0] for s in subs]
    pos_m=[]
    for i in pos:
        pos_m.append(int(i))
    type_m=[]
    for i in types:
        type_m.append(i)
    s1=[0]*3004
    s2=ref
    for i in range(len(pos_m)):
        pointer = pos_m[i] #-1
        s1[pointer]=1
        s2 = s2[:pointer] + type_m[i][1] + s2[pointer + 1:]
    s0=''.join(map(str,s1))
    return(len(subs),s0, s2)   ## len(subs):number of mutations; s0: binary; s2:character sequence

def getConsensus(ssMuts, REF):
    ref = open(REF,"r").readlines()[1]
    
    types=[re.findall(r'[A-Z]+',s) for s in ssMuts]
    l = len(types)
    subs = []
    #### deletions ############
    D = []
    I = []
    for i in range(l):
        if types[i]==["D"]:
            D.append(ssMuts[i])
        elif types[i]==["I"]:
            I.append(ssMuts[i])
    ## parse mutations
    subs = [ele for ele in ssMuts if ele not in D and ele not in I]
#         if len(subs)>= numMuts:
    matrixALL = parseSubs(subs, ref)
    ## whether there is long deletions
    D_long = parseDel(D,maxLen=200)
    if D_long is None or "LongDels" not in D_long:
        longDel =0
    else:
        longDel =1
    return(matrixALL,longDel)


#########################################################################################
## run main scripts
#########################################################################################
if __name__ == "__main__":
    """
    generate consensus seqs
    """
    SAMPLE_NAME = sys.argv[1] ## 16N
    INBAM = sys.argv[2]+SAMPLE_NAME+"/" ## /data/zhaolian/LineageTracing/SMALT/1.CRC/
    INUC = sys.argv[3]+SAMPLE_NAME+"/" ## /data/zhaolian/LineageTracing/SMALT/1.CRC/
    OUTDIR = sys.argv[4] ## /data/zhaolian/LineageTracing/SMALT/1.CRC/
    numUC=3                  ## number of ccs reads in one consensus cluster
    mutFreq=0.6              ## export the mutation if mutation frequency >= mutFreq 
    quality_threshold = 50   ## mapping quality
    REF="/data/zhaolian/LineageTracing/SMALT/reference/3k_HMFonly.fasta"
    
    os.chdir(INUC)
    txtfiles = []
    for file in glob.glob("*.uc"):
        txtfiles.append(file)
        
    out_char=OUTDIR+SAMPLE_NAME+"_"+str(numUC)+"_"+str(mutFreq)+".fa"           ## output 1. fasta file
    out_phy=OUTDIR+SAMPLE_NAME+"_"+str(numUC)+"_"+str(mutFreq)+"_temp.phy"      ## output 2. phy file
    out_mut=OUTDIR+SAMPLE_NAME+"_"+str(numUC)+"_"+str(mutFreq)+"_mutations.txt" ## output 3. numUMI,numMut file
    out_nonexist=OUTDIR+SAMPLE_NAME+"_"+str(numUC)+"_"+str(mutFreq)+"_no_bam.txt" ## output 4. numUMI,numMut file
    with open (out_char, "w") as f_char, open(out_phy, "w") as f_phy, open(out_mut, "w") as f_mut, open(out_nonexist, "w") as f_nonexist:
        f_phy.write(str("ref")+" "+''.join(map(str,[0]*3004))+"\n")
        for file in txtfiles:
            umi = file.replace(".uc","")
            bamfile=INBAM+umi+".bam"
            if os.path.isfile(bamfile):
                ucfile=INUC+umi+".uc"
                dqnames=domiCluster(ucfile)
                if len(dqnames) >= numUC:
                    ucbam=dictBam(bamfile,dqnames)
    #                 sub=ucbam
                    ssMuts = ccSeq(ucbam, mutFreq, quality_threshold)
                    cc = getConsensus(ssMuts, REF)
                    f_char.write(">"+str(umi)+"\n"+str(cc[0][2])+"\n")
                    f_phy.write(str(umi)+" "+str(cc[0][1])+"\n")
                    f_mut.write(str(umi)+"\t"+str(len(dqnames))+"\t"+str(cc[0][0])+"\t"+str(cc[1])+"\n") ##umi numReads numMuts LongDel(1)
                else:
                    print(str(umi)+" less than "+str(numUC)+" reads")
            else:
                f_nonexist.write(str(umi)+"\n") ##umi numReads numMuts LongDel(1)
    phy=pd.read_csv(out_phy,sep = " ",header = None,names=["num","3004"])
    numR=phy.shape[0]
    phy=phy.rename({'num': str(numR)}, axis='columns')
    phy.to_csv(OUTDIR+SAMPLE_NAME+"_"+str(numUC)+"_"+str(mutFreq)+".phy",sep = " ",index=False)  ## output 2. phy file, add header
    os.remove(out_phy) 
