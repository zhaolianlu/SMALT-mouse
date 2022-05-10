#!/bin/sh
#PBS -N iqtree_4N
#PBS -l mem=50gb,walltime=240:00:00,nodes=1:ppn=4
#PBS -o /data/zhaolian/LineageTracing/DSS/PacBio/scripts/8.run_4N.out
#PBS -e /data/zhaolian/LineageTracing/DSS/PacBio/scripts/8.run_4N.error
#PBS -V 
#PBS -S /bin/bash


PIPEDIR='/data/zhaolian/LineageTracing/DSS/PacBio/scripts'
WORKDIR='/data/zhaolian/LineageTracing/DSS/PacBio/8.tree'
cd $WORKDIR


echo PBS: current array id is $PBS_ARRAYID 
starttime=`date +'%Y-%m-%d %H:%M:%S'` 
### 执行程序
iqtree -nt AUTO -s 4N.phy 
###运行结束 

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
echo "本次运行时间:" $((end_seconds-start_seconds))s  
echo $starttime 
echo $endtime 
