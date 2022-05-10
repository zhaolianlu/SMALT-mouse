#!/usr/bin/perl -w
use strict;
if (@ARGV != 3) {
        die ".pl mem ppn path fastq_name sample_name\n";
}

#####################################################import parameters##################################################
#my $mem=$ARGV[0];
my $ppn=$ARGV[0];
my $node=$ARGV[1];
my $SAMPLENAME=$ARGV[2];
##########################Shell file#######################################################################################
open OUT,">8.run_iqtree_${SAMPLENAME}_2.sh";
my $std;
$std.= "#!/bin/sh\n";
$std.= "#PBS -N iqtree_$SAMPLENAME\n";
$std.= "#PBS -l walltime=240:00:00,nodes=$node:ppn=$ppn\n"; #mem=$mem,
#$std.= "#HSCHED -s Projects+software+Species\n";
#$std.= "#PPN limit $ppn\n";
$std.= "#PBS -o /data/zhaolian/LineageTracing/DSS/PacBio/scripts/8.run_$SAMPLENAME.out\n";
$std.= "#PBS -e /data/zhaolian/LineageTracing/DSS/PacBio/scripts/8.run_$SAMPLENAME.error\n";
#$std.= "#PBS -t \$a \n";
$std.= "#PBS -V \n";
$std.= "#PBS -S /bin/bash\n";
$std.= "\n";
$std.= "\n";

$std.= "PIPEDIR='/data/zhaolian/LineageTracing/DSS/PacBio/scripts'\n";
$std.= "WORKDIR='/data/zhaolian/LineageTracing/DSS/PacBio/8.tree_2'\n";
#$std.= "mkdir -p $WORKDIR";
$std.= "cd \$WORKDIR\n\n\n";

$std.= "echo PBS: current array id is \$PBS_ARRAYID \n";
$std.= "starttime=`date +'%Y-%m-%d %H:%M:%S'` \n";
	
$std.= "### 执行程序\n";
$std.= "iqtree -nt AUTO -s $SAMPLENAME.phy -o ref \n";
$std.= "###运行结束 \n\n";
$std.= "endtime=`date +'%Y-%m-%d %H:%M:%S'`\n";
$std.= "start_seconds=\$(date --date=\"\$starttime\" +%s);\n";
$std.= "end_seconds=\$(date --date=\"\$endtime\" +%s);\n";
$std.= "echo \"本次运行时间:\" \$((end_seconds-start_seconds))s  \n";
$std.= "echo \$starttime \n";
$std.= "echo \$endtime \n";
print OUT $std;
close OUT;

system "qsub 8.run_iqtree_${SAMPLENAME}_2.sh ";

#usage perl 8.iqtree_submit.pl 4 1 4N

