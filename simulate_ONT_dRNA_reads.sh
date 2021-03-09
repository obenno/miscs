#! /usr/bin/env bash


## This script is to use pbsim2 get simulated reads.
##
## profiles were studied from real reads, transcript
## depth were estimated separately and used for simulation

## Deps: bedtools>=2.28, minimap2, pbsim2, parallel

## Tested with bedtools v2.29.2
##             minimap2 2.17-r941
##             pbsim2   github e71f789

realReads=$1 # fastq format
refTrans=$2 # fasta format
threads=$3
outReads=$4 # fastq format

## First round run pbsim2, get profile from whole dataset
if [[ -e sample_profile_profile1.fastq ]];
then
    rm sample_profile_profile1.*
fi

outdir1=$(mktemp -d -p ./)
pbsim --depth 20 --prefix $outdir1/profiling \
      --sample-fastq $realReads \
      --sample-profile-id profile1 $refTrans

## Map real reads to reference transcriptome
minimap2 -t $threads -ax map-ont $refTrans $realReads |
    samtools sort -@ $threads > tmp.aln.bam

## Get depth estimation
bedtools genomecov -d -split -ibam tmp.aln.bam |
    bedtools groupby -g 1 -c 2,3 -o max,sum |
    awk '{print $1"\t"$3/$2}' > tmp.depth

## Run pbsim parallelly with distinct depth per transcript
outdir2=$(mktemp -d -p ./)
awk '$2>0{print "pbsim --depth "$2" --prefix '$outdir2/'"$1" --id-prefix "$1"_ --sample-profile-id profile1 <(fasta_formatter -t -i Gmax_a2v1.transcriptome.Chr01.fa | grep "$1" | awk '\''{print \">\"$1\"\\n\"$2}'\'')"}' tmp.depth |
    parallel -P $threads

cat $outdir2/*.fastq > $outReads
rm -rf $outdir1 $outdir2 tmp.depth tmp.aln.bam
