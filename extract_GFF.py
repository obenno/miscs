#! /usr/bin/env python

import os, sys
import argparse
import gffutils

parser = argparse.ArgumentParser(description='Extract Records from GFF/GTF')
parser.add_argument('-i', '--input',
                    help='input transcripts list')
parser.add_argument('-g', '--gff',
                    help='input GFF/GTF file')
parser.add_argument('-o', '--output',
                     help='extracted GFF/GTF file')
## or Output to Stdout
args = parser.parse_args()

## Create sqlite db from input GFF/GTF
with open(args.input, "r") as f:
    transcripts = f.read().splitlines()

db = gffutils.create_db(args.gff, ':memory:', force=True,
                        keep_order=True,
                        merge_strategy='merge',
                        sort_attribute_values=True)

out = open(args.output, "w")
SubFeatures = ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR')

outRecords = []

## Extract records according to transcript
##print("##gff-version 3\n##extracted from " + args.gff)
##print("##ID list filename: " + args.input)
outRecords.append("##gff-version 3\n##extracted from " + args.gff)
outRecords.append("##ID list filename: " + args.input)
for t in transcripts:
    ## Get gene record
    for gene in db.parents(t):
        ##print(gene)
        outRecords.append(gene)
    ## Get transcript itself
    outRecords.append(db[t])
    ## Get all children for transcript
    for i in db.children(t, featuretype=SubFeatures, order_by='start'):
        ##print(i)
        outRecords.append(i)

## Unique method is from stackoverflow:
## https://stackoverflow.com/questions/12897374/get-unique-values-from-a-list-in-python
outRecords_u = list(dict.fromkeys(outRecords))
for r in outRecords_u:
    print(r, file=out)

out.close()