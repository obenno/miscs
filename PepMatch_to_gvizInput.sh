#! /bin/bash

# Input is the output of mapPep_to_CDS.py
# Convert its output to Gviz input

input = $1
output = $2

gffcompare $input

awk '
    BEGIN{
        print "chromosome\tstart\tend\twidth\tstrand\tfeature\tgene\ttranscript\tsymbol"
    }
    ' > $output

awk -F"\t" '
    NR==FNR&&$3=="transcript"{
        split($9,tmp," ");
        a[tmp[2]]=tmp[6]
    }
    NR>FNR&&$3=="exon"{
        split($9,tmp," ");
        tmp[2]=a[tmp[2]];
        print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"tmp[1]" "tmp[2]" "tmp[3]" "tmp[4]" "tmp[5]" "tmp[6]
    }
    ' gffcmp.combined.gtf gffcmp.combined.gtf |
    awk '
        {
            print $1"\t"$4"\t"$5"\t"$5-$4+1"\t"$7"\tCDS\t"$12"\t"$10"\t"substr($10,1,index($10,"_")-1)
        }
        ' |
    sed 's/\"//g; s/;//g' >> $output

## Change feature to "protein_coding"
sed -i 's/CDS/protein_coding/' $output

