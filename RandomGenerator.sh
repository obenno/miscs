#! /bin/bash

## Input is Reads_1.fq and Reads_2.fq
## Output is out_1.fq and out_2.fq

## Usage:
## RandomGenerator.sh Reads_1.fq.gz Reads_2.fq.gz 100 out_1.fq out_2.fq
## Extracted reads could also be 0~1, indicating percentage

num=`zcat $1 | wc -l`
reads=$(( $num/4 ))
if [[ $2 =~ "fq" ]] || [[ $2 =~ "fastq" ]];
then
    time=$3
    input1=$1
    input2=$2
    out1=$4
    out2=$5
else
    input1=$1
    time=$2
    out1=$3
fi

#rm -rf $3
zcat $input1 |
    awk -v reads="$reads" -v time="$time" '
    BEGIN{
        print "reads: "reads > "/dev/stderr"
        print "repeat or percentage: "time > "/dev/stderr"
        if(time > 1){
            rep=time
        }else if(time>0 && time<=1){
            rep=reads*time
        }else{
            print "Please provide extracted reads num or percentage" > "/dev/stderr"
            exit 1
        }
        for(i=1;i<=rep;i++){
            a[int(rand()*reads)*4-3]
        }
    }
    {k=NR; if(k in a){print; getline; print; getline; print; getline; print}}' > $out1

if [[ -n $out2 ]];
then
    zcat $input2 | awk '
    BEGIN{
        print "Pair end mode..." > "/dev/stderr"
    }
    NR==FNR{
        if(FNR%4==1){
            if($1~/\/1$/){
                a[substr($1,1,length($1)-2)]
            }else if($2!=""){
                split($2, tmp, ":")
                if(length(tmp)==4){
                    a[$1]
                }
            }else if($1~/^@SRR/){
                a[$1]
            }else{
                a[$1]
            }
        }
    }
    NR>FNR{
        if(FNR%4==1){
            k=0
            if($1~/\/2$/){
                if(substr($1,1,length($1)-2) in a){
                    k=1
                }
            }else if($1 in a){
                  k=1
            }
            if(k==1){
                print; getline; print; getline; print; getline; print
            }
        }
    }
    ' $out1 - > $out2
fi
