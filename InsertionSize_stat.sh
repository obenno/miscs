#! /bin/bash


## What the insert size used here is column 9 of sam input,
## which is actually fragment length (Oberved Template length)

## The calculation will based on first 5000 pairs/10000 reads

samtools view $1 | head -10000 |
    awk '
        BEGIN{n=0}
        $9>0{n+=1;a[n]=$9}
        END{
            sum=0
            sum_squares=0;
            for(i=1;i<=length(a);i++){sum+=a[i]};
            mean=sum/length(a)
            for(i=1;i<=length(a);i++){sum_squares+=(a[i]-mean)*(a[i]-mean)}
            div=sqrt(sum_squares/length(a))
            print "Mean: "mean
            print "Div: "div
        }
        '
