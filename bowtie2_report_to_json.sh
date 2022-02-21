#! /usr/bin/env bash

set -euo pipefail

## This script is to convert bowtie2 mapping report to json format
## require gawk, jq installed

sampleName=$1
inputReport=$2
## use stdout for result

totalFragments=$(awk 'NR==1{print $1}' $inputReport)
concordant_unique=$(awk 'NR==4{print $1}' $inputReport)
concordant_unique_rate=$(awk 'NR==4{print $2}' $inputReport| sed 's/(//; s/)//')
concordant_multiple=$(awk 'NR==5{print $1}' $inputReport)
concordant_multiple_rate=$(awk 'NR==5{print $2}' $inputReport| sed 's/(//; s/)//')
discordant_unique=$(awk 'NR==8{print $1}' $inputReport)
non_concordant_unique=$(awk 'NR==13{print $1}' $inputReport)
non_concordant_multiple=$(awk 'NR==14{print $1}' $inputReport)
totalRate=$(awk 'NR==15{print $1}' $inputReport)

jq -n \
   --arg sampleName "$sampleName" \
   --arg totalFragments "$totalFragments" \
   --arg concordant_unique "$concordant_unique" \
   --arg concordant_unique_rate "$concordant_unique_rate" \
   --arg discordant_unique "$discordant_unique" \
   --arg non_concordant_unique "$non_concordant_unique" \
   --arg non_concordant_multiple "$non_concordant_multiple" \
   --arg totalRate "$totalRate" \
   '{sampleName: $sampleName, totalFragments: $totalFragments, uniqueReadPair: $concordant_unique, uniqueMappingRate: $concordant_unique_rate, discordant_unique: $discordant_unique, non_concordant_unique: $non_concordant_unique, non_concordant_multiple: $non_concordant_multiple, totalMappingRate: $totalRate}'
