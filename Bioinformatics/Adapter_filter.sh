#!/bin/env bash

TRIMMOMATIC=$(find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null | head -1)

java -jar $TRIMMOMATIC SE -threads 8 -phred33 $1.R1.fq.gz $1.Ra.fq.gz HEADCROP:7 &> $1.trim1.log &

zcat $1.R2.fq.gz | seqtk seq -r | mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/TCCGATCT/) {match($0, "TCCGATCT"); if (RSTART < 20 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/TCCGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }' | 
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/.CCGATCT/) {match($0, ".CCGATCT"); if (RSTART < 20 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/.CCGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }' |
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/T.CGATCT/) {match($0, "T.CGATCT"); if (RSTART < 20 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/T.CGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }' |
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/TC.GATCT/) {match($0, "TC.GATCT"); if (RSTART < 20 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/TC.GATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }' | 
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/TCC.ATCT/) {match($0, "TCC.ATCT"); if (RSTART < 20 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/TCC.ATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }'  > $1.temp2

mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/^CCGATCT/) {match($0, "^CCGATCT"); if (length($0) > 100) {print substr($0,RSTART+14); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/^CCGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+14);j=0} }' $1.temp2 | 
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/^CGATCT/) {match($0, "^CGATCT"); if (length($0) > 100) {print substr($0,RSTART+13); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/^CGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+13);j=0} }' |
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/^GATCT/) {match($0, "^GATCT"); if (length($0) > 100) {print substr($0,RSTART+12); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/^GATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+12);j=0} }' |  seqtk seq -r | gzip > $1.R2a.fq.gz 

rm $1.temp2

wait 
wait

zcat $1.Ra.fq.gz |  mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/TCCGATCT/) {match($0, "TCCGATCT"); if (RSTART < 10 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/TCCGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }' | 
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/.CCGATCT/) {match($0, ".CCGATCT"); if (RSTART < 10 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/.CCGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }' |
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/T.CGATCT/) {match($0, "T.CGATCT"); if (RSTART < 10 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/T.CGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }' |
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/TC.GATCT/) {match($0, "TC.GATCT"); if (RSTART < 10 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/TC.GATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }' | 
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/TCC.ATCT/) {match($0, "TCC.ATCT"); if (RSTART < 10 && RSTART +15 < length($0)) {print substr($0,RSTART+15); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/TCC.ATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+15);j=0} }'  > $1.temp1

mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/^CCGATCT/) {match($0, "^CCGATCT"); if (length($0) > 100) {print substr($0,RSTART+14); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/^CCGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+14);j=0} }' $1.temp1 | 
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/^CGATCT/) {match($0, "^CGATCT"); if (length($0) > 100) {print substr($0,RSTART+13); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/^CGATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+13);j=0} }' |
mawk ' BEGIN {j=0}{if (NR%2 == 1) print $0 
else if (NR %2 == 0 && NR%4 != 0 && $0 ~/^GATCT/) {match($0, "^GATCT"); if (length($0) > 100) {print substr($0,RSTART+12); j=1} else print $0}
else if (NR %2 == 0 && NR%4 != 0 && $0 !~/^GATCT/) {j=0; print $0}
else if (NR %4 == 0 && j == 0) print $0
else if (NR %4 == 0 && j == 1) {print substr($0,RSTART+12);j=0} }' |  gzip > $1.R1a.fq.gz 

rm $1.temp1
