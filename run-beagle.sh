#!/bin/bash
# Written 2013 by Peter Ralph and Graham Coop
# 
# contact: petrel.harp@gmail.com
#
#     To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty. 
# 
# You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 
# 
#
# Run beagle.

# remove previous lock file
if [[ -e run-beagle-finished ]]
then
    rm run-beagle-finished
fi

DATAFILE=$1
if ! [[ -e $DATAFILE ]]
then
    echo "File" $DATAFILE "does not exist."
    echo "Usage: ./run-beagle.sh DATAFILE [NRUNS]"
    exit 1
fi
if [[ $# -gt 1 ]]
then
    NRUNS=$2
else
    NRUNS=1
fi
if [[ -e /home/tmp ]]
then
    JAVA="java -Djava.io.tmpdir=/home/tmp"
else
    JAVA="java"
fi
echo "using" $JAVA

for k in $(seq $NRUNS)
do
    BFIX=$RANDOM
    nice -19 $JAVA -Xmx1500m -jar /usr/local/beagle/beagle.jar \
        unphased=$DATAFILE missing=0 fastibd=true seed=$BFIX gprobs=false \
        out=$BFIX &>$BFIX-beagle.log
done

touch run-beagle-finished

exit 0
