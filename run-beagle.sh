#!/bin/bash
#    Copyright 2012, 2013, Peter Ralph and Graham Coop
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
