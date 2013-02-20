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
# Usage:
#    mixup-by-regions.sh (chromosome number)
# Numbers of individuals we have:
#      Denmark         Latvia       Slovakia       Slovenia       Bulgaria 
#             1              1              1              1              2 
#       Finland         Norway        Ukraine        Albania         Cyprus 
#             2              2              2              3              4 
#        Russia         Greece         Turkey Czech Republic         Sweden 
#             5              6              7             12             12 
#       Austria        Romania    Netherlands        Hungary         Poland 
#            13             16             17             20             20 
#       Belgium        Ireland     Yugoslavia        Germany         France 
#            36             60             70             75             98 
#         Spain       Portugal          Italy United Kingdom    Switzerland 
#           128            129            193            394           1050 
#
# ... restrict to Hungary and above.

if [ -d /home/ibd/data ] 
then
    BASEDIR=/home/ibd/data
else 
    BASEDIR=$HOME/projects/ibd/data
fi
SCRIPTDIR=/home/peter/projects/genome

# Note must be in the same order as in Euro-samples-info.tsv.
# COUNTRIES='"France" "Spain" "Poland" "Hungary" "Yugoslavia" "Switzerland" "Italy" "Germany" "Portugal" "United.Kingdom" "Belgium" "Ireland"'
# COUNTRYPATT='\(France\)\|\(Spain\)\|\(Poland\)\|\(Hungary\)\|\(Yugoslavia\)\|\(Switzerland\)\|\(Italy\)\|\(Germany\)\|\(Portugal\)\|\(United Kingdom\)\|\(Belgium\)\|\(Ireland\)'
# COUNTRY_CSPLIT="/France/ /Spain/ /Poland/ /Hungary/ /Yugoslavia/ /Switzerland/ /Italy/ /Germany/ /Portugal/ /United.Kingdom/ /Belgium/ /Ireland/"
COUNTRIES='"Belgium" "France" "Germany" "Hungary" "Ireland" "Italy" "Poland" "Portugal" "Spain" "Switzerland" "United.Kingdom" "Yugoslavia"'
COUNTRYPATT='\(Belgium\)\|\(France\)\|\(Germany\)\|\(Hungary\)\|\(Ireland\)\|\(Italy\)\|\(Poland\)\|\(Portugal\)\|\(Spain\)\|\(Switzerland\)\|\(United Kingdom\)\|\(Yugoslavia\)'
COUNTRY_CSPLIT='/Belgium/ /France/ /Germany/ /Hungary/ /Ireland/ /Italy/ /Poland/ /Portugal/ /Spain/ /Switzerland/ /United.Kingdom/ /Yugoslavia/'

# make files with samples corresponding to countries
cat $BASEDIR/POPRES/european_labels/Euro-samples-info.tsv | cut -f 1,2 | grep "$COUNTRYPATT" | csplit -f country -z - $COUNTRY_CSPLIT 
CFILES=$(ls country??)

# get order of columns from input file
CHROM=$1
BFILE=POPRES_chr$CHROM.beagle
BDIR=$BASEDIR/POPRES/beagle-input
MAPFILE=$BASEDIR/genetic_maps/marker.genetic$CHROM.gmap
zcat $BDIR/$BFILE.gz | head -n 1 | tr ' ' '\n' > columns_all
for CFILE in $CFILES
do
    COLNUMS=$(cat $CFILE | cut -f 1 | grep -n -x -f - columns_all | tr ":" "\t" | cut -f 1 | tr "\n" "," | sed -e "s/,$//")
    echo $COLNUMS > colnums_$CFILE
    zcat $BDIR/$BFILE.gz | cut -f 1,2,$COLNUMS -d ' ' | python $SCRIPTDIR/mixup-genomes.py -i - -m $MAPFILE -g reorder.$CFILE.$BFILE.gz -o - | cut -f 3- -d ' ' >mixedup.$CFILE.$BFILE
done
zcat $BDIR/$BFILE.gz | cut -f 1,2 -d ' ' | paste -d ' ' - $(ls mixedup.*.$BFILE) | gzip -c > RANDOM.REGIONAL.$BFILE.gz
# save log files and whatnot.
tar -cvzf reorder.REGIONAL.$BFILE.tar.gz reorder.*.$BFILE.gz colnums_* columns_all $CFILES 

# clean up
for CFILE in $CFILES
do
    rm colnums_$CFILE
    rm $CFILE
    rm mixedup.$CFILE.$BFILE
    rm reorder.$CFILE.$BFILE.gz
done
rm columns_all
