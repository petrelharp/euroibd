# Written 2013 by Peter Ralph and Graham Coop
# 
# contact: petrel.harp@gmail.com
#
#     To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty. 
# 
# You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 
# 
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# usage:
#  make-hapmap-plus-popres.sh CHR

chr=$1
PFIX=$RANDOM

BASEDIR=/home/peter/projects/ibd/data
SCRIPTDIR=/home/peter/projects/genome

# pull out 60 Swiss:
COUNTRIES='"Belgium" "France" "Germany" "Hungary" "Ireland" "Italy" "Poland" "Portugal" "Spain" "Switzerland" "United.Kingdom" "Yugoslavia"'
# COUNTRIES=$(echo $COUNTRIES | tr " " "\n" | sort)
COUNTRYPATT=$(echo $COUNTRIES | sed -e 's/"\</\\(/g' | sed -e 's/\>"/\\)/g' | sed -e 's/ /\\|/g' | tr '.' ' ')
COUNTRY_CSPLIT=$(echo $COUNTRIES | tr '"' '/')

# make files with samples corresponding to countries
cat $BASEDIR/POPRES/european_labels/Euro-samples-info.tsv | cut -f 1,2 | grep "$COUNTRYPATT" | sort -k 2 | csplit -f country -z - $COUNTRY_CSPLIT 
CFILES=$(ls country??)

# country09 is the swiss
head -n 60 country09 | cut -f 2 | sed -e 's/\(.*\)/^\1$/' >remove.inds
KEEPCOLS=$( zcat ${BASEDIR}/POPRES/beagle-input/POPRES_chr${chr}.beagle.gz | head -n 1 | tr " " "\n" | grep -v -n -f remove.inds | cut -f 1 -d ':' | tr "\n" "," | sed -e "s/,$//")

# make temporary file without those
zcat ${BASEDIR}/POPRES/beagle-input/POPRES_chr${chr}.beagle.gz | cut -f $KEEPCOLS -d ' ' | gzip -c >temp.chr${chr}.beagle.gz

python ${SCRIPTDIR}/insert-ibd-segments.py -i ${BASEDIR}/HAPMAP/HAPMAP_PHASED_chr${chr}_CEU_r21_nr_all.beagle.gz -j temp.chr${chr}.beagle.gz -m ${BASEDIR}/genetic_maps/marker.genetic${chr}.gmap -o HAPMAP_POPRES_${PFIX}_chr${chr}.beagle.gz -g HAPMAP_POPRES_${PFIX}_chr${chr}.TRUE.ibd.gz $2
