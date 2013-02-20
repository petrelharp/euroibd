#!/usr/bin/python
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
description = '''Randomly mix up the input genotypes.'''

import gzip
import random
import sys
from optparse import OptionParser
from math import isnan

def fileopt(fname,opts):
    '''Return the file referred to by fname, open with options opts;
    if fname is "-" return stdin/stdout; if fname ends with .gz run it through gzip.
    '''
    if fname == "-":
        if opts == "r":
            fobj = sys.stdin
        elif opts == "w":
            fobj = sys.stdout
        else:
            print "Something not right here."
    elif fname[len(fname)-3:len(fname)]==".gz":
        fobj = gzip.open(fname,opts)
    else:
        fobj = open(fname,opts)
    return fobj

parser = OptionParser(description=description)
parser.add_option("-l","--maxlen",dest="maxlen",help="maximum length (in cM) of segments (default=0.2)",default=0.2)
parser.add_option("-i","--snpfile",dest="snpfile",help="name of input file (or '-' for stdin)",default="-")
parser.add_option("-o","--outfile",dest="outfile",help="name of output file (or '-' for stdout)",default="-")
parser.add_option("-m","--mapfile",dest="mapfile",help="name of map file (or '-' for stdin)")
parser.add_option("-g","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="mixup-genomes.log")
(options,args) =  parser.parse_args()

maxlen = float(options.maxlen)
minlen = maxlen/2

# For BP instead of cM,
# maxbp = int( 1e6 * maxlen )
# minbp = int( 0.5 * 1e6 * maxlen )
# change 'mapsite[0]' to 'mapsite[2]' below
# and replace 'uniform' with 'randint'

snpfile = fileopt(options.snpfile, "r")
outfile = fileopt(options.outfile, "w")
mapfile = fileopt(options.mapfile, "r")
logfile = fileopt(options.logfile, "w")
# # TEST
# snpfile = fileopt("POPRES_chr22.beagle.gz", "r")
# outfile = fileopt("RANDOM_chr22.beagle.gz", "w")
# mapfile = fileopt("marker.genetic22.gmap", "r")

# Read in mapfile
mapheader = mapfile.readline().split()
mapd = {}
for line in mapfile:
    # of the form CHROM NAME MAPPOS BPPOS
    chrom,name,mappos,bp = line.split()
    if mappos == "NA":
        mappos = "NaN"  # translate to python
    mapd[ name ] = float(mappos), chrom, int(bp)

# set things up
header = snpfile.readline().split()
outfile.write(" ".join(header)+"\n")
# Number of samples
nsamps = len(header)/2-1
logfile.write("Working on " + str(nsamps) + " samples, cutting between " + str(minlen) + " and " + str(maxlen) + ".\n")

# The current order we are writing things out in
order = range(nsamps)
snpsite = snpfile.readline().split()
code = snpsite[0]
name = snpsite[1]
snpsite = [ snpsite[2*k:2*k+2] for k in range(1,len(snpsite)/2) ]
# Flush out the initial run of unknown positions
while isnan(mapd[ name ][0]):
    logfile.write("Skipping over missing position site " + str(name) + "\n")
    outfile.write(" ".join( [code,name] + [" ".join(x) for x in snpsite] )+"\n")
    snpsite = snpfile.readline().split()
    code = snpsite[0]
    name = snpsite[1]
    snpsite = [ snpsite[2*k:2*k+2] for k in range(1,len(snpsite)/2) ]
# Lines are of the form:
#    I id 474 474 584 584 608 608
#    M SNP_A-2174774 A G G G G A
mapsite = mapd[ name ]
# the list of future (bp) positions at which each site switches
switches = [ mapsite[0] + random.uniform(minlen,maxlen) for x in xrange(nsamps) ]
outfile.write(" ".join( [code,name] + [" ".join(x) for x in snpsite] )+"\n")

for line in snpfile:
    snpsite = line.split()
    code = snpsite[0]
    name = snpsite[1]
    snpsite = [ snpsite[2*k:2*k+2] for k in range(1,len(snpsite)/2) ]
    # Current position
    mapsite = mapd[ name ]
    switchnow = [ x<=mapsite[0] for x in switches ]
    nswitch = sum(switchnow)
    if nswitch==1:
        # choose another random one to switch
        switchnow[ random.choice( [ n for n,x in zip(xrange(nsamps),switchnow) if not x ] ) ] = True
    if nswitch>=1:
        # reorder these (by cyclic permutation)
        reordered = [ y for x,y in zip(switchnow,order) if x ]
        logfile.write("Switching " + str(nswitch) + " at " + str(mapsite[0]) + ": " + " ".join( map(str,reordered) ) + "\n")
        reordered = [ reordered.pop() ] + reordered
        j = 0
        for i in xrange(nsamps):
            if switchnow[i]:
                order[i] = reordered[j]
                j = j+1
    # rearrange and write out the site 
    snpsite = [ snpsite[i] for i in order ]
    outfile.write(" ".join( [code,name] + [" ".join(x) for x in snpsite] )+"\n")
    # record order everywhere?
    # logfile.write(str(mapsite[0]) + " " + " ".join( map(str,order) ) + "\n" )
    # choose the next set of places to switch
    switches = [ x if not y else mapsite[0]+random.uniform(minlen,maxlen) for x,y in zip(switches,switchnow) ]

snpfile.close()
mapfile.close()
outfile.close()
logfile.close()
