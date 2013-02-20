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
description = '''Merge two sets of genomes, randomly creating ibd segments within one set.
    Output only on lines shared by both sets.
    '''

import gzip
import random
import sys
from optparse import OptionParser
from math import isnan, exp, log
import pdb

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


def parsesite(snpsite):
    '''Read in a line from a snp file; 
    keeping the two haplotypes of an individual together.
    '''
    snpsite = snpsite.split()
    code = snpsite[0]
    name = snpsite[1]
    snpsite = [ snpsite[2*k:2*k+2] for k in range(1,len(snpsite)/2) ]
    return code,name,snpsite

def pairfiles(sf,tf,mapd):
    '''Read in sites from two files until their next matching site.
    Note that mapd enties should be in the form:
    mapd[ name ] = float(mappos), chrom, int(bp)
    '''
    sline = sf.readline()
    tline = tf.readline()
    while sline != '' and tline != '':
        scode,sname,ssite = parsesite(sline)
        tcode,tname,tsite = parsesite(tline)
        smap = mapd[sname][0]
        tmap = mapd[tname][0]
        sbp = mapd[sname][2]
        tbp = mapd[tname][2]
        if isnan( sbp ) or isnan( smap ):
            # print "skipping ", sname
            sline = sf.readline()
        elif isnan( tbp ) or isnan( tmap ):
            # print "skipping ", tname
            tline = tf.readline()
        elif sbp < tbp:
            # print "skipping ", sname, sbp, tbp
            sline = sf.readline()
        elif tbp < sbp:
            # print "skipping ", tname, sbp, tbp
            tline = tf.readline()
        else:
            # print "reading ", sname, tname, sbp, tbp
            yield scode, sname, ssite, tsite, mapd[sname]
            sline = sf.readline()
            tline = tf.readline()

parser = OptionParser(description=description)
parser.add_option("-i","--ibdfile",dest="ibdfile",help="name of input file to copy segments within (or '-' for stdin)",default="-")
parser.add_option("-j","--otherfile",dest="otherfile",help="name of other input file to merge with (or '-' for stdin)")
parser.add_option("-o","--outfile",dest="outfile",help="name of output file (or '-' for stdout)",default="-")
parser.add_option("-m","--mapfile",dest="mapfile",help="name of map file (or '-' for stdin)")
parser.add_option("-g","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="insert-ibd-segments.log")
parser.add_option("-a","--minwait",dest="minwait",help="Minimum length to wait before copying another segment (default 2cM)", default=2.0)
parser.add_option("-b","--maxwait",dest="maxwait",help="Maximum length to wait before copying another segment (default 3cM)", default=3.0)
parser.add_option("-x","--minlen",dest="minlen",help="Minimum length of segment (default 0.5cM)", default=0.5)
parser.add_option("-y","--maxlen",dest="maxlen",help="Maximum length of segment (default 10cM)", default=10.0)
parser.add_option("-u","--mutprob",dest="mutprob",help="Mutation probability per site (default=.002)", default=0.002)
parser.add_option("-s","--missprob",dest="missprob",help="Missing probability per site (default=.023)", default=0.023)
# parser.add_option("-v","--mutable",dest="mutable",help="Proportion of individuals to apply 5x higher mutation rates to.", default=.1)
(options,args) =  parser.parse_args()

ibdfile = fileopt(options.ibdfile, "r")
otherfile = fileopt(options.otherfile, "r")
outfile = fileopt(options.outfile, "w")
mapfile = fileopt(options.mapfile, "r")
logfile = fileopt(options.logfile, "w")
minwait = float(options.minwait)
maxwait = float(options.maxwait)
minlen = float(options.minlen)
maxlen = float(options.maxlen)
mutprob = float(options.mutprob)
missprob = float(options.missprob)
# mutable = float(options.mutable)

def random_seglen(maxlen=maxlen,minlen=minlen):
    '''random segment length to be copied -- uniform on a log scale'''
    return exp( log(maxlen/minlen)*random.uniform(0,1)+log(minlen) )


# Missing code
missing_code = '0'

# # TEST
# ibdfile = fileopt("short1_chr22.beagle", "r")
# otherfile = fileopt("short2_chr22.beagle", "r")
# mapfile = fileopt("/home/peter/projects/ibd/data/genetic_maps/marker.genetic22.gmap", "r")
# outfile = fileopt("INSERTED-test-output", "w")
# logfile = fileopt("-","w")
# minwait = 1.0
# maxwait = 5.0

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
ibdheader = ibdfile.readline().split()
otherheader = otherfile.readline().split()
# record the ids of the individuals 
ibdids = [ x for i,x in enumerate(ibdheader[2:]) if i%2 ]
otherids = [ x for i,x in enumerate(otherheader[2:]) if i%2 ]
outheader = ibdheader + otherheader[2:]
outfile.write(" ".join(outheader)+"\n")
# Number of samples
nibd = len(ibdheader)/2-1
nother = len(otherheader)/2-1
nsamps = nibd + nother
logfile.write("# Working on " + str(nsamps) + " samples.\n")
logfile.write("id1 id2 start end\n")

# create the iterator
stfiles = pairfiles(ibdfile, otherfile, mapd)
# Read until first matching site
code, name, ssite, tsite, mapsite = stfiles.next()
# Maintain lists of when each individual is next allowed to copy
#  and initialize them offset from the first position
swait = [ mapsite[0]+random.uniform(0,1) for x in ssite ]
# Also maintain the list of who is copying from where:
#  ...this is either an integer in range(nibd) or -1, meaning "not copying"
#  ... and a 0 or 1 for which haplotype
order = [ (-1,-1,-1) for x in range(nibd) ]
# and the list of how long to copy for
copyuntil = [ 0 for x in swait ]
# and keep list of segments in a dict
segments = {}
# ... and the set of individuals last copied to
lastcopied = [ [-1,-1] for x in range(nibd) ]
# output first line
outfile.write(" ".join( [code,name] + [" ".join(x) for x in ssite+tsite] )+"\n")
# for error checking
allnums = []

for code, name, ssite, tsite, mapsite in stfiles:
    # print str(mapsite[0]) + ":" + str(swait)
    # terminate segments
    for j in range(nibd):
        if copyuntil[j] <= mapsite[0]:
            i = order[j][0] 
            if (i,j) in segments:
                seg = segments.pop( (i,j) )
                # record site we began on and last site we wrote to file as start and end of segment
                logfile.write(" ".join( [ ibdids[i], ibdids[j], str(seg[0]), str(mapd[lastsite][0]) ] ) + "\n" )
            order[j] = -1,-1,-1
    # choose new sources and targets
    ibdable = [ j for j,x in enumerate(swait) if x <= mapsite[0] ]
    ibdable = random.sample( ibdable, len(ibdable)-len(ibdable)%2 )
    sourceable = ibdable[:len(ibdable)/2]
    targetable = ibdable[len(ibdable)/2:2*(len(ibdable)/2)]
    newseglens = [ random_seglen() for i in sourceable ]
    # print str(zip( sourceable, targetable, newseglens ))
    for i,j,x in zip( sourceable, targetable, newseglens ):
        # don't copy if we've recently copied these individuals
        if not i in lastcopied[j] and not j in lastcopied[i]:
            # copy from i to j for x length on randomly chosen haplotype
            order[j] = [i] + random.sample([0,1],2)
            lastcopied[i].append(j); lastcopied[i].pop(0)
            lastcopied[j].append(i); lastcopied[j].pop(0)
            copyuntil[j] = x+mapsite[0]
            swait[i] = swait[j] = copyuntil[j] + random.uniform(minwait,maxwait)
            # record chosen segments
            segments[i,j] = mapsite[0],copyuntil[j]
            # logfile.write(" ".join( [ ibdids[i], ibdids[j], str(mapsite[0]), str(copyuntil[j]) ] ) + "\n" )
    # find set of alleles
    alleles = set( sum(ssite+tsite,[]) )
    alleles.discard(missing_code)
    # record numbers of alleles (for error checking)
    if alleles:
        allnums.append(len(alleles))
    # don't output if it's monomorphic
    if alleles and len(alleles)>1:
        # copy over segments
        ssite = [ ( [x,y] if (i == -1) else ( [x,ssite[i][u]] if v else [ssite[i][u],y] ) ) for (i,u,v),(x,y) in zip(order, ssite) ]
        # for visualization debugging:
        # ssite = [ ( [x,y] if (i == -1) else ( [x,'*'] if v else ['*',y] ) ) for (i,u,v),(x,y) in zip(order, ssite) ]
        # insert mutations
        mutate = [ random.uniform(0,1) < mutprob for k in ssite ]
        ssite = [ ( [x,y] if not mut else ( [alleles.difference([x]).pop(),y] if random.choice([0,1]) else [x,alleles.difference([y]).pop()] ) ) for (x,y),mut in zip(ssite,mutate) ]
        # insert missings
        missing = [ random.uniform(0,1) < missprob for k in ssite ]
        ssite = [ ( [x,y] if not miss else [missing_code,missing_code] ) for (x,y),miss in zip(ssite, missing) ]
        # output
        # print "--", code, name, " mutated: ", sum(mutate), " missing: ", sum(missing), ssite, tsite
        outfile.write(" ".join( [code,name] + [" ".join(x) for x in ssite+tsite] )+"\n")
        lastsite = name
    # else:
        # print "Skipping monomorphic site ", mapsite[0]

# output remaining segments falling off the ends
for (i,j) in segments:
    seg = segments[ (i,j) ]
    logfile.write(" ".join( [ ibdids[i], ibdids[j], str(seg[0]), str(mapd[lastsite][0]) ] ) + "\n" )

for j in range(nibd):
    if copyuntil[j] <= mapsite[0]:
        i = order[j][0] 
        if (i,j) in segments:
            seg = segments[ i,j ]
            # record site we began on and last site we wrote to file as start and end of segment
        order[j] = -1,-1,-1

if any( [ x>2 for x in allnums ] ):
    alltab = [ allnums.count(x) for x in set(allnums) ]
    print "WARNING:"
    print "More than three alleles at some sites:", alltab

ibdfile.close()
otherfile.close()
outfile.close()
logfile.close()


