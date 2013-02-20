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
description = ''' 
  Turn hapmap legend/phased files into input required by beagle:
  Each panel has 3 files per chromosome/region. For example, for chr22 in the CEU panel there are the files
        
  /phased/genotypes_chr22_CEU_r21_nr_fwd_legend.txt
  /phased/genotypes_chr22_CEU_r21_nr_fwd_sample.txt
  /phased/genotypes_chr22_CEU_r21_nr_fwd_phased
   
  The _legend.txt file contains a legend detailing the rs id, base pair position, the allele coded 0
  and the allele coded 1 for each of the segregating SNPs e.g.
  
  NOTE: we have changed the rs to SNP_ID!!

  rs              position        0       1
  rs11089130      14431347        C       G
  rs738829        14432618        A       G
  rs915674        14433624        A       G
   
  The _phased contains the haplotypes arranged one haplotype per row.
   
  For the CEU and YRI the haplotypes are arranged as follows:
   
  row 1 - trio 1 parent 1 transmitted haplotype
  row 2 - trio 1 parent 1 untransmitted haplotype
  row 3 - trio 1 parent 2 transmitted haplotype
  row 4 - trio 1 parent 2 untransmitted haplotype
  row 5 - trio 2 parent 1 transmitted haplotype
  row 6 - trio 2 parent 1 untransmitted haplotype
  row 7 - trio 2 parent 2 transmitted haplotype
  row 8 - trio 2 parent 2 untransmitted haplotype
  .
  .
 The _sample.txt file contains an ordered list of individual ids that corresponds to the _phased file.  
 Individuals are arranged in the same order in the all/, phased/ and consensus/ datasets. 
 For the CEU and YRI panels the individual id's are arranged as follows but the _phased files do not 
 contain the haplotypes of the children as these can be infered from the parents haplotypes
  
 trio 1 parent 1
 trio 1 parent 2
 trio 2 parent 1
 trio 2 parent 2
 .
 .
 .
 trio 30 parent 1
 trio 30 parent 2
 trio 1 child
 trio 2 child
 .
 .
 .
 trio 30 child
 And output is one line per snp, like follows:

 I id 474 474 584 584 608 608 660 660 843 843 1311 1311 1725 1725 1890 1890 1923 1923
 M SNP_A-2314782 C C C C C C C C C C C C C C C C C C
 M SNP_A-1941632 G G C G G G C G G G C G G G C G G G
'''

from plr_utils import *
from optparse import OptionParser
# import pdb

parser = OptionParser(description=description)
parser.add_option("-l","--legendfile",dest="legendfile",help="name of input _legend file (or '-' for stdin)")
parser.add_option("-s","--samplefile",dest="samplefile",help="name of input _sample file (or '-' for stdin)")
parser.add_option("-p","--phasedfile",dest="phasedfile",help="name of input _phased file (or '-' for stdin)")
parser.add_option("-m","--mapfile",dest="mapfile",help="name of input map file to say which snps to keep (or '-' for stdin)")
parser.add_option("-o","--outfile",dest="outfile",help="name of output file (or '-' for stdin)")
(options,args) =  parser.parse_args()

legendfile = fileopt(options.legendfile, "r")
samplefile = fileopt(options.samplefile, "r")
phasedfile = fileopt(options.phasedfile, "r")
mapfile = fileopt(options.mapfile, "r")
outfile = fileopt(options.outfile, "w")

# Read in legendfile
legheader = legendfile.readline().split()
legend = [ line.split() for line in legendfile ]

# Read in mapfile
mapheader = mapfile.readline().split()
mapd = {}
for line in mapfile:
    # of the form CHROM NAME MAPPOS BPPOS
    chrom,name,mappos,bp = line.split()
    if mappos == "NA":
        mappos = "NaN"  # translate to python
    mapd[ name ] = float(mappos), chrom, int(bp)

# Read in samplefile
samples = []
for indid in samplefile:
    indid = indid.split()[0]
    samples += [ indid, indid ]

# Read in phased file
haplotypes = [ line.split() for line in phasedfile ]
nsnps = len(haplotypes[0])
if not nsnps == len(legend):
    raise ValueError("Something wrong?  Number of snps is not consistent.")
if len(haplotypes)<len(samples):
    print "Warning:", len(samples)-len(haplotypes), "samples not represented in haplotype file."
    samples = samples[:len(haplotypes)]
if not len(haplotypes) == len(samples):
    raise ValueError("Something wrong?  More haplotypes than samples.")

# begin output file
outfile.write("I id " + " ".join( samples ) + "\n")

for k in range(nsnps):
    snpid,bp,a0,a1 = legend[k]
    # recode alleles
    alleles = map( lambda i: [a0,a1][int(i)], [ x[k] for x in haplotypes ] )
    # check that this is in the collection of snps to output
    if mapd.get(snpid):
        # if any( [ bp < x[1] for x in legend[:k] ] ):
        #     pdb.set_trace()
        outfile.write( "M " + snpid + " " + " ".join(alleles) + "\n" )

outfile.close()
legendfile.close()
phasedfile.close()
samplefile.close()
