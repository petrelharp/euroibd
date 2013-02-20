# usage: python process_fibd.py score_threshold gap_threshold fibd_file1.gz fibd_file2.gz ... > outfile

# modified by plr from original by Browning
#  from http://faculty.washington.edu/sguy/beagle/ibd_and_hbd/ibd_and_hbd.html

import sys, gzip
from math import isnan
# import pdb

if len(sys.argv)<3:
  for i in range(50): print "*",
  print
  print >> sys.stderr, "This program takes Beagle fastIBD files from multiple seeds and combines overlapping IBD tracts for each pair of individuals."; print
  print >> sys.stderr, "The score_threshold parameter is provided so that one can apply a more stringent threshold than that used when generating the fastIBD files.  A recommended value is 1e-10."; print
  print >> sys.stderr, "Also, merge any adjacent blocks that are separated by no more than the distance gap_threshold."
  print >> sys.stderr, "Usage: python process_fibd.py score_threshold gap_threshold mapfile file1.fibd.gz file2.fibd.gz ... (pipe output to selected file)"; print
  print >> sys.stderr, "Output is one line per tract per pair of individuals, with index of first marker and last marker of tract (starting from 0, so first marker on chromosome is 0)"
  for i in range(50): print "*",
  print
  print >> sys.stderr, "After this in output are the number of segments that have been merged, the maximum score among those, and the minimum score among those."
  print
  raise SystemExit

nibsfiles = len(sys.argv) - 4
thresh = float(sys.argv[1])
# size of maximum gap to merge
gapthresh = float(sys.argv[2])

# Merge gaps:
# ... need to use map not snp number
mapfile = file(sys.argv[3],"r")
mapheader = mapfile.readline().split()
mapd = {}
j = 0 
for line in mapfile:
    # of the form CHROM NAME MAPPOS BPPOS
    chrom,name,mappos,bp = line.split()
    if mappos == "NA":
        mappos = "NaN"  # translate to python
    mapd[ j ] = float(mappos), chrom, int(bp), name
    j += 1
# For our purposes here replace nan's with min/max locations
minpos = min( [ x[0] for x in mapd.itervalues() if not isnan(x[0]) ] )
maxpos = max( [ x[0] for x in mapd.itervalues() if not isnan(x[0]) ] )
midbp = sum( [ x[2] for x in mapd.itervalues() ] ) / len(mapd)
for j in mapd:
    pos,chrom,bp,name = mapd[j]
    if isnan(pos): 
        if bp < midbp:
            pos = minpos
        if bp > midbp:
            pos = maxpos
        mapd[j] = pos,chrom,bp,name

results = {}
for i in range(nibsfiles):
  infile = gzip.open(sys.argv[i+4])
  for line in infile:
    # note: what is read is NOT the "length" but it is different (by one position) than the output.  
    # Beagle reports the starting (inclusive) and ending (exclusive) positions.
    (id1, id2, start, length, score) = line.split()
    start = int(start); end = int(length)-1; score = float(score)
    startpos = mapd[start][0]
    endpos = mapd[end][0]
    minscore = maxscore = score
    if score > thresh:
       continue
    if (id1,id2) in results:
       currentlist = results[(id1,id2)]
       markremove = [False for x in currentlist]
       for j in range(len(currentlist)):
                    x = currentlist[j]
                    ovlap = x[0]<=end and start<=x[1]
                    # is a gap?
                    gap = mapd[x[0]][0]<=endpos+gapthresh and startpos<=mapd[x[1]][0]+gapthresh
                    # is shorter than an adjacent segment?
                    # gaplen is max( start1-end2, start2-end1 ) since the min is negative
                    gap = gap and ( max( mapd[x[0]][0]-endpos, startpos-mapd[x[1]][0] ) < max( endpos-startpos, mapd[x[1]][0]-mapd[x[0]][0] ) )
                    if ovlap or gap:
                        # overlap
                        markremove[j] = True
                        start = min(x[0],start)
                        end = max(x[1],end)
                        startpos = mapd[start][0]
                        endpos = mapd[end][0]
                        minscore = min(minscore,x[3])
                        maxscore = max(maxscore,x[4])
       nsegs = 1 + sum( [ x[2] for x,r in zip(currentlist,markremove) if r] )
       if sum(markremove):
           results[(id1,id2)] = [x for x,r in zip(currentlist,markremove) if not r]
       results[(id1,id2)].append([start,end,nsegs,minscore,maxscore])
    else:
       results[(id1,id2)] = [[start,end,1,minscore,maxscore]]

# Note: changed 'minscore' to 'score'.
print "id1 id2 start end nsegs score maxscore"
for x in results:
  id1 = x[0]; id2 = x[1]
  for y in results[x]:
    print " ".join(map(str,x) + map(str,y))


