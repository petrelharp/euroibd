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
source("/home/peter/projects/genome/ibd-blocks-fns.R")
source("/home/peter/projects/genome/laplace-inversion-fns.R")  # for adjust.hist

load("all-blocks-winnowed.Rdata")

#######
# Compare endpoint distributions to the genetic map 
chrmap <- getchrmap(1:22)
chrmap$bp <- chrmap$bp + .chrbpstarts[chrmap$chr]
chrmap$map <- chrmap$map + .chrstarts[chrmap$chr]
# Rearrange so map position increases
chrmap <- chrmap[ order( chrmap$map, chrmap$bp ), ]
all( diff(chrmap$map)>=0 ) # TRUE

###########
# Compute number of overlapping blocks in different length categories along the genome
#  as well as looking up the number of nearby snps
lens <- c(1,2.5,4,6,8,10,Inf) 
windowsize <- 3 # cM
chrmap <- getchrmap(1:22)
chrmap$bp <- chrmap$bp + .chrbpstarts[chrmap$chr]
chrmap$map <- chrmap$map + .chrstarts[chrmap$chr]
# Rearrange so map position increases
chrmap <- chrmap[ order( chrmap$map, chrmap$bp ), ]
chrmap <- na.omit(chrmap)
all( diff(chrmap$map)>=0 ) # TRUE

olaps <- lapply( 1:(length(lens)-1), function (k) {
        z <- with( subset(blocks, maplen>lens[k] & maplen<=lens[k+1]), overlap( mapstart+.chrstarts[chrom], mapend+.chrstarts[chrom] ) )
        colnames(z) <- c("map","nblocks")
        z <- as.data.frame(z)
        z$chrom <- findInterval( z$map, .chrstarts )
        # how many nearby snps?
        z$nsnps <- NA
        for (chr in 1:22) {
            z$nsnps[z$chrom==chr] <- with(subset(z,chrom==chr), findInterval( map+windowsize, chrmap[chrmap$chr==chr,"map"] ) - findInterval( map-windowsize, chrmap[chrmap$chr==chr,"map"] ) )
            # normalize by length of window (affects blocks near the ends)
            z$nsnps[z$chrom==chr] <- z$nsnps[z$chrom==chr] / with(subset(z,chrom==chr), ( pmin(.chrstarts[chrom+1],map+windowsize) - pmax(.chrstarts[chrom],map-windowsize) ) ) 
        }
        z <- na.omit(z)
        resid.nblocks <- with(z, resid( loess( nblocks ~ nsnps, statistics="approximate", trace.hat="approximate", span=.2 )) )
        z$z <- scale(resid.nblocks)
        return( z )
    } )
names(olaps) <- paste(lens[1:(length(lens)-1)],"--",lens[-1],"cM",sep="")
# ylims <- range( unlist( lapply( olaps, function (x) quantile( (x[,"nblocks"]-mean(x[,"nblocks"]))/sqrt(var(x[,"nblocks"])), c(.001,.999) ) ) ) )
ylims <- range( unlist( lapply( olaps, function (x) quantile( x$z, c(.001,.999) ) ) ) )

#############
# Locate centromeres etc
centromeres <- read.table(paste(.mapdir,"centromeres-hg18",sep=""),header=TRUE)
centromeres$bpstart <- with(centromeres, bpstart + .chrbpstarts[chrom])
centromeres$bpend <- with(centromeres, bpend + .chrbpstarts[chrom])
centromeres$mapstart <- approx( x=chrmap$bp, y=chrmap$map, xout=centromeres$bpstart )$y
centromeres$mapend <- approx( x=chrmap$bp, y=chrmap$map, xout=centromeres$bpend )$y
if (FALSE) {  # Add all the inversions
    inversions <- read.csv(paste(.mapdir,"inversions.csv",sep=""),header=TRUE)
}
# interesting regions
#  for use on chrmap made to increase
#  bp2map <- function (x,chr) { with( chrmap[order(chrmap$bp),], map[ findInterval(x+.chrbpstarts[chr],bp) ] ) }
#   MHC: chromosome 6, from MOG to COL11A2, 29632000 to 33145000 =  43.04498 to 44.78392
#   Inversion on chr17 is from 44750000 for 900Kb = 63.64734 to 65.63818
#   Inversion p11q24 on chr8 is ??
#   Inversion: chr8:126845071-129872217 = 124.7991 130.6848
#   Inversions:  chr15:18870124-20099321 = NA-NA, chr15:28524207-30669954=13.24916-15.87232
#   Inversions:  chr16:1184606-1276384=NA-NA, chr16:14882692-15031240= 28.13764-28.21072, chr16:21485316-22620002=37.97415-38.41155
hilight <- list( MHC=c(6,43.04498,44.78392), "17q21.31"=c(17,63.64734, 65.63818), "8p"=c(8,13.59273,18.96307), "9p"=c(9,22.31770,22.49888), "15p"=c(15,13.28827,16.05908) ) 
#, "8q24.13"=c(8,124.4754,130.6775), "15q"=c(15,13.24916,15.87232), "16p"=c(16,28.13764,28.21072), "16q"=c(16,37.97415,38.41155)  )
# add centromeres
centros <- 1:22; names(centros) <- paste("centro",1:22)
hilight <- c( hilight, with(centromeres, lapply( centros, function (k) { c( chrom[k], mapstart[k]-.chrstarts[chrom[k]], mapend[k]-.chrstarts[chrom[k]]) } ) ) )
# # add inversions
# hilight <- c( hilight, with(inversions, lapply( 1:nrow(inversions), function (k) { c( chr[k], mapstart[k], mapend[k] ) } ) ) )

save(olaps, centros, centromeres, hilight, chrmap, lens, windowsize, ylims, file="overlaps.RData")
