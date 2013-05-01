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
# Various helper functions for working with ibd block files.
require(plyr)
require(colorspace)
# require(fields)

# base directory: takes the first one of these that it finds.
.basedir <- suppressWarnings( system("ls -d /home/ibd/ /home/peter/projects/ibd/", intern=TRUE, ignore.stderr=TRUE)[1] )
#.basedir <- "/home/peter/projects/ibd/"
# where are the genetic map files?
.mapdir <- paste(.basedir,"data/genetic_maps/",sep="")
# which set of maps to use?  All end in CHRNUM.gmap .
.mapbase <- "marker.genetic"
# where are the .fibd files?
.blocksdir <- paste(.basedir,"data/POPRES/ibdblocks",sep="")
# where is the list of individual sample information?
.geogdir <- paste(.basedir,"data/POPRES/european_labels/",sep="")
# where is the PC info
.pcadir <- paste(.basedir,"data/POPRES/pca_euro/",sep="")

# return the genetic map file
getchrmap <- function (chr, mapdir=.mapdir, mapbase=.mapbase) {
    do.call( rbind, lapply( chr, function (chrom) {
        read.table(paste(mapdir,mapbase,chrom,".gmap",sep=""), header=TRUE)
    } ) )
}

# Chromosome map lengths, from:
# .chrlens <- sapply(1:22, function (k) { chrmap <- getchrmap(k); c( max(chrmap$map,na.rm=TRUE) ) } )
.chrlens <- c( 262.00830, 244.03647, 209.97332, 197.48005, 189.58412, 174.70526, 172.35610, 156.52810, 146.28861, 160.40336, 143.88007, 156.40816, 119.34065, 104.83476, 110.93523, 118.68929, 119.61966, 105.74726, 92.57134, 83.62057, 54.79314, 55.59191 )
.chrstarts <- c(0, cumsum(.chrlens) )
# and in bp
# .chrbps <- sapply(1:22, function (k) { chrmap <- getchrmap(k); c( max(chrmap$bp,na.rm=TRUE) ) } ) 
.chrbps <- as.numeric( c( 247135059, 242663303, 199318155, 191167888, 180625439, 170747902, 158798338, 146264218, 140185941, 135284541, 134449982, 132276874, 114092980, 106356482, 100210760, 88684276, 78605474, 76115554, 63785051, 62376958, 46902240, 49522492 ) )
.chrbpstarts <- c(0,cumsum(.chrbps))

# Abbreviations
country.abbrev.list <- do.call( rbind, list( 
        c("Belgium", "BE"), 
        c("France", "FR"), 
        c("Austria", "AT"), 
        c("Bulgaria", "BG"), 
        c("Italy", "IT"), 
        c("Poland", "PL"), 
        c("Czech Republic", "CZ"), 
        c("Cyprus", "CY"), 
        c("Portugal", "PT"), 
        c("Denmark", "DK"), 
        c("Latvia", "LV"), 
        c("Romania", "RO"), 
        c("Germany", "DE"), 
        # c("Lithuania", "LT"), 
        c("Slovenia", "SI"), 
        # c("Estonia", "EE"), 
        # c("Luxembourg", "LU"), 
        c("Slovakia", "SK"), 
        c("Ireland", "IE"), 
        c("Hungary", "HU"), 
        c("Finland", "FI"), 
        c("Greece", "EL"), 
        # c("Malta", "MT"), 
        c("Sweden", "SE"), 
        c("Spain", "ES"), 
        c("Netherlands", "NL"), 
        c("United Kingdom", "UK"), 
        # c("Iceland", "IS"), 
        c("Norway", "NO"), 
        # c("Liechtenstein", "LI"), 
        c("Switzerland", "CH"), 
        c("Croatia", "HR"), 
        c("Montenegro", "ME"), 
        c("Turkey", "TR"), 
        c("Albania", "AL"), 
        c("Serbia", "RS"), 
        # c("Bosnia and Herzegovina", "BA"), 
        # c("Armenia", "AM"), 
        # c("Belarus", "BY"), 
        # c("Georgia", "GE"), 
        # c("Azerbaijan", "AZ"), 
        # c("Moldova", "MD"), 
        c("Ukraine", "UA"), 
        # c("Algeria", "DZ"), 
        # c("Lebanon", "LB"), 
        # c("Syria", "SY"), 
        c("Russia", "RU"), 
        c("Swiss French", "CHf"), 
        c("Swiss German", "CHd"), 
        c("Yugoslavia", "YU"),
        c("Bosnia", "BO"),     
        c("Croatia", "CR"),    
        c("England", "EN"),    
        c("Kosovo", "KO"),
        c("Macedonia", "MA"),
        c("Montenegro", "MO"),
        c("Scotland", "SC"),
        c("Serbia", "SR")
    ) )
country.abbrev.list <- country.abbrev.list[order(country.abbrev.list[,1]),]
countryabbrevs <- country.abbrev.list[,2]
names(countryabbrevs) <- country.abbrev.list[,1]
# Colors for countries
countrycols <- rainbow( length(countryabbrevs) )
names(countrycols) <- names(countryabbrevs)


# read in a (non-combined) fibd file
getblocks <- function (chr, prefix, blocksdir=.blocksdir, mapdir=.mapdir, mapbase=.mapbase, geogdir=.geogdir, namebase="POPRES_chr", endblocks=TRUE, remove.qc=TRUE, resolution="fine") {
    # If prefix matches "combined" or "winnowed", read in those files, respectively;
    # otherwise look for prefix.POPRES_chrXX, picking one at random if prefix is not specified.
    if (missing(prefix)) {
        # choose one at random
        filenames <- list.files(blocksdir, paste(".*\\.", namebase, chr, "\\..*fibd\\.gz", sep=""))
        prefix <- sub("\\.POPRES_chr.*","",sample(filenames,1))
        ftype <- NA
    } else {
        # ftype will be NA if it doesn't match
        ftype <- c("combined","winnowed")[pmatch(prefix, c("combined","winnowed"))]
    }
    if (is.na(ftype)) {
        # fibd as output by beagle
        filenames <- list.files(blocksdir, paste(prefix, "\\.", namebase, chr, "\\..*fibd\\.gz", sep=""), full.names=TRUE)
        if (length(filenames)!=1) { stop("Ambiguous/nonexistant block information: ", prefix, " and ", chr) }
        header <- FALSE
    } else {
        if (ftype == "winnowed") { ftype <- "combined.winnowed" }
        # combined as output by us
        filenames <- list.files(blocksdir, paste(namebase, chr, "\\..*", ftype, "\\.fibd\\.gz", sep=""), full.names=TRUE)
        if (length(filenames)!=1) { stop("Ambiguous/nonexistant block information: ", prefix, " and ", chr, " : ", filenames) }
        header=TRUE
    }
    readblocks( filenames[1], chr, mapdir=mapdir, mapbase=mapbase, geogdir=geogdir, header=header, endblocks=endblocks, remove.qc=remove.qc, resolution=resolution )
}

readblocks <- function (filename, chr, mapdir=.mapdir, mapbase=.mapbase, geogdir=.geogdir, header=TRUE, endblocks=TRUE, remove.qc=FALSE, resolution="coarse") {
    # Actually read in the block file
    # can use this for nonstandard filenames
    print(paste("Reading from", filename))
    if (header) {
        # ours have headers
        blocks <- read.table(filename, header=TRUE)
    } else {
        # as beagle outputs it
        blocks <- read.table(filename, col.names=c("id1","id2","start","end","score"))
    }
    blocks$chrom <- chr
    blocks <- addpos( blocks, chr, mapdir=mapdir, mapbase=mapbase, endblocks=endblocks )
    blocks <- addgeog( blocks, geogdir=geogdir, remove.qc=remove.qc, resolution=resolution )
    return(blocks)
}

getsampleinfo <- function(geogdir=.geogdir, resolution="fine", samplefile=switch(resolution,"coarse"="Euro-samples-info.tsv","Euro-samples-info-fine.tsv"), remove.qc=FALSE) {
    z <- read.table(paste(geogdir, samplefile, sep=""), header=TRUE)
    if (remove.qc) {
        z <- z[z$YESOK,]
        z$COUNTRY_SELF <- factor(z$COUNTRY_SELF)
    }
    return(z)
}

getcitydists <- function(geogdir=.geogdir) {
    as.matrix( read.table(paste(geogdir, "citypair.dists.tsv", sep=""), header=TRUE, row.names=1, check.names=FALSE) ) 
}

addlatlong <- function(blocks, ccols=setdiff(colnames(blocks)[grep("^country",colnames(blocks))],"countrypair"), geogdir=.geogdir, sampleinfo=getsampleinfo(geogdir=geogdir)) {
    # add lat and long info to data frame blocks
    suffixes <- gsub("^country","",ccols)
    ll <- ddply( sampleinfo, "COUNTRY_SELF", summarize, lat=mean(lat), long=mean(long) )
    for ( k in seq_along(ccols) ) {
        blocks[,paste("lat",suffixes[k],sep="")] <- ll$lat[match(blocks[,ccols[k]],ll$COUNTRY_SELF)]
        blocks[,paste("long",suffixes[k],sep="")] <- ll$long[match(blocks[,ccols[k]],ll$COUNTRY_SELF)]
    }
    return(blocks)
}

addgeog <- function(blocks, sampleinfo, citydists, idcols=grep("^id",colnames(blocks)), geogdir=.geogdir, remove.qc=FALSE, resolution="coarse") {
    # Add geographical information: country1, country2, city1, city2, and gdist
    #  ... optionally removing "mixed" or other ambiguous individuals
    # And also PCA info

    # expand columns idcols
    # and append these labels
    idcol.names <- gsub( "^id", "", colnames(blocks)[idcols] )

    if (missing(sampleinfo)) { sampleinfo <- getsampleinfo(geogdir=geogdir,resolution=resolution) }
    if (missing(citydists)) { citydists <- getcitydists(geogdir=geogdir) }

    # no longer a problem (?)
    # if (is.factor(blocks[,idcols[1]])) { stop("Uh-oh: id1 is a factor?") }

    # Add country information and city information
    for (k in 1:length(idcols)) {
        blocks[,paste("country",idcol.names[k],sep='')] <- sampleinfo$COUNTRY_SELF[match(blocks[,idcols[k]],sampleinfo$SUBJID)]
        blocks[,paste("city",idcol.names[k],sep='')] <- sampleinfo$CITY_SELF[match(blocks[,idcols[k]],sampleinfo$SUBJID)]
        # and "geographic id"
        blocks[,paste("gid",idcol.names[k],sep='')] <- sampleinfo$GEOGID[match(blocks[,idcols[k]],sampleinfo$SUBJID)]
    }
    countrycols <- match( paste("country",idcol.names,sep=""), colnames(blocks) )
    citycols <- match( paste("city",idcol.names,sep=""), colnames(blocks) )
    # add a "countrypair" variable
    if (length(idcols)==2) {
        countrypairs <- gsub(' ', '.', outer(levels(blocks$country1),levels(blocks$country2),paste,sep="-") )
        countrypairs[row(countrypairs)>col(countrypairs)] <- t(countrypairs)[row(countrypairs)>col(countrypairs)]
        blocks$countrypair <- factor( countrypairs[cbind(as.numeric(blocks$country1),as.numeric(blocks$country2))] )
    
        # Add geographic distance
        havelocation <- (blocks$country1 %in% dimnames(citydists)[[1]]) & (blocks$country2 %in% dimnames(citydists)[[1]])
        blocks$gdist <- NA
        blocks$gdist[havelocation] <- apply( as.matrix(blocks[havelocation,c("country1","country2")]), 1, 
                function (x) { citydists[x[1],x[2]] } 
            )
    
        # and PC distance
        pc1s <- sampleinfo[match(blocks$id1,sampleinfo$SUBJID),c("PC1","PC2")]
        pc2s <- sampleinfo[match(blocks$id2,sampleinfo$SUBJID),c("PC1","PC2")]
        blocks$pcdist <- sqrt( rowSums( (pc1s-pc2s)^2 ) )
    }

    if (remove.qc) {
        # removed mixed etc
        removeindivs <- with(sampleinfo, SUBJID[ !YESOK ])
        blocks <- subset( blocks, ! ( apply( blocks[,idcols], 1, function (x) any( x %in% removeindivs ) ) ) ) 
    }

    # make sure country{1,2} and city{1,2} have the same levels and drop unused ones
    countries <- unique( unlist( lapply( blocks[,countrycols], function (x) { levels(x)[table(x)>0] } ) ) )
    for (k in countrycols) {
        blocks[,k] <- factor(blocks[,k], levels=countries)
    }
    cities <- unique( unlist( lapply( blocks[,citycols], function (x) { levels(x)[table(x)>0] } ) ) )
    for (k in citycols) {
        blocks[,k] <- factor(blocks[,k], levels=cities)
    }

    return(blocks)
}

# Look at PCs
# helper function for scatterplotting:
# blah <- function (x,y,minlen=5,maxlen=Inf,alpha=.5) { with( subset(blocks.w.pcs, maplen>minlen & minlen<maxlen & countrypair==gsub(" ",".",paste(sort(c(x,y)),collapse="-"))), plot( pc1, pc2, col=adjustcolor( "black", alpha ), pch=20, xlim=quantile(pc1,c(.01,.99),na.rm=TRUE), ylim=quantile(pc2,c(.01,.99),na.rm=TRUE) ) ) }
# helper function for histogramming:
histpcs <- function (x,y,minlen1=0,maxlen1=2,minlen2=maxlen1,maxlen2=Inf) { 
    opar <- par(mfrow=c(1,2)); 
    cp <- gsub(" ",".",paste(sort(c(x,y)),collapse="-"))
    with( subset(blocks.w.pcs, maplen>minlen1 & maplen<maxlen1 & countrypair==cp), hist(pc1,breaks=100,freq=FALSE,main=paste(x,y)) )
    with( subset(blocks.w.pcs, maplen>minlen2 & maplen<maxlen2 & countrypair==cp), hist(pc1,breaks=100,freq=FALSE,add=TRUE,col=adjustcolor("red",0.5)) )
    with( subset(blocks.w.pcs, maplen>minlen1 & maplen<maxlen1 & countrypair==cp), hist(pc2,breaks=100,freq=FALSE,main=paste(x,y)) )
    with( subset(blocks.w.pcs, maplen>minlen2 & maplen<maxlen2 & countrypair==cp), hist(pc2,breaks=100,freq=FALSE,add=TRUE,col=adjustcolor("red",0.5)) )
    par(opar) 
} 

addpcs <- function (blocks, sampleinfo=getsampleinfo(remove.qc=TRUE), pccols=colnames(blocks)[grep("^pc[1-9]*\\..*",colnames(blocks))]) {
    # convert PC information, if observed PC info is (X + Y)/2, where X is the signal and Y is the mean for that country.
    #  assume columns containing pc info are named of the form pcX.idY
    pcnums <- unique( substr(gsub("\\..*","",pccols),start=3,stop=100) )
    for (pcnum in pcnums) {
        pcnumcols <- pccols[ grep(paste("^pc",pcnum,"\\.",sep=""),pccols) ]
        pcids <- gsub(".*\\.id","",pcnumcols)
        cmeans <- tapply( sampleinfo[,paste("PC",pcnum,sep='')], sampleinfo$COUNTRY_SELF, mean )
        blocks[,paste("pc",pcnum,sep='')] <- 0
        for (pcid in pcids) {
            pccol <- paste("pc",pcnum,".id",pcid,sep="")
            countrycol <- paste("country",pcid,sep="")
            blocks[,paste("pc",pcnum,sep='')] <- blocks[,paste("pc",pcnum,sep='')] + ( 2*blocks[,pccol] - cmeans[match(blocks[,countrycol],names(cmeans))] ) / length(pcids)
        }
    }
    return(blocks)
}

permgeog <- function ( blocks, sampleinfo, idcols=grep("^id",colnames(blocks)) ) {
    # Randomly permute (only) the countries indivs are assigned to
    #   and return the permuted country* columns
    # sampleinfo is as provided by getsampleinfo(remove.qc=TRUE)
    idcol.names <- gsub( "^id", "", idcols)
    # permute countries
    sampleinfo$COUNTRY_SELF <- sample( sampleinfo$COUNTRY_SELF)
    permcountries <- lapply( 1:length(idcols), function (k) {
            sampleinfo$COUNTRY_SELF[match(blocks[,idcols[k]],sampleinfo$SUBJID)]
        } )
    names(permcountries) <- paste("country",idcol.names,sep='')
    permcountries <- data.frame( permcountries )
    return(permcountries)
}

addpos <- function(blocks, chr, chrmap, mapdir=.mapdir, mapbase=.mapbase, endblocks=TRUE, colname="") {
    # Add map and bp position, length, etc to the blocks dataframe
    # If endblocks then include positions for blocks ending in unmapped areas
    if (!missing(chr) & is.data.frame(chr)) { chrmap <- chr }
    if (missing(chrmap)) { chrmap <- getchrmap(chr, mapdir=mapdir, mapbase=mapbase) }
    if (endblocks) {
        # replace telomere positions by maximal positions we have
        upt <- is.na(chrmap$map) & ( 1:nrow(chrmap) < nrow(chrmap)/2 )
        dnt <- is.na(chrmap$map) & ( 1:nrow(chrmap) > nrow(chrmap)/2 )
        chrmap$map[upt] <- min(chrmap$map,na.rm=TRUE)
        chrmap$map[dnt] <- max(chrmap$map,na.rm=TRUE)
    }
    # Note that beagle outputs a 0-based indexing system,
    #  so that marker 0 is the first one in the map file.
    blocks[,paste("mapstart",colname,sep='')] <- with(chrmap, map[blocks[,paste("start",colname,sep='')]+1] )
    blocks[,paste("mapend",colname,sep='')] <- with(chrmap, map[blocks[,paste("end",colname,sep='')]+1] )
    blocks[,paste("mapmid",colname,sep='')] <- (blocks[,paste("mapstart",colname,sep='')]+blocks[,paste("mapend",colname,sep='')])/2
    blocks[,paste("maplen",colname,sep='')] <- blocks[,paste("mapend",colname,sep='')]-blocks[,paste("mapstart",colname,sep='')]
    # And position in basepairs
    blocks[,paste("bpstart",colname,sep='')] <- with(chrmap, bp[blocks[,paste("start",colname,sep='')]+1] )
    blocks[,paste("bpend",colname,sep='')] <- with(chrmap, bp[blocks[,paste("end",colname,sep='')]+1] )
    blocks[,paste("bpmid",colname,sep='')] <- (blocks[,paste("bpstart",colname,sep='')]+blocks[,paste("bpend",colname,sep='')])/2
    blocks[,paste("bplen",colname,sep='')] <- blocks[,paste("bpend",colname,sep='')]-blocks[,paste("bpstart",colname,sep='')]
    # Length in number of markers
    blocks[,paste("nsnps",colname,sep='')] <- blocks[,paste("end",colname,sep='')]-blocks[,paste("start",colname,sep='')]+1
    return(blocks)
}

countpairs <- function (idnums, blocks, sampleinfo=getsampleinfo(remove.qc=TRUE)) {
    # count numbers of pairs between countries in idnums (or in blocks)
    if (missing(idnums)) {
        if (!missing(blocks)) {
            idnums <- unique( c(blocks$id1,blocks$id2) )
        } else {
            idnums <- sampleinfo$SUBJID
        }
    }
    sampleinfo <- droplevels( sampleinfo[ sampleinfo$SUBJID %in% idnums, ] )
    ncountries <- table( sampleinfo$COUNTRY_SELF )
    npairs <- outer( ncountries, ncountries, "*" )
    diag(npairs) <- choose( ncountries, 2 )
    pairnames <- gsub(" ", ".", outer(names(ncountries), names(ncountries), paste, sep="-"), fixed=TRUE )[upper.tri(npairs,diag=TRUE)]
    npairs <- npairs[upper.tri(npairs,diag=TRUE)]
    names(npairs) <- pairnames
    return(npairs)
}

subset.blocks <- function (blocks, idnums, only=FALSE, reorder.ids=FALSE) {
    # subset out blocks involving individuals given in idnums
    # if only is TRUE then only return blocks for which both are in idnums.
    # if reorder.ids then make sure first column is always something in idnums.
    if (is.data.frame(idnums)) { idnums <- as.matrix(idnums) }
    in.id1 <- blocks$id1 %in% idnums
    in.id2 <- blocks$id2 %in% idnums
    if (only) {
        usethese <- in.id1 & in.id2
    } else {
        usethese <-  in.id1 | in.id2 
        if (reorder.ids) {
            # make sure id1 and acommpanying info corresponds to something in idnums
            for (u in c("id","country","city","gid")) {
                u1 <- paste(u,"1",sep='')
                u2 <- paste(u,"2",sep='')
                x <- blocks[ !in.id1, u1 ]
                blocks[ !in.id1, u1 ] <- blocks[ !in.id1, u2 ]
                blocks[ !in.id1, u2 ] <- x
            }
        }
    }
    return( blocks[usethese,] )
}

getinds <- function (blocks) {
    # return list of individuals appearing in blocks
    if (is.factor(blocks$id1)) {
        inds <- unique( levels(blocks$id1), levels(blocks$id2) )
    } else {
        inds <- unique( blocks$id1, blocks$id2 )
    }
    return( inds )
}

match.blocks <- function (b1, b2) {
    # Find blocks in b1 and b2 that overlap,
    # both individual-wise and map-wise.
    # Return matrix of pairs (i,j) such that b1[i,] overlaps b2[j,]
    #  along with the *overlapping* length
    if (is.factor(b1$id1)) { 
        b1$id1 <- levels(b1$id1)[b1$id1]
        b1$id2 <- levels(b1$id2)[b1$id2]
    }
    matches <- ldply( 1:nrow(b1), function (k) {
                x <- b1[k,]
                mmm <- which( (b2$chrom==x$chrom) & (b2$id1 %in% c(x$id1,x$id2)) & (b2$id2 %in% c(x$id1,x$id2)) & (x$mapstart <= b2$mapend) & (x$mapend >= b2$mapstart) )
                if ( length(mmm)>0 ) {
                    lll <- pmin(b2$mapend[mmm],x$mapend)-pmax(b2$mapstart[mmm],x$mapstart)
                    return( data.frame(cbind(k,mmm,lll)) )
                } else { return(data.frame()) }
            } ) 
    if (!is.null(matches)) {
        names(matches)<-c("b1","b2","omaplen")
    }
    return(matches)
}

overlap <- function (starts, ends, weights=rep(1,length(starts))) {
    # Compute how much IBD overlaps each position on the genome
    # Returns an (nx2) matrix with 
    #    jumps[,1] = position
    #    jumps[,2] = number of blocks overlapping that position (until the next position)
    # Note that the total number of pairs of overlapping blocks is
    #   the sum over starting positions of the number of other blocks covering that position, i.e.
    #    sum( pmax(0, jumps[ c(TRUE,diff(jumps[,2])>0), 2 ] - 1 ) )
    locs <- sort(unique(c(starts,ends)))
    ups <- tapply( weights, factor(starts,levels=locs), sum )
    ups[is.na(ups)] <- 0
    downs <- tapply( weights, factor(ends,levels=locs), sum )
    downs[is.na(downs)] <- 0
    jumps <- cbind( locs, ups - downs )
    jumps[,2] <- cumsum(jumps[,2])
    return( jumps )
    # return( stepfun( jumps[,1], jumps[,2] ) )
}



get.gaps <- function (blocks) {
    # Return the set of gaps between blocks in the same indivs on the same chromsome
    gaps <- NULL
    for (chr in unique(blocks$chrom)) {
        verbose <- nrow(blocks)>1000
        if (verbose) { print(chr) }
        gblocks <- subset(blocks, chrom==chr)
        newgaps <- ddply( gblocks, c("id1","id2"), function(x) {
                if ( nrow(x)>1 ){
                    x <- x[order(x$start),]
                    x.id1 <- x$id1[1]
                    x.id2 <- x$id2[1]
                    x.chrom <- x$chrom[1]
                    x.lens <- x$maplen
                    x.scores <- x$score
                    nsnps <- 1 + x$start[-1] - x$end[-nrow(x)]
                    x <- as.vector( t(x[,c("mapstart","mapend")]) )
                    x <- x[-c(1,length(x))]
                    dim(x) <- c(2,length(x)/2)
                    x <- data.frame( t(x) )
                    names(x) <- c("mapstart","mapend")
                    x$id1 <- x.id1
                    x$id2 <- x.id2
                    x$chrom <- x.chrom
                    x$leftmaplen <- x.lens[-length(x.lens)]
                    x$rightmaplen <- x.lens[-1]
                    x$leftscore <- x.scores[-length(x.scores)]
                    x$rightscore <- x.scores[-1]
                    x$nsnps <- nsnps
                    return( x )
                } else { NULL }
            }, .progress=ifelse( verbose, "text", "none" ) )
        if (nrow(newgaps) > 0) { gaps <- rbind(gaps,newgaps) }
    }
    gaps$maplen <- gaps$mapend - gaps$mapstart
    gaps$mapmid <- gaps$mapstart + gaps$maplen/2
    gaps$leftmapstart <- gaps$mapstart - gaps$leftmaplen
    gaps$leftmapend <- gaps$mapstart
    gaps$rightmapstart <- gaps$mapend
    gaps$rightmapend <- gaps$mapend + gaps$rightmaplen
    return( gaps )
}


coverage <- function (blocks) {
    # return stepfun that is coverage of ibd blocks along genome
    startends <- rbind(
            cbind( blocks$mapstart + .chrstarts[blocks$chrom], +1 ),
            cbind( blocks$mapend + .chrstarts[blocks$chrom], -1 )
        )
    startends <- startends[ order(startends[,1]), ]
    startends[,2] <- cumsum(startends[,2])
    return(stepfun(x=c(startends[,1],sum(.chrlens[1:max(blocks$chrom)])),y=c(0,startends[,2],0)))
}



## PLOTTING FUNCTIONS

blocks.to.map <- function (blocks,chrstarts=.chrstarts) {
    # take a set of blocks and return a list of (start,end) vectors
    # so that they can be plotted on one line with segments()
    if (! "chrom" %in% names(blocks)) { stop("chrom not defined") }
    return( list( x0=blocks$mapstart+chrstarts[blocks$chrom], x1=blocks$mapend+chrstarts[blocks$chrom] ) )
}

plotblocks <- function(blocks, yvals, chroms=1:22, add=FALSE, lwd=3, yadj=0, xlab="chromosome", ylab="", chrspace=5, scale.lines=FALSE, xaxis=TRUE, xlim, ylim, ...) {
    # helper function:
    # plots all blocks in blocks,
    # each at position yvals
    # optionally shifted by yadj
    # putting chrspace between each chromosome
    chrstarts <- .chrstarts
    chrstarts[ -chroms ] <- NA
    chrstarts[ c(chroms,23) ] <- c(0,cumsum(.chrlens[chroms]+chrspace))
    chrends <- chrstarts-chrspace
    chrmids <- chrstarts[chroms] + .chrlens[chroms]/2
    segs <- blocks.to.map(blocks,chrstarts)
    # default for yvals is id1
    if (missing(yvals)) { 
        segs$y0 <- match( blocks$id1, unique(blocks$id1) )
    } else { 
        segs$y0 <- yvals 
    }
    segs$y0 <- segs$y0 + yadj
    if (!add) {
        if (missing(xlim)) { xlim <- range(chrstarts,na.rm=TRUE) }
        if (missing(ylim)) { ylim <- range(segs$y0,na.rm=TRUE) }
        plot( 0, type="n", xlab=xlab, xlim=xlim, ylim=ylim, xaxt="n", ylab=ylab, ... )
        rect( chrends, par("usr")[3], chrstarts, par("usr")[4], col=grey(.75) )
        if (scale.lines) {
            # add vertical lines for scale
            xscale <- ifelse( diff(xrange)>500, 100, 10 )
            abline(v=seq(xrange[1],xrange[2],by=xscale), col="grey", lty=2)
        }
        abline(v=c(chrstarts,chrends), col=grey(.8))
        if (xaxis) axis(side=1, at=chrmids, labels=chroms, tick=FALSE)
    }
    do.call("segments", c(segs, lwd=lwd, lend=1, list(...)) )
    return(invisible(segs))
}

plotpairs <- function (idnums, blocks, ...) {
    # Plot blocks for pairs of individuals, one per line.
    # each row gives a pair of individuals
    if (is.null(dim(idnums))) { dim(idnums) <- c(1,length(idnums)) }
    idnums <- unique(idnums)
    if (is.data.frame(idnums)) { idnums <- as.matrix(idnums) }
    blocks <- ldply( 1:nrow(idnums), function(k) { 
            x <- subset.blocks(blocks, idnums[k,], only=TRUE) 
            if (nrow(x)>0) { return( cbind(k, x) ) } else { return( NULL ) }  
        } )
    plotblocks( blocks, yvals=blocks[,1], ...  )
    return( invisible(blocks) )
}

plotmany <- function (idnums, blocks, ...) {
    # plot all blocks corresponding to idnums,
    # one per line, in color.
    blocks <- ldply( idnums, function (idnum) { subset.blocks(blocks, idnum, reorder.ids=TRUE) } )
    if( is.null(blocks) ) { stop("No such identifiers", idnums) }
    others <- unique( c(blocks$gid1, blocks$gid2 ) )
    nothers <- length(others)
    cols <- adjustcolor( rainbow_hcl(nothers, c=100, l=50), 0.75 )
    plotblocks( blocks, yvals=jitter(match( blocks$gid1, others ), amount=0.1), col=cols[match( blocks$gid2, others )], ... )
    if (length(idnums)<42) {
        abline( h=(1:(length(idnums)-1))+0.5, lty=2, col="grey" )
    }
}

plotindivs <- function(blocks, add=FALSE, offset=TRUE, allids=sort(unique(c(blocks$id1,blocks$id2))), ...) {
    # Plot each block twice, one row for each individual
    nids <- length(allids)
    id1 <- factor(blocks$id1,levels=allids)
    id2 <- factor(blocks$id2,levels=allids)
    y1 <- as.numeric(id1) + if (offset) { (as.numeric(id2)/(nids+1)-0.5)*.6 } else { 0 }
    y2 <- as.numeric(id2) + if (offset) { (as.numeric(id1)/(nids+1)-0.5)*.6 } else { 0 }
    plotblocks(blocks, yvals=y2, col=rainbow_hcl(nids,c=90,l=80)[as.numeric(id1)], ylim=c(.35,nids+.65), add=add, ... )
    plotblocks(blocks, yvals=y1, col=rainbow_hcl(nids,c=90,l=80)[as.numeric(id2)], add=TRUE, ...)
}

plotindiv <- function (idnum, blocks, cols, chroms=1, ...) {
    # Plot all segements shared with a given individual,
    # with one other country per line,
    # colored by individual by default
    op <- par(mar=c(5,8,4,2)+.1)
    blocks.id <- subset.blocks(blocks, idnum, reorder.ids=TRUE)
    if (missing(cols)) {
        cols <- colorize(blocks.id$id2)
    }
    plotblocks(blocks.id, jitter(as.numeric(blocks.id$country2)), chroms=chroms,yaxt="n",main=paste(blocks.id$country1[1],idnum),col=cols,ylab="",...)
    axis(2,at=1:nlevels(blocks.id$country2),labels=levels(blocks.id$country2),las=1)
    par(op)
    return(invisible(blocks.id))
}

## For plotting gaps:

plotgapsegs <- function (x, ...) {
    names(x)[names(x)%in% c("mapstart","mapend")] <- c("gapmapstart","gapmapend")
    names(x)[names(x)%in% c("leftmapstart","leftmapend")] <- c("mapstart","mapend")
    plotblocks(x, col="red", ...)
    names(x)[names(x)%in% c("mapstart","mapend")] <- c("leftmapstart","leftmapend")
    names(x)[names(x)%in% c("rightmapstart","rightmapend")] <- c("mapstart","mapend")
    plotblocks(x, col="blue", add=TRUE, ...)
    names(x)[names(x)%in% c("mapstart","mapend")] <- c("rightmapstart","rightmapend")
    names(x)[names(x)%in% c("gapmapstart","gapmapend")] <- c("mapstart","mapend")
    plotblocks(x,add=TRUE, ...)
}



##
# Computation

cor.counts <- function (blocks, distbins, countries=levels(blocks$country1), indivs=sort(unique(c(blocks$id1,blocks$id2))), countrymatch="block" ) {
    # For each individual, each country in countries, and each distance d in distbins,
    # compute the number of pairs of blocks within distance d of eachother
    # that are both of that country, 
    # and the total number within distance d of eachother of any country.
    #
    # countrymatch is "block" or "indiv" to look for blocks nearby of the same country as either
    #   the block in question or the individual in question
    #
    # Returns a table whose columns are
    # indiv country1 distance total [ ... countries ... ]
    #   that tabulates the counts of pairs in each category.
    # distbins <- exp( seq( log(.5), log(32), len=30 ) )
    # countries <- names(sort(table(blocks$country1),decreasing=TRUE))[1:8]
    # below we end up getting country codes not character strings...
    if (is.character(countries)) { 
        countrycodes <- match( countries, levels(blocks$country1) ) 
    } else if (is.numeric(countries)) {
        countrycodes <- countries
    } else {
        stop("I don't understand countries: ", countries)
    }
    # Get the big and little bins.  Note that cut() uses left-open intervals,
    # so including 0 as the smallest cutpoint will turn distance 0 into NA.
    if (all(is.finite(distbins[distbins>0]))) distbins <- c(distbins, Inf)
    if (all(distbins >= 0)) distbins <- c(-Inf,distbins)
    allcounts <- sapply( indivs, function (idnum) {
                x <- subset.blocks(blocks, idnum, reorder.ids=TRUE)   # 35 seconds
                # Note if two blocks are (a,b) and (c,d) then:
                #   the gap between them is either (c-b) or (a-d), whichever is positive,
                #   or they overlap if both are negative.
                # on the same chromosome?
                # matrix of block-block distances
                bbD <- outer(x$mapstart,x$mapend,"-")
                # on the same chromosome?
                bbD[ !outer(x$chrom,x$chrom,"==") ] <- NA
                bbD <- pmax(bbD,t(bbD))  # negative numbers is extent of overlap
                bbD <- cut( bbD[upper.tri(bbD,diag=FALSE)], distbins )  # discretize
                # do countries match?
                if (pmatch(countrymatch,c("block","individual"))==1) {
                    M <- outer(x$country2,x$country2,function(x,y) ifelse(x==y, x, NA))
                } else {
                    M <- outer(x$country2,x$country1,function(x,y) ifelse(x==y, x, NA))
                }
                return(  cbind( total=table(bbD), table( bbD, factor( M[upper.tri(M,diag=FALSE)], levels=countrycodes ) ) ) )
            } )
    dim(allcounts) <- c( length(distbins)-1, length(countries)+1, length(indivs) )
    dimnames(allcounts) <- list( distbins[-1], c("totals",countries), indivs )
    return(allcounts)
}

# Triples


polarized <- function (triples, varname="rate", pattern=c("aab","baa")[1], minsamples=20) {
    # From the data.frame as above,
    # return the matrix yy such that
    # either
    #  pattern="aab" => yy[a,b] = number of a,a,b triple blocks / number of triples a,a,b  
    # or
    #  pattern="baa" => yy[a,b] = number of b,a,a triple blocks / number of triples b,a,a  
    # ... so yy[a,b] increases with migration from a->b
    if (pattern=="aab") {
        yy <- apply( expand.grid(names(nsamples)[nsamples>minsamples],names(nsamples)[nsamples>minsamples]), 1, 
            function (x) {
                ind <- sort(x)   # sorts into decreasing order
                triples[ triples$countryA==x[1] & triples$countryB==ind[1] & triples$countryC==ind[2], varname ]
            } )
    } else if (pattern=="baa") {
        yy <- apply( expand.grid(names(nsamples)[nsamples>minsamples],names(nsamples)[nsamples>minsamples]), 1, 
            function (x) {
                triples[ triples$countryA==x[2] & triples$countryB==x[1] & triples$countryC==x[1], varname ]
            } )
    } else {
        stop("Pattern", pattern, "not recognized.")
    }
    dim(yy) <- rep(sum(nsamples>minsamples),2)
    dimnames(yy) <- list(names(nsamples)[nsamples>minsamples],names(nsamples)[nsamples>minsamples])
    return(yy)
}

three.rate <- function (triples, varname="rate", minsamples=20) {
        yy <- sapply( names(nsamples)[nsamples>minsamples],
            function (x) {
                triples[ triples$countryA==x & triples$countryB==x & triples$countryC==x, varname ]
            } )
        names(yy) <- names(nsamples)[nsamples>minsamples]
        return(yy)
}


# Computation

conditional.means <- function (x,n,tabx=table(x)) {
    # given a vector of obsevations x,
    # return a vector of length n
    # whose (k+1)th element is the mean of a k-size-biased pick from x, 
    # i.e. E[ X(X-1)...(X-k) ]/E[ X(X-1)...(X-k+1) ]
    # and is therefore constant in expectation if x is Poisson
    cmeans <- sapply( setdiff(n,0), function (k) sum( apply(outer(as.numeric(names(tabx)),0:k,"-"),1,prod)*tabx )/sum( apply(outer(as.numeric(names(tabx)),0:(k-1),"-"),1,prod)*tabx ) ) 
    if (0 %in% n) { cmeans <- c( sum(tabx*as.numeric(names(tabx)))/sum(tabx), cmeans ) }
    return(cmeans)
}


# Sample without undesired behavior at n=0
tsample <- function (x,n) { if(length(x)==1 & n>0) { x } else { sample(x,min(length(x),n)) } }

## On maps
require(rgdal)

euplot <- function (x,scale=15,lab,themap,cols=countrycols,mincex=.25,
        longs=indivinfo$long[match(names(x),indivinfo$COUNTRY_SELF)],
        lats=indivinfo$lat[match(names(x),indivinfo$COUNTRY_SELF)],
        legend=FALSE,legendloc=c(35,60),legendunits="",
        ...) {
    # map("world",xlim=c(-10,38), ylim=c(35,61),proj="globular",col=grey(.50),mar=c(1,1,1,1),resolution=0, ...)
    # xy <- mapproject( longs, lats )  # note: uses previous projection, passing this in messes it up.
    xylims <- project( cbind(long=c(-10,38),lat=c(35,61)), proj=proj4string(themap) )
    plot( themap, xlim=xylims[,1], ylim=xylims[,2], mar=c(1,1,1,1), border=grey(.5) )
    xy <- project( cbind(longs,lats), proj=proj4string(themap) )
    points( xy[,1], xy[,2], cex=pmax(mincex,sqrt(scale*abs(x))), pch=21, col=ifelse(x>0,"black","red"), bg=adjustcolor(cols[names(x)],.5), ... )
    if (!missing(lab)) { text( par("usr")[1:2]%*%c(.95,.05), par("usr")[4:3]%*%c(.95,.05), lab, pos=4, cex=1.2 ) }
    if (legend) {
        topleft <- project( cbind(legendloc[1],legendloc[2]), proj=proj4string(themap) )
        legsize <- 10^(floor(log10(12^2/scale)))
        points( topleft[,1], topleft[,2], cex=sqrt(scale*legsize), pch=21 )
        text( topleft[,1], topleft[,2], labels=paste(legsize,legendunits), pos=4, cex=.75, offset=0 )
    }
    return( invisible( xy ) )
}

## Plotting overlaps

plotseg <- function (olaps, chroms=1:22, chrstarts=.chrstarts, xlims, ylims=range(unlist(lapply(olaps,function(x)quantile(x$z,c(.001,.999),na.rm=TRUE)))), hilight=NULL, snps=FALSE, posvar="map", ylab="normalized # of segments", do.xlab=TRUE, cols=adjustcolor(rainbow(length(olaps)),.6), mar=c(2,4,1,2)+.1, legend=TRUE, chrspace=5, ... ) {
    # Pass mar=NULL if you don't want it to mess with the margins.
    chrst <- chrstarts
    chrst[ -chroms ] <- NA
    chrst[ c(chroms,length(chrst)) ] <- c(0,cumsum(diff(chrstarts)[chroms]+chrspace))
    chrends <- chrst-chrspace
    chrmids <- chrst[chroms] + diff(chrstarts)[chroms]/2
    if (missing(xlims)) xlims <- range(chrst,na.rm=TRUE)
    if (!is.null(mar)) { opar <- par(mar=mar) }
    plot( 0, type='n', xaxt='n', ylab=ylab, xlab='', xlim=xlims, ylim=ylims )
    # axis(1,labels=FALSE)
    # interesting regions and centromeres
    for ( k in seq_along(hilight) ) {
        x <- hilight[[k]]
        lims <- c(chrst[x[1]]+x[2], chrst[x[1]]+x[3])
        if ( x[1] %in% chroms & ( (lims[1] - xlims[2]) * (lims[2] - xlims[1]) < 0 ) ) {
            label <- names(hilight)[k]
            border <- adjustcolor("black",.5)
            if (diff(lims)>.5) { lwd <- 0; lty <- 1 } else { lwd <- 1; lty <- 5 }
            if (substring(names(hilight)[k],1,6)=="centro") { lwd <- 3; lty <- 1; label <- "c"; border <- grey(.75) }
            rect( xleft=lims[1], xright=lims[2], ybottom=ylims[1], ytop=ylims[2], col=grey(.90), border=border, lty=lty, lwd=lwd )
            if (nchar(label)>0) { text( lims[1], ylims[2], labels=label, adj=c(0,1), col='red' ) }
        }
    }
    for (k in 1:length(olaps)) {
        chroffsets <- chrst[olaps[[k]]$chrom] - chrstarts[olaps[[k]]$chrom]
        usethese <-  (olaps[[k]][,"chrom"] %in% chroms) & (olaps[[k]][,posvar]+chroffsets >= xlims[1]) & (olaps[[k]][,posvar]+chroffsets <= xlims[2])
        lines( olaps[[k]][usethese,posvar]+chroffsets[usethese], olaps[[k]][usethese,"z"], col=cols[k] )
        # SNP density
        if (k==1 & snps) {
            lines( olaps[[1]][usethese,posvar]+chroffsets[usethese], ylims[1] + (olaps[[1]][usethese,"nsnps"])/quantile(olaps[[1]][,"nsnps"],.95), col=adjustcolor("black",0.5)  )
            abline(h=ylims[1]+c(0,1),col=grey(.2))
        }
    }
    rect( chrends, par("usr")[3], chrst, par("usr")[4], col=grey(.75) )
    # abline(v=chrstarts,lwd=3,col=grey(.3))
    abline(h=0)
    axis(1, at=chrst, labels=FALSE)
    if (do.xlab) axis(1, at=chrmids, labels=paste("chromosome",seq_along(chrst))[chroms], tick=FALSE)
    if (legend) { legend("topright",lty=1,col=adjustcolor(cols,1),legend=names(olaps),...) }
    if (!is.null(mar)) { par(opar) }
    return( invisible( list( cols=cols, names=names(olaps) ) ) )
}


## Other

legend.svg <- function ( x, y, labels, ... ) {
    # ... can include things like col=, pch=, etc.
    require("RSVGTipsDevice")
    oop <- options("stringsAsFactors"=FALSE)  # otherwise this messes up passing in col=...
    pargs <- data.frame(x=x,y=y)
    pargs <- do.call( cbind, c( list(pargs), list(...) ) )
    if ("cex" %in% names(pargs)) {
        # reorder so biggest circles are underneath of the others
        reord <- order( pargs$cex, decreasing=TRUE )
        pargs <- pargs[reord,]
        labels <- labels[reord]
    }
    for (i in 1:length(labels)) {
        setSVGShapeToolTip(desc=labels[i])
        do.call(points, as.list(pargs[i,]))
    }
    options(oop)
}

colorize <- function (x, nc=32, colfn=function (n) rainbow_hcl(n,c=100,l=50), zero=FALSE, trim=0) {
    if (is.numeric(x) & trim>0) {
        x[ x<quantile(x,trim,na.rm=TRUE) ] <- quantile(x,trim,na.rm=TRUE)
        x[ x>quantile(x,1-trim,na.rm=TRUE) ] <- quantile(x,1-trim,na.rm=TRUE)
    }
    if (is.numeric(x)) {
        if (zero) {
            breaks <- seq( (-1)*max(abs(x),na.rm=TRUE), max(abs(x),na.rm=TRUE), length.out=nc )
        } else {
            breaks <- seq( min(x,na.rm=TRUE), max(x,na.rm=TRUE), length.out=nc )
        }
        x <- cut(x,breaks=breaks,include.lowest=TRUE)
    } else {
        x <- factor(x)
    }
    return( colfn(nlevels(x))[as.numeric(x)] )
}

pointify <- function (x, nc=10, trim=0) {
    as.numeric( cut( as.numeric(x), breaks=nc ) )
}

sample.df <- function (x,size,replace=FALSE,prob=NULL) {
    x[ sample.int(nrow(x),size,replace=replace,prob=prob),]
}

darken <- function (color,x) { y <- color; y[is.na(y)] <- "#FFFFFF"; ifelse( is.na(color), NA, hex(mixcolor(x,HSV(0,1,0),hex2RGB(substring(y,1,7)))) ) }
lighten <- function (color,x) { y <- color; y[is.na(y)] <- "#FFFFFF"; ifelse( is.na(color), NA, hex(mixcolor(x,HSV(0,1,1),hex2RGB(substring(y,1,7)))) ) }

rowplot <- function (rmatrix, sdmatrix, clist=colnames(rmatrix), rlist=clist, x, y, xlim, abbrevs=sapply(rlist,identity), cols=countrycols[rlist], oneline=(length(rlist)>10), log='', symmetrize=FALSE, ylab="long", ylabs=if (ylab=="long") paste(clist," (",abbrevs[clist],")",sep="") else abbrevs[clist],...) {
    # Plot the matrix rmatrix[rlist, clist], using labels,
    #   with each *column* of the matrix grouped together.
    # The vector clist gives the groupings on the y-axis,
    #   the vector rlist gives the entries displayed in each.
    # Standard deviations will be added if sdmatrix is present,
    #    and oneline controls whether to put dotted lines for each entry. 
    # Note do not need to pass in a matrix if x and y are given;
    #    passing in x and y as factors works, too.
    if (length(dim(rmatrix))==2) rmatrix <- rmatrix[rlist,clist]
    if (missing(x)) x <- col(rmatrix)
    if (missing(y)) y <- row(rmatrix)
    if (is.factor(x)) x <- match(x,clist)
    if (is.factor(y)) y <- match(y,rlist)
    if (symmetrize & is.null(dim(rmatrix))) {
        nondup <- x!=y
        xx <- c(x,y[nondup])
        yy <- c(y,x[nondup])
        rr <- c(rmatrix,rmatrix[nondup])
        x <- xx; y <- yy; rmatrix <- rr
    }
    if (log != '')  {
        rmatrix[rmatrix<=0] <- NA
    }
    ypos <- x+(y-min(y,na.rm=TRUE))/(1.2*diff(range(y,finite=TRUE)))-1/(2*1.2)
    if (ylab=="long") {
        opar <- par(mar=c( 4, 10, 1, 1 ) )
    } else if (ylab=="none") {
        opar <- par(mar=c( 4, 0, 1, 1 ) )
    } else {
        opar <- par(mar=c( 4, 5, 1, 1 ) )
    }
    if (missing(xlim)) xlim <- range(rmatrix,finite=TRUE)
    plot( 0, type='n', xlim=xlim, ylim=range(ypos,finite=TRUE), yaxt='n', ylab='', log=log, ... )
    abline(h=0.5+0:length(clist), lwd=2, col=grey(.5))
    if (ylab!="none") { axis( 2, at=1:length(clist), labels=ylabs, las=2 ) }
    text( rmatrix, ypos, labels=abbrevs[rlist[y]], col=darken(cols[y],0.9) )
    if (!oneline) {
        abline(h=ypos,lty=3,col=grey(.8))
    }
    if (!missing(sdmatrix)) {
        if (length(dim(sdmatrix))==2) { sdmatrix <- sdmatrix[rlist,clist][cbind(as.vector(x),as.vector(y))] }
        arrows( x0=rmatrix, x1=rmatrix+sdmatrix, y0=ypos, col=adjustcolor(cols[y],.2), lwd=2, angle=90, length=0 )
        arrows( x0=rmatrix, x1=rmatrix-sdmatrix, y0=ypos, col=adjustcolor(cols[y],.2), lwd=2, angle=90, length=0 )
    }
    par(opar)
}

plotses <- function (x, y, yse, ...) {
    # Draw standard errors
    #   avoiding "zero-length arrow" warning
    nonz <- !is.na(yse) & (yse>.001*diff(par("usr")[c(3,4)])/par("fin")[2])
    arrows(x0=x[nonz], y0=(y-2*yse)[nonz], y1=(y+2*yse)[nonz], angle=90, code=3, length=.1, ...)
}

hcolor <- function (z,alpha=.75,cols=adjustcolor(heat_hcl(64,h=c(40,360),l=70,c.=c(70,100)),alpha),nc=length(cols),...) {
    # coloring function for hplot
    cols[as.numeric(cut(pmin(1,pmax(-1,z)),breaks=seq(-1,1,length.out=nc+1),include.lowest=TRUE))]
}
hplot <- function (z,x=as.vector(col(z)),y=as.vector(row(z)),scale=1,max.cex=4,labs,xlabs=labs,ylabs=labs,alpha=.75,...) {
    # Like heatmap but with circles...
    plot( x=x, y=y, pch=20, cex=max.cex*sqrt(abs(z)/scale), col=hcolor(z/scale,alpha=alpha), xaxt='n', yaxt='n', xlab="", ylab="", ... )
    if (!missing(labs) | !missing(xlabs)) 
        axis(1, at=1:length(xlabs), labels=xlabs, las=2) 
    if (!missing(labs) | !missing(ylabs)) 
        axis(2, at=1:length(ylabs), labels=ylabs, las=2)
}

shrinkarrows <- function( x0, y0, x1=x0, y1=y0, amount=par("cxy")[1], ... ) {
    # Shrink each arrow by amount on both ends
    dx <- x1-x0
    dy <- y1-y0
    lens <- sqrt( dx^2 + dy^2 )
    dx <- dx / lens  # unit vectorize
    dy <- dy / lens
    x0 <- x0 + dx*amount
    y0 <- y0 + dy*amount
    x1 <- x1 - dx*amount
    y1 <- y1 - dy*amount
    arrows(x0,y0,x1,y1,...)
}

textlab <- function(x,...,nudge=.05,nudgex=nudge,nudgey=nudge) {
    # like text() but allow e.g. "topright" as in legend()
    if (is.character(x)) {
        auto <- match.arg(x, c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"))
        nx <- 0
        usr <- par("usr")
        if (par("xlog")) usr[1L] <- 10^usr[1L]
        if (par("ylog")) usr[2L] <- 10^usr[2L]
        x <- usr[1L:2L] %*% switch( auto, 
                bottomright=, topright=, right=c(nudgex,1-nudgex),
                bottomleft=, topleft=, left=c(1-nudgex,nudgex),
                top=, bottom=, center=c(.5,.5)
            )
        y <- usr[3L:4L] %*% switch( auto, 
                bottomright=, bottomleft=, bottom=c(1-nudgey,nudgey),
                topleft=, topright=, top=c(nudgey,1-nudgey),
                left=, right=, center=c(.5,.5)
            )
        pos <- switch( auto, 
                bottomright=, topright=, right=2,
                bottomleft=, topleft=, left=4,
                top=1, bottom=3, center=NULL
            )
        text( x, y, pos=pos, ... )
    } else { text(x,...) }
}

colaxis <- function ( side, at, labels=TRUE, col=par("fg"), line=par("mgp")[2], cex=1, ... ) {
    # make the axis labels in color
    labs <- axis( side, at=at, labels=FALSE, ... )
    if (length(labels==1) && labels==TRUE) { labels <- labs } 
    mtext( text=labels, side=side, line=line, at=at, col=col, cex=par("cex.axis")*par("cex")*cex, ... )
}

axislab <- function ( side, at, labels=TRUE, eps=par("cxy")[2], ... ) {
    # add *nonoverlapping* axis labels
    ord <- order(at)
    at <- at[ord]
    if (length(labels)>1) {
        labels <- labels[ord]
    }
    labat <- cumsum( c(0, pmax(eps,diff(at)) ) )
    labat <- min(at) + (max(at)-min(at))*labat/max(labat)
    axis( side, at, tick=TRUE, labels=FALSE, ... )
    axis( side, labat, tick=FALSE, labels=labels, ... )
}
