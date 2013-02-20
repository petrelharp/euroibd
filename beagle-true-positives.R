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

if (file.exists("true-pos-everything.Rdata")) {
    load("true-pos-everything.Rdata")
} else {
    # basename1 <- "HAPMAP_POPRES_NOERROR_23456_chr1"
    # basename22 <- "HAPMAP_POPRES_NOERROR_7729_chr22"
    basenames <- list(
            c( 1, "HAPMAP_POPRES_29562_chr1"),
            c( 2, "HAPMAP_POPRES_15766_chr2"),
            c( 8, "HAPMAP_POPRES_20862_chr8"),
            c( 9, "HAPMAP_POPRES_29549_chr9"),
            c( 10, "HAPMAP_POPRES_26294_chr10"),
            c( 11, "HAPMAP_POPRES_12869_chr11"),
            c( 12, "HAPMAP_POPRES_27981_chr12"),
            # c( 18, "HAPMAP_POPRES_7106_chr18"),  # none.
            c( 20, "HAPMAP_POPRES_9010_chr20"),
            c( 21, "HAPMAP_POPRES_1041_chr21"),
            c( 22, "HAPMAP_POPRES_8531_chr22")
        )
    tblocks <- do.call( rbind, lapply( basenames, function (x) {
                z <- read.table( paste(x[2], ".TRUE.ibd.gz",sep=""), header=TRUE )
                z$chrom <- as.numeric(x[1])
                return(z)
            } ) )

    # True blocks -- recorded as map position
    names(tblocks)[ names(tblocks) == "start" ] <- "mapstart"
    names(tblocks)[ names(tblocks) == "end" ] <- "mapend"
    tblocks$maplen <- tblocks$mapend - tblocks$mapstart
    hapinds <- levels(tblocks$id1)
    tblocks$id2 <- factor( tblocks$id2, levels=hapinds )

    # Get only blocks seen between those from hapmap
    # rawblocks <- readblocks("HAPMAP_POPRES_22411_chr1.combined.winnowed.fibd.gz",1)
    # rawblocks <- rbind( rawblocks, readblocks("HAPMAP_POPRES_26528_chr22.combined.winnowed.fibd.gz",22) )
    rawblocks <- do.call( rbind, lapply( basenames, function (x) {
                readblocks( paste(x[2], ".combined.nogaps.fibd.gz",sep=""), chr=as.numeric(x[1]), mapbase="hapmap.popres.genetic" )
            } ) )
    rawblocks <- subset.blocks(rawblocks, hapinds, only=TRUE)
    rawblocks$id1 <- factor( rawblocks$id1, levels=hapinds )
    rawblocks$id2 <- factor( rawblocks$id2, levels=hapinds )

    # Match up blocks
    matches <- match.blocks( tblocks, rawblocks )
    matches$tmaplen <- tblocks[ matches[,1], "maplen" ]
    matches$imaplen <- rawblocks[ matches[,2], "maplen" ]
    matches$score <- rawblocks[ matches[,2], "score" ]
    # does the 'true' block have multiple overlapping inferred blocks?
    matches$tmulti <- sapply(matches[,1], function (k) { sum( matches[,1] == k ) } )
    matches$imulti <- sapply(matches[,2], function (k) { sum( matches[,2] == k ) } )
    # Add to the collection of true blocks:
    #  total length of inferred overlapping segments
    tblocks$itotal <- sapply( 1:nrow(tblocks), function (k) { sum( matches$imaplen[matches[,1]==k] ) } )
    #  total overlapping length of inferred segments
    tblocks$ototal <- sapply( 1:nrow(tblocks), function (k) { sum( matches$omaplen[matches[,1]==k] ) } )
    #  number of overlapping inferred segments
    tblocks$ninferred <- sapply( 1:nrow(tblocks), function (k) { sum( matches[,1]==k ) } )
    # minimum score of overlapping inferred segments
    tblocks$score <- sapply(sapply(1:nrow(tblocks), function (k) matches$score[matches[,1]==k]),
            function(x) if(length(x)>0) { min(x) } else {NA} )

    # Find gaps in the inferred blocks
    allgaps <- get.gaps(rawblocks)
    missedgaps  <- ddply( matches[,c("b1","b2")], "b1", function (x) {
            if (nrow(x)<=1) {
                return( data.frame() )
            } else { return( get.gaps(rawblocks[x$b2,]) ) } 
        } )
    missedgaps$tmaplen <- tblocks[missedgaps$b1,"maplen"]
    tblocks$gaplen <- sapply( 1:nrow(tblocks), function (k) { sum( missedgaps$maplen[missedgaps$b1==k] ) } )

    # Above we looked at the raw output.  Here we evaluate our winnowing strategy.
    obsblocks <- do.call( rbind, lapply( basenames, function (x) {
                readblocks(paste(x[2],".combined.winnowed.fibd.gz",sep=""),as.numeric(x[1]),mapbase="hapmap.popres.genetic")
            } ) )
    obsblocks <- subset.blocks(obsblocks, hapinds, only=TRUE)
    obsblocks$id1 <- factor( obsblocks$id1, levels=hapinds )
    obsblocks$id2 <- factor( obsblocks$id2, levels=hapinds )

    obsmatches <- match.blocks( tblocks, obsblocks )
    obsmatches$tmaplen <- tblocks[ obsmatches[,1], "maplen" ]
    obsmatches$imaplen <- obsblocks[ obsmatches[,2], "maplen" ]
    obsmatches$score <- obsblocks[ obsmatches[,2], "score" ]
    # does the 'true' block have multiple overlapping inferred blocks?
    obsmatches$tmulti <- sapply(obsmatches[,1], function (k) { sum( obsmatches[,1] == k ) } )
    obsmatches$imulti <- sapply(obsmatches[,2], function (k) { sum( obsmatches[,2] == k ) } )
    # Add to the collection of true blocks:
    #  total length of inferred overlapping segments
    tblocks$obstotal <- sapply( 1:nrow(tblocks), function (k) { sum( obsmatches$imaplen[obsmatches[,1]==k] ) } )
    #  number of overlapping inferred segments
    tblocks$obsninferred <- sapply( 1:nrow(tblocks), function (k) { sum( obsmatches[,1]==k ) } )

    # Save!
    save(rawblocks, obsblocks, tblocks, hapinds, matches, allgaps, missedgaps, file="true-pos-everything.Rdata")
}

if (FALSE) {

png(file="alltrue-inferred-blocks.png", width=10*144, height=8*144, res=144)
# Plot all the true blocks
plotblocks(tblocks, chrom=c(2), yvals=as.numeric(tblocks$id2), col=rainbow(60)[as.numeric(tblocks$id1)])
plotblocks(tblocks, chrom=c(2), yvals=as.numeric(tblocks$id1), col=rainbow(60)[as.numeric(tblocks$id2)], add=TRUE)
# And the inferred blocks
plotblocks(rawblocks, chrom=c(2), yvals=as.numeric(rawblocks$id2), col=adjustcolor(rainbow(60)[as.numeric(rawblocks$id1)],.8), add=TRUE, yadj=.4)
plotblocks(rawblocks, chrom=c(2), yvals=as.numeric(rawblocks$id1), col=adjustcolor(rainbow(60)[as.numeric(rawblocks$id2)],.8), add=TRUE, yadj=.4)
dev.off()

# A random sample of pairs
pdf(file="false-pos-pairs.pdf", width=7, height=5)
spairs <- tblocks[ sample.int(nrow(tblocks),20),c("id1","id2")]
plotpairs( spairs, tblocks, chrom=c(1,2,22), col="black", yadj=-.2)
plotpairs( spairs, rawblocks, chrom=c(1,2,22), add=TRUE, yadj=+.2, col=c("red","green")[1+(rawblocks$score<1e-10)] )
legend( "topleft", col=c("black","green","red"), lty=1, legend=c("true", "inferred, score <1e-10", "inferred, score > 1e-10"))
dev.off()

# Look at true-versus-inferred lengths
with(matches, plot( tmaplen, imaplen ) )
abline(0,1)

# true versus inferred-overlapping
png(file="true-inferred-maplen.png", width=7*144, height=6*144, res=144)
nicols <- adjustcolor(c("black","red","green","purple","orange","brown"),.5)[1:length(unique(tblocks$ninferred))]
with(tblocks, plot( maplen, itotal+gaplen, col=nicols[ninferred+1], xlab="true length", ylab="total inferred length", pch=ifelse(is.na(score)|score>1e-10, 1, 20), main="True IBD blocks" ) )
with(tblocks[tblocks$gaplen>0,], segments( x0=maplen, y0=itotal+gaplen, x1=maplen, y1=itotal ) )
with(tblocks, lines(lowess( maplen[ninferred>0], (itotal+gaplen)[ninferred>0] ), lwd=2) )
legend("topleft", col=c(nicols,"black","black"), pch=c(rep(20,length(nicols)),1,NA), lty=c(rep(NA,length(nicols)+1),1), legend=c(paste(1:length(nicols)-1," overlapping blocks" ), "score>1e-10", "gap length") )
abline(0,1,lty=2)
dev.off()

# true versus overlapping inferred
with(tblocks, plot( maplen, ototal, col=nicols[ninferred+1] ) )
legend("topleft", col=nicols, pch=1, legend=paste(1:length(nicols)-1," overlapping blocks") )
abline(0,1)

# Length distribution of gaps
png(file="inferred-gap-hist.png", width=4*144, height=3*144, res=144)
hist(allgaps$maplen, xlim=c(0,5), breaks=(0:3000)/10, ylim=c(0,40), col="grey", main="", xlab="maplength")
par(new=TRUE)
hist(missedgaps$maplen, breaks=(0:40)/10, col="red", ylim=c(0,40), xlim=c(0,5), main="", xlab="")
dev.off()

# distribution of ratio of gaplen to segement length?
with(allgaps, hist(maplen/(leftmaplen+rightmaplen), breaks=(0:600)/4, col="grey", xlim=c(0,20), ylim=c(0,500), main="", xlab="gap/maplength") )
par(new=TRUE)
with( missedgaps, hist(maplen/(leftmaplen+rightmaplen), breaks=(0:600)/4, col="red", xlim=c(0,20), ylim=c(0,500), main="", xlab="") )

# something like this ratio against true block length
with(missedgaps, plot(tmaplen, maplen ) )
with(missedgaps, plot(tmaplen, maplen/(leftmaplen+rightmaplen) ) )

with( missedgaps, plot( maplen, leftmaplen+rightmaplen ) )
with( missedgaps, plot( maplen, tmaplen ) )

# Importance of snp density?
with( missedgaps, plot( maplen, nsnps, type="n") )
# what does snp density in segments look like?
with( rawblocks, points( maplen, nsnps, pch=2, cex=0.5, col="orange" ) )
with( allgaps, points( maplen, nsnps, pch=16, col="red") )
with( missedgaps, points( maplen, nsnps, pch=16 ) )
with( missedgaps[missedgaps$maplen<.4,], abline( lm( nsnps~maplen )$coef ) )
legend("topleft", col=c("black","red","black"), pch=c(16,16,2), legend=c("true gap","all gaps","inferred IBD"))



# look at a subset of gaps
gpairs <- missedgaps[ sample.int(nrow(missedgaps),20), c("id1","id2") ]
plotpairs( gpairs, tblocks, chrom=c(1,2,22), col="black", yadj=-.2)
plotpairs( gpairs, rawblocks, chrom=c(1,2,22), add=TRUE, yadj=+.2, col=c("red","green")[1+(rawblocks$score<1e-10)] )
legend( "topleft", col=c("black","green","red"), lty=1, legend=c("true", "inferred, score <1e-10", "inferred, score > 1e-10"))

# compare to all gaps
png(file="gaplen-distrn.png", width=5*144, height=4*144, res=144)
with( allgaps, plot( maplen, leftmaplen+rightmaplen, xlim=c(0,5), ylab="adjacent seg lengths", xlab="gap length" ) )
with( missedgaps, points( maplen, leftmaplen+rightmaplen, col="red", pch=20 ) )
legend("topright", col=c("black","red"), pch=c(1,20), legend=c("gap, not ibd", "gap, should be ibd") )
dev.off()

# How does gappiness relate to length of segment?
png(file="ninferred-maplen-gaps.png", width=5*144, height=4*144, res=144, pointsize=10)
with(tblocks, plot(maplen, jitter(ninferred), col=rainbow(20)[ifelse(is.na(score),1,pmin(1-floor(log10(score)),20))]), xlab="length of ibd segment", ylab="number of distinct inferred blocks")
legend("topleft",pch=1,col=rainbow(20),legend=1:20, cex=0.5)
dev.off()

with(tblocks, plot(gaplen, jitter(ninferred), col=rainbow(20)[ifelse(is.na(score),1,pmin(1-floor(log10(score)),20))]))
legend("topleft",pch=1,col=rainbow(20),legend=1:20)

with(tblocks, plot(maplen, tmaplen, col=rainbow(20)[ifelse(is.na(score),1,pmin(1-floor(log10(score)),20))]) )
with(tblocks, plot(maplen, tmaplen/maplen, col=rainbow(20)[ifelse(is.na(score),1,pmin(1-floor(log10(score)),20))]) )


#######
# What is our power?

# estimate power in bins
tb.quants <- quantile( tblocks$maplen, probs=seq(0,1,length.out=40) )
tb.mids <- tb.quants[-1] - diff(tb.quants)/2
tb.counts <- with(tblocks, table( cut(maplen, breaks=tb.quants), (ninferred>0 & score<1e-9) ) )
tb.power <- tb.counts[,"TRUE"]/apply(tb.counts,1,sum)
powerlenscoreq <- glm( (ninferred>0 & score<1e-9) ~ maplen + I(sqrt(maplen)), data=tblocks[tblocks$maplen>0.5,], family=binomial)

plcols <- c("green","black")[ 2-is.na(tblocks$score) ]
plcols[ (!is.na(tblocks$score) & tblocks$score>1e-9) ] <- "red"
mlvals <- (1:100)/10; 

png(file="powerplot.png", width=8*144, height=5*144, res=144)
with( tblocks, plot( maplen, jitter((ninferred>0 & score<1e-9)*1.0, amount=.05), col=plcols, ylab="prob found by beagle", xlim=c(0,10) ) )
lines( tb.mids, tb.power, lwd=2 )
lines(mlvals, predict(powerlenscoreq,newdata=data.frame(maplen=mlvals), type="response"), lty=2, col="red", lwd=2 )
abline(v=0.5, lty=2)
legend("bottomright", pch=c(1,1,1,NA,NA), lty=c(NA,NA,NA,1,2), title="true blocks:", col=c("black","red","green"), legend=c("score <1e-9", "score >1e-9", "not found", "percent found", "glm fit") )
dev.off()

# for talk
pdf(file="powerplot.pdf", width=3.5, height=3, pointsize=10)
par(mai=par("mai")*c(.7,.8,.1,.1))
with( tblocks, plot( maplen, jitter((ninferred>0 & score<1e-9)*1.0, amount=.05), pch=20, col=adjustcolor("black",.2), xlab="map length (cM)", ylab="prob found by beagle", xlim=c(0,10), cex.lab=1.2, mgp=c(2.3,1,0) ) )
lines( tb.mids, tb.power, lwd=2 )
lines(mlvals, predict(powerlenscoreq,newdata=data.frame(maplen=mlvals), type="response"), lty=2, col="red", lwd=2 )
abline(v=0.5, lty=2)
legend("bottomright", lty=1, legend=c("percent found", "glm fit"), col=c("black","red"), lwd=3 )
dev.off()

# Look at logit scores for fitting:
#  (i.e. why we put in a sqrt(maplen) term
logit <- function (x) { ifelse( x<0 | x>1, 0, log( x/(1-x) ) ) }
plot( tb.mids, logit( tb.power ) )
powerlenscore <- glm( (ninferred>0 & score<1e-9) ~ maplen, data=tblocks[tblocks$maplen>0.5,], family=binomial)
lines(mlvals, logit( predict(powerlenscore,newdata=data.frame(maplen=mlvals), type="response") ), lty=2, col="red" )
powerlenscoreq <- glm( (ninferred>0 & score<1e-9) ~ maplen + I(sqrt(maplen)), data=tblocks[tblocks$maplen>0.5,], family=binomial)
lines(mlvals, logit( predict(powerlenscoreq,newdata=data.frame(maplen=mlvals), type="response") ), lty=2 )

pdf("power-and-trueinferred-maplen.pdf", height=6, width=7.5, title="Power simulations", pointsize=14)
pmai <- par("mai")
par(mai=pmai*c(0,1,1/2,1))
layout((1:2),heights=c(1,2))

# What is our power?
powerlenscore <- glm( (ninferred>0 & score<1e-9) ~ maplen, data=tblocks, family=binomial)
with( tblocks, plot( maplen, jitter((ninferred>0 & score<1e-9)*1.0, amount=.05), xaxt="n", xlab="", ylab="Power", xlim=c(0,20), col="red", pch=20, cex=.25, yaxt="n" ) )
axis(side=2,at=c(0,0.5,1))
mlvals <- (1:100)/10; 
lines(mlvals, predict(powerlenscore,newdata=data.frame(maplen=mlvals), type="response"), lwd=2 )
text(20,0,labels="A",cex=3,pos=3)

# true versus inferred-overlapping
par(mai=pmai*c(1,1,0,1))
nicols <- adjustcolor(c("black","red","green","purple","orange","brown"),.5)[1:length(unique(tblocks$ninferred))]
with(tblocks, plot( maplen, (itotal+gaplen)*(!is.na(score)&score<1e-9), col=nicols[(!is.na(score)&score<1e-9)*ninferred+1], xlim=c(0,20), xlab="true length (cM)", ylab="total inferred length (cM)", pch=20 ) )
# with(tblocks[tblocks$gaplen>0,], segments( x0=maplen, y0=itotal+gaplen, x1=maplen, y1=itotal ) )
with(tblocks, lines(lowess( maplen[ninferred>0], (itotal+gaplen)[ninferred>0] ), lwd=2) )
legend("topleft", col=c(nicols), pch=rep(20,length(nicols)), legend=paste(1:length(nicols)-1," overlapping blocks" ) )
# legend("topleft", col=c(nicols,"black"), pch=c(1,rep(20,3)), legend=c("0 overlapping blocks","1 overlapping blocks","2 overlapping blocks",">2 overlapping blocks") )
abline(0,1,lty=2)
text(20,0,labels="B",cex=3,pos=3)


dev.off()

}
