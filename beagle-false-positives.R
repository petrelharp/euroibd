# Written 2013 by Peter Ralph and Graham Coop
# 
# contact: petrel.harp@gmail.com
#
#     To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty. 
# 
# You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 
# 
# Usage:
# Find false positive rate from beagle runs on randomly mixed up chromosomes.

source("ibd-blocks-fns.R")

# Load in false positives
if (!file.exists("false-pos-everything.Rdata")) {
    chromlist <- c(1,13,22)
    fpblocks <- do.call( rbind, lapply( chromlist, function (chrom) {
            fp <- readblocks(paste("RANDOM.REGIONAL.POPRES_chr",chrom,".combined.fibd.gz",sep=""), chr=chrom, header=TRUE, endblocks=FALSE)
            fp$chrom <- chrom
            return(fp)
        } ) )
    save(fpblocks, chromlist, file="false-pos-everything.Rdata")
} else {
    load("false-pos-everything.Rdata")
}

# Observed blocks
blocks <- do.call( rbind, lapply( chromlist, function (chrom) getblocks(chrom,"combined") ) )


# Restrict to those seen in at least two runs:
fpgb <- (fpblocks$nsegs > 1) & (fpblocks$score <= 1e-9) & (!is.na(fpblocks$maplen))
gb <- (blocks$nsegs > 1) & (blocks$score <= 1e-9)

require(hexbin)
with(fpblocks, plot( hexbin( log10(score) ~ log10(maplen) ), style="centroids" ) )

# How many FP at different score cutoffs and for different numbers of runs
png(file="chr22-nfp.png",width=5*72,height=5*72,res=72)
scuts <- exp( seq(log(1e-10),log(1e-15),len=100) )
nfp <- sapply(scuts, function (cutoff) sum( fpblocks$score<=cutoff, na.rm=TRUE ) )
npairs <- choose( length( unique( c(fpblocks$id1, fpblocks$id2) ) ), 2 )
nfp <- nfp/npairs
plot(scuts, nfp, log="x", type="l", xlab="score cutoff", ylab="number of false positives per pair", ylim=c(0,max(nfp)))
# by number of runs cutoff
for (k in 1:5) {
    nfp <- sapply(scuts, function (cutoff) sum( (fpblocks$nsegs>k & fpblocks$score<=cutoff)/npairs, na.rm=TRUE ) )
    lines(scuts, nfp, col=rainbow(5)[k], lty=2)
}
# How many FP at different score cutoffs and for different lengths
lcuts <- (1:5)-0.5
# by number of runs cutoff
for (k in 1:5) {
    nfp <- sapply(scuts, function (cutoff) sum( (fpblocks$maplen>lcuts[k] & fpblocks$score<=cutoff & fpblocks$nsegs>1)/npairs, na.rm=TRUE ) )
    lines(scuts, nfp, col=rainbow(5)[k], lty=1)
}
legend("topleft", col=rainbow(5), lty=c(rep(2,5),0,rep(1,5)), legend=c(paste("seen in at least", 1:5, "runs"),"seen in at least 2 runs and...", paste("at least", lcuts, "cM long")), cex=0.5) 
dev.off()



# Tabulate this.
fptab <- with(fpblocks, table( pmin(15,floor(-log10(score))), nsegs, floor(maplen) ) )
names(dimnames(fptab)) <- c("log10score","nsegs","maplenCM")
ftable(addmargins(fptab))

# Small table
smfptab <- with(fpblocks, table( score <= 1e-10, nsegs>1, pmin(4,floor(maplen)) ) )
dimnames(smfptab) <- list(score=c(">1e-10","<=1e-10"),nsegs=c("1",">1"),maplenMb=c(paste(c(0:3),"Mb"),">4Mb") )
ftable(addmargins(smfptab))

# Plot these:
png(file="chr22-pos-vs-score-with-fp.png", width=5*72, height=5*72, res=72)
with(blocks[gb,], plot( maplen, score, log="xy", pch=20, col=adjustcolor("black",0.2), ylim=c(1e-23,1e-10), xlim=c(.5,10) ) )
with(fpblocks[fpgb,], points( maplen, score, pch=20, col=adjustcolor("red",0.2) ) )
legend("bottomleft",col=c("black","red"),legend=c("observed","falsepos"),pch=20)
dev.off()

# Plot proportions
png(file="chr22-observed-rates.png", width=8*72, height=5*72, res=72)
with(blocks[gb,], plot( ecdf(log10(maplen)) ) )
with(fpblocks[fpgb,], lines( ecdf(log10(maplen)), col="red" ) )
lines( seq(-2,1,len=50), sapply( 10^seq(-2,1,len=50), function (x) { sum( fpgb & fpblocks$maplen>x, na.rm=TRUE ) / sum( gb & blocks$maplen>x, na.rm=TRUE ) } ), col="green" )
legend("topleft", col=c("black","red","green"), legend=c("observed ecdf","falsepos ecdf","(num fp >x)/(num obs >x)"), lty=1)
dev.off()

# Examine position along the chromosome
png(file="chr22-fp-locations.png", width=8*72, height=6*72, res=72)
with( blocks[gb,], plot( mapmid, maplen, pch=20, col=adjustcolor("black",.3)) )
with( fpblocks[fpgb,], points( mapmid, maplen, pch=20, col=adjustcolor("red",.3) ) )
rug( chrmap$map, col=adjustcolor("black"))
legend("topleft",col=c("black","red"), pch=20, legend=c("observed","false pos"))
dev.off()

# Plot some examples
png(file="300-examples.png", width=7*72, height=5*72, res=72)
plotmany( unique( sample(fpblocks$id1,300)), fpblocks[fpgb,], chrom=22, lwd=3 )
dev.off()

#########
# FP rate by country pair
# False positives
if (!file.exists("false-pos-winnowed.Rdata")) {
    chromlist <- c(1,13,22)
    fpblocks <- do.call( rbind, lapply( chromlist, function (chrom) {
            fp <- readblocks(paste("RANDOM.REGIONAL.POPRES_chr",chrom,".combined.winnowed.fibd.gz",sep=""), chr=chrom, header=TRUE, endblocks=FALSE)
            fp$chrom <- chrom
            return(fp)
        } ) )
    indivinfo <- getsampleinfo(remove.qc=TRUE)
    fpblocks <- subset(fpblocks, (!is.na(maplen)) & indivinfo$YESOK[match(id1,indivinfo$SUBJID)] & indivinfo$YESOK[match(id2,indivinfo$SUBJID)] )  # remove qc'ed indivs
    fp.ids <- unique( c(fpblocks$id1, fpblocks$id2) )
    indivinfo <- droplevels( subset( indivinfo, SUBJID %in% fp.ids ) )
    # We rearranged genomes before splitting out England, Wales, etc from UK
    # check this by untarring e.g. reorder.REGIONAL.POPRES_chr13.beagle.tar.gz to test and running:
    #    tmpindivinfo <- getsampleinfo(remove.qc=FALSE)
    #    fnames <- list.files("test","^country.*",full.names=TRUE); z <- do.call( rbind, lapply( fnames, function (x) read.table(x,header=FALSE,stringsAsFactors=FALSE) ) );
    #    names(z) <- c("SUBJID","COUNTRY_SELF"); z$other <- as.character( indivinfo$COUNTRY_SELF[ match(z$SUBJID,indivinfo$SUBJID) ] );
    #    z <- subset(z,tmpindivinfo$YESOK[match(SUBJID,tmpindivinfo$SUBJID)]); 
    #    subset(z,COUNTRY_SELF!=other); all( sort(unique(z$SUBJID)) == sort(fp.ids) )  # TRUE
    for (x in list( 
            list( "United Kingdom", c("England","Scotland","Wales") ),
            list( "Switzerland", c("Swiss French","Swiss German") ),
            list( "Yugoslavia", c("Albania","Bosnia","Croatia","Serbia","Macedonia","Kosovo") )
            ) ) {
        levels(indivinfo$COUNTRY_SELF)[levels(indivinfo$COUNTRY_SELF) %in% x[[2]]] <- x[[1]]
    }
    indivinfo$COUNTRY_SELF <- ordered(indivinfo$COUNTRY_SELF,levels=sort(levels(indivinfo$COUNTRY_SELF)))
    # subset( as.data.frame(with(subset(indivinfo,SUBJID %in% fp.ids), table(COUNTRY_SELF,ORIG_COUNTRYSELF))), Freq>0 )
    fpblocks$country1 <- ordered( indivinfo$COUNTRY_SELF[ match(fpblocks$id1,indivinfo$SUBJID) ], levels=(levels(indivinfo$COUNTRY_SELF)) )
    fpblocks$country2 <- ordered( indivinfo$COUNTRY_SELF[ match(fpblocks$id2,indivinfo$SUBJID) ], levels=(levels(indivinfo$COUNTRY_SELF)) )
    fpblocks$countrypair <- gsub(" ",".",with(fpblocks, paste( levels(fpblocks$country1)[ifelse(country1<country2,country1,country2)], levels(fpblocks$country1)[ifelse(country1<country2,country2,country1)], sep="-" ) ) ) 
    # Our observed data
    obsblocks <- subset( do.call( rbind, lapply( chromlist, function (chrom) getblocks(chrom,"winnowed") ) ), (id1 %in% fp.ids) & (id2 %in% fp.ids) )
    obsblocks$country1 <- ordered( indivinfo$COUNTRY_SELF[ match(obsblocks$id1,indivinfo$SUBJID) ], levels=(levels(indivinfo$COUNTRY_SELF)) )
    obsblocks$country2 <- ordered( indivinfo$COUNTRY_SELF[ match(obsblocks$id2,indivinfo$SUBJID) ], levels=(levels(indivinfo$COUNTRY_SELF)) )
    obsblocks$countrypair <- gsub(" ",".",with(obsblocks, paste( levels(obsblocks$country1)[ifelse(country1<country2,country1,country2)], levels(obsblocks$country1)[ifelse(country1<country2,country2,country1)], sep="-" ) ) ) 
    save(chromlist, indivinfo, fpblocks, obsblocks, file="false-pos-winnowed.Rdata")
} else {
    load("false-pos-winnowed.Rdata")
}

# number of pairs
npairs <- countpairs(sampleinfo=indivinfo)
totalpairs <- sum(npairs)

# FP rate by length
lenbins <- quantile( obsblocks$maplen[obsblocks$maplen>0.0], seq(0,1,length.out=10), na.rm=TRUE )
midbins <- lenbins[-1]-diff(lenbins)/2; midbins[length(midbins)] <- midbins[length(midbins)-1]+3
fp.counts <- with(fpblocks, table( cut( maplen[maplen>0.0], breaks=lenbins ) ) )
true.counts <- with(obsblocks, table( cut( maplen[maplen>0.0], breaks=lenbins ) ) )
## Now look by country
fp.c.counts <- with(fpblocks, table( cut( maplen[maplen>0.0], breaks=lenbins ), countrypair[maplen>0.0] ) )
fp.c.counts <- fp.c.counts[,match(names(npairs),dimnames(fp.c.counts)[[2]])]  # make order agree with npairs
# fp.c.counts[is.na(fp.c.counts)] <- 0
true.c.counts <- with(obsblocks, table( cut( maplen[maplen>0.0], breaks=lenbins ), countrypair[maplen>0.0] ) )
true.c.counts <- true.c.counts[,match(names(npairs),dimnames(true.c.counts)[[2]])]  # make order agree with npairs
# true.c.counts[is.na(true.c.counts)] <- 0
### only those we have at least 300 fpblocks for
manyblocks <- names(npairs)[ apply(fp.c.counts,2,sum,na.rm=TRUE)>300 ]
paircols <- rainbow(length(manyblocks))
names(paircols) <- manyblocks

## Absolute rate
# png(file="false-pos-by-country.png", width=7*144, height=5*144, res=144)
pdf(file="false-pos-by-country.pdf", width=8, height=5, pointsize=10)
plot(midbins, as.numeric(fp.counts/totalpairs)/diff(lenbins), type='l', xlim=c(0,10), lwd=2, ylim=c(0,.035), main="false positive rate")
abline(h=c(0,1))
invisible( sapply( manyblocks, function (k) lines( midbins, (fp.c.counts[,k]/npairs[k])/diff(lenbins), col=paircols[k] ) ) )
lines(midbins, exp(-7.108-1.879*midbins+3.705*sqrt(midbins)), lwd=2, lty=3)
legend("topright",legend=manyblocks, lty=1, col=paircols)
text( midbins[1], (fp.c.counts[1,manyblocks]/npairs[manyblocks])/diff(lenbins)[1], labels=manyblocks)
dev.off()

# Huh, look at PT-PT blocks compared to others
# pdf(file="example-false-pos.pdf", width=10,height=8,pointsize=10)
png(file="example-false-pos.png", width=10*144,height=8*144,res=144)
layout(1:3)
par(mar=c(2,4,2,2)+.1)
plotindivs(subset(fpblocks,(country1==country2)&(country2=="Portugal")),chroms=chromlist,main=paste("false positives, Portugal, n=",nsamples["Portugal"]))
theseones <- sample(subset(indivinfo,COUNTRY_SELF=="Switzerland")$SUBJID,nsamples["Portugal"])
plotindivs(subset(fpblocks,(id1%in%theseones)&(id2%in%theseones)),chroms=chromlist,main=paste("false positives, Switzerland, n=",nsamples["Portugal"]))
theseones <- sample(subset(indivinfo,COUNTRY_SELF=="Italy")$SUBJID,nsamples["Portugal"])
plotindivs(subset(fpblocks,(id1%in%theseones)&(id2%in%theseones)),chroms=chromlist,main=paste("false positives, Italy, n=",nsamples["Portugal"]))
dev.off()
# Yup; no doubt.

## Relative to observed rate 
png(file="false-pos-relative.png", width=7*144, height=5*144, res=144)
plot(midbins, as.numeric(fp.counts/totalpairs)/as.numeric(true.counts/totalpairs), type='l', xlim=c(0,10), ylim=c(0,1))
abline(h=c(0,1), col='red')
abline(h=.1, lty=3, col='red')
dev.off()

# by country
pdf(file="false-pos-relative-by-country.pdf", width=7, height=5, pointsize=10)
par(mar=c(4,10,2,4)+.1)
matplot(midbins,(fp.c.counts/true.c.counts)[,manyblocks],type='l',col=rainbow(ncol(fp.c.counts)),yaxt='n', xlab="block length (cM)", ylab='', ylim=c(0,1), main="false positive rate relative to observed" )
axis(side=4); 
axislab(side=2,at=(fp.c.counts/true.c.counts)[1,manyblocks], labels=manyblocks, las=2, cex.axis=8/12)
dev.off()

## Compared to observed rate
# png(file="false-pos-and-true-rate.png", width=7*144, height=5*144, res=144)
pdf(file="false-pos-and-true-rate.pdf", width=7, height=5, pointsize=10)
plot(midbins, sum(.chrlens)*as.numeric(true.counts/totalpairs)/(diff(lenbins)*sum(.chrlens[chromlist])), type='l', xlim=c(0,6), xlab="map length (cM)", ylab="rate per pair", main="False positive rates")
invisible( sapply( manyblocks, function (k) lines( midbins, sum(.chrlens)*(fp.c.counts[,k]/npairs[k])/(diff(lenbins)*sum(.chrlens[chromlist])), col=paircols[k] ) ) )
label.these <- (fp.c.counts[1,]/npairs)>.01
text( midbins[1], sum(.chrlens)*(fp.c.counts[1,]/npairs)[label.these]/(diff(lenbins)[1]*sum(.chrlens[chromlist])), labels=colnames(fp.c.counts)[label.these], cex=8/12, col=adjustcolor("black",.5))
lines(midbins, sum(.chrlens)*as.numeric(fp.counts/totalpairs)/(diff(lenbins)*sum(.chrlens[chromlist])), lty=2)
lines(midbins, sum(.chrlens)*exp(-13.704-2.095*midbins+4.381*sqrt(midbins)), lwd=2, lty=3)
legend("topright", legend=c("Observed mean IBD rate","Estimated mean false positive rate","Fitted mean fp rate","Estimated mean fp rate, by country-pair"), lty=c(1,2,3,1),col=c("black","black","black","red"))
# abline(h=c(0,1), col='red')
# abline(h=.1, lty=3, col='red')
dev.off()

# Look at short ones
makecp <- function(x,y) gsub(" ",".",paste(sort(c(x,y)),collapse="-"))
clist <- c("Italy","Spain","Portugal","Swiss.German","Swiss.French","United.Kingom","Germany","France","Ireland")
tmp <- sapply( clist, function (x) {
        c( with( subset(fpblocks,countrypair==makecp(x,x)), sum( maplen>.4 & maplen<.5 ) ) / npairs[makecp(x,x)],
           with( subset(obsblocks,countrypair==makecp(x,x)), sum( maplen>.4 & maplen<.5 ) ) / npairs[makecp(x,x)],
           with( subset(fpblocks,countrypair==makecp(x,"Spain")), sum( maplen>.4 & maplen<.5 ) ) / npairs[makecp(x,"Spain")],
           with( subset(obsblocks,countrypair==makecp(x,"Spain")), sum( maplen>.4 & maplen<.5 ) ) / npairs[makecp(x,"Spain")],
           with( subset(fpblocks,countrypair==makecp(x,"Portugal")), sum( maplen>.4 & maplen<.5 ) ) / npairs[makecp(x,"Portugal")],
           with( subset(obsblocks,countrypair==makecp(x,"Portugal")), sum( maplen>.4 & maplen<.5 ) ) / npairs[makecp(x,"Portugal")]
           )
   } )
barplot( data.frame( self=tmp[2,]-tmp[1,], ES=tmp[4,]-tmp[3,], PT=tmp[6,]-tmp[5,] ), beside=TRUE )

# Fit something to the mean block density per (total) maplen
lenbins <- quantile( obsblocks$maplen[obsblocks$maplen>0.5], seq(0,1,length.out=100), na.rm=TRUE )
midbins <- lenbins[-1]-diff(lenbins)/2
fp.counts <- with(fpblocks, table( cut( maplen[maplen>0.5], breaks=lenbins ) ) )
true.counts <- with(obsblocks, table( cut( maplen[maplen>0.5], breaks=lenbins ) ) )
plot(midbins, as.numeric(fp.counts/totalpairs)/as.numeric(true.counts/totalpairs), type='l', xlim=c(0,10), log='y')

fplm <-  lm( log(as.numeric(fp.counts/totalpairs)/(diff(lenbins)*sum(.chrlens[chromlist]))) ~ midbins, subset=midbins>1 )
fplmsq <-  lm( log(as.numeric(fp.counts/totalpairs)/(diff(lenbins)*sum(.chrlens[chromlist]))) ~ midbins + I(sqrt(midbins)), subset=midbins>1 & midbins<6)
#    (Intercept)           midbins  I(sqrt(midbins))  
#      -7.108            -1.879             3.705  
fplmsqsq <-  lm( log(as.numeric(fp.counts/totalpairs)/(diff(lenbins)*sum(.chrlens[chromlist]))) ~ midbins + I(sqrt(midbins))+ I(midbins^(-1)), subset=midbins>1 )
#     (Intercept)           midbins  I(sqrt(midbins))      I(midbins^2)  
#       -14.1369           -2.5989            5.2891            0.0402  

par(mfrow=c(1,2))
plot(midbins, as.numeric(fp.counts/totalpairs)/(diff(lenbins)*sum(.chrlens[chromlist])), type='b', xlim=c(0,6), log='y')
lines( midbins, exp(predict( fplm, newdata=list(midbins=midbins) )), lty=3, col='red' )
lines( midbins, exp(predict( fplmsq, newdata=list(midbins=midbins) )), lty=2, col='green' )
lines( midbins, exp(predict( fplmsqsq, newdata=list(midbins=midbins) )), lty=4, col='purple' )
plot(midbins, as.numeric(fp.counts/totalpairs)/(diff(lenbins)*sum(.chrlens[chromlist])), type='b', xlim=c(0,300) )
lines( seq(0,300,length.out=100), exp(predict( fplm, newdata=list(midbins=seq(0,300,length.out=100)) )), lty=3, col='red' )
lines( seq(0,300,length.out=100), exp(predict( fplmsq, newdata=list(midbins=seq(0,300,length.out=100)) )), lty=2, col='green' )
lines( seq(0,300,length.out=100), exp(predict( fplmsqsq, newdata=list(midbins=seq(0,300,length.out=100)) )), lty=4, col='purple' )
