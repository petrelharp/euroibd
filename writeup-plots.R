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
source("/home/peter/projects/genome/laplace-inversion-fns.R")

##########
# Power and false positives
load("../true-pos/true-pos-everything.Rdata")
load("../false-pos/false-pos-winnowed.Rdata")

# Actual blocks and other metainformation
load("../ibdblocks/all-blocks-winnowed-fine.Rdata")
load("../ibdblocks/eda-data-fine.Rdata")

# make a variable which is country-pair
blocks$countrypair <- countrypairs[cbind(as.numeric(blocks$country1),as.numeric(blocks$country2))]
# and individual-pair
blocks$indivpair <- factor( paste( pmin(blocks$id1,blocks$id2), pmax(blocks$id1,blocks$id2), sep="-" ) )

#####
# Categories

easterns <- c( "Albania", "Kosovo", "Croatia", "Bosnia", "Montenegro", "Serbia", "Yugoslavia", "Slovenia", "Macedonia", "Greece", "Bulgaria", "Romania", "Poland", "Austria", "Hungary", "Slovakia", "Czech Republic", "Russia", "Ukraine" )
northerns <- c("Latvia", "Finland", "Sweden", "Norway", "Denmark")
mideasterns <- c("Cyprus","Turkey")
southerns <- c("Italy","Portugal","Spain")
westerns <- c("France","United Kingdom","England", "Scotland", "Ireland","Belgium","Netherlands", "Switzerland", "Swiss French","Swiss German", "Germany")

newcats <- list( 
        I=southerns,
        W=westerns,
        N=northerns,
        E=easterns,
        TC=mideasterns
    )

#######
# False pos, power, and various observed distributions

# divide pairs by distance
citydists <- getcitydists()
npairs <- countpairs(blocks=blocks)
npairs <- cbind( data.frame( do.call(rbind,strsplit(names(npairs),"-")) ), npairs )
names(npairs) <- c("country1","country2","npairs")
npairs$gdist <- citydists[ cbind( match(npairs$country1,gsub(" ",".",rownames(citydists))), match(npairs$country2,gsub(" ",".",colnames(citydists))) ) ]
allpairs <- sum(npairs$npairs)
neighborpairs <- sum( npairs$npairs[ npairs$gdist<1000 & npairs$country1 != npairs$country2 ] )
nonneighborpairs <- sum( npairs$npairs[ npairs$gdist>=1000 & npairs$country1 != npairs$country2 ] )

## False pos stuff
fphist <- hist( fpblocks$maplen, breaks=100, plot=FALSE )
fppairs <- sum( countpairs(blocks=fpblocks) )


## True pos stuff
tb.quants <- quantile( tblocks$maplen, probs=seq(0,1,length.out=40) )
tb.mids <- tb.quants[-1] - diff(tb.quants)/2
tb.counts <- with(tblocks, table( cut(maplen, breaks=tb.quants), (ninferred>0 & score<1e-9) ) )
tb.power <- tb.counts[,"TRUE"]/apply(tb.counts,1,sum)
powerlenscoreq <- glm( (ninferred>0 & score<1e-9) ~ maplen + I(sqrt(maplen)), data=tblocks[tblocks$maplen>0.5,], family=binomial)

plcols <- c("green","black")[ 2-is.na(tblocks$score) ]
plcols[ (!is.na(tblocks$score) & tblocks$score>1e-9) ] <- "red"
mlvals <- (1:100)/10; 


# Power, false positive rates, and observed IBD rates.
pdf(file="power-and-fp.pdf", width=7, height=4, pointsize=10)
layout( cbind( c(1,1), c(2,3) ), heights=c(1.8,1) )

nicols <- adjustcolor(c("black","red","green","purple","orange","brown"),.5)[1:length(unique(tblocks$ninferred))]
with(tblocks, plot( maplen, ifelse(is.na(score)|score>1e-9,0,itotal+gaplen), col=nicols[ifelse(is.na(score)|score>1e-9,0,ninferred)+1], xlab="true length", ylab="total inferred length", pch=20, main="(A)" ) )
with(tblocks[tblocks$gaplen>0,], segments( x0=maplen, y0=itotal+gaplen, x1=maplen, y1=itotal ) )
with(tblocks, lines(lowess( maplen[ninferred>0], (itotal+gaplen)[ninferred>0] ), lwd=2) )
legend("topleft", col=c(nicols,"black"), pch=c(rep(20,length(nicols)),NA), lty=c(rep(NA,length(nicols)),1), legend=c(paste(1:length(nicols)-1," overlapping blocks" ), "gap length") )
abline(0,1,lty=2)

opar <- par(mai=c(.05,par("mai")[2:4]))
with(withinhist, plot( mids, (counts/diff(breaks))/sum(choose(nsamples,2)), type='l', xlim=c(0,10), ylab="rate per pair", lwd=2, col='black', xaxt='n', xlab='', main="(B)") )
with(fphist, lines( mids, (sum(.chrlens)/sum(.chrlens[chromlist]))*(counts/diff(breaks))/fppairs, lty=2 ) )
with(neighborhist, lines( mids, (counts/diff(breaks))/neighborpairs, lwd=2, col="blue") )
with(nonneighborhist, lines( mids, (counts/diff(breaks))/nonneighborpairs, lwd=2, col="purple") )
legend("topright", legend=c("Within countries","Nearby countries","Distant countries","False positives"), lty=c(1,1,1,rep(3,length(clist)),2), col=c("black","blue","purple",'black') )

par(mai=c(opar$mai[1:2],.05,par("mai")[4]))
plot( tb.mids, tb.power, lwd=2, type='l', xlab="map length (cM)", ylab="power", xlim=c(0,10), yaxt='n' )
axis(2,at=c(0,0.5,1))
lines(mlvals, predict(powerlenscoreq,newdata=data.frame(maplen=mlvals), type="response"), lty=2, col="red", lwd=2 )

dev.off()



###
# Barplot of IBD rates

# ... or: plot the matrices --
## Indices in poppairs to pull out a (symmetric) matrix
clist <- c("United Kingdom", "Ireland", "Germany", "Swiss German", "Swiss French", "France", "Spain", "Portugal", "Italy", "Yugoslavia", "Hungary", "Poland")
rlist <- intersect( gsub("."," ", c(mideasterns,southerns,westerns,northerns,easterns), fixed=TRUE ), names(nsamples) )
ind.matrix <- apply( expand.grid(rlist,rlist), 1, sort )
ind.matrix <- gsub(" ",".", paste( ind.matrix[1,], ind.matrix[2,], sep="-" ) )
ind.matrix <- match(ind.matrix, poppairs$countrypair)
dim(ind.matrix) <- c(length(rlist),length(rlist))
dimnames(ind.matrix) <- list( rlist, rlist )
rate.matrix.1cM <- matrix( poppairs$nblocks1cM[ind.matrix]/poppairs$npairs[ind.matrix], nrow=length(rlist), dimnames=dimnames(ind.matrix))
rate.matrix.5cM <- matrix( poppairs$nblocks5cM[ind.matrix]/poppairs$npairs[ind.matrix], nrow=length(rlist), dimnames=dimnames(ind.matrix))
rate.matrix.10cM <- matrix( poppairs$nblocks10cM[ind.matrix]/poppairs$npairs[ind.matrix], nrow=length(rlist), dimnames=dimnames(ind.matrix))
npair.matrix <- matrix( poppairs$npairs[ind.matrix], nrow=length(rlist), dimnames=dimnames(ind.matrix))

rlist <- intersect( gsub("."," ", c(mideasterns,southerns,westerns,northerns,rev(easterns)), fixed=TRUE ), names(nsamples) )
pdf(file="sharing-rates-dotchart-long.pdf", width=8, height=20, pointsize=10)
    rmatrix <- (rate.matrix.1cM + rate.matrix.5cM + rate.matrix.10cM)
    vmatrix <- apply( nblock.vars[rlist,rlist,], c(1,2), sum, na.rm=TRUE )
    rowplot( rmatrix, sqrt(vmatrix), clist=rlist, rlist=rlist, log='x', xlab="rate per pair", oneline=TRUE, abbrevs=countryabbrevs, ylab=TRUE )
    abline(v=c(.25,.5,1,2),col=adjustcolor("black",.25))
dev.off()

require(maps)
mplot <- function (x,scale=15,label=x,cols=adjustcolor(countrycols[countries],.5),...) {
    map("world",xlim=c(-10,38), ylim=c(35,61),proj="globular",col=grey(.50),mar=c(1,1,1,1),resolution=0)
    countries <- levels(poppairs$country1)
    ii <- match(x,countries)
    longs <- indivinfo$long[match(countries,indivinfo$COUNTRY_SELF)]
    lats <- indivinfo$lat[match(countries,indivinfo$COUNTRY_SELF)]
    xy <- mapproject( longs, lats, proj="globular")
    z <- with( subset(poppairs,country1%in%x|country2%in%x), tapply((nblocks1cM+nblocks5cM+nblocks10cM)/npairs, ifelse(country1%in%x,country2,country1), sum) )
    # points( xy$x, xy$y, cex=sqrt(scale*rowSums(rmatrix[,x,drop=FALSE])), pch=21, bg=adjustcolor(countrycols[rownames(rmatrix)],.5), ... )
    points( xy$x, xy$y, cex=sqrt(scale*z), pch=21, bg=cols, ... )
    text( xy$x[ii], xy$y[ii], "*", cex=6 )
    text( par("usr")[1:2]%*%c(.95,.05), par("usr")[4:3]%*%c(.95,.05), label, pos=4, cex=1.2, font=2 )
}


####
# Plot example blocks and block density along the genome
load("overlaps.RData")

# for paper:
nrandoms <- 24
subset.countries <- list(c("Albania","Kosovo"),"Portugal","Spain","Swiss German","Italy","United Kingdom","Ireland")
random.ids <- unlist( lapply(subset.countries, function (country) sample( indivinfo$SUBJID[indivinfo$COUNTRY_SELF %in% country], min(nrandoms,sum(nsamples[country])) ) ) )
random.subset <- subset.blocks( blocks, idnums=random.ids, reorder.ids=TRUE, only=TRUE )

chroms <- c(1,8)
pdf( file="random-indivs-and-overlap-chr1-8.pdf", width=7, height=4, pointsize=10 )
layout( 1:2 )
opar <- par(mar=c(0,4,1,1)+.1)
plotblocks( random.subset, chroms=chroms, col=adjustcolor(countrycols,.75)[random.subset$country2], yvals=(match( random.subset$id1, random.ids )), yaxt='n', lwd=6, xlab="", xaxis=FALSE )
plotblocks( random.subset, chroms=chroms, col=adjustcolor(countrycols,.75)[random.subset$country1], yvals=(match( random.subset$id2, random.ids )), yaxt='n', lwd=6, add=TRUE )
abline(h=nrandoms*(0:length(subset.countries))+.5, lty=2)
axis(side=2, at=nrandoms*((1:length(subset.countries))-.5), labels=countryabbrevs[sapply(subset.countries,function(x)x[1])], las=2)
segments( x0=rep((par("usr")[1])/2,length(subset.countries)), y0=nrandoms*((1:length(subset.countries))-.75), y1=nrandoms*((1:length(subset.countries))-.25), col=countrycols[sapply(subset.countries,function(x)x[1])], lwd=10, lend=1 )
# 
par(mar=c(3,4,0,1)+.1)
plotseg( olaps, chrom=chroms, hilight=hilight, snps=TRUE, mar=NULL, bg=adjustcolor("white",.9), cex=6/10, ylab="norm. # blocks" )
dev.off()

pdf(file="overlap-rates-all-chr-onepage.pdf", width=6.5, height=8.5, pointsize=10)
labpos <- rep("topleft",22)
labpos[c(8,9,15)] <- "topright"
par(mar=c(1,3,0,1),mfcol=c(11,2))
for (chrom in 1:11) {
    figreg <- c(0, (3/4)*.chrlens[chrom]/.chrlens[1], (11-chrom)/11, (11-chrom+1)/11  )
    par(fig=figreg,new=TRUE)
    tmp <- plotseg( olaps, chrom=chrom, hilight=hilight[setdiff(names(hilight),c("17q21.31","9p","15p"))], snps=TRUE, mar=NULL, bg=adjustcolor("white",.9), cex=6/10, ylab="", do.xlab=FALSE, legend=FALSE )
    textlab(labpos[chrom], paste("chr", chrom))
    figreg[1] <- 1-(3/4)*.chrlens[22-chrom+1]/.chrlens[1]
    figreg[2] <- 1
    par(fig=figreg,new=TRUE)
    plotseg( olaps, chrom=22-chrom+1, hilight=hilight[setdiff(names(hilight),c("17q21.31","9p","15p"))], snps=TRUE, mar=NULL, bg=adjustcolor("white",.9), cex=6/10, ylab="", do.xlab=FALSE, legend=FALSE )
    textlab(labpos[22-chrom+1], paste("chr", 22-chrom+1))
}
par( fig=c( .47, .67, 2/11, 4/11 ), new=TRUE, ann=FALSE )
plot( rep(0,length(tmp$names)), 1:length(tmp$names), type='n', xlim=c(0,4), ylim=c(0,length(tmp$names)+1), xaxt='n', yaxt='n' )
mtext("Legend", side=3)
for ( k in 1:length(tmp$names) ) {
    lines( c(0,1), c(k,k), col=tmp$cols[k], lwd=3 )
    text( 1, k, labels=tmp$names[k], pos=4 )
}
dev.off()

###
# Decay with distance
load("lentab.Rdata")

newcat <- rep(names(newcats),times=sapply(newcats,length)); names(newcat) <- unlist(newcats)
newcat <- newcat[levels(lentab$countryX)]
tmp.X <- newcat[lentab$countryX]
tmp.Y <- newcat[lentab$countryY]
lentab$smcat <- as.factor( paste( ifelse(tmp.X<tmp.Y,tmp.X,tmp.Y), ifelse(tmp.X<tmp.Y,tmp.Y,tmp.X), sep="-" ) )
levels( lentab$smcat ) <- c(
        "E-E"="E-E", 
        "N-N"="N-N", 
        "W-W"="W-W",
        "E-N"="between E,N,W", 
        "E-W"="between E,N,W", 
        "N-W"="between E,N,W", 
        "E-I"="I-(I,E,N,W)", 
        "I-I"="I-(I,E,N,W)", 
        "I-N"="I-(I,E,N,W)", 
        "I-W"="I-(I,E,N,W)", 
        "E-TC"="TC-any", 
        "I-TC"="TC-any", 
        "N-TC"="TC-any", 
        "TC-TC"="TC-any", 
        "TC-W"="TC-any"
    )[ levels(lentab$smcat) ]
# smcat.cols <- c( "E-E"="#66C2A5", "I-(E,N,W)"="#FC8D62", "between E,N,W"="#8DA0CB", "TC-any"="#E78AC3", "N-N"="#A6D854", "W-W"="#FFD92F" )
smcat.cols <- rainbow_hcl(nlevels(lentab$smcat), c=90); names(smcat.cols) <- levels(lentab$smcat)
catcols <- c(E="#66C2A5", W="#FC8D62", I="#8DA0CB", TC="#E78AC3", N="#A6D854")
# brewer.pal(n=5,name="Set2")
# catcols <- c(E="red", W="green", I="orange", TC="magenta", N="blue")
ccatcols <- catcols[newcat]; names(ccatcols) <- names(newcat)

if (FALSE) {
    # old categories
    lentab$smcat <- lentab$catXY
    levels( lentab$smcat ) <- c( "E-E"="within NEastern", "E-M"="to Turkey", "E-N"="within NEastern", "E-S"="to Italy except Turkey", "E-W"="Western-NEastern", "M-M"="to Turkey", "M-N"="to Turkey", "M-S"="to Turkey", "M-W"="to Turkey", "N-N"="within NEastern", "N-S"="to Italy except Turkey", "N-W"="Western-NEastern", "S-S"="to Italy except Turkey", "S-W"="to Italy except Turkey", "W-W"="within Western")[ levels(lentab$catXY) ]
    smcat.cols <- c( "within Eastern"="red", "to Turkey"="magenta", "to Italy except Turkey"="orange", "within Western and to Eastern"="green" )[levels(lentab$smcat)]
}

# tmp <- with( subset(lentab, countryX<countryY), data.frame( countryX=countryX, countryY=countryY, maplen=maplen, ibd=ibd, gdist=gdist, smcat=smcat, catXY=catXY, npairs=npairs, cex=cex ) )
# levels( tmp$maplen ) <- unlist(  list( "(1,2]"="(1,3]", "(2,2.5]"="(1,3]", "(2.5,3]"="(1,3]", "(3,3.5]"="(3,5]", "(3.5,4]"="(3,5]", "(4,5]"="(3,5]", "(5,6]"=">5", "(6,8]"=">5", "(8,100]"=">5" )[levels(tmp$maplen)] )
# tmp <- ddply( tmp, c("countryX","countryY","maplen"), summarise, countryX=countryX[1], countryY=countryY[1], ibd=sum(ibd), gdist=gdist[1], smcat=smcat[1], catXY=catXY[1], maplen=maplen[1], npairs=npairs[1], cex=cex[1] )
# xx <- seq(0,4000,length.out=101)

tmp <- with( subset(lentab, countryX<countryY), data.frame( countryX=countryX, countryY=countryY, maplen=maplen, ibd=ibd, gdist=gdist, smcat=smcat, npairs=npairs, cex=cex ) )
# levels( tmp$maplen ) <- unlist(  list( "(1,2]"="(1,3]", "(2,2.5]"="(1,3]", "(2.5,3]"="(1,3]", "(3,3.5]"="(3,5]", "(3.5,4]"="(3,5]", "(4,5]"="(3,5]", "(5,6]"=">5", "(6,8]"=">5", "(8,100]"=">5" )[levels(tmp$maplen)] )
tmp <- ddply( tmp, c("countryX","countryY","maplen"), summarise, countryX=countryX[1], countryY=countryY[1], ibd=sum(ibd), gdist=gdist[1], smcat=smcat[1], maplen=maplen[1], npairs=npairs[1], cex=cex[1] )
xx <- seq(0,4000,length.out=101)

# COMBINED FOR PAPER
pdf(file="sharing-rates-and-maps.pdf", width=7, height=7, pointsize=10)
layout( matrix(1:9,nrow=3,byrow=TRUE), heights=c(1,1,2), widths=c(1.2,1,1) )
# layout( rbind( cbind( matrix(rep(1:6,each=3),nrow=2,byrow=T), matrix(7,nrow=2,ncol=3) ), rep(8:10,each=4) ), heights=c(1,1,2), widths=c(1.2,1,1) )
# easter eggs
opar <- par(mar=c(0,0,0,0)+.1)
theseones <- c("Ireland","Sweden","Poland","France","Italy","Serbia")
for (k in seq_along(theseones)) {
    mplot( theseones[k], scale=15, label=paste( LETTERS[k], ") ", theseones[k], sep=""), cols=adjustcolor(ccatcols,.75) )
}
# decay w/ distance
for (k in 1:nlevels(tmp$maplen)) {
    par( mar=c( 4, ifelse(k%%3==1,4,0), 2, 1 )+.1 )
    glmfit <- glm( ibd ~ gdist * smcat, family='poisson', weights=npairs, data=tmp, subset=as.numeric(maplen)==k )
    plot( ibd ~ gdist, data=tmp, subset=(as.numeric(maplen)==k), log='y', pch=21, cex=pmin(cex,3), col=adjustcolor(smcat.cols[smcat],.8), bg=adjustcolor(smcat.cols[smcat],.4), ylim=c(.003,2), yaxt=ifelse(k%%3==1,"s","n"), xlab=ifelse(k%%3==2,"dist (km)",""), ylab=ifelse(k%%3==1,"# blocks per pair",""), main=paste(LETTERS[6+k],") ",levels(maplen)[k],"cM",sep="") )
    for (cxy in (1:nlevels(tmp$smcat))[table(tmp$smcat)>50] ) { 
        lines( xx, exp( predict( glmfit, newdata=list(gdist=xx,smcat=rep(levels(tmp$smcat)[cxy],length(xx))) ) ), col=smcat.cols[cxy] )
    }
    if (k%%3==0) { legend("topright",legend=levels(tmp$smcat),col=smcat.cols,pt.bg=adjustcolor(smcat.cols,.3),pch=21,cex=10/10,pt.cex=2) }
}
par(opar)
dev.off()

# SVG VERSION
require("RSVGTipsDevice")
devSVGTips(file="sharing-rates.svg", width=10, height=7.5, toolTipMode=1, title="Number of blocks shared, by geographic distance")
layout( matrix(1:3,nrow=1,byrow=TRUE), widths=c(1.2,1,1) )
# decay w/ distance
for (k in 1:nlevels(tmp$maplen)) {
    par( mar=c( 4, ifelse(k%%3==1,4,0), 2, 1 )+.1 )
    glmfit <- glm( ibd ~ gdist * smcat, family='poisson', weights=npairs, data=tmp, subset=as.numeric(maplen)==k )
    plot( ibd ~ gdist, data=tmp, subset=(as.numeric(maplen)==k), log='y', pch=21, cex=pmin(cex,3), col=adjustcolor(smcat.cols[smcat],.8), bg=adjustcolor(smcat.cols[smcat],.4), ylim=c(.003,2), yaxt=ifelse(k%%3==1,"s","n"), xlab=ifelse(k%%3==2,"dist (km)",""), ylab=ifelse(k%%3==1,"# blocks per pair",""), main=paste(LETTERS[6+k],") ",levels(maplen)[k],"cM",sep=""), type='n' )
    with(subset(tmp,(as.numeric(maplen)==k)), legend.svg( gdist, ibd, labels=paste(countryX,countryY,sep="-"), pch=21, cex=pmin(cex,3), col=adjustcolor(smcat.cols[smcat],.8), bg=adjustcolor(smcat.cols[smcat],.4) ) )
    for (cxy in (1:nlevels(tmp$smcat))[table(tmp$smcat)>50] ) { 
        lines( xx, exp( predict( glmfit, newdata=list(gdist=xx,smcat=rep(levels(tmp$smcat)[cxy],length(xx))) ) ), col=smcat.cols[cxy] )
    }
    if (k%%3==0) { legend("topright",legend=levels(tmp$smcat),col=smcat.cols,pt.bg=adjustcolor(smcat.cols,.3),pch=21,cex=10/10,pt.cex=2) }
}
dev.off()


# Old version
pdf(file="nblock-by-gdist.pdf",width=7,height=4,pointsize=10)
layout(matrix(c(rep(1:3,each=2),4,5),nrow=2),widths=c(1.2,1.2,1,1))
# correlation plot
opar <- par(mar=c(2,2,2,1)+.1)
with( subset(lentab,maplen==levels(lentab$maplen)[2]), hplot( x=match(countryY,unlist(catlist)), y=match(countryX,unlist(catlist)), z=ccor, max.cex=2, main="correlations" ) )
for (k in 1:length(catlist)) {
    axis(1, at=which(unlist(catlist)%in%catlist[[k]]), labels=countryabbrevs[catlist[[k]]], las=2, cex.axis=8/10, col.axis=ccatcols[catlist[[k]][1]])
    axis(2, at=which(unlist(catlist)%in%catlist[[k]]), labels=countryabbrevs[catlist[[k]]], las=2, cex.axis=8/10, col.axis=ccatcols[catlist[[k]][1]])
}
# map of categories 
euplot( x=sapply(newcat,is.character), cols=ccatcols )
legend("left",pch=21,pt.cex=2,pt.bg=adjustcolor(catcols,.5),legend=names(catcols),cex=8/10)
for (k in 1:nlevels(tmp$maplen)) {
    par( mar=c( 4, ifelse(k%%3==1,4,0), 2, 1 )+.1 )
    glmfit <- glm( ibd ~ gdist * smcat, family='poisson', weights=npairs, data=tmp, subset=as.numeric(maplen)==k )
    plot( ibd ~ gdist, data=tmp, subset=(as.numeric(maplen)==k), log='y', pch=21, cex=pmin(cex,3), col=adjustcolor(smcat.cols[smcat],.5), bg=adjustcolor(smcat.cols[smcat],.2), ylim=c(.003,2), yaxt=ifelse(k%%3==1,"s","n"), xlab=ifelse(k%%3==2,"dist (km)",""), ylab=ifelse(k%%3==1,"# blocks per pair",""), main=paste(levels(maplen)[k],"cM",sep="") )
    # checkthese <- c("Sweden")
    # with( subset(tmp, (as.numeric(maplen)==k) & ( (countryX%in%checkthese) | (countryY%in%checkthese) )), text( gdist, ibd, labels=countryabbrevs[ifelse(countryX%in%checkthese,countryY,countryX)], cex=.5 ) )
    for (cxy in 1:nlevels(tmp$smcat)) { 
        lines( xx, exp( predict( glmfit, newdata=list(gdist=xx,smcat=rep(levels(tmp$smcat)[cxy],length(xx))) ) ), col=smcat.cols[cxy] )
    }
    if (k%%3==0) { legend("topright",legend=levels(tmp$smcat),col=smcat.cols,pt.bg=adjustcolor(smcat.cols,.3),pch=21,cex=8/10,pt.cex=2) }
    # put on scale label
    segments( x0=3500, y0=.003, y1=.003*exp(-1000*(coef(glmfit)["gdist"]+mean(coef(glmfit)[grep("gdist:",names(coef(glmfit)))]))), lwd=2 )
    text( x=3500, y=.0035, labels="1000km", pos=2 )
}
par(opar)
dev.off()


# For grant
tmp <- with( subset(lentab, countryX<countryY), data.frame( countryX=countryX, countryY=countryY, maplen=maplen, ibd=ibd, gdist=gdist, smcat=smcat, catXY=catXY, npairs=npairs, cex=cex ) )
levels( tmp$maplen ) <- unlist(  list( "(1,2]"="(1,5]", "(2,2.5]"="(1,5]", "(2.5,3]"="(1,5]", "(3,3.5]"="(1,5]", "(3.5,4]"="(1,5]", "(4,5]"="(1,5]", "(5,6]"=">5", "(6,8]"=">5", "(8,100]"=">5" )[levels(tmp$maplen)] )
tmp <- ddply( tmp, c("countryX","countryY","maplen"), summarise, countryX=countryX[1], countryY=countryY[1], ibd=sum(ibd), gdist=gdist[1], catXY=catXY[1], maplen=maplen[1], npairs=npairs[1], cex=cex[1] )
tmp$smcat <- tmp$catXY
levels( tmp$smcat ) <- c( "E-E"="within NEastern", "E-M"="to Southern", "E-N"="within NEastern", "E-S"="to Southern", "E-W"="Other", "M-M"="to Southern", "M-N"="to Southern", "M-S"="to Southern", "M-W"="to Southern", "N-N"="within NEastern", "N-S"="to Southern", "N-W"="Other", "S-S"="to Southern", "S-W"="to Southern", "W-W"="Other")[ levels(tmp$catXY) ]
smcat.cols <- c( "within NEastern"="magenta", "to Southern"="blue", "Other"="green" )[levels(tmp$smcat)]
xx <- seq(0,4000,length.out=101)

pdf(file="nblock-by-gdist-grant.pdf",width=4,height=2.5,pointsize=10)
layout(t(1:2),widths=c(1.25,1))
for (k in 1:nlevels(tmp$maplen)) {
    par( mar=c( 3.1, ifelse(k%%nlevels(tmp$maplen)==1,2.8,0), 2, 1 )+.1, cex.axis=8/10, cex.lab=8/10, cex.main=10/10, mgp=c(2,.8,0) )
    glmfit <- glm( ibd ~ gdist * smcat, family='poisson', weights=npairs, data=tmp, subset=as.numeric(maplen)==k )
    plot( ibd ~ gdist, data=tmp, subset=(as.numeric(maplen)==k), log='y', pch=21, cex=pmin(cex,3)/1.7, col=adjustcolor(smcat.cols[smcat],.5), bg=adjustcolor(smcat.cols[smcat],.2), ylim=c(.002,2), yaxt=ifelse(k%%nlevels(tmp$maplen)==1,"s","n"), xlab="dist (km)", ylab=ifelse(k%%nlevels(tmp$maplen)==1,"# blocks per pair",""), main=paste(levels(maplen)[k],"cM",sep="") )
    for (cxy in 1:nlevels(tmp$smcat)) { 
        lines( xx, exp( predict( glmfit, newdata=list(gdist=xx,smcat=rep(levels(tmp$smcat)[cxy],length(xx))) ) ), col=smcat.cols[cxy] )
    }
    if (k%%nlevels(tmp$maplen)==0) { legend("topright",legend=levels(tmp$smcat),col=smcat.cols,pt.bg=adjustcolor(smcat.cols,.3),pch=21,cex=8/10,pt.cex=2) }
}
par(opar)
dev.off()


###
# Substructure

f <- function (a,b,x,y=Inf,margins=TRUE) { 
    # make a table of numbers of blocks shared
    x <- with( subset(blocks, maplen>x & maplen<=y & countrypair==gsub(' ','.',paste(sort(c(a,b)),collapse='-'))),  {
        if ( a==b ) {
            table( factor(id1,levels=indivinfo$SUBJID[indivinfo$COUNTRY_SELF==a]), factor(id2,levels=indivinfo$SUBJID[indivinfo$COUNTRY_SELF==b]) ) + table( factor(id2,levels=indivinfo$SUBJID[indivinfo$COUNTRY_SELF==a]), factor(id1,levels=indivinfo$SUBJID[indivinfo$COUNTRY_SELF==b]) )
        } else {
            table( factor(ifelse(country1==a,id1,id2),levels=indivinfo$SUBJID[indivinfo$COUNTRY_SELF==a]), factor(ifelse(country1==a,id2,id1),levels=indivinfo$SUBJID[indivinfo$COUNTRY_SELF==b]) )
        } } )
    x <- x[order(rowSums(x)),order(colSums(x)),drop=FALSE]; 
    if (margins) { addmargins(x) } else { x } 
}


ff <- function (a,b,x=2,y=Inf,n=30,plotit=FALSE,rescale=FALSE,...) {
    # plot histograms of the margins of that table, with poisson expectations alongside
    fzz <- f(a,b,x,y)
    zz <- rev(rev(sort(fzz[,"Sum"]))[-1])
    fna1 <- hist(zz,breaks=seq(0,length.out=n,by=ceiling(1.2*max(zz)/n)),plot=FALSE)
    zp1 <- with(fna1, nsamples[a]*diff( ppois(q=c(-1,breaks[-1]),lambda=mean(zz)) ) )
    zz <- rev(rev(sort(fzz["Sum",]))[-1])
    fna2 <- hist(zz,breaks=seq(0,length.out=n,by=ceiling(1.2*max(zz)/n)),plot=FALSE)
    zp2 <- with(fna2, nsamples[b]*diff( ppois(q=c(-1,breaks[-1]),lambda=mean(zz)) ) )
    if (rescale) {
        # make x-axis "rate per pair"
        for (f2 in c("mids","breaks")) {
            fna1[[f2]] <- fna1[[f2]]/nsamples[b]
            fna2[[f2]] <- fna2[[f2]]/nsamples[a]
        }
    }
    if (plotit) {
        plot(fna1,main=paste("rate of",b,"IBD across",a)); polygon( fna1$mids, zp1, col=adjustcolor('red',.25) )
        if (a!=b) plot(fna2,main=paste("rate of",a,"IBD across",b)); polygon( fna2$mids, zp2, col=adjustcolor('red',.25) )
    }
    return(invisible(list(hist1=fna1,pois1=zp1,hist2=fna2,pois2=zp2)))
}

# get p-values for differences between (b blocks in a) and (d blocks in c)
compare.means <- function (a,b,c,d,x=2,y=Inf) {
    ab <- f(a,b,x,y,margins=FALSE)  # rows are indivs of a
    if (a==b) { ab <- ab[upper.tri(ab,diag=FALSE)] }
    ab <- as.vector(ab)
    cd <- f(c,d,x,y,margins=FALSE)
    if (c==d) { cd <- cd[upper.tri(cd,diag=FALSE)] }
    cd <- as.vector(cd)
    ans <- c( ab=mean(ab), cd=mean(cd), se=sqrt( (var(ab)+var(cd))/(length(ab)+length(cd)) ), pvalue=wilcox.test(ab,cd)$p.value)
    return(ans)
}


pdf(file="various-substructure.pdf", width=7, height=4, pointsize=10)
layout( matrix(c(1,2,3,3,4,5,6,6),nrow=2), widths=c(1.5,2,1.5,2) )
par(mar=c(4,4,2,1), mgp=c(3,1,0))
with( ff("Italy","Swiss French",x=2,y=6,n=60), {  
        plot(hist1, main="CHf blocks in IT", col=adjustcolor("blue",0.5),ylim=c(0,max(pois1)), xlab="# blocks", ylab='');
        polygon( hist1$mids, pois1, col=adjustcolor('red',0.25) )
    } )
with( ff("Italy","United Kingdom",x=2,y=6,n=60), {  
        plot(hist1, main="UK blocks in IT", col=adjustcolor("blue",0.5),ylim=c(0,max(pois1)), xlab="# blocks", ylab='');
        polygon( hist1$mids, pois1, col=adjustcolor('red',0.25) )
    } )
# plist <- c("France","Switzerland","Italy","Greece","Turkey","Cyprus")
plist <- c("France", "Italy", "Greece", "Turkey", "Cyprus")
pcols <- countrycols; pcols[plist] <- rainbow_hcl(length(plist),l=60,c=70)
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( United.Kingdom.blocks ~ Swiss.French.blocks, pch=21, xlab="# blocks with French-speaking Swiss (CHf)", ylab="# blocks with United Kingdom (UK)", cex=1.5, xlim=c(150,700), ylim=c(55,250), subset=(COUNTRY_SELF=="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), col=NA, main="(A)",  xpd=NA )
        points( United.Kingdom.blocks ~ Swiss.French.blocks, cex=1.5, subset=(COUNTRY_SELF!="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), pch=ifelse(COUNTRY_SELF%in%c("Switzerland","France"),22,23), col=adjustcolor("black",.25) )
    } )
legend("topright", legend=countryabbrevs[plist], pch=c(22,21,23,23,23), pt.cex=1.5, bg="white", pt.bg=pcols[plist], col=ifelse(plist=="Italy",NA,"black") )
clist <- c("Ireland","Germany","Italy","Poland","Slovakia")
tmp <- lapply( clist, function (b)  ff("United Kingdom",b,x=0,y=Inf,n=60) ); names(tmp) <- clist
with( tmp[["Ireland"]], {  
        plot(hist1, main="IE blocks in UK", col=adjustcolor("blue",0.5),ylim=c(0,max(pois1)), xlab="# blocks", ylab='');
        polygon( hist1$mids, pois1, col=adjustcolor('red',0.25) )
    } )
with( tmp[["Germany"]], {  
        plot(hist1, main="DE blocks in UK", col=adjustcolor("blue",0.5),ylim=c(0,max(pois1)), xlab="# blocks", ylab='');
        polygon( hist1$mids, pois1, col=adjustcolor('red',0.25) )
    } )
# plist <- c("Poland","Germany","United Kingdom","Ireland")
plist <- c("Germany","United Kingdom","Ireland")
pcols <- countrycols; pcols[plist] <- rainbow_hcl(length(plist),l=60,c=70)
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( Germany.blocks ~ Ireland.blocks, pch=21, xlab="# blocks with Ireland (IE)", ylab="# blocks with Germany (DE)", bg=adjustcolor(pcols[COUNTRY_SELF],.75), cex=1.5, col=NA, subset=(COUNTRY_SELF=="United Kingdom"), main="(B)" ) 
        points( Germany.blocks ~ Ireland.blocks, pch=ifelse(COUNTRY_SELF%in%c("Germany","Poland"),22,23), bg=adjustcolor(pcols[COUNTRY_SELF],.75), col=adjustcolor("black",.25), cex=1.5, subset=(COUNTRY_SELF!="United Kingdom") ) 
    } )
legend("topright",pch=c(22,21,23),pt.bg=pcols[plist],legend=countryabbrevs[plist],pt.cex=1.5,col=ifelse(plist=="United Kingdom",NA,"black") )
# with( subset(indivinfo,COUNTRY_SELF=="United Kingdom"&Italy.blocks>110), {
#         text( Ireland.blocks, Germany.blocks, labels="*", cex=2.5 )
#         # points( Ireland.blocks, Germany.blocks, col=adjustcolor(pcols["Italy"],.5), pch=20, cex=1.5 )
#         text( mean(Ireland.blocks), mean(Germany.blocks)-13, labels=">10 blocks with IT", pos=4, cex=.75 )
#     } )
with( subset(indivinfo,COUNTRY_SELF=="United Kingdom"&Slovakia.blocks>20), {
        # points( Ireland.blocks, Germany.blocks, col=adjustcolor(pcols["Slovakia"],.5), cex=1.5, pch=20)
        text( mean(Ireland.blocks), mean(Germany.blocks)+2, labels=countryabbrevs["Slovakia"], col=countrycols["Slovakia"], cex=8/10 )
    } )
dev.off()

with( subset(indivinfo,COUNTRY_SELF%in%c("Italy","Greece","Turkey","Cyprus","France","Switzerland","Albania","Kosovo")), {
        plot( United.Kingdom.blocks ~ France.blocks, pch=20, xlab="# blocks with CHf", ylab="# blocks with UK", cex=2, col=adjustcolor(countrycols[COUNTRY_SELF],.25), subset=COUNTRY_SELF=="Italy" ) 
        points( United.Kingdom.blocks ~ France.blocks, pch=20, cex=2, col=adjustcolor(countrycols[COUNTRY_SELF],.25), subset=COUNTRY_SELF!="Italy" ) 
        # identify( France.blocks, United.Kingdom.blocks, labels=COUNTRY_SELF )
    } )
with( subset(indivinfo,COUNTRY_SELF%in%c("Albania","Kosovo","Yugoslavia","Montenegro","Bosnia","Serbia")), {
        points( United.Kingdom.blocks ~ France.blocks, pch=19, cex=2, col=adjustcolor(countrycols[COUNTRY_SELF],.25) ) 
        identify( France.blocks, United.Kingdom.blocks, labels=COUNTRY_SELF )
    } )

# for grant
clist <- c("Ireland","Germany")
tmp <- lapply( clist, function (b)  ff("United Kingdom",b,x=0,y=Inf,n=60) ); names(tmp) <- clist
pdf(file="various-substructure-grant.pdf", width=3.5, height=3, pointsize=10)
layout( matrix(c(1,2,3,3),nrow=2), widths=c(1.5,2) )
par(mar=c(2,2,2,1), mgp=c(1,1,0))
with( tmp[["Ireland"]], {  
        plot(hist1, main="IE blocks in UK", col=adjustcolor("blue",0.5),ylim=c(0,max(pois1)), xlab="", ylab='');
        polygon( hist1$mids, pois1, col=adjustcolor('red',0.25) )
    } )
with( tmp[["Germany"]], {  
        plot(hist1, main="DE blocks in UK", col=adjustcolor("blue",0.5),ylim=c(0,max(pois1)), xlab="", ylab='');
        polygon( hist1$mids, pois1, col=adjustcolor('red',0.25) )
    } )
plist <- c("Germany","United Kingdom","Ireland")
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( Germany.blocks ~ Ireland.blocks, pch=21, xlab="# blocks with Ireland", ylab="# blocks with Germany", xaxt='n', yaxt='n', bg=adjustcolor(countrycols[COUNTRY_SELF],.75), cex=1.5, col=NA, subset=(COUNTRY_SELF=="United Kingdom") ) 
        points( Germany.blocks ~ Ireland.blocks, pch=ifelse(COUNTRY_SELF%in%c("Germany","Poland"),22,23), bg=adjustcolor(countrycols[COUNTRY_SELF],.75), col=adjustcolor("black",.25), cex=1.5, subset=(COUNTRY_SELF!="United Kingdom") ) 
    } )
legend("topright",pch=c(22,22,21,23),pt.bg=countrycols[plist],legend=countryabbrevs[plist],pt.cex=1.5,col=ifelse(plist=="United Kingdom",NA,"black") )
dev.off()


# for talk: self rates
selfrates <- with( subset(blocks, country1==country2), table(country1)/choose(nsamples[levels(country1)],2) )
pdf(file="selfrates.pdf", width=7, height=7, pointsize=10)
par(mar=c(0,0,1.5,0));
euplot( x=selfrates, legend=TRUE )
dev.off()

# Significant substructure?  Do a permutation test.
if (!file.exists("signif.RData")) {
    # temporary computations
    minlen <- 1
    zzz <- with( subset(blocks,maplen>minlen), table( factor(id1,levels=(indivinfo$SUBJID)), country2 ) )
    zzz <- zzz + with( subset(blocks,maplen>minlen), table( factor(id2,levels=(indivinfo$SUBJID)), country1 ) )
    nbl <- ( tapply( 1:dim(zzz)[1], indivinfo$COUNTRY_SELF, function (kk) { colSums( zzz[kk,,drop=FALSE] ) } ) )  # sum of numbers of blocks shared
    ynames <- list( names(nbl[[1]]), names(nbl) )
    nbl2 <- ( tapply( 1:dim(zzz)[1], indivinfo$COUNTRY_SELF, function (kk) { colSums( zzz[kk,,drop=FALSE]^2 ) } ) ) # sum of squared numbers
    signif <- unlist(nbl); dim(signif) <- sapply(ynames,length); dimnames(signif) <- ynames
    signif <- as.data.frame.table( signif ); colnames(signif) <- c("countryY","countryX","sumbl")
    signif <- signif[ c("countryX","countryY","sumbl")]
    signif$countryX <- ordered(signif$countryX,levels=levels(blocks$country1))
    signif$countryY <- ordered(signif$countryY,levels=levels(blocks$country1))
    signif$sumbl2 <- unlist(nbl2)
    signif$npairs <- with(signif, ifelse( countryX==countryY, choose(nsamples[countryX],2), nsamples[countryX]*nsamples[countryY] ) )
    signif$nblocks <- with(signif, ifelse(countryX==countryY,1/2,1) * sumbl )
    signif$ibd <- with(signif, nblocks / npairs )
    # Mean and SD of numbers of blocks that X-indivs share with country Y
    #  ... divide these by nsamples[countryY] to get naturally comparable things between countries.
    signif$meanblocks <- with(signif, sumbl / nsamples[countryX] )
    signif$sdblocks <- with(signif, sqrt( (1/(nsamples[countryX]-1)) * ( sumbl2 - sumbl^2 / nsamples[countryX] ) ) )
    signif$sdibd <- signif$sdblocks / nsamples[signif$countryY]
    # Significance for substructure
    signif$p.sd <- NA  # a p-value
    signif$z.sd <- NA  # a z-score
    for (c1 in levels(blocks$country1)) {
        all.ids <- indivinfo$SUBJID[indivinfo$COUNTRY_SELF==c1]
        if (length(all.ids)>1) {
            subblocks <- with( subset( blocks, (country1==c1 | country2==c1) & country1!=country2 & maplen>minlen ),
                    data.frame( idA=ifelse(country1==c1,id1,id2), 
                        idB=ifelse(country1==c1,id2,id1), 
                        countryB=factor( levels(blocks$country1)[ifelse(country1==c1,country2,country1)], levels=levels(blocks$country1) )
                        )
                    )
            boots <- replicate( 1000, {
                    subblocks$idA <- sample( all.ids, nrow(subblocks), replace=TRUE )
                    zzz <- with( subblocks, table( list( factor(idA,levels=all.ids), countryB ) ) )
                    ( colSums( zzz^2, dims=1 ) ) # dimensions of this are (countryY,maplen)
                } )
            signif$p.sd[ signif$countryX==c1 ] <- rowMeans( boots>=signif$sumbl2[signif$countryX==c1] )  # right tail p-value
            signif$z.sd[ signif$countryX==c1 ] <- (signif$sumbl2[signif$countryX==c1] - rowMeans( boots ) ) / sqrt( rowMeans( sweep(boots,1,rowMeans(boots),"-")^2 ) )  # right tail "z-value"
        }
    }
    attr(signif,"minlen") <- minlen
    save(signif,file="signif.RData")
} else {
    load("signif.RData")
}

# Summary of how much substructure there is -- for paper
pdf(file="substructure_summaries.pdf", width=6.5, height=6.5, pointsize=10)
# layout(t(1:2), widths=c(1.8,1))
layout(1:2, heights=c(1,3))
par(mar=c(5,4,1,1)+.1)
with( subset(signif,nsamples[countryY]>10), hist( p.sd, breaks=30, xlab='p-value', main="A) Significance of substructure") )
with( signif, rowplot( z.sd, rlist=levels(countryX)[nsamples[levels(countryX)]>3], clist=levels(countryX)[nsamples[levels(countryX)]>10], x=countryX, y=countryY, log="x", abbrevs=countryabbrevs, oneline=TRUE, xlab="z score", main="B) Degree of substructure" ) )
abline(v=5, col="grey")
# with( signif, rowplot( p.sd, rlist=levels(countryX)[nsamples[levels(countryX)]>3], clist=levels(countryX)[nsamples[levels(countryX)]>24], x=countryX, y=countryY, log="x", abbrevs=countryabbrevs, oneline=TRUE, ylab="none", xlab="p-value", main="B" ) )
dev.off()

# Just the correlations, for the paper (supplement)
pdf( file="sharing-correlations.pdf", width=6.5, height=8, pointsize=10 )
layout(cbind( matrix(1:6,nrow=3,byrow=TRUE),rep(7,3) ), widths=c(3,3,1) )
for (ml in c(1,3,5,6,7,8)) {
    # pdf( file=paste("sharing-correlations",gsub("[][(),.]","-",levels(lentab$maplen)[ml]),".pdf",sep=""), width=5, height=4, pointsize=10 )
    # layout(t(1:2),widths=c(3,1))
    par(mar=c(3,3,3,1)+.1)
    with( subset(lentab,maplen==levels(lentab$maplen)[ml]), hplot( x=match(countryY,unlist(newcats)), y=match(countryX,unlist(newcats)), z=ccor, max.cex=1, main=paste("correlations,",levels(lentab$maplen)[ml], "cM"), alpha=1 ) )
    for (k in 1:length(newcats)) {
        axis(1, at=which(unlist(newcats)%in%newcats[[k]]), labels=countryabbrevs[newcats[[k]]], las=2, cex.axis=8/10, col.axis=ccatcols[newcats[[k]][1]])
        axis(2, at=which(unlist(newcats)%in%newcats[[k]]), labels=countryabbrevs[newcats[[k]]], las=2, cex.axis=8/10, col.axis=ccatcols[newcats[[k]][1]])
    }
    abline(h=0.5+cumsum(sapply(newcats,length))[-length(newcats)],v=0.5+cumsum(sapply(newcats,length))[-length(newcats)])
    # plot( rep(0,17), seq(-1,1,length.out=17), col=hcolor(seq(-1,1,length.out=17), alpha=1), xaxt='n', ylab='correlation', pch=20, cex=4*sqrt(abs(seq(-1,1,length.out=17))) )
    # dev.off()
}
plot( rep(0,33), seq(-1,1,length.out=33), col=hcolor(seq(-1,1,length.out=33), alpha=1), xaxt='n', ylab='correlation', pch=20, cex=4*sqrt(abs(seq(-1,1,length.out=33))) )
dev.off()


####
# Triples
load("poptriples.Rdata")

# Plot ALL triples
clist <- unlist(newcats)[nsamples[unlist(newcats)]>1]
pdf(file="all-triples.pdf", width=6, height=12, pointsize=10)
for (x in clist) {
    rlist <- with(subset(triples,countryA==x&countryB%in%clist&countryC%in%clist), { trate <- rate-pairrate; c( countryB[!is.na(trate)&trate>0], countryC[!is.na(trate)&trate>0] ) } )
    if (length(rlist)>0) {
        with(subset(triples,countryA==x), rowplot( x=c(match(countryB,clist),match(countryC,clist)), y=c(match(countryC,clist),match(countryB,clist)), rmatrix=rep(rate-pairrate,2), rlist=clist, clist=clist, abbrevs=countryabbrevs[clist], main=x, oneline=TRUE, log='x', xlim=c(1e-4,.033) ) )
        abline(v=c(5e-4,1e-3,2e-3),col=adjustcolor("black",.4))
    }
}
dev.off()

# Compare to normalization
with( triples, plot( rate, (poppairs$nblocks/poppairs$npairs)[ match(countrypair1,poppairs$countrypair) ] + (poppairs$nblocks/poppairs$npairs)[ match(countrypair2,poppairs$countrypair) ] + (poppairs$nblocks/poppairs$npairs)[ match(countrypair3,poppairs$countrypair) ] ) )

aab.corrected <- (1/2)*( (polarized(triples,pattern="aab",minsamples=5)-(1/2)*polarized(triples,pattern="aab",varname="pairrate",minsamples=5)) + (polarized(triples,pattern="baa",minsamples=5)-(1/2)*polarized(triples,pattern="baa",varname="pairrate",minsamples=5)) )

pdf(file="triple-rates.pdf", width=7, height=4, pointsize=10)
# Not normalized
rowplot(t(aab.corrected), clist=clist, log='x', xlab="rate per xxy triple",oneline=TRUE,abbrevs=countryabbrevs)
dev.off()

rowplot(t(aab.corrected)-aab.corrected, clist=clist, log='', xlab="rate difference xxy-yyx per triple",oneline=TRUE,abbrevs=countryabbrevs)

# Normalized by AAA rate
rowplot( t( (polarized(triples,pattern="aab",minsamples=5)-(1/2)*polarized(triples,pattern="aab",varname="pairrate",minsamples=5))/(three.rate(triples,minsamples=5)-(1/2)*three.rate(triples,varname="pairrate",minsamples=5)) ), clist=clist, log='x', xlab="rate per triple",oneline=TRUE,abbrevs=countryabbrevs)
rowplot( t( (polarized(triples,pattern="baa",minsamples=5)-(1/2)*polarized(triples,pattern="baa",varname="pairrate",minsamples=5))/(three.rate(triples,minsamples=5)-(1/2)*three.rate(triples,varname="pairrate",minsamples=5)) ), clist=clist, log='x', xlab="rate per triple",oneline=TRUE,abbrevs=countryabbrevs)


#####
# Table of rates
tmp <- with( indivinfo, data.frame( 
        abbrev=countryabbrevs[levels(COUNTRY_SELF)], 
        nsamples=as.numeric(table(COUNTRY_SELF)), 
        self.ibd=with(subset(poppairs,country1==country2),(nblocks/npairs)[match(levels(COUNTRY_SELF),country1)]),   # self rate
        other.ibd= sapply( levels(COUNTRY_SELF), function (x) with(subset(poppairs,(country1==x | country2==x) & country1!=country2), mean(nblocks/npairs) ) )  # mean rate across other pops
        # total.ibd=as.numeric(table(blocks$country1)+table(blocks$country2))/as.numeric(table(COUNTRY_SELF))  # total number per indiv
    ) )
newcats <- list( 
        E=c( "Slovakia", "Greece", "Yugoslavia", "Albania", "Bosnia", "Montenegro", "Macedonia", "Kosovo", "Serbia", "Bulgaria", "Romania", "Poland", "Hungary", "Czech Republic", "Russia", "Slovenia", "Ukraine", "Croatia", "Austria"),
        TC=c("Turkey","Cyprus"),
        N=c( "Sweden", "Norway", "Denmark", "Latvia", "Finland" ),
        W=c("France", "United Kingdom", "Scotland", "England", "Ireland", "Swiss German", "Swiss French", "Switzerland", "Belgium", "Netherlands", "Germany" ),
        I=c("Italy","Spain","Portugal")
    )
tmp <- tmp[ match(unlist(newcats),rownames(tmp)), ]

#######
## sample sizes

longs=indivinfo$long[match(levels(indivinfo$COUNTRY_SELF),indivinfo$COUNTRY_SELF)]
lats=indivinfo$lat[match(levels(indivinfo$COUNTRY_SELF),indivinfo$COUNTRY_SELF)]
orig.xy <- xy <- mapproject( longs, lats )  # note: uses previous projection, passing this in messes it up.
eps <- 1.2e-2
for (k in 1:100) {
    dists <- sqrt(outer(xy$x,xy$x,"-")^2 + outer(xy$y,xy$y,"-")^2 )
    diag(dists) <- Inf
    if (all(dists>=.9*eps)) { break }
    ij <- arrayInd(which.min(dists),.dim=dim(dists))
    x0 <- mean(xy$x[ij])
    dx <- diff(xy$x[ij])
    y0 <- mean(xy$y[ij])
    dy <- diff(xy$y[ij])
    if (dx^2+dy^2==0) { dy <- 1 }
    xy$x[ij] <- x0 + (1/2)*c(1,-1)*eps*dx/sqrt(dx^2+dy^2)
    xy$y[ij] <- y0 + (1/2)*c(1,-1)*eps*dy/sqrt(dx^2+dy^2)
    cat(ij, k, dists[ij[1],ij[2]], sqrt(diff(xy$x[ij])^2+diff(xy$y[ij])^2), "\n" )
}

pdf(file="sample-sizes.pdf",width=3, height=3)
par(mai=c(0,0,0,0))
map("world",xlim=c(-10,38), ylim=c(35,61),proj="globular",col=grey(.50),mar=c(1,1,1,1),resolution=0 )
text( xy$x, xy$y, labels=nsamples, cex=.5, col='red' )
# arrows( x0=xy$x, x1=orig.xy$x, y0=xy$y, y1=orig.xy$y, length=0 )
dev.off()
