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
source("ibd-blocks-fns.R")
source("laplace-inversion-fns.R")


# Actual blocks and other metainformation
load("all-blocks-winnowed-fine.Rdata")
load("eda-data-fine.Rdata")
# make a variable which is country-pair
blocks$countrypair <- countrypairs[cbind(as.numeric(blocks$country1),as.numeric(blocks$country2))]
# and individual-pair
blocks$indivpair <- factor( paste( pmin(blocks$id1,blocks$id2), pmax(blocks$id1,blocks$id2), sep="-" ) )

####
# PCs of individual-by-individual sharing matrix
require(spam)
require(vegan)

# positions of means
long.lat <- with(indivinfo, cbind( tapply(long,COUNTRY_SELF,mean), tapply(lat,COUNTRY_SELF,mean))  ) 

na.zero <- function (x) { x[is.na(x)] <- 0; x }

tsample <- function (x,n) { if(length(x)==1 & n>0) { x } else { sample(x,min(length(x),n)) } }

get.pcs <- function (countrylist, rotate=NULL, rescale=!is.null(rotate), nsubsamp=1000, renorm="none", recenter=FALSE, eps=0) {
    subsamp.indivs <- unlist( lapply( countrylist, function (country) tsample(indivinfo$SUBJID[indivinfo$COUNTRY_SELF==country], nsubsamp) ) )
    subsamp.indivinfo <- droplevels( indivinfo[indivinfo$SUBJID%in%subsamp.indivs,] )
    subsamp.indpairs <- subset( indpairs, (id1 %in% subsamp.indivinfo$SUBJID) & (id2 %in% subsamp.indivinfo$SUBJID) )
    # long blocks
    ind.by.ind <- with(subsamp.indpairs, list( values=nblocks1cM+nblocks5cM+nblocks10cM, i=factor(id1,levels=subsamp.indivinfo$SUBJID), j=factor(id2,levels=subsamp.indivinfo$SUBJID) ) )
    if (renorm=="indiv") {
        # Normalize by each individual's marginal amount of sharing
        ind.marg <- with(ind.by.ind,  na.zero(tapply( values, i, sum )) +  na.zero(tapply( values, j, sum )) )
        # geometric
        ind.by.ind$values <- log(1+ with(ind.by.ind, values / sqrt( ind.marg[as.numeric(i)] * ind.marg[as.numeric(j)] ) ) )
        # arithmetic
        # ind.by.ind$values <- with(ind.by.ind, values / ( ind.marg[as.numeric(i)] + ind.marg[as.numeric(j)] ) )
    } else if (renorm=="country") {
        # Normalize so each country gets weight one (using eps to not upweight singles too much)
        countrysizes <- pmax(eps, with( subsamp.indivinfo, nsamples[ COUNTRY_SELF ] ))
        ind.by.ind$values <- with( ind.by.ind, values / sqrt( countrysizes[as.numeric(i)] * countrysizes[as.numeric(j)] ) )
    }
    ind.by.ind <- spam( x=ind.by.ind, nrow=nrow(subsamp.indivinfo), ncol=nrow(subsamp.indivinfo) )
    ind.by.ind <- ind.by.ind + t(ind.by.ind)
    if (recenter) {
        # center the matrix
        ind.by.ind <- ind.by.ind - apply(ind.by.ind,1,sum)[row(ind.by.ind)]/nrow(ind.by.ind) - apply(ind.by.ind,2,sum)[col(ind.by.ind)]/ncol(ind.by.ind) + sum(ind.by.ind)/prod(dim(ind.by.ind))
    }
    # dimnames(ind.by.ind) <- list( subsamp.indivinfo$SUBJID, subsamp.indivinfo$SUBJID )
    # ibi.pca <- svd(exp(-ind.by.ind))
    ibi.pca <- svd(ind.by.ind)
    xy <- ibi.pca$v[,1:2]
    mean.pca <- do.call( cbind, lapply( 1:10, function (k) tapply( ibi.pca$v[,k], subsamp.indivinfo$COUNTRY_SELF, mean ) ) )
    long.lat <- with(subsamp.indivinfo, cbind( tapply(long,COUNTRY_SELF,mean), tapply(lat,COUNTRY_SELF,mean))  ) 
    countries <- with(subsamp.indivinfo, levels(COUNTRY_SELF)[as.numeric(COUNTRY_SELF)] )
    if (rescale) {
        # match the means
        # rotate <- procrustes( X=long.lat[rownames(mean.pca),], Y=mean.pca[,1:2], scale=TRUE )
        # match the individuals
        rotate <- procrustes( X=long.lat[countries,], Y=xy, scale=TRUE )
        xy <- predict(rotate, newdata=xy)
        mean.pca <- predict(rotate, mean.pca)
    } else if (is.matrix(rotate)) {
        xy <- (xy)%*%rotate
        mean.pca <- (mean.pca)%*%rotate
    } else { rotate=rotate }
    return( list( xy=xy, countries=countries, mean.pca=mean.pca, rotate=rotate, v=ibi.pca$v[,1:10] ) )
}

# PCA of sharing matrix
# pcs <- get.pcs( setdiff(levels(indivinfo$COUNTRY_SELF),c("Yugoslavia","Albania","Spain","Portugal")), rescale=FALSE, rotate=matrix(c(0,-1,-1,0),nrow=2), renorm=TRUE, nsubsamp=10 )
# mean rates of sharing for countries within distance x
nearby <- function (x) sapply( names(countrycols), function(ccc) with( subset(poppairs, (country1==ccc | country2==ccc) & gdist <= x ), mean( (nblocks1cM+nblocks5cM+nblocks10cM)/npairs ) ) )

pclist <- lapply( list( all=levels(indivinfo$COUNTRY_SELF), noAL=setdiff(levels(indivinfo$COUNTRY_SELF),c("Albania","Kosovo")) ), 
    function (ccc) {
        lapply( c(none="none",country="country"), function (renorm) { get.pcs( ccc, rescale=FALSE, renorm=renorm, nsubsamp=100, eps=0 ) } ) 
    } )

# look at the first few pcs
npcs <- 3
pdf(file="tmp-pcs.pdf",width=2.5*choose(npcs,2),height=10,pointsize=10)
layout( matrix(1:(4*choose(npcs,2)),nrow=4,byrow=TRUE) )
par(mar=c(4,4,0,0))
for ( usethese in c("all","noAL") ) {
    with( pclist[[usethese]][["country"]], {
            for (i in 1:(npcs-1)) for (j in (i+1):npcs) {
                plot( v[,j], v[,i], col=adjustcolor(countrycols[countries],.5), xlab=paste("PC",j), ylab=paste("PC",i), pch=20, cex=1.5 )
                text( mean.pca[,j], mean.pca[,i], labels=countryabbrevs[rownames(mean.pca)], )
                if (i==1 & j==2) textlab("topleft",usethese)
            }
    } ) 
    with( pclist[[usethese]][["none"]], {
            for (i in 1:(npcs-1)) for (j in (i+1):npcs) {
                plot( v[,j], v[,i], col=adjustcolor(countrycols[countries],.5), xlab=paste("PC",j), ylab=paste("PC",i), pch=20, cex=1.5 )
                text( mean.pca[,j], mean.pca[,i], labels=countryabbrevs[rownames(mean.pca)], )
                if (i==1 & j==2) textlab("topleft",usethese)
            }
    } ) 
}
dev.off()

pdf(file="indiv-sharing-map-of-europe.pdf", width=7, height=3, pointsize=10)
layout( matrix(1:3,nrow=1) )
par(mar=c(2,2,1,0)+.1,mgp=c(.5,.5,.5))
with( pclist[["all"]][["country"]], {
        i <- 1; j <- 2
        plot( v[,j], v[,i], col=adjustcolor(countrycols[countries],.5), xlab=paste("PC",j), ylab=paste("PC",i), pch=20, cex=1.5, xaxt='n', yaxt='n' )
        text( mean.pca[,j], mean.pca[,i], labels=countryabbrevs[rownames(mean.pca)], )
        textlab( "topleft", "All countries" )
    } )
with( pclist[["noAL"]][["country"]], {
        i <- 1; j <- 2
        plot( v[,j], v[,i], col=adjustcolor(countrycols[countries],.5), xlab=paste("PC*",j), ylab=paste("PC*",i), pch=20, cex=1.5, xaxt='n', yaxt='n' )
        text( mean.pca[,j], mean.pca[,i], labels=countryabbrevs[rownames(mean.pca)], )
        textlab("topleft", "without AL/KO" )
        i <- 1; j <- 3
        plot( v[,j], v[,i], col=adjustcolor(countrycols[countries],.5), xlab=paste("PC*",j), ylab=paste("PC*",i), pch=20, cex=1.5, xaxt='n', yaxt='n' )
        text( mean.pca[,j], mean.pca[,i], labels=countryabbrevs[rownames(mean.pca)], )
        textlab("topright", "without AL/KO" )
    } )
# map of self rates
# xy <- euplot( nearby( 10 ), scale=6, lab="Self IBD" )
# text( xy, labels=countryabbrevs, cex=8/10 )
dev.off()


