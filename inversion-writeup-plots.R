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

load("eda-data-fine.Rdata")

# minimum length to use
minlen <- 2  # cM

# Get the overall length distribution
load("all-blocks-winnowed-fine.Rdata")
# Exclude chromosome 19 and short blocks
minblocklen <- 1  # note NOT the minimum length to use in discretization
omitchroms <- c( )# omit chr19?
blocks <- subset( blocks, !is.na(blocks$maplen) & blocks$maplen>minblocklen & !(chrom%in%omitchroms) )  
# with(blocks, hist(maplen, breaks=100, xlim=c(0,20)) )
# get a good discretization
nbins <- 100
maxbinsize <- 1
lenbins <- quantile( blocks$maplen, probs=seq(0,1,length.out=nbins+1) )
lenbins <- c( lenbins[c(diff(lenbins)<maxbinsize,FALSE)], seq( lenbins[max(which(diff(lenbins)<maxbinsize))+1], max(lenbins), maxbinsize ) )
lenbins <- lenbins[-length(lenbins)]
lenbins[1] <- minblocklen
binsizes <- diff(c(lenbins,max(.chrlens)))
midbins <- lenbins[-1]-diff(lenbins)/2

lendist <- hist( blocks$maplen, breaks=c(lenbins,100), plot=FALSE )$counts

if (!file.exists("new-disc-trans.Rdata")) {
    gens <- cumsum( rep(1:36, each=10) )
    system.time( L <- error.disc.trans( lenbins=lenbins, gens=gens, chrlens=.chrlens[setdiff(seq_along(.chrlens),omitchroms)] ) )
    # user   system  elapsed 
    # 4947.929    1.728 5038.586 
    fp <- disc.fp(lenbins, chrlens=.chrlens)
    fp <- runmed(fp,k=15)
    save(L, fp, file="new-disc-trans.Rdata")
} else {
    load(file="new-disc-trans.Rdata")
    gens <- attr(L,"gens")
}

if (any(attr(L,"lenbins")!=lenbins) | any(attr(L,"gens")!=gens)) {
    stop("Mismatching discretizations, oops!")
}

# get initial initial estimate
S <- as.vector( smooth( hist( blocks$maplen, breaks=c(lenbins,100), plot=FALSE )$counts ) )
S <- S/mean(S)
# tail of S becomes 0
S[lenbins>20] <- exp( predict( lm( log(S) ~ lenbins + I(lenbins^2), subset=(lenbins>20) & (S>0) ), newdata=list(lenbins=lenbins[lenbins>20]) ) )
est.Sigma <- function ( np, X, alpha=100 ) {
    # estimate Sigma by combining the population average with the smoothed, observed data
    # where Sigma is the covariance matrix of X/np
    smX <- as.vector( smooth(X) )
    ppp <- (alpha*np/(2543640+alpha*np))  # 2543640 = total # pairs
    Sigma <- ppp*smX + (1-ppp)*sum(smX)*S/sum(S)
    return( Sigma/np^2 )
}
Ksvd <- svd( 1/sqrt(as.vector(S)) * L )
total.npairs <- choose( length(unique(blocks$id1)), 2 )
linverse <- linv( X=lendist/total.npairs, S=as.vector(S), Sigma=est.Sigma(np=total.npairs,X=lendist), Ksvd=Ksvd, maxk=30, fp=fp )
sminverse <- smoothed.nonneg( lS=linverse, Ksvd=Ksvd, bestk=15, lam=1 )
ans <- sinv( lendist, sminverse, L, fp, npairs=total.npairs, lam=0, gamma=1e6 )

toplot.fn <- function (plotthese) { lapply( plotthese, function (subc) {
                subind1 <- with( subset(indivinfo,COUNTRY_SELF%in%subc[[1]]), unlist( tapply( SUBJID, COUNTRY_SELF, function (x) tsample(x,1000) ) ) )
                subind2 <- with( subset(indivinfo,COUNTRY_SELF%in%subc[[2]]), unlist( tapply( SUBJID, COUNTRY_SELF, function (x) tsample(x,1000) ) ) )
                npairs <- length(subind1) * length(subind2) - choose( length(intersect(subind1,subind2)), 2 ) - length(intersect(subind1,subind2))
                lendist <- with( subset(blocks, (id1%in%subind1 & id2%in%subind2) | (id1%in%subind2 & id2%in%subind1) ), 
                    hist( maplen, breaks=c(lenbins,100), plot=FALSE )$counts )
                initans <- sinv( lendist, xinit=ans$par+1e-7, L, fp, npairs=npairs, minlen=minlen, gamma=10, lam=10 )
                list(
                        pointy = initans,
                        smooth = squoosh( lendist, xinit=initans, L, fp, npairs=initans$npairs, minlen=minlen, gamma=100, lam=0, fitfn="loglik", relthresh=2/initans$npairs, twid=.05 )
                    )
            } )
    }
lentimes <- function ( x, genbins, cumulative=TRUE ) {
    # return matrix of cumulative proportion of blocks at each gen coming from each bin
    props <- cbind( fp, sapply( seq_along(genbins)[-1], function (ell) L %*% ifelse(gens>genbins[ell-1]&gens<=genbins[ell],x,0) ) ) / diff(c(lenbins,100))
    # props <- sweep( props, 1, rowSums(props), "/" ) 
    if (cumulative) props <- rbind( 1e-8, apply(props,1,cumsum) )
    return(props)
}
plot.lentimes <- function ( ans, genbins, legend=TRUE ) {
    pcols <- rainbow_hcl(length(genbins),start=0,end=.75*360)
    with( ans,  {
        plot( lenbins, (lendist/npairs)/diff(c(lenbins,100)), type='n', xlab="", ylab="IBD rate", xlim=c(minlen,10), ylim=c(0,1) )
        lt <- lentimes( par, genbins )
        for (k in 1:(nrow(lt)-1)) { polygon( c(lenbins,rev(lenbins)), c(lt[k,],rev(lt[k+1,])), col=pcols[k] ) }
        lines( lenbins, (lendist/npairs)/diff(c(lenbins,100)), lwd=1 )
        if (legend) legend("topright", legend=c("false pos",paste((30/2)*genbins[-length(genbins)],(30/2)*genbins[-1],sep="-")), title="Years ago",fill=pcols, cex=8/10)
    } )
}



vikings <- c("Denmark","Norway","Sweden","Finland","Latvia")
balkans <- c("Yugoslavia", "Bulgaria", "Romania","Croatia","Bosnia","Montenegro","Macedonia", "Serbia", "Slovenia")
all.plotthese <- list( 
    UK= list( 
            "UK-Ireland"=list("United Kingdom","Ireland"),
            "UK-France"=list("United Kingdom","France"),
            "UK-Italy"=list("United Kingdom","Italy"),
            "UK-Turkey/Cyprus"=list("United Kingdom",c("Turkey","Cyprus"))
        ),
    Balkans = list( 
            "Balkans"=list(balkans,balkans),
            "Balkans-Albanian"=list(balkans,c("Albania","Kosovo")),
            "Balkans-Italy"=list(balkans,"Italy"),
            "Balkans-France"=list(balkans,"France")
        ),
    Albanian = list(
            "Albanian"=list(c("Albania","Kosovo"),c("Albania","Kosovo")),
            "Alb-Balk"=list(c("Albania","Kosovo"),balkans),
            "Alb-IT"=list(c("Albania","Kosovo"),"Italy"),
            "Alb-EL"=list(c("Albania","Kosovo"),"Greece"),
            "Alb-PL"=list(c("Albania","Kosovo"),"Poland")
        )
    )
all.toplot <- lapply( all.plotthese, toplot.fn )

pdf(file="inversion-distributions.pdf", width=7, height=7, pointsize=10)
layout( matrix( 1:10,nrow=5), heights=c(1.3,1,1,1.3,2) )
whichone <- 0
for ( y in c("Balkans","UK") ) {
    toplot <- all.toplot[[y]]
    cutgens <- c(0,36,100,169,289,529)
    genbins <- lapply(seq_along(cutgens)[-1],function(k)cutgens[c(k-1,k)]+c(1,0))
    pcols <- rainbow_hcl(length(cutgens),start=0,end=.75*360)
    opar <- par(mar=c(.25,3,3.25,.25)+.1,mgp=c(2,1,0))
    for ( x in names(toplot) ) {
        do.xlab.above <- x==names(toplot)[1] 
        do.xlab.below <- x==names(toplot)[length(toplot)] 
        do.ylab <- x==names(toplot)[2] 
        do.legends <- x==names(toplot)[1]
        if ( do.xlab.below ) { par(mar=par("mar")+c(3,0,0,0)) }
        whichone <- whichone+1
        with( toplot[[x]], {
                # Plot various smoothings
                plot( gens*30/2, coal.to.anc(pointy$par,gens), xlim=c(0,4000), ylim=c(0,8), type='n', ylab="", xlab="", xaxt="n" )
                polygon( c(gens,rev(gens))*30/2, c(coal.to.anc(pointy$par,gens),rep(0,length(gens))), col=adjustcolor("black",.4) )
                polygon( c(gens,rev(gens))*30/2, c(coal.to.anc(smooth$par,gens),rep(0,length(gens))), col=adjustcolor("red",.4) )
                abline(v=cutgens[-1]*30/2,col=adjustcolor("black",.25))
                textlab( "topleft", paste( LETTERS[whichone], ") ", x, sep=""), nudgey=.1 )
                # if (do.legends) legend("topright", title="penalization", legend=c("none","roughness"), fill=adjustcolor(c("black","red"),.4) )
                if (do.xlab.above) { genticks <- pretty((1/2)*gens[gens*30/2<4000]);axis(3,at=genticks*30,labels=genticks); mtext("generations ago", side=3, line=2, cex=8/12 ) }
                if (do.xlab.below) { axis(1); mtext("years ago", side=1, line=2, cex=8/12 ) }
                if (do.ylab) { mtext("# genetic ancestors", side=2, line=2, cex=8/12) }
        } )
        if ( do.xlab.above ) { par(mar=par("mar")+c(0,0,-3,0)) }
    }
    par(opar)
    # layout(cbind(matrix(1:4,nrow=2,byrow=TRUE),c(5,5)), heights=c(rep(1,length(toplot)-1),1.4), widths=c(1,1,.5) )
    opar <- par(mar=c(3,7.2,.25,3)+.1,mgp=c(2,1,0))
    for ( x in names(toplot[ifelse(y=="UK",2,1)]) ) {
        whichone <- whichone+1
        do.xlab <- TRUE
        do.ylab <- TRUE
        do.legends <- y=="UK"
        with( toplot[[x]][["smooth"]], {
                lt <- lentimes( par, cutgens )
                plot( lenbins, (lendist/npairs)/diff(c(lenbins,100)), type='n', ylab="", xlab="", xlim=c(minlen,6), xaxt="n", yaxt='n'  )
                for (k in 1:(nrow(lt)-1)) { polygon( c(lenbins,rev(lenbins)), c(lt[k,],rev(lt[k+1,])), col=pcols[k], ) }
                lines( lenbins, (lendist/npairs)/diff(c(lenbins,100)), lwd=1, col='red' )
                if (do.legends) {
                    legend("topright", legend=c("false pos",sapply(genbins,function(z)paste(z*30/2,collapse="-"))), title="years ago", fill=pcols, cex=1.2)
                }
                labind <- which.min(abs(lenbins-par("usr")[1]))
                lablocs <- lt[,labind][-1] - diff(lt[,labind])/2
                uselocs <- diff(lt[,labind])>.01
                axis(2, at=c(lablocs[uselocs],par("usr")[4]), labels=c(c("false pos",sapply(genbins,function(z)paste(z,collapse="-")))[uselocs],"generations ago"), las=2 )
                if (do.xlab) { axis(1); mtext("block length (cM)", side=1, line=2, cex=10/12 ) }
                if (do.ylab) { axis(4); mtext("IBD rate", side=4, line=2, cex=10/12) }
                textlab("topleft", paste(LETTERS[whichone], ") ", x, sep=""), cex=1.2, nudgex=.25, nudgey=.1 )
            } )
    }
}
dev.off()

pdf(file="inversion-distributions-coalrate.pdf", width=7, height=7, pointsize=10)
layout( matrix( 1:10,nrow=5), heights=c(1.3,1,1,1.3,2) )
whichone <- 0
for ( y in c("Balkans","UK") ) {
    toplot <- all.toplot[[y]]
    cutgens <- c(0,36,100,169,289,529)
    genbins <- lapply(seq_along(cutgens)[-1],function(k)cutgens[c(k-1,k)]+c(1,0))
    pcols <- rainbow_hcl(length(cutgens),start=0,end=.75*360)
    opar <- par(mar=c(.25,3,3.25,.25)+.1,mgp=c(2,1,0))
    for ( x in names(toplot) ) {
        do.xlab.above <- x==names(toplot)[1] 
        do.xlab.below <- x==names(toplot)[length(toplot)] 
        do.ylab <- x==names(toplot)[2] 
        do.legends <- x==names(toplot)[1]
        if ( do.xlab.below ) { par(mar=par("mar")+c(3,0,0,0)) }
        whichone <- whichone+1
        with( toplot[[x]], {
                # Plot various smoothings
                plot( gens*30/2, pointy$par, xlim=c(0,4000), ylim=c(0,1.9e-5), type='n', ylab="", xlab="", xaxt="n" )
                polygon( c(gens,rev(gens))*30/2, c((pointy$par),rep(0,length(gens))), col=adjustcolor("black",.4) )
                polygon( c(gens,rev(gens))*30/2, c((smooth$par),rep(0,length(gens))), col=adjustcolor("red",.4) )
                abline(v=cutgens[-1]*30/2,col=adjustcolor("black",.25))
                textlab( "topleft", paste( LETTERS[whichone], ") ", x, sep=""), nudgey=.1 )
                # if (do.legends) legend("topright", title="penalization", legend=c("none","roughness"), fill=adjustcolor(c("black","red"),.4) )
                if (do.xlab.above) { genticks <- pretty((1/2)*gens[gens*30/2<4000]);axis(3,at=genticks*30,labels=genticks); mtext("generations ago", side=3, line=2, cex=8/12 ) }
                if (do.xlab.below) { axis(1); mtext("years ago", side=1, line=2, cex=8/12 ) }
                if (do.ylab) { mtext("coalescent rate", side=2, line=2, cex=8/12) }
        } )
        if ( do.xlab.above ) { par(mar=par("mar")+c(0,0,-3,0)) }
    }
    par(opar)
    # layout(cbind(matrix(1:4,nrow=2,byrow=TRUE),c(5,5)), heights=c(rep(1,length(toplot)-1),1.4), widths=c(1,1,.5) )
    opar <- par(mar=c(3,7.2,.25,3)+.1,mgp=c(2,1,0))
    for ( x in names(toplot[ifelse(y=="UK",2,1)]) ) {
        whichone <- whichone+1
        do.xlab <- TRUE
        do.ylab <- TRUE
        do.legends <- y=="UK"
        with( toplot[[x]][["smooth"]], {
                lt <- lentimes( par, cutgens )
                plot( lenbins, (lendist/npairs)/diff(c(lenbins,100)), type='n', ylab="", xlab="", xlim=c(minlen,6), xaxt="n", yaxt='n'  )
                for (k in 1:(nrow(lt)-1)) { polygon( c(lenbins,rev(lenbins)), c(lt[k,],rev(lt[k+1,])), col=pcols[k], ) }
                lines( lenbins, (lendist/npairs)/diff(c(lenbins,100)), lwd=1, col='red' )
                if (do.legends) {
                    legend("topright", legend=c("false pos",sapply(genbins,function(z)paste(z*30/2,collapse="-"))), title="years ago", fill=pcols, cex=1.2)
                }
                labind <- which.min(abs(lenbins-par("usr")[1]))
                lablocs <- lt[,labind][-1] - diff(lt[,labind])/2
                uselocs <- diff(lt[,labind])>.01
                axis(2, at=c(lablocs[uselocs],par("usr")[4]), labels=c(c("false pos",sapply(genbins,function(z)paste(z,collapse="-")))[uselocs],"generations ago"), las=2 )
                if (do.xlab) { axis(1); mtext("block length (cM)", side=1, line=2, cex=10/12 ) }
                if (do.ylab) { axis(4); mtext("IBD rate", side=4, line=2, cex=10/12) }
                textlab("topleft", paste(LETTERS[whichone], ") ", x, sep=""), cex=1.2, nudgex=.25, nudgey=.1 )
            } )
    }
}
dev.off()

######
# boxplots

if (!file.exists("bplot-inversions.Rdata")) {
    thesecs <- list( "AL"=c("Albania","Kosovo"), "S-C"=c("Bosnia","Croatia","Serbia","Montenegro","Yugoslavia"), "R-B"=c("Romania","Bulgaria"), PL="Poland", HU="Hungary", UK=c("United Kingdom","England","Scotland","Wales"), IE="Ireland", DE="Germany", CHd="Swiss German", CHf="Swiss French", FR="France", IT="Italy", Iber=c("Spain","Portugal"), Bel=c("Belgium", "Netherlands"), Bal=c("Latvia", "Finland", "Sweden", "Norway", "Denmark") )
    theselatlongs <- list( 
            lat=sapply( thesecs, function (x) mean( subset(indivinfo,COUNTRY_SELF%in%x)$lat ) ),
            long=sapply( thesecs, function (x) mean( subset(indivinfo,COUNTRY_SELF%in%x)$long ) )
        )
    genbins <- list(c(0,36), c(37,100), c(101,169), c(170,289) )
    bplot <- list()
    for (x in names(thesecs)) {
        bplot[[x]] <- list()
        for (y in names(thesecs)) {
            if (y %in% names(bplot) & y!=x) {
                bplot[[x]][[y]] <- bplot[[y]][[x]]
            } else {
                lendist <- with( subset(blocks, (country1%in%thesecs[[x]] & country2%in%thesecs[[y]]) | (country1%in%thesecs[[y]] & country2%in%thesecs[[x]]) ), hist( maplen, breaks=c(lenbins,100), plot=FALSE )$counts )
                npairs <- if (x==y) { choose( sum(indivinfo$COUNTRY_SELF %in% thesecs[[x]]), 2 ) } else { sum(indivinfo$COUNTRY_SELF %in% thesecs[[x]]) * sum(indivinfo$COUNTRY_SELF %in% thesecs[[y]]) }
                bplot[[x]][[y]] <- get.coal.bounds(lendist,npairs,genbins=genbins, xinit=sminverse+1e-7, L=L, fp=fp, minlen=2, fitfn="loglik" )
            }
        }
    }
    save( thesecs, genbins, bplot, theselatlongs, file="bplot-inversions.Rdata" )
} else {
    load("bplot-inversions.Rdata")
}

get.bplotnums <- function (usethese,andthese=usethese) {
    genbins <- bplot[[1]][[1]]$genbins
    bplotnums <- do.call( rbind, lapply( seq_along(genbins), function (k) {
            x <- expand.grid( country1=usethese, country2=andthese )
            x$lower <- as.vector( sapply( bplot[usethese], function (z1) sapply(z1[andthese],function (z2) sum( ( diff(c(0,gens)) * coal.to.anc(z2$lower.ans[[k]]$par,gens) )[ gens>genbins[[k]][1] & gens<=genbins[[k]][2] ] ) ) ) )
            x$pointy <- as.vector( sapply( bplot[usethese], function (z1) sapply(z1[andthese],function (z2) sum( ( diff(c(0,gens)) * coal.to.anc(z2$ans$par,gens) )[ gens>genbins[[k]][1] & gens<=genbins[[k]][2] ] ) ) ) )
            x$smooth <- as.vector( sapply( bplot[usethese], function (z1) sapply(z1[andthese],function (z2) sum( ( diff(c(0,gens)) * coal.to.anc(z2$sq.ans$par,gens) )[ gens>genbins[[k]][1] & gens<=genbins[[k]][2] ] ) ) ) )
            x$upper <- as.vector( sapply( bplot[usethese], function (z1) sapply(z1[andthese],function (z2) sum( ( diff(c(0,gens)) * coal.to.anc(z2$upper.ans[[k]]$par,gens) )[ gens>genbins[[k]][1] & gens<=genbins[[k]][2] ] ) ) ) )
            x$lower.coal <- as.vector( sapply( bplot[usethese], function (z1) sapply(z1[andthese],function (z2) mean( ( z2$lower.ans[[k]]$par )[ gens>genbins[[k]][1] & gens<=genbins[[k]][2] ] ) ) ) )
            x$pointy.coal <- as.vector( sapply( bplot[usethese], function (z1) sapply(z1[andthese],function (z2) mean( ( z2$ans$par )[ gens>genbins[[k]][1] & gens<=genbins[[k]][2] ] ) ) ) )
            x$smooth.coal <- as.vector( sapply( bplot[usethese], function (z1) sapply(z1[andthese],function (z2) mean( ( z2$sq.ans$par )[ gens>genbins[[k]][1] & gens<=genbins[[k]][2] ] ) ) ) )
            x$upper.coal <- as.vector( sapply( bplot[usethese], function (z1) sapply(z1[andthese],function (z2) mean( ( z2$upper.ans[[k]]$par )[ gens>genbins[[k]][1] & gens<=genbins[[k]][2] ] ) ) ) )
            x$xval <- with(x, as.numeric(country1) + as.numeric(country2)/(1.2*nlevels(country2)) )
            x$ya <- paste(genbins[[k]]*30/2,collapse="-")
            return(x)
        } ) )
    bplotnums$ya <- factor(bplotnums$ya,levels=sapply(genbins,function(x)paste(x*30/2,collapse="-")),ordered=TRUE)
    return(bplotnums)
}

## boxplots
##  grouped by date range

for (all.boxplots in c(TRUE,FALSE)) {

if (all.boxplots)  {
    # ALL
    usethese <- names(thesecs)
    bplotnums <- get.bplotnums(usethese)
    pdf(file="inversion-boxplots-long.pdf",width=20,height=10,pointsize=10)
    ylimlist <- with( subset(bplotnums, pointy>1e-6 & smooth>1e-6), list( tapply( pmax(pointy,smooth), ya, function (x) c(0,max(x)) ), tapply( pmax(pointy,smooth), ya, function (x) c(0,quantile(x,.85)) ) ) )
    layout( (1:(2*nlevels(bplotnums$ya))), heights=rep(c(rep(1,nlevels(bplotnums$ya)-1),1.4),2) )

} else {
    # Short version for paper
    usethese <- c("S-C","PL","R-B","DE","UK","IT","Iber")
    bplotnums <- get.bplotnums(usethese)
    pdf(file="inversion-boxplots.pdf",width=7,height=5,pointsize=10)
    ylimlist <- with( subset(bplotnums, pointy>1e-6 & smooth>1e-6), list( tapply( pmax(pointy,smooth), ya, function (x) c(0,max(x)) ) ) )
    layout( (1:nlevels(bplotnums$ya)), heights=c(rep(1,nlevels(bplotnums$ya)-1),1.4) )
}

for (ylims in ylimlist) {
    opar <- par(mar=c(0,4,1,2)+.1)
    for (k in 1:nlevels(bplotnums$ya)) {
        do.xlab <- (k==nlevels(bplotnums$ya))
        if ( do.xlab ) {
            par(mar=c(4,4,1,2)+.1)
        }
        with( subset(bplotnums,ya==levels(ya)[k]), {
                cols <- countrycols[sapply(thesecs[levels(country1)],function(x)x[1])]
                names(cols) <- levels(country1)
                plot( 0, 0, type='n', xlim=c(1,nlevels(country1)+1), ylim=ylims[[k]], xlab="", xaxt='n', ylab="# genetic ancestors", )
                segments( x0=xval, y0=pmin(lower,pointy,smooth), y1=pmax(upper,pointy,smooth), col="grey" )
                rect( xleft=xval-1/(3*nlevels(country2)), xright=xval+1/(3*nlevels(country2)), ybottom=pointy, ytop=smooth, col=cols[as.numeric(country2)], border=cols[as.numeric(country2)] )
                # text( x=xval, y=0, labels=ifelse(lower>1e-7,"*","") )
                abline( v=(2:nlevels(country2)), lty=2 )
                textlab( "topright",  labels=paste(levels(ya)[k],"years ago"), cex=1.2, nudgey=.1, nudgex=.03 )
                textlab( "topright",  labels=paste("(",paste(floor(genbins[[k]]/2),collapse="-")," generations ago)",sep=""), cex=1.2, nudgey=.23, nudgex=.03 )
                if ( do.xlab ) {
                    colaxis( 1, at=xval, labels=(country2), las=3, col=cols[levels(country2)[country2]] )
                    # omgp <- par(mgp=c(3,3,0))
                    # axis( 1, at=(1:nlevels(country1))+0.5, tick=FALSE, labels=levels(country1) )
                    colaxis( 1, at=(1:nlevels(country1))+0.5, tick=FALSE, labels=levels(country1), line=3, col=cols[levels(country1)] )
                    # par(omgp)
                }
                if (k==1) {
                    colaxis( 3, at=(1:nlevels(country1))+0.5, tick=FALSE, labels=levels(country1), line=0, col=cols[levels(country1)], cex=1.5 )
                }
            } )
    }
    par(opar)
}
dev.off()
}


# In coalescent units
usethese <- names(thesecs)
bplotnums <- get.bplotnums(usethese)
pdf(file="inversion-boxplots-long-coal.pdf",width=20,height=10,pointsize=10)
ylimlist <- with( subset(bplotnums, pointy.coal>1e-8 & smooth.coal>1e-8), list( tapply( pmax(pointy.coal,smooth.coal), ya, function (x) c(0,max(x)) ), tapply( pmax(pointy.coal,smooth.coal), ya, function (x) c(0,quantile(x,.85)) ) ) )
layout( (1:(2*nlevels(bplotnums$ya))), heights=rep(c(rep(1,nlevels(bplotnums$ya)-1),1.4),2) )
opar <- par(mar=c(0,4,1,2)+.1)
for (ylims in ylimlist) {
    opar <- par(mar=c(0,4,1,2)+.1)
    for (k in 1:nlevels(bplotnums$ya)) {
        do.xlab <- (k==nlevels(bplotnums$ya))
        if ( do.xlab ) {
            par(mar=c(4,4,1,2)+.1)
        }
        with( subset(bplotnums,ya==levels(ya)[k]), {
                cols <- countrycols[sapply(thesecs[levels(country1)],function(x)x[1])]
                names(cols) <- levels(country1)
                plot( 0, 0, type='n', xlim=c(1,nlevels(country1)+1), ylim=ylims[[k]], xlab="", xaxt='n', ylab="coalescence rate" )
                segments( x0=xval, y0=pmin(lower,pointy.coal,smooth.coal), y1=pmax(upper.coal,pointy.coal,smooth.coal), col="grey" )
                rect( xleft=xval-1/(3*nlevels(country2)), xright=xval+1/(3*nlevels(country2)), ybottom=pointy.coal, ytop=smooth.coal, col=cols[as.numeric(country2)], border=cols[as.numeric(country2)] )
                # text( x=xval, y=0, labels=ifelse(lower>1e-7,"*","") )
                abline( v=(2:nlevels(country2)), lty=2 )
                textlab( "topright",  labels=paste(levels(ya)[k],"years ago"), cex=1.2, nudgey=.1, nudgex=.03 )
                textlab( "topright",  labels=paste("(",paste(floor(genbins[[k]]/2),collapse="-")," generations ago)",sep=""), cex=1.2, nudgey=.23, nudgex=.03 )
                if ( do.xlab ) {
                    colaxis( 1, at=xval, labels=(country2), las=3, col=cols[levels(country2)[country2]] )
                    # omgp <- par(mgp=c(3,3,0))
                    # axis( 1, at=(1:nlevels(country1))+0.5, tick=FALSE, labels=levels(country1) )
                    colaxis( 1, at=(1:nlevels(country1))+0.5, tick=FALSE, labels=levels(country1), line=3, col=cols[levels(country1)] )
                    # par(omgp)
                }
                if (k==1) {
                    colaxis( 3, at=(1:nlevels(country1))+0.5, tick=FALSE, labels=levels(country1), line=0, col=cols[levels(country1)], cex=1.5 )
                }
            } )
    }
    par(opar)
}
dev.off()


# look at these more closely
plot.all.inversions <- function (bplot) {
    # layout(matrix(1:9,nrow=3))
    # layout( cbind( rbind( matrix(1,nrow=4,ncol=5), c(2:6) ), c(7:11)  ) )
    layout( rbind( rep(1,5), c(2:6), c(7:11) ), heights=c(3,1,1) )
    par(mar=c(3,4,0,0)+.1)
    for ( y in names(bplot) ) {
        for (x in names(bplot[[y]]) ) {
            thisone <- with(bplot[[y]][[x]], c( list(pointy=ans,smooth=sq.ans), lower.ans, upper.ans ) )
            names(thisone)[3:length(thisone)] <- with(bplot[[y]][[x]], c( paste( sapply(genbins,paste,collapse="-"), "lower" ), paste( sapply(genbins,paste,collapse="-"), "upper" ) ) )
            plot( 0, 0, type='n', xlim=c(0,4000), ylim=c(0,2.5e-5), xlab="ya", ylab="coal rate" )
            for (k in seq_along(thisone)) { 
                lines( gens*30/2, thisone[[k]]$par, col=rainbow(10)[k] ) 
                polygon( c(gens,rev(gens))*30/2, c(thisone[[k]]$par,rep(0,length(gens))), col=adjustcolor(rainbow(10)[k],.4) )
            }
            legend("topright",title=paste(x,y),legend=names(thisone),fill=adjustcolor(rainbow(10),.4)[1:length(thisone)])
            for (k in seq_along(thisone)) { 
                plot.lentimes( thisone[[k]], c(0,sapply(bplot[[y]][[x]]$genbins,function(x)x[2]),1000) )
                textlab("topleft",names(thisone)[k])
            }
        }
    }
}

pdf(file="boxplotted-inversions.pdf",width=10,height=7,pointsize=10)
plot.all.inversions(bplot)
dev.off()

# pdf(file="selected-inversions.pdf",width=8,height=10,pointsize=15)
# layout( matrix(1:6,nrow=3) )
# par(mar=c(4,4,0,0)+.1)
# for (xy in list( c("AL","AL"), c("AL","PL"), c("AL","IT"), c("IE","IE"), c("IE","UK"), c("IE","AL") ) ) {
#     x <- xy[1]; y <- xy[2]
#     thisone <- rev( with(bplot[[y]][[x]], c( list(pointy=ans,smooth=sq.ans), lower.ans, upper.ans ) ) )
#     names(thisone)[3:length(thisone)] <- with(bplot[[y]][[x]], c( paste( sapply(genbins,paste,collapse="-"), "lower" ), paste( sapply(genbins,paste,collapse="-"), "upper" ) ) )
#     plot( 0, 0, type='n', xlim=c(0,4000), ylim=c(0,2.5e-5), xlab='generations ago', ylab='# genetic common anc.', cex.axis=1, cex.lab=1 )
#     abline(v=500*(1:7),col="grey")
#     for (k in seq_along(thisone)) { 
#         lines( gens*30/2, thisone[[k]]$par, col=rainbow(10)[k] ) 
#         polygon( c(gens,rev(gens))*30/2, c(thisone[[k]]$par,rep(0,length(gens))), col=adjustcolor(rainbow(10)[k],.4) )
#     }
#     textlab("topright",paste(xy,collapse="-"),cex=2)
# }
# dev.off()

pdf(file="example-inversion-bounds.pdf", width=6.5, height=8, pointsize=10)
thisone <- with( bplot[['PL']][['DE']], c( list(pointy=ans,smooth=sq.ans), lower.ans, upper.ans) )
genbins <- bplot[['PL']][['DE']]$genbins
names(thisone)[3:length(thisone)] <- with(bplot[[y]][[x]], c( paste( sapply(genbins,paste,collapse="-"), "lower" ), paste( sapply(genbins,paste,collapse="-"), "upper" ) ) )
thisone <- thisone[c(1,2,2+rep(1:length(genbins),each=2)+rep(c(0,length(genbins)),length(genbins)))]
pcols <- adjustcolor(c("black","red",rep(rainbow(length(genbins),start=.3,end=.7),each=2)),.5)
layout(1:length(thisone),heights=c(1.8,rep(1,length(thisone)-2),1.8))
for (k in seq_along(thisone)) { 
    par(mar=c(ifelse(k==length(thisone),4,0),4,ifelse(k==1,5,0),1)+.1)
    whichbin <- (k - 1 )%/% 2
    plot( 0, 0, type='n', xlim=c(0,4000), ylim=c(0,2.7e-5), xlab="", ylab="coal rate", xaxt='n', yaxt='n' )
    axis(1, labels=(k==length(thisone)) )
    axis(2, labels=(k%%2==1))
    if (k==1) {
        title(main="Germany-Poland consistent histories",line=4)
        genticks <- pretty((1/2)*gens[gens*30/2<4000])
        axis(3,at=genticks*30,labels=genticks); mtext("generations ago", side=3, line=2, cex=8/12 )
    }
    if (k==length(thisone)) mtext("years ago", side=1, line=3, cex=par("cex.axis")*par("cex") )
    if (whichbin) polygon( rep(genbins[[whichbin]],each=2)*30/2, c(0,1e-4,1e-4,0), col=grey(.8) )
    lines( gens*30/2, thisone[[k]]$par, col=rainbow(10)[k] ) 
    polygon( c(gens,rev(gens))*30/2, c(thisone[[k]]$par,rep(0,length(gens))), col=pcols[k] )
    textlab("topright",gsub("pointy","MLE",names(thisone)[k]),font=2,nudgey=.1)
}
dev.off()

pdf(file="tinyinversions.pdf", width=18, height=18, pointsize=8)
layout( rbind( c(0,1:length(bplot)), cbind(length(bplot)+(1:length(bplot)), 2*length(bplot)+matrix(1:(length(bplot)^2),nrow=length(bplot)) ) ) )
par(mar=c(0,0,0,0))
for (k in 1:2) { for (y in names(bplot)) {
    plot(0,0,type='n',xaxt='n',yaxt='n',xlim=c(0,2),ylim=c(0,2))
    text(1,1,y,cex=4)
} }
for (y in names(bplot)) {
    for (x in names(bplot[[y]])) {
        thisone <- with(bplot[[y]][[x]], list(pointy=ans,smooth=sq.ans) )
        # plot lower/upper bounds also:
        # thisone <- rev( with(bplot[[y]][[x]], c( list(pointy=ans,smooth=sq.ans), lower.ans, upper.ans ) ) )
        # names(thisone)[3:length(thisone)] <- with(bplot[[y]][[x]], c( paste( sapply(genbins,paste,collapse="-"), "lower" ), paste( sapply(genbins,paste,collapse="-"), "upper" ) ) )
        plot( 0, 0, type='n', xlim=c(0,4000), ylim=c(0,2.5e-5), xaxt='n', yaxt='n' )
        abline(v=500*(1:7),col="grey")
        for (k in seq_along(thisone)) { 
            lines( gens*30/2, thisone[[k]]$par, col=rainbow(10)[k] ) 
            polygon( c(gens,rev(gens))*30/2, c(thisone[[k]]$par,rep(0,length(gens))), col=adjustcolor(c("black","red",rainbow(10))[k],.4) )
        }
    }
}
dev.off()

