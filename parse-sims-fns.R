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


coalprob <- function ( gens, opts, nesize=opts$nesize, migprob=opts$migprob, ploidy=2 ) {
    # compute coalescent probs from nesize and migprob, recursively
    # note that gens is in *generations* =  2*meioses
    # 
    if (is.null(dim(migprob(1)))) {
        migprobfn <- function (t) { x <- migprob(t); dim(x) <- c(1,length(x)); return(x) }
    } else {
        migprobfn <- migprob
    }
    if (is.null(colnames(migprobfn(1)))) {
        states <- 1:ncol(migprobfn(1))
    } else {
        states <- colnames(migprobfn(1))
    }
    coal <- array(0,dim=c(max(gens),length(states),length(states)))
    dimnames(coal) <- list( NULL, states, states )
    # this indexes where coalescences can occur
    samestates <- c( as.vector( diag(length(states))>0 ), FALSE )
    # transition matrix for pairs: last state is the graveyard (coalescence)
    transprobs <- cbind( rbind( diag(length(states)) %x% diag(length(states)), 0 ), 0 )
    for (t in 1:max(gens)) {
        if (!missing(opts) & t>opts$ngens) {
            # only ran the simulation so long
            coal[t,,] <- 0
        } else {
            # nesize works in *diploid* individuals, so actual coalescence is smaller by a factor of 1/ploidy
            C <- as.vector(1/(ploidy*nesize(t)))
            M <- cbind( rbind( migprobfn(t) %x% migprobfn(t), 0 ), 0 )
            M[ samestates, ] <- (1-C)*M[ samestates, ] 
            M[ samestates, ncol(M) ] <- C
            transprobs <- transprobs %*% M
            coal[t,,] <- transprobs[-nrow(M),ncol(M)]
        }
    }
    return(coal[match(gens,1:max(gens)),,,drop=FALSE])
}

mergerate <- function ( x, L, opts, minminlen=0.1, gaplen=5 ) {
    # approx for the proportion of falsely merged blocks
    thesegens <- 1:(max(attr(L,"gens"))/2)
    coalgens <- coalprob(thesegens,opts)
    f <- function (t,x) { ( gaplen*t/100*(1+x*t/100) - 2*(1+x*t/100+(1/2)*(x*t/100)^2) )*exp(-x*t/100) }
    sum( sqrt(2)*sum(.chrlens-x)*(2*thesegens/100)*( f(2*thesegens,x)-f(2*thesegens,2*minminlen) ) * coalgens )
}

predict.blocks <- function ( L, opts, times=FALSE ) {
    # L gives number of blocks per constant rate of coalescence in windows of numbers of *meioses*;
    #  coalprob works in 2*meioses (generations) here; so half of these are zero;
    #  here we put in those zeros.
    coalgens <- coalprob(1:(max(attr(L,"gens"))/2),opts)
    coal <- numeric( 2*prod(dim(coalgens)) )
    coal[ 2*(1:prod(dim(coalgens))) ] <- coalgens
    dim(coal) <- c( 2*dim(coalgens)[1], dim(coalgens)[-1] )
    coal <- apply( coal, c(2,3), function (x) diff(c(0,cumsum(x)[attr(L,"gens")]))/diff(c(0,attr(L,"gens"))) )
    sampsize <- opts$sampsize
    npairs <- outer(sampsize,sampsize,"*")
    diag(npairs) <- choose(sampsize,2)
    coal <- coal * rep(npairs, each=dim(coal)[1])
    dim(coal) <- c( dim(coal)[1], length(coal)/dim(coal)[1] )
    coal <- coal[,upper.tri(npairs,diag=TRUE),drop=FALSE]
    if (times) {
         blocklens <- L * rep(coal,each=nrow(L))
         dim(blocklens) <- c(dim(L),dim(coal)[2])
         dimnames(blocklens) <- list( NULL, NULL, paste( rownames(npairs)[row(npairs)], colnames(npairs)[col(npairs)], sep="-" )[upper.tri(npairs,diag=TRUE)] )
    } else {
        blocklens <- L%*%coal
        colnames(blocklens) <- paste( rownames(npairs)[row(npairs)], colnames(npairs)[col(npairs)], sep="-" )[upper.tri(npairs,diag=TRUE)] 
    }
    return( (blocklens) )
}

# make stuff easier functions

loadblocks <- function ( opts, ploidy=2 ) {
    blocks <- read.table(opts$ibdfile,header=TRUE,sep=" ")
    # blocks$id1.orig <- blocks$id1
    # blocks$id2.orig <- blocks$id2
    blocks$id1 <- blocks$id1 %/% ploidy
    blocks$id2 <- blocks$id2 %/% ploidy
    blocks$chrom <- findInterval( blocks$start, .chrstarts/100 )
    blocks$mapstart <- blocks$start*100 - .chrstarts[blocks$chrom]  # these sims have recorded in cumulative distance
    blocks$mapend <- blocks$end*100 - .chrstarts[blocks$chrom]
    blocks$maplen <- with(blocks, mapend-mapstart)
    blocks$country1 <- names(opts$sampsizes)[ findInterval( blocks$id1, cumsum( c(0, unlist( opts$sampsizes ) ) ) ) ]
    blocks$country2 <- names(opts$sampsizes)[ findInterval( blocks$id2, cumsum( c(0, unlist( opts$sampsizes ) ) ) ) ]
    blocks$countrypair <- with( blocks, ifelse(country1<country2,paste(country1,country2,sep="-"),paste(country2,country1,sep="-")) )
    return(blocks)
}

init.sinv <- function( lendist, lenbins, npairs, L, fp ) {
    S <- as.vector(smooth(lendist))
    S <- pmax(1e-4, S/mean(S) )
    # tail of S becomes 0
    est.Sigma <- function ( np, X, alpha=100 ) {
        # estimate Sigma by combining the population average with the smoothed, observed data
        # where Sigma is the covariance matrix of X/np
        smX <- as.vector( smooth(X) )
        ppp <- (alpha*np/(npairs+alpha*np))  # 2543640 = total # pairs
        Sigma <- ppp*smX + (1-ppp)*sum(smX)*S/sum(S)
        return( Sigma/np^2 )
    }
    Ksvd <- svd( 1/sqrt(as.vector(S)) * L )
    linverse <- linv( X=lendist/npairs, S=as.vector(S), Sigma=est.Sigma(np=npairs,X=lendist), Ksvd=Ksvd, maxk=30, fp=fp )
    return( smoothed.nonneg( lS=linverse, Ksvd=Ksvd, bestk=15, lam=1 ) )
}


plot.ans <- function (anslist,opts,thispair,L,genscale=TRUE,coalrate=TRUE,dothese,main,legend1=TRUE,legend2=FALSE, tcols, plots=c("coal","spectrum"), spectrum.xlim, spectrum.ylim, ...) {
    # plot nice stuff about a named list of sinv objects
    if (missing(opts) && 'opts'%in%names(anslist)) { opts <- anslist[['opts']] }
    if ( missing(L) && 'L'%in%names(anslist) ) { L <- anslist[['L']] }
    if ('anslist'%in%names(anslist)) { anslist <- anslist[['anslist']] }
    if ( (!missing(thispair)) && (thispair %in% names(anslist))) { anslist <- anslist[[thispair]] }
    if (missing(dothese)) { dothese <- names(anslist)[sapply(anslist,class)=="sinv"] }
    if (missing(tcols)) { tcols <- rainbow_hcl(length(dothese)) }
    if (is.null(names(dothese))) { names(dothese) <- dothese }
    if (missing(main) && !missing(thispair)) {
        main <- thispair
    } else if (missing(main)) {
        main <- ""
    } else if (length(main)==1) { main <- rep(main,2) }
    if (genscale) { timescale <- 1/2 } else { timescale <- 29/2 }
    npairs <- outer(opts$sampsizes,opts$sampsizes,"*")
    diag(npairs) <- choose(opts$sampsizes,2)
    gens <- attr(L,"gens")
    lenbins <- attr(L,"lenbins")
    midbins <- c(lenbins,max(.chrlens))[-1]-diff(c(lenbins,max(.chrlens)))/2
    binsizes <- diff(c(lenbins,max(.chrlens)))
    coal <- coalprob(gens%/%2,opts) 
    cdn <- dimnames(coal)
    dim(coal) <- c(dim(coal)[1],prod(dim(coal)[-1]))
    dimnames(coal) <- list( NULL, outer(cdn[[2]], cdn[[3]], paste, sep="-") )
    coal <- coal[,upper.tri(npairs,diag=TRUE),drop=FALSE]
    # account for coal working in generations
    coal <- coal/2
    npairs <- npairs[upper.tri(npairs,diag=TRUE)]
    predicted <- predict.blocks(L,opts)
    if (missing(thispair)) {
        coal <- rowSums(coal*npairs[col(coal)]/sum(npairs))
        predicted <- rowSums(predicted*npairs[col(predicted)]/sum(npairs))
    } else {
        coal <- coal[,thispair]
        predicted <- predicted[,thispair]
    }
    coalvals <- data.frame(c( list(theoretical=coal), lapply(dothese, function (x) anslist[[match(x,names(anslist))]]$par) ))
    # plot coalescent rates or number of ancestors?
    if (coalrate) {
        coallab <- "coal rate"
    } else {
        coalvals <- data.frame( lapply( coalvals, coal.to.anc, gens ) )
        coallab <- "# shared ancestors"
    }
    spectra <- data.frame( c( list("theoretical"=predicted/binsizes, "observed"=anslist[[1]]$lendist/binsizes), 
                lapply(dothese, function (x) with(anslist[[x]], (npairs * (L%*%par))/binsizes ) )
            ) )
    ### BEGIN PLOTTING
    if ("coal" %in% plots) {
        plot( gens*timescale, rowMeans(coalvals), type='n', ylab=coallab, xlab=if(genscale){"generations ago"} else {"years ago"}, main=main[1], ... )
        for (k in seq_along(dothese)) {
            polygon( c(gens,rev(gens))*timescale, c(coalvals[[names(dothese)[k]]],rep(0,length(gens))), col=adjustcolor(tcols[k],.4) )
        }
        lines( gens*timescale, coalvals$theoretical, col='green', lwd=2 )  # generations
        if (is.numeric(opts$ngens)) { abline(v=opts$ngens*timescale*2, col='grey', lty=2) }
        if (legend1) { legend("topright", fill=c(NA,tcols[seq_along(dothese)]), border=c("green",rep("black",length(dothese))), legend=c("theoretical",names(dothese)) ) }
    }
    # predicted and observed blocks
    if ("spectrum" %in% plots) {
        if (missing(spectrum.xlim)) { spectrum.xlim <- range(midbins[is.finite(log(midbins))]) }
        if (missing(spectrum.ylim)) { spectrum.ylim <- range((spectra$theoretical/npairs)[is.finite(log(spectra$theoretical/npairs))]) }
        plot( midbins, spectra$theoretical/npairs, type='l', col='green', lwd=2, log='xy', xlab="block length (cM)", ylab="density", main=main[2], xlim=spectrum.xlim, ylim=spectrum.ylim )
        lines( midbins, spectra$observed/npairs )
        for (k in match(dothese,names(anslist))) {
            lines( midbins, spectra[[names(dothese)[k]]]/npairs, col=tcols[k] )
        }
        if (legend2) { legend("topright", lty=1, lwd=c(2,1,rep(1,length(anslist))), col=c("green","black",tcols[seq_along(dothese)]), legend=c("theoretical","observed",names(anslist))) }
    }
    return( invisible( list(gens=gens, coal=coalvals, midbins=midbins, spectra=spectra, npairs=npairs, lenbins=lenbins ) ) )
}

do.everything <- function (prefix,minblocklen,fprate=function(x)0,maxgen,genbins=NULL) {
    opts <- simopts[[prefix]]  # of note is nesize and ngens
    blocks <- loadblocks(opts)
    if (missing(minblocklen)) { minblocklen <- opts$minlen }
    if (missing(maxgen)) { maxgen <- opts$ngens }

    indivinfo <- data.frame( SUBJID=1:sum(opts$sampsizes)-1, COUNTRY_SELF=names(opts$sampsizes)[unlist(lapply(1:length(opts$sampsizes),function(k) rep(k,opts$sampsizes[k])))] )
    nsamples <- table(indivinfo$COUNTRY_SELF)
    npairs <- outer(nsamples,nsamples,"*")
    diag(npairs) <- choose(nsamples,2)
    names(npairs) <- paste( names(nsamples)[row(npairs)], names(nsamples)[col(npairs)], sep="-" )
    npairs <- npairs[upper.tri(npairs,diag=TRUE)]

    # get a good discretization
    nbins <- 100
    maxbinsize <- 1
    lenbins <- with( subset(blocks,maplen>minblocklen), quantile( maplen, probs=seq(0,1,length.out=nbins+1) ) )
    lenbins <- c( lenbins[c(diff(lenbins)<maxbinsize,FALSE)], seq( lenbins[max(which(diff(lenbins)<maxbinsize))+1], max(lenbins), maxbinsize ) )
    lenbins <- lenbins[-length(lenbins)]
    lenbins[1] <- minblocklen
    binsizes <- diff(c(lenbins,max(.chrlens)))
    midbins <- c( lenbins[-1]-diff(lenbins)/2, lenbins[length(lenbins)] )
    maxbin <- 5*ceiling(max(blocks$maplen)/5)

    lendist <- do.call( cbind, with( subset(blocks,maplen>minblocklen), tapply( maplen, countrypair, function (x) hist(x, breaks=c(lenbins,maxbin), plot=FALSE)$counts ) ) )

    # gens <- cumsum( rep(1:36, each=10) )
    gens <- 1:(2*maxgen)
    L <- theoretical.operator( lenbins=lenbins, gens=gens, chrlens=.chrlens )
    fp <- rep(0, length(lenbins))

    # add in unaccounted-for false positives
    predicted.blocklens <- predict.blocks( L, opts )
    oops <- rpois( length(lendist), fprate(midbins)*predicted.blocklens )
    true.lendist <- lendist
    lendist <- lendist + oops

    sminverse <- init.sinv(rowSums(lendist), lenbins, sum(npairs), L, fp)
    ans <- sinv( rowSums(lendist), xinit=sminverse, L, fp, npairs=sum(npairs), lam=0, gamma=1e6 )

    bothans <- lapply( colnames(lendist), function (x) {
            noconstraints <- sinv( lendist[,x], xinit=ans, L, fp, npairs=npairs[x], gamma=10, lam=10 )
            list(
                pointy = noconstraints,
                smooth = squoosh( lendist[,x], xinit=noconstraints, L, fp, npairs=npairs[x], gamma=100, lam=0, relthresh=2/npairs[x] ),
                smoother = squoosh( lendist[,x], xinit=noconstraints, L, fp, npairs=npairs[x], gamma=100, lam=0, relthresh=4/npairs[x] ),
                lower.bounds <- lapply( genbins, function (gb) squoosh( lendist[,x], npairs[x], xinit=noconstraints, L, fp, lam=0, gamma=0, del=1, delvec=( (gens>=gb[1]) & (gens<gb[2]) ) ) ),
                upper.bounds <- lapply( genbins, function (gb) squoosh( lendist[,x], npairs[x], xinit=noconstraints, L, fp, lam=0, gamma=0, del=1, delvec=( (gens<gb[1]) | (gens>=gb[2]) ) ) )
            ) }
        )
    names(bothans) <- colnames(lendist)
    return( list(opts=opts, anslist=bothans, L=L, fp=fp, npairs=npairs ) )
}




## automatically get stuff out of python dict text format

parsedict <- function(x) {
    # Parse the text string corresponding to a python dict
    x <- gsub("\\<None\\>", "NA", gsub(":", "=", strsplit( gsub( "}.*$", "", gsub( "^.*\\{", "", x ) ), "," )[[1]] ) )
    x <- paste( "list(", paste( x, collapse=", " ), ")" )
    x <- gsub("'([0-9.]{1,})'", "\\1", x )
    return( eval(parse(text=x)) )
}

getoptions <- function (prefix, logfilename ) {
    # read in command-line options from the log file
    logfilenames <- list.files( ".", paste("^", prefix, ".*\\.log$",sep="") )
    if (length(logfilenames)>1) {
        stop("Ambiguous log file prefix.")
    } else {
        logfilename <- logfilenames[1]
    }
    opts <- parsedict( system( paste( "grep '^options'", logfilename ), intern=TRUE ) )
    names(opts)[names(opts)=="sampsizes"] <- "opts.sampsizes"
    opts$sampsizes <- parsedict( system( paste( "grep '^sampsizes:'", logfilename ), intern=TRUE ) )
    return( opts )
}
