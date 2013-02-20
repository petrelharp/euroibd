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
require(MASS)      # for ??/
require(mgcv)      # for mroot
require(quadprog)  # for solve.QP


sinv <- function ( lendist, xinit, L, fp, npairs, mingen=min(attr(L,"gens")), maxgen=max(attr(L,"gens")), minlen=min(attr(L,"lenbins")), maxlen=max(attr(L,"lenbins")), lam=0, gamma=0, del=0, lb=1e-10, control=list(maxit=10000), retlik=FALSE, ... ) {
    # Invert by maximum likelihood plus a penalty sum( lam*x^2 ) + gamma*|Dx|^2 + sum( del*x )^2
    #  and constraining x>= lb and support(x) in [mingen,maxgen]
    # Note lam and del can be vectors, to constrain size only at certain areas.
    #
    # Strictly positive lower bound may be needed to avoid getting stuck at 0.
    gens <- attr(L,"gens")
    lenbins <- attr(L,"lenbins")
    if (class(xinit)=="sinv") {
        xinit <- xinit$par
    }
    if (length(lb)==1) { lb <- rep(lb,length(xinit)) }
    dgens <- diff(gens[gens>=mingen & gens<=maxgen])
    mu <- (lendist/npairs)[ lenbins>=minlen & lenbins<=maxlen ]
    fp <- fp[ lenbins>=minlen & lenbins<=maxlen ]
    L <- L[ lenbins>=minlen & lenbins<=maxlen, gens>=mingen & gens<=maxgen ]
    L1 <- colSums(L)
    # badness avoided by taking lb > 0?
    # fn <- function (x) { z <- as.vector(L%*%x+fp); L1%*%x - mu %*% ifelse(z>0,log(ifelse(z>0,z,1)),0) + (1/npairs)*( sum(lam*x^2) + gamma*sum((diff(x)/dgens)^2) + sum(del*x)^2 ) }
    fn <- function (x) { z <- as.vector(L%*%x+fp); L1%*%x - mu %*% ifelse(mu>0,log(z),0) + (1/npairs)*( sum(lam*x^2) + gamma*sum((diff(x)/dgens)^2) + sum(del*x)^2 ) }
    gr <- function (x) { z <- as.vector(L%*%x+fp); L1 - ifelse(z>0, mu/z, 0 )%*%L  + (1/npairs)*(2*lam*x + gamma*2*( c(0,diff(x)/dgens) - c(diff(x)/dgens,0) ) + 2*del*sum(del*x) ) }
    # hess <- function (x) { z <- L%*%x+fp; t(L) %*% diag(ifelse(z>0, mu / z^2, 0 )) %*% L + lam*2 }  # missing gamma term
    ans <- optim( par=xinit[gens>=mingen & gens<=maxgen], fn=fn, gr=gr, lower=lb[gens>=mingen & gens<=maxgen], method="L-BFGS-B", control=control )
    if (ans$convergence != 0) { warning("Convergence ", ans$convergence,": ",ans$message) }
    ans$par <- c( rep(0,sum(gens<mingen)), ans$par, rep(0,sum(gens>maxgen)) )
    ans <- c(ans, list( lendist=lendist, xinit=xinit, mingen=mingen, maxgen=maxgen, minlen=minlen, maxlen=maxlen, lam=lam, gamma=gamma, del=del, lb=lb, npairs=npairs ) )
    if (retlik) { ans$loglik <- fn }
    class(ans) <- "sinv"
    return(ans)
}

loglik <- function (xx,L,fp) {
    # get negative log likelihood of a sinv object
    # note that it omits constant factors (e.g. lfactorial(sum(lendist))/npairs)
    lenbins <- attr(L,"lenbins")
    gens <- attr(L,"gens")
    with( xx, {
        mu <- lendist/npairs
        x <- par
        uselens <- (lenbins>=minlen) & (lenbins<=maxlen) & (mu>0)
        usegens <-  (gens>=mingen) & (gens<=maxgen)
        L1 <- colSums(L[uselens,usegens])
        z <- as.vector(L[uselens,usegens]%*%x[usegens]+fp[uselens])
        return( L1%*%x[usegens] - mu[uselens] %*% log(z) )
    } )
}

liklims <- function (x,L,fp,probs=.95,fitfn=function(y)loglik(y,L,fp)) {
    # Estimate via bootstrapping reasonable limits for the likelihood
    # z <- replicate( 1000, {
    #         with(x, rpois( length(lenbins), npairs * (L%*%par + fp) ) )
    #     } )
    # logliks <- apply(z,2,function(zz) {
    #         x$lendist <- zz; fitfn(x) } )
    logliks <- replicate( 1000, {
            x$lendist <- with(x, rpois( length(lenbins), npairs * (L%*%par + fp) ) )
            fitfn(x)
        } )
    return( quantile(logliks,probs=probs) )
}

# Helper plotting fn
plot.sinv <- function (more.ans,ans,genbins=2*c(0,10,30,120,2000),title=NA) {
    if ( missing(ans) & "par"%in%names(more.ans) ) { ans <- more.ans; more.ans <- NULL }
    else if ( missing(ans) ) { ans <- more.ans[[1]] }
    layout(matrix(c(1,2,3,3),nrow=2))
    oe <- observed.error( ans, L=L, fp=fp, p=c(.25,.75), twid=.1 )
    # Estimated coal rate
    plot( 29*gens/2, ans$par, type='l', xlim=c(0,4000), ylim=range(ans$par,ans$xinit,finite=TRUE), col='blue', ylab="coal rate", lwd=2 )
    rug(29*genbins/2, lwd=3)
    lines( 29*gens/2, ans$xinit, col=adjustcolor('red',0.5), lty=1, lwd=1 )
    for (k in seq_along(more.ans)) {
        col <- rainbow(length(more.ans)+10)[k]
        lines( 29*gens/2, more.ans[[k]]$par, col=col, lty=1, lwd=2 )
        abline(v=29*c(more.ans[[k]]$mingen/2,more.ans[[k]]$maxgen/2),lwd=4,col=adjustcolor(col,.5))
    }
    # Residuals
    plot( lenbins, ( L%*% ans$par + fp - ans$lendist/ans$npairs )/apply(oe,1,diff), type='l', log='x', col=adjustcolor('black',.75), lwd=4, ylab="normalized resids" )
    abline(v=c(ans$minlen,ans$maxlen),lwd=4,col=adjustcolor("black",.5))
    abline(h=0)
    for (k in seq_along(more.ans)[-1]) {
        col <- rainbow(length(more.ans)+10)[k]
        lines( lenbins, ( L%*% more.ans[[k]]$par + fp - more.ans[[k]]$lendist/more.ans[[k]]$npairs )/apply(oe,1,diff), col=adjustcolor(col,.75), lwd=2 )
        abline(v=c(more.ans[[k]]$minlen,more.ans[[k]]$maxlen),lwd=4,col=adjustcolor(col,.5))
    }
    # Block len distrns
    plot( lenbins, (ans$lendist/ans$npairs)/diff(c(lenbins,100)), type='l', log='xy', col=adjustcolor('black',.75), lwd=4, ylab="length distrn" )
    abline(v=c(ans$minlen,ans$maxlen),lwd=4,col=adjustcolor("black",.5))
    lines( lenbins, ( L%*% ans$par + fp )/diff(c(lenbins,100)), col='blue')
    lines( lenbins, ( L%*% ans$xinit + fp )/diff(c(lenbins,100)), col='red')
    for (k in seq_along(more.ans)) {
        col <- rainbow(length(more.ans)+10)[k]
        lines( lenbins, ( L%*% more.ans[[k]]$par + fp )/diff(c(lenbins,100)), col=adjustcolor(col,.75), lwd=2, )
        abline(v=c(more.ans[[k]]$minlen,more.ans[[k]]$maxlen),lwd=4,col=adjustcolor(col,.5))
        for (ell in seq_along(genbins)[-1]) {
            lines( lenbins, ( L%*% ifelse(gens>genbins[ell-1]&gens<=genbins[ell],more.ans[[k]]$par,0) )/diff(c(lenbins,100)), col=adjustcolor(col,.75), lwd=2, lty=2+ell )
        }
    }
    polygon( x=c(lenbins,rev(lenbins)), y=pmax(1e-10,c(oe[,1]/diff(c(lenbins,100)),rev(oe[,2]/diff(c(lenbins,100))))), border=NA, col=adjustcolor("blue",.5) )
    legend("topright",lty=c(rep(1,3+length(names(more.ans))),2+seq_along(genbins)[-1]),lwd=2, col=c("black","blue","red",rainbow(length(more.ans)+10)[seq_along(names(more.ans))],rep("black",max(0,length(genbins)-1))), legend=c("truth","estimated","initial guess",names(more.ans),paste("--",genbins[-1]/2)), title=title )
    layout(1)
}

coal.to.anc <- function (x,gens) {
    # convert coalescent rate to numbers of genetic ancestors
    x*(sum(.chrlens)*gens+length(.chrlens))
}

observed.error <- function (ans, npairs=ans$npairs, L, fp, p=c(.05,.95), eps=1e-10, twid=0 ) {
    # return lower and upper bounds on the actual length distribution
    lenbins <- attr(L,"lenbins")
    mu <- pmax(eps, L%*%ans$par + fp )
    oe <- sapply( p, function (p) {
            x <- qpois( p, lambda=mu*npairs )/npairs
        } )
    oe <- oe * rep(c(1-twid,1+twid),each=dim(oe)[1])
    return(oe)
}

squoosh <- function ( lendist, npairs, xinit, L, fp, lam=1, gamma=1e2, del=0, lamvec, delvec, stepsize=2, maxsteps=50, fitfn=function(x)loglik(x,L,fp), thresh, relthresh=2/npairs, minlen=2, twid=.05, control=list(trace=0,maxit=10000), debug=FALSE, ... ) {
    # Penalize by lamvec * weight, with weight increasing in multiplicative increments of stepsize,
    #   until the maximum fitfn() is greater than thresh.
    # Move by factors of 10 and then by factors of stepsize.
    if (class(xinit)=="sinv") {
        ans <- list( par=xinit$par, npairs=npairs, lendist=lendist, minlen=xinit$minlen, maxlen=xinit$maxlen, mingen=xinit$mingen, maxgen=xinit$maxgen, convergence=0 )
        # we want to only use the estimated coalescent rates from this solution; not others.
        if (missing(minlen)) { minlen <- ans$minlen } else { ans$minlen <- minlen }
        if (missing(lamvec)) { lamvec <- rep(1,length(ans$par)) }
        if (missing(delvec)) { delvec <- rep(1,length(ans$par)) }
        ans$npairs <- npairs
        ans$lendist <- lendist
    } else {
        if (missing(lamvec)) { lamvec <- rep(1,length(xinit)) }
        if (missing(delvec)) { delvec <- rep(1,length(xinit)) }
        ans <- sinv( lendist, xinit, L, fp, npairs=npairs, lam=lam*lamvec, gamma=gamma, del=del*delvec, minlen=minlen, control=control, ... )
    }
    initpar <- ans$par
    if (is.character(fitfn) && fitfn=="loglik") { fitfn <- function (x) loglik(x,L,fp) }
    if (is.character(fitfn) && fitfn=="resids") {
        # alternative: fit by residuals.
        fitfn <- function (x) {
                oe <- observed.error( x, L=L, fp=fp, twid=twid, p=c(.25,.75), eps=2e-7 )
                stderrs <- ifelse( oe[,1]==oe[,2], Inf, apply(oe,1,diff) )
                resids <- ( L%*% x$par + fp - x$lendist/x$npairs )/stderrs
                return( sum(abs(resids[lenbins>x$minlen])) )
            }
    }
    # if (missing(thresh)) { thresh <- liklims(ans,L,fp,probs=1-relthresh,fitfn=fitfn) }
    if (missing(thresh)) { thresh <- fitfn(ans)+relthresh }
    if ( ans$convergence!=0 || fitfn(ans)>thresh ) { warning("Initial penalization exceeded threshold or didn't converge."); return(ans) }
    # Aim is to get convergent answer with fitfn below thresh, with penalizations as large as possible.
    # To do this, increase penalizations by step until the criteria fail;
    #  then undo the last increase;
    #  and decrease the step size and repeat, until the minimum resolution of stepsize is reached.
    step <- 16*stepsize
    naccepted <- 0
    for (k in 1:maxsteps) {
        cat(".")
        lam <- lam*step
        gamma <- gamma*step
        del <- del*step
        ansp <- sinv( lendist, initpar, L, fp, npairs=npairs, lam=lam*lamvec, gamma=gamma, del=del*delvec, minlen=minlen, control=control, ... )
        if ( ( ansp$convergence!=0 || fitfn(ansp)>thresh ) && step<=stepsize ) { 
            # we're at minimum resolution, and that step didn't work. stop and keep previous answer.
            break; 
        } else if( ansp$convergence!=0 || fitfn(ansp)>thresh ) { 
            # that step didn't work; decrease resolution and try again from previous answer
            gamma <- gamma/step; lam <- lam/step; del <- del/step; step <- step/2
        } else { 
            # ok, that step worked. keep it.
            naccepted <- naccepted + 1
            cat("+")
            ans <- ansp 
        }
    }
    cat("\n")
    ans$fitfn <- fitfn
    if ( naccepted == 0 | ans$convergence != 0 ) { warning("Did not finally achieve convergence."); if (debug) {browser()}; return(ansp) }
    return( ans )
}

coal.bounds <- function (countries, genbins=list( c(0,36), c(37,100), c(101,169)), ...) {
    # genbins <- list( c(0,36), c(37,100), c(101,169), c(170,289) )  # roughly 500,1500,2500, and 4000 years ago
    # genbins <- list( c(0,36), c(37,100), c(101,169) )  # roughly 500,1500,2500 years ago
    if (length(grep("-",countries))==0) { 
        cp <- gsub(" ",".",paste(sort(countries),collapse="-")) 
    } else {
        cp <- countries
    }
    lendist <- with( subset(blocks, countrypair==cp), hist( maplen, breaks=c(lenbins,100), plot=FALSE )$counts )
    # npairs <- ifelse(countries[1]==countries[2],choose(nsamples[countries[1]],2),nsamples[countries[1]]*nsamples[countries[2]])
    npairs <- poppairs$npairs[ match(cp,poppairs$countrypair) ]
    get.coal.bounds( lendist, npairs, genbins, ... )
}

get.coal.bounds <- function ( lendist, npairs, genbins, xinit, L, fp, ... ) {
    ans <- sinv( lendist, xinit=xinit, L=L, fp=fp, npairs=npairs, lam=0, gamma=10, del=0, ... )
    sq.ans <- squoosh( lendist, npairs, xinit=ans, L, fp, lam=1, gamma=100, ... )
    lower.ans <- lapply( genbins, function (gb) squoosh( lendist, npairs, xinit=ans, L, fp, lam=0, gamma=0, del=1, delvec=( (gens>=gb[1]) & (gens<gb[2]) ), ... ) )
    upper.ans <- lapply( genbins, function (gb) squoosh( lendist, npairs, xinit=ans, L, fp, lam=0, gamma=0, del=1, delvec=( (gens<gb[1]) | (gens>=gb[2]) ), ... ) )
    return( list( ans=ans, sq.ans=sq.ans, lower.ans=lower.ans, upper.ans=upper.ans, genbins=genbins ) )
}


########

# truncated exponential generation
rtexp <- function(n, rate, trunc=Inf) { -log(1-runif(n)*(1-exp(-abs(trunc)*rate)))/rate } 

fpfun <- function (maplen) {
    # False positive rate, per pair and per cM, estimated for blocks with score <1e-9
    exp(-13.704-2.095*maplen+4.381*sqrt(maplen))
}

disc.fp <- function (lenbins, chrlens=100, r=1/100) {
    # Returns discretized false positive rates per pair.
    f <- fpfun
    fprates <- numeric(length(lenbins))
    for (G in chrlens) {
        fprates <- fprates + G * sapply( seq_along(lenbins), function (k) integrate( f, lenbins[k], c(pmin(lenbins,G),G)[k+1] )$value )
    }
    return(fprates)
}

powerfun <- function (maplen) {
    # Power function estimated by logistic regression, for blocks with score <1e-9
    # z <- -7.2109449 - 0.5446777*maplen + 5.7685511*sqrt(maplen)
    # return( exp(z)/(1+exp(z)) )
    # power function estimated in fit-error-model.R
    warning("Using old power function?  Use error.density.")
    return( 1 - 1 / ( 1 + 0.05718153 * maplen * exp( 0.95613164 * maplen ) ) )
}

simulate.blocks <- function (t, npairs, mincutoff=1, chrlens=100, true.lengths=FALSE) {
    # True lengths:
    #  simulate from distribution proportional to (t/100)*((t/100)*(G-x)+1)exp(-xt/100)
    #  by rejection sampling 
    npairs <- 1000
    cprob <- 1-(1/(chrlens[1]+200/t))*(100/t)*(1-(1+chrlens[1]*t/100)*exp(-chrlens[1]*t/100)) # mean prob of being not rejected
    # denominator: \int_0^G (t/100)*((t/100)*(G-x)+2)exp(-xt/100) dx
    denom <- (chrlens[1]*t/100+2)*(1-exp(-chrlens[1]*t/100)) - (1-(chrlens[1]*t/100)*exp(-chrlens[1]*t/100)) 
    xx <- rtexp( floor(npairs* denom / cprob), rate=t/100, trunc=chrlens[1] )
    xx <- xx[ (rbinom(length(xx),1,1-xx/(chrlens[1]+200/t))==1) ] 
    if (true.lengths) {
        return(xx)
    }
    # Add errors:
    cens <- rbinom(length(xx),1,prob=prob.censor(xx))
    signs <- 1-2*rbinom(length(xx),1,prob=prob.down(xx))
    yy <- ifelse( cens, NA, xx + ifelse( signs<0, -rtexp(length(xx),rate=down.rate(xx),trunc=xx-mincutoff), pmax(0,mincutoff-xx)+rtexp(length(xx),rate=up.rate(x),trunc=chrlens[1]) ) )
    return(yy)
}

adjust.hist <- function (x) {
    # Adjust a histogram for fp rates and power
    # i.e. observed = true * power + fp
    # so tru = (observed-fp)/power
    for (z in c("counts","intensities","density")) {
        x[[z]] <- pmax(0, (x[[z]] - fpfun(x$mids))/powerfun(x$mids) )
    }
    return( x )
}

theoretical <- function (x0,t,x1,chrlen) {
    # mean number per pair in theory of length x
    # either density or integrated over an interval
    #    vectorized in x and t, not chrlen.
    # note factor of 4 is ploidy squared.
    if (missing(x1)) {
        # density per cM
        4 * ( (t/100)*((t/100)*pmax(0,chrlen-x0)+2*(x0<chrlen))*exp(-x0*t/100) )
    } else {
        # total numbers
        4 * ( (pmax(0,chrlen-x0)*t/100 + (x0<chrlen))*exp(-x0*t/100) - (pmax(0,chrlen-x1)*t/100 + (x1<chrlen))*exp(-pmin(x1,chrlen)*t/100) )
    }
}

theoretical.gaps <- function (x0,t,gaplen,minminlen,x1,chrlen) {
    # mean number of gaps no more than gaplen long due to single-generation correlations
    f <- function (x,t) { 
        4 * (pmax(0,x-2*minminlen)^2-pmax(0,x-2*minminlen-gaplen)^2)/100^2 * ( (t-2)*((t-2)*pmax(0,chrlen-x)/100+2*(x<chrlen))*exp(-x*t/100+4*minminlen/100) ) +
            (t-2) * (pmax(0,x-2*minminlen)^2-pmax(0,x-2*minminlen-gaplen)^2)/100^2 * ( (t-1)*((t-1)*pmax(0,chrlen-x)/100+2*(x<chrlen))*exp(-x*t/100+2*minminlen/100) ) 
    }
    if (missing(x1)) {
        # density per cM
        f(x0,t)
    } else {
        if (length(t)==1){ t <- rep(t,length(x0)) }
        sapply( seq_along(x0), function (k) integrate(f,lower=x0[k],upper=x1[k],t=t[k])$value )
    }
}

theoretical.operator <- function( lenbins, mingens, maxgens, gens=mingens:maxgens, chrlens=100, fn=theoretical, ... ) {
    # Work out the operator without error
    L <- matrix(0, nrow=length(lenbins), ncol=length(gens) )
    for (chrlen in chrlens) {
        L <- L + outer( 1:length(lenbins), gens, function (k,t) { fn( x0=lenbins[k], t=t, x1=c(lenbins[-1],Inf)[k], chrlen=chrlen, ... ) } )
    }
    # Adjust for bin length
    gendiff <- diff(gens)
    L <- L %*% diag( c(gendiff, floor(mean(rev(gendiff)[1:5]))) )
    attr(L,"lenbins") <- lenbins
    attr(L,"gens") <- gens
    return( L )
}

meanrate <- function (x, coalprobs, y) {
    # Return the mean rate of IBD of length at least x
    #  without error
    # if the coalescent rate is coalprobs
    coaldist <- - diff( cumprod( 1-c(0,coalprobs) ) )
    ans <- sum( coaldist * rowSums( sapply( .chrlens, function (chrlen) theoretical(x0=x,t=seq_along(coaldist),x1=y,chrlen=chrlen) ) ) )
    return( ans )
}

prob.censor <- function (x,mincutoff=1) {
    # fitted probability of being unobserved given true length of x
    1/(1+0.0772355*x^2*exp(0.5423082*x))
}
prob.down <- function (x, mincutoff=1) {
    # probability the observed block is shorter than the true block
    ( 1 - 1/(1+0.5066205*pmax(0,x-mincutoff)*exp(0.6761991*pmax(x-mincutoff))))*0.3419458
}

# for testing:
# prob.down <- function (x,mincutoff=1) { ifelse( x>mincutoff, 1, 0 ) }
# prob.censor <- function (x,mincutoff=1) { ifelse( x>mincutoff, 0, 1 ) }

up.rate <- function (x) {
    # parameter for (conditioned) exponential distr'n of observed-true length given true length of x if observed > true
    # 1.399283 + 1/(3487.896989*x)
    1.399283
}
down.rate <- function (x) {
    # parameter for (conditioned) exponential distr'n of observed-true length given true length of x if observed < true
    pmin( 12, (0.4009342+1/(0.1816122*x)) )
}

error.density <- function (x,y0,y1,mincutoff=1) {
    # For each block of length x, the mean number of blocks of length y observed
    # From fit-error-model.R
    if (missing(y1)) {
        ifelse( is.na(y0)|is.na(x), prob.censor(x,mincutoff),
            (1 - prob.censor(x,mincutoff)) * ifelse( y0<x, 
                prob.down(x,mincutoff)*down.rate(x)*exp(-down.rate(x)*(x-y0))/(1-exp(-down.rate(x)*(x-mincutoff))),  # down -- note x>y0>mincutoff
                (1-prob.down(x,mincutoff))*up.rate(x)*exp(-up.rate(x)*(y0-pmax(x,mincutoff)))  # up -- note y>max(x,mincutoff)
                ) )
    } else {
        if (any(is.na(y0) | is.na(y1))) {
            stop("Can't have both y0 and y1 with some limits NA.")
        } else {
            (1 - prob.censor(x)) * (
                prob.down(x,mincutoff)*abs( exp(-down.rate(x)*( abs(pmax(0,x-y0)) )) - exp(-down.rate(x)*( abs(pmax(0,x-y1)) )) )/(1-exp(-down.rate(x)*x))  #down 
                + (1-prob.down(x,mincutoff))*abs( exp(-up.rate(x)*( abs(pmin(0,pmax(x,mincutoff)-y0)) )) - exp(-up.rate(x)*( abs(pmin(0,pmax(x,mincutoff)-y1)) )) )   # up
                )
        }
    }
}

theoretical.error <- function (y0,t,y1,chrlen, ...) {
    # Combine theoretical and error density to get mean number of blocks observed between y0 and y1 from time t
    xscale <- log(10)/max(0.4698545,t/100)
    if (missing(y1)) {
        xlocation <- y0
        fscale <- max(0.0001, theoretical(x0=xlocation,t=t,chrlen=chrlen) * error.density(x=xlocation,y0=y0) )
        lower <- max( 0, xlocation-xscale )
        upper <- min( xlocation+xscale, chrlen )
        res <- integrate( f=function (x) { ( theoretical(x0=x*xscale+xlocation,t=t,chrlen=chrlen) * error.density(x*xscale+xlocation,y0=y0) )/fscale }, lower=(lower-(xlocation))/xscale, upper=(upper-(xlocation))/xscale, ... )
    } else {
        y1 <- min(y1, chrlen)
        xlocation <- (y0+y1)/2
        fscale <- max(0.0001, theoretical(x0=xlocation,t=t,chrlen=chrlen) * max( error.density(x=xlocation,y0=y0), error.density(x=xlocation,y0=y1) ) )
        lower <- max( 0, xlocation-xscale )
        upper <- min( xlocation+xscale, chrlen )
        xscale <- max( xscale, y1-y0 )
        res <- integrate( f=function (x) { ( theoretical(x0=x*xscale+xlocation,t=t,chrlen=chrlen) * error.density(x=x*xscale+xlocation,y0=y0,y1=y1) )/fscale }, lower=(lower-(xlocation))/xscale, upper=(upper-(xlocation))/xscale, ... )
        # tryCatch( res <- integrate( f=function (x) { ( theoretical(x0=x*xscale+xlocation,t=t,chrlen=chrlen) * error.density(x=x*xscale+xlocation,y0=y0,y1=y1) )/fscale }, lower=(lower-(xlocation))/xscale, upper=(upper-(xlocation))/xscale, ... ), error=function (e) recover(), finally=cat(".") )
        # integrate( f=function (x) { ( theoretical(x0=x,t=t,chrlen=chrlen) * error.density(x=x,y0=y0,y1=y1) ) }, lower=lower, upper=upper )
    }
    res$value <- res$value*fscale*xscale
    res$abs.error <- res$abs.error*fscale*xscale
    return(res)
}

error.disc.trans <- function (lenbins, mingens, maxgens, gens=mingens:maxgens, chrlens=100) {
    # Return lenbins x mingens:maxgens matrix from theoretical.error
    #   giving the mean number of blocks falling in each interval of lenbins
    #   that come from a constant unit of coalescence rate across an interval of gens (in meioses)
    L <- matrix(0, nrow=length(lenbins), ncol=length(gens) )
    for (chrlen in chrlens) {
        for (k in seq_along(lenbins)) {
            for (ell in seq_along(gens)) {
                L[ k, ell ] <- L[ k, ell ] + theoretical.error( y0=lenbins[k], t=gens[ell], y1=c(lenbins[-1],Inf)[k], chrlen=chrlen )$value
            }
        }
    }
    # Adjust for bin length
    gendiff <- diff(gens)
    L <- L %*% diag( c(gendiff, floor(mean(rev(gendiff)[1:5]))) )
    attr(L,"lenbins") <- lenbins
    attr(L,"gens") <- gens
    return( L )
}

disc.trans <- function (lenbins, mingens=1, empirical.gen=0, maxgens=100, minbin=1, chrlens=100, nbins=500, r=1/100, lenfns, integrals=TRUE) {
    # Returns the discretized operator
    # chrlens is a vector of chromosome lengths
    # lenbins is vector of lower enpoints of bins of block length
    # r is recombination rate.
    if (missing(lenbins)) { 
        lenbins <- exp( seq(log(minbin),log(max(chrlens)),length.out=nbins+1) )[-(nbins+1)] 
    } else {
        nbins <- length(lenbins)
    }
    L <- matrix( 0, nrow=length(lenbins), ncol=maxgens-mingens+1 )
    if (empirical.gen>=mingens) {
        #  Want lenfns[[chrom]][[t]](x) to give the mean number of blocks of size x 
        f <- function (x,t,chrom) { lenfns[[chrom]][[t]](x) * powerfun(x) }
        L[,(mingens:empirical.gen)-mingens+1] <- integrate.disc.trans( lenbins, mingens, empirical.gen, chrlens, f, integrals=integrals )
    }
    if (empirical.gen<maxgens) {
        f <- function (x,t,chrom) { r*t*(r*t*(.chrlens[[chrom]]-x)+2)*exp(-r*t*x) * powerfun(x) }
        L[,(max(mingens,(empirical.gen+1)):maxgens)-mingens+1] <- integrate.disc.trans( lenbins, max(mingens,empirical.gen+1), maxgens, chrlens, f, integrals=integrals )
    }
    # Factor of 2!
    L <- 2*L
    attr(L,"lenbins") <- lenbins
    return( L )
}

integrate.disc.trans <- function (lenbins, mingens, maxgens, chrlens, f, integrals=TRUE ) {
    # Do the integration over bins
    #   f should be a function from (x,t,G) to mean number of length x blocks after t meioses on chromosome of length G
    #   sum over chromosomes since this is mean numbers
    L <- matrix( 0, nrow=length(lenbins), ncol=maxgens-mingens+1 )
    for (chrom in 1:length(chrlens)) {
        # f <- function (x,t) { (G-x)*lenfn(t,x) * powerfun(x) + G*fpfun(x) }
        G <- chrlens[chrom]
        nzeros <- sum(lenbins>=G)
        binlens <- diff(c(lenbins[lenbins<G],G))
        midbins <- lenbins[lenbins<G] + binlens/2
        if (integrals) {
            L <- L + sapply( mingens:maxgens, function (t) {
                    sapply( seq_along(lenbins), function (k) {
                            val <- if (lenbins[k]<G) {
                                integrate( f, lenbins[k], c(pmin(lenbins,G),G)[k+1], t=t, chrom=chrom )$value
                            } else { 0 }
                            # include the atom for the entire chromosome
                            if (lenbins[k]<=G & c(lenbins,Inf)[k+1]>G) { val <- val + f(G) }
                            return(val)
                        } )
                } )
        } else {
            L <- L + sapply( mingens:maxgens, function (t) {
                    c( f(x=midbins,t=t,chrom=chrom)*binlens, rep(0,nzeros) )
                } )
        }
    }
    return(L)
}


svd.project <- function ( x, v, k=NCOL(v), nonneg=FALSE, normalize=FALSE, px=t(v[,1:max(k),drop=FALSE])%*%x, lowerbound=0 ) {
    # Project x into the span of the first k columns of v
    #  which are assumed to be orthonormal.
    #  k can be a vector
    # px gives the coefficients.
    #
    # If nonneg, then require that the projection is nonnegative
    # i.e. find coefficients b1, ..., bn
    #   such that sum( (x-b*v)^2 ) is minimized
    #   with the constraint that b*v >= 0
    # ...so if d[j] = x^T v[,j]
    #    and if D[i,j] = v[,i]^T v[,j]
    #    then we minimize (1/2) b^T D b - d^T b
    #    subject to the constraint that v b >= 0 .
    maxk <- max(k)
    if (length(px)<maxk) { px <- c(px,rep(0,maxk-length(px))) }
    if ( "v" %in% names(v) ) { v <- v$v }  # allow svd() objects
    v <- v[,1:maxk,drop=FALSE]
    if (any( !all.equal(t(v)%*%v, diag(maxk)) )) { warning("Matrix v is not orthogonal.") }
    if (nonneg) {
        qpk <- lapply( k, function (kk) {
                qp <- solve.QP( diag(kk), px[1:kk], t(v[,1:kk]), bvec=rep(lowerbound,nrow(v)) )[c("solution", "value", "unconstrained.solution", "iterations", "iact")]
            } )
        ppx <- sapply( qpk, function (qp) c( qp$solution, rep(0,maxk-length(qp$solution)) ) )
    } else {
        ppx <- sapply( k, function (j) c( px[1:j], rep(0,maxk-j) ) )
    }
    ans <- v %*% ppx
    if (normalize) {
        # ensure results have same L1 norm as x
        # ... note if the projection is not close to x, this may be a BAD idea.
        normf <- apply(abs(ans), 2, sum)/sum(abs(x))
        ans <- sweep( ans, 2, normf, "/" )
        pps <- sweep( ppx, 2, normf, "/" )
    }
    attr(ans, "ppx") <- ppx
    return( ans )
}

nearest.nonneg <- function ( ppx, v, k=length(ppx), maxk=ncol(v), x=v[,1:k,drop=FALSE]%*%ppx[1:k] ) {
    # Find the smallest-norm element in the column space of v
    #  having projection ppx in the first k coordinates
    # Here x is the function to be projected;
    #  if both x and ppx are present then x should be in the span of the first k columns of v.
    if ( "v" %in% names(v) ) { v <- v$v }  # allow svd() objects
    if (missing(ppx)) {
        x <- svd.project( x, v, k=k )
        ppx <- attr(x,'ppx')
    }
    ppx <- c( ppx[1:k], rep(0, maxk-k) )
    ppx[(k+1):maxk] <- solve.QP( diag(maxk-k), dvec=rep(0,maxk-k), Amat=t(v[,(k+1):maxk]), bvec=-x )$solution 
    ans <- v[,1:maxk] %*% ppx
    attr(ans, "ppx") <- ppx
    return( ans )
}

avg.nonneg <- function ( lS, Ksvd=lS$Ksvd, maxk=lS$maxk, bestk=lS$noisek, nsamp=200 ) {
    # Sample from noise about U,
    #  project into nonnegativeness,
    #  and take the average.
    sqCovU <- mroot(lS$covU[1:bestk,1:bestk])
    Z <- rnorm(nsamp*bestk)
    dim(Z) <- c(bestk,nsamp)
    Z <- (lS$U[1:bestk] + apply(Z, 2, function (z) (sqCovU%*%z)) )/Ksvd$d[1:bestk]
    PZ <- apply( Z, 2, function (z) svd.project( v=Ksvd$v, k=maxk, px=z, nonneg=TRUE ) )
    return( apply(PZ,1,mean) )
    # return( PZ )
}

smoothed.nonneg <- function (U, lS, Ksvd=lS$Ksvd, maxk=ncol(Ksvd$v), bestk=lS$noisek, lam, lowerbound=0 ) {
    # Solve in V-space,
    # minimizing |U'-U|^2 + lam |Df|^2
    if (missing(U)) { U <- lS$U }
    if (length(U)<maxk) { U <- c( U, rep(0, maxk-length(U) ) ) }
    deriv <- apply( Ksvd$v[,1:maxk], 2, diff )
    covUinv <- ginv(lS$covU)
    pmat <- outer( Ksvd$d[1:maxk], Ksvd$d[1:maxk] ) * covUinv
    # pmat <- diag( Ksvd$d[1:maxk] )  # transforms V's to U's
    if (missing(lam)) { lam <- sum( pmat^2 )/sum( deriv^2 ) }
    # normalize a bit to help the algorithm out
    normf <- pmat[1] + lam*sum(deriv[,1]^2)
    # measure closeness in U, so with (b = smoothed V) and (c = observed/estimated V)
    # solve.QP solves: min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0
    # note: 1/2 (b-c)^T D (b-c) = 1/2 b^T D b - c^T D b + (const)
    ppx <- solve.QP( Dmat=(pmat + lam*t(deriv)%*%deriv)/normf, dvec=(Ksvd$d[1:maxk]*covUinv%*%U[1:maxk])/normf, Amat=t(Ksvd$v[,1:maxk]), bvec=rep(lowerbound,nrow(Ksvd$v)) )$solution
    ans <- Ksvd$v[,1:maxk] %*% ppx
    attr(ans, "ppx") <- ppx   # coordinates in V-space
    attr(ans, "weights") <- c( sum((pmat%*%ppx-U[1:maxk])^2), lam*sum((deriv%*%ppx)^2) )  # weight on fitting versus smoothing
    attr(ans, "lam") <- lam
    return( ans )
}

linv <- function (X, S=rep(1,length(X)), Sigma=S, maxk, K, Ksvd, fp) {
    # Compute things relevant to inversion with covariance matrix Sigma
    #   for discretization K
    #   and applied to data X
    # "Find f", where
    #   E[X] = Kf+g  and cov(X) = Sigma
    # with g the (assumed known) false positive rate
    # and solve this via
    #   Y = sqrt(1/S) X
    #   cov(Y) = sqrt(1/S) Sigma sqrt(1/S)
    # (where S is chosen to be close to Sigma)
    # and Ksvd gives the decomposition
    #   sqrt(1/S) K = u diag(d) v'
    # so
    #   v diag(1/d) u' E[Y] = f
    # and let
    #   U = u' Y
    #   E[U] = diag(d) v' f
    #   cov(U) = u' sqrt(1/S) Sigma sqrt(1/S) u
    # Returns:
    # svd = (u,d,v) where sqrt(1/S) K = u diag(d) v'
    #   U = u' sqrt(1/S) X
    # approxes[k] = sum_{i=1}^k v[,i] (U[i]/d[i])
    # Note: should never use dimensions with d < sqrt(.Machine$double.eps)
    #  The sqrt arises because we want, for instance,
    #    ginv(K) %*% K %*% x = x     ... maybe???
    if (is.vector(S) || min(dim(S))==1) { 
        sqSinv <- diag( 1/sqrt(as.vector(S)) )
        S <- diag(as.vector(S)) 
    } else {
        sqSinv <- ginv( mroot(S) )
    }
    if (is.vector(Sigma) || min(dim(Sigma))==1) { Sigma <- diag(as.vector(Sigma)) }
    # may want to use a common svd for many linv's,
    # and if it is large may not want to return it every time
    if (missing(Ksvd)) { 
        Ksvd <- svd( sqSinv %*% K ) 
        return.Ksvd <- TRUE 
    } else {
        return.Ksvd <- FALSE
    }
    if (missing(maxk)) { maxk <- min( length(Ksvd$d), 2*max( which( Ksvd$d > sqrt(.Machine$double.eps) ) ) ) }
    U <- with(Ksvd, t(u) %*% sqSinv %*% (X-fp) )
    covU <- with(Ksvd, t(u) %*% sqSinv %*% Sigma %*% t(sqSinv) %*% u)
    # guess at where the approximation is breaking down
    noisek <- min( maxk, max( which( abs(U/sqrt(diag(covU)))>5 ) ) )
    # Subsequent approximations
    approxes <- with(Ksvd, sapply(1:maxk, function (k) v[,1:k,drop=FALSE] %*% (U[1:k]/d[1:k]) ) )
    linvobj <-  list( U=U, covU=covU, approxes=approxes, maxk=maxk, noisek=noisek )
    if (return.Ksvd) { linvobj <- c( linvobj, list( Ksvd=Ksvd ) ) }
    class( linvobj ) <- "linv"
    return(linvobj)
}

plot.linv <- function (lS, maxk=lS$maxk, noisek=lS$noisek, ylims=range(lS$approxes[,1:noisek]), cols=rainbow(maxk+5), ...) {
    # Plot the successive approximate inverses
    par(mfrow=c(2,1))
    plot( 1, 1, xlim=c(1,nrow(lS$approxes)), ylim=ylims, type="n", xlab="time", ylab="inverse" )
    plotk <- min(noisek+4,maxk)
    legend("topright", legend=1:plotk, lty=1, col=cols[1:plotk])
    for (k in 1:plotk) {
        lines( lS$approxes[,k], col=cols[k], ... )
    }
    with(lS, plot( U/sqrt(diag(covU)), col=c("black","red")[1+(seq_along(U)<=noisek)] ) )
    abline(h=c(-5,5))
    return( invisible( lS$approxes[,1:plotk] ) )
}

print.linv <- function (lS) { 
    with(lS, { print( paste("First", noisek, "terms:") ); print(U[1:noisek]) } )
}

predict.linv <- function (lS, x, projx, bestk, Ksvd=lS$Ksvd, tol=.1, nonneg=TRUE, ...) {
    #  Returns an estimate of x'f = Z
    #  by truncating the sum
    #   x'f = sum_k (v'x)_k (v'f)_k
    #  at the maximal k such that
    #    sqrt( var(Z) ) / |Z| < tol
    # where w is of the form (1,1,...,1,0,0,...,0)
    #         Z = (w v'x)' diag(1/d) U
    #    var(Z) = (w v'x/d)' cov(U) (w v'x/d)
    # Here lS should be a 'linv' object.
    #   ... need to pass in Ksvd if not carried with lS.
    # Constrain results to be nonnegative if nonneg is TRUE.
    if ( is.null(Ksvd) ) { stop("Missing the svd decomposition.") }
    if (missing(x)) { x <- NULL }
    if (missing(projx)) {
        projx <- project( x, Ksvd$v, 1:lS$maxk, nonneg=nonneg )
    }
    ppx <- attr(projx, "ppx")  # the coefficients of the projections
    Z <- t(ppx/Ksvd$d[1:lS$maxk]) %*% lS$U[1:lS$maxk]
    sdZ <- sqrt( diag( t(ppx/Ksvd$d[1:lS$maxk]) %*% lS$covU[1:lS$maxk,1:lS$maxk] %*% (ppx/Ksvd$d[1:lS$maxk]) ) )
    if (missing(bestk)) { 
        if (any(abs(sdZ/Z)<tol)) { 
            bestk <- max(which(abs(sdZ/Z)<tol)) 
        } else {
            bestk <- 1
        }
    }
    out <- list( ans=Z[bestk], Z=Z, sdZ=sdZ, x=x, ppx=ppx, bestk=bestk, maxk=lS$maxk, Ksvd=lS$Ksvd ) 
    class(out) <- "lapprox"
    return(out)
}

plot.lapprox <- function ( lA, Ksvd=lA$Ksvd, maxk=lA$maxk, bestk=lA$bestk, newplot=TRUE, ... ) {
    # visualize the approximation to x'f
    if (newplot & !is.null(Ksvd)) { opar <- par(mfrow=c(1,2)) }
    plot( lA$Z[1:maxk], ylim=range( lA$Z[1:min(bestk,maxk)] ), col=c("black","red")[1+(1:maxk<=bestk)], xlab="number of terms", ylab="projection" )
    plotses( 1:maxk, lA$Z[1:maxk], lA$sdZ[1:maxk] )
    abline(h=0, lty=2)
    if (!is.null(Ksvd)) {
        xapx <- Ksvd$v[,1:lA$maxk] %*% lA$ppx[,1:lA$bestk]
        plot( xapx[,1], type="n", ylim=range(c(lA$x,xapx)), xlab="", ylab="value", ... )
        if (!is.null(lA$x)) { lines( lA$x ) }
        for (k in 1:ncol(xapx)) { lines(xapx[,k], col=rainbow(ncol(xapx)+5 )[k]) }
        legend("topright", legend=1:ncol(xapx), col=rainbow(ncol(xapx)+5)[1:ncol(xapx)], lty=1)
    }
    if (newplot & !is.null(Ksvd)) { par(opar) }
    return( invisible( list( y=lA$Z[1:maxk], xapx=xapx ) ) )
}

print.lapprox <- function (lA) { 
    print(paste("Approximation using", lA$bestk, "terms:", lA$ans))
    print(paste(" standard error:", lA$sdZ[lA$bestk]))
}

