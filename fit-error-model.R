# Written 2013 by Peter Ralph and Graham Coop
# 
# contact: petrel.harp@gmail.com
#
#     To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty. 
# 
# You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 
# 
#
source("ibd-blocks-fns.R")

load("true-pos-everything.Rdata")

# produce plots?
produce.plots <- TRUE
if (produce.plots) {
    pdf(file="error-fit-plots.pdf", width=7, height=7, pointsize=10)
    layout((1:3))
}

# Estimate the error kernel, that gives, for a block of length x, the probability we observe it at length y.
mincutoff <- 1

# create truncated vectors
tblocks$xx <- tblocks$maplen
tblocks$yy <- ifelse(tblocks$obstotal>mincutoff, tblocks$obstotal, NA)

if (FALSE) {
    tblocks$yy <- tblocks$obstotal
    resid.lm <- with(subset(tblocks,yy>.5), nls( xx/yy ~ 1/(1+aa/yy) , start=list(aa=.1) ) )
    shrink.yy <- with(subset(tblocks, yy>mincutoff), loess( xx ~ I(yy/(1+coef(resid.lm)/yy)), trace.hat="approx", family="symm", span=.4) )
    ytest <- seq(0.5,15,length.out=100)
    with( tblocks, plot( xx ~ I(yy / (1 + 0.1297/yy)) ) )
    abline(0,1,col='red')
    lines( ytest, predict(shrink.yy,newdata=ytest), col='blue', lwd=2 )
    #     aa     bb     cc 
    # 0.0000 0.7610 0.2001 
    tblocks$yy <- with(tblocks, yy * predict( resid.lm, newdata=list(yy=yy) ) )
    ytest <- seq(0.5,15,length.out=100)
    shrink.yy <- with(subset(tblocks, yy>mincutoff), loess( xx/yy ~ yy, trace.hat="approx", family="symm", span=.7) )
    # tblocks$yy <- with(tblocks, yy * 1-(1)*(1-xx/predict(shrink.yy,xx)) )
    plot( xx/yy ~ yy, data=tblocks, subset=yy>0.1, ylim=c(0,2) )
    abline(h=1,col='red')
    lines( ytest, (predict( resid.lm, newdata=list(yy=ytest) )), col='red', lwd=2 )
    lines( ytest, predict(shrink.yy,newdata=ytest), col='blue', lwd=2 )
    plot( yy*predict(resid.lm,newdata=list(yy=yy)) ~ xx, data=tblocks )
    abline(0,1, col='red')
}

lenbins <- c(0, exp( seq(log(mincutoff),log(ceiling(max(tblocks$xx[!is.na(tblocks$yy)]))+.1),length.out=10) ) )
midbins <- lenbins[-1] - diff(lenbins)/2
obsbins <- c(0, exp( seq(log(mincutoff),log(ceiling(max(tblocks$yy[!is.na(tblocks$yy)]))+.1),length.out=10) ) )
obsmids <- obsbins[-1] - diff(obsbins)/2
tblocks$tbins <- with( tblocks, cut( maplen,breaks=lenbins,right=FALSE) )

# Fit proportions of y>x, c<y<x, and y==c
fit.cens <- with(tblocks, optim( par=c(.5,.5,.5,.5,.5), f=function(r) { 
        if (r[1]<0 | r[1]>1 | r[3]>1) { Inf } else { 
            # pneg <- r[3]
            # pcens <- (1-r[1])/(1+r[1]*(exp(r[2]*xx)-1))
            pneg <- ( 1 - 1/(1+r[3]*pmax(0,xx-mincutoff)*exp(r[4]*pmax(xx-mincutoff))))*r[5]
            pcens <- 1/(1+r[1]*xx^2*exp(r[2]*xx))
            mean( -log( ifelse( is.na(yy), pcens, ifelse( yy<xx, (1-pcens)*pneg, (1-pcens)*(1-pneg) ) ) ) )
        } } ) )
prob.censor <- function (x,mincutoff=1) {
    # fitted probability of being unobserved given true length of x
    1/(1+fit.cens$par[1]*x^2*exp(fit.cens$par[2]*x))
}
prob.down <- function (x, mincutoff=1) {
    # probability the observed block is shorter than the true block
    ( 1 - 1/(1+fit.cens$par[3]*pmax(0,x-mincutoff)*exp(fit.cens$par[4]*pmax(x-mincutoff))))*fit.cens$par[5]
}
hist.xx <- with( tblocks, hist( xx, breaks=100, plot=FALSE ) )
prop.cens <- with( tblocks, tapply( factor(as.numeric(yy<xx),levels=c(0,1)), cut(xx,breaks=hist.xx$breaks), function (y) { table(y,useNA="always")/length(y) } ) )
prop.cens <- do.call( rbind, lapply( prop.cens, function (z) if (length(z)==0) { c(NA,NA,NA) } else { z } ) )
with(hist.xx, plot( mids, prop.cens[,1]/(1-prop.cens[,3]), xlab="block length (cM)", ylab="probability", main=expression(paste("fitted forms for ",c(x)," and ",gamma(x))) ) )  # down
with(hist.xx, points( mids, prop.cens[,2]/(1-prop.cens[,3]), col='red' )  ) # up
with(hist.xx, points( mids, prop.cens[,3], col='green' ) )
legend("topright", legend=c("censored","y<x","y>x"), pch=1, col=c('green','red','black') )
with(hist.xx, lines( mids, prob.down(mids), col='red' ) )
with(hist.xx, lines( mids, 1-prob.down(mids), col='black' ) )
with(hist.xx, lines( mids, prob.censor(mids), col='green' ) )
fit.cens$par
# [1] 0.0772355 0.5423082 0.5066205 0.6761991 0.3419458

# Fit just those with y>x
fit.up <- with(tblocks, optim( par=c(.3,.05), f=function(r) { 
        rate.up <- r[1] + 1/(r[2]*xx)
        if (any(rate.up<0,na.rm=TRUE)) { Inf } else { 
            mean( -log(rate.up[yy>xx]) + (rate.up*(yy-xx))[yy>xx], na.rm=TRUE ) 
        } } ) )
up.rate <- function (x) {
    # parameter for (conditioned) exponential distr'n of observed-true length given true length of x if observed > true
    fit.up$par[1] + 1/(fit.up$par[2]*x)
}
if (!produce.plots) { layout(t(1:2)) }
hist.up <- with(tblocks, hist( (yy-xx)[yy>xx], breaks=100, plot=FALSE ) )
with(hist.up, plot( density ~ mids, log='y', xlab="block length (cM)", main=expression(paste("fitted form for ",lambda['+'] (x))) ) )
with(tblocks, lines( hist.up$mids, sapply(hist.up$mids, function (y) { mean( up.rate(xx)*exp(-up.rate(xx)*y) ) } ) ) )
kfits <- sapply( 1:nlevels(tblocks$tbins), function (k) {
    ttt <- subset( tblocks, tbins==levels(tbins)[k] )
    hist.up <- with(ttt, hist( (yy-xx)[yy>xx], breaks=100, plot=FALSE ) )
    with(hist.up, points( density ~ mids, col=rainbow(nlevels(tblocks$tbins))[as.numeric(ttt$tbins[1])] ) )
    with(ttt, 1/mean((yy-xx)[yy>xx],na.rm=TRUE) )
} )
legend("topright",legend=levels(tblocks$tbins),col=rainbow(nlevels(tblocks$tbins)),pch=1)
if (!produce.plots) {
    plot( midbins, kfits )
    lines( midbins, up.rate(midbins) )
}
fit.up$par
# [1]    1.399283 3487.896989


# and those with y<x
fit.down <- with(tblocks, optim( par=c(.3,.05), f=function(r) { 
            rate.down <- r[1] + 1/(r[2]*xx)
            if (any(rate.down<0,na.rm=TRUE)) { Inf } else { 
                mean( -log(rate.down[yy<xx&!is.na(yy)]) + (rate.down*(xx-yy))[yy<xx&!is.na(yy)], na.rm=TRUE ) 
            } } ) )
down.rate <- function (x) {
    # parameter for (conditioned) exponential distr'n of observed-true length given true length of x if observed < true
    pmin( 12, fit.down$par[1] + 1/(fit.down$par[2]*x) )
}
hist.down <- with(tblocks, hist( (xx-yy)[yy<xx&!is.na(yy)], breaks=100, plot=FALSE ) )
with(hist.down, plot( density ~ mids, log='y', ylim=c(.01, 30), xlab="block length (cM)", main=expression(paste("fitted form for ",lambda['-'](x))) ) )
with(tblocks, lines( hist.down$mids, sapply(hist.down$mids, function (y) { mean(( down.rate(xx)*exp(-down.rate(xx)*y)/(1-exp(-down.rate(xx)*xx)) )[xx>y] ) } ) ) )
kfits <- sapply( 1:nlevels(tblocks$tbins), function (k) {
        if (with(tblocks, sum( tbins==levels(tbins)[k] & yy<xx & !is.na(yy) ) >10)) {
            ttt <- subset( tblocks, tbins==levels(tbins)[k] )
            hist.down <- with(ttt, hist( (xx-yy)[yy<xx&!is.na(yy)], breaks=100, plot=FALSE ) )
            with(hist.down, points( density ~ mids, col=rainbow(nlevels(tblocks$tbins))[as.numeric(ttt$tbins[1])] ) )
            with(ttt, 1/mean((xx-yy)[yy<xx&!is.na(yy)],na.rm=TRUE))
        } else { NA }
} )
legend("topright",legend=levels(tblocks$tbins),col=rainbow(nlevels(tblocks$tbins)),pch=1)
if (!produce.plots) {
    plot( midbins, kfits )
    lines( midbins, down.rate(midbins) )
}
fit.down$par
# [1] 0.4009342 0.1816122

# Done with plot for paper
if (produce.plots) { dev.off() }
