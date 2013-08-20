# Written 2013 by Peter Ralph and Graham Coop
# 
# contact: petrel.harp@gmail.com
#
#     To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty. 
# 
# You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 
# 
#
#!/usr/bin/R

###
# Do sensitivity analysis with simulation results (produced by coalpedigree)
source("laplace-inversion-fns.R")
source("ibd-blocks-fns.R")
source("sim-params.R")  # gets simopts and coalprob()
source("parse-sims-fns.R")  # gets simopts and coalprob()

ploidy <- 2
prefixes <- c(
    "fixed"="28176-growing-migration-1",    # fixed size
    "old"="16062-demographics-expansion", # expanding 100 gens ago
    "young"="26758-demographics-expansion",  # expanding 50 gens ago
    "complex"="27854-demographics-complex"   # complex
    #"28860-migration-3"            # two pops, one expands...
)

# ( for x in 28176-growing-migration-1 16062-demographics-expansion 26758-demographics-expansion 27854-demographics-complex; do python ~/projects/coalpedigree/winnow.py -b $x.fibd.gz -l ${x}-winnowed.log -o ${x}-winnowed.fibd.gz -g .05 -n 0.005 -m 0.0; done ) &

if (!file.exists("various-inversions.Rdata")) {
    minminlen <- 0.5
    full.inversions <- lapply(prefixes,do.everything,minblocklen=minminlen)
    # moretime.inversions <- lapply(prefixes,do.everything,maxgen=1000)
    # midtime.inversions <- lapply(prefixes,do.everything,maxgen=500)
    long.inversions <- lapply(prefixes,do.everything,minblocklen=2)
    # long.moretime.inversions <- lapply(prefixes,do.everything,minblocklen=2,maxgen=1000)
    # long.midtime.inversions <- lapply(prefixes,do.everything,minblocklen=2,maxgen=500)
    fp.inversions <- lapply(prefixes,do.everything,minblocklen=2,fprate=function(x) (1)*exp(-(x-2)))
    little.fp.inversions <- lapply(prefixes,do.everything,minblocklen=2,fprate=function(x) (.1)*exp(-(x-2)))
    longer.fp.inversions <- lapply(prefixes,do.everything,minblocklen=2,fprate=function(x) (.1)*exp(-(x-2)/4))

    # save(prefixes,full.inversions,moretime.inversions,midtime.inversions,long.inversions,long.moretime.inversions,long.midtime.inversions,fp.inversions,little.fp.inversions,longer.fp.inversions,file="various-inversions.Rdata")
    save(prefixes,minminlen,full.inversions,long.inversions,fp.inversions,little.fp.inversions,longer.fp.inversions,file="various-inversions.Rdata")
} else {
    load("various-inversions.Rdata")
}

# or load("various-inversions-nonwinnowed.Rdata")

pdf(file=paste("sensitivity-results.pdf",sep="-"), width=12, height=3*3.5, pointsize=10)
ylims <- lapply(full.inversions, function (x) range(x$anslist[['a-a']][[1]]$par) )
layout(matrix(1:(4*sum(sapply(full.inversions,function(x)length(x$anslist)))),nrow=4))
par(mar=c(3,3,1,1)+.1)
lapply( seq_along(full.inversions), function (k) {
        x <- full.inversions[[k]]
        y <- long.inversions[[k]]
        for (thispair in names(x$anslist)) {
            plot.ans(x$anslist[[thispair]], x$opts, thispair, x$L, ylim=ylims[[k]], main=paste(prefixes[k],"all") )
            plot.ans(y$anslist[[thispair]], y$opts, thispair, y$L, ylim=ylims[[k]], main=">2cM" )
        }
    } )
layout(matrix(1:(6*sum(sapply(full.inversions,function(x)length(x$anslist)))),nrow=6))
if (exists("moretime.inversions")) {
    lapply( seq_along(full.inversions), function (k) {
            x <- midtime.inversions[[k]]
            y <- long.midtime.inversions[[k]]
            z <- moretime.inversions[[k]]
            for (thispair in names(x$anslist)) {
                plot.ans(x$anslist[[thispair]], x$opts, thispair, x$L, ylim=ylims[[k]], main=paste(prefixes[k],"all, to 500 gens") )
                plot.ans(y$anslist[[thispair]], y$opts, thispair, y$L, ylim=ylims[[k]], main=">2cM to 500 gens" )
                plot.ans(z$anslist[[thispair]], z$opts, thispair, z$L, ylim=ylims[[k]], main="all, to 1000 gens" )
            }
        } )
}
lapply( seq_along(full.inversions), function (k) {
        x <- fp.inversions[[k]]
        y <- little.fp.inversions[[k]]
        z <- longer.fp.inversions[[k]]
        for (thispair in names(x$anslist)) {
            plot.ans(x$anslist[[thispair]], x$opts, thispair, x$L, ylim=ylims[[k]], main=paste(prefixes[k],"2x, 2-3cM") )
            plot.ans(y$anslist[[thispair]], y$opts, thispair, y$L, ylim=ylims[[k]], main="1.1x, 2--3cM" )
            plot.ans(z$anslist[[thispair]], y$opts, thispair, y$L, ylim=ylims[[k]], main="1.1x, 2--6cM" )
        }
    } )
dev.off()


# for paper:

pdf(file="spectra-comparisons.pdf",width=6.5,height=3,pointsize=10)
tmp <- lapply( c("fixed","old","young","complex"), function (x)  plot.ans( full.inversions[[x]], thispair='a-a', plots="none" ) )
cols <- rainbow_hcl(length(tmp))
layout(t(1:2))
par(mar=c(4,4,1,1)+.1)
ylims <- range( unlist( sapply( tmp, function (x) x$coal$theoretical ) ), na.rm=TRUE ) 
plot( 0, 0, type='n', xlab="generations ago", ylab="coalescent rate", xlim=c(0,300), ylim=ylims, log='y' )
lapply( seq_along(tmp), function (k) with(tmp[[k]], { 
            lines( gens, coal$theoretical, col=cols[k], lwd=2 ) 
        } ) )
ylims <- range( unlist( sapply( tmp, function (x) x$spectra$theoretical/x$npairs ) ) ) 
plot( 0, 0, type='n', xlab="block length (cM)", ylab="density", xlim=c(0.5,45), ylim=ylims, log='xy' )
lapply( seq_along(tmp), function (k) with(tmp[[k]], { 
            lines( midbins, (spectra$observed/npairs), col=cols[k], lwd=2 )
            lines( midbins, (spectra$theoretical/npairs), col=cols[k], lty=2, lwd=2 ) 
        } ) )
legend("topright",legend=LETTERS[1:4],lwd=2,col=cols)
dev.off()

# age distributions of different block lengths
pdf(file="age-distributions.pdf", width=6.5, height=3, pointsize=10)
plotlens <- c(2,5,8)
plotlims <- c(600,300,150)
dists <- lapply( full.inversions, function (x) {
        tmp <- with( x, predict.blocks( L, opts, times=TRUE ) )[,,1]
        lenbins <- attr(x$L,"lenbins")
        gens <- attr(x$L,"gens")
        usegens <- which(gens%%2==0)
        plotbins <- findInterval(plotlens, lenbins)
        p <- sweep( tmp[plotbins,usegens], 1, rowSums(tmp[plotbins,usegens]), "/" )
        # p <- tmp[plotbins,usegens]
        return( list( t=gens[usegens], p=p ) )
    })
layout(t(1:length(plotlens)))
for (k in seq_along(plotlens)) {
    plot( 0, 0, type='n', xlab="generations ago", ylab="density", xlim=c(0,plotlims[k]), ylim=range(unlist(sapply(dists,function(x)x$p[k,]))), main=paste("age distribution,", plotlens[k],"cM") )
    abline(v=50/plotlens[k],lty=2)
    invisible( lapply( seq_along(dists), function (j) {
                lines( dists[[j]]$t, dists[[j]]$p[k,], col=cols[j] )
            } ) )
}
legend("topright",legend=LETTERS[1:4],lty=1,col=cols)
dev.off()

pdf(file="full-inversions.pdf",width=6.5,height=5,pointsize=10)
layout(matrix(1:6,nrow=2))
plot.ans( full.inversions[["fixed"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main=c("(A)","(B)"), ylim=c(0,2e-5) )
plot.ans( full.inversions[["old"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main=c("(C)","(D)"), ylim=c(0,2e-5), legend1=FALSE )
plot.ans( full.inversions[["young"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main=c("(E)","(F)"), ylim=c(0,2e-5), legend1=FALSE )
dev.off()


pdf(file="long-inversions.pdf",width=6.5,height=5,pointsize=10)
layout(matrix(1:6,nrow=2))
plot.ans( long.inversions[["fixed"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main=c("(A)","(B)"), ylim=c(0,2e-5) )
plot.ans( long.inversions[["old"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main=c("(C)","(D)"), ylim=c(0,2e-5), legend1=FALSE )
plot.ans( long.inversions[["young"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main=c("(E)","(F)"), ylim=c(0,2e-5), legend1=FALSE )
dev.off()


pdf(file="complex-inversions.pdf", width=6.5, height=5, pointsize=10)
layout(matrix(1:6,nrow=2,byrow=TRUE))
plot.ans( full.inversions[["complex"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main="A", ylim=c(0,2e-5), plots="coal", legend1=FALSE )
plot.ans( full.inversions[["complex"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main="B", xlim=c(0,80), ylim=c(0,2e-5), plots="coal" )
plot.ans( full.inversions[["complex"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main="C", plots="spectrum" )
plot.ans( long.inversions[["complex"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main="D", ylim=c(0,2e-5), plots="coal", legend1=FALSE )
plot.ans( long.inversions[["complex"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main="E", xlim=c(0,80), ylim=c(0,2e-5), plots="coal" )
plot.ans( long.inversions[["complex"]], thispair='a-a', dothese=c("MLE"="pointy","smooth"="smooth"), main="F", plots="spectrum" )
dev.off()

pdf(file="fp-inversions.pdf",width=6.5,height=5,pointsize=10)
layout(matrix(1:9,nrow=3))
par(mar=c(4,4,2,1)+.1)
ylims <- lapply(full.inversions, function (x) range(x$anslist[['a-a']][[1]]$par) )
lapply( seq_along(full.inversions), function (k) {
        x <- fp.inversions[[k]]
        y <- little.fp.inversions[[k]]
        z <- longer.fp.inversions[[k]]
        for (thispair in names(x$anslist)) {
            plot.ans(y$anslist[[thispair]], y$opts, thispair, y$L, ylim=ylims[[k]], plots="coal", main=LETTERS[k], dothese=c("MLE"="pointy","smooth"="smooth"), legend1=(k==1) )
            plot.ans(z$anslist[[thispair]], y$opts, thispair, y$L, ylim=ylims[[k]], plots="coal", main=LETTERS[3+k], dothese=c("MLE"="pointy","smooth"="smooth"), legend1=FALSE  )
            plot.ans(x$anslist[[thispair]], x$opts, thispair, x$L, ylim=ylims[[k]], plots="coal", main=LETTERS[6+k], dothese=c("MLE"="pointy","smooth"="smooth"), legend1=FALSE )
        }
    } )
dev.off()

if (FALSE) {

    pdf(file="moretime-inversions.pdf",width=6.5,height=5,pointsize=10)
    layout(matrix(1:6,nrow=2))
    plot.ans( long.moretime.inversions[['fixed']], thispair='a-a', plots='coal', dothese=c("MLE"="pointy","smooth"="smooth"), main="(A)", ylim=c(0,5e-5), legend1=FALSE )
    plot.ans( long.moretime.inversions[['fixed']], thispair='a-a', plots='coal', dothese=c("MLE"="pointy","smooth"="smooth"), main="(E)", xlim=c(0,250), ylim=c(0,1e-5), legend1=FALSE )
    plot.ans( long.moretime.inversions[['old']], thispair='a-a', plots='coal', dothese=c("MLE"="pointy","smooth"="smooth"), main="(B)", ylim=c(0,2e-4), legend1=FALSE )
    plot.ans( long.moretime.inversions[['old']], thispair='a-a', plots='coal', dothese=c("MLE"="pointy","smooth"="smooth"), main="(F)" , xlim=c(0,250), ylim=c(0,5e-5) )
    plot.ans( long.moretime.inversions[['young']], thispair='a-a', plots='coal', dothese=c("MLE"="pointy","smooth"="smooth"), main="(C)", ylim=c(0,2e-4), legend1=FALSE )
    plot.ans( long.moretime.inversions[['young']], thispair='a-a', plots='coal', dothese=c("MLE"="pointy","smooth"="smooth"), main="(G)" , xlim=c(0,250), ylim=c(0,1.5e-5), legend1=FALSE )
    dev.off()

##
# box plots
genbins <- list( c(0,36), c(37,100), c(101,169), c(170,289), c(290,666) )  # roughly 500,1500,2500, 4000, and 10000 years ago
names(genbins) <- sapply( genbins, paste, collapse="-" )
bin.ans <- function (ans,gens,genbins) {
    # return mean coalescent rates across each of genbins
    sapply(genbins, function (gb) {
            mean( ans$par[ gens>gb[1] & gens<=gb[2] ] )
        } )
}

getbins <- function (everything,genbins) {
    L <- everything$L
    gens <- attr(L,"gens")
    lapply( everything$anslist, function (x) {
            binned <- cbind(
                    sapply( x[ setdiff(names(x),c("lower.bounds","upper.bounds")) ], bin.ans, gens=gens, genbins=genbins ),
                    "lower" = sapply( seq_along(genbins), function (k) bin.ans( x$lower.bounds[[k]], gens, genbins )[k] ),
                    "upper" = sapply( seq_along(genbins), function (k) bin.ans( x$upper.bounds[[k]], gens, genbins )[k] )
                )
                rownames(binned) <- sapply(genbins,paste,collapse="-")
                return(binned)
        } )
}

}
