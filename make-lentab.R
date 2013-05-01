source("ibd-blocks-fns.R")
load("all-blocks-winnowed-fine.Rdata")
load("eda-data-fine.Rdata")

if (!file.exists("lentab.Rdata")) {
    lenbreaks <- c(1,3,5,100)
    blocks$country1 <- as.ordered(blocks$country1)
    blocks$country2 <- as.ordered(blocks$country2)

    all( indivinfo$YESOK ) # TRUE ? else remove.qc=TRUE
    # temporary computations
    zzz <- with( subset(blocks,maplen>min(lenbreaks)), table( factor(id1,levels=(indivinfo$SUBJID)), country2, cut(maplen,lenbreaks) ) )
    zzz <- zzz + with( subset(blocks,maplen>min(lenbreaks)), table( factor(id2,levels=(indivinfo$SUBJID)), country1, cut(maplen,lenbreaks) ) )
    nbl <- ( tapply( 1:dim(zzz)[1], indivinfo$COUNTRY_SELF, function (kk) { colSums( zzz[kk,,,drop=FALSE] ) } ) )  # sum of numbers of blocks shared
    ynames <- list( dimnames(nbl[[1]])[[1]], dimnames(nbl[[1]])[[2]], names(nbl) )
    nbl2 <- ( tapply( 1:dim(zzz)[1], indivinfo$COUNTRY_SELF, function (kk) { colSums( zzz[kk,,,drop=FALSE]^2 ) } ) ) # sum of squared numbers
    zzz <- with( subset(blocks,maplen>min(lenbreaks)), tapply( maplen, list( factor(id1,levels=(indivinfo$SUBJID)), country2, cut(maplen,lenbreaks) ), sum ) )  
    # zzz has NA if no such blocks
    zzz[is.na(zzz)] <- 0
    tzzz <- with( subset(blocks,maplen>min(lenbreaks)), tapply( maplen, list( factor(id2,levels=(indivinfo$SUBJID)), country1, cut(maplen,lenbreaks) ), sum ) )
    tzzz[is.na(tzzz)] <- 0
    zzz <- zzz + tzzz
    totbl <- ( tapply( 1:dim(zzz)[1], indivinfo$COUNTRY_SELF, function (kk) { colSums( zzz[kk,,,drop=FALSE] ) } ) )  # sum of lengths of blocks shared
    znames <- list( dimnames(totbl[[1]])[[1]], dimnames(totbl[[1]])[[2]], names(totbl) )
    all(unlist(ynames)==unlist(znames)) # TRUE
    totbl2 <- ( tapply( 1:dim(zzz)[1], indivinfo$COUNTRY_SELF, function (kk) { colSums( zzz[kk,,,drop=FALSE]^2 ) } ) ) # sum of squared lengths
    # put into df
    lentab <- unlist(nbl); dim(lentab) <- sapply(ynames,length); dimnames(lentab) <- ynames
    lentab <- as.data.frame.table( lentab ); colnames(lentab) <- c("countryY","maplen","countryX","sumbl")
    lentab <- lentab[ c("countryX","countryY","maplen","sumbl")]
    levels(lentab$maplen) <- levels(lentab$maplen)
    lentab$countryX <- ordered(lentab$countryX,levels=levels(blocks$country1))
    lentab$countryY <- ordered(lentab$countryY,levels=levels(blocks$country1))
    lentab$countrypair <- factor( with(lentab, gsub(" ",".", paste( levels(lentab$countryX)[ifelse(countryX<countryY,countryX,countryY)], levels(lentab$countryX)[ifelse(countryX<countryY,countryY,countryX)], sep="-" ) ) ) )
    lentab$sumbl2 <- unlist(nbl2)
    lentab$tmaplen <- unlist(totbl)
    lentab$tmaplen2 <- unlist(totbl2)
    lentab$npairs <- with(lentab, ifelse( countryX==countryY, choose(nsamples[countryX],2), nsamples[countryX]*nsamples[countryY] ) )
    lentab$nblocks <- with(lentab, ifelse(countryX==countryY,1/2,1) * sumbl )
    lentab$ibd <- with(lentab, nblocks / npairs )
    # Mean and SD of total length per other individual
    lentab$meanlen <- with(lentab, ifelse(countryX==countryY,1/2,1) * tmaplen / npairs )
    lentab$sdlen <- with( lentab, sqrt( (1/(nsamples[countryX]-1)) * ( tmaplen2 - tmaplen^2 / nsamples[countryX] ) ) )
    # Mean and SD of numbers of blocks that X-indivs share with country Y
    #  ... divide these by nsamples[countryY] to get naturally comparable things between countries.
    lentab$meanblocks <- with(lentab, sumbl / nsamples[countryX] )
    lentab$sdblocks <- with(lentab, sqrt( (1/(nsamples[countryX]-1)) * ( sumbl2 - sumbl^2 / nsamples[countryX] ) ) )
    lentab$sdibd <- lentab$sdblocks / nsamples[lentab$countryY]
    # Correlations in sharing rates
    lentab$ccor <- NA
    cmx <- as.array( xtabs( ibd ~ countryX + countryY + maplen, data=lentab ) )  # equivalently, xtabs( meanblocks /nsamples[countryY] ~ countryX + countryY + maplen, data=lentab )
    cmx[ rep( outer( dimnames(cmx)[[1]] , dimnames(cmx)[[2]], "==" ), dim(cmx)[3] ) ] <- NA
    for (x in levels(lentab$countryX)) { 
        for (y in levels(lentab$countryX)) {
            for ( k in 1:nlevels(lentab$maplen) ) {
                lentab$ccor[ lentab$countryX==x & lentab$countryY==y & as.numeric(lentab$maplen)==k ] <- cor( (cmx[x,,k]), (cmx[y,,k]), use="complete.obs", method="pearson" )
            }
        }
    }

    lentab$gdist <- poppairs$gdist[ match(lentab$countrypair,poppairs$countrypair) ] 
    lentab$cex <- poppairs$cex[ match(lentab$countrypair,poppairs$countrypair) ] 

    # Substructure:
    # significance for sdibd?  Do a permutation test.
    lentab$p.sd <- NA  # a p-value
    lentab$z.sd <- NA  # a z-score
    for (c1 in levels(blocks$country1)) {
        all.ids <- indivinfo$SUBJID[indivinfo$COUNTRY_SELF==c1]
        if (length(all.ids)>1) {
            subblocks <- with( subset( blocks, (country1==c1 | country2==c1) & country1!=country2 & maplen>min(lenbreaks)), 
                    data.frame( idA=ifelse(country1==c1,id1,id2), 
                        idB=ifelse(country1==c1,id2,id1), 
                        countryB=factor( levels(blocks$country1)[ifelse(country1==c1,country2,country1)], levels=levels(blocks$country1) ), 
                        maplen=cut(maplen,lenbreaks)
                        )
                    )
            boots <- replicate( 1000, {
                    subblocks$idA <- sample( all.ids, nrow(subblocks), replace=TRUE )
                    zzz <- with( subblocks, table( list( factor(idA,levels=all.ids), countryB, maplen ) ) )
                    ( colSums( zzz^2, dims=1 ) ) # dimensions of this are (countryY,maplen)
                } )
            lentab$p.sd[ lentab$countryX==c1 ] <- rowMeans( boots>=lentab$sumbl2[lentab$countryX==c1], dims=2 )  # right tail p-value
            lentab$z.sd[ lentab$countryX==c1 ] <- (lentab$sumbl2[lentab$countryX==c1] - rowMeans( boots, dims=2 ) ) / sqrt( rowMeans( sweep(boots,c(1,2),rowMeans(boots,dims=2),"-")^2, dims=2 ) )  # right tail "z-value"
        }
    }

    newcats <- list( 
            I=c("Italy","Spain","Portugal"),
            W=c("France", "United Kingdom", "Scotland", "England", "Ireland", "Swiss German", "Swiss French", "Switzerland", "Belgium", "Netherlands", "Germany" ),
            N=c( "Sweden", "Norway", "Denmark", "Latvia", "Finland" ),
            E=c( "Slovakia", "Greece", "Yugoslavia", "Albania", "Bosnia", "Montenegro", "Macedonia", "Kosovo", "Serbia", "Bulgaria", "Romania", "Poland", "Hungary", "Czech Republic", "Russia", "Slovenia", "Ukraine", "Croatia", "Austria"),
            TC=c("Turkey","Cyprus")
        )
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

    ## SAVE:
    save(lentab,lenbreaks,newcat,newcats,catcols,smcat.cols,ccatcols,file="lentab.Rdata")
} else {
    load("lentab.Rdata")
}

