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
#!/usr/bin/R
# Do some EDA on the final set of blocks

source("ibd-blocks-fns.R")

if (!file.exists("all-blocks-winnowed-fine.Rdata")) {
    blocks <- lapply(1:22, function (chrom) { getblocks(chrom, "winnowed", resolution="fine", remove.qc=TRUE) } )
    blocks <- do.call(rbind, blocks)
    save(blocks, file="all-blocks-winnowed-fine.Rdata")
} else {
    load(file="all-blocks-winnowed-fine.Rdata")
}

if (!file.exists("all-blocks-winnowed-with-rellys.Rdata")) {
    sampleinfo <- getsampleinfo()
    okbutrellys <- with(sampleinfo, SUBJID[ KEEP_EURO & !DROPPED_IN_PCA & !MIXED ] )
    blocks <- lapply(1:22, function (chrom) { 
            bl <- getblocks(chrom, "winnowed", remove.qc=FALSE) 
            subset( bl, ( id1 %in% okbutrellys & id2 %in% okbutrellys ) )
        } )
    blocks <- do.call(rbind, blocks)
    save(blocks, file="all-blocks-winnowed-with-rellys.Rdata")
}

## NOTES:
# what to do about short (length-zero) blocks?

edadata.filename <- "eda-data-fine.Rdata"
if (file.exists(edadata.filename)) {
    load(edadata.filename)
    # make a variable which is country-pair
    blocks$countrypair <- countrypairs[cbind(as.numeric(blocks$country1),as.numeric(blocks$country2))]
    # and individual-pair
    blocks$indivpair <- factor( paste( pmin(blocks$id1,blocks$id2), pmax(blocks$id1,blocks$id2), sep="-" ) )
} else {
    # Get info on the individuals, and drop the excluded ones
    indivs <- unique( c(blocks$id1,blocks$id2) )
    indivinfo <- getsampleinfo(resolution="fine",remove.qc=TRUE)
    # check:
    sum( ! indivs %in% indivinfo$SUBJID )  # == 0?
    sum( !indivinfo$YESOK[match(indivs, indivinfo$SUBJID)] ) # == 0?

    # Find top 6 regions (XXX may want to pick a different set of regions?)
    indivinfo$COUNTRY_SELF <- factor(indivinfo$COUNTRY_SELF, levels=levels(blocks$country1))
    nsamples <- table(indivinfo$COUNTRY_SELF)
    countries <- levels(indivinfo$COUNTRY_SELF)
    bigcountries <- names(sort(nsamples,decreasing=TRUE))[1:8] # "Switzerland" "United Kingdom" "Italy" "Spain" "Portugal" "France" 
    bigcountrypairs <- outer(bigcountries,bigcountries,function(x,y) { 
            u <- gsub(" ",".",x,fixed=TRUE)
            v <- gsub(" ",".",y,fixed=TRUE)
            paste(pmin(u,v),pmax(u,v),sep="-") 
        } )[upper.tri(diag(bigcountries),diag=TRUE)]

    # Colors for countries
    countrycols <- rainbow( nlevels(blocks$country1) )
    names(countrycols) <- levels(blocks$country1)
    # Abbreviations
    country.abbrev.list <- do.call( rbind, list( 
            c("Belgium", "BE"), 
            c("France", "FR"), 
            c("Austria", "AT"), 
            c("Bulgaria", "BG"), 
            c("Italy", "IT"), 
            c("Poland", "PL"), 
            c("Czech Republic", "CZ"), 
            c("Cyprus", "CY"), 
            c("Portugal", "PT"), 
            c("Denmark", "DK"), 
            c("Latvia", "LV"), 
            c("Romania", "RO"), 
            c("Germany", "DE"), 
            c("Lithuania", "LT"), 
            c("Slovenia", "SI"), 
            c("Estonia", "EE"), 
            c("Luxembourg", "LU"), 
            c("Slovakia", "SK"), 
            c("Ireland", "IE"), 
            c("Hungary", "HU"), 
            c("Finland", "FI"), 
            c("Greece", "EL"), 
            c("Malta", "MT"), 
            c("Sweden", "SE"), 
            c("Spain", "ES"), 
            c("Netherlands", "NL"), 
            c("United Kingdom", "UK"), 
            c("Iceland", "IS"), 
            c("Norway", "NO"), 
            c("Liechtenstein", "LI"), 
            c("Switzerland", "CH"), 
            c("Croatia", "HR"), 
            c("Montenegro", "ME"), 
            c("Turkey", "TR"), 
            c("Albania", "AL"), 
            c("Serbia", "RS"), 
            c("Bosnia and Herzegovina", "BA"), 
            c("Armenia", "AM"), 
            c("Belarus", "BY"), 
            c("Georgia", "GE"), 
            c("Azerbaijan", "AZ"), 
            c("Moldova", "MD"), 
            c("Ukraine", "UA"), 
            c("Algeria", "DZ"), 
            c("Lebanon", "LB"), 
            c("Syria", "SY"), 
            c("Russia", "RU"), 
            c("Swiss French", "CHf"), 
            c("Swiss German", "CHd"), 
            c("Yugoslavia", "YU"),
            c("Bosnia", "BO"),     
            c("Croatia", "CR"),    
            c("England", "EN"),    
            c("Kosovo", "KO"),
            c("Macedonia", "MA"),
            c("Montenegro", "MO"),
            c("Scotland", "SC"),
            c("Serbia", "SR")
        ) )
    countryabbrevs <- country.abbrev.list[ match(levels(blocks$country1),country.abbrev.list[,1]), 2 ]
    names(countryabbrevs) <- levels(blocks$country1)

    # Categories
    easterns <- c("Albania", "Kosovo", "Slovakia", "Greece", "Yugoslavia", "Bulgaria", "Romania", "Poland", "Hungary", "Czech Republic", "Russia", "Slovenia","Ukraine", "Croatia","Bosnia","Montenegro","Macedonia", "Serbia" )
    northerns <- c("Latvia", "Finland", "Sweden", "Norway", "Denmark")
    mideasterns <- c("Cyprus","Turkey")
    southerns <- c("Italy","Portugal","Spain")
    westerns <- c("France","United Kingdom","England", "Scotland", "Ireland","Belgium","Netherlands", "Switzerland", "Swiss French","Swiss German", "Germany", "Austria")
    catlist <- list( E=easterns, N=northerns, TC=mideasterns, I=southerns, W=westerns )
    # catlist <- list( SW=southerns, NW=westerns, NE=northerns, CE=easterns, SE=mideasterns )
    # countrycats <- unlist( lapply( names(catlist), function (x) rep(x,length(catlist[[x]])) ) )
    # names(countrycats) <- unlist(catlist)

    # catlist <- list( 
    #         I=c("Italy","Spain","Portugal"),
    #         W=c("France", "United Kingdom", "Ireland", "Swiss German", "Swiss French", "Switzerland", "Belgium", "Netherlands", "Germany", "Denmark", "Norway", "Sweden"), 
    #         E=c("Finland", "Austria", "Slovenia", "Latvia", "Poland", "Russia", "Czech Republic", "Hungary", "Romania", "Yugoslavia", "Albania", "Greece", "Ukraine", "Bulgaria", "Slovakia"),
    #         TC=c("Turkey","Cyprus")
    #     )
    countrycats <- rep(names(catlist),times=sapply(catlist,length)); names(countrycats) <- unlist(catlist)
    countrycat <- function (countryX) {
        if (is.factor(countryX)) { countryX <- levels(countryX)[as.numeric(countryX)] }
        countrycats[countryX]
    }
    catpair <- function (countryX,countryY,collapse=TRUE) {
        tmp.X <- countrycat(countryX)
        tmp.Y <- countrycat(countryY)
        smcat <- as.factor( paste( ifelse(tmp.X<tmp.Y,tmp.X,tmp.Y), ifelse(tmp.X<tmp.Y,tmp.Y,tmp.X), sep="-" ) )
        if (collapse) {
            levels( smcat ) <- c("E-E"="E-E", "E-I"="I", "E-TC"="TC", "E-W"="E-W", "E-N"="E-E", "N-N"="E-E", "I-N"="I", "N-TC"="TC", "N-W"="E-W", "I-I"="I", "I-TC"="TC", "I-W"="I", "TC-TC"="TC", "TC-W"="TC", "W-W"="W-W")[ levels(smcat) ]
        }
        return( smcat )
    }
    catcols <- c(E="red", W="green", I="orange", TC="magenta")
    ccatcols <- catcols[countrycats]; names(ccatcols) <- names(countrycats)

    # which blocks are in big countries
    bigblocks <- (blocks$country1 %in% bigcountries) & (blocks$country2 %in% bigcountries)

    countrypairs <- gsub(' ', '.', outer(levels(blocks$country1),levels(blocks$country2),paste,sep="-") )
    countrypairs[row(countrypairs)>col(countrypairs)] <- t(countrypairs)[row(countrypairs)>col(countrypairs)]
    # make a variable which is country-pair
    blocks$countrypair <- countrypairs[cbind(as.numeric(blocks$country1),as.numeric(blocks$country2))]
    # and individual-pair
    blocks$indivpair <- factor( paste( pmin(blocks$id1,blocks$id2), pmax(blocks$id1,blocks$id2), sep="-" ) )
    blocks$country1 <- as.ordered(blocks$country1)
    blocks$country2 <- as.ordered(blocks$country2)

    # compute summaries by country pair
    poppairs <- ddply( blocks, "countrypair", summarise,
                country1=min(country1[1],country2[1]), 
                country2=max(country1[1],country2[1]), 
                # country1=country1[1], country2=country2[1],
                nblocks=length(maplen),                            # number of blocks 
                nblocks1cM=sum(maplen>=1&maplen<5),                          # number of blocks 1cM < * < 5cM
                nblocks5cM=sum(maplen>=5&maplen<10),                          # number of blocks 5cM < * < 10cM
                nblocks10cM=sum(maplen>=10),                        # number of blocks >10cM
                nsharepairs=nrow(unique(data.frame(id1,id2))),   # number of pairs sharing anything
                nchroms=nrow(unique(data.frame(id1,id2,chrom))), # number of chromosomes shared
                totalmaplen=sum(maplen),                           # total map length shared
                totalmaplen.1cM=sum(maplen[maplen>1]),             # total map length shared in blocks >1cM
                gdist=gdist[1],                                    # geographic distance
                npairs=ifelse( country1[1]==country2[1], choose(nsamples[country1[1]],2), nsamples[country1[1]]*nsamples[country2[1]] ),          # number of pairs (i.e. the denominator)
            .progress="text"
        )
    poppairs$cex <- pmax(sqrt(poppairs$npairs)/50,.25)   #  point sizes reflecting sample sizes

    # And summaries by individual pair
    indpairs <- with(blocks, tapply( maplen, indivpair, length ) )
    indpairs <- data.frame( list( nblocks=indpairs ) )
    indpairs$id1 <- with(blocks, tapply(id1, indivpair, function(x) x[1] ) )
    indpairs$id2 <- with(blocks, tapply(id2, indivpair, function(x) x[1] ) )
    indpairs$country1 <- factor( with(blocks, tapply(as.numeric(country1), indivpair, function(x) x[1] ) ) )
    levels(indpairs$country1) <- levels(blocks$country1)
    indpairs$country2 <- factor( with(blocks, tapply(as.numeric(country2), indivpair, function(x) x[1] ) ) )
    levels(indpairs$country2) <- levels(blocks$country2)
    indpairs$nblocks1cM <- with(blocks, tapply( maplen, indivpair, function (x) sum(x>=1 & x<5) ) )
    indpairs$nblocks5cM <- with(blocks, tapply( maplen, indivpair, function (x) sum(x>=5 & x<10) ) )
    indpairs$nblocks10cM <- with(blocks, tapply( maplen, indivpair, function (x) sum(x>=10) ) )
    indpairs$nchroms <- with(blocks, tapply( chrom, indivpair, function (x) length(unique(x)) ) )
    indpairs$totalmaplen <- with(blocks, tapply( maplen, indivpair, sum ) )
    indpairs$totalmaplen.1cM <- with(blocks, tapply( maplen, indivpair, function (x) sum(x[x>1]) ) )
    indpairs$gdist <- with(blocks, tapply( gdist, indivpair, function(x) x[1] ) )
    indpairs$countrypair <- factor( countrypairs[cbind(as.numeric(indpairs$country1),as.numeric(indpairs$country2))] )

    # Also by individual:
    #  total number of blocks
    for ( cname in c("nblocks","nblocks1cM","nblocks5cM","nblocks10cM","totalmaplen") ) {
        x <- tapply( indpairs[,cname], indpairs[,"id1"], sum ) 
        y <- tapply( indpairs[,cname], indpairs[,"id2"], sum )
        indivinfo[,cname] <- 0
        indivinfo[,cname][na.omit(match(names(x),indivinfo$SUBJID))] <- x 
        indivinfo[,cname][na.omit(match(names(y),indivinfo$SUBJID))] <- y + indivinfo[,cname][na.omit(match(names(y),indivinfo$SUBJID))]
    }
    # total number of blocks by country
    indiv.countryblocks1 <- do.call(rbind, with(blocks, tapply(country2, id1, table) ) )
    indiv.countryblocks2 <- do.call(rbind, with(blocks, tapply(country1, id2, table) ) )
    indiv.countryblock.names <- paste(gsub(" ",".",levels(blocks$country1),fixed=TRUE),"blocks",sep=".")
    for ( dn in indiv.countryblock.names ) { indivinfo[,dn] <- 0 }
    indivinfo[ na.omit(match(dimnames(indiv.countryblocks1)[[1]],indivinfo$SUBJID)), paste(gsub(" ",".",dimnames(indiv.countryblocks1)[[2]]),"blocks",sep=".") ] <- indiv.countryblocks1
    indivinfo[ na.omit(match(dimnames(indiv.countryblocks2)[[1]],indivinfo$SUBJID)), paste(gsub(" ",".",dimnames(indiv.countryblocks2)[[2]]),"blocks",sep=".") ] <- indiv.countryblocks2 + indivinfo[ na.omit(match(dimnames(indiv.countryblocks2)[[1]],indivinfo$SUBJID)), paste(gsub(" ",".",dimnames(indiv.countryblocks2)[[2]]),"blocks",sep=".") ] 
    lims <- c(1,5,10,Inf)
    for ( k in 1:3 ) {
        usethese <- (blocks$maplen > lims[k]) & (blocks$maplen<=lims[k+1])
        indiv.countryblocks1 <- do.call(rbind, with(blocks[usethese,], tapply(country2, id1, table) ) )
        indiv.countryblocks2 <- do.call(rbind, with(blocks[usethese,], tapply(country1, id2, table) ) )
        indiv.countryblock.names <- paste(gsub(" ",".",levels(blocks$country1),fixed=TRUE),".blocks.",lims[k],"cM",sep="")
        for ( dn in indiv.countryblock.names ) { indivinfo[,dn] <- 0 }
        indivinfo[ na.omit(match(dimnames(indiv.countryblocks1)[[1]],indivinfo$SUBJID)), paste(gsub(" ",".",dimnames(indiv.countryblocks1)[[2]]),".blocks.",lims[k],"cM",sep="") ] <- indiv.countryblocks1
        indivinfo[ na.omit(match(dimnames(indiv.countryblocks2)[[1]],indivinfo$SUBJID)), paste(gsub(" ",".",dimnames(indiv.countryblocks2)[[2]]),".blocks.",lims[k],"cM",sep="") ] <- indiv.countryblocks2 + indivinfo[ na.omit(match(dimnames(indiv.countryblocks2)[[1]],indivinfo$SUBJID)), paste(gsub(" ",".",dimnames(indiv.countryblocks2)[[2]]),".blocks.",lims[k],"cM",sep="") ] 
    }


    # [x,y] gives the variance of the number of x blocks, across individuals of country y.
    #   note -- not symmetric.
    nblock.vars <- c( 
            sapply( levels(blocks$country1), function (ccc) {
                cname <- gsub(" ",".",paste(ccc,".blocks.1cM",sep=''))
                tapply( indivinfo[,cname], indivinfo[,"COUNTRY_SELF"], var )/nsamples[ccc]^2
            } ),
            sapply( levels(blocks$country1), function (ccc) {
                cname <- gsub(" ",".",paste(ccc,".blocks.5cM",sep=''))
                tapply( indivinfo[,cname], indivinfo[,"COUNTRY_SELF"], var )/nsamples[ccc]^2
            } ),
            sapply( levels(blocks$country1), function (ccc) {
                cname <- gsub(" ",".",paste(ccc,".blocks.10cM",sep=''))
                tapply( indivinfo[,cname], indivinfo[,"COUNTRY_SELF"], var )/nsamples[ccc]^2
            } )
        )
    dim(nblock.vars) <- c( nlevels(blocks$country1), nlevels(blocks$country1), 3 )
    dimnames(nblock.vars) <- list( levels(blocks$country1), levels(blocks$country1), c("1cM","5cM","10cM") )

    # which pairs of individuals are in big countries
    bigindpairs <- (indpairs$country1 %in% bigcountries) & (indpairs$country2 %in% bigcountries)

    save(poppairs, indpairs, indivinfo, nblock.vars, nsamples, countrycols, countryabbrevs, catlist, countrycats, countrypairs, bigblocks, bigindpairs, bigcountrypairs, file="eda-data-fine.Rdata")
}
