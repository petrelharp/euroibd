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
# Process Europeans.out,
# making Euro-samples-info.csv
# for our use.

require(maps)
require(mapproj)
require(fields)  # for rdist.earth
require(plyr)

.basedir <- suppressWarnings( system("ls -d /home/ibd/ /home/peter/projects/ibd/", intern=TRUE, ignore.stderr=TRUE)[1] )
.pcadir <- paste(.basedir,"data/POPRES/pca_euro/",sep="")

# Information about of european samples
euro <- read.table("Europeans.out")
# List of countries, with 0 or 1 if in europe or not
in.Europe <- read.table("dbgap_In_Europe.out",as.is=T)
europe.countries<-in.Europe[in.Europe$V1==1,2]

# Combine certain countries/locations together.
combine.country<-list()
# a finer grouping
combine.country <- c( combine.country, list(c("Russia","USSR")) )
combine.country <- c( combine.country, list(c("Netherlands","Holland")) )

# Make new COUNTRY_SELF that agrees with grandparents
# record original information
euro$ORIG_COUNTRYSELF<-euro$COUNTRY_SELF
gfolx <- as.matrix( euro[,c("COUNTRY_MGM","COUNTRY_MGF","COUNTRY_PGM","COUNTRY_PGF")] )
##For individuals with no granf
par.country <- as.matrix( euro[,c("COUNTRY_FATHER","COUNTRY_MOTHER")] )

euro$MIXEDGFOLX <- apply(gfolx, 1, function (x) length(unique(x))>1 )
num.grandpar<-apply(gfolx,1,function(grand){4-sum(grand=="")}) 


# Is the individual "mixed"?
euro$MIXED <- euro$GROUPING_PCA_LABEL1 == "Mix" | euro$COUNTRY_SELF=="Europe" | euro$MIXEDGFOLX
euro$COUNTRY_SELF <- NA
euro$COUNTRY_GFOLX <- NA
euro$COUNTRY_SELF[!euro$MIXED] <- euro$COUNTRY_GFOLX[!euro$MIXED] <- gfolx[!euro$MIXED,1]
## If no grandpar info. available use country self
euro$COUNTRY_SELF[num.grandpar==0]<-as.character(euro$ORIG_COUNTRYSELF[num.grandpar==0]) 

# combine countries
for (i in 1:length(combine.country)) {
    euro$COUNTRY_SELF[euro$COUNTRY_SELF %in% combine.country[[i]]]<-combine.country[[i]][1]
}

# split countries by language
# Switzerland
euro$COUNTRY_SELF[ euro$COUNTRY_SELF=="Switzerland" & euro$PRIMARY_LANGUAGE=="French" ] <- "Swiss French"
euro$COUNTRY_SELF[ euro$COUNTRY_SELF=="Switzerland" & euro$PRIMARY_LANGUAGE=="German" ] <- "Swiss German"
# Balkans:
# with( subset(euro,COUNTRY_SELF %in% c("Yugoslavia","Albania","Serbia", "Macedonia", "Montenegro", "Croatia", "Kosovo", "Bosnia" )), table( COUNTRY_SELF, droplevels(PRIMARY_LANGUAGE) ) )
# COUNTRY_SELF    Albanian Bosnian Croatian French Hungarian Kosovan Macedonian Romanian Serbian Serbo-Croatian Yugoslavian
#   Albania     0        3       0        0      0         0       0          0        0       0              0           0
#   Bosnia      0        0       4        0      0         0       0          0        0       1              4           0
#   Croatia     0        0       0        7      0         0       0          0        0       0              1           0
#   Kosovo      1       11       0        0      0         0       2          0        0       0              2           1
#   Macedonia   0        0       0        0      0         0       0          4        0       0              0           0
#   Montenegro  0        0       0        0      0         0       0          0        0       1              0           0
#   Serbia      0        1       0        0      0         1       0          0        0       6              4           0
#   Yugoslavia  1        5       0        1      1         0       1          0        1       2              3           4
balkans <- c( "Albania", "Bosnia", "Croatia", "Kosovo", "Macedonia", "Montenegro", "Serbia", "Yugoslavia" )
balkan.langs <- c( "Albanian", "Bosnian", "Croatian", "Kosovan", "Macedonian", "Serbian", "Serbo-Croatian" )
# notes: "kosovan" probably = "albanian"
#  "serbo-croatian" includes serbian, croatian, bosnian and probably = "yugoslavian"
#   macedonian "forms a continuum" of south slavic languages with bulgarian and serbo-croatian
euro$COUNTRY_SELF[ euro$COUNTRY_SELF%in%c("Yugoslavia","Serbia") & euro$PRIMARY_LANGUAGE=="Albanian" ] <- "Albania"
euro$COUNTRY_SELF[ euro$COUNTRY_SELF%in%c("Yugoslavia") & euro$PRIMARY_LANGUAGE=="Croatian" ] <- "Croatia"
euro$COUNTRY_SELF[ euro$COUNTRY_SELF%in%c("Yugoslavia") & euro$PRIMARY_LANGUAGE=="Kosovan" ] <- "Kosovo"
euro$COUNTRY_SELF[ euro$COUNTRY_SELF%in%c("Yugoslavia") & euro$PRIMARY_LANGUAGE=="Serbian" ] <- "Serbia"


# Passed QC and all reported ancestors as European?
in.Europe<-read.table("dbgap_In_Europe.out",as.is=T)
europe.countries<-in.Europe[in.Europe$V1==1,2]
reduced.euro<-euro[euro$COUNTRY_SELF %in% europe.countries,]
ancs.in.euro<- apply(reduced.euro[,c("COUNTRY_FATHER","COUNTRY_MOTHER","COUNTRY_MGF","COUNTRY_MGM","COUNTRY_PGF","COUNTRY_PGM")],1,function(kin){ all(kin %in% c(europe.countries,"")) } )
keep.euro <- subset( reduced.euro, STATUS_PASSED_QC2=="Y" & ancs.in.euro )$SUBJID
euro$KEEP_EURO <- euro$SUBJID %in% keep.euro$V1

## read in geographic information
data(world.cities)
names(world.cities) <- c("CITY_SELF", "COUNTRY_SELF", "CITY_POP", "lat", "long", "capital")
# choose largest population city in each country
largest <- ddply(world.cities, "COUNTRY_SELF", function (x) x[which.max(x$CITY_POP),])

# UK -> United Kingdom
largest[largest$CITY_SELF=="London" & largest$COUNTRY_SELF=="UK","COUNTRY_SELF"] <- "United Kingdom"
# "Switzerland" city -> Bern
largest[largest$COUNTRY_SELF=="Switzerland",] <- world.cities[ world.cities$COUNTRY_SELF=="Switzerland" & world.cities$CITY_SELF=="Bern", ]
swfrench <- world.cities[world.cities$CITY_SELF =="Geneva", ]
swfrench$COUNTRY_SELF <- "Swiss French"
swgerman <- world.cities[world.cities$CITY_SELF =="Zurich", ]
swgerman$COUNTRY_SELF <- "Swiss German"
largest <- rbind(largest,swfrench,swgerman)
# Balkans
addthese <- list( 
        c( "Serbia", "Belgrade" ),
        c( "Bosnia", "Sarajevo" ),
        c( "Kosovo", "Pristina" ),
        c( "Montenegro", "Podgorica" ),
        c( "Yugoslavia", "Belgrade" ),
        c( "England", "London" ),
        c( "Wales", "Cardiff" ),
        c( "Scotland", "Glasgow" )
    )
for (x in addthese) {
    addthis <- world.cities[ world.cities$CITY_SELF==x[2], ]
    addthis <- addthis[which.max(addthis$CITY_POP),]
    addthis$COUNTRY_SELF <- x[1]
    largest <- rbind( largest, addthis )
}

## Attach geographic information to euro
euro <- merge( euro, largest[,1:5], by="COUNTRY_SELF", all.x=TRUE, sort=FALSE )

# Add in PCA info
pcas<-read.table(paste(.pcadir,"euro_nooutlier.pcavec",sep=""),skip=1,as.is=TRUE)
colnames(pcas)[-c(1,12)]<-paste("PC",1:10,sep="")
pcas$labels<-sapply(pcas$V1,function(x){strsplit(x,split="\\:")[[1]][1]})
euro <- merge(euro, pcas[,-match(c("V1","V12"),names(pcas))], by.x="SUBJID", by.y="labels", all.x=TRUE, sort=FALSE)

#  Related individuals according to kinship
related<-read.table(paste(.pcadir,"kinship.above.cutoff",sep=""),as.is=TRUE,head=TRUE)  #added by G
euro$CLOSE_REL<-FALSE
euro$CLOSE_REL[euro$SUBJID %in% related$IID1]<-TRUE  ##mark one of the pair to drop

# Dropped in PCA analysis 
euro$DROPPED_IN_PCA <- ! euro$SUBJID %in% pcas$labels

# Flag for "use these ones"
euro$YESOK <- euro$KEEP_EURO & !euro$DROPPED_IN_PCA & !euro$MIXED & !euro$CLOSE_REL

# Add new ID that groups indivs by country
euro$GEOGID <- do.call( order, euro[,c("COUNTRY_MGM","COUNTRY_MGF","COUNTRY_PGM","COUNTRY_PGF")] )

##Add % missing data
missing.data<-read.table(paste(.pcadir,"europeans.imiss",sep=""),head=TRUE,as.is=TRUE)
euro$F_MISS<-missing.data$F_MISS[match(euro$SUBJID,missing.data$IID)]

# Sort by COUNTRY_SELF
euro <- euro[ order(euro$COUNTRY_SELF), ]

# create table for paper
if (FALSE) {
    tmp <- with( subset(euro,YESOK), table( ORIG_COUNTRYSELF, COUNTRY_GFOLX, PRIMARY_LANGUAGE, COUNTRY_SELF ) )
    tmp <- as.data.frame(tmp)
    tmp <- subset(tmp,Freq>0)
    tmp
}

# Write out
write.table(euro, "Euro-samples-info-fine.tsv", sep="\t", row.names=FALSE)

# compute pairwise geographic distances
largest <- subset(largest, largest$COUNTRY_SELF %in% euro$COUNTRY_SELF)
citydists <- rdist.earth(largest[,c("long","lat")],miles=FALSE)
dimnames(citydists) <- list( largest$COUNTRY_SELF, largest$COUNTRY_SELF )
write.table(citydists, file="citypair.dists.tsv", sep="\t")

# This is not efficient to access precomputed.
# # and pairwise PC distances
# pca.pos <- pcas[,c("PC1","PC2")]
# pca.dist <- rdist(pca.pos,pca.pos)
# rownames(pca.dist) <- pcas$labels
# colnames(pca.dist) <- pcas$labels
# write.table(pca.dist, file=paste(.pcadir,"pca.dists.tsv",sep=""), sep="\t")

