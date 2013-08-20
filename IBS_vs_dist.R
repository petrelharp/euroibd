# Written 2013 by Peter Ralph and Graham Coop
# 
# contact: petrel.harp@gmail.com
#
#     To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty. 
# 
# You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 
# 
#


##find rels

source("ibd-blocks-fns.R")

#indivinfo<-read.table("/home/ibd/data/POPRES/european_labels/Euro-samples-info.tsv",head=TRUE,as.is=TRUE)
indivinfo<-getsampleinfo(remove.qc=TRUE) 

if(FALSE){
kin<-read.table(paste(.pcadir,"kin.genome.gz",sep=""),as.is=T,head=T)

kin<-kin[kin$IID1 %in% indivinfo$SUBJID & kin$IID2 %in% indivinfo$SUBJID,]

country1<-rep(NA,nrow(kin))
country2<-rep(NA,nrow(kin))

for(i in 1:nrow(indivinfo)){
	country<-indivinfo$COUNTRY_SELF[i] ##need to watch out for levels
	ind<-indivinfo$SUBJID[i]
	country1[kin$IID1==ind]<-country
	country2[kin$IID2==ind]<-country
	
}

mean.IBS<-tapply(kin$DST,kin[,c("country1","country2")],mean)
num.pairs.IBS<-tapply(kin$DST,kin[,c("country1","country2")],length)
num.pairs.IBS<-apply(num.pairs.IBS,1,function(x){x[is.na(x)]<-0;x})
mean.IBS<-apply(mean.IBS,1,function(x){x[is.na(x)]<-0;x})


tot.IBS<-(mean.IBS*num.pairs.IBS+t(mean.IBS*num.pairs.IBS))/(num.pairs.IBS+t(num.pairs.IBS))
##not symetric
save(mean.IBS,num.pairs.IBS,tot.IBS,file=paste(.pcadir,"mean_IBS.Robj",sep=""))
}

load(file=paste(.pcadir,"mean_IBS.Robj",sep=""))
load("/home/ibd/data/POPRES/ibdblocks/eda-data.Rdata")

countries<-rownames(tot.IBS)
countryB<-character()
countryA<-character()
column.tot.IBS<-numeric()

for(country in countries){
	countryA<-c(countryA,rep(country,length(countries)))
	countryB<-c(countryB,countries)
	column.tot.IBS<-c(column.tot.IBS,tot.IBS[country,])
}


countrypairs.IBS<-cbind(countryA,countryB,column.tot.IBS)

gdists<-apply(countrypairs.IBS,1,function(my.pair){
	this.pair<-(my.pair[1]==poppairs$country1 & my.pair[2]==poppairs$country2) | (my.pair[2]==poppairs$country1 & my.pair[1]==poppairs$country2);
	if(sum(this.pair)){
		return(poppairs$gdist[this.pair])
		}else{
		return(NA)	
		}
	})

countrypairs.IBS<-cbind(countrypairs.IBS,as.numeric(gdists))



newcats <- list( 
        I=c("Italy","Spain","Portugal"),
        W=c("France", "United Kingdom", "Scotland", "England", "Ireland", "Swiss German", "Swiss French", "Switzerland", "Belgium", "Netherlands", "Germany" ),
        N=c( "Sweden", "Norway", "Denmark", "Latvia", "Finland" ),
        E=c( "Slovakia", "Greece", "Yugoslavia", "Albania", "Bosnia", "Montenegro", "Macedonia", "Kosovo", "Serbia", "Bulgaria", "Romania", "Poland", "Hungary", "Czech Republic", "Russia", "Slovenia", "Ukraine", "Croatia", "Austria"),
        TC=c("Turkey","Cyprus")
    )
newcat <- rep(names(newcats),times=sapply(newcats,length)); names(newcat) <- unlist(newcats)
tmp.X <- newcat[countryA]
tmp.Y <- newcat[countryB]
country.pairs <- as.factor( paste( ifelse(tmp.X<tmp.Y,tmp.X,tmp.Y), ifelse(tmp.X<tmp.Y,tmp.Y,tmp.X), sep="-" ) )
levels( country.pairs ) <- c(
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
    )[ levels(country.pairs) ]
# smcat.cols <- c( "E-E"="#66C2A5", "I-(E,N,W)"="#FC8D62", "between E,N,W"="#8DA0CB", "TC-any"="#E78AC3", "N-N"="#A6D854", "W-W"="#FFD92F" )
smcat.cols <- rainbow_hcl(nlevels(country.pairs), c=90); names(smcat.cols) <- levels(country.pairs)
country.cols <- smcat.cols[country.pairs]

nsamples <- table(indivinfo$COUNTRY_SELF)
npairs <- ifelse( countryA==countryB, choose(nsamples[countryA],2), nsamples[countryA]*nsamples[countryB] ) 
country.cex <- pmin(3,pmax(sqrt(npairs)/50,.25))   #  point sizes reflecting sample sizes

pdf(file="IBS_vs_dist.pdf")
plot(gdists,column.tot.IBS,bg=adjustcolor(country.cols,.4),col=adjustcolor(country.cols,.8),cex=country.cex,pch=21,xlab='geographic distance (km)',ylab='mean probability of IBS')
for(comp.col in unique(country.cols)) abline(lm(column.tot.IBS[country.cols==comp.col]~gdists[country.cols==comp.col]),col=comp.col)
dev.off()

require("RSVGTipsDevice")
devSVGTips(file="IBS_vs_dist.svg", width=10, height=7.5, toolTipMode=1, title="Number of blocks shared, by geographic distance")
plot(gdists,column.tot.IBS,col=country.cols,type='n',xlab='geographic distance (km)',ylab='mean probability of IBS')
legend.svg(gdists,column.tot.IBS,labels=paste(countryA,countryB,sep="-"),bg=adjustcolor(country.cols,.4),col=adjustcolor(country.cols,.8),cex=country.cex,pch=21)
dev.off()
# fix bug in RSVGTipsDevice
system("sed -i -e 's/client\\([XY]\\)/page\\1/' IBS_vs_dist.svg")

save(file=paste(.pcadir,"IBS_by_dist.Robj",sep=""),gdists,column.tot.IBS,countryA,countryB,country.cols)
load("IBS_by_dist.Robj")
