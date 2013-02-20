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
source("~/projects/genome/ibd-blocks-fns.R")

load("eda-data-fine.Rdata")

library("colorspace")

#####
pdf(file="Italy_PCA_vs_IBD.pdf")
plist <- c("France", "Italy", "Greece", "Turkey", "Cyprus")
pcols <- countrycols; pcols[plist] <- rainbow_hcl(length(plist),l=60,c=70)
layout(matrix(1:4,nrow=2))
par(mar=c(4,4,1,1))
#xlim=c(150,700), ylim=c(55,250),
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( United.Kingdom.blocks ~ PC1, pch=21,ylim=c(55,250),xlim=c(-0.07,.01), xlab="PC1", ylab="# blocks with UK", cex=1.5,  
            subset=(COUNTRY_SELF=="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), col=NA, )
        points( United.Kingdom.blocks ~ PC1, cex=1.5, subset=(COUNTRY_SELF!="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), 
            pch=ifelse(COUNTRY_SELF%in%c("Switzerland","France"),22,23), col=adjustcolor("black",.25) )
    } )
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( Swiss.French.blocks ~ PC1, pch=21,ylim=c(150,700),xlim=c(-0.07,.01), xlab="PC1", ylab="# blocks with CHf", cex=1.5,  
            subset=(COUNTRY_SELF=="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), col=NA, main="" )
        points( Swiss.French.blocks ~ PC1, cex=1.5, subset=(COUNTRY_SELF!="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), 
            pch=ifelse(COUNTRY_SELF%in%c("Switzerland","France"),22,23), col=adjustcolor("black",.25) )
    } )
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( United.Kingdom.blocks ~ PC2, pch=21,xlim=c(-0.03,0.05),ylim=c(55,250), xlab="PC2", ylab="# blocks with UK", cex=1.5,  
            subset=(COUNTRY_SELF=="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), col=NA, main="" )
        points( United.Kingdom.blocks ~ PC2, cex=1.5, subset=(COUNTRY_SELF!="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), 
            pch=ifelse(COUNTRY_SELF%in%c("Switzerland","France"),22,23), col=adjustcolor("black",.25) )
    } )
legend("topright", legend=countryabbrevs[plist], pch=c(22,21,23,23,23), pt.cex=1.5, bg="white", pt.bg=pcols[plist], col=ifelse(plist=="Italy",NA,"black") )
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( Swiss.French.blocks ~ PC2, pch=21,xlim=c(-0.03,0.05),ylim=c(150,700),, xlab="PC2", ylab="# blocks with CHf", cex=1.5,  
            subset=(COUNTRY_SELF=="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), col=NA, main="" )
        points( Swiss.French.blocks ~ PC2, cex=1.5, subset=(COUNTRY_SELF!="Italy"), bg=adjustcolor(pcols[COUNTRY_SELF],.75), 
            pch=ifelse(COUNTRY_SELF%in%c("Switzerland","France"),22,23), col=adjustcolor("black",.25) )
    } )
dev.off()

####UK
pdf(file="UK_PCA_vs_IBD.pdf")
layout(matrix(1:4,nrow=2))
par(mar=c(4,4,1,1))
plist <- c("Germany","United Kingdom","Ireland")
pcols <- countrycols; pcols[plist] <- rainbow_hcl(length(plist),l=60,c=70)
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( Germany.blocks ~ PC1, pch=21, xlab="PC1", ylab="# blocks with Germany", bg=adjustcolor(pcols[COUNTRY_SELF],.75), cex=1.5, col=NA, 
            subset=(COUNTRY_SELF=="United Kingdom"), main="") 
        points( Germany.blocks ~ PC1, pch=ifelse(COUNTRY_SELF%in%c("Germany","Poland"),22,23), bg=adjustcolor(pcols[COUNTRY_SELF],.75), 
            col=adjustcolor("black",.25), cex=1.5, subset=(COUNTRY_SELF!="United Kingdom") ) 
    } )
legend("topleft",pch=c(22,21,23),pt.bg=pcols[plist],legend=countryabbrevs[plist],pt.cex=1.5,col=ifelse(plist=="United Kingdom",NA,"black") )
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( Ireland.blocks ~ PC1, pch=21, xlab="PC1", ylab="# blocks with Ireland", bg=adjustcolor(pcols[COUNTRY_SELF],.75), cex=1.5, col=NA, 
            subset=(COUNTRY_SELF=="United Kingdom"), main="") 
        points( Ireland.blocks ~ PC1, pch=ifelse(COUNTRY_SELF%in%c("Germany","Poland"),22,23), bg=adjustcolor(pcols[COUNTRY_SELF],.75), 
            col=adjustcolor("black",.25), cex=1.5, subset=(COUNTRY_SELF!="United Kingdom") ) 
    } )
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( Germany.blocks ~ PC2, pch=21, xlab="PC2", ylab="# blocks with Germany", bg=adjustcolor(pcols[COUNTRY_SELF],.75), cex=1.5, col=NA, 
            subset=(COUNTRY_SELF=="United Kingdom"), main="") 
        points( Germany.blocks ~ PC2, pch=ifelse(COUNTRY_SELF%in%c("Germany","Poland"),22,23), bg=adjustcolor(pcols[COUNTRY_SELF],.75), 
            col=adjustcolor("black",.25), cex=1.5, subset=(COUNTRY_SELF!="United Kingdom") ) 
    } )
with( subset(indivinfo,COUNTRY_SELF%in%plist), {
        plot( Ireland.blocks ~ PC2, pch=21, xlab="PC2", ylab="# blocks with Ireland", bg=adjustcolor(pcols[COUNTRY_SELF],.75), cex=1.5, col=NA, 
            subset=(COUNTRY_SELF=="United Kingdom"), main="") 
        points( Ireland.blocks ~ PC2, pch=ifelse(COUNTRY_SELF%in%c("Germany","Poland"),22,23), bg=adjustcolor(pcols[COUNTRY_SELF],.75), 
            col=adjustcolor("black",.25), cex=1.5, subset=(COUNTRY_SELF!="United Kingdom") ) 
    } )
dev.off()



country.symbol <- rep(c(21,22,23,24,25),8)
names(country.symbol) <- levels(indivinfo$COUNTRY_SELF)
nsamples <- table(indivinfo$COUNTRY_SELF)

pdf("pca-map.pdf",width=6.5,height=6.5,pointsize=10)
par(mar=c(5,4,1,1)+.1)
pcols <- countrycols; 
with(indivinfo,plot(PC1,PC2,bg=adjustcolor(pcols[COUNTRY_SELF],0.5),pch=country.symbol[COUNTRY_SELF],col=NA,xlim=c(-0.07,0.05)))
country.means <- do.call(rbind, with(indivinfo, tapply( 1:nrow(indivinfo), COUNTRY_SELF, function (k) { colMeans( indivinfo[k,c("PC1","PC2")] ) } ) ) )
text( country.means, labels=countryabbrevs )
points( country.means, pch=21, col=NA, bg=adjustcolor(countrycols,.25), cex=4 )
legend( "topright",pt.bg=countrycols,legend=countryabbrevs,pch=country.symbol,col=NA,cex=8/12,pt.cex=0.9)
dev.off()
