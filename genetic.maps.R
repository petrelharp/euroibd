# Written 2013 by Peter Ralph and Graham Coop
# 
# contact: petrel.harp@gmail.com
#
#     To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty. 
# 
# You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>. 
# 
#
female.map<-read.table("/home/ibd/data/genetic_maps/female.gmap.gz",as.is=T,head=T)
male.map<-read.table("/home/ibd/data/genetic_maps/male.gmap.gz",as.is=T,head=T)

# sex.avg<-cbind( female.map[female.map$chr!="chrX",],male.map$cM)
# colnames(sex.avg)<-c("chr","rs","pos","cM.female","cM.male")
# sex.avg$cM.avg<-(sex.avg$cM.female + sex.avg$cM.male)/2

sex.avg <- merge( female.map, male.map, by=c("chr","snp","pos"), all=TRUE, sort=FALSE )
sex.avg$cM.avg <- rowMeans( sex.avg[,c("cM.x","cM.y")], na.rm=TRUE )

sex.avg$chr[ sex.avg$chr=="chrX" ] <- "chr23"

marker.map<-read.table("/home/ibd/data/POPRES/original/phg000027.POPRES.genotype-calls.Affy500K.p1.MULTI.marker-info/POPRES_Snps_QC2.txt.gz",as.is=T)
colnames(marker.map)<-c("snp.id","rs","alleles","strand")

chrs<-1:23
for(chr in chrs){
    chr.markers<-read.table(paste("/home/ibd/data/POPRES/dbgap_chr/POPRES_chr",chr,".map",sep=""),as.is=T)
    colnames(chr.markers)<-c("chr","snp.id","other","pos")

    chr.map<-merge(chr.markers,marker.map,by.x="snp.id",by.y="snp.id",sort=FALSE)

    genetic.chr <- sex.avg[sex.avg$chr ==paste("chr",chr,sep=""),]

    genetic.chr$genetic.map<-rep(NA,nrow(genetic.chr))
    genetic.chr <- genetic.chr[ order(genetic.chr$pos), ]
    genetic.chr$genetic.map[!is.na(genetic.chr$cM.avg)]<-cumsum(genetic.chr$cM.avg[!is.na(genetic.chr$cM.avg)])

    # chr.map$genetic.map<-sapply(chr.map$pos,interpolate.map,genetic.chr)
    chr.map$genetic.map <- with(genetic.chr, approx( x=pos, y=genetic.map, xout=chr.map$pos )$y )
    new.chr<-chr.map[order(chr.map$pos),]
    diffs<-diff(new.chr$genetic,na.rm=T)
    write.table(file=paste("/home/ibd/data/genetic_maps/marker.genetic",chr,".gmap",sep=""), new.chr[,c("chr","snp.id","genetic.map","pos")], quote=FALSE)
    cat(chr,all(diffs[!is.na(diffs)]>=0),max(new.chr$genetic,na.rm=T),"\n")
}


# interpolate.map<-function(pos,genetic.chr){
#     # this does the same thing as approx() but not vectorized
# 
#     where<-genetic.chr$pos <= pos
#     l<-max(which(where))
#     r<-min(which(!where))
# 
#     l.pos<-genetic.chr$pos[l]
#     r.pos<-genetic.chr$pos[r]
#     l.g.pos<-genetic.chr$genetic.map[l]
#     r.g.pos<-genetic.chr$genetic.map[r]
# 
#     if(is.na(l.g.pos) | is.na(r.g.pos)) return(NA)
#     grad<-(r.g.pos-l.g.pos)/(r.pos-l.pos)
# 
#     genetic.pos<- l.g.pos + grad*(pos-l.pos)
#     if(genetic.pos <= r.g.pos & genetic.pos >=l.g.pos) return(genetic.pos)
# 
#     cat("messed up",pos,l.g.pos,genetic.pos,r.g.pos,"\n" )
# }
# 
